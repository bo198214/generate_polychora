using System;
using System.Collections.Generic;
using System.Linq;

namespace D4BB.Geometry
{
    public class TrueConvexHull4D
    {
        public List<double[]> Vertices;
        public List<int[]> Edges = new List<int[]>();
        public List<int[]> Faces = new List<int[]>();
        public List<int[]> Cells = new List<int[]>();
        /// <summary>Outward unit normals, one per cell (parallel to Cells list).</summary>
        public List<double[]> Normals = new List<double[]>();

        public static TrueConvexHull4D Compute(List<double[]> points, double eps = 1e-7,
                                               Action<int> onCellDiscovered = null,
                                               Action<int, int> onRidgeProcessed = null)
        {
            var hull = new TrueConvexHull4D { Vertices = points };
            if (points.Count < 5) return hull;

            // Centroid computed once — enables O(1) normal orientation in Pivot.
            var centroid = new double[4];
            foreach (var p in points) for (int c = 0; c < 4; c++) centroid[c] += p[c] / points.Count;

            var discoveredCells = new List<List<int>>();
            var discoveredNormals = new List<double[]>();
            var firstNormal = FindInitialNormal(points, centroid, eps);
            var firstCell = GetSupportFace(points, firstNormal, eps);

            // If the initial supporting hyperplane only touches 1-3 vertices it is
            // too degenerate to seed the gift-wrapping.  Rebuild a proper 4-vertex
            // simplex around the extreme vertex via Gram-Schmidt and re-derive the normal.
            // If the initial cell is too small to seed the algorithm, search for a
            // proper face normal via all triplets among the K nearest neighbours of
            // the extreme vertex.  O(K^3 · V) — fast for any V with K=30.
            if (firstCell.Count < 4)
            {
                int p0idx = firstCell.Count > 0 ? firstCell[0] : 0;
                int K = Math.Min(50, points.Count - 1);
                var nearby = Enumerable.Range(0, points.Count)
                    .Where(i => i != p0idx)
                    .OrderBy(i => { double s=0; for(int k=0;k<4;k++){double dx=points[i][k]-points[p0idx][k]; s+=dx*dx;} return s; })
                    .Take(K).ToList();

                // Use a tighter eps so near-degenerate hyperplanes don't pass IsExtreme.
                double searchEps = Math.Max(eps / 100.0, 1e-9);
                bool found = false;
                var p0dir = Normalize(Sub(points[p0idx], centroid)); // for potential future sorting
                for (int a = 0; a < nearby.Count && !found; a++)
                for (int b = a+1; b < nearby.Count && !found; b++)
                for (int c = b+1; c < nearby.Count && !found; c++)
                {
                    var raw = Cross4DRaw(Sub(points[nearby[a]], points[p0idx]),
                                        Sub(points[nearby[b]], points[p0idx]),
                                        Sub(points[nearby[c]], points[p0idx]));
                    if (Mag(raw) < 1e-10) continue;
                    var n = Normalize(raw);
                    double d = Dot(points[p0idx], n);
                    if (Dot(centroid, n) > d) { n = Scale(n, -1); d = -d; }
                    if (!IsExtreme(points, n, d, searchEps)) continue;
                    // Use GetSupportFace with eps (normal eps) so the cell is consistent
                    // with the rest of the algorithm which also uses eps.
                    var cand = GetSupportFace(points, n, eps);
                    if (cand.Count >= 4) { firstNormal = n; firstCell = cand; found = true; }
                }
            }

            discoveredCells.Add(firstCell);
            discoveredNormals.Add(firstNormal);
            onCellDiscovered?.Invoke(1);   // signal: initialisation done, ridge loop starting

            var ridgeQueue = new Queue<(List<int> ridge, double[] currentNormal)>();
            foreach (var f in FindMaximalFacets(points, firstCell, 3, eps)) ridgeQueue.Enqueue((f, firstNormal));

            var cellKeys = new HashSet<string> { GetKey(firstCell) };
            var ridgeKeys = new HashSet<string>(ridgeQueue.Select(r => GetKey(r.ridge)));
            var discoveredFaces = ridgeQueue.Select(r => r.ridge).ToList();

            int ridgesDone = 0;
            while (ridgeQueue.Count > 0)
            {
                var (ridge, prevNormal) = ridgeQueue.Dequeue();
                onRidgeProcessed?.Invoke(++ridgesDone, ridgeQueue.Count);
                var nextNormal = Pivot(points, ridge, prevNormal, centroid, eps);
                if (nextNormal == null) continue;
                var nextCell = GetSupportFace(points, nextNormal, eps);
                if (cellKeys.Add(GetKey(nextCell)))
                {
                    discoveredCells.Add(nextCell);
                    discoveredNormals.Add(nextNormal);
                    onCellDiscovered?.Invoke(discoveredCells.Count);
                    foreach (var subFacet in FindMaximalFacets(points, nextCell, 3, eps))
                        if (ridgeKeys.Add(GetKey(subFacet))) { discoveredFaces.Add(subFacet); ridgeQueue.Enqueue((subFacet, nextNormal)); }
                }
            }

            var edgeSet = new HashSet<string>();
            foreach (var face in discoveredFaces)
                foreach (var edge in FindMaximalFacets(points, face, 2, eps))
                    if (edgeSet.Add(GetKey(edge))) hull.Edges.Add(edge.ToArray());

            hull.Cells = discoveredCells.Select(c => c.ToArray()).ToList();
            hull.Faces = discoveredFaces.Select(f => f.ToArray()).ToList();
            hull.Normals = discoveredNormals;
            return hull;
        }

        static List<List<int>> FindMaximalFacets(List<double[]> allPts, List<int> verts, int dim, double eps)
        {
            if (dim == 2) return FindPolygonEdges(allPts, verts, eps);
            var basis = GetOrthoBasis(allPts, verts, dim);
            var v0 = allPts[verts[0]];
            var proj = verts.Select(v => Project(Sub(allPts[v], v0), basis)).ToList();
            var facets = new List<List<int>>();
            var seen = new HashSet<string>();
            for (int i = 0; i < proj.Count; i++)
            for (int j = i + 1; j < proj.Count; j++)
            for (int k = j + 1; k < proj.Count; k++)
            {
                var n = Normalize(Cross3D(Sub(proj[j], proj[i]), Sub(proj[k], proj[i])));
                if (Mag(n) < 1e-10) continue;
                double d = Dot(proj[i], n);
                if (IsExtreme(proj, n, d, eps)) {
                    var f = new List<int>();
                    for (int l = 0; l < proj.Count; l++) if (Math.Abs(Dot(proj[l], n) - d) < eps) f.Add(verts[l]);
                    if (seen.Add(GetKey(f))) facets.Add(f);
                }
            }
            // Keep only maximal facets — discard any whose vertex set is a proper subset of a larger facet
            var facetSets = facets.Select(f => new HashSet<int>(f)).ToList();
            var maximal = new List<List<int>>();
            for (int i = 0; i < facets.Count; i++)
            {
                bool isSubset = false;
                for (int j = 0; j < facets.Count; j++)
                {
                    if (facets[j].Count <= facets[i].Count) continue;
                    if (facets[i].All(v => facetSets[j].Contains(v))) { isSubset = true; break; }
                }
                if (!isSubset) maximal.Add(facets[i]);
            }
            return maximal;
        }

        static List<List<int>> FindPolygonEdges(List<double[]> allPts, List<int> verts, double eps)
        {
            var basis = GetOrthoBasis(allPts, verts, 2);
            var v0 = allPts[verts[0]];
            var proj = verts.Select(v => Project(Sub(allPts[v], v0), basis)).ToList();
            var edges = new List<List<int>>();
            for (int i = 0; i < proj.Count; i++)
            for (int j = i + 1; j < proj.Count; j++)
            {
                var vec = Sub(proj[j], proj[i]);
                var n = Normalize(new[] { -vec[1], vec[0] });
                double d = Dot(proj[i], n);
                if (IsExtreme(proj, n, d, eps)) edges.Add(new List<int> { verts[i], verts[j] });
            }
            return edges;
        }

        static double[] Pivot(List<double[]> pts, List<int> ridge, double[] prevN, double[] centroid, double eps)
        {
            // For large polytopes the minimum dihedral angle scales as ~1/√V.
            // Shrink the "same-cell" guard proportionally so tiny-angle neighbours
            // are not mistaken for the current cell.
            double thetaEps = Math.Max(1e-14, 1e-9 / Math.Sqrt(pts.Count));

            var v0 = pts[ridge[0]];
            var e1 = Sub(pts[ridge[1]], v0);
            var e2 = Sub(pts[ridge[2]], v0);

            // Build a second basis vector for the 2D rotation plane (complement of ridge span).
            // p2 = unit vector perpendicular to e1, e2, prevN — gives signed rotation angle.
            var p2raw = Cross4DRaw(e1, e2, prevN);
            var p2 = Mag(p2raw) > 1e-9 ? Normalize(p2raw) : null;

            var ridgeSet = new HashSet<int>(ridge);
            double minTheta = double.MaxValue;
            double[] bestN = null;
            for (int i = 0; i < pts.Count; i++)
            {
                if (ridgeSet.Contains(i)) continue;
                var raw = Cross4DRaw(e1, e2, Sub(pts[i], v0));
                if (Mag(raw) < 1e-10) continue;
                var n = Normalize(raw);
                double d = Dot(n, v0);
                // All ridge vertices must lie on this hyperplane (required for valid pivot)
                bool ridgeOnPlane = true;
                for (int r = 3; r < ridge.Count; r++) {
                    if (Math.Abs(Dot(n, pts[ridge[r]]) - d) > eps) { ridgeOnPlane = false; break; }
                }
                if (!ridgeOnPlane) continue;
                // O(1) orientation: centroid must lie on the negative side.
                if (Dot(n, centroid) - d > 0) n = Scale(n, -1);
                // Signed rotation angle from prevN to n in the (prevN, p2) plane.
                // Maps to (0, 2π]: skip zero (same cell) and pick smallest positive angle.
                double cosT = Dot(n, prevN);
                double sinT = p2 != null ? -Dot(n, p2) : 0;  // negative = correct CCW direction
                double theta = Math.Atan2(sinT, cosT);
                if (theta <= thetaEps) theta += 2 * Math.PI;
                if (theta > 2 * Math.PI - thetaEps) continue;
                if (theta < minTheta) { minTheta = theta; bestN = n; }
            }
            // Verify the best candidate is actually a supporting hyperplane (numerical safety net).
            if (bestN != null) {
                double d = Dot(bestN, v0);
                if (!IsExtreme(pts, bestN, d, eps)) bestN = null;
            }
            return bestN;
        }

        static bool IsExtreme(List<double[]> pts, double[] n, double d, double eps)
        {
            int side = 0;
            foreach (var p in pts) {
                double v = Dot(n, p) - d; if (Math.Abs(v) < eps) continue;
                if (side == 0) side = v > 0 ? 1 : -1; else if ((v > 0 && side == -1) || (v < 0 && side == 1)) return false;
            }
            return true;
        }

        static double[] FindInitialNormal(List<double[]> pts, double[] centroid, double eps)
        {
            int p0 = 0; for (int i = 1; i < pts.Count; i++) if (pts[i][0] < pts[p0][0]) p0 = i;

            // Greedy simplex (O(V) per step): the original approach, correct for all sizes.
            var basis = new List<double[]>();
            var chosen = new List<int> { p0 };
            for (int step = 0; step < 3; step++) {
                double maxDist = 0; int best = -1;
                for (int i = 0; i < pts.Count; i++) {
                    if (chosen.Contains(i)) continue;
                    var v = Sub(pts[i], pts[p0]);
                    foreach (var b in basis) v = Sub(v, Scale(b, Dot(v, b)));
                    double d = Mag(v);
                    if (d > maxDist) { maxDist = d; best = i; }
                }
                if (best < 0 || maxDist < 1e-9) break;
                var nb = Sub(pts[best], pts[p0]);
                foreach (var b in basis) nb = Sub(nb, Scale(b, Dot(nb, b)));
                basis.Add(Normalize(nb)); chosen.Add(best);
            }
            if (chosen.Count == 4) {
                var raw = Cross4DRaw(Sub(pts[chosen[1]], pts[p0]),
                                    Sub(pts[chosen[2]], pts[p0]),
                                    Sub(pts[chosen[3]], pts[p0]));
                if (Mag(raw) > 1e-10) {
                    var n = Normalize(raw);
                    double d = Dot(pts[p0], n);
                    if (Dot(centroid, n) > d) { n = Scale(n, -1); d = -d; }
                    if (IsExtreme(pts, n, d, eps)) return n;
                }
            }

            // Fallback for large near-spherical polytopes where the greedy simplex
            // hyperplane cuts through the interior: use normalize(p0−centroid) which
            // is always outward and typically extreme for symmetric polytopes.  O(V).
            for (int i = 0; i < pts.Count; i++) {
                double dist = 0;
                for (int k = 0; k < 4; k++) { double dx = pts[i][k]-centroid[k]; dist += dx*dx; }
                if (dist > 0) { double d0 = 0; for (int k=0;k<4;k++) d0+=(pts[p0][k]-centroid[k])*(pts[p0][k]-centroid[k]); if (dist > d0) p0 = i; }
            }
            { var n = Normalize(Sub(pts[p0], centroid));
              double d = Dot(pts[p0], n);
              if (IsExtreme(pts, n, d, eps)) return n; }

            // Last resort: axis-aligned directions (always valid supporting hyperplanes).
            for (int ax = 0; ax < 4; ax++) for (int sg = -1; sg <= 1; sg += 2) {
                var na = new double[4]; na[ax] = sg;
                double da = pts.Max(p => Dot(p, na));
                if (IsExtreme(pts, na, da, eps)) return na;
            }
            return new[] { 1.0, 0, 0, 0 };
        }

        static List<int> GetSupportFace(List<double[]> pts, double[] n, double eps)
        {
            double max = pts.Max(p => Dot(p, n));
            var l = new List<int>();
            for (int i = 0; i < pts.Count; i++) if (Math.Abs(Dot(pts[i], n) - max) < eps) l.Add(i);
            return l;
        }

        static double Dot(double[] a, double[] b) { double s = 0; for (int i = 0; i < a.Length; i++) s += a[i]*b[i]; return s; }
        static double[] Sub(double[] a, double[] b) { var r = new double[a.Length]; for (int i = 0; i < a.Length; i++) r[i] = a[i]-b[i]; return r; }
        static double Mag(double[] a) => Math.Sqrt(Dot(a, a));
        static double[] Normalize(double[] a) { double m = Mag(a); return m < 1e-15 ? a : Scale(a, 1/m); }
        static double[] Scale(double[] a, double s) { var r = new double[a.Length]; for (int i = 0; i < a.Length; i++) r[i] = a[i]*s; return r; }
        static double[] Cross3D(double[] u, double[] v) => new[] { u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0] };
        static double[] Cross4DRaw(double[] u, double[] v, double[] w) {
            var r = new double[4];
            r[0] =  (u[1]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[1]*w[3]-v[3]*w[1]) + u[3]*(v[1]*w[2]-v[2]*w[1]));
            r[1] = -(u[0]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[2]-v[2]*w[0]));
            r[2] =  (u[0]*(v[1]*w[3]-v[3]*w[1]) - u[1]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[1]-v[1]*w[0]));
            r[3] = -(u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0]));
            return r;
        }
        static double[][] GetOrthoBasis(List<double[]> pts, List<int> verts, int dim) {
            var b = new double[dim][]; var v0 = pts[verts[0]]; int f = 0;
            for (int i = 1; i < verts.Count && f < dim; i++) {
                var v = Sub(pts[verts[i]], v0);
                for (int j = 0; j < f; j++) v = Sub(v, Scale(b[j], Dot(v, b[j])));
                if (Mag(v) > 1e-9) b[f++] = Normalize(v);
            }
            // Fill remaining slots with axis vectors so b never contains null.
            // Degenerate cells (dim < actual span) will project all points onto the
            // same subplane → IsExtreme will find no faces, which is correct behaviour.
            for (; f < dim; f++)
                for (int ax = 0; ax < 4; ax++) {
                    var c = new double[4]; c[ax] = 1.0;
                    for (int j = 0; j < f; j++) c = Sub(c, Scale(b[j], Dot(c, b[j])));
                    if (Mag(c) > 1e-9) { b[f] = Normalize(c); break; }
                }
            return b;
        }
        static double[] Project(double[] p, double[][] b) { var r = new double[b.Length]; for (int i = 0; i < b.Length; i++) r[i] = Dot(p, b[i]); return r; }
        static string GetKey(IEnumerable<int> idx) => string.Join(",", idx.OrderBy(x => x));
    }
}
