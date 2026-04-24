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

        public static TrueConvexHull4D Compute(List<double[]> points, double eps = 1e-7)
        {
            var hull = new TrueConvexHull4D { Vertices = points };
            if (points.Count < 5) return hull;

            var discoveredCells = new List<List<int>>();
            var firstNormal = FindInitialNormal(points, eps);
            var firstCell = GetSupportFace(points, firstNormal, eps);
            discoveredCells.Add(firstCell);

            var ridgeQueue = new Queue<(List<int> ridge, double[] currentNormal)>();
            foreach (var f in FindMaximalFacets(points, firstCell, 3, eps)) ridgeQueue.Enqueue((f, firstNormal));

            var cellKeys = new HashSet<string> { GetKey(firstCell) };
            var ridgeKeys = new HashSet<string>(ridgeQueue.Select(r => GetKey(r.ridge)));
            var discoveredFaces = ridgeQueue.Select(r => r.ridge).ToList();

            while (ridgeQueue.Count > 0)
            {
                var (ridge, prevNormal) = ridgeQueue.Dequeue();
                var nextNormal = Pivot(points, ridge, prevNormal, eps);
                if (nextNormal == null) continue;
                var nextCell = GetSupportFace(points, nextNormal, eps);
                if (cellKeys.Add(GetKey(nextCell)))
                {
                    discoveredCells.Add(nextCell);
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
            return facets;
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

        static double[] Pivot(List<double[]> pts, List<int> ridge, double[] prevN, double eps)
        {
            var v0 = pts[ridge[0]];
            var e1 = Sub(pts[ridge[1]], v0);
            var e2 = Sub(pts[ridge[2]], v0);
            double minA = double.MaxValue;
            double[] bestN = null;
            for (int i = 0; i < pts.Count; i++)
            {
                if (ridge.Contains(i)) continue;
                var n = Cross4D(e1, e2, Sub(pts[i], v0));
                if (Mag(n) < 1e-10) continue;
                if (Dot(n, prevN) > 1.0 - 1e-9) continue; 
                double d = Dot(n, pts[i]);
                if (IsExtreme(pts, n, d, eps)) {
                    double a = Math.Acos(Math.Max(-1, Math.Min(1, Dot(n, prevN))));
                    if (a < minA) { minA = a; bestN = n; }
                }
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

        static double[] FindInitialNormal(List<double[]> pts, double eps)
        {
            int p0 = 0; for (int i = 1; i < pts.Count; i++) if (pts[i][0] < pts[p0][0]) p0 = i;
            for (int i = 0; i < pts.Count; i++)
            for (int j = i+1; j < pts.Count; j++)
            for (int k = j+1; k < pts.Count; k++)
            {
                var n = Cross4D(Sub(pts[i], pts[p0]), Sub(pts[j], pts[p0]), Sub(pts[k], pts[p0]));
                if (Mag(n) < 1e-10) continue;
                if (IsExtreme(pts, n, Dot(pts[p0], n), eps)) return n;
                if (IsExtreme(pts, Scale(n, -1), Dot(pts[p0], Scale(n, -1)), eps)) return Scale(n, -1);
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
        static double[] Cross4D(double[] u, double[] v, double[] w) {
            var r = new double[4];
            r[0] =  (u[1]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[1]*w[3]-v[3]*w[1]) + u[3]*(v[1]*w[2]-v[2]*w[1]));
            r[1] = -(u[0]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[2]-v[2]*w[0]));
            r[2] =  (u[0]*(v[1]*w[3]-v[3]*w[1]) - u[1]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[1]-v[1]*w[0]));
            r[3] = -(u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0]));
            return Normalize(r);
        }
        static double[][] GetOrthoBasis(List<double[]> pts, List<int> verts, int dim) {
            var b = new double[dim][]; var v0 = pts[verts[0]]; int f = 0;
            for (int i = 1; i < verts.Count && f < dim; i++) {
                var v = Sub(pts[verts[i]], v0);
                for (int j = 0; j < f; j++) v = Sub(v, Scale(b[j], Dot(v, b[j])));
                if (Mag(v) > 1e-9) b[f++] = Normalize(v);
            }
            return b;
        }
        static double[] Project(double[] p, double[][] b) { var r = new double[b.Length]; for (int i = 0; i < b.Length; i++) r[i] = Dot(p, b[i]); return r; }
        static string GetKey(IEnumerable<int> idx) => string.Join(",", idx.OrderBy(x => x));
    }
}
