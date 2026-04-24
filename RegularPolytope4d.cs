using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MIConvexHull;

namespace D4BB.Geometry
{
    public class RegularPolytope4d
    {
        public List<double[]> vertices;
        public List<(int a, int b)> edges;
        // Each cell: list of triangle index-triples [v0,v1,v2] (triangulated)
        public List<List<int[]>> cells;
        // Outward 4D unit normal per cell
        public List<double[]> cellNormals;
        public double edgeLength;

        // -----------------------------------------------------------------------
        // Build everything from a vertex list.
        // -----------------------------------------------------------------------
        public static RegularPolytope4d FromVertices(List<double[]> verts, double hullTolerance = 1e-9)
        {
            var p = new RegularPolytope4d { vertices = verts };
            p.edgeLength = ComputeEdgeLength(verts);
            p.edges      = ComputeEdges(verts, p.edgeLength);
            (p.cells, p.cellNormals) = ComputeCells(verts, hullTolerance);
            return p;
        }

        // -----------------------------------------------------------------------
        // Edges: all pairs at minimum pairwise distance
        // -----------------------------------------------------------------------
        static double ComputeEdgeLength(List<double[]> verts)
        {
            double min = double.MaxValue;
            for (int i = 0; i < verts.Count; i++)
            for (int j = i + 1; j < verts.Count; j++)
            {
                double d = Dist(verts[i], verts[j]);
                if (d < min - 1e-6) min = d;
            }
            return min;
        }

        static List<(int, int)> ComputeEdges(List<double[]> verts, double edgeLen)
        {
            var edges = new List<(int, int)>();
            for (int i = 0; i < verts.Count; i++)
            for (int j = i + 1; j < verts.Count; j++)
                if (Math.Abs(Dist(verts[i], verts[j]) - edgeLen) < 1e-4)
                    edges.Add((i, j));
            return edges;
        }

        // -----------------------------------------------------------------------
        // Cells: use MIConvexHull in 4D.
        // Each returned ConvexFace has 4 vertices (tetrahedron) and a Normal.
        // For non-simplicial cells (cubes, octahedra, dodecahedra) the hull
        // triangulates them; we group coplanar triangles into one cell.
        // -----------------------------------------------------------------------
        static (List<List<int[]>>, List<double[]>) ComputeCells(List<double[]> verts, double tolerance = 1e-9)
        {
            var hullVerts = verts.Select((v, i) => new HullVertex(v, i)).ToList();
            var hull = ConvexHull.Create<HullVertex, HullFace>(hullVerts, tolerance);
            if (hull?.Faces == null) throw new Exception("ConvexHull failed");

            // Deduplicate hull faces: MIConvexHull sometimes returns the same tetrahedron
            // multiple times (ridge faces in nearly-degenerate configurations like the 120-cell).
            // Keep only one instance per unique sorted vertex set.
            var allRawFaces = hull.Faces.ToArray();
            {
                var seen = new HashSet<string>();
                var unique = new List<HullFace>();
                foreach (var f in allRawFaces)
                {
                    var key = string.Join(",", f.Vertices.Select(v => ((HullVertex)v).Index).OrderBy(x => x));
                    if (seen.Add(key)) unique.Add(f);
                }
                allRawFaces = unique.ToArray();
            }
            var rawFaces = allRawFaces;
            int nFaces = rawFaces.Length;

            // Precompute (outward normal, hyperplane offset d = n·v0) per face.
            var ns = new double[nFaces][];
            var ds = new double[nFaces];
            for (int i = 0; i < nFaces; i++)
            {
                ns[i] = rawFaces[i].Normal;
                double d = 0;
                var v0 = rawFaces[i].Vertices[0].Position;
                for (int k = 0; k < 4; k++) d += ns[i][k] * v0[k];
                ds[i] = d;
            }

            // Union-Find: merge faces that share the same supporting hyperplane.
            // Two faces are coplanar iff their normals are parallel (dot ≈ 1) and
            // their offsets match — avoids brittle string-rounding of near-unit vectors.
            var parent = Enumerable.Range(0, nFaces).ToArray();
            for (int i = 0; i < nFaces; i++)
            for (int j = i + 1; j < nFaces; j++)
            {
                double dot = 0;
                for (int k = 0; k < 4; k++) dot += ns[i][k] * ns[j][k];
                if (dot > 0.99)
                    UFUnion(parent, i, j);
            }

            // Rescue pass: iteratively merge small groups into large ones by vertex-inclusion.
            // Uses consensus (averaged) normal per large group — robust vs. sliver outliers.
            // A vertex is "in" a hyperplane (n, d) iff |n·v - d| < vertexEps.
            const double vertexEps = 0.15;
            bool anyMerge = true;
            while (anyMerge)
            {
                anyMerge = false;

                // Recount groups.
                var rootCount = new Dictionary<int, int>();
                for (int i = 0; i < nFaces; i++)
                {
                    int r = UFFind(parent, i);
                    rootCount[r] = rootCount.TryGetValue(r, out int c) ? c + 1 : 1;
                }
                if (rootCount.Count == 0) break;

                // "Rescue threshold": any group with >= rescueMin faces is a valid host.
                // Using median would exclude host groups that are themselves below-average
                // due to a split (120-cell: split host has size 23 < median 25).
                // A low absolute floor (4) is safe: the vertex-inclusion test rejects wrong hosts.
                int totalFaces = rootCount.Values.Sum();
                int rescueMin  = Math.Max(4, totalFaces / rootCount.Count / 3);

                // Large groups: potential merge targets.
                var largeGroots = rootCount.Where(kv => kv.Value >= rescueMin).Select(kv => kv.Key).ToArray();
                if (largeGroots.Length == 0) break;

                // Build CONSENSUS (n, d) for each large group (average of all face normals).
                var groupNormals = new Dictionary<int, double[]>();
                var groupOffsets = new Dictionary<int, double>();
                foreach (int gr in largeGroots)
                {
                    double[] avg = new double[4];
                    int cnt = 0;
                    for (int i = 0; i < nFaces; i++)
                        if (UFFind(parent, i) == gr)
                        { for (int k = 0; k < 4; k++) avg[k] += ns[i][k]; cnt++; }
                    double len = Math.Sqrt(avg.Sum(x => x * x));
                    if (len > 1e-12) for (int k = 0; k < 4; k++) avg[k] /= len;
                    groupNormals[gr] = avg;
                    double dSum = 0;
                    for (int i = 0; i < nFaces; i++)
                        if (UFFind(parent, i) == gr)
                        {
                            var v0 = rawFaces[i].Vertices[0].Position;
                            for (int k = 0; k < 4; k++) dSum += avg[k] * v0[k];
                        }
                    groupOffsets[gr] = dSum / cnt;
                }

                // Try to absorb every small group into the closest large group.
                for (int i = 0; i < nFaces; i++)
                {
                    int ri = UFFind(parent, i);
                    if (rootCount.TryGetValue(ri, out int ric) && ric >= rescueMin) continue;

                    var verts4 = rawFaces[i].Vertices;
                    int bestGroup = -1;
                    double bestMaxErr = double.MaxValue;

                    foreach (int gr in largeGroots)
                    {
                        if (UFFind(parent, i) == gr) continue; // already same group
                        double[] n = groupNormals[gr];
                        double d  = groupOffsets[gr];
                        double maxErr = 0;
                        foreach (var v in verts4)
                        {
                            double dot2 = 0;
                            for (int k = 0; k < 4; k++) dot2 += n[k] * v.Position[k];
                            maxErr = Math.Max(maxErr, Math.Abs(dot2 - d));
                        }
                        if (maxErr < vertexEps && maxErr < bestMaxErr)
                        { bestMaxErr = maxErr; bestGroup = gr; }
                    }

                    if (bestGroup >= 0 && UFFind(parent, i) != UFFind(parent, bestGroup))
                    {
                        UFUnion(parent, i, bestGroup);
                        anyMerge = true;
                    }
                }
            }

            // Second pass: merge groups with similar CONSENSUS normals.
            // Consensus normal = average of all face normals in the group (robust vs. sliver outliers).
            {
                var groupFaces = new Dictionary<int, List<int>>();
                for (int i = 0; i < nFaces; i++)
                {
                    int r = UFFind(parent, i);
                    if (!groupFaces.ContainsKey(r)) groupFaces[r] = new List<int>();
                    groupFaces[r].Add(i);
                }

                var roots = groupFaces.Keys.ToArray();
                var consensusN = new Dictionary<int, double[]>();
                var consensusD = new Dictionary<int, double>();
                foreach (int r in roots)
                {
                    double[] avg = new double[4];
                    foreach (int fi in groupFaces[r])
                        for (int k = 0; k < 4; k++) avg[k] += ns[fi][k];
                    double len = Math.Sqrt(avg.Sum(x => x * x));
                    if (len > 1e-12) for (int k = 0; k < 4; k++) avg[k] /= len;
                    consensusN[r] = avg;
                    // Offset: average of n·v0 over all faces using the consensus normal.
                    double dSum = 0;
                    foreach (int fi in groupFaces[r])
                    {
                        var v0 = rawFaces[fi].Vertices[0].Position;
                        for (int k = 0; k < 4; k++) dSum += avg[k] * v0[k];
                    }
                    consensusD[r] = dSum / groupFaces[r].Count;
                }

                // Merge groups whose consensus normals are nearly parallel with same offset.
                for (int i = 0; i < roots.Length; i++)
                for (int j = i + 1; j < roots.Length; j++)
                {
                    int ri = UFFind(parent, roots[i]), rj = UFFind(parent, roots[j]);
                    if (ri == rj) continue;
                    // Use consensus normals of the current roots (might have changed via prior merges).
                    double[] ni = consensusN[ri], nj = consensusN[rj];
                    double dot2 = 0;
                    for (int k = 0; k < 4; k++) dot2 += ni[k] * nj[k];
                    if (dot2 > 0.999 && Math.Abs(consensusD[ri] - consensusD[rj]) < 0.1)
                    {
                        UFUnion(parent, roots[i], roots[j]);
                        // Update consensus for new root.
                        int newRoot = UFFind(parent, roots[i]);
                        double[] merged = new double[4];
                        for (int k = 0; k < 4; k++) merged[k] = (ni[k] + nj[k]) / 2;
                        double mlen = Math.Sqrt(merged.Sum(x => x * x));
                        if (mlen > 1e-12) for (int k = 0; k < 4; k++) merged[k] /= mlen;
                        consensusN[newRoot] = merged;
                        consensusD[newRoot] = (consensusD[ri] + consensusD[rj]) / 2;
                    }
                }
            }

            var groups = new Dictionary<int, (List<int[]> tris, double[] norm, double apothem)>();
            for (int i = 0; i < nFaces; i++)
            {
                int root = UFFind(parent, i);
                if (!groups.ContainsKey(root))
                    groups[root] = (new List<int[]>(), ns[i], ds[i]);
                groups[root].tris.Add(rawFaces[i].Vertices.Select(v => ((HullVertex)v).Index).ToArray());
            }

            // Apothem filter: for regular polytopes all cells share the same apothem.
            // Spurious hull faces (MIConvexHull numerical artifacts) have a different apothem.
            // Discard groups whose apothem deviates by > 3 % from the median.
            if (groups.Count > 1)
            {
                var apothems = groups.Values.Select(g => g.apothem).OrderBy(x => x).ToArray();
                double medianApothem = apothems[apothems.Length / 2];
                double tol = Math.Abs(medianApothem) * 0.03 + 1e-6;
                var validRoots = groups.Keys
                    .Where(r => Math.Abs(groups[r].apothem - medianApothem) <= tol)
                    .ToList();
                foreach (var k in groups.Keys.Except(validRoots).ToList())
                    groups.Remove(k);
            }

            return (groups.Values.Select(g => g.tris).ToList(),
                    groups.Values.Select(g => g.norm).ToList());
        }

        static int UFFind(int[] parent, int i) =>
            parent[i] == i ? i : (parent[i] = UFFind(parent, parent[i]));

        static void UFUnion(int[] parent, int i, int j)
        {
            parent[UFFind(parent, i)] = UFFind(parent, j);
        }

        // -----------------------------------------------------------------------
        // Back-face culling: returns indices of visible cells for a given eye pos
        // -----------------------------------------------------------------------
        public List<int> VisibleCells(double[] eye)
        {
            var visible = new List<int>();
            // polytope center (all regular polytopes used here are centered at origin)
            double[] center = new double[4];
            for (int i = 0; i < cells.Count; i++)
            {
                // cell center = average of all its vertex positions
                double[] cellCenter = CellCenter(i);
                double dot = 0;
                for (int k = 0; k < 4; k++) dot += cellNormals[i][k] * (eye[k] - cellCenter[k]);
                if (dot > 0) visible.Add(i);
            }
            return visible;
        }

        double[] CellCenter(int cellIdx)
        {
            var tris = cells[cellIdx];
            var used = new HashSet<int>();
            foreach (var tri in tris) foreach (var v in tri) used.Add(v);
            double[] c = new double[4];
            foreach (int vi in used) for (int k = 0; k < 4; k++) c[k] += vertices[vi][k];
            for (int k = 0; k < 4; k++) c[k] /= used.Count;
            return c;
        }

        // -----------------------------------------------------------------------
        // Helpers
        // -----------------------------------------------------------------------
        static double Dist(double[] a, double[] b)
        {
            double s = 0;
            for (int i = 0; i < 4; i++) s += (a[i] - b[i]) * (a[i] - b[i]);
            return Math.Sqrt(s);
        }

        // -----------------------------------------------------------------------
        // MIConvexHull adapter types
        // -----------------------------------------------------------------------
        class HullVertex : IVertex
        {
            public double[] Position { get; }
            public int Index { get; }
            public HullVertex(double[] pos, int index) { Position = pos; Index = index; }
        }

        class HullFace : ConvexFace<HullVertex, HullFace> { }
    
        // -----------------------------------------------------------------------
        // Generator for 47 uniform convex polychora
        // -----------------------------------------------------------------------
        public static RegularPolytope4d Generate(string group, int mask)
        {
            var verts = PolychoraGenerator.GenerateVertices(group, mask);
            return FromVertices(verts);
        }

        public class PolychoraGenerator
        {
            public static List<double[]> GenerateVertices(string groupName, int activeNodes)
            {
                double[,] matrix = GetCoxeterMatrix(groupName);
                double[][] mirrors = GetMirrorNormals(matrix);
                double[] p = SolveGeneratorPoint(mirrors, activeNodes);
                return GenerateOrbit(p, mirrors);
            }

            private static double[,] GetCoxeterMatrix(string group)
            {
                switch (group.ToUpper())
                {
                    case "A4": return new double[,] { {1,3,2,2}, {3,1,3,2}, {2,3,1,3}, {2,2,3,1} };
                    case "B4": return new double[,] { {1,4,2,2}, {4,1,3,2}, {2,3,1,3}, {2,2,3,1} };
                    case "F4": return new double[,] { {1,3,2,2}, {3,1,4,2}, {2,4,1,3}, {2,2,3,1} };
                    case "H4": return new double[,] { {1,5,2,2}, {5,1,3,2}, {2,3,1,3}, {2,2,3,1} };
                    default: throw new ArgumentException("Unknown group");
                }
            }

            private static double[][] GetMirrorNormals(double[,] matrix)
            {
                int n = matrix.GetLength(0);
                double[][] normals = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    normals[i] = new double[n];
                    for (int j = 0; j < i; j++)
                    {
                        double dot = -Math.Cos(Math.PI / matrix[i, j]);
                        double sum = 0;
                        for (int k = 0; k < j; k++) sum += normals[i][k] * normals[j][k];
                        normals[i][j] = (dot - sum) / (normals[j][j] == 0 ? 1 : normals[j][j]);
                    }
                    double ssq = 0;
                    for (int k = 0; k < i; k++) ssq += normals[i][k] * normals[i][k];
                    normals[i][i] = Math.Sqrt(Math.Max(0, 1.0 - ssq));
                }
                return normals;
            }

            private static double[] SolveGeneratorPoint(double[][] mirrors, int activeNodes)
            {
                int n = mirrors.Length;
                double[] d = new double[n];
                for (int i = 0; i < n; i++)
                    d[i] = ((activeNodes >> i) & 1) != 0 ? 1.0 : 0.0;
                double[] p = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++) sum += mirrors[i][j] * p[j];
                    p[i] = (d[i] - sum) / (mirrors[i][i] == 0 ? 1 : mirrors[i][i]);
                }
                return p;
            }

            private static List<double[]> GenerateOrbit(double[] start, double[][] mirrors)
            {
                var vertices = new List<double[]>();
                var stack = new Stack<double[]>();
                stack.Push(start);
                vertices.Add(start);
                var seen = new HashSet<string>();
                seen.Add(VKey(start));
                while (stack.Count > 0)
                {
                    double[] v = stack.Pop();
                    for (int i = 0; i < mirrors.Length; i++)
                    {
                        double dot = 0;
                        for (int k = 0; k < 4; k++) dot += v[k] * mirrors[i][k];
                        double[] next = new double[4];
                        for (int k = 0; k < 4; k++) next[k] = v[k] - 2 * dot * mirrors[i][k];
                        string key = VKey(next);
                        if (!seen.Contains(key))
                        {
                            seen.Add(key);
                            vertices.Add(next);
                            stack.Push(next);
                        }
                    }
                }
                return vertices;
            }

            private static string VKey(double[] v)
            {
                return string.Join(",", v.Select(x => Math.Round(x, 8).ToString("F8")));
            }
        }
}
}
