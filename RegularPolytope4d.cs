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
        public List<List<int[]>> cells;
        public List<double[]> cellNormals;
        public double edgeLength;

        public static RegularPolytope4d FromVertices(List<double[]> verts, double hullTolerance = 1e-9)
        {
            var p = new RegularPolytope4d { vertices = verts };
            p.edgeLength = ComputeEdgeLength(verts);
            p.edges      = ComputeEdges(verts, p.edgeLength);
            (p.cells, p.cellNormals) = ComputeCells(verts, hullTolerance);
            return p;
        }

        static double ComputeEdgeLength(List<double[]> verts)
        {
            double min = double.MaxValue;
            int n = Math.Min(verts.Count, 100);
            for (int i = 0; i < n; i++)
            for (int j = i + 1; j < verts.Count; j++)
            {
                double d = Dist(verts[i], verts[j]);
                if (d > 1e-6 && d < min - 1e-7) min = d;
            }
            return min;
        }

        static List<(int, int)> ComputeEdges(List<double[]> verts, double edgeLen)
        {
            var edges = new List<(int, int)>();
            for (int i = 0; i < verts.Count; i++)
            for (int j = i + 1; j < verts.Count; j++)
                if (Math.Abs(Dist(verts[i], verts[j]) - edgeLen) < 1e-3)
                    edges.Add((i, j));
            return edges;
        }

        static (List<List<int[]>>, List<double[]>) ComputeCells(List<double[]> verts, double tolerance = 1e-9)
        {
            var hullVerts = verts.Select((v, i) => new HullVertex(v, i)).ToList();
            var hull = ConvexHull.Create<HullVertex, HullFace>(hullVerts, tolerance);
            if (hull?.Faces == null) throw new Exception("ConvexHull failed");

            var rawFaces = hull.Faces.ToArray();
            int nFaces = rawFaces.Length;

            // Grouping by adjacency and coplanarity
            var parent = Enumerable.Range(0, nFaces).ToArray();
            
            // Build adjacency map (faces sharing a 2D ridge = 3 vertices)
            var ridgeToFaces = new Dictionary<string, List<int>>();
            for (int i = 0; i < nFaces; i++)
            {
                int[] v = rawFaces[i].Vertices.Select(vh => ((HullVertex)vh).Index).OrderBy(x => x).ToArray();
                // 4 vertices per face (tetrahedron), 4 ridges (triangles)
                int[][] ridges = { 
                    new[] { v[0], v[1], v[2] }, new[] { v[0], v[1], v[3] },
                    new[] { v[0], v[2], v[3] }, new[] { v[1], v[2], v[3] }
                };
                foreach (var r in ridges) {
                    string key = string.Join(",", r);
                    if (!ridgeToFaces.ContainsKey(key)) ridgeToFaces[key] = new List<int>();
                    ridgeToFaces[key].Add(i);
                }
            }

            for (int i = 0; i < nFaces; i++)
            {
                var ni = rawFaces[i].Normal;
                double di = 0; 
                for (int k = 0; k < 4; k++) di += ni[k] * rawFaces[i].Vertices[0].Position[k];

                // For each ridge of face i, check the other face sharing it
                int[] vi = rawFaces[i].Vertices.Select(vh => ((HullVertex)vh).Index).OrderBy(x => x).ToArray();
                int[][] ridges = { 
                    new[] { vi[0], vi[1], vi[2] }, new[] { vi[0], vi[1], vi[3] },
                    new[] { vi[0], vi[2], vi[3] }, new[] { vi[1], vi[2], vi[3] }
                };
                foreach (var r in ridges) {
                    string key = string.Join(",", r);
                    foreach (int j in ridgeToFaces[key]) {
                        if (i == j) continue;
                        var nj = rawFaces[j].Normal;
                        double dot = 0; for (int k = 0; k < 4; k++) dot += ni[k] * nj[k];
                        if (dot > 0.99) UFUnion(parent, i, j);
                    }
                }
            }

            var groups = new Dictionary<int, (List<int[]> tris, double[] norm)>();
            for (int i = 0; i < nFaces; i++)
            {
                int root = UFFind(parent, i);
                if (!groups.ContainsKey(root))
                {
                    double[] avgN = new double[4];
                    for (int k = 0; k < 4; k++) avgN[k] = rawFaces[root].Normal[k];
                    groups[root] = (new List<int[]>(), avgN);
                }
                groups[root].tris.Add(rawFaces[i].Vertices.Select(vh => ((HullVertex)vh).Index).ToArray());
            }

            return (groups.Values.Select(g => g.tris).ToList(),
                    groups.Values.Select(g => g.norm).ToList());
        }

        static int UFFind(int[] parent, int i) => parent[i] == i ? i : (parent[i] = UFFind(parent, parent[i]));
        static void UFUnion(int[] parent, int i, int j) { parent[UFFind(parent, i)] = UFFind(parent, j); }

        public List<int> VisibleCells(double[] eye)
        {
            var visible = new List<int>();
            for (int i = 0; i < cells.Count; i++)
            {
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

        static double Dist(double[] a, double[] b)
        {
            double s = 0;
            for (int i = 0; i < 4; i++) s += (a[i] - b[i]) * (a[i] - b[i]);
            return Math.Sqrt(s);
        }

        class HullVertex : IVertex {
            public double[] Position { get; }
            public int Index { get; }
            public HullVertex(double[] pos, int index) { Position = pos; Index = index; }
        }
        class HullFace : ConvexFace<HullVertex, HullFace> { }

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
                for (int i = 0; i < n; i++) d[i] = ((activeNodes >> i) & 1) != 0 ? 1.0 : 0.0;
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
                        if (Math.Abs(dot) < 1e-8) continue;
                        double[] next = new double[4];
                        for (int k = 0; k < 4; k++) next[k] = v[k] - 2 * dot * mirrors[i][k];
                        string key = VKey(next);
                        if (!seen.Contains(key)) { seen.Add(key); vertices.Add(next); stack.Push(next); }
                    }
                }
                return vertices;
            }

            private static string VKey(double[] v) { return string.Join(",", v.Select(x => Math.Round(x, 8).ToString("F8"))); }
        }
    }
}
