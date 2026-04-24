using System;
using System.Collections.Generic;
using System.Linq;
using MIConvexHull;

namespace D4BB.Geometry
{
    public class TrueConvexHull4D
    {
        public List<double[]> Vertices;
        public List<int[]> Edges;
        public List<int[]> Faces;
        public List<int[]> Cells;

        public static TrueConvexHull4D Compute(List<double[]> inputVerts, double epsilon = 1e-7)
        {
            var hull = new TrueConvexHull4D();
            hull.Vertices = inputVerts;

            var vertexData = inputVerts.Select((v, i) => new HullVertex(v, i)).ToList();
            var simplicialHull = ConvexHull.Create<HullVertex, HullFace>(vertexData, epsilon);
            if (simplicialHull?.Faces == null) throw new Exception("Simplicial Hull failed");

            // 1. Group tetrahedra into 3D cells via unique hyperplanes
            var cellVertices = new Dictionary<string, List<int>>();
            foreach (var f in simplicialHull.Faces)
            {
                var n = f.Normal;
                double d = 0; for (int k = 0; k < 4; k++) d += n[k] * f.Vertices[0].Position[k];
                string key = $"{Math.Round(n[0], 6):F6},{Math.Round(n[1], 6):F6},{Math.Round(n[2], 6):F6},{Math.Round(n[3], 6):F6}|{Math.Round(d, 5):F5}";
                if (!cellVertices.ContainsKey(key)) cellVertices[key] = new List<int>();
                foreach (var v in f.Vertices) {
                    int idx = ((HullVertex)v).Index;
                    if (!cellVertices[key].Contains(idx)) cellVertices[key].Add(idx);
                }
            }
            hull.Cells = cellVertices.Values.Select(l => l.ToArray()).ToList();

            // 2. Identify 1D edges (distance based for uniform polychora)
            double minD = double.MaxValue;
            for (int i = 0; i < Math.Min(inputVerts.Count, 100); i++)
                for (int j = i + 1; j < inputVerts.Count; j++) {
                    double d = Dist(inputVerts[i], inputVerts[j]);
                    if (d > 1e-6 && d < minD - 1e-7) minD = d;
                }
            
            var edgeList = new List<int[]>();
            for (int i = 0; i < inputVerts.Count; i++)
                for (int j = i + 1; j < inputVerts.Count; j++)
                    if (Math.Abs(Dist(inputVerts[i], inputVerts[j]) - minD) < 1e-4)
                        edgeList.Add(new[] { i, j });
            hull.Edges = edgeList;

            // 3. Identify 2D faces as unique intersections of two 3D cells
            var facesSet = new HashSet<string>();
            var cellsArr = hull.Cells.ToArray();
            for (int i = 0; i < cellsArr.Length; i++)
                for (int j = i + 1; j < cellsArr.Length; j++)
                {
                    var intersection = cellsArr[i].Intersect(cellsArr[j]).OrderBy(x => x).ToList();
                    if (intersection.Count >= 3)
                        facesSet.Add(string.Join(",", intersection));
                }
            hull.Faces = facesSet.Select(s => s.Split(',').Select(int.Parse).ToArray()).ToList();

            return hull;
        }

        static double Dist(double[] a, double[] b) {
            double s = 0; for (int i = 0; i < 4; i++) s += (a[i]-b[i])*(a[i]-b[i]);
            return Math.Sqrt(s);
        }

        class HullVertex : IVertex {
            public double[] Position { get; }
            public int Index { get; }
            public HullVertex(double[] pos, int index) { Position = pos; Index = index; }
        }
        class HullFace : ConvexFace<HullVertex, HullFace> { }
    }
}
