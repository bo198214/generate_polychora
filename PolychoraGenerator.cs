using System;
using System.Collections.Generic;
using System.Linq;

namespace D4BB.Geometry
{
    public class PolychoraGenerator
    {
        public static List<double[]> GenerateVertices(string groupName, int activeNodes)
        {
            double[,] matrix = GetCoxeterMatrix(groupName);
            int n = matrix.GetLength(0);
            double[][] mirrors = GetMirrorNormals(matrix);

            // Generator point P such that P . mirrors[i] = 1 if bit i is set, else 0
            double[] p = SolveGeneratorPoint(mirrors, activeNodes);

            return GenerateOrbit(p, mirrors);
        }

        private static double[,] GetCoxeterMatrix(string group)
        {
            switch (group.ToUpper())
            {
                case "A4": return new double[,] { {1,3,2,2}, {3,1,3,2}, {2,3,1,3}, {2,2,3,1} };
                case "B4": return new double[,] { {1,4,2,2}, {4,1,3,2}, {2,3,1,3}, {2,2,3,1} };
                // B4R = B4 with nodes reversed (4-bond at positions 2-3)
                case "B4R": return new double[,] { {1,3,2,2}, {3,1,3,2}, {2,3,1,4}, {2,2,4,1} };
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

            // Solve N * P = D
            // Since N is lower triangular (from GetMirrorNormals), use forward substitution
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
                    
                    // Reflect: v' = v - 2(v.n)n
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
            // +0.0 normalizes IEEE 754 negative zero (-0.0) to avoid duplicate vertices
            return string.Join(",", v.Select(x => (Math.Round(x, 8) + 0.0).ToString("F8")));
        }
    }
}
