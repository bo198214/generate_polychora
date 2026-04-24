using System;
using System.Collections.Generic;
using System.Linq;

namespace D4BB.Geometry
{
    /// <summary>Generates vertices for non-Wythoffian convex uniform 4-polytopes.</summary>
    public static class SnubGenerator
    {
        /// <summary>
        /// Snub 24-cell (snic/sadi): 96 vertices at all even permutations of (0, ±1, ±φ, ±φ²)
        /// where φ = (1+√5)/2 is the golden ratio. Edge length 2.
        /// </summary>
        public static List<double[]> SnubIcositetrachoron()
        {
            double phi  = (1.0 + Math.Sqrt(5.0)) / 2.0; // φ ≈ 1.618
            double phi2 = phi * phi;                      // φ² = φ+1 ≈ 2.618

            // Build all even permutations of (0, ε1, ε2*φ, ε3*φ²) with ε ∈ {+1,-1}
            var verts = new List<double[]>();
            double[] vals = { 0, 1, phi, phi2 };

            // All 4! = 24 permutations of indices {0,1,2,3}
            int[] idx = { 0, 1, 2, 3 };
            foreach (var perm in Permutations(idx))
            {
                if (!IsEvenPermutation(perm)) continue; // keep only even permutations (12 of them)

                // 0 has no sign, the others can be ±
                for (int s1 = -1; s1 <= 1; s1 += 2)
                for (int s2 = -1; s2 <= 1; s2 += 2)
                for (int s3 = -1; s3 <= 1; s3 += 2)
                {
                    double[] signs = { 1, s1, s2, s3 }; // position of 0 in perm gets sign 1 (×0)
                    var v = new double[4];
                    for (int k = 0; k < 4; k++)
                        v[k] = signs[perm[k]] * vals[perm[k]];
                    verts.Add(v);
                }
            }
            return verts; // 12 × 8 = 96 vertices
        }

        /// <summary>
        /// Grand antiprism: 600-cell (H4/an=8) minus two antipodal decagonal rings of 10 vertices each.
        /// The two rings lie in mutually orthogonal planes; their removal yields the 100-vertex grand antiprism.
        /// </summary>
        public static List<double[]> GrandAntiprism()
        {
            // Grand antiprism = 600-cell (H4/an=8) minus two orthogonal decagonal rings.
            // Ring A: 10 vertices with highest z²+w² fraction (most "ZW-polar")
            // Ring B: 10 vertices with highest x²+y² fraction (most "XY-polar")
            var v600 = PolychoraGenerator.GenerateVertices("H4", 8); // 120 vertices

            double ZWfrac(double[] v) => v[2]*v[2] + v[3]*v[3]; // z²+w²
            double XYfrac(double[] v) => v[0]*v[0] + v[1]*v[1]; // x²+y²

            var ringZW = Enumerable.Range(0, v600.Count)
                .OrderByDescending(i => ZWfrac(v600[i])).Take(10).ToHashSet();
            var ringXY = Enumerable.Range(0, v600.Count)
                .Where(i => !ringZW.Contains(i))
                .OrderByDescending(i => XYfrac(v600[i])).Take(10).ToHashSet();

            var remove = new HashSet<int>(ringZW);
            remove.UnionWith(ringXY);  // 10 + 10 = 20 removed

            return v600.Where((_, i) => !remove.Contains(i)).ToList(); // 120 - 20 = 100
        }

        // ── helpers ──────────────────────────────────────────────────────────

        static IEnumerable<int[]> Permutations(int[] arr)
        {
            if (arr.Length == 1) { yield return (int[])arr.Clone(); yield break; }
            for (int i = 0; i < arr.Length; i++)
            {
                int[] rest = new int[arr.Length - 1];
                int ri = 0;
                for (int j = 0; j < arr.Length; j++) if (j != i) rest[ri++] = arr[j];
                foreach (var sub in Permutations(rest))
                {
                    var p = new int[arr.Length];
                    p[0] = arr[i];
                    Array.Copy(sub, 0, p, 1, sub.Length);
                    yield return p;
                }
            }
        }

        /// <summary>Counts inversions; even permutation iff inversion count is even.</summary>
        static bool IsEvenPermutation(int[] p)
        {
            int inv = 0;
            for (int i = 0; i < p.Length; i++)
                for (int j = i + 1; j < p.Length; j++)
                    if (p[i] > p[j]) inv++;
            return inv % 2 == 0;
        }
    }
}
