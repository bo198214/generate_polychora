using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using D4BB.Geometry;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    /// <summary>
    /// Finds the correct (group, activeNodes) for remaining failing polytopes
    /// and regenerates their JSON files.
    /// </summary>
    [TestFixture]
    public class ExtendedSearchTest
    {
        static readonly string OutDir = Path.GetFullPath(Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "output"));

        // All expected polytope topologies
        static readonly (string name, int v, int e, int f, int c)[] AllExpected = {
            // B4 family
            ("tah",   64, 128,  80, 16), ("rico",  32,  96,  80, 16),
            ("tico",  48, 120,  88, 16), ("thic", 192, 384, 224, 32),
            ("xic",   64, 192, 160, 32), ("scic", 192, 480, 320, 32), ("gic", 384, 768, 448, 64),
            // F4 family
            ("tico_f",96,288,216,24),("frico",288,864,672,96),
            ("grico",576,1152,816,240),("drico",576,1440,1008,144),("trico",1152,2304,1344,192),
            // H4 family - small enough to compute
            ("thi",2400,4800,3120,720),("tex",1440,4320,3600,720),
        };

        [Test, CancelAfter(600000)]
        public void FindAndRegenerateH4()
        {
            // H4: check an=3 (V=2400=thi?) and an=12 (V=1440=tex?)
            // Plus an=6,9,10,11,13,14,15 for the others
            var h4Candidates = new[] { 3, 6, 9, 10, 11, 12, 13, 14, 15 };
            var byVEFC = AllExpected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);

            foreach (int an in h4Candidates)
            {
                var pts = PolychoraGenerator.GenerateVertices("H4", an);
                Console.WriteLine($"H4/an={an}: V={pts.Count}");
                if (pts.Count > 4000) { Console.WriteLine("  (too large, skipping)"); continue; }
                var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                var key = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                string match = byVEFC.TryGetValue(key, out var n) ? $" <-- {n} ✓" : "";
                Console.WriteLine($"  hull ({key.Item1},{key.Item2},{key.Item3},{key.Item4}){match}");
                if (match.Length > 0) {
                    SaveJson(n, $"H4/{an}", pts);
                    Console.WriteLine($"  Saved {n}.json");
                }
            }
            Assert.Pass();
        }

        [Test, CancelAfter(300000)]
        public void FindF4LargeActiveNodes()
        {
            // Compute topology for F4 an=7,11,13,14,15 (V=576 or 1152)
            var byVEFC = AllExpected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);
            var candidates = new[] { 5, 7, 10, 11, 13, 14, 15 };
            Console.WriteLine("=== F4 large activeNodes ===");
            foreach (int an in candidates) {
                var pts = PolychoraGenerator.GenerateVertices("F4", an);
                Console.WriteLine($"F4/an={an}: V={pts.Count}");
                if (pts.Count > 1200) { Console.WriteLine("  (too large)"); continue; }
                var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                var key = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                string match = byVEFC.TryGetValue(key, out var n) ? $" <-- {n} ✓" : "";
                Console.WriteLine($"  hull ({key.Item1},{key.Item2},{key.Item3},{key.Item4}){match}");
                if (match.Length > 0) { SaveJson(n, $"F4/{an}", pts); Console.WriteLine($"  Saved {n}.json"); }
            }
            Assert.Pass();
        }

        [Test, CancelAfter(120000)]
        public void FindB4F4Remaining()
        {
            // Check ALL B4 and F4 activeNodes for remaining matches
            var byVEFC = AllExpected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);
            foreach (var (grp, maxV) in new[] { ("B4", 500), ("F4", 600) })
            {
                Console.WriteLine($"\n=== {grp} ===");
                for (int an = 1; an <= 15; an++) {
                    var pts = PolychoraGenerator.GenerateVertices(grp, an);
                    if (pts.Count > maxV) { continue; }
                    var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                    var key = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                    if (byVEFC.TryGetValue(key, out var n)) {
                        Console.WriteLine($"  {grp}/an={an}: ({key.Item1},{key.Item2},{key.Item3},{key.Item4}) <-- {n} ✓");
                        SaveJson(n, $"{grp}/{an}", pts);
                    }
                }
            }
            Assert.Pass();
        }

        void SaveJson(string name, string source, List<double[]> verts) {
            var json = JsonSerializer.Serialize(new { name, source, vertices = verts },
                new JsonSerializerOptions { WriteIndented = true });
            File.WriteAllText(Path.Combine(OutDir, $"{name}.json"), json);
        }
    }
}
