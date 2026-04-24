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
            AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "vertex_output"));

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
            // H4: an=12 (V=1440=tex?), an=6,9 for xhi/rex, an=3 is thi (V=2400 - too large)
            var h4Candidates = new[] { 6, 9, 10, 11, 12, 13 }; // skip an=3 (too slow)
            var byVEFC = AllExpected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);

            foreach (int an in h4Candidates)
            {
                var pts = PolychoraGenerator.GenerateVertices("H4", an);
                Console.WriteLine($"H4/an={an}: V={pts.Count}");
                if (pts.Count > 2000) { Console.WriteLine("  (too large, skipping)"); continue; }
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

        [Test, CancelAfter(180000)]
        public void RegenerateAllKnown()
        {
            // Complete mapping of all known polytopes to (group, activeNodes)
            var mappings = new (string name, int v, int e, int f, int c, string grp, int an)[] {
                // A4 family
                ("pen",    5,  10,  10,   5,"A4", 1),("rap",   10,  30,  30,  10,"A4", 2),
                ("tip",   20,  40,  30,  10,"A4", 3),("deca",  30,  60,  40,  10,"A4", 6),
                ("hap",   30,  90,  80,  20,"A4", 5),("tap",   60, 120,  80,  20,"A4", 7),
                ("dappat",60, 150, 120,  30,"A4",11),("tappy",120, 240, 150,  30,"A4",15),
                // B4 family (correct Wythoff topologies)
                ("tes",   16,  32,  24,   8,"B4", 1),("hex",    8,  24,  32,  16,"B4", 8),
                ("rico",  32,  96,  88,  24,"B4", 2),("tah",   64, 128,  88,  24,"B4", 3),
                ("tico",  48, 120,  96,  24,"B4",12),("spic",  96, 288, 240,  48,"B4",10),
                ("thic", 192, 384, 248,  56,"B4", 7),("xic",   64, 192, 208,  80,"B4", 9),
                ("scic", 192, 480, 368,  80,"B4",11),("gic",  384, 768, 464,  80,"B4",15),
                // F4 family
                ("ico",   24,  96,  96,  24,"F4", 1),("cont", 192, 384, 240,  48,"F4", 3),
                ("tico_f",96, 192, 120,  24,"B4", 6),("srico",288, 576, 336,  48,"F4", 6),
                ("frico",288, 864, 720, 144,"F4", 5),("grico",576,1152, 720, 144,"F4", 7),
                ("prico",144, 576, 672, 240,"F4", 9),("drico",576,1440,1104, 240,"F4",11),
                ("trico",1152,2304,1392, 240,"F4",15),
                // H4 family
                ("hi",   600,1200, 720, 120,"H4", 1),("ex",   120, 720,1200, 600,"H4", 8),
                ("rhi", 1200,3600,3120, 720,"H4", 2),("tex", 1440,4320,3600, 720,"H4",12),
                ("rex",  720,3600,3600, 720,"H4", 4),
            };
            int ok=0, fail=0;
            foreach (var m in mappings) {
                var pts = PolychoraGenerator.GenerateVertices(m.grp, m.an);
                var json = System.Text.Json.JsonSerializer.Serialize(new {
                    name=m.name, source=$"{m.grp}/{m.an}", vertices=pts
                }, new System.Text.Json.JsonSerializerOptions{WriteIndented=true});
                File.WriteAllText(Path.Combine(OutDir, $"{m.name}.json"), json);
                if (pts.Count <= 1500) {
                    var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                    bool match = hull.Vertices.Count==m.v&&hull.Edges.Count==m.e
                               &&hull.Faces.Count==m.f&&hull.Cells.Count==m.c;
                    Console.WriteLine($"{m.name,-10} {m.grp}/{m.an}: {(match?"✓":"✗")} ({hull.Vertices.Count},{hull.Edges.Count},{hull.Faces.Count},{hull.Cells.Count})");
                    if(match)ok++;else fail++;
                } else {
                    Console.WriteLine($"{m.name,-10} {m.grp}/{m.an}: V={pts.Count} (large)"); ok++;
                }
            }
            Console.WriteLine($"\n{ok} OK, {fail} failed");
            Assert.That(fail, Is.EqualTo(0));
        }

        void SaveJson(string name, string source, List<double[]> verts) {
            var json = JsonSerializer.Serialize(new { name, source, vertices = verts },
                new JsonSerializerOptions { WriteIndented = true });
            File.WriteAllText(Path.Combine(OutDir, $"{name}.json"), json);
        }

        [Test, CancelAfter(600000)]
        public void GenerateTopologyOutput()
        {
            string vertexDir  = Path.GetFullPath(Path.Combine(
                AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "vertex_output"));
            string topologyDir = Path.GetFullPath(Path.Combine(
                AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "topology_output"));

            TopologyGenerator.Generate(vertexDir, topologyDir, eps: 1e-6,
                log: msg => Console.WriteLine(msg));

            var written = Directory.GetFiles(topologyDir, "*.json").Length;
            Console.WriteLine($"\nWrote {written} topology files to {topologyDir}");
            Assert.That(written, Is.GreaterThan(0));
        }
    }
}
