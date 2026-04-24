using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using D4BB.Geometry;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    [TestFixture]
    public class DiagTest
    {
        [Test]
        public void DiagBitruncated()
        {
            // Bitruncated tesseract vertices: permutations of (+-1, +-1, +-1, +-(1+sqrt2))
            // with positive product (even number of negative signs)
            double b = 1.0 + Math.Sqrt(2);
            var pts = new List<double[]>();
            int[] positions = {0,1,2,3};
            foreach (int bPos in positions) {
                var others = positions.Where(i => i != bPos).ToArray();
                foreach (double bSign in new[]{1.0,-1.0}) {
                    for (int mask = 0; mask < 8; mask++) {
                        double s0 = (mask&1) != 0 ? -1 : 1;
                        double s1 = (mask&2) != 0 ? -1 : 1;
                        double s2 = (mask&4) != 0 ? -1 : 1;
                        if (bSign * s0 * s1 * s2 > 0) { // even parity
                            var v = new double[4];
                            v[bPos] = bSign * b;
                            v[others[0]] = s0; v[others[1]] = s1; v[others[2]] = s2;
                            pts.Add(v);
                        }
                    }
                }
            }
            Console.WriteLine($"V={pts.Count}");
            var hull = TrueConvexHull4D.Compute(pts, 1e-6);
            Console.WriteLine($"hull: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
            Console.WriteLine($"Expected: V=32 E=96 F=80 C=16");
            Assert.Pass();
        }

        [Test]
        public void DiagReversedB4()
        {
            // Reversed B4 matrix: 4-bond at nodes 2-3 instead of 0-1
            // Predicts: an=4 gives V=32 (bitruncated tesseract = rico?)
            var reversedB4 = "B4r";
            // We implement inline for now
            // Expected: reversed an=8=tes, an=1=hex, an=4=rico...
            // Actually let's just run the generator with reversed node mapping
            // by trying reversed activeNodes for B4:
            // If original uses matrix[0,1]=4, reversed has matrix[2,3]=4
            // We'll use the F4 matrix trick: check if F4 gives rico topology
            var expectedRico = (32, 96, 80, 16);
            for (int an = 1; an <= 15; an++) {
                var pts = PolychoraGenerator.GenerateVertices("B4", an);
                if (pts.Count != 32) continue;
                var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                Console.WriteLine($"B4 an={an}: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
            }
            Assert.Pass();
        }

        [Test]
        public void DiagB4Orbits()
        {
            for (int an = 1; an <= 15; an++)
            {
                var pts = PolychoraGenerator.GenerateVertices("B4", an);
                Console.WriteLine($"B4 an={an,2}: {pts.Count,5} vertices, first={string.Join(",", pts[0].Select(x=>x.ToString("F2")))}");
            }
            Assert.Pass();
        }

        [Test]
        public void DiagRico()
        {
            var json = File.ReadAllText(@"P:\workspace\generate_polychora\output\rico.json");
            using var doc = JsonDocument.Parse(json);
            var pts = doc.RootElement.GetProperty("vertices")
                .EnumerateArray()
                .Select(v => v.EnumerateArray().Select(x => x.GetDouble()).ToArray())
                .ToList();

            var hull = TrueConvexHull4D.Compute(pts, 1e-6);
            Console.WriteLine($"rico: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
            Console.WriteLine($"Euler: {hull.Vertices.Count - hull.Edges.Count + hull.Faces.Count - hull.Cells.Count}");
            var cellSizes = hull.Cells.GroupBy(c => c.Length).OrderBy(g => g.Key);
            foreach (var g in cellSizes) Console.WriteLine($"  cells with {g.Key} verts: {g.Count()}");
            Assert.Pass();
        }

        [Test]
        public void RegenerateKnownPolychora()
        {
            // Known mappings: (name, expectedV, expectedE, expectedF, expectedC, group, activeNodes)
            var mappings = new (string name, int v, int e, int f, int c, string group, int an)[] {
                ("pen",    5,  10,  10,   5, "A4",  1), ("rap",   10,  30,  30,  10, "A4",  2),
                ("tip",   20,  40,  30,  10, "A4",  3), ("deca",  30,  60,  40,  10, "A4",  6),
                ("hap",   30,  90,  80,  20, "A4",  5), ("tap",   60, 120,  80,  20, "A4",  7),
                ("dappat",60, 150, 120,  30, "A4", 11), ("tappy",120, 240, 150,  30, "A4", 15),
                ("tes",   16,  32,  24,   8, "B4",  1), ("hex",    8,  24,  32,  16, "B4",  8),
                ("spic",  96, 288, 240,  48, "B4", 10),
                ("ico",   24,  96,  96,  24, "F4",  1), ("cont", 192, 384, 240,  48, "F4",  3),
                ("srico",288, 576, 336,  48, "F4",  6), ("prico",144, 576, 672, 240, "F4",  9),
                ("hi",   600,1200, 720, 120, "H4",  1), ("ex",   120, 720,1200, 600, "H4",  8),
                ("rhi", 1200,3600,3120, 720, "H4",  2), ("rex",  720,3600,3600, 720, "H4",  4),
            };
            string outDir = Path.GetFullPath(Path.Combine(
                AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "output"));
            int ok = 0, fail = 0;
            foreach (var m in mappings) {
                var pts = PolychoraGenerator.GenerateVertices(m.group, m.an);
                var json = JsonSerializer.Serialize(new {
                    name = m.name, group = $"{m.group}/{m.an}", vertices = pts
                }, new JsonSerializerOptions { WriteIndented = true });
                File.WriteAllText(Path.Combine(outDir, $"{m.name}.json"), json);
                if (pts.Count <= 1500) {
                    var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                    bool match = hull.Vertices.Count==m.v && hull.Edges.Count==m.e
                               && hull.Faces.Count==m.f && hull.Cells.Count==m.c;
                    Console.WriteLine($"{m.name,-10} {m.group}/{m.an}: {(match?"✓":"✗")} got ({hull.Vertices.Count},{hull.Edges.Count},{hull.Faces.Count},{hull.Cells.Count})");
                    if (match) ok++; else fail++;
                } else {
                    Console.WriteLine($"{m.name,-10} {m.group}/{m.an}: V={pts.Count} (large)");
                    ok++;
                }
            }
            Console.WriteLine($"\n{ok} regenerated OK, {fail} failed");
            Assert.That(fail, Is.EqualTo(0), "Some regenerated polytopes don't match expected topology");
        }
    }
}
