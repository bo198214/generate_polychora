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
    public class RegularPolytope4dTests
    {
        // Theoretical Topologies for all 47 non-prismatic uniform polychora
        // (name, v, e, f, c)
        static readonly (string name, int v, int e, int f, int c)[] All47 = {
            // A4 family
            ("pen", 5, 10, 10, 5), ("tip", 20, 40, 30, 10), ("rap", 10, 30, 30, 10), ("deca", 30, 60, 40, 10),
            ("hap", 30, 90, 80, 20), ("tap", 60, 120, 80, 20), ("sadi", 20, 60, 60, 20), ("dappat", 60, 150, 120, 30),
            ("tappy", 120, 240, 150, 30),
            // B4 family (Wikipedia B4 topologies: tat, rit, thex, grit, sidpith, proh, gidpith)
            ("tes", 16, 32, 24, 8), ("hex", 8, 24, 32, 16), ("tah", 64, 128, 88, 24), ("rico", 32, 96, 88, 24),
            ("tico", 48, 120, 96, 24), ("spic", 96, 288, 240, 48), ("thic", 192, 384, 248, 56), ("xic", 64, 192, 208, 80),
            ("scic", 192, 480, 368, 80), ("gic", 384, 768, 464, 80),
            ("srit", 96, 288, 248, 56), ("prit", 192, 480, 368, 80),
            // F4 family
            ("ico", 24, 96, 96, 24), ("cont", 192, 384, 240, 48), ("tico_f", 96, 192, 120, 24), ("srico", 288, 576, 336, 48),
            ("frico", 288, 864, 720, 144), ("grico", 576, 1152, 720, 144), ("prico", 144, 576, 672, 240), ("drico", 576, 1440, 1104, 240),
            ("trico", 1152, 2304, 1392, 240),
            // H4 family
            ("hi", 600, 1200, 720, 120), ("ex", 120, 720, 1200, 600), ("rhi", 1200, 3600, 3120, 720), ("thi", 2400, 4800, 3120, 720),
            ("rex", 720, 3600, 3600, 720), ("tex", 1440, 4320, 3600, 720), ("xhi", 3600, 7200, 4320, 720), ("tphi", 7200, 14400, 7920, 720),
            ("tpi", 7200, 14400, 8640, 1440), ("sphi", 3600, 10800, 8640, 1440), ("spex", 3600, 10800, 9600, 2400), ("xex", 2400, 7200, 6000, 1200),
            ("dphi", 7200, 18000, 12240, 1440), ("dex", 7200, 18000, 13200, 2400), ("gishi", 14400, 28800, 15840, 1440),
            // Non-Wythoffian
            ("gap", 100, 500, 720, 320), ("snic", 96, 432, 480, 144)
        };

        [TestCaseSource(nameof(All47))]
        public void TrueConvexHullVerification((string name, int v, int e, int f, int c) expected)
        {
            var verts = LoadJson(expected.name);
            var hull = TrueConvexHull4D.Compute(verts, 1e-6);

            Assert.That(hull.Vertices.Count, Is.EqualTo(expected.v), $"{expected.name}: V mismatch");
            Assert.That(hull.Edges.Count, Is.EqualTo(expected.e), $"{expected.name}: E mismatch");
            Assert.That(hull.Faces.Count, Is.EqualTo(expected.f), $"{expected.name}: F mismatch");
            Assert.That(hull.Cells.Count, Is.EqualTo(expected.c), $"{expected.name}: C mismatch");
            
            // Euler check: V - E + F - C = 0
            Assert.That(hull.Vertices.Count - hull.Edges.Count + hull.Faces.Count - hull.Cells.Count, Is.EqualTo(0), $"{expected.name}: Euler characteristic violated");
        }

        static List<double[]> LoadJson(string name)
        {
            var candidates = new[] {
                Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "Assets", "uniform_polychora", $"{name}.json"),
                Path.Combine("P:\\workspace\\generate_polychora\\vertex_output", $"{name}.json")
            };
            foreach (var path in candidates)
                if (File.Exists(path)) return ParseVertices(File.ReadAllText(path));
            throw new FileNotFoundException($"{name}.json not found");
        }

        static List<double[]> ParseVertices(string json)
        {
            using var doc = JsonDocument.Parse(json);
            return doc.RootElement.GetProperty("vertices")
                .EnumerateArray()
                .Select(v => v.EnumerateArray().Select(x => x.GetDouble()).ToArray())
                .ToList();
        }
    }
}
