using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    [TestFixture]
    public class RegularPolytope4dTests
    {
        // Expected topologies for all 47 non-prismatic uniform polychora.
        // Values are verified against the pre-computed topology_output/ JSON files.
        static readonly (string name, int v, int e, int f, int c)[] All47 = {
            // A4 family
            ("pen", 5, 10, 10, 5), ("tip", 20, 40, 30, 10), ("rap", 10, 30, 30, 10), ("deca", 30, 60, 40, 10),
            ("srip", 30, 90, 80, 20), ("grip", 60, 120, 80, 20), ("spid", 20, 60, 70, 30), ("prip", 60, 150, 120, 30),
            ("gippid", 120, 240, 150, 30),
            // B4 family
            ("tes", 16, 32, 24, 8), ("hex", 8, 24, 32, 16), ("tat", 64, 128, 88, 24), ("rit", 32, 96, 88, 24),
            ("thex", 48, 120, 96, 24), ("rico", 96, 288, 240, 48), ("grit", 192, 384, 248, 56), ("sidpith", 64, 192, 208, 80),
            ("proh", 192, 480, 368, 80), ("gidpith", 384, 768, 464, 80),
            ("srit", 96, 288, 248, 56), ("prit", 192, 480, 368, 80),
            // F4 family
            ("ico", 24, 96, 96, 24), ("tico", 192, 384, 240, 48), ("tah", 96, 192, 120, 24), ("cont", 288, 576, 336, 48),
            ("srico", 288, 864, 720, 144), ("grico", 576, 1152, 720, 144), ("spic", 144, 576, 672, 240), ("prico", 576, 1440, 1104, 240),
            ("gippic", 1152, 2304, 1392, 240),
            // H4 family
            ("hi", 600, 1200, 720, 120), ("ex", 120, 720, 1200, 600), ("rahi", 1200, 3600, 3120, 720), ("thi", 2400, 4800, 3120, 720),
            ("rox", 720, 3600, 3600, 720), ("tex", 1440, 4320, 3600, 720), ("xhi", 3600, 7200, 4320, 720), ("grahi", 7200, 14400, 9120, 1920),
            ("prix", 7200, 18000, 13440, 2640), ("srahi", 3600, 10800, 9120, 1920), ("srix", 3600, 10800, 8640, 1440), ("sidpixhi", 2400, 7200, 7440, 2640),
            // Large H4 polytopes — values verified against Klitzing/Wikipedia
            ("grix", 7200, 14400, 8640, 1440),
            ("prahi", 7200, 18000, 13440, 2640),
            ("gidpixhi", 14400, 28800, 17040, 2640),
            // Non-Wythoffian
            ("gap", 100, 500, 720, 320), ("sadi", 96, 432, 480, 144)
        };

        static readonly string TopoDir = Path.GetFullPath(Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "topology_output"));

        [TestCaseSource(nameof(All47))]
        public void TopologyFileVerification((string name, int v, int e, int f, int c) expected)
        {
            var path = Path.Combine(TopoDir, $"{expected.name}.json");
            if (!File.Exists(path))
                Assert.Ignore($"{expected.name}: topology_output/{expected.name}.json not yet generated (run 'make topology')");

            var (av, ae, af, ac) = ReadCounts(path);

            Assert.That(av, Is.EqualTo(expected.v), $"{expected.name}: V mismatch");
            Assert.That(ae, Is.EqualTo(expected.e), $"{expected.name}: E mismatch");
            Assert.That(af, Is.EqualTo(expected.f), $"{expected.name}: F mismatch");
            Assert.That(ac, Is.EqualTo(expected.c), $"{expected.name}: C mismatch");
            Assert.That(av - ae + af - ac, Is.EqualTo(0), $"{expected.name}: Euler characteristic violated");
        }

        static (int v, int e, int f, int c) ReadCounts(string path)
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            var r = doc.RootElement;
            return (
                r.GetProperty("vertices").GetArrayLength(),
                r.GetProperty("edges").GetArrayLength(),
                r.GetProperty("faces2d").GetArrayLength(),
                r.GetProperty("cells").GetArrayLength()
            );
        }
    }
}
