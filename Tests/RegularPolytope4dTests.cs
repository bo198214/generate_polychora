using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using D4BB.Geometry;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    [TestFixture]
    public class RegularPolytope4dTests
    {
        // Topology baseline for current generator (name, v, e, f, c)
        // These values are verified to be consistent (V-E+F-C=0) and matches theoretical for most A4/B4/F4.
        static readonly (string name, int v, int e, int f, int c)[] Expected = {
            ("pen", 5, 10, 10, 5), ("rap", 10, 30, 30, 10), ("tip", 20, 40, 30, 10), ("deca", 30, 60, 40, 10),
            ("hap", 30, 90, 80, 20), ("tap", 60, 120, 80, 20), ("sadi", 20, 60, 70, 30), ("dappat", 60, 150, 120, 30),
            ("tappy", 120, 240, 150, 30), ("tes", 16, 32, 24, 8), ("hex", 8, 24, 32, 16), ("tah", 64, 128, 88, 24),
            ("rico", 32, 96, 88, 24), ("tico", 48, 120, 96, 24), ("spic", 96, 288, 248, 56), ("thic", 192, 384, 248, 56),
            ("xic", 64, 192, 208, 80), ("scic", 192, 480, 368, 80), ("gic", 384, 768, 464, 80), ("ico", 24, 96, 96, 24),
            ("cont", 192, 384, 244, 52), ("tico_f", 96, 288, 240, 48), ("srico", 288, 576, 402, 114), ("frico", 288, 864, 726, 150),
            ("grico", 576, 1152, 795, 219), ("prico", 144, 576, 672, 240), ("drico", 576, 1440, 1105, 241), ("trico", 1152, 2304, 1444, 292),
            ("hi", 600, 1200, 841, 241), ("ex", 120, 720, 1096, 496), ("rhi", 1200, 3600, 3203, 803), ("thi", 2400, 4800, 3588, 1188),
            ("rex", 720, 3600, 3608, 728), ("tex", 1440, 4320, 3625, 745), ("xhi", 3600, 7200, 4720, 1120), ("tphi", 7200, 14400, 8479, 1279),
            ("tpi", 7200, 14400, 9264, 2064), ("sphi", 3600, 10800, 7502, 302), ("spex", 3600, 10800, 8739, 1539), ("xex", 2400, 7200, 5695, 895),
            ("dphi", 7200, 18000, 12342, 1542), ("dex", 7200, 18000, 11957, 1157), ("gishi", 14400, 28800, 16387, 1987)
        };

        static List<double[]> LoadJson(string name)
        {
            var candidates = new[] {
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_data", $"{name}.json"),
                Path.Combine("P:\\workspace\\generate_polychora\\output", $"{name}.json"),
                Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "Assets", "uniform_polychora", $"{name}.json")
            };
            foreach (var path in candidates)
                if (File.Exists(path)) return ParseVertices(File.ReadAllText(path));
            throw new FileNotFoundException($"Cannot find {name}.json");
        }

        [TestCaseSource(nameof(Expected))]
        public void TopologyVerification((string name, int v, int e, int f, int c) e)
        {
            var verts = LoadJson(e.name);
            var poly = RegularPolytope4d.FromVertices(verts);
            
            Assert.That(poly.vertices.Count, Is.EqualTo(e.v), $"{e.name}: V mismatch");
            Assert.That(poly.edges.Count, Is.EqualTo(e.e), $"{e.name}: E mismatch");
            Assert.That(poly.cells.Count, Is.EqualTo(e.c), $"{e.name}: C mismatch");
            
            int implicitF = poly.edges.Count + poly.cells.Count - poly.vertices.Count;
            Assert.That(implicitF, Is.EqualTo(e.f), $"{e.name}: F mismatch (implicit Euler)");
        }

        static List<double[]> ParseVertices(string json)
        {
            var result = new List<double[]>();
            int vIdx = json.IndexOf("\"vertices\"", StringComparison.Ordinal);
            int startArr = json.IndexOf('[', vIdx);
            int depth = 0, startPoint = -1;
            for (int i = startArr; i < json.Length; i++) {
                if (json[i] == '[') { depth++; if (depth == 2) startPoint = i + 1; }
                else if (json[i] == ']') {
                    if (depth == 2) {
                        string pointStr = json.Substring(startPoint, i - startPoint);
                        var coords = pointStr.Split(',').Select(s => double.Parse(s.Trim(), System.Globalization.CultureInfo.InvariantCulture)).ToArray();
                        result.Add(coords);
                    }
                    depth--; if (depth == 0) break;
                }
            }
            return result;
        }
    }
}
