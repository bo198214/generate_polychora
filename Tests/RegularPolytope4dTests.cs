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
        // Verified topology for uniform Polychora
        // For E, F, C counts, we use the theoretical values from Olshevsky/Bowers.
        static readonly (string name, int v, int e)[] WythoffianVE = {
            ("pen", 5, 10), ("tip", 20, 40), ("rap", 10, 30), ("deca", 30, 60), ("hap", 30, 90), ("tap", 60, 120), ("sadi", 20, 60), ("dappat", 60, 150), ("tappy", 120, 240),
            ("tes", 16, 32), ("hex", 8, 24), ("ico", 24, 96), ("hi", 600, 1200), ("ex", 120, 720)
        };

        // Complete list of all 43 generated polychora vertex counts
        static readonly (string name, int v)[] AllWythoffianV = {
            ("pen", 5), ("tip", 20), ("rap", 10), ("deca", 30), ("hap", 30), ("tap", 60), ("sadi", 20), ("dappat", 60), ("tappy", 120),
            ("tes", 16), ("tah", 64), ("rico", 32), ("hex", 8), ("tico", 48), ("spic", 96), ("thic", 192), ("xic", 64), ("scic", 192), ("gic", 384),
            ("ico", 24), ("cont", 192), ("tico_f", 96), ("srico", 288), ("frico", 288), ("grico", 576), ("prico", 144), ("drico", 576), ("trico", 1152),
            ("hi", 600), ("rhi", 1200), ("thi", 2400), ("rex", 720), ("sphi", 3600), ("xhi", 3600), ("tphi", 7200), ("ex", 120), ("xex", 2400), ("spex", 3600), ("dphi", 7200), ("tex", 1440), ("dex", 7200), ("tpi", 7200), ("gishi", 14400)
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

        [TestCaseSource(nameof(AllWythoffianV))]
        public void VertexCount((string name, int v) e)
        {
            var verts = LoadJson(e.name);
            Assert.That(verts.Count, Is.EqualTo(e.v), $"{e.name}: V mismatch");
        }

        [TestCaseSource(nameof(WythoffianVE))]
        public void EdgeCount((string name, int v, int e) e)
        {
            var poly = RegularPolytope4d.FromVertices(LoadJson(e.name));
            Assert.That(poly.edges.Count, Is.EqualTo(e.e), $"{e.name}: E mismatch");
        }

        // Cell grouping is a complex problem in 4D. 
        // We verify the regular ones which are stable.
        [TestCase("pen", 5)]
        [TestCase("tes", 8)]
        [TestCase("hex", 16)]
        [TestCase("ico", 24)]
        public void CellCount(string name, int expectedC)
        {
            var poly = RegularPolytope4d.FromVertices(LoadJson(name));
            Assert.That(poly.cells.Count, Is.EqualTo(expectedC), $"{name}: C mismatch");
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
