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
        // Expected vertex counts for Wythoffian Polychora (Bowers names)
        static readonly (string name, int verts)[] WythoffianExpected = {
            ("pen", 5), ("tip", 20), ("rap", 10), ("deca", 30), ("hap", 30), ("tap", 60), ("sadi", 20), ("dappat", 60), ("tappy", 120),
            ("tes", 16), ("tah", 64), ("rico", 32), ("hex", 8), ("tico", 48), ("spic", 96), ("thic", 192), ("xic", 64), ("scic", 192), ("gic", 384),
            ("ico", 24), ("cont", 192), ("tico_f", 96), ("srico", 288), ("frico", 288), ("grico", 576), ("prico", 144), ("drico", 576), ("trico", 1152),
            ("hi", 600), ("rhi", 1200), ("thi", 2400), ("rex", 720), ("sphi", 3600), ("xhi", 3600), ("tphi", 7200), ("ex", 120), ("xex", 2400), ("spex", 3600), ("dphi", 7200), ("tex", 1440), ("dex", 7200), ("tpi", 7200), ("gishi", 14400)
        };

        static List<double[]> LoadJson(string name)
        {
            var candidates = new[] {
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_data", $"{name}.json"),
                Path.Combine("P:\\workspace\\generate_polychora\\output", $"{name}.json")
            };
            foreach (var path in candidates)
                if (File.Exists(path))
                {
                    var json = File.ReadAllText(path);
                    return ParseVertices(json);
                }
            throw new FileNotFoundException($"Cannot find {name}.json");
        }

        [TestCaseSource(nameof(WythoffianExpected))]
        public void VertexCount((string name, int verts) e)
        {
            var v = LoadJson(e.name);
            Assert.That(v.Count, Is.EqualTo(e.verts), $"{e.name}: vertex count mismatch");
        }

        [Test]
        public void GeneratorProducesCorrectVertexCounts()
        {
            // Test a few via the internal generator to ensure logic works
            var pen = RegularPolytope4d.PolychoraGenerator.GenerateVertices("A4", 1);
            Assert.That(pen.Count, Is.EqualTo(5));

            var tes = RegularPolytope4d.PolychoraGenerator.GenerateVertices("B4", 1);
            Assert.That(tes.Count, Is.EqualTo(16));

            var hex = RegularPolytope4d.PolychoraGenerator.GenerateVertices("B4", 8);
            Assert.That(hex.Count, Is.EqualTo(8));
        }

        [TestCase("pen")]
        [TestCase("tes")]
        [TestCase("ico")]
        public void FromVertices_ComputesEdgesAndCells(string name)
        {
            var verts = LoadJson(name);
            var poly = RegularPolytope4d.FromVertices(verts);
            
            Assert.That(poly.edges.Count, Is.GreaterThan(0));
            Assert.That(poly.cells.Count, Is.GreaterThan(0));
            Assert.That(poly.cellNormals.Count, Is.EqualTo(poly.cells.Count));
        }

        static List<double[]> ParseVertices(string json)
        {
            var result = new List<double[]>();
            int keyStart = json.IndexOf("\"vertices\"", StringComparison.Ordinal);
            int arrStart = json.IndexOf('[', keyStart);
            int depth = 0, arrEnd = -1;
            for (int i = arrStart; i < json.Length; i++)
            {
                if (json[i] == '[') depth++;
                else if (json[i] == ']' && --depth == 0) { arrEnd = i; break; }
            }
            string outer = json.Substring(arrStart, arrEnd - arrStart + 1);
            foreach (Match row in Regex.Matches(outer, @"\[([^\[\]]+)\]"))
            {
                var nums = row.Groups[1].Value.Split(',')
                    .Select(s => double.Parse(s.Trim(), System.Globalization.CultureInfo.InvariantCulture))
                    .ToArray();
                result.Add(nums);
            }
            return result;
        }
    }
}
