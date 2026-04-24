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
    }
}
