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
