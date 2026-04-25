using System;
using System.IO;
using System.Linq;
using System.Text.Json;
using D4BB.Geometry;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    /// <summary>
    /// Extended search tests for finding correct (group, activeNodes) mappings.
    /// Generation of vertex and topology files is handled by Program.cs (dotnet run).
    /// </summary>
    [TestFixture]
    public class ExtendedSearchTest
    {
        static readonly string OutDir = Path.GetFullPath(Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "vertex_output"));

        [Test, CancelAfter(300000)]
        public void FindF4LargeActiveNodes()
        {
            var expected = new (string name, int v, int e, int f, int c)[] {
                ("frico",288,864,720,144), ("grico",576,1152,720,144),
                ("drico",576,1440,1104,240), ("trico",1152,2304,1392,240),
            };
            var byVEFC = expected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);
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
            }
            Assert.Pass();
        }

        [Test, CancelAfter(600000)]
        public void FindAndRegenerateH4()
        {
            var expected = new (string name, int v, int e, int f, int c)[] {
                ("tex",1440,4320,3600,720),
            };
            var byVEFC = expected.ToDictionary(x => (x.v, x.e, x.f, x.c), x => x.name);
            var h4Candidates = new[] { 6, 9, 10, 11, 12, 13 };
            Console.WriteLine("=== H4 activeNodes ===");
            foreach (int an in h4Candidates) {
                var pts = PolychoraGenerator.GenerateVertices("H4", an);
                Console.WriteLine($"H4/an={an}: V={pts.Count}");
                if (pts.Count > 2000) { Console.WriteLine("  (too large, skipping)"); continue; }
                var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                var key = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                string match = byVEFC.TryGetValue(key, out var n) ? $" <-- {n} ✓" : "";
                Console.WriteLine($"  hull ({key.Item1},{key.Item2},{key.Item3},{key.Item4}){match}");
            }
            Assert.Pass();
        }
    }
}
