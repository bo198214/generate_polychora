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
    // Known topology for all 6 regular 4D polytopes
    static readonly (string name, int verts, int edges, int faces2d, int cells3d)[] Expected = {
        ("5-cell",   5,   10,  10,   5),
        ("8-cell",  16,   32,  24,   8),
        ("16-cell",  8,   24,  32,  16),
        ("24-cell", 24,   96,  96,  24),
        ("120-cell", 600, 1200, 720, 120),
        ("600-cell", 120, 720, 1200, 600), // slow — enable for full verification
        ("bitruncated_5_cell",30,90,20,20)

    };

    static List<double[]> LoadJson(string name)
    {
        // Look relative to the test project; fall back to unity-projects path
        var candidates = new[] {
            Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_data", $"{name}.json"),
            Path.Combine("p:\\unity-projects\\tesserian-6.3\\Assets\\uniform_polychora", $"{name}.json")
        };
        foreach (var path in candidates)
            if (File.Exists(path))
            {
                var json = File.ReadAllText(path);
                return ParseVertices(json);
            }
        throw new FileNotFoundException($"Cannot find {name}.json");
    }

    [TestCaseSource(nameof(Expected))]
    public void VertexCount((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var v = LoadJson(e.name);
        Assert.That(v.Count, Is.EqualTo(e.verts), $"{e.name}: vertex count");
    }

    [TestCaseSource(nameof(Expected))]
    public void AllEdgesEqualLength((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var verts = LoadJson(e.name);
        var poly  = RegularPolytope4d.FromVertices(verts);
        foreach (var (a, b) in poly.edges)
        {
            double d = Dist(verts[a], verts[b]);
            Assert.That(d, Is.EqualTo(poly.edgeLength).Within(1e-4), $"{e.name}: edge ({a},{b}) length");
        }
    }

    [TestCaseSource(nameof(Expected))]
    public void EdgeCount((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var poly = RegularPolytope4d.FromVertices(LoadJson(e.name));
        Assert.That(poly.edges.Count, Is.EqualTo(e.edges), $"{e.name}: edge count");
    }

    [TestCaseSource(nameof(Expected))]
    public void CellCount((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var poly = RegularPolytope4d.FromVertices(LoadJson(e.name));
        Assert.That(poly.cells.Count, Is.EqualTo(e.cells3d), $"{e.name}: cell count");
    }

    [TestCaseSource(nameof(Expected))]
    public void CellNormalsAreUnitVectors((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var poly = RegularPolytope4d.FromVertices(LoadJson(e.name));
        foreach (var n in poly.cellNormals)
        {
            double len = Math.Sqrt(n.Sum(x => x * x));
            Assert.That(len, Is.EqualTo(1.0).Within(1e-4), $"{e.name}: normal not unit");
        }
    }

    [TestCaseSource(nameof(Expected))]
    public void BackFaceCullingHalvesVisibleCells((string name, int verts, int edges, int faces2d, int cells3d) e)
    {
        var poly = RegularPolytope4d.FromVertices(LoadJson(e.name));
        // Eye far along w-axis — should see roughly half the cells
        var eye  = new double[] { 0, 0, 0, 10 };
        var vis  = poly.VisibleCells(eye);
        Assert.That(vis.Count, Is.GreaterThan(0).And.LessThan(poly.cells.Count),
            $"{e.name}: back-face culling should remove some cells");
    }

    static List<double[]> ParseVertices(string json)
    {
        var result = new List<double[]>();
        int keyStart = json.IndexOf("\"vertices\"", StringComparison.Ordinal);
        if (keyStart < 0) throw new FormatException("No 'vertices' key found in JSON");
        int arrStart = json.IndexOf('[', keyStart);
        if (arrStart < 0) throw new FormatException("No 'vertices' array in JSON");
        int depth = 0, arrEnd = -1;
        for (int i = arrStart; i < json.Length; i++)
        {
            if (json[i] == '[') depth++;
            else if (json[i] == ']' && --depth == 0) { arrEnd = i; break; }
        }
        if (arrEnd < 0) throw new FormatException("Unmatched '[' in 'vertices' array");
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

    static double Dist(double[] a, double[] b)
    {
        double s = 0;
        for (int i = 0; i < 4; i++) s += (a[i]-b[i])*(a[i]-b[i]);
        return Math.Sqrt(s);
    }
}
}
