using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace D4BB.Geometry
{
    /// <summary>
    /// Reads vertex data from vertex_output/, computes convex hull topology,
    /// and writes the result to topology_output/ as JSON with fields:
    /// vertices, edges, faces2d, cells, normals.
    /// </summary>
    public static class TopologyGenerator
    {
        public static void Generate(
            string vertexDir,
            string topologyDir,
            double eps = 0,      // 0 = auto-scale with circumradius
            Action<string> log = null)
        {
            Directory.CreateDirectory(topologyDir);
            // Sort by vertex count (ascending) so small polytopes are processed first.
            var files = Directory.GetFiles(vertexDir, "*.json")
                .Select(f => (path: f, count: ReadVertexCount(f)))
                .OrderBy(x => x.count)
                .Select(x => x.path)
                .ToArray();
            log?.Invoke($"Generating topology for {files.Length} polytopes...");

            foreach (var file in files)
            {
                var name = Path.GetFileNameWithoutExtension(file);
                var outPath = Path.Combine(topologyDir, $"{name}.json");
                if (File.Exists(outPath)) { log?.Invoke($"  {name}: already exists, skipping"); continue; }
                try
                {
                    var (vertices, description) = ReadVerticesAndDesc(file);
                    double e = eps <= 0 ? AutoEps(vertices) : eps;
                    var hull = TrueConvexHull4D.Compute(vertices, e);

                    WriteTopology(outPath, name, description, hull);
                    log?.Invoke($"  {name}: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
                }
                catch (Exception ex)
                {
                    log?.Invoke($"  {name}: ERROR – {ex.Message}");
                }
            }
        }

        public static void GenerateOne(
            string vertexDir, string topologyDir, string name,
            double eps = 0, Action<string> log = null)
        {
            Directory.CreateDirectory(topologyDir);
            var file = Path.Combine(vertexDir, $"{name}.json");
            if (!File.Exists(file)) { log?.Invoke($"  {name}: vertex file not found"); return; }
            var (vertices, description) = ReadVerticesAndDesc(file);
            if (eps <= 0) eps = AutoEps(vertices);

            log?.Invoke($"  {name}: V={vertices.Count} eps={eps:e1} — initialising...");
            Console.Write("    ");
            var hull = TrueConvexHull4D.Compute(vertices, eps,
                onCellDiscovered: n => {
                    if (n == 1) { Console.WriteLine("done"); Console.Write("    "); return; }
                    Console.Write('#');
                    if (n % 100 == 0) { Console.WriteLine($"  {n}"); Console.Write("    "); }
                });
            Console.WriteLine();
            var outPath = Path.Combine(topologyDir, $"{name}.json");
            WriteTopology(outPath, name, description, hull);
            log?.Invoke($"  {name}: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
        }

        /// <summary>
        /// Scales the base eps (1e-6) by the mean circumradius of the vertices,
        /// so the tolerance is appropriate regardless of the polytope's coordinate scale.
        /// </summary>
        static double AutoEps(List<double[]> verts, double baseEps = 1e-6)
        {
            double sum = 0;
            foreach (var v in verts) { double r = 0; foreach (var x in v) r += x*x; sum += Math.Sqrt(r); }
            double radius = sum / verts.Count;
            return baseEps * Math.Max(1.0, radius);
        }

        static int ReadVertexCount(string path)
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            return doc.RootElement.GetProperty("vertices").GetArrayLength();
        }

        static (List<double[]> vertices, string description) ReadVerticesAndDesc(string path)
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            var root = doc.RootElement;
            var desc = root.TryGetProperty("description", out var d) ? d.GetString() ?? "" : "";
            var verts = root.GetProperty("vertices")
                .EnumerateArray()
                .Select(v => v.EnumerateArray().Select(x => x.GetDouble()).ToArray())
                .ToList();
            return (verts, desc);
        }

        static void WriteTopology(string outPath, string name, string description, TrueConvexHull4D hull)
        {
            using var stream = File.Create(outPath);
            using var writer = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true });

            writer.WriteStartObject();
            writer.WriteString("name", name);
            writer.WriteString("description", description);

            writer.WritePropertyName("vertices");
            writer.WriteStartArray();
            foreach (var v in hull.Vertices) { writer.WriteStartArray(); foreach (var x in v) writer.WriteNumberValue(Math.Round(x, 8)); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WritePropertyName("edges");
            writer.WriteStartArray();
            foreach (var e in hull.Edges) { writer.WriteStartArray(); foreach (var i in e) writer.WriteNumberValue(i); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WritePropertyName("faces2d");
            writer.WriteStartArray();
            foreach (var f in hull.Faces) { writer.WriteStartArray(); foreach (var i in f) writer.WriteNumberValue(i); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WritePropertyName("cells");
            writer.WriteStartArray();
            foreach (var c in hull.Cells) { writer.WriteStartArray(); foreach (var i in c) writer.WriteNumberValue(i); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WritePropertyName("cell_faces");
            writer.WriteStartArray();
            foreach (var cf in hull.CellFaces) { writer.WriteStartArray(); foreach (var i in cf) writer.WriteNumberValue(i); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WritePropertyName("normals");
            writer.WriteStartArray();
            foreach (var n in hull.Normals) { writer.WriteStartArray(); foreach (var x in n) writer.WriteNumberValue(Math.Round(x, 8)); writer.WriteEndArray(); }
            writer.WriteEndArray();

            writer.WriteEndObject();
        }
    }
}
