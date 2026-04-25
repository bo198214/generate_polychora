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
            double eps = 1e-6,
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
                    var hull = TrueConvexHull4D.Compute(vertices, eps);

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
            double eps = 1e-6, Action<string> log = null)
        {
            Directory.CreateDirectory(topologyDir);
            var file = Path.Combine(vertexDir, $"{name}.json");
            if (!File.Exists(file)) { log?.Invoke($"  {name}: vertex file not found"); return; }
            var (vertices, description) = ReadVerticesAndDesc(file);
            var hull = TrueConvexHull4D.Compute(vertices, eps);
            var outPath = Path.Combine(topologyDir, $"{name}.json");
            WriteTopology(outPath, name, description, hull);
            log?.Invoke($"  {name}: V={hull.Vertices.Count} E={hull.Edges.Count} F={hull.Faces.Count} C={hull.Cells.Count}");
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
            // Write to a temp file first; rename atomically to avoid corrupt partials on abort.
            var tmp = outPath + ".tmp";
            using var stream = File.Create(tmp);
            using var writer = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true });

            writer.WriteStartObject();
            writer.WriteString("name", name);
            writer.WriteString("description", description);

            // vertices: array of 4D coordinate arrays
            writer.WritePropertyName("vertices");
            writer.WriteStartArray();
            foreach (var v in hull.Vertices)
            {
                writer.WriteStartArray();
                foreach (var x in v) writer.WriteNumberValue(Math.Round(x, 8));
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            // edges: array of [i, j] index pairs
            writer.WritePropertyName("edges");
            writer.WriteStartArray();
            foreach (var e in hull.Edges)
            {
                writer.WriteStartArray();
                foreach (var i in e) writer.WriteNumberValue(i);
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            // faces2d: array of vertex-index lists (2D faces / ridges)
            writer.WritePropertyName("faces2d");
            writer.WriteStartArray();
            foreach (var f in hull.Faces)
            {
                writer.WriteStartArray();
                foreach (var i in f) writer.WriteNumberValue(i);
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            // cells: array of vertex-index lists (3D cells)
            writer.WritePropertyName("cells");
            writer.WriteStartArray();
            foreach (var c in hull.Cells)
            {
                writer.WriteStartArray();
                foreach (var i in c) writer.WriteNumberValue(i);
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            // normals: outward unit normals, one 4D vector per cell
            writer.WritePropertyName("normals");
            writer.WriteStartArray();
            foreach (var n in hull.Normals)
            {
                writer.WriteStartArray();
                foreach (var x in n) writer.WriteNumberValue(Math.Round(x, 8));
                writer.WriteEndArray();
            }
            writer.WriteEndArray();

            writer.WriteEndObject();
            writer.Flush();
            stream.Close();

            // Atomic rename: only a completed file appears at the final path.
            if (File.Exists(outPath)) File.Delete(outPath);
            File.Move(tmp, outPath);
        }
    }
}
