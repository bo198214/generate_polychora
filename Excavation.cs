using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;

namespace D4BB.Geometry
{
    /// <summary>
    /// Cell excavation — the elementary nonconvex CSG operation on a polychoron boundary:
    /// replace one cell by the lateral cells of a pyramid whose apex points INTO the solid
    /// (the mirror image of the CRF community's "augmentation"; in 3D this is Bonnie
    /// Stewart's excavation, the move that builds the Stewart toroids).  The result is an
    /// embedded, nonconvex polychoron whose faces are all still regular polygons.
    ///
    /// Requirements checked at runtime:
    ///  - the excavated cell must have simplicial (triangular) 2-faces, so all new faces
    ///    are triangles (no cyclic-ordering questions) and all new cells are pyramids with
    ///    unit edges: apex distance = edge length e requires circumradius R &lt; e, giving
    ///    depth h = sqrt(e² − R²);
    ///  - the apex must lie strictly inside every ORIGINAL cell hyperplane (otherwise the
    ///    dent would punch through the far side — this rules out e.g. the 16-cell, whose
    ///    pentachoron dent of depth 1.118e exceeds its 1.0e thickness);
    ///  - for multiple excavations the pyramids must not overlap (only checked implicitly:
    ///    each excavation re-validates against the current boundary).
    /// The closedness invariant (every face in exactly 2 cells) and Euler characteristic 0
    /// are re-validated after every operation.
    /// </summary>
    public static class Excavation
    {
        public class Topo
        {
            public string Name = "", Description = "";
            public List<double[]> Vertices = new();
            public List<int[]> Edges = new();
            public List<int[]> Faces = new();
            public List<int[]> Cells = new();
            public List<int[]> CellFaces = new();
            public List<double[]> Normals = new();
        }

        /// <summary>Generates the curated set of excavated polychora into <paramref name="outDir"/>.</summary>
        public static void GenerateCuratedSet(string topoDir, string outDir, Action<string> log)
        {
            Directory.CreateDirectory(outDir);

            // 600-cell with one pentachoron dimple.
            var ex1 = Load(Path.Combine(topoDir, "ex.json"));
            ExcavateCell(ex1, 0, log);
            Save(ex1, Path.Combine(outDir, "excavated-ex.json"), "excavated-ex",
                "600-cell with one tetrahedral cell excavated by a pentachoron (nonconvex, regular-faced)");
            log($"  excavated-ex: V={ex1.Vertices.Count} E={ex1.Edges.Count} F={ex1.Faces.Count} C={ex1.Cells.Count}");

            // 600-cell with two antipodal dimples.
            var ex2 = Load(Path.Combine(topoDir, "ex.json"));
            int first = 0;
            ExcavateCell(ex2, first, log);
            int anti = FindAntipodalCell(ex2, first);
            ExcavateCell(ex2, anti, log);
            Save(ex2, Path.Combine(outDir, "bi-excavated-ex.json"), "bi-excavated-ex",
                "600-cell with two antipodal tetrahedral cells excavated by pentachora (nonconvex, regular-faced)");
            log($"  bi-excavated-ex: V={ex2.Vertices.Count} E={ex2.Edges.Count} F={ex2.Faces.Count} C={ex2.Cells.Count}");

            // 24-cell with one octahedral-pyramid dent — the pyramid depth equals the
            // inradius, so the dent's apex sits exactly at the polytope's center.
            var ico = Load(Path.Combine(topoDir, "ico.json"));
            ExcavateCell(ico, 0, log);
            Save(ico, Path.Combine(outDir, "excavated-ico.json"), "excavated-ico",
                "24-cell with one octahedral cell excavated by an octahedral pyramid reaching the exact center (nonconvex, regular-faced)");
            log($"  excavated-ico: V={ico.Vertices.Count} E={ico.Edges.Count} F={ico.Faces.Count} C={ico.Cells.Count}");

            // Snub 24-cell with one icosahedral-pyramid dimple (shallow: h ≈ 0.31e).
            var sadi = Load(Path.Combine(topoDir, "sadi.json"));
            int ike = Enumerable.Range(0, sadi.Cells.Count).First(i => sadi.Cells[i].Length == 12);
            ExcavateCell(sadi, ike, log);
            Save(sadi, Path.Combine(outDir, "excavated-sadi.json"), "excavated-sadi",
                "Snub 24-cell with one icosahedral cell excavated by an icosahedral pyramid (nonconvex, regular-faced)");
            log($"  excavated-sadi: V={sadi.Vertices.Count} E={sadi.Edges.Count} F={sadi.Faces.Count} C={sadi.Cells.Count}");
        }

        public static void ExcavateCell(Topo t, int cellIdx, Action<string> log)
        {
            var cellVerts = t.Cells[cellIdx];
            var cellFaces = t.CellFaces[cellIdx];
            foreach (var fi in cellFaces)
                if (t.Faces[fi].Length != 3)
                    throw new InvalidOperationException(
                        $"Excavation requires simplicial cell faces; face {fi} has {t.Faces[fi].Length} vertices.");

            // Circumcenter/edge/circumradius of the (regular) cell.
            var center = new double[4];
            foreach (var v in cellVerts) for (int k = 0; k < 4; k++) center[k] += t.Vertices[v][k] / cellVerts.Length;
            double edge = double.MaxValue;
            for (int i = 0; i < cellVerts.Length; i++)
                for (int j = i + 1; j < cellVerts.Length; j++)
                    edge = Math.Min(edge, Dist(t.Vertices[cellVerts[i]], t.Vertices[cellVerts[j]]));
            double R = Dist(center, t.Vertices[cellVerts[0]]);
            double h2 = edge * edge - R * R;
            if (h2 <= 0)
                throw new InvalidOperationException(
                    $"Cell {cellIdx}: circumradius {R:f4} ≥ edge {edge:f4} — no unit-edged pyramid exists (cf. dodecahedron).");
            double h = Math.Sqrt(h2);

            var n0 = t.Normals[cellIdx];
            var apex = new double[4];
            for (int k = 0; k < 4; k++) apex[k] = center[k] - h * n0[k];

            // The apex must stay strictly inside every current cell hyperplane (punch-through guard).
            for (int j = 0; j < t.Cells.Count; j++)
            {
                if (j == cellIdx) continue;
                double d = Dot(t.Normals[j], t.Vertices[t.Cells[j][0]]);
                double margin = d - Dot(t.Normals[j], apex);
                if (margin < 1e-6 * edge)
                    throw new InvalidOperationException(
                        $"Cell {cellIdx}: dent apex violates hyperplane of cell {j} (margin {margin:e2}) — dent would punch through.");
            }

            int apexIdx = t.Vertices.Count;
            t.Vertices.Add(apex);

            // Cell edge set (from its triangular faces).
            var cellEdges = new HashSet<(int, int)>();
            foreach (var fi in cellFaces)
            {
                var f = t.Faces[fi];
                for (int i = 0; i < 3; i++)
                {
                    int a = f[i], b = f[(i + 1) % 3];
                    cellEdges.Add(a < b ? (a, b) : (b, a));
                }
            }

            foreach (var v in cellVerts) t.Edges.Add(new[] { v, apexIdx });

            var edgeToNewFace = new Dictionary<(int, int), int>();
            foreach (var (a, b) in cellEdges)
            {
                edgeToNewFace[(a, b)] = t.Faces.Count;
                t.Faces.Add(new[] { a, b, apexIdx });
            }

            // One lateral pyramid cell per original face of the excavated cell.
            var newCells = new List<int[]>();
            var newCellFaces = new List<int[]>();
            var newNormals = new List<double[]>();
            foreach (var fi in cellFaces)
            {
                var f = t.Faces[fi];
                var cv = f.Append(apexIdx).OrderBy(x => x).ToArray();
                var cf = new List<int> { fi };
                for (int i = 0; i < 3; i++)
                {
                    int a = f[i], b = f[(i + 1) % 3];
                    cf.Add(edgeToNewFace[a < b ? (a, b) : (b, a)]);
                }
                cf.Sort();

                // Outward normal of the dent cell = pointing into the excavated cavity,
                // i.e. towards the removed cell's circumcenter.
                var u1 = Sub(t.Vertices[f[1]], t.Vertices[f[0]]);
                var u2 = Sub(t.Vertices[f[2]], t.Vertices[f[0]]);
                var u3 = Sub(apex, t.Vertices[f[0]]);
                var nn = Normalize(Cross4(u1, u2, u3));
                var cellCentroid = new double[4];
                foreach (var v in cv) for (int k = 0; k < 4; k++) cellCentroid[k] += t.Vertices[v][k] / cv.Length;
                if (Dot(nn, Sub(center, cellCentroid)) < 0) for (int k = 0; k < 4; k++) nn[k] = -nn[k];

                newCells.Add(cv);
                newCellFaces.Add(cf.ToArray());
                newNormals.Add(nn);
            }

            t.Cells.RemoveAt(cellIdx);
            t.CellFaces.RemoveAt(cellIdx);
            t.Normals.RemoveAt(cellIdx);
            t.Cells.AddRange(newCells);
            t.CellFaces.AddRange(newCellFaces);
            t.Normals.AddRange(newNormals);

            Validate(t, edge);
        }

        static int FindAntipodalCell(Topo t, int _)
        {
            // After the first excavation the original cell 0 is gone; the antipode of its
            // circumcenter is found among the remaining cells by centroid opposition.
            // (Works because the 600-cell is centrally symmetric.)
            var apex = t.Vertices[^1];                       // last added vertex = first dent apex
            double best = double.MaxValue; int bestIdx = -1;
            for (int i = 0; i < t.Cells.Count; i++)
            {
                var c = new double[4];
                foreach (var v in t.Cells[i]) for (int k = 0; k < 4; k++) c[k] += t.Vertices[v][k] / t.Cells[i].Length;
                double s = 0; for (int k = 0; k < 4; k++) { double d = c[k] + apex[k] / Mag(apex) * Mag(c); s += d * d; }
                if (s < best) { best = s; bestIdx = i; }
            }
            return bestIdx;
        }

        static void Validate(Topo t, double edge)
        {
            var inc = new int[t.Faces.Count];
            foreach (var cf in t.CellFaces) foreach (var fi in cf) inc[fi]++;
            int dangling = inc.Count(x => x != 2);
            if (dangling > 0)
                throw new InvalidOperationException($"Excavation broke closedness: {dangling} faces not shared by exactly 2 cells.");
            int euler = t.Vertices.Count - t.Edges.Count + t.Faces.Count - t.Cells.Count;
            if (euler != 0)
                throw new InvalidOperationException($"Excavation broke Euler characteristic: {euler} != 0.");
            // All new cells sit flat on their hyperplane and all pyramid edges are unit.
            var apex = t.Vertices[^1];
            foreach (var e in t.Edges.Where(e => e[0] == t.Vertices.Count - 1 || e[1] == t.Vertices.Count - 1))
            {
                double d = Dist(t.Vertices[e[0]], t.Vertices[e[1]]);
                if (Math.Abs(d - edge) > 1e-6 * edge)
                    throw new InvalidOperationException($"Pyramid edge length {d:f6} deviates from cell edge {edge:f6}.");
            }
        }

        // ── JSON I/O (same schema as topology_output) ────────────────────────

        public static Topo Load(string path)
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            var r = doc.RootElement;
            var t = new Topo
            {
                Name = r.GetProperty("name").GetString() ?? "",
                Description = r.TryGetProperty("description", out var d) ? d.GetString() ?? "" : "",
            };
            foreach (var v in r.GetProperty("vertices").EnumerateArray())
                t.Vertices.Add(v.EnumerateArray().Select(x => x.GetDouble()).ToArray());
            foreach (var e in r.GetProperty("edges").EnumerateArray())
                t.Edges.Add(e.EnumerateArray().Select(x => x.GetInt32()).ToArray());
            foreach (var f in r.GetProperty("faces2d").EnumerateArray())
                t.Faces.Add(f.EnumerateArray().Select(x => x.GetInt32()).ToArray());
            foreach (var c in r.GetProperty("cells").EnumerateArray())
                t.Cells.Add(c.EnumerateArray().Select(x => x.GetInt32()).ToArray());
            foreach (var cf in r.GetProperty("cell_faces").EnumerateArray())
                t.CellFaces.Add(cf.EnumerateArray().Select(x => x.GetInt32()).ToArray());
            foreach (var n in r.GetProperty("normals").EnumerateArray())
                t.Normals.Add(n.EnumerateArray().Select(x => x.GetDouble()).ToArray());
            return t;
        }

        public static void Save(Topo t, string path, string name, string description)
        {
            using var stream = File.Create(path);
            using var w = new Utf8JsonWriter(stream, new JsonWriterOptions { Indented = true });
            w.WriteStartObject();
            w.WriteString("name", name);
            w.WriteString("description", description);
            void IntArrays(string key, List<int[]> rows)
            {
                w.WritePropertyName(key); w.WriteStartArray();
                foreach (var r in rows) { w.WriteStartArray(); foreach (var i in r) w.WriteNumberValue(i); w.WriteEndArray(); }
                w.WriteEndArray();
            }
            void DblArrays(string key, List<double[]> rows)
            {
                w.WritePropertyName(key); w.WriteStartArray();
                foreach (var r in rows) { w.WriteStartArray(); foreach (var x in r) w.WriteNumberValue(Math.Round(x, 8)); w.WriteEndArray(); }
                w.WriteEndArray();
            }
            DblArrays("vertices", t.Vertices);
            IntArrays("edges", t.Edges);
            IntArrays("faces2d", t.Faces);
            IntArrays("cells", t.Cells);
            IntArrays("cell_faces", t.CellFaces);
            DblArrays("normals", t.Normals);
            w.WriteEndObject();
        }

        static double Dot(double[] a, double[] b) { double s = 0; for (int i = 0; i < 4; i++) s += a[i] * b[i]; return s; }
        static double[] Sub(double[] a, double[] b) { var r = new double[4]; for (int i = 0; i < 4; i++) r[i] = a[i] - b[i]; return r; }
        static double Mag(double[] a) => Math.Sqrt(Dot(a, a));
        static double Dist(double[] a, double[] b) => Mag(Sub(a, b));
        static double[] Normalize(double[] a) { double m = Mag(a); var r = new double[4]; for (int i = 0; i < 4; i++) r[i] = a[i] / m; return r; }
        static double[] Cross4(double[] u, double[] v, double[] w)
        {
            var r = new double[4];
            r[0] =  (u[1]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[1]*w[3]-v[3]*w[1]) + u[3]*(v[1]*w[2]-v[2]*w[1]));
            r[1] = -(u[0]*(v[2]*w[3]-v[3]*w[2]) - u[2]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[2]-v[2]*w[0]));
            r[2] =  (u[0]*(v[1]*w[3]-v[3]*w[1]) - u[1]*(v[0]*w[3]-v[3]*w[0]) + u[3]*(v[0]*w[1]-v[1]*w[0]));
            r[3] = -(u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0]));
            return r;
        }
    }
}
