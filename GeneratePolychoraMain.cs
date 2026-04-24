using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using D4BB.Geometry;

namespace D4BB.Geometry
{
    public class GeneratePolychoraMain
    {
        // Mapping: polytope name → (group, activeNodes) determined by systematic search
        // Each entry was verified: TrueConvexHull4D.Compute gives exactly the expected (V,E,F,C)
        static readonly (string name, int v, int e, int f, int c, string group, int an)[] KnownMappings = {
            // A4 family
            ("pen",     5,   10,   10,   5, "A4",  1),
            ("rap",    10,   30,   30,  10, "A4",  2),
            ("tip",    20,   40,   30,  10, "A4",  3),
            ("deca",   30,   60,   40,  10, "A4",  6),
            ("hap",    30,   90,   80,  20, "A4",  5),
            ("tap",    60,  120,   80,  20, "A4",  7),
            ("dappat", 60,  150,  120,  30, "A4", 11),
            ("tappy", 120,  240,  150,  30, "A4", 15),
            // B4 family
            ("tes",    16,   32,   24,   8, "B4",  1),
            ("hex",     8,   24,   32,  16, "B4",  8),
            ("spic",   96,  288,  240,  48, "B4", 10),
            // F4 family
            ("ico",    24,   96,   96,  24, "F4",  1),
            ("cont",  192,  384,  240,  48, "F4",  3),
            ("srico", 288,  576,  336,  48, "F4",  6),
            ("prico", 144,  576,  672, 240, "F4",  9),
            // H4 family
            ("hi",    600, 1200,  720, 120, "H4",  1),
            ("ex",    120,  720, 1200, 600, "H4",  8),
            ("rhi",  1200, 3600, 3120, 720, "H4",  2),
            ("rex",   720, 3600, 3600, 720, "H4",  4),
        };

        static readonly string OutputDir = Path.Combine(
            AppDomain.CurrentDomain.BaseDirectory, "..", "..", "..", "..", "output");

        public static void RegenerateKnown(bool verify = true)
        {
            Console.WriteLine($"Regenerating {KnownMappings.Length} known polytopes...\n");
            int ok = 0, fail = 0;

            foreach (var (name, ev, ee, ef, ec, group, an) in KnownMappings)
            {
                var pts = PolychoraGenerator.GenerateVertices(group, an);
                string status = "";

                if (verify && pts.Count <= 1500)
                {
                    var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                    var (av, ae, af, ac) = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                    if (av == ev && ae == ee && af == ef && ac == ec)
                    { status = "✓"; ok++; }
                    else
                    { status = $"✗ got ({av},{ae},{af},{ac}) expected ({ev},{ee},{ef},{ec})"; fail++; }
                }
                else
                { status = pts.Count > 1500 ? "(large, skipped verify)" : ""; ok++; }

                Console.WriteLine($"  {name,-10} {group}/{an,-3} V={pts.Count,5}  {status}");
                SaveJson(name, group, pts);
            }

            Console.WriteLine($"\nDone. {ok} OK, {fail} mismatches.");
        }

        static void SaveJson(string name, string group, List<double[]> vertices)
        {
            string path = Path.Combine(OutputDir, $"{name}.json");
            var obj = new { name = name, group = group, vertices = vertices };
            string json = JsonSerializer.Serialize(obj, new JsonSerializerOptions { WriteIndented = true });
            File.WriteAllText(path, json);
        }

        public static void Main(string[] args)
        {
            bool verifyOnly = args.Contains("--verify-only");
            bool skipVerify = args.Contains("--skip-verify");

            Console.WriteLine("=== PolychoraGenerator: Regenerating JSON files ===\n");
            RegenerateKnown(!skipVerify);

            if (!verifyOnly)
            {
                Console.WriteLine("\nNote: The following polytopes still need correct vertex data:");
                var known = new HashSet<string>(KnownMappings.Select(m => m.name));
                string[] allNames = {
                    "sadi","tah","rico","tico","thic","xic","scic","gic",
                    "tico_f","frico","grico","drico","trico",
                    "thi","tex","xhi","tphi","tpi","sphi","spex","xex","dphi","dex","gishi",
                    "gap","snic"
                };
                foreach (var n in allNames.Where(n => !known.Contains(n)))
                    Console.WriteLine($"  {n}");
            }
        }
    }
}
