using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using D4BB.Geometry;

// Usage:
//   dotnet run -- vertices          → regenerate all vertex JSONs in vertex_output/
//   dotnet run -- topology          → compute convex hull topology for all vertex JSONs
//   dotnet run -- all               → both steps in sequence

// When run via "dotnet run", the working directory is the project directory.
// AppContext.BaseDirectory points to bin/Debug/net8.0/ — we go 3 levels up to find the project.
var root = Directory.GetCurrentDirectory();
if (!Directory.Exists(Path.Combine(root, "vertex_output")))
{
    root = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", ".."));
    if (!Directory.Exists(Path.Combine(root, "vertex_output")))
        root = AppContext.BaseDirectory;
}

var vertexDir   = Path.Combine(root, "vertex_output");
var topologyDir = Path.Combine(root, "topology_output");

string cmd = args.Length > 0 ? args[0] : "help";
switch (cmd)
{
    case "vertices":
        RunVertexGeneration(vertexDir);
        break;
    case "topology":
        RunTopologyGeneration(vertexDir, topologyDir);
        break;
    case "all":
        RunVertexGeneration(vertexDir);
        RunTopologyGeneration(vertexDir, topologyDir);
        break;
    default:
        Console.WriteLine("Usage: dotnet run -- <vertices|topology|all>");
        Console.WriteLine("  vertices  – regenerate vertex JSONs in vertex_output/");
        Console.WriteLine("  topology  – compute topology JSONs in topology_output/");
        Console.WriteLine("  all       – both steps");
        break;
}

static void RunVertexGeneration(string vertexDir)
{
    Console.WriteLine($"Generating vertex files → {vertexDir}");
    Directory.CreateDirectory(vertexDir);

    // (group, activeNodes) → polytope name, verified by TrueConvexHull4D
    var mappings = new (string name, string grp, int an)[] {
        ("pen","A4",1), ("rap","A4",2), ("tip","A4",3), ("deca","A4",6),
        ("hap","A4",5), ("tap","A4",7), ("dappat","A4",11), ("tappy","A4",15),
        ("tes","B4",1), ("hex","B4",8), ("rico","B4",2), ("tah","B4",3),
        ("tico","B4",12), ("spic","B4",10), ("thic","B4",7), ("xic","B4",9),
        ("scic","B4",11), ("gic","B4",15),
        ("ico","F4",1), ("cont","F4",3), ("tico_f","B4",6), ("srico","F4",6),
        ("frico","F4",5), ("grico","F4",7), ("prico","F4",9),
        ("drico","F4",11), ("trico","F4",15),
        ("hi","H4",1), ("ex","H4",8), ("rhi","H4",2), ("tex","H4",12), ("rex","H4",4),
        // Missing B4 forms
        ("srit","B4",5), ("prit","B4",13),
    };

    int ok = 0;
    foreach (var (name, grp, an) in mappings)
    {
        var pts = PolychoraGenerator.GenerateVertices(grp, an);
        var json = System.Text.Json.JsonSerializer.Serialize(
            new { name, source = $"{grp}/{an}", vertices = pts },
            new System.Text.Json.JsonSerializerOptions { WriteIndented = true });
        File.WriteAllText(Path.Combine(vertexDir, $"{name}.json"), json);
        Console.WriteLine($"  {name,-10} {grp}/{an}: {pts.Count} vertices");
        ok++;
    }

    // Non-Wythoffian polytopes
    WriteNonWythoffian(vertexDir, "snic", "snub-24-cell", SnubGenerator.SnubIcositetrachoron());
    WriteNonWythoffian(vertexDir, "gap",  "grand-antiprism", SnubGenerator.GrandAntiprism());
    ok += 2;

    Console.WriteLine($"Done: {ok} vertex files written.\n");
}

static void WriteNonWythoffian(string dir, string name, string source, List<double[]> pts)
{
    var json = System.Text.Json.JsonSerializer.Serialize(
        new { name, source, vertices = pts },
        new System.Text.Json.JsonSerializerOptions { WriteIndented = true });
    File.WriteAllText(Path.Combine(dir, $"{name}.json"), json);
    Console.WriteLine($"  {name,-10} {source}: {pts.Count} vertices");
}

static void RunTopologyGeneration(string vertexDir, string topologyDir)
{
    Console.WriteLine($"Computing topology → {topologyDir}");
    TopologyGenerator.Generate(vertexDir, topologyDir, eps: 1e-6,
        log: msg => Console.WriteLine(msg));
    Console.WriteLine("Done.\n");
}
