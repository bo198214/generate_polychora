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

void RunVertexGeneration(string vertexDir)
{
    Console.WriteLine($"Generating vertex files → {vertexDir}");
    Directory.CreateDirectory(vertexDir);

    // (name, group, activeNodes, description)
    var mappings = new (string name, string grp, int an, string desc)[] {
        // A4 family (5-cell / pentachoron symmetry)
        ("pen",    "A4", 1,  "Pentachoron (5-cell)"),
        ("rap",    "A4", 2,  "Rectified pentachoron"),
        ("tip",    "A4", 3,  "Truncated pentachoron"),
        ("deca",   "A4", 6,  "Decachoron (bitruncated 5-cell)"),
        ("hap",    "A4", 5,  "Small rhombated pentachoron"),
        ("tap",    "A4", 7,  "Great rhombated pentachoron"),
        ("dappat", "A4", 11, "Prismatorhombated pentachoron"),
        ("tappy",  "A4", 15, "Great prismatodecachoron (omnitruncated 5-cell)"),
        // B4 family (tesseract / 16-cell symmetry)
        ("tes",    "B4", 1,  "Tesseract (8-cell)"),
        ("hex",    "B4", 8,  "Hexadecachoron (16-cell)"),
        ("rico",   "B4", 2,  "Rectified tesseract"),
        ("tah",    "B4", 3,  "Truncated tesseract"),
        ("tico",   "B4", 12, "Truncated hexadecachoron (truncated 16-cell)"),
        ("spic",   "B4", 10, "Rectified icositetrachoron (rectified 24-cell, B4 form)"),
        ("thic",   "B4", 7,  "Cantitruncated tesseract"),
        ("xic",    "B4", 9,  "Runcinated tesseract"),
        ("scic",   "B4", 11, "Runcitruncated tesseract"),
        ("gic",    "B4", 15, "Omnitruncated tesseract"),
        ("srit",   "B4", 5,  "Cantellated tesseract (small rhombated tesseract)"),
        ("prit",   "B4", 13, "Runcitruncated hexadecachoron"),
        // F4 family (24-cell / icositetrachoron symmetry)
        ("ico",    "F4", 1,  "Icositetrachoron (24-cell)"),
        ("cont",   "F4", 3,  "Truncated icositetrachoron (truncated 24-cell)"),
        ("tico_f", "B4", 6,  "Bitruncated tesseract"),
        ("srico",  "F4", 6,  "Small rhombated icositetrachoron"),
        ("frico",  "F4", 5,  "Cantellated 24-cell"),
        ("grico",  "F4", 7,  "Cantitruncated 24-cell"),
        ("prico",  "F4", 9,  "Prismatorhombated icositetrachoron"),
        ("drico",  "F4", 11, "Runcitruncated icositetrachoron"),
        ("trico",  "F4", 15, "Omnitruncated icositetrachoron"),
        // H4 family (120-cell / 600-cell symmetry)
        ("hi",     "H4", 1,  "Hecatonicosachoron (120-cell)"),
        ("ex",     "H4", 8,  "Hexacosichoron (600-cell)"),
        ("rhi",    "H4", 2,  "Rectified 120-cell"),
        ("tex",    "H4", 12, "Truncated 600-cell"),
        ("rex",    "H4", 4,  "Rectified 600-cell"),
    };

    int ok = 0;
    var opts = new System.Text.Json.JsonSerializerOptions { WriteIndented = true };
    foreach (var (name, grp, an, desc) in mappings)
    {
        var pts = PolychoraGenerator.GenerateVertices(grp, an);
        var json = System.Text.Json.JsonSerializer.Serialize(
            new { name, description = desc, source = $"{grp}/{an}", vertices = pts }, opts);
        File.WriteAllText(Path.Combine(vertexDir, $"{name}.json"), json);
        Console.WriteLine($"  {name,-10} {grp}/{an}: {pts.Count} vertices  ({desc})");
        ok++;
    }

    // Non-Wythoffian polytopes
    WriteNonWythoffian(vertexDir, "snic", "Snub icositetrachoron (snub 24-cell)", "snub-24-cell",    SnubGenerator.SnubIcositetrachoron());
    WriteNonWythoffian(vertexDir, "gap",  "Grand antiprism",                      "grand-antiprism", SnubGenerator.GrandAntiprism());
    ok += 2;

    Console.WriteLine($"Done: {ok} vertex files written.\n");
}

void WriteNonWythoffian(string dir, string name, string description, string source, List<double[]> pts)
{
    var json = System.Text.Json.JsonSerializer.Serialize(
        new { name, description, source, vertices = pts },
        new System.Text.Json.JsonSerializerOptions { WriteIndented = true });
    File.WriteAllText(Path.Combine(dir, $"{name}.json"), json);
    Console.WriteLine($"  {name,-10} {source}: {pts.Count} vertices  ({description})");
}

void RunTopologyGeneration(string vertexDir, string topologyDir)
{
    Console.WriteLine($"Computing topology → {topologyDir}");
    TopologyGenerator.Generate(vertexDir, topologyDir, eps: 1e-6,
        log: msg => Console.WriteLine(msg));
    Console.WriteLine("Done.\n");
}
