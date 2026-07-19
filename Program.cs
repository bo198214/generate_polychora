using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text.Json;
using D4BB.Geometry;

// Always use en-US so decimal points are '.' regardless of system locale.
CultureInfo.DefaultThreadCurrentCulture   = CultureInfo.GetCultureInfo("en-US");
CultureInfo.DefaultThreadCurrentUICulture = CultureInfo.GetCultureInfo("en-US");

// Replace Console.Out with an auto-flushing writer so progress dots appear
// immediately even when stdout is piped through PowerShell.
Console.SetOut(new StreamWriter(Console.OpenStandardOutput(), leaveOpen: true) { AutoFlush = true });

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
    case "excavate":
        // Nonconvex regular-faced polychora via cell excavation (see Excavation.cs)
        Excavation.GenerateCuratedSet(topologyDir, Path.Combine(root, "crf_output"),
            msg => Console.WriteLine(msg));
        break;
    case "topology-one" when args.Length >= 2:
        // Generate topology for a single polytope by name (used by Makefile)
        TopologyGenerator.GenerateOne(vertexDir, topologyDir, args[1], eps: 0,  // 0 = auto-scale with circumradius
            log: msg => Console.WriteLine(msg));
        break;
    default:
        Console.WriteLine("Usage: dotnet run -- <vertices|topology|all|topology-one NAME>");
        Console.WriteLine("  vertices         – regenerate all vertex JSONs");
        Console.WriteLine("  topology         – compute topology for missing files");
        Console.WriteLine("  all              – both steps");
        Console.WriteLine("  topology-one N   – (re)compute topology for one polytope");
        Console.WriteLine("  excavate         – generate curated nonconvex excavated polychora → crf_output/");
        break;
}

void RunVertexGeneration(string vertexDir)
{
    Console.WriteLine($"Generating vertex files → {vertexDir}");
    Directory.CreateDirectory(vertexDir);

    // (name, group, activeNodes, description)
    // Names are the Bowers acronyms as listed on
    // https://en.wikipedia.org/wiki/Uniform_4-polytope (aligned 2026-07-19).
    var mappings = new (string name, string grp, int an, string desc)[] {
        // A4 family (5-cell / pentachoron symmetry)
        ("pen",    "A4", 1,  "Pentachoron (5-cell)"),
        ("rap",    "A4", 2,  "Rectified pentachoron"),
        ("tip",    "A4", 3,  "Truncated pentachoron"),
        ("deca",   "A4", 6,  "Decachoron (bitruncated 5-cell)"),
        ("srip",   "A4", 5,  "Cantellated 5-cell (small rhombated pentachoron)"),
        ("grip",   "A4", 7,  "Cantitruncated 5-cell (great rhombated pentachoron)"),
        ("spid",   "A4", 9,  "Runcinated 5-cell (spid)"),
        ("prip",   "A4", 11, "Runcitruncated 5-cell (prismatorhombated pentachoron)"),
        ("gippid", "A4", 15, "Omnitruncated 5-cell (great prismatodecachoron)"),
        // B4 family (tesseract / 16-cell symmetry)
        ("tes",    "B4", 1,  "Tesseract (8-cell)"),
        ("hex",    "B4", 8,  "Hexadecachoron (16-cell)"),
        ("rit",    "B4", 2,  "Rectified tesseract"),
        ("tat",    "B4", 3,  "Truncated tesseract"),
        ("thex",   "B4", 12, "Truncated 16-cell (truncated hexadecachoron)"),
        ("rico",   "B4", 10, "Rectified 24-cell (rectified icositetrachoron)"),
        ("grit",   "B4", 7,  "Cantitruncated tesseract (great rhombated tesseract)"),
        ("sidpith","B4", 9,  "Runcinated tesseract (small disprismatotesseractihexadecachoron)"),
        ("proh",   "B4", 11, "Runcitruncated tesseract (prismatorhombated hexadecachoron)"),
        ("gidpith","B4", 15, "Omnitruncated tesseract (great disprismatotesseractihexadecachoron)"),
        ("srit",   "B4", 5,  "Cantellated tesseract (small rhombated tesseract)"),
        ("prit",   "B4", 13, "Runcitruncated 16-cell (prismatorhombated tesseract)"),
        ("tah",    "B4", 6,  "Bitruncated tesseract"),
        // F4 family (24-cell / icositetrachoron symmetry)
        ("ico",    "F4", 1,  "Icositetrachoron (24-cell)"),
        ("tico",   "F4", 3,  "Truncated 24-cell (truncated icositetrachoron)"),
        ("cont",   "F4", 6,  "Bitruncated 24-cell (tetracontoctachoron)"),
        ("srico",  "F4", 5,  "Cantellated 24-cell (small rhombated icositetrachoron)"),
        ("grico",  "F4", 7,  "Cantitruncated 24-cell"),
        ("spic",   "F4", 9,  "Runcinated 24-cell (small prismatotetracontoctachoron)"),
        ("prico",  "F4", 11, "Runcitruncated 24-cell (prismatorhombated icositetrachoron)"),
        ("gippic", "F4", 15, "Omnitruncated 24-cell (great prismatotetracontoctachoron)"),
        // H4 family (120-cell / 600-cell symmetry)
        ("hi",     "H4", 1,  "Hecatonicosachoron (120-cell)"),
        ("ex",     "H4", 8,  "Hexacosichoron (600-cell)"),
        ("rahi",   "H4", 2,  "Rectified 120-cell (rectified hecatonicosachoron)"),
        ("thi",    "H4", 3,  "Truncated 120-cell (thi)"),
        ("tex",    "H4", 12, "Truncated 600-cell"),
        ("rox",    "H4", 4,  "Rectified 600-cell (rectified hexacosichoron)"),
        ("sidpixhi","H4", 9,  "Runcinated 120-cell (sidpixhi)"),
        ("srahi",  "H4", 5,  "Cantellated 120-cell (srahi)"),
        ("xhi",    "H4", 6,  "Bitruncated 120-cell (xhi)"),
        ("srix",   "H4", 10, "Cantellated 600-cell (srix)"),
        ("grahi",  "H4", 7,  "Cantitruncated 120-cell (grahi)"),
        ("prix",   "H4", 11, "Runcitruncated 120-cell (prix)"),
        ("prahi",  "H4", 13, "Runcitruncated 600-cell (prismatorhombated hecatonicosachoron)"),
        ("grix",   "H4", 14, "Cantitruncated 600-cell (grix)"),
        ("gidpixhi","H4",15, "Great disprismatohexacosihecatonicosachoron (gidpixhi)"),
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
    WriteNonWythoffian(vertexDir, "sadi", "Snub 24-cell (snub disicositetrachoron)", "snub-24-cell",    SnubGenerator.SnubIcositetrachoron());
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
    TopologyGenerator.Generate(vertexDir, topologyDir, eps: 0,  // 0 = auto-scale with circumradius
        log: msg => Console.WriteLine(msg));
    Console.WriteLine("Done.\n");
}
