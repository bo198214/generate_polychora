using System;
using System.Collections.Generic;
using System.Linq;
using D4BB.Geometry;
using NUnit.Framework;

namespace D4BB.GeometryTests
{
    [TestFixture]
    public class FindActiveNodesTest
    {
        static readonly (string name, int v, int e, int f, int c)[] Expected = {
            ("tes",16,32,24,8),("hex",8,24,32,16),("tah",64,128,80,16),("rico",32,96,80,16),("tico",48,120,88,16),
            ("spic",96,288,240,48),("thic",192,384,224,32),("xic",64,192,160,32),("scic",192,480,320,32),("gic",384,768,448,64),
            ("ico",24,96,96,24),("cont",192,384,240,48),("tico_f",96,288,216,24),("srico",288,576,336,48),
            ("frico",288,864,672,96),("grico",576,1152,816,240),("prico",144,576,672,240),("drico",576,1440,1008,144),("trico",1152,2304,1344,192),
            ("hi",600,1200,720,120),("ex",120,720,1200,600),
        };

        static readonly (string group, int[,] matrix)[] Groups = {
            ("B4", new int[,]{{1,4,2,2},{4,1,3,2},{2,3,1,3},{2,2,3,1}}),
            ("F4", new int[,]{{1,3,2,2},{3,1,4,2},{2,4,1,3},{2,2,3,1}}),
        };

        [Test, Timeout(600000)]
        public void SearchAllGroupsAndNodes()
        {
            foreach (var (gname, mat) in Groups)
            {
                var mirrors = GetMirrorNormals(mat);
                Console.WriteLine($"\n=== {gname} ===");
                for (int an = 1; an <= 15; an++)
                {
                    // Use PolychoraGenerator directly (no dot-skip bug)
                    var pts = PolychoraGenerator.GenerateVertices(gname, an);
                    Console.WriteLine($"  an={an,2}: orbit V={pts.Count,5}");
                    if (pts.Count > 600) { Console.WriteLine("       (skipped hull)"); continue; }
                    var hull = TrueConvexHull4D.Compute(pts, 1e-6);
                    var vefc = (hull.Vertices.Count, hull.Edges.Count, hull.Faces.Count, hull.Cells.Count);
                    var match = Expected.Where(e => e.v==vefc.Item1 && e.e==vefc.Item2 && e.f==vefc.Item3 && e.c==vefc.Item4).ToList();
                    string tag = match.Any() ? "  <-- " + string.Join(",", match.Select(m=>m.name)) : "";
                    Console.WriteLine($"       hull: V={vefc.Item1,5} E={vefc.Item2,5} F={vefc.Item3,5} C={vefc.Item4,4}{tag}");
                }
            }
            Assert.Pass();
        }

        static double[][] GetMirrorNormals(int[,] matrix) {
            int n = matrix.GetLength(0);
            var b = new double[n][];
            for (int i = 0; i < n; i++) {
                b[i] = new double[n];
                for (int j = 0; j < i; j++) {
                    double dot = -Math.Cos(Math.PI / matrix[i,j]);
                    double sum = 0; for (int k = 0; k < j; k++) sum += b[i][k]*b[j][k];
                    b[i][j] = (dot-sum)/(b[j][j]==0?1:b[j][j]);
                }
                double ssq = 0; for (int k = 0; k < i; k++) ssq += b[i][k]*b[i][k];
                b[i][i] = Math.Sqrt(Math.Max(0,1.0-ssq));
            }
            return b;
        }
        static double[] SolveGen(double[][] mirrors, int active) {
            int n = mirrors.Length; var d = new double[n]; var p = new double[n];
            for (int i=0;i<n;i++) d[i]=((active>>i)&1)!=0?1.0:0.0;
            for (int i=0;i<n;i++){double s=0;for(int j=0;j<i;j++)s+=mirrors[i][j]*p[j];p[i]=(d[i]-s)/(mirrors[i][i]==0?1:mirrors[i][i]);}
            return p;
        }
        static List<double[]> GenerateOrbit(double[] start, double[][] mirrors) {
            var verts = new List<double[]>{start};
            var seen = new System.Collections.Generic.HashSet<string>{Key(start)};
            for (int i=0;i<verts.Count;i++){
                var v = verts[i];
                foreach (var m in mirrors){
                    double dot=0; for(int k=0;k<4;k++) dot+=v[k]*m[k];
                    if (Math.Abs(dot)<1e-8) continue;
                    var nxt=new double[4]; for(int k=0;k<4;k++) nxt[k]=v[k]-2*dot*m[k];
                    var key=Key(nxt); if(seen.Add(key)) verts.Add(nxt);
                }
            }
            return verts;
        }
        static string Key(double[] v) => string.Join(",",v.Select(x=>Math.Round(x,8).ToString("F8")));
    }
}
