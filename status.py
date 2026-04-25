import os
import json

vdir = 'vertex_output'
tdir = 'topology_output'
v=set(os.listdir(vdir))
t=set(os.listdir(tdir))
print(f"vertex_output: {len(v)} files")
print(f"topology_output: {len(t)} files")

m=sorted(v-t)
for fname in m:
    if not fname.endswith(".json"):
        continue
    path = os.path.join(vdir, fname)
    name = fname[:-5]

    with open(path, encoding="utf-8-sig") as f:
        d = json.load(f)

    v  = len(d.get("vertices",  []))
    desc = d.get("description", "")

    print(f"{name:10}  V={v:5} {desc[:50]}")


