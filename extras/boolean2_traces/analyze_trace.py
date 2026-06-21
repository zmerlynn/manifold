#!/usr/bin/env python3
# Summarize the trace and locate the open walk: cluster filtered_output segment
# endpoints by eps-proximity and report any cluster with in-degree != out-degree
# (the vertex where the boundary walk cannot close). That centroid is the ideal
# zoom center.
import json
import sys

trace = json.load(open(sys.argv[1]))
eps = trace["eps"]
print(f"eps = {eps:.6e}   rule = {trace.get('rule')}")
print(f"{len(trace['phases'])} phases:")
for i, ph in enumerate(trace["phases"]):
    print(f"  {i}: {ph['name']:<24} pts={len(ph.get('points',[]))} "
          f"segs={len(ph.get('segments',[]))} polys={len(ph.get('polygons',[]))}")


def cluster_endpoints(segs, tol):
    pts = []
    for s in segs:
        pts.append(tuple(s["a"]))
        pts.append(tuple(s["b"]))
    uniq = []
    for p in pts:
        for c in uniq:
            if abs(p[0] - c[0]) <= tol and abs(p[1] - c[1]) <= tol:
                break
        else:
            uniq.append(p)
    # net degree per cluster: +1 for tail(a), -1 for head(b)
    deg = [0] * len(uniq)
    cin = [0] * len(uniq)
    cout = [0] * len(uniq)
    def find(p):
        for j, c in enumerate(uniq):
            if abs(p[0] - c[0]) <= tol and abs(p[1] - c[1]) <= tol:
                return j
        return -1
    for s in segs:
        cout[find(tuple(s["a"]))] += 1
        cin[find(tuple(s["b"]))] += 1
    return uniq, cin, cout


def poly_area(verts):
    a = 0.0
    n = len(verts)
    for i in range(n):
        x0, y0 = verts[i]
        x1, y1 = verts[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return a / 2.0

for ph in trace["phases"]:
    if ph["name"] == "final_polygons":
        for k, poly in enumerate(ph.get("polygons", [])):
            print(f"\nfinal_polygons[{k}]: {len(poly['verts'])} verts, "
                  f"area={poly_area(poly['verts']):.4f}")

for name in ("filtered_output_edges",):
    ph = next((p for p in trace["phases"] if p["name"] == name), None)
    if not ph or not ph.get("segments"):
        continue
    print(f"\n=== {name}: degree balance (tol=eps) ===")
    uniq, cin, cout = cluster_endpoints(ph["segments"], eps)
    imbalanced = []
    for j, c in enumerate(uniq):
        if cin[j] != cout[j]:
            imbalanced.append(c)
            print(f"  IMBALANCED v=({c[0]:.12f},{c[1]:.12f}) in={cin[j]} out={cout[j]}")
    if imbalanced:
        cx = sum(c[0] for c in imbalanced) / len(imbalanced)
        cy = sum(c[1] for c in imbalanced) / len(imbalanced)
        print(f"  --> zoom center ({cx:.12f},{cy:.12f})  ({len(imbalanced)} imbalanced)")
        # pairwise separations among imbalanced verts (in eps units)
        for a in range(len(imbalanced)):
            for b in range(a + 1, len(imbalanced)):
                dx = imbalanced[a][0] - imbalanced[b][0]
                dy = imbalanced[a][1] - imbalanced[b][1]
                d = (dx * dx + dy * dy) ** 0.5
                print(f"     sep[{a},{b}] = {d:.4e} = {d/eps:.3f} eps")
    else:
        print("  (no imbalance found at tol=eps)")
