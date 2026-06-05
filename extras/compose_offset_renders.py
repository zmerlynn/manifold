#!/usr/bin/env python3
# Copyright 2026 The Manifold Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Diagnostic only (not shipped). Merge clipper2 + boolean2 scenario JSON
# (from render_offset_scenarios) into one self-contained HTML: 3 columns per
# scenario (Clipper2 fill | vertex diff | Boolean2 fill), inline SVG, no deps.
# The two backends agree on the offset *shape* on well-conditioned inputs, so
# the differences are vertex-level; the middle column draws each output vertex
# (C2 hollow blue ring, B2 filled orange dot) to expose count/placement diffs.
#
# Usage: compose_offset_renders.py clipper2.json boolean2.json out.html
import json
import sys

c2 = json.load(open(sys.argv[1]))
b2 = json.load(open(sys.argv[2]))
out_path = sys.argv[3]

PAD = 0.08          # fraction of bbox to pad
BLUE = "#2c7fb8"    # clipper2
ORANGE = "#e6550d"  # boolean2
GRAY = "#999"       # input outline


def bounds(rings_lists):
    xs, ys = [], []
    for rings in rings_lists:
        for r in rings:
            for x, y in r:
                xs.append(x)
                ys.append(-y)  # flip y for SVG
    if not xs:
        return (0, 0, 1, 1)
    minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)
    w, h = maxx - minx, maxy - miny
    if w == 0:
        w = 1
    if h == 0:
        h = 1
    pad = PAD * max(w, h)
    return (minx - pad, miny - pad, w + 2 * pad, h + 2 * pad)


def path_el(rings, fill, stroke, fa, sw):
    parts = []
    for r in rings:
        if len(r) < 2:
            continue
        parts.append("M" + " ".join(f"{x:.6g},{-y:.6g}" for x, y in r) + "Z")
    d = "".join(parts)
    if not d:
        return ""
    fillattr = (f'fill="{fill}" fill-opacity="{fa}" fill-rule="evenodd"'
                if fill else 'fill="none"')
    return (f'<path d="{d}" {fillattr} stroke="{stroke}" stroke-width="{sw}" '
            'vector-effect="non-scaling-stroke"/>')


def dots_el(rings, r, fill, stroke, sw):
    out = []
    fa = f'fill="{fill}"' if fill else 'fill="none"'
    sa = (f'stroke="{stroke}" stroke-width="{sw}" vector-effect="non-scaling-stroke"'
          if stroke else "")
    for ring in rings:
        for x, y in ring:
            out.append(f'<circle cx="{x:.6g}" cy="{-y:.6g}" r="{r:.6g}" {fa} {sa}/>')
    return "".join(out)


def svg(vb, inner, size=300):
    minx, miny, w, h = vb
    return (f'<svg viewBox="{minx:.6g} {miny:.6g} {w:.6g} {h:.6g}" '
            f'width="{size}" height="{size}" '
            f'preserveAspectRatio="xMidYMid meet">{inner}</svg>')


def cap(tag, e, other):
    nv, ar, rg = e["numVert"], e["area"], e["rings"]
    dv = "" if nv == other["numVert"] else f' <b>(Δ{nv-other["numVert"]:+d})</b>'
    da = ("" if abs(e["area"] - other["area"]) < 1e-6 * max(1, abs(e["area"]))
          else ' <b>area≠</b>')
    return f'<div class="cap">{tag}: verts={nv}{dv} area={ar:.4g}{da} rings={rg}</div>'


rows = []
for sa, sb in zip(c2["scenarios"], b2["scenarios"]):
    vb = bounds([sa["input"], sa["output"], sb["output"]])
    rdot = 0.012 * max(vb[2], vb[3])
    inlayer = path_el(sa["input"], None, GRAY, 0, 1.0)
    panel_c2 = svg(vb, inlayer + path_el(sa["output"], BLUE, BLUE, 0.55, 1.3))
    panel_b2 = svg(vb, inlayer + path_el(sb["output"], ORANGE, ORANGE, 0.55, 1.3))
    # Middle: faint shared boundary, then C2 hollow rings + B2 filled dots.
    panel_ov = svg(vb,
                   path_el(sa["output"], None, "#ccc", 0, 1.2)
                   + dots_el(sa["output"], rdot, None, BLUE, 1.4)
                   + dots_el(sb["output"], rdot * 0.55, ORANGE, None, 0))
    rows.append(f'''<section>
  <h2>{sa["name"]}</h2>
  <div class="grid">
    <figure><figcaption>Clipper2</figcaption>{panel_c2}{cap("C2", sa, sb)}</figure>
    <figure><figcaption>vertices</figcaption>{panel_ov}
       <div class="cap"><span style="color:{BLUE}">&#9711; C2</span> &nbsp; <span style="color:{ORANGE}">&#9679; B2</span> &nbsp; lone marker = backend-only</div></figure>
    <figure><figcaption>Boolean2</figcaption>{panel_b2}{cap("B2", sb, sa)}</figure>
  </div>
</section>''')

html = f'''<!doctype html><meta charset="utf-8">
<title>Clipper2 vs Boolean2 - offset joins</title>
<style>
 body{{font:14px/1.4 system-ui,sans-serif;margin:24px;color:#222;background:#fafafa}}
 h1{{font-size:20px}} h2{{font-size:15px;margin:18px 0 6px}}
 .grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;max-width:1000px}}
 figure{{margin:0;background:#fff;border:1px solid #ddd;border-radius:6px;padding:8px;text-align:center}}
 figcaption{{font-weight:600;color:#555;margin-bottom:4px}}
 svg{{width:100%;height:auto;background:#fff}}
 .cap{{font-size:12px;color:#444;margin-top:6px}}
 section{{margin-bottom:8px}}
</style>
<h1>Clipper2 vs Boolean2 &mdash; CrossSection offset joins</h1>
<p>Same input (gray outline) offset through each backend. Left/right show the
   filled result. The backends agree on shape here, so the middle column draws
   every output vertex (<span style="color:{BLUE}">&#9711; C2 ring</span> /
   <span style="color:{ORANGE}">&#9679; B2 dot</span>): a lone marker is a vertex
   only one backend emits.</p>
{"".join(rows)}
'''
open(out_path, "w").write(html)
print(f"wrote {out_path} ({len(c2['scenarios'])} scenarios)")
