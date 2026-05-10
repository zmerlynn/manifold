#!/usr/bin/env python3
"""Convert the JTS robust/overlay XML test corpus to a flat text format.

Source:
  github.com/locationtech/jts → modules/tests/src/test/resources/testxml/
  robust/overlay/*.xml — <case> blocks with hex-WKB <a> and optional <b>,
  plus a <test><op name="overlayAreaTest" arg1="A" arg2="B"/> ...</test>
  block per case. Some files use unionArea with a numeric oracle.

Output: a single text file consumed by overlap2d_proto's `jtscorpus`
mode. One case per block, separated by blank lines:

  CASE <n> <op> <source_file>#<idx> [<expected_value>]
  A <num_rings>
  <num_pts> x y x y ...        (one ring per line; coords space-separated)
  ...
  B <num_rings>                 (omit block entirely if no B)
  <num_pts> x y x y ...
  <blank line>

Only `overlayAreaTest`, `overlayareatest`, and `unionArea` ops are
emitted (intersection-as-WKB, unionLength, isValid skipped). The WKB
decoder handles ISO and EWKB Polygon (3) / MultiPolygon (6) /
GeometryCollection (7) / MultiPoint (4) / MultiLineString (5), including
SRID/Z/M flags; non-polygon geometries contribute no rings. Inputs that
fail to parse are skipped with a stderr note (no abort).

Usage:
  jts_to_corpus.py <jts_overlay_dir> <output.txt>
"""

import os
import struct
import sys
import xml.etree.ElementTree as ET
from glob import glob


def hex_to_bytes(text):
    return bytes.fromhex("".join(text.split()))


def parse_wkt_rings(text):
    """Extract all rings from a WKT geometry string.

    A 'ring' is the content of an innermost parenthesized block (one with
    no nested parens before its closing ')'). This naturally handles
    POLYGON, MULTIPOLYGON, and GEOMETRYCOLLECTION without needing to
    distinguish them — POINT/LINESTRING contribute leaf parens that yield
    rings with fewer than 3 points, which we drop.

    Coordinates within each ring are comma-separated `x y [z [m]]`
    triples; we use the first two components and ignore Z/M.
    """
    rings = []
    n = len(text)
    i = 0
    while i < n:
        if text[i] == '(':
            # Scan to the matching ')' or first nested '('.
            j = i + 1
            has_inner = False
            while j < n and text[j] != ')':
                if text[j] == '(':
                    has_inner = True
                    break
                j += 1
            if not has_inner and j < n:
                ring = []
                for piece in text[i + 1:j].split(','):
                    nums = piece.split()
                    if len(nums) >= 2:
                        try:
                            ring.append((float(nums[0]), float(nums[1])))
                        except ValueError:
                            pass
                if len(ring) >= 2 and ring[0] == ring[-1]:
                    ring.pop()
                if len(ring) >= 3:
                    rings.append(ring)
                i = j  # land on the ')'
        i += 1
    return rings


def parse_geometry(text):
    """Decode WKT or hex-WKB to a list of rings. Whitespace-trimmed input
    starting with '(' or an alphabetic char is treated as WKT; anything
    else as hex."""
    s = text.strip()
    if not s:
        return []
    if s[0].isalpha() or s[0] == '(':
        return parse_wkt_rings(s)
    return parse_wkb(hex_to_bytes(s))[0]


def parse_wkb(b, off=0):
    """Return (rings, new_offset). Rings come from any Polygon/MultiPolygon
    embedded in the WKB stream; other geometry types are skipped (their
    bytes are still consumed)."""
    order = b[off]
    off += 1
    if order not in (0, 1):
        raise ValueError(f"bad WKB byte order {order} at offset {off-1}")
    fmt_u32 = "<I" if order == 1 else ">I"
    fmt_d = "<d" if order == 1 else ">d"
    t = struct.unpack_from(fmt_u32, b, off)[0]
    off += 4

    has_z = bool(t & 0x80000000)
    has_m = bool(t & 0x40000000)
    has_srid = bool(t & 0x20000000)
    if has_srid:
        off += 4
    iso = t & 0x1FFFFFFF  # strip EWKB flag bits
    base = iso & 0xFF
    if iso >= 1000:
        if iso >= 3000:
            has_z = True
            has_m = True
            base = iso - 3000
        elif iso >= 2000:
            has_m = True
            base = iso - 2000
        else:
            has_z = True
            base = iso - 1000

    coord_extra = (8 if has_z else 0) + (8 if has_m else 0)
    rings = []

    if base == 1:  # Point
        off += 16 + coord_extra
    elif base == 2:  # LineString
        n = struct.unpack_from(fmt_u32, b, off)[0]
        off += 4
        off += n * (16 + coord_extra)
    elif base == 3:  # Polygon
        nr = struct.unpack_from(fmt_u32, b, off)[0]
        off += 4
        for _ in range(nr):
            np = struct.unpack_from(fmt_u32, b, off)[0]
            off += 4
            ring = []
            for _ in range(np):
                x = struct.unpack_from(fmt_d, b, off)[0]
                off += 8
                y = struct.unpack_from(fmt_d, b, off)[0]
                off += 8
                if has_z:
                    off += 8
                if has_m:
                    off += 8
                ring.append((x, y))
            if len(ring) >= 2 and ring[0] == ring[-1]:
                ring.pop()
            if len(ring) >= 3:
                rings.append(ring)
    elif base in (4, 5, 6, 7):  # MultiPoint / MultiLineString / MultiPolygon / GeometryCollection
        nc = struct.unpack_from(fmt_u32, b, off)[0]
        off += 4
        for _ in range(nc):
            sub_rings, off = parse_wkb(b, off)
            rings.extend(sub_rings)
    else:
        raise ValueError(f"unknown WKB type {t} (base {base})")

    return rings, off


def emit_geom(name, rings, out):
    out.write(f"{name} {len(rings)}\n")
    for ring in rings:
        out.write(str(len(ring)))
        for x, y in ring:
            out.write(f" {x!r} {y!r}")
        out.write("\n")


def main(jts_dir, output):
    files = sorted(glob(os.path.join(jts_dir, "*.xml")))
    if not files:
        print(f"no XML files found in {jts_dir}", file=sys.stderr)
        return 1
    case_n = 0
    skipped = 0
    with open(output, "w") as out:
        out.write(
            "# JTS overlay corpus, converted from JTS testxml/robust/overlay.\n"
            "# Source: github.com/locationtech/jts (EDL/EPL dual-licensed; EDL\n"
            "# is BSD-3-Clause-equivalent and Apache-2.0 inbound-compatible).\n"
            "#\n"
            "# Format per case:\n"
            "#   CASE <n> <op> <source_file>#<idx> [<expected_value>]\n"
            "#   A <num_rings>\n"
            "#   <num_pts> x y x y ...    (one ring per line)\n"
            "#   B <num_rings>            (omit if no B)\n"
            "#   <num_pts> x y x y ...\n"
            "#   <blank line>\n"
            "#\n"
            "# Operations:\n"
            "#   overlayAreaTest — invariant check: |A∪B| + |A∩B| = |A| + |B|\n"
            "#   unionArea       — single-input regularization, <expected_value>\n"
            "#                     is the JTS-stored union area\n\n"
        )
        for f in files:
            base = os.path.basename(f)
            try:
                root = ET.parse(f).getroot()
            except ET.ParseError as e:
                print(f"skip {base}: parse error {e}", file=sys.stderr)
                continue
            for idx, case in enumerate(root.findall(".//case")):
                op_el = case.find(".//op")
                if op_el is None:
                    continue
                op = op_el.attrib.get("name", "")
                if op not in ("overlayAreaTest", "overlayareatest", "unionArea"):
                    continue
                a_el = case.find("a")
                if a_el is None or not a_el.text or not a_el.text.strip():
                    skipped += 1
                    continue
                try:
                    a_rings = parse_geometry(a_el.text)
                except Exception as e:
                    print(f"skip {base}#{idx} (A): {e}", file=sys.stderr)
                    skipped += 1
                    continue
                if not a_rings:
                    skipped += 1
                    continue
                b_rings = []
                b_el = case.find("b")
                if b_el is not None and b_el.text and b_el.text.strip():
                    try:
                        b_rings = parse_geometry(b_el.text)
                    except Exception as e:
                        print(f"skip {base}#{idx} (B): {e}", file=sys.stderr)
                        skipped += 1
                        continue
                expected = ""
                if op_el.text:
                    val_text = op_el.text.strip()
                    if val_text:
                        try:
                            float(val_text)
                            expected = " " + val_text
                        except ValueError:
                            pass
                case_n += 1
                out.write(f"CASE {case_n} {op} {base}#{idx}{expected}\n")
                emit_geom("A", a_rings, out)
                if b_rings:
                    emit_geom("B", b_rings, out)
                out.write("\n")
    print(f"wrote {case_n} cases ({skipped} skipped)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "usage: jts_to_corpus.py <jts_overlay_dir> <output.txt>",
            file=sys.stderr,
        )
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2]))
