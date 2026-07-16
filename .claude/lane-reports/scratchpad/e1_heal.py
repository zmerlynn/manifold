#!/usr/bin/env python3
# E1 de-risk: soup-level T-junction healing. Split each triangle edge at any
# welded vertex lying strictly on it, re-weld, measure opens. Tests whether
# SUBDIVISION coordination (T-junction consistency) closes the openscad opens,
# separated from RETENTION.
#
# usage: python3 e1_heal.py [PROD|<soupfile>]
#   PROD        -> heal the production emit (from regen_dump.txt)
#   <soupfile>  -> heal a true_boundary.py soup (9 hex floats/line = one tri)
import sys, math
import importlib.util
spec = importlib.util.spec_from_file_location("oracle", "/tmp/oracle-asm/oracle.py")
oracle = importlib.util.module_from_spec(spec); spec.loader.exec_module(oracle)
parse=oracle.parse; build_impl=oracle.build_impl; open_edges=oracle.open_edges

comps=parse(); comp=max(comps,key=lambda c:len(c["tris"])); eps=comp["eps"]

def load_soup(path):
    emit=[]
    with open(path) as fp:
        for ln in fp:
            p=ln.split()
            if len(p)<9: continue
            cs=[float.fromhex(x) for x in p[:9]]
            emit.append((-2,(cs[0],cs[1],cs[2]),(cs[3],cs[4],cs[5]),(cs[6],cs[7],cs[8])))
    return emit

arg = sys.argv[1] if len(sys.argv)>1 else "PROD"
emit = comp["emit"] if arg=="PROD" else load_soup(arg)

def report(tag, verts, tv):
    opens,_=open_edges(verts,tv)
    fan={}
    for o in opens: fan[len(o["hes"])]=fan.get(len(o["hes"]),0)+1
    print("%-14s verts=%d tris=%d OPENS=%d fan=%s"%(tag,len(verts),len(tv),len(opens),dict(sorted(fan.items()))))
    return opens

verts, tv, tvprov, tvsrc = build_impl(emit, eps)
report("baseline", verts, tv)

# T-junction healing on the welded soup. A vertex v lies strictly on edge (a,b)
# iff collinear (perp dist < tolC*len) and strictly between (0<t<1 with param
# margin). Split the triangle carrying (a,b) into (a,v,c)+(v,b,c). Iterate.
def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def cross(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def nrm(a): return math.sqrt(dot(a,a))

# grid of verts for fast "on edge" query: bucket by eps cell
def cellOf(p): return (math.floor(p[0]/eps),math.floor(p[1]/eps),math.floor(p[2]/eps))
grid={}
for i,p in enumerate(verts): grid.setdefault(cellOf(p),[]).append(i)

def verts_near_seg(a,b):
    # candidate vertex indices within a few eps of the segment bbox
    lo=[min(a[k],b[k]) for k in range(3)]; hi=[max(a[k],b[k]) for k in range(3)]
    c0=cellOf(tuple(lo[k]-2*eps for k in range(3))); c1=cellOf(tuple(hi[k]+2*eps for k in range(3)))
    out=[]
    for cx in range(c0[0],c1[0]+1):
        for cy in range(c0[1],c1[1]+1):
            for cz in range(c0[2],c1[2]+1):
                out.extend(grid.get((cx,cy,cz),[]))
    return out

TOLC = 1e-9   # perp-distance tol for "on line" (between eps and true feature sep)
def on_edge(vi, ai, bi):
    a=verts[ai]; b=verts[bi]; v=verts[vi]
    ab=sub(b,a); L=nrm(ab)
    if L==0: return None
    t=dot(sub(v,a),ab)/(L*L)
    if not (eps/L < t < 1.0-eps/L): return None
    proj=(a[0]+t*ab[0],a[1]+t*ab[1],a[2]+t*ab[2])
    if nrm(sub(v,proj)) > TOLC: return None
    return t

# heal to fixpoint (bounded passes)
tv=list(tv)
for _pass in range(6):
    newtv=[]; splits=0
    for tri in tv:
        a,b,c=tri
        # find a splitter on any edge; split ONE edge per triangle per pass
        done=False
        for (x,y,z) in ((a,b,c),(b,c,a),(c,a,b)):
            cands=verts_near_seg(verts[x],verts[y])
            best=None; bestt=None
            for vi in cands:
                if vi==x or vi==y or vi==z: continue
                t=on_edge(vi,x,y)
                if t is None: continue
                if best is None or t<bestt: best=vi; bestt=t
            if best is not None:
                newtv.append((x,best,z)); newtv.append((best,y,z)); splits+=1; done=True; break
        if not done: newtv.append(tri)
    tv=newtv
    if splits==0: break
report("healed", verts, tv)
