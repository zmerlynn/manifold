#!/usr/bin/env python3
# E1 de-risk (fast): for each OPEN edge classify T-JUNCTION (a welded vertex lies
# strictly on the open edge -> subdivision inconsistency) vs GENUINE HOLE (no
# vertex on it -> retention/absent-partner). Separates the two failure facets.
# usage: python3 e1_diag.py [PROD|<soupfile>]
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
verts, tv, tvprov, tvsrc = build_impl(emit, eps)
opens,_=open_edges(verts,tv)
fan={}
for o in opens: fan[len(o["hes"])]=fan.get(len(o["hes"]),0)+1
print("%s verts=%d tris=%d OPENS=%d fan=%s"%(arg,len(verts),len(tv),len(opens),dict(sorted(fan.items()))))

def sub(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def nrm(a): return math.sqrt(dot(a,a))

# perp distance tol for "on line": between eps and the model min feature sep.
TOLC=1e-8
def on_edge(vi,ai,bi):
    a=verts[ai]; b=verts[bi]; v=verts[vi]
    ab=sub(b,a); L=nrm(ab)
    if L==0: return None
    t=dot(sub(v,a),ab)/(L*L)
    if not (eps/L < t < 1.0-eps/L): return None   # strictly interior in param
    proj=(a[0]+t*ab[0],a[1]+t*ab[1],a[2]+t*ab[2])
    if nrm(sub(v,proj))>TOLC: return None
    return t

nV=len(verts)
tjunc=0; hole=0; tjex=[]
for ki,o in enumerate(opens):
    e=o["edge"]; ai,bi=e
    found=None
    for vi in range(nV):            # 18 * ~1500 verts = trivial
        if vi==ai or vi==bi: continue
        t=on_edge(vi,ai,bi)
        if t is not None: found=(vi,t); break
    if found is not None:
        tjunc+=1; tjex.append((ki,found[0],round(found[1],4)))
    else: hole+=1
print("  T-JUNCTION (vertex on open edge, subdivision) = %d  %s"%(tjunc,tjex[:8]))
print("  GENUINE HOLE (no vertex on edge, retention)    = %d"%hole)
