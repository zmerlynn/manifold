#!/usr/bin/env python3
# PER-EDGE FAN TRUTH for the C++ engine's opens (final18 circular-fan method):
# around each open edge, find the TRUE {w>=1} boundary carriers by exact
# winding over the radial sectors of all carrier input faces, and compare with
# the sheets the engine actually emitted.  Names the missing-partner class.
# usage: python3 e1_fan.py [N=40] [R=1e-7]
import sys, math
from fractions import Fraction as Fr
import importlib.util
spec = importlib.util.spec_from_file_location("oracle", "/tmp/oracle-asm/oracle.py")
o = importlib.util.module_from_spec(spec); spec.loader.exec_module(o)

comps=o.parse(); comp=max(comps,key=lambda c:len(c["tris"])); eps=comp["eps"]
Ff=o.build_faces_f(comp)
V=[o.R3(v) for v in comp["verts"]]; T=comp["tris"]; nT=len(T)

# exact plane per input face + engine-equivalent geometric plane groups
plane=[]
fkey=[None]*nT
groups={}
for f in range(nT):
    v0,v1,v2,_=T[f]
    A,B,C=V[v0],V[v1],V[v2]
    n=o.cross(o.vsub(B,A),o.vsub(C,A)); d=o.dot(n,A)
    plane.append((n,d,A,B,C))
    lead=[i for i in range(3) if n[i]!=0]
    if lead:
        q=n[lead[0]]
        fkey[f]=(n[0]/q,n[1]/q,n[2]/q,d/q)
        groups.setdefault(fkey[f],len(groups))
gid=[groups.get(fkey[f],-1) if fkey[f] is not None else -1 for f in range(nT)]

# engine soup with per-tri gid
gids=[]; emit=[]
for ln in open("/tmp/oracle-asm/e1_soup.txt"):
    p=ln.split(); gids.append(int(p[0]))
    cs=[float.fromhex(x) for x in p[1:10]]
    emit.append((-2,(cs[0],cs[1],cs[2]),(cs[3],cs[4],cs[5]),(cs[6],cs[7],cs[8])))
verts,tv,tvprov,tsrc=o.build_impl(emit,eps)
opens,_=o.open_edges(verts,tv)
print("engine soup: verts=%d tris=%d opens=%d"%(len(verts),len(tv),len(opens)))

NS=40; R=1e-7
for a in sys.argv[1:]:
    if a.startswith("N="): NS=int(a.split("=")[1])
    if a.startswith("R="): R=float(a.split("=")[1])

def norm(a):
    l=math.sqrt(sum(x*x for x in a)); return tuple(x/l for x in a)
def subf(a,b): return (a[0]-b[0],a[1]-b[1],a[2]-b[2])
def dotf(a,b): return sum(x*y for x,y in zip(a,b))
def crossf(a,b): return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])

def in_tri3(P,A,B,C,margin):
    # 3D point-in-tri with relative margin (P assumed near the plane)
    n=crossf(subf(B,A),subf(C,A)); nn=dotf(n,n)
    if nn==0: return False
    m=margin*nn
    s0=dotf(n,crossf(subf(B,A),subf(P,A)))
    s1=dotf(n,crossf(subf(C,B),subf(P,B)))
    s2=dotf(n,crossf(subf(A,C),subf(P,C)))
    return s0>=-m and s1>=-m and s2>=-m

import random
random.seed(11)
sample = list(range(len(opens)))
random.shuffle(sample)
sample = sample[:NS]

from collections import Counter
classes=Counter(); missdetail=Counter()
for oi in sample:
    op=opens[oi]
    e=op["edge"]; pa,pb=verts[e[0]],verts[e[1]]
    M=tuple((pa[k]+pb[k])/2 for k in range(3))
    d=norm(subf(pb,pa))
    # perpendicular basis
    ref=(1.0,0,0) if abs(d[0])<0.9 else (0,1.0,0)
    u=norm(crossf(d,ref)); v=crossf(d,u)
    # carrier input faces: plane contains both endpoints (tol), tri spans M
    carriers=[]
    for f in range(nT):
        n,dd,A,B,C=plane[f]
        nf=(float(n[0]),float(n[1]),float(n[2])); df=float(dd)
        nl=math.sqrt(dotf(nf,nf))
        if nl==0: continue
        if abs(dotf(nf,pa)-df)/nl>1e-6 or abs(dotf(nf,pb)-df)/nl>1e-6: continue
        Af=tuple(float(x) for x in A); Bf=tuple(float(x) for x in B); Cf=tuple(float(x) for x in C)
        if not in_tri3(M,Af,Bf,Cf,1e-9): continue
        carriers.append(f)
    # crossing angles of each carrier plane with the circle
    angs=[]
    for f in carriers:
        n,dd,A,B,C=plane[f]
        nf=(float(n[0]),float(n[1]),float(n[2])); df=float(dd)
        Aa=R*dotf(nf,u); Bb=R*dotf(nf,v); Cc=df-dotf(nf,M)
        Rr=math.hypot(Aa,Bb)
        if Rr==0: continue
        cc=max(-1.0,min(1.0,Cc/Rr))
        phi=math.atan2(Bb,Aa)
        for s in (1,-1):
            th=phi+s*math.acos(cc)
            P=tuple(M[k]+R*(math.cos(th)*u[k]+math.sin(th)*v[k]) for k in range(3))
            Af=tuple(float(x) for x in A); Bf=tuple(float(x) for x in B); Cf=tuple(float(x) for x in C)
            if in_tri3(P,Af,Bf,Cf,1e-9):
                angs.append((th%(2*math.pi),f))
    if not angs:
        classes["no-carrier-crossings"]+=1; continue
    angs.sort()
    # sector windings (exact) at interval midpoints
    ths=[a for a,_ in angs]
    mids=[]
    for k in range(len(ths)):
        t0=ths[k]; t1=ths[(k+1)%len(ths)]
        if k==len(ths)-1: t1+=2*math.pi
        mids.append(((t0+t1)/2)%(2*math.pi))
    ws=[]
    okw=True
    for th in mids:
        P=tuple(M[k]+R*(math.cos(th)*u[k]+math.sin(th)*v[k]) for k in range(3))
        w,dg=o.wSf(o.R3(P),Ff)
        if dg: okw=False
        ws.append(w)
    if not okw:
        classes["degen-probe"]+=1; continue
    # true boundary sheets: crossing k separates sector k-1 (mid[k-1]) and
    # sector k (mid[k]); boundary iff the {w>=1} class flips there
    true_sheets=[]   # (angle, face)
    for k in range(len(ths)):
        wprev=ws[(k-1)%len(ws)]; wnext=ws[k]
        if (wprev>=1)!=(wnext>=1):
            true_sheets.append((ths[k],angs[k][1]))
    # emitted sheets at this edge: incident tris' third-vertex angles + gid
    emitted=[]
    for h in op["hes"]:
        t=tv[h//3]; g=gids[tsrc[h//3]]
        third=[x for x in t if x not in e]
        if not third: continue
        tp=verts[third[0]]
        w3=subf(tp,M)
        th=math.atan2(dotf(w3,v),dotf(w3,u))%(2*math.pi)
        emitted.append((th,g))
    # classify: for each true sheet, is an emitted sheet within 0.15 rad?
    missing=[]
    for (th,f) in true_sheets:
        near=any(abs(((th-te+math.pi)%(2*math.pi))-math.pi)<0.15 for te,_ in emitted)
        if not near: missing.append((th,f))
    extra=[]
    for (te,g) in emitted:
        near=any(abs(((te-th+math.pi)%(2*math.pi))-math.pi)<0.15 for th,_ in true_sheets)
        if not near: extra.append((te,g))
    egids=set(g for _,g in emitted)
    if not true_sheets:
        classes["no-true-boundary (engine over-emits)"]+=1
    elif missing and not extra:
        # in-plane (same group as an emitting sheet) vs cross-plane partner
        for th,f in missing:
            if gid[f] in egids: missdetail["missing IN-PLANE (same group %s)"%""]+=1
            else: missdetail["missing CROSS-PLANE group"]+=1
        classes["missing-partner"]+=1
    elif missing and extra:
        classes["missing+extra (shifted)"]+=1
    elif extra:
        classes["extra-only"]+=1
    else:
        classes["all-matched (unweld/T-junction)"]+=1
print("classes over %d sampled opens:"%len(sample))
for k,c in classes.most_common(): print("  %-42s %d"%(k,c))
for k,c in missdetail.most_common(): print("  detail: %-34s %d"%(k,c))
