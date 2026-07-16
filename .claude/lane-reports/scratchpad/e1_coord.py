#!/usr/bin/env python3
# E1 COORDINATED ARRANGEMENT + IDENTITY-CARRIED EMISSION - offline validation.
#
# The coordinated recipe (what the C++ engine will implement):
#   - PER-PLANE exact 2D arrangements (exactly-coplanar faces grouped -> one
#     arrangement; fold multiplicity via per-cell covering jump), NOT per-face.
#   - Chords = genuine tri-tri seam segments (plane^plane line clipped to both
#     tris, exact rationals). Subdivision consistency across planes is BY VALUE:
#     exact arithmetic makes every shared vertex the identical rational point
#     (identity = exact value; the C++ port replicates this with canonical
#     once-only constructions).
#   - Per-cell classify: w_above via exact ray-cast winding at an interior point
#     offset by an ADAPTIVE sub-feature normal step (below the min distance to
#     any relevant foreign plane - the "infinitesimal" done soundly);
#     w_below = w_above + jump, CERTIFIED by probing both sides
#     (measured dw == combinatorial jump = arrangement-completeness check).
#   - Emit boundary cells ({w>=1} transition) as IDENTITY-INDEXED triangles
#     (no position weld). Opens computed on the index topology. Target: 0.
#   - Round each identity ONCE to double; identity-carried topology unchanged;
#     report min identity separation / degenerate tris; then the position-weld
#     comparison (build_impl) = the offline identity-weld mutation anchor.
#
# usage: python3 e1_coord.py [SAMPLE=n] [NOGWN]
import sys, math, time
from fractions import Fraction as Fr
import importlib.util
spec = importlib.util.spec_from_file_location("oracle", "/tmp/oracle-asm/oracle.py")
oracle = importlib.util.module_from_spec(spec); spec.loader.exec_module(oracle)
parse=oracle.parse; R3=oracle.R3; vsub=oracle.vsub; vadd=oracle.vadd
smul=oracle.smul; cross=oracle.cross; dot=oracle.dot
build_faces_f=oracle.build_faces_f; wSf=oracle.wSf
build_impl=oracle.build_impl; open_edges=oracle.open_edges

Z=Fr(0); ONE=Fr(1)
t_start=time.time()
def log(msg):
    sys.stderr.write("[%6.1fs] %s\n"%(time.time()-t_start,msg)); sys.stderr.flush()

comps=parse(); comp=max(comps,key=lambda c:len(c["tris"])); eps=comp["eps"]
eps_fr=Fr(eps)
V=[R3(v) for v in comp["verts"]]; T=comp["tris"]; nT=len(T)
Ff=build_faces_f(comp)
log("parsed: nT=%d eps=%.4g"%(nT,eps))

# ---------- exact planes + canonical plane groups ----------
plane_n=[None]*nT; plane_d=[None]*nT
fkey=[None]*nT; fsgn=[0]*nT
groups={}   # key -> dict(N=(n0,n1,n2) canonical, d, members=[(f,sgn)])
ndegen=0
for f in range(nT):
    v0,v1,v2,_=T[f]; A,B,C=V[v0],V[v1],V[v2]
    n=cross(vsub(B,A),vsub(C,A)); d=dot(n,A)
    plane_n[f]=n; plane_d[f]=d
    lead=None
    for i in range(3):
        if n[i]!=0: lead=i; break
    if lead is None: ndegen+=1; continue
    # canonical GEOMETRIC plane: divide by the SIGNED leading component so
    # anti-oriented coplanar faces share one group (leading comp exactly +1);
    # sgn = the face normal's direction relative to the canonical N.
    q=n[lead]
    s=1 if q>0 else -1
    key=(n[0]/q,n[1]/q,n[2]/q,d/q)
    fkey[f]=key; fsgn[f]=s
    g=groups.setdefault(key,dict(N=(key[0],key[1],key[2]),d=key[3],members=[]))
    g["members"].append((f,s))
log("plane groups=%d degenerate faces=%d"%(len(groups),ndegen))

# ---------- float bboxes ----------
fbb=[]
for f in range(nT):
    v0,v1,v2,_=T[f]
    xs=[comp["verts"][v][0] for v in (v0,v1,v2)]
    ys=[comp["verts"][v][1] for v in (v0,v1,v2)]
    zs=[comp["verts"][v][2] for v in (v0,v1,v2)]
    fbb.append((min(xs),max(xs),min(ys),max(ys),min(zs),max(zs)))
def bb_overlap(i,j,m=1e-7):
    a=fbb[i]; b=fbb[j]
    return (a[0]-m<=b[1] and b[0]-m<=a[1] and a[2]-m<=b[3] and b[2]-m<=a[3]
            and a[4]-m<=b[5] and b[4]-m<=a[5])

# ---------- chords: genuine tri-tri seam segments, exact ----------
def line_tri_tparams(L0,dirv,triA,triB,triC):
    n=cross(vsub(triB,triA),vsub(triC,triA))
    an=(abs(n[0]),abs(n[1]),abs(n[2])); drop=an.index(max(an))
    ax=[i for i in range(3) if i!=drop]
    def pr(P): return (P[ax[0]],P[ax[1]])
    a2,b2,c2=pr(triA),pr(triB),pr(triC); l0=pr(L0); d2=(dirv[ax[0]],dirv[ax[1]])
    if d2==(Z,Z): return None
    hits=[]
    for (p,q) in ((a2,b2),(b2,c2),(c2,a2)):
        r=(q[0]-p[0],q[1]-p[1]); den=d2[0]*r[1]-d2[1]*r[0]
        if den==0: continue
        qp=(p[0]-l0[0],p[1]-l0[1])
        t=(qp[0]*r[1]-qp[1]*r[0])/den; u=(qp[0]*d2[1]-qp[1]*d2[0])/den
        if 0<=u<=1: hits.append(t)
    if len(hits)<2: return None
    return (min(hits),max(hits))

group_chords={k:[] for k in groups}   # key -> list of (Pa3,Pb3)
nchord=0
for f in range(nT):
    kf=fkey[f]
    if kf is None: continue
    n=plane_n[f]; d=plane_d[f]
    A=V[T[f][0]]; B=V[T[f][1]]; C=V[T[f][2]]
    for g in range(nT):
        if g==f or fkey[g] is None or fkey[g]==kf: continue
        if not bb_overlap(f,g): continue
        ng=plane_n[g]; dg=plane_d[g]
        dirn=cross(n,ng)
        if dirn==(Z,Z,Z): continue          # parallel distinct plane
        adirs=[abs(float(dirn[0])),abs(float(dirn[1])),abs(float(dirn[2]))]
        k=adirs.index(max(adirs))            # fix coord k: det = +-dirn[k] != 0
        i0,i1=[i for i in range(3) if i!=k]
        a11,a12,r1=n[i0],n[i1],d-n[k]*A[k]
        a21,a22,r2=ng[i0],ng[i1],dg-ng[k]*A[k]
        det=a11*a22-a12*a21
        if det==0: continue                  # cannot happen for k=argmax|dirn|
        X=(r1*a22-a12*r2)/det; Y=(a11*r2-r1*a21)/det
        P0=[Z,Z,Z]; P0[k]=A[k]; P0[i0]=X; P0[i1]=Y; P0=tuple(P0)
        tf=line_tri_tparams(P0,dirn,A,B,C)
        tg=line_tri_tparams(P0,dirn,V[T[g][0]],V[T[g][1]],V[T[g][2]])
        if not tf or not tg: continue
        lo=max(tf[0],tg[0]); hi=min(tf[1],tg[1])
        if lo>=hi: continue
        Pa=tuple(P0[j]+lo*dirn[j] for j in range(3))
        Pb=tuple(P0[j]+hi*dirn[j] for j in range(3))
        group_chords[kf].append((Pa,Pb)); nchord+=1
log("chords=%d"%nchord)

# ---------- per-plane exact arrangement ----------
import functools
def build_arrangement(segs):
    # segs: list of ((x,y),(x,y)) exact 2D. Returns verts2 (exact), adj (sets).
    S=len(segs)
    pts_on=[[(Z,segs[i][0]),(ONE,segs[i][1])] for i in range(S)]
    for i in range(S):
        p,q=segs[i]; r=(q[0]-p[0],q[1]-p[1])
        for j in range(i+1,S):
            a,b=segs[j]; s=(b[0]-a[0],b[1]-a[1])
            den=r[0]*s[1]-r[1]*s[0]
            qp=(a[0]-p[0],a[1]-p[1])
            if den!=0:
                t=(qp[0]*s[1]-qp[1]*s[0])/den; u=(qp[0]*r[1]-qp[1]*r[0])/den
                if 0<=t<=1 and 0<=u<=1:
                    X=(p[0]+t*r[0],p[1]+t*r[1])
                    pts_on[i].append((t,X)); pts_on[j].append((u,X))
            else:
                # parallel: collinear overlap -> mutual endpoint insertion
                if qp[0]*r[1]-qp[1]*r[0]!=0: continue   # parallel, not collinear
                rr=r[0]*r[0]+r[1]*r[1]
                if rr==0: continue
                for E in (a,b):
                    t=((E[0]-p[0])*r[0]+(E[1]-p[1])*r[1])/rr
                    if 0<t<1: pts_on[i].append((t,E))
                ss=s[0]*s[0]+s[1]*s[1]
                if ss==0: continue
                for E in (p,q):
                    u=((E[0]-a[0])*s[0]+(E[1]-a[1])*s[1])/ss
                    if 0<u<1: pts_on[j].append((u,E))
    vid={}; verts=[]
    def vof(p):
        if p not in vid: vid[p]=len(verts); verts.append(p)
        return vid[p]
    adj={}
    for i in range(S):
        seq=[]
        for _,p in sorted(pts_on[i],key=lambda x:x[0]):
            v=vof(p)
            if not seq or seq[-1]!=v: seq.append(v)
        for k in range(len(seq)-1):
            a,b=seq[k],seq[k+1]
            adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)
    return verts,adj

def trace_cells(verts,adj):
    # exact CCW angular sort -> rotation-system face trace. CCW (s>0) loops.
    def cmp_dir(v):
        pv=verts[v]
        def cmp(a,b):
            pa=verts[a]; pb=verts[b]
            ax,ay=pa[0]-pv[0],pa[1]-pv[1]; bx,by=pb[0]-pv[0],pb[1]-pv[1]
            ha=0 if (ay>0 or (ay==0 and ax>0)) else 1
            hb=0 if (by>0 or (by==0 and bx>0)) else 1
            if ha!=hb: return -1 if ha<hb else 1
            cr=ax*by-ay*bx
            if cr>0: return -1
            if cr<0: return 1
            return 0
        return cmp
    cwprev={}
    for v in adj:
        nb=sorted(adj[v],key=functools.cmp_to_key(cmp_dir(v)))
        cwprev[v]={nb[k]:nb[k-1] for k in range(len(nb))}
    used=set(); cells=[]; negloops=[]; bad=0
    for a in list(adj.keys()):
        for b in list(adj[a]):
            if (a,b) in used: continue
            loop=[]; ca,cb=a,b; ok=True
            for _ in range(200000):
                used.add((ca,cb)); loop.append(ca)
                nb=cwprev[cb].get(ca)
                if nb is None: ok=False; break
                ca,cb=cb,nb
                if (ca,cb)==(a,b): break
            else: ok=False
            if not ok: bad+=1; continue
            # excise spurs (out-and-back: ..., u, X, u, ... -> ..., u, ...)
            changed=True
            while changed and len(loop)>=3:
                changed=False
                m=len(loop)
                for k in range(m):
                    if loop[(k-1)%m]==loop[(k+1)%m]:
                        # remove k and one duplicate neighbour instance
                        hi=max((k+1)%m,k); lo=min((k+1)%m,k)
                        del loop[hi]; del loop[lo]
                        changed=True; break
            if len(loop)<3: continue
            s=Z
            for k in range(len(loop)):
                x1,y1=verts[loop[k]]; x2,y2=verts[loop[(k+1)%len(loop)]]
                s+=x1*y2-x2*y1
            # KEEP the raw (possibly pinched) loop: positive loops are cells
            # (their pinch-connected negative lobes = hole boundaries, handled
            # by the weakly-simple triangulation); negative loops are
            # disconnected island hole-rings (assigned by containment below).
            if s>0: cells.append(loop)
            elif s<0: negloops.append(loop)
    return cells,negloops,bad

# ---------- exact polygon triangulation (earclip + diagonal fallback) ----------
def area2(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])

def seg_proper_cross(p,q,a,b):
    d1=area2(p,q,a); d2=area2(p,q,b); d3=area2(a,b,p); d4=area2(a,b,q)
    return ((d1>0)!=(d2>0)) and ((d3>0)!=(d4>0)) and d1!=0 and d2!=0 and d3!=0 and d4!=0

def on_open_seg(p,q,v):
    if area2(p,q,v)!=0: return False
    lo0,hi0=(p[0],q[0]) if p[0]<=q[0] else (q[0],p[0])
    lo1,hi1=(p[1],q[1]) if p[1]<=q[1] else (q[1],p[1])
    return (lo0<v[0]<hi0) or (lo1<v[1]<hi1)

def diag_split(poly,ids):
    # recursive exact diagonal-split triangulation of a (weakly-)simple CCW
    # polygon (collinear runs + keyhole-duplicated vertices allowed).
    n=len(poly)
    if n<3: return []
    if n==3:
        return [tuple(ids)] if area2(*poly)>0 else []
    for i in range(n):
        a,b,c=poly[i-1],poly[i],poly[(i+1)%n]
        for jo in range(2,n-1):
            j=(i+jo)%n
            d=poly[j]
            if d==b: continue   # duplicate position (keyhole twin)
            # in-cone at i (collinear prev/next treated as convex half-plane)
            if area2(a,b,c)>=0:
                if not(area2(b,c,d)>0 and area2(b,d,a)>0): continue
            else:
                if (area2(b,c,d)<=0 and area2(b,d,a)<=0): continue
            ok=True
            for k in range(n):
                k2=(k+1)%n
                if k in (i,j) or k2 in (i,j): continue
                if seg_proper_cross(b,d,poly[k],poly[k2]): ok=False; break
            if not ok: continue
            for k in range(n):
                if k in (i,j): continue
                if on_open_seg(b,d,poly[k]): ok=False; break
            if not ok: continue
            # split at (i,j)
            p1=[];i1=[];k=i
            while True:
                p1.append(poly[k]); i1.append(ids[k])
                if k==j: break
                k=(k+1)%n
            p2=[];i2=[];k=j
            while True:
                p2.append(poly[k]); i2.append(ids[k])
                if k==i: break
                k=(k+1)%n
            return diag_split(p1,i1)+diag_split(p2,i2)
    return None   # no valid diagonal (should not happen for simple polygons)

def triangulate(poly):
    # poly: CCW list of exact 2D points (simple; may have collinear runs).
    # Earclip; on stall, exact diagonal-split fallback on the remainder.
    n=len(poly); idx=list(range(n)); tris=[]
    stall=False
    while len(idx)>3:
        found=False
        m=len(idx)
        for k in range(m):
            a,b,c=idx[(k-1)%m],idx[k],idx[(k+1)%m]
            if area2(poly[a],poly[b],poly[c])<=0: continue
            ok=True
            for j in idx:
                if j in (a,b,c): continue
                d1=area2(poly[a],poly[b],poly[j])
                d2=area2(poly[b],poly[c],poly[j])
                d3=area2(poly[c],poly[a],poly[j])
                if d1>=0 and d2>=0 and d3>=0: ok=False; break
            if ok:
                tris.append((a,b,c)); idx.pop(k); found=True; break
        if not found:
            rem=diag_split([poly[i] for i in idx],idx)
            if rem is None: stall=True
            else: tris.extend(rem)
            return tris,stall
    if len(idx)==3:
        a,b,c=idx
        if area2(poly[a],poly[b],poly[c])!=0: tris.append((a,b,c))
    return tris,stall

def loop_area(loop,verts):
    s=Z
    for k in range(len(loop)):
        x1,y1=verts[loop[k]]; x2,y2=verts[loop[(k+1)%len(loop)]]
        s+=x1*y2-x2*y1
    return s

_SLOPES=[Fr(1,7919),Fr(3,104729),Fr(-5,1299709),Fr(7,15485863),Fr(-11,32452843)]
def point_in_loop(p,loop,verts,stats):
    # exact parity along ray p + t*(1,r), retry slopes on degeneracy
    for r in _SLOPES:
        cnt=0; ok=True
        for k in range(len(loop)):
            A=verts[loop[k]]; B=verts[loop[(k+1)%len(loop)]]
            dx=B[0]-A[0]; dy=B[1]-A[1]
            det=dy-r*dx
            if det==0:
                # parallel to ray; degenerate only if collinear with it
                if (A[1]-p[1])==r*(A[0]-p[0]): ok=False; break
                continue
            u=(-(r)*(p[0]-A[0])+(p[1]-A[1]))/det
            t=(dx*(p[1]-A[1])-dy*(p[0]-A[0]))/det
            if u==0 or u==1 or t==0: ok=False; break
            if 0<u<1 and t>0: cnt+=1
        if ok: return cnt%2==1
        stats["parityretry"]+=1
    return False

def attach_holes(cells,negloops,verts,stats):
    # assign each disconnected hole ring to its innermost containing cell,
    # then keyhole-splice it in (bridge = exact-visibility diagonal).
    holeof={}
    for nl in negloops:
        p=verts[nl[0]]
        best=None; bestA=None
        for ci,loop in enumerate(cells):
            if point_in_loop(p,loop,verts,stats):
                a=abs(loop_area(loop,verts))
                if bestA is None or a<bestA: bestA=a; best=ci
        if best is None: stats["orphanhole"]+=1; continue
        holeof.setdefault(best,[]).append(nl)
    out=[]
    for ci,loop in enumerate(cells):
        merged=list(loop)
        pend=[list(h) for h in holeof.get(ci,())]
        while pend:
            spliced=False
            for hi,hole in enumerate(pend):
                # all ring edges/vertices to respect: merged + every pending hole
                rings=[merged]+pend
                done=False
                for i in range(len(merged)):
                    b=verts[merged[i]]
                    for j in range(len(hole)):
                        d=verts[hole[j]]
                        if b==d: continue
                        ok=True
                        for ring in rings:
                            m=len(ring)
                            for k in range(m):
                                P=verts[ring[k]]; Q=verts[ring[(k+1)%m]]
                                if P==b or P==d or Q==b or Q==d: continue
                                if seg_proper_cross(b,d,P,Q): ok=False; break
                            if not ok: break
                        if ok:
                            for ring in rings:
                                for k in ring:
                                    v=verts[k]
                                    if v==b or v==d: continue
                                    if on_open_seg(b,d,v): ok=False; break
                                if not ok: break
                        if ok:
                            rot=hole[j:]+hole[:j]
                            merged=merged[:i+1]+rot+[hole[j]]+merged[i:]
                            done=True; break
                    if done: break
                if done:
                    pend.pop(hi); spliced=True; break
            if not spliced:
                stats["splicefail"]+=1
                break
        out.append(merged)
    return out

# ---------- classification helpers ----------
def point_in_tri2(p,a,b,c):
    d1=area2(a,b,p); d2=area2(b,c,p); d3=area2(c,a,p)
    neg=(d1<0)or(d2<0)or(d3<0); pos=(d1>0)or(d2>0)or(d3>0)
    return not(neg and pos)

SAMPLE=None; DO_GWN=True; DUMPAT=None; DUMPR=4e-4
for a in sys.argv[1:]:
    if a.startswith("SAMPLE="): SAMPLE=int(a.split("=")[1])
    if a=="NOGWN": DO_GWN=False
    if a.startswith("DUMPAT="): DUMPAT=tuple(float(x) for x in a.split("=")[1].split(","))
    if a.startswith("DUMPR="): DUMPR=float(a.split("=")[1])

idmap={}      # exact 3D tuple -> global id
idpos=[]      # global id -> exact 3D
def gid_of(p3):
    i=idmap.get(p3)
    if i is None:
        i=len(idmap); idmap[p3]=i; idpos.append(p3)
    return i

emitted=[]    # (i0,i1,i2) identity triples, outward-oriented
emitprov=[]   # parallel: emitting plane-group index
cellrec=[]    # (probeA, probeB, wA, wB) for GWN spot check
stats=dict(cells=0,jumpcells=0,emitcells=0,tristall=0,badloop=0,dangling=0,
           jumpviol=0,degenprobe=0,zeroarea=0,duptri=0,negloops=0,orphanhole=0,
           splicefail=0,areamis=0,parityretry=0)
minalt=1e300; nsliver=0

gkeys=list(groups.keys())
log("classifying %d plane groups..."%len(gkeys))
seen_tris=set()
for gi,key in enumerate(gkeys):
    G=groups[key]; N=G["N"]; dpl=G["d"]; members=G["members"]
    anf=(abs(float(N[0])),abs(float(N[1])),abs(float(N[2])))
    drop=anf.index(max(anf)); ax=[i for i in range(3) if i!=drop]
    def pr(P): return (P[ax[0]],P[ax[1]])
    nk=N[drop]
    def lift(p2):
        P=[Z,Z,Z]; P[ax[0]]=p2[0]; P[ax[1]]=p2[1]
        P[drop]=(dpl-N[ax[0]]*p2[0]-N[ax[1]]*p2[1])/nk
        return tuple(P)
    segs=[]
    mem2d=[]
    for (f,s) in members:
        A=V[T[f][0]]; B=V[T[f][1]]; C=V[T[f][2]]
        a2,b2,c2=pr(A),pr(B),pr(C)
        mem2d.append((a2,b2,c2,s))
        segs.append((a2,b2)); segs.append((b2,c2)); segs.append((c2,a2))
    for (Pa,Pb) in group_chords[key]:
        segs.append((pr(Pa),pr(Pb)))
    verts2,adj=build_arrangement(segs)
    for v in adj:
        if len(adj[v])==1:
            stats["dangling"]+=1
            p3=lift(verts2[v])
            log("DANGLING group=%d nmem=%d pos=(%.9g,%.9g,%.9g)"%(
                gi,len(members),float(p3[0]),float(p3[1]),float(p3[2])))
    cells,negloops,bad=trace_cells(verts2,adj)
    stats["badloop"]+=bad
    stats["cells"]+=len(cells)
    stats["negloops"]+=len(negloops)
    if negloops:
        cells=attach_holes(cells,negloops,verts2,stats)
    # L1 norm of canonical N for the adaptive offset
    L1N=abs(N[0])+abs(N[1])+abs(N[2])
    for loop in cells:
        poly=[verts2[i] for i in loop]
        # covering jump at an interior point (largest sub-tri centroid)
        tris,stall=triangulate(poly)
        if stall:
            stats["tristall"]+=1
            c0=lift(poly[0])
            log("TRISTALL group=%d n=%d at=(%.9g,%.9g,%.9g)"%(
                gi,len(poly),float(c0[0]),float(c0[1]),float(c0[2])))
        # AREA CERTIFICATE: exact triangulation coverage (holes respected)
        sarea=loop_area(loop,verts2)
        tarea=Z
        for (ia,ib,ic) in tris: tarea+=area2(poly[ia],poly[ib],poly[ic])
        if tarea!=sarea:
            stats["areamis"]+=1
            c0=lift(poly[0])
            log("AREAMIS group=%d n=%d tris=%d at=(%.9g,%.9g,%.9g)"%(
                gi,len(poly),len(tris),float(c0[0]),float(c0[1]),float(c0[2])))
        if not tris: continue
        best=None; bestA=Z
        for (ia,ib,ic) in tris:
            A2=area2(poly[ia],poly[ib],poly[ic])
            if A2>bestA: bestA=A2; best=(ia,ib,ic)
        if best is None: continue
        c2=tuple((poly[best[0]][j]+poly[best[1]][j]+poly[best[2]][j])/3
                 for j in range(2))
        jump=0
        for (a2,b2,cc2,s) in mem2d:
            if point_in_tri2(c2,a2,b2,cc2): jump+=s
        if jump==0: continue
        stats["jumpcells"]+=1
        c3=lift(c2)
        c3f=(float(c3[0]),float(c3[1]),float(c3[2]))
        # adaptive offset: below min dist to any relevant foreign plane
        tmax=eps_fr/(100*L1N)
        for g in range(nT):
            if fkey[g] is None or fkey[g]==key: continue
            bbg=fbb[g]; m=1e-6
            if not (bbg[0]-m<=c3f[0]<=bbg[1]+m and bbg[2]-m<=c3f[1]<=bbg[3]+m
                    and bbg[4]-m<=c3f[2]<=bbg[5]+m): continue
            ng=plane_n[g]
            num=dot(ng,c3)-plane_d[g]
            if num==0: continue      # plane through the point, sheet elsewhere
            if num<0: num=-num
            L1g=abs(ng[0])+abs(ng[1])+abs(ng[2])
            tg=num/(2*L1g*L1N)
            if tg<tmax: tmax=tg
        off=smul(tmax,N)
        pA=vadd(c3,off); pB=vsub(c3,off)
        wA,dgA=wSf(pA,Ff); wB,dgB=wSf(pB,Ff)
        if dgA or dgB: stats["degenprobe"]+=1
        if wB-wA!=jump:
            stats["jumpviol"]+=1
            log("JUMPVIOL group=%d jump=%d wA=%d wB=%d n=%d at=(%.9g,%.9g,%.9g)"%(
                gi,jump,wA,wB,len(poly),c3f[0],c3f[1],c3f[2]))
            # adjudicate with the fully-exact winding (slow; rare)
            global _FEXACT
            try: _FEXACT
            except NameError: _FEXACT=oracle.build_faces(comp)
            eA,dA2=oracle.winding(pA,_FEXACT); eB,dB2=oracle.winding(pB,_FEXACT)
            log("JUMPVIOL-EXACT wA=%d(%s) wB=%d(%s) -> fast %s"%(
                -eA,dA2,-eB,dB2,
                "WRONG" if (-eA!=wA or -eB!=wB) else "CONFIRMED (jump wrong?)"))
            wA,wB=-eA,-eB
        aboveIn=(wA>=1); belowIn=(wB>=1)
        if DUMPAT is not None and math.dist(c3f,DUMPAT)<DUMPR:
            lv=[lift(p) for p in poly]
            log("CELLDUMP g=%d jump=%d wA=%d wB=%d emit=%s n=%d tmax=%.3g loop=%s"%(
                gi,jump,wA,wB,aboveIn!=belowIn,len(poly),float(tmax),
                " ".join("(%.10g,%.10g,%.10g)"%(float(p[0]),float(p[1]),float(p[2]))
                         for p in lv)))
        if aboveIn==belowIn:
            cellrec.append((pA,pB,wA,wB))
            continue
        orient=1 if belowIn else -1
        stats["emitcells"]+=1
        cellrec.append((pA,pB,wA,wB))
        # identity-indexed emission
        l3=[lift(p) for p in poly]
        gids=[gid_of(p) for p in l3]
        for (ia,ib,ic) in tris:
            Pa,Pb,Pc=l3[ia],l3[ib],l3[ic]
            nrm=cross(vsub(Pb,Pa),vsub(Pc,Pa))
            sgn=dot(nrm,N)
            if sgn==0: stats["zeroarea"]+=1; continue
            i0,i1,i2=gids[ia],gids[ib],gids[ic]
            if (sgn>0)!=(orient>0): i1,i2=i2,i1
            keyt=tuple(sorted((i0,i1,i2)))
            if keyt in seen_tris: stats["duptri"]+=1
            seen_tris.add(keyt)
            emitted.append((i0,i1,i2)); emitprov.append(gi)
            # sliver stats on rounded coords
            da=[float(Pa[j]) for j in range(3)]; db=[float(Pb[j]) for j in range(3)]
            dc=[float(Pc[j]) for j in range(3)]
            e1=[db[j]-da[j] for j in range(3)]; e2=[dc[j]-da[j] for j in range(3)]
            cr=[e1[1]*e2[2]-e1[2]*e2[1],e1[2]*e2[0]-e1[0]*e2[2],
                e1[0]*e2[1]-e1[1]*e2[0]]
            ar=0.5*math.sqrt(sum(x*x for x in cr))
            lmax=max(math.dist(da,db),math.dist(db,dc),math.dist(dc,da))
            if lmax>0:
                alt=2*ar/lmax
                if alt>0 and alt<minalt: minalt=alt
                if ar<eps: nsliver+=1
    if gi%200==0:
        log("  group %d/%d cells=%d emit=%d ids=%d"%(gi,len(gkeys),
            stats["cells"],len(emitted),len(idmap)))

log("emitted tris=%d identities=%d"%(len(emitted),len(idmap)))
log("stats=%s minAlt=%.3e slivers(area<eps)=%d"%(stats,minalt,nsliver))

# ---------- THE NUMBER: opens on the identity-indexed topology ----------
he={}
for (a,b,c) in emitted:
    for (u,w) in ((a,b),(b,c),(c,a)):
        k=(min(u,w),max(u,w))
        he.setdefault(k,[0,0])
        he[k][0 if u<w else 1]+=1
opens=[k for k,v in he.items() if v[0]!=v[1]]
fan={}
for k in opens:
    n=sum(he[k]); fan[n]=fan.get(n,0)+1
print("IDENTITY-INDEXED OPENS = %d  fan=%s"%(len(opens),dict(sorted(fan.items()))))
# open-edge provenance dump
tri_at={}
for ti,(a,b,c) in enumerate(emitted):
    for (u,w) in ((a,b),(b,c),(c,a)):
        tri_at.setdefault((min(u,w),max(u,w)),[]).append(ti)
for k in opens[:40]:
    pa=idpos[k[0]]; pb=idpos[k[1]]
    gs=sorted(set(emitprov[t] for t in tri_at.get(k,[])))
    print("  OPEN id(%d,%d) groups=%s p=(%.9g,%.9g,%.9g)->(%.9g,%.9g,%.9g) len=%.3g"
          %(k[0],k[1],gs,float(pa[0]),float(pa[1]),float(pa[2]),
            float(pb[0]),float(pb[1]),float(pb[2]),
            math.dist([float(x) for x in pa],[float(x) for x in pb])))
print("certificates: jumpviol=%d tristall=%d badloop=%d dangling=%d degenprobe=%d zeroarea=%d duptri=%d"
      %(stats["jumpviol"],stats["tristall"],stats["badloop"],stats["dangling"],
        stats["degenprobe"],stats["zeroarea"],stats["duptri"]))

# ---------- exact volume ----------
vol6=Z
for (a,b,c) in emitted:
    vol6+=dot(idpos[a],cross(idpos[b],idpos[c]))
print("EXACT VOLUME = %.15g"%(float(vol6)/6.0))

# ---------- rounding: identity-carried ----------
rpos=[(float(p[0]),float(p[1]),float(p[2])) for p in idpos]
# min separation among distinct identities (report in eps units)
used=set()
for (a,b,c) in emitted: used.update((a,b,c))
usedl=sorted(used)
minsep=1e300
pairs=[]
byx=sorted(usedl,key=lambda i:rpos[i][0])
for ii in range(len(byx)):
    for jj in range(ii+1,len(byx)):
        a=byx[ii]; b=byx[jj]
        if rpos[b][0]-rpos[a][0]>10*eps: break
        dd=math.dist(rpos[a],rpos[b])
        if dd<minsep: minsep=dd
        if dd<eps: pairs.append((dd,a,b))
print("MIN IDENTITY SEPARATION = %.4g = %.3f eps  (identities used=%d)"
      %(minsep,minsep/eps,len(used)))
pairs.sort()
print("SUB-EPS identity pairs = %d; smallest:"%len(pairs))
for (dd,a,b) in pairs[:6]:
    print("  %.4g (%.3f eps) id%d=(%.17g,%.17g,%.17g) id%d=(%.17g,%.17g,%.17g)"
          %(dd,dd/eps,a,*rpos[a],b,*rpos[b]))
rdegen=sum(1 for (a,b,c) in emitted if rpos[a]==rpos[b] or rpos[b]==rpos[c] or rpos[a]==rpos[c])
print("rounded coincident-corner tris = %d"%rdegen)

# ---------- position-weld comparison (the identity-weld mutation anchor) ----------
soup=[(-2,rpos[a],rpos[b],rpos[c]) for (a,b,c) in emitted]
pv,ptv,_,_=build_impl(soup,eps)
popens,_=open_edges(pv,ptv)
pfan={}
for o in popens: pfan[len(o["hes"])]=pfan.get(len(o["hes"]),0)+1
print("POSITION-WELD (build_impl eps grid) OPENS = %d fan=%s  (verts %d vs identities %d)"
      %(len(popens),dict(sorted(pfan.items())),len(pv),len(used)))

# ---------- GWN spot check vs the input oracle ----------
if DO_GWN:
    comp2=dict(verts=rpos,tris=[(a,b,c,0) for (a,b,c) in emitted])
    Ff2=build_faces_f(comp2)
    import random
    random.seed(7)
    recs=cellrec if SAMPLE is None else random.sample(cellrec,min(SAMPLE,len(cellrec)))
    mism=0; ntest=0
    for (pA,pB,wA,wB) in recs:
        for (P,w) in ((pA,wA),(pB,wB)):
            we,dg=wSf(P,Ff2)
            if dg: continue
            ntest+=1
            if (1 if w>=1 else 0)!=we: mism+=1
    print("GWN SPOT CHECK: %d probes, mismatches=%d"%(ntest,mism))

# save
with open("/tmp/oracle-asm/e1_coord_soup.txt","w") as fp:
    for (a,b,c) in emitted:
        A,B,C=rpos[a],rpos[b],rpos[c]
        fp.write("%d %d %d %s %s %s %s %s %s %s %s %s\n"%(a,b,c,
            A[0].hex(),A[1].hex(),A[2].hex(),B[0].hex(),B[1].hex(),B[2].hex(),
            C[0].hex(),C[1].hex(),C[2].hex()))
log("wrote /tmp/oracle-asm/e1_coord_soup.txt")
