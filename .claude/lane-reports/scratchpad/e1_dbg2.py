#!/usr/bin/env python3
# Locate the missing in-plane crossing inside group 742 (plane of face 1276):
# probe-walk w_above from just-off the open edge A-B (742b side) toward the
# containing cell's centroid; find the jump; name the face(s) whose sheet
# crosses there and why the chord is absent.
import sys, math
from fractions import Fraction as Fr
import importlib.util
spec = importlib.util.spec_from_file_location("oracle", "/tmp/oracle-asm/oracle.py")
oracle = importlib.util.module_from_spec(spec); spec.loader.exec_module(oracle)
parse=oracle.parse; R3=oracle.R3; vsub=oracle.vsub; vadd=oracle.vadd
smul=oracle.smul; cross=oracle.cross; dot=oracle.dot
build_faces_f=oracle.build_faces_f; wSf=oracle.wSf
Z=Fr(0); ONE=Fr(1)

comps=parse(); comp=max(comps,key=lambda c:len(c["tris"])); eps=comp["eps"]
V=[R3(v) for v in comp["verts"]]; T=comp["tris"]; nT=len(T)
Ff=build_faces_f(comp)

plane_n=[None]*nT; plane_d=[None]*nT; fkey=[None]*nT
groups={}; gorder=[]
for f in range(nT):
    v0,v1,v2,_=T[f]; A,B,C=V[v0],V[v1],V[v2]
    n=cross(vsub(B,A),vsub(C,A)); d=dot(n,A)
    plane_n[f]=n; plane_d[f]=d
    lead=[i for i in range(3) if n[i]!=0]
    if not lead: continue
    q=n[lead[0]]
    key=(n[0]/q,n[1]/q,n[2]/q,d/q)
    fkey[f]=key
    if key not in groups:
        groups[key]=dict(N=(key[0],key[1],key[2]),d=key[3],members=[])
        gorder.append(key)
    groups[key]["members"].append((f,1 if q>0 else -1))

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

GKEY=gorder[742]
G=groups[GKEY]; N=G["N"]; dpl=G["d"]; members=G["members"]
print("g742 members:",[f for f,_ in members])
anf=(abs(float(N[0])),abs(float(N[1])),abs(float(N[2])))
drop=anf.index(max(anf)); ax=[i for i in range(3) if i!=drop]
def pr(P): return (P[ax[0]],P[ax[1]])
nk=N[drop]
def lift(p2):
    P=[Z,Z,Z]; P[ax[0]]=p2[0]; P[ax[1]]=p2[1]
    P[drop]=(dpl-N[ax[0]]*p2[0]-N[ax[1]]*p2[1])/nk
    return tuple(P)
L1N=abs(N[0])+abs(N[1])+abs(N[2])

def probe(c3):
    c3f=(float(c3[0]),float(c3[1]),float(c3[2]))
    tmax=Fr(eps)/(100*L1N)
    for g in range(nT):
        if fkey[g] is None or fkey[g]==GKEY: continue
        bbg=fbb[g]; m=1e-6
        if not (bbg[0]-m<=c3f[0]<=bbg[1]+m and bbg[2]-m<=c3f[1]<=bbg[3]+m
                and bbg[4]-m<=c3f[2]<=bbg[5]+m): continue
        ng=plane_n[g]
        num=dot(ng,c3)-plane_d[g]
        if num==0: continue
        if num<0: num=-num
        L1g=abs(ng[0])+abs(ng[1])+abs(ng[2])
        tg=num/(2*L1g*L1N)
        if tg<tmax: tmax=tg
    off=smul(tmax,N)
    wA,dgA=wSf(vadd(c3,off),Ff); wB,dgB=wSf(vsub(c3,off),Ff)
    return wA,wB,dgA or dgB

# A and B exact: A = triple {1263,1302,1276}; B = chord(1263x742) clip endpoint.
# Reconstruct A: intersection of planes 1263, 1302, 742-canonical (exact Cramer).
def plane3(f1,f2):  # intersect plane f1, plane f2, canonical G plane
    n1,d1=plane_n[f1],plane_d[f1]; n2,d2=plane_n[f2],plane_d[f2]
    n3,d3=N,dpl
    M=[[n1[0],n1[1],n1[2]],[n2[0],n2[1],n2[2]],[n3[0],n3[1],n3[2]]]
    rhs=[d1,d2,d3]
    det=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
        -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
        +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    if det==0: return None
    def col(k):
        MM=[row[:] for row in M]
        for r in range(3): MM[r][k]=rhs[r]
        return (MM[0][0]*(MM[1][1]*MM[2][2]-MM[1][2]*MM[2][1])
               -MM[0][1]*(MM[1][0]*MM[2][2]-MM[1][2]*MM[2][0])
               +MM[0][2]*(MM[1][0]*MM[2][1]-MM[1][1]*MM[2][0]))
    return (col(0)/det,col(1)/det,col(2)/det)

A3=plane3(1263,1302)
print("A =",[float(x) for x in A3])
# B: endpoint of chord(1263 x g742) clipped by tri 1263's edge shared with 1322:
# B is on planes 1263, 742, 1322.
B3=plane3(1263,1322)
print("B =",[float(x) for x in B3])

A2=pr(A3); B2=pr(B3)
M2=tuple((A2[j]+B2[j])/2 for j in range(2))
dirv=(B2[0]-A2[0],B2[1]-A2[1])
perp=(-dirv[1],dirv[0])
# side of C (the (B,A,C) cell): C = triple {1302,1322,742}
C3=plane3(1302,1322)
C2=pr(C3)
sC=(C2[0]-M2[0])*perp[0]+(C2[1]-M2[1])*perp[1]
sgn=-1 if sC>0 else 1   # opposite side from C
dl=Fr(1,10**7)
ln=math.sqrt(float(perp[0])**2+float(perp[1])**2)
step=dl/Fr(ln)  # rational scale approx
P2=(M2[0]+sgn*step*perp[0],M2[1]+sgn*step*perp[1])
P3=lift(P2)
wA,wB,dg=probe(P3)
print("just-off-edge (742b side, 1e-7): wA=%d wB=%d dg=%s"%(wA,wB,dg))

# walk outward: sample w at increasing distances along the same perpendicular
prev=(wA,wB)
for expo in range(6,1,-1):
    dl2=Fr(1,10**expo)
    st=dl2/Fr(ln)
    Q2=(M2[0]+sgn*st*perp[0],M2[1]+sgn*st*perp[1])
    Q3=lift(Q2)
    w1,w2,dg=probe(Q3)
    print("  dist=1e-%d: wA=%d wB=%d dg=%s"%(expo,w1,w2,dg))

# scan for sheets crossing between the near point and dist=1e-2 along perp
print("scan for crossing sheets along the perpendicular:")
far=Fr(1,100)/Fr(ln)
Q2far=(M2[0]+sgn*far*perp[0],M2[1]+sgn*far*perp[1])
Pn=lift(P2); Pf=lift(Q2far)
for g in range(nT):
    if fkey[g] is None or fkey[g]==GKEY: continue
    ng=plane_n[g]; dgq=plane_d[g]
    s0=dot(ng,Pn)-dgq; s1=dot(ng,Pf)-dgq
    if s0==0 or s1==0 or (s0>0)==(s1>0): continue
    t=s0/(s0-s1)
    X=tuple(Pn[j]+t*(Pf[j]-Pn[j]) for j in range(3))
    # inside tri g?
    P,Q,R=V[T[g][0]],V[T[g][1]],V[T[g][2]]
    n=plane_n[g]
    an=(abs(n[0]),abs(n[1]),abs(n[2])); dr=an.index(max(an))
    axg=[i for i in range(3) if i!=dr]
    def prg(PP): return (PP[axg[0]],PP[axg[1]])
    def area2(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    p2=prg(X); a2t,b2t,c2t=prg(P),prg(Q),prg(R)
    d1=area2(a2t,b2t,p2); d2=area2(b2t,c2t,p2); d3=area2(c2t,a2t,p2)
    neg=(d1<0)or(d2<0)or(d3<0); pos=(d1>0)or(d2>0)or(d3>0)
    if neg and pos: continue
    print("  SHEET f=%d crosses at t=%.4g pos=(%.9g,%.9g,%.9g) bbov(1276)=%s bbov(1240)=%s"%(
        g,float(t),float(X[0]),float(X[1]),float(X[2]),
        bb_overlap(1276,g),bb_overlap(1240,g)))
    # does the chord (member,g) exist / where does it clip?
    for (fm,_) in members:
        n1=plane_n[fm]; d1q=plane_d[fm]
        dirn=cross(n1,plane_n[g])
        if dirn==(Z,Z,Z): print("    member %d: parallel"%fm); continue
        if not bb_overlap(fm,g): print("    member %d: bb_overlap FALSE"%fm); continue
        AA=V[T[fm][0]]
        adirs=[abs(float(dirn[0])),abs(float(dirn[1])),abs(float(dirn[2]))]
        k=adirs.index(max(adirs))
        i0,i1=[i for i in range(3) if i!=k]
        a11,a12,r1=n1[i0],n1[i1],d1q-n1[k]*AA[k]
        a21,a22,r2=plane_n[g][i0],plane_n[g][i1],plane_d[g]-plane_n[g][k]*AA[k]
        det=a11*a22-a12*a21
        if det==0: print("    member %d: det=0"%fm); continue
        Xc=(r1*a22-a12*r2)/det; Yc=(a11*r2-r1*a21)/det
        P0=[Z,Z,Z]; P0[k]=AA[k]; P0[i0]=Xc; P0[i1]=Yc; P0=tuple(P0)
        tf=line_tri_tparams(P0,dirn,V[T[fm][0]],V[T[fm][1]],V[T[fm][2]])
        tg=line_tri_tparams(P0,dirn,V[T[g][0]],V[T[g][1]],V[T[g][2]])
        if not tf or not tg:
            print("    member %d: clip tf=%s tg=%s"%(fm,tf is not None,tg is not None)); continue
        lo=max(tf[0],tg[0]); hi=min(tf[1],tg[1])
        print("    member %d: interval lo=%.6g hi=%.6g %s"%(
            fm,float(lo),float(hi),"EMPTY" if lo>=hi else "chord exists"))

# ---------- trace g742's arrangement; find cells adjacent to the corner ----------
import functools
segs=[]; sprov=[]
for (f,s) in members:
    Am=V[T[f][0]]; Bm=V[T[f][1]]; Cm=V[T[f][2]]
    segs.append((pr(Am),pr(Bm))); sprov.append(("edge",f))
    segs.append((pr(Bm),pr(Cm))); sprov.append(("edge",f))
    segs.append((pr(Cm),pr(Am))); sprov.append(("edge",f))
    n1=plane_n[f]; d1q=plane_d[f]
    for g in range(nT):
        if g==f or fkey[g] is None or fkey[g]==GKEY: continue
        if not bb_overlap(f,g): continue
        ng=plane_n[g]; dgq=plane_d[g]
        dirn=cross(n1,ng)
        if dirn==(Z,Z,Z): continue
        adirs=[abs(float(dirn[0])),abs(float(dirn[1])),abs(float(dirn[2]))]
        k=adirs.index(max(adirs))
        i0,i1=[i for i in range(3) if i!=k]
        a11,a12,r1=n1[i0],n1[i1],d1q-n1[k]*Am[k]
        a21,a22,r2=ng[i0],ng[i1],dgq-ng[k]*Am[k]
        det=a11*a22-a12*a21
        if det==0: continue
        Xc=(r1*a22-a12*r2)/det; Yc=(a11*r2-r1*a21)/det
        P0=[Z,Z,Z]; P0[k]=Am[k]; P0[i0]=Xc; P0[i1]=Yc; P0=tuple(P0)
        tf=line_tri_tparams(P0,dirn,Am,Bm,Cm)
        tg=line_tri_tparams(P0,dirn,V[T[g][0]],V[T[g][1]],V[T[g][2]])
        if not tf or not tg: continue
        lo=max(tf[0],tg[0]); hi=min(tf[1],tg[1])
        if lo>=hi: continue
        Pa=tuple(P0[j]+lo*dirn[j] for j in range(3))
        Pb=tuple(P0[j]+hi*dirn[j] for j in range(3))
        segs.append((pr(Pa),pr(Pb))); sprov.append(("chord",f,g))

S=len(segs)
pts_on=[[(Z,segs[i][0]),(ONE,segs[i][1])] for i in range(S)]
for i in range(S):
    p,q=segs[i]; r=(q[0]-p[0],q[1]-p[1])
    for j in range(i+1,S):
        a,b=segs[j]; s2=(b[0]-a[0],b[1]-a[1])
        den=r[0]*s2[1]-r[1]*s2[0]
        qp=(a[0]-p[0],a[1]-p[1])
        if den!=0:
            t=(qp[0]*s2[1]-qp[1]*s2[0])/den; u=(qp[0]*r[1]-qp[1]*r[0])/den
            if 0<=t<=1 and 0<=u<=1:
                X=(p[0]+t*r[0],p[1]+t*r[1])
                pts_on[i].append((t,X)); pts_on[j].append((u,X))
        else:
            if qp[0]*r[1]-qp[1]*r[0]!=0: continue
            rr=r[0]*r[0]+r[1]*r[1]
            if rr==0: continue
            for E in (a,b):
                t=((E[0]-p[0])*r[0]+(E[1]-p[1])*r[1])/rr
                if 0<t<1: pts_on[i].append((t,E))
            ss=s2[0]*s2[0]+s2[1]*s2[1]
            if ss==0: continue
            for E in (p,q):
                u=((E[0]-a[0])*s2[0]+(E[1]-a[1])*s2[1])/ss
                if 0<u<1: pts_on[j].append((u,E))
vid={}; verts2=[]
def vof(p):
    if p not in vid: vid[p]=len(verts2); verts2.append(p)
    return vid[p]
adj={}
for i in range(S):
    seq=[]
    for tpar,p in sorted(pts_on[i],key=lambda x:x[0]):
        v=vof(p)
        if not seq or seq[-1]!=v: seq.append(v)
    for k in range(len(seq)-1):
        a,b=seq[k],seq[k+1]
        adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)

vA=vid.get(A2); vB=vid.get(B2)
print("A2 in arrangement:",vA,"B2:",vB,"A-B adjacent:",
      vB in adj.get(vA,set()) if vA is not None else None)

def cmp_dir(v):
    pv=verts2[v]
    def cmp(a,b):
        pa=verts2[a]; pb=verts2[b]
        ax_,ay=pa[0]-pv[0],pa[1]-pv[1]; bx,by=pb[0]-pv[0],pb[1]-pv[1]
        ha=0 if (ay>0 or (ay==0 and ax_>0)) else 1
        hb=0 if (by>0 or (by==0 and bx>0)) else 1
        if ha!=hb: return -1 if ha<hb else 1
        cr=ax_*by-ay*bx
        if cr>0: return -1
        if cr<0: return 1
        return 0
    return cmp
cwprev={}
for v in adj:
    nb=sorted(adj[v],key=functools.cmp_to_key(cmp_dir(v)))
    cwprev[v]={nb[k]:nb[k-1] for k in range(len(nb))}
used=set()
print("cells through A or B:")
for a0 in list(adj.keys()):
    for b0 in list(adj[a0]):
        if (a0,b0) in used: continue
        loop=[]; ca,cb=a0,b0; ok=True
        for _ in range(200000):
            used.add((ca,cb)); loop.append(ca)
            nb=cwprev[cb].get(ca)
            if nb is None: ok=False; break
            ca,cb=cb,nb
            if (ca,cb)==(a0,b0): break
        else: ok=False
        if not ok or len(loop)<3: continue
        if (vA in loop) or (vB in loop):
            s=Z
            for k in range(len(loop)):
                x1,y1=verts2[loop[k]]; x2,y2=verts2[loop[(k+1)%len(loop)]]
                s+=x1*y2-x2*y1
            # classify at a strictly-interior point: centroid of first ear-ish
            print("  loop n=%d s>0=%s hasA=%s hasB=%s AB-consec=%s"%(
                len(loop),s>0,vA in loop,vB in loop,
                any((loop[k]==vA and loop[(k+1)%len(loop)]==vB) or
                    (loop[k]==vB and loop[(k+1)%len(loop)]==vA)
                    for k in range(len(loop)))))
            if s>0 and len(loop)<600:
                pts=[verts2[i] for i in loop]
                cx=sum(p[0] for p in pts)/len(pts); cy=sum(p[1] for p in pts)/len(pts)
                print("    verts near corner:",
                      [("%.8g,%.8g,%.8g"%tuple(float(x) for x in lift(verts2[i])))
                       for i in loop if math.dist(
                           [float(x) for x in lift(verts2[i])],
                           [-16.7107814,7.52182514,-207.42406])<5e-4])

# ---------- classify the big cell at its largest-subtri centroid, fast vs exact ----------
def area2(a,b,c): return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
# retrace to capture the big cell's loop
used2=set(); bigloop=None
for a0 in list(adj.keys()):
    for b0 in list(adj[a0]):
        if (a0,b0) in used2: continue
        loop=[]; ca,cb=a0,b0; ok=True
        for _ in range(200000):
            used2.add((ca,cb)); loop.append(ca)
            nb=cwprev[cb].get(ca)
            if nb is None: ok=False; break
            ca,cb=cb,nb
            if (ca,cb)==(a0,b0): break
        else: ok=False
        if not ok or len(loop)<3: continue
        if vA in loop and vB in loop and len(loop)>10:
            s=Z
            for k in range(len(loop)):
                x1,y1=verts2[loop[k]]; x2,y2=verts2[loop[(k+1)%len(loop)]]
                s+=x1*y2-x2*y1
            if s>0: bigloop=loop
print("bigloop n=%s"%(len(bigloop) if bigloop else None))
if bigloop:
    poly=[verts2[i] for i in bigloop]
    # earclip (same as e1_coord)
    n=len(poly); idx=list(range(n)); tris=[]
    while len(idx)>3:
        found=False; m=len(idx)
        for k in range(m):
            a,b,c=idx[(k-1)%m],idx[k],idx[(k+1)%m]
            if area2(poly[a],poly[b],poly[c])<=0: continue
            ok=True
            for j in idx:
                if j in (a,b,c): continue
                d1=area2(poly[a],poly[b],poly[j]); d2=area2(poly[b],poly[c],poly[j]); d3=area2(poly[c],poly[a],poly[j])
                if d1>=0 and d2>=0 and d3>=0: ok=False; break
            if ok: tris.append((a,b,c)); idx.pop(k); found=True; break
        if not found: break
    if len(idx)==3: tris.append(tuple(idx))
    best=None; bestA=Z
    for t3 in tris:
        A2q=area2(poly[t3[0]],poly[t3[1]],poly[t3[2]])
        if A2q>bestA: bestA=A2q; best=t3
    c2=tuple((poly[best[0]][j]+poly[best[1]][j]+poly[best[2]][j])/3 for j in range(2))
    c3=lift(c2)
    print("big cell centroid=(%.9g,%.9g,%.9g) area2max=%.3g ntris=%d"%(
        float(c3[0]),float(c3[1]),float(c3[2]),float(bestA),len(tris)))
    w1,w2,dg=probe(c3)
    print("FAST probes at centroid: wA=%d wB=%d dg=%s"%(w1,w2,dg))
    # exact winding (slow) at the same probe points
    F=oracle.build_faces(comp)
    c3f=(float(c3[0]),float(c3[1]),float(c3[2]))
    tmax=Fr(eps)/(100*L1N)
    for g in range(nT):
        if fkey[g] is None or fkey[g]==GKEY: continue
        bbg=fbb[g]; mm=1e-6
        if not (bbg[0]-mm<=c3f[0]<=bbg[1]+mm and bbg[2]-mm<=c3f[1]<=bbg[3]+mm
                and bbg[4]-mm<=c3f[2]<=bbg[5]+mm): continue
        ng=plane_n[g]
        num=dot(ng,c3)-plane_d[g]
        if num==0: continue
        if num<0: num=-num
        L1g=abs(ng[0])+abs(ng[1])+abs(ng[2])
        tg=num/(2*L1g*L1N)
        if tg<tmax: tmax=tg
    off=(tmax*N[0],tmax*N[1],tmax*N[2])
    pAq=vadd(c3,off); pBq=vsub(c3,off)
    ea,da=oracle.winding(pAq,F); eb,db=oracle.winding(pBq,F)
    print("EXACT winding at centroid probes: wA=%d(%s) wB=%d(%s)"%(-ea,da,-eb,db))

# ---------- pinch analysis of the big loop ----------
print("bigloop vertex ids:",bigloop)
from collections import Counter
cnt=Counter(bigloop)
rep=[v for v,c in cnt.items() if c>1]
print("repeated vertices:",rep)
def split_pinch(loop):
    seen={}
    for k,v in enumerate(loop):
        if v in seen:
            i=seen[v]
            l1=loop[i:k]; l2=loop[:i]+loop[k:]
            return split_pinch(l1)+split_pinch(l2)
        seen[v]=k
    return [loop] if len(loop)>=3 else []
subs=split_pinch(bigloop)
for sub in subs:
    s=Z
    for k in range(len(sub)):
        x1,y1=verts2[sub[k]]; x2,y2=verts2[sub[(k+1)%len(sub)]]
        s+=x1*y2-x2*y1
    print("  sub n=%d s=%.4g s>0=%s hasA=%s hasB=%s"%(
        len(sub),float(s),s>0,vA in sub,vB in sub))
