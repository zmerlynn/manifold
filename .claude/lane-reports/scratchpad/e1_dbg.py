#!/usr/bin/env python3
# Debug the residual micro-corner opens: rebuild the involved plane groups'
# arrangements and dump the neighborhood structure + cross-plane identity
# presence + sub-edge presence per open edge.
import sys, math, time
from fractions import Fraction as Fr
import importlib.util
spec = importlib.util.spec_from_file_location("oracle", "/tmp/oracle-asm/oracle.py")
oracle = importlib.util.module_from_spec(spec); spec.loader.exec_module(oracle)
parse=oracle.parse; R3=oracle.R3; vsub=oracle.vsub; vadd=oracle.vadd
smul=oracle.smul; cross=oracle.cross; dot=oracle.dot
Z=Fr(0); ONE=Fr(1)

comps=parse(); comp=max(comps,key=lambda c:len(c["tris"])); eps=comp["eps"]
V=[R3(v) for v in comp["verts"]]; T=comp["tris"]; nT=len(T)

plane_n=[None]*nT; plane_d=[None]*nT; fkey=[None]*nT; fsgn=[0]*nT
groups={}
gorder=[]
for f in range(nT):
    v0,v1,v2,_=T[f]; A,B,C=V[v0],V[v1],V[v2]
    n=cross(vsub(B,A),vsub(C,A)); d=dot(n,A)
    plane_n[f]=n; plane_d[f]=d
    lead=None
    for i in range(3):
        if n[i]!=0: lead=i; break
    if lead is None: continue
    q=n[lead]; s=1 if q>0 else -1
    key=(n[0]/q,n[1]/q,n[2]/q,d/q)
    fkey[f]=key; fsgn[f]=s
    if key not in groups:
        groups[key]=dict(N=(key[0],key[1],key[2]),d=key[3],members=[])
        gorder.append(key)
    groups[key]["members"].append((f,s))

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

# args: group indices + focus point
GIDX=[int(x) for x in sys.argv[1].split(",")]
FP=tuple(float(x) for x in sys.argv[2].split(","))
RAD=float(sys.argv[3]) if len(sys.argv)>3 else 5e-4

gkeys=[gorder[i] for i in GIDX]
print("focus groups:")
for gi,key in zip(GIDX,gkeys):
    G=groups[key]
    print(" g%d nmem=%d N=(%.6g,%.6g,%.6g) d=%.9g members=%s"%(
        gi,len(G["members"]),float(key[0]),float(key[1]),float(key[2]),
        float(key[3]),[f for f,_ in G["members"]][:8]))

# pairwise plane distances at the focus point
print("plane signed dist at focus (in eps):")
for gi,key in zip(GIDX,gkeys):
    n=key[:3]; d=key[3]
    num=float(n[0])*FP[0]+float(n[1])*FP[1]+float(n[2])*FP[2]-float(d)
    nl=math.sqrt(float(n[0])**2+float(n[1])**2+float(n[2])**2)
    print("  g%d dist=%.4g eps"%(gi,num/nl/eps))

def build_group(key):
    G=groups[key]; N=G["N"]; dpl=G["d"]; members=G["members"]
    anf=(abs(float(N[0])),abs(float(N[1])),abs(float(N[2])))
    drop=anf.index(max(anf)); ax=[i for i in range(3) if i!=drop]
    def pr(P): return (P[ax[0]],P[ax[1]])
    nk=N[drop]
    def lift(p2):
        P=[Z,Z,Z]; P[ax[0]]=p2[0]; P[ax[1]]=p2[1]
        P[drop]=(dpl-N[ax[0]]*p2[0]-N[ax[1]]*p2[1])/nk
        return tuple(P)
    segs=[]; prov=[]
    for (f,s) in members:
        A=V[T[f][0]]; B=V[T[f][1]]; C=V[T[f][2]]
        segs.append((pr(A),pr(B))); prov.append(("edge",f))
        segs.append((pr(B),pr(C))); prov.append(("edge",f))
        segs.append((pr(C),pr(A))); prov.append(("edge",f))
        n=plane_n[f]; d=plane_d[f]
        for g in range(nT):
            if g==f or fkey[g] is None or fkey[g]==key: continue
            if not bb_overlap(f,g): continue
            ng=plane_n[g]; dg=plane_d[g]
            dirn=cross(n,ng)
            if dirn==(Z,Z,Z): continue
            adirs=[abs(float(dirn[0])),abs(float(dirn[1])),abs(float(dirn[2]))]
            k=adirs.index(max(adirs))
            i0,i1=[i for i in range(3) if i!=k]
            a11,a12,r1=n[i0],n[i1],d-n[k]*A[k]
            a21,a22,r2=ng[i0],ng[i1],dg-ng[k]*A[k]
            det=a11*a22-a12*a21
            if det==0: continue
            X=(r1*a22-a12*r2)/det; Y=(a11*r2-r1*a21)/det
            P0=[Z,Z,Z]; P0[k]=A[k]; P0[i0]=X; P0[i1]=Y; P0=tuple(P0)
            tf=line_tri_tparams(P0,dirn,A,B,C)
            tg=line_tri_tparams(P0,dirn,V[T[g][0]],V[T[g][1]],V[T[g][2]])
            if not tf or not tg: continue
            lo=max(tf[0],tg[0]); hi=min(tf[1],tg[1])
            if lo>=hi: continue
            Pa=tuple(P0[j]+lo*dirn[j] for j in range(3))
            Pb=tuple(P0[j]+hi*dirn[j] for j in range(3))
            segs.append((pr(Pa),pr(Pb))); prov.append(("chord",f,g))
    # arrangement with per-vertex provenance of the segments through it
    S=len(segs)
    pts_on=[[(Z,segs[i][0]),(ONE,segs[i][1])] for i in range(S)]
    for i in range(S):
        p,q=segs[i]; r=(q[0]-p[0],q[1]-p[1])
        for j in range(i+1,S):
            a,b=segs[j]; s2=(b[0]-a[0],b[1]-a[1])
            den=r[0]*s2[1]-r[1]*s2[0]
            qp=(a[0]-p[0],a[1]-p[1])
            if den!=0:
                t=(qp[0]*s2[1]-qp[1]*s2[0])/den; u=(qp[0]*r[0]+0,0) # placeholder
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
    vid={}; verts=[]
    def vof(p):
        if p not in vid: vid[p]=len(verts); verts.append(p)
        return vid[p]
    adj={}; vsegs={}
    for i in range(S):
        seq=[]
        for tpar,p in sorted(pts_on[i],key=lambda x:x[0]):
            v=vof(p)
            if not seq or seq[-1]!=v: seq.append(v)
        for v in seq: vsegs.setdefault(v,set()).add(i)
        for k in range(len(seq)-1):
            a,b=seq[k],seq[k+1]
            adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)
    return dict(segs=segs,prov=prov,verts=verts,adj=adj,vsegs=vsegs,
                lift=lift,pr=pr)

built={gi:build_group(k) for gi,k in zip(GIDX,gkeys)}

# neighborhood dump
print("\nneighborhood (r=%.1e) around (%.9g,%.9g,%.9g):"%(RAD,*FP))
id2g={}
for gi in GIDX:
    B=built[gi]
    for vi,p2 in enumerate(B["verts"]):
        p3=B["lift"](p2)
        p3f=(float(p3[0]),float(p3[1]),float(p3[2]))
        if math.dist(p3f,FP)>RAD: continue
        provs=[B["prov"][s] for s in B["vsegs"].get(vi,())]
        id2g.setdefault(p3,[]).append((gi,vi))
        print(" g%d v%d (%.10g,%.10g,%.10g) deg=%d segs=%s"%(
            gi,vi,p3f[0],p3f[1],p3f[2],len(B["adj"].get(vi,())),provs[:6]))
print("\ncross-plane identity presence (same exact 3D value):")
for p3,lst in id2g.items():
    if len(lst)>1:
        print("  shared by %s at (%.10g,%.10g,%.10g)"%(lst,float(p3[0]),float(p3[1]),float(p3[2])))
