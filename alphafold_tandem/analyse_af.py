import json,glob,os,sys
import numpy as np
ROOT='/tmp/claude-1000/-home-robson-PetaChem/05591632-d63e-4446-a0e3-108a439ad793/scratchpad/afruns/AF_runs'
AU=2.5417464519
def ang(u,v): return np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1)))
def kabsch(P,Q):
    pc,qc=P.mean(0),Q.mean(0); H=(P-pc).T@(Q-qc)
    U,_,Vt=np.linalg.svd(H); d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T; return R,qc-R@pc
# crystal reference
ca={};cr2={}
for l in open('/home/robson/PetaChem/venus_dimer.pdb'):
    if not l.startswith(('ATOM','HETATM')): continue
    c=l[21];nm=l[12:16].strip();x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    if nm=='CA': ca.setdefault(c,{})[int(l[22:26])]=x
    if l[17:20].strip()=='CR2': cr2.setdefault(c,{})[nm]=x
IMID=['N2','CA2','C2','N3','O2']
def laxis(d): return np.mean([d[a] for a in IMID],axis=0)-d['OH']
def read(fn):
    at={};caa={}
    for l in open(fn):
        if l.startswith('ATOM'):
            p=l.split(); k=(int(p[8]),p[3]); v=np.array([float(p[10]),float(p[11]),float(p[12])])
            at[k]=v
            if p[3]=='CA': caa[int(p[8])]=v
    return at,caa
rows=[]
for d in sorted(glob.glob(f'{ROOT}/x_*')):
    jr=glob.glob(f'{d}/*job_request.json')
    r=json.load(open(jr[0])); r=r[0] if isinstance(r,list) else r
    n=len(r['sequences'][0]['proteinChain']['sequence']); seed=r['modelSeeds'][0]
    tyr2 = 331 if n==495 else 339
    b2off = 265 if n==495 else 273
    for k in range(5):
        cif=f'{d}/{os.path.basename(d)[2:]}_model_{k}.cif'.replace('x_','')
        cif=glob.glob(f'{d}/*_model_{k}.cif')
        if not cif: continue
        at,caa=read(cif[0])
        # delta, Cusick: 90 - theta/2 on OH->CB of the two Tyr precursors
        u=at[(67,'CB')]-at[(67,'OH')]; v=at[(tyr2,'CB')]-at[(tyr2,'OH')]
        theta=ang(u,v); delta=90-theta/2
        # grafted alpha + register RMSD
        out=[]
        for off,ch in ((1,'A'),(b2off,'B')):
            pr=[(a,a-off) for a in caa if 0<=a-off<=230 and (a-off) in ca[ch]]
            P=np.array([caa[a] for a,_ in pr]); Q=np.array([ca[ch][b] for _,b in pr])
            R,t=kabsch(Q,P); rms=float(np.sqrt(((Q@R.T+t-P)**2).sum(1).mean()))
            out.append((rms,laxis({kk:(w@R.T+t) for kk,w in cr2[ch].items()}),R,t,P,Q))
        alpha=ang(out[0][1],out[1][1])
        # register: fit barrel1, carry barrel2
        R1,t1=kabsch(out[0][5],out[0][4])
        cross=float(np.sqrt(((out[1][5]@R1.T+t1-out[1][4])**2).sum(1).mean()))
        # PAE, structured cores only
        fd=glob.glob(f'{d}/*full_data_{k}.json')
        pae=np.nan
        if fd:
            P_=np.array(json.load(open(fd[0]))['pae'],float)
            A=slice(0,231); B=slice(b2off-1,b2off+230)
            pae=float(np.concatenate([P_[A,B].ravel(),P_[B,A].ravel()]).mean())
        rows.append(dict(n=n,seed=seed,model=k,delta=delta,theta=theta,alpha=alpha,
                         reg=cross,fit=out[0][0],pae=pae))
np.save('/tmp/claude-1000/-home-robson-PetaChem/05591632-d63e-4446-a0e3-108a439ad793/scratchpad/af_rows.npy',rows,allow_pickle=True)
import collections
for n in (495,511):
    r=[x for x in rows if x['n']==n]
    if not r: continue
    D=np.array([x['delta'] for x in r]); A=np.array([x['alpha'] for x in r])
    G=np.array([x['reg'] for x in r]); PA=np.array([x['pae'] for x in r])
    F=np.array([x['fit'] for x in r])
    seeds=sorted(set(x['seed'] for x in r))
    print(f"=== {n}-aa construct: {len(seeds)} seeds x 5 models = {len(r)} structures ===")
    print(f"   barrel fold RMSD vs crystal : {F.mean():.2f} A")
    print(f"   inter-domain PAE (cores)    : {PA.mean():.1f} A")
    print(f"   delta (Cusick)              : {D.mean():.1f} +- {D.std():.1f}   range {D.min():.1f}-{D.max():.1f}")
    print(f"   grafted alpha               : {A.mean():.1f} +- {A.std():.1f}   range {A.min():.1f}-{A.max():.1f}")
    print(f"   register RMSD vs crystal    : {G.mean():.1f} +- {G.std():.1f}   range {G.min():.1f}-{G.max():.1f}")
    # between vs within seed
    bys=collections.defaultdict(list)
    for x in r: bys[x['seed']].append(x['delta'])
    mu=np.array([np.mean(v) for v in bys.values()]); wi=np.array([np.std(v) for v in bys.values()])
    print(f"   delta: between-seed sd {mu.std():.1f}   mean within-seed sd {wi.mean():.1f}")
    print()
