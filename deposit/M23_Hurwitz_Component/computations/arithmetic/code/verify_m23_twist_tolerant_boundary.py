#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip, hashlib, json, math, time
from collections import defaultdict
from pathlib import Path

def cycles(p):
    seen=bytearray(len(p)); out=[]
    for a in range(len(p)):
        if seen[a]: continue
        c=[]; x=a
        while not seen[x]: seen[x]=1; c.append(x); x=p[x]
        out.append(c)
    return out
def parity(p): return 1 if (len(p)-len(cycles(p)))%2==0 else -1
def odd_commuter(s,cs,forbidden=frozenset()):
    by=defaultdict(list); forbidden=set(forbidden)
    for c in cs:
        if len(c)>1 and len(c)%2 and not forbidden.intersection(c): by[len(c)].append(c)
    for l in sorted(by):
        if len(by[l])>=2:
            a,b=by[l][:2]; q=list(range(len(s)))
            for x,y in zip(a,b): q[x]=y; q[y]=x
            return tuple(q)
    raise RuntimeError('no odd commuting cycle swap')
def build_base(s,cs,m):
    nmap=list(range(len(s))); t=list(range(len(s)))
    for c in cs:
        l=len(c); mm=m%l
        for k,x in enumerate(c): t[x]=c[(k+mm)%l]; nmap[x]=c[(k*mm)%l]
    return tuple(nmap),tuple(t)
def compose_left(a,c): return tuple(a[c[x]] for x in range(len(a)))

def verify(input_path: Path) -> dict:
    with gzip.open(input_path,'rt') as f: data=json.load(f)
    ops=[tuple(x) for x in data['forward_operators']]; css=[cycles(s) for s in ops]
    fixed=[[c[0] for c in cs if len(c)==1] for cs in css]
    odd=[odd_commuter(s,cs) for s,cs in zip(ops,css)]
    pairs={'B13_fixed_pair':(1,fixed[1][0],fixed[1][1]),'B14_fixed_pair':(2,fixed[2][0],fixed[2][1])}
    units=[m for m in range(1,9240) if math.gcd(m,9240)==1]
    start=time.time(); count=0; samples=[]
    for m in units:
        base=[]; targets=[]
        for idx,(s,cs) in enumerate(zip(ops,css)):
            a,t=build_base(s,cs,m)
            if parity(a)<0: a=compose_left(a,odd[idx])
            assert parity(a)>0 and all(a[s[x]]==t[a[x]] for x in range(len(s)))
            base.append(a); targets.append(t)
        for name,(j,A,B) in pairs.items():
            for source,target in ((A,B),(B,A)):
                ns=list(base); a=list(base[j]); a[source],a[target]=a[target],a[source]; a=tuple(a)
                if parity(a)<0: a=compose_left(a,odd[j])
                assert parity(a)>0 and a[source]==target
                assert all(a[ops[j][x]]==targets[j][a[x]] for x in range(len(a)))
                ns[j]=a
                if len(samples)<4:
                    h=hashlib.sha256()
                    for q in ns:
                        for x in q: h.update(x.to_bytes(2,'little'))
                    samples.append({'m':m,'pair':name,'source':source+1,'target':target+1,'sha256':h.hexdigest()})
                count+=1
    return {
      'certificate':'M23_TWIST_TOLERANT_BOUNDARY_NONEXCLUSION','geometric_monodromy':'A_980','normalizer_in_S980':'S_980',
      'modulus':9240,'unit_count':len(units),'pairs_one_based':{k:[a+1,b+1] for k,(j,a,b) in pairs.items()},
      'verified_nonempty_transporter_cases':count,'expected_cases':len(units)*4,'phi_choice':'identity','all_n_i_even':True,
      'equation_checked':'n_i o s_i = s_i^m o n_i, equivalent to s_i = n_i^-1 s_i^m n_i',
      'sample_witness_hashes':samples,'conclusion':'All multiplier/pair/direction transporter sets are nonempty; peripheral power-conjugacy alone cannot separate either width-one boundary pair.',
      'elapsed_seconds':time.time()-start,'pass':count==len(units)*4,
    }
def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser(); ap.add_argument('--input',type=Path,default=root/'data/orbit_certificate.json.gz'); ap.add_argument('--json',type=Path,default=root/'results/M23_twist_tolerant_boundary_audit.json'); args=ap.parse_args()
    r=verify(args.input); args.json.parent.mkdir(parents=True,exist_ok=True); args.json.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
    print(r['certificate']); print('PASS',r['pass']); print('UNIT_COUNT',r['unit_count']); print('PAIRS',r['pairs_one_based']); print('VERIFIED_NONEMPTY_TRANSPORTER_CASES',r['verified_nonempty_transporter_cases']); print('ALL_N_I_EVEN',r['all_n_i_even']); print('PHI_CHOICE',r['phi_choice'])
    return 0 if r['pass'] else 1
if __name__=='__main__': raise SystemExit(main())
