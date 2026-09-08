#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, itertools
from pathlib import Path
import m23_hensel_core as h
P=h.P
SCRIPT=Path(__file__).resolve()

def script_sha256():
    return hashlib.sha256(SCRIPT.read_bytes()).hexdigest()

# Exact stage-four obstruction interpolation on the total-degree <=3 Newton grid.
def exponents(nvars,degree):
    out=[]
    def rec(prefix,left,k):
        if k==1:out.append(tuple(prefix+[left]));return
        for v in range(left+1):rec(prefix+[v],left-v,k-1)
    for total in range(degree+1):rec([],total,nvars)
    return out
exps=exponents(4,3);assert len(exps)==35
pts=exps[:]
def evmon(pt,e):
    z=1
    for x,k in zip(pt,e):z=z*pow(x,k,P)%P
    return z
A=[[evmon(pt,e) for e in exps] for pt in pts]
assert len(h.rref(A)[1])==35
vals=[h.stage4_obstruction(pt) for pt in pts]
coeffs=[]
for j in range(12):
    sol,_,_=h.solve_affine(A,[v[j] for v in vals]);assert sol is not None;coeffs.append(sol)
# Factor F4 = 11 L1 L2 L3.
roots=[(16,92,93,51),(22,83,104,50),(134,126,28,172)]
def L(pt,r):return (pt[3]-r[0]-r[1]*pt[0]-r[2]*pt[1]-r[3]*pt[2])%P
def F4(pt):return 11*L(pt,roots[0])*L(pt,roots[1])*L(pt,roots[2])%P
# Determine and certify coordinate multipliers from interpolation nodes.
idx=next(i for i,pt in enumerate(pts) if F4(pt))
fnode=F4(pts[idx]);scalars=[vals[idx][j]*pow(fnode,-1,P)%P for j in range(12)]
for pt,v in zip(pts,vals):assert v==[s*F4(pt)%P for s in scalars]
for pt in [(5,7,11,13),(190,62,203,107),(127,108,33,179),(0,0,1,0),(210,210,210,210)]:
    assert h.stage4_obstruction(pt)==[s*F4(pt)%P for s in scalars]
assert scalars==[112,207,36,189,9,12,119,84,19,82,116,70]
# Pair intersections coincide; count union.
# Differences of root forms are scalar multiples, while any two full hyperplane equations are independent.
d1=[(roots[0][j]-roots[1][j])%P for j in range(4)]
d2=[(roots[0][j]-roots[2][j])%P for j in range(4)]
assert len(h.rref([d1[1:],d2[1:]])[1])==1
p4_count=3*P**3-2*P**2
assert p4_count==28092751
# Exact affine sensitivity at the p4 -> p5 gate. Taylor degree is affine; five spanning points determine it.
span_points=[(0,0,0,16),(0,0,0,22),(1,0,0,108),(0,1,0,109),(0,0,1,67)]
assert len(h.rref([[1,*pt] for pt in span_points])[1])==5
def sensitivity(pt):
    d=h.lift_to_p4(pt);assert d is not None
    u2,z3=d['u2'],d['z3']
    c,_,_=h.next_obstruction_with_kernel(u2,3,z3,[0]*5)
    cols=[]
    for i in range(5):
        g=[0]*5;g[i]=1
        v,_,_=h.next_obstruction_with_kernel(u2,3,z3,g)
        cols.append([(v[j]-c[j])%P for j in range(12)])
    return c,cols
for pt in span_points:
    _,cols=sensitivity(pt);assert all(x==0 for col in cols for x in col)
# Resolve proposed points.
proposed={
 'direct_127_211':(127,108,33,179),
 'roman_carry':(0,0,1,0),
 'original_selected':(190,62,203,107),
}
proposed_residuals={k:h.stage4_obstruction(v) for k,v in proposed.items()}
assert any(proposed_residuals['direct_127_211'])
assert any(proposed_residuals['roman_carry'])
assert proposed_residuals['original_selected']==[0]*12
# Read the exhaustive point list and log.
parser=argparse.ArgumentParser()
ROOT=Path(__file__).resolve().parents[2]
parser.add_argument('--points',default=str(ROOT/'results'/'family_211'/'stage5_zero_points.txt'))
parser.add_argument('--exhaustive-log',default=str(ROOT/'results'/'family_211'/'exhaustive_stage5.log'))
parser.add_argument('--json',default=str(ROOT/'results'/'family_211'/'M23_complete_211_family.json'))
args=parser.parse_args()
points_path=Path(args.points);log_path=Path(args.exhaustive_log)
points=[tuple(map(int,line.split())) for line in points_path.read_text().splitlines() if line.strip()]
assert len(points)==len(set(points))==212
# Every listed point is a p4 survivor and has zero p5 residual. Continue all branches.
records=[]
for pt in points:
    assert h.stage4_obstruction(pt)==[0]*12
    d=h.lift_to_p4(pt);assert d is not None
    FF=h.phi(d['u3']);assert all(v%(P**4)==0 for v in FF)
    c5=h.obs([(v//(P**4))%P for v in FF]);assert c5==[0]*12
    rec=h.adaptive_branch(pt,max_n=6)
    assert rec['highest_precision']==6
    assert rec.get('failure_precision')==7
    assert rec['failure']=='lookahead'
    assert all(item['A_rank']==0 for item in rec['history'])
    assert any(rec['failure_c'])
    records.append({'point':pt,'terminal_residual':rec['failure_c']})
# Affine span and the two normalized equations.
rows=[[1,*pt] for pt in points]
rank=len(h.rref(rows)[1]);assert rank==3
_,basis,_=h.solve_affine(rows,[0]*len(rows));assert len(basis)==2
expected_eqs=[[205,9,200,1,0],[100,156,190,0,1]]
for eq in expected_eqs:
    assert all(sum(eq[j]*row[j] for j in range(5))%P==0 for row in rows)
# Exhaustive log aggregation is an independently generated traversal statement.
log=log_path.read_text()
assert 'TRIPLES 9393931 ROOTS 28092751 ZEROES 212 BAD 0' in log
sha_points=hashlib.sha256(points_path.read_bytes()).hexdigest()
sha_log=hashlib.sha256(log_path.read_bytes()).hexdigest()
terminal_blob=json.dumps(records,sort_keys=True,separators=(',',':')).encode()
result={
 'certificate':'M23_COMPLETE_211_FAMILY',
 'environment':{'script_sha256':script_sha256()},
 'prime':P,
 'stage_two_hyperplane':'alpha0 + 14 alpha1 + 124 alpha2 + 155 alpha3 + 91 alpha4 = 113',
 'stage_four_obstruction':{
   'factorization':'11*L1*L2*L3',
   'linear_factors':[
    'alpha4-(16+92 alpha1+93 alpha2+51 alpha3)',
    'alpha4-(22+83 alpha1+104 alpha2+50 alpha3)',
    'alpha4-(134+126 alpha1+28 alpha2+172 alpha3)'],
   'coordinate_multipliers':scalars,
   'survivor_count':p4_count,
 },
 'stage_five':{
   'kernel_sensitivity_zero_on_complete_survivor_locus':True,
   'survivor_count':len(points),
   'affine_span_equations':[
    '205 + 9 alpha1 + 200 alpha2 + alpha3 = 0',
    '100 + 156 alpha1 + 190 alpha2 + alpha4 = 0'],
 },
 'terminal_lifting':{
   'all_212_lift_through_power':6,
   'all_212_fail_at_power':7,
   'terminal_sensitivity_rank':0,
   'terminal_residuals_nonzero':True,
   'terminal_records_sha256':hashlib.sha256(terminal_blob).hexdigest(),
 },
 'proposed_directions':{k:{'u4':v,'alpha':h.alpha_from_u(v),'stage_four_residual':proposed_residuals[k]} for k,v in proposed.items()},
 'files':{'points_sha256':sha_points,'exhaustive_log_sha256':sha_log},
 'scope':'The complete nondegenerate released m3=179 family has no compatible all-order 211-adic lift. This does not itself construct or exclude rational points on other residue branches or a central idempotent in a separate genuine interior fiber algebra.',
 'pass':True,
}
Path(args.json).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(result['certificate'])
print('PASS',result['pass'])
print('SCRIPT_SHA256',result['environment']['script_sha256'])
print('P4_SURVIVORS',p4_count)
print('P5_SURVIVORS',len(points))
print('ALL_LIFT_THROUGH',6)
print('ALL_FAIL_AT',7)
print('POINTS_SHA256',sha_points)
