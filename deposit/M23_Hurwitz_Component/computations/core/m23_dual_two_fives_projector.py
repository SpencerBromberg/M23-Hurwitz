#!/usr/bin/env python3
from pathlib import Path
import json
import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup
ROOT=Path(__file__).resolve().parent
D=[5619,5707,10252,11406,13499,15797,17032,18695,20064,21190]
R=[1,2,4,6,7,8,5,0,9,3]
K=[8,9,3,2,6,7,4,5,0,1]
cycle_p0=[21190,11406,17032,15797,20064]
cycle_p1=[5707,10252,13499,18695,5619]
pairs=list(zip(cycle_p0,cycle_p1))

def comp(p,q): return [p[q[i]] for i in range(len(p))]
def inv(p):
    z=[0]*len(p)
    for i,j in enumerate(p): z[j]=i
    return z
def ppow(p,e):
    z=list(range(len(p))); b=p[:]
    while e:
        if e&1: z=comp(b,z)
        b=comp(b,b); e//=2
    return z
def order(p):
    z=list(range(len(p)))
    for n in range(1,100):
        z=comp(p,z)
        if z==list(range(len(p))): return n
    raise RuntimeError('order bound exceeded')
def matrix_of_perm(p):
    M=sp.zeros(len(p))
    for i,j in enumerate(p): M[j,i]=1
    return M

assert ppow(R,5)==list(range(10))
assert ppow(K,2)==list(range(10)) and all(K[i]!=i for i in range(10))
assert comp(R,K)==comp(K,R)
coord={}
for n,s in enumerate(cycle_p0): coord[s]=(n,0)
for n,s in enumerate(cycle_p1): coord[s]=(n,1)
assert set(coord)==set(D)
pos={s:i for i,s in enumerate(D)}
for s,(n,p) in coord.items():
    assert coord[D[R[pos[s]]]]==((n+1)%5,p)
    assert coord[D[K[pos[s]]]]==(n,1-p)
fibers={n:sorted(s for s,c in coord.items() if c[0]==n) for n in range(5)}
assert [len(fibers[n]) for n in range(5)]==[2,2,2,2,2]
I=sp.eye(10); MK=matrix_of_perm(K); MR=matrix_of_perm(R)
Pp=(I+MK)/2; Pm=(I-MK)/2
assert Pp*Pp==Pp and Pm*Pm==Pm and Pp*Pm==sp.zeros(10) and Pp+Pm==I
assert Pp.rank()==5 and Pm.rank()==5 and MR*MK==MK*MR
X=[None]*10
for s,(n,p) in coord.items():
    target=next(t for t,c in coord.items() if c==((2*n)%5,p))
    X[pos[s]]=pos[target]
assert order(X)==4 and comp(X,K)==comp(K,X)
assert comp(X,comp(R,inv(X)))==ppow(R,2)
G40=PermutationGroup([Permutation(R),Permutation(K),Permutation(X)])
G20=PermutationGroup([Permutation([(i+1)%5 for i in range(5)]),Permutation([(2*i)%5 for i in range(5)])])
assert G40.order()==40 and G20.order()==20

t,y,v=sp.symbols('t y v')
beta10=sp.expand(1-(1-t**2)**5); beta5=sp.expand(1-(1-v)**5)
F=sp.expand(beta10-y)
assert sp.degree(F,t)==10 and sp.expand(beta10.subs(t,-t)-beta10)==0
assert sp.expand(beta10-beta5.subs(v,t**2))==0
disc=sp.factor(sp.discriminant(F,t))
assert disc==10_000_000_000*y*(y-1)**8
res={
 'certificate':'M23_DUAL_TWO_FIVES_PROJECTOR',
 'D10_states_zero_based':D,
 'two_five_cycles':[cycle_p0,cycle_p1],
 'projected_key_fibers':fibers,
 'fiber_profile':[2,2,2,2,2],
 'kappa_pairs':pairs,
 'ledger_action':{'q1':'(n,p)->(n+1,p)','kappa':'(n,p)->(n,p+1)','kappa_free':True,'commute':True},
 'projectors':{'P_plus_rank':5,'P_minus_rank':5,'idempotent':True,'orthogonal':True},
 'multiplier2':{'X_order':4,'X_commutes_kappa':True,'X_q1_Xinv':'q1^2','ledger_group_order':40,'quotient_order_mod_kappa':20,'quotient_type':'F20'},
 'beta10_model':{'beta10':str(beta10),'kappa':'t -> -t','quotient_coordinate':'v=t^2','beta5':str(beta5),'discriminant':str(disc),'etale_open':'y*(y-1) != 0'},
 'scope':{'rigidified_two_five_ledger_certified':True}
}
(ROOT/'M23_DUAL_TWO_FIVES_PROJECTOR.json').write_text(json.dumps(res,indent=2,sort_keys=True)+'\n')
print('D10_SIZE',len(D))
print('Q1_CYCLE_TYPE 5^2')
print('PROJECTED_KEY_MULTIPLICITIES 2 2 2 2 2')
print('KAPPA_FREE True')
print('KAPPA_COMMUTES_Q1 True')
print('P_PLUS_RANK',Pp.rank())
print('P_MINUS_RANK',Pm.rank())
print('LEDGER_GROUP_ORDER',G40.order())
print('QUOTIENT_GROUP_ORDER',G20.order())
print('BETA10',beta10)
print('BETA10_DISCRIMINANT',disc)
print('ALL_DUAL_TWO_FIVES_PROJECTOR_CHECKS_PASS True')
