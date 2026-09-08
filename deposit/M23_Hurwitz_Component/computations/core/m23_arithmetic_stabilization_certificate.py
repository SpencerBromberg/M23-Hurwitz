#!/usr/bin/env python3
"""Exact arithmetic-stabilization certificate for the M23 ten-involution layer.

This certificate proves only the finite Nielsen/frame statements printed in the
article.  The certificate records the finite marking structure; interior rational-point activation belongs to the five-branch
Hurwitz surface or a global ten-sheet Hurwitz subcover.
"""
import os, io, contextlib, importlib.util, math
from sympy.combinatorics import Permutation, PermutationGroup
HERE=os.path.dirname(os.path.abspath(__file__))
BASE=os.path.join(HERE,'base')

def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    with contextlib.redirect_stdout(io.StringIO()): s.loader.exec_module(m)
    return m
h=load(os.path.join(BASE,'m23_huang_ten_2a.py'),'huang')
m=load(os.path.join(BASE,'m23_ten_2a_explicit.py'),'target')
N=23

def hm(T,i,d):
    T=list(T);a,b=T[i],T[i+1]
    if d==1: T[i],T[i+1]=b,h.compose(h.compose(h.inv(b),a),b)
    else: T[i],T[i+1]=h.compose(h.compose(a,b),h.inv(a)),a
    return tuple(T)

def order(p): return int(Permutation(list(p)).order())
def canon(T):
    R=len(T); sig=[sum((1<<i) for i,p in enumerate(T) if p[v]==v) for v in range(N)]; mn=min(sig); best=None
    for root in [v for v,z in enumerate(sig) if z==mn]:
        lab=[-1]*N; seq=[root]; lab[root]=0; k=0
        while k<len(seq):
            u=seq[k];k+=1
            for p in T:
                w=p[u]
                if lab[w]<0: lab[w]=len(seq);seq.append(w)
        if len(seq)!=N: continue
        key=tuple(lab[p[u]] for p in T for u in seq)
        if best is None or key<best: best=key
    return best

def dec(key,R): return tuple(tuple(key[i*N:(i+1)*N]) for i in range(R))
def cmove(key,i,d,R): return canon(hm(dec(key,R),i,d))
def sig_h(x): return h.conj(h.h,h.inv(x))
def closure(gens):
    e=h.I(); gs=list(gens)+[h.inv(g) for g in gens]; S={e}; todo=[e]
    while todo:
        a=todo.pop()
        for g in gs:
            b=h.compose(a,g)
            if b not in S:S.add(b);todo.append(b)
    return S

def tp(P): return tuple(P(i) for i in range(N))
def conj_inv(p,x): return h.compose(h.compose(h.inv(x),p),x) # x^-1 p x

# A. A short local B4 word on the released 23A factorization.
A=tuple(h.T[2:6])
local_word=[(2,-1),(1,-1),(0,1),(0,1),(1,-1),(0,1),(2,1),(2,1),(1,-1)]
for i,d in local_word: A=hm(A,i,d)
assert h.prod(A)==h.g2
B=tuple(sig_h(x) for x in reversed(A))
R=(h.T[0],h.T[1])+A+B
assert h.prod(R)==h.I() and h.prod(R[:2])==h.g1 and h.prod(R[2:6])==h.g2 and h.prod(R[6:])==h.g3
iota=[0,1,9,8,7,6,5,4,3,2]
assert all(sig_h(R[i])==R[iota[i]] for i in range(10))

# B. The already released twelve-move interleaving, followed by pair fusion.
inter=[(8,-1),(7,1),(6,-1),(5,-1),(4,1),(3,-1),(8,-1),(7,-1),(6,-1),(5,-1),(8,-1),(7,1)]
U=R
for i,d in inter: U=hm(U,i,d)
P=tuple(h.compose(U[2*i],U[2*i+1]) for i in range(5))
orders=[order(x) for x in P]
assert orders==[2,2,3,2,2] and h.prod(P)==h.I()
Gorder=PermutationGroup([h.sp(x) for x in P]).order(); assert Gorder==10200960

# C. Exact nine-move B5 membership word from the printed target seed to P.
h5_word=[(3,1),(2,-1),(3,1),(2,-1),(1,1),(0,-1),(1,1),(0,1),(2,1)]
s=canon(tuple(m.seed))
for i,d in h5_word:s=cmove(s,i,d,5)
assert s==canon(P)

# D. The true Huang descent involution and a D10<F20 marking on the SAME refinement.
hA=[order(h.compose(h.h,a)) for a in A]
assert hA==[4,3,5,3]
a=A[2]
r=h.compose(h.h,a)
assert order(r)==5 and order(a)==2 and order(h.h)==2
# D10 inversion by h: h^-1 r h = r^-1 (h=h^-1).
assert conj_inv(r,h.h)==h.inv(r)

G=PermutationGroup([h.sp(h.a),h.sp(h.b)])
C=G.centralizer(h.sp(h.h)); assert C.order()==2688
r2=h.compose(r,r)
solutions=[]
for X in C.generate_schreier_sims():
    x=tp(X)
    if h.compose(x,x)==h.h and conj_inv(r,x)==r2:
        solutions.append(x)
assert len(solutions)==3
x=solutions[0]
assert order(x)==4 and h.compose(x,x)==h.h and conj_inv(h.h,x)==h.h and conj_inv(r,x)==r2
assert len(closure([r,x]))==20

# E. The five rotations plus five reflected markings form the intrinsic D10 packet.
packet=set(); rk=h.I()
for _ in range(5):
    packet.add(rk);packet.add(h.compose(rk,h.h));rk=h.compose(rk,r)
assert len(packet)==10
assert all(conj_inv(y,x) in packet for y in packet)
# right multiplication by h is the fixed-point-free dual pairing
assert all(h.compose(y,h.h)!=y for y in packet)

# F. Arithmetic multiplier decomposition across the stabilized 5A/3A data.
units=[1,7,11,13,17,19,23,29]
C4=[];z=1
for _ in range(4): C4.append(z);z=(z*7)%30
assert C4==[1,7,19,13] and (11*11)%30==1
assert sorted(set(C4+[(11*z)%30 for z in C4]))==units
assert 7%5==2 and 7%3==1
assert 11%5==1 and 11%3==2
# On the unique order-3 fused pair, swapping the two involution factors inverts the product.
a3,b3=U[4],U[5]
assert order(a3)==order(b3)==2
assert h.compose(a3,b3)==P[2]
assert h.compose(b3,a3)==h.inv(P[2])
assert h.compose(P[2],P[2])==h.inv(P[2]) # 11 == 2 mod 3
# On all order-2 fused entries, the 11th power is the entry itself.
for j in (0,1,3,4):
    assert order(P[j])==2

# G. Raw duals still collapse after forgetting the rigidified markings.
Rdual=tuple(sig_h(z) for z in R)
assert Rdual==tuple(h.conj(h.h,z) for z in R)

print('LOCAL_B4_WORD_LENGTH',len(local_word))
print('LOCAL_B4_WORD',' '.join('q%d%s'%(i+1,'^-1' if d<0 else '') for i,d in local_word))
print('Q_SYMMETRIC_REFINEMENT',True)
print('PAIR_FUSION_TYPES',orders)
print('PAIR_FUSION_GENERATED_ORDER',Gorder)
print('TARGET_H5_MEMBERSHIP_WORD_LENGTH',len(h5_word))
print('TARGET_H5_MEMBERSHIP_WORD',' '.join('q%d%s'%(i+1,'^-1' if d<0 else '') for i,d in h5_word))
print('HUANG_HA_ORDERS',hA)
print('D10_LEAF_INDEX_1BASED',3)
print('D10_R_ORDER',order(r))
print('CENTRALIZER_H_ORDER',C.order())
print('MULTIPLIER2_LIFTS',len(solutions))
print('F20_X_ORDER',order(x))
print('F20_X_SQUARE_EQUALS_H',h.compose(x,x)==h.h)
print('F20_R_TO_R2',conj_inv(r,x)==r2)
print('F20_SUBGROUP_ORDER',len(closure([r,x])))
print('LOCAL_FIVE_PLUS_DUAL_PACKET_SIZE',len(packet))
print('UNITS_MOD_30',units)
print('C4_GENERATOR_7_ORBIT',C4)
print('MULTIPLIER_7_MOD5_MOD3',(7%5,7%3))
print('MULTIPLIER_11_MOD5_MOD3',(11%5,11%3))
print('UNIQUE_3A_PAIR_SWAP_GIVES_INVERSE',h.compose(b3,a3)==h.inv(P[2]))
print('RAW_DUAL_COLLAPSE_AFTER_FORGETTING_MARKINGS',True)
print('ALL_ARITHMETIC_STABILIZATION_CHECKS_PASS',True)
