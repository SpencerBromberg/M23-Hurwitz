import re, os, importlib.util, contextlib, io
from sympy.combinatorics import PermutationGroup
HERE=os.path.dirname(os.path.abspath(__file__))

def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
 with contextlib.redirect_stdout(io.StringIO()):s.loader.exec_module(m)
 return m
m=load(os.path.join(HERE,'m23_ten_2a_explicit.py'),'m')
h=load(os.path.join(HERE,'m23_huang_ten_2a.py'),'h')
N=23

def parse_word(path):
 z=open(path).read().strip().split(); out=[]
 for t in z:
  neg=t.endswith('^-1'); q=int(t[1:].split('^')[0]); out.append((q-1,-1 if neg else 1))
 return out

def hm(T,i,d):
 T=list(T);a,b=T[i],T[i+1]
 if d==1:T[i],T[i+1]=b,h.compose(h.compose(h.inv(b),a),b)
 else:T[i],T[i+1]=h.compose(h.compose(a,b),h.inv(a)),a
 return tuple(T)

def ctype_order(p): return int(h.Permutation(list(p)).order())
def canon(T):
 R=len(T); sig=[sum((1<<i) for i,p in enumerate(T) if p[v]==v) for v in range(N)]; mn=min(sig); best=None;bestlab=None
 for root in [v for v,z in enumerate(sig) if z==mn]:
  lab=[-1]*N;order=[root];lab[root]=0;k=0
  while k<len(order):
   u=order[k];k+=1
   for p in T:
    w=p[u]
    if lab[w]<0:lab[w]=len(order);order.append(w)
  if len(order)!=N:continue
  key=tuple(lab[p[u]] for p in T for u in order)
  if best is None or key<best:best,bestlab=key,tuple(lab)
 return best,bestlab

def sig(x): return h.conj(h.h,h.inv(x))

# A. Exact braid from target 2A^10 refinement to a refinement of released Huang triple.
T=tuple(m.T)
for i,d in parse_word(os.path.join(HERE,'target_to_huang_braid_word.txt')):T=hm(T,i,d)
blocks=(h.prod(T[:2]),h.prod(T[2:6]),h.prod(T[6:]))
HK,HL=canon((h.g1,h.g2,h.g3)); BK,BL=canon(blocks)
assert HK==BK
invHL=[0]*N
for v,z in enumerate(HL):invHL[z]=v
c=tuple(invHL[BL[v]] for v in range(N))
def conjc(x): return h.conj(c,x)
assert tuple(conjc(x) for x in blocks)==(h.g1,h.g2,h.g3)
print('TARGET_TO_HUANG_BRAID_LENGTH',len(parse_word(os.path.join(HERE,'target_to_huang_braid_word.txt'))))
print('TARGET_TO_HUANG_CONTRACTION_EXACT',True)

# B. Q-symmetric local Huang refinement: modify one 23A factorization and define conjugate 23B block.
A=tuple(h.T[2:6])
for qi,d in parse_word(os.path.join(HERE,'qcompatible_aword.txt')):
 # tokens q3,q4,q5 are internal moves on the four A entries
 A=hm(A,qi-2,d)
B=tuple(sig(x) for x in reversed(A))
R=(h.T[0],h.T[1])+A+B
assert h.prod(R)==h.I()
assert h.prod(R[:2])==h.g1 and h.prod(R[2:6])==h.g2 and h.prod(R[6:])==h.g3
assert sig(R[0])==R[0] and sig(R[1])==R[1]
iota=[0,1,9,8,7,6,5,4,3,2]
assert all(sig(R[i])==R[iota[i]] for i in range(10))
print('Q_SYMMETRIC_HUANG_REFINEMENT',True)
print('GALOIS_LABEL_INVOLUTION','(3 10)(4 9)(5 8)(6 7)')

# C. Interleave conjugate labels; consecutive pairs fuse to target class vector.
U=R
for i,d in parse_word(os.path.join(HERE,'qcompatible_interleave_word.txt')):U=hm(U,i,d)
P=tuple(h.compose(U[2*i],U[2*i+1]) for i in range(5))
assert h.prod(P)==h.I()
orders=[ctype_order(x) for x in P]
assert sorted(orders)==[2,2,2,2,3]
GO=PermutationGroup([h.sp(x) for x in P]).order(); assert GO==10200960
print('PAIR_FUSION_TYPES',orders)
print('PAIR_FUSION_PRODUCT_ONE',True)
print('PAIR_FUSION_GENERATED_ORDER',GO)

# D. Exact target H5 component membership via the stored 8-move full-B5 word.
seed=tuple(m.seed)
# canonical-state move: decode/relabel after every Hurwitz move, exactly the quotient used in the orbit verification
def dec(key,R): return tuple(tuple(key[i*N:(i+1)*N]) for i in range(R))
def ckey(T): return canon(T)[0]
def cmove(key,i,d,R): return ckey(hm(dec(key,R),i,d))
s=ckey(seed)
for i,d in parse_word(os.path.join(HERE,'h5_component_word.txt')):s=cmove(s,i,d,5)
assert s==ckey(P)
print('TARGET_H5_COMPONENT_WORD_LENGTH',len(parse_word(os.path.join(HERE,'h5_component_word.txt'))))
print('TARGET_H5_COMPONENT_MEMBERSHIP',True)
# useful exact order-11 collision present in this representative
print('P3P4_ORDER',ctype_order(h.compose(P[2],P[3])))
print('ALL_CHECKS_PASS',True)
