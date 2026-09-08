#!/usr/bin/env python3
"""Exact algebraic and finite certificate for the two six-branch collision tails."""
from __future__ import annotations
from collections import Counter, deque
from itertools import permutations
from pathlib import Path
import re
import sympy as sp

N=23

def cyc(s):
 p=list(range(N))
 for grp in re.findall(r'\(([^()]*)\)',s):
  c=[int(x)-1 for x in re.split(r'[, ]+',grp.strip()) if x]
  for a,b in zip(c,c[1:]+c[:1]): p[a]=b
 return tuple(p)
def comp(a,b): return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
 r=[0]*len(a)
 for i,j in enumerate(a): r[j]=i
 return tuple(r)
def orbits(gens,n):
 unseen=set(range(n)); out=[]
 while unseen:
  root=min(unseen); seen={root}; q=deque([root])
  while q:
   x=q.popleft()
   for g in gens:
    y=g[x]
    if y not in seen: seen.add(y);q.append(y)
  out.append(tuple(sorted(seen))); unseen-=seen
 return sorted(out,key=lambda z:(-len(z),z))
def profile_on(p,orb):
 loc={x:i for i,x in enumerate(orb)}; q=tuple(loc[p[x]] for x in orb)
 seen=set(); L=[]
 for i in range(len(q)):
  if i in seen: continue
  x=i;l=0
  while x not in seen: seen.add(x);l+=1;x=q[x]
  L.append(l)
 return tuple(sorted(L,reverse=True))
def restricted_triples(a,b):
 c=inv(comp(a,b)); out=[]
 for o in orbits([a,b],N): out.append((len(o),profile_on(a,o),profile_on(b,o),profile_on(c,o)))
 return out

u=cyc('(1 15 6)(2 18 16)(4 20 9)(5 8 22)(7 12 21)(10 14 19)')
v=cyc('(1 6 16)(2 15 18)(4 7 21)(8 22 13)(9 20 12)(10 17 19)')
h4=cyc('(2 7)(4 13)(6 8)(9 22)(10 11)(14 17)(15 20)(16 18)')
h5=cyc('(1 20)(3 7)(4 18)(6 19)(9 22)(10 23)(13 17)(14 16)')
R33=restricted_triples(u,v); R22=restricted_triples(h4,h5)
print('3A3A_ORBIT_DEGREES',[r[0] for r in R33])
print('2A2A_ORBIT_DEGREES',[r[0] for r in R22])
print('3A3A_RESTRICTED_TRIPLES',R33)
print('2A2A_RESTRICTED_TRIPLES',R22)

x=sp.symbols('x')
models={
 'alpha6':((12*x**2+12*x-1)**3,(12*x**2-12*x-1)**3),
 'alpha4':(4*x*(x-2)**3,(2*x**2+10*x-1)**2),
 'beta6':((x**3+1)**2,4*x**3),
 'beta3':(4*x**3-3*x+1,sp.Integer(2)),
 'beta2':(x**2,x**2+1),
}
def root_mults(poly,d):
 p=sp.Poly(sp.expand(poly),x); out=[]
 if p.is_zero: raise ValueError('zero fiber polynomial')
 for factor,e in sp.factor_list(p.as_expr())[1]:
  out += [int(e)]*sp.Poly(factor,x).degree()
 out += [d-p.degree()] if p.degree()<d else []
 return tuple(sorted(out,reverse=True))
def passport(num,den):
 num=sp.expand(num);den=sp.expand(den);d=max(sp.degree(num,x),sp.degree(den,x))
 return (root_mults(num,d),root_mults(num-den,d),root_mults(den,d))
passports={k:passport(*v) for k,v in models.items()}
for k,vv in passports.items(): print(k.upper()+'_PASSPORT',vv)
expected={
 'alpha6':((3,3),(2,2,1,1),(3,3)),
 'alpha4':((3,1),(3,1),(2,2)),
 'beta6':((2,2,2),(2,2,2),(3,3)),
 'beta3':((2,1),(2,1),(3,)),
 'beta2':((2,),(2,),(1,1)),
}
# Branch-value order may permute for alpha models; compare multisets of profiles.
model_ok=all(sorted(passports[k])==sorted(expected[k]) for k in expected)

# Small symmetric-group uniqueness up to simultaneous conjugacy.
def cprof(p):
 seen=set();L=[]
 for i in range(len(p)):
  if i in seen: continue
  x=i;l=0
  while x not in seen:seen.add(x);l+=1;x=p[x]
  L.append(l)
 return tuple(sorted(L,reverse=True))
def transitive3(a,b,c):
 n=len(a);seen={0};q=deque([0]);gens=(a,b,c,inv(a),inv(b),inv(c))
 while q:
  x=q.popleft()
  for g in gens:
   y=g[x]
   if y not in seen:seen.add(y);q.append(y)
 return len(seen)==n
def conj(p,h): return comp(comp(h,p),inv(h))
def canon3(T,alls):
 return min(tuple(v for p in (conj(T[0],h),conj(T[1],h),conj(T[2],h)) for v in p) for h in alls)
def uniqueness(passp):
 n=sum(passp[0]); alls=[tuple(p) for p in permutations(range(n))]
 pools=[[p for p in alls if cprof(p)==prof] for prof in passp]
 keys=set();sol=0
 for a in pools[0]:
  for b in pools[1]:
   c=inv(comp(a,b))
   if cprof(c)!=passp[2] or not transitive3(a,b,c): continue
   sol+=1;keys.add(canon3((a,b,c),alls))
 return sol,len(keys)
unique={k:uniqueness(expected[k]) for k in expected}
for k,(sol,classes) in unique.items(): print(k.upper()+'_TRANSITIVE_TRIPLES',sol,'SIMULTANEOUS_CONJUGACY_CLASSES',classes)

# Compare actual restricted triples with the model passport classes up to branch-value permutation.
def triple_signature(T): return sorted([T[0],T[1],T[2]])
actual33=[r for r in R33 if r[0]>1]; actual22=[r for r in R22 if r[0]>1]
want33=[sorted(expected['alpha6'])]*2+[sorted(expected['alpha4'])]*2
want22=[sorted(expected['beta6'])]+[sorted(expected['beta3'])]*4+[sorted(expected['beta2'])]
actual33sig=sorted([sorted([r[1],r[2],r[3]]) for r in actual33])
actual22sig=sorted([sorted([r[1],r[2],r[3]]) for r in actual22])
match33=actual33sig==sorted(want33)
match22=actual22sig==sorted(want22)
ok=([r[0] for r in R33]==[6,6,4,4,1,1,1] and [r[0] for r in R22]==[6,3,3,3,3,2,1,1,1]
    and model_ok and all(c==1 for _,c in unique.values()) and match33 and match22)
print('ACTUAL_3A3A_MODELS_MATCH',match33)
print('ACTUAL_2A2A_MODELS_MATCH',match22)
print('ALL_COLLISION_TAIL_CHECKS_PASS',ok)
raise SystemExit(0 if ok else 1)
