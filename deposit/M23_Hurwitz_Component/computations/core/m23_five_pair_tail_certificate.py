#!/usr/bin/env python3
"""Exact five-pair local-tail and genus-ledger certificate."""
from __future__ import annotations
from collections import deque
from pathlib import Path
import contextlib, importlib.util, io
import sympy as sp
HERE=Path(__file__).resolve().parent
BASE=HERE/'base'
def load(path,name):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
 with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(m)
 return m
c=load(BASE/'m23_common_component_certificate.py','common')
h=c.h
N=23

def sheet_orbits(a,b):
 gens=(a,b,h.inv(a),h.inv(b));unseen=set(range(N));out=[]
 while unseen:
  root=min(unseen);seen={root};q=deque([root])
  while q:
   x=q.popleft()
   for g in gens:
    y=g[x]
    if y not in seen:seen.add(y);q.append(y)
  out.append(len(seen));unseen-=seen
 return sorted(out,reverse=True)
profiles=[sheet_orbits(c.U[2*i],c.U[2*i+1]) for i in range(5)]
orders=[c.ctype_order(p) for p in c.P]
print('PAIR_PRODUCT_ORDERS',orders)
for i,p in enumerate(profiles,1): print(f'PAIR{i}_SHEET_ORBITS',p)
size84=[4,4,2,2,2,2,2,2,1,1,1]
size14=[4,4,4,4,1,1,1,1,1,1,1]
size3=[6,3,3,3,3,2,1,1,1]
profile_ok=(profiles[0]==size84 and profiles[1]==size14 and profiles[2]==size3 and profiles[3]==size14 and profiles[4]==size84)

x=sp.symbols('x')
nu4_num=(x**2+1)**2; nu4_den=4*x**2
nu2_num=x**2;nu2_den=x**2+1

def root_mults(poly,d):
 p=sp.Poly(sp.expand(poly),x);out=[]
 for fac,e in sp.factor_list(p.as_expr())[1]:out += [int(e)]*sp.Poly(fac,x).degree()
 if p.degree()<d:out.append(d-p.degree())
 return tuple(sorted(out,reverse=True))
def passport(num,den):
 d=max(sp.degree(num,x),sp.degree(den,x));return (root_mults(num,d),root_mults(num-den,d),root_mults(den,d))
p4=passport(nu4_num,nu4_den);p2=passport(nu2_num,nu2_den)
print('NU4_PASSPORT',p4);print('NU2_PASSPORT',p2)
formula_ok=sorted(p4)==sorted(((2,2),(2,2),(2,2))) and sorted(p2)==sorted(((2,),(2,),(1,1)))
components=[len(p) for p in profiles]
branches=[15,15,11,15,15]
contrib=[b-k for b,k in zip(branches,components)]
print('TAIL_COMPONENT_COUNTS',components)
print('NODE_BRANCH_COUNTS',branches)
print('GENUS_CONTRIBUTIONS',contrib)
print('GENUS_SUM',sum(contrib))
ok=(orders==[2,2,3,2,2] and profile_ok and formula_ok and contrib==[4,4,2,4,4] and sum(contrib)==18)
print('ALL_FIVE_PAIR_TAIL_CHECKS_PASS',ok)
raise SystemExit(0 if ok else 1)
