#!/usr/bin/env python3
"""Exact six-branch M23 constructions.

Produces two verified product-one generating tuples:
(1) a two-point stabilization of the four-branch (2A,2A,3A,5A) tuple;
(2) a pure Harbater--Mumford three-pair tuple from the standard ATLAS generators.
"""
from __future__ import annotations
import json
from math import gcd
from pathlib import Path
from sympy.combinatorics import Permutation, PermutationGroup

N=23
ID=Permutation(list(range(N)))
def perm_from_cycles(cycles):
    arr=list(range(N))
    for cyc in cycles:
        c=[x-1 for x in cyc]
        for i,x in enumerate(c): arr[x]=c[(i+1)%len(c)]
    return Permutation(arr)
def product(seq):
    p=ID
    for g in seq: p=p*g
    return p
def cycles_with_fixed(p):
    return p.cyclic_form+[[i] for i in range(N) if p(i)==i]
def index(p): return N-len(cycles_with_fixed(p))
def cycle_type(p):
    d={}
    for c in cycles_with_fixed(p): d[len(c)]=d.get(len(c),0)+1
    return {str(k):v for k,v in sorted(d.items())}
def lcm(a,b): return a*b//gcd(a,b)
def audit(seq,labels,name):
    inds=[index(x) for x in seq]
    two_g=2-2*N+sum(inds)
    assert two_g%2==0
    adj=[int(product([seq[i],seq[i+1]]).order()) for i in range(len(seq)-1)]
    branch_mod=1
    for x in seq: branch_mod=lcm(branch_mod,int(x.order()))
    return {
      'name':name,
      'class_vector':labels,
      'orders':[int(x.order()) for x in seq],
      'cycle_types':[cycle_type(x) for x in seq],
      'indices':inds,
      'index_sum':sum(inds),
      'source_genus':two_g//2,
      'product_one':bool(product(seq).is_Identity),
      'generated_group_order':int(PermutationGroup(seq).order()),
      'generated_group_is_M23':int(PermutationGroup(seq).order())==10200960,
      'adjacent_product_orders':adj,
      'branch_class_power_modulus':branch_mod,
      'width_one_adjacent_pairs':[i+1 for i,o in enumerate(adj) if o==1],
    }

# Printed four-branch tuple.
g1=perm_from_cycles([(1,7),(2,20),(4,11),(5,10),(6,18),(9,21),(16,23),(17,22)])
g2=perm_from_cycles([(1,20),(2,7),(5,10),(6,9),(8,14),(13,19),(17,22),(18,21)])
g3=perm_from_cycles([(1,21,2),(3,23,14),(4,16,18),(5,8,9),(6,12,22),(13,15,20)])
g4=perm_from_cycles([(2,6,22,12,21),(3,8,5,18,23),(4,9,14,16,11),(7,20,15,19,13)])
stabilized=[g1,g2,g3,g4,g1,g1]

# Standard ATLAS generators, with ord(a)=2, ord(b)=4, ord(ab)=23.
a=perm_from_cycles([(1,2),(3,4),(7,8),(9,10),(13,14),(15,16),(19,20),(21,22)])
b=perm_from_cycles([(1,16,11,3),(2,9,21,12),(4,5,8,23),(6,22,14,18),(13,20),(15,17)])
c=a*b
hm=[a,a**-1,b,b**-1,c,c**-1]

out={
 'stabilized_r6':audit(stabilized,['2A','2A','3A','5A','2A','2A'],'four-branch tuple plus inverse 2A pair'),
 'pure_HM_r6':audit(hm,['2A','2A','4A','4A','23A','23B'],'three inverse pairs (a,a^-1,b,b^-1,ab,(ab)^-1)'),
}
out['stabilized_r6']['checks_pass']=all([
 out['stabilized_r6']['product_one'],out['stabilized_r6']['generated_group_is_M23'],
 out['stabilized_r6']['orders']==[2,2,3,5,2,2],out['stabilized_r6']['source_genus']==8,
 out['stabilized_r6']['width_one_adjacent_pairs']==[5]])
out['pure_HM_r6']['checks_pass']=all([
 out['pure_HM_r6']['product_one'],out['pure_HM_r6']['generated_group_is_M23'],
 out['pure_HM_r6']['orders']==[2,2,4,4,23,23],out['pure_HM_r6']['source_genus']==22,
 out['pure_HM_r6']['width_one_adjacent_pairs']==[1,3,5]])
out['overall_pass']=out['stabilized_r6']['checks_pass'] and out['pure_HM_r6']['checks_pass']
Path(__file__).with_name('r6_exact_certificate.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
if not out['overall_pass']: raise SystemExit(1)
