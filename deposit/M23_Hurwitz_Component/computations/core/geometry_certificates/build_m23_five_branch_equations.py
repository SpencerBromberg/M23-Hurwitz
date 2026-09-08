#!/usr/bin/env python3
"""Expand the 72 exact coefficient equations for the five-branch reconstruction target."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sympy as sp
OUT=Path(__file__).resolve().parent
x,lam,mu=sp.symbols('x lambda mu')

def coeffs(prefix,deg): return list(sp.symbols(' '.join(f'{prefix}{i}' for i in range(deg+1))))
def conv(a,b):
 out=[sp.Integer(0)]*(len(a)+len(b)-1)
 for i,ai in enumerate(a):
  for j,bj in enumerate(b): out[i+j]+=ai*bj
 return out
def power(a,n):
 out=[sp.Integer(1)]
 for _ in range(n): out=conv(out,a)
 return out
def add(a,b,scale=1):
 n=max(len(a),len(b));return [(a[i] if i<len(a) else 0)+scale*(b[i] if i<len(b) else 0) for i in range(n)]
def polystr(a):
 return ' + '.join(f'({sp.sstr(c)})*x^{i}' for i,c in enumerate(a) if c!=0) or '0'
A0=coeffs('a0_',8); L0=coeffs('l0_',7)
Ai=coeffs('ai_',8); Li=coeffs('li_',7)
A1=coeffs('a1_',8); L1=coeffs('l1_',7)
Al=coeffs('al_',8); Ll=coeffs('ll_',7)
C=coeffs('c_',6); M=coeffs('m_',5)
P=conv(power(A0,2),L0); Q=conv(power(Ai,2),Li)
F1=conv(power(A1,2),L1); Fl=conv(power(Al,2),Ll); Fc=conv(power(C,3),M)
assert all(len(z)==24 for z in [P,Q,F1,Fl,Fc])
E1=add(add(P,Q,-1),F1,-1)
E2=add(add(P,Q,-lam),Fl,-1)
E3=add(add(P,Q,-mu),Fc,-1)
eqs=E1+E2+E3
variables=A0+L0+Ai+Li+A1+L1+Al+Ll+C+M+[lam,mu]
assert len(eqs)==72 and len(variables)==83
payload={'coefficient_variable_count_before_gauge':83,'equation_count':72,
 'equation_blocks':{'P-Q-A1^2L1':24,'P-lambda Q-Alambda^2Llambda':24,'P-mu Q-C^3M':24},
 'degrees':{'A0':8,'L0':7,'Ainf':8,'Linf':7,'A1':8,'L1':7,'Alambda':8,'Llambda':7,'C':6,'M':5},
 'variables':[str(v) for v in variables],
 'equations':[sp.sstr(e) for e in eqs],
 'scope':'Exact five-branch ramification reconstruction target; Nielsen data select the M23 component.'}
j=OUT/'M23_five_branch_equations.json'; t=OUT/'M23_five_branch_equations.txt'
j.write_text(json.dumps(payload,indent=2)+'\n')
with t.open('w') as f:
 f.write('M23 FIVE-BRANCH COEFFICIENT MODEL\n'+'='*60+'\n\n')
 for name,a in [('A0',A0),('L0',L0),('Ainf',Ai),('Linf',Li),('A1',A1),('L1',L1),('Alambda',Al),('Llambda',Ll),('C',C),('M',M)]: f.write(f'{name}(x) = {polystr(a)}\n')
 f.write('\nP=A0^2 L0; Q=Ainf^2 Linf\n')
 for base,block,label in [(0,E1,'P-Q-A1^2 L1'),(24,E2,'P-lambda Q-Alambda^2 Llambda'),(48,E3,'P-mu Q-C^3 M')]:
  f.write(f'\n{label}\n')
  for k,e in enumerate(block,base): f.write(f'E{k:02d} = {sp.sstr(e)}\n')
for path in [j,t]: print(path.name,hashlib.sha256(path.read_bytes()).hexdigest())
print('PASS equation_count=72 variable_count_before_gauge=83 degrees=23')
