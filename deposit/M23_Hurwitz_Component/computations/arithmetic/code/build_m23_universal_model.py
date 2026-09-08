#!/usr/bin/env python3
"""Build exact coefficient equations for the ordered M23 Hurwitz curve."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parents[1] / 'data'
OUT.mkdir(parents=True, exist_ok=True)
x, lam = sp.symbols('x lambda')

def conv(a, b):
    out = [sp.Integer(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i+j] += ai * bj
    return out

def power(a, n):
    out = [sp.Integer(1)]
    for _ in range(n):
        out = conv(out, a)
    return out

def add(a, b, scale_b=1):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + scale_b*(b[i] if i < len(b) else 0) for i in range(n)]

def deriv(a):
    return [sp.Integer(i+1)*a[i+1] for i in range(len(a)-1)]

def poly_string(coeffs):
    terms=[]
    for i,c in enumerate(coeffs):
        if c == 0: continue
        terms.append(f'({sp.sstr(c)})*x^{i}')
    return ' + '.join(terms) if terms else '0'

# 48 variables in the rational-coefficient gauge.
a4,a3,a2,a1,a0=sp.symbols('a4 a3 a2 a1 a0')
aa4,aa3,aa2,aa1,aa0=sp.symbols('aa4 aa3 aa2 aa1 aa0')
q2,q1,q0=sp.symbols('q2 q1 q0')
r3,r2,r1,r0=sp.symbols('r3 r2 r1 r0')
b6,b5,b4,b3,b2,b1,b0=sp.symbols('b6 b5 b4 b3 b2 b1 b0')
s7,s6,s5,s4,s3,s2,s1,s0=sp.symbols('s7 s6 s5 s4 s3 s2 s1 s0')
c7,c6,c5,c4,c3,c2,c1,c0=sp.symbols('c7 c6 c5 c4 c3 c2 c1 c0')
t7,t6,t5,t4,t3,t2,t1,t0=sp.symbols('t7 t6 t5 t4 t3 t2 t1 t0')
variables=[a4,a3,a2,a1,a0,aa4,aa3,aa2,aa1,aa0,q2,q1,q0,r3,r2,r1,r0,
 b6,b5,b4,b3,b2,b1,b0,s7,s6,s5,s4,s3,s2,s1,s0,
 c7,c6,c5,c4,c3,c2,c1,c0,t7,t6,t5,t4,t3,t2,t1,t0]
assert len(variables)==48

# Coefficients are listed from x^0 upward.
A3=[a0,a1,a2,a3,a4,0,1]
A1=[aa0,aa1,aa2,aa3,aa4,1]
Q5=[q0,q1,q2,0,1]
Q1=[r0,r1,r2,r3]
B2=[b0,b1,b2,b3,b4,b5,b6,0,1]
B1=[s0,s1,s2,s3,s4,s5,s6,s7]
C2=[c0,c1,c2,c3,c4,c5,c6,c7,1]
C1=[t0,t1,t2,t3,t4,t5,t6,t7]

p=conv(power(A3,3),A1)
q=conv(power(Q5,5),Q1)
Bfiber=conv(power(B2,2),B1)
Cfiber=conv(power(C2,2),C1)
assert len(p)==24 and len(q)==24 and len(Bfiber)==24 and len(Cfiber)==24
E1=add(add(p,q,-1),Bfiber,-1)
E2=add(add(p,q,-lam),Cfiber,-1)
equations=E1+E2
assert len(equations)==48

# Compact exact Wronskian coefficient vector, without full multivariate expansion.
W=add(conv(deriv(p),q),conv(p,deriv(q)),-1)
while W and W[-1] == 0:
    W.pop()
forced=conv(conv(conv(power(A3,2),B2),C2),power(Q5,4))
assert len(W)==45 and len(forced)==45

forms={'A3':A3,'A1':A1,'Q5':Q5,'Q1':Q1,'B2':B2,'B1':B1,'C2':C2,'C1':C1,'p':p,'q':q}
ledger={
 'base_coordinate':'lambda',
 'coefficient_variable_count':48,
 'equation_count':48,
 'equation_blocks':{'p-q-B2^2*B1':24,'p-lambda*q-C2^2*C1':24},
 'degrees':{k:len(v)-1 for k,v in forms.items()},
 'variables':[str(v) for v in variables],
 'map':'beta(lambda,u)=lambda',
 'plane_elimination_target':'H(lambda,u)=0 on the normalized M23 component',
 'wronskian':'p_x*q-p*q_x=c*A3^2*B2*C2*Q5^4, c != 0',
 'scope':'Exact universal coefficient model. Expanded H requires saturation, M23-component selection, quotient, and elimination.'
}

payload={**ledger,
 'forms_coefficients_low_to_high':{k:[sp.sstr(c) for c in v] for k,v in forms.items()},
 'equations':[sp.sstr(e) for e in equations],
 'wronskian_coefficients':[sp.sstr(e) for e in W],
 'forced_wronskian_factor_coefficients':[sp.sstr(e) for e in forced]}
json_path=OUT/'M23_universal_equations.json'
json_path.write_text(json.dumps(payload,indent=2),encoding='utf-8')
text_path=OUT/'M23_universal_equations.txt'
with text_path.open('w',encoding='utf-8') as f:
 f.write('M23 UNIVERSAL COEFFICIENT MODEL\n'+'='*72+'\n\n')
 for name in ['A3','A1','Q5','Q1','B2','B1','C2','C1']:
  f.write(f'{name}(x) = {poly_string(forms[name])}\n')
 f.write('\np(x)=A3(x)^3 A1(x)\nq(x)=Q5(x)^5 Q1(x)\n')
 f.write('\nE_0,...,E_23 = coefficients of p-q-B2^2 B1\n')
 for i,e in enumerate(E1): f.write(f'E{i:02d} = {sp.sstr(e)}\n')
 f.write('\nE_24,...,E_47 = coefficients of p-lambda*q-C2^2 C1\n')
 for i,e in enumerate(E2,24): f.write(f'E{i:02d} = {sp.sstr(e)}\n')
 f.write('\nClosed map on the normalized M23 component: beta(lambda,u)=lambda.\n')
summary=OUT/'BUILD_SUMMARY.json'
summary.write_text(json.dumps(ledger,indent=2),encoding='utf-8')
for path in [json_path,text_path,summary]:
 print(path.name,hashlib.sha256(path.read_bytes()).hexdigest())
print('PASS equation_count=48 variable_count=48 degrees=23,23')
