#!/usr/bin/env python3
"""Exact rank-one idempotent extraction from a primitive-element algebra.

Input: a square-free polynomial P(Z) in Q[Z] defining A = Q[Z]/(P).
For each rational root r, this constructs
    e_r = (P(Z)/(Z-r)) / (P'(r))  mod P,
which is the primitive central idempotent of the Q-factor Z=r. It then sets
U_r = 1-e_r and verifies U_r^27=U_r and Tr(1-U_r^26)=1.

This is the decisive exact computation once a genuine interior primitive
polynomial P_s(Z) has been exported from the saturated Hurwitz fiber algebra.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import sympy as sp

Z = sp.Symbol('Z')

def coeff_list(poly: sp.Poly):
    d = poly.degree()
    return [str(sp.Rational(poly.nth(i))) for i in range(d + 1)]

def mod_poly(expr, P):
    return sp.Poly(expr, Z, domain=sp.QQ).rem(P).as_expr()

def pow_mod(expr, n, P):
    result = sp.Integer(1)
    base = mod_poly(expr, P)
    while n:
        if n & 1:
            result = mod_poly(result*base, P)
        n >>= 1
        if n:
            base = mod_poly(base*base, P)
    return result

def mul_matrix(expr, P):
    n = P.degree()
    cols=[]
    for j in range(n):
        y=sp.Poly(mod_poly(expr*Z**j,P),Z,domain=sp.QQ)
        cols.append([sp.Rational(y.nth(i)) for i in range(n)])
    return sp.Matrix(n,n,lambda i,j: cols[j][i])

def parse_polynomial(args):
    if args.expr:
        return sp.Poly(sp.sympify(args.expr, locals={'Z':Z}),Z,domain=sp.QQ)
    data=json.loads(Path(args.json).read_text())
    if 'expr' in data:
        return sp.Poly(sp.sympify(data['expr'],locals={'Z':Z}),Z,domain=sp.QQ)
    coeffs=[sp.Rational(c) for c in data['coefficients_low_to_high']]
    return sp.Poly(sum(c*Z**i for i,c in enumerate(coeffs)),Z,domain=sp.QQ)

def main():
    ap=argparse.ArgumentParser()
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--expr',help="Polynomial expression, e.g. 'Z**27-Z'")
    g.add_argument('--json',help='JSON with expr or coefficients_low_to_high')
    ap.add_argument('--output',default='idempotent_certificate.json')
    args=ap.parse_args()
    P=parse_polynomial(args)
    if P.degree() < 1:
        raise SystemExit('polynomial must have positive degree')
    if sp.gcd(P,P.diff()).degree()!=0:
        raise SystemExit('polynomial is not square-free; algebra is not etale')
    fac=sp.factor_list(P.as_expr(), modulus=None)
    rational=[]
    for f,m in fac[1]:
        fp=sp.Poly(f,Z,domain=sp.QQ)
        if fp.degree()==1:
            a=fp.nth(1); b=fp.nth(0)
            rational.append(-b/a)
    certs=[]
    for r in sorted(set(rational), key=lambda q: (float(q),str(q))):
        Q=sp.div(P,sp.Poly(Z-r,Z,domain=sp.QQ),domain=sp.QQ)[0]
        den=sp.Rational(Q.eval(r))
        if den==0: raise AssertionError('square-free factor gave zero denominator')
        e=mod_poly(Q.as_expr()/den,P)
        U=mod_poly(1-e,P)
        e2=mod_poly(e*e,P)
        U27=pow_mod(U,27,P)
        one_minus=mod_poly(1-pow_mod(U,26,P),P)
        Me=mul_matrix(e,P)
        cert={
          'rational_root':str(r),
          'idempotent_e_expr':str(sp.expand(e)),
          'U_expr':str(sp.expand(U)),
          'e_coefficients_low_to_high':coeff_list(sp.Poly(e,Z,domain=sp.QQ)),
          'U_coefficients_low_to_high':coeff_list(sp.Poly(U,Z,domain=sp.QQ)),
          'e_squared_equals_e':sp.expand(e2-e)==0,
          'U27_equals_U':sp.expand(U27-U)==0,
          'one_minus_U26_equals_e':sp.expand(one_minus-e)==0,
          'multiplication_trace_e':str(Me.trace()),
          'multiplication_rank_e':int(Me.rank()),
          'rank_one_trace_one':Me.trace()==1 and Me.rank()==1,
        }
        cert['checks_pass']=all([cert['e_squared_equals_e'],cert['U27_equals_U'],cert['one_minus_U26_equals_e'],cert['rank_one_trace_one']])
        certs.append(cert)
    out={
      'polynomial':str(P.as_expr()),
      'degree':P.degree(),
      'square_free':True,
      'factorization':str(sp.factor(P.as_expr())),
      'rational_linear_factors':len(certs),
      'certificates':certs,
      'overall_pass':bool(certs) and all(c['checks_pass'] for c in certs),
    }
    Path(args.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,indent=2,sort_keys=True))
    if not out['overall_pass']:
        raise SystemExit(2)

if __name__=='__main__': main()
