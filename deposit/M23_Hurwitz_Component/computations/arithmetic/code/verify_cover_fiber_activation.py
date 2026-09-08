#!/usr/bin/env python3
"""Exact finite-etale fiber polynomial and rational-root activation verifier.

Input is a rational multiplication matrix M_z for a candidate primitive element z
and the coordinate vector of 1 in the chosen Q-basis. The script:
  1. verifies that 1,z,...,z^(n-1) span (Krylov determinant nonzero),
  2. computes F_a(T)=det(TI-M_z),
  3. verifies F_a is square-free,
  4. finds every rational root r,
  5. constructs e_r=G_r(z)/G_r(r), F_a=(T-r)G_r,
  6. verifies e_r^2=e_r, rank(m_e)=Tr(m_e)=1,
  7. sets U_r=1-e_r and verifies U_r^27=U_r and 1-U_r^26=e_r.

All arithmetic is exact over Q.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

T = sp.Symbol("T")


def q(x):
    return sp.Rational(str(x))


def load_input(path: Path):
    data = json.loads(path.read_text())
    M = sp.Matrix([[q(v) for v in row] for row in data["multiplication_matrix"]])
    one = sp.Matrix([q(v) for v in data["unit_vector"]])
    if M.rows != M.cols:
        raise ValueError("multiplication matrix must be square")
    if one.rows != M.rows or one.cols != 1:
        raise ValueError("unit vector has wrong dimension")
    return data, M, one


def mat_pow_apply(M, v, n):
    out = v
    for _ in range(n):
        out = M * out
    return out


def eval_poly_matrix(poly: sp.Poly, M: sp.Matrix):
    n = M.rows
    out = sp.zeros(n)
    power = sp.eye(n)
    for i in range(poly.degree() + 1):
        out += poly.nth(i) * power
        power = power * M
    return sp.simplify(out)


def matrix_to_json(M):
    return [[str(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]


def vector_to_json(v):
    return [str(v[i, 0]) for i in range(v.rows)]


def rational_roots(P: sp.Poly):
    roots = []
    for f, mult in sp.factor_list(P.as_expr())[1]:
        fp = sp.Poly(f, T, domain=sp.QQ)
        if fp.degree() == 1:
            roots.extend([-fp.nth(0) / fp.nth(1)] * mult)
    return sorted(set(roots), key=lambda x: (float(x), str(x)))


def sha256_text(text: str):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    input_path = Path(args.input)
    data, M, one = load_input(input_path)
    n = M.rows

    krylov_cols = []
    v = one
    for _ in range(n):
        krylov_cols.append(v)
        v = M * v
    K = sp.Matrix.hstack(*krylov_cols)
    krylov_rank = K.rank()
    primitive = krylov_rank == n

    cp = M.charpoly(T).as_poly(domain=sp.QQ)
    square_free = sp.gcd(cp, cp.diff()).degree() == 0
    roots = rational_roots(cp)

    certs = []
    I = sp.eye(n)
    for r in roots:
        G = sp.div(cp, sp.Poly(T-r, T, domain=sp.QQ), domain=sp.QQ)[0]
        den = sp.Rational(G.eval(r))
        if den == 0:
            raise AssertionError("square-free rational factor has zero derivative")
        eM = eval_poly_matrix(sp.Poly(G.as_expr()/den, T, domain=sp.QQ), M)
        UM = I - eM
        e_idem = eM*eM == eM
        rank = eM.rank()
        tr = sp.trace(eM)
        U27 = UM**27 == UM
        one_minus_U26 = I - UM**26 == eM
        eval_vector = eM * one
        certs.append({
            "root": str(r),
            "G_polynomial": str(G.as_expr()),
            "G_at_root": str(den),
            "idempotent_matrix": matrix_to_json(eM),
            "idempotent_on_unit": vector_to_json(eval_vector),
            "idempotent_squared_equals_idempotent": bool(e_idem),
            "multiplication_rank": int(rank),
            "multiplication_trace": str(tr),
            "rank_one_trace_one": bool(rank == 1 and tr == 1),
            "U27_equals_U": bool(U27),
            "one_minus_U26_equals_idempotent": bool(one_minus_U26),
            "checks_pass": bool(e_idem and rank == 1 and tr == 1 and U27 and one_minus_U26),
        })

    input_text = input_path.read_text()
    out = {
        "input_label": data.get("label", ""),
        "input_sha256": sha256_text(input_text),
        "dimension": n,
        "krylov_rank": int(krylov_rank),
        "primitive_element_verified": bool(primitive),
        "fiber_polynomial": str(cp.as_expr()),
        "fiber_polynomial_degree": int(cp.degree()),
        "fiber_polynomial_square_free": bool(square_free),
        "fiber_polynomial_factorization": str(sp.factor(cp.as_expr())),
        "rational_roots": [str(r) for r in roots],
        "rational_root_count": len(roots),
        "root_certificates": certs,
        "clause_ii_cover_to_Fa_pass": bool(primitive and cp.degree() == n and square_free),
        "clause_iii_rational_root_activation_pass": bool(certs and all(c["checks_pass"] for c in certs)),
    }
    out["overall_pass"] = bool(out["clause_ii_cover_to_Fa_pass"] and out["clause_iii_rational_root_activation_pass"])
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out, indent=2, sort_keys=True))
    raise SystemExit(0 if out["overall_pass"] else 2)


if __name__ == "__main__":
    main()
