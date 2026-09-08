#!/usr/bin/env python3
"""Exact five-state cyclic carrier attached to the certified pair fusion."""
from __future__ import annotations
import contextlib, importlib.util, io
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / 'base'
spec = importlib.util.spec_from_file_location('basecert', BASE / 'm23_common_component_certificate.py')
c = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(c)

P = tuple(c.P)

def delta(T):
    T = tuple(T)
    for i in range(4):
        T = c.hm(T, i, 1)
    return T

keys = []
T = P
for _ in range(6):
    keys.append(c.ckey(T))
    T = delta(T)

checks = {
    'FIRST_FIVE_DISTINCT': len(set(keys[:5])) == 5,
    'FIFTH_ITERATE_RETURNS': keys[5] == keys[0],
    'FUSED_PRODUCT_ONE': c.h.prod(P) == c.h.I(),
    'FUSED_ORDERS': [c.ctype_order(x) for x in P] == [2, 2, 3, 2, 2],
}
for k, v in checks.items():
    print(k, bool(v))
print('DELTA_ORBIT_LENGTH', 5 if checks['FIRST_FIVE_DISTINCT'] and checks['FIFTH_ITERATE_RETURNS'] else 0)
print('BETA5_PASSPORT', '1^5,5,5')
print('ALL_BETA5_FINITE_CARRIER_CHECKS_PASS', all(checks.values()))
raise SystemExit(0 if all(checks.values()) else 1)
