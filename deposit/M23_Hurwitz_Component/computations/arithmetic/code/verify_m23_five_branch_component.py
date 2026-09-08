#!/usr/bin/env python3
"""Reconstruct the valid five-branch M23 seed component.

The computation starts from an explicit product-one generating tuple of type
(2A,2A,2A,2A,3A), closes its colored braid orbit under q1,q2,q3,q4^2,
and records the four boundary cycle histograms.  It makes no arithmetic
multiplier or rational-point claim.
"""
from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path
from typing import Iterable, Sequence
import argparse
import json
import platform
import sys
import time

import sympy
from sympy.combinatorics import Permutation, PermutationGroup

try:
    import resource
except ImportError:  # pragma: no cover
    resource = None  # type: ignore[assignment]

N = 23
M23_ORDER = 10_200_960
ID = tuple(range(N))
Perm = tuple[int, ...]
Tuple5 = tuple[Perm, Perm, Perm, Perm, Perm]


def cyc(*cycles: Iterable[int]) -> Perm:
    p = list(range(N))
    for cycle in cycles:
        c = [x - 1 for x in cycle]
        for a, b in zip(c, c[1:] + c[:1], strict=True):
            p[a] = b
    return tuple(p)


def mul(p: Perm, q: Perm) -> Perm:
    """Composition p after q: (p*q)(i)=p(q(i))."""
    return tuple(p[q[i]] for i in range(N))


@lru_cache(maxsize=None)
def inv(p: Perm) -> Perm:
    out = [0] * N
    for i, j in enumerate(p):
        out[j] = i
    return tuple(out)


def product(tup: Sequence[Perm]) -> Perm:
    out = ID
    for g in tup:
        out = mul(out, g)
    return out


def order(p: Perm) -> int:
    result = 1
    for length in cycle_lengths(p):
        result = lcm(result, length)
    return result


def cycle_lengths(p: Sequence[int]) -> list[int]:
    seen = [False] * len(p)
    out: list[int] = []
    for i in range(len(p)):
        if seen[i]:
            continue
        x = i
        length = 0
        while not seen[x]:
            seen[x] = True
            length += 1
            x = p[x]
        out.append(length)
    return sorted(out)


def cycle_type(p: Sequence[int]) -> dict[int, int]:
    return dict(sorted(Counter(cycle_lengths(p)).items()))


def index(p: Sequence[int]) -> int:
    return sum(length - 1 for length in cycle_lengths(p))


def relabel_perm(g: Perm, lab: Sequence[int]) -> Perm:
    out = [0] * N
    for old in range(N):
        out[lab[old]] = lab[g[old]]
    return tuple(out)


@lru_cache(maxsize=None)
def canonical(tup: Tuple5) -> Tuple5:
    generators = list(tup) + [inv(g) for g in tup]
    best_key: tuple[int, ...] | None = None
    best_tuple: Tuple5 | None = None
    for root in range(N):
        labels: list[int | None] = [None] * N
        labels[root] = 0
        queue = deque([root])
        next_label = 1
        while queue:
            x = queue.popleft()
            for g in generators:
                y = g[x]
                if labels[y] is None:
                    labels[y] = next_label
                    next_label += 1
                    queue.append(y)
        if next_label != N:
            raise RuntimeError("canonicalization requires a transitive tuple")
        lab = [int(x) for x in labels]
        relabeled = tuple(relabel_perm(g, lab) for g in tup)
        key = tuple(value for g in relabeled for value in g)
        if best_key is None or key < best_key:
            best_key = key
            best_tuple = relabeled  # type: ignore[assignment]
    assert best_tuple is not None
    return best_tuple


def qmove(tup: Tuple5, i: int, sign: int = 1) -> Tuple5:
    out = list(tup)
    if sign == 1:
        a, b = out[i], out[i + 1]
        out[i] = mul(mul(a, b), inv(a))
        out[i + 1] = a
    elif sign == -1:
        a, b = out[i], out[i + 1]
        out[i] = b
        out[i + 1] = mul(mul(inv(b), a), b)
    else:
        raise ValueError("sign must be +/-1")
    return tuple(out)  # type: ignore[return-value]


def apply_generator(tup: Tuple5, name: str, inverse_flag: bool = False) -> Tuple5:
    sign = -1 if inverse_flag else 1
    if name == "q1":
        return qmove(tup, 0, sign)
    if name == "q2":
        return qmove(tup, 1, sign)
    if name == "q3":
        return qmove(tup, 2, sign)
    if name == "q4sq":
        return qmove(qmove(tup, 3, sign), 3, sign)
    raise ValueError(name)


@lru_cache(maxsize=None)
def action_can(tup: Tuple5, name: str, inverse_flag: bool = False) -> Tuple5:
    return canonical(apply_generator(tup, name, inverse_flag))


def build_orbit(seed: Tuple5, names: Sequence[str]) -> tuple[list[Tuple5], dict[Tuple5, int]]:
    start = canonical(seed)
    states = [start]
    position = {start: 0}
    queue = deque([start])
    while queue:
        state = queue.popleft()
        for name in names:
            # Forward closure is enough: on a finite set, the semigroup generated
            # by bijections contains the inverse powers as well.
            nxt = action_can(state, name, False)
            if nxt not in position:
                position[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
    return states, position


def build_operator(states: Sequence[Tuple5], position: dict[Tuple5, int], name: str) -> list[int]:
    return [position[action_can(state, name, False)] for state in states]


def permutation_histogram(p: Sequence[int]) -> dict[str, object]:
    lengths = cycle_lengths(p)
    return {
        "num_cycles": len(lengths),
        "fixed_points": sum(length == 1 for length in lengths),
        "max_cycle_length": max(lengths),
        "histogram": dict(sorted(Counter(lengths).items())),
        "index": sum(length - 1 for length in lengths),
    }


def fingerprint_tuple(tup: Tuple5) -> str:
    payload = bytes(value + 1 for g in tup for value in g)
    return sha256(payload).hexdigest()


SEED: Tuple5 = (
    cyc((2, 6), (4, 12), (5, 8), (10, 17), (13, 22), (14, 19), (15, 16), (20, 21)),
    cyc((1, 15), (2, 19), (3, 8), (4, 17), (6, 7), (11, 23), (12, 21), (14, 18)),
    cyc((2, 7), (4, 13), (6, 8), (9, 22), (10, 11), (14, 17), (15, 20), (16, 18)),
    cyc((1, 20), (3, 7), (4, 18), (6, 19), (9, 22), (10, 23), (13, 17), (14, 16)),
    cyc((3, 19, 14), (4, 21, 15), (5, 7, 8), (10, 17, 11), (12, 16, 20), (13, 22, 18)),
)

GENERATOR_NAMES = ("q1", "q2", "q3", "q4sq")
EXPECTED_HISTOGRAMS = {
    "q1": {1: 24, 2: 292, 3: 796, 4: 1248, 5: 980, 6: 1428},
    "q2": {1: 24, 2: 292, 3: 796, 4: 1248, 5: 980, 6: 1428},
    "q3": {1: 24, 2: 292, 3: 796, 4: 1248, 5: 980, 6: 1428},
    "q4sq": {1: 18, 2: 468, 3: 582, 4: 984, 5: 900, 7: 840, 8: 192, 11: 264},
}


def verify() -> dict[str, object]:
    started = time.perf_counter()

    seed_group = PermutationGroup([Permutation(list(g)) for g in SEED])
    generated_order = int(seed_group.order())

    # The fixed ATLAS copy used by the four-branch certificate.
    from verify_m23_hurwitz_certificate import ATLAS_A, ATLAS_B, to_sympy
    atlas_group = PermutationGroup([to_sympy(ATLAS_A), to_sympy(ATLAS_B)])
    membership = [bool(atlas_group.contains(Permutation(list(g)))) for g in SEED]

    orders = [order(g) for g in SEED]
    types = [cycle_type(g) for g in SEED]
    indices = [index(g) for g in SEED]
    source_genus = 1 - N + sum(indices) // 2
    adjacent_orders = [order(mul(SEED[i], SEED[i + 1])) for i in range(4)]

    states, position = build_orbit(SEED, GENERATOR_NAMES)
    operators = {name: build_operator(states, position, name) for name in GENERATOR_NAMES}
    histograms = {name: permutation_histogram(p) for name, p in operators.items()}

    checks = {
        "seed_membership_in_atlas_copy": all(membership),
        "seed_group_order": generated_order == M23_ORDER,
        "seed_product_one": product(SEED) == ID,
        "seed_orders": orders == [2, 2, 2, 2, 3],
        "seed_indices": indices == [8, 8, 8, 8, 12],
        "source_genus_zero": source_genus == 0,
        "adjacent_orders": adjacent_orders == [6, 6, 3, 11],
        "orbit_size": len(states) == 21_456,
        "operator_bijections": all(sorted(p) == list(range(len(states))) for p in operators.values()),
        "histograms": all(
            histograms[name]["histogram"] == EXPECTED_HISTOGRAMS[name]
            for name in GENERATOR_NAMES
        ),
    }

    elapsed = time.perf_counter() - started
    peak_rss_kb = None
    if resource is not None:
        peak_rss_kb = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    return {
        "title": "M23 five-branch colored braid seed component",
        "scope": "geometric finite computation only; no multiplier, Q-cusp, or rational-point conclusion",
        "conventions": {
            "permutation_composition": "p after q",
            "product_one": "g1*g2*g3*g4*g5=identity in this convention",
            "colored_braid_generators": list(GENERATOR_NAMES),
        },
        "seed": {
            "fingerprint_sha256": fingerprint_tuple(SEED),
            "membership_in_atlas_copy": membership,
            "generated_group_order": generated_order,
            "orders": orders,
            "cycle_types": types,
            "indices": indices,
            "product_one": product(SEED) == ID,
            "source_genus": source_genus,
            "adjacent_product_orders": adjacent_orders,
        },
        "component": {
            "size": len(states),
            "operator_histograms": histograms,
        },
        "checks": checks,
        "pass": all(checks.values()),
        "environment": {
            "python": platform.python_version(),
            "sympy": sympy.__version__,
            "platform": platform.platform(),
        },
        "performance": {"elapsed_seconds": elapsed, "peak_rss_kb": peak_rss_kb},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    result = verify()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print("M23_FIVE_BRANCH_COMPONENT_CERTIFICATE")
    print("PASS", result["pass"])
    print("GENERATED_GROUP_ORDER", result["seed"]["generated_group_order"])
    print("PRODUCT_ONE", result["seed"]["product_one"])
    print("ORDERS", result["seed"]["orders"])
    print("INDICES", result["seed"]["indices"])
    print("SOURCE_GENUS", result["seed"]["source_genus"])
    print("ADJACENT_PRODUCT_ORDERS", result["seed"]["adjacent_product_orders"])
    print("COMPONENT_SIZE", result["component"]["size"])
    for name in GENERATOR_NAMES:
        h = result["component"]["operator_histograms"][name]
        print(name.upper(), "CYCLES", h["num_cycles"], "FIXED", h["fixed_points"], "MAX", h["max_cycle_length"], "HIST", h["histogram"])
    failed = [name for name, ok in result["checks"].items() if not ok]
    print("FAILED_CHECKS", failed)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
