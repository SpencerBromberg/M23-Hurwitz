#!/usr/bin/env python3
"""Exact, deterministic certificate for an ordered M23 Hurwitz orbit.

The program starts from two explicit permutations generating the standard
23-point ATLAS copy of M23 and from an explicit product-one tuple of class type
(2A,2A,3A,5A).  It then performs the following exact checks:

* group order, standard-generator orders, transitivity and 2-transitivity;
* tuple membership, product one, generation, cycle types and source genus;
* the Hurwitz-move inverse and Artin relations on every enumerated state;
* canonicalization modulo simultaneous relabeling of the 23 sheets;
* breadth-first closure under three pure point-push words and their inverses;
* construction of the three peripheral permutations and the equal-slot half twist;
* operator inverses, product-one relation, connectedness, cycle partitions;
* Riemann--Hurwitz genus of the compactified ordered Hurwitz cover.

No floating-point arithmetic and no stored degree-980 permutation arrays are
used.  The only non-standard dependency is SymPy.

Usage:
    python verify_m23_ordered_hurwitz.py --json result.json
    python verify_m23_ordered_hurwitz.py --json result.json \
        --dump-data orbit_certificate.json.gz
"""

from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence
import argparse
import gzip
import hashlib
import json
import platform
import struct
import sys
import time

import sympy
from sympy.combinatorics import Permutation, PermutationGroup

try:  # Available on Unix; the verifier remains portable without it.
    import resource
except ImportError:  # pragma: no cover - Windows fallback
    resource = None  # type: ignore[assignment]

N = 23
M23_ORDER = 10_200_960
M22_ORDER = 443_520
ID = tuple(range(N + 1))

Perm = tuple[int, ...]
Tuple4 = tuple[Perm, Perm, Perm, Perm]
Word = tuple[tuple[int, int], ...]


def perm_from_cycles(cycles: Iterable[Iterable[int]]) -> Perm:
    p = list(range(N + 1))
    for cycle in cycles:
        c = list(cycle)
        if len(c) < 2:
            continue
        for i, x in enumerate(c):
            p[x] = c[(i + 1) % len(c)]
    return tuple(p)


def compose(a: Perm, b: Perm) -> Perm:
    """Right-action product: i^(a b) = (i^a)^b."""
    return tuple([0] + [b[a[i]] for i in range(1, N + 1)])


@lru_cache(maxsize=None)
def inverse(a: Perm) -> Perm:
    out = [0] * (N + 1)
    for i in range(1, N + 1):
        out[a[i]] = i
    return tuple(out)


def power(a: Perm, exponent: int) -> Perm:
    if exponent < 0:
        return power(inverse(a), -exponent)
    result = ID
    base = a
    e = exponent
    while e:
        if e & 1:
            result = compose(result, base)
        base = compose(base, base)
        e >>= 1
    return result


def conjugate_standard(g: Perm, h: Perm) -> Perm:
    """Return h^-1 g h."""
    return compose(compose(inverse(h), g), h)


def conjugate_by_left(h: Perm, g: Perm) -> Perm:
    """Return h g h^-1."""
    return compose(compose(h, g), inverse(h))


def tuple_product(t: Tuple4) -> Perm:
    out = ID
    for g in t:
        out = compose(out, g)
    return out


def tuple_conjugate(t: Tuple4, h: Perm) -> Tuple4:
    return tuple(conjugate_standard(g, h) for g in t)  # type: ignore[return-value]


def cycle_lengths(p: Sequence[int], start: int = 1) -> list[int]:
    used = [False] * len(p)
    lengths: list[int] = []
    for i in range(start, len(p)):
        if used[i]:
            continue
        x = i
        length = 0
        while not used[x]:
            used[x] = True
            length += 1
            x = p[x]
        lengths.append(length)
    return sorted(lengths)


def cycle_type(lengths: Sequence[int]) -> dict[int, int]:
    return dict(sorted(Counter(lengths).items()))


def permutation_order(p: Perm) -> int:
    order = 1
    for length in cycle_lengths(p):
        order = order * length // gcd(order, length)
    return order


def permutation_index(p: Sequence[int], start: int = 1) -> int:
    return sum(length - 1 for length in cycle_lengths(p, start=start))


def to_sympy(p: Perm) -> Permutation:
    return Permutation([p[i] - 1 for i in range(1, N + 1)])


# Standard ATLAS degree-23 generators, written on {1,...,23}.
# Provenance: ATLAS of Finite Groups / ATLAS online M23 natural
# permutation representation.  The certificate independently checks
# ord(a)=2, ord(b)=4, ord(ab)=23 and |<a,b>|=10,200,960.
ATLAS_A = perm_from_cycles(
    [(1, 2), (3, 4), (7, 8), (9, 10),
     (13, 14), (15, 16), (19, 20), (21, 22)]
)
ATLAS_B = perm_from_cycles(
    [(1, 16, 11, 3), (2, 9, 21, 12), (4, 5, 8, 23),
     (6, 22, 14, 18), (13, 20), (15, 17)]
)

# Explicit product-one tuple of class type (2A,2A,3A,5A).
TUPLE: Tuple4 = (
    perm_from_cycles(
        [(1, 7), (2, 20), (4, 11), (5, 10), (6, 18),
         (9, 21), (16, 23), (17, 22)]
    ),
    perm_from_cycles(
        [(1, 20), (2, 7), (5, 10), (6, 9), (8, 14),
         (13, 19), (17, 22), (18, 21)]
    ),
    perm_from_cycles(
        [(1, 21, 2), (3, 23, 14), (4, 16, 18),
         (5, 8, 9), (6, 12, 22), (13, 15, 20)]
    ),
    perm_from_cycles(
        [(2, 6, 22, 12, 21), (3, 8, 5, 18, 23),
         (4, 9, 14, 16, 11), (7, 20, 15, 19, 13)]
    ),
)

EXPECTED_TUPLE_TYPES = [
    {1: 7, 2: 8},
    {1: 7, 2: 8},
    {1: 5, 3: 6},
    {1: 3, 5: 4},
]
EXPECTED_OPERATOR_TYPES = [
    {1: 24, 2: 136, 3: 178, 5: 30},
    {1: 2, 2: 21, 3: 16, 4: 41, 5: 46, 7: 38, 8: 12, 11: 12},
    {1: 2, 2: 21, 3: 16, 4: 41, 5: 46, 7: 38, 8: 12, 11: 12},
]


def qmove(t: Tuple4, i: int, sign: int) -> Tuple4:
    """Hurwitz move q_i or q_i^-1, with zero-based i in {0,1,2}."""
    values = list(t)
    a, b = values[i], values[i + 1]
    if sign == 1:
        values[i] = compose(compose(a, b), inverse(a))
        values[i + 1] = a
    elif sign == -1:
        values[i] = b
        values[i + 1] = compose(compose(inverse(b), a), b)
    else:
        raise ValueError("sign must be +1 or -1")
    return tuple(values)  # type: ignore[return-value]


def apply_word(t: Tuple4, word: Word) -> Tuple4:
    out = t
    for i, exponent in word:
        sign = 1 if exponent > 0 else -1
        for _ in range(abs(exponent)):
            out = qmove(out, i, sign)
    return out


Q1: Word = ((0, 1),)
Q2: Word = ((1, 1),)
Q3: Word = ((2, 1),)
Q1_INV: Word = ((0, -1),)
Q2_INV: Word = ((1, -1),)
Q3_INV: Word = ((2, -1),)

B12: Word = ((0, 2),)
B12_INV: Word = ((0, -2),)
B13: Word = ((1, 1), (0, 2), (1, -1))
B13_INV: Word = ((1, 1), (0, -2), (1, -1))
B14: Word = ((2, 1), (1, 1), (0, 2), (1, -1), (2, -1))
B14_INV: Word = ((2, 1), (1, 1), (0, -2), (1, -1), (2, -1))

FORWARD_WORDS = (B12, B13, B14)
INVERSE_WORDS = (B12_INV, B13_INV, B14_INV)
ALL_WORDS = (B12, B12_INV, B13, B13_INV, B14, B14_INV)
WORD_NAMES = ("B12", "B13", "B14")


@lru_cache(maxsize=None)
def canonical_tuple_and_key(t: Tuple4) -> tuple[Tuple4, tuple[int, ...]]:
    """Canonical form under simultaneous conjugacy in S_23.

    For every possible root sheet, the connected colored Schreier graph is
    labeled by deterministic breadth-first traversal with edge colors
    g1,...,g4,g1^-1,...,g4^-1.  The lexicographically least resulting image
    vector is returned.  For transitive tuples this is a complete invariant of
    simultaneous sheet relabeling; the paper contains the proof.
    """
    generators = list(t)
    traversal_generators = generators + [inverse(g) for g in generators]
    best_key: tuple[int, ...] | None = None
    best_tuple: Tuple4 | None = None

    for first_sheet in range(1, N + 1):
        old_to_new = [0] * (N + 1)
        old_to_new[first_sheet] = 1
        queue = [first_sheet]
        next_label = 2
        cursor = 0

        while cursor < len(queue):
            x = queue[cursor]
            cursor += 1
            for g in traversal_generators:
                y = g[x]
                if old_to_new[y] == 0:
                    old_to_new[y] = next_label
                    next_label += 1
                    queue.append(y)

        if next_label != N + 1:
            continue

        images: list[Perm] = []
        flat: list[int] = []
        for g in generators:
            relabeled = [0] * (N + 1)
            for old in range(1, N + 1):
                relabeled[old_to_new[old]] = old_to_new[g[old]]
            relabeled_tuple = tuple(relabeled)
            images.append(relabeled_tuple)
            flat.extend(relabeled[1:])

        key = tuple(flat)
        if best_key is None or key < best_key:
            best_key = key
            best_tuple = tuple(images)  # type: ignore[assignment]

    if best_key is None or best_tuple is None:
        raise RuntimeError("canonicalization failed; tuple is not transitive")
    return best_tuple, best_key


def build_orbit(seed: Tuple4) -> tuple[list[Tuple4], dict[tuple[int, ...], int]]:
    canonical_seed, seed_key = canonical_tuple_and_key(seed)
    orbit = [canonical_seed]
    index = {seed_key: 0}
    cursor = 0
    while cursor < len(orbit):
        t = orbit[cursor]
        cursor += 1
        for word in ALL_WORDS:
            canonical, key = canonical_tuple_and_key(apply_word(t, word))
            if key not in index:
                index[key] = len(orbit)
                orbit.append(canonical)
    return orbit, index


def build_operator(
    orbit: Sequence[Tuple4],
    index: dict[tuple[int, ...], int],
    word: Word,
) -> tuple[int, ...]:
    images: list[int] = []
    for t in orbit:
        _, key = canonical_tuple_and_key(apply_word(t, word))
        images.append(index[key])
    return tuple(images)


def compose_zero_based(a: Sequence[int], b: Sequence[int]) -> tuple[int, ...]:
    """Right-action product on {0,...,n-1}."""
    return tuple(b[a[i]] for i in range(len(a)))


def inverse_zero_based(p: Sequence[int]) -> tuple[int, ...]:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def orbit_components(operators: Sequence[Sequence[int]]) -> list[list[int]]:
    n = len(operators[0])
    unseen = set(range(n))
    components: list[list[int]] = []
    generators = list(operators) + [inverse_zero_based(op) for op in operators]
    while unseen:
        start = min(unseen)
        seen = {start}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for op in generators:
                y = op[x]
                if y not in seen:
                    seen.add(y)
                    queue.append(y)
        components.append(sorted(seen))
        unseen.difference_update(seen)
    return sorted(components, key=lambda c: (-len(c), c[0]))


def cycle_lengths_zero_based(p: Sequence[int]) -> list[int]:
    used = [False] * len(p)
    lengths: list[int] = []
    for i in range(len(p)):
        if used[i]:
            continue
        x = i
        length = 0
        while not used[x]:
            used[x] = True
            length += 1
            x = p[x]
        lengths.append(length)
    return sorted(lengths)


def fingerprint_tuple(t: Tuple4) -> str:
    payload = json.dumps(
        [[g[i] for i in range(1, N + 1)] for g in t],
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def hash_orbit_keys(keys: Iterable[tuple[int, ...]]) -> str:
    h = hashlib.sha256()
    for key in sorted(keys):
        h.update(bytes(key))  # entries lie in {1,...,23}; key length is fixed.
    return h.hexdigest()


def hash_operator(operator: Sequence[int]) -> str:
    h = hashlib.sha256()
    for image in operator:
        h.update(struct.pack(">H", image))
    return h.hexdigest()


def hash_operators(operators: Sequence[Sequence[int]]) -> str:
    h = hashlib.sha256()
    for operator in operators:
        for image in operator:
            h.update(struct.pack(">H", image))
    return h.hexdigest()


def braid_relation_checks(orbit: Sequence[Tuple4]) -> dict[str, bool]:
    q1q2q1: Word = ((0, 1), (1, 1), (0, 1))
    q2q1q2: Word = ((1, 1), (0, 1), (1, 1))
    q2q3q2: Word = ((1, 1), (2, 1), (1, 1))
    q3q2q3: Word = ((2, 1), (1, 1), (2, 1))
    q1q3: Word = ((0, 1), (2, 1))
    q3q1: Word = ((2, 1), (0, 1))

    inverse_ok = True
    artin_12_ok = True
    artin_23_ok = True
    far_ok = True
    product_ok = True
    class_vector_ok = True
    peripheral_relation_ok = True

    for t in orbit:
        for i in range(3):
            inverse_ok &= qmove(qmove(t, i, 1), i, -1) == t
            inverse_ok &= qmove(qmove(t, i, -1), i, 1) == t
        artin_12_ok &= apply_word(t, q1q2q1) == apply_word(t, q2q1q2)
        artin_23_ok &= apply_word(t, q2q3q2) == apply_word(t, q3q2q3)
        far_ok &= apply_word(t, q1q3) == apply_word(t, q3q1)
        product_ok &= tuple_product(t) == ID
        class_vector_ok &= [cycle_type(cycle_lengths(g)) for g in t] == EXPECTED_TUPLE_TYPES

        lhs = apply_word(apply_word(apply_word(t, B12), B13), B14)
        g1 = t[0]
        rhs: Tuple4 = tuple(conjugate_by_left(g1, g) for g in t)  # type: ignore[assignment]
        peripheral_relation_ok &= lhs == rhs

    return {
        "qmove_inverse_relations_all_states": inverse_ok,
        "artin_q1q2q1_all_states": artin_12_ok,
        "artin_q2q3q2_all_states": artin_23_ok,
        "far_commutativity_q1q3_all_states": far_ok,
        "product_one_all_states": product_ok,
        "ordered_class_vector_all_states": class_vector_ok,
        "peripheral_product_is_inner_conjugation_all_states": peripheral_relation_ok,
    }


def canonicalization_self_check(seed: Tuple4) -> bool:
    conjugators = [
        ATLAS_A,
        ATLAS_B,
        compose(ATLAS_A, ATLAS_B),
        compose(ATLAS_B, ATLAS_A),
        power(compose(ATLAS_A, ATLAS_B), 7),
        seed[0],
        seed[2],
        compose(seed[1], seed[3]),
    ]
    _, seed_key = canonical_tuple_and_key(seed)
    return all(canonical_tuple_and_key(tuple_conjugate(seed, h))[1] == seed_key for h in conjugators)


def script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def verify() -> tuple[dict[str, object], dict[str, object]]:
    started = time.perf_counter()

    atlas_group = PermutationGroup([to_sympy(ATLAS_A), to_sympy(ATLAS_B)])
    tuple_group = PermutationGroup([to_sympy(g) for g in TUPLE])

    atlas_order = int(atlas_group.order())
    tuple_group_order = int(tuple_group.order())
    atlas_generator_orders = [permutation_order(ATLAS_A), permutation_order(ATLAS_B)]
    atlas_ab_order = permutation_order(compose(ATLAS_A, ATLAS_B))
    atlas_transitive = bool(atlas_group.is_transitive())
    point_stabilizer = atlas_group.stabilizer(0)
    point_stabilizer_order = int(point_stabilizer.order())
    point_stabilizer_orbits = sorted(len(orbit) for orbit in point_stabilizer.orbits())
    atlas_two_transitive = point_stabilizer_orbits == [1, N - 1]

    stabilizer_chain_orders: list[int] = []
    stabilizer_chain_orbit_sizes: list[list[int]] = []
    chain_group = atlas_group
    for fixed_count in range(5):
        stabilizer_chain_orders.append(int(chain_group.order()))
        stabilizer_chain_orbit_sizes.append(sorted(len(o) for o in chain_group.orbits()))
        if fixed_count < 4:
            chain_group = chain_group.stabilizer(fixed_count)
    atlas_four_transitive = stabilizer_chain_orbit_sizes[3] == [1, 1, 1, N - 3]
    atlas_five_transitive = stabilizer_chain_orbit_sizes[4] == [1, 1, 1, 1, N - 4]

    tuple_membership = [bool(atlas_group.contains(to_sympy(g))) for g in TUPLE]
    tuple_types = [cycle_type(cycle_lengths(g)) for g in TUPLE]
    tuple_orders = [permutation_order(g) for g in TUPLE]
    tuple_indices = [permutation_index(g) for g in TUPLE]
    source_genus = 1 - N + sum(tuple_indices) // 2

    orbit, index = build_orbit(TUPLE)
    forward_operators = [build_operator(orbit, index, word) for word in FORWARD_WORDS]
    inverse_operators = [build_operator(orbit, index, word) for word in INVERSE_WORDS]
    q1_operator = build_operator(orbit, index, Q1)
    q1_lengths = cycle_lengths_zero_based(q1_operator)
    q1_type = cycle_type(q1_lengths)
    q1_order = 1
    for length in q1_lengths:
        q1_order = q1_order * length // gcd(q1_order, length)
    q1_square_is_b12 = compose_zero_based(q1_operator, q1_operator) == forward_operators[0]

    _, seed_key = canonical_tuple_and_key(TUPLE)
    seed_index_zero_based = index[seed_key]
    seed_cycle_lengths: list[int] = []
    for operator in forward_operators:
        x = seed_index_zero_based
        length = 1
        y = operator[x]
        while y != x:
            length += 1
            y = operator[y]
        seed_cycle_lengths.append(length)

    identity = tuple(range(len(orbit)))

    operator_inverse_checks = [
        compose_zero_based(forward, backward) == identity
        and compose_zero_based(backward, forward) == identity
        for forward, backward in zip(forward_operators, inverse_operators)
    ]
    product = compose_zero_based(
        compose_zero_based(forward_operators[0], forward_operators[1]),
        forward_operators[2],
    )
    components = orbit_components(forward_operators)
    operator_lengths = [cycle_lengths_zero_based(op) for op in forward_operators]
    operator_types = [cycle_type(lengths) for lengths in operator_lengths]
    operator_indices = [sum(length - 1 for length in lengths) for lengths in operator_lengths]
    hurwitz_genus = 1 - len(orbit) + sum(operator_indices) // 2

    relation_checks = braid_relation_checks(orbit)
    canonicalization_check = canonicalization_self_check(TUPLE)

    closure_under_six_words = True
    for t in orbit:
        for word in ALL_WORDS:
            _, key = canonical_tuple_and_key(apply_word(t, word))
            if key not in index:
                closure_under_six_words = False
                break
        if not closure_under_six_words:
            break

    checks = {
        "atlas_order": atlas_order == M23_ORDER,
        "atlas_generator_orders_2_4": atlas_generator_orders == [2, 4],
        "atlas_product_order_23": atlas_ab_order == 23,
        "atlas_transitive": atlas_transitive,
        "atlas_two_transitive": atlas_two_transitive,
        "atlas_four_transitive": atlas_four_transitive,
        "atlas_not_five_transitive": not atlas_five_transitive,
        "point_stabilizer_order": point_stabilizer_order == M22_ORDER,
        "tuple_membership": tuple_membership == [True, True, True, True],
        "tuple_generated_group_order": tuple_group_order == M23_ORDER,
        "tuple_product_one": tuple_product(TUPLE) == ID,
        "tuple_orders": tuple_orders == [2, 2, 3, 5],
        "tuple_cycle_types": tuple_types == EXPECTED_TUPLE_TYPES,
        "tuple_indices": tuple_indices == [8, 8, 12, 16],
        "source_genus": source_genus == 0,
        "canonicalization_conjugacy_self_check": canonicalization_check,
        "orbit_size": len(orbit) == 980,
        "orbit_closed_under_six_point_pushes": closure_under_six_words,
        "operator_inverse_checks": operator_inverse_checks == [True, True, True],
        "operator_product_one_on_inner_classes": product == identity,
        "operator_connectedness": [len(c) for c in components] == [980],
        "operator_cycle_types": operator_types == EXPECTED_OPERATOR_TYPES,
        "operator_indices": operator_indices == [612, 792, 792],
        "equal_slot_half_twist_cycle_type": q1_type == {2: 12, 3: 38, 4: 68, 5: 30, 6: 70},
        "equal_slot_half_twist_order": q1_order == 60,
        "equal_slot_half_twist_square_is_b12": q1_square_is_b12,
        "printed_seed_index_one_based": seed_index_zero_based + 1 == 1,
        "printed_seed_peripheral_cycle_lengths": seed_cycle_lengths == [1, 11, 11],
        "compactified_genus": hurwitz_genus == 119,
        **relation_checks,
    }
    passed = all(checks.values())

    elapsed = time.perf_counter() - started
    peak_rss_kb: int | None = None
    if resource is not None:
        peak_rss_kb = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    invariant_data = {
        "conventions": {
            "permutation_action": "right",
            "product": "i^(xy)=(i^x)^y",
            "hurwitz_move": "q_i(a,b)=(a*b*a^-1,a)",
            "base": "ordered M_0,4 = P^1 minus {0,1,infinity}",
            "equivalence": "simultaneous S_23-conjugacy; equal to inner equivalence for this action",
        },
        "atlas_copy": {
            "degree": N,
            "order": atlas_order,
            "generator_orders": atlas_generator_orders,
            "product_order": atlas_ab_order,
            "transitive": atlas_transitive,
            "two_transitive": atlas_two_transitive,
            "four_transitive": atlas_four_transitive,
            "five_transitive": atlas_five_transitive,
            "point_stabilizer_order": point_stabilizer_order,
            "point_stabilizer_orbit_sizes": point_stabilizer_orbits,
            "stabilizer_chain_orders": stabilizer_chain_orders,
            "stabilizer_chain_orbit_sizes": stabilizer_chain_orbit_sizes,
        },
        "degree_23_tuple": {
            "fingerprint_sha256": fingerprint_tuple(TUPLE),
            "membership_in_atlas_copy": tuple_membership,
            "generated_group_order": tuple_group_order,
            "product_one": tuple_product(TUPLE) == ID,
            "orders": tuple_orders,
            "cycle_types": tuple_types,
            "indices": tuple_indices,
            "source_genus": source_genus,
        },
        "ordered_hurwitz_orbit": {
            "size": len(orbit),
            "operator_names": list(WORD_NAMES),
            "operator_product_one_on_inner_classes": product == identity,
            "operator_inverse_checks": operator_inverse_checks,
            "component_sizes": [len(c) for c in components],
            "cycle_types": operator_types,
            "indices": operator_indices,
            "total_index": sum(operator_indices),
            "compactified_genus": hurwitz_genus,
            "orbit_keys_sha256": hash_orbit_keys(index.keys()),
            "operator_sha256": [hash_operator(op) for op in forward_operators],
            "operators_combined_sha256": hash_operators(forward_operators),
            "equal_slot_half_twist": {
                "word": "q1",
                "cycle_type": q1_type,
                "order": q1_order,
                "square_is_b12": q1_square_is_b12,
            },
            "printed_seed": {
                "index_one_based": seed_index_zero_based + 1,
                "peripheral_cycle_lengths": dict(zip(WORD_NAMES, seed_cycle_lengths)),
            },
        },
        "checks": checks,
        "pass": passed,
    }

    run_data = {
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "sympy": sympy.__version__,
            "platform": platform.platform(),
            "script_sha256": script_sha256(),
        },
        "performance": {
            "elapsed_seconds": round(elapsed, 6),
            "peak_rss_kb": peak_rss_kb,
        },
    }

    dump_data = {
        "metadata": {
            "orbit_size": len(orbit),
            "key_length": 4 * N,
            "operator_names": list(WORD_NAMES),
            "indexing": "zero-based for orbit operators; one-based for degree-23 permutations",
        },
        "sorted_canonical_keys": [list(key) for key in sorted(index.keys())],
        "forward_operators": [list(op) for op in forward_operators],
        "inverse_operators": [list(op) for op in inverse_operators],
    }

    return {**invariant_data, **run_data}, dump_data


def format_cycle_type(value: dict[int, int]) -> str:
    return " ".join(f"{length}^{multiplicity}" for length, multiplicity in value.items())


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_gzip_json(path: Path, value: object) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="write full verification result as JSON")
    parser.add_argument(
        "--dump-data",
        type=Path,
        help="write canonical keys and degree-980 operators as compressed JSON",
    )
    args = parser.parse_args()

    result, dump_data = verify()
    degree = result["degree_23_tuple"]
    hurwitz = result["ordered_hurwitz_orbit"]
    atlas = result["atlas_copy"]

    print("M23_ORDERED_HURWITZ_CERTIFICATE")
    print("PASS", result["pass"])
    print("PYTHON", result["environment"]["python"])
    print("SYMPY", result["environment"]["sympy"])
    print("SCRIPT_SHA256", result["environment"]["script_sha256"])
    print("ATLAS_GROUP_ORDER", atlas["order"])
    print("ATLAS_GENERATOR_ORDERS", atlas["generator_orders"])
    print("ATLAS_AB_ORDER", atlas["product_order"])
    print("ATLAS_TWO_TRANSITIVE", atlas["two_transitive"])
    print("ATLAS_FOUR_TRANSITIVE", atlas["four_transitive"])
    print("ATLAS_FIVE_TRANSITIVE", atlas["five_transitive"])
    print("STABILIZER_CHAIN_ORDERS", atlas["stabilizer_chain_orders"])
    print("POINT_STABILIZER_ORDER", atlas["point_stabilizer_order"])
    print("TUPLE_MEMBERSHIP", degree["membership_in_atlas_copy"])
    print("TUPLE_GENERATED_GROUP_ORDER", degree["generated_group_order"])
    print("TUPLE_PRODUCT_ONE", degree["product_one"])
    print("TUPLE_ORDERS", degree["orders"])
    print("TUPLE_INDICES", degree["indices"])
    print("SOURCE_GENUS", degree["source_genus"])
    print("ORDERED_HURWITZ_ORBIT_SIZE", hurwitz["size"])
    print("OPERATOR_PRODUCT_ONE", hurwitz["operator_product_one_on_inner_classes"])
    print("OPERATOR_INVERSE_CHECKS", hurwitz["operator_inverse_checks"])
    print("OPERATOR_COMPONENT_SIZES", hurwitz["component_sizes"])
    print("EQUAL_SLOT_HALF_TWIST_CYCLE_TYPE", format_cycle_type(hurwitz["equal_slot_half_twist"]["cycle_type"]))
    print("EQUAL_SLOT_HALF_TWIST_ORDER", hurwitz["equal_slot_half_twist"]["order"])
    print("EQUAL_SLOT_HALF_TWIST_SQUARE_IS_B12", hurwitz["equal_slot_half_twist"]["square_is_b12"])
    print("PRINTED_SEED_INDEX_ONE_BASED", hurwitz["printed_seed"]["index_one_based"])
    print("PRINTED_SEED_POINT_PUSH_LENGTHS", hurwitz["printed_seed"]["peripheral_cycle_lengths"])
    for name, value, digest in zip(
        hurwitz["operator_names"],
        hurwitz["cycle_types"],
        hurwitz["operator_sha256"],
    ):
        print(name + "_CYCLE_TYPE", format_cycle_type(value))
        print(name + "_SHA256", digest)
    print("ORBIT_KEYS_SHA256", hurwitz["orbit_keys_sha256"])
    print("OPERATORS_COMBINED_SHA256", hurwitz["operators_combined_sha256"])
    print("OPERATOR_INDICES", hurwitz["indices"])
    print("TOTAL_INDEX", hurwitz["total_index"])
    print("COMPACTIFIED_HURWITZ_GENUS", hurwitz["compactified_genus"])
    print("TUPLE_FINGERPRINT_SHA256", degree["fingerprint_sha256"])
    print("ELAPSED_SECONDS", result["performance"]["elapsed_seconds"])
    print("PEAK_RSS_KB", result["performance"]["peak_rss_kb"])

    failed_checks = [name for name, ok in result["checks"].items() if not ok]
    print("FAILED_CHECKS", failed_checks)

    if args.json:
        write_json(args.json, result)
        print("JSON", args.json)
    if args.dump_data:
        write_gzip_json(args.dump_data, dump_data)
        print("DUMP_DATA", args.dump_data)

    return 0 if result["pass"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
