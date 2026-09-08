#!/usr/bin/env python3
"""Exact exhaustiveness verification for the ordered (2A,2A,3A,5A) M23 Nielsen class.

The program starts from an explicit degree-23 copy of M23 and the explicit
product-one tuple printed in the accompanying article.  It then:

1. constructs the complete 2A and 3A conjugacy classes by conjugation closure;
2. fixes the displayed element d in 5A;
3. exhaustively tests every ordered pair (x,y) in 2A x 2A, solving
       c = (x y)^(-1) d^(-1),
   and counts the solutions with c in 3A;
4. counts which solutions are transitive on the 23 sheets;
5. computes |C_M23(d)| and the resulting number of generating inner classes.

All group operations and membership tests are exact.  NumPy is used only to
batch the finite permutation-array comparisons; two independent exact byte-key
membership methods are compared on every batch.

Tested with Python 3.13.5, NumPy 2.3.1 and SymPy 1.14.0.
"""

from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence
import argparse
import hashlib
import json
import platform
import struct
import sys
import time

import numpy as np
import sympy
from sympy.combinatorics import Permutation, PermutationGroup

try:
    import resource
except ImportError:  # pragma: no cover
    resource = None  # type: ignore[assignment]

N = 23
M23_ORDER = 10_200_960
ID = tuple(range(N + 1))
Perm = tuple[int, ...]
Tuple4 = tuple[Perm, Perm, Perm, Perm]


def perm_from_cycles(cycles: Iterable[Iterable[int]]) -> Perm:
    p = list(range(N + 1))
    for raw in cycles:
        cycle = list(raw)
        if len(cycle) < 2:
            continue
        for i, x in enumerate(cycle):
            p[x] = cycle[(i + 1) % len(cycle)]
    return tuple(p)


def compose(a: Perm, b: Perm) -> Perm:
    """Right-action product: i^(a b)=(i^a)^b."""
    return tuple([0] + [b[a[i]] for i in range(1, N + 1)])


@lru_cache(maxsize=None)
def inverse(a: Perm) -> Perm:
    out = [0] * (N + 1)
    for i in range(1, N + 1):
        out[a[i]] = i
    return tuple(out)


def conjugate(g: Perm, h: Perm) -> Perm:
    """Return h^(-1) g h."""
    return compose(compose(inverse(h), g), h)


def tuple_product(t: Tuple4) -> Perm:
    out = ID
    for g in t:
        out = compose(out, g)
    return out


def cycle_lengths(p: Sequence[int]) -> list[int]:
    used = [False] * len(p)
    lengths: list[int] = []
    for start in range(1, len(p)):
        if used[start]:
            continue
        x = start
        length = 0
        while not used[x]:
            used[x] = True
            x = p[x]
            length += 1
        lengths.append(length)
    return sorted(lengths)


def cycle_type(p: Sequence[int]) -> dict[int, int]:
    return dict(sorted(Counter(cycle_lengths(p)).items()))


def permutation_order(p: Perm) -> int:
    value = 1
    for length in cycle_lengths(p):
        value = value * length // gcd(value, length)
    return value


def to_sympy(p: Perm) -> Permutation:
    return Permutation([p[i] - 1 for i in range(1, N + 1)])


# Provenance: ATLAS of Finite Groups / ATLAS online M23 natural
# permutation representation; the companion ordered-Hurwitz
# certificate checks the standard-generator signatures.
ATLAS_A = perm_from_cycles(
    [(1, 2), (3, 4), (7, 8), (9, 10),
     (13, 14), (15, 16), (19, 20), (21, 22)]
)
ATLAS_B = perm_from_cycles(
    [(1, 16, 11, 3), (2, 9, 21, 12), (4, 5, 8, 23),
     (6, 22, 14, 18), (13, 20), (15, 17)]
)

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

EXPECTED_TYPES = [
    {1: 7, 2: 8},
    {1: 7, 2: 8},
    {1: 5, 3: 6},
    {1: 3, 5: 4},
]


def conjugacy_class(rep: Perm) -> list[Perm]:
    """Conjugacy closure under generators of the fixed M23 copy."""
    steps = (ATLAS_A, ATLAS_B, inverse(ATLAS_A), inverse(ATLAS_B))
    seen = {rep}
    queue = deque([rep])
    while queue:
        g = queue.popleft()
        for s in steps:
            h = conjugate(g, s)
            if h not in seen:
                seen.add(h)
                queue.append(h)
    return sorted(seen)


def as_numpy(permutations: Sequence[Perm]) -> np.ndarray:
    """Convert one-based tuple permutations to zero-based uint8 image arrays."""
    return np.asarray(
        [[p[i] - 1 for i in range(1, N + 1)] for p in permutations],
        dtype=np.uint8,
    )


def transitive_on_23(x: np.ndarray, y: np.ndarray, d: np.ndarray, dinv: np.ndarray) -> bool:
    all_bits = (1 << N) - 1
    seen = 1
    stack = [0]
    while stack:
        point = stack.pop()
        for generator in (x, y, d, dinv):
            image = int(generator[point])
            bit = 1 << image
            if not (seen & bit):
                seen |= bit
                stack.append(image)
    return seen == all_bits


def hash_pairs(pairs: Sequence[tuple[int, int]]) -> str:
    h = hashlib.sha256()
    for i, j in pairs:
        h.update(struct.pack(">HH", i, j))
    return h.hexdigest()


def script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def verify(batch_size: int = 32) -> dict[str, object]:
    started = time.perf_counter()

    group = PermutationGroup([to_sympy(ATLAS_A), to_sympy(ATLAS_B)])
    tuple_group = PermutationGroup([to_sympy(g) for g in TUPLE])
    group_order = int(group.order())
    tuple_group_order = int(tuple_group.order())

    class_2a = conjugacy_class(TUPLE[0])
    class_3a = conjugacy_class(TUPLE[2])
    class_2a_size = len(class_2a)
    class_3a_size = len(class_3a)

    x_array = as_numpy(class_2a)
    c3_array = as_numpy(class_3a)
    d_array = as_numpy([TUPLE[3]])[0]
    d_inverse = np.argsort(d_array).astype(np.uint8)

    byte_key_dtype = np.dtype((np.void, N))
    c3_keys = np.ascontiguousarray(c3_array).view(byte_key_dtype).ravel()
    c3_keys_sorted = np.sort(c3_keys)

    # Confirm the solved formula on the printed seed before the exhaustive pass.
    seed_x = np.asarray([TUPLE[0][i] - 1 for i in range(1, N + 1)], dtype=np.uint8)
    seed_y = np.asarray([TUPLE[1][i] - 1 for i in range(1, N + 1)], dtype=np.uint8)
    seed_c = np.asarray([TUPLE[2][i] - 1 for i in range(1, N + 1)], dtype=np.uint8)
    seed_xy = seed_y[seed_x]
    seed_candidate = d_inverse[np.argsort(seed_xy).astype(np.uint8)]
    seed_formula_ok = bool(np.array_equal(seed_candidate, seed_c))

    solution_pairs: list[tuple[int, int]] = []
    dual_membership_methods_agree = True
    m = class_2a_size

    # There are 3795^2 = 14,402,025 ordered pairs.  For a batch of x values,
    # xy is obtained exactly as y[x] under the stated right-action convention.
    for start in range(0, m, batch_size):
        x_batch = x_array[start:start + batch_size]
        xy = x_array[:, x_batch].transpose(1, 0, 2)
        inverse_xy = np.argsort(xy, axis=2).astype(np.uint8)
        candidates = d_inverse[inverse_xy]
        candidate_keys = (
            np.ascontiguousarray(candidates)
            .reshape(-1, N)
            .view(byte_key_dtype)
            .ravel()
        )

        mask_isin = np.isin(candidate_keys, c3_keys, assume_unique=False)
        positions = np.searchsorted(c3_keys_sorted, candidate_keys)
        mask_search = np.zeros(candidate_keys.shape[0], dtype=bool)
        in_range = positions < c3_keys_sorted.shape[0]
        mask_search[in_range] = c3_keys_sorted[positions[in_range]] == candidate_keys[in_range]
        dual_membership_methods_agree &= bool(np.array_equal(mask_isin, mask_search))

        mask = mask_isin.reshape(len(x_batch), m)
        rows, columns = np.nonzero(mask)
        solution_pairs.extend((int(row + start), int(column)) for row, column in zip(rows, columns))

    transitive_pairs: list[tuple[int, int]] = []
    for i, j in solution_pairs:
        if transitive_on_23(x_array[i], x_array[j], d_array, d_inverse):
            transitive_pairs.append((i, j))

    centralizer = group.centralizer(PermutationGroup([to_sympy(TUPLE[3])]))
    centralizer_order = int(centralizer.order())
    normalized_transitive_count = (
        len(transitive_pairs) // centralizer_order
        if centralizer_order and len(transitive_pairs) % centralizer_order == 0
        else None
    )

    checks = {
        "atlas_group_order": group_order == M23_ORDER,
        "printed_tuple_generates_m23": tuple_group_order == M23_ORDER,
        "printed_tuple_product_one": tuple_product(TUPLE) == ID,
        "printed_tuple_cycle_types": [cycle_type(g) for g in TUPLE] == EXPECTED_TYPES,
        "seed_solution_formula": seed_formula_ok,
        "class_2a_size": class_2a_size == 3_795,
        "class_3a_size": class_3a_size == 56_672,
        "class_2a_unique": len(set(class_2a)) == class_2a_size,
        "class_3a_unique": len(set(class_3a)) == class_3a_size,
        "dual_membership_methods_agree": dual_membership_methods_agree,
        "fixed_d_product_one_solutions": len(solution_pairs) == 127_200,
        "fixed_d_transitive_solutions": len(transitive_pairs) == 14_700,
        "fixed_d_intransitive_solutions": len(solution_pairs) - len(transitive_pairs) == 112_500,
        "centralizer_order": centralizer_order == 15,
        "normalized_transitive_count": normalized_transitive_count == 980,
    }

    elapsed = time.perf_counter() - started
    peak_rss_kb: int | None = None
    if resource is not None:
        peak_rss_kb = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    return {
        "pass": all(checks.values()),
        "checks": checks,
        "conventions": {
            "permutation_action": "right",
            "product": "i^(xy)=(i^x)^y",
            "fixed_fourth_entry": "the printed 5A element d",
            "solved_entry": "c=(xy)^(-1)d^(-1)",
        },
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "sympy": sympy.__version__,
            "platform": platform.platform(),
            "script_sha256": script_sha256(),
        },
        "group": {
            "degree": N,
            "order": group_order,
            "printed_tuple_generated_group_order": tuple_group_order,
            "printed_tuple_orders": [permutation_order(g) for g in TUPLE],
            "printed_tuple_cycle_types": [cycle_type(g) for g in TUPLE],
        },
        "class_enumeration": {
            "class_2a_size": class_2a_size,
            "class_3a_size": class_3a_size,
            "ordered_pairs_tested": class_2a_size * class_2a_size,
            "fixed_d_product_one_solutions": len(solution_pairs),
            "fixed_d_transitive_solutions": len(transitive_pairs),
            "fixed_d_intransitive_solutions": len(solution_pairs) - len(transitive_pairs),
            "centralizer_order": centralizer_order,
            "normalized_transitive_count": normalized_transitive_count,
            "solution_pair_sha256": hash_pairs(solution_pairs),
            "transitive_pair_sha256": hash_pairs(transitive_pairs),
        },
        "performance": {
            "batch_size": batch_size,
            "elapsed_seconds": round(elapsed, 6),
            "peak_rss_kb": peak_rss_kb,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    if args.batch_size <= 0:
        parser.error("--batch-size must be positive")

    result = verify(batch_size=args.batch_size)
    counts = result["class_enumeration"]
    environment = result["environment"]

    print("M23_NIELSEN_EXHAUSTIVENESS_CERTIFICATE")
    print("PASS", result["pass"])
    print("PYTHON", environment["python"])
    print("NUMPY", environment["numpy"])
    print("SYMPY", environment["sympy"])
    print("SCRIPT_SHA256", environment["script_sha256"])
    print("CLASS_2A_SIZE", counts["class_2a_size"])
    print("CLASS_3A_SIZE", counts["class_3a_size"])
    print("ORDERED_PAIRS_TESTED", counts["ordered_pairs_tested"])
    print("FIXED_D_PRODUCT_ONE_SOLUTIONS", counts["fixed_d_product_one_solutions"])
    print("FIXED_D_TRANSITIVE_SOLUTIONS", counts["fixed_d_transitive_solutions"])
    print("FIXED_D_INTRANSITIVE_SOLUTIONS", counts["fixed_d_intransitive_solutions"])
    print("CENTRALIZER_ORDER", counts["centralizer_order"])
    print("NORMALIZED_TRANSITIVE_COUNT", counts["normalized_transitive_count"])
    print("SOLUTION_PAIR_SHA256", counts["solution_pair_sha256"])
    print("TRANSITIVE_PAIR_SHA256", counts["transitive_pair_sha256"])
    print("ELAPSED_SECONDS", result["performance"]["elapsed_seconds"])
    print("PEAK_RSS_KB", result["performance"]["peak_rss_kb"])
    failed = [name for name, ok in result["checks"].items() if not ok]
    print("FAILED_CHECKS", failed)

    if args.json:
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print("JSON", args.json)

    return 0 if result["pass"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
