#!/usr/bin/env python3
"""Construct the full reduced j-line monodromy for the M23 class multiset.

The construction starts from the explicit (2A,2A,3A,5A) tuple used by the
primary certificate. It closes the full Hurwitz-braid orbit across all twelve
orderings of the class multiset, quotients by the standard Klein four group
Q'' = <q1 q3^-1, sh^2>, and constructs the reduced generators

    gamma_1 = sh,
    gamma_infinity = q2,
    gamma_0 = (gamma_1 gamma_infinity)^-1.

For the right-action convention used in the bundle these become

    gamma_0 = q1 q2,
    gamma_1 = q1 q2 q1,
    gamma_infinity = q2

on the Q'' quotient. The script also identifies the 980-state equal-class
symmetrized system as one of the three ordering-color slices of the resulting
2940-state reduced action.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
from typing import Sequence

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[1]
CODE = SCRIPT.parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import sympy
import verify_m23_ordered_hurwitz as base

Tuple4 = base.Tuple4
Perm = tuple[int, ...]

Q_GENERATORS = [
    base.Q1, base.Q1_INV,
    base.Q2, base.Q2_INV,
    base.Q3, base.Q3_INV,
]


def compose(a: Sequence[int], b: Sequence[int]) -> Perm:
    """Right-action composition: apply a, then b."""
    return tuple(b[a[i]] for i in range(len(a)))


def power(p: Sequence[int], exponent: int) -> Perm:
    n = len(p)
    if exponent < 0:
        inv = [0] * n
        for i, image in enumerate(p):
            inv[image] = i
        return power(tuple(inv), -exponent)
    result: Perm = tuple(range(n))
    base_p: Perm = tuple(p)
    e = exponent
    while e:
        if e & 1:
            result = compose(result, base_p)
        base_p = compose(base_p, base_p)
        e >>= 1
    return result


def cycle_lengths(p: Sequence[int]) -> list[int]:
    seen = [False] * len(p)
    out: list[int] = []
    for start in range(len(p)):
        if seen[start]:
            continue
        x = start
        length = 0
        while not seen[x]:
            seen[x] = True
            length += 1
            x = p[x]
        out.append(length)
    return sorted(out)


def cycle_type(p: Sequence[int]) -> dict[int, int]:
    return dict(sorted(Counter(cycle_lengths(p)).items()))


def permutation_order(p: Sequence[int]) -> int:
    return math.lcm(*cycle_lengths(p))


def operator_sha256(p: Sequence[int]) -> str:
    digest = hashlib.sha256()
    for value in p:
        digest.update(int(value).to_bytes(4, "little"))
    return digest.hexdigest()


def script_sha256() -> str:
    return hashlib.sha256(SCRIPT.read_bytes()).hexdigest()


def ordering_representatives() -> dict[tuple[int, int, int, int], base.Word]:
    """Choose one adjacent-braid word for each ordering of (2,2,3,5)."""
    start_ordering = (2, 2, 3, 5)
    words = [base.Q1, base.Q2, base.Q3]
    representatives: dict[tuple[int, int, int, int], base.Word] = {
        start_ordering: (),
    }
    queue = deque([start_ordering])
    while queue:
        ordering = queue.popleft()
        for i, generator in enumerate(words):
            target_list = list(ordering)
            target_list[i], target_list[i + 1] = target_list[i + 1], target_list[i]
            target = tuple(target_list)
            if target not in representatives:
                representatives[target] = representatives[ordering] + generator
                queue.append(target)
    if len(representatives) != 12:
        raise RuntimeError("expected twelve orderings of the class multiset")
    return representatives


def build_full_orbit() -> tuple[list[Tuple4], dict[tuple[int, ...], int]]:
    """Build all 11,760 states from the complete 980-state ordered class.

    The pure-braid orbit for the initial ordering is already the complete
    ordered generating Nielsen class.  A braid word carrying the initial
    ordering to each of the other eleven orderings gives a bijection of
    Nielsen classes.  Transporting all 980 states through these twelve words
    therefore constructs the full braid orbit without an expensive global
    breadth-first search.
    """
    ordered_orbit, _ = base.build_orbit(base.TUPLE)
    base.canonical_tuple_and_key.cache_clear()
    base.inverse.cache_clear()

    orbit: list[Tuple4] = []
    index: dict[tuple[int, ...], int] = {}
    for ordering, word in ordering_representatives().items():
        start = len(orbit)
        for state in ordered_orbit:
            canonical, key = base.canonical_tuple_and_key(base.apply_word(state, word))
            if key in index:
                raise RuntimeError(
                    f"duplicate canonical state while constructing ordering {ordering}"
                )
            index[key] = len(orbit)
            orbit.append(canonical)
        if len(orbit) - start != 980:
            raise RuntimeError(f"ordering {ordering} does not contain 980 states")
        base.canonical_tuple_and_key.cache_clear()
        base.inverse.cache_clear()

    if len(orbit) != 11760:
        raise RuntimeError("full all-orderings orbit does not have 11,760 states")
    return orbit, index


def build_q_operators(
    orbit: list[Tuple4],
    index: dict[tuple[int, ...], int],
) -> list[list[int]]:
    operators = [[], [], []]
    words = [base.Q1, base.Q2, base.Q3]
    for i, tuple_state in enumerate(orbit):
        for operator_index, word in enumerate(words):
            transformed = base.apply_word(tuple_state, word)
            _, key = base.canonical_tuple_and_key(transformed)
            try:
                operators[operator_index].append(index[key])
            except KeyError as exc:
                raise RuntimeError("full all-orderings set is not braid-stable") from exc
        if (i + 1) % 500 == 0:
            base.canonical_tuple_and_key.cache_clear()
            base.inverse.cache_clear()
    base.canonical_tuple_and_key.cache_clear()
    base.inverse.cache_clear()
    return operators

def orbit_partition(n: int, operators: Sequence[Sequence[int]]) -> tuple[list[list[int]], list[int]]:
    inverses: list[Perm] = []
    for operator in operators:
        inverse = [0] * n
        for i, image in enumerate(operator):
            inverse[image] = i
        inverses.append(tuple(inverse))

    all_operators = [tuple(operator) for operator in operators] + inverses
    seen = [False] * n
    blocks: list[list[int]] = []
    block_of = [-1] * n

    for start in range(n):
        if seen[start]:
            continue
        queue = deque([start])
        seen[start] = True
        block: list[int] = []
        while queue:
            x = queue.popleft()
            block.append(x)
            for operator in all_operators:
                y = operator[x]
                if not seen[y]:
                    seen[y] = True
                    queue.append(y)
        block.sort()
        block_number = len(blocks)
        for x in block:
            block_of[x] = block_number
        blocks.append(block)

    return blocks, block_of


def quotient_operator(
    operator: Sequence[int],
    blocks: Sequence[Sequence[int]],
    block_of: Sequence[int],
) -> Perm:
    images: list[int] = []
    for block in blocks:
        target_blocks = {block_of[operator[x]] for x in block}
        if len(target_blocks) != 1:
            raise RuntimeError("operator does not descend to the Q'' quotient")
        images.append(next(iter(target_blocks)))
    return tuple(images)


def class_ordering(tuple_state: Tuple4) -> tuple[int, int, int, int]:
    return tuple(base.permutation_order(g) for g in tuple_state)  # type: ignore[return-value]


def finish_certificate(
    orbit: list[Tuple4],
    index: dict[tuple[int, ...], int],
    q_operators: list[list[int]],
) -> dict[str, object]:
    q1, q2, q3 = map(tuple, q_operators)
    full_degree = len(orbit)

    braid_checks = {
        "q1_q3_commute": compose(q1, q3) == compose(q3, q1),
        "q1_q2_q1_equals_q2_q1_q2": (
            compose(compose(q1, q2), q1) == compose(compose(q2, q1), q2)
        ),
        "q2_q3_q2_equals_q3_q2_q3": (
            compose(compose(q2, q3), q2) == compose(compose(q3, q2), q3)
        ),
    }

    q_double_prime_1 = compose(q1, power(q3, -1))
    shift = compose(compose(q1, q2), q3)
    q_double_prime_2 = power(shift, 2)

    blocks, block_of = orbit_partition(
        full_degree,
        [q_double_prime_1, q_double_prime_2],
    )

    reduced_q1 = quotient_operator(q1, blocks, block_of)
    reduced_q2 = quotient_operator(q2, blocks, block_of)
    reduced_q3 = quotient_operator(q3, blocks, block_of)

    gamma0 = compose(reduced_q1, reduced_q2)
    gamma1 = compose(compose(reduced_q1, reduced_q2), reduced_q1)
    gamma_infinity = reduced_q2

    reduced_degree = len(blocks)
    identity = tuple(range(reduced_degree))
    reduced_components, _ = orbit_partition(
        reduced_degree,
        [gamma0, gamma1, gamma_infinity],
    )

    cycle_counts = [
        len(cycle_lengths(gamma0)),
        len(cycle_lengths(gamma1)),
        len(cycle_lengths(gamma_infinity)),
    ]
    indices = [reduced_degree - count for count in cycle_counts]
    genus = 1 - reduced_degree + sum(indices) // 2

    ordering_counts = Counter(class_ordering(tuple_state) for tuple_state in orbit)
    block_ordering_sets = [
        tuple(sorted(class_ordering(orbit[x]) for x in block))
        for block in blocks
    ]
    color_keys = sorted(set(block_ordering_sets))
    color_index = {key: i for i, key in enumerate(color_keys)}
    colors = [color_index[key] for key in block_ordering_sets]

    color_actions: dict[str, dict[str, int]] = {}
    for name, operator in (
        ("gamma0", gamma0),
        ("gamma1", gamma1),
        ("gamma_infinity", gamma_infinity),
    ):
        action: dict[str, int] = {}
        for color in range(len(color_keys)):
            targets = {
                colors[operator[i]]
                for i, source_color in enumerate(colors)
                if source_color == color
            }
            if len(targets) != 1:
                raise RuntimeError("reduced generator does not preserve ordering colors")
            action[str(color)] = next(iter(targets))
        color_actions[name] = action

    target_ordering = (3, 2, 2, 5)
    target_color = next(
        color for color, key in enumerate(color_keys) if target_ordering in key
    )
    slice_blocks = [i for i, color in enumerate(colors) if color == target_color]
    local_index = {block: i for i, block in enumerate(slice_blocks)}

    full_beta13 = compose(compose(q2, power(q1, 2)), power(q2, -1))
    reduced_beta13 = quotient_operator(full_beta13, blocks, block_of)
    reduced_eta23 = power(compose(reduced_beta13, reduced_q2), -1)

    def restrict_to_slice(operator: Sequence[int]) -> Perm:
        if any(colors[operator[block]] != target_color for block in slice_blocks):
            raise RuntimeError("symmetrized operator does not preserve its color slice")
        return tuple(local_index[operator[block]] for block in slice_blocks)

    slice_beta13 = restrict_to_slice(reduced_beta13)
    slice_beta3 = restrict_to_slice(reduced_q2)
    slice_eta23 = restrict_to_slice(reduced_eta23)
    slice_identity = tuple(range(len(slice_blocks)))
    symmetrized_slice = {
        "size": len(slice_blocks),
        "color": target_color,
        "ordering": list(target_ordering),
        "beta13_cycle_type": cycle_type(slice_beta13),
        "beta3_cycle_type": cycle_type(slice_beta3),
        "eta23_cycle_type": cycle_type(slice_eta23),
        "product_one": compose(
            compose(slice_beta13, slice_beta3), slice_eta23
        ) == slice_identity,
    }

    gamma_types = [
        cycle_type(gamma0),
        cycle_type(gamma1),
        cycle_type(gamma_infinity),
    ]
    gamma_orders = [
        permutation_order(gamma0),
        permutation_order(gamma1),
        permutation_order(gamma_infinity),
    ]

    checks = {
        **braid_checks,
        "full_orbit_size_11760": full_degree == 11760,
        "twelve_class_orderings": len(ordering_counts) == 12,
        "980_states_per_ordering": set(ordering_counts.values()) == {980},
        "q_double_prime_generators_are_involutions": (
            permutation_order(q_double_prime_1) == 2
            and permutation_order(q_double_prime_2) == 2
        ),
        "q_double_prime_generators_commute": (
            compose(q_double_prime_1, q_double_prime_2)
            == compose(q_double_prime_2, q_double_prime_1)
        ),
        "q_double_prime_action_free": all(len(block) == 4 for block in blocks),
        "reduced_degree_2940": reduced_degree == 2940,
        "gamma_orders": gamma_orders == [3, 2, 18480],
        "gamma_cycle_types": gamma_types == [
            {3: 980},
            {2: 1470},
            {2: 14, 3: 38, 4: 89, 5: 30, 6: 86,
             8: 41, 10: 46, 14: 38, 16: 12, 22: 12},
        ],
        "gamma_product_one": compose(compose(gamma0, gamma1), gamma_infinity) == identity,
        "reduced_action_transitive": [len(component) for component in reduced_components] == [2940],
        "reduced_genus_43": genus == 43,
        "three_ordering_colors": Counter(colors) == Counter({0: 980, 1: 980, 2: 980}),
        "symmetrized_slice_identification": (
            symmetrized_slice["size"] == 980
            and symmetrized_slice["beta13_cycle_type"]
                == {1: 2, 2: 21, 3: 16, 4: 41, 5: 46, 7: 38, 8: 12, 11: 12}
            and symmetrized_slice["beta3_cycle_type"]
                == {2: 12, 3: 38, 4: 68, 5: 30, 6: 70}
            and symmetrized_slice["eta23_cycle_type"] == {2: 490}
            and bool(symmetrized_slice["product_one"])
        ),
    }

    return {
        "certificate": "M23_FULL_REDUCED_JLINE_CERTIFICATE",
        "scope": (
            "full Hurwitz-braid orbit across all class orderings, quotient by "
            "Q''=<q1 q3^-1, sh^2>, and the resulting reduced j-line monodromy"
        ),
        "environment": {
            "python": platform.python_version(),
            "sympy": sympy.__version__,
            "platform": platform.platform(),
            "script_sha256": script_sha256(),
        },
        "full_all_orderings_orbit_size": full_degree,
        "ordering_count": len(ordering_counts),
        "states_per_ordering": {
            str(key): value for key, value in sorted(ordering_counts.items())
        },
        "q_double_prime": {
            "generators": ["q1 q3^-1", "sh^2"],
            "generator_orders": [
                permutation_order(q_double_prime_1),
                permutation_order(q_double_prime_2),
            ],
            "commute": (
                compose(q_double_prime_1, q_double_prime_2)
                == compose(q_double_prime_2, q_double_prime_1)
            ),
            "orbit_count": len(blocks),
            "orbit_sizes": dict(Counter(len(block) for block in blocks)),
            "free": all(len(block) == 4 for block in blocks),
        },
        "reduced_degree": reduced_degree,
        "gamma_names": [
            "gamma0=q1 q2",
            "gamma1=sh=q1 q2 q1",
            "gamma_infinity=q2",
        ],
        "gamma_orders": gamma_orders,
        "gamma_cycle_types": gamma_types,
        "gamma_cycle_counts": cycle_counts,
        "gamma_indices": indices,
        "gamma_product_one": compose(compose(gamma0, gamma1), gamma_infinity) == identity,
        "component_sizes": [len(component) for component in reduced_components],
        "genus": genus,
        "elliptic_fibers": {
            "normalized_j_0_points": 980,
            "normalized_j_0_ramification_index": 3,
            "classical_j_1728_points": 1470,
            "classical_j_1728_ramification_index": 2,
        },
        "cusp_count": cycle_counts[2],
        "minimum_cusp_width": min(cycle_lengths(gamma_infinity)),
        "ordering_colors": {
            "count": len(color_keys),
            "sizes": dict(Counter(colors)),
            "representatives": [
                [list(ordering) for ordering in key]
                for key in color_keys
            ],
            "actions": color_actions,
        },
        "symmetrized_slice": symmetrized_slice,
        "operator_sha256": [
            operator_sha256(gamma0),
            operator_sha256(gamma1),
            operator_sha256(gamma_infinity),
        ],
        "checks": checks,
        "pass": all(checks.values()),
    }


def print_summary(result: dict[str, object]) -> None:
    print(result["certificate"])
    print("PASS", result["pass"])
    print("FULL_ALL_ORDERINGS_ORBIT_SIZE", result["full_all_orderings_orbit_size"])
    print("REDUCED_DEGREE", result["reduced_degree"])
    print("GAMMA_ORDERS", result["gamma_orders"])
    print("GAMMA_CYCLE_TYPES", result["gamma_cycle_types"])
    print("GAMMA_CYCLE_COUNTS", result["gamma_cycle_counts"])
    print("GAMMA_INDICES", result["gamma_indices"])
    print("GAMMA_PRODUCT_ONE", result["gamma_product_one"])
    print("COMPONENT_SIZES", result["component_sizes"])
    print("REDUCED_GENUS", result["genus"])
    print("ELLIPTIC_FIBERS", result["elliptic_fibers"])
    print("CUSP_COUNT", result["cusp_count"])
    print("MINIMUM_CUSP_WIDTH", result["minimum_cusp_width"])
    print("SYMMETRIZED_SLICE", result["symmetrized_slice"])
    failed = [name for name, ok in result["checks"].items() if not ok]  # type: ignore[union-attr]
    print("FAILED_CHECKS", failed)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json",
        type=Path,
        default=ROOT / "results" / "full_reduced_jline.json",
    )
    args = parser.parse_args()

    orbit, index = build_full_orbit()
    q_operators = build_q_operators(orbit, index)
    result = finish_certificate(orbit, index, q_operators)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print_summary(result)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
