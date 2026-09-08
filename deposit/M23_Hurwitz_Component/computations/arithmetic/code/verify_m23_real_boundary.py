#!/usr/bin/env python3
"""Exact real-structure and stable-signature certificate for the M23 Hurwitz component.

The script reconstructs the ordered 980-state Nielsen orbit from the degree-23
seed tuple.  It then:

* constructs the unique fiber permutation realizing complex conjugation for a
  real base point lambda_0 < 0 and the standard symmetric peripheral paths;
* verifies the exact conjugation relations with B12, B13, and B14;
* determines the arithmetic monodromy group from the parity of this operator;
* computes the induced involutions on all three compactified boundary fibers;
* counts real boundary points by cusp width;
* attaches an intrinsic stable-component signature to every boundary point;
* identifies every point whose signature is unique in its full geometric
  boundary fiber, thereby certifying individual Q-rationality.

All calculations use exact permutations and finite group algorithms.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Iterable, Sequence

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[1]
CODE = SCRIPT.parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))

import sympy
from sympy.combinatorics import PermutationGroup
import verify_m23_hurwitz_certificate as base
import verify_m23_order11_boundary as boundary

Perm0 = tuple[int, ...]


def script_sha256() -> str:
    return hashlib.sha256(SCRIPT.read_bytes()).hexdigest()


def inverse_operator(p: Sequence[int]) -> Perm0:
    out = [0] * len(p)
    for i, image in enumerate(p):
        out[image] = i
    return tuple(out)


def compose_operator(a: Sequence[int], b: Sequence[int]) -> Perm0:
    """Right-action product: x^(a b)=(x^a)^b."""
    return tuple(b[a[i]] for i in range(len(a)))


def operator_cycles(p: Sequence[int]) -> list[list[int]]:
    seen = [False] * len(p)
    result: list[list[int]] = []
    for start in range(len(p)):
        if seen[start]:
            continue
        cycle: list[int] = []
        x = start
        while not seen[x]:
            seen[x] = True
            cycle.append(x)
            x = p[x]
        result.append(cycle)
    return result


def permutation_parity(p: Sequence[int]) -> str:
    cycles = len(operator_cycles(p))
    return "even" if (len(p) - cycles) % 2 == 0 else "odd"


def solve_colored_intertwiner(
    source_generators: Sequence[Sequence[int]],
    target_generators: Sequence[Sequence[int]],
    root: int = 0,
) -> list[Perm0]:
    """Enumerate colored graph isomorphisms from source to target generators.

    A solution c satisfies c(s_i(x))=t_i(c(x)) for every color i and point x.
    The source action is transitive, so a root image determines at most one
    solution.
    """
    n = len(source_generators[0])
    source = list(source_generators) + [inverse_operator(p) for p in source_generators]
    target = list(target_generators) + [inverse_operator(p) for p in target_generators]
    solutions: list[Perm0] = []

    for root_image in range(n):
        image = [-1] * n
        preimage = [-1] * n
        image[root] = root_image
        preimage[root_image] = root
        queue = deque([root])
        valid = True

        while queue and valid:
            x = queue.popleft()
            y = image[x]
            for s, t in zip(source, target):
                sx = s[x]
                ty = t[y]
                if image[sx] == -1:
                    if preimage[ty] != -1:
                        valid = False
                        break
                    image[sx] = ty
                    preimage[ty] = sx
                    queue.append(sx)
                elif image[sx] != ty:
                    valid = False
                    break

        if valid and all(value >= 0 for value in image):
            candidate = tuple(image)
            if all(
                candidate[s[x]] == t[candidate[x]]
                for s, t in zip(source, target)
                for x in range(n)
            ):
                solutions.append(candidate)
    return solutions


def group_invariants(generators: Iterable[base.Perm]) -> dict[str, object]:
    group = PermutationGroup([base.to_sympy(g) for g in generators])
    return {
        "order": int(group.order()),
        "orbit_sizes": sorted(len(orbit) for orbit in group.orbits()),
    }


def stable_signature(
    state: base.Tuple4,
    width: int,
) -> dict[str, object]:
    a, b, c, d = state
    h = base.compose(a, b)
    h_inverse = base.inverse(h)
    left = group_invariants((a, b, h_inverse))
    right = group_invariants((h, c, d))
    return {
        "width": width,
        "node_order": base.permutation_order(h),
        "left_component": left,
        "right_component": right,
    }


def signature_key(signature: dict[str, object]) -> tuple[object, ...]:
    left = signature["left_component"]
    right = signature["right_component"]
    assert isinstance(left, dict) and isinstance(right, dict)
    return (
        int(signature["width"]),
        int(signature["node_order"]),
        int(left["order"]),
        tuple(int(x) for x in left["orbit_sizes"]),
        int(right["order"]),
        tuple(int(x) for x in right["orbit_sizes"]),
    )


def boundary_action(
    operator: Sequence[int],
    point_map: Sequence[int],
) -> tuple[list[list[int]], list[int]]:
    cycles = operator_cycles(operator)
    cycle_index = {point: i for i, cycle in enumerate(cycles) for point in cycle}
    induced: list[int] = []
    for cycle in cycles:
        targets = {cycle_index[point_map[x]] for x in cycle}
        if len(targets) != 1:
            raise RuntimeError("point map does not descend to the boundary cycle quotient")
        induced.append(next(iter(targets)))
    return cycles, induced


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--json",
        type=Path,
        default=ROOT / "results" / "M23_real_boundary_certificate.json",
    )
    args = parser.parse_args()

    orbit, index = boundary.build_raw_ordered_orbit(base.TUPLE)
    operators = [
        boundary.build_operator(orbit, index, word)
        for word in (base.B12, base.B13, base.B14)
    ]
    b12, b13, b14 = operators
    n = len(orbit)
    identity = tuple(range(n))
    b12_inverse, b13_inverse, b14_inverse = map(inverse_operator, operators)

    # With a real base point lambda_0<0 and symmetric peripheral paths:
    #   kappa(gamma_0)=gamma_0^{-1},
    #   kappa(gamma_1)=gamma_0 gamma_1^{-1} gamma_0^{-1},
    #   kappa(gamma_infinity)=gamma_infinity^{-1}.
    target_0 = b12_inverse
    target_1 = compose_operator(
        compose_operator(b12, b13_inverse),
        b12_inverse,
    )
    target_infinity = b14_inverse

    solutions = solve_colored_intertwiner(
        (b12, b13, b14),
        (target_0, target_1, target_infinity),
    )
    if len(solutions) != 1:
        raise RuntimeError(f"expected one complex-conjugation operator, found {len(solutions)}")
    conjugation = solutions[0]

    # The local path correction at lambda=1 converts the conjugated inertia
    # subgroup back to the chosen B13 generator.
    local_maps = {
        "lambda_0": conjugation,
        "lambda_1": compose_operator(conjugation, b12),
        "lambda_infinity": conjugation,
    }
    prefixes = {
        "lambda_0": (),
        "lambda_1": base.Q2,
        "lambda_infinity": base.Q3 + base.Q2,
    }
    operator_by_fiber = {
        "lambda_0": b12,
        "lambda_1": b13,
        "lambda_infinity": b14,
    }

    boundary_results: dict[str, object] = {}
    total_boundary_points = 0
    total_real_boundary_points = 0
    rational_points: list[dict[str, object]] = []

    for fiber in ("lambda_0", "lambda_1", "lambda_infinity"):
        operator = operator_by_fiber[fiber]
        point_map = local_maps[fiber]
        cycles, induced = boundary_action(operator, point_map)
        total_boundary_points += len(cycles)
        fixed_cycle_indices = [i for i, image in enumerate(induced) if image == i]
        total_real_boundary_points += len(fixed_cycle_indices)

        signatures: list[dict[str, object]] = []
        for cycle in cycles:
            representative = min(cycle)
            state = orbit[representative]
            prefix = prefixes[fiber]
            if prefix:
                state = base.apply_word(state, prefix)
            signature = stable_signature(state, len(cycle))
            signatures.append(signature)

        multiplicities = Counter(signature_key(signature) for signature in signatures)
        fixed_by_width = Counter(len(cycles[i]) for i in fixed_cycle_indices)
        fixed_by_node_order = Counter(
            int(signatures[i]["node_order"]) for i in fixed_cycle_indices
        )

        singleton_rows: list[dict[str, object]] = []
        for i, (cycle, signature) in enumerate(zip(cycles, signatures)):
            if multiplicities[signature_key(signature)] != 1:
                continue
            row = {
                "fiber": fiber,
                "representative_one_based": min(cycle) + 1,
                "cycle_one_based": [x + 1 for x in cycle],
                "real": induced[i] == i,
                "signature": signature,
            }
            singleton_rows.append(row)
            rational_points.append(row)

        boundary_results[fiber] = {
            "geometric_point_count": len(cycles),
            "real_point_count": len(fixed_cycle_indices),
            "nonreal_point_count": len(cycles) - len(fixed_cycle_indices),
            "real_points_by_width": dict(sorted(fixed_by_width.items())),
            "real_points_by_node_order": dict(sorted(fixed_by_node_order.items())),
            "fixed_cycle_representatives_one_based": [
                min(cycles[i]) + 1 for i in fixed_cycle_indices
            ],
            "singleton_signature_count": len(singleton_rows),
            "singleton_signatures": singleton_rows,
        }

    conjugation_cycles = Counter(len(cycle) for cycle in operator_cycles(conjugation))
    local_1 = local_maps["lambda_1"]
    local_1_cycles = Counter(len(cycle) for cycle in operator_cycles(local_1))
    interval_maps = {
        "(-infinity,0)": conjugation,
        "(0,1)": compose_operator(conjugation, b12),
        "(1,infinity)": compose_operator(b14, conjugation),
    }
    interval_fixed_counts = {
        name: sum(1 for i, image in enumerate(op) if image == i)
        for name, op in interval_maps.items()
    }

    checks = {
        "ordered_orbit_size_980": n == 980,
        "peripheral_product_one": compose_operator(compose_operator(b12, b13), b14) == identity,
        "unique_complex_conjugation_operator": len(solutions) == 1,
        "complex_conjugation_is_involution": compose_operator(conjugation, conjugation) == identity,
        "complex_conjugation_relation_B12": all(
            conjugation[b12[x]] == b12_inverse[conjugation[x]] for x in range(n)
        ),
        "complex_conjugation_relation_B13": all(
            conjugation[b13[x]] == target_1[conjugation[x]] for x in range(n)
        ),
        "complex_conjugation_relation_B14": all(
            conjugation[b14[x]] == b14_inverse[conjugation[x]] for x in range(n)
        ),
        "lambda1_local_map_inverts_B13": all(
            local_1[b13[x]] == b13_inverse[local_1[x]] for x in range(n)
        ),
        "complex_conjugation_cycle_type_1_26_2_477": conjugation_cycles == Counter({1: 26, 2: 477}),
        "complex_conjugation_odd": permutation_parity(conjugation) == "odd",
        "lambda1_local_map_cycle_type_1_26_2_477": local_1_cycles == Counter({1: 26, 2: 477}),
        "real_interval_sheet_counts_26_26_30": list(interval_fixed_counts.values()) == [26, 26, 30],
        "boundary_counts_368_188_188": [
            boundary_results[f]["geometric_point_count"]
            for f in ("lambda_0", "lambda_1", "lambda_infinity")
        ] == [368, 188, 188],
        "real_boundary_counts_26_28_28": [
            boundary_results[f]["real_point_count"]
            for f in ("lambda_0", "lambda_1", "lambda_infinity")
        ] == [26, 28, 28],
        "total_boundary_points_744": total_boundary_points == 744,
        "total_real_boundary_points_82": total_real_boundary_points == 82,
        "total_nonreal_boundary_points_662": total_boundary_points - total_real_boundary_points == 662,
        "width_one_pairs_exchanged": (
            boundary_results["lambda_1"]["real_points_by_width"].get(1, 0) == 0
            and boundary_results["lambda_infinity"]["real_points_by_width"].get(1, 0) == 0
        ),
        "eight_singleton_stable_signatures": len(rational_points) == 8,
        "all_singleton_signatures_real": all(bool(row["real"]) for row in rational_points),
    }

    result = {
        "certificate": "M23_REAL_STRUCTURE_AND_BOUNDARY_SIGNATURE_CERTIFICATE",
        "scope": (
            "complex conjugation on the ordered 980-sheet cover, arithmetic monodromy, "
            "real boundary census, and Q-rationality from singleton stable signatures"
        ),
        "environment": {
            "python": platform.python_version(),
            "sympy": sympy.__version__,
            "script_sha256": script_sha256(),
        },
        "ordered_orbit_size": n,
        "complex_conjugation": {
            "root_image_one_based": conjugation[0] + 1,
            "cycle_type": dict(sorted(conjugation_cycles.items())),
            "fixed_nielsen_sheets": conjugation_cycles[1],
            "parity": permutation_parity(conjugation),
            "arithmetic_monodromy": "S_980",
        },
        "real_interval_sheet_counts": interval_fixed_counts,
        "boundary": boundary_results,
        "boundary_totals": {
            "geometric_points": total_boundary_points,
            "real_points": total_real_boundary_points,
            "nonreal_points": total_boundary_points - total_real_boundary_points,
            "certified_Q_rational_points": len(rational_points),
        },
        "certified_Q_rational_boundary_points": rational_points,
        "width_one_arithmetic": {
            "lambda_1": "one imaginary quadratic closed point",
            "lambda_infinity": "one imaginary quadratic closed point",
            "field_discriminants": "undetermined by the present certificate",
        },
        "checks": checks,
        "pass": all(checks.values()),
    }

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(result["certificate"])
    print("PASS", result["pass"])
    print("COMPLEX_CONJUGATION_CYCLE_TYPE", result["complex_conjugation"]["cycle_type"])
    print("COMPLEX_CONJUGATION_PARITY", result["complex_conjugation"]["parity"])
    print("ARITHMETIC_MONODROMY", result["complex_conjugation"]["arithmetic_monodromy"])
    print("REAL_INTERVAL_SHEET_COUNTS", interval_fixed_counts)
    print("REAL_BOUNDARY_COUNTS", [boundary_results[f]["real_point_count"] for f in ("lambda_0", "lambda_1", "lambda_infinity")])
    print("BOUNDARY_TOTALS", result["boundary_totals"])
    print("CERTIFIED_Q_RATIONAL_BOUNDARY_POINTS", len(rational_points))
    print("FAILED_CHECKS", [name for name, ok in checks.items() if not ok])
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
