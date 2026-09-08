#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import construct_m23_reduced_jline as reduced
import verify_m23_hurwitz_certificate as base


def inverse_reversal_operator(
    orbit: list[base.Tuple4], index: dict[tuple[int, ...], int]
) -> tuple[int, ...]:
    values: list[int] = []
    for row, state in enumerate(orbit, start=1):
        transformed = tuple(base.inverse(g) for g in reversed(state))
        _, key = base.canonical_tuple_and_key(transformed)
        if key not in index:
            raise RuntimeError("inverse reversal leaves the full orbit")
        values.append(index[key])
        if row % 500 == 0:
            base.canonical_tuple_and_key.cache_clear()
            base.inverse.cache_clear()
    base.canonical_tuple_and_key.cache_clear()
    base.inverse.cache_clear()
    return tuple(values)


def fixed_count(p: tuple[int, ...]) -> int:
    return sum(i == image for i, image in enumerate(p))


def verify() -> dict[str, object]:
    orbit, index = reduced.build_full_orbit()
    q1, q2, q3 = map(tuple, reduced.build_q_operators(orbit, index))
    identity = tuple(range(len(orbit)))

    a = reduced.compose(q1, reduced.power(q3, -1))
    sh = reduced.compose(reduced.compose(q1, q2), q3)
    b = reduced.power(sh, 2)
    c = inverse_reversal_operator(orbit, index)

    ab = reduced.compose(a, b)
    ac = reduced.compose(a, c)
    bc = reduced.compose(b, c)
    abc = reduced.compose(ab, c)

    elements = {
        "1": identity,
        "a": a,
        "b": b,
        "c": c,
        "ab": ab,
        "ac": ac,
        "bc": bc,
        "abc": abc,
    }

    fixed = {name: fixed_count(p) for name, p in elements.items()}
    group_order = len(set(elements.values()))
    blocks_e, _ = reduced.orbit_partition(len(orbit), [a, b, c])
    e_orbit_sizes = dict(sorted(Counter(len(block) for block in blocks_e).items()))

    blocks_q, block_of_q = reduced.orbit_partition(len(orbit), [a, b])
    c_quotient = reduced.quotient_operator(c, blocks_q, block_of_q)
    c_quotient_type = reduced.cycle_type(c_quotient)

    character_dimensions: dict[str, int] = {}
    for epsilon in (1, -1):
        for delta in (1, -1):
            for kappa in (1, -1):
                key = f"{epsilon:+d},{delta:+d},{kappa:+d}"
                character_dimensions[key] = (11760 + 328 * epsilon * delta * kappa) // 8

    checks = {
        "full_orbit_size_11760": len(orbit) == 11760,
        "a_b_c_involutions": all(
            reduced.compose(p, p) == identity for p in (a, b, c)
        ),
        "a_b_c_commute": all(
            reduced.compose(x, y) == reduced.compose(y, x)
            for x, y in ((a, b), (a, c), (b, c))
        ),
        "group_order_8": group_order == 8,
        "inverse_reversal_braid_relations": (
            reduced.compose(reduced.compose(c, q1), c) == reduced.power(q3, -1)
            and reduced.compose(reduced.compose(c, q2), c) == reduced.power(q2, -1)
            and reduced.compose(reduced.compose(c, q3), c) == reduced.power(q1, -1)
        ),
        "orbit_sizes_82x4_1429x8": e_orbit_sizes == {4: 82, 8: 1429},
        "fixed_point_census": fixed
        == {"1": 11760, "a": 0, "b": 0, "c": 0, "ab": 0, "ac": 0, "bc": 0, "abc": 328},
        "character_dimensions_1511_1429": Counter(character_dimensions.values())
        == Counter({1511: 4, 1429: 4}),
        "q_double_prime_2940_blocks": len(blocks_q) == 2940
        and Counter(len(block) for block in blocks_q) == Counter({4: 2940}),
        "reduced_c_cycle_type_1_82_2_1429": c_quotient_type == {1: 82, 2: 1429},
    }

    return {
        "certificate": "M23_RANK3_HARMONIC_EXTENSION",
        "pass": all(checks.values()),
        "checks": checks,
        "full_orbit_size": len(orbit),
        "group_order": group_order,
        "group_orbit_sizes": {str(k): v for k, v in e_orbit_sizes.items()},
        "fixed_point_counts": fixed,
        "character_dimensions": character_dimensions,
        "character_dimension_multiplicities": {
            str(k): v for k, v in sorted(Counter(character_dimensions.values()).items())
        },
        "q_double_prime_block_count": len(blocks_q),
        "inverse_reversal_on_reduced_blocks_cycle_type": {
            str(k): v for k, v in c_quotient_type.items()
        },
        "braid_intertwining": {
            "c_q1_c": "q3^-1",
            "c_q2_c": "q2^-1",
            "c_q3_c": "q1^-1",
        },
        "conclusion": (
            "Inverse reversal extends the free Klein-four reduction to an exact "
            "(C2)^3 action. The eight character spaces have dimensions 1511 "
            "and 1429 by parity, and the induced involution on the 2940 "
            "Klein blocks has cycle type 1^82 2^1429."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    result = verify()
    if args.json:
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("M23_RANK3_HARMONIC_EXTENSION")
    print("PASS", result["pass"])
    print("FULL_ORBIT_SIZE", result["full_orbit_size"])
    print("GROUP_ORDER", result["group_order"])
    print("GROUP_ORBIT_SIZES", result["group_orbit_sizes"])
    print("FIXED_POINT_COUNTS", result["fixed_point_counts"])
    print("CHARACTER_DIMENSION_MULTIPLICITIES", result["character_dimension_multiplicities"])
    print("REDUCED_CYCLE_TYPE", result["inverse_reversal_on_reduced_blocks_cycle_type"])
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
