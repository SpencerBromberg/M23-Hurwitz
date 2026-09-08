#!/usr/bin/env python3
"""Exact verification of a strict power-intertwiner model on the 980-point orbit.

This program proves a finite permutation-graph statement. Let s_12,s_13,s_14
be the three exported peripheral permutations. For every unit m modulo the
least common multiple of their orders, it tests whether a bijection Phi exists
with

    Phi s_ij = s_ij^m Phi

for all three colors, while interchanging either of the two actual fixed-point
pairs of s_13 or s_14. The common modulus is 9240, and all 1920 units are
checked in both directions for both pairs.

Scope: this is a complete certificate for the displayed strict semilinear
model. It also computes the simultaneous centralizer of the three colored
operators. Arithmetic Galois transport requires the separate
normalizer--power-conjugacy condition stated in the paper.
"""

from __future__ import annotations

from collections import Counter, deque
from pathlib import Path
from typing import Sequence
import argparse
import gzip
import hashlib
import json
import math
import platform
import sys
import time

Perm = tuple[int, ...]

EXPECTED_TYPES = [
    {1: 24, 2: 136, 3: 178, 5: 30},
    {1: 2, 2: 21, 3: 16, 4: 41, 5: 46, 7: 38, 8: 12, 11: 12},
    {1: 2, 2: 21, 3: 16, 4: 41, 5: 46, 7: 38, 8: 12, 11: 12},
]
EXPECTED_ORDERS = [30, 9240, 9240]
EXPECTED_FIXED = [
    [1, 96, 118, 162, 236, 249, 299, 381, 398, 506, 508, 535,
     542, 554, 562, 579, 581, 600, 602, 625, 785, 788, 804, 956],
    [454, 914],
    [634, 974],
]


def load_dump(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def script_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def cycle_lengths(p: Sequence[int]) -> list[int]:
    seen = [False] * len(p)
    lengths: list[int] = []
    for start in range(len(p)):
        if seen[start]:
            continue
        x = start
        length = 0
        while not seen[x]:
            seen[x] = True
            length += 1
            x = p[x]
        lengths.append(length)
    return sorted(lengths)


def cycle_type(p: Sequence[int]) -> dict[int, int]:
    return dict(sorted(Counter(cycle_lengths(p)).items()))


def permutation_order(p: Sequence[int]) -> int:
    return math.lcm(*cycle_lengths(p))


def fixed_points_one_based(p: Sequence[int]) -> list[int]:
    return [i + 1 for i, image in enumerate(p) if image == i]


def compose(a: Sequence[int], b: Sequence[int]) -> Perm:
    return tuple(b[a[i]] for i in range(len(a)))


def power(p: Sequence[int], exponent: int) -> Perm:
    n = len(p)
    result: Perm = tuple(range(n))
    base: Perm = tuple(p)
    e = exponent
    while e:
        if e & 1:
            result = compose(result, base)
        base = compose(base, base)
        e >>= 1
    return result


def strict_intertwiner_exists(
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
    source: int,
    target: int,
) -> bool:
    """Propagate the unique rooted color-preserving candidate map."""
    n = len(left[0])
    image = [-1] * n
    preimage = [-1] * n
    image[source] = target
    preimage[target] = source
    queue: deque[int] = deque([source])

    while queue:
        x = queue.popleft()
        y = image[x]
        for source_op, target_op in zip(left, right):
            next_x = source_op[x]
            next_y = target_op[y]
            current = image[next_x]
            if current == -1:
                if preimage[next_y] != -1:
                    return False
                image[next_x] = next_y
                preimage[next_y] = next_x
                queue.append(next_x)
            elif current != next_y:
                return False

    return all(value >= 0 for value in image)


def verify(input_path: Path) -> dict[str, object]:
    started = time.perf_counter()
    data = load_dump(input_path)
    operators = [tuple(op) for op in data["forward_operators"]]  # type: ignore[index]
    metadata = data["metadata"]  # type: ignore[index]
    n = int(metadata["orbit_size"])  # type: ignore[index]
    names = list(metadata["operator_names"])  # type: ignore[index]

    types = [cycle_type(op) for op in operators]
    orders = [permutation_order(op) for op in operators]
    fixed = [fixed_points_one_based(op) for op in operators]
    modulus = math.lcm(*orders)
    units = [m for m in range(1, modulus) if math.gcd(m, modulus) == 1]

    # The two natural two-point fixed divisors in the exported peripheral data.
    pairs = {
        "B13_fixed_pair": [point - 1 for point in fixed[1]],
        "B14_fixed_pair": [point - 1 for point in fixed[2]],
    }

    # The simultaneous centralizer is the color-preserving automorphism group
    # of the connected three-colored directed graph. A root image determines
    # at most one automorphism, so all n target roots give an exact count.
    centralizer_root_images = [
        target + 1
        for target in range(n)
        if strict_intertwiner_exists(operators, operators, 0, target)
    ]

    x_y_commute = compose(operators[0], operators[1]) == compose(operators[1], operators[0])

    successes: list[dict[str, object]] = []
    per_pair_tests = {name: 0 for name in pairs}
    for m in units:
        twisted = [power(op, m) for op in operators]
        for pair_name, pair in pairs.items():
            if len(pair) != 2:
                raise RuntimeError(f"{pair_name} is not a two-point pair")
            a, b = pair
            for direction, source, target in (
                ("forward", a, b),
                ("reverse", b, a),
            ):
                per_pair_tests[pair_name] += 1
                if strict_intertwiner_exists(operators, twisted, source, target):
                    successes.append(
                        {
                            "pair": pair_name,
                            "m": m,
                            "direction": direction,
                            "source_one_based": source + 1,
                            "target_one_based": target + 1,
                        }
                    )

    checks = {
        "orbit_size": n == 980,
        "operator_names": names == ["B12", "B13", "B14"],
        "cycle_types": types == EXPECTED_TYPES,
        "operator_orders": orders == EXPECTED_ORDERS,
        "fixed_points": fixed == EXPECTED_FIXED,
        "common_modulus": modulus == 9240,
        "unit_count": len(units) == 1920,
        "tests_per_pair": all(value == 3840 for value in per_pair_tests.values()),
        "total_tests": sum(per_pair_tests.values()) == 7680,
        "simultaneous_centralizer_order": len(centralizer_root_images) == 1,
        "simultaneous_centralizer_root_images": centralizer_root_images == [1],
        "first_two_generators_noncommuting": not x_y_commute,
        "no_strict_swap": len(successes) == 0,
    }
    passed = all(checks.values())

    return {
        "certificate": "M23_STRICT_POWER_MODEL_CERTIFICATE",
        "scope": (
            "finite rooted colored-digraph theorem for the strict relations "
            "Phi*s_ij=s_ij^m*Phi, with exact simultaneous centralizer; "
            "arithmetic transport uses the separate normalizer-power-conjugacy criterion"
        ),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "script_sha256": script_sha256(),
            "input_sha256": file_sha256(input_path),
        },
        "input": str(input_path),
        "orbit_size": n,
        "operator_names": names,
        "cycle_types": types,
        "operator_orders": orders,
        "fixed_points_one_based": fixed,
        "common_power_modulus": modulus,
        "unit_count": len(units),
        "audited_pairs_one_based": {
            name: [point + 1 for point in pair] for name, pair in pairs.items()
        },
        "tests_per_pair": per_pair_tests,
        "total_tests": sum(per_pair_tests.values()),
        "simultaneous_centralizer_order": len(centralizer_root_images),
        "simultaneous_centralizer_root_images_one_based": centralizer_root_images,
        "first_two_generators_commute": x_y_commute,
        "successes": successes,
        "checks": checks,
        "pass": passed,
        "elapsed_seconds": time.perf_counter() - started,
    }


def print_summary(result: dict[str, object]) -> None:
    print(result["certificate"])
    print("PASS", result["pass"])
    print("SCOPE", result["scope"])
    print("ORBIT_SIZE", result["orbit_size"])
    print("OPERATOR_ORDERS", result["operator_orders"])
    print("FIXED_POINTS_ONE_BASED", result["fixed_points_one_based"])
    print("COMMON_POWER_MODULUS", result["common_power_modulus"])
    print("UNIT_COUNT", result["unit_count"])
    print("AUDITED_PAIRS_ONE_BASED", result["audited_pairs_one_based"])
    print("TESTS_PER_PAIR", result["tests_per_pair"])
    print("TOTAL_TESTS", result["total_tests"])
    print("SIMULTANEOUS_CENTRALIZER_ORDER", result["simultaneous_centralizer_order"])
    print("SIMULTANEOUS_CENTRALIZER_ROOT_IMAGES_ONE_BASED", result["simultaneous_centralizer_root_images_one_based"])
    print("FIRST_TWO_GENERATORS_COMMUTE", result["first_two_generators_commute"])
    print("TOTAL_SUCCESSES", len(result["successes"]))  # type: ignore[arg-type]
    failed = [name for name, ok in result["checks"].items() if not ok]  # type: ignore[union-attr]
    print("FAILED_CHECKS", failed)


def main() -> int:
    default_input = Path(__file__).resolve().parents[1] / "data" / "orbit_certificate.json.gz"
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    result = verify(args.input)
    print_summary(result)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
