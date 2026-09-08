#!/usr/bin/env python3
"""Exact trace-closure verification for the M23 article.

This certificate computes every trace used in the article:
* inverse reversal and its products with the reduced j-line generators;
* setwise action on the elliptic and cusp cycles;
* the split quadratic X^2-X-1 over F_211 and its 127-coordinate change;
* the Frobenius trace of the released stage-five obstruction;
* the C4 action on the 3x3 reference box.

The script distinguishes permutation trace, field trace, and a central
idempotent. It does not infer a rational point from a Reynolds projector.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import construct_m23_reduced_jline as reduced
import verify_m23_hurwitz_certificate as base
import verify_m23_rank3_harmonic_extension as rank3

P = 211
C5 = (202,108,168,30,147,178,91,84,32,31,129,66)
M23_ELEMENT_ORDERS = (1,2,3,4,5,6,7,8,11,14,15,23)
ADDRESS_SET_C = (2,3,5,7,13,17,19,31,61,127)


def fixed_count(p: tuple[int, ...]) -> int:
    return sum(i == image for i, image in enumerate(p))


def cycles(p: tuple[int, ...]) -> list[tuple[int, ...]]:
    seen = [False] * len(p)
    out: list[tuple[int, ...]] = []
    for start in range(len(p)):
        if seen[start]:
            continue
        cycle: list[int] = []
        x = start
        while not seen[x]:
            seen[x] = True
            cycle.append(x)
            x = p[x]
        out.append(tuple(cycle))
    return out


def cycle_action_data(p: tuple[int, ...], c: tuple[int, ...]) -> dict[str, object]:
    cs = cycles(p)
    cycle_of = {x: j for j, cyc in enumerate(cs) for x in cyc}
    action: list[int] = []
    fixed_sizes: Counter[int] = Counter()
    pointwise = 0
    for j, cyc in enumerate(cs):
        targets = {cycle_of[c[x]] for x in cyc}
        if len(targets) != 1:
            raise AssertionError("involution does not act on cycle set")
        target = next(iter(targets))
        action.append(target)
        if target == j:
            fixed_sizes[len(cyc)] += 1
            if all(c[x] == x for x in cyc):
                pointwise += 1
    return {
        "cycle_count": len(cs),
        "action_cycle_type": reduced.cycle_type(tuple(action)),
        "setwise_fixed_cycles": sum(fixed_sizes.values()),
        "pointwise_fixed_cycles": pointwise,
        "fixed_cycle_sizes": {str(k): v for k, v in sorted(fixed_sizes.items())},
    }


def main_result() -> dict[str, object]:
    orbit, index = reduced.build_full_orbit()
    q1, q2, q3 = map(tuple, reduced.build_q_operators(orbit, index))
    a = reduced.compose(q1, reduced.power(q3, -1))
    shift = reduced.compose(reduced.compose(q1, q2), q3)
    b = reduced.power(shift, 2)
    blocks, block_of = reduced.orbit_partition(len(orbit), [a, b])

    rq1 = reduced.quotient_operator(q1, blocks, block_of)
    rq2 = reduced.quotient_operator(q2, blocks, block_of)
    gamma0 = reduced.compose(rq1, rq2)
    gamma1 = reduced.compose(reduced.compose(rq1, rq2), rq1)
    gamma_inf = rq2

    c_full = rank3.inverse_reversal_operator(orbit, index)
    c = reduced.quotient_operator(c_full, blocks, block_of)
    c_gamma1 = reduced.compose(c, gamma1)
    c_gamma_inf = reduced.compose(c, gamma_inf)

    elliptic = cycle_action_data(gamma1, c)
    cusp = cycle_action_data(gamma_inf, c)

    # Split finite-field opening.
    roots = (33, 179)
    q_values = tuple((x*x - x - 1) % P for x in roots)
    inv127 = pow(127, -1, P)
    inverse_roots = tuple((inv127*x) % P for x in roots)
    q127_values = tuple((n*n - inv127*n - inv127*inv127) % P for n in inverse_roots)
    disc = 5 % P
    disc127 = (disc * inv127 * inv127) % P

    # Frobenius and field trace of the released obstruction.
    frobenius_c5 = tuple(pow(x, P, P) for x in C5)
    trace_c5 = tuple((x+y) % P for x,y in zip(C5, frobenius_c5))
    c5_times_127 = tuple((127*x) % P for x in C5)
    c5_times_108 = tuple((inv127*x) % P for x in C5)

    # C4 reference box.
    cells = tuple((r,c0) for r in (-1,0,1) for c0 in (-1,0,1))
    cell_index = {x:i for i,x in enumerate(cells)}
    rot = tuple(cell_index[(-c0,r)] for r,c0 in cells)
    rot2 = reduced.compose(rot, rot)
    rot3 = reduced.compose(rot2, rot)
    identity9 = tuple(range(9))
    box_orbits, _ = reduced.orbit_partition(9, [rot])
    box_projector_rank = len(box_orbits)

    checks = {
        "reduced_degree_2940": len(blocks) == 2940,
        "c_cycle_type": reduced.cycle_type(c) == {1:82, 2:1429},
        "c_commutes_gamma1": reduced.compose(c,gamma1) == reduced.compose(gamma1,c),
        "c_conjugates_gamma_inf_to_inverse": (
            reduced.compose(reduced.compose(c,gamma_inf),c) == reduced.power(gamma_inf,-1)
        ),
        "trace_c_82": fixed_count(c) == 82,
        "trace_c_gamma1_6": fixed_count(c_gamma1) == 6,
        "trace_c_gamma_inf_6": fixed_count(c_gamma_inf) == 6,
        "elliptic_cycle_action": elliptic["action_cycle_type"] == {1:44,2:713},
        "elliptic_fixed_cycle_split": elliptic["setwise_fixed_cycles"] == 44 and elliptic["pointwise_fixed_cycles"] == 41,
        "cusp_cycle_action": cusp["action_cycle_type"] == {1:44,2:181},
        "cusp_fixed_sizes": cusp["fixed_cycle_sizes"] == {"3":6,"4":15,"6":8,"8":7,"10":8},
        "finite_field_roots": q_values == (0,0),
        "finite_field_discriminant_square": pow(65,2,P) == disc,
        "inverse_127": inv127 == 108 and (127*108) % P == 1,
        "inverse_chart_roots": inverse_roots == (188,131) and q127_values == (0,0),
        "inverse_chart_discriminant": disc127 == 84 and pow(57,2,P) == 84,
        "c5_frobenius_fixed": frobenius_c5 == C5,
        "c5_trace_nonzero": trace_c5 == (193,5,125,60,83,145,182,168,64,62,47,132),
        "c5_scalings_nonzero": any(c5_times_127) and any(c5_times_108),
        "box_rotation_order_four": reduced.compose(rot3,rot) == identity9,
        "box_unique_fixed_cell": fixed_count(rot) == 1 and cells[next(i for i,x in enumerate(rot) if i==x)] == (0,0),
        "box_orbits_1_4_4": Counter(len(o) for o in box_orbits) == Counter({4:2,1:1}),
        "box_reynolds_rank_three": box_projector_rank == 3,
        "order_spectrum_and_address_set_distinct": set(M23_ELEMENT_ORDERS) != set(ADDRESS_SET_C),
        "address_only_values_not_group_orders": set(ADDRESS_SET_C)-set(M23_ELEMENT_ORDERS) == {13,17,19,31,61,127},
    }

    return {
        "certificate": "M23_TRACE_CLOSURE",
        "pass": all(checks.values()),
        "checks": checks,
        "reduced_trace_table": {
            "tr(c)": fixed_count(c),
            "tr(c gamma_1)": fixed_count(c_gamma1),
            "tr(c gamma_infinity)": fixed_count(c_gamma_inf),
            "cycle_type(c gamma_1)": {str(k):v for k,v in reduced.cycle_type(c_gamma1).items()},
            "cycle_type(c gamma_infinity)": {str(k):v for k,v in reduced.cycle_type(c_gamma_inf).items()},
        },
        "elliptic_cycle_action": elliptic,
        "cusp_cycle_action": cusp,
        "finite_field_quadratic": {
            "polynomial": "X^2-X-1",
            "roots_mod_211": list(roots),
            "discriminant": disc,
            "sqrt_discriminant": 65,
            "inverse_127": inv127,
            "inverse_chart_roots": list(inverse_roots),
            "inverse_chart_discriminant": disc127,
        },
        "selected_obstruction": {
            "c5": list(C5),
            "frobenius_c5": list(frobenius_c5),
            "field_trace_c5": list(trace_c5),
            "127_c5": list(c5_times_127),
            "108_c5": list(c5_times_108),
        },
        "reference_box": {
            "rotation_cycle_type": {str(k):v for k,v in reduced.cycle_type(rot).items()},
            "orbit_sizes": sorted(len(o) for o in box_orbits),
            "reynolds_projector_rank": box_projector_rank,
            "unique_fixed_cell": [0,0],
        },
        "sets": {
            "M23_element_orders_natural_action": list(M23_ELEMENT_ORDERS),
            "external_address_set_C": list(ADDRESS_SET_C),
        },
        "conclusion": (
            "All supplied trace combinations are computed. The Reynolds/field-trace "
            "projector has an invariant line, but the selected F_211 obstruction is "
            "Frobenius-fixed and has nonzero field trace. The 127 chart change preserves "
            "that nonvanishing. A rank-one invariant linear projector is therefore kept "
            "distinct from a rank-one central idempotent of an interior fiber algebra."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = main_result()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(result["certificate"])
    print("PASS", result["pass"])
    print("REDUCED_TRACE_TABLE", result["reduced_trace_table"])
    print("FIELD_TRACE_C5", result["selected_obstruction"]["field_trace_c5"])
    print("REFERENCE_BOX", result["reference_box"])
    return 0 if result["pass"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
