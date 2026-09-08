#!/usr/bin/env python3
"""Exact mixed-boundary stabilizer pullback and dual-involution certificate."""
from __future__ import annotations

import argparse
import json
from collections import deque, Counter
from pathlib import Path

from sympy.combinatorics import Permutation, PermutationGroup

HERE = Path(__file__).resolve().parent


def parse_vector(lines: list[str], key: str) -> list[int]:
    prefix = key + " "
    for line in lines:
        if line.startswith(prefix):
            return [int(x) for x in line[len(prefix):].split()]
    raise RuntimeError(f"missing {key}")


def compose(p: list[int], q: list[int]) -> list[int]:
    return [p[q[i]] for i in range(len(p))]


def inverse(p: list[int]) -> list[int]:
    out = [0] * len(p)
    for i, j in enumerate(p):
        out[j] = i
    return out


def subset_image(s: tuple[int, ...], p: list[int]) -> tuple[int, ...]:
    return tuple(sorted(p[i] for i in s))


def cycle_hist(p: list[int]) -> dict[int, int]:
    seen = [False] * len(p)
    h: Counter[int] = Counter()
    for i in range(len(p)):
        if seen[i]:
            continue
        j = i
        n = 0
        while not seen[j]:
            seen[j] = True
            n += 1
            j = p[j]
        h[n] += 1
    return dict(sorted(h.items()))


def cycle_count(p: list[int]) -> int:
    return sum(cycle_hist(p).values())


def genus_from_three(perms: list[list[int]]) -> int:
    d = len(perms[0])
    idx = sum(d - cycle_count(p) for p in perms)
    num = 2 - 2 * d + idx
    assert num % 2 == 0
    return num // 2


def subset_orbit(root: tuple[int, ...], generators: list[list[int]]):
    invs = [inverse(g) for g in generators]
    steps = generators + invs
    identity = list(range(len(generators[0])))
    reps: dict[tuple[int, ...], list[int]] = {root: identity}
    order = [root]
    q = deque([root])
    schreier: list[list[int]] = []
    while q:
        s = q.popleft()
        r = reps[s]
        for g in steps:
            s2 = subset_image(s, g)
            gr = compose(g, r)
            if s2 not in reps:
                reps[s2] = gr
                order.append(s2)
                q.append(s2)
            else:
                rp = reps[s2]
                h = compose(inverse(rp), gr)
                assert subset_image(root, h) == root
                schreier.append(h)
    return order, reps, schreier


def action_on_subsets(subsets: list[tuple[int, ...]], g: list[int]) -> list[int]:
    pos = {s: i for i, s in enumerate(subsets)}
    return [pos[subset_image(s, g)] for s in subsets]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clutching-log", default=str(HERE / "geometry_certificates" / "m23_five_branch_clutching.out"))
    ap.add_argument("--projector-json", default=str(HERE / "M23_DUAL_TWO_FIVES_PROJECTOR.json"))
    ap.add_argument("--dual-log", default=str(HERE / "m23_dual_conjugacy_collapse.out"))
    ap.add_argument("--json", default=str(HERE / "M23_RECURSIVE_PULLBACK_CERTIFICATE.json"))
    args = ap.parse_args()

    lines = Path(args.clutching_log).read_text().splitlines()
    a = parse_vector(lines, "MIXED_Q1_LOCAL")
    b = parse_vector(lines, "MIXED_Q2_LOCAL")
    k18 = parse_vector(lines, "MIXED_CENTRAL_INVOLUTION")
    r10 = tuple(parse_vector(lines, "R10_LOCAL"))
    r10_global = parse_vector(lines, "R10_GLOBAL")
    assert len(a) == len(b) == len(k18) == 18
    assert len(r10) == len(r10_global) == 10

    A, B, K = Permutation(a), Permutation(b), Permutation(k18)
    G = PermutationGroup([A, B])
    g_order = int(G.order())
    assert g_order == 92_897_280
    assert G.is_transitive()
    assert K.order() == 2 and all(k18[i] != i for i in range(18))
    assert A * K == K * A and B * K == K * B
    assert G.contains(K)
    center = G.center()
    assert int(center.order()) == 2 and center.contains(K)

    assert subset_image(r10, a) == r10
    assert subset_image(r10, b) != r10
    assert subset_image(r10, k18) == r10

    subsets, reps, schreier = subset_orbit(r10, [a, b])
    subset_degree = len(subsets)
    assert subset_degree == 126
    h_order = g_order // subset_degree
    assert h_order == 737_280 and g_order % subset_degree == 0

    rpos = {x: i for i, x in enumerate(r10)}
    induced = []
    for h in schreier:
        induced.append(Permutation([rpos[h[x]] for x in r10]))
    H10 = PermutationGroup(induced)
    assert int(H10.order()) == 3840 and H10.is_transitive()

    k10_list = [rpos[k18[x]] for x in r10]
    K10 = Permutation(k10_list)
    assert K10.order() == 2
    assert all(k10_list[i] != i for i in range(10))
    assert all(K10 * g == g * K10 for g in H10.generators)

    projector = json.loads(Path(args.projector_json).read_text())
    assert projector["D10_states_zero_based"] == r10_global
    expected_pairs = {tuple(sorted(x)) for x in projector["kappa_pairs"]}
    actual_pairs = {
        tuple(sorted((r10_global[i], r10_global[k10_list[i]])))
        for i in range(10)
        if i < k10_list[i]
    }
    assert actual_pairs == expected_pairs

    blocks = []
    seen = set()
    for i in range(10):
        if i in seen:
            continue
        j = k10_list[i]
        block = tuple(sorted((i, j)))
        blocks.append(block)
        seen.update(block)
    assert len(blocks) == 5
    block_pos = {p: i for i, p in enumerate(blocks)}

    def block_perm(g: Permutation) -> Permutation:
        out = []
        for p in blocks:
            q = tuple(sorted((g(p[0]), g(p[1]))))
            assert q in block_pos
            out.append(block_pos[q])
        return Permutation(out)

    Q5 = PermutationGroup([block_perm(g) for g in H10.generators])
    assert int(Q5.order()) == 120 and Q5.is_transitive()

    a10 = Permutation([rpos[a[x]] for x in r10])
    C10 = PermutationGroup([a10, K10])
    assert int(C10.order()) == 10 and C10.is_transitive()
    C5 = PermutationGroup([block_perm(a10)])
    assert int(C5.order()) == 5 and C5.is_transitive()

    # Compactified branch data over the mixed-boundary M_0,4 base.
    ab = compose(b, a)
    c = inverse(ab)
    sub_a = action_on_subsets(subsets, a)
    sub_b = action_on_subsets(subsets, b)
    sub_c = action_on_subsets(subsets, c)
    assert cycle_hist(sub_a) == {1: 1, 5: 4, 15: 7}
    assert cycle_hist(sub_b) == {1: 1, 5: 4, 15: 7}
    assert cycle_hist(sub_c) == {3: 42}
    assert genus_from_three([sub_a, sub_b, sub_c]) == 31

    subset_pos = {s: i for i, s in enumerate(subsets)}
    lpts = []
    lid = {}
    for si, s in enumerate(subsets):
        for x in s:
            lid[(si, x)] = len(lpts)
            lpts.append((si, x))

    def total_perm(g: list[int]) -> list[int]:
        out = []
        for si, x in lpts:
            s2 = subset_image(subsets[si], g)
            out.append(lid[(subset_pos[s2], g[x])])
        return out

    la, lb, lc = total_perm(a), total_perm(b), total_perm(c)
    assert cycle_hist(la) == {5: 20, 10: 2, 15: 68, 30: 4}
    assert cycle_hist(lb) == {5: 20, 10: 2, 15: 68, 30: 4}
    assert cycle_hist(lc) == {3: 420}
    assert genus_from_three([la, lb, lc]) == 327

    qpts = []
    qid = {}
    for si, s in enumerate(subsets):
        used = set()
        for x in s:
            if x in used:
                continue
            y = k18[x]
            assert y in s
            pair = tuple(sorted((x, y)))
            used.update(pair)
            qid[(si, pair)] = len(qpts)
            qpts.append((si, pair))

    def quotient_perm(g: list[int]) -> list[int]:
        out = []
        for si, pair in qpts:
            s2 = subset_image(subsets[si], g)
            p2 = tuple(sorted((g[pair[0]], g[pair[1]])))
            out.append(qid[(subset_pos[s2], p2)])
        return out

    qa, qb, qc = quotient_perm(a), quotient_perm(b), quotient_perm(c)
    assert cycle_hist(qa) == {5: 12, 15: 38}
    assert cycle_hist(qb) == {5: 12, 15: 38}
    assert cycle_hist(qc) == {3: 210}
    assert genus_from_three([qa, qb, qc]) == 161

    dual_text = Path(args.dual_log).read_text()
    raw_dual_same_inner = "P_INNER_KEYS_EQUAL True" in dual_text and "ALL_DUAL_COLLAPSE_CHECKS_PASS True" in dual_text
    assert raw_dual_same_inner
    boundary_pairs_distinct = all(x != y for x, y in actual_pairs)
    assert boundary_pairs_distinct

    result = {
        "certificate": "M23_RECURSIVE_PULLBACK_CERTIFICATE",
        "mixed_boundary": {
            "degree": 18,
            "monodromy_order": g_order,
            "center_order": 2,
            "central_involution_free": True,
        },
        "ten_state_set": {
            "global_states": r10_global,
            "subset_orbit_size": subset_degree,
            "set_stabilizer_order": h_order,
            "induced_stabilizer_order": int(H10.order()),
            "induced_stabilizer_transitive": True,
            "central_involution_restriction_equals_certified_kappa": True,
            "kappa_pairs": sorted([list(x) for x in actual_pairs]),
            "cyclic_subgroup_q1_kappa_order": int(C10.order()),
            "cyclic_subgroup_transitive": True,
        },
        "five_pair_quotient": {
            "induced_stabilizer_quotient_order": int(Q5.order()),
            "induced_stabilizer_quotient_is_S5": True,
            "q1_quotient_order": int(C5.order()),
        },
        "compactified_covers_over_mixed_boundary_base": {
            "set_stabilizer_base": {
                "degree": 126,
                "branch_cycle_histograms": [cycle_hist(sub_a), cycle_hist(sub_b), cycle_hist(sub_c)],
                "genus": 31,
            },
            "ten_sheet_pullback_component": {
                "degree_over_original_base": 1260,
                "branch_cycle_histograms": [cycle_hist(la), cycle_hist(lb), cycle_hist(lc)],
                "genus": 327,
            },
            "five_sheet_kappa_quotient": {
                "degree_over_original_base": 630,
                "branch_cycle_histograms": [cycle_hist(qa), cycle_hist(qb), cycle_hist(qc)],
                "genus": 161,
            },
        },
        "typed_duality": {
            "boundary_kappa_pairs_distinct_hurwitz_states": boundary_pairs_distinct,
            "raw_h_conjugate_dual_representatives_same_inner_key": raw_dual_same_inner,
            "two_actions_are_distinct_typed_operations": True,
        },
        "scope": {
            "canonical_set_stabilizer_base_change_geometrizes_ten_state_cell": True,
            "global_comparison_source": "five_sheet_kappa_quotient_Q_H",
            "global_comparison_target": "mixed_boundary_quotient_D45_mod_kappa18",
            "cyclic_beta5_carrier": "distinguished_C5_subcarrier_on_ten_state_fiber",
        },
        "pass": True,
    }
    Path(args.json).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print("MIXED_BOUNDARY_MONODROMY_ORDER", g_order)
    print("MIXED_BOUNDARY_CENTER_ORDER", center.order())
    print("R10_SUBSET_ORBIT_SIZE", subset_degree)
    print("R10_SET_STABILIZER_ORDER", h_order)
    print("R10_INDUCED_STABILIZER_ORDER", H10.order())
    print("R10_INDUCED_STABILIZER_TRANSITIVE", H10.is_transitive())
    print("BOUNDARY_KAPPA_EQUALS_LEDGER_KAPPA", actual_pairs == expected_pairs)
    print("KAPPA_QUOTIENT_ACTION_ORDER", Q5.order())
    print("KAPPA_QUOTIENT_ACTION_IS_S5", Q5.order() == 120)
    print("CYCLIC_Q1_KAPPA_SUBGROUP_ORDER", C10.order())
    print("STABILIZER_BASE_GENUS", 31)
    print("TEN_SHEET_PULLBACK_GENUS", 327)
    print("FIVE_SHEET_QUOTIENT_GENUS", 161)
    print("BOUNDARY_KAPPA_PAIRS_DISTINCT_HURWITZ_STATES", boundary_pairs_distinct)
    print("RAW_DUAL_REPRESENTATIVES_SAME_INNER_KEY", raw_dual_same_inner)
    print("ALL_RECURSIVE_PULLBACK_CHECKS_PASS", True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
