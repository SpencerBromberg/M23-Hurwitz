#!/usr/bin/env bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
ROOT="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-python3}"
CXX="${CXX:-g++}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp -R "$ROOT/computations" "$TMP/computations"

CORE="$TMP/computations/core"
ARITH="$TMP/computations/arithmetic"
COMPARE="$ROOT/computations/compare_recorded_json.py"
compare_json() {
  "$PYTHON" -B "$COMPARE" "$1" "$2"
}

printf '%s\n' '[core] ordered Nielsen geometry'
cd "$CORE"
"$PYTHON" -B geometry_certificates/verify_m23_ordered_hurwitz.py --json "$TMP/ordered.json"
compare_json "$TMP/ordered.json" "$ROOT/computations/core/geometry_certificates/verify_m23_ordered_hurwitz.json"
"$PYTHON" -B geometry_certificates/verify_m23_nielsen_exhaustiveness.py --json "$TMP/exhaustiveness.json"
compare_json "$TMP/exhaustiveness.json" "$ROOT/computations/core/geometry_certificates/verify_m23_nielsen_exhaustiveness.json"
"$CXX" -O3 -std=c++20 geometry_certificates/m23_five_branch_clutching.cpp -o "$TMP/m23_five_branch"
"$TMP/m23_five_branch" | tee "$TMP/m23_five_branch_clutching.out"
"$PYTHON" -B geometry_certificates/build_m23_four_branch_equations.py
"$PYTHON" -B geometry_certificates/build_m23_five_branch_equations.py

printf '%s\n' '[core] refinements, collision tails, and cyclic carriers'
"$PYTHON" -B base/m23_ten_2a_explicit.py
"$PYTHON" -B base/m23_huang_ten_2a.py
"$PYTHON" -B base/m23_huang_descent_involution.py
"$PYTHON" -B base/m23_huang_3a_split_descent.py
"$PYTHON" -B base/m23_orientation_exact.py > "$TMP/m23_orientation_exact.out"
cmp "$TMP/m23_orientation_exact.out" "$ROOT/computations/core/base/m23_orientation_exact.out"
"$PYTHON" -B base/m23_2a2a_to_2a_orbits.py
"$PYTHON" -B base/m23_common_component_certificate.py
"$PYTHON" -B m23_six_bridge_certificate.py
"$PYTHON" -B m23_beta5_finite_carrier_certificate.py | tee "$TMP/beta5.log"
"$PYTHON" -B m23_collision_tail_certificate.py
"$PYTHON" -B m23_five_pair_tail_certificate.py
"$PYTHON" -B m23_dual_conjugacy_collapse.py | tee "$TMP/dual_conjugacy.out"
"$PYTHON" -B m23_fixed_dual_fiber_certificate.py
"$PYTHON" -B m23_arithmetic_stabilization_certificate.py
"$PYTHON" -B m23_dual_two_fives_projector.py
"$PYTHON" -B m23_recursive_pullback_certificate.py \
  --clutching-log "$TMP/m23_five_branch_clutching.out" \
  --projector-json "$CORE/M23_DUAL_TWO_FIVES_PROJECTOR.json" \
  --dual-log "$TMP/dual_conjugacy.out" \
  --json "$TMP/recursive_pullback.json"
"$PYTHON" -B m23_comparison_morphism_certificate.py \
  --pullback-json "$TMP/recursive_pullback.json" \
  --beta5-log "$TMP/beta5.log" \
  --json "$TMP/comparison_morphism.json"
"$CXX" -O3 -std=c++20 m23_dual_two_fives_certificate.cpp -o "$TMP/m23_dual_two_fives"
"$TMP/m23_dual_two_fives"
"$CXX" -O3 -std=c++17 m23_ten_width11_completion_certificate.cpp -o "$TMP/m23_ten_width11_completion"
"$TMP/m23_ten_width11_completion" | tee "$TMP/m23_ten_width11_completion.out"
"$CXX" -O3 -std=c++17 width11_boundary_resolution/verify_m23_width11_boundary_resolution.cpp -o "$TMP/m23_width11_boundary_resolution"
"$TMP/m23_width11_boundary_resolution" > "$TMP/m23_width11_boundary_resolution.out"
cmp "$TMP/m23_width11_boundary_resolution.out" width11_boundary_resolution/verification.log
"$CXX" -O3 -std=c++20 m23_raw_dickson_separation_certificate.cpp -o "$TMP/m23_raw_dickson"
"$TMP/m23_raw_dickson"

printf '%s\n' '[arithmetic] target algebra and finite arithmetic checks'
cd "$ARITH"
"$PYTHON" -B code/verify_cover_fiber_activation.py --input data/b27_multiplication_input.json --output "$TMP/b27_activation.json"
compare_json "$TMP/b27_activation.json" "$ROOT/computations/arithmetic/results/M23_b27_cover_fiber_activation.json"
"$PYTHON" -B code/extract_rank_one_idempotent.py --expr 'Z^27-Z' --output "$TMP/b27_idempotent.json"
compare_json "$TMP/b27_idempotent.json" "$ROOT/computations/arithmetic/results/M23_b27_idempotent.json"
(cd code && "$PYTHON" -B build_r6_certificate.py) > "$TMP/r6.log"
(cd code && "$PYTHON" -B build_m23_universal_model.py) > "$TMP/universal_model.log"
cmp "$ARITH/data/M23_universal_equations.json" "$ROOT/computations/arithmetic/data/M23_universal_equations.json"
cmp "$ARITH/data/M23_universal_equations.txt" "$ROOT/computations/arithmetic/data/M23_universal_equations.txt"
cmp "$ARITH/data/BUILD_SUMMARY.json" "$ROOT/computations/arithmetic/data/BUILD_SUMMARY.json"
cmp "$TMP/universal_model.log" "$ROOT/computations/arithmetic/results/M23_universal_model.log"
"$CXX" -O3 -std=c++17 code/compute_color_diagonal_pairs.cpp -o "$TMP/color_diagonal"
"$TMP/color_diagonal" > "$TMP/color_diagonal.log"
"$CXX" -O3 -std=c++17 code/compute_point_push.cpp -o "$TMP/point_push"
"$TMP/point_push" > "$TMP/point_push.log"

"$PYTHON" -B code/verify_m23_strict_power_model.py --input data/orbit_certificate.json.gz --json "$TMP/strict_power.json"
compare_json "$TMP/strict_power.json" "$ROOT/computations/arithmetic/results/M23_strict_power_model.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_hurwitz_certificate.py --json "$TMP/primary.json"
compare_json "$TMP/primary.json" "$ROOT/computations/arithmetic/results/M23_primary.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_nielsen_exhaustiveness.py --json "$TMP/arithmetic_exhaustiveness.json"
compare_json "$TMP/arithmetic_exhaustiveness.json" "$ROOT/computations/arithmetic/results/M23_exhaustiveness.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_real_boundary.py --json "$TMP/real_boundary.json"
compare_json "$TMP/real_boundary.json" "$ROOT/computations/arithmetic/results/M23_real_boundary.json"
"$PYTHON" -B code/verify_m23_twist_tolerant_boundary.py --input data/orbit_certificate.json.gz --json "$TMP/frame_free_transport.json"
compare_json "$TMP/frame_free_transport.json" "$ROOT/computations/arithmetic/results/M23_frame_free_transport.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_five_branch_component.py --json "$TMP/five_branch_component.json"
compare_json "$TMP/five_branch_component.json" "$ROOT/computations/arithmetic/results/M23_five_branch_component.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_trace_closure.py --json "$TMP/trace_closure.json"
compare_json "$TMP/trace_closure.json" "$ROOT/computations/arithmetic/results/M23_trace_closure.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_rank3_harmonic_extension.py --json "$TMP/rank3_harmonic_extension.json"
compare_json "$TMP/rank3_harmonic_extension.json" "$ROOT/computations/arithmetic/results/M23_rank3_harmonic_extension.json"
PYTHONPATH=code "$PYTHON" -B code/verify_m23_order11_boundary.py --json "$TMP/order11_boundary.json"
compare_json "$TMP/order11_boundary.json" "$ROOT/computations/arithmetic/results/M23_order11_boundary_certificate.json"
PYTHONPATH=code/family_211 "$PYTHON" -B code/family_211/verify_m23_complete_211_family.py \
  --points results/family_211/stage5_zero_points.txt \
  --exhaustive-log results/family_211/exhaustive_stage5.log \
  --json "$TMP/family_211.json"
compare_json "$TMP/family_211.json" "$ROOT/computations/arithmetic/results/M23_complete_211_family.json"

if [[ "${QUICK:-0}" != "1" ]]; then
  printf '%s\n' '[heavy] reduced j-line reconstructions'
  cd "$CORE"
  "$PYTHON" -B geometry_certificates/verify_m23_reduced_jline.py --json "$TMP/reduced_jline.json"
  compare_json "$TMP/reduced_jline.json" "$ROOT/computations/core/geometry_certificates/verify_m23_reduced_jline.json"
  cd "$ARITH"
  PYTHONPATH=code "$PYTHON" -B code/construct_m23_reduced_jline.py --json "$TMP/arithmetic_reduced_jline.json"
  compare_json "$TMP/arithmetic_reduced_jline.json" "$ROOT/computations/arithmetic/results/M23_reduced_jline.json"
fi

printf '%s\n' '[integrity] embedded script provenance'
"$PYTHON" -B "$ROOT/computations/verify_embedded_script_hashes.py"

if command -v gap >/dev/null 2>&1; then
  printf '%s\n' '[independent] GAP/CTblLib Frobenius cross-check'
  gap -q "$ROOT/computations/INDEPENDENT_GAP_CROSSCHECK.g"
else
  printf '%s\n' '[independent] GAP/CTblLib cross-check skipped (gap not installed)'
fi

printf '%s\n' 'ALL_PUBLIC_CERTIFICATE_CHECKS_PASS'
