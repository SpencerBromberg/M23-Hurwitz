M23 five-branch width-11 boundary resolution certificate

Input/convention:
  Reconstructs the 21,456-state colored Nielsen component for
  (2A,2A,2A,2A,3A) from the explicit degree-23 seed used in the manuscript.
  The mixed collision operator is q4^2.

Certified output:
  * q4^2 has exactly 264 cycles of width 11.
  * q1 and q2 induce four transitive components, each of degree 66.
  * On each component the reduced branch-cycle triple has passports
      gamma0: 3^22
      gamma1: 2^33
      gammainf: 3^1 4^3 5^3 6^6
    and Riemann-Hurwitz gives genus 0.
  * kappa11=(q1 q2 q3)c in the recorded right-action convention is an
    involution, conjugates q4^2 to (q4^2)^(-1), has no fixed width-11 cusp,
    and exchanges components 0<->2 and 1<->3.

Compile/run:
  g++ -O3 -std=c++17 verify_m23_width11_boundary_resolution.cpp -o verify
  ./verify

The verification.log included here is the output of that command in the
current container.
