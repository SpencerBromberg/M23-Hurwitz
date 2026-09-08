# Independent GAP/CTblLib cross-check

`INDEPENDENT_GAP_CROSSCHECK.g` uses only GAP's ATLAS-derived CTblLib table for
`M23` and the Frobenius character formula. It is independent of the Python and
C++ braid-orbit implementations in this release.

It checks the fixed-target class-product counts 98, 2688, and 127200. The first
two independently confirm the totals underlying the recorded local orbit
splittings `84+14` and `1344+896+224+224`. The third independently confirms the
full fixed-5A product-one count used by the four-branch exhaustiveness
certificate before the generating/intransitive split. 
Run with a GAP installation containing CTblLib:

```text
gap -q computations/INDEPENDENT_GAP_CROSSCHECK.g
```

The master replay runs this cross-check automatically when GAP with CTblLib is
available and otherwise records a clean skip, so the Python/C++ suite remains
self-contained in the environment documented in `ENVIRONMENT.md`.
