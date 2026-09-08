# Galois Covers of M23 over Q: A Positive-Dimensional Approach

Spencer Harkness Bromberg and Alexander Xavier Bromberg  
Revised: August 19, 2026 (June 23, 2026)

Zenodo DOI: `10.5281/zenodo.22010546`

This archive accompanies the article *Galois Covers of M23 over Q: A Positive-Dimensional Approach - Braid-Orbit Reconstruction, Exact Clutching, Algebraic Collision Tails, and Six-Branch Refinements*.

The root contains the final article, its LaTeX source, and a one-page standalone abstract. `computations/core/` contains the exact finite Nielsen, clutching, refinement, collision-tail, cyclic-carrier, mixed-boundary stabilizer-pullback, native ten-state return, width-11 boundary, and comparison-morphism certificates. `computations/arithmetic/` contains the real-boundary, reduced-j-line arithmetic, degree-27 target-algebra, boundary, transport, and coefficient-model certificates. `STRUCTURE.md` gives the complete package map. The subfolder `computations/core/width11_boundary_resolution/` preserves the supplied width-11 source, recorded verification output, and local SHA-256 checksums.

## Verification

Run from the archive root:

```bash
bash RUN_ALL_CERTIFICATES.sh
```

The script creates a temporary working copy, sets `PYTHONDONTWRITEBYTECODE=1`, compiles the C++ certificates in the temporary directory, and replays the public source-level verification suite. Fresh JSON outputs are compared against the recorded mathematical payloads, and embedded `script_sha256` values are checked against the exact source files in the archive. The exhaustive fixed-5A search and reduced-j-line reconstruction are the largest memory stages.

Python dependencies and the tested toolchain are recorded in `requirements.txt` and `ENVIRONMENT.md`. The runner honors `PYTHON` and `CXX` environment variables and defaults to `python3` and `g++`.

For an independent character-theoretic check, `computations/INDEPENDENT_GAP_CROSSCHECK.g` applies the Frobenius formula to GAP's ATLAS-derived CTblLib table for M23; see `computations/INDEPENDENT_GAP_CROSSCHECK.md`.

Verify the shipped files with:

```bash
sha256sum -c SHA256SUMS.txt
```

## Licenses

The article and explanatory documentation use CC-BY-NC-ND-4.0. Source code uses GPL-3.0-only. Computational data and recorded outputs use CC-BY-NC-4.0. See `LICENSES.md`.
