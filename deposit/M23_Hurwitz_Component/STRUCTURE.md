# Package structure

This deposit is organized by mathematical role.

| Path | Contents |
|---|---|
| `M23_Hurwitz_Component.pdf` | Final article for citation and reading |
| `M23_Hurwitz_Component.tex` | Final LaTeX source |
| `ABSTRACT.pdf` | One-page standalone abstract |
| `ABSTRACT.tex`, `ABSTRACT.txt` | Standalone abstract source and plain text |
| `computations/core/` | Four-branch and five-branch Nielsen geometry, clutching, collision tails, common refinements, cyclic carriers, the canonical ten-state stabilizer pullback, the native ten-state return, the width-11 boundary decomposition, the comparison-morphism certificate, and arithmetic stabilization |
| `computations/core/width11_boundary_resolution/` | Supplied width-11 C++ verifier, recorded verification output, README, and local SHA-256 checksums |
| `computations/arithmetic/` | Real-boundary arithmetic, reduced j-line arithmetic, degree-27 target algebra, boundary and transport checks, coefficient-model data, and 2/11-family certificates |
| `RUN_ALL_CERTIFICATES.sh` | Source-level verification from a temporary working copy, including fresh-vs-recorded JSON comparison |
| `computations/INDEPENDENT_GAP_CROSSCHECK.g` | Optional independent Frobenius character-formula check using GAP/CTblLib |
| `requirements.txt`, `ENVIRONMENT.md` | Pinned Python dependencies and tested toolchain/runtime notes |
| `SHA256SUMS.txt` | Cryptographic manifest for the deposit |
| `CITATION.cff` | Citation metadata and Zenodo DOI |
| `LICENSES.md` | License scope for manuscript, code, and data |
| `GPL-3.0-only.txt` | Full GPL-3.0-only text for source code |

The public archive contains the final article and the mathematical certificate materials required to verify the stated computations. The directory names describe mathematical content and expose a single publication layer.
