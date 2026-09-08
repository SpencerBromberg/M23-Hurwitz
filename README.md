# Galois Covers of M23 over Q: A Positive-Dimensional Approach

[![Replay certificates](https://github.com/SpencerBromberg/M23-Hurwitz/actions/workflows/replay.yml/badge.svg)](https://github.com/SpencerBromberg/M23-Hurwitz/actions/workflows/replay.yml)

Spencer Harkness Bromberg and Alexander Xavier Bromberg

**Version of record:** Zenodo, DOI [10.5281/zenodo.22010546](https://doi.org/10.5281/zenodo.22010546)

This repository is a convenience mirror. `deposit/M23_Hurwitz_Component/` is a byte-identical copy of the Zenodo deposit: the article, its LaTeX source, the standalone abstract, and the full computational certificate suite. Nothing in that directory is edited here; changes are made on Zenodo as new versions and then mirrored.

| Git tag | Zenodo version |
|---|---|
| `zenodo-v1` | Version 1 (2026-08-19) |

## Verify the mirror against the deposit

```bash
cd deposit/M23_Hurwitz_Component
sha256sum -c SHA256SUMS.txt
```

## Replay the certificates

```bash
cd deposit/M23_Hurwitz_Component
pip install -r requirements.txt
bash RUN_ALL_CERTIFICATES.sh        # QUICK=1 skips the two heaviest stages
```

The runner compiles the C++ certificates in a temporary directory, replays the Python and C++ verification suite, compares fresh JSON output against the recorded payloads, and checks embedded `script_sha256` values against the source files. Toolchain and runtime figures are in `ENVIRONMENT.md`. A GAP/CTblLib character-formula cross-check runs automatically if `gap` is on the path.

The same replay runs in CI on every push (see the badge above), and can be triggered manually from the Actions tab.

`computations/core/width11_boundary_resolution/` is an independently generated second implementation of the width-11 boundary check, kept alongside the primary certificate `m23_ten_width11_completion_certificate.cpp` so the two can be compared.

## Citation

Use the DOI above, or the "Cite this repository" button (`CITATION.cff`).

## Licenses

- Article, abstract, and explanatory documentation: CC-BY-NC-ND-4.0
- Source code in `computations/` and `RUN_ALL_CERTIFICATES.sh`: GPL-3.0-only
- Computational data and recorded certificate outputs: CC-BY-NC-4.0

Full scope in `LICENSE.md`; license texts in `LICENSES/`.
