# Reproducibility environment

The release was tested on Linux x86-64 with Python 3.13.5, NumPy 2.3.5,
SymPy 1.14.0, and g++ 14.2.0. The replay wrapper accepts `PYTHON` and `CXX`
environment variables and defaults to `python3` and `g++`.

The C++ sources are compiled under the language level required by each
certificate (`-std=c++17` or `-std=c++20` as recorded in
`RUN_ALL_CERTIFICATES.sh`). Exact Python dependencies are pinned in the root
`requirements.txt`.

Representative measurements on the release-build container are:

| Certificate stage | Wall time | Peak RSS |
| --- | ---: | ---: |
| Ordered 980-state reconstruction | 2.2 s | 169 MB |
| Fixed-5A exhaustive pair count | 12.3 s | 232 MB |
| Five-branch 21,456-state component | 28.5 s | 368 MB |
| Trace closure | 20.1 s | 173 MB |
| Rank-three harmonic extension | 20.2 s | 174 MB |

Runtime is hardware dependent. An independent audit run of the heavier
reconstruction stages reported peak memory below 3 GB.

`computations/INDEPENDENT_GAP_CROSSCHECK.g` is an optional independent
character-theoretic check and requires GAP with CTblLib. GAP is intentionally
not a dependency of the main Python/C++ replay.
