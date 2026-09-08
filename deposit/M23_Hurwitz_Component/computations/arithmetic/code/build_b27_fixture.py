#!/usr/bin/env python3
"""Build the exact multiplication-by-Z matrix for Q[Z]/(Z^27-Z)."""
import json
from pathlib import Path
n=27
M=[[0 for _ in range(n)] for _ in range(n)]
# Z * Z^j = Z^(j+1), j=0,...,25
for j in range(26):
    M[j+1][j]=1
# Z * Z^26 = Z^27 = Z
M[1][26]=1
one=[1]+[0]*26
out={
  "label":"B27 exact finite-etale fixture Q[Z]/(Z^27-Z)",
  "basis":[f"Z^{i}" for i in range(n)],
  "unit_vector":one,
  "multiplication_matrix":M,
  "defining_relation":"Z^27-Z"
}
Path("b27_multiplication_input.json").write_text(json.dumps(out,indent=2)+"\n")
