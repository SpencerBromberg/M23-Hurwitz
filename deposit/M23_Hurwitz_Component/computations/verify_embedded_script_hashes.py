#!/usr/bin/env python3
"""Verify embedded script_sha256 fields against source files in this release."""
from __future__ import annotations
import hashlib, json, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def hashes_in(x):
    if isinstance(x, dict):
        for k,v in x.items():
            if k == "script_sha256" and isinstance(v,str):
                yield v
            yield from hashes_in(v)
    elif isinstance(x,list):
        for v in x: yield from hashes_in(v)

def main() -> int:
    sources={}
    for p in ROOT.rglob("*"):
        if p.is_file() and p.suffix in {".py",".cpp",".sh"}:
            h=hashlib.sha256(p.read_bytes()).hexdigest()
            sources.setdefault(h,[]).append(p.relative_to(ROOT).as_posix())
    bad=[]; checked=0
    for p in ROOT.rglob("*.json"):
        try: data=json.loads(p.read_text(encoding="utf-8"))
        except Exception: continue
        for h in hashes_in(data):
            checked += 1
            if h not in sources:
                bad.append((p.relative_to(ROOT).as_posix(),h))
    print("EMBEDDED_SCRIPT_HASHES_CHECKED",checked)
    if bad:
        for p,h in bad: print("STALE_SCRIPT_SHA256",p,h)
        return 1
    print("ALL_EMBEDDED_SCRIPT_HASHES_RESOLVE True")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
