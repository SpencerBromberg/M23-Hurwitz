#!/usr/bin/env python3
"""Compare freshly generated certificate JSON with the recorded certificate.

Only volatile runtime/environment metadata is ignored.  Mathematical content,
PASS flags, finite counts, cycle data, hashes outside the environment block,
and witness data must agree exactly.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

VOLATILE_KEYS = {
    "environment", "performance", "runtime_seconds", "elapsed_seconds",
    "wall_seconds", "wall_time_seconds", "peak_memory_bytes",
    "peak_memory_kb", "peak_rss_kb", "timestamp", "generated_at",
}

def normalize(x):
    if isinstance(x, dict):
        return {k: normalize(v) for k, v in sorted(x.items()) if k not in VOLATILE_KEYS}
    if isinstance(x, list):
        return [normalize(v) for v in x]
    return x

def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("fresh", type=Path)
    ap.add_argument("recorded", type=Path)
    args=ap.parse_args()
    a=normalize(json.loads(args.fresh.read_text(encoding="utf-8")))
    b=normalize(json.loads(args.recorded.read_text(encoding="utf-8")))
    if a != b:
        print("RECORDED_JSON_MISMATCH", args.recorded)
        return 1
    print("RECORDED_JSON_MATCH", args.recorded)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
