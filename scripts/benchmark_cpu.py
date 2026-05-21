#!/usr/bin/env python3
"""Rebuild and time the CPU MEMD implementation at several problem sizes.

The MEMD code uses compile-time constants in inc/MEMD_Config.h, so this script
patches those values, rebuilds, generates matching synthetic data, and runs the
binary for each case. The original config is restored when the script exits.
"""
import argparse
import csv
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "inc" / "MEMD_Config.h"

DEFAULT_CASES = [
    # n, ndim, ndir. Keep these modest enough for a laptop CPU.
    (1001, 16, 64),
    (2001, 16, 64),
    (5001, 16, 64),
    (1001, 32, 64),
    (1001, 16, 128),
]


def run(cmd, *, stdout=None):
    return subprocess.run(cmd, cwd=ROOT, check=True, text=True, stdout=stdout)


def patch_define(text: str, name: str, value: int) -> str:
    return re.sub(rf"^#define\s+{name}\s+\d+", f"#define {name} {value}", text, flags=re.M)


def set_config(n: int, ndim: int, ndir: int, original: str) -> None:
    text = original
    text = patch_define(text, "N", n)
    text = patch_define(text, "NDIM", ndim)
    text = patch_define(text, "NDIR", ndir)
    CONFIG.write_text(text)


def parse_cases(values):
    if not values:
        return DEFAULT_CASES
    cases = []
    for value in values:
        try:
            n, ndim, ndir = [int(x) for x in value.lower().split("x")]
        except ValueError as exc:
            raise SystemExit(f"Invalid case '{value}'. Use format NxNDIMxNDIR, e.g. 5001x16x64") from exc
        cases.append((n, ndim, ndir))
    return cases


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", help="benchmark case as NxNDIMxNDIR; may be repeated")
    parser.add_argument("--results", default="output/benchmark_cpu.csv")
    parser.add_argument("--keep-outputs", action="store_true", help="do not delete generated IMF/residue files")
    args = parser.parse_args()

    cases = parse_cases(args.case)
    original = CONFIG.read_text()
    results_path = ROOT / args.results
    results_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    try:
        for n, ndim, ndir in cases:
            label = f"N{n}_D{ndim}_R{ndir}"
            data_path = ROOT / "data" / f"bench_{label}.csv"
            output_prefix = ROOT / "output" / f"bench_{label}"

            print(f"\n=== {label} ===", flush=True)
            set_config(n, ndim, ndir, original)

            t0 = time.perf_counter()
            run(["make", "clean"], stdout=subprocess.DEVNULL)
            run(["make"], stdout=subprocess.DEVNULL)
            build_sec = time.perf_counter() - t0

            run(["python3", "scripts/generate_demo_data.py", "--n", str(n), "--ndim", str(ndim), "--out", str(data_path)])

            with (ROOT / "output" / f"{label}.log").open("w") as log:
                t0 = time.perf_counter()
                run(["./memd", str(data_path), str(output_prefix)], stdout=log)
                run_sec = time.perf_counter() - t0

            imfs = len(list((ROOT / "output").glob(f"bench_{label}_imf*.txt")))
            rows.append({
                "N": n,
                "NDIM": ndim,
                "NDIR": ndir,
                "values": n * ndim,
                "projection_points": n * ndim * ndir,
                "build_sec": f"{build_sec:.3f}",
                "run_sec": f"{run_sec:.3f}",
                "imfs": imfs,
            })
            print(f"run_sec={run_sec:.3f}, imfs={imfs}", flush=True)

            if not args.keep_outputs:
                for path in (ROOT / "output").glob(f"bench_{label}_*.txt"):
                    path.unlink()
    finally:
        CONFIG.write_text(original)
        # Rebuild the default binary after restoring config.
        run(["make"], stdout=subprocess.DEVNULL)

    with results_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nwrote {results_path}")


if __name__ == "__main__":
    main()
