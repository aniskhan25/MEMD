#!/usr/bin/env python3
"""Generate synthetic multivariate signals for MEMD demos/benchmarks."""
import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=1001, help="number of samples")
    parser.add_argument("--ndim", type=int, default=16, help="number of channels/dimensions")
    parser.add_argument("--out", default="data/demo_signal.csv", help="output CSV path")
    args = parser.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="") as f:
        writer = csv.writer(f)
        for i in range(args.n):
            t = i / (args.n - 1)
            row = []
            for d in range(args.ndim):
                phase = 2.0 * math.pi * d / args.ndim
                # Shared oscillatory modes with dimension-specific phase/amplitude,
                # plus a slow trend. This is suitable for exercising MEMD.
                value = (
                    math.sin(2 * math.pi * 7 * t + phase)
                    + 0.45 * math.sin(2 * math.pi * 23 * t + 0.5 * phase)
                    + 0.20 * math.sin(2 * math.pi * 61 * t + 1.7 * phase)
                    + 0.30 * (t - 0.5) * math.cos(phase)
                )
                row.append(f"{value:.10f}")
            writer.writerow(row)

    print(f"wrote {out} ({args.n} rows x {args.ndim} columns)")


if __name__ == "__main__":
    main()
