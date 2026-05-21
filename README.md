# MEMD C++ demo

This repository contains an old C++ implementation of Multivariate Empirical Mode Decomposition (MEMD).

## Build

```sh
make
```

## Demo data

The code is configured for `N=1001` samples and `NDIM=16` channels in `inc/MEMD_Config.h`.
Generate a compatible synthetic multichannel signal:

```sh
make demo-data
```

This writes `data/demo_signal.csv` with 1001 rows and 16 comma-separated values per row.

## Run

```sh
mkdir -p output
./memd data/demo_signal.csv output/memd
```

or simply:

```sh
make run
```

Outputs are tab-separated matrices:

- `output/memd_imf01.txt`, `output/memd_imf02.txt`, ...
- `output/memd_residue.txt`

## Input format

CSV with exactly `N * NDIM` numeric values, normally arranged as `N` rows by `NDIM` columns.
To use a different shape, update `N` and `NDIM` in `inc/MEMD_Config.h` and rebuild.

## CPU scaling benchmark

Because `N`, `NDIM`, and `NDIR` are compile-time constants, the benchmark script temporarily patches `inc/MEMD_Config.h`, rebuilds, generates matching data, and times the run.

```sh
make benchmark
```

Custom cases use `NxNDIMxNDIR` format:

```sh
python3 scripts/benchmark_cpu.py --case 1001x16x64 --case 5001x16x64 --case 1001x32x128
```

Results are written to `output/benchmark_cpu.csv`. The original config is restored afterwards.
