#!/usr/bin/env python3
"""Figures for the Temperature 12k multi-method composite.

  reconstruction_gmst.png    -- consensus median + 5-95% band vs age (ka BP)
  reconstruction_methods.png -- per-method medians overlaid
  sample_depth.png           -- pooled ensemble members contributing per bin

Every PNG in the figures dir is surfaced on the static-Pages tile UI.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make(in_csv: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(in_csv)
    if "age_bp" in df.columns:
        kyr = df["age_bp"].to_numpy(dtype=float) / 1000.0
    else:
        kyr = (1950 - df["year"].to_numpy(dtype=float)) / 1000.0

    # 1) consensus GMST with uncertainty band
    fig, ax = plt.subplots(figsize=(9, 4.5))
    if {"lo_95", "hi_95"}.issubset(df.columns):
        ax.fill_between(kyr, df["lo_95"], df["hi_95"], color="0.8",
                        label="5-95% (multi-method)", linewidth=0)
    ax.plot(kyr, df["mean"], color="firebrick", lw=2, label="consensus median")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))      # 12 ka left, present right
    ax.set_xlabel("Age (ka BP)")
    ax.set_ylabel("Temperature anomaly (°C, vs 1800-1900)")
    ax.set_title("Temperature 12k multi-method composite GMST")
    ax.legend(loc="lower left", fontsize=8, frameon=False)
    ax.grid(alpha=0.3, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(out_dir / "reconstruction_gmst.png", dpi=130)
    plt.close(fig)

    # 2) per-method medians overlay
    meth_cols = [c for c in df.columns if c.endswith("_median")]
    if meth_cols:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        for c in meth_cols:
            ax.plot(kyr, df[c], lw=1.3, label=c.replace("_median", "").upper())
        ax.plot(kyr, df["mean"], color="k", lw=2.2, label="consensus")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))
        ax.set_xlabel("Age (ka BP)")
        ax.set_ylabel("Temperature anomaly (°C)")
        ax.set_title("Per-method medians")
        ax.legend(loc="lower left", fontsize=8, ncol=2, frameon=False)
        ax.grid(alpha=0.3, linewidth=0.5)
        fig.tight_layout()
        fig.savefig(out_dir / "reconstruction_methods.png", dpi=130)
        plt.close(fig)

    # 3) sample depth
    if "n_members" in df.columns:
        fig, ax = plt.subplots(figsize=(9, 2.8))
        ax.fill_between(kyr, 0, df["n_members"], color="steelblue", alpha=0.7)
        ax.set_xlim(np.nanmax(kyr), np.nanmin(kyr))
        ax.set_xlabel("Age (ka BP)")
        ax.set_ylabel("ensemble members")
        ax.set_title("Pooled ensemble depth")
        fig.tight_layout()
        fig.savefig(out_dir / "sample_depth.png", dpi=130)
        plt.close(fig)

    print(f"[make_figures] wrote figures to {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    make(args.in_csv, args.out_dir)


if __name__ == "__main__":
    main()
