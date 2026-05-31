#!/usr/bin/env python3
"""Aggregate revision-run outputs into a table-ready CSV."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

from revision_experiment import aggregate_result_dirs, write_aggregate_csv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Aggregate metadata.json and summary.json files from revision runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "roots",
        nargs="*",
        default=["results/revision"],
        help="Run directories or parent directories containing revision runs.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output CSV path. Defaults to stdout.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    rows = aggregate_result_dirs([Path(root) for root in args.roots])
    write_aggregate_csv(rows, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
