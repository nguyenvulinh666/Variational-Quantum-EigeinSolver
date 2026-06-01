#!/usr/bin/env python3
"""Generate plain shell command manifests for revision sweeps."""

from __future__ import annotations

import argparse
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional


DETERMINISTIC_METHODS = {"psr", "fd", "cobyla", "qnbda_psr", "qnbda_fd"}
STOCHASTIC_METHODS = {"spsa", "qnbda_spsa", "qnspsa_psr", "qnspsa_fd", "qnspsa_spsa"}
DEFAULT_METHODS = [
    "psr",
    "fd",
    "spsa",
    "cobyla",
    "qnbda_psr",
    "qnbda_fd",
    "qnbda_spsa",
    "qnspsa_psr",
    "qnspsa_fd",
    "qnspsa_spsa",
]


@dataclass(frozen=True)
class ManifestJob:
    method: str
    num_qubits: int
    h: float
    seed: int


def parse_ints(value: str) -> List[int]:
    values: List[int] = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        if ".." in token:
            start_s, end_s = token.split("..", 1)
            start = int(start_s)
            end = int(end_s)
            step = 1 if end >= start else -1
            values.extend(range(start, end + step, step))
        else:
            values.append(int(token))
    return values


def parse_floats(value: str) -> List[float]:
    return [float(token.strip()) for token in value.split(",") if token.strip()]


def parse_methods(value: str) -> List[str]:
    if value == "default":
        return DEFAULT_METHODS.copy()
    return [token.strip() for token in value.split(",") if token.strip()]


def jobs_for(
    methods: Iterable[str],
    qubits: Iterable[int],
    fields: Iterable[float],
    deterministic_seeds: Iterable[int],
    stochastic_seeds: Iterable[int],
) -> List[ManifestJob]:
    jobs: List[ManifestJob] = []
    deterministic_seed_values = list(deterministic_seeds)
    stochastic_seed_values = list(stochastic_seeds)
    for method in methods:
        seeds = (
            stochastic_seed_values
            if method in STOCHASTIC_METHODS
            else deterministic_seed_values
        )
        for num_qubits in qubits:
            for field in fields:
                for seed in seeds:
                    jobs.append(ManifestJob(method, num_qubits, field, seed))
    return jobs


def command_for(args: argparse.Namespace, job: ManifestJob) -> str:
    env_parts = []
    if args.workers:
        env_parts.append(f"VQE_WORKERS={shlex.quote(str(args.workers))}")
    if args.parallel_backend:
        env_parts.append(f"VQE_PARALLEL_BACKEND={shlex.quote(args.parallel_backend)}")

    output_root = Path(args.results_root) / args.label
    command = [
        args.python,
        "run_revision_sweep.py",
        "--method",
        job.method,
        "--n",
        str(job.num_qubits),
        "--h",
        str(job.h),
        "--ansatz",
        args.ansatz,
        "--reps",
        str(args.reps),
        "--entanglement",
        args.entanglement,
        "--iterations",
        str(args.iterations),
        "--learning-rate",
        str(args.learning_rate),
        "--seed",
        str(job.seed),
        "--shots",
        args.shots,
        "--out",
        str(output_root),
        "--run-label",
        args.label,
        "--overwrite",
    ]
    if args.exact_max_qubits is not None:
        command.extend(["--exact-max-qubits", str(args.exact_max_qubits)])

    return " ".join(env_parts + [shlex.join(command)])


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate shell commands for VQE revision cluster sweeps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--preset", choices=["smoke", "full"], default="smoke")
    parser.add_argument("--label", default=None)
    parser.add_argument("--results-root", default="results/revision")
    parser.add_argument("--python", default="python")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--parallel-backend", default=None)
    parser.add_argument("--methods", default="default")
    parser.add_argument("--qubits", default=None)
    parser.add_argument("--fields", default=None)
    parser.add_argument("--deterministic-seeds", default=None)
    parser.add_argument("--stochastic-seeds", default=None)
    parser.add_argument("--iterations", type=int, default=None)
    parser.add_argument("--ansatz", default="RealAmplitudes")
    parser.add_argument("--reps", type=int, default=2)
    parser.add_argument("--entanglement", default="reverse_linear")
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--shots", default="none")
    parser.add_argument("--exact-max-qubits", type=int, default=12)
    parser.add_argument("--out", type=Path, default=None)
    return parser


def defaults_for_preset(args: argparse.Namespace) -> None:
    if args.preset == "smoke":
        args.label = args.label or "revision_smoke"
        args.qubits = args.qubits or "6"
        args.fields = args.fields or "2.0"
        args.deterministic_seeds = args.deterministic_seeds or "0"
        args.stochastic_seeds = args.stochastic_seeds or "0,1"
        args.iterations = args.iterations or 20
        return

    args.label = args.label or "revision_full"
    args.qubits = args.qubits or "6..12"
    args.fields = args.fields or "2.0"
    args.deterministic_seeds = args.deterministic_seeds or "0"
    args.stochastic_seeds = args.stochastic_seeds or "0..6"
    args.iterations = args.iterations or 500


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    defaults_for_preset(args)

    methods = parse_methods(args.methods)
    qubits = parse_ints(args.qubits)
    fields = parse_floats(args.fields)
    deterministic_seeds = parse_ints(args.deterministic_seeds)
    stochastic_seeds = parse_ints(args.stochastic_seeds)

    commands = [
        command_for(args, job)
        for job in jobs_for(
            methods,
            qubits,
            fields,
            deterministic_seeds,
            stochastic_seeds,
        )
    ]

    if args.out is None:
        for command in commands:
            print(command)
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(commands) + "\n", encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
