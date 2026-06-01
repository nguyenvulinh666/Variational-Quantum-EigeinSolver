#!/usr/bin/env python3
"""Create table-ready summaries and reviewer-facing plots from revision runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


METHOD_ORDER = [
    "COBYLA",
    "SPSA",
    "PSR",
    "FD",
    "QN-BDA+SPSA",
    "QN-BDA+FD",
    "QN-BDA+PSR",
    "QN-SPSA+SPSA",
    "QN-SPSA+FD",
    "QN-SPSA+PSR",
]

METHOD_STYLES = {
    "COBYLA": {"color": "tab:blue", "linestyle": "-"},
    "SPSA": {"color": "tab:orange", "linestyle": "-"},
    "PSR": {"color": "tab:red", "linestyle": "-"},
    "FD": {"color": "tab:green", "linestyle": "--"},
    "QN-BDA+SPSA": {"color": "tab:purple", "linestyle": "-"},
    "QN-BDA+FD": {"color": "tab:olive", "linestyle": "--"},
    "QN-BDA+PSR": {"color": "tab:brown", "linestyle": "-"},
    "QN-SPSA+SPSA": {"color": "tab:pink", "linestyle": "-"},
    "QN-SPSA+FD": {"color": "tab:cyan", "linestyle": "--"},
    "QN-SPSA+PSR": {"color": "tab:gray", "linestyle": "-"},
}


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def result_dirs(roots: Sequence[Path]) -> List[Path]:
    paths = []
    for root in roots:
        if (root / "metadata.json").exists() and (root / "summary.json").exists():
            paths.append(root)
        else:
            paths.extend(
                path.parent
                for path in sorted(root.rglob("metadata.json"))
                if (path.parent / "summary.json").exists()
            )
    return sorted(set(paths))


def parse_float(value: str) -> Optional[float]:
    if value == "" or value is None:
        return None
    parsed = float(value)
    if math.isnan(parsed):
        return None
    return parsed


def trajectory_for(run_dir: Path, metadata: Dict[str, Any]) -> List[Dict[str, Any]]:
    cost = metadata.get("nominal_quantum_eval_cost", {})
    evals_per_update = cost.get("total_evals_per_update")
    charged_records = 0
    rows = []
    with (run_dir / "trajectory.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            phase = row["phase"]
            if phase in {"update", "objective"} and evals_per_update is not None:
                charged_records += 1
                quantum_evals = charged_records * evals_per_update
            elif evals_per_update is None:
                quantum_evals = None
            else:
                quantum_evals = 0
            rows.append(
                {
                    "step": int(row["step"]),
                    "phase": phase,
                    "relative_error": parse_float(row["relative_error"]),
                    "quantum_evals": quantum_evals,
                }
            )
    return rows


def load_runs(roots: Sequence[Path]) -> List[Dict[str, Any]]:
    runs = []
    for run_dir in result_dirs(roots):
        metadata = read_json(run_dir / "metadata.json")
        summary = read_json(run_dir / "summary.json")
        runs.append(
            {
                "run_dir": run_dir,
                "metadata": metadata,
                "summary": summary,
                "trajectory": trajectory_for(run_dir, metadata),
            }
        )
    return runs


def mean(values: Iterable[Optional[float]]) -> Optional[float]:
    clean = [value for value in values if value is not None]
    if not clean:
        return None
    return statistics.fmean(clean)


def stdev(values: Iterable[Optional[float]]) -> Optional[float]:
    clean = [value for value in values if value is not None]
    if len(clean) < 2:
        return 0.0 if clean else None
    return statistics.stdev(clean)


def write_summary_table(runs: Sequence[Dict[str, Any]], path: Path) -> None:
    groups: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = defaultdict(list)
    for run in runs:
        metadata = run["metadata"]
        key = (
            metadata.get("method"),
            metadata.get("method_display_name"),
            metadata.get("num_qubits"),
            metadata.get("h"),
            metadata.get("ansatz"),
            metadata.get("reps"),
            metadata.get("entanglement"),
        )
        groups[key].append(run)

    fieldnames = [
        "method",
        "method_display_name",
        "num_qubits",
        "h",
        "ansatz",
        "reps",
        "entanglement",
        "runs",
        "final_relative_error_mean",
        "final_relative_error_std",
        "best_relative_error_mean",
        "best_relative_error_std",
        "best_step_mean",
        "nominal_total_quantum_evals_mean",
        "final_energy_mean",
        "best_energy_mean",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for key, group in sorted(groups.items(), key=lambda item: item[0]):
            method, display, num_qubits, h_value, ansatz, reps, entanglement = key
            summaries = [run["summary"] for run in group]
            writer.writerow(
                {
                    "method": method,
                    "method_display_name": display,
                    "num_qubits": num_qubits,
                    "h": h_value,
                    "ansatz": ansatz,
                    "reps": reps,
                    "entanglement": entanglement,
                    "runs": len(group),
                    "final_relative_error_mean": mean(
                        summary.get("final_relative_error") for summary in summaries
                    ),
                    "final_relative_error_std": stdev(
                        summary.get("final_relative_error") for summary in summaries
                    ),
                    "best_relative_error_mean": mean(
                        summary.get("best_relative_error") for summary in summaries
                    ),
                    "best_relative_error_std": stdev(
                        summary.get("best_relative_error") for summary in summaries
                    ),
                    "best_step_mean": mean(
                        summary.get("best_step") for summary in summaries
                    ),
                    "nominal_total_quantum_evals_mean": mean(
                        summary.get("nominal_total_quantum_evals")
                        for summary in summaries
                    ),
                    "final_energy_mean": mean(
                        summary.get("final_energy") for summary in summaries
                    ),
                    "best_energy_mean": mean(
                        summary.get("best_energy") for summary in summaries
                    ),
                }
            )


def method_sort_key(method_name: str) -> Tuple[int, str]:
    try:
        return (METHOD_ORDER.index(method_name), method_name)
    except ValueError:
        return (len(METHOD_ORDER), method_name)


def averaged_series(
    group: Sequence[Dict[str, Any]],
    x_key: str,
) -> List[Tuple[float, float]]:
    by_x: Dict[float, List[float]] = defaultdict(list)
    for run in group:
        for row in run["trajectory"]:
            x_value = row.get(x_key)
            y_value = row.get("relative_error")
            if x_value is None or y_value is None:
                continue
            by_x[float(x_value)].append(float(y_value))
    return [(x_value, statistics.fmean(values)) for x_value, values in sorted(by_x.items())]


def series_stats(
    group: Sequence[Dict[str, Any]],
    x_key: str,
) -> List[Tuple[float, float, float, float, int]]:
    by_x: Dict[float, List[float]] = defaultdict(list)
    for run in group:
        for row in run["trajectory"]:
            x_value = row.get(x_key)
            y_value = row.get("relative_error")
            if x_value is None or y_value is None:
                continue
            by_x[float(x_value)].append(float(y_value))

    stats = []
    for x_value, values in sorted(by_x.items()):
        center = statistics.fmean(values)
        if len(values) > 1:
            lower = percentile(values, 15.865)
            upper = percentile(values, 84.135)
        else:
            lower = center
            upper = center
        stats.append((x_value, center, lower, upper, len(values)))
    return stats


def percentile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("percentile requires at least one value")
    if len(ordered) == 1:
        return ordered[0]

    position = (len(ordered) - 1) * q / 100.0
    lower_index = int(math.floor(position))
    upper_index = int(math.ceil(position))
    if lower_index == upper_index:
        return ordered[lower_index]

    fraction = position - lower_index
    return (
        ordered[lower_index] * (1.0 - fraction)
        + ordered[upper_index] * fraction
    )


def plot_group(
    runs: Sequence[Dict[str, Any]],
    output: Path,
    title: str,
    x_key: str,
    x_label: str,
) -> None:
    import matplotlib.pyplot as plt

    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for run in runs:
        grouped[run["metadata"].get("method_display_name", run["metadata"]["method"])].append(
            run
        )

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for method_name in sorted(grouped, key=method_sort_key):
        stats = series_stats(grouped[method_name], x_key)
        if not stats:
            continue
        x_values = [item[0] for item in stats]
        y_values = [item[1] for item in stats]
        lower_values = [max(item[2], 1e-14) for item in stats]
        upper_values = [max(item[3], 1e-14) for item in stats]
        counts = [item[4] for item in stats]
        style = METHOD_STYLES.get(method_name, {})
        color = style.get("color")

        if max(counts) > 1 and lower_values != upper_values:
            ax.fill_between(
                x_values,
                lower_values,
                upper_values,
                color=color,
                alpha=0.16,
                linewidth=0,
            )

        ax.plot(
            x_values,
            y_values,
            linewidth=1.8,
            label=method_name,
            color=color,
            linestyle=style.get("linestyle", "-"),
        )

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Relative error")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def plot_results(runs: Sequence[Dict[str, Any]], output_dir: Path) -> List[Path]:
    groups: Dict[Tuple[Any, ...], List[Dict[str, Any]]] = defaultdict(list)
    for run in runs:
        metadata = run["metadata"]
        key = (
            metadata.get("num_qubits"),
            metadata.get("h"),
            metadata.get("ansatz"),
            metadata.get("reps"),
            metadata.get("entanglement"),
        )
        groups[key].append(run)

    outputs: List[Path] = []
    for key, group in sorted(groups.items()):
        num_qubits, h_value, ansatz, reps, entanglement = key
        slug = f"n{num_qubits}_h{str(h_value).replace('.', 'p')}_{ansatz}_reps{reps}_{entanglement}"
        title = f"{ansatz}, n={num_qubits}, h={h_value}, reps={reps}"
        iter_path = output_dir / f"relative_error_vs_iteration_{slug}.pdf"
        eval_path = output_dir / f"relative_error_vs_quantum_evals_{slug}.pdf"
        plot_group(group, iter_path, title, "step", "Parameter-update step")
        plot_group(group, eval_path, title, "quantum_evals", "Nominal quantum evaluations")
        outputs.extend([iter_path, eval_path])
    return outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create summary tables and plots from revision VQE runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("roots", nargs="+", help="Run directories or parent directories.")
    parser.add_argument("--out-dir", type=Path, default=Path("results/revision_plots"))
    parser.add_argument("--table", type=Path, default=None)
    parser.add_argument("--no-plots", action="store_true")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    runs = load_runs([Path(root) for root in args.roots])
    if not runs:
        parser.error("No completed revision runs found.")

    table_path = args.table or (args.out_dir / "revision_method_summary.csv")
    write_summary_table(runs, table_path)
    print(table_path)

    if not args.no_plots:
        outputs = plot_results(runs, args.out_dir)
        for output in outputs:
            print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
