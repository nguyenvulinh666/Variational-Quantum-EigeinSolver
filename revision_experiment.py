"""Structured experiment layer for manuscript-revision VQE runs.

This module intentionally treats ``CoreVQEModified.py`` as the legacy numerical
backend. The classes here organize configuration, optimizer dispatch, cost
accounting, result writing, and aggregation without rewriting the optimizer
math.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Type

import numpy as np


def parse_shots(value: str) -> Optional[int]:
    if value.lower() in {"none", "null", "statevector", "exact"}:
        return None
    shots = int(value)
    if shots <= 0:
        raise argparse.ArgumentTypeError("--shots must be positive or 'none'")
    return shots


def safe_float_token(value: float) -> str:
    return str(value).replace("-", "m").replace(".", "p")


def git_commit() -> Optional[str]:
    try:
        output = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return output.strip()
    except Exception:
        return None


@dataclass(frozen=True)
class ExperimentConfig:
    method: str
    num_qubits: int
    h: float
    j: float = 1.0
    ansatz: str = "RealAmplitudes"
    reps: int = 1
    entanglement: str = "reverse_linear"
    iterations: int = 500
    learning_rate: float = 0.01
    shots: Optional[int] = None
    seed: int = 0
    initial_mode: str = "constant"
    initial_value: float = -0.5
    initial_scale: float = 0.5
    initial_point_file: Optional[str] = None
    initial_point_delimiter: Optional[str] = None
    history_length: int = 5
    exact_max_qubits: int = 12
    save_fubini_history: bool = False
    out: Path = Path("results/revision")
    run_label: Optional[str] = None
    overwrite: bool = False
    command: str = ""

    @classmethod
    def from_args(
        cls,
        args: argparse.Namespace,
        method: str,
        command: Optional[Sequence[str]] = None,
    ) -> "ExperimentConfig":
        command_text = " ".join(command if command is not None else sys.argv)
        return cls(
            method=method,
            num_qubits=args.num_qubits,
            h=args.h,
            j=args.j,
            ansatz=args.ansatz,
            reps=args.reps,
            entanglement=args.entanglement,
            iterations=args.iterations,
            learning_rate=args.learning_rate,
            shots=args.shots,
            seed=args.seed,
            initial_mode=args.initial_mode,
            initial_value=args.initial_value,
            initial_scale=args.initial_scale,
            initial_point_file=args.initial_point_file,
            initial_point_delimiter=args.initial_point_delimiter,
            history_length=args.history_length,
            exact_max_qubits=args.exact_max_qubits,
            save_fubini_history=args.save_fubini_history,
            out=Path(args.out),
            run_label=args.run_label,
            overwrite=args.overwrite,
            command=command_text,
        )

    def run_dir(self) -> Path:
        label = self.run_label
        if label is None:
            label = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        slug = (
            f"{label}__method-{self.method}__ansatz-{self.ansatz}"
            f"__n-{self.num_qubits}__h-{safe_float_token(self.h)}"
            f"__seed-{self.seed}"
        )
        return self.out / slug


@dataclass
class TrajectoryRecord:
    step: int
    phase: str
    energy: float
    elapsed_s: float
    parameters: np.ndarray

    def to_csv_row(self, exact_energy: Optional[float]) -> Dict[str, Any]:
        if exact_energy is not None and abs(exact_energy) > 0:
            rel_error: Any = abs((float(self.energy) - exact_energy) / exact_energy)
        else:
            rel_error = ""

        return {
            "step": self.step,
            "phase": self.phase,
            "energy": self.energy,
            "relative_error": rel_error,
            "elapsed_s": self.elapsed_s,
            "parameters_json": json.dumps(
                np.asarray(self.parameters, dtype=float).tolist(),
                separators=(",", ":"),
            ),
        }


@dataclass
class RunResult:
    run_dir: Path
    metadata: Dict[str, Any]
    trajectory: List[TrajectoryRecord]
    summary: Dict[str, Any]
    fubini_history: List[np.ndarray] = field(default_factory=list)


@dataclass(frozen=True)
class CostBreakdown:
    gradient_evals_per_update: Optional[int]
    metric_evals_per_update: Optional[int]
    total_evals_per_update: Optional[int]
    qnbda_parameter_layers: Optional[int]

    def to_dict(self) -> Dict[str, Optional[int]]:
        return {
            "gradient_evals_per_update": self.gradient_evals_per_update,
            "metric_evals_per_update": self.metric_evals_per_update,
            "total_evals_per_update": self.total_evals_per_update,
            "qnbda_parameter_layers": self.qnbda_parameter_layers,
        }


class CostModel:
    """Nominal quantum-evaluation accounting used for reviewer plots/tables."""

    @staticmethod
    def count_qnbda_layers(core_module: Any, ansatz: Any) -> Optional[int]:
        try:
            layers = core_module.Separate_Circuit_Apart(ansatz)
        except Exception:
            return None
        return sum(1 for layer in layers if getattr(layer, "num_parameters", 0) > 0)

    @classmethod
    def for_method(cls, core_module: Any, method: str, ansatz: Any) -> CostBreakdown:
        p = int(ansatz.num_parameters)
        qnbda_layers = cls.count_qnbda_layers(core_module, ansatz)
        metric_cost = {
            "psr": 0,
            "fd": 0,
            "spsa": 0,
            "qnbda_psr": qnbda_layers,
            "qnbda_spsa": qnbda_layers,
            "qnspsa_psr": 4,
            "qnspsa_spsa": 4,
            "qnspsa_psr_mc": None,
            "qnspsa_spsa_mc": None,
            "cobyla": 0,
        }[method]
        gradient_cost = {
            "psr": 2 * p,
            "fd": 2 * p,
            "spsa": 2,
            "qnbda_psr": 2 * p,
            "qnbda_spsa": 2,
            "qnspsa_psr": 2 * p,
            "qnspsa_spsa": 2,
            "qnspsa_psr_mc": 2 * p,
            "qnspsa_spsa_mc": 2,
            "cobyla": 1,
        }[method]

        if metric_cost is None or gradient_cost is None:
            total_cost = None
        else:
            total_cost = int(metric_cost + gradient_cost)

        return CostBreakdown(
            gradient_evals_per_update=gradient_cost,
            metric_evals_per_update=metric_cost,
            total_evals_per_update=total_cost,
            qnbda_parameter_layers=qnbda_layers,
        )


class ResultWriter:
    FIELDNAMES = [
        "step",
        "phase",
        "energy",
        "relative_error",
        "elapsed_s",
        "parameters_json",
    ]

    @staticmethod
    def write_json(path: Path, payload: Dict[str, Any]) -> None:
        with path.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")

    def write(self, result: RunResult) -> None:
        result.run_dir.mkdir(parents=True, exist_ok=True)
        self.write_json(result.run_dir / "metadata.json", result.metadata)
        self.write_json(result.run_dir / "summary.json", result.summary)

        parameters = np.asarray(
            [record.parameters for record in result.trajectory],
            dtype=float,
        )
        np.save(result.run_dir / "parameters.npy", parameters)

        exact_energy = result.metadata.get("exact_ground_energy")
        with (result.run_dir / "trajectory.csv").open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=self.FIELDNAMES)
            writer.writeheader()
            for record in result.trajectory:
                writer.writerow(record.to_csv_row(exact_energy))

        if result.fubini_history:
            np.save(
                result.run_dir / "fubini_history.npy",
                np.asarray(result.fubini_history, dtype=float),
            )


def bind_ansatz(ansatz: Any, parameters: np.ndarray) -> Any:
    return ansatz.bind_parameters(
        {theta: parameters[i] for i, theta in enumerate(ansatz.parameters)}
    )


def build_ansatz(ansatz_name: str, num_qubits: int, entanglement: str, reps: int) -> Any:
    from qiskit.circuit.library import EfficientSU2, RealAmplitudes

    if ansatz_name == "RealAmplitudes":
        return RealAmplitudes(
            num_qubits,
            entanglement=entanglement,
            reps=reps,
            insert_barriers=True,
        ).decompose()
    if ansatz_name == "EfficientSU2":
        return EfficientSU2(
            num_qubits,
            entanglement=entanglement,
            reps=reps,
            insert_barriers=True,
        ).decompose()
    raise ValueError(f"Unknown ansatz '{ansatz_name}'")


def initial_parameters(config: ExperimentConfig, num_parameters: int) -> np.ndarray:
    if config.initial_point_file:
        path = Path(config.initial_point_file)
        if path.suffix == ".npy":
            values = np.load(path)
        else:
            values = np.loadtxt(path, delimiter=config.initial_point_delimiter)
        values = np.asarray(values, dtype=float).reshape(-1)
        if values.size != num_parameters:
            raise ValueError(
                f"Initial point has {values.size} values, expected {num_parameters}"
            )
        return values

    if config.initial_mode == "uniform":
        return np.random.uniform(
            low=-config.initial_scale,
            high=config.initial_scale,
            size=num_parameters,
        )

    return np.zeros(num_parameters, dtype=float) + config.initial_value


def exact_ground_energy(
    hamiltonian: Any,
    num_qubits: int,
    exact_max_qubits: int,
) -> Optional[float]:
    if num_qubits > exact_max_qubits:
        return None
    try:
        matrix = np.asarray(hamiltonian.to_matrix(), dtype=complex)
        return float(np.linalg.eigvalsh(matrix).min().real)
    except Exception:
        return None


@dataclass
class ExperimentContext:
    config: ExperimentConfig
    core_module: Any
    hamiltonian: Any
    ansatz: Any
    initial_point: np.ndarray
    sampler: Any
    trajectory: List[TrajectoryRecord]
    fubini_history: List[np.ndarray]
    started: float

    def elapsed(self) -> float:
        return time.time() - self.started

    def record(self, phase: str, energy: float, parameters: Iterable[float]) -> None:
        self.trajectory.append(
            TrajectoryRecord(
                step=len(self.trajectory),
                phase=phase,
                energy=float(energy),
                elapsed_s=self.elapsed(),
                parameters=np.asarray(parameters, dtype=float).copy(),
            )
        )


class Optimizer:
    method_key: str = ""
    display_name: str = ""

    def run(self, context: ExperimentContext) -> None:
        raise NotImplementedError


class LegacyGradientOptimizer(Optimizer):
    core_function_name: str = ""
    uses_qnspsa_state: bool = False

    def _initial_energy(self, context: ExperimentContext) -> float:
        return float(
            context.core_module.Transverse_Ising_Measurement(
                context.hamiltonian,
                bind_ansatz(context.ansatz, context.initial_point),
                context.config.shots,
                context.sampler,
            )
        )

    def run(self, context: ExperimentContext) -> None:
        initial_energy = self._initial_energy(context)
        context.record("initial", initial_energy, context.initial_point)

        def callback(
            parameters: Iterable[float],
            energy: float,
            fubini_matrix_previous: Optional[np.ndarray] = None,
        ) -> None:
            context.record("update", float(energy), parameters)
            if (
                context.config.save_fubini_history
                and fubini_matrix_previous is not None
            ):
                context.fubini_history.append(
                    np.asarray(fubini_matrix_previous, dtype=float)
                )

        function = getattr(context.core_module, self.core_function_name)
        if self.uses_qnspsa_state:
            previous_fubini_matrix = np.eye(context.ansatz.num_parameters, dtype=float)
            last_n_steps = np.zeros(context.config.history_length, dtype=float)
            last_n_steps[0] = initial_energy
            function(
                context.hamiltonian,
                context.initial_point,
                context.config.learning_rate,
                context.ansatz,
                context.config.iterations,
                0,
                context.config.shots,
                callback,
                context.sampler,
                previous_fubini_matrix,
                last_n_steps,
            )
            return

        function(
            context.hamiltonian,
            context.initial_point,
            context.config.learning_rate,
            context.ansatz,
            context.config.iterations,
            context.config.shots,
            callback,
            context.sampler,
        )


class PSROptimizer(LegacyGradientOptimizer):
    method_key = "psr"
    display_name = "PSR"
    core_function_name = "Customize_Parameter_Shift_Rule"


class FDOptimizer(LegacyGradientOptimizer):
    method_key = "fd"
    display_name = "FD"
    core_function_name = "Customize_Finite_Difference"


class SPSAOptimizer(LegacyGradientOptimizer):
    method_key = "spsa"
    display_name = "SPSA"
    core_function_name = "Customize_SPSA"


class QNBDAPSROptimizer(LegacyGradientOptimizer):
    method_key = "qnbda_psr"
    display_name = "QN-BDA+PSR"
    core_function_name = "Customize_Quantum_Natural_Gradient_Descent"


class QNBDASPSAOptimizer(LegacyGradientOptimizer):
    method_key = "qnbda_spsa"
    display_name = "QN-BDA+SPSA"
    core_function_name = "Customize_QN_SPSA_blocking"


class QNSPSAPSROptimizer(LegacyGradientOptimizer):
    method_key = "qnspsa_psr"
    display_name = "QN-SPSA+PSR"
    core_function_name = "Customize_QNSPSA_PRS_blocking"
    uses_qnspsa_state = True


class QNSPSASPSAOptimizer(LegacyGradientOptimizer):
    method_key = "qnspsa_spsa"
    display_name = "QN-SPSA+SPSA"
    core_function_name = "Customize_QNSPSA_SPSA_blocking"
    uses_qnspsa_state = True


class QNSPSAPSRMonteCarloOptimizer(LegacyGradientOptimizer):
    method_key = "qnspsa_psr_mc"
    display_name = "QN-SPSA+PSR-MC"
    core_function_name = "Customize_QNSPSA_PRS_blocking_MonteCarlo"
    uses_qnspsa_state = True


class QNSPSASPSAMonteCarloOptimizer(LegacyGradientOptimizer):
    method_key = "qnspsa_spsa_mc"
    display_name = "QN-SPSA+SPSA-MC"
    core_function_name = "Customize_QNSPSA_SPSA_MonteCarlo"
    uses_qnspsa_state = True


class COBYLAOptimizer(Optimizer):
    method_key = "cobyla"
    display_name = "COBYLA"

    def run(self, context: ExperimentContext) -> None:
        try:
            from qiskit.algorithms.optimizers import COBYLA
        except Exception as exc:
            raise RuntimeError("Qiskit COBYLA is unavailable in this environment") from exc

        optimizer = COBYLA(maxiter=context.config.iterations)

        def objective(parameters: Iterable[float]) -> float:
            params = np.asarray(parameters, dtype=float)
            energy = float(
                context.core_module.Transverse_Ising_Measurement(
                    context.hamiltonian,
                    bind_ansatz(context.ansatz, params),
                    context.config.shots,
                    context.sampler,
                )
            )
            context.record("objective", energy, params)
            return energy

        optimizer.minimize(fun=objective, x0=context.initial_point)


OPTIMIZER_CLASSES: List[Type[Optimizer]] = [
    PSROptimizer,
    FDOptimizer,
    SPSAOptimizer,
    QNBDAPSROptimizer,
    QNBDASPSAOptimizer,
    QNSPSAPSROptimizer,
    QNSPSASPSAOptimizer,
    QNSPSAPSRMonteCarloOptimizer,
    QNSPSASPSAMonteCarloOptimizer,
    COBYLAOptimizer,
]

OPTIMIZER_REGISTRY: Dict[str, Type[Optimizer]] = {
    optimizer.method_key: optimizer for optimizer in OPTIMIZER_CLASSES
}

METHOD_DISPLAY_NAMES: Dict[str, str] = {
    optimizer.method_key: optimizer.display_name for optimizer in OPTIMIZER_CLASSES
}

METHOD_ALIASES: Dict[str, str] = {
    "parameter_shift_rule": "psr",
    "parameter-shift-rule": "psr",
    "constant_psr": "psr",
    "constant+psr": "psr",
    "finite_difference": "fd",
    "finite-difference": "fd",
    "constant_fd": "fd",
    "constant+fd": "fd",
    "constant_spsa": "spsa",
    "constant+spsa": "spsa",
    "qn_bda_psr": "qnbda_psr",
    "qn-bda+psr": "qnbda_psr",
    "qnbda+psr": "qnbda_psr",
    "qng": "qnbda_psr",
    "qng_psr": "qnbda_psr",
    "qn_bda_spsa": "qnbda_spsa",
    "qn-bda+spsa": "qnbda_spsa",
    "qnbda+spsa": "qnbda_spsa",
    "qn_spsa_psr": "qnspsa_psr",
    "qn-spsa+psr": "qnspsa_psr",
    "qnspsa+psr": "qnspsa_psr",
    "qn_spsa_spsa": "qnspsa_spsa",
    "qn-spsa+spsa": "qnspsa_spsa",
    "qnspsa+spsa": "qnspsa_spsa",
}


def normalize_method(name: str) -> str:
    cleaned = name.strip().lower().replace("-", "_").replace("+", "_")
    cleaned = cleaned.replace(" ", "_")
    if cleaned in METHOD_DISPLAY_NAMES:
        return cleaned
    if cleaned in METHOD_ALIASES:
        return METHOD_ALIASES[cleaned]
    original_alias = name.strip().lower()
    if original_alias in METHOD_ALIASES:
        return METHOD_ALIASES[original_alias]
    raise ValueError(
        f"Unknown method '{name}'. Available methods: {', '.join(METHOD_DISPLAY_NAMES)}"
    )


def parse_methods(raw_methods: Optional[List[str]]) -> List[str]:
    if not raw_methods:
        return ["qnspsa_psr"]

    expanded: List[str] = []
    for item in raw_methods:
        for token in item.split(","):
            token = token.strip()
            if not token:
                continue
            if token in {"all", "all-current", "all_current"}:
                expanded.extend(METHOD_DISPLAY_NAMES.keys())
            else:
                expanded.append(normalize_method(token))

    seen = set()
    ordered: List[str] = []
    for method in expanded:
        if method not in seen:
            ordered.append(method)
            seen.add(method)
    return ordered


def summarize_records(
    records: List[TrajectoryRecord],
    exact_energy: Optional[float],
    cost: CostBreakdown,
) -> Dict[str, Any]:
    energies = np.asarray([record.energy for record in records], dtype=float)
    best_idx = int(np.argmin(energies))
    final_energy = float(energies[-1])
    best_energy = float(energies[best_idx])

    summary: Dict[str, Any] = {
        "num_records": len(records),
        "final_energy": final_energy,
        "best_energy": best_energy,
        "best_step": int(records[best_idx].step),
        "exact_ground_energy": exact_energy,
        "nominal_total_quantum_evals": None,
    }

    if cost.total_evals_per_update is not None:
        charged_records = sum(
            1 for record in records if record.phase in {"update", "objective"}
        )
        summary["nominal_total_quantum_evals"] = int(
            cost.total_evals_per_update * charged_records
        )

    if exact_energy is not None and abs(exact_energy) > 0:
        summary["final_relative_error"] = float(
            abs((final_energy - exact_energy) / exact_energy)
        )
        summary["best_relative_error"] = float(
            abs((best_energy - exact_energy) / exact_energy)
        )
        summary["final_absolute_error"] = float(abs(final_energy - exact_energy))
        summary["best_absolute_error"] = float(abs(best_energy - exact_energy))
    else:
        summary["final_relative_error"] = None
        summary["best_relative_error"] = None
        summary["final_absolute_error"] = None
        summary["best_absolute_error"] = None

    return summary


def build_metadata(
    config: ExperimentConfig,
    optimizer: Optimizer,
    ansatz: Any,
    exact_energy: Optional[float],
    cost: CostBreakdown,
    qiskit_module: Any,
) -> Dict[str, Any]:
    return {
        "method": config.method,
        "method_display_name": optimizer.display_name,
        "num_qubits": config.num_qubits,
        "j": config.j,
        "h": config.h,
        "ansatz": config.ansatz,
        "ansatz_name": ansatz.name,
        "ansatz_num_parameters": int(ansatz.num_parameters),
        "ansatz_num_qubits": int(ansatz.num_qubits),
        "reps": config.reps,
        "entanglement": config.entanglement,
        "iterations": config.iterations,
        "learning_rate": config.learning_rate,
        "shots": config.shots,
        "seed": config.seed,
        "initial_mode": config.initial_mode,
        "initial_value": config.initial_value,
        "initial_scale": config.initial_scale,
        "exact_ground_energy": exact_energy,
        "nominal_quantum_eval_cost": cost.to_dict(),
        "git_commit": git_commit(),
        "qiskit_version": getattr(qiskit_module, "__version__", None),
        "command": config.command,
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }


def run_experiment(config: ExperimentConfig) -> RunResult:
    import CoreVQEModified as core_module
    import qiskit
    from qiskit.primitives import Sampler

    run_dir = config.run_dir()
    if run_dir.exists() and not config.overwrite:
        raise FileExistsError(
            f"{run_dir} already exists; pass --overwrite or choose --run-label"
        )

    np.random.seed(config.seed)
    random.seed(config.seed)

    optimizer = OPTIMIZER_REGISTRY[config.method]()
    sampler = Sampler()
    hamiltonian = core_module.Ising_hamiltonian(config.num_qubits, config.j, config.h)
    ansatz = build_ansatz(
        config.ansatz,
        config.num_qubits,
        config.entanglement,
        config.reps,
    )
    initial_point = initial_parameters(config, ansatz.num_parameters)
    exact_energy = exact_ground_energy(
        hamiltonian,
        config.num_qubits,
        config.exact_max_qubits,
    )
    cost = CostModel.for_method(core_module, config.method, ansatz)
    metadata = build_metadata(config, optimizer, ansatz, exact_energy, cost, qiskit)

    trajectory: List[TrajectoryRecord] = []
    fubini_history: List[np.ndarray] = []
    started = time.time()
    context = ExperimentContext(
        config=config,
        core_module=core_module,
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        initial_point=initial_point,
        sampler=sampler,
        trajectory=trajectory,
        fubini_history=fubini_history,
        started=started,
    )

    optimizer.run(context)

    elapsed = time.time() - started
    metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    metadata["elapsed_s"] = elapsed
    summary = summarize_records(trajectory, exact_energy, cost)
    summary["elapsed_s"] = elapsed

    result = RunResult(
        run_dir=run_dir,
        metadata=metadata,
        trajectory=trajectory,
        summary=summary,
        fubini_history=fubini_history,
    )
    ResultWriter().write(result)
    return result


AGGREGATE_FIELDS = [
    "run_dir",
    "method",
    "method_display_name",
    "num_qubits",
    "h",
    "j",
    "ansatz",
    "reps",
    "entanglement",
    "iterations",
    "learning_rate",
    "shots",
    "seed",
    "final_energy",
    "best_energy",
    "final_absolute_error",
    "best_absolute_error",
    "final_relative_error",
    "best_relative_error",
    "best_step",
    "exact_ground_energy",
    "nominal_total_quantum_evals",
    "elapsed_s",
    "git_commit",
    "started_utc",
    "finished_utc",
]


def aggregate_result_dirs(roots: Sequence[Path]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    metadata_paths = []
    for root in roots:
        if (root / "metadata.json").exists():
            metadata_paths.append(root / "metadata.json")
        else:
            metadata_paths.extend(sorted(root.rglob("metadata.json")))

    for metadata_path in sorted(set(metadata_paths)):
        run_dir = metadata_path.parent
        summary_path = run_dir / "summary.json"
        if not summary_path.exists():
            continue
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        row = {
            "run_dir": str(run_dir),
            "method": metadata.get("method"),
            "method_display_name": metadata.get("method_display_name"),
            "num_qubits": metadata.get("num_qubits"),
            "h": metadata.get("h"),
            "j": metadata.get("j"),
            "ansatz": metadata.get("ansatz"),
            "reps": metadata.get("reps"),
            "entanglement": metadata.get("entanglement"),
            "iterations": metadata.get("iterations"),
            "learning_rate": metadata.get("learning_rate"),
            "shots": metadata.get("shots"),
            "seed": metadata.get("seed"),
            "final_energy": summary.get("final_energy"),
            "best_energy": summary.get("best_energy"),
            "final_absolute_error": summary.get("final_absolute_error"),
            "best_absolute_error": summary.get("best_absolute_error"),
            "final_relative_error": summary.get("final_relative_error"),
            "best_relative_error": summary.get("best_relative_error"),
            "best_step": summary.get("best_step"),
            "exact_ground_energy": summary.get("exact_ground_energy"),
            "nominal_total_quantum_evals": summary.get(
                "nominal_total_quantum_evals"
            ),
            "elapsed_s": summary.get("elapsed_s"),
            "git_commit": metadata.get("git_commit"),
            "started_utc": metadata.get("started_utc"),
            "finished_utc": metadata.get("finished_utc"),
        }
        rows.append(row)
    return rows


def write_aggregate_csv(rows: Sequence[Dict[str, Any]], output: Optional[Path]) -> None:
    if output is None:
        handle = sys.stdout
        close_handle = False
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        handle = output.open("w", encoding="utf-8", newline="")
        close_handle = True

    try:
        writer = csv.DictWriter(handle, fieldnames=AGGREGATE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in AGGREGATE_FIELDS})
    finally:
        if close_handle:
            handle.close()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run clean manuscript-revision VQE sweeps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--method",
        action="append",
        help=(
            "Optimizer method. Repeat or comma-separate values. Use all-current "
            "for the methods currently implemented by the source repo."
        ),
    )
    parser.add_argument("--num-qubits", "--n", dest="num_qubits", type=int, required=True)
    parser.add_argument("--h", type=float, required=True, help="Transverse field strength")
    parser.add_argument("--j", type=float, default=1.0, help="Ising coupling strength")
    parser.add_argument(
        "--ansatz",
        choices=["RealAmplitudes", "EfficientSU2"],
        default="RealAmplitudes",
    )
    parser.add_argument("--reps", type=int, default=1)
    parser.add_argument("--entanglement", default="reverse_linear")
    parser.add_argument("--iterations", type=int, default=500)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--shots", type=parse_shots, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--initial-mode",
        choices=["constant", "uniform"],
        default="constant",
    )
    parser.add_argument("--initial-value", type=float, default=-0.5)
    parser.add_argument("--initial-scale", type=float, default=0.5)
    parser.add_argument("--initial-point-file")
    parser.add_argument(
        "--initial-point-delimiter",
        default=None,
        help="Delimiter for text initial-point files; default accepts whitespace.",
    )
    parser.add_argument("--history-length", type=int, default=5)
    parser.add_argument("--exact-max-qubits", type=int, default=12)
    parser.add_argument("--save-fubini-history", action="store_true")
    parser.add_argument("--out", default="results/revision")
    parser.add_argument("--run-label")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned run directories without executing optimizers.",
    )
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.iterations <= 0:
        parser.error("--iterations must be positive")
    if args.num_qubits <= 0:
        parser.error("--num-qubits must be positive")
    if args.reps < 0:
        parser.error("--reps must be non-negative")
    if args.history_length <= 0:
        parser.error("--history-length must be positive")


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    methods = parse_methods(args.method)
    command = [Path(sys.argv[0]).name] + (argv if argv is not None else sys.argv[1:])

    configs = [ExperimentConfig.from_args(args, method, command) for method in methods]
    if args.dry_run:
        for config in configs:
            print(config.run_dir())
        return 0

    completed: List[Path] = []
    for config in configs:
        completed.append(run_experiment(config).run_dir)

    for path in completed:
        print(path)
    return 0
