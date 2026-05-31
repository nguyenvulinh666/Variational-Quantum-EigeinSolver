# Variational Quantum Eigensolver: A Comparative Analysis of Classical and Quantum Optimization Methods
<strong>Abtract: </strong><br>
In this study, we delved into several optimization methods, both classical and quantum, and ana-
lyzed the quantum advantage that each of these methods offered, and then we proposed a new combinatorial
optimization scheme, deemed as QN-SPSA+PSR which combines calculating approximately Fubini-study metric
(QN-SPSA) and the exact evaluation of gradient by Parameter-Shift Rule (PSR). The QN-SPSA+PSR method
integrates the QN-SPSA computational efficiency with the precise gradient computation of the PSR, improving
both stability and convergence speed while maintaining low computational consumption. Our results provide
a new potential quantum supremacy in the VQE’s optimization subroutine and enhance viable paths toward
efficient quantum simulations on Noisy Intermediate-Scale Quantum Computing (NISQ) devices. Additionally,
we also conducted a detailed study of quantum circuit ansatz structures in order to find the one that would work
best with the Ising model and NISQ, in which we utilized the symmetry of the investigated model. <br>

<strong>Instruction for users: </strong> <br>
- All of the main source codes are contained within the CoreVQEModified.py file, which is used and executed by the run_VQE_modified_on_HPC_sample.py script to generate the dataset.
- All data visualizations were created by the Plot_data.ipynb notebook. <br>
- A list of required versions are provided in the Requirements.txt file.

## Revision sweep runner

Use `run_revision_sweep.py` for new manuscript-revision experiments. It is a
thin CLI wrapper around `revision_experiment.py`, which contains the structured
experiment config, optimizer adapters, cost model, result writer, and aggregate
helpers. The legacy optimizer implementations remain in `CoreVQEModified.py`.
New runs stay separate from the legacy `energy/`, `parameter/`, and
`fubini_matrix_previous/` text dumps by writing one self-contained directory per
configuration under `results/revision/`.

Example single run:

```bash
python run_revision_sweep.py \
  --method qnspsa_psr \
  --n 12 \
  --h 2.0 \
  --ansatz RealAmplitudes \
  --reps 2 \
  --entanglement reverse_linear \
  --iterations 500 \
  --seed 17 \
  --shots none \
  --out results/revision
```

Useful method names are `psr`, `fd`, `spsa`, `qnbda_psr`, `qnbda_spsa`,
`qnspsa_psr`, `qnspsa_spsa`, `qnspsa_psr_mc`, `qnspsa_spsa_mc`, and `cobyla`.
Use `--method all-current` to run every method currently exposed by the source
repo. Repeat `--method` or comma-separate methods to run a smaller set.

Each run directory contains:

- `metadata.json`: configuration, ansatz size, random seed, Qiskit version, git
  commit, and nominal quantum-evaluation cost per optimizer update.
- `trajectory.csv`: step, phase, energy, relative error when exact diagonalizing
  is feasible, elapsed time, and the parameter vector as JSON.
- `parameters.npy`: dense parameter-history array for post-processing.
- `summary.json`: final/best energy, final/best absolute and relative errors,
  and nominal total quantum-evaluation count.

For cluster use, prefer one scheduler task per `(method, n, h, seed)` rather
than wrapping a large Python `ProcessPoolExecutor` around the sweep. The
optimizer implementations can parallelize PSR/FD gradient calls internally; set
`VQE_WORKERS=32` or another scheduler-appropriate value to cap per-run workers.
Use `VQE_PARALLEL_BACKEND=serial` for debugging and local smoke tests.

Dry-run example:

```bash
python run_revision_sweep.py --method qnspsa_psr,cobyla --n 12 --h 2.0 \
  --run-label smoke --dry-run
```

Aggregate completed runs into a table-ready CSV:

```bash
python aggregate_revision_results.py results/revision \
  --out results/revision_summary.csv
```

Generate cluster command manifests:

```bash
python generate_revision_manifest.py --preset smoke \
  --label revision_smoke_20260531 \
  --out results/revision_smoke_20260531_commands.txt

python generate_revision_manifest.py --preset full \
  --label revision_full_20260531 \
  --out results/revision_full_20260531_commands.txt
```

Create seed-aggregated tables and plots:

```bash
python plot_revision_results.py results/revision/revision_full_20260531 \
  --out-dir results/revision_full_20260531_plots
```
