"""Campaign launcher: fan out the Table-2 (method, seed) jobs across cores / VMs.

Each (method, seed) is an independent job that reuses the native harness'
checkpoint/resume logic, so re-running is safe and resumes partial work.

Examples
--------
List the full plan (runs nothing):
    python run_campaign.py --list

Run everything on this VM, 6 jobs at a time:
    python run_campaign.py --jobs 6

Split across 4 VMs (run this on each, changing k = 0,1,2,3):
    python run_campaign.py --shard 0/4 --jobs 6
    python run_campaign.py --shard 1/4 --jobs 6
    ...

Run only the cheap stochastic methods, many at once:
    python run_campaign.py --methods spsa,qnbda_spsa,qnspsa_spsa --jobs 24

Notes
-----
Each custom-optimizer job runs in its own python subprocess; inside it the
optimizer may spawn an inner ProcessPoolExecutor of up to p (=36) workers for
gradient methods. Pick --jobs so jobs * 36 stays within your core count
(e.g. 255 cores -> --jobs 6 for gradient-heavy methods, higher for SPSA-based
methods which use no inner pool).
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HARNESS_MODULE = "run_VQE_modified_on_HPC_sample_full"  # full-entanglement harness

# method key -> (CoreVQE function name, n_seeds, harness)
JOBS_SPEC = {
    # deterministic: 1 seed
    "qnbda_psr":   ("Customize_Quantum_Natural_Gradient_Descent", 1, "hpc"),
    "fd":          ("Customize_Finite_Difference",                 1, "hpc"),
    "qnbda_fd":    ("Customize_QN_FD",                             1, "hpc"),
    "psr":         ("Customize_Parameter_Shift_Rule",             1, "hpc"),
    # stochastic: 7 seeds
    "spsa":        ("Customize_SPSA",                              7, "hpc"),
    "qnbda_spsa":  ("Customize_QN_SPSA_blocking",                 7, "hpc"),
    "qnspsa_psr":  ("Customize_QNSPSA_PRS_blocking",              7, "hpc"),
    "qnspsa_fd":   ("Customize_QNSPSA_FD_blocking",               7, "hpc"),
    "qnspsa_spsa": ("Customize_QNSPSA_SPSA_blocking",             7, "hpc"),
    # COBYLA: separate script, single run
    "cobyla":      ("COBYLA",                                      1, "cobyla"),
}

METHOD_ORDER = list(JOBS_SPEC.keys())


def build_jobs(methods, n, h, reps, iterations, seeds_override=None):
    jobs = []
    for key in methods:
        func, n_seeds, harness = JOBS_SPEC[key]
        if seeds_override is not None:
            n_seeds = seeds_override
        for seed in range(n_seeds):
            jobs.append({"key": key, "func": func, "seed": seed, "harness": harness,
                         "n": n, "h": h, "reps": reps, "iterations": iterations})
    return jobs


def energy_filename(job):
    # mirrors the harness naming for --list / verification
    n, h, reps, it = job["n"], job["h"], job["reps"], job["iterations"]
    return (f"{job['func']} - LR 0.01 - shots None - interation {it} - "
            f"RealAmplitudes({n},{reps}) - full - J1h{h} - Energy - {job['seed']+1}.txt")


def job_command(job):
    if job["harness"] == "cobyla":
        # NOTE: COBYLA reps/entanglement are set inside VQE_qiskit.py, not by --reps.
        return [sys.executable, "-c", "from VQE_qiskit import main; main()"]
    code = (
        f"import {HARNESS_MODULE} as r; "
        f"r.main(({job['n']}, {float(job['h'])}, r.{job['func']}, {job['seed']}, {job['reps']}))"
    )
    return [sys.executable, "-c", code]


def run_job(job):
    start = time.time()
    cmd = job_command(job)
    env = dict(
        os.environ,
        VQE_GRAD_WORKERS=str(job["grad_workers"]),
        VQE_ITERATIONS=str(job["iterations"]),
        VQE_REPS=str(job["reps"]),  # used by VQE_qiskit (COBYLA); harness uses the tuple
    )
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return job, proc.returncode, time.time() - start, proc.stderr[-500:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methods", default="all",
                    help="comma-separated method keys, or 'all'")
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--reps", type=int, default=2, help="ansatz reps (1=legacy table, 2=reps-2 table)")
    ap.add_argument("--iterations", type=int, default=500, help="optimizer iterations per run")
    ap.add_argument("--seeds", type=int, default=None,
                    help="override seed count for ALL methods (e.g. 7 random restarts)")
    ap.add_argument("--random-init", action="store_true",
                    help="random per-seed init (uniform[-1,1]) instead of constant -0.5")
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--jobs", type=int, default=6, help="concurrent jobs on this VM")
    ap.add_argument("--grad-workers", default="auto",
                    help="inner gradient workers per job: 'auto' (cores/jobs), or an int "
                         "(1 = serial inner). Total cores used ~= jobs * grad-workers.")
    ap.add_argument("--shard", default=None,
                    help="k/N: run only jobs where index %% N == k (for multi-VM)")
    ap.add_argument("--list", action="store_true", help="print plan and exit")
    args = ap.parse_args()

    methods = METHOD_ORDER if args.methods == "all" else args.methods.split(",")
    for m in methods:
        if m not in JOBS_SPEC:
            ap.error(f"unknown method '{m}'. valid: {', '.join(METHOD_ORDER)}")

    if args.random_init:
        os.environ["VQE_RANDOM_INIT"] = "1"  # inherited by job subprocesses
    jobs = build_jobs(methods, args.n, args.h, args.reps, args.iterations, args.seeds)

    if args.shard:
        k, total = (int(x) for x in args.shard.split("/"))
        jobs = [j for i, j in enumerate(jobs) if i % total == k]
        shard_note = f" (shard {k}/{total})"
    else:
        shard_note = ""

    # Lever #1: size the inner gradient pool so (concurrent jobs) x (workers) ~ cores.
    concurrent_jobs = max(1, min(args.jobs, len(jobs)))
    if args.grad_workers == "auto":
        grad_workers = max(1, (os.cpu_count() or 1) // concurrent_jobs)
    else:
        grad_workers = max(1, int(args.grad_workers))
    for j in jobs:
        j["grad_workers"] = grad_workers

    print(f"Campaign: {len(jobs)} jobs{shard_note}, --jobs={args.jobs}, "
          f"grad_workers={grad_workers} (~{concurrent_jobs*grad_workers} cores), "
          f"n={args.n}, reps={args.reps}, h={args.h}")
    for j in jobs:
        print(f"  [{j['key']:<12} seed {j['seed']}]  -> energy/{energy_filename(j)}")

    if args.list:
        return 0

    print("\nLaunching...\n")
    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_job, j) for j in jobs]
        for fut in as_completed(futures):
            job, rc, elapsed, err = fut.result()
            status = "ok" if rc == 0 else f"FAIL(rc={rc})"
            print(f"  {status:<12} {job['key']:<12} seed {job['seed']}  {elapsed:.0f}s")
            if rc != 0:
                failures.append((job, err))

    if failures:
        print(f"\n{len(failures)} job(s) failed:")
        for job, err in failures:
            print(f"  {job['key']} seed {job['seed']}:\n    {err.strip()}")
        return 1
    print("\nAll jobs complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
