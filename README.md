# Variational Quantum Eigensolver: A Comparative Analysis of Classical and Quantum Optimization Methods
<strong>Abtract: </strong><br>
In this study https://arxiv.org/abs/2412.19176, we delved into several optimization methods, both classical and quantum, and ana-
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
- All of the main source codes are contained within the CoreVQEModified.py file, which is used and executed by the run_VQE_modified_on_HPC_sample.py script to generate the dataset. run_VQE_modified_on_HPC_sample_v3.py is a modified version (assisted by Claude for better parallelization) of CoreVQEModified.py that is faster implementations on VMs/clusters. 
- All data visualizations were created by the Plot_data.ipynb notebook. <br>
- A list of required versions are provided in the Requirements.txt file.

## Figure 3 and Figure 4 campaigns

The revised Figure 3 and Figure 4 data are generated from the native legacy
trajectory format under `energy/`. The existing legacy rows remain in `energy/`; the
new rows added for the complete 3 x 3 optimizer grid are PSR, QN-BDA+FD, and
QN-SPSA+FD for both reverse-linear and full entanglement. The file
`figures/figure3_4_new_rows_manifest_20260603.txt` lists the exact raw
trajectory files added for those rows.

The launchers used on the VM were:

- `run_campaign.py`: reverse-linear entanglement, using
  `run_VQE_modified_on_HPC_sample_v3.py`.
- `run_campaign_full.py`: full entanglement, using
  `run_VQE_modified_on_HPC_sample_full.py`.

Both launchers write directly to `energy/` using the same filename convention as
the legacy data. The `--list` option prints the planned output filenames without
running the campaign.

Figure 3 uses the 12-qubit, reps=1, h=2.0, 500-iteration reverse-linear data.
The command used to generate the three new raw-data rows was:

```bash
python run_campaign.py \
  --methods psr,qnbda_fd,qnspsa_fd \
  --reps 1 \
  --iterations 500 \
  --jobs 7 \
  --grad-workers 36 \
  --h 2.0
```

Figure 4 compares the same 12-qubit, reps=1, h=2.0, 500-iteration trajectories
for reverse-linear and full entanglement. The reverse-linear data are generated
with the same command above. The full-entanglement rows were generated with:

```bash
python run_campaign_full.py \
  --methods psr,qnbda_fd,qnspsa_fd \
  --reps 1 \
  --iterations 500 \
  --jobs 7 \
  --grad-workers 36 \
  --h 2.0
```

These commands regenerate the raw trajectories behind the paper Figure 3 and
Figure 4 panels. The manuscript plotting scripts used for final figure styling
are kept with the paper source and are not committed in this repository.


Useful method keys are `cobyla`, `spsa`, `fd`, `psr`, `qnbda_spsa`,
`qnbda_fd`, `qnbda_psr`, `qnspsa_spsa`, `qnspsa_fd`, and `qnspsa_psr`.
Deterministic methods run one seed by default; stochastic methods run seven
seeds. This default is the setting used for the Figure 3/4 consistency check:
PSR and QN-BDA+FD each produce one raw trajectory, while QN-SPSA+FD produces
seven stochastic trajectories.

For cluster use, choose `--jobs` and `--grad-workers` so that
`jobs * grad-workers` fits the available cores. Re-running a campaign is safe:
the harnesses reuse the native checkpoint/resume behavior and overwrite only the
corresponding method/seed trajectory files.
