<strong>Abtract</strong>: <br>
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

<strong>Instruction for users: <\strong><br>
All of the main source codes are contained within the CoreVQEModified.py file, which is used and executed by the run_VQE_modified_on_HPC_sample.py script to generate the dataset. Finally, all data visualizations were created by the Plot_data.ipynb notebook. <br>
A list of required versions are provided in the Requirements.txt file.
