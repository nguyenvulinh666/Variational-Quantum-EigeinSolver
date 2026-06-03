import numpy as np
from qiskit.execute_function import execute
from qiskit import BasicAer
backend = BasicAer.get_backend('qasm_simulator')
from qiskit import QuantumCircuit
import qiskit
from qiskit import Aer
from qiskit.primitives import Sampler
from numpy.linalg import eig
from qiskit.opflow import Z, I, X, Y
from qiskit.quantum_info import Statevector, Operator, Pauli
from qiskit.circuit import Parameter
import time



# import model
from CoreVQEModified_v3 import Ising_hamiltonian
# import ansatz 
from qiskit.circuit.library import RealAmplitudes, EfficientSU2
# import optimize
from CoreVQEModified_v3 import Customize_Finite_Difference, Customize_Parameter_Shift_Rule, Customize_Quantum_Natural_Gradient_Descent, Customize_SPSA, Customize_QNSPSA_PRS_blocking, Customize_QN_SPSA_blocking, Customize_QNSPSA_SPSA_blocking, Customize_QNSPSA_PRS_blocking_MonteCarlo, Customize_QNSPSA_SPSA_MonteCarlo, Customize_QNSPSA_FD_blocking, Customize_QN_FD

# import measurement
from CoreVQEModified_v3 import Transverse_Ising_Measurement

from qiskit.primitives import Sampler

# Import Cobyla from Qiskit
import numpy as np
from qiskit.algorithms.optimizers import COBYLA
from qiskit.providers.basicaer import StatevectorSimulatorPy  # local simulator
from qiskit.algorithms import VQE

import os

from qiskit.primitives import Estimator, Sampler, BaseEstimator, BackendEstimator



def main(params):
    num_qubits, h, optimize, run_time, reps = params
    # num_qubits = 7
    J = 1
    # h = 2
    sampler = Sampler()
    
    
    intermediate_info = {
        'parameters': [],
        'energy': [],
    }
    
    
    history_length = 5
    
    entanglement = 'reverse_linear'
    
    hamiltonian = Ising_hamiltonian(num_qubits, J, h)
    
    ansatz = RealAmplitudes(num_qubits, entanglement, reps, insert_barriers=True).decompose()
    
    #ansatz = EfficientSU2(num_qubits, entanglement=entanglement, reps=reps, insert_barriers=True).decompose()
    
    interation = int(os.environ.get("VQE_ITERATIONS", 500))

    eta = 0.01
    shots = None


    file_name_parameters = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz.name}({num_qubits},{reps}) - {entanglement} - J{J}h{h} - Parameters - {run_time+1}'
    file_name_energy = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz.name}({num_qubits},{reps}) - {entanglement} - J{J}h{h} - Energy - {run_time+1}'

    if str(optimize.__name__)[:16] == 'Customize_QNSPSA':
        file_name_fubini_matrix_previous = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz.name}({num_qubits},{reps}) - J{J}h{h} - PreviousFubiniMatrix - {run_time+1}'

    
    # Call back function
    def callback(parameters, energy, fubini_matrix_previous=None):
        if str(optimize.__name__)[:16] == 'Customize_QNSPSA':
            #print(file_name_fubini_matrix_previous)
            file_fubini_matrix_previous = open(f'./fubini_matrix_previous/{file_name_fubini_matrix_previous}.txt', "a+")
            file_fubini_matrix_previous.write(f'{str((fubini_matrix_previous).tolist())} \n')
            file_fubini_matrix_previous.close()
        
        file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
        file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')

        file_parameters.write(f'{str(list(parameters))} \n')
        file_energy.write(f'{str(energy)} \n')

        file_parameters.close()
        file_energy.close()
   
    file_length = 0
    
    # Initialized file
    if os.path.isfile(f'./parameter/{file_name_parameters}.txt'):
        data_file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'r').readlines()
        initial_point = eval(data_file_parameters[-1])
        
    else:
        if os.environ.get("VQE_RANDOM_INIT"):
            # Random restart: per-run_time reproducible random init (uniform[-1,1]),
            # kept near the trainable scale of the legacy -0.5 start to avoid the
            # barren plateau of full-range [-pi,pi] at high depth.
            np.random.seed(run_time)
            initial_point = np.random.uniform(-1.0, 1.0, ansatz.num_parameters)
        else:
            initial_point = np.zeros(ansatz.num_parameters) - 0.5 #0.5
        internal_energy = Transverse_Ising_Measurement(hamiltonian, ansatz.bind_parameters({theta: initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        if callback is not None:
            file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
            file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')  
            
            file_parameters.write('Parameters \n \n')
            file_energy.write('Energy \n \n')
            
            file_parameters.close()
            file_energy.close()
            
            if str(optimize.__name__)[:16] == 'Customize_QNSPSA':
                #previous_fubini_matrix = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
                previous_fubini_matrix = np.eye(ansatz.num_parameters, dtype = float)
                file_fubini_matrix_previous = open(f'./fubini_matrix_previous/{file_name_fubini_matrix_previous}.txt', "a+")
                file_fubini_matrix_previous.write('Fubini-study metric previous \n \n')
                file_fubini_matrix_previous.close()
                callback(initial_point, internal_energy, previous_fubini_matrix)
                
            else: 
                callback(initial_point, internal_energy)
    
  
    data_file_parameter = open(f'./parameter/{file_name_parameters}.txt', 'r').readlines()
    file_length = len(data_file_parameter)

    start_time = time.time()

    # Real interation left
    modified_interation = interation + 2 - file_length
    
    if modified_interation <= 0:
      return
    
    # Optimizer
    if str(optimize.__name__)[:16] != 'Customize_QNSPSA':
        energy = optimize(hamiltonian, initial_point, eta, ansatz, modified_interation, shots, callback, sampler)        
    else:
        if file_length > 3:
            step = file_length - 3
            
            data_file_fubini_study_previous = open(f'./fubini_matrix_previous/{file_name_fubini_matrix_previous}.txt', 'r').readlines()
            data_file_energy = open(f'./energy/{file_name_energy}.txt', 'r').readlines()
            #print((data_file_fubini_study_previous[-1]))
            if step > history_length:
              last_n_steps = np.flip([float(data_file_energy[-i].split(' ')[0]) for i in range(1, history_length+1)])
            else:
              #print(np.flip([(data_file_energy[-i].split(' ')[0]) for i in range(1, step+2)]))
              last_n_steps = np.array(list(np.flip([float(data_file_energy[-i].split(' ')[0]) for i in range(1, step+2)])) + [0]*(history_length- step))
            
            previous_fubini_matrix = np.array(eval(data_file_fubini_study_previous[-1]))
            
            
        else:
            last_n_steps = np.zeros(history_length)
            internal_energy = Transverse_Ising_Measurement(hamiltonian, ansatz.bind_parameters({theta: initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
            last_n_steps[0] = internal_energy
            #previous_fubini_matrix = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
            previous_fubini_matrix = np.eye(ansatz.num_parameters, dtype = float)
            step = 0
        energy = optimize(hamiltonian, initial_point, eta, ansatz, modified_interation, step, shots, callback, sampler, previous_fubini_matrix, last_n_steps)      
    end_time = time.time()


    file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
    file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')

    file_parameters.write(f'\nTime \n')
    file_energy.write(f'\n Time \n')

    file_parameters.write(f'{str(end_time-start_time)} \n')
    file_energy.write(f'{str(end_time-start_time)} \n')

    file_parameters.close()
    file_energy.close()

    if str(optimize.__name__)[:16] == 'Customize_QNSPSA':
        file_fubini_matrix_previous = open(f'./fubini_matrix_previous/{file_name_fubini_matrix_previous}.txt', "a+")
        file_fubini_matrix_previous.write(f'\n Time \n')
        file_fubini_matrix_previous.write(f'{str(end_time-start_time)} \n')
        file_fubini_matrix_previous.close()

    return


if __name__ == '__main__':
    params = []
    num_qubits = [12] #[3,4,5,6,7,8,9,10,11]

    # One optimizer per launch. Pick the method by index and re-run the script
    # once per method (see the index table below). The seed count is chosen
    # automatically: deterministic methods run once, stochastic methods run 7x.
    #
    #   [0] Customize_Quantum_Natural_Gradient_Descent  (QN-BDA+PSR)   deterministic
    #   [1] Customize_QNSPSA_PRS_blocking               (QN-SPSA+PSR)  stochastic
    #   [2] Customize_QN_SPSA_blocking                  (QN-BDA+SPSA)  stochastic
    #   [3] Customize_QNSPSA_SPSA_blocking              (QN-SPSA+SPSA) stochastic
    #   [4] Customize_Finite_Difference                 (FD)           deterministic
    #   [5] Customize_SPSA                              (SPSA)         stochastic
    #   [6] Customize_QNSPSA_FD_blocking                (QN-SPSA+FD)   stochastic
    #   [7] Customize_QN_FD                             (QN-BDA+FD)    deterministic
    #   [8] Customize_Parameter_Shift_Rule              (PSR)          deterministic
    optimizes = [Customize_Quantum_Natural_Gradient_Descent, Customize_QNSPSA_PRS_blocking, Customize_QN_SPSA_blocking, Customize_QNSPSA_SPSA_blocking, Customize_Finite_Difference, Customize_SPSA, Customize_QNSPSA_FD_blocking, Customize_QN_FD, Customize_Parameter_Shift_Rule]

    DETERMINISTIC = {
        Customize_Quantum_Natural_Gradient_Descent,
        Customize_Finite_Difference,
        Customize_QN_FD,
        Customize_Parameter_Shift_Rule,
    }

    # >>> Edit this index to select the method for this launch <<<
    selected_index = 7
    selected = optimizes[selected_index]

    # >>> Edit reps for the target table (1 = legacy reps=1, 2 = reps=2) <<<
    reps = 2

    run_time = [0] if selected in DETERMINISTIC else [0, 1, 2, 3, 4, 5, 6]

    #external_field = np.linspace(0,1.8, 10) + 0.0
    external_field = np.array([2]) + 0.0
    for m in range(len(num_qubits)):
      for i in range(len(external_field)):
          for k in range(len(run_time)):
              params.append((num_qubits[m], np.round(external_field[i],1), selected, run_time[k], reps))
    import concurrent.futures
    executor = concurrent.futures.ProcessPoolExecutor()
    executor.map(main, params)
    #for i in range(len(params)):
    #  main(params[i])