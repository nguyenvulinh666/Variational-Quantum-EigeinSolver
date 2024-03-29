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
from CoreVQEModified import Ising_hamiltonian
# import ansatz 
from CoreVQEModified import Customize_RealAmplidues, Customize_EfficientSU2
# import optimize
from CoreVQEModified import Customize_Finite_Difference, Customize_Parameter_Shift_Rule, Customize_Quantum_Natural_Gradient_Descent, Customize_SPSA, Customize_QNSPSA_PRS_blocking, Customize_QN_SPSA_blocking, Customize_QNSPSA_SPSA_blocking, Customize_HESPSA_SPSA_blocking, Customize_HE_PRS

from qiskit.primitives import Sampler

# Import Cobyla from Qiskit
import numpy as np
from qiskit.algorithms.optimizers import COBYLA
from qiskit.providers.basicaer import StatevectorSimulatorPy  # local simulator
from qiskit.algorithms import VQE


def main():
    num_qubits = 7
    reps = 1
    J = 1
    h = 2
    sampler = Sampler()
    
    
    intermediate_info = {
        'parameters': [],
        'energy': [],
    }


    hamiltonian = Ising_hamiltonian(num_qubits, J, h)
    
    ansatz = Customize_RealAmplidues(num_qubits, reps)
    ansatz_name = 'RealAmplidues'

    #ansatz = Customize_EfficientSU2(num_qubits, reps)
    #ansatz_name = 'EffcientSU2'

    interation = 500
    initial_point = np.zeros(ansatz.num_parameters) - 0.5

    eta = 0.01
    shots = None

    number_of_times = 7
    
    

    for i in range(number_of_times):
        seed = i    
        #optimize = Customize_Finite_Difference
        #optimize = Customize_Parameter_Shift_Rule
        #optimize = Customize_Quantum_Natural_Gradient_Descent
        #optimize = Customize_SPSA
        optimize = Customize_QNSPSA_PRS_blocking
        #optimize = Customize_QN_SPSA_blocking
        #optimize = Customize_QNSPSA_SPSA_blocking
        # optimize = Customize_HESPSA_SPSA_blocking
        # optimize = Customize_HE_PRS

        #optimize = COBYLA
        #optimize_run = COBYLA(maxiter=500)

        file_name_parameters = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz_name}({num_qubits},{reps}) - J{J}h{h} - seed{seed} - Parameters - {i+1}'
        file_name_energy = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz_name}({num_qubits},{reps}) - J{J}h{h} - seed{seed} - Energy - {i+1}'

        file_parameters = open(f'{file_name_parameters}.txt', 'a+')
        file_energy = open(f'{file_name_energy}.txt', 'a+')

        file_parameters.write('Parameters \n \n')
        file_energy.write('Energy \n \n')

        file_parameters.close()
        file_energy.close()

        def callback(parameters, energy):
        #def callback(nfev, parameters, energy, stddev):
            #intermediate_info['parameters'].append(parameters)
            #intermediate_info['energy'].append(energy)
            file_parameters = open(f'{file_name_parameters}.txt', 'a+')
            file_energy = open(f'{file_name_energy}.txt', 'a+')

            file_parameters.write(f'{str(parameters)} \n')
            file_energy.write(f'{str(energy)} \n')

            file_parameters.close()
            file_energy.close()
        
        np.random.seed(seed=seed)
        start_time = time.time()
        
        #energy = Customize_Finite_Difference(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)
        #energy = Customize_Parameter_Shift_Rule(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)       
        #energy = Customize_Quantum_Natural_Gradient_Descent(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)     
        #energy = Customize_SPSA(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)
        energy = Customize_QNSPSA_PRS_blocking(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)        
        #energy = Customize_QN_SPSA_blocking(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)
        #energy = Customize_QNSPSA_SPSA_blocking(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)      
        # energy = Customize_HESPSA_SPSA_blocking(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)
        # energy = Customize_HE_PRS(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)

        #local_vqe = VQE(ansatz=ansatz,
        #        optimizer=optimize_run,
        #        initial_point=initial_point,
        #        quantum_instance=StatevectorSimulatorPy(),
        #        callback=callback)

        #energy = local_vqe.compute_minimum_eigenvalue(hamiltonian)
        
        end_time = time.time()


        file_parameters = open(f'{file_name_parameters}.txt', 'a+')
        file_energy = open(f'{file_name_energy}.txt', 'a+')

        file_parameters.write(f'\n Time \n')
        file_energy.write(f'\n Time \n')

        file_parameters.write(f'{str(end_time-start_time)} \n')
        file_energy.write(f'{str(end_time-start_time)} \n')

        file_parameters.close()
        file_energy.close()


if __name__ == '__main__':
    main()