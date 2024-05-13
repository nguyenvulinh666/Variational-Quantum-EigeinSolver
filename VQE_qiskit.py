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
from qiskit.circuit.library import RealAmplitudes, EfficientSU2
from qiskit.primitives import Sampler

# Import Cobyla from Qiskit
import numpy as np
from qiskit.algorithms.optimizers import COBYLA
from qiskit.providers.basicaer import StatevectorSimulatorPy  # local simulator
from qiskit.algorithms import VQE


def main():
    num_qubits = 12
    reps = 1
    J = 1
    h = 2.0
    sampler = Sampler()
    
    
    intermediate_info = {
        'parameters': [],
        'energy': [],
    }

    entanglement = 'reverse_linear'
    
    hamiltonian = Ising_hamiltonian(num_qubits, J, h)
    
    ansatz = RealAmplitudes(num_qubits, entanglement, reps, insert_barriers=True).decompose()
    #ansatz = EfficientSU2(num_qubits, entanglement=entanglement, reps=reps, insert_barriers=True).decompose()
    
    interation = 500
    initial_point = np.zeros(ansatz.num_parameters) - 0.5

    eta = 0.01
    shots = None

    number_of_times = 7
    
    
    for i in range(number_of_times):
        optimize = COBYLA  
        optimize_run = COBYLA(maxiter=interation)

        file_name_parameters = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz.name}({num_qubits},{reps}) - {entanglement} - J{J}h{h} - Parameters - {i+1}'
        file_name_energy = f'{str(optimize.__name__)} - LR {eta} - shots {shots} - interation {interation} - {ansatz.name}({num_qubits},{reps}) - {entanglement} - J{J}h{h} - Energy - {i+1}'
        
        print("bug")
        print(file_name_parameters)
        print("bug")
        file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
        file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')

        file_parameters.write('Parameters \n \n')
        file_energy.write('Energy \n \n')

        file_parameters.close()
        file_energy.close()

        #def callback(parameters, energy):
        def callback(nfev, parameters, energy, stddev):
            #intermediate_info['parameters'].append(parameters)
            #intermediate_info['energy'].append(energy)
            file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
            file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')
            
            print(energy)
            
            file_parameters.write(f'{str(list(parameters))} \n')
            file_energy.write(f'{str(energy)} \n')

            file_parameters.close()
            file_energy.close()
        
        start_time = time.time()
        
        #energy = optimize(hamiltonian, initial_point, eta, ansatz, interation, shots, callback, sampler)

        local_vqe = VQE(ansatz=ansatz,
                optimizer=optimize_run,
                initial_point=initial_point,
                quantum_instance=StatevectorSimulatorPy(),
                callback=callback)
        
        np.random.seed(i)
        
        energy = local_vqe.compute_minimum_eigenvalue(hamiltonian)
        
        end_time = time.time()


        file_parameters = open(f'./parameter/{file_name_parameters}.txt', 'a+')
        file_energy = open(f'./energy/{file_name_energy}.txt', 'a+')

        file_parameters.write(f'\n Time \n')
        file_energy.write(f'\n Time \n')

        file_parameters.write(f'{str(end_time-start_time)} \n')
        file_energy.write(f'{str(end_time-start_time)} \n')

        file_parameters.close()
        file_energy.close()


if __name__ == '__main__':
    main()