from qiskit.circuit import ParameterVector
import numpy as np
from qiskit.execute_function import execute
from qiskit import BasicAer
backend = BasicAer.get_backend('qasm_simulator')
from qiskit import QuantumCircuit
import qiskit
from qiskit import Aer
# from qiskit.primitives import Sampler
from numpy.linalg import eig
from qiskit.opflow import Z, I, X, Y
from qiskit.quantum_info import Statevector, Operator, Pauli
from qiskit.circuit import Parameter
import multiprocessing 
import time
from qiskit.primitives  import Sampler

from qiskit.primitives import Estimator, Sampler, BaseEstimator, BackendEstimator

# Separate quantum circuit to layer
from qiskit.converters import circuit_to_dag, dag_to_circuit

sampler = Sampler()

def mini_derivate(param):
    theta_position, h, divide, initial_point, operator, ansatz, shots, sampler = param
    plus_parameter = initial_point.copy()
    plus_parameter[theta_position] += h
    minus_parameter = initial_point.copy()
    minus_parameter[theta_position] -= h
    result = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/(2*divide)
    return result

def SwapTest(circ1, circ2):
    """
    circ1: for the parameter circuit of VQD
    circ2: for the parameter circuit of VQE, this circuit have the infromation from the VQE running, or the k-1 of VQD running
    shots: number of shots
    backend: backend for running
    """

    qc = circ1.compose(circ2.inverse())
    qc.measure_all()
    job = sampler.run(qc, shots = None)
    result = (job.result().quasi_dists[0].binary_probabilities())
    
    result = result['0'*qc.num_qubits] if '0'*qc.num_qubits in result.keys() else 0
    
    return result
    
# def SwapTest(circ1, circ2):
#     """
#     circ1: for the parameter circuit of VQD
#     circ2: for the parameter circuit of VQE, this circuit have the infromation from the VQE running, or the k-1 of VQD running
#     shots: number of shots
#     backend: backend for running
#     """

#     # Perpare the circuit
#     num_qubits = circ1.num_qubits
#     circ = QuantumCircuit(2*num_qubits+1) # Plus one for the ancilla qubit

#     circ = circ.compose(circ1, [i for i in range(1, num_qubits+1)])
#     circ = circ.compose(circ2, [i for i in range(num_qubits+1, num_qubits*2+1)])
#     circ.barrier()

#     circ.h(0)

#     for i in range(1, num_qubits+1):
#         circ.cswap(0, i, i+num_qubits)

#     circ.h(0)  

#     overlap_value = 1 - 2*Statevector(circ).probabilities([0])[1]
    
    return overlap_value, circ

def matrix_power(matrix, coeff):
    eigen_value, eigen_vector = eig(matrix)

    diagonalize_matrix = np.linalg.pinv(eigen_vector) @ matrix @ eigen_vector

    diagonalize_matrix = np.diag(diagonalize_matrix)

    # print(len(diagonalize_matrix))
    # print(diagonalize_matrix)

    record_matrix = []

    # print(diagonalize_matrix)

    for i in range(len(eigen_vector)):
        # if diagonalize_matrix[i] == 0:
        #     record_matrix.append(0)
        #     continue
        # print(np.complex(diagonalize_matrix[i])**(0.1)) 
        
        record_matrix.append((diagonalize_matrix[i]+0j)**coeff)

    # print(f'record_matrix: {record_matrix}')

    diagonalize_matrix = np.diag(record_matrix)

    return np.array(eigen_vector @ diagonalize_matrix @ np.linalg.pinv(eigen_vector))

# Ising hamiltonian
def Ising_hamiltonian(num_qubits, J, h):
    hamiltonian = 0

    if num_qubits == 1:
        hamiltonian = -J*Z - h*X
        return hamiltonian.reduce()

    if num_qubits == 2:
        hamiltonian = -J*Z^Z
        hamiltonian -= h*X^I
        hamiltonian -= h*I^X
        return hamiltonian.reduce()


    
    hamiltonian += -J*Z^(I^(num_qubits-2))^Z
    for i in range(num_qubits-1):
        hamiltonian += (-J)*((I^(i))^(Z^Z)^(I^(num_qubits-2-i)))
    
    for i in range(num_qubits):
        hamiltonian -= h*((I^(i))^(X)^(I^(num_qubits-1-i)))

    return hamiltonian.reduce()

# def Transverse_Ising_Measurement(hamiltonian, quantum_circuit, shots, sampler):
#     """
#     hamiltonian: the Pauli matrix we want to measure for the Quantum Nature Gradient Descent with writing in this from (pauli_matrix+is_position). For exammple: (Z1Z2) or (X3)
#     quantum_circuit: the circuit for measurement 
#     shots: the number of shots that we want to execute, if we give the value to shot we will change from statevector_simulator to qasm_simulator
#     """
#     number_of_qubits = quantum_circuit.num_qubits

#     # if len(hamiltonian) > (number_of_qubits):
#     #     raise Exception("The position of Pauli matrix exceeds the number of qubits")

#     Z_list = {}
#     X_list = {}

#     for i in range(len(hamiltonian)):
#         if str(hamiltonian.primitive._pauli_list[i]).__contains__('Z'):
#             Z_list[str(hamiltonian.primitive._pauli_list[i])] = hamiltonian.coeffs.real[i]
#         if str(hamiltonian.primitive._pauli_list[i]).__contains__('X'):
#             X_list[str(hamiltonian.primitive._pauli_list[i])] = hamiltonian.coeffs.real[i]    


    
#     expectation_value_z = 0
#     quantum_circuit_z = quantum_circuit.copy()
#     quantum_circuit_z.measure_all()

#     job_z = sampler.run(quantum_circuit_z, shots = shots)

#     result_z = (job_z.result().quasi_dists[0].binary_probabilities())

#     # print(Z_list)

#     for k in range(len(Z_list)):
#         position_of_non_I_gate = np.where(np.array(list(str(list(Z_list.keys())[k]))) != 'I')[0]
#         # print(list(Z_list.keys())[k])
#         # print((list(Z_list.keys()))[k])
#         # print(position_of_non_I_gate)
#         for i in range(len(result_z)):
#             extra_minus = 1
#             # print( (list(result_z.keys()))[i])
#             for j in range(len(position_of_non_I_gate)):
#                 if (list(result_z.keys()))[i][number_of_qubits - position_of_non_I_gate[j]-1] == '1':
#                     extra_minus *= -1
#             expectation_value_z += result_z[list(result_z.keys())[i]]*extra_minus*list(Z_list.values())[k]

#     quantum_circuit_x_to_meausurement = QuantumCircuit(quantum_circuit.num_qubits)
#     quantum_circuit_x_to_meausurement.h(i for i in range(quantum_circuit.num_qubits))

#     expectation_value_x = 0
#     quantum_circuit_x = quantum_circuit.compose(quantum_circuit_x_to_meausurement)
#     quantum_circuit_x.measure_all()

#     # print(quantum_circuit_x)
    

#     job_x = sampler.run(quantum_circuit_x, shots = shots)

#     result_x = (job_x.result().quasi_dists[0].binary_probabilities())



#     for k in range(len(X_list)):
#         position_of_non_I_gate = np.where(np.array(list(str(list(X_list.keys())[k]))) != 'I')[0]
#         # print(list(X_list.keys())[k])
#         # print((list(Z_list.keys()))[k])
#         # print(position_of_non_I_gate)
#         for i in range(len(result_x)):
#             extra_minus = 1
#             # print( (list(result_x.keys()))[i])
#             for j in range(len(position_of_non_I_gate)):
#                 if (list(result_x.keys()))[i][number_of_qubits - position_of_non_I_gate[j]-1] == '1':
#                     extra_minus *= -1
#             expectation_value_x += result_x[list(result_x.keys())[i]]*extra_minus*list(X_list.values())[k]

#     return expectation_value_z + expectation_value_x

def Transverse_Ising_Measurement(hamiltonian, quantum_circuit, shots, sampler):
    return Estimator().run(quantum_circuit, hamiltonian).result().values[0]



# Parameter shift rule
def Customize_Parameter_Shift_Rule(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: Hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: The number of interation
        shots: The number of shots
        callback: A callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """
    
    internal_initial_point = initial_point.copy()
    
    # Write Interation Energy to variable Energy
    if callback is None:
      energy = []
      internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
      energy.append(internal_energy)
      
    for i in range(interation):
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       

        # grad = np.zeros(ansatz.num_parameters)      
        # for i in range(ansatz.num_parameters):
        #     plus_parameter = internal_initial_point.copy()
        #     plus_parameter[i] += np.pi/2
        #     minus_parameter = internal_initial_point.copy()
        #     minus_parameter[i] -= np.pi/2
        #     grad[i] = (learning_rate)*(Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/2
        
        params = []
        grad = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, np.pi/2, 1, internal_initial_point, operator, ansatz, shots, sampler))

        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        grad = list(executor.map(mini_derivate, params))

        internal_initial_point =  np.subtract(internal_initial_point, learning_rate*grad)

        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')

        if callback is not None:
            callback(internal_initial_point, internal_energy)
    
    if callback is None:
        return energy


# Gradient
def Customize_Finite_Difference(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: Hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: The number of interation
        shots: The number of shots
        callback: A callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """
    
    
    internal_initial_point = initial_point.copy()
    h = 1e-2 # 0.05
    
    # Write Interation Energy to variable Energy
    if callback is None:
      energy = []
      internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
      energy.append(internal_energy)

    for i in range(interation):
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        
        # Evaluate the gradient of the cost function (non-parallel way)
        # grad = np.zeros(ansatz.num_parameters)      
        # for i in range(ansatz.num_parameters):
        #     plus_parameter = internal_initial_point.copy()
        #     plus_parameter[i] += h
        #     minus_parameter = internal_initial_point.copy()
        #     minus_parameter[i] -= h
        #     grad[i] = (learning_rate)*(Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/(2*h)
        
        # Evaluate the gradient of the cost function
        params = []
        grad = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, h, h, internal_initial_point, operator, ansatz, shots, sampler))

        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        grad = list(executor.map(mini_derivate, params))

        internal_initial_point =  np.subtract(internal_initial_point, learning_rate*np.array(grad))
        
        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        
        # Write Interation Energy to variable Energy
        if callback is None:
          energy.append(internal_energy)
          
        # Write Interation Data to file
        if callback is not None:
            callback(internal_initial_point, internal_energy)

    if callback is None:
        return energy


def Separate_Circuit_Apart(ansatz):

    super_circuit = []

    dag = circuit_to_dag(ansatz)

    divide_to_layer_circuit = []
    for layer in dag.layers():
        divide_to_layer_circuit.append(dag_to_circuit(layer['graph']))

    num_layer = 0
    for layer in range(len(divide_to_layer_circuit)):
        layer += num_layer
        if layer == len(divide_to_layer_circuit):
            break
        
        internal_quantum_circuit = divide_to_layer_circuit[layer]
        if internal_quantum_circuit[0].operation.params:
            super_circuit.append(internal_quantum_circuit)
        
        if internal_quantum_circuit[0].operation.name == "barrier":
            sub_circuit = QuantumCircuit(ansatz.num_qubits)
            num_layer += 1
            internal_layer = layer + num_layer

            while divide_to_layer_circuit[internal_layer].data[0].operation.name != "barrier":      
                sub_circuit = sub_circuit.compose(divide_to_layer_circuit[internal_layer])
                internal_layer = layer
                num_layer += 1
                internal_layer += num_layer
            super_circuit.append(sub_circuit)

    return super_circuit



def Customize_Quantum_Natural_Gradient_Descent(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: the number of interation
        shots: The number of shots
        callback: a callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """
      
    def calculate_block_diagonal_fubini_metric(super_circuit, ansatz, parameter, shots, sampler):
        def Measure_element_of_Fubini_Study_metric(circuit, circuit_for_measurement, i, j, shots, sampler):
            if i != j:
                term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)
                term3 = ['I']*len(circuit)
                # print(circuit)
                # Change rotation gate to pauli gate
                # Term1
                if circuit[i].operation.name == 'rx':
                    term1[i] = 'X'
                    term2[i] = 'X'
                if circuit[i].operation.name == 'ry':
                    term1[i] = 'Y'
                    term2[i] = 'Y'
                if circuit[i].operation.name == 'rz':
                    term1[i] = 'Z'
                    term2[i] = 'Z'
                
                if circuit[j].operation.name == 'rx':
                    term1[j] = 'X'
                    term3[j] = 'X'
                if circuit[j].operation.name == 'ry':
                    term1[j] = 'Y'
                    term3[j] = 'Y'
                if circuit[j].operation.name == 'rz':
                    term1[j] = 'Z'
                    term3[j] = 'Z'
    
                # # Term 2
                # if circuit[i].operation.name == 'rx':
                #     term2[i] = 'X'
                # if circuit[i].operation.name == 'ry':
                #     term2[i] = 'Y'
                # if circuit[i].operation.name == 'rz':
                #     term2[i] = 'Z'
    
                # # Term 3
                # if circuit[j].operation.name == 'rx':
                #     term3[j] = 'X'
                # if circuit[j].operation.name == 'ry':
                #     term3[j] = 'Y'
                # if circuit[j].operation.name == 'rz':
                #     term3[j] = 'Z'
                
    
                term1 = ''.join(term1[::-1])
                term2 = ''.join(term2[::-1])
                term3 = ''.join(term3[::-1])
    
                return  (Transverse_Ising_Measurement(term1, circuit_for_measurement, shots, sampler) - Transverse_Ising_Measurement(term2, circuit_for_measurement, shots, sampler)*Transverse_Ising_Measurement(term3, circuit_for_measurement, shots, sampler))/4
            else:
                # term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)
    
                if circuit[i].operation.name == 'rx':
                    term2[i] = 'X'
                if circuit[i].operation.name == 'ry':
                    term2[i] = 'Y'
                if circuit[i].operation.name == 'rz':
                    term2[i] = 'Z'
    
                term2 = ''.join(term2[::-1])
    
                # new_gate = circuit_for_measurement.compose(new_gate)
    
                # return (Measurement(term1, circuit_for_measurement, shots, backend)[0] - Measurement(term2, circuit_for_measurement, shots, backend)[0].real**2)/4
                return (1 - Transverse_Ising_Measurement(term2, circuit_for_measurement, shots, sampler).real**2)/4
  
      
        # Caculate the Fubini_study_metric
        fubini_study_metric = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
        initial_point = parameter.copy()
        # num_parameter = 0
    
    
        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                parameter_previous = 0
                internal_circuit = QuantumCircuit(ansatz.num_qubits)
    
                for j in range(i):
                    parameter_previous += super_circuit[j].num_parameters
                    internal_circuit = internal_circuit.compose(super_circuit[j])
                    
                # print(internal_circuit)
                internal_circuit = internal_circuit.bind_parameters({theta: initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                # num_parameter += super_circuit[i].num_parameters
    
                # print(num_parameter)
                
                for l in range(super_circuit[i].num_parameters):
                    for m in range(l+1, super_circuit[i].num_parameters):
                        fubini_study_metric[l + parameter_previous][m + parameter_previous] = Measure_element_of_Fubini_Study_metric(super_circuit[i], internal_circuit, l, m, shots, sampler)
                        fubini_study_metric[m + parameter_previous][l + parameter_previous] = fubini_study_metric[l + parameter_previous][m + parameter_previous]
                for l in range(super_circuit[i].num_parameters):
                        fubini_study_metric[l + parameter_previous][l + parameter_previous] = Measure_element_of_Fubini_Study_metric   (super_circuit[i], internal_circuit, l, l, shots, sampler)
        return fubini_study_metric
                        
    energy = []
    internal_initial_point = initial_point.copy()
    super_circuit = Separate_Circuit_Apart(ansatz)
    
    for i in range(interation):
        
        fubini_study_metric = np.array(calculate_block_diagonal_fubini_metric(super_circuit, ansatz, internal_initial_point, shots, sampler))
        
        # Exact Gradient of the cost function
        # grad = np.zeros(ansatz.num_parameters)
        # for i in range(ansatz.num_parameters):
        #     plus_parameter = internal_initial_point.copy()
        #     plus_parameter[i] += np.pi/2
        #     minus_parameter = internal_initial_point.copy()
        #     minus_parameter[i] -= np.pi/2
        #     grad[i] = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/2

    
        params = []
        grad = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, np.pi/2, 1, internal_initial_point, operator, ansatz, shots, sampler))


        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        grad = list(executor.map(mini_derivate, params))


        #print(fubini_study_metric)
        FS_metric_inv = np.linalg.pinv(fubini_study_metric)
        # print(FS_metric_inv)
        # print(grad)

        combine = learning_rate*FS_metric_inv.dot(np.array(grad))

        internal_initial_point =  np.subtract(internal_initial_point, combine)
    
        # Measure the expectation of our hamiltonian
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})   
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        energy.append(internal_energy)
        print(internal_energy)
        #print(f'{internal_initial_point} ---------')

        if callback is not None:
            callback(internal_initial_point, internal_energy)
    
    if callback is None:
        return energy


def Customize_SPSA(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: the number of interation
        shots: The number of shots
        callback: a callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """
    
    internal_initial_point = initial_point.copy()

    
    # Write Interation Energy to variable Energy
    if callback is None:
      energy = []
      internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
      energy.append(internal_energy)

    for k in range(interation):        
        internal_energy =  Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        #print(f'{internal_initial_point} ---------')
        #print(internal_energy)
        

        # Optimization part
        grad = np.zeros(ansatz.num_parameters)
        # ak = a/(1+k+A)**alpha
        # ck = c/(1+k)**gamma
        ak = learning_rate
        ck = 0.01
        
        # print(f'ak_{k}: {ak}')
        # print(f'ck_{k}: {ck}')
        random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
        plus_parameter = np.array(internal_initial_point.copy())
        plus_parameter = np.add(plus_parameter,random*ck)
        minus_parameter = np.array(internal_initial_point.copy())
        minus_parameter = np.subtract(minus_parameter,random*ck)
        grad_func = ak*(Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)) / (2*ck)
        
        grad = np.add(grad, grad_func*random)

        internal_initial_point = np.subtract(internal_initial_point, grad)
        
        if callback is None:
            energy.append(internal_energy)
        
        if callback is not None:
            callback(internal_initial_point, internal_energy)
    if callback is None:
        return energy



    
    
def Customize_QNSPSA_PRS_blocking(operator, initial_point, learning_rate, ansatz, interation, step, shots, callback, sampler, previous_fubini_matrix, last_n_steps):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: the number of interation
        shots: The number of shots
        callback: a callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """

    beta = 0.001

    internal_initial_point = initial_point.copy()

    energy = []

    regularized_fubini_matrix_previous = previous_fubini_matrix.copy()
    grad = np.zeros(ansatz.num_parameters)
    
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    
    energy.append(internal_energy)
    
    #print(last_n_steps)
  
    for k in range(interation):  
        k += step + 1
        #print(k)
        next_energy = 0

        # gradPRS = np.zeros(ansatz.num_parameters)
        # for i in range(ansatz.num_parameters):
        #     plus_parameter = internal_initial_point.copy()
        #     plus_parameter[i] += np.pi/2
        #     minus_parameter = internal_initial_point.copy()
        #     minus_parameter[i] -= np.pi/2
        #     gradPRS[i] = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/2
        params = []
        gradPRS = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, np.pi/2, 1, internal_initial_point, operator, ansatz, shots, sampler))

        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        gradPRS = list(executor.map(mini_derivate, params))

        # print('gradPRS: ', gradPRS)

        while True:    
            grad = np.zeros(ansatz.num_parameters)
            print(last_n_steps)
            ck = 0.01
            ak = learning_rate
            # Natural Gradient Part
            random1 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            random2 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            #random1 = np.array([1.0 for _ in range(ansatz.num_parameters)])
            #random2 = np.array([-1.0 for _ in range(ansatz.num_parameters)])
            initial_plus1_plus2 = np.add(internal_initial_point, np.add(ck*random1, ck*random2))
            initial_plus1 = np.add(internal_initial_point, ck*random1)
            initial_minus1_plus2 = np.subtract(internal_initial_point, np.subtract(ck*random1, ck*random2))
            initial_minus1 = np.subtract(internal_initial_point, ck*random1)

            ansatz_initial = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1_plus2 = ansatz.bind_parameters({theta: initial_plus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1 = ansatz.bind_parameters({theta: initial_plus1[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1_plus2 = ansatz.bind_parameters({theta: initial_minus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1 =  ansatz.bind_parameters({theta: initial_minus1[k] for k, theta in enumerate(ansatz.parameters)})  

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2) - SwapTest(ansatz_initial, ansatz_plus1) - SwapTest(ansatz_initial, ansatz_minus1_plus2) + SwapTest(ansatz_initial, ansatz_minus1)


            # print(deltaF)
            # Fubini_matrix
            fubini_matrix = -1/2*(deltaF/(2*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))/2
            
            # exponentially_smoothed_fubini  = fubini_matrix

            # Data of previous regularized study metric
            exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*fubini_matrix.copy()    
            regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

            # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()

            grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradPRS)
            
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            tolerance = 2 * last_n_steps.std() if (k-1 > len(last_n_steps)) else 2 * last_n_steps[:k].std()
            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            #print(next_energy)
            #print(f'tolerance: {tolerance}')
            # print('gradPRS: ', gradPRS)
            #print(next_energy)

            if next_energy <= internal_energy + tolerance:
                break
            
            #print(next_energy)
            # if True:
                # break            

        
        regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()
        
        # ghi data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')
        # print(f'gradSPSA {gradSPSA}')
        # print(f'regularized_fubini_matrix: {np.linalg.pinv(regularized_fubini_matrix)}')

        if callback is not None:
            callback(internal_initial_point, internal_energy, regularized_fubini_matrix_previous)

        last_n_steps[(k) % len(last_n_steps)] = internal_energy
    return energy


def Customize_QNSPSA_PRS_blocking_MonteCarlo(operator, initial_point, learning_rate, ansatz, interation, step, shots, callback, sampler, previous_fubini_matrix, last_n_steps):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: the number of interation
        shots: The number of shots
        callback: a callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """

    beta = 0.001

    internal_initial_point = initial_point.copy()


    energy = []

    regularized_fubini_matrix_previous = previous_fubini_matrix.copy()
    grad = np.zeros(ansatz.num_parameters)
    
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)
    #print(last_n_steps)
    #print(internal_energy)
    for k in range(interation):
        #print(last_n_steps)
        k += step + 1
        next_energy = 0

        # gradPRS = np.zeros(ansatz.num_parameters)
        # for i in range(ansatz.num_parameters):
        #     plus_parameter = internal_initial_point.copy()
        #     plus_parameter[i] += np.pi/2
        #     minus_parameter = internal_initial_point.copy()
        #     minus_parameter[i] -= np.pi/2
        #     gradPRS[i] = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/2
        params = []
        gradPRS = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, np.pi/2, 1, internal_initial_point, operator, ansatz, shots, sampler))

        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        gradPRS = list(executor.map(mini_derivate, params))

        # print('gradPRS: ', gradPRS)
        
        monte_carlo_matrix = []

        while True:
            if len(monte_carlo_matrix) > int(ansatz.num_parameters*ansatz.num_qubits/4):
                break
                
            #print(k)
            #print("True")  
            grad = np.zeros(ansatz.num_parameters)
            ck = 0.01
            ak = learning_rate
            # Natural Gradient Part
            random1 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            random2 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            #random1 = np.array([1.0 for _ in range(ansatz.num_parameters)])
            #random2 = np.array([-1.0 for _ in range(ansatz.num_parameters)])
            initial_plus1_plus2 = np.add(internal_initial_point, np.add(ck*random1, ck*random2))
            initial_plus1 = np.add(internal_initial_point, ck*random1)
            initial_minus1_plus2 = np.subtract(internal_initial_point, np.subtract(ck*random1, ck*random2))
            initial_minus1 = np.subtract(internal_initial_point, ck*random1)

            ansatz_initial = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1_plus2 = ansatz.bind_parameters({theta: initial_plus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1 = ansatz.bind_parameters({theta: initial_plus1[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1_plus2 = ansatz.bind_parameters({theta: initial_minus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1 =  ansatz.bind_parameters({theta: initial_minus1[k] for k, theta in enumerate(ansatz.parameters)})  

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2) - SwapTest(ansatz_initial, ansatz_plus1) - SwapTest(ansatz_initial, ansatz_minus1_plus2) + SwapTest(ansatz_initial, ansatz_minus1)


            # print(deltaF)
            # Fubini_matrix
            fubini_matrix = -1/2*(deltaF/(2*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))/2
            
            monte_carlo_matrix.append(fubini_matrix.copy())

            # Data of previous regularized study metric
            exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*np.array(np.mean(monte_carlo_matrix, axis=0))    
            regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

            # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()

            grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradPRS)

            
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            tolerance = 2 * last_n_steps.std() if (k-1 > len(last_n_steps) ) else 2 * last_n_steps[:k].std()
            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            #print(f'tolerance: {tolerance}')
            # print('gradPRS: ', gradPRS)
            #(next_energy)

            if next_energy <= internal_energy + tolerance:
                break

            # if True:
                # break            

        
        regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()
        
        # ghi data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')
        # print(f'gradSPSA {gradSPSA}')
        # print(f'regularized_fubini_matrix: {np.linalg.pinv(regularized_fubini_matrix)}')

        if callback is not None:
            callback(internal_initial_point, internal_energy, regularized_fubini_matrix_previous)

        last_n_steps[(k) % len(last_n_steps)] = internal_energy
    return energy
    


def Customize_QN_SPSA_blocking(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """
    
    def calculate_block_diagonal_fubini_metric(super_circuit, ansatz, parameter, shots, sampler):
        def Measure_element_of_Fubini_Study_metric(circuit, circuit_for_measurement, i, j, shots, sampler):
            if i != j:
                term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)
                term3 = ['I']*len(circuit)
                # print(circuit)
                # Change rotation gate to pauli gate
                # Term1
                if circuit[i].operation.name == 'rx':
                    term1[i] = 'X'
                    term2[i] = 'X'
                if circuit[i].operation.name == 'ry':
                    term1[i] = 'Y'
                    term2[i] = 'Y'
                if circuit[i].operation.name == 'rz':
                    term1[i] = 'Z'
                    term2[i] = 'Z'
                
                if circuit[j].operation.name == 'rx':
                    term1[j] = 'X'
                    term3[j] = 'X'
                if circuit[j].operation.name == 'ry':
                    term1[j] = 'Y'
                    term3[j] = 'Y'
                if circuit[j].operation.name == 'rz':
                    term1[j] = 'Z'
                    term3[j] = 'Z'
    
                # # Term 2
                # if circuit[i].operation.name == 'rx':
                #     term2[i] = 'X'
                # if circuit[i].operation.name == 'ry':
                #     term2[i] = 'Y'
                # if circuit[i].operation.name == 'rz':
                #     term2[i] = 'Z'
    
                # # Term 3
                # if circuit[j].operation.name == 'rx':
                #     term3[j] = 'X'
                # if circuit[j].operation.name == 'ry':
                #     term3[j] = 'Y'
                # if circuit[j].operation.name == 'rz':
                #     term3[j] = 'Z'
                
    
                term1 = ''.join(term1[::-1])
                term2 = ''.join(term2[::-1])
                term3 = ''.join(term3[::-1])
    
                return  (Transverse_Ising_Measurement(term1, circuit_for_measurement, shots, sampler) - Transverse_Ising_Measurement(term2, circuit_for_measurement, shots, sampler)*Transverse_Ising_Measurement(term3, circuit_for_measurement, shots, sampler))/4
            else:
                # term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)
    
                if circuit[i].operation.name == 'rx':
                    term2[i] = 'X'
                if circuit[i].operation.name == 'ry':
                    term2[i] = 'Y'
                if circuit[i].operation.name == 'rz':
                    term2[i] = 'Z'
    
                term2 = ''.join(term2[::-1])
    
                # new_gate = circuit_for_measurement.compose(new_gate)
    
                # return (Measurement(term1, circuit_for_measurement, shots, backend)[0] - Measurement(term2, circuit_for_measurement, shots, backend)[0].real**2)/4
                return (1 - Transverse_Ising_Measurement(term2, circuit_for_measurement, shots, sampler).real**2)/4
  
      
        # Caculate the Fubini_study_metric
        fubini_study_metric = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
        initial_point = parameter.copy()
        # num_parameter = 0
    
    
        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                parameter_previous = 0
                internal_circuit = QuantumCircuit(ansatz.num_qubits)
    
                for j in range(i):
                    parameter_previous += super_circuit[j].num_parameters
                    internal_circuit = internal_circuit.compose(super_circuit[j])
                    
                # print(internal_circuit)
                internal_circuit = internal_circuit.bind_parameters({theta: initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                # num_parameter += super_circuit[i].num_parameters
    
                # print(num_parameter)
                
                for l in range(super_circuit[i].num_parameters):
                    for m in range(l+1, super_circuit[i].num_parameters):
                        fubini_study_metric[l + parameter_previous][m + parameter_previous] = Measure_element_of_Fubini_Study_metric(super_circuit[i], internal_circuit, l, m, shots, sampler)
                        fubini_study_metric[m + parameter_previous][l + parameter_previous] = fubini_study_metric[l + parameter_previous][m + parameter_previous]
                for l in range(super_circuit[i].num_parameters):
                        fubini_study_metric[l + parameter_previous][l + parameter_previous] = Measure_element_of_Fubini_Study_metric   (super_circuit[i], internal_circuit, l, l, shots, sampler)
        return fubini_study_metric


    
    internal_initial_point = initial_point.copy()


    energy = []


    super_circuit = Separate_Circuit_Apart(ansatz)

    grad = np.zeros(ansatz.num_parameters)

    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)

    for k in range(interation):         
        fubini_study_metric = np.array(calculate_block_diagonal_fubini_metric(super_circuit, ansatz, internal_initial_point, shots, sampler))
  
        grad = np.zeros(ansatz.num_parameters)
        # SPSA part
        gradSPSA = np.zeros(ansatz.num_parameters)
        # ak = a/(1+k+A)**alpha
        # ck = c/(1+k)**gamma
        ck = 0.01
        ak = learning_rate
        
        # print(f'ak_{k}: {ak}')
        # print(f'ck_{k}: {ck}')


        # SPSA
        random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
        plus_parameter = np.array(internal_initial_point.copy())
        plus_parameter = np.add(plus_parameter,random*ck)
        minus_parameter = np.array(internal_initial_point.copy())
        minus_parameter = np.subtract(minus_parameter,random*ck)
        gradSPSA = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/(2*(ck))*random
        
        
        
        grad = learning_rate*np.linalg.pinv(fubini_study_metric).dot(gradSPSA)            


        
        # write data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        energy.append(internal_energy)
 

        if callback is not None:
            callback(internal_initial_point, internal_energy)

    return energy

    
def Customize_QNSPSA_SPSA_blocking(operator, initial_point, learning_rate, ansatz, interation, step, shots, callback, sampler, previous_fubini_matrix, last_n_steps):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: Hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: The number of interation
        shots: The number of shots
        callback: A callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """

    beta = 0.001
    
    internal_initial_point = initial_point.copy()

    regularized_fubini_matrix_previous = previous_fubini_matrix.copy()
    grad = np.zeros(ansatz.num_parameters)
    
    
    energy = []
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)

    for k in range(interation):   
        k += step + 1
        next_energy = 0
        
        while True:    
            grad = np.zeros(ansatz.num_parameters)
            # SPSA part
            gradSPSA = np.zeros(ansatz.num_parameters)
            # ak = a/(1+k+A)**alpha
            # ck = c/(1+k)**gamma
            ck = 0.01
            ak = learning_rate
            
            random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            plus_parameter = np.array(internal_initial_point.copy())
            plus_parameter = np.add(plus_parameter,random*ck)
            minus_parameter = np.array(internal_initial_point.copy())
            minus_parameter = np.subtract(minus_parameter,random*ck)
            gradSPSA = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/(2*(ck))*random
            
            # print(f'gradSPSA {gradSPSA}')

            
            # Natural Gradient Part
            random1 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            random2 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            # random1 = np.array([1.0 for _ in range(ansatz.num_parameters)])
            # random2 = np.array([-1.0 for _ in range(ansatz.num_parameters)])
            initial_plus1_plus2 = np.add(internal_initial_point, np.add(ck*random1, ck*random2))
            initial_plus1 = np.add(internal_initial_point, ck*random1)
            initial_minus1_plus2 = np.subtract(internal_initial_point, np.subtract(ck*random1, ck*random2))
            initial_minus1 = np.subtract(internal_initial_point, ck*random1)

            ansatz_initial = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1_plus2 = ansatz.bind_parameters({theta: initial_plus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1 = ansatz.bind_parameters({theta: initial_plus1[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1_plus2 = ansatz.bind_parameters({theta: initial_minus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1 =  ansatz.bind_parameters({theta: initial_minus1[k] for k, theta in enumerate(ansatz.parameters)})  

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2) - SwapTest(ansatz_initial, ansatz_plus1) -SwapTest(ansatz_initial, ansatz_minus1_plus2) + SwapTest(ansatz_initial, ansatz_minus1)


            # print(deltaF)
            # Fubini_matrix
            fubini_matrix = -1/2*(deltaF/(2*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))/2
            
            # Data of previous regularized study metric
            exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*fubini_matrix.copy()    
            regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

            # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()

            grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradSPSA)

            
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            
            tolerance = 2 * last_n_steps.std() if (k-1 > len(last_n_steps)) else 2 * last_n_steps[:k].std()

            #print(f'tolerance: {tolerance}')

            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            if next_energy <= internal_energy + tolerance:
                break            
            
        
        regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()
        
        # ghi data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        if callback is None:
          energy.append(internal_energy)

        if callback is not None:
            callback(internal_initial_point, internal_energy, regularized_fubini_matrix_previous)

        last_n_steps[(k) % len(last_n_steps)] = internal_energy
    if callback is None:
      return energy
      
      
      
def Customize_QNSPSA_SPSA_MonteCarlo(operator, initial_point, learning_rate, ansatz, interation, step, shots, callback, sampler, previous_fubini_matrix, last_n_steps):
    """
    Args:
        operator: PauliOp
        initial_point: An initial parameter values for the optimizer. 
        learning_rate: Hyper-parameter used to govern the pace at which an algorithm updates or learns the values of a parameter estimate
        ansatz: A parameterized circuit used as Ansatz for the wave function.
        interation: The number of interation
        shots: The number of shots
        callback: A callback that can access the intermediate data during the optimization.
        sampler: Sampler instance.
    
    Returns:
        Interation energy 
    """

    beta = 0.001
    
    internal_initial_point = initial_point.copy()

    regularized_fubini_matrix_previous = previous_fubini_matrix.copy()
    grad = np.zeros(ansatz.num_parameters)
    
    
    energy = []
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)

    for k in range(interation):   
        #(last_n_steps)
        k += step + 1
        next_energy = 0
        
        monte_carlo_matrix = []
        monte_carlo_SPSA = []
        
        while True:    
            grad = np.zeros(ansatz.num_parameters)
            # SPSA part
            gradSPSA = np.zeros(ansatz.num_parameters)
            # ak = a/(1+k+A)**alpha
            # ck = c/(1+k)**gamma
            ck = 0.01
            ak = learning_rate
            
            random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            plus_parameter = np.array(internal_initial_point.copy())
            plus_parameter = np.add(plus_parameter,random*ck)
            minus_parameter = np.array(internal_initial_point.copy())
            minus_parameter = np.subtract(minus_parameter,random*ck)
            gradSPSA = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/(2*(ck))*random
            
            monte_carlo_SPSA.append(gradSPSA.copy())
            # print(f'gradSPSA {gradSPSA}')

            
            # Natural Gradient Part
            random1 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            random2 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            # random1 = np.array([1.0 for _ in range(ansatz.num_parameters)])
            # random2 = np.array([-1.0 for _ in range(ansatz.num_parameters)])
            initial_plus1_plus2 = np.add(internal_initial_point, np.add(ck*random1, ck*random2))
            initial_plus1 = np.add(internal_initial_point, ck*random1)
            initial_minus1_plus2 = np.subtract(internal_initial_point, np.subtract(ck*random1, ck*random2))
            initial_minus1 = np.subtract(internal_initial_point, ck*random1)

            ansatz_initial = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1_plus2 = ansatz.bind_parameters({theta: initial_plus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_plus1 = ansatz.bind_parameters({theta: initial_plus1[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1_plus2 = ansatz.bind_parameters({theta: initial_minus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
            ansatz_minus1 =  ansatz.bind_parameters({theta: initial_minus1[k] for k, theta in enumerate(ansatz.parameters)})  

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2) - SwapTest(ansatz_initial, ansatz_plus1) -SwapTest(ansatz_initial, ansatz_minus1_plus2) + SwapTest(ansatz_initial, ansatz_minus1)


            # print(deltaF)
            # Fubini_matrix
            fubini_matrix = -1/2*(deltaF/(2*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))/2
            
            monte_carlo_matrix.append(fubini_matrix.copy())
            
            # Data of previous regularized study metric
            exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*np.mean(monte_carlo_matrix, axis=0)   
            regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

            # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()

            grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradSPSA)

            monte_carlo_SPSA.append(gradSPSA.copy())
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            
            tolerance = 2 * last_n_steps.std() if (k-1 > len(last_n_steps)) else 2 * last_n_steps[:k].std()

            #print(f'tolerance: {tolerance}')

            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            if next_energy <= internal_energy + tolerance:
                break            
            
        
        regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()
        
        # ghi data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        if callback is None:
          energy.append(internal_energy)

        if callback is not None:
            callback(internal_initial_point, internal_energy, regularized_fubini_matrix_previous)

        last_n_steps[(k) % len(last_n_steps)] = internal_energy
    if callback is None:
      return energy
    


