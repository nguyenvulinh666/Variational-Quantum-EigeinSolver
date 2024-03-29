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
    
    return result, qc
    
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

def Transverse_Ising_Measurement(hamiltonian, quantum_circuit, shots, sampler):
    """
    hamiltonian: the Pauli matrix we want to measure for the Quantum Nature Gradient Descent with writing in this from (pauli_matrix+is_position). For exammple: (Z1Z2) or (X3)
    quantum_circuit: the circuit for measurement 
    shots: the number of shots that we want to execute, if we give the value to shot we will change from statevector_simulator to qasm_simulator
    """
    number_of_qubits = quantum_circuit.num_qubits

    # if len(hamiltonian) > (number_of_qubits):
    #     raise Exception("The position of Pauli matrix exceeds the number of qubits")

    Z_list = {}
    X_list = {}

    for i in range(len(hamiltonian)):
        if str(hamiltonian.primitive._pauli_list[i]).__contains__('Z'):
            Z_list[str(hamiltonian.primitive._pauli_list[i])] = hamiltonian.coeffs.real[i]
        if str(hamiltonian.primitive._pauli_list[i]).__contains__('X'):
            X_list[str(hamiltonian.primitive._pauli_list[i])] = hamiltonian.coeffs.real[i]    


    
    expectation_value_z = 0
    quantum_circuit_z = quantum_circuit.copy()
    quantum_circuit_z.measure_all()

    job_z = sampler.run(quantum_circuit_z, shots = shots)

    result_z = (job_z.result().quasi_dists[0].binary_probabilities())

    # print(Z_list)

    for k in range(len(Z_list)):
        position_of_non_I_gate = np.where(np.array(list(str(list(Z_list.keys())[k]))) != 'I')[0]
        # print(list(Z_list.keys())[k])
        # print((list(Z_list.keys()))[k])
        # print(position_of_non_I_gate)
        for i in range(len(result_z)):
            extra_minus = 1
            # print( (list(result_z.keys()))[i])
            for j in range(len(position_of_non_I_gate)):
                if (list(result_z.keys()))[i][number_of_qubits - position_of_non_I_gate[j]-1] == '1':
                    extra_minus *= -1
            expectation_value_z += result_z[list(result_z.keys())[i]]*extra_minus*list(Z_list.values())[k]

    quantum_circuit_x_to_meausurement = QuantumCircuit(quantum_circuit.num_qubits)
    quantum_circuit_x_to_meausurement.h(i for i in range(quantum_circuit.num_qubits))

    expectation_value_x = 0
    quantum_circuit_x = quantum_circuit.compose(quantum_circuit_x_to_meausurement)
    quantum_circuit_x.measure_all()

    # print(quantum_circuit_x)
    

    job_x = sampler.run(quantum_circuit_x, shots = shots)

    result_x = (job_x.result().quasi_dists[0].binary_probabilities())



    for k in range(len(X_list)):
        position_of_non_I_gate = np.where(np.array(list(str(list(X_list.keys())[k]))) != 'I')[0]
        # print(list(X_list.keys())[k])
        # print((list(Z_list.keys()))[k])
        # print(position_of_non_I_gate)
        for i in range(len(result_x)):
            extra_minus = 1
            # print( (list(result_x.keys()))[i])
            for j in range(len(position_of_non_I_gate)):
                if (list(result_x.keys()))[i][number_of_qubits - position_of_non_I_gate[j]-1] == '1':
                    extra_minus *= -1
            expectation_value_x += result_x[list(result_x.keys())[i]]*extra_minus*list(X_list.values())[k]

    return expectation_value_z + expectation_value_x


def Customize_EfficientSU2(number_qubits, number_of_subcircuit):
    """
    number_qubits: The amounts of qubits in out system
    number_of_subcircuit: The amounts of subcuirt for parameterized our system
    su2gate: is the gate we want to parameterized with
    insert_barriers: add barriers in our circuit for nice looking
    The function will return the circuit having the parameter and we can update these parameter in our code
    """
    circuit = QuantumCircuit(number_qubits)

    theta = ParameterVector(r'$\theta$', number_qubits*2 + number_of_subcircuit*number_qubits*2)


    # Use for the 1-qubit case
    if number_qubits == 1:
        theta1 = Parameter(theta[0])
        circuit.ry(theta1, 0)
        theta2 = Parameter(theta[1])
        circuit.rz(theta2, 0)
        circuit.barrier()

        for i in range(number_of_subcircuit):
            theta1 = Parameter(theta[2*i+2])
            circuit.ry(theta1, 0)
            theta2 = Parameter(theta[2*i+3])
            circuit.rz(theta2, 0)

            if i != number_of_subcircuit - 1:
                circuit.barrier()

        return circuit


    def add_subcircuit(circuit, stop_barrier):
        sub_circuit = QuantumCircuit(number_qubits)
        number_parameter = circuit.num_parameters

        sub_circuit.barrier()

        # sub_circuit.cx(number_qubits-1, 0)
        for i in (range(number_qubits-1)):
            sub_circuit.cx(i, i+1)

        # sub_circuit.cx(number_qubits-1, 0)

        sub_circuit.barrier()


        # hmm cause of lack of my knowledge, I will work with the RealAmplitudes ansat, circular entanglement
        for i in range(0, number_qubits):
            # theta = Parameter(theta[number_parameter+i])
            sub_circuit.ry(theta[number_parameter+i], i)
        for i in range(0, number_qubits):
            # theta = Parameter(theta[number_qubits+number_parameter+i])
            sub_circuit.rz(theta[number_qubits+number_parameter+i], i)



        # if stop_barrier != number_of_subcircuit-1:
        #     sub_circuit.barrier()

        # sub_circuit.draw('mpl', style = 'iqx')

        return sub_circuit


    for i in range(0, number_qubits):
        # theta = Parameter(theta[i])
        circuit.ry(theta[i], i)

    for i in range(0, number_qubits):
        # theta = Parameter(theta[number_qubits+i])
        circuit.rz(theta[number_qubits+i], i)

    for i in range(number_of_subcircuit):
        circuit = circuit.compose(add_subcircuit(circuit, i))

    return circuit


def Customize_RealAmplidues(number_qubits, number_of_subcircuit):
    """
    number_qubits: The amounts of qubits in out system
    number_of_subcircuit: The amounts of subcuirt for parameterized our system
    su2gate: is the gate we want to parameterized with
    insert_barriers: add barriers in our circuit for nice looking
    The function will return the circuit having the parameter and we can update these parameter in our code
    """
    circuit = QuantumCircuit(number_qubits)


    # Use for the 1-qubit case
    if number_qubits == 1:

        for i in range(0, number_qubits):
                theta = Parameter(r'$\theta[{}]$'.format(i))
                circuit.ry(theta, i)
                circuit.barrier()
        for i in range(number_of_subcircuit):
            theta1 = Parameter(r'$\theta[{}]$'.format(1+i))
            circuit.ry(theta1, 0)
            if i != number_of_subcircuit - 1:
                circuit.barrier()


        return circuit


    def add_subcircuit(circuit, stop_barrier):
        sub_circuit = QuantumCircuit(number_qubits)
        number_parameter = circuit.num_parameters

        sub_circuit.barrier()

        # sub_circuit.cx(number_qubits-1, 0)
        for i in (range(number_qubits-1)):
            sub_circuit.cx(i, i+1)

        sub_circuit.barrier()


        # hmm cause of lack of my knowledge, I will work with the RealAmplitudes ansat, circular entanglement
        for i in range(0, number_qubits):
            theta = Parameter(r'$\theta[{}]$'.format(number_parameter+i))
            sub_circuit.ry(theta, i)


        # if stop_barrier != number_of_subcircuit-1:
        #     sub_circuit.barrier()

        # sub_circuit.draw('mpl', style = 'iqx')

        return sub_circuit


    for i in range(0, number_qubits):
        theta = Parameter(r'$\theta[{}]$'.format(i))
        circuit.ry(theta, i)

    # for i in range(0, number_qubits):
    #     theta = Parameter(r'$\theta[{}]$'.format(number_qubits+i))
    #     circuit.rz(theta, i)

    for i in range(number_of_subcircuit):
        circuit = circuit.compose(add_subcircuit(circuit, i))

    return circuit


# Parameter shift rule
def Customize_Parameter_Shift_Rule(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    learning_rate: learning rate
    """
    energy = []
    internal_initial_point = initial_point.copy()
    

    for i in range(interation):
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')

        if callback is not None:
            callback(internal_initial_point, internal_energy)

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
    
    if callback is None:
        return energy


# Gradient
def Customize_Finite_Difference(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    learning_rate: learning rate
    """
    energy = []
    internal_initial_point = initial_point.copy()
    h = 1e-2 # 0.05

    

    for i in range(interation):
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')

        if callback is not None:
            callback(internal_initial_point, internal_energy)

        
        params = []
        grad = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            params.append((i, h, h, internal_initial_point, operator, ansatz, shots, sampler))

        import concurrent.futures
        executor = concurrent.futures.ProcessPoolExecutor()
        grad = list(executor.map(mini_derivate, params))

        internal_initial_point =  np.subtract(internal_initial_point, learning_rate*grad)
    
    if callback is None:
        return energy


def Separate_Circuit_Apart(ansatz):
    # Divide the circuit to subcircuit of parameter circuit and non parameter circuit
    """
    Note that Separate circuit ansatz just works with the custom ansatz, and not work well with the qiskit ansatz
    """
    super_circuit = []
    no_name = 0

    ansatz_barrier = 0 

    for i in range(len(ansatz)):
        if ansatz[i].operation.name == 'barrier':
            ansatz_barrier += 1


    while no_name < (ansatz.size() + ansatz_barrier):
        if ansatz[no_name].operation.params:
            sub_circuit = QuantumCircuit(ansatz.num_qubits)
            for i in range(no_name, ansatz.num_qubits + no_name):
                sub_circuit.append(ansatz[i])
            super_circuit.append(sub_circuit)
            no_name += ansatz.num_qubits 
            
        elif ansatz[no_name].operation.name == 'barrier':
            no_name += 1
            
        else:
            sub_circuit = QuantumCircuit(ansatz.num_qubits)
            
            while not ansatz[no_name].operation.params:
                if ansatz[no_name].operation.name == 'barrier':
                    no_name += 1
                    break
                sub_circuit.append(ansatz[no_name])
                no_name += 1 

                if no_name == (ansatz.size() + ansatz_barrier - 1):
                    break

            super_circuit.append(sub_circuit)
    return super_circuit




def Customize_Quantum_Natural_Gradient_Descent(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    
    
    energy = []
    internal_initial_point = initial_point.copy()

    super_circuit = Separate_Circuit_Apart(ansatz)
    

    for i in range(interation):
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')

        if callback is not None:
            callback(internal_initial_point, internal_energy)


        fubini_study_metric = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
        # Measure the fubini-study metric
        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                g_internal = [[[] for _ in range(super_circuit[i].num_parameters)] for _ in range(super_circuit[i].num_parameters)]
                internal_circuit = QuantumCircuit(super_circuit[i].num_qubits)

                parameter_previous = 0
                for j in range(i):
                    parameter_previous += super_circuit[j].num_parameters
                    internal_circuit = internal_circuit.compose(super_circuit[j])
                # print(parameter_previous)
                

                internal_circuit = internal_circuit.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                
                # Add string
                string = ['I']*ansatz.num_qubits

                for l in range(super_circuit[i].num_parameters):
                    if super_circuit[i][l].operation.name == 'rx':
                        string[l] = 'X'
                    if super_circuit[i][l].operation.name == 'ry':
                        string[l] = 'Y'
                    if super_circuit[i][l].operation.name == 'rz':
                        string[l] = 'Z'

                quantum_circuit_to_meausurement = QuantumCircuit(ansatz.num_qubits)

                for l in range(len(string)):
                    if string[l] == "X":
                        quantum_circuit_to_meausurement.h(l)
                    if string[l] == "Y":
                        quantum_circuit_to_meausurement.rx(np.pi/2, l)
                    
                
                quantum_circuit = internal_circuit.compose(quantum_circuit_to_meausurement)
                quantum_circuit.measure_all()
                # print('string: ', string)
                # print(quantum_circuit)
                job = sampler.run(quantum_circuit, shots = shots)

                result = (job.result().quasi_dists[0].binary_probabilities())
                # print(result)
                # term 2 in fubini-study
                term2 = np.zeros(len(string))
                for l in range(super_circuit[i].num_parameters):
                    expectation_value = 0
                    position_of_non_I_gate = [(ansatz.num_qubits - 1) - l]
                    for m in range(len(result)):
                        extra_minus = 1
                        for n in range(len(position_of_non_I_gate)):
                            if (list(result.keys()))[m][position_of_non_I_gate[n]] == '1':
                                extra_minus *= -1
                        expectation_value += result[list(result.keys())[m]]*extra_minus
                    term2[l] = expectation_value

                # print('term2:', term2)
                
                
                for l in range(super_circuit[i].num_parameters):
                    for m in range(l, super_circuit[i].num_parameters):
                        if l == m:
                            fubini_study_metric[l + parameter_previous, l + parameter_previous] = (1 - term2[l]**2)/4
                        else:

                            term1_lm = 0
                            position_of_non_I_gate = [(ansatz.num_qubits - 1) - l, (ansatz.num_qubits - 1) - m]
                            for h in range(len(result)):
                                extra_minus = 1
                                for n in range(len(position_of_non_I_gate)):
                                    if (list(result.keys()))[h][position_of_non_I_gate[n]] == '1':
                                        extra_minus *= -1
                                term1_lm += result[list(result.keys())[h]]*extra_minus
                            
                            fubini_study_metric[l + parameter_previous, m + parameter_previous] = (term1_lm - term2[l]*term2[m])/4
                            fubini_study_metric[m + parameter_previous, l + parameter_previous] = fubini_study_metric[l + parameter_previous, m + parameter_previous]
                # print(fubini_study_metric)
                # print(internal_circuit)
                        


    
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

        combine = learning_rate*FS_metric_inv.dot(grad)

        internal_initial_point =  np.subtract(internal_initial_point, combine)
    
    if callback is None:
        return energy


def Customize_SPSA(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """
    
    internal_initial_point = initial_point.copy()

    energy = []

    for k in range(interation):        
        internal_energy =  Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        #print(f'{internal_initial_point} ---------')
        #print(internal_energy)
        energy.append(internal_energy)
        
        if callback is not None:
            callback(internal_initial_point, internal_energy)

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
    
    return energy



def Customize_QNSPSA_PRS_blocking(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    beta = 0.001

    internal_initial_point = initial_point.copy()

    history_length = 5

    last_n_steps = np.zeros(history_length)

    energy = []

    regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
    grad = np.zeros(ansatz.num_parameters)
    # first measurement outside
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)
    #print(internal_energy)
    #print(f'{internal_initial_point} ---------')
    if callback is not None:
        callback(internal_initial_point, internal_energy)
    
    last_n_steps[0] = internal_energy

    for k in range(interation-1):   
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
            ck = 0.01
            ak = learning_rate
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

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2)[0] - SwapTest(ansatz_initial, ansatz_plus1)[0] - SwapTest(ansatz_initial, ansatz_minus1_plus2)[0] + SwapTest(ansatz_initial, ansatz_minus1)[0]


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
            tolerance = 2 * last_n_steps.std() if (k > history_length) else 2 * last_n_steps[:k+1].std()
            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            #print(f'tolerance: {tolerance}')
            # print('gradPRS: ', gradPRS)
            #print(next_energy)

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
            callback(internal_initial_point, internal_energy)

        last_n_steps[(k+1) % history_length] = internal_energy
    return energy


def Customize_QN_SPSA_blocking(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    beta = 0.001
    
    internal_initial_point = initial_point.copy()

    history_length = 5

    last_n_steps = np.zeros(history_length)

    energy = []

    

    super_circuit = Separate_Circuit_Apart(ansatz)

    regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
    grad = np.zeros(ansatz.num_parameters)
    # first measurement outside
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)
    #print(internal_energy)
    #print(f'{internal_initial_point} ---------')
    if callback is not None:
        callback(internal_initial_point, internal_energy)
    
    last_n_steps[0] = internal_energy

    for k in range(interation-1):   
        next_energy = 0
        fubini_study_metric = np.zeros((ansatz.num_parameters, ansatz.num_parameters))

        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                g_internal = [[[] for _ in range(super_circuit[i].num_parameters)] for _ in range(super_circuit[i].num_parameters)]
                internal_circuit = QuantumCircuit(super_circuit[i].num_qubits)

                parameter_previous = 0
                for j in range(i):
                    parameter_previous += super_circuit[j].num_parameters
                    internal_circuit = internal_circuit.compose(super_circuit[j])
                

                internal_circuit = internal_circuit.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                
                # Add string
                string = ['I']*ansatz.num_qubits

                for l in range(super_circuit[i].num_parameters):
                    if super_circuit[i][l].operation.name == 'rx':
                        string[l] = 'X'
                    if super_circuit[i][l].operation.name == 'ry':
                        string[l] = 'Y'
                    if super_circuit[i][l].operation.name == 'rz':
                        string[l] = 'Z'

                quantum_circuit_to_meausurement = QuantumCircuit(ansatz.num_qubits)

                for l in range(len(string)):
                    if string[l] == "X":
                        quantum_circuit_to_meausurement.h(l)
                    if string[l] == "Y":
                        quantum_circuit_to_meausurement.rx(np.pi/2, l)
                    
                
                quantum_circuit = internal_circuit.compose(quantum_circuit_to_meausurement)
                quantum_circuit.measure_all()
                # print('string: ', string)
                # print(quantum_circuit)
                job = sampler.run(quantum_circuit, shots = shots)

                result = (job.result().quasi_dists[0].binary_probabilities())
                # print(result)
                # term 2 in fubini-study
                term2 = np.zeros(len(string))
                for l in range(super_circuit[i].num_parameters):
                    expectation_value = 0
                    position_of_non_I_gate = [(ansatz.num_qubits - 1) - l]
                    for m in range(len(result)):
                        extra_minus = 1
                        for n in range(len(position_of_non_I_gate)):
                            if (list(result.keys()))[m][position_of_non_I_gate[n]] == '1':
                                extra_minus *= -1
                        expectation_value += result[list(result.keys())[m]]*extra_minus
                    term2[l] = expectation_value

                # print('term2:', term2)
                
                
                for l in range(super_circuit[i].num_parameters):
                    for m in range(l, super_circuit[i].num_parameters):
                        if l == m:
                            fubini_study_metric[l + parameter_previous, l + parameter_previous] = (1 - term2[l]**2)/4
                        else:

                            term1_lm = 0
                            position_of_non_I_gate = [(ansatz.num_qubits - 1) - l, (ansatz.num_qubits - 1) - m]
                            for h in range(len(result)):
                                extra_minus = 1
                                for n in range(len(position_of_non_I_gate)):
                                    if (list(result.keys()))[h][position_of_non_I_gate[n]] == '1':
                                        extra_minus *= -1
                                term1_lm += result[list(result.keys())[h]]*extra_minus
                            
                            fubini_study_metric[l + parameter_previous, m + parameter_previous] = (term1_lm - term2[l]*term2[m])/4
                            fubini_study_metric[m + parameter_previous, l + parameter_previous] = fubini_study_metric[l + parameter_previous, m + parameter_previous]
       
        while True:    
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

            
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            
            tolerance = 2 * last_n_steps.std() if (k > history_length) else 2 * last_n_steps[:k+1].std()

            

            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            #print(f'tolerance: {tolerance}')
            #print(next_energy)

            if next_energy <= internal_energy + tolerance:
                break
            # if True:
            #     break            

        
        # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()
        
        # ghi data
        internal_initial_point = np.subtract(internal_initial_point, grad)
        internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
        energy.append(internal_energy)
        #print(internal_energy)
        #print(f'{internal_initial_point} ---------')
        # print(f'gradSPSA {gradSPSA}')
        # print(f'regularized_fubini_matrix: {np.linalg.pinv(regularized_fubini_matrix)}')

        if callback is not None:
            callback(internal_initial_point, internal_energy)

        last_n_steps[(k+1) % history_length] = internal_energy
    return energy

def Customize_QNSPSA_SPSA_blocking(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    beta = 0.001
    
    internal_initial_point = initial_point.copy()

    history_length = 5

    last_n_steps = np.zeros(history_length)

    energy = []

    regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
    grad = np.zeros(ansatz.num_parameters)
    # first measurement outside
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)
    #print(internal_energy)
    #print(f'{internal_initial_point} ---------')
    if callback is not None:
        callback(internal_initial_point, internal_energy)
    
    last_n_steps[0] = internal_energy

    for k in range(interation-1):   
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

            deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2)[0] - SwapTest(ansatz_initial, ansatz_plus1)[0] -SwapTest(ansatz_initial, ansatz_minus1_plus2)[0] + SwapTest(ansatz_initial, ansatz_minus1)[0]


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
            
            tolerance = 2 * last_n_steps.std() if (k > history_length) else 2 * last_n_steps[:k+1].std()

            #print(f'tolerance: {tolerance}')

            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            if next_energy < internal_energy + tolerance:
                break            

        
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
            callback(internal_initial_point, internal_energy)

        last_n_steps[(k+1) % history_length] = internal_energy
    return energy


def Customize_HESPSA_SPSA_blocking(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    beta = 0.001
    
    internal_initial_point = initial_point.copy()

    history_length = 5

    last_n_steps = np.zeros(history_length)

    energy = []

    regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))
    grad = np.zeros(ansatz.num_parameters)
    # first measurement outside
    internal_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)
    energy.append(internal_energy)
    print(internal_energy)
    print(f'{internal_initial_point} ---------')
    if callback is not None:
        callback(internal_initial_point, internal_energy)
    
    last_n_steps[0] = internal_energy

    for k in range(interation-1):   
        next_energy = 0

        while True:    
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

            deltaF = Transverse_Ising_Measurement(operator, ansatz_plus1_plus2, shots, sampler) - Transverse_Ising_Measurement(operator, ansatz_plus1, shots, sampler) - Transverse_Ising_Measurement(operator, ansatz_minus1_plus2, shots, sampler) + Transverse_Ising_Measurement(operator, ansatz_minus1, shots, sampler)


            # print(deltaF)
            # Fubini_matrix
            fubini_matrix = (deltaF/(4*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))
            
            # Data of previous regularized study metric
            exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*fubini_matrix.copy()    
            regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

            # regularized_fubini_matrix_previous = regularized_fubini_matrix.copy()

            grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradSPSA)

            
            # cập nhật tham số
            internal_initial_point_while = np.subtract(internal_initial_point, grad)
            
            tolerance = 2 * last_n_steps.std() if (k > history_length) else 2 * last_n_steps[:k+1].std()

            #print(f'tolerance: {tolerance}')

            next_energy = Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: internal_initial_point_while[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler)

            if next_energy < internal_energy + tolerance:
                break            

            # if True:
            #     break

        
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
            callback(internal_initial_point, internal_energy)

        last_n_steps[(k+1) % history_length] = internal_energy
    return energy

def Customize_HE_PRS(operator, initial_point, learning_rate, ansatz, interation, shots, callback, sampler):
    """
    operator: the pauli matrix
    interation: number of interations 
    ansatz: the ansatz for perparing the parameterized circuit
    eta: learning rate
    initial_points: random sample of number uses at the beggining of running
    shots: number of shots
    backend: the backend
    """

    energy = []

    internal_initial_point = initial_point.copy()
    # regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))

    for k in range(interation):
        print(f'{internal_initial_point} ---------')
        internal_ansatz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = Transverse_Ising_Measurement(operator, internal_ansatz, shots, sampler)
        
        energy.append(internal_energy)
        print(internal_energy)


        if callback is not None:
            callback(internal_initial_point, internal_energy)

        # Update the parameter points
        grad = np.zeros(ansatz.num_parameters)
        Hessian_Matrix = np.zeros((ansatz.num_parameters, ansatz.num_parameters))

        for i in range(ansatz.num_parameters):
            plus_parameter = internal_initial_point.copy()
            plus_parameter[i] += np.pi/2
            minus_parameter = internal_initial_point.copy()
            minus_parameter[i] -= np.pi/2
            grad[i] = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler))/2
            Hessian_Matrix[i,i] = (Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) + Transverse_Ising_Measurement(operator, ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, sampler) - 2*internal_energy)/2
        
        # exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*Hessian_Matrix.copy()    
        # #     regularized_fubini_matrix = np.add(matrix_power(np.dot(exponentially_smoothed_fubini,exponentially_smoothed_fubini), 1/2).real, beta*np.identity(ansatz.num_parameters))

        # regularized_fubini_matrix_previous = exponentially_smoothed_fubini.copy()

        a = learning_rate*np.linalg.pinv(Hessian_Matrix).dot(grad)


        internal_initial_point =  np.subtract(internal_initial_point,a)
        # print(f'a: {a}')

        # print(internal_energy)

    if callback is None:
        return energy

