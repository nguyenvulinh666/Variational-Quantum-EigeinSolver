import numpy as np
from qiskit.execute_function import execute
from qiskit import BasicAer
import qiskit
from qiskit.circuit import Parameter
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, Aer, execute
import re
import Magicfile
from qiskit.opflow import Z, I, X



def Measurement(operator, quantum_circuit, shots, backend):
    """
    operator: the Pauli matrix we want to measure for the Quantum Nature Gradient Descent with writing in this from (pauli_matrix+is_position). For exammple: (Z1Z2) or (X3)
    quantum_circuit: the circuit for measurement 
    shots: the number of shots that we want to executee
    !!!COMMENT: In my way of measurement, I reverse the qiskit-way of reprent a arbitrary state for more easily handle in my code, if I don't my code would be meesy and changing alot, and more I want measure the same as the we repsent them like XZZ to be X o Z o Z not be reversed like the qiskit does hmmmmmmmmmmmm
    """
    number_of_qubits = quantum_circuit.num_qubits
    

    if len(operator) > (number_of_qubits):
        raise Exception("The position of Pauli matrix exceeds the number of qubits")

    quantum_circuit_to_measurement = QuantumCircuit(number_of_qubits)


    expectation_value = 0

    for i in range(len(operator)):
        if operator[i] == 'X':
            quantum_circuit_to_measurement.ry(-np.pi/2, i)
        if operator[i] == 'Y':
            quantum_circuit_to_measurement.rx(np.pi/2, i)

    new_quantum_circuit_to_measurement = quantum_circuit.compose(quantum_circuit_to_measurement)
    new_quantum_circuit_to_measurement.measure_all()


    job = execute(new_quantum_circuit_to_measurement, backend, shots=shots)
    result = job.result()
    counts = result.get_counts()

    # Reversed the way qiskit represnts state
    counts_change =  {(key)[::-1] : value for key, value in counts.items()}


    position_of_non_I_gate = np.where(np.array(list(operator)) != 'I')[0]

    for i in range(len(counts_change)):
        extra_minus = 1
        for j in range(len(position_of_non_I_gate)):
            if list(counts_change.keys())[i][position_of_non_I_gate[j]] == '1':
                extra_minus *= -1       
        expectation_value += counts_change[list(counts_change.keys())[i]]*extra_minus
        

    expectation_value = expectation_value/shots

    return expectation_value, counts_change, new_quantum_circuit_to_measurement


def Customize_EfficientSU2(number_qubits, number_of_subcircuit):
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
        theta1 = Parameter(r'$\theta[{}]$'.format(0))
        circuit.ry(theta1, 0)
        theta2 = Parameter(r'$\theta[{}]$'.format(1))
        circuit.rz(theta2, 0)
        circuit.barrier()

        for i in range(number_of_subcircuit):
            theta1 = Parameter(r'$\theta[{}]$'.format(2*i+2))
            circuit.ry(theta1, 0)
            theta2 = Parameter(r'$\theta[{}]$'.format(2*i+3))
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
        sub_circuit.barrier()


        # hmm cause of lack of my knowledge, I will work with the RealAmplitudes ansat, circular entanglement
        for i in range(0, number_qubits):
            theta = Parameter(r'$\theta[{}]$'.format(number_parameter+i))
            sub_circuit.ry(theta, i)
        for i in range(0, number_qubits):
            theta = Parameter(r'$\theta[{}]$'.format(number_qubits+number_parameter+i))
            sub_circuit.rz(theta, i)

        return sub_circuit


    for i in range(0, number_qubits):
        theta = Parameter(r'$\theta[{}]$'.format(i))
        circuit.ry(theta, i)

    for i in range(0, number_qubits):
        theta = Parameter(r'$\theta[{}]$'.format(number_qubits+i))
        circuit.rz(theta, i)

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

        sub_circuit.cx(number_qubits-1, 0)
        for i in (range(number_qubits-1)):
            sub_circuit.cx(i, i+1)

        sub_circuit.barrier()


        # hmm cause of lack of my knowledge, I will work with the RealAmplitudes ansat, circular entanglement
        for i in range(0, number_qubits):
            theta = Parameter(r'$\theta[{}]$'.format(number_parameter+i))
            sub_circuit.ry(theta, i)

        return sub_circuit


    for i in range(0, number_qubits):
        theta = Parameter(r'$\theta[{}]$'.format(i))
        circuit.ry(theta, i)

    for i in range(number_of_subcircuit):
        circuit = circuit.compose(add_subcircuit(circuit, i))

    return circuit


import numpy as np
from scipy.linalg import block_diag

# Parameter shift rule
def Customize_Parameter_Shift_Rule(operator, parameter, eta, ansatz, shots, backend, interation):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    eta: learning rate
    """
    number_of_parameters = len(parameter)
    grad = np.zeros(number_of_parameters)
    

    for i in range(len(operator)):
        grad_func = np.zeros(number_of_parameters)
        for j in range(number_of_parameters):
            plus_parameter = parameter.copy()
            plus_parameter[j] += np.pi/2
            minus_parameter = parameter.copy()
            minus_parameter[j] -= np.pi/2
            grad_func[j] = (eta)*operator.coeffs.real[i]*(Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0] - Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0])/2
        grad = np.add(grad, grad_func)
    
    return np.array(grad)


# Gradient
def Customize_Gradient(operator, parameter, eta, ansatz, shots, backend, interation):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    eta: learning rate
    """
    number_of_parameters = len(parameter)
    grad = np.zeros(number_of_parameters)
    h = .01

    for i in range(len(operator)):
        grad_func = np.zeros(number_of_parameters)
        for j in range(number_of_parameters):
            plus_parameter = parameter.copy()
            plus_parameter[j] += h
            minus_parameter = parameter.copy()
            minus_parameter[j] -= h
            grad_func[j] = (eta)*operator.coeffs.real[i]*(Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: plus_parameter[i] for i,
            theta in enumerate(ansatz.parameters)}), shots, backend)[0] - Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0])/(2*h)
        grad = np.add(grad, grad_func)
    
    return np.array(grad)

def Seperate_Circuit_Apart(ansantz):
    # Divide the circuit to subcircuit of parameter circuit and non parameter circuit
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




def Customize_Quantum_Natural_Graident_Descent(operator, parameter, eta, ansatz, shots, backend, interation, super_circuit):
    """
    operator: The pauli operator
    parameter: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """

    number_of_parameters = ansatz.num_parameters

    # return parameter_shift_rules(operator, parameter, ansatz)

    # Next, I will calculate the Fubini-Study metric, but first I would like to learn the way of extracting information from the ansatz, and calculate the Fubini-Study metric
    def fubini_stydy_metric(ansatz, parameter, shots, backend):
        def Measure_element_of_Fubini_Study_metric(circuit, circuit_for_measurement, i, j, shots, backend):
            if i != j:
                term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)
                term3 = ['I']*len(circuit)
                
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
                

                term1 = ''.join(term1)
                term2 = ''.join(term2)
                term3 = ''.join(term3)

                return  (Measurement(term1, circuit_for_measurement, shots, backend)[0] - Measurement(term2, circuit_for_measurement, shots, backend)[0]*Measurement(term3, circuit_for_measurement, shots, backend)[0])/4
            else:
                # term1 = ['I']*len(circuit)
                term2 = ['I']*len(circuit)

                if circuit[i].operation.name == 'rx':
                    term2[i] = 'X'
                if circuit[i].operation.name == 'ry':
                    term2[i] = 'Y'
                if circuit[i].operation.name == 'rz':
                    term2[i] = 'Z'

                term2 = ''.join(term2)

                # new_gate = circuit_for_measurement.compose(new_gate)

                # return (Measurement(term1, circuit_for_measurement, shots, backend)[0] - Measurement(term2, circuit_for_measurement, shots, backend)[0].real**2)/4
                return (1 - Measurement(term2, circuit_for_measurement, shots, backend)[0].real**2)/4


        # Caculate the Fubini-Study metric in 1-qubit
        if ansatz.num_qubits == 1:
            initial_point = parameter
            # num_parameter = 0
            g = np.zeros((number_of_parameters, number_of_parameters))

            for i in range(len(super_circuit)):
                if super_circuit[i][0].operation.params:
                    internal_circuit = QuantumCircuit(super_circuit[i].num_qubits)
                    for j in range(i+1):
                        internal_circuit = internal_circuit.compose(super_circuit[j])

                    internal_circuit = internal_circuit.bind_parameters({theta: initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                    # num_parameter += super_circuit[i].num_parameters

                    g[i][i] = Measure_element_of_Fubini_Study_metric(super_circuit[i], internal_circuit, 0, 0, shots, backend)
                    # print(np.round(g_internal,2))
            return g

        
        # Caculate the Fubini_study_metric
        g = []
        initial_point = parameter
        # num_parameter = 0


        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                g_internal = [[[] for _ in range(super_circuit[i].num_parameters)] for _ in range(super_circuit[i].num_parameters)]
                internal_circuit = QuantumCircuit(super_circuit[i].num_qubits)
                for j in range(i+1):
                    internal_circuit = internal_circuit.compose(super_circuit[j])

                internal_circuit = internal_circuit.bind_parameters({theta: initial_point[i] for i, theta in enumerate(internal_circuit.parameters)})
                # num_parameter += super_circuit[i].num_parameters

                # print(num_parameter)
                
                for l in range(super_circuit[i].num_parameters):
                    for m in range(super_circuit[i].num_parameters):
                        g_internal[l][m] = Measure_element_of_Fubini_Study_metric(super_circuit[i], internal_circuit, l, m, shots, backend)
                # print(np.round(g_internal,2))
                if g != []:
                    g = block_diag(g, g_internal)
                else:
                    g = block_diag(g_internal)
        # print(g)
        return g
    

    # Gradient of the loss function
    grad = np.zeros(len(parameter))
    
    grad = Customize_Parameter_Shift_Rule(operator, parameter, eta, ansatz, shots, backend, interation)

    
    # grad = gradient(operator, parameter, ansatz)

    # Inverse of Fubini study metric
    FS_metric = np.linalg.pinv(fubini_stydy_metric(ansatz, parameter, shots, backend)).real
    # print(fubini_stydy_metric(ansatz, parameter, shots, backend))

    # print(FS_metric)

    combine = FS_metric.dot(grad)

    # print(f'Fubini-Study metric: {FS_metric}')

    # print(f'combine: {combine}')
    return combine




def  Customize_imaginary_time_evolution(operator, interation, ansatz, eta, initial_point, shots, backend, callback):
    """
    operator: the pauli matrix
    interation: number of interations 
    ansatz: the ansatz for perparing the parameterized circuit
    eta: learning rate but in this case is the time interval of the imaginary time evolution
    initial_points: random sample of number uses at the beggining of running
    shots: number of shots
    backend: the backend
    """

    circ_U = [[[] for _ in range(ansatz.num_parameters)] for _ in range(ansatz.num_parameters)]

    super_circuit = Magicfile.Seperate_Circuit_Apart(ansatz)

    # return the circuit for perpare the evaluate the value of matrix U  
    def circ_for_U_part(ansatz, l, m):
        layer_of_l = 0
        position_of_l = 0
        layer_of_m = 0
        position_of_m = 0


        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                # data of the parrameter layer
                a = str(super_circuit[i].data[-1].operation.params[0])
                b = int(re.findall(r'\d+', a)[0])

                if l <= b:
                    layer_of_l = i
                    position_of_l = ansatz.num_qubits - 1 - (b - l)
                    break
                    

        for i in range(len(super_circuit)):
            if super_circuit[i][0].operation.params:
                # data of the parrameter layer
                a = str(super_circuit[i].data[-1].operation.params[0])
                b = int(re.findall(r'\d+', a)[0])

                if m <= b:
                    layer_of_m = i
                    position_of_m = ansatz.num_qubits - 1 - (b - m)
                    break
        


        circ = QuantumCircuit(ansatz.num_qubits+1)
        ancilla_circ = QuantumCircuit(1)
        ancilla_circ.h(0)
        ancilla_circ.x(0)
        circ = circ.compose(ancilla_circ)


        h = 0
        position_h = 0
            
        if l > m:
            h = l
            position_h = position_of_l 
            l = m
            position_of_l = position_of_m
            m = h
            position_of_m = position_h
            

        if l <= m:
            internal_circuit = QuantumCircuit(ansatz.num_qubits)
            for i in range(layer_of_l):
                internal_circuit = internal_circuit.compose(super_circuit[i])
            
            circ = circ.compose(internal_circuit, [i for i in range(1, ansatz.num_qubits+1)])

            circ.barrier()

            if super_circuit[layer_of_l].data[position_of_l].operation.name == 'rx':
                circ.cx(0, position_of_l+1)
            if super_circuit[layer_of_l].data[position_of_l].operation.name == 'ry':
                circ.cy(0, position_of_l+1)
            if super_circuit[layer_of_l].data[position_of_l].operation.name == 'rz':
                circ.cz(0, position_of_l+1)

            circ.barrier()

            circ.x(0)

            internal_circuit = QuantumCircuit(ansatz.num_qubits)
            for i in range(layer_of_l, layer_of_m):
                internal_circuit = internal_circuit.compose(super_circuit[i])

            circ = circ.compose(internal_circuit, [i for i in range(1, ansatz.num_qubits+1)])
            
            circ.barrier()

            if super_circuit[layer_of_m].data[position_of_m].operation.name == 'rx':
                circ.cx(0, position_of_m+1)
            if super_circuit[layer_of_m].data[position_of_m].operation.name == 'ry':
                circ.cy(0, position_of_m+1)
            if super_circuit[layer_of_m].data[position_of_m].operation.name == 'rz':
                circ.cz(0, position_of_m+1)

            circ.barrier()

            circ.h(0)   

            measurement = ClassicalRegister(1)
            circ.add_register(measurement)
            circ.measure(0, 0)

        return circ


    def circ_for_C_part(ansatz, l):
        layer_of_l = 0
        position_of_l = 0

        for i in range(len(super_circuit)):   
            if super_circuit[i][0].operation.params:
                # data of the parrameter layer
                a = str(super_circuit[i].data[-1].operation.params[0])
                # b = 0

                # for num in range(len(a)):
                #     if a[num].isnumeric():
                #         b = int(a[num])
                b = int(re.findall(r'\d+', a)[0])

                if l <= b:
                    layer_of_l = i
                    position_of_l = ansatz.num_qubits - 1 - (b - l)
                    # print(f'layer_of_m: {layer_of_m}')
                    # print(f'position_of_m: {position_of_m}')
                    break
            
        circ = QuantumCircuit(ansatz.num_qubits+1)
        ancilla_circ = QuantumCircuit(1)
        ancilla_circ.h(0)
        ancilla_circ.s(0)
        ancilla_circ.x(0)
        circ = circ.compose(ancilla_circ)

        internal_circuit = QuantumCircuit(ansatz.num_qubits)
        for i in range(layer_of_l):
            internal_circuit = internal_circuit.compose(super_circuit[i])
        
        circ = circ.compose(internal_circuit, [i for i in range(1, ansatz.num_qubits+1)])

        circ.barrier()

        if super_circuit[layer_of_l].data[position_of_l].operation.name == 'rx':
            circ.cx(0, position_of_l+1)
        if super_circuit[layer_of_l].data[position_of_l].operation.name == 'ry':
            circ.cy(0, position_of_l+1)
        if super_circuit[layer_of_l].data[position_of_l].operation.name == 'rz':
            circ.cz(0, position_of_l+1)

        circ.barrier()

        circ.x(0)

        internal_circuit = QuantumCircuit(ansatz.num_qubits)
        for i in range(layer_of_l, len(super_circuit)):
            internal_circuit = internal_circuit.compose(super_circuit[i])

        circ = circ.compose(internal_circuit, [i for i in range(1, ansatz.num_qubits+1)])
        
        circ.barrier()
        return circ

    energy = []

    internal_initial_point = initial_point.copy()


    for i in range(ansatz.num_parameters):
        circ_U[i][i] = circ_for_U_part(ansatz, i, i)


    for i in range(ansatz.num_parameters):
        for j in range(i+1, ansatz.num_parameters):
            circ_U[i][j] = circ_for_U_part(ansatz, i, j)
            circ_U[j][i] = circ_for_U_part(ansatz, i, j)

    for i in range(interation):
        print(f'{internal_initial_point} ---------')
        internal_anstaz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = 0
        for j in range(len(operator)):
            internal_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
        
        energy.append(internal_energy)
        print(internal_energy)


        

        # The U_matrix
        U_matrix = np.zeros((ansatz.num_parameters, ansatz.num_parameters))

        for element in range(ansatz.num_parameters):
            U_matrix[element][element] = 1

        
        for element in range(ansatz.num_parameters):
            for j in range(element+1, ansatz.num_parameters):
                internal_ansatz = circ_U[element][j].copy()
                internal_ansatz = internal_ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(internal_ansatz.parameters)})  
                
                result = execute(internal_ansatz, backend, shots=shots).result()

                counts = result.get_counts()
                overlap_value = 0
                
                if '0' in counts:
                    # Key exists, do something with data['value']
                    overlap_value = counts['0']
                else:
                    # Key doesn't exist, return to default value (0 in this case)
                    overlap_value = 0
                
                U_matrix[element][j] = overlap_value/shots - 1
                U_matrix[j][element] = overlap_value/shots - 1


        # Matrix U
        U_matrix = np.array(U_matrix/4)
        grad_vector_C = np.zeros(ansatz.num_parameters)
    

        for element in range(len(operator)):
            grad_func = np.zeros(ansatz.num_parameters)
            for j in range(ansatz.num_parameters):
                plus_parameter = internal_initial_point.copy()
                plus_parameter[j] += np.pi/2
                minus_parameter = internal_initial_point.copy()
                minus_parameter[j] -= np.pi/2
                grad_func[j] = (1)*operator.coeffs.real[element]*(Measurement(str(operator.primitive._pauli_list[element]), ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0] - Measurement(str(operator.primitive._pauli_list[element]), ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0])/2
            grad_vector_C = np.subtract(grad_vector_C, grad_func)      
        
        print(f'grad_vector_C: {grad_vector_C}')
        print(f'inverse of U matrix: {np.linalg.pinv(U_matrix)} ')
        internal_initial_point = np.add(internal_initial_point, eta*np.dot(np.linalg.pinv(U_matrix), grad_vector_C))
        # print(eta*np.dot(np.linalg.pinv(U_matrix), grad_vector_C))
        # internal_initial_point = np.add(internal_initial_point, eta*grad_vector_C)
    if callback is None:
        return energy

        

def Customize_Nelder_Mead(operator, interation, ansatz, initial_point, shots, backend, callback):
    """
    operator: the pauli matrix
    interation: number of interations 
    ansatz: the ansatz for perparing the parameterized circuit
    initial_points: random sample of number uses at the beggining of running
    shots: number of shots
    backend: the backend
    """
    def bubbleSort(arr):
        n = len(arr)
        
        # Traverse through all array elements
        for i in range(n):
            swapped = False
    
            # Last i elements are already in place
            for j in range(0, n-i-1):
    
                # Traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if arr[j] > arr[j+1]:
                    arr[j], arr[j+1] = arr[j+1], arr[j]
                    swapped = True
            if (swapped == False):
                break
        return arr

    energy = []

    data_points = [[] for _ in range(ansatz.num_parameters+1)]
    data_energy = [[] for _ in range(ansatz.num_parameters+1)]

    data_points[0] = initial_point

    # data_points[0] = np.array([1.4, 5.5])
    # data_points[1] = np.array([1.6, 5.5])
    # data_points[2] = np.array([1.5, 6])

    for i in range(1, ansatz.num_parameters+1):
        # np.random.seed(i)

        # random = np.random.rand(ansatz.num_parameters)*2*np.pi
        # random = np.ones(ansatz.num_parameters)
        
        # data_points[i] = np.subtract(initial_point, random)
        # data_points[i][i-1]  = np.pi*2
        if i % 2 == 0:
            random = np.ones(ansatz.num_parameters)*2*np.pi*i/ansatz.num_parameters
            data_points[i] = np.subtract(initial_point, random)
        else:
            random = np.ones(ansatz.num_parameters)*(-2)*np.pi*i/ansatz.num_parameters
            data_points[i] = np.subtract(initial_point, random)
        
        

    
    for i in range(ansatz.num_parameters+1):
        # Measurement energy
        internal_anstaz = ansatz.bind_parameters({theta: data_points[i][k] for k, theta in enumerate(ansatz.parameters)})       
        internal_energy = 0
        for j in range(len(operator)):
            internal_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
        data_energy[i] = internal_energy
    

    for i in range(interation):
        
        # Sort
        data_energy_sort = bubbleSort(data_energy.copy())
        data_points_sort = [[] for _ in range(ansatz.num_parameters+1)]

        for j in range(len(data_energy_sort)):
            data_points_sort[j] = data_points[np.where(data_energy == data_energy_sort[j])[0][0]]

        data_energy = data_energy_sort
        data_points = data_points_sort
    
        avarage_energy = sum(data_energy[i] for i in range(len(data_energy)))/len(data_energy)
        energy.append(avarage_energy)
        
        print(avarage_energy)
        print(f'data_points: {data_points}')
        print(data_energy)

        # Nelder Mead part
        # Calculate the centroid point
        centroid = np.zeros(ansatz.num_parameters)

        for i in range(ansatz.num_parameters):
            for j in range(ansatz.num_parameters):
                centroid[i] += data_points[j][i]

        centroid = centroid/ansatz.num_parameters 
        # print(centroid)

        # Reflection
        reflection_coefficient = 1

        reflection_point = np.zeros(ansatz.num_parameters)

        reflection_point = centroid + reflection_coefficient*(centroid - np.array(data_points[-1]))
        # print(f'reflection_point: {reflection_point}')
        internal_anstaz = ansatz.bind_parameters({theta: reflection_point[k] for k, theta in enumerate(ansatz.parameters)})       
        
        # Measure the reflection energy
        reflection_energy = 0
        for j in range(len(operator)):
            reflection_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
        
        if data_energy[0] <= reflection_energy and reflection_energy < data_energy[-2]:
            print('reflection')
            data_points[-1] = reflection_point
            data_energy[-1] = reflection_energy
            continue
        

        # Expansion
        expansion_coefficient = 1

        if reflection_energy < data_energy[0]:
            print('expansion')
            expansion_point = np.zeros(ansatz.num_parameters)

            for j in range(ansatz.num_parameters):
                expansion_point[j] = centroid[j] + expansion_coefficient*(reflection_point[j] - centroid[j])

            internal_anstaz = ansatz.bind_parameters({theta: expansion_point[k] for k, theta in enumerate(ansatz.parameters)})       
            # Measure the expansion energy
            expansion_energy = 0
            for j in range(len(operator)):
                expansion_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]

            if expansion_energy < reflection_energy:
                data_points[-1] = expansion_point
                data_energy[-1] = expansion_energy
                continue
            
            else:
                data_points[-1] = reflection_point
                data_energy[-1] = reflection_energy
                continue
        
        # Contraction
        contraction_coefficient = 0.5

        if reflection_energy < data_energy[-1]:
            contraction_point = np.zeros(ansatz.num_parameters)

            for j in range(ansatz.num_parameters):
                contraction_point[j] = centroid[j] + contraction_coefficient*(reflection_point[j] - centroid[j])

            internal_anstaz = ansatz.bind_parameters({theta: contraction_point[k] for k, theta in enumerate(ansatz.parameters)})       
            # Measure the contraction energy
            contraction_energy = 0
            for j in range(len(operator)):
                contraction_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]

            if contraction_energy < reflection_energy:
                print('contraction1')
                data_points[-1] = contraction_point
                data_energy[-1] = contraction_energy
                continue
        
        else:
            contraction_point = np.zeros(ansatz.num_parameters)

            for j in range(ansatz.num_parameters):
                contraction_point[j] = centroid[j] + contraction_coefficient*(data_points[-1][j] - centroid[j])

            internal_anstaz = ansatz.bind_parameters({theta: contraction_point[k] for k, theta in enumerate(ansatz.parameters)})       
            # Measure the contraction energy
            contraction_energy = 0
            for j in range(len(operator)):
                contraction_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]

            if contraction_energy < data_energy[-1]:
                print('contraction2')
                data_points[-1] = contraction_point
                data_energy[-1] = contraction_energy
                continue
        
        # Shrink
        shrink_coefficient = 0.5
        print('shrink')
        
        for k in range(1, ansatz.num_parameters+1):
            shrink_point = np.zeros(ansatz.num_parameters)
            
            for j in range(ansatz.num_parameters):    
                shrink_point[j] = data_points[0][j] + shrink_coefficient*(data_points[k][j] - data_points[0][j])

            internal_anstaz = ansatz.bind_parameters({theta: shrink_point[k] for k, theta in enumerate(ansatz.parameters)})       
            
            internal_energy = 0
            # shink energy
            for j in range(len(operator)):
                internal_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
            
            data_points[k] = shrink_point
            data_energy[k] = internal_energy



    return energy


def pauli_measurement(operator, ansatz, initial_point, shots, backend):
    """ 
    operator: indicate the sum of serries of Tensor Product Pauli Operator (Hamiltonian)
    ansatz: the parameter circuit
    initial_point: the parameter for our parameter circuit
    shots: number of shots
    backend: the quantum simulator
    retun: the Hamiltonian value
    """
    internal_anstaz = ansatz.bind_parameters({theta: initial_point[k] for k, theta in enumerate(ansatz.parameters)})  

    internal_energy = 0
    for j in range(len(operator)):
        internal_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
    return internal_energy



def Customize_SPSA(operator, interation, ansatz, initial_point, shots, backend, callback):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """


    # Scaling parameter for the gradient step size.
    # c = np.array([pauli_measurement(operator, ansatz, initial_point, shots, backend) for _ in range(10)]).std()
    c = 1

    # Stability constant for the step size 10% of number of interations
    A = 0.1*interation

    # Scaling exponent for the step size
    alpha = 0.602

    # Scaling exponent for the gradient step size
    gamma = 0.101

    #  Smallest change of magnitude desired in the elements of θ.
    delTheta = 0.05
    
    
    magnitude_g0 = np.abs(Customize_Parameter_Shift_Rule(operator, initial_point, c, ansatz, shots, backend, interation=None)).mean()
    # Scaling parameter for the step size
    a = delTheta*((A+1)**alpha)/magnitude_g0
    
    internal_initial_point = initial_point.copy()

    energy = []

    for k in range(interation):        
        internal_energy = pauli_measurement(operator, ansatz, internal_initial_point, shots, backend)
        print(internal_energy)
        print(internal_initial_point)
        energy.append(internal_energy)


        # Optimization part
        grad = np.zeros(ansatz.num_parameters)
        ak = a/(1+k+A)**alpha
        ck = c/(1+k)**gamma
        
        print(f'ak_{k}: {ak}')
        print(f'ck_{k}: {ck}')
        for i in range(len(operator)):
            random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            # print(random)
            plus_parameter = np.array(internal_initial_point.copy())
            plus_parameter = np.add(plus_parameter,random*ck)
            minus_parameter = np.array(internal_initial_point.copy())
            minus_parameter = np.subtract(minus_parameter,random*ck)
            grad_func = ak*operator.coeffs.real[i]*(Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0] - Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0])/(2*ck)
            
            grad = np.add(grad, grad_func*random)
        print(f'grad: {grad}')
        internal_initial_point = np.subtract(internal_initial_point, grad)
    
    return energy


def SwapTest(circ1, circ2, shots, backend):
    """
    circ1: for the parameter circuit of VQD
    circ2: for the parameter circuit of VQE, this circuit have the infromation from the VQE running, or the k-1 of VQD running
    shots: number of shots
    backend: backend for running
    """

    # Perpare the circuit
    num_qubits = circ1.num_qubits
    circ = QuantumCircuit(2*num_qubits+1) # Plus one for the ancilla qubit

    circ = circ.compose(circ1, [i for i in range(1, num_qubits+1)])
    circ = circ.compose(circ2, [i for i in range(num_qubits+1, num_qubits*2+1)])
    circ.barrier()

    circ.h(0)

    for i in range(1, num_qubits+1):
        circ.cswap(0, i, i+num_qubits)

    circ.h(0) 


    # Make measurement
    measurement = ClassicalRegister(1)
    circ.add_register(measurement)
    circ.measure(0, 0)

    result = execute(circ, backend, shots=shots).result()

    counts = result.get_counts()

    overlap_value = 0

    if '1' in counts:
        # Key exists, do something with data['value']
        overlap_value = counts['1']
    else:
        # Key doesn't exist, return to default value (0 in this case)
        overlap_value = 0
    
    overlap_value = 1 - 2*overlap_value/shots

    

    
    return overlap_value, circ

def Customize_QN_SPSA(operator, interation, ansatz, initial_point, shots, backend, callback):
    """
    operator: The pauli operator
    interation: number of interation
    initial_point: the initial point that we will update until we end up with the desired point
    ansatz: the parameterized circuit that we want to update 
    """


    # Scaling parameter for the gradient step size.
    # c = np.array([pauli_measurement(operator, ansatz, initial_point, shots, backend) for _ in range(10)]).std()
    # c = 0.08
    c = 1

    # Stability constant for the step size 10% of number of interations
    A = 0.1*interation

    # Scaling exponent for the step size
    alpha = 0.602

    # Scaling exponent for the gradient step size
    gamma = 0.101

    #  Smallest change of magnitude desired in the elements of θ.
    # delTheta = 0.001
    delTheta = 0.05
    
    
    magnitude_g0 = np.abs(Customize_Parameter_Shift_Rule(operator, initial_point, c, ansatz, shots, backend, interation=None)).mean()
    # Scaling parameter for the step size
    a = delTheta*((A+1)**alpha)/magnitude_g0
    # a = 0.001

    #  regularization coefficient
    beta = 0.05
    
    internal_initial_point = initial_point.copy()

    energy = []

    regularized_fubini_matrix_previous = np.zeros((ansatz.num_parameters, ansatz.num_parameters))

    for k in range(interation):        
        internal_energy = pauli_measurement(operator, ansatz, internal_initial_point, shots, backend)
        print(internal_energy)
        print(internal_initial_point)
        energy.append(internal_energy)

        grad = np.zeros(ansatz.num_parameters)
        # SPSA part
        gradSPSA = np.zeros(ansatz.num_parameters)
        ak = a/(1+k+A)**alpha
        # ak = 0.001
        ck = c/(1+k)**gamma
        # ck = 0.05
        
        print(f'ak_{k}: {ak}')
        print(f'ck_{k}: {ck}')
        for i in range(len(operator)):
            random = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
            # print(random)
            plus_parameter = np.array(internal_initial_point.copy())
            plus_parameter = np.add(plus_parameter,random*ck)
            minus_parameter = np.array(internal_initial_point.copy())
            minus_parameter = np.subtract(minus_parameter,random*ck)
            grad_func = operator.coeffs.real[i]*(Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: plus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0] - Measurement(str(operator.primitive._pauli_list[i]), ansatz.bind_parameters({theta: minus_parameter[i] for i, theta in enumerate(ansatz.parameters)}), shots, backend)[0])/(2*ck)
            
            gradSPSA = np.add(gradSPSA, grad_func*random)
        
        print(f'gradSPSA {gradSPSA}')

        
        # Natural Gradient Part
        random1 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
        random2 = np.array([np.random.choice([-1,1]) for _ in range(ansatz.num_parameters)])
        initial_plus1_plus2 = np.add(initial_point, np.add(ck*random1, ck*random2))
        initial_plus1 = np.add(initial_point, ck*random1)
        initial_minus1_plus2 = np.subtract(initial_point, np.subtract(ck*random1, ck*random2))
        initial_minus1 = np.subtract(initial_point, ck*random1)

        ansatz_initial = ansatz.bind_parameters({theta: initial_point[k] for k, theta in enumerate(ansatz.parameters)})  
        ansatz_plus1_plus2 = ansatz.bind_parameters({theta: initial_plus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
        ansatz_plus1 = ansatz.bind_parameters({theta: initial_plus1[k] for k, theta in enumerate(ansatz.parameters)})  
        ansatz_minus1_plus2 = ansatz.bind_parameters({theta: initial_minus1_plus2[k] for k, theta in enumerate(ansatz.parameters)})  
        ansatz_minus1 =  ansatz.bind_parameters({theta: initial_minus1[k] for k, theta in enumerate(ansatz.parameters)})  

        deltaF = SwapTest(ansatz_initial, ansatz_plus1_plus2, shots, backend)[0] - SwapTest(ansatz_initial, ansatz_plus1, shots, backend)[0] -SwapTest(ansatz_initial, ansatz_minus1_plus2, shots, backend)[0] + SwapTest(ansatz_initial, ansatz_minus1, shots, backend)[0]

        # Fubini_matrix
        fubini_matrix = -1/2*(deltaF/(2*ck**2))*(np.array(np.array([random1]).T*random2) + np.array(np.array([random2]).T*random1))/2
        
        # Data of previous regularized study metric
        exponentially_smoothed_fubini = k/(k+1)*regularized_fubini_matrix_previous + 1/(k+1)*fubini_matrix.copy()

    
        regularized_fubini_matrix = np.add(np.power(exponentially_smoothed_fubini*exponentially_smoothed_fubini, 1/2), beta*np.identity(ansatz.num_parameters))

        regularized_fubini_matrix_previous = regularized_fubini_matrix

        # print(f'regularized_fubini_matrix---{regularized_fubini_matrix}')

        grad = ak*np.linalg.pinv(regularized_fubini_matrix).dot(gradSPSA)

        print(f'grad: {grad/ak}')
        print(f'parameter {internal_initial_point}')
        internal_initial_point = np.subtract(internal_initial_point, grad)
    
    return energy


def Customize_Variational_Quantum_EigenSolver(operator, optimize, interation, ansatz, eta, initial_point, shots, backend, callback):
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
    
    b = np.zeros(ansatz.num_parameters)

    for i in range(interation):
        print(f'{internal_initial_point} ---------')
        internal_anstaz = ansatz.bind_parameters({theta: internal_initial_point[k] for k, theta in enumerate(ansatz.parameters)})       
        # Measure the expectation of our hamiltonian
        internal_energy = 0
        for j in range(len(operator)):
            internal_energy += operator.coeffs.real[j] * Measurement(str(operator.primitive._pauli_list[j]), internal_anstaz, shots, backend)[0]
            # print( operator.coeffs.real[j])
            # print(operator.primitive._pauli_list[j])
            # print(internal_energy)
            
            # measure = Measurement('X', internal_anstaz, 10000, backend=Aer.get_backend('qasm_simulator'))[0]*10
            # print(f'///{measure}///')
        
        # print(f'{internal_initial_point} hmm')
        
        energy.append(internal_energy)
        print(internal_energy)


        if callback is not None:
            callback(internal_initial_point, internal_energy)

        # Update the parameter points
        if optimize == Customize_Quantum_Natural_Graident_Descent:
            super_circuit = Seperate_Circuit_Apart(ansatz)
            a = np.add(optimize(operator, internal_initial_point, eta, ansatz, shots, backend, interation, super_circuit), -0.9*b)
            # b = a.copy()
            internal_initial_point =  np.subtract(internal_initial_point,a)
            # print(internal_energy)

            continue

        a = optimize(operator, internal_initial_point, eta, ansatz, shots, backend, interation)
        # a = np.add(optimize(operator, internal_initial_point, eta, ansatz, shots, backend, interation), -0.9*b)
        # b = a.copy()
        internal_initial_point =  np.subtract(internal_initial_point,a)
        print(f'a: {a}')

        # print(internal_energy)

    if callback is None:
        return energy




def Ising_hamiltonian(num_qubits, h):
    hamiltonian = 0

    if num_qubits == 1:
        hamiltonian = Z + h*X
        return hamiltonian.reduce()

    if num_qubits == 2:
        hamiltonian = Z^Z
        hamiltonian += h*X^I
        hamiltonian += h*I^X
        return hamiltonian.reduce()


    
    hamiltonian += Z^(I^(num_qubits-2))^Z
    for i in range(num_qubits-1):
        hamiltonian += (I^(i))^(Z^Z)^(I^(num_qubits-2-i))
    
    for i in range(num_qubits):
        hamiltonian += h*((I^(i))^(X)^(I^(num_qubits-1-i)))

    return hamiltonian.reduce()



