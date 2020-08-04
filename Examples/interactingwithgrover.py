from qiskit.aqua.algorithms import Grover
from qiskit.aqua import QuantumInstance
from qiskit import Aer
from qiskit import BasicAer
from qiskit import execute
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
from qiskit.visualization import plot_histogram
from pytket.circuit import Circuit
from pytket.qiskit import tk_to_qiskit
import ALALI

#In this example, we will be setting up Grover's search on top of ALALI for Guanine. 
# The results of this algorithm could be used to determine the differences in polarity, diffusional response, electron hole
# behavior, and more. Note: No particular use case is given here. Rather, this tutorial is to show you how simple it is to connect
# ALALI to IBM Q libraries. Please read the wiki to see how to best understand results for real use cases.


#First, we will go ahead and get a graph for glycine.

def getAtomicGraph(moleculeInfo):
    glycineRep = ALALI.molecule(moleculeInfo, 'smi')
    ALALI.print_circuit(glycineRep.circuit)
    print(glycineRep.graph.nodes(data=True))
    return glycineRep.graph 

#now, we will set up our oracle function. This sets up the qubits so that they have different qubit frequencies and wavelengths
#that can be used to represent the different characteristic that Grover will search for. In this case, 
#we are going to set up a function that sets up every qubit's gate frequency of the number of valance electrons in each element
#to help set up grover to help qubits represent electronic weight differences between atoms.
def construct_oracle(graph):
    """Generate the initial, unrouted circuit with the Rz gates based on electrostatics, given a networkx graph

    Parameters
    ----------
    mol_graph : networkx graph
        Networkx graph with nodes:atoms and edges:bonds

    Returns
    -------
    pytket Circuit
        simple pytket circuit with a qubit for each atom

    Notes
    -----
    The qubits are labelled with integers that match the integers found in
    the inputted networkx graph node labels"""
    # - generate list of edges along with list of inverted edges
    # - combine the inverted edges with non-inverted ones

    edges = list(graph.edges)
    inv_edges = [(edge[1], edge[0]) for edge in edges]
    combinededges = [val for pair in zip(edges, inv_edges) for val in pair]

    # Create circuit and unpack each edge, applying it as a our new circuit
    constraint_test_circuit = Circuit(graph.number_of_nodes())
    for qubit in range(graph.number_of_nodes()):
        constraint_test_circuit.H(qubit)

    # Go through every node and add the Rz gate that tries to represent electronic differences between the different atoms
    for index in range(0, graph.number_of_nodes()):
        if (graph.nodes[index]['symbol'] == 'O'):
            constraint_test_circuit.Rz(6*(3.14)/2, index)
        if (graph.nodes[index]['symbol'] == 'H'):
            constraint_test_circuit.Rz(1*(3.14)/2, index)
        if (graph.nodes[index]['symbol'] == 'N'):
            constraint_test_circuit.Rz(5*(3.14)/2, index)
        if (graph.nodes[index]['symbol'] == 'C'):
            constraint_test_circuit.Rz(4*(3.14)/2, index)

    return constraint_test_circuit

#Now, we will actually set up our Grover diffusion algorithm. We will be following the same steps as laid out by IBM Qiskit's tutorials. (https://qiskit.org/textbook/ch-algorithms/grover.html#3.-Example:-3-Qubits-)
#The main difference between their guide and ours is that we expanded upon their 3 qubit setup to create cells of such activity across all qubits, as 9 (number of nodes here) is a multiple of 3. 
#But, it would make sense to use the CX and 2 gate setup as shown for any multiple of two. 

def runGrovers(oracle):
    """Apply inversion about the average step of Grover's algorithm."""

    qubits = oracle.qubits
    nqubits = len(qubits)

    for q in range(nqubits):
        oracle.H(q)
        oracle.X(q)
    
      # Do controlled-Z
    for q in range(nqubits-2):
        oracle.H(q+2)
        oracle.CCX(q,q+1,q+2)
        oracle.H(q+2)

    for q in range(nqubits):
        oracle.X(q)
        oracle.H(q)
    ALALI.print_circuit(oracle)
    return oracle


glycineSmiles = 'C(C(=O)O)N'
oracle = construct_oracle(getAtomicGraph(glycineSmiles))
#need to convert the pytket created circuit into one that IBM can understand
grover_circuit = tk_to_qiskit(runGrovers(oracle))
#need to measure all results
grover_circuit.measure_all()
#Now, we can setup an actual device to try out our search algorithm, yay!
backend = Aer.get_backend('qasm_simulator')
shots = 1024
results = execute(grover_circuit, backend=backend, shots=shots).result()
answer = results.get_counts()
plot_histogram(answer)

#Okay, now you should see a pretty graph showing you your different qubit states and their respective probabilities. 
#As you can see, this algorithm does show the species within the chemical system that have the most polarity. However, this could 
#be further refined to best actually fit the respecive noises seen by each element involved.



