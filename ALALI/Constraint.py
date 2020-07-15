import networkx as nx
from pytket.circuit import Circuit
from pytket.qiskit import tk_to_qiskit
from pytket.routing import Architecture, route
from pytket.passes import SequencePass, DecomposeMultiQubitsIBM, DecomposeSingleQubitsIBM
from pytket.predicates import CompilationUnit

# Load graphml file and convert to NetworkX Graph
def load_graphml(graphml_file):
    return nx.read_graphml(graphml_file)

# Generate the paired CXs representing interactions between atoms
def cx_gen(mol_graph):

    # generate list of tuples containing the edges in the integer-labeled graph
    edges = list(mol_graph.edges) #EDITED TO BE MOL_GRAPH INSTEAD OF MOLE_GRAPH_INTS, SINCE NETWORKX GRAPHS ARE DIFFERENT
    # create a list of INVERTED edges from the original edge list
    inv_edges = [(edge[1], edge[0]) for edge in edges]
    # - Pair together an edge with its inverted counterpart
    # - Store them in one large list
    return [val for pair in zip(edges, inv_edges) for val in pair]

# generate the initial, unrouted circuit with the CXs
def circuit_gen(graph):
    # - generate list of edges along with list of inverted edges
    # - combine the inverted edges with non-inverted ones
    edges = list(graph.edges)
    inv_edges = [(edge[1], edge[0]) for edge in edges]
    cxs = [val for pair in zip(edges, inv_edges) for val in pair]

    # Create circuit and unpack each edge, applying it as a CX
    constraint_test_circuit = Circuit(graph.number_of_nodes())
    for edge in cxs:
        constraint_test_circuit.CX(*edge)

    # return the mapping from atoms to integers and the generated constraint circuit
    return constraint_test_circuit

# convert circuit to IBM Qiskit compatible verison
def circuit_to_ibm(tk_circuit):
    return tk_to_qiskit(tk_circuit)

# print the circuit (implicitly convert to IBM Qiskit supported gates)
def print_circuit(tk_circuit):
    print(tk_to_qiskit(tk_circuit))

# Route the circuit to a given architecture map
def route_circuit(tk_circuit, architecture_map):
    architecture = Architecture(architecture_map)
    routed_circuit = route(tk_circuit, architecture)
    return routed_circuit

# If a subsequent print of the routed circuit fails...
# Invoke the decomposition function with varying "levels" of decomposition
def decompose_circuit(tk_circuit, decomposition_lvl):
    cu = CompilationUnit(tk_circuit)
    # 1 to decompose BRIDGE gates, 2 to decompose as far as possible (includes acceptable SWAPs!)
    if decomposition_lvl == 1:
        seqpass = SequencePass([DecomposeMultiQubitsIBM()])
    elif decomposition_lvl == 2: 
        seqpass = SequencePass([DecomposeMultiQubitsIBM(), DecomposeSingleQubitsIBM()])
    # apply decomposition and return the resulting circuit
    seqpass.apply(cu)
    return cu.circuit
