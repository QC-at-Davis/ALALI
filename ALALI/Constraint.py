import networkx as nx
from pytket.circuit import Circuit
from pytket.qiskit import tk_to_qiskit
from pytket.routing import Architecture, route
from pytket.passes import SequencePass, DecomposeMultiQubitsIBM, DecomposeSingleQubitsIBM
from pytket.predicates import CompilationUnit

def load_graphml(graphml_file):
    """Load graphml file and convert to Networkx graph
    
    Networkx graphs can be exported as .graphml files, and this function serves to re-import
    those networkx graphs. It does so while retaining all of the data.

    Parameters
    ----------
    graphml_file : str
        path to .graphml file containing a molecule's networkx graph

    Returns
    -------
    networkx Graph"""
    return nx.read_graphml(graphml_file)

def cx_gen(mol_graph):
    """Generate the paired CXs representing interactions between atoms, given networkx graph.

    Parameters
    ----------
    mol_graph : networkx graph
        Networkx graph with nodes:atoms and edges:bonds

    Returns
    -------
    list
        list of graph edges paired together"""
    # generate list of tuples containing the edges in the integer-labeled graph
    edges = list(mol_graph.edges)
    # create a list of INVERTED edges from the original edge list
    inv_edges = [(edge[1], edge[0]) for edge in edges]
    # - Pair together an edge with its inverted counterpart
    # - Store them in one large list
    return [val for pair in zip(edges, inv_edges) for val in pair]


def circuit_gen(graph):
    """Generate the initial, unrouted circuit with the CXs, given a networkx graph

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
    cxs = [val for pair in zip(edges, inv_edges) for val in pair]

    # Create circuit and unpack each edge, applying it as a CX
    constraint_test_circuit = Circuit(graph.number_of_nodes())
    for edge in cxs:
        constraint_test_circuit.CX(*edge)

    # return the mapping from atoms to integers and the generated constraint circuit
    return constraint_test_circuit

def circuit_to_ibm(tk_circuit):
    """Convert circuit to IBM Qiskit compatible verison, given pytket circuit

    Parameters
    ----------
    tk_circuit : pytket Circuit
        A pytket circuit

    Returns
    -------
    pytket Circuit
        An IBM Qiskit compatible pytket circuit"""
    return tk_to_qiskit(tk_circuit)

def print_circuit(tk_circuit):
    """Convert circuit to IBM Qiskit supported gates, and print that circuit

    Parameters
    ----------
    tk_circuit : pytket Circuit
        A pytket circuit"""
    print(tk_to_qiskit(tk_circuit))

def route_circuit(tk_circuit, architecture_map):
    """Route the circuit to a given architecture map

    Parameters
    ----------
    tk_circuit : pytket Circuit
        A pytket circuit
    architecture_map : list
        A list of qubit pairings, which are in the form of tuples

    Returns
    -------
    pytket Circuit
        A pytket circuit routed for the given architecture

    Notes
    -----
    IBMQ architectures are provided in the IBMLayouts module, accessible with ALALI.ibm"""
    architecture = Architecture(architecture_map)
    routed_circuit = route(tk_circuit, architecture)
    return routed_circuit

# If a subsequent print of the routed circuit fails...
# Invoke the decomposition function with varying "levels" of decomposition
def decompose_circuit(tk_circuit, decomposition_lvl):
    """Function to decompose circuit, with decomposition_lvl 1 decomposing BRIDGEs and 2 decomposing as much as possible (SWAPs).

    While normal pytket Circuits may contain bridge gates, these are not compatible with IBM Qiskit architectures. Therefore,
    one might have to decompose the circuit before converting it to an IBM Qiskit compatible version.

    Parameters
    ----------
    tk_circuit : pytket Circuit
        A pytket circuit
    decomposition_lvl : integer
        use 1 to decompose BRIDGEs, 2 to decompose as far as possible including SWAPs

    Returns
    -------
    pytket Circuit
        A pytket circuit decomposed
    """
    cu = CompilationUnit(tk_circuit)
    # 1 to decompose BRIDGE gates, 2 to decompose as far as possible (includes acceptable SWAPs!)
    if decomposition_lvl == 1:
        seqpass = SequencePass([DecomposeMultiQubitsIBM()])
    elif decomposition_lvl == 2: 
        seqpass = SequencePass([DecomposeMultiQubitsIBM(), DecomposeSingleQubitsIBM()])
    # apply decomposition and return the resulting circuit
    seqpass.apply(cu)
    return cu.circuit
