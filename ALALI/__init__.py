"""ALALI: Automated Loading assuming Atomic Level Interactions

ALALI helps map molecules to quantum computing circuits. With one line, a molecule can be converted to a pytket circuit, RDKit molecule, and networkx graph. These objects can be easily accessed as class attributes, and can be outputted as files with a simple output() function. Furthermore, we provide functions to route pytket circuits to IBMQ machines.

ALALI consists of one class, molecule, and several other functions for working with pytket circuits. They can all be accessed easily with:
>>>ALALI.$FUNCTION_OR_CLASS_NAME

Classes
-------
molecule : the main class used to input data and create an RDKit/networkx/pytket representation. View help(ALALI.molecule) for more information.

Functions
---------
load_graphml : load a graphml file as a networkx graph
cx_gen : generate paired CXs for bonded atoms in a molecule
circuit_gen : generate an unrouted pytket circuit from a networkx graph
circuit_to_ibm : convert circuit to IBM Qiskit compatible version, from pytket circuit
print_circuit : Convert circuit to IBM Qiskit compatible version, and print
route_circuit : given a circuit and architecture map, route a circuit
decompose_circuit : given a circuit and decomposition level, decompose a circuit.
View help(ALALI.$FUNCTION_NAME) for more information on the above functions.

Example
-------
import ALALI
proline = ALALI.molecule("C1CC(NC1)C(=O)O", "smi")
RDKit_molecule = proline.mol
networkx_graph = proline.graph
pytket_circuit = proline.circuit
routed_circuit = ALALI.route_circuit(proline.circuit, ALALI.ibm.ibmq_20_tokyo)
decomposed_circuit = ALALI.decompose_circuit(routed_circuit, 1)
ibm_circuit = ALALI.circuit_to_ibm(decomposed_circuit)
print(ibm_circuit)
"""
from ALALI.Chemistry import *
from ALALI.IBMLayouts import *
from ALALI.Constraint import *
__version__ = "1.0.0"



