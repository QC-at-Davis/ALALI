from rdkit.Chem import AllChem as Chem # AllChem has most of Chem, and also includes advanced algorithms
import networkx as nx
import sys
from pytket.circuit import Circuit

class molecule:
	"""Class for working with a molecule and converting to an RDKit/Networkx/Pytket representation.

	Parameters
	----------
	data : str
		molecule data, either a SMILES string or a .mol file path string.
	filetype : str
		type of data input, 'smi' for SMILES or 'mol' for .mol.
	addHs : boolean, optional
		set to False if Hydrogens should not be added to the molecule.
	threeD : boolean, optional
		set to False if 3D conformation does not need to be calculated.
	doall : boolean, optional
		set to False if RDKit/Networkx/Pytket representations should not be created initially.

	Attributes
	----------
	data : str
		molecule data, either a SMILES string or a .mol file path string.
	filetype : str
		type of data input, 'smi' for SMILES or 'mol' for .mol.
	addHs : boolean, optional
		set to False if Hydrogens should not be added to the molecule.
	threeD : boolean, optional
		set to False if 3D conformation does not need to be calculated.
	doall : boolean, optional
		set to False if RDKit/Networkx/Pytket representations should not be created initially.
	mol : RDKit molecule
		RDKit molecule directly created from the data.
	graph : networkx graph
		Networkx graph with nodes:atoms and edges:bonds. 
		Nodes have the attribute 'symbol' which gives the element name (str).
		Edges have the attribute 'type' which gives the number of bonds (1, 1.5, 2, etc).
	circuit : pytket circuit. 
		Qubits:atoms and CXs:bonds, where the qubits are labelled with integers 
		that match the networkx graph node labels and RDKit molecule atomIdx
	"""

	def __init__(self, data, filetype, addHs=True, threeD=True, doall=True):
		"""When a molecule is created, an RDKit molecule is automatically created from the input data
		and filetype and stored in self.mol"""
		# initial data
		self.data = data
		self.filetype = filetype
		# since initial data is required, we can automatically create an RDKit molecule
		self.inputMol()
		# optional parameters for whether to add hydrogens or a 3D conformation
		self.addHs = addHs
		self.threeD = threeD
		# initialize graph and circuit as None, so that we can ignore them if they haven't been created yet
		self.graph = None
		self.circuit = False # Use False as default, since that != operator doesn't work on pytket circuit object
		# automatically create graph and circuit if one chooses to
		if doall:
			self.doall()

	def inputMol(self):
		"""Function for inputting molecule data and creating an RDKit molecule."""
		# import from smiles string
		if self.filetype == "smi":
			try:
				self.mol = Chem.MolFromSmiles(self.data) # generate RDKit molecule
			# raise an error if RDKit molecule creation goes wrong
			except:
				raise Exception("Molecule could not be created. Please check that the data and filetype inputs are correct.")
				sys.exit(1)
		# import from .mol file
		elif filetype == "mol":
			try:
				self.mol = Chem.MolFromMolFile(self.data) # generate RDKit molecule
			# raise an error if RDKit molecule creation goes wrong
			except:
				raise Exception("Molecule could not be created. Please check that the data and filetype inputs are correct.")
				sys.exit(1)
		# if filetype doesn't fall into those categories, raise an error
		else:
			raise Exception("Accepted data inputs: SMILES string or MOL filepath\nAccepted filetype inputs: 'smi' or 'mol'\nTry again.")
			sys.exit(1)

	def doall(self):
		"""Adds molecule features, creates a networkx graph and pytket circuit by calling class functions""" 
		self.addToMol() # Add hydrogens and conformation data, if applicable
		self.graphMol()	# Generate a networkx graph
		self.circuitMol() # Generate a pytket circuit

	def addToMol(self):
		"""Add hydrogens to molecule and/or calculate 3D conformation of molecule.

		Notes
		-----
		Whether or not to add hydrogens or generate a conformation can be set with the `addHs` and 
		`threeD` arguments when the molecule instance is created. You can also edit these parameters
		with:

		>>> molecule=(example_data, example_filetype) #create molecule instance
		>>> molecule.addHs = True #if you want hydrogens
		>>> molecule.addHs = False #if you do not want hydrogens

		and the same goes for `threeD`"""
		if self.addHs:
			self.mol = Chem.AddHs(self.mol)
		if self.threeD:
			Chem.EmbedMolecule(self.mol, randomSeed=0xf00d)

	def graphMol(self):
		"""Create a networkx graph for molecule, with nodes:atoms and edges:bonds. 

		Notes
		-----
		Graph nodes are labelled with unique integers. The atom element is stored in the node attribute 
		"symbol" as a string.
		The bond type (single, double, etc) is stored in the edge attribute "type" as a double."""
		# create networkx graph instance
		self.graph = nx.Graph()
		# add atoms as nodes, labelled with the RDKit atom index and symbol attribute with the element type (string)
		for atom in self.mol.GetAtoms():
			self.graph.add_node(atom.GetIdx(), symbol = atom.GetSymbol())
		# add bonds as edges, with type attribute with the bond type (double)
		for bond in self.mol.GetBonds():
			self.graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), type = bond.GetBondTypeAsDouble())

	def circuitMol(self):
		"""a function for creating a pytket circuit for the molecule, with qubits:atoms and CXs:bonds.
		
		Notes
		----------
		Pytket circuit qubits are labelled with unique index integers that align to the 
		RDKit molecule and networkx graph.
		If there doesn't yet exist a graph to circuit, a graph will be created automatically."""
		# If graph doesn't exist yet, create one
		if self.graph == None:
			self.graphMol()
		# Generate a pytket circuit
	    # - generate list of edges along with list of inverted edges
	    # - combine the inverted edges with non-inverted ones
		edges = list(self.graph.edges)
		inv_edges = [(edge[1], edge[0]) for edge in edges]
		cxs = [val for pair in zip(edges, inv_edges) for val in pair]
	    # Create circuit and unpack each edge, applying it as a CX
		self.circuit = Circuit(self.graph.number_of_nodes())
		for edge in cxs:
			self.circuit.CX(*edge)

	def output(self, name):
		"""Outputs RDKit mol, networkx graph, and/or circuit commands as files.
		
		Parameters
		----------
		name : string
			name used for output files"""
		# save RDKit molecule + conformation data (if applicable) as .mol file
		if self.mol != None:
			molBlock = Chem.MolToMolBlock(self.mol)
			print(molBlock,file=open(name+".mol",'w+'))
		# save networkx graph as .graphml file
		if self.graph != None:
			nx.write_graphml(self.graph, name+".graphml")
		# save pytket circuit commands as .txt file
		if self.circuit:
			print(self.circuit.get_commands(),file=open(name+".txt", 'w+'))












