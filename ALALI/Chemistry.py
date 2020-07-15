from rdkit.Chem import AllChem as Chem # AllChem has most of Chem, and also includes advanced algorithms
import networkx as nx
import matplotlib.pyplot as plt
import sys

class molecule:
	# Initialize input data and variables, run the self.main() function
	def __init__(self, data, filetype, name="molecule", AddHs=True, ThreeD=True, doall=True):
		self.data = data
		self.filetype = filetype
		self.name = name
		self.AddHs = AddHs
		self.ThreeD = ThreeD
		self.mol = None
		self.graph = nx.Graph()
		if doall:
			self.doall()

	# Start with SMILES or MOL format, create RDKit molecule, add Hydrogens, add 3D force field using ETKDG
	def inputMol(self):
		if self.filetype == "smi":
			self.mol = Chem.MolFromSmiles(self.data)
		elif self.filetype == "mol":
			self.mol = Chem.MolFromMolFile(self.data)
		else:
			raise Exception("Accepted data inputs: SMILES string or MOL filepath\nAccepted filetype inputs: 'smi' or 'mol'\nTry again.")
			sys.exit(1)

	# Add hydrogens, make 3D
	def addToMol(self):
		if self.mol == None:
			raise Exception("Molecule was incorrectly created. Cannot add Hydrogens or 3D coordinates.")
			sys.exit(1)
		else:
			if self.AddHs:
				self.mol = Chem.AddHs(self.mol)
			if self.ThreeD:
				Chem.EmbedMolecule(self.mol, randomSeed=0xf00d)

	# Create graph of nodes=atoms and edges=bonds.
	# TODO: Add geometry data, charges
	def graphMol(self):
		if self.mol == None:
			raise Exception("Molecule was incorrectly created. Graph couldn't be created.")
			sys.exit(1)
		else:
			for atom in self.mol.GetAtoms():
				self.graph.add_node(atom.GetIdx(), symbol = atom.GetSymbol())
			for bond in self.mol.GetBonds():
				self.graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), type = bond.GetBondTypeAsDouble())

	# Output in the form of files (.mol file, networkx graph as .graphml) and/or print (.mol data, networkx graph visualized)
	def output(self,file=True, printing=True):
		molBlock = Chem.MolToMolBlock(self.mol)
		if file:
			print(molBlock,file=open(self.name+".mol",'w+'))
			nx.write_graphml(self.graph, self.name+".graphml")
		if printing:
			print(molBlock)
			nx.draw(self.graph, with_labels=True)
			plt.show()

	# input molecule, add Hs/3D coordinates, and create graph
	def doall(self): #does all of the functions, to be run in __init__()
		self.inputMol()
		self.addToMol()
		self.graphMol()
