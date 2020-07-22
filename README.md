# ALALI - Automated Loading assuming Atomic Level Interactions

ALALI is an object-oriented solution to mapping a molecule to a quantum computer. It maps atoms to qubits and bonds to CX gates. All you need is a SMILES representation of your molecule, and a few lines of code!

## Features
* Starting from a SMILES string, access:
	* **Pytket circuit** (atoms:qubits, bonds:CXs)
	* **Networkx graph** (atoms:nodes, bonds:edges) + visualization
	* **RDKit molecule**
* Add hydrogens to the molecule  
* Work with 3D molecules using ETKDG conformation generation
* Access bond types (single, double, aromatic) and atom elements as graph attributes
* Easily output graph and molecule data to .graphml and .mol files
* Map to several IBMQ topologies with our built-in IBMQ topology library
* Deconstruct circuits to avoid BRIDGE and SWAP gates, if necessary

## Getting Started
### Prerequisites
RDKit is a prerequisite for ALALI. It can be installed by following [these instructions](https://www.rdkit.org/docs/Install.html).  
It can be easily installed with: 
```
conda install -c conda-forge rdkit
```
Pytket is also a prerequisite, but it will be automatically installed with the pip command shown below.
It is important to note that Pytket currently does not support Windows, only Linux/MacOS. 
Therefore, ALALI also only supports Linux/MacOS.
### Installation
Once RDKit is installed, ALALI can be installed with:

```
pip install ALALI
```

### Example
ALALI can be easily imported with:  

	import ALALI

From there, a molecule instance can be created:

```
example = ALALI.molecule(data, filetype, addHs=True, threeD=True, doall=True)
```

* `data` = a SMILES string or a .mol file path
* `filetype`: either `'smi'` for SMILES data, or `'mol'` for a .mol file
* `addHs`: for adding Hydrogens to a molecule
* `threeD`: for generating a 3D conformation using RDKit's ETKDG method
* `doall`: for automatically generating an RDKit molecule and networkx graph when you create the instance

RDKit molecule, networkx graph, and pytket circuit can be accessed as class attributes

```
networkx_graph = example.graph
rdkit_molecule = example.mol
pytket_circuit = example.circuit
```

Now, we can output a .mol file (RDKit molecule), .graphml file (Networkx graph), and .txt file (Pytket circuit commands)
```
example.output("example_molecule")
```

That's it! For more in-depth examples and documentation, see the [Github wiki](https://github.com/QC-at-Davis/ALALI/wiki)
or use the help() command to see what individual functions do. We've provided several convenient pytket functions to route
circuits to certain toplogies, including IBMQ topologies (again, this can be found in the documentation).

## Support
If there is something wrong, or a feature you would like to see, please feel free to contact Jack Goon at jbgoon@ucdavis.edu. Furthermore, feel free to submit issues on [Github](https://github.com/QC-at-Davis/biopython-modeling/issues). 