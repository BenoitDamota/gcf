from rdkit.Chem.rdmolops import GetAdjacencyMatrix, SanitizeMol
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Kekulize, FastFindRings

import rdkit.Chem as Chem
import networkx as nx
import numpy as np

def get_bridge_bonds_matrix(adj_mat, mol):
    """
    Returning a boolean matrix of size (n_defined_atoms, n_defined_atoms) representing whether bonds of the
    molecular graph are bridges.
    """

    # Converting the molecular graph to a NetworkX object
    nx_mol_graph = nx.from_numpy_array(adj_mat)

    # Initialization of the output matrix of bridge bonds
    output_bridges_matrix = np.full((mol.GetNumAtoms(), mol.GetNumAtoms()), False)

    # Extracting the list of bridges in the molecular simple graph
    bridges_list = list(nx.bridges(nx_mol_graph))

    for bridge in bridges_list:
        output_bridges_matrix[bridge[0], bridge[1]] = True
        output_bridges_matrix[bridge[1], bridge[0]] = True

    return output_bridges_matrix

def update_mol_representation(mol_graph):
    """
    Updating internal RDKIT representation of the molecular graph
    """

    # Sanitizing mol if needed
    SanitizeMol(mol_graph)

    # Clear computed properties
    mol_graph.ClearComputedProps()

    # Kekulization of the molecular graph
    Kekulize(mol_graph)

    # Setting all atoms to non aromatics
    for i in range(mol_graph.GetNumAtoms()):
        mol_graph.GetAtomWithIdx(i).SetIsAromatic(False)

    # Setting all bonds to non aromatics
    for i in range(mol_graph.GetNumAtoms()):
        for j in range(mol_graph.GetNumAtoms()):
            bond = mol_graph.GetBondBetweenAtoms(i, j)
            if bond is not None:
                bond.SetIsAromatic(False)

    # Updating the property cache of atoms
    for i in range(mol_graph.GetNumAtoms()):
        mol_graph.GetAtomWithIdx(i).UpdatePropertyCache()

    # Updating RDKit representation
    mol_graph.UpdatePropertyCache()
    FastFindRings(mol_graph)

def compute_generic_scaffold_with_insat(smi):
    """
    Compute a generic cyclic feature (with insaturations) of a cyclic subgraph. 

    For that, it computes the Murcko scaffold on a subgraph (that includes a cycle), then all atoms 
    with a coordination number of 4 or less are converted to carbon atoms. Since hypervalent carbon 
    produce RDkit errors, hypervalent atoms are left unchanged. 
    """
    mol = Chem.MolFromSmiles(smi)
    smi_scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smi),includeChirality=False)
    mol_scaf = Chem.MolFromSmiles(smi_scaf)
    update_mol_representation(mol_scaf)

    for i in range(mol_scaf.GetNumAtoms()):
        if mol_scaf.GetAtomWithIdx(i).GetTotalValence() <= 4 :
            # Changing atomic number
            mol_scaf.GetAtomWithIdx(i).SetAtomicNum(6)
            # Setting formal charge to 0
            mol_scaf.GetAtomWithIdx(i).SetFormalCharge(0)
        # else keep the atom type (eg. P with valence = 5)

    update_mol_representation(mol_scaf)
    smi_scaf = Chem.MolToSmiles(mol_scaf)
    return smi_scaf


def compute_generic_cyclic_features_with_insat(smi):
    """
    Compute all generic cyclic features (with insaturations) of a molecules. 

    For that, it breaks all the vertices (bonds) that do not belong to a cycle are deleted (called 
    bridge in graph vocabulary). Then, compute_generic_scaffold_with_insat is called on all 
    remaining cyclic subgraphs.

    """
    # build molecular graph
    mol = Chem.MolFromSmiles(smi)

    #build RWMol object
    rwmol = RWMol(mol)
    update_mol_representation(rwmol)
    
    # compute adjacency matrix, bridges, then coordinates of bridges
    adjmat = GetAdjacencyMatrix(rwmol)
    mol_bridges = get_bridge_bonds_matrix(adjmat, rwmol)
    x, y = np.where(mol_bridges == True)
    
    # delete all bridges
    for i in range(len(x)):
        rwmol.RemoveBond(int(x[i]), int(y[i]))
        
    # update representation
    update_mol_representation(rwmol)
    
    # computes smiles of rings features and remove unique atoms
    features = [s for s in Chem.MolToSmiles(rwmol).split(".") if len(s) > 4]
    smi_features = [compute_generic_scaffold_with_insat(f) for f in features]
            
    return smi_features
