
#pip install biopython

from Bio import pairwise2
from Bio.PDB import PDBParser, Superimposer
import numpy as np

def parse_pdb(filename):
    parser = PDBParser()
    structure = parser.get_structure(filename, filename)
    model = structure[0]  # Assuming there is only one model in the structure
    return model

def calculate_rmsd(model1, model2):
    atoms1 = [atom.get_coord() for atom in model1.get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom.get_coord() for atom in model2.get_atoms() if atom.get_id() == 'CA']
    
    if len(atoms1) != len(atoms2):
        raise ValueError("Number of CA atoms in the two models are different")
    
    return np.sqrt(np.mean(np.sum((np.array(atoms1) - np.array(atoms2))**2, axis=1)))

def calculate_tmscore(model1, model2):
    super_imposer = Superimposer()
    atoms1 = [atom for atom in model1.get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom for atom in model2.get_atoms() if atom.get_id() == 'CA']
    
    if len(atoms1) != len(atoms2):
        raise ValueError("Number of CA atoms in the two models are different")
    
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(model2.get_atoms())

    return super_imposer.rms, super_imposer.tm_score

if __name__ == "__main__":
    pdb_file1 = "C:\Users\GAFFNEM2\Documents\gitrepos\tmscoring\tmscoring\tests\pdb1.pdb"
    pdb_file2 = "C:\Users\GAFFNEM2\Documents\gitrepos\tmscoring\tmscoring\tests\pdb2.pdb"
    
    # Parse PDB files
    structure1 = parse_pdb(pdb_file1)
    structure2 = parse_pdb(pdb_file2)
    
    # Calculate RMSD
    rmsd = calculate_rmsd(structure1, structure2)
    
    # Calculate TM-score
    rms, tmscore = calculate_tmscore(structure1, structure2)
    
    print(f"RMSD between the structures: {rmsd:.2f} Å")
    print(f"TM-score between the structures: {tmscore:.4f}")
