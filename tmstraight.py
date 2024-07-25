from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

ref_structure = "C:/Users/GAFFNEM2/Downloads/1k43.pdb"
model_structure = "C:/Users/GAFFNEM2/Downloads/1k43.pdb"

def calculate_tm_score(ref_structure, model_structure):
    # Load structures
    parser = PDBParser()
    ref_struct = parser.get_structure('reference', ref_structure)
    model_struct = parser.get_structure('model', model_structure)
    
    # Extract atoms for superimposition
    ref_atoms = []
    model_atoms = []
    
    for ref_model, model_model in zip(ref_struct.get_models(), model_struct.get_models()):
        for ref_chain, model_chain in zip(ref_model.get_chains(), model_model.get_chains()):
            ref_atoms.extend(ref_chain.get_atoms())
            model_atoms.extend(model_chain.get_atoms())
    
   # Convert atom coordinates to numpy arrays
    ref_coords = np.array([atom.coord for atom in ref_atoms])
    model_coords = np.array([atom.coord for atom in model_atoms])

    # Superimpose structures
    sup = SVDSuperimposer()
    sup.set(ref_coords, model_coords)
    #sup.apply_transformation(model_coords)
    #sup.apply(model_coords)

    #Apply transformation to model coordinates
    transformed_model_coords = sup.apply(model_coords)
    
    # Calculate TM-score
    N = len(ref_atoms)
    d = np.sqrt(sum(sup.rms ** 2 for atom in model_atoms) / N)
    TM_score = 1 / (1 + d ** 2)
    
    return TM_score

# Example usage
#ref_structure = "C:/Users/GAFFNEM2/Downloads/1k43.pdb"
#model_structure = "C:/Users/GAFFNEM2/Downloads/1k43.pdb"
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Calculate TM-score between two protein structures.")
    parser.add_argument("--ref_structure", required=True, help="Path to reference structure (PDB file).")
    parser.add_argument("--model_structure", required=True, help="Path to model structure (PDB file).")
    args = parser.parse_args()

    ref_structure = args.ref_structure
    model_structure = args.model_structure
    
    tm_score = calculate_tm_score(ref_structure, model_structure)
    print(f"TM-score between the structures: {tm_score:.4f}")

tm_score = calculate_tm_score(ref_structure, model_structure)
print(f"TM-score between the structures: {tm_score:.4f}")
