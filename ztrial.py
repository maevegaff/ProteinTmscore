from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

def calculate_tm_score(ref_structure, model_structure):
    # Load structures
    parser = PDBParser(QUIET=True)
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
    try:
        sup.set(ref_coords, model_coords)
        transformed_model_coords = sup.apply(model_coords)
    except Exception as e:
        print(f"Superimposition failed: {e}")
        return None
    
    # Calculate TM-score
    N = len(ref_atoms)  # Number of atoms
    try:
        d = np.sqrt(np.sum(sup.rms ** 2) / N)
        TM_score = 1 / (1 + d ** 2)
    except Exception as e:
        print(f"TM-score calculation failed: {e}")
        return None
    
    return TM_score

# Example usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Calculate TM-score between two protein structures.")
    parser.add_argument("--ref_structure", required=True, help="Path to reference structure (PDB file).")
    parser.add_argument("--model_structure", required=True, help="Path to model structure (PDB file).")
    args = parser.parse_args()

    ref_structure = args.ref_structure
    model_structure = args.model_structure
    
    tm_score = calculate_tm_score(ref_structure, model_structure)
    if tm_score is not None:
        print(f"TM-score between the structures: {tm_score:.4f}")
    else:
        print("TM-score calculation failed.")
