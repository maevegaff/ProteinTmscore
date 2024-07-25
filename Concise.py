#!/usr/bin/env python3

# The MIT License
#
# Copyright (c) 2010-2016 Anders S. Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

#imports
from __future__ import division
import argparse
import math
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, Superimposer

def calculate_rmsd(ref_structure, model_structure, start_residue, end_residue):
    # Load structures
    parser = PDB.PDBParser(QUIET=True)
    ref_struct = parser.get_structure('reference', ref_structure)
    model_struct = parser.get_structure('model', model_structure)
    
    # Select what residues numbers you wish to align
    # and put them in a list
    atoms_to_be_aligned = range(start_residue, end_residue + 1)

    # Get the first model in the structures
    ref_model = ref_struct[0]
    model_model = model_struct[0]

    # Make a list of the atoms (in the structures) you wish to align.
    ref_atoms = []
    model_atoms = []

    # Iterate over all chains in the model to find all residues
    for ref_chain in ref_model:
        for ref_res in ref_chain:
            # Check if residue number is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                ref_atoms.append(ref_res['CA'])  # Use 'CA' atom for alignment

    # Do the same for the model structure
    for model_chain in model_model:
        for model_res in model_chain:
            if model_res.get_id()[1] in atoms_to_be_aligned:
                model_atoms.append(model_res['CA'])  # Use 'CA' atom for alignment

    # Initialize the superimposer
    super_imposer = PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, model_atoms)
    super_imposer.apply(model_model.get_atoms())

    # Calculate RMSD
    rmsd = super_imposer.rms

    # Save RMSD into an output file
    with open("Testrmsd.txt", 'w') as rmsd_out:
        rmsd_out.write(str(rmsd))

    print(f"RMSD: {rmsd:.4f}")

    return rmsd

###this is the problem!!!!
def calculate_tm_score(ref_structure, model_structure):
    # Implement your TM-score calculation here
    # This part is still missing in your provided code
    parser = PDBParser()
    ref_struct = parser.get_structure('reference', ref_structure)
    model_struct = parser.get_structure('model', model_structure)
    
    # Perform superimposition
    superimposer = Superimposer()
    # Perform your superimposition steps here, e.g., based on residue selection
    
    # Apply transformation to model structure
    superimposer.set_atoms(ref_struct.get_atoms(), model_struct.get_atoms())
    superimposer.apply(model_structure)
    
    # Calculate TM-score (replace with actual TM-score calculation method)
    tm_score = superimposer.tm

    #tm_score = 0.0  # Replace with actual TM-score calculation

    # Save TM-score into an output file
    with open("tmscore.txt", 'w') as tm_out:
        tm_out.write(f"TM-score: {tm_score}")

    print(f"TM-score: {tm_score}")

    return tm_score

def main():
    parser = argparse.ArgumentParser(
        description='Calculate RMSD and TM-score between two protein structures.')
    parser.add_argument('--start_residue', type=int, required=True,
                        help='Start residue number for alignment')
    parser.add_argument('--end_residue', type=int, required=True,
                        help='End residue number for alignment')
    parser.add_argument('--ref_structure', required=True,
                        help='Path to reference structure (PDB file)')
    parser.add_argument('--model_structure', required=True,
                        help='Path to model structure (PDB file)')
    args = parser.parse_args()

    # Calculate RMSD
    rmsd = calculate_rmsd(args.ref_structure, args.model_structure,
                          args.start_residue, args.end_residue)

    # Calculate TM-score
    tm_score = calculate_tm_score(args.ref_structure, args.model_structure)

if __name__ == "__main__":
    main()
