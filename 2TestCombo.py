from __future__ import division
import subprocess
import sys
import os
import argparse
import numpy as np
from numpy.testing import assert_almost_equal, TestCase
from nose.plugins.skip import SkipTest
from shutil import which
from Bio import PDB
import tmscoring

def __main__():
    parser = argparse.ArgumentParser(
        description='Residues to be aligned')

    parser.add_argument(
                        '--start_residue', default=None,
                        help='start residue')
    parser.add_argument(
                        '--end_residue', default=None,
                        help='end residue')
    parser.add_argument(
                        '--ref_structure', default=None,
                        help='reference structure')
    parser.add_argument(
                        '--model', default=None,
                        help='model structure')
    parser.add_argument(
                        '--aligned_structure', default=None,
                        help='aligned structure')
    parser.add_argument(
                        '--rmsd', default=None,
                        help='rmsd')
    args = parser.parse_args()
   
    #Hard coding file path 
    pdb1="C:\\Users\\GAFFNEM2\\Downloads\\7rvk (2).pdb"
    pdb2="C:\\Users\\GAFFNEM2\\Downloads\\7rvk (2).pdb"

     # Print statements to debug
    print(f"pdb1 path: {pdb1}")
    print(f"pdb2 path: {pdb2}")


   #following existing logic
    sc = tmscoring.TMscoring(pdb1, pdb2)
    _, tm, rmsd = sc.optimise()

 # Select what residues numbers you wish to align
    # and put them in a list
    start_id = int(args.start_residue)
    end_id = int(args.end_residue)
    atoms_to_be_aligned = range(start_id, end_id + 1)

    # Start the parser
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)

# Get the structures
    ref_structure = pdb_parser.get_structure("reference", args.ref_structure)
    sample_structure = pdb_parser.get_structure("sample", args.model)

 # Use the first model in the pdb-files for alignment
    # Change the number 0 if you want to align to another structure
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]

    # Make a list of the atoms (in the structures) you wish to align.
    # In this case we use CA atoms whose index is in the specified range
    ref_atoms = []
    sample_atoms = []

    # Iterate of all chains in the model in order to find all residues
    for ref_chain in ref_model:
        # Iterate of all residues in each model in order to find proper atoms
        for ref_res in ref_chain:
            # Check if residue number ( .get_id() ) is in the list
            if ref_res.get_id()[1] in atoms_to_be_aligned:
                # Append CA atom to list
                ref_atoms.append(ref_res['CA'])

    # Do the same for the sample structure
    for sample_chain in sample_model:
        for sample_res in sample_chain:
            if sample_res.get_id()[1] in atoms_to_be_aligned:
                sample_atoms.append(sample_res['CA'])


class TestAligningBase(TestCase):
    def test_matrix(self):
        align_object = tmscoring.Aligning('pdb1.pdb', 'pdb2.pdb')
        np.random.seed(124)
        for _ in range(100):
            theta, phi, psi = 2 * np.pi * np.random.random(3)
            dx, dy, dz = 10 * np.random.random(3)

            matrix = align_object.get_matrix(theta, phi, psi, dx, dy, dz)
            rotation = matrix[:3, :3]
            assert_almost_equal(1, np.linalg.det(rotation), decimal=6)
            assert_almost_equal(1, np.linalg.det(matrix), decimal=6)

    def test_tm_valuex(self):
        align_object = tmscoring.Aligning('pdb1.pdb', 'pdb2.pdb')
        np.random.seed(124)
        for _ in range(100):
            theta, phi, psi = 2 * np.pi * np.random.random(3)
            dx, dy, dz = 10 * np.random.random(3)

            tm = align_object._tm(theta, phi, psi, dx, dy, dz)

            assert np.all(0 <= -tm / align_object.N)

    def test_load_data_alignment(self):
        align_object = tmscoring.Aligning('pdb1.pdb', 'pdb2.pdb', mode='align')
        assert align_object.coord1.shape[0] == 4
        assert align_object.coord2.shape[0] == 4
        assert align_object.coord1.shape == align_object.coord2.shape

    def test_load_data_index(self):
        align_object = tmscoring.Aligning('pdb1.pdb', 'pdb2.pdb', mode='index')
        assert align_object.coord1.shape[0] == 4
        assert align_object.coord2.shape[0] == 4
        assert align_object.coord1.shape == align_object.coord2.shape

def test_identity():
    sc = tmscoring.TMscoring('pdb1.pdb', 'pdb1.pdb')
    assert sc.tmscore(0, 0, 0, 0, 0, 0) == 1

    sc = tmscoring.RMSDscoring('pdb1.pdb', 'pdb1.pdb')
    assert sc.rmsd(0, 0, 0, 0, 0, 0) == 0.0


def test_tm_output():
    if which("TMscore") is None:
        raise SkipTest('TMscore is not installed in the system.')

    pdb1, pdb2 = 'pdb1.pdb', 'pdb2.pdb'
    sc = tmscoring.TMscoring(pdb1, pdb2)
    _, tm, rmsd = sc.optimise()

    p = subprocess.Popen('TMscore {} {} | grep TM-score | grep d0'.format(pdb1, pdb2), stdout=subprocess.PIPE, shell=True)
    ref_tm = float(p.communicate()[0].decode('utf-8').split('=')[1].split('(')[0])
    assert_almost_equal(ref_tm, tm, decimal=2)

    p = subprocess.Popen('TMscore {} {} | grep RMSD | grep common'.format(pdb1, pdb2),
                         stdout=subprocess.PIPE, shell=True)
    ref_rmsd = float(p.communicate()[0].decode('utf-8').split('=')[1])
    assert abs(ref_rmsd - rmsd) < 0.1

def test_repeated():
    pdb1, pdb2 = 'pdbrep_1.pdb', 'pdbrep_2.pdb'
    sc = tmscoring.TMscoring(pdb1, pdb2)
    _, tm, rmsd = sc.optimise()

    assert_almost_equal(tm, 0.27426501120343644)
    assert_almost_equal(rmsd, 15.940038528551929)

if __name__ == "__main__":
    __main__()
