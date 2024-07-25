##Debuging trial 

##Import Statments
from __future__ import division
import math
import numpy as np
import iminuit
from Bio import PDB
from Bio.Align import PairwiseAligner
from tmscoring import get_tm, get_rmsd

##constant-data type
DTYPE = np.float64

 #Alignment Class initalisation
class Aligning(object):
    def __init__(self, pdb_1, pdb_2, mode='index', chain_1='A', chain_2='A', d0s=5.):
        self.pdb1 = pdb_1
        self.pdb2 = pdb_2
        self.chain_1 = chain_1
        self.chain_2 = chain_2

        if mode == 'align':
            self._load_data_alignment(chain_1, chain_2)
        elif mode == 'index':
            self._load_data_index(chain_1, chain_2)
        else:
            raise ValueError('Unrecognised mode {}'.format(mode))

        d0 = 1.24 * (self.N - 15) ** (1.0 / 3.0) - 1.8
        self.d02 = d0 ** 2
        self.d0s2 = d0s ** 2

        self._values = dict(dx=0, dy=0, dz=0, theta=0, phi=0, psi=0)

#including in tmscore but causing errors at the moment: intented to set up parameters 

    def __call__(self, theta, phi, psi, dx, dy, dz):
        #raise NotImplementedError('This method should be overriden by subclasses')
        return self.tmscore(theta, phi, psi, dx, dy, dz)

    def get_default_values(self):
        out = dict(dx=0, dy=0, dz=0, theta=0, phi=0, psi=0)
        dx, dy, dz, _ = np.mean(self.coord1 - self.coord2, axis=1)
        out['dx'] = dx
        out['dy'] = dy
        out['dz'] = dz

        vec1 = self.coord1[:-1, 1] - self.coord1[:-1, -1]
        vec2 = self.coord2[:-1, 1] - self.coord2[:-1, -1]
        vec1 /= np.linalg.norm(vec1)
        vec2 /= np.linalg.norm(vec2)

        v = np.cross(vec1, vec2)
        s = np.linalg.norm(v) + np.finfo(DTYPE).eps

        c = vec1.dot(vec2)
        vx = np.array([[0, -v[2], v[1]],
                       [v[2], 0, -v[0]],
                       [-v[1], v[0], 0]], dtype=DTYPE)
        rotation_matrix = np.eye(3) + vx + vx.dot(vx) * (1 - c) / (s * s)

        out['theta'] = math.atan2(rotation_matrix[2, 1], rotation_matrix[2, 2])
        out['phi'] = math.atan2(-rotation_matrix[2, 0],
                                math.hypot(rotation_matrix[2, 1],
                                           rotation_matrix[2, 2]))
        out['psi'] = math.atan2(rotation_matrix[1, 0], rotation_matrix[0, 0])

        return out

    @staticmethod
    def get_matrix(theta, phi, psi, dx, dy, dz,
                   matrix=np.zeros((4, 4), dtype=DTYPE),
                   angles=np.zeros(3, dtype=DTYPE)):
        angles[0] = theta
        angles[1] = phi
        angles[2] = psi

        cx, cy, cz = np.cos(angles)
        sx, sy, sz = np.sin(angles)

        rotation = matrix[:3, :3]
        rotation.flat = (cx * cz - sx * cy * sz,
                         cx * sz + sx * cy * cz, sx * sy,
                         -sx * cz - cx * cy * sz,
                         -sx * sz + cx * cy * cz, cx * sy,
                         sy * sz,
                         -sy * cz, cy)

        matrix[:3, 3] = dx, dy, dz
        matrix[3, 3] = 1.
        return matrix

    def optimise(self, restart=True):
        if restart:
            default = self.get_default_values()
        else:
            default = self.get_current_values()

# Create a dictionary of parameter names and their initial values
        parameters = {
            'theta': default['theta'],
            'phi': default['phi'],
            'psi': default['psi'],
            'dx': default['dx'],
            'dy': default['dy'],
            'dz': default['dz']
        }


        m = iminuit.Minuit(self, print_level=0, pedantic=False, errordef=self.errordef(),
                           **default)

       
       
        m.migrad()

        _values = m.values
        self._values = _values
        return _values, self.tmscore(**_values), self.rmsd(**_values)

    def get_current_values(self):
        return self._values



##Defining functions for calculations

    def _tm(self, theta, phi, psi, dx, dy, dz):
        matrix = self.get_matrix(theta, phi, psi, dx, dy, dz)
        coord = matrix.dot(self.coord2)
        dist = coord - self.coord1

        d_i2 = (dist * dist).sum(axis=0)
        tm = -(1 / (1 + (d_i2 / self.d02)))

        return tm

    def _s(self, theta, phi, psi, dx, dy, dz):
        matrix = self.get_matrix(theta, phi, psi, dx, dy, dz)
        coord = matrix.dot(self.coord2)
        dist = coord - self.coord1

        d_i2 = (dist * dist).sum(axis=0)
        tm = -(1 / (1 + (d_i2 / self.d0s2)))

        return tm

    def _rmsd(self, theta, phi, psi, dx, dy, dz):
        matrix = self.get_matrix(theta, phi, psi, dx, dy, dz)
        coord = matrix.dot(self.coord2)
        dist = coord - self.coord1
        return (dist * dist)

    def tmscore(self, theta, phi, psi, dx, dy, dz):
        return -np.mean(self._tm(theta, phi, psi, dx, dy, dz))

    def sscore(self, theta, phi, psi, dx, dy, dz):
        return -np.mean(self._s(theta, phi, psi, dx, dy, dz))

    def rmsd(self, theta, phi, psi, dx, dy, dz):
        return np.sqrt(np.mean(self._rmsd(theta, phi, psi, dx, dy, dz)))

    def write(self, outputfile='out.pdb', appended=False):
        matrix = self.get_matrix(**self.get_current_values())

        out = open(outputfile, 'w')
        atomid = 1
        if appended:
            for line in open(self.pdb1):
                if not line.startswith('ATOM') or (line[21] != self.chain_1 and line[21] != ' '):
                    continue
                out.write(line[:7])
                out.write('{: >4}'.format(atomid))
                atomid += 1
                out.write(line[11:21])
                out.write('A')
                out.write(line[22:])

        for line in open(self.pdb2):
            if not line.startswith('ATOM') or (line[21] != self.chain_2 and line[21] != ' '):
                continue

            x = float(line[32:38])
            y = float(line[39:46])
            z = float(line[48:54])
            vec = np.array([x, y, z, 1])
            x, y, z, _ = matrix.dot(vec)

            out.write(line[:7])
            out.write('{: >4}'.format(atomid))
            atomid += 1
            out.write(line[11:21])
            out.write('B')
            out.write(line[22:30])
            out.write('{:>8.3f}{:>8.3f}{:>8.3f}'.format(x, y, z))
            out.write(line[54:])

        out.close()

    def _load_data_alignment(self, chain1, chain2):
        parser = PDB.PDBParser(QUIET=True)
        ppb = PDB.PPBuilder()

        structure1 = parser.get_structure(chain1, self.pdb1)
        structure2 = parser.get_structure(chain2, self.pdb2)

        seq1 = str(ppb.build_peptides(structure1)[0].get_sequence())
        seq2 = str(ppb.build_peptides(structure2)[0].get_sequence())

        aligner = PairwiseAligner()
        alignments = aligner.align(seq1, seq2)

        indexes = set(i for i, (s1, s2) in enumerate(zip(alignments[0].seqA, alignments[0].seqB))
                      if s1 != '-' and s2 != '-')
        coord1 = np.hstack([np.concatenate((r['CA'].get_coord(), (1,)))[:, None]
                            for i, r in enumerate(structure1.get_residues())
                            if i in indexes and 'CA' in r]).astype(DTYPE, copy=False)
        coord2 = np.hstack([np.concatenate((r['CA'].get_coord(), (1,)))[:, None]
                            for i, r in enumerate(structure2.get_residues())
                            if i in indexes and 'CA' in r]).astype(DTYPE, copy=False)

        self.coord1 = coord1
        self.coord2 = coord2
        self.N = len(seq1)

    def _load_data_index(self, chain1, chain2):
        parser = PDB.PDB

##testing

pdb1 = "C:\\Users\\GAFFNEM2\\Downloads\\7rvk (2).pdb"
pdb2 = "C:\\Users\\GAFFNEM2\\Downloads\\7rvk (2).pdb"

tm_score = get_tm(pdb1, pdb2)
rmsd_score = get_rmsd(pdb1, pdb2)

print(f"TM score: {tm_score}")
print(f"RMSD: {rmsd_score}")
