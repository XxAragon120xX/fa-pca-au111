import numpy as np
from Bio.PDB import PDBParser

lig = PDBParser(QUIET=True).get_structure("lig", "ligando_cristal.pdb")
xyz = np.array([a.coord for a in lig.get_atoms()])
print("Átomos del ligando:", len(xyz))
print("Centro de la caja (Å):", xyz.mean(axis=0).round(2))
print("Extensión (Å):", (xyz.max(axis=0) - xyz.min(axis=0)).round(1))