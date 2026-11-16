from Bio.PDB import PDBList, PDBParser, PDBIO, Select
import numpy as np

PDB_ID = "4LRH"   # verifica con la salida que sea FRα + ácido fólico

PDBList().retrieve_pdb_file(PDB_ID, pdir=".", file_format="pdb")
s = PDBParser(QUIET=True).get_structure(PDB_ID, f"pdb{PDB_ID.lower()}.ent")

# 1) Qué hetero-moléculas hay
print("Heteroátomos por cadena:")
for ch in s[0]:
    het = {}
    for r in ch:
        if r.id[0].startswith("H") and r.resname != "HOH":
            het[r.resname] = het.get(r.resname, 0) + 1
    print(f"  cadena {ch.id}: {het}")

LIG = "FOL"       # nombre del ligando en el PDB (ajústalo si la salida dice otro)
CHAIN = "A"

class Prot(Select):
    def accept_chain(self, c): return c.id == CHAIN
    def accept_residue(self, r): return r.id[0] == " "      # solo aminoácidos
class Lig(Select):
    def accept_chain(self, c): return c.id == CHAIN
    def accept_residue(self, r): return r.resname == LIG

io = PDBIO(); io.set_structure(s)
io.save("proteina.pdb", Prot())
io.save("ligando_cristal.pdb", Lig())

xyz = np.array([a.coord for a in s[0][CHAIN] if a.get_parent().resname == LIG])
print("Átomos del ligando:", len(xyz))
print("Centro de la caja (Å):", xyz.mean(axis=0).round(2))