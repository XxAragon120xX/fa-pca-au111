"""Genera la geometría 3D de un ligando y la información de sus fragmentos.
Uso: python make_ligand.py NOMBRE SMILES salida.sdf salida.json
- RDKit ETKDG (semilla 42, 60 confórmeros) + MMFF94.
- Si el ligando tiene γ-carboxilato (FA), se elige el confórmero más extendido
  (máxima distancia pterín–γ-C), necesario para la configuración "de pie sobre la cola".
  Si no, el de menor energía MMFF."""
import sys, json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

name, smi, out_sdf, out_json = sys.argv[1:5]
m = Chem.AddHs(Chem.MolFromSmiles(smi))
cids = AllChem.EmbedMultipleConfs(m, numConfs=60, randomSeed=42)
res = AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=5000)

def match(smarts):
    return list(m.GetSubstructMatches(Chem.MolFromSmarts(smarts)))

pter = match('[#7;H2]~[#6]1~[#7]~[#6]2~[#7]~[#6]~[#6]~[#7]~[#6]2~[#6](~[#8])~[#7]1')
assert len(pter) == 1, 'pterín no encontrado'
pter = list(pter[0])
amino, o4, n3 = pter[0], pter[10], pter[11]
gamma = match('[CH2]-[CX3](=[OX1])-[OX1-]')          # γ-carboxilato de FA: C(=O)O- unido a CH2
alpha = match('[CH1](-[NX3])-[CX3](=[OX1])-[OX1-]')   # α-carboxilato de FA
coo6 = match('c-[CX3](=[OX1])-[OX1-]')               # carboxilato del PCA (unido al anillo)
paba = match('c1ccc(cc1)C(=O)N')                     # anillo del pABA
info = {'name': name, 'smiles': smi, 'charge': Chem.GetFormalCharge(m),
        'pterin': pter, 'amino_N': amino, 'O4': o4, 'N3': n3,
        'heavy': [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1],
        'gamma_COO': list(gamma[0][1:]) if gamma else [],
        'alpha_COO': list(alpha[0][2:]) if alpha else [],
        'ring_COO': list(coo6[0][1:]) if coo6 else [],
        'pABA_ring': list(paba[0][:6]) if paba else []}

E = [e for _, e in res]
if info['gamma_COO']:
    gC = info['gamma_COO'][0]
    def ext(cid):
        X = m.GetConformer(cid).GetPositions()
        return np.linalg.norm(X[pter].mean(0) - X[gC])
    best = max(cids, key=ext)
    print(f'{name}: confórmero más extendido, pterín–γC = {ext(best):.1f} Å')
else:
    best = int(np.argmin(E))
    print(f'{name}: confórmero de menor energía MMFF')
mol = Chem.Mol(m, confId=best)
c = Chem.Conformer(m.GetConformer(best)); mol.RemoveAllConformers(); mol.AddConformer(c, assignId=True)
mol.SetProp('_Name', name)
Chem.MolToMolFile(mol, out_sdf)
json.dump(info, open(out_json, 'w'), indent=1)
print(f'{name}: carga {info["charge"]}, {mol.GetNumAtoms()} átomos, pterín {len(pter)} átomos, '
      f'γ-COO {info["gamma_COO"]}, α-COO {info["alpha_COO"]}, COO anillo {info["ring_COO"]}')
