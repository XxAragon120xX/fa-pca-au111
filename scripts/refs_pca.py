"""Coordenadas de referencia (marco del cristal 4LRH) de los 12 átomos del pterín,
en este orden de etiquetas cristalográficas (ligando FOL)."""
import numpy as np

LAB = ['NA2', 'C2', 'N1', 'C8A', 'N8', 'C7', 'C6', 'N5', 'C4A', 'C4', 'O4', 'N3']
# índices dentro del residuo PCA (orden del SDF del docking = orden de la topología)
PCA_PT = [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14]
PCA_COO = [7, 8, 9]              # C y los dos O del carboxilato
PCA_ELEM = 'NCNCNCCCOONCCON'     # elementos esperados de los 15 átomos pesados
# nombres de los mismos 12 átomos en la topología de FA (FOA), mismo orden que LAB
FOA_PT = ['N1', 'C1', 'N2', 'C2', 'N3', 'C3', 'C4', 'N6', 'C18', 'C19', 'O6', 'N7']

CRY = np.array([
    [0.248, 22.065, 9.241], [1.119, 22.308, 8.239], [0.639, 22.832, 7.090], [1.497, 23.087, 6.020],
    [1.058, 23.601, 4.854], [1.911, 23.826, 3.827], [3.291, 23.522, 3.974], [3.750, 23.015, 5.134],
    [2.901, 22.790, 6.159], [3.340, 22.233, 7.433], [4.507, 21.969, 7.599], [2.415, 22.022, 8.410]])
NAT = np.array([   # pose 6 del docking (tipo nativa)
    [0.651, 22.328, 9.365], [1.386, 22.523, 8.242], [0.819, 22.908, 7.149], [1.641, 23.105, 6.026],
    [0.997, 23.538, 4.923], [1.791, 23.740, 3.832], [3.174, 23.532, 3.813], [3.815, 23.092, 4.945],
    [3.018, 22.892, 6.022], [3.610, 22.428, 7.295], [4.794, 22.156, 7.463], [2.724, 22.302, 8.345]])
TOP = np.array([   # pose 1 del docking (carboxilato en el sitio del O4)
    [0.076, 25.368, 0.942], [0.745, 24.784, 1.968], [0.276, 24.838, 3.168], [1.012, 24.200, 4.181],
    [0.456, 24.264, 5.409], [1.161, 23.638, 6.396], [2.371, 22.965, 6.196], [2.923, 22.909, 4.939],
    [2.217, 23.535, 3.967], [2.719, 23.525, 2.577], [3.776, 23.021, 2.214], [1.905, 24.147, 1.653]])
O4_SITE = CRY[LAB.index('O4')]

def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))

def crystal_fit(u, cry_ca_map):
    """Devuelve (ca, keep, ref, rc) para superponer Cα de cada frame sobre el cristal (mismo número de residuo)."""
    ca = u.select_atoms('protein and name CA')
    keep = [i for i, a in enumerate(ca) if a.resid in cry_ca_map]
    ref = np.array([cry_ca_map[ca[i].resid] for i in keep])
    return ca, keep, ref, ref.mean(0)
