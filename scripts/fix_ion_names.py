"""genion escribe los iones con nombre de átomo 'Na'/'Cl'; la topología usa 'Na+'/'Cl-'. Uso: python fix_ion_names.py sys.gro"""
import sys
f = sys.argv[1]; G = open(f).read().split('\n'); n = int(G[1]); k = 0
for i in range(2, 2 + n):
    rn = G[i][5:10].strip()
    if rn in ('Na+', 'Cl-'):
        G[i] = G[i][:10] + f'{rn:>5}' + G[i][15:]; k += 1
open(f, 'w').write('\n'.join(G)); print(f'   {k} iones renombrados')
