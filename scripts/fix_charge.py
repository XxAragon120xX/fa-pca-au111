"""Reparte el error de redondeo de las cargas de un mol2 para que sumen exactamente un entero.
Uso: python fix_charge.py archivo.mol2 carga_entera"""
import sys
f, Q = sys.argv[1], int(sys.argv[2])
L = open(f).read().split('\n')
i0 = L.index('@<TRIPOS>ATOM')
i1 = next(i for i in range(i0 + 1, len(L)) if L[i].startswith('@<TRIPOS>'))
q = [float(l.split()[8]) for l in L[i0 + 1:i1]]
d = (Q - sum(q)) / len(q)
new = [round(x + d, 6) for x in q]
new[-1] = round(new[-1] + (Q - sum(new)), 6)
for k, x in zip(range(i0 + 1, i1), new):
    c = L[k].split()
    L[k] = f'{int(c[0]):>7} {c[1]:<8}{float(c[2]):>10.4f}{float(c[3]):>10.4f}{float(c[4]):>10.4f} {c[5]:<6}{int(c[6]):>4} {c[7]:<8}{x:>10.6f}'
open(f, 'w').write('\n'.join(L))
print(f'{f}: carga {sum(q):+.4f} -> {sum(new):+.4f}')
