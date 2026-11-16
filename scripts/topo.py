"""Arma la topología del sistema Au(111) + ligando + agua a partir de base.top (ParmEd, GAFF2/TIP3P/iones JC).
Uso: python topo.py base.top LIG n_Au n_WAT salida.top
Agrega el tipo AU (INTERFACE-FF, Heinz et al. 2008: sigma 0.2629 nm, epsilon 22.145 kJ/mol, q = 0),
restricciones de posición del ligando (#ifdef POSRES_LIG) y del Au (#ifdef POSRES_AU)."""
import re, sys

base, lig, nau, nw, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
L = open(base).read().split('\n')

def hdr(l):
    m = re.match(r'\s*\[\s*(\w+)\s*\]', l.split(';')[0]); return m.group(1).lower() if m else None

# localizar secciones
idx = [(i, hdr(l)) for i, l in enumerate(L) if hdr(l)]
first_mol = next(i for i, h in idx if h == 'moleculetype')
sys_i = next(i for i, h in idx if h == 'system')
head = L[:first_mol]
# añadir el tipo AU al final de [ atomtypes ]
at = next(i for i, h in idx if h == 'atomtypes')
end_at = next(i for i, h in idx if i > at)
j = end_at
while j > at and not L[j - 1].strip():
    j -= 1
head = L[:j] + ['AU            79  196.966570  0.00000000  A        0.2629      22.145   ; Au, INTERFACE-FF (Heinz 2008)'] + L[j:first_mol]

# bloques de moleculetype
mols, starts = {}, [i for i, h in idx if h == 'moleculetype'] + [sys_i]
for a, b in zip(starts[:-1], starts[1:]):
    k = a + 1
    while not L[k].split(';')[0].strip():
        k += 1
    mols[L[k].split()[0]] = L[a:b]
assert lig in mols, f'no encuentro el moleculetype {lig}: {list(mols)}'
wat = next(n for n in mols if n.upper() in ('WAT', 'HOH', 'SOL'))

# restricciones del ligando: átomos pesados (masa > 2)
heavy, sec = [], None
for l in mols[lig]:
    s = l.split(';')[0].strip()
    if hdr(l):
        sec = hdr(l); continue
    if sec == 'atoms' and s:
        c = s.split()
        if float(c[7]) > 2.0:
            heavy.append(int(c[0]))
posre = ['', '#ifdef POSRES_LIG', '[ position_restraints ]'] + [f'{i:6d}  1  1000  1000  1000' for i in heavy] + ['#endif', '']
au = ['[ moleculetype ]', '; Name  nrexcl', 'AU      1', '', '[ atoms ]',
      '    1   AU   1   AU   AU   1   0.000000   196.966570', '',
      '#ifdef POSRES_AU', '[ position_restraints ]', '    1  1  1000  1000  1000', '#endif', '']

body = head[:]
body += mols[lig] + posre
for n in mols:
    if n != lig:
        body += mols[n]
body += au
body += ['[ system ]', f'Au(111) + {lig} en agua', '', '[ molecules ]', '; Compound       #mols',
         f'AU      {nau}', f'{lig}      1', f'{wat}      {nw}', '']
open(out, 'w').write('\n'.join(body))
print(f'{out}: moleculetypes {list(mols) + ["AU"]}, {len(heavy)} átomos pesados del ligando restringibles, agua = {wat}')
