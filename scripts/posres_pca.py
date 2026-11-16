"""Crea posre_prot.itp, posre_lig.itp y sys_posres.top a partir de sys.top (ParmEd).
Restringe los átomos pesados (masa > 2) de la proteína (primer moleculetype) y de PCA,
1000 kJ/mol/nm2, bajo #ifdef POSRES (igual que en FRa_FA). También cuenta los puentes S–S."""
import re

L = open('sys.top').read().split('\n')
mols = []  # (nombre, línea de inicio)
sec = None
for i, l in enumerate(L):
    s = l.split(';')[0].strip()
    m = re.match(r'\[\s*(\w+)\s*\]', s)
    if m:
        sec = m.group(1).lower()
        if sec == 'moleculetype':
            j = i + 1
            while not L[j].split(';')[0].strip():
                j += 1
            mols.append((L[j].split()[0], i))
        continue

def block(k):
    a = mols[k][1]
    b = mols[k + 1][1] if k + 1 < len(mols) else next(i for i, l in enumerate(L) if re.match(r'\s*\[\s*system\s*\]', l))
    return a, b

def heavy_and_ss(k):
    a, b = block(k)
    sec, heavy, names, ss = None, [], {}, 0
    for l in L[a:b]:
        s = l.split(';')[0].strip()
        if not s or s.startswith('#'):
            continue
        m = re.match(r'\[\s*(\w+)\s*\]', s)
        if m:
            sec = m.group(1).lower(); continue
        c = s.split()
        if sec == 'atoms':
            names[int(c[0])] = c[4]
            if float(c[7]) > 2.0:
                heavy.append(int(c[0]))
        elif sec == 'bonds' and names.get(int(c[0])) == 'SG' and names.get(int(c[1])) == 'SG':
            ss += 1
    return heavy, ss

kp = 0
kl = next(k for k, (n, _) in enumerate(mols) if n.upper() == 'PCA')
hp, ss = heavy_and_ss(kp)
hl, _ = heavy_and_ss(kl)

def itp(fn, idx):
    with open(fn, 'w') as f:
        f.write('[ position_restraints ]\n; atom  type  fx  fy  fz\n')
        for i in idx:
            f.write(f'{i:>6}  1  1000  1000  1000\n')

itp('posre_prot.itp', hp)
itp('posre_lig.itp', hl)

ins = {block(kp)[1]: 'posre_prot.itp', block(kl)[1]: 'posre_lig.itp'}
out = []
for i, l in enumerate(L):
    if i in ins:
        out += ['#ifdef POSRES', f'#include "{ins[i]}"', '#endif', '']
    out.append(l)
open('sys_posres.top', 'w').write('\n'.join(out))
print(f'moleculetypes: {[n for n, _ in mols]}')
print(f'proteína: {len(hp)} átomos pesados restringidos | PCA: {len(hl)} | puentes S-S: {ss}')
