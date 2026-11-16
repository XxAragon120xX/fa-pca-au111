#!/usr/bin/env python3
"""make_scf.py — arma los SCF de producción A (complejo), B (slab congelado) y C (molécula congelada)
a partir de la relajación de cada configuración del FA (FA_fix/<conf>/relax.out).

Uso (en ~/Au/IN/FA_fix):
    python3 make_scf.py                      # todas las carpetas FA_* con relax.out
    python3 make_scf.py FA_V-Z-Pt FA_V-X-B   # solo esas
Crea ~/Au/SCF_final/FAfix_<conf>/{A.in,B.in,C.in,run_ABC.sh} y copia post_electronic.sh.
Protocolo idéntico al del PCA: 60/600 Ry, 4x4x1, MV 0.01 Ry, D3, conv_thr 1e-8.
"""
import re, sys, os, glob, shutil
HOME = os.path.expanduser('~')
DEST = os.environ.get('SCF_ROOT', os.path.join(HOME, 'Au', 'SCF_final'))
POST = next((p for p in (os.path.join(HOME, 'Au', 'post_electronic.sh'), os.path.join(HOME, 'Au', 'SCF_final', 'post_electronic.sh'),
                          os.path.join(HOME, 'post_electronic.sh')) if os.path.exists(p)), None)
PSEUDO = {'Au': ('196.97', 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'), 'C': ('12.011', 'C.pbe-n-rrkjus_psl.1.0.0.UPF'),
          'O': ('15.999', 'O.pbe-n-rrkjus_psl.1.0.0.UPF'), 'N': ('14.007', 'N.pbe-n-rrkjus_psl.1.0.0.UPF'),
          'H': ('1.008', 'H.pbe-rrkjus_psl.1.0.0.UPF')}

def cell_from_in(fn):
    L = open(fn).read().splitlines()
    for i, l in enumerate(L):
        if l.strip().upper().startswith('CELL_PARAMETERS'):
            return [L[i + k].split()[:3] for k in (1, 2, 3)]
    sys.exit(f'sin CELL_PARAMETERS en {fn}')

def final_positions(out, nat):
    s = open(out).read()
    conv = 'bfgs converged' in s
    blocks = [m.start() for m in re.finditer(r'ATOMIC_POSITIONS', s)]
    if not blocks:
        return None, conv
    L = s[blocks[-1]:].splitlines()
    assert 'angstrom' in L[0].lower(), 'se esperaban posiciones en angstrom: ' + L[0]
    pos = []
    for l in L[1:]:
        t = l.split()
        if len(t) >= 4 and re.match(r'^[A-Za-z]{1,2}$', t[0]):
            pos.append((t[0], *map(float, t[1:4])))
            if len(pos) == nat: break
        elif pos: break
    return pos, conv

def write_in(fn, prefix, pos, cell, local_tf=False):
    sp = [e for e in ('Au', 'C', 'O', 'N', 'H') if any(p[0] == e for p in pos)]
    mix = ("  mixing_beta      = 0.2\n  mixing_mode      = 'local-TF'\n" if local_tf else
           "  mixing_beta      = 0.05\n  mixing_mode      = 'plain'\n")
    txt = f"""&CONTROL
  calculation   = 'scf'
  restart_mode  = 'from_scratch'
  prefix        = '{prefix}'
  outdir        = './tmp'
  pseudo_dir    = '/home/dagiraldo/Au/Pseudos/'
  tprnfor       = .true.
  tstress       = .true.
/
&SYSTEM
  ibrav         = 0
  nat           = {len(pos)}
  ntyp          = {len(sp)}
  ecutwfc       = 60
  ecutrho       = 600
  occupations   = 'smearing'
  smearing      = 'mv'
  degauss       = 0.01
  vdw_corr      = 'grimme-d3'
  dftd3_version = 4
  nosym         = .true.
/
&ELECTRONS
  conv_thr         = 1.0d-8
{mix}  mixing_ndim      = 16
  electron_maxstep = 800
  diagonalization  = 'david'
  startingwfc      = 'atomic+random'
/
ATOMIC_SPECIES
""" + ''.join(f'  {e:2s}  {PSEUDO[e][0]:>7s}  {PSEUDO[e][1]}\n' for e in sp)
    txt += '\nK_POINTS automatic\n  4 4 1  0 0 0\n\nCELL_PARAMETERS angstrom\n'
    txt += ''.join('  %s  %s  %s\n' % tuple(v) for v in cell)
    txt += '\nATOMIC_POSITIONS angstrom\n' + ''.join('%-3s %16.10f %16.10f %16.10f\n' % p for p in pos)
    open(fn, 'w').write(txt)

RUN = """#!/bin/bash
#SBATCH --job-name=ABC_{tag}
#SBATCH --partition=multi
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=2
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH --output=ABC.log
#SBATCH --error=ABC.err
module load quantum-espresso/7.2
ulimit -s unlimited
export PW_TIMEOUT=140000 STEP_TIMEOUT=7200 GRACE=600
bash post_electronic.sh
"""

confs = sys.argv[1:] or sorted(d for d in glob.glob('FA_*') if os.path.isdir(d))
for c in confs:
    rin, rout = f'{c}/relax.in', f'{c}/relax.out'
    if not os.path.exists(rout):
        print(f'{c:13s} sin relax.out, se omite'); continue
    nat = int(re.search(r'nat\s*=\s*(\d+)', open(rin).read()).group(1))
    pos, conv = final_positions(rout, nat)
    if not pos or len(pos) != nat:
        print(f'{c:13s} no se pudieron leer {nat} posiciones, se omite'); continue
    if not conv:
        print(f'{c:13s} AVISO: la relajación NO dice "bfgs converged"; no se generan SCF (reiniciar primero)'); continue
    cell = cell_from_in(rin)
    au = [p for p in pos if p[0] == 'Au']; mol = [p for p in pos if p[0] != 'Au']
    d = os.path.join(DEST, 'FAfix_' + c.replace('FA_', '')); os.makedirs(d, exist_ok=True)
    write_in(f'{d}/A.in', 'A', pos, cell)
    write_in(f'{d}/B.in', 'B', au, cell)
    write_in(f'{d}/C.in', 'C', mol, cell, local_tf=True)
    open(f'{d}/run_ABC.sh', 'w').write(RUN.format(tag=c.replace('FA_', '')))
    if POST: shutil.copy(POST, d)
    zt = max(p[3] for p in au); hmin = min(p[3] for p in mol) - zt
    print(f'{c:13s} OK -> {d}   nAu={len(au)} nmol={len(mol)}  átomo más bajo de FA a {hmin:.2f} Å sobre el Au'
          + ('' if POST else '   (copiar post_electronic.sh a mano)'))
