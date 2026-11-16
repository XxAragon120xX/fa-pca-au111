#!/bin/bash
#SBATCH --partition=multi
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=2
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH --output=pipeline.log
#SBATCH --open-mode=append
# =============================================================================
# pipeline.sh — relax -> SCF A/B/C -> post-procesado, en un solo trabajo que
# se reenvía solo hasta terminar. Cada etapa se salta si ya está hecha.
#
# Uso (dentro de ~/Au/IN/FA_fix/<conf>, p. ej. FA_V-Z-Pt):
#     sbatch -J P_V-Z-Pt pipeline.sh
# Para que arranque cuando termine una relajación que ya está corriendo:
#     sbatch -J P_V-Z-Pt --dependency=afterany:<JOBID> pipeline.sh
# Estado:   cat pipeline.log      Detener: touch STOP (no se reenvía más)
# =============================================================================
set -u
module load quantum-espresso/7.2 2>/dev/null
ulimit -s unlimited
CONF=$(basename "$PWD")                       # FA_V-Z-Pt
TAG=${CONF#FA_}                               # V-Z-Pt
FIX=$(dirname "$PWD")                         # ~/Au/IN/FA_fix
SCF=${SCF_ROOT:-$HOME/Au/SCF_final}/FAfix_$TAG
MPI=${MPI-"srun -n ${SLURM_NTASKS:-64}"}
PW=${PW:-pw.x}
MARGIN=${MARGIN:-2400}                        # s que se reservan antes del límite
POST_MIN=${POST_MIN:-21600}                   # s mínimos para lanzar el post-procesado (6 h)
MAXRESUB=${MAXRESUB:-6}
T0=$(date +%s)
LIMIT=${LIMIT:-172800}                        # 2 días
if [ -n "${SLURM_JOB_ID:-}" ]; then           # tiempo real que le queda al trabajo
  L=$(squeue -h -j "$SLURM_JOB_ID" -o %L 2>/dev/null)
  if [[ $L =~ ^(([0-9]+)-)?(([0-9]+):)?([0-9]+):([0-9]+)$ ]]; then
    LIMIT=$(( ${BASH_REMATCH[2]:-0}*86400 + 10#${BASH_REMATCH[4]:-0}*3600 + 10#${BASH_REMATCH[5]}*60 + 10#${BASH_REMATCH[6]} ))
  fi
fi
left(){ echo $(( LIMIT - ($(date +%s) - T0) )); }
log(){ echo "[$(date '+%F %T')] [$CONF] $*"; }
done_ok(){ grep -q "JOB DONE" "$1" 2>/dev/null && grep -q "$2" "$1" 2>/dev/null; }
set_ctrl(){ # set_ctrl archivo clave valor  -> fija/añade una variable en &CONTROL
  python3 - "$1" "$2" "$3" <<'PY'
import re,sys
f,k,v=sys.argv[1:]; s=open(f).read()
pat=re.compile(r'^\s*%s\s*=.*$'%k,re.M|re.I)
s=pat.sub('  %s = %s'%(k,v),s) if pat.search(s) else re.sub(r'(&CONTROL[^\n]*\n)',r'\1  %s = %s\n'%(k,v),s,count=1,flags=re.I)
open(f,'w').write(s)
PY
}
resubmit(){
  [ -f STOP ] && { log "STOP presente: no se reenvía"; exit 0; }
  N=$(cat "$FIX/$CONF/.resub" 2>/dev/null || echo 0); N=$((N+1)); echo $N > "$FIX/$CONF/.resub"
  [ $N -gt $MAXRESUB ] && { log "ERROR: $MAXRESUB reenvíos sin terminar; revisar a mano"; exit 1; }
  if [ -n "${SLURM_JOB_ID:-}" ]; then
    J=$(sbatch --parsable -J "P_$TAG" --chdir="$FIX/$CONF" "$FIX/$CONF/pipeline.sh"); log "reenviado (#$N) como trabajo $J"
  else log "(modo prueba) aquí se reenviaría (#$N)"; fi
  exit 0
}
run_pw(){ # run_pw input output -> corre pw.x con max_seconds para cortar limpio antes del límite
  local secs=$(( $(left) - MARGIN ))
  [ $secs -lt 1800 ] && return 2
  set_ctrl "$1" max_seconds "$secs"
  log "pw.x $1 -> $2 (max_seconds=$secs)"
  $MPI $PW -npool 4 -in "$1" > "$2" 2>&1
}

log "inicio | tiempo disponible $(( LIMIT/3600 )) h | SCF en $SCF"

# ---------------------------------------------------------------- 1. relajación
if ! grep -q "bfgs converged" relax.out 2>/dev/null; then
  INP=relax.in
  if [ -s relax.out ]; then                   # hubo una corrida que no terminó: seguir desde la última geometría
    python3 - <<'PY'
import re
s=open('relax.out').read(); b=[m.start() for m in re.finditer(r'ATOMIC_POSITIONS',s)]
inp=open('relax.in').read()
if b:
    blk=s[b[-1]:].splitlines(); pos=[]
    for l in blk[1:]:
        t=l.split()
        if len(t)>=4 and re.match(r'^[A-Za-z]{1,2}$',t[0]): pos.append(l)
        elif pos: break
    head,tail=re.split(r'ATOMIC_POSITIONS[^\n]*\n',inp,maxsplit=1)
    old=[l for l in tail.splitlines() if l.strip()]
    nat=len(pos); rest='\n'.join(old[nat:])
    inp=head+'ATOMIC_POSITIONS angstrom\n'+'\n'.join(pos)+'\n'+rest+'\n'
open('relax_cont.in','w').write(re.sub(r"restart_mode\s*=\s*'[^']*'","restart_mode  = 'from_scratch'",inp))
PY
    INP=relax_cont.in
    cat relax.out >> relax_hist.out; log "relajación incompleta: se continúa desde la última geometría"
  fi
  run_pw $INP relax.out; rc=$?
  [ $rc = 2 ] && resubmit
  grep -q "bfgs converged" relax.out || { log "relajación sin terminar en este trabajo"; resubmit; }
  log "relajación convergida: $(grep 'bfgs converged' relax.out | tail -1)"
else log "relajación ya convergida: se salta"; fi

# ---------------------------------------------------------------- 2. inputs A/B/C
if [ ! -f "$SCF/A.in" ]; then
  ( cd "$FIX" && python3 make_scf.py "$CONF" ) || { log "ERROR en make_scf.py"; exit 1; }
fi
cd "$SCF" || { log "ERROR: no existe $SCF"; exit 1; }

# ---------------------------------------------------------------- 3. SCF A, B, C
for s in A B C; do
  if done_ok $s.out "convergence has been achieved"; then log "SCF $s ya convergido: se salta"; continue; fi
  if [ -s $s.out ] && grep -q "Maximum CPU time exceeded" $s.out; then
    set_ctrl $s.in restart_mode "'restart'"; cp $s.out $s.out.parte_$(date +%m%d%H%M); log "SCF $s se reanuda (restart)"
  fi
  run_pw $s.in $s.out; rc=$?
  [ $rc = 2 ] && resubmit
  done_ok $s.out "convergence has been achieved" || { grep -q "Maximum CPU time exceeded" $s.out && resubmit; log "ERROR: SCF $s no convergió; revisar $SCF/$s.out"; exit 1; }
  log "SCF $s convergido: $(grep '^!' $s.out | tail -1)"
done

# ---------------------------------------------------------------- 4. post-procesado
if [ ! -f .post_done ]; then
  [ $(left) -lt $POST_MIN ] && { log "poco tiempo para el post-procesado"; resubmit; }
  RUN_PW=0 GRACE=600 bash post_electronic.sh && touch .post_done
fi
log "TERMINADO: E(A),E(B),E(C) y post-procesado listos en $SCF"
