#!/bin/bash
# ============================================================================
# post_electronic.sh  —  post-procesado electrónico QE para adsorbato/Au(111)
# Metodología Sección 2.4: DOS (dos.x), PDOS + spilling (projwfc.x),
# densidad de valencia (pp.x, cube) + Bader, Δρ = ρ(A) − ρ(B) − ρ(C),
# Δρ(z) y potencial electrostático promediados en z (average.x).
#
# Uso (dentro de la carpeta del sistema):
#   bash post_electronic.sh            -> A B C y Δρ (salta lo ya hecho)
#   bash post_electronic.sh A          -> solo el sistema A
#   bash post_electronic.sh C drho     -> C y luego Δρ
# Variables de entorno opcionales: NP, MPI, SERIAL, BIN, BADER
#   (watch_post.sh las usa para lanzar todo dentro de un trabajo que ya corre)
# ============================================================================
set -u

# ----------------------------- AJUSTAR --------------------------------------
COMPLEX=A          # complejo PCA/Au(111)
SLAB=B             # slab limpio congelado
MOL=C              # molécula aislada congelada
NP=${NP:-${SLURM_NTASKS:-16}}
MPI=${MPI:-"srun -n ${NP}"}       # lanzador paralelo
SERIAL=${SERIAL:-""}              # lanzador serial (average.x, bader)
BIN=${BIN:-""}                    # ruta a binarios QE; vacío = PATH
BADER=${BADER:-bader}
DE=0.02            # paso de energía DOS/PDOS (eV)
AWIN=4.45          # ventana del promedio macroscópico (bohr) ≈ d interplanar Au(111)
# ----------------------------------------------------------------------------

STEPS=${*:-"$COMPLEX $SLAB $MOL drho"}
log(){ echo "[$(date '+%F %T')] $*"; }
getstr(){ grep -oiE "$1[[:space:]]*=[[:space:]]*['\"][^'\"]*['\"]" "$2" | head -1 \
          | sed -E "s/.*['\"]([^'\"]*)['\"]$/\1/"; }
GRACE=${GRACE:-600}               # s de espera tras "JOB DONE" antes de forzar el cierre
STEP_TIMEOUT=${STEP_TIMEOUT:-3600} # s máximos por paso
ok(){ grep -q "JOB DONE" "$1" 2>/dev/null; }
run(){
  # Lanza un paso y lo vigila: si el programa ya escribió "JOB DONE" pero el proceso
  # (srun/MPI) no termina, lo cierra tras GRACE s. Evita cuelgues en MPI_Finalize.
  local lg=$1; shift; log "-> $* (log: $lg)"
  "$@" > "$lg" 2>&1 &
  local pid=$! t0=$(date +%s) td=0 rc
  while kill -0 $pid 2>/dev/null; do
    sleep 10
    [ $td = 0 ] && ok "$lg" && td=$(date +%s)
    if [ $td != 0 ] && [ $(( $(date +%s) - td )) -gt $GRACE ]; then
      log "   $lg ya dice JOB DONE pero el proceso no termina: se fuerza el cierre"
      kill -TERM $pid 2>/dev/null; sleep 15; kill -KILL $pid 2>/dev/null; break
    fi
    if [ $(( $(date +%s) - t0 )) -gt $STEP_TIMEOUT ]; then
      log "   TIMEOUT de paso (${STEP_TIMEOUT}s) en $lg"
      kill -TERM $pid 2>/dev/null; sleep 15; kill -KILL $pid 2>/dev/null; break
    fi
  done
  wait $pid 2>/dev/null; rc=$?
  if ok "$lg" || [ $rc = 0 ]; then return 0; fi
  log "   ERROR (código $rc), revisar $lg"; return 1
}

log "Post-procesado en $(pwd) | pasos: $STEPS | MPI='$MPI'"

declare -A PREFIX OUTDIR OK
for s in $COMPLEX $SLAB $MOL; do
  OK[$s]=0; PREFIX[$s]=""; OUTDIR[$s]=""
  [ -f "$s.in" ] || { log "ERROR: falta $s.in"; continue; }
  p=$(getstr prefix "$s.in"); o=$(getstr outdir "$s.in")
  PREFIX[$s]=${p:-pwscf}; OUTDIR[$s]=${o:-${ESPRESSO_TMPDIR:-./}}
  if grep -q "JOB DONE" "$s.out" 2>/dev/null && grep -q "convergence has been achieved" "$s.out"; then
    OK[$s]=1
  fi
done

# ---- SCF faltantes: si un sistema no tiene SCF convergido y existe su .in, se corre pw.x aquí
RUN_PW=${RUN_PW:-1}
PW_TIMEOUT=${PW_TIMEOUT:-10800}
if [ "$RUN_PW" = 1 ]; then
  for s in $COMPLEX $SLAB $MOL; do
    [ "${OK[$s]}" = 1 ] && continue
    [ -f "$s.in" ] || continue
    mkdir -p "${OUTDIR[$s]}"
    [ -f "$s.out" ] && mv "$s.out" "$s.out.fallido_$(date +%m%d%H%M)"
    log "===== SCF de $s no convergido o ausente: corriendo pw.x ====="
    STEP_TIMEOUT=$PW_TIMEOUT run $s.out $MPI ${BIN}pw.x -npool 4 -in $s.in
    if grep -q "JOB DONE" "$s.out" 2>/dev/null && grep -q "convergence has been achieved" "$s.out"; then
      OK[$s]=1; log "   SCF $s convergido: $(grep '^!' $s.out | tail -1)"
    else
      log "   ERROR: SCF $s no convergió; revisar $s.out"
    fi
  done
fi

for a in $COMPLEX $SLAB $MOL; do for b in $COMPLEX $SLAB $MOL; do
  if [ "$a" \< "$b" ] && [ -n "${PREFIX[$a]}" ] && \
     [ "${OUTDIR[$a]%/}/${PREFIX[$a]}" = "${OUTDIR[$b]%/}/${PREFIX[$b]}" ]; then
    log "ERROR: $a y $b comparten outdir/prefix; abortando."; exit 1
  fi
done; done

command -v "$BADER" >/dev/null 2>&1 || { log "AVISO: '$BADER' no está en el PATH; se omite Bader"; BADER=""; }
NR3=$(grep -m1 "FFT dimensions" "$COMPLEX.out" 2>/dev/null | sed -E 's/.*\([[:space:]]*[0-9]+,[[:space:]]*[0-9]+,[[:space:]]*([0-9]+)\).*/\1/')
NR3=${NR3:-720}

post_system(){
  local s=$1 p=${PREFIX[$1]} o=${OUTDIR[$1]}
  if [ -f ".post_done_$s" ]; then log "$s ya post-procesado (.post_done_$s); se salta"; return; fi
  if [ "${OK[$s]}" != 1 ]; then log "AVISO: $s.out sin SCF convergido/terminado; se omite $s"; return; fi
  log "===== Sistema $s (prefix=$p, outdir=$o) ====="

  cat > dos_$s.in <<EOF
&dos
  prefix='$p', outdir='$o', fildos='$s.dos', DeltaE=$DE
/
EOF
  if [ -f "$s.dos" ] && ok dos_$s.out; then log "   dos $s ya hecho"; else
  run dos_$s.out $MPI ${BIN}dos.x -inp dos_$s.in; fi

  if ls "${o%/}/$p.save"/wfc* >/dev/null 2>&1; then
    mkdir -p pdos_$s
    cat > projwfc_$s.in <<EOF
&projwfc
  prefix='$p', outdir='$o', filpdos='pdos_$s/$s', DeltaE=$DE
/
EOF
    if ok projwfc_$s.out; then log "   projwfc $s ya hecho"; else
    run projwfc_$s.out $MPI ${BIN}projwfc.x -inp projwfc_$s.in; fi
    log "   $(grep -i 'spilling' projwfc_$s.out | tail -1)"
  else
    log "AVISO: no hay wfc* en ${o%/}/$p.save; se omite projwfc para $s"
  fi

  cat > pp_rho_$s.in <<EOF
&inputpp
  prefix='$p', outdir='$o', filplot='$s.rho', plot_num=0
/
&plot
  iflag=3, output_format=6, fileout='${s}_rho.cube'
/
EOF
  if [ -f "${s}_rho.cube" ] && [ -f "$s.rho" ] && ok pp_rho_$s.out; then log "   densidad $s ya hecha"; else
  run pp_rho_$s.out $MPI ${BIN}pp.x -inp pp_rho_$s.in; fi

  cat > pp_v_$s.in <<EOF
&inputpp
  prefix='$p', outdir='$o', filplot='$s.v', plot_num=11
/
EOF
  if [ -f "$s.v" ] && ok pp_v_$s.out; then log "   potencial $s ya hecho"; else
  run pp_v_$s.out $MPI ${BIN}pp.x -inp pp_v_$s.in; fi

  if [ -f "$s.v" ]; then
    printf "1\n%s\n1.0\n%s\n3\n%s\n" "$s.v" "$NR3" "$AWIN" > avg_v_$s.in
    run avg_v_$s.out $SERIAL ${BIN}average.x -inp avg_v_$s.in
    [ -f avg.dat ] && mv avg.dat ${s}_vavg.dat
  fi

  if [ -n "$BADER" ] && [ -f "${s}_rho.cube" ] && [ ! -s "bader_$s/ACF.dat" ]; then
    mkdir -p bader_$s
    ( cd bader_$s && $SERIAL $BADER ../${s}_rho.cube > bader.out 2>&1 ) \
      && log "   Bader $s listo" || log "   ERROR en Bader $s"
  fi

  if [ -f "$s.dos" ] && [ -f "${s}_rho.cube" ] && [ -f "${s}_vavg.dat" ] && \
     ok projwfc_$s.out; then
    touch ".post_done_$s"; log "   $s completo"
  else
    log "   $s incompleto; se repetirá en la próxima ejecución"
  fi
}

post_drho(){
  if [ -f drho.cube ] && [ -f drho_zavg.dat ]; then log "Δρ ya calculado; se salta"; return; fi
  for f in $COMPLEX.rho $SLAB.rho $MOL.rho; do
    [ -f "$f" ] || { log "AVISO: falta $f; no se calcula Δρ"; return; }
  done
  log "===== Δρ = ρ($COMPLEX) − ρ($SLAB) − ρ($MOL) ====="
  cat > pp_drho.in <<EOF
&inputpp
/
&plot
  nfile=3,
  filepp(1)='$COMPLEX.rho', weight(1)=1.0,
  filepp(2)='$SLAB.rho',    weight(2)=-1.0,
  filepp(3)='$MOL.rho',     weight(3)=-1.0,
  iflag=3, output_format=6, fileout='drho.cube'
/
EOF
  if [ -f drho.cube ] && ok pp_drho.out; then log "   drho.cube ya hecho"; else
  run pp_drho.out $MPI ${BIN}pp.x -inp pp_drho.in; fi
  printf "3\n%s\n1.0\n%s\n-1.0\n%s\n-1.0\n%s\n3\n%s\n" \
    "$COMPLEX.rho" "$SLAB.rho" "$MOL.rho" "$NR3" "$AWIN" > avg_drho.in
  run avg_drho.out $SERIAL ${BIN}average.x -inp avg_drho.in
  [ -f avg.dat ] && mv avg.dat drho_zavg.dat
}

for step in $STEPS; do
  case $step in
    drho) post_drho ;;
    *)    post_system "$step" ;;
  esac
done

log "===== RESUMEN ====="
for s in $COMPLEX $SLAB $MOL; do
  [ "${OK[$s]}" = 1 ] || continue
  ef=$(grep -m1 "the Fermi energy is" "$s.out" | awk '{print $5}')
  sp=$(grep -i 'spilling' projwfc_$s.out 2>/dev/null | tail -1 | awk '{print $NF}')
  log "$s: E_F=${ef:-?} eV  spilling=${sp:-?}  hecho=$([ -f .post_done_$s ] && echo sí || echo no)"
done
log "Fin"
