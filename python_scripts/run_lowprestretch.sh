#!/usr/bin/env bash
# Lower prestretch: refit the pruned final configuration (B1) with a lower upper bound on lambda0.
# Run from python_scripts/ on the VM (needs the --bound option of optimise_laos_v2.py).
#
# B1 optimum (final2_B1, J 0.0534): lambda0 1.0506, s_taut 1.0951, koff0 0.0630, phi 0.9997, dxs 3.20 nm.
# lambda0 stays free, only its upper bound is lowered (1.035 or 1.03); s_taut stays free in [1.05, 1.20].
# At fixed s_taut, a lower lambda0 lowers the rest force and the rest stiffness of the filaments
# (-19 % at 1.035, -25 % at 1.03; the taut stretch is unchanged); a lower s_taut brings the filament
# closer to taut and restores the stiffness, but moves the taut point to smaller strains (gamma ~0.16 instead
# of ~0.18). The starts use the s_taut that gives the B1 rest stiffness at the start lambda0
# (1.081 at lambda0 1.03, 1.077 at 1.025); keq (free in U2, U4) raises the stiffness through the bound content.
#   U1  lambda0 <= 1.035                 U2  lambda0 <= 1.035, keq free
#   U3  lambda0 <= 1.03                  U4  lambda0 <= 1.03,  keq free
# Same settings as B1: Fortran engine, factor 6, dt 0.025 s, 7 cycles, 7 sweep amplitudes, sweep weight 1.5,
# dxc = 0.1 nm, eta = 1 fixed; keq = 14 unless free.
#
# WORKERS processes per fit (default 6), MAXJOBS = min(cores, MAXCORES) / WORKERS fits at a time
# (MAXCORES default 24: never more than 24 busy processes by this queue; lower it when another queue runs).
# Fits with a result.json are skipped. Outputs optimisation_runs/final2_U<k>, logs final2_U<k>.log,
# queue log stageLP.log.
set -e
# The script keeps running as the job queue: re-launch itself with nohup in the background, so that it
# survives the end of the SSH session (plain `bash run_lowprestretch.sh` is enough).
if [ -z "$STAGELP_DETACHED" ]; then
  STAGELP_DETACHED=1 nohup bash "$0" "$@" > stageLP.log 2>&1 &
  echo "Low-prestretch queue running in the background (pid $!); progress in stageLP.log, fits in final2_U*.log"
  exit 0
fi
WORKERS=${WORKERS:-6}
MAXCORES=${MAXCORES:-24}   # upper limit on the cores used by all fits together (MAXJOBS * WORKERS)
NCORES=$(nproc)   # before the export below: GNU nproc returns OMP_NUM_THREADS when it is set
[ "$NCORES" -le "$MAXCORES" ] || NCORES=$MAXCORES
MAXJOBS=$(( NCORES / WORKERS ))
[ "$MAXJOBS" -ge 1 ] || MAXJOBS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
FACTOR=6
FDIR=../10_CL_DIFF_CS_SFR_FI_TC
echo "$(date '+%F %T') WORKERS=$WORKERS MAXJOBS=$MAXJOBS: at most $(( MAXJOBS * WORKERS )) busy processes (MAXCORES=$MAXCORES, $(nproc --all) cores)"

# --- build the closed-system driver for this factor if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f$FACTOR ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f$FACTOR)" ]; then
  bash build_laos_closed.sh $FACTOR
fi

BASE="--consistent-kinetics --engine fortran --dt 0.025 --factor $FACTOR --cycles 7 --x0 keq=14.0 --x0 eta=1.0"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158"
PRUNED="--fix dxc_nm=0.1"
REST="--x0 koff0=0.063 --x0 phi=0.999 --x0 dxs_nm=3.2"
S1035="--bound lambda0=1.0,1.035 --x0 lambda0=1.03 --x0 s_taut=1.081 $REST"
S1030="--bound lambda0=1.0,1.03 --x0 lambda0=1.025 --x0 s_taut=1.077 $REST"

launch() {  # name, options; keeps at most MAXJOBS fits running; skips fits that already have a result
  if [ -f optimisation_runs/final2_$1/result.json ]; then echo "$(date '+%F %T') skip final2_$1 (done)"; return; fi
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n || true; done
  nohup python -B optimise_laos_v2.py $2 --workers $WORKERS --out optimisation_runs/final2_$1 > final2_$1.log 2>&1 &
  echo "$(date '+%F %T') started final2_$1 (pid $!)"
}

launch U1 "$BASE $G7 --sweep-weight 1.5 --max-nfev 60 $PRUNED $S1035"
launch U2 "$BASE $G7 --sweep-weight 1.5 --max-nfev 60 $PRUNED --free keq $S1035"
launch U3 "$BASE $G7 --sweep-weight 1.5 --max-nfev 60 $PRUNED $S1030"
launch U4 "$BASE $G7 --sweep-weight 1.5 --max-nfev 60 $PRUNED --free keq $S1030"
wait
echo "$(date '+%F %T') low-prestretch fits finished"
