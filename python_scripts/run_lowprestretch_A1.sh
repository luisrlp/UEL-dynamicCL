#!/usr/bin/env bash
# Lower prestretch with the A1 reservoir (constant-tension reservoir, force target; laos_closed_f6 built from folder 10).
# Run from python_scripts/ (locally or on the VM).
#
# A1 optimum (final2_B1_A1, J7 0.04917): lambda0 1.0514, s_taut 1.0958, koff0 0.0640, phi 0.9991, dxs 2.62 nm.
# The 8-parameter check (final2_C1_A1) gives corr(lambda0, keq) = -0.92: a larger keq lowers lambda0 at almost no
# cost, as in the legacy K_eq profile (final2_Pkeq_45: lambda0 1.0463 at keq 45, J7 -1.2 % vs final2_B1).
#   A1_Pkeq_45   keq fixed at 45, lambda0 free           (A1 version of final2_Pkeq_45)
#   A1_L1044     lambda0 fixed at 1.044, keq free in [3, 100]
#   A1_L1040     lambda0 fixed at 1.040, keq free in [3, 100]
#   A1_L1042     lambda0 fixed at 1.042, keq free in [3, 100]   (fills in the profile)
# The fixed-lambda0 fits give the cost profile J7(lambda0) with the other parameters (and keq) refitted, so the lowest
# lambda0 compatible with the data can be read off it.
# Starts: legacy Pkeq_45 optimum shifted by the B1 -> B1_A1 change (lambda0 +0.0008, koff0 +0.001, dxs -0.6 nm).
# Same settings as final2_B1_A1: Fortran engine, factor 6, dt 0.025 s, 7 cycles, 7 sweep amplitudes, sweep weight 1.5,
# dxc = 0.1 nm and eta = 1 fixed.
#
# WORKERS processes per fit (default 8), MAXJOBS = min(cores, MAXCORES) / WORKERS fits at a time (MAXCORES default 8:
# one fit at a time locally; on the VM e.g. MAXCORES=24 WORKERS=6 runs all four at once).
# Fits with a result.json are skipped. Outputs optimisation_runs/final2_A1_*, logs final2_A1_*.log, queue log stageLP_A1.log.
set -e
if [ -z "$STAGELP_DETACHED" ]; then
  STAGELP_DETACHED=1 nohup bash "$0" "$@" > stageLP_A1.log 2>&1 &
  echo "A1 low-prestretch queue running in the background (pid $!); progress in stageLP_A1.log, fits in final2_A1_*.log"
  exit 0
fi
WORKERS=${WORKERS:-8}
MAXCORES=${MAXCORES:-8}   # upper limit on the cores used by all fits together (MAXJOBS * WORKERS)
NCORES=$(nproc)   # before the export below: GNU nproc returns OMP_NUM_THREADS when it is set
[ "$NCORES" -le "$MAXCORES" ] || NCORES=$MAXCORES
MAXJOBS=$(( NCORES / WORKERS ))
[ "$MAXJOBS" -ge 1 ] || MAXJOBS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ -z "$PY" ]; then if [ -x ../venv/bin/python ]; then PY=../venv/bin/python; else PY=python; fi; fi
FACTOR=6
FDIR=../10_CL_DIFF_CS_SFR_FI_TC
echo "$(date '+%F %T') WORKERS=$WORKERS MAXJOBS=$MAXJOBS: at most $(( MAXJOBS * WORKERS )) busy processes (MAXCORES=$MAXCORES, $(nproc --all) cores), $PY"

# --- build the closed-system driver for this factor if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f$FACTOR ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f$FACTOR)" ]; then
  bash build_laos_closed.sh $FACTOR
fi

BASE="--consistent-kinetics --engine fortran --dt 0.025 --factor $FACTOR --cycles 7 --x0 eta=1.0"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 --sweep-weight 1.5 --max-nfev 40"
PRUNED="--fix dxc_nm=0.1"
REST="--x0 s_taut=1.0953 --x0 koff0=0.0613 --x0 phi=0.999 --x0 dxs_nm=2.95"
KEQFREE="--free keq --bound keq=3,100 --x0 keq=60"

launch() {  # name, options; keeps at most MAXJOBS fits running; skips fits that already have a result
  if [ -f optimisation_runs/final2_$1/result.json ]; then echo "$(date '+%F %T') skip final2_$1 (done)"; return; fi
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n || true; done
  nohup $PY -B optimise_laos_v2.py $2 --workers $WORKERS --out optimisation_runs/final2_$1 > final2_$1.log 2>&1 &
  echo "$(date '+%F %T') started final2_$1 (pid $!)"
}

launch A1_Pkeq_45 "$BASE $G7 $PRUNED --x0 keq=45 --x0 lambda0=1.0471 $REST"
launch A1_L1044   "$BASE $G7 $PRUNED $KEQFREE --fix lambda0=1.044 $REST"
launch A1_L1040   "$BASE $G7 $PRUNED $KEQFREE --fix lambda0=1.040 $REST"
launch A1_L1042   "$BASE $G7 $PRUNED $KEQFREE --fix lambda0=1.042 $REST"
wait
echo "$(date '+%F %T') A1 low-prestretch fits finished"
