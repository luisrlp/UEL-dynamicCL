#!/usr/bin/env bash
# Refit with the Langmuir (1 - theta_b) kinetics (material_TC.py, kineticsFunc.f90; branch tc-langmuir-kinetics).
# Run from python_scripts/ (on the VM).
#
# Starts: the A1 (old-law) optima with koff0 rescaled by 1/(1 - theta_b0), i.e. the same rest turnover:
#   L1044:   lambda0 1.044 fixed, s_taut 1.09416, koff0 0.06062 -> 0.0995, phi 0.99967, dxs 3.085 nm, keq 64.45 (free)
#   Pkeq_45: lambda0 1.04705 free, s_taut 1.0952, koff0 0.06108 -> 0.0997, phi 0.99985, dxs 3.417 nm, keq 45 fixed
# final3_L1044 is the candidate final set; final3_Pkeq_45 gives the near-minimum J, to check that lambda0 = 1.044
# stays within the 95 % limit (112 dJ/Jmin <= 3.84, i.e. +3.4 %; old law: +1.9 %).
# Same settings as the A1 fits: Fortran engine, factor 6, dt 0.025 s, 7 cycles, 7 sweep amplitudes, sweep weight 1.5,
# dxc = 0.1 nm and eta = 1 fixed.
#
# WORKERS processes per fit (default 8), MAXJOBS = min(cores, MAXCORES) / WORKERS fits at a time
# (e.g. MAXCORES=16 WORKERS=8 runs both at once). Fits with a result.json are skipped.
# Outputs optimisation_runs/final3_*, logs final3_*.log, queue log stageLG.log (records the git commit).
set -e
if [ -z "$STAGELG_DETACHED" ]; then
  STAGELG_DETACHED=1 nohup bash "$0" "$@" > stageLG.log 2>&1 &
  echo "Langmuir refit queue running in the background (pid $!); progress in stageLG.log, fits in final3_*.log"
  exit 0
fi
WORKERS=${WORKERS:-8}
MAXCORES=${MAXCORES:-8}   # upper limit on the cores used by all fits together (MAXJOBS * WORKERS)
NCORES=$(nproc)
[ "$NCORES" -le "$MAXCORES" ] || NCORES=$MAXCORES
MAXJOBS=$(( NCORES / WORKERS ))
[ "$MAXJOBS" -ge 1 ] || MAXJOBS=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ -z "$PY" ]; then if [ -x ../venv/bin/python ]; then PY=../venv/bin/python; else PY=python; fi; fi
FACTOR=6
FDIR=../10_CL_DIFF_CS_SFR_FI_TC
echo "$(date '+%F %T') WORKERS=$WORKERS MAXJOBS=$MAXJOBS (MAXCORES=$MAXCORES, $(nproc --all) cores), $PY"
echo "$(date '+%F %T') code: $(git rev-parse --abbrev-ref HEAD) $(git rev-parse HEAD)"
[ -z "$(git status --porcelain --untracked-files=no -- $FDIR/src $FDIR/global.f90 $FDIR/main_closed.f90 \
        material_TC.py optimise_laos_v2.py)" ] || echo "$(date '+%F %T') WARNING: uncommitted changes in the model sources"

# --- build the closed-system driver for this factor if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f$FACTOR ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f$FACTOR)" ]; then
  bash build_laos_closed.sh $FACTOR
fi

BASE="--consistent-kinetics --engine fortran --dt 0.025 --factor $FACTOR --cycles 7 --x0 eta=1.0"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 --sweep-weight 1.5 --max-nfev 60"
PRUNED="--fix dxc_nm=0.1"

launch() {  # name, options; keeps at most MAXJOBS fits running; skips fits that already have a result
  if [ -f optimisation_runs/final3_$1/result.json ]; then echo "$(date '+%F %T') skip final3_$1 (done)"; return; fi
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n || true; done
  nohup $PY -B optimise_laos_v2.py $2 --workers $WORKERS --out optimisation_runs/final3_$1 > final3_$1.log 2>&1 &
  echo "$(date '+%F %T') started final3_$1 (pid $!)"
}

launch L1044   "$BASE $G7 $PRUNED --fix lambda0=1.044 --free keq --bound keq=3,100 --x0 keq=64.45 \
                --x0 s_taut=1.09416 --x0 koff0=0.0995 --x0 phi=0.99967 --x0 dxs_nm=3.085"
launch Pkeq_45 "$BASE $G7 $PRUNED --x0 keq=45 --x0 lambda0=1.04705 \
                --x0 s_taut=1.0952 --x0 koff0=0.0997 --x0 phi=0.99985 --x0 dxs_nm=3.417"
wait
echo "$(date '+%F %T') Langmuir refits finished"