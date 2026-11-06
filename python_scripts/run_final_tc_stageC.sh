#!/usr/bin/env bash
# Final C4 optimisation, Stage C (final polish, robustness and profiles), see paper/paper_plan.md.
# Run from python_scripts/ on the VM after Stage B (needs the same files as run_final_tc_stageB.sh).
#
# Stage B outcome: one basin; B1 (5 free) J 0.0552 = C1 (8 free) 0.0551, C2 (s4 start) merges into it;
# C3 (Guo catch, dxc 1.5 nm) J 0.0690. phi -> 0.999 (near the bound) with |corr(phi, dxs)| 0.96-0.99: only the
# onset force of the slip arm is identified. Last-cycle change 0.010-0.015 Pa with 5 cycles; B2-B4 at max-nfev.
#
# Fixed throughout: keq 14, eta 1, dxc 0.1 nm (unless profiled), beta 0.438, mu0stretch 38600 pN.
#   F1  / F2   B1 / B2 polish from their optima: 7 cycles, max-nfev 100, 5 free           (final set candidates)
#   F1p / F2p  same with phi = 0.999 fixed (4 free; dxs carries the slip onset force)    (removes the degeneracy)
#   M01-M10    multi-start, B1 settings (5 cycles, 7 amplitudes, w 1.5), Latin hypercube in the pruned space
#              (seed 20261105; lambda0 1.00-1.05, s_taut 1.05-1.20, koff0 0.03-0.3 log-uniform, phi 0.50-0.99,
#               dxs 0.5-4 nm), max-nfev 60                                                 (is the B1 basin global?)
#   Pdxc_*     profile of dxc (0.2 0.3 0.5 0.75 1.0 2.0 nm; B1 = 0.1, C3 = 1.5), B1 settings, start B1 (dxc <= 0.5)
#              or C3 optimum (dxc >= 0.75)                                                 (catch distance allowed)
#   Pkeq_*     profile of keq (5 10 20 30 45; B1 = 14), B1 settings, start B1 optimum      (keq not identified)
# Objective J = loop_RMS^2 + w^2 sweep_RMS_log10^2; J of M, Pdxc, Pkeq compares directly with B1 / C3.
#
# 25 fits x 3 workers, at most MAXJOBS fits at a time (default 7 -> 21 cores); about 10-40 min per fit.
set -e
# The script keeps running as the job queue: re-launch itself with nohup in the background, so that it
# survives the end of the SSH session (plain `bash run_final_tc_stageC.sh` is enough; log in stageC.log).
if [ -z "$STAGEC_DETACHED" ]; then
  STAGEC_DETACHED=1 nohup bash "$0" "$@" > stageC.log 2>&1 &
  echo "Stage C queue running in the background (pid $!); progress in stageC.log, fits in final_*.log"
  exit 0
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKERS=${WORKERS:-3}
MAXJOBS=${MAXJOBS:-7}
FDIR=../10_CL_DIFF_CS_SFR_FI_TC

# --- build the closed-system driver for factor 3 (laos_closed_f3) if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f3 ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f3)" ]; then
  bash build_laos_closed.sh 3
fi

BASE="--consistent-kinetics --engine fortran --dt 0.1 --factor 3 --x0 keq=14.0 --x0 eta=1.0"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158"
B1SET="$BASE $G7 --sweep-weight 1.5 --cycles 5"

# start vectors (exact optima from optimisation_runs/final_B1, _B2, _C3 / Latin hypercube)
B1X="--x0 lambda0=1.0522343393590514 --x0 s_taut=1.0954514815541176 --x0 koff0=0.06296073654516689 --x0 phi=0.9992549604520296 --x0 dxs_nm=2.5982450050334887"
B2X="--x0 lambda0=1.042127644572663 --x0 s_taut=1.0922412237244665 --x0 koff0=0.06280411786336111 --x0 phi=0.9994024517190676 --x0 dxs_nm=2.3268457527245388"
B1X_NOPHI="--x0 lambda0=1.0522343393590514 --x0 s_taut=1.0954514815541176 --x0 koff0=0.06296073654516689 --x0 dxs_nm=2.5982450050334887"
B2X_NOPHI="--x0 lambda0=1.042127644572663 --x0 s_taut=1.0922412237244665 --x0 koff0=0.06280411786336111 --x0 dxs_nm=2.3268457527245388"
C3X="--x0 lambda0=1.0547979541948291 --x0 s_taut=1.098914413510173 --x0 koff0=0.06742489042731893 --x0 phi=0.6980068565292581 --x0 dxs_nm=0.7532110932812232"
M01="--x0 lambda0=1.046189 --x0 s_taut=1.110601 --x0 koff0=0.191916 --x0 phi=0.718459 --x0 dxs_nm=2.040256"
M02="--x0 lambda0=1.024088 --x0 s_taut=1.109205 --x0 koff0=0.045093 --x0 phi=0.904903 --x0 dxs_nm=2.651349"
M03="--x0 lambda0=1.007943 --x0 s_taut=1.12723 --x0 koff0=0.066817 --x0 phi=0.862872 --x0 dxs_nm=0.722873"
M04="--x0 lambda0=1.013317 --x0 s_taut=1.07426 --x0 koff0=0.112623 --x0 phi=0.567181 --x0 dxs_nm=2.311446"
M05="--x0 lambda0=1.035686 --x0 s_taut=1.157178 --x0 koff0=0.037572 --x0 phi=0.525133 --x0 dxs_nm=3.054334"
M06="--x0 lambda0=1.034944 --x0 s_taut=1.141366 --x0 koff0=0.184554 --x0 phi=0.676589 --x0 dxs_nm=3.548874"
M07="--x0 lambda0=1.044337 --x0 s_taut=1.054116 --x0 koff0=0.135623 --x0 phi=0.626509 --x0 dxs_nm=1.128775"
M08="--x0 lambda0=1.025658 --x0 s_taut=1.088805 --x0 koff0=0.054959 --x0 phi=0.958263 --x0 dxs_nm=1.213403"
M09="--x0 lambda0=1.018961 --x0 s_taut=1.181953 --x0 koff0=0.078777 --x0 phi=0.752004 --x0 dxs_nm=1.868239"
M10="--x0 lambda0=1.001155 --x0 s_taut=1.185379 --x0 koff0=0.260518 --x0 phi=0.839141 --x0 dxs_nm=3.791113"

launch() {  # name, options; keeps at most MAXJOBS fits running
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n || true; done
  nohup python -B optimise_laos_v2.py $2 --workers $WORKERS --out optimisation_runs/final_$1 > final_$1.log 2>&1 &
  echo "$(date '+%F %T') started final_$1 (pid $!)"
}

# final polish
launch F1  "$BASE $G7 --sweep-weight 1.5 --cycles 7 --max-nfev 100 --fix dxc_nm=0.1 $B1X"
launch F2  "$BASE $G7 --sweep-weight 1.0 --cycles 7 --max-nfev 100 --fix dxc_nm=0.1 $B2X"
launch F1p "$BASE $G7 --sweep-weight 1.5 --cycles 7 --max-nfev 100 --fix dxc_nm=0.1 --fix phi=0.999 $B1X_NOPHI"
launch F2p "$BASE $G7 --sweep-weight 1.0 --cycles 7 --max-nfev 100 --fix dxc_nm=0.1 --fix phi=0.999 $B2X_NOPHI"
# multi-start
for k in 01 02 03 04 05 06 07 08 09 10; do
  eval X=\$M$k
  launch M$k "$B1SET --max-nfev 60 --fix dxc_nm=0.1 $X"
done
# profile of the catch distance
for v in 0.2 0.3 0.5; do launch Pdxc_$v "$B1SET --max-nfev 40 --fix dxc_nm=$v $B1X"; done
for v in 0.75 1.0 2.0; do launch Pdxc_$v "$B1SET --max-nfev 40 --fix dxc_nm=$v $C3X"; done
# profile of keq (--x0 keq later on the line overrides the BASE value)
for v in 5 10 20 30 45; do launch Pkeq_$v "$B1SET --max-nfev 40 --fix dxc_nm=0.1 $B1X --x0 keq=$v"; done
wait
echo "$(date '+%F %T') Stage C finished"
