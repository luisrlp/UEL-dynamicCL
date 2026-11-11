#!/usr/bin/env bash
# Final C4 optimisation, second pass, Stages B and C together (production fits, controls, multi-start, profiles).
# Run from python_scripts/ on the VM after the second-pass Stage A (run_final_tc_stageA.sh, final2_A_s1..s4).
#
# Second-pass Stage A (2026-11-11; factor 6, dt 0.025, 5 cycles, 8 free):
#   J = loop_RMS^2 + 1.5^2 sweep_RMS_log10^2:  s2 0.0600 | s1 0.0629 | s3 0.0663 (max-nfev) | s4 0.0682 (max-nfev)
#   s1, s2: main basin (phi 0.90-0.97, koff0 ~0.06, a_ratio 1.13-1.16); s3, s4: phi 0.55-0.73.
#   Not identified, as in the first pass: keq (std error 300-950, |corr| with lambda0 / koff0 0.95-0.98),
#   eta (|corr| with s_taut 0.96-0.99), dxc_nm (0.18-0.46 nm, std error 0.3-2.3 nm).
#   The pruned final set final_F1p_adaptive_f6 (keq 14, eta 1, phi 0.999, dxc 0.1) is 8 % better than every
#   Stage A optimum at the same settings (J 0.0546 vs 0.0597 at 7 cycles): the 8-parameter fits stall in the
#   flat valley of keq and eta. -> same pruning as the first pass: keq = 14, eta = 1, dxc = 0.1 nm fixed;
#   free: lambda0, s_taut, koff0, phi, dxs_nm.
#
# All fits at production settings: Fortran engine with adaptive steps, factor 6, output spacing dt 0.025 s,
# 7 cycles (fitted cycle steady). Stages B and C run together; the profiles and the stage-C starts use the
# current final set (final_F1p_adaptive_f6) instead of the B1 optimum, so nothing waits for Stage B.
#   B1         main configuration: 7 amplitudes, sweep weight 1.5, start P1 (s1 optimum, keq/eta/dxc fixed)
#   B1b        same from P2 (s2 optimum, projected)                                       (one basin?)
#   B2         7 amplitudes, weight 1        B3  8 amplitudes (+0.2), weight 1.5       B4  8 amplitudes, weight 1
#   C1         all 8 free, start from the final set                                       (cost of the pruning)
#   C2         all 8 free, start from the s4 optimum                                      (does the second basin survive?)
#   C3         dxc = 1.5 nm fixed (Guo & Guilford 2006 catch arm), start P1               (single-molecule catch law)
#   M01-M10    multi-start in the pruned space, Latin hypercube of the first pass (seed 20261105), max-nfev 60
#   Pdxc_*     profile of dxc (0 0.2 0.3 0.5 from the final set; 0.75 1.0 2.0 from the first-pass C3 optimum)
#   Pkeq_*     profile of keq (5 10 20 30 45), start from the final set
# Objective J = loop_RMS^2 + w^2 sweep_RMS_log10^2 (w = sweep weight); the least-squares cost is N_data * J / 2.
#
# 30 fits. Cost: about 13 CPU-min per evaluated point (loop + 7 sweep amplitudes), 30-150 points per fit.
# WORKERS processes per fit (default 6), at most MAXJOBS fits at a time (default: cores / WORKERS).
# Needs on the VM: the current 10_CL_DIFF_CS_SFR_FI_TC folder, optimise_laos_v2.py, build_laos_closed.sh,
# main.py, gfortran. Outputs optimisation_runs/final2_<name>, logs final2_<name>.log, queue log stage2BC.log.
set -e
# The script keeps running as the job queue: re-launch itself with nohup in the background, so that it
# survives the end of the SSH session (plain `bash run_final2_tc_stageBC.sh` is enough).
if [ -z "$STAGE2BC_DETACHED" ]; then
  STAGE2BC_DETACHED=1 nohup bash "$0" "$@" > stage2BC.log 2>&1 &
  echo "Stage B+C queue running in the background (pid $!); progress in stage2BC.log, fits in final2_*.log"
  exit 0
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKERS=${WORKERS:-6}
MAXJOBS=${MAXJOBS:-$(( $(nproc) / WORKERS ))}
[ "$MAXJOBS" -ge 1 ] || MAXJOBS=1
FACTOR=6
FDIR=../10_CL_DIFF_CS_SFR_FI_TC
echo "$(date '+%F %T') WORKERS=$WORKERS MAXJOBS=$MAXJOBS ($(nproc) cores)"

# --- build the closed-system driver for this factor if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f$FACTOR ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f$FACTOR)" ]; then
  bash build_laos_closed.sh $FACTOR
fi

BASE="--consistent-kinetics --engine fortran --dt 0.025 --factor $FACTOR --cycles 7 --x0 keq=14.0 --x0 eta=1.0"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158"
G8="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 0.2"
PRUNED="--fix dxc_nm=0.1"
FULL="--free keq --free eta"

# start vectors (exact optima: optimisation_runs/final2_A_s1, _s2, _s4, final_F1p_adaptive_f6, final_C3)
P1="--x0 lambda0=1.0539884641290787 --x0 s_taut=1.0959300920640591 --x0 koff0=0.06155289237084141 --x0 phi=0.8993017937100841 --x0 dxs_nm=1.1326559769270155"
P2="--x0 lambda0=1.0448200819830096 --x0 s_taut=1.0780625522518184 --x0 koff0=0.05831440323823828 --x0 phi=0.9735847576279267 --x0 dxs_nm=1.5874661183877783"
FSET="--x0 lambda0=1.0507017430729444 --x0 s_taut=1.094998268604536 --x0 koff0=0.06204682345303625 --x0 phi=0.999 --x0 dxs_nm=2.799063441642941"
FSET8="--x0 lambda0=1.0507017430729444 --x0 s_taut=1.094998268604536 --x0 koff0=0.06204682345303625 --x0 phi=0.999 --x0 dxc_nm=0.1 --x0 dxs_nm=2.799063441642941 --x0 keq=14.0 --x0 eta=1.0"
S4="--x0 lambda0=1.0477341391200912 --x0 s_taut=1.0897205825407814 --x0 koff0=0.04955263975872563 --x0 phi=0.5519885432994929 --x0 dxc_nm=0.4224608176465315 --x0 dxs_nm=0.6533241278269428 --x0 keq=44.56974331830002 --x0 eta=0.9848933673479993"
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
  nohup python -B optimise_laos_v2.py $2 --workers $WORKERS --out optimisation_runs/final2_$1 > final2_$1.log 2>&1 &
  echo "$(date '+%F %T') started final2_$1 (pid $!)"
}

# Stage B: production fits and controls (most important first)
launch B1  "$BASE $G7 --sweep-weight 1.5 --max-nfev 100 $PRUNED $P1"
launch C1  "$BASE $G7 --sweep-weight 1.5 --max-nfev 100 $FULL $FSET8"
launch B1b "$BASE $G7 --sweep-weight 1.5 --max-nfev 100 $PRUNED $P2"
launch C2  "$BASE $G7 --sweep-weight 1.5 --max-nfev 100 $FULL $S4"
launch C3  "$BASE $G7 --sweep-weight 1.5 --max-nfev 100 --fix dxc_nm=1.5 $P1"
launch B2  "$BASE $G7 --sweep-weight 1.0 --max-nfev 100 $PRUNED $P1"
launch B3  "$BASE $G8 --sweep-weight 1.5 --max-nfev 100 $PRUNED $P1"
launch B4  "$BASE $G8 --sweep-weight 1.0 --max-nfev 100 $PRUNED $P1"
# Stage C: multi-start in the pruned space
for k in 01 02 03 04 05 06 07 08 09 10; do
  eval X=\$M$k
  launch M$k "$BASE $G7 --sweep-weight 1.5 --max-nfev 60 $PRUNED $X"
done
# profile of the catch distance
for v in 0 0.2 0.3 0.5; do launch Pdxc_$v "$BASE $G7 --sweep-weight 1.5 --max-nfev 40 --fix dxc_nm=$v $FSET"; done
for v in 0.75 1.0 2.0;  do launch Pdxc_$v "$BASE $G7 --sweep-weight 1.5 --max-nfev 40 --fix dxc_nm=$v $C3X"; done
# profile of keq (--x0 keq later on the line overrides the BASE value; keq stays fixed by default)
for v in 5 10 20 30 45; do launch Pkeq_$v "$BASE $G7 --sweep-weight 1.5 --max-nfev 40 $PRUNED $FSET --x0 keq=$v"; done
wait
echo "$(date '+%F %T') Stages B+C finished"
