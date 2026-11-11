#!/usr/bin/env bash
# Final C4 optimisation, Stage A (multi-start screening), see paper/paper_plan.md. Run from python_scripts/ on the VM.
# Second pass (2026-11-11) with the Fortran engine and adaptive time stepping (the first pass, final_A_s1..s4 of
# 30 Oct, used the Python engine with fixed steps of 0.1 s and the old sub-stepping, 16-125 h per fit).
# Engine: closed-system Fortran driver 10_CL_DIFF_CS_SFR_FI_TC/main_closed.f90; its time steps are controlled by the
# PNEWDT of the network routine (Backward-Euler error estimate of the kinetics), --dt is the output spacing and the
# maximum step. Needs on the VM: the current 10_CL_DIFF_CS_SFR_FI_TC folder (main_closed.f90, src/material_AFFCL),
# optimise_laos_v2.py, build_laos_closed.sh, main.py and gfortran. The driver is rebuilt below if missing or older
# than its sources (optimise_laos_v2 stops otherwise). Four independent fits, WORKERS workers each (default 5).
#   C4 kinetics, point stress, factor FACTOR (default 3), output spacing DT (default 0.05 s), 5 cycles,
#   8 free parameters: lambda0, s_taut, koff0, phi, dxc_nm, dxs_nm, keq, eta (beta 0.438 and mu0stretch 38600 fixed),
#   LAOS loop (Fig. 4a) + amplitude sweep (Fig. 4b, 7 amplitudes up to 0.158), sweep weight 1.5, max-nfev 40.
# Same four starts as the first pass, so that the two can be compared:
# Start 1 = T1 set (c5_w15_tc_T1; dxc, keq at their literature values, eta = 1 as in all earlier fits);
# starts 2-4 = Latin hypercube in the parameter ranges of paper_plan.md (seed 20261101, k_off0 log-uniform).
# Bounds (optimise_laos_v2.PARAMS): lambda0 1.00-1.07 (all starts <= 1.05), dxc/dxs >= 0.1 nm, eta 0.6-1.0.
# Outputs: optimisation_runs/final2_A_s1..s4 (the first-pass folders final_A_s* are not touched).
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKERS=${WORKERS:-5}
FACTOR=${FACTOR:-3}
DT=${DT:-0.05}
FDIR=../10_CL_DIFF_CS_SFR_FI_TC

# --- build the closed-system driver for this factor if needed (see build_laos_closed.sh) ---
if [ ! -x $FDIR/laos_closed_f$FACTOR ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed_f$FACTOR)" ]; then
  bash build_laos_closed.sh $FACTOR
fi

COMMON="--consistent-kinetics --engine fortran --free keq --free eta \
 --sweep-weight 1.5 --sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 \
 --dt $DT --factor $FACTOR --cycles 5 --max-nfev 40"

S1="--x0 lambda0=1.0499999534924342 --x0 s_taut=1.1007813434271365 --x0 koff0=0.07876008986637978 --x0 phi=0.907 \
 --x0 dxc_nm=1.5 --x0 dxs_nm=1.3546358818091297 --x0 keq=14.0 --x0 eta=1.0"
S2="--x0 lambda0=1.02347 --x0 s_taut=1.0821 --x0 koff0=0.0445805 --x0 phi=0.665527 --x0 dxc_nm=2.86859 \
 --x0 dxs_nm=1.36601 --x0 keq=25.6808 --x0 eta=0.964595"
S3="--x0 lambda0=1.04279 --x0 s_taut=1.12831 --x0 koff0=0.114946 --x0 phi=0.614186 --x0 dxc_nm=0.76099 \
 --x0 dxs_nm=0.33344 --x0 keq=5.84984 --x0 eta=0.652369"
S4="--x0 lambda0=1.00317 --x0 s_taut=1.19452 --x0 koff0=0.179968 --x0 phi=0.859971 --x0 dxc_nm=1.43864 \
 --x0 dxs_nm=1.5438 --x0 keq=17.4516 --x0 eta=0.842001"

for k in 1 2 3 4; do
  eval X0=\$S$k
  nohup python -B optimise_laos_v2.py $COMMON $X0 --workers $WORKERS \
    --out optimisation_runs/final2_A_s$k > final2_A_s$k.log 2>&1 &
  echo "started final2_A_s$k (pid $!)"
done
