#!/usr/bin/env bash
# Final C4 optimisation, Stage A (multi-start screening at production settings), see paper/paper_plan.md.
# Run from python_scripts/ on the VM. Needs the current optimise_laos_v2.py (--consistent-kinetics, LSQ tol 1e-5,
# eta parameter), amplitude_sweep_runs.py and material_TC.py. Four independent fits, 5 workers each (adapt WORKERS).
#   C4 kinetics, point stress, dt 0.1 s, factor 3, 5 cycles,
#   8 free parameters: lambda0, s_taut, koff0, phi, dxc_nm, dxs_nm, keq, eta (beta 0.438 and mu0stretch 38600 fixed),
#   LAOS loop (Fig. 4a) + amplitude sweep (Fig. 4b, 7 amplitudes up to 0.158), sweep weight 1.5, max-nfev 40.
# Start 1 = T1 set (c5_w15_tc_T1; dxc, keq at their literature values, eta = 1 as in all earlier fits);
# starts 2-4 = Latin hypercube in the parameter ranges of paper_plan.md (seed 20261101, k_off0 log-uniform).
# Bounds (optimise_laos_v2.PARAMS): lambda0 1.00-1.07 (all starts <= 1.05), dxc/dxs >= 0.1 nm, eta 0.6-1.0.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKERS=${WORKERS:-5}
COMMON="--consistent-kinetics --free keq --free eta \
 --sweep-weight 1.5 --sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 \
 --dt 0.1 --factor 3 --cycles 5 --max-nfev 40"

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
    --out optimisation_runs/final_A_s$k > final_A_s$k.log 2>&1 &
  echo "started final_A_s$k (pid $!)"
done
