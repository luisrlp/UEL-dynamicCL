#!/usr/bin/env bash
# Final C4 optimisation, Stage B (production polish), see paper/paper_plan.md. Run from python_scripts/ on the VM.
#
# Stage A (final_A_s1..s4, Python engine, 30 Oct) re-evaluated with the Fortran engine (current directions and g-hat):
#   J = loop_RMS^2 + 1.5^2 sweep_RMS_log10^2:  s2 0.0581 | s3 0.0585 | s1 0.0619 | s4 0.0665 (second basin, phi 0.66)
#   s1-s3 share one basin: lambda0 1.045-1.051, a_ratio 1.12-1.15, koff0 0.056-0.059, phi 0.93-0.99, dxs 1.3-2.0 nm.
#   Not identified (pruning rule): keq (std error 170-950, |corr| with lambda0 / koff0 0.95-0.98),
#   eta (|corr| with s_taut 0.97-0.998), dxc_nm (at the lower bound 0.1 nm in s2, near it in s1, s3).
#   -> main configurations: keq = 14 (Lieleg 2009), eta = 1 (rigor HMM), dxc = 0.1 nm fixed;
#      free: lambda0, s_taut, koff0, phi, dxs_nm. Start P1 = s1 optimum with these three values.
#
#   B1 (main)  7 sweep amplitudes (to 0.158), weight 1.5     B2  7 amplitudes, weight 1
#   B3         8 amplitudes (+ 0.2), weight 1.5              B4  8 amplitudes, weight 1
#   C1 (control) B1 settings, all 8 parameters free from the s2 optimum (cost of the pruning)
#   C2 (control) B1 settings, all 8 parameters free from the s4 optimum (does the second basin survive?)
#   C3 (control) B1 with dxc = 1.5 nm fixed (Guo & Guilford 2006 catch arm): cost of imposing the single-molecule
#                catch-slip law (fitted laws are nearly force-independent below ~5 pN; LAOS forces mostly < 2 pN)
#
# Objective J = loop_RMS^2 + w^2 * sweep_RMS_log10^2 (w = sweep weight); the least-squares cost is N_data * J / 2.
# Engine: closed-system Fortran driver (10_CL_DIFF_CS_SFR_FI_TC/main_closed.f90), ~26 s per loop + 7-amplitude sweep
# on one core (Python: ~45 min). Production settings: dt 0.1, factor 3, 5 cycles, max-nfev 40 (LSQ tol 1e-5).
# Needs on the VM: optimise_laos_v2.py, amplitude_sweep_runs.py, material_TC.py, dfd.py (current versions)
# and the folder 10_CL_DIFF_CS_SFR_FI_TC next to python_scripts/. The driver is rebuilt below if missing or
# older than its sources (optimise_laos_v2 stops otherwise). Seven fits x 3 workers = 21 cores (adapt WORKERS).
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
WORKERS=${WORKERS:-3}
FDIR=../10_CL_DIFF_CS_SFR_FI_TC

# --- build the closed-system driver (laos_closed) if needed ---
if [ ! -x $FDIR/laos_closed ] || [ -n "$(find $FDIR/global.f90 $FDIR/main_closed.f90 $FDIR/src/material_AFFCL -name '*.f90' -newer $FDIR/laos_closed)" ]; then
  ( cd $FDIR && cat global.f90 > umat.f90 && \
    find . -type f -path '*src/material_AFFCL*' -name '*.f90' -exec cat {} >> umat.f90 \; && \
    gfortran -O2 -o laos_closed umat.f90 main_closed.f90 )
  echo "built $FDIR/laos_closed"
fi

BASE="--consistent-kinetics --engine fortran --dt 0.1 --factor 3 --cycles 5 --max-nfev 40"
G7="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158"
G8="--sweep-gammas 0.0126 0.025 0.05 0.079 0.1 0.126 0.158 0.2"

# pruned start (s1 optimum; keq, eta, dxc fixed)
P1="--x0 lambda0=1.0514636806123838 --x0 s_taut=1.0935008421037167 --x0 koff0=0.05630012135707198 \
 --x0 phi=0.9268128731506545 --x0 dxs_nm=1.292674573055505 --x0 keq=14.0 --x0 eta=1.0 --fix dxc_nm=0.1"
# full starts (Stage A optima, all 8 free)
S2="--free keq --free eta --x0 lambda0=1.045152730091424 --x0 s_taut=1.0736940328978668 --x0 koff0=0.05775814120056595 \
 --x0 phi=0.9895433042217543 --x0 dxc_nm=0.10177148593294981 --x0 dxs_nm=1.9571443271556497 \
 --x0 keq=28.245092386699753 --x0 eta=0.8537799699316309"
S4="--free keq --free eta --x0 lambda0=1.046957108365689 --x0 s_taut=1.0904593895004608 --x0 koff0=0.050958735276755 \
 --x0 phi=0.6584549951413714 --x0 dxc_nm=0.43179870490835515 --x0 dxs_nm=0.7915242641584863 \
 --x0 keq=45.31323760328048 --x0 eta=0.9900828401040771"

launch() {  # name, options
  nohup python -B optimise_laos_v2.py $BASE $2 --workers $WORKERS --out optimisation_runs/final_$1 > final_$1.log 2>&1 &
  echo "started final_$1 (pid $!)"
}
launch B1 "$G7 --sweep-weight 1.5 $P1"
launch B2 "$G7 --sweep-weight 1.0 $P1"
launch B3 "$G8 --sweep-weight 1.5 $P1"
launch B4 "$G8 --sweep-weight 1.0 $P1"
launch C1 "$G7 --sweep-weight 1.5 $S2"
launch C2 "$G7 --sweep-weight 1.5 $S4"
launch C3 "$G7 --sweep-weight 1.5 ${P1/dxc_nm=0.1/dxc_nm=1.5}"
