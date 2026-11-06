#!/usr/bin/env bash
# Build the closed-system LAOS driver (10_CL_DIFF_CS_SFR_FI_TC/main_closed.f90 + the C4 UMAT routines) once per
# integration factor. FACTOR is a compile-time constant (global.f90, NDIR = 10 FACTOR^2), so each factor gets its own
# executable, 10_CL_DIFF_CS_SFR_FI_TC/laos_closed_f<factor>, compiled from a temporary copy of global.f90 with that
# FACTOR. The repository files (global.f90, umat.f90) are not modified, so the Abaqus UMAT/UEL keep their own FACTOR.
#
#   bash build_laos_closed.sh          # factors 3 6 9
#   bash build_laos_closed.sh 3 6      # selected factors
#
# Rebuild after any change of global.f90, main_closed.f90 or src/material_AFFCL/*.f90 (optimise_laos_v2.py stops
# if an executable is older than its sources).
set -e
FDIR=$(cd "$(dirname "$0")/../10_CL_DIFF_CS_SFR_FI_TC" && pwd)
FACTORS=${*:-3 6 9}
for k in $FACTORS; do
  TMP=$(mktemp -d)
  sed -E "s/PARAMETER \(FACTOR=[0-9]+\)/PARAMETER (FACTOR=$k)/" "$FDIR/global.f90" > "$TMP/global.f90"
  grep -q "PARAMETER (FACTOR=$k)" "$TMP/global.f90" || { echo "could not set FACTOR=$k in global.f90"; exit 1; }
  { cat "$TMP/global.f90"; find "$FDIR/src/material_AFFCL" -type f -name '*.f90' -exec cat {} \; ; } > "$TMP/umat.f90"
  ( cd "$TMP" && gfortran -O2 -I"$FDIR" -J"$TMP" -o "$FDIR/laos_closed_f$k" umat.f90 "$FDIR/main_closed.f90" )
  rm -rf "$TMP"
  echo "built $FDIR/laos_closed_f$k (NDIR = $((10 * k * k)))"
done
