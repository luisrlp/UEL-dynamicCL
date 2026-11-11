#!/usr/bin/env bash
# Build the sensitivity-analysis material-point driver (10_CL_DIFF_CS_SFR_FI_TC/main_sa.f90) against the UMAT
# routines of one model folder, once per integration factor:
#   10 = 10_CL_DIFF_CS_SFR_FI_TC (consistent kinetics, TC)
#    9 = 9_CL_DIFF_CS_SFR_FI     (non-consistent kinetics); xmech.f90 of folder 10 is added (outputs only)
# FACTOR is a compile-time constant (global.f90, NDIR = 10 FACTOR^2): each factor gets its own executable,
# <folder>/sa_point_f<factor>, compiled from a temporary copy of global.f90 with that FACTOR. The repository files
# (global.f90, umat.f90) are not modified.
#
#   bash build_sa.sh 10          # TC, factor 6
#   bash build_sa.sh 9 3 6       # non-TC, factors 3 and 6
#
# Rebuild after any change of global.f90, main_sa.f90 or src/material_AFFCL/*.f90.
set -e
REPO=$(cd "$(dirname "$0")/.." && pwd)
case "$1" in
  10) FDIR="$REPO/10_CL_DIFF_CS_SFR_FI_TC" ;;
  9)  FDIR="$REPO/9_CL_DIFF_CS_SFR_FI" ;;
  *)  echo "usage: bash build_sa.sh <9|10> [factors]"; exit 1 ;;
esac
shift
FACTORS=${*:-6}
DRIVER="$REPO/10_CL_DIFF_CS_SFR_FI_TC/main_sa.f90"
XMECH="$REPO/10_CL_DIFF_CS_SFR_FI_TC/src/material_AFFCL/xmech.f90"
for k in $FACTORS; do
  TMP=$(mktemp -d)
  sed -E "s/PARAMETER \(FACTOR=[0-9]+\)/PARAMETER (FACTOR=$k)/" "$FDIR/global.f90" > "$TMP/global.f90"
  grep -q "PARAMETER (FACTOR=$k)" "$TMP/global.f90" || { echo "could not set FACTOR=$k in global.f90"; exit 1; }
  { cat "$TMP/global.f90"; find "$FDIR/src/material_AFFCL" -type f -name '*.f90' -exec cat {} \; ;
    [ -f "$FDIR/src/material_AFFCL/xmech.f90" ] || cat "$XMECH"; } > "$TMP/umat.f90"
  ( cd "$TMP" && gfortran -O2 -I"$FDIR" -J"$TMP" -o "$FDIR/sa_point_f$k" umat.f90 "$DRIVER" )
  rm -rf "$TMP"
  echo "built $FDIR/sa_point_f$k (NDIR = $((10 * k * k)))"
done
