cat global.f90 > umat.f90
find . -type f -path '*src/material_AFFCL*' -name '*.f90' -exec cat {} >> umat.f90 \;
gfortran -o test_diffcl umat.f90 main_umat.f90
