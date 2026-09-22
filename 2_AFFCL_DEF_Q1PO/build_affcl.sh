cat global.f90 > umat_affcl.f90
find . -type f -path '*src/material_AFFCL/*' -name '*.f90' -exec cat {} +> umat_affcl.f90 
gfortran -o affcl.o umat_affcl.f90 main_umat.f90
