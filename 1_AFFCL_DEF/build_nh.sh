find . -type f -path '*src/material/*' -name '*.f90' -exec cat {} +> umat_nh.f90 
gfortran -o nh.o umat_nh.f90 main_umat.f90
