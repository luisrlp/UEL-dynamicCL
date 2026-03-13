
      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err
      INTEGER NWP,NELEM,NCH,NSDV,NTERM,FACTOR,NDIR,NGP
      DOUBLE PRECISION  ONE, TWO, THREE, FOUR, SIX, ZERO
      DOUBLE PRECISION HALF,THIRD
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=1)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=1000)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER (NTERM=60) ! 60
      PARAMETER (FACTOR=6) ! 6
      PARAMETER (NDIR=20 * FACTOR**2)
      PARAMETER (NGP=8)
      PARAMETER(NELEM=1, NSDV=7)
      PARAMETER(ZERO=0.D0, ONE=1.0D0,TWO=2.0D0)
      PARAMETER(THREE=3.0D0,FOUR=4.0D0,SIX=6.0D0)
      PARAMETER(HALF=0.5d0,THIRD=1.d0/3.d0)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CHARACTER(256) DIR2, DIR3
      PARAMETER (DIR2='prefdir.inp')
      PARAMETER (DIR3='etadir.inp')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8, allocatable :: globalSdv(:,:,:)

      end module global
