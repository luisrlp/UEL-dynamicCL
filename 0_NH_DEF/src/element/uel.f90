!************************************************************************
!
! User element for transient fluid permeation, and large 
!  elastic deformation in 2D or 3D.  This is for plane strain,
!  axisymetric, and 3D.
!
! Solution variables (or nodal variables) are the displacements and the
!  chemical potential.
! 
! This subroutine is for the following element types
!  > two-dimensional 4 node isoparametric element as shown below
!       with 1pt (reduced) or 4pt (full) gauss integration.
!  > three-dimensional 8 node isoparametric element as shown below
!       with 1pt (reduced) or 8pt (full) gauss integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
!  Mechanical, traction- and pressure-type boundary conditions 
!   may be applied to the dummy mesh using the Abaqus built-in 
!   commands *Dload or *Dsload.
!
! Surface flux boundary conditions are supported in the following
!  elements.  Based on our convention, the face on which the fliud
!  flux is applied is the "label", i.e.
!  - U1,U2,U3,U4,... refer to fluid fluxes applied to faces 
!                     1,2,3,4,... respectively,
!
!     
!              A eta (=xi_2)
!  4-node      |
!   quad       |Face 3
!        4-----------3
!        |     |     |
!        |     |     |
!  Face 4|     ------|---> xi (=xi_1)
!        |           | Face2
!        |           |
!        1-----------2
!          Face 1
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!     Face numbering follows:
!       Face 1 = nodes 1,2,3,4
!       Face 2 = nodes 5,8,7,6
!       Face 3 = nodes 1,5,6,2
!       Face 4 = nodes 2,6,7,3
!       Face 5 = nodes 3,7,8,4
!       Face 6 = nodes 4,8,5,1
!
! Shawn A. Chester, December 2010 -- as used in my prior publications
! Shawn A. Chester, December 2013 -- modified for public distribution
!
!***********************************************************************
!
! User element statement in the input file (set ? values as needed):
!
!  2D elements
!  *User Element,Nodes=4,Type=U?,Iproperties=2,Properties=9,Coordinates=2,Variables=?,Unsymm
!  1,2,11
!
!  3D elements
!  *User Element,Nodes=8,Type=U3,Iproperties=2,Properties=9,Coordinates=3,Variables=?,Unsymm
!  1,2,3,11
!
!
!     State Variables
!     --------------------------------------------------------------
!     Global SDV's (used for visualization)
!       1) polymer volume fraction (phi)
!
!     Local SDV's (used for the solution procedure)
!       j = 0
!       do k = 1,nInttPt
!          svars(1+j) = phi ---- polymer volume fraction at integ pt k
!          j = j + nlSdv
!       end loop over k
!
!     In the input file, set 'User output variables'= number of global SDV's
!
!     In the input file, set 'ngSdv'= number of global SDV's
!
!     In the input file, set 'nlSdv'= number of local SDV's
!
!     In the input file, set 'varibles'=(nlSdv*nInttPt)
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!     Kbulk  = props(1) ! Bulk modulus
!     C10    = props(2) ! Chi parameter
!     C01     = props(3) ! Coefficient of permeability
!     nlSdv  = jprops(1) ! Number of local sdv's per integ pt
!     ngSdv  = jprops(2) ! Number of global sdv's per integ pt
!
!***********************************************************************

subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
            PROPS, NPROPS, coords, MCRD, NNODE, Uall, DUall, Vel, &
            Accn, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, &
            NDLOAD, JDLTYP, ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, &
            DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD)
  use global
  implicit none

  ! Variables defined in UEL, passed back to Abaqus
  real(8), intent(out) :: RHS(MLVARX, 1), AMATRX(NDOFEL, NDOFEL), &
                    SVARS(NSVARS), ENERGY(8)
  
  ! Variables passed into UEL
  real(8), intent(in) :: PROPS(NPROPS), coords(MCRD, NNODE), &
                   Uall(NDOFEL), DUall(MLVARX, 1), Vel(NDOFEL), &
                   Accn(NDOFEL), TIME(2), DTIME, PARAMS(1), &
                   ADLMAG(MDLOAD, 1), PREDEF(2, NPREDF, NNODE), &
                   DDLMAG(MDLOAD, 1), PNEWDT, PERIOD
  integer, intent(in) :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, &
                   JTYPE, KSTEP, KINC, JELEM, NDLOAD, JDLTYP(MDLOAD, 1), &
                   NPREDF, LFLAGS(4), MLVARX, MDLOAD, JPROPS(NJPROP), &
                   NJPROP

  integer :: lenJobName, lenOutDir, nDim, nIntt, nInttS
  character(len=256) :: jobName, outDir, fileName

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nIntt = 8  ! number of volume integration points
  nInttS = 1 ! number of surface integration points
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !----------------------------------------------------------------
  ! Perform initial checks
  ! Open the debug/error message file
  ! call getJobName(jobName, lenJobName)
  ! call getOutDir(outDir, lenOutDir)
  ! fileName = trim(outDir(1:lenOutDir)) // '\aaMSGS_' // &
  !            trim(jobName(1:lenJobName)) // '.dat'
  ! open(unit=80, file=fileName, status='unknown')

  ! Check the procedure type, this should be a coupled displacement
  ! which are any of the following (64, 65, 72, 73)
  if ((lflags(1) == 1) .or. (lflags(1) == 65) .or. & !64,65
     (lflags(1) == 72) .or. (lflags(1) == 73)) then
   ! all is good
  else
   write(*,*) 'Abaqus does not have the right procedure'
   write(*,*) 'go back and check the procedure type'
   write(*,*) 'lflags(1)=', lflags(1)
   ! write(80,*) 'Abaqus does not have the right procedure'
   ! write(80,*) 'go back and check the procedure type'
   ! write(80,*) 'lflags(1)=', lflags(1)
   call exit
  endif

  ! Make sure Abaqus knows you are doing a large deformation problem
  ! I think this only matters when it comes to output in viewer
  if (lflags(2) == 0) then
   ! lflags(2)=0 -> small disp.
   ! lflags(2)=1 -> large disp.
   write(*,*) 'Abaqus thinks you are doing'
   write(*,*) 'a small displacement analysis'
   write(*,*) 'go in and set nlgeom=yes'
   ! write(80,*) 'Abaqus thinks you are doing'
   ! write(80,*) 'a small displacement analysis'
   ! write(80,*) 'go in and set nlgeom=yes'
   call exit
  endif

  ! Check to see if you are doing a general step or a linear perturbation step
  if (lflags(4) == 1) then
   ! lflags(4)=0 -> general step
   ! lflags(4)=1 -> linear perturbation step
   write(*,*) 'Abaqus thinks you are doing'
   write(*,*) 'a linear perturbation step'
   ! write(80,*) 'Abaqus thinks you are doing'
   ! write(80,*) 'a linear perturbation step'
   call exit
  endif

  ! Do nothing if a ``dummy'' step
  if (dtime == 0.0) return

  ! Done with initial checks
  !----------------------------------------------------------------
  if (jtype == 3) then
   ! This is a 3D analysis
   nDim = 3
   call U3D8(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
           PROPS, NPROPS, coords, MCRD, NNODE, Uall, DUall, Vel, &
           Accn, JTYPE, TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, &
           NDLOAD, JDLTYP, ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, &
           DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD, &
           nDim, nIntt, nInttS)
  else
   ! We have a problem...
   write(*,*) 'Element type not supported, jtype=', jtype
   ! write(80,*) 'Element type not supported, jtype=', jtype
   call exit
  endif

  ! Done with this element, RHS and AMATRX already returned
  ! as output from the specific element routine called
  !----------------------------------------------------------------
  return
end subroutine UEL






