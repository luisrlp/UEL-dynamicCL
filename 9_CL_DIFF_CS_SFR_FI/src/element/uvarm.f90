      subroutine UVARM(UVAR, DIRECT, T, TIME, DTIME, CMNAME, ORNAME, &
                               NUVARM, NOEL, NPT, LAYER, KSPT, KSTEP, KINC, NDI, NSHR, COORD, &
                               JMAC, JMATYP, MATLAYO, LACCFLA)

        ! This subroutine is used to transfer SDV's from the UEL
        ! onto the dummy mesh for viewing. Note that an offset of
        ! ElemOffset is used between the real mesh and the dummy mesh.
        ! If your model has more than ElemOffset UEL elements, then
        ! this will need to be modified.

        use global

        include 'aba_param.inc'

        character(len=80) :: CMNAME, ORNAME
        character(len=3) :: FLGRAY(15)
        DOUBLE PRECISION :: UVAR(NUVARM), DIRECT(3,3), T(3,3), TIME(2)
        integer :: ARRAY(15), JARRAY(15), JMAC(*), JMATYP(*), COORD(*)
        integer :: i1

        ! The dimensions of the variables FLGRAY, ARRAY and JARRAY
        ! must be set equal to or greater than 15.

        ! uvar(1) = globalSdv(noel-ElemOffset,npt,1)
        ! for example
        ! uvar(2) = globalSdv(noel-ElemOffset,npt,2)

        ! Bounds of globalSdv (allocated by the UEL; see the check in U3D8)
        if (.not. allocated(globalSdv)) then
            UVAR = 0.0d0
            return
        end if
        if ((noel-ElemOffset < 1) .or. (noel-ElemOffset > size(globalSdv,1)) .or. &
            (npt < 1) .or. (npt > size(globalSdv,2)) .or. (nuvarm < nsdv)) then
            write(*,*) 'ERROR in UVARM: globalSdv index out of bounds'
            write(*,*) '  noel =', noel, '  ElemOffset =', ElemOffset, '  numElem =', size(globalSdv,1)
            write(*,*) '  npt =', npt, '  nuvarm =', nuvarm, '  nsdv =', nsdv
            call exit
        end if

        do i1 = 1, nsdv
            UVAR(i1) = globalSdv(noel-ElemOffset, npt, i1)
        end do
        ! uvar(3) = globalSdv(noel-ElemOffset,npt,3)
        ! uvar(4) = globalSdv(noel-ElemOffset,npt,4)

        return
      end subroutine UVARM
      
