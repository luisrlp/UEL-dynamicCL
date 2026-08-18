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
        do i1 = 1, nsdv
            UVAR(i1) = globalSdv(noel-ElemOffset, npt, i1)
        end do
        ! uvar(3) = globalSdv(noel-ElemOffset,npt,3)
        ! uvar(4) = globalSdv(noel-ElemOffset,npt,4)

        return
      end subroutine UVARM
