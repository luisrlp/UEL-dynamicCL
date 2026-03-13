subroutine resetdfgrd(dfgrd, ndi)
      use global
      implicit none

      integer, intent(in) :: ndi
      real(8), intent(out) :: dfgrd(ndi, ndi)

      dfgrd = 0.0
      dfgrd(1,1) = 1.0
      dfgrd(2,2) = 1.0
      dfgrd(3,3) = 1.0

end subroutine resetdfgrd
