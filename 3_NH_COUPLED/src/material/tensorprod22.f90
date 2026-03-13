subroutine tensorprod2(a, b, c, ndi)
       use global
       implicit none
       integer, intent(in) :: ndi
       real(8), intent(in) :: a(ndi, ndi), b(ndi, ndi)
       real(8), intent(out) :: c(ndi, ndi, ndi, ndi)
       integer :: i, j, k, l

       do i = 1, ndi
              do j = 1, ndi
                     do k = 1, ndi
                            do l = 1, ndi
                                   c(i, j, k, l) = a(i, j) * b(k, l)
                            end do
                     end do
              end do
       end do

end subroutine tensorprod2
