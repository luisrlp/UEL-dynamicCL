subroutine AssembleElement(nDim, nNode, ndofel, &
                           Ru, Rc, &
                           Kuu, Kuc, Kcu, Kcc, & 
                           rhs, amatrx)

   ! Subroutine to assemble the local elements residual and tangent

   implicit none

   integer :: i, j, nDim, nNode, ndofel, nDofN
   integer :: A11, A12, B11, B12

   real(8), intent(in) :: Ru(nDim*nNode, 1), Rc(nNode, 1)
   real(8), intent(in) :: Kuu(nDim*nNode, nDim*nNode)
   real(8), intent(in) :: Kuc(nDim*nNode, nNode), Kcu(nNode, nDim*nNode)
   real(8), intent(in) :: Kcc(nNode, nNode)
   real(8), intent(out) :: rhs(ndofel, 1), amatrx(ndofel, ndofel)

   ! Total number of degrees of freedom per node
   nDofN = ndofel / nNode

   ! Initialize
   rhs = 0.d0
   amatrx = 0.d0

   !!!!!! 2D !!!!!!
   if (nDim == 2) then
      !
      ! Assemble the element level residual
      !
      do i = 1, nNode
         A11 = nDofN * (i - 1) + 1
         A12 = nDim * (i - 1) + 1
         !
         ! Displacement
         !
         rhs(A11, 1) = Ru(A12, 1)
         rhs(A11 + 1, 1) = Ru(A12 + 1, 1)
         !
         ! Chemical Potential
         !
         rhs(A11 + 2, 1) = Rc(i, 1)
      end do
      !
      ! Assemble the element level tangent matrix
      !
      do i = 1, nNode
         do j = 1, nNode
            A11 = nDofN * (i - 1) + 1
            A12 = nDim * (i - 1) + 1
            B11 = nDofN * (j - 1) + 1
            B12 = nDim * (j - 1) + 1
            !
            ! Displacement
            !
            amatrx(A11, B11) = Kuu(A12, B12)
            amatrx(A11, B11 + 1) = Kuu(A12, B12 + 1)
            amatrx(A11 + 1, B11) = Kuu(A12 + 1, B12)
            amatrx(A11 + 1, B11 + 1) = Kuu(A12 + 1, B12 + 1)
            !
            ! Chemical Potential
            !
            amatrx(A11 + 2, B11 + 2) = Kcc(i,j)
            !
            ! Displacement - Chemical Potential
            !
            amatrx(A11, B11 + 2) = Kuc(A12, j)
            amatrx(A11 + 1, B11 + 2) = Kuc(A12 + 1, j)
            !
            ! Chemical Potential - Displacement
            !
            amatrx(A11 + 2, B11) = Kcu(i, B12)
            amatrx(A11 + 2, B11 + 1) = Kcu(i, B12 + 1)
         end do
      end do

   !!!!!! 3D !!!!!!
   elseif (nDim == 3) then
      !
      ! Assemble the element level residual
      !
      do i = 1, nNode
         A11 = nDofN * (i - 1) + 1
         A12 = nDim * (i - 1) + 1
         !
         ! Displacement
         !
         rhs(A11, 1) = Ru(A12, 1)
         rhs(A11 + 1, 1) = Ru(A12 + 1, 1)
         rhs(A11 + 2, 1) = Ru(A12 + 2, 1)
         !
         ! Chemical Potential
         !
         rhs(A11 + 3, 1) = Rc(i, 1)
      end do
      !
      ! Assemble the element level tangent matrix
      !
      do i = 1, nNode
         do j = 1, nNode
            A11 = nDofN * (i - 1) + 1
            A12 = nDim * (i - 1) + 1
            B11 = nDofN * (j - 1) + 1
            B12 = nDim * (j - 1) + 1
            !
            ! Displacement
            !
            amatrx(A11, B11) = Kuu(A12, B12)
            amatrx(A11, B11 + 1) = Kuu(A12, B12 + 1)
            amatrx(A11, B11 + 2) = Kuu(A12, B12 + 2)
            amatrx(A11 + 1, B11) = Kuu(A12 + 1, B12)
            amatrx(A11 + 1, B11 + 1) = Kuu(A12 + 1, B12 + 1)
            amatrx(A11 + 1, B11 + 2) = Kuu(A12 + 1, B12 + 2)
            amatrx(A11 + 2, B11) = Kuu(A12 + 2, B12)
            amatrx(A11 + 2, B11 + 1) = Kuu(A12 + 2, B12 + 1)
            amatrx(A11 + 2, B11 + 2) = Kuu(A12 + 2, B12 + 2)
            !
            ! Chemical Potential
            !
            amatrx(A11 + 3, B11 + 3) = Kcc(i,j)
            !
            ! Displacement - Chemical Potential
            !
            amatrx(A11, B11 + 3) = Kuc(A12, j)
            amatrx(A11 + 1, B11 + 3) = Kuc(A12 + 1, j)
            amatrx(A11 + 2, B11 + 3) = Kuc(A12 + 2, j)
            !
            ! Chemical Potential - Displacement
            !
            amatrx(A11 + 3, B11) = Kcu(i, B12)
            amatrx(A11 + 3, B11 + 1) = Kcu(i, B12 + 1)
            amatrx(A11 + 3, B11 + 2) = Kcu(i, B12 + 2)
         end do
      end do

   else
      write(*, *) 'How did you get nDim=', nDim
      call exit
   end if

   return
end subroutine AssembleElement

