subroutine rnd_gennor (mu, sd, phrase, n, array)

  ! Based on subroutine test_gennor 
  ! See test_gennor (or other test subroutines in main.f90) to verify statistics of the generated distribution
  
  implicit none
    
  integer ( kind = 4 ) n
  integer ( kind = 4 ) i
  integer ( kind = 4 ) seed1
  integer ( kind = 4 ) seed2
  real ( kind = 4 ) gennor
  real ( kind = 4 ), intent(out) :: array(n)
  real ( kind = 4 ), intent(in)  :: mu, sd
  character ( len = * ) phrase

    
  !  Initialize the generators.
  call initialize_gen ( )   
  
  !  Set the seeds based on the phrase.
  call phrtsd ( phrase, seed1, seed2 )
  
  !  Initialize all generators.
  call set_initial_seed ( seed1, seed2 )

  !  Generate N samples.
  do i = 1, n
    array(i) = gennor ( mu, sd )
    !write(*,*) array(i)
  end do

  return
end
