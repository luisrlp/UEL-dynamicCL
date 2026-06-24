subroutine thetafFunc(thetaf, f, df, args, nargs)
    ! This subroutine serves as the function we would like to solve for                                         
    ! the free crosslinker volume fraction (thetaf = cf/cfmax)                                                  
    ! by finding thetaf such that f = 0                                                                         
                                                                                                                
    implicit none                                                                                               
                                                                                                                
    integer, intent(in)              :: nargs                                                                   
    DOUBLE PRECISION, intent(in out) :: thetaf                                                                   
    DOUBLE PRECISION, intent(out)    :: f, df                                                                    
    DOUBLE PRECISION, intent(in)     :: args(nargs)                                                              
                                                                                                                
    DOUBLE PRECISION                 :: mu, mu0, Rgas, theta, chi, Vmol, Kbulk                                                           
    DOUBLE PRECISION                 :: detF, RT, Jc, Je, cb, cfmax
    DOUBLE PRECISION, parameter      :: zero  = 0.0d0                                                                         
    DOUBLE PRECISION, parameter      :: one   = 1.0d0                                                                         
    DOUBLE PRECISION, parameter      :: two   = 2.0d0                                                                         
                                                                                                                
    ! Obtain relevant quantities                                                                                
    mu    = args(1)                                                                                             
    mu0   = args(2)                                                                                             
    Rgas  = args(3)                                                                                             
    theta = args(4)                                                                                             
    chi   = args(5)                                                                                             
    Vmol  = args(6)                                                                                             
    Kbulk = args(7)                                                                                             
    detF  = args(8)                                                                                             
    cb    = args(9)                                                         
    cfmax = args(10)                                                                    
                                                                                                                
    ! Compute the useful quantity                                                                               
    RT = Rgas * theta                                                                                           
                                                                                                                
    ! Compute the swelling ratio J^c
    Jc = one + Vmol * cb + Vmol * cfmax * thetaf
    
    ! Compute Elastic Volume Ratio J^e                                                                          
    Je = detF / Jc                                                       
                                                                                                                
    ! Compute the residual f(thetaf) = 0                                                                        
    f = (mu0 - mu) / RT &                                                                                       
        + log(thetaf / (one - thetaf)) &                                                                          
        + chi * (one - two * thetaf) &                                                                            
        - ((Kbulk * Vmol) / RT) * (log(Je) / Jc)
                                                                                                                
    ! Compute the exact analytical tangent df/dthetaf                                                           
    df = (one / thetaf) + (one / (one - thetaf)) &                                                              
        - two * chi &                                                                                            
        + ((Kbulk * Vmol) / RT) * (Vmol * cfmax) * (one + log(Je)) / (Jc * Jc)  

end subroutine thetafFunc

