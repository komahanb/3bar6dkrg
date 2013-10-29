  subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct)
    use dimkrig,only:fctindx
    implicit none

    integer  :: ndvart,fct
    double precision :: X(ndvart),fobj,dfdD(ndvart),x3
    double precision ::  rho, L, sigmay, pi, p, E, Fs  
    
    fctindx=fct

    call calcf(x,ndvart,12,fobj)

    call calcdf(x,ndvart,12,dfdD)   

    return
  end subroutine CalcstuffBFGS
