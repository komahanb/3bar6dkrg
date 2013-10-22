program krig
  use dimKrig
      implicit none
      integer,parameter:: N=2
      integer::i
      double precision F, X(N),sigmax(N),fmeantmp,fvartmp,fmeanprimetmp(n),fvarprimetmp(n)
      double precision DAT(5)
      integer IDAT(5)
      integer kprob,NMC
      integer IERR
      double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
      double precision  rho, L, sigmay, pi, p, E, Fs 

      !if (id_proc.eq.0) print *,'Calculate Objective',X

      NMC=100000

      Fs=DAT(1)
      kprob=IDAT(1)
      do i=1,N
         sigmax(i)=DAT(i+1)
         Xsave(i)=X(i)
      end do     

!!$      gtol=1e-4
!!$
!!$      low(1:N-1)=X(1:N-1)-sigmax(1:N-1)
!!$      up(1:N-1)=X(1:N-1)+sigmax(1:N-1)
!!$
!!$      call optimize(N-1,X,N,fmax,gradmax,low,up,gtol,.true.,.false.,1)


!---- MEAN and VARIANCE OF worst OBJECTIVE FUNCTION

      call Krigingestimate(1,X,N,sigmax,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,NMC,0)
    end program krig
