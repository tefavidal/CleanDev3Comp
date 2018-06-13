      subroutine ic(t,Nx,Ny,beta,gamma,ro)
      
      implicit none

      double precision t
      integer Nx, Ny

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0


      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision gamma0,beta0,ro0
      double precision dke(Nx,Ny),dsigma(Nx,Ny)
      integer nfix, i, j
      double precision factor


!     %%%%%Without Initial State
         nfix=1
         gamma0=0.006
         beta0=0.3167
         ro0=0.8953

!         gamma0(1)=0.5
!         beta0(1)=0.5
!         ro0(1)=0.5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do i=1,Nx
       do j=1,Ny

            ro(i,j)=ro0
            beta(i,j)=beta0
            gamma(i,j)=gamma0

       enddo
      enddo




         gamma01=gamma0
         beta01=beta0
         ro01=ro0


      return

      end

