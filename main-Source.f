      program main
      
      implicit none


      INTEGER, PARAMETER :: Nx=500
      INTEGER, PARAMETER :: Ny=80

      double precision dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      common /const/ dL1,dL2,dk,dc,dalpha,depsilon,depsilonp,
     .               dlambda1,dlambda2,s1,s2,vd,tend,tout,dt,tE,
     .               dx,dy,tol,isf,itstart,pi,amplit,prob

      double precision gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0

      common /param/ gamma01,beta01,ro01,Diffgamma,dke0,dk1,dsigma0
      double precision beta(Nx,Ny),gamma(Nx,Ny),ro(Nx,Ny)
      double precision dke(Nx,Ny),dsigma(Nx,Ny),TS(Nx,Ny)
      double precision t, iPeriod, it, it1
      double precision romin, romax
      integer counter
      character(len=4) ct1
      character(len=17) ct2
      character(len=21) ct3

      double precision gamma0(10),ro0(10),beta0(10)
      integer ier,pgbeg, i, j
      t=0.d0
      counter=0;
      iPeriod=0
      call anfang(t,Nx,Ny,beta,gamma,ro,TS)
          open(10,file ='OutputData2D/data   0'
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,gamma,ro,beta,TS)
          close(10)
      ct2='OutputData2D/data'

 5    continue


      call ODE(t,Nx,Ny,beta,gamma,ro,TS,dke,dsigma)

      write(6,*) 'real t= '
      write(6,'(F6.2)') t/dk1
      counter=counter+1
      write(ct1,'(I4)') counter
      ct3 = ct2 // ct1

      write(6,*) ct3

          open(10,file =ct3
     .     ,status = 'unknown',form = 'formatted')
              call out(t,Nx,Ny,gamma,ro,beta,TS)
          close(10)


      if (t+dt .lt. tend) then
         goto 5

      endif

      close(10)

!

!     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     WRITES FINAL STATE
      open(42,file ='OutputData2D/Final-State'
     . ,status = 'unknown',form = 'formatted')
      call outFinal(t,Nx,Ny,gamma,ro,beta,TS)
      close(42)

      end
      


