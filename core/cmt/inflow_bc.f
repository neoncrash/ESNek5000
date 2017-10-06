C> @file Dirichlet states for inflow boundary conditions
      subroutine inflow(nvar,f,e,facew,wbc)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTBCDATA'
      integer nvar,f,e
      real facew(nx1,nz1,2*ldim,nelt,nvar)
      real wbc(nx1,nz1,2*ldim,nelt,nvar)

! JH021717 compare
!     call inflow_rflu(nvar,f,e,facew,wbc)
      call inflow_inviscid(nvar,f,e,facew,wbc)

      return
      end
