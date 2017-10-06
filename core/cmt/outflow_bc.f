C> @file outflow_bc.f Dirichlet states for outflow boundary conditions
      subroutine outflow(nvar,f,e,facew,wbc) ! don't really need nvar anymore
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'

      integer  nvar,f,e
      real facew(nx1,nz1,2*ldim,nelt,nvar)
      real wbc(nx1,nz1,2*ldim,nelt,nvar)

      call outflow_rflu(nvar,f,e,facew,wbc)

      return
      end
