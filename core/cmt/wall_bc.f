      subroutine wallbc_inviscid(nstate,f,e,facew,wbc)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      include 'CMTBCDATA'

      integer nstate,f,e
      real    facew(nx1*nz1,2*ndim,nelt,nstate)
      real    wbc(nx1*nz1,2*ndim,nelt,nstate) 

! JH102016
! rind state for inviscid fluxes is different from viscous fluxes
      call reflect_rind(nstate,f,e,facew,wbc)

      return
      end
