C> @file surface_fluxes.f Routines for surface terms on RHS.

C> \ingroup isurf
C> @{
C> Restrict and copy face data and compute inviscid numerical flux 
C> \f$\oint \mathbf{H}^{c\ast}\cdot\mathbf{n}dA\f$ on face points
      subroutine fluxes_full_field
!-----------------------------------------------------------------------
! JH060314 First, compute face fluxes now that we have the primitive variables
! JH091514 renamed from "surface_fluxes_inviscid" since it handles all fluxes
!          that we compute from variables stored for the whole field (as
!          opposed to one element at a time).
!-----------------------------------------------------------------------
      include 'SIZE'
      include 'DG'
      include 'SOLN'
      include 'CMTDATA'
      include 'INPUT'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*3*lfq) ! might not need ldim
! JH070214 OK getting different answers whether or not the variables are
!          declared locally or in common blocks. switching to a different
!          method of memory management that is more transparent to me.
      common /CMTSURFLX/ fatface(heresize),notyet(hdsize)
      real fatface,notyet
      integer eq
      character*32 cname
      nfq=nx1*nz1*2*ndim*nelt
      nstate = nqq
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq

      i_cvars=(iu1-1)*nfq+1
      do eq=1,toteq
         call faceu(eq,fatface(i_cvars))
! JH080317 at least get the product rule right until we figure out how
!          we want the governing equations to look
         call invcol2(fatface(i_cvars),fatface(iwm+nfq*(iph-1)),nfq)
         i_cvars=i_cvars+nfq
      enddo

      call face_state_commo(fatface(iwm),fatface(iwp),nfq,nstate
     >                     ,dg_hndl)

      call InviscidBC(fatface(iwm),fatface(iwp),nstate)

      call InviscidFlux(fatface(iwm),fatface(iwp),fatface(iflx)
     >                 ,nstate,toteq)

!     call face_flux_commo(fatface(iflx),fatface(iflx),ndg_face,toteq,
!    >                     flux_hndl) ! for non-symmetric gs_op someday

C> @}

      return
      end

!-----------------------------------------------------------------------

      subroutine faceu(ivar,yourface)
! get faces of conserved variables stored contiguously
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'
      integer e
      real yourface(nx1,nz1,2*ldim,nelt)

      do e=1,nelt
         call full2face_cmt(1,nx1,ny1,nz1,iface_flux(1,e),
     >                      yourface(1,1,1,e),u(1,1,1,ivar,e))
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine fillq(ivar,field,wminus,yourface)
      include 'SIZE'
      include 'DG'

      integer ivar! intent(in)
      real field(nx1,ny1,nz1,nelt)! intent(in)
!     real, intent(out)wminus(7,nx1*nz1*2*ldim*nelt) ! gs_op no worky
      real wminus(nx1*nz1*2*ndim*nelt,*)! intent(out)
      real yourface(nx1,nz1,2*ndim,*)
      integer e,f

      nxz  =nx1*nz1
      nface=2*ndim

      call full2face_cmt(nelt,nx1,ny1,nz1,iface_flux,yourface,field)

      do i=1,ndg_face
         wminus(i,ivar)=yourface(i,1,1,1)
      enddo

      return
      end

!-----------------------------------------------------------------------

      subroutine face_state_commo(mine,yours,nf,nstate,handle)

! JH060414 if we ever want to be more intelligent about who gets what,
!          who gives what and who does what, this is the place where all
!          that is done. At the very least, gs_op may need the transpose
!          flag set to 1. Who knows. Everybody duplicates everything for
!          now.
! JH070714 figure out gs_op_fields, many, vec, whatever (and the
!          corresponding setup) to get this done for the transposed
!          ordering of state variables. I want variable innermost, not
!          grid point.

      integer handle,nf,nstate ! intent(in)
      real yours(*),mine(*)

      ntot=nf*nstate
      call copy(yours,mine,ntot)
!-----------------------------------------------------------------------
! operation flag is second-to-last arg, an integer
!                                                1 ==> +
      call gs_op_fields(handle,yours,nf,nstate,1,1,0)
      call sub2 (yours,mine,ntot)
      return
      end

!-----------------------------------------------------------------------

      subroutine face_flux_commo(flux1,flux2,nf,neq,handle)
! JH060514 asymmetric transposed gs_op, gs_unique magic may be needed if
!          we ever decide to avoid redundancy. For now, this routine
!          doesn't need to do anything.
      integer ntot,handle
      real flux1(*),flux2(*)
! JH061814 It doesn't need to do anything, but a sanity check would be
!          wise.
      return
      end

!-------------------------------------------------------------------------------

      subroutine surface_integral_full(vol,flux)
! Integrate surface fluxes for an entire field. Add contribution of flux
! to volume according to add_face2full_cmt
      include 'SIZE'
      include 'GEOM'
      include 'DG'
      include 'CMTDATA'
      real vol(nx1*ny1*nz1*nelt),flux(*)
      integer e,f

! weak form until we get the time loop rewritten
!     onem=-1.0
!     ntot=nx1*nz1*2*ndim*nelt
!     call cmult(flux,onem,ntot)
! weak form until we get the time loop rewritten
      call add_face2full_cmt(nelt,nx1,ny1,nz1,iface_flux,vol,flux)

      return
      end

!-------------------------------------------------------------------------------

      subroutine diffh2graduf(e,eq,graduf)
! peels off diffusiveH into contiguous face storage via restriction operator R
! for now, stores {{gradU}} for igu
      include  'SIZE'
      include  'DG' ! iface
      include  'CMTDATA'
      include  'GEOM'
      integer e,eq
      real graduf(nx1*nz1*2*ndim,nelt,toteq)
      common /scrns/ hface(lx1*lz1,2*ldim)
     >              ,normal(lx1*ly1,2*ldim)
      real hface, normal

      integer f

      nf    = nx1*nz1*2*ndim*nelt
      nfaces=2*ndim
      nxz   =nx1*nz1
      nxzf  =nxz*nfaces
      nxyz  = nx1*ny1*nz1

      call rzero(graduf(1,e,eq),nxzf) !   . dot nhat -> overwrites beginning of flxscr
      do j =1,ndim
         if (j .eq. 1) call copy(normal,unx(1,1,1,e),nxzf)
         if (j .eq. 2) call copy(normal,uny(1,1,1,e),nxzf)
         if (j .eq. 3) call copy(normal,unz(1,1,1,e),nxzf)
         call full2face_cmt(1,nx1,ny1,nz1,iface_flux,hface,diffh(1,j)) 
         call addcol3(graduf(1,e,eq),hface,normal,nxzf)
      enddo
      call col2(graduf(1,e,eq),area(1,1,1,e),nxzf)

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_cmt(flxscr,gdudxk,wminus)
! Hij^{d*}
      include 'SIZE'
      include 'CMTDATA'
      include 'DG'

      real flxscr(nx1*nz1*2*ndim*nelt,toteq)
      real gdudxk(nx1*nz1*2*ndim,nelt,toteq)
      real wminus(nx1*nz1,2*ndim,nelt,nqq)
      real const
      integer e,eq,f

      nxz = nx1*nz1
      nfaces=2*ndim
      nxzf=nxz*nfaces
      nfq =nx1*nz1*nfaces*nelt
      ntot=nfq*toteq

      call copy (flxscr,gdudxk,ntot) ! save AgradU.n
      const = 0.5
      call cmult(gdudxk,const,ntot)
!-----------------------------------------------------------------------
! supa huge gs_op to get {{AgradU}}
! operation flag is second-to-last arg, an integer
!                                                   1 ==> +
      call gs_op_fields(dg_hndl,gdudxk,nfq,toteq,1,1,0)
!-----------------------------------------------------------------------
      call sub2  (flxscr,gdudxk,ntot) ! overwrite flxscr with
                                      !           -
                                      ! (AgradU.n)  - {{AgradU.n}}
! [v] changes character on boundaries, so we actually need
! 1. (AgradU.n)- on Dirichlet boundaries
      call igu_dirichlet(flxscr,gdudxk)
! 2. (Fbc.n)- on Neumann boundaries
      call bcflux(flxscr,gdudxk,wminus)
      call chsign(flxscr,ntot) ! needs to change with sign changes

      return
      end

!-----------------------------------------------------------------------

      subroutine igu_dirichlet(flux,agradu)
! Acts on ALL boundary faces because I'm lazy. SET NEUMANN BC AFTER THIS
! CALL. BCFLUX IS PICKIER ABOUT THE BOUNDARY FACES IT ACTS ON.
      include 'SIZE'
      include 'TOTAL'
      integer e,eq,f
      real flux(nx1*nz1,2*ndim,nelt,toteq)
      real agradu(nx1*nz1,2*ndim,nelt,toteq)
      character*3 cb2

      nxz=nx1*nz1
      nfaces=2*ndim

      ifield=1
      do e=1,nelt
         do f=1,nfaces
            cb2=cbc(f, e, ifield)
            if (cb2.ne.'E  '.and.cb2.ne.'P  ') then ! cbc bndy.
! all Dirichlet conditions result in IGU being
! strictly one-sided, so we undo 0.5*QQT
! UNDER THE ASSUMPTIONS THAT
! 1. agradu's actual argument is really gdudxk AND
! 2. IT HAS ALREADY BEEN MULTIPLIED BY 0.5
! 3. gs_op has not changed it at all.
! overwriting flux with it and and multiplying it 2.0 should do the trick
               do eq=1,toteq
!                  call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
!! in fact, that copy should not be necessary at all. TEST WITHOUT IT
!                  call cmult(flux(1,f,e,eq),2.0,nxz)
! JH112216 This may be better because agradu (without the factor of 1/2) is
!          needed for some Neumann conditions (like adiabatic walls)
                   call cmult(agradu(1,f,e,eq),2.0,nxz)
                   call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
               enddo
            endif
         enddo
      enddo

      return
      end
