C> @file diffusive_cmt.f routines for diffusive fluxes.
C> Some surface. Some volume. All pain. Jacobians and other factorizations.

C> \ingroup vsurf
C> @{
C> ummcu = \f$\mathbf{U}^--\{\{\mathbf{U}\}\}\f$
      subroutine imqqtu(ummcu,uminus,uplus)
! Computes (I-0.5*QQT)U for all five conserved variables.
! See call in compute_rhs_and_dt for important documentation
!                                     -
! Spoiler: for SEM this comes out to U -{{U}}, which when
!          spoken is "UMinus Minus the Central flux of U" which I
!          then abbreviate as ummcu
      include 'SIZE'

      real ummcu (nx1*nz1*2*ndim*nelt,toteq) ! intent(out)
      real uminus(nx1*nz1*2*ndim*nelt,toteq) ! intent(in)
      real uplus (nx1*nz1*2*ndim*nelt,toteq) ! intent(in)
      integer ivar

      nf = nx1*nz1*2*ndim*nelt
      const=-0.5

! U-{{U}} on interior faces. first just do it on all faces.
      do ivar=1,toteq
         call add3(ummcu(1,ivar),uminus(1,ivar),uplus(1,ivar),nf)
         call cmult(ummcu(1,ivar),const,nf)        !         -
         call add2(ummcu(1,ivar),uminus(1,ivar),nf)!ummcu = U -{{U}}
      enddo

C> @}
      return
      end

!-----------------------------------------------------------------------

C> \ingroup bcond
C> @{
C> umubc = \f$\mathbf{U}^--\mathbf{U}^D\f$
      subroutine imqqtu_dirichlet(umubc,wminus,wplus)
! v+ undefined on boundary faces, so (I-0.5QQ^T) degenerates to 
! [[U]] with respect to the Dirichlet boundary state
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'TSTEP' ! for ifield
      include 'CMTDATA'
      real umubc (nx1*nz1,2*ndim,nelt,toteq) ! intent(out)
      real wminus(nx1*nz1,2*ndim,nelt,nqq),
     >     wplus (nx1*nz1,2*ndim,nelt,nqq)
      real nTol
      integer e,f
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      nxz = nx1*nz1
      nface=2*ndim
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

      do e=1,nelt
      do f=1,nface

         cb=cbc(f,e,ifield)
         if (cb.ne.'E  '.and.cb.ne.'P  ') then ! cbc bndy. this routine
                                               ! had better not touch any
                                               ! interior face
! JH031315 flux added to argument list. BC routines preserve wminus for
!          obvious reasons and fill wplus with good stuff for everybody:
!          imposed states for Dirichlet conditions, and important things
!          for viscous numerical fluxes.
! JH060215 added SYM bc. Just use it as a slip wall hopefully.
! JH111416 This may look like lazy duplication, but there is a good chance
!          that this may replace/consolidate BC calls in InviscidFlus.
!          basically, nothing in wplus is trustworthy, so we are going to
!          recompute and overwrite iu1 through iu5 in wplus with stuff we
!          do trust (UBC, to be exact).
!           if (cb.eq.'v  ' .or. cb .eq. 'V  ') then
!             call inflow2(nqq,f,e,wminus,wplus)
            if (cb .eq. 'W  ' .or. cb .eq.'I  ')then
              call wallbc2(nqq,f,e,wminus,wplus)
            endif

!  -
! U - UBC
            do ivar=1,toteq
               call sub3(umubc(1,f,e,ivar),wminus(1,f,e,iu1+ivar-1),
     >                                      wplus(1,f,e,iu1+ivar-1),nxz)
            enddo

         endif 
      enddo
      enddo

C> @}
      return
      end

!-----------------------------------------------------------------------

      subroutine half_iku_cmt(res,diffh,e)
      include 'SIZE'
      include 'MASS'
! diffh has D AgradU. half_iku_cmt applies D^T BM1 to it and increments
! the residual res with the result
      integer e ! lopsided. routine for one element must reference bm1
      real res(nx1,ny1,nz1),diffh(nx1*ny1*nz1,ndim)

      n=nx1*ny1*nz1

      do j=1,ndim
         call col2(diffh(1,j),bm1(1,1,1,e),n)
!        call col2(diffh(1,j),phig(1,1,1,e),n) ! still no idea where phi goes
      enddo

!     const=-1.0 ! I0
      const=1.0  ! *-1 in time march
      call gradm11_t(res,diffh,const,e)

      return
      end
