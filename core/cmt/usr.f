      subroutine evaluate_dealiased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'
     
      integer  e,eq

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2

      n=nxd*nyd*nzd

      if (eq .eq. 1) then ! convective flux of mass=rho u_j=U_{j+1}

         do j=1,ndim
            call intp_rstd(convh(1,j),u(1,1,1,eq+j,e),nx1,nxd,if3d,0)
         enddo

      else

c To be consistent with momentum equation, for mass balance flux vector is 
c computed by multiplying rho by u_j
         call intp_rstd(ju1,phig(1,1,1,e),nx1,nxd,if3d,0)
         call intp_rstd(ju2,pr(1,1,1,e),nx1,nxd,if3d,0)

         if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
            call add2col2(convh(1,eq-1),ju1,ju2,n)

         elseif (eq .eq. 5) then

            call intp_rstd(convh(1,1),u(1,1,1,eq,e),nx1,nxd,if3d,0)
            call add2col2(convh(1,1),ju1,ju2,n)
            do j=2,ndim
               call copy(convh(1,j),convh(1,1),n)
            enddo
            call col2(convh(1,1),vxd(1,1,1,e),n)
            call col2(convh(1,2),vyd(1,1,1,e),n)
            call col2(convh(1,3),vzd(1,1,1,e),n)

         else
            if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
            if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
            call exitt
         endif

      endif
     
      return
      end

      subroutine evaluate_aliased_conv_h(e,eq)
! computed as products between primitive variables and conserved variables.
! if you want to write rho u_i u_j as (rho u_i) (rho u_j) (rho^{-1}), this
! is the place to do it
      include  'SIZE'
      include  'SOLN'
      include  'DEALIAS'
      include  'CMTDATA'
      include  'INPUT'

      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ju1(ldd),ju2(ldd)!,ur(ldd),us(ldd),ud(ldd),tu(ldd)
      real ju1,ju2
      integer  e,eq

      n=nxd*nyd*nzd

      call copy(ju1,phig(1,1,1,e),n)
      call copy(ju2,pr(1,1,1,e),n)

      if (eq .lt. 5) then ! self-advection of rho u_i by rho u_i u_j

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         if (if3d) call col2(convh(1,3),vzd(1,1,1,e),n)
         if(eq. gt. 1) call add2col2(convh(1,eq-1),ju1,ju2,n)

      elseif (eq .eq. 5) then

         call copy(convh(1,1),u(1,1,1,eq,e),n)
         call add2col2(convh(1,1),ju1,ju2,n)
         do j=2,ndim
            call copy(convh(1,j),convh(1,1),n)
         enddo
         call col2(convh(1,1),vxd(1,1,1,e),n)
         call col2(convh(1,2),vyd(1,1,1,e),n)
         call col2(convh(1,3),vzd(1,1,1,e),n)

      else
         if(nio.eq.0) write(6,*) 'eq=',eq,'really must be <= 5'
         if(nio.eq.0) write(6,*) 'aborting in evaluate_conv_h'
         call exitt
      endif

      return
      end
C> @}
!-----------------------------------------------------------------------

      subroutine cmt_flow_ics
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'

      integer e
      nxyz1 = nx1*ny1*nz1
      n     = nxyz1*lelt*toteq
      if (ifrestart)then
         do e=1,nelt
            call copy(U(1,1,1,2,e),vx(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,3,e),vy(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,4,e),vz(1,1,1,e),nxyz1) 
            call copy(U(1,1,1,5,e),t(1,1,1,e,1),nxyz1) 
            call copy(U(1,1,1,1,e),pr(1,1,1,e),nxyz1) 
         enddo
         call copy(tlag(1,1,1,1,1,2),t(1,1,1,1,2),nxyz1*nelt) ! s_{n-1}
         call copy(tlag(1,1,1,1,2,1),t(1,1,1,1,3),nxyz1*nelt) ! s_n
      endif
      call rzero(res1,n)
!     call copy(res2,t(1,1,1,1,5),n) ! art visc hardcoding. old entropy resid
      call rzero(res2,n) ! Actually,...
      return
      end

!-----------------------------------------------------------------------
      subroutine compute_forcing(e,eq_num)
      include  'SIZE'
      include  'INPUT'
      include  'GEOM'
      include  'MASS'
      include  'SOLN'
      include  'CMTDATA'
      include  'DEALIAS'
      
      integer e,eq_num
      parameter (ldd=lxd*lyd*lzd)
c     common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ur(lx1,ly1,lz1),us(lx1,ly1,lz1),ut(lx1,ly1,lz1),
     >    rdumz(lx1,ly1,lz1)

      nxyz=nx1*ny1*nz1
      if(eq_num.ne.1.and.eq_num.ne.5)then

        call gradl_rst(ur(1,1,1),us(1,1,1),ut(1,1,1),
     >                                        phig(1,1,1,e),lx1,if3d)
        if(if3d) then ! 3d
          if(eq_num.eq.2) then
            do i=1,nxyz
              rdumz(i,1,1) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e) +
     >              ut(i,1,1)*TXM1(i,1,1,e))
            enddo
          elseif(eq_num.eq.3) then
            do i=1,nxyz
              rdumz(i,1,1) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e) +
     >              ut(i,1,1)*TYM1(i,1,1,e))
            enddo
          elseif(eq_num.eq.4) then
            do i=1,nxyz
              rdumz(i,1,1) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RZM1(i,1,1,e) +
     >              us(i,1,1)*SZM1(i,1,1,e) +
     >              ut(i,1,1)*TZM1(i,1,1,e))
            enddo
          endif
        else ! end 3d, 2d
          if(eq_num.eq.2) then
            do i=1,nxyz
              rdumz(i,1,1) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RXM1(i,1,1,e) +
     >              us(i,1,1)*SXM1(i,1,1,e))
            enddo
          elseif(eq_num.eq.3) then
            do i=1,nxyz
              rdumz(i,1,1) = 1.0d+0/JACM1(i,1,1,e)*
     >             (ur(i,1,1)*RYM1(i,1,1,e) +
     >              us(i,1,1)*SYM1(i,1,1,e))
            enddo
          endif ! end 2d
      endif ! eqn nums 2-4

c     multiply by pressure
      call col2(rdumz,pr(1,1,1,e),nxyz)

c     gradphi.tauij !NOPE
!     do i=1,nxyz !   ????????????SiGN??????
!        rdumz(i,1,1)=rdumz(i,1,1)+diffh(i,1)*
!              (1.0d+0/JACM1(i,1,1,e)*
!    >             (ur(i,1,1)*RXM1(i,1,1,e) +
!    >              us(i,1,1)*SXM1(i,1,1,e) +
!    >              ut(i,1,1)*TXM1(i,1,1,e)))
!     enddo
!     do i=1,nxyz !   ????????????SiGN??????
!        rdumz(i,1,1)=rdumz(i,1,1)+diffh(i,2)*
!              (1.0d+0/JACM1(i,1,1,e)*
!    >             (ur(i,1,1)*RYM1(i,1,1,e) +
!    >              us(i,1,1)*SYM1(i,1,1,e) +
!    >              ut(i,1,1)*TYM1(i,1,1,e)))
!     enddo
!     do i=1,nxyz !   ????????????SiGN??????
!        rdumz(i,1,1)=rdumz(i,1,1)+diffh(i,3)*
!              (1.0d+0/JACM1(i,1,1,e)*
!    >             (ur(i,1,1)*RZM1(i,1,1,e) +
!    >              us(i,1,1)*SZM1(i,1,1,e) +
!    >              ut(i,1,1)*TZM1(i,1,1,e)))
!     enddo

        if (eq_num.eq.4.and.ldim.eq.2)then

        else
           call subcol3(res1(1,1,1,e,eq_num),rdumz(1,1,1)
     >                  ,bm1(1,1,1,e),nxyz)
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
        endif
      elseif(eq_num.eq.5)then
           call subcol3(res1(1,1,1,e,eq_num),usrf(1,1,1,eq_num)
     $                  ,bm1(1,1,1,e),nxyz) 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cmtusrf(e)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'TSTEP'
      include 'PARALLEL'

      integer e,eg

      if(istep.eq.1)then
        n = nx1*ny1*nz1*5
        call rzero(usrf,n)
      endif
      eg = lglel(e)
      do k=1,nz1
         do j=1,ny1
            do i=1,nx1
               call NEKASGN(i,j,k,e)
               call userf(i,j,k,eg)
               usrf(i,j,k,2) = FFX
               usrf(i,j,k,3) = FFY
               usrf(i,j,k,4) = FFZ
               usrf(i,j,k,5) = (U(i,j,k,2,e)*FFX + U(i,j,k,3,e)*FFY
     &                       +  U(i,j,k,4,e)*FFZ)/ U(i,j,k,1,e)
            enddo
         enddo
      enddo

      return
      end 

!-----------------------------------------------------------------------

C> \ingroup vfjac
C> @{
C> flux = \f$\mathscr{A}\f$ dU = \f$\left(\mathscr{A}^{\mbox{NS}}+\mathscr{A}^{\mbox{EVM}}\right) \f$dU 
      subroutine agradu(flux,du,e,eq)
      include 'SIZE'
      include 'CMTDATA'
!JH122716 Apply viscous flux jacobian \mathscr{A} to a notional gradient
!         of the unknowns U (gradU in viscous_cmt and iku, U-{{U}} and
!         U-U_D in igtu_cmt, [U] in ihu (CHECK LI'S NOTES AGAIN))
! Yes, I know theoretically that EVM and NSE have different diffusive
! fluxes for different physical reasons. But I am combining them because
! 1. s_{ij} is shared by both and cuts down on code redundancy
! 2. we may wish to run EVM AND NSE together (I know this is frowned upon)

! but let's consider the normal use case of either NSE OR EVM
! Compressible Navier-Stokes      |
! mu=mu(T)                        |
! lambda=-2/3mu                   |
! nu_s=0                          |
! kappa=kappa(T)                  |  All of this 
!                                 |  should be done in
! Entropy visosity method (EVM)   |  uservp. very error-prone
! mu=mu_s(R_s,|u|+c,h)            |  requires deliberate user attention
! lambda=0 for EVM                |
! nu_s=nu_s(mu_s)                 |
! kappa=0                         |
! I need a flag, ifevm, for controlling calls (entropy_residual
! and longitudinal viscous fluxes, mostly due to grad rho) SOMEDAY
!-----------------------------------------------------------------------
! constructive feedback is always welcome
! flux is zero on entry
!-----------------------------------------------------------------------
      integer e, eq
      real flux(nx1*ny1*nz1,ndim),du(nx1*ny1*nz1,toteq,ndim)

C> \f$\tau_{ij}\f$ and \f$u_j \tau_{ij}\f$.  \f$\lambda=0\f$ and \f$\kappa=0\f$
C> for EVM
      call fluxj_ns (flux,du,e,eq)
C> \f$\nu_s \nabla \rho\f$, \f$\nu_s \left(\nabla \rho \right) \otimes \mathbf{u}\f$
C> and \f$\nu_s \nabla \left(\rho e\right)\f$.  \f$\nu_s=0\f$ for Navier-Stokes
      call fluxj_evm(flux,du,e,eq)

! no idea where phi goes. put it out front
!     call col2(flux,phig(1,1,1,e),nx1*ny1*nz1)

C> @}
      return
      end

!-----------------------------------------------------------------------

C> \f$ \tau_{ij}=2 \mu\sigma_{ij} + \lambda \Delta \delta_{ij}\f$
C> Navier-Stokes, so no mass diffusion. uservp provides properties.
C> Implemented via maxima-generated code
      subroutine fluxj_ns(flux,gradu,e,eq)
! viscous flux jacobian for compressible Navier-Stokes equations (NS)
! SOLN and CMTDATA are indexed, assuming vdiff has been filled by uservp
! somehow. In serious need of debugging and replacement.
      include 'SIZE'
      include 'INPUT'! TRIAGE?
      include 'SOLN' ! TRIAGE?

      parameter (ldd=lx1*ly1*lz1)
      common /ctmp1/ viscscr(lx1,ly1,lz1)
      real viscscr

      integer e,eq
      real flux(nx1*ny1*nz1,ndim),gradu(nx1*ny1*nz1,toteq,ndim)
      integer eijk3(3,3,3)
!     data eijk2 / 0, -1, 1, 0/
      data eijk3
     >/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,0,0/

      n=nx1*ny1*nz1

! This is a disaster that I might want to program less cleverly
      if (eq .lt. toteq) then ! TRIAGE. CAN'T GET AGRADU_NS to WORK
                              ! for ENERGY EQUATION. MAXIMA ROUTINES
                              ! BELOW
!!      do j=1,ndim
!!         do k=1,ndim
!!            ieijk=0
!!!           if (eq .lt. toteq .and. eq .gt. 1) ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
!!            if (eq.gt.1)ieijk=eijk3(eq-1,j,k) ! does this work in 2D?
!!
!!            if (ieijk .eq. 0) then
!!              call agradu_ns(flux(1,j),gradu(1,1,k),viscscr,e,
!!     >                           eq,j,k)
!!            endif
!!         enddo
!!      enddo
! JH110716 Maxima routines added for every viscous flux.
!          agradu_ns has failed all verification checks for homentropic vortex
!          initialization.
!          start over
        if (eq.eq.2) then
           call A21kldUldxk(flux(1,1),gradu,e)
           call A22kldUldxk(flux(1,2),gradu,e)
           call A23kldUldxk(flux(1,3),gradu,e)
        elseif (eq.eq.3) then
           call A31kldUldxk(flux(1,1),gradu,e)
           call A32kldUldxk(flux(1,2),gradu,e)
           call A33kldUldxk(flux(1,3),gradu,e)
        elseif (eq.eq.4) then
           call A41kldUldxk(flux(1,1),gradu,e)
           call A42kldUldxk(flux(1,2),gradu,e)
           call A43kldUldxk(flux(1,3),gradu,e)
        endif

      else ! Energy equation courtesy of thoroughly-checked maxima
           ! until I can get agradu_ns working correctly
         if (if3d) then
            call a53kldUldxk(flux(1,3),gradu,e)
         else
            call rzero(gradu(1,1,3),nx1*ny1*nz1*toteq)
            call rzero(vz(1,1,1,e),nx1*ny1*nz1)
         endif
         call a51kldUldxk(flux(1,1),gradu,e)
         call a52kldUldxk(flux(1,2),gradu,e)
      endif

      return
      end

!-----------------------------------------------------------------------
!     subroutine agradu_ns(gijklu,dut,visco,e,eq,jflux,kdir) ! in junkyard...
!-----------------------------------------------------------------------

C> viscous flux jacobian for entropy viscosity Euler regularization of
C> Guermond and Popov (2014) SIAM JAM 74(2) that do NOT overlap with
C> the compressible Navier-Stokes equations (NS).
      subroutine fluxj_evm(flux,du,e,eq)
! SOLN and CMTDATA are indexed, assuming vdiff has been filled by uservp
! somehow. In serious need of debugging and replacement.
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'! TRIAGE?
      include 'SOLN' ! TRIAGE?
      include 'CMTDATA'

      parameter (ldd=lx1*ly1*lz1)
      common /ctmp1/ viscscr(lx1,ly1,lz1) ! I'ma keep this
      real viscscr

      integer e,eq,eq2
      real flux(nx1*ny1*nz1,ndim),du(nx1*ny1*nz1,toteq,ndim)

      n=nx1*ny1*nz1

! diffusion due to grad rho
      if (eq .eq. 1) then
         do j=1,ndim ! flux+= viscscr*nu_s*grad (rho)
            call addcol3(flux(1,j),vdiff(1,1,1,e,inus),du(1,1,j),n)
         enddo
      else
         if (eq.lt.toteq) then
            call copy(viscscr,du(1,1,eq-1),n)
            call col2(viscscr,vdiff(1,1,1,e,inus),n)
            call addcol3(flux(1,1),viscscr,vx(1,1,1,e),n)
            call addcol3(flux(1,2),viscscr,vy(1,1,1,e),n)
            if (if3d) call addcol3(flux(1,3),viscscr,vz(1,1,1,e),n)

         else ! energy equation

            if(if3d) then ! mass diffusion term
               call vdot3(viscscr,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),
     >                            vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),n)
            else
               call vdot2(viscscr,vx(1,1,1,e),vy(1,1,1,e),
     >                            vx(1,1,1,e),vy(1,1,1,e),n)
            endif
            call col2(viscscr,vdiff(1,1,1,e,inus),n)
            do j=1,ndim
               call addcol3(flux(1,j),du(1,1,j),viscscr,n)
            enddo

            do j=1,ndim
               do eq2=2,ndim+1
                  call col4(viscscr,du(1,eq2,j),u(1,1,1,eq2,e),
     >                           vdiff(1,1,1,e,inus),n)
                  call invcol2(viscscr,vtrans(1,1,1,e,irho),n) ! scr=nu_s*U/rho
                  call sub2(flux(1,j),viscscr,n)
               enddo
               call addcol3(flux(1,j),du(1,toteq,j),vdiff(1,1,1,e,inus),
     >                      n)
            enddo
         endif ! eq<toteq?

      endif ! eq==1?

      return
      end

!-----------------------------------------------------------------------

      subroutine compute_transport_props
! get vdiff props
! viscosity in imu
! second viscosity in ilam; second viscosity is usually -2/3mu
! but we refuse to assume Stokes' hypothesis for the user
! second viscosity=0 in the EVM for Euler gas dynamics
! thermal conductivity in iknd;
! mass diffusivity for EVM in inus
! via nekasgn
      include 'SIZE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'SOLN'
      include 'CMTDATA'

      integer   e

      do e=1,nelt
         ieg=lglel(e)
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            call nekasgn(i,j,k,e)
            call cmtasgn(i,j,k,e)
            call uservp(i,j,k,ieg)
            vdiff(i,j,k,e,imu)  = mu   ! NEKUSE
            vdiff(i,j,k,e,ilam) = lambda!NEKUSE
            vdiff(i,j,k,e,iknd) = udiff! NEKUSE
            vdiff(i,j,k,e,inus) = nu_s ! CMTDATA
         enddo
         enddo
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------
! TRIAGE BELOW UNTIL I CAN FIX AGRADU_NS
!-----------------------------------------------------------------------
      subroutine a51kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5x+cv*lambda*u1*dU4z-kmcvmu*u3*dU4x+cv*lambda*u1*dU3y
     1   -kmcvmu*u2*dU3x+cv*mu*u3*dU2z+cv*mu*u2*dU2y+(cv*lambda-
     2   K+2*cv*mu)*u1*dU2x-cv*lambdamu*u1*u3*dU1z-cv*lambdamu
     3   *u1*u2*dU1y+(K*u3**2-cv*mu*u3**2+K*u2**2-cv*mu*u2**2-cv*la
     4   mbda*u1**2+K*u1**2-2*cv*mu*u1**2-E*K)*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a52kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5y+cv*lambda*u2*dU4z-kmcvmu*u3*dU4y+cv*mu*u3*dU3z+(cv
     1   *lambda-K+2*cv*mu)*u2*dU3y+cv*mu*u1*dU3x-kmcvmu*u1*dU2y+
     2   cv*lambda*u2*dU2x-cv*lambdamu*u2*u3*dU1z+(K*u3**2-cv*mu
     3   *u3**2-cv*lambda*u2**2+K*u2**2-2*cv*mu*u2**2+K*u1**2-cv*mu*
     4   u1**2-E*K)*dU1y-cv*lambdamu*u1*u2*dU1x)/(cv*rho)
      enddo
      return
      end
      subroutine a53kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU3x=dU(i,3,1)
         dU4x=dU(i,4,1)
         dU5x=dU(i,5,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         dU3y=dU(i,3,2)
         dU4y=dU(i,4,2)
         dU5y=dU(i,5,2)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         dU3z=dU(i,3,3)
         dU4z=dU(i,4,3)
         dU5z=dU(i,5,3)
         rho   =vtrans(i,1,1,ie,irho)
         cv    =vtrans(i,1,1,ie,icv)/rho
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         K     =vdiff(i,1,1,ie,iknd)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         E     =U(i,1,1,toteq,ie)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*(dU5z-E*dU1z)+cv*u3*(lambda*dU4z+2*mu*dU4z+lambda*dU3y+lambda
     1   *dU2x)-K*u3*dU4z+cv*mu*u2*(dU4y+dU3z)+cv*mu*u1*(dU4x+dU2z)-
     2   K*u2*dU3z-K*u1*dU2z-cv*(lambda+2*mu)*u3**2*dU1z+K*u3**2*dU1z+
     3   K*u2**2*dU1z-cv*mu*u2**2*dU1z+K*u1**2*dU1z-cv*mu*u1**2*dU1z-c
     4   v*(lambda+mu)*u2*u3*dU1y-cv*(lambda+mu)*u1*u3*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine A21kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=
     >(lambda*(dU4z+dU3y-u3*dU1z-u2*dU1y)+lambdamu*(dU2x-u1*dU1x))/rho
      enddo
      return
      end
      subroutine A22kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU3x=dU(i,3,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         flux(i)=mu*(dU3x+dU2y-u1*dU1y-u2*dU1x)/rho
      enddo
      return
      end
      subroutine A23kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU4x=dU(i,4,1)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4x+dU2z-u1*dU1z-u3*dU1x)/rho
      enddo
      return
      end

      subroutine A31kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU3x=dU(i,3,1)
         dU1y=dU(i,1,2)
         dU2y=dU(i,2,2)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         flux(i)=mu*(dU3x+dU2y-u1*dU1y-u2*dU1x)/rho
      enddo
      return
      end
      subroutine A32kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=(lambda*(dU4z+dU2x-u3*dU1z-u1*dU1x)+
     >   lambdamu*(dU3y-u2*dU1y))/rho
      enddo
      return
      end
      subroutine A33kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1y=dU(i,1,2)
         dU4y=dU(i,4,2)
         dU1z=dU(i,1,3)
         dU3z=dU(i,3,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4y+dU3z-u2*dU1z-u3*dU1y)/rho
      enddo
      return
      end

      subroutine A41kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU4x=dU(i,4,1)
         dU1z=dU(i,1,3)
         dU2z=dU(i,2,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4x+dU2z-u1*dU1z-u3*dU1x)/rho
      enddo
      return
      end
      subroutine A42kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1y=dU(i,1,2)
         dU4y=dU(i,4,2)
         dU1z=dU(i,1,3)
         dU3z=dU(i,3,3)
         rho   =vtrans(i,1,1,ie,irho)
         mu    =vdiff(i,1,1,ie,imu)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         flux(i)=mu*(dU4y+dU3z-u2*dU1z-u3*dU1y)/rho
      enddo
      return
      end
      subroutine A43kldUldxk(flux,dU,ie)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA' ! gradu lurks
      real K,E,kmcvmu,lambdamu
      real dU(nx1*ny1*nz1,toteq,3)
      real flux(nx1*ny1*nz1)
      npt=lx1*ly1*lz1
      do i=1,npt
         dU1x=dU(i,1,1)
         dU2x=dU(i,2,1)
         dU1y=dU(i,1,2)
         dU3y=dU(i,3,2)
         dU1z=dU(i,1,3)
         dU4z=dU(i,4,3)
         rho   =vtrans(i,1,1,ie,irho)
         lambda=vdiff(i,1,1,ie,ilam)
         mu    =vdiff(i,1,1,ie,imu)
         u1    =vx(i,1,1,ie)
         u2    =vy(i,1,1,ie)
         u3    =vz(i,1,1,ie)
         lambdamu=lambda+2.0*mu
         flux(i)=(lambda*(dU3y+dU2x-u2*dU1y-u1*dU1x)+
     >lambdamu*(dU4z-u3*dU1z))/rho
      enddo
      return
      end

!-----------------------------------------------------------------------
      subroutine setup_cmt_param
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CMTDATA'
      INCLUDE 'CMTBCDATA'

      real  MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      
      external MixtPerf_R_CpG, MixtPerf_T_DPR, MixtPerf_C_GRT
     >                 ,MixtPerf_Ho_CpTUVW,MixtPerf_Cp_CvR,MixtPerf_R_M
     >                 ,MixtPerf_G_CpR      

      cip_adhoc=10.0
      cvgref     = param(104)
c     gmaref     = param(105)
      molmass    = param(106)
      muref      = param(107)
      coeflambda = param(108)
      suthcoef   = param(109)
      reftemp    = param(110)
      prlam      = param(111)
      pinfty     = param(112)
      rgasref    = MixtPerf_R_M(molmass,dum)
      cpgref     = MixtPerf_Cp_CvR(cvgref,rgasref)
      gmaref     = MixtPerf_G_CpR(cpgref,rgasref) 
! put these in rea file someday
      ifsip = .false.
      return
      end
c------------------------------------------------------------------------
C> @file driver3_cmt.f routines for primitive variables, usr-file interfaces
C> and properties

C> Compute primitive variables (velocity, thermodynamic state) from 
C> conserved unknowns U
      subroutine compute_primitive_vars
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'CMTDATA'
      include 'SOLN'
      include 'DEALIAS' ! until we are comfortable with setup_convect

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ energy(lx1,ly1,lz1),scr(lx1,ly1,lz1)
      integer e, eq

      nxyz= nx1*ny1*nz1
      ntot=nxyz*nelt

      do e=1,nelt
         call invcol3(vx(1,1,1,e),u(1,1,1,irpu,e),u(1,1,1,irg,e),nxyz)
         call invcol3(vy(1,1,1,e),u(1,1,1,irpv,e),u(1,1,1,irg,e),nxyz)
!        if (if3d)
         call invcol3(vz(1,1,1,e),u(1,1,1,irpw,e),u(1,1,1,irg,e),nxyz)
! first kinetic energy
         if (if3d) then
            call vdot3(scr,
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
         else
            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
         endif
         call invcol2(scr,u(1,1,1,irg,e),nxyz)
         call cmult(scr,0.5,nxyz)
! then to internal energy
         call sub3(energy,u(1,1,1,iret,e),scr,nxyz)
! now mass-specific
         call invcol2(energy,u(1,1,1,irg,e),nxyz)
! don't forget to get density where it belongs
         call invcol3(vtrans(1,1,1,e,irho),u(1,1,1,irg,e),phig(1,1,1,e),
     >                nxyz)
         call tdstate(e,energy)
      enddo

! setup_convect has the desired effect
! if IFPART=F
! if IFCHAR=F
! if IFCONS=T
! if igeom .ne. 1
! if param(99) .ge. 0
!-----------------------------------------------------------------------
!     call setup_convect(0)
!-----------------------------------------------------------------------
! to make life easier until we master this stuff and harmonize even better with nek,
! I'm including 'DEALIAS' and calling set_convect_cons here
      if (nxd.gt.nx1) then
         call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
      else
         call copy(vxd,vx,ntot) 
         call copy(vyd,vy,ntot) 
         call copy(vzd,vz,ntot) 
      endif

      return
      end
!-----------------------------------------------------------------------

C> Compute thermodynamic state for element e from internal energy.
C> usr file.
      subroutine tdstate(e,energy)!,energy)
c compute the gas properties. We will have option to add real gas models
c We have perfect gas law. Cvg is stored full field
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'
      include 'PARALLEL'
      include 'NEKUSE'
      integer   e,eg
      real energy(nx1,ny1,nz1)

      eg = lglel(e)
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         call nekasgn(i,j,k,e)
         call cmtasgn(i,j,k,e)
         e_internal=energy(i,j,k) !nekasgn should do this, but can't
         call cmt_userEOS(i,j,k,eg)
         vtrans(i,j,k,e,icp)= cp*rho
         vtrans(i,j,k,e,icv)= cv*rho
         t(i,j,k,e,1)       = temp
         pr(i,j,k,e)        = pres
         csound(i,j,k,e)    = asnd
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine cmtasgn (ix,iy,iz,e)
      include 'SIZE'
      include 'SOLN'
      include 'CMTDATA'
      include 'NEKUSE'

      integer e,eqnum

      do eqnum=1,toteq
         varsic(eqnum)=u(ix,iy,iz,eqnum,e)  
      enddo
      phi  = phig  (ix,iy,iz,e)
      rho  = vtrans(ix,iy,iz,e,irho)
      pres = pr    (ix,iy,iz,e)
      if (rho.ne.0) then
         cv   = vtrans(ix,iy,iz,e,icv)/rho
         cp   = vtrans(ix,iy,iz,e,icp)/rho
      endif
      asnd = csound(ix,iy,iz,e)
      mu     = vdiff(ix,iy,iz,e,imu)
      udiff  = vdiff(ix,iy,iz,e,iknd)
! MAKE SURE WE''RE NOT USING UTRANS FOR ANYTHING IN pre-v16 code!!
      lambda = vdiff(ix,iy,iz,e,ilam)

      return
      end

!-----------------------------------------------------------------------

      subroutine cmt_ics
! overlaps with setics. -DCMT will require IFDG as well
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      nxyz2=nx2*ny2*nz2       ! Initialize all fields:
      ntot2=nxyz2*nelv
      nxyz1=nx1*ny1*nz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
      ltott=lelt*nxyz1
      ntotcv=lelt*nxyz1*toteq
      call rzero(phig,ltott)
      call rzero(csound,ltott)
      call rzero(vtrans,ltott*ldimt1)
      call rzero(vdiff ,ltott*ldimt1)
      call rzero(u,ntotcv)
      call cmtuic
      if(ifrestart) call my_full_restart !  Check restart files. soon...

C print min values
      xxmax = glmin(xm1,ntott)
      yymax = glmin(ym1,ntott)
      zzmax = glmin(zm1,ntott)

      vxmax = glmin(vx,ntotv)
      vymax = glmin(vy,ntotv)
      vzmax = glmin(vz,ntotv)
      prmax = glmin(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmin(t ,ntott)

      if (nio.eq.0) then
         write(6,19) xxmax,yymax,zzmax
   19    format('Cxyz min  ',5g25.18)
      endif
      if (nio.eq.0) then
         write(6,20) vxmax,vymax,vzmax,prmax,ttmax
   20    format('Cuvwpt min',5g25.18)
      endif

c print max values
      xxmax = glmax(xm1,ntott)
      yymax = glmax(ym1,ntott)
      zzmax = glmax(zm1,ntott)

      vxmax = glmax(vx,ntotv)
      vymax = glmax(vy,ntotv)
      vzmax = glmax(vz,ntotv)
      prmax = glmax(pr,ntot2)

      ntot = nxyz1*nelt
      ttmax = glmax(t ,ntott)

      if (nio.eq.0) then
         write(6,16) xxmax,yymax,zzmax
   16    format('Cxyz max  ',5g25.18)
      endif

      if (nio.eq.0) then
         write(6,17) vxmax,vymax,vzmax,prmax,ttmax
   17    format('Cuvwpt max',5g25.18)
      endif

c     ! save velocity on fine mesh for dealiasing
!     call setup_convect(2) ! check what this does again. might be a good
!                           ! idea, or it might be counterproductive
      if(nio.eq.0) then
        write(6,*) 'done :: set initial conditions, CMT-nek'
        write(6,*) ' '
      endif
      return
      end

!-----------------------------------------------------------------------

      subroutine cmtuic
! overlaps with setics. -DCMT will require IFDG as well
      include 'SIZE'
      include 'SOLN'
      include 'PARALLEL'
      include 'CMTDATA'
      include 'NEKUSE'
      integer e,eg
      do e=1,nelt
         eg = lglel(e)
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1           
            call nekasgn (i,j,k,e)
            call cmtasgn (i,j,k,e)
            call useric  (i,j,k,eg)
            vx(i,j,k,e) = ux
            vy(i,j,k,e) = uy
            vz(i,j,k,e) = uz
            vtrans(i,j,k,e,irho)  = rho
            vtrans(i,j,k,e,icv)= rho*cv
            vtrans(i,j,k,e,icp)= rho*cp
            phig(i,j,k,e)  = phi
            pr(i,j,k,e)    = pres
            u(i,j,k,irg,e) = phi*rho
            u(i,j,k,irpu,e)= phi*rho*ux
            u(i,j,k,irpv,e)= phi*rho*uy
            u(i,j,k,irpw,e)= phi*rho*uz
            u(i,j,k,iret,e)= phi*rho*(cv*temp+0.5*(ux**2+uy**2+uz**2))
            vdiff(i,j,k,e,imu) = mu
            vdiff(i,j,k,e,iknd)= udiff
            vdiff(i,j,k,e,ilam)= lambda
            t(i,j,k,e,1) = temp
         enddo
         enddo
         enddo
      enddo
      return
      end
C> @file ausm.f Riemann solvers and other rocflu miscellany
! ******************************************************************************
!
! Purpose: Compute convective fluxes using AUSM+ scheme.
!
! Description: None.
!
! Input: 
!   nx          x-component of face normal
!   ny          y-component of face normal
!   nz          z-component of face normal
!   nm          Magnitude of face normal
!   fs          Face speed
!   rl          Density of left state
!   ul          x-component of velocity of left state
!   vl          y-component of velocity of left state
!   wl          z-component of velocity of left state   
!   Hl		Total enthalpy of left state
!   al		Speed of sound of left state
!   pl          Pressure of left state
!   rr          Density of right state
!   ur          x-component of velocity of right state
!   vr          y-component of velocity of right state
!   wr          z-component of velocity of right state  
!   pr          Pressure of right state
!   Hr		Total enthalpy of right state
!   ar		Speed of sound of right state
!
! Output: 
!   flx         Fluxes
!   vf          Face velocities ! NOT USED IN CMT-NEK YET
!
! Notes: 
!   1. Liou M.-S., Progress towards an improved CFD method: AUSM+, AIAA Paper
!      95-1701, 1995
!   2. Do not use computation of face speed of sound which leads to exact 
!      capturing of isolated normal shock waves because of robustness problems
!      for unsteady flows and because that formulation is not applicable to 
!      anything but calorically and thermally perfect gases.
!
! ******************************************************************************

C> \ingroup isurf
C> @{
C> Computes inviscid numerical surface flux from AUSM+ Riemann solver
      SUBROUTINE AUSM_FluxFunction(ntot,nx,ny,nz,nm,fs,rl,ul,vl,wl,pl,
     >                         al,tl,rr,ur,vr,wr,pr,ar,tr,flx,cpl,cpr)

!     IMPLICIT NONE ! HAHAHHAHHAHA
! ******************************************************************************
! Definitions and declarations
! ******************************************************************************
      real MixtPerf_Ho_CpTUVW
      external MixtPerf_Ho_CpTUVW

! ==============================================================================
! Arguments
! ==============================================================================
      integer ntot
      REAL al(ntot),ar(ntot),fs(ntot),nm(ntot),nx(ntot),ny(ntot),
     >     nz(ntot),pl(ntot),pr(ntot),rl(ntot),rr(ntot),ul(ntot),
     >     ur(ntot),vl(ntot),vr(ntot),wl(ntot),wr(ntot),cpl(ntot),
     >     cpr(ntot),tl(ntot),tr(ntot)! INTENT(IN) ::
      REAL flx(ntot,5)!,vf(3) ! INTENT(OUT) ::

! ==============================================================================
! Locals
! ==============================================================================

      REAL af,mf,mfa,mfm,mfp,ml,mla,mlp,mr,mra,mrm,pf,ql,qr,vml,vmr,
     >        wtl,wtr,Hl,Hr

! ******************************************************************************
! Start, compute face state
! ******************************************************************************
      call invcol2(cpl,rl,ntot)
      call invcol2(cpr,rr,ntot)

      do i=1,ntot
         Hl = MixtPerf_Ho_CpTUVW(cpl(i),tl(i),ul(i),vl(i),wl(i))
         Hr = MixtPerf_Ho_CpTUVW(cpr(i),tr(i),ur(i),vr(i),wr(i))

         ql = ul(i)*nx(i) + vl(i)*ny(i) + wl(i)*nz(i) - fs(i)
         qr = ur(i)*nx(i) + vr(i)*ny(i) + wr(i)*nz(i) - fs(i)

         af = 0.5*(al(i)+ar(i)) ! NOTE not using original formulation, see note
         ml  = ql/af
         mla = ABS(ml)

         mr  = qr/af
         mra = ABS(mr)    

         IF ( mla .le. 1.0 ) THEN 
            mlp = 0.25*(ml+1.0)*(ml+1.0) + 0.125*(ml*ml-1.0)*(ml*ml-1.0)
            wtl = 0.25*(ml+1.0)*(ml+1.0)*(2.0-ml) +
     >            0.1875*ml*(ml*ml-1.0)*(ml*ml-1.0)
         ELSE
            mlp = 0.5*(ml+mla)
            wtl = 0.5*(1.0+ml/mla)
         END IF ! mla

         IF ( mra .le. 1.0 ) THEN 
            mrm = -0.25*(mr-1.0)*(mr-1.0)-0.125*(mr*mr-1.0)*(mr*mr-1.0)
            wtr = 0.25*(mr-1.0)*(mr-1.0)*(2.0+mr) -
     >            0.1875*mr*(mr*mr-1.0)*(mr*mr-1.0)
         ELSE
            mrm = 0.5*(mr-mra)
            wtr = 0.5*(1.0-mr/mra)
         END IF ! mla

         mf  = mlp + mrm
         mfa = ABS(mf)
         mfp = 0.5*(mf+mfa)
         mfm = 0.5*(mf-mfa)

         pf = wtl*pl(i) + wtr*pr(i)

! ******************************************************************************
! Compute fluxes
! ******************************************************************************

!        vf(1) = mfp*ul + mfm*ur ! I'm sure we'll need this someday
!        vf(2) = mfp*vl + mfm*vr
!        vf(3) = mfp*wl + mfm*wr

         flx(i,1)=(af*(mfp*rl(i)      +mfm*rr(i)   )        )*nm(i)
         flx(i,2)=(af*(mfp*rl(i)*ul(i)+mfm*rr(i)*ur(i))+pf*nx(i))*
     >            nm(i)
         flx(i,3)=(af*(mfp*rl(i)*vl(i)+mfm*rr(i)*vr(i))+pf*ny(i))*
     >            nm(i)
         flx(i,4)=(af*(mfp*rl(i)*wl(i)+mfm*rr(i)*wr(i))+pf*nz(i))*
     >            nm(i)
         flx(i,5)=(af*(mfp*rl(i)*Hl   +mfm*rr(i)*Hr) + pf*fs(i))*
     >            nm(i)
      enddo
C> @}
      return
      END

!-----------------------------------------------------------------------
! NOT LONG FOR THIS WORLD

      SUBROUTINE CentralInviscid_FluxFunction(ntot,nx,ny,nz,fs,ul,pl,
     >                                     ur,pr,flx)
! JH081915 More general, more obvious
! JH111815 HEY GENIUS THIS MAY BE SECOND ORDER AND THUS KILLING YOUR
!          CONVERGENCE. REPLACE WITH AUSM AND SHITCAN IT
! JH112015 This isn't why walls aren't converging. There's something
!          inherently second-order about your wall pressure. Think!
      real nx(ntot),ny(ntot),nz(ntot),fs(ntot),ul(ntot,5),pl(ntot),
     >     ur(ntot,5),pr(ntot) ! intent(in)
      real flx(ntot,5)! intent(out),dimension(5) ::

      do i=1,ntot
         rl =ul(i,1)
         rul=ul(i,2)
         rvl=ul(i,3)
         rwl=ul(i,4)
         rel=ul(i,5)

         rr =ur(i,1)
         rur=ur(i,2)
         rvr=ur(i,3)
         rwr=ur(i,4)
         rer=ur(i,5)

         ql = (rul*nx(i) + rvl*ny(i) + rwl*nz(i))/rl - fs(i)
         qr = (rur*nx(i) + rvr*ny(i) + rwr*nz(i))/rr - fs(i)

         flx(i,1) = 0.5*(ql* rl+ qr*rr               )
         flx(i,2) = 0.5*(ql* rul+pl(i)*nx(i) + qr* rur     +pr(i)*nx(i))
         flx(i,3) = 0.5*(ql* rvl+pl(i)*ny(i) + qr* rvr     +pr(i)*ny(i))
         flx(i,4) = 0.5*(ql* rwl+pl(i)*nz(i) + qr* rwr     +pr(i)*nz(i))
         flx(i,5) = 0.5*(ql*(rel+pl(i))+pl(i)*fs(i)+qr*(rer+pr(i))+
     >               pr(i)*fs(i))
      enddo

      return
      end
!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: Collect relations for static and total speed of sound for perfect
!   gases.
!
! Description: None.
!
! Notes: None.
!
!******************************************************************************
!
! $Id: MixtPerf.f,v 1.5 2015/07/17 15:58:14 mrugeshs Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

      FUNCTION MixtPerf_C_Co2GUVW(Co2,G,U,V,W)
      IMPLICIT NONE
      REAL Co2,G,U,V,W ! INTENT(IN) W
      REAL MixtPerf_C_Co2GUVW
      MixtPerf_C_Co2GUVW = SQRT(Co2 - 0.5*(G - 1.0)*(U*U + V*V + W*W))
      END

      FUNCTION MixtPerf_C_DGP(D,G,P)
      IMPLICIT NONE
      REAL D,G,P! INTENT(IN) 
      REAL MixtPerf_C_DGP
      MixtPerf_C_DGP = SQRT(G*P/D)
      END

      FUNCTION MixtPerf_C_GHoVm2(G,Ho,Vm2)
      IMPLICIT NONE
      REAL G,Ho,Vm2! INTENT(IN) 
      REAL MixtPerf_C_GHoVm2
      MixtPerf_C_GHoVm2 = SQRT((G - 1.0)*(Ho - 0.5*Vm2))
      END

      FUNCTION MixtPerf_C_GRT(G,R,T)
      IMPLICIT NONE
      REAL G,R,T! INTENT(IN) 
      REAL MixtPerf_C_GRT
      MixtPerf_C_GRT = SQRT(G*R*T)
      END

      FUNCTION MixtPerf_C2_GRT(G,R,T)
      IMPLICIT NONE
      REAL G,R,T! INTENT(IN) 
      REAL MixtPerf_C2_GRT
      MixtPerf_C2_GRT = G*R*T
      END

      FUNCTION MixtPerf_Co2_CGUVW(C,G,U,V,W)
      IMPLICIT NONE
      REAL C,G,U,V,W! INTENT(IN) 
      REAL MixtPerf_Co2_CGUVW
      MixtPerf_Co2_CGUVW = C*C + 0.5*(G - 1.0)*(U*U + V*V + W*W)
      END

      FUNCTION MixtPerf_Cv_CpR(Cp,R)
      IMPLICIT NONE
      REAL Cp,R! INTENT(IN) 
      REAL MixtPerf_Cv_CpR
      MixtPerf_Cv_CpR = Cp - R  
      END

      FUNCTION MixtPerf_Cp_CvR(Cv,R)
      IMPLICIT NONE
      REAL Cv,R! INTENT(IN) 
      REAL MixtPerf_Cp_CvR
      MixtPerf_Cp_CvR = Cv + R  
      END

      FUNCTION MixtPerf_D_CGP(C,G,P)
      IMPLICIT NONE
      REAL C,G,P! INTENT(IN) 
      REAL MixtPerf_D_CGP
      MixtPerf_D_CGP = G*P/(C*C)
      END

      FUNCTION MixtPerf_D_DoGMa(D,G,Ma)
      IMPLICIT NONE
      REAL D,G,Ma! INTENT(IN) 
      REAL MixtPerf_D_DoGMa
      MixtPerf_D_DoGMa = D/ (1.0 + 0.5*(G-1.0)*Ma*Ma)**(1.0/(G-1.0))
      END

      FUNCTION MixtPerf_D_PRT(P,R,T)
      IMPLICIT NONE
      REAL P,R,T ! INTENT(IN) 
      REAL MixtPerf_D_PRT
      MixtPerf_D_PRT = P/(R*T)
      END

      FUNCTION MixtPerf_Eo_DGPUVW(D,G,P,U,V,W)
      IMPLICIT NONE
      REAL D,G,P,U,V,W! INTENT(IN) 
      REAL MixtPerf_Eo_DGPUVW
      MixtPerf_Eo_DGPUVW = P/(D*(G - 1.0)) + 0.5*(U*U + V*V + W*W)
      END

      FUNCTION MixtPerf_Eo_DGPVm(D,G,P,Vm)
      IMPLICIT NONE
      REAL D,G,P,Vm! INTENT(IN) 
      REAL MixtPerf_Eo_DGPVm
      MixtPerf_Eo_DGPVm = P/(D*(G - 1.0)) + 0.5*Vm*Vm
      END

      FUNCTION MixtPerf_Eo_GRTUVW(G,R,T,U,V,W)
      IMPLICIT NONE
      REAL  G,R,T,U,V,W
      REAL  MixtPerf_Eo_GRTUVW
      MixtPerf_Eo_GRTUVW = R*T/(G - 1.0) + 0.5*(U*U + V*V + W*W)
      END

      FUNCTION MixtPerf_G_CpR(Cp,R)
      IMPLICIT NONE
      REAL  Cp,R
      REAL  MixtPerf_G_CpR
      MixtPerf_G_CpR = Cp/(Cp - R)
      END

      FUNCTION MixtPerf_Ho_CpTUVW(Cp,T,U,V,W)
      IMPLICIT NONE
      REAL  Cp,T,U,V,W
      REAL  MixtPerf_Ho_CpTUVW
      MixtPerf_Ho_CpTUVW = Cp*T + 0.5*(U*U + V*V + W*W)
      END

      FUNCTION MixtPerf_M_R(R)
      IMPLICIT NONE
      REAL  R
      REAL  MixtPerf_M_R
      MixtPerf_M_R = 8314.3/R
      END

      FUNCTION MixtPerf_P_DEoGVm2(D,Eo,G,Vm2)
      IMPLICIT NONE
      REAL  D,Eo,G,Vm2
      REAL  MixtPerf_P_DEoGVm2
      MixtPerf_P_DEoGVm2 = (G - 1.0)*D*(Eo - 0.5*Vm2)
      END

      FUNCTION MixtPerf_P_DRT(D,R,T)
      IMPLICIT NONE
      REAL  D,R,T
      REAL  MixtPerf_P_DRT
      MixtPerf_P_DRT = D*R*T
      END

      FUNCTION MixtPerf_P_GMaPo(G,Ma,Po)
      IMPLICIT NONE
      REAL  G,Ma,Po
      REAL  MixtPerf_P_GMaPo
      MixtPerf_P_GMaPo = Po/((1.0 + 0.5*(G - 1.0)*Ma*Ma)**(G/(G - 1.0)))
      END

      FUNCTION MixtPerf_P_DDoGPo(D,Do,G,Po)
      IMPLICIT NONE
      REAL  D,Do,G,Po
      REAL  MixtPerf_P_DDoGPo
      MixtPerf_P_DDoGPo = Po*(D/Do)**G
      END

      FUNCTION MixtPerf_P_GPoTTo(G,Po,T,To)
      IMPLICIT NONE
      REAL  G,Po,T,To
      REAL  MixtPerf_P_GPoTTo
      MixtPerf_P_GPoTTo = Po*(T/To)**(G/(G - 1.0))
      END

      FUNCTION MixtPerf_Po_GPTTo(G,P,T,To)
      IMPLICIT NONE
      REAL  G,P,T,To
      REAL  MixtPerf_Po_GPTTo
      MixtPerf_Po_GPTTo = P*(To/T)**(G/(G - 1.0))
      END

      FUNCTION MixtPerf_Po_CGPUVW(C,G,P,U,V,W)
      IMPLICIT NONE
      REAL  C,G,P,U,V,W
      REAL  MixtPerf_Po_CGPUVW
      MixtPerf_Po_CGPUVW =
     >        P*(1.0 + 0.5*(G - 1.0)*(U*U+V*V+W*W)/(C*C))**(G/(G - 1.0))
      END

      FUNCTION MixtPerf_R_CpG(Cp,G)
      IMPLICIT NONE
      REAL  Cp,G
      REAL  MixtPerf_R_CpG
      MixtPerf_R_CpG = Cp - Cp/G  
      END

      FUNCTION MixtPerf_R_M(M,whatev)
      IMPLICIT NONE
      REAL  M,whatev
      REAL  MixtPerf_R_M
      MixtPerf_R_M = 8314.3/M
      END

      FUNCTION MixtPerf_T_CGR(C,G,R)
      IMPLICIT NONE
      REAL  C,G,R
      REAL  MixtPerf_T_CGR
      MixtPerf_T_CGR = C*C/(G*R)
      END

      FUNCTION MixtPerf_T_CpHoVm2(Cp,Ho,Vm2)
      IMPLICIT NONE
      REAL  Cp,Ho,Vm2
      REAL  MixtPerf_T_CpHoVm2
      MixtPerf_T_CpHoVm2 = (Ho-0.5*Vm2)/Cp
      END

      FUNCTION MixtPerf_T_CvEoVm2(Cv,Eo,Vm2)
      IMPLICIT NONE
      REAL  Cv,Eo,Vm2
      REAL  MixtPerf_T_CvEoVm2
      MixtPerf_T_CvEoVm2 = (Eo-0.5*Vm2)/Cv
      END

      FUNCTION MixtPerf_T_DPR(D,P,R)
      IMPLICIT NONE
      REAL  D,P,R
      REAL  MixtPerf_T_DPR
      MixtPerf_T_DPR = P/(D*R)
      END

      FUNCTION MixtPerf_HM_T_DGMP(D,G,M,P)
      IMPLICIT NONE
      REAL  D,G,M,P
      REAL  MixtPerf_HM_T_DGMP
      MixtPerf_HM_T_DGMP = (G*M*M*P + 1.0)/D 
      END

      FUNCTION MixtPerf_T_GMaTo(G,Ma,To)
      IMPLICIT NONE
      REAL  G,Ma,To
      REAL  MixtPerf_T_GMaTo
      MixtPerf_T_GMaTo = To/(1.0 + 0.5*(G - 1.0)*Ma*Ma)
      END

      FUNCTION MixtPerf_To_CpTUVW(Cp,T,U,V,W)
      IMPLICIT NONE
      REAL  Cp,T,U,V,W
      REAL  MixtPerf_To_CpTUVW
      MixtPerf_To_CpTUVW = T + 0.5*(U*U + V*V + W*W)/Cp
      END

      FUNCTION MixtPerf_Vm_C2Co2G(C2,Co2,G)
      IMPLICIT NONE
      REAL  C2,Co2,G
      REAL  MixtPerf_Vm_C2Co2G
      IF ( Co2 .gt. C2 ) THEN  
         MixtPerf_Vm_C2Co2G = SQRT(2.0/(G - 1.0)*(Co2 - C2))
      ELSE 
         MixtPerf_Vm_C2Co2G = 0.0
      END IF ! Co2
      END

! JH060614 stitched bloodily into cmt-nek. Before that,
!******************************************************************************
!
! RCS Revision history:
!
! $Log: MixtPerf.f,v $
! Revision 1.5  2015/07/17 15:58:14  mrugeshs
!  - MS making the code fortran 77 compatible -- not there yet. regression
!       testing going on
!
! Revision 1.4  2015/07/17 15:21:36  jhackl
! JH071715 more fortran 77 goodness
!
! Revision 1.3  2015/03/03 21:09:25  jhackl
! compiles gfortran 4.8.2
!
! Revision 1.2  2015/02/27 22:13:27  mrugeshs
!  - MS022715 - Added function to compute Cp when Cv and R are known
!
! Revision 1.1  2014/06/30 16:42:11  mrugeshs
! - Add CMT souce code
!
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:48:47  haselbac
! Initial revision after changing case
!
! Revision 1.2  2002/05/28 13:44:44  haselbac
! Added new functions
!
! Revision 1.1  2002/05/04 16:16:52  haselbac
! Initial revision
!
!******************************************************************************

      subroutine compute_entropy(s)
! computes entropy at istep and pushes the stack down for previous
! steps needed to compute ds/dt via finite difference (for now).
! hardcoded for Burgers equation. More later when folded into CMT-nek
! for Burgers, s=energy=1/2 U^2
      include 'SIZE'
      include 'TOTAL'  ! tlag is lurking. be careful
      include 'CMTDATA'
! I've always seen lorder=3, but I still need three steps
!          s(:,               1,       1)  ! entropy at current step
!          s(:,               2,       1)  ! entropy at step n-1
!          s(:,               1,       2)  ! entropy at step n-2
      real s(lx1*ly1*lz1*lelt,lorder-1,*)
      real ntol
      integer e
      data icalld /0/
      save icalld

      n=nx1*ny1*nz1
      ntot=n*nelt
      ntol=1.0e-10

      if (icalld .eq. 0) then
         write(6,*) 'zeroing out entropy stack',istep
         icalld=1
         call rzero(s,ntot)
         call rzero(s(1,1,2),ntot) ! s_{n-1}
         call rzero(s(1,2,1),ntot) ! s_n
      endif

! compute the current entropy. This actually needs to go back in the
! usr file because it's EOS-dependent
      rgam=rgasref/(gmaref-1.0)
      do i=1,ntot
         rho=max(vtrans(i,1,1,1,irho),ntol)
         s(i,1,1)=rgam*rho*log(pr(i,1,1,1)/(rho**gmaref))
      enddo

      if (stage .eq. 1) then
! push the stack
         call copy(s(1,1,2),s(1,2,1),ntot) ! s_{n-1}=s_n
         call copy(s(1,2,1),s(1,1,1),ntot) ! fill s_n
      endif

      return
      end

!-----------------------------------------------------------------------
! This may not actually stay in the usr file after further refactoring
C> \ingroup bcond
C> @{
C> Determining IGU contribution to boundary flux. 0 for artificial
C> viscosity, and strictly interior for physical viscosity.
      subroutine bcflux(flux,agradu,qminus)
! Need proper indexing and nekasgn & cmtasgn calls
      include 'SIZE'
      include 'INPUT'
      include 'DG'
!     include 'NEKUSE'
      include 'TSTEP' ! wait how do we know what ifield is?
      integer e,eq,f
      real flux  (nx1*nz1,2*ndim,nelt,toteq)
      real agradu(nx1*nz1,2*ndim,nelt,toteq)
!     real qminus(nx1*nz1,2*ndim,nelt,nqq) ! include CMTDATA?
      real qminus(*) ! 'scuse me. comin' through
      common /nekcb/ cb
      character*3 cb

      nfaces=2*ndim
      nxz=nx1*nz1
      ifield=1

      do e=1,nelt
         do f=1,nfaces
            if (cbc(f, e, ifield).ne.'E  '.and.
     >          cbc(f, e, ifield).ne.'P  ') then ! cbc bndy
               cb=cbc(f,e,ifield)
               if (cb .eq. 'I  ') then ! NEUMANN CONDITIONS GO HERE
!-------------------------------------------------------------
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
                  call rzero(flux(1,f,e,1),nxz)
                  do eq=2,ndim+1
                     call copy(flux(1,f,e,eq),agradu(1,f,e,eq),nxz)
                  enddo
! METHOD "B", ADIABATIC NO-SLIP
                  call rzero(flux(1,f,e,toteq),nxz)
! METHOD "A", ADIABATIC NO-SLIP augments with viscous work. triage below
! because, predictably, NOW I need to computate AgradU on surfaces and I don't
! have general code for that.
                  call a5adiabatic_wall(flux(1,1,1,toteq),f,e,agradu,
     >                                  qminus)
! JH112216 HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
!                 cbu=cb
!                 do eq=1,toteq
!                    call userflux(flux(1,f,e,eq)) ! replace this with userbc
!                 enddo
               else  ! if (cb .eq. 'SYM') then ! NEED PHYSICAL VISC TEST
! JH031617 But this code block basically guarantees that artificial viscosity
!          does not contribute to viscous fluxes at boundaries.
                  do eq=1,toteq
                     call rzero(flux(1,f,e,eq),nxz)
                  enddo
               endif
            endif
         enddo
      enddo

C> @}
      return
      end

      subroutine a5adiabatic_wall(eflx,f,e,dU,wstate)
      include 'SIZE'
      include 'INPUT'
      include 'GEOM' ! for UNX under ADIABATIC WALL METHOD "A"
      include 'CMTDATA'
      real eflx  (nx1*nz1,2*ndim,nelt) ! better be zero on entry
      real dU    (nx1*nz1,2*ndim,nelt,toteq)
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      common /scrns/ flxscr(lx1*lz1)
      real flxscr
      integer e,f

      nxz=nx1*nz1

      call rzero(eflx(1,f,e),nxz)
      call rzero(hface,nxz)

      call a51dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,unx(1,1,f,e),nxz)
      call a52dUadia(flxscr,f,e,dU,wstate)
      call add2col2(eflx(1,f,e),flxsxcr,uny(1,1,f,e),nxz)
      if (if3d) then
         call a53dUadia(flxscr,f,e,dU,wstate)
         call add2col2(eflx(1,f,e),flxsxcr,unz(1,1,f,e),nxz)
      endif
      return
      end

      subroutine a51dUadia(flux,f,ie,dU,wstate)
! same as A51 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5x+cv*lambda*u1*dU4z-kmcvmu*u3*dU4x+cv*lambda*u1*dU3y
     1   -kmcvmu*u2*dU3x+cv*mu*u3*dU2z+cv*mu*u2*dU2y+(cv*lambda-
     2   K+2*cv*mu)*u1*dU2x-cv*lambdamu*u1*u3*dU1z-cv*lambdamu
     3   *u1*u2*dU1y+(K*u3**2-cv*mu*u3**2+K*u2**2-cv*mu*u2**2-cv*la
     4   mbda*u1**2+K*u1**2-2*cv*mu*u1**2-E*K)*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a52dUadia(flux,f,ie,dU,wstate)
! same as A52 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*dU5y+cv*lambda*u2*dU4z-kmcvmu*u3*dU4y+cv*mu*u3*dU3z+(cv
     1   *lambda-K+2*cv*mu)*u2*dU3y+cv*mu*u1*dU3x-kmcvmu*u1*dU2y+
     2   cv*lambda*u2*dU2x-cv*lambdamu*u2*u3*dU1z+(K*u3**2-cv*mu
     3   *u3**2-cv*lambda*u2**2+K*u2**2-2*cv*mu*u2**2+K*u1**2-cv*mu*
     4   u1**2-E*K)*dU1y-cv*lambdamu*u1*u2*dU1x)/(cv*rho)
      enddo
      return
      end

      subroutine a53dUadia(flux,f,ie,dU,wstate)
! same as A53 for volume flux, but
! 1. uses surface storage of quantities in wstate <-qminus (intent(in))
! 2. SETS K=0. ADIABATIC WALLS HAVE VISCOUS HEATING, BUT DON'T CONDUCT
      include 'SIZE'
      include 'CMTDATA'
      real wstate(nx1*nz1,2*ndim,nelt,nqq)
      real dU    (nx1*nz1,2*ndim,nelt,toteq,3)
      real flux  (nx1*ny1*nz1)
      real K,E,kmcvmu,lambdamu
      integer f
      npt=nx1*nz1

      do i=1,npt
         dU1x=dU(i,f,ie,1,1)
         dU2x=dU(i,f,ie,2,1)
         dU3x=dU(i,f,ie,3,1)
         dU4x=dU(i,f,ie,4,1)
         dU5x=dU(i,f,ie,5,1)
         dU1y=dU(i,f,ie,1,2)
         dU2y=dU(i,f,ie,2,2)
         dU3y=dU(i,f,ie,3,2)
         dU4y=dU(i,f,ie,4,2)
         dU5y=dU(i,f,ie,5,2)
         dU1z=dU(i,f,ie,1,3)
         dU2z=dU(i,f,ie,2,3)
         dU3z=dU(i,f,ie,3,3)
         dU4z=dU(i,f,ie,4,3)
         dU5z=dU(i,f,ie,5,3)
         rho   =wstate(i,f,ie,irho)
         cv    =wstate(i,f,ie,icvf)/rho
         lambda=wstate(i,f,ie,ilamf)
         mu    =wstate(i,f,ie,imuf)
         K     =0.0 ! ADIABATIC HARDCODING
         u1    =wstate(i,f,ie,iux)
         u2    =wstate(i,f,ie,iuy)
         u3    =wstate(i,f,ie,iuz)
         E     =wstate(i,f,ie,iu5)/rho
         lambdamu=lambda+mu
         kmcvmu=K-cv*mu
         flux(i)=
     >(K*(dU5z-E*dU1z)+cv*u3*(lambda*dU4z+2*mu*dU4z+lambda*dU3y+lambda
     1   *dU2x)-K*u3*dU4z+cv*mu*u2*(dU4y+dU3z)+cv*mu*u1*(dU4x+dU2z)-
     2   K*u2*dU3z-K*u1*dU2z-cv*(lambda+2*mu)*u3**2*dU1z+K*u3**2*dU1z+
     3   K*u2**2*dU1z-cv*mu*u2**2*dU1z+K*u1**2*dU1z-cv*mu*u1**2*dU1z-c
     4   _v*(lambda+mu)*u2*u3*dU1y-cv*(lambda+mu)*u1*u3*dU1x)/(cv*rho)
      enddo
      return
      end

!-----------------------------------------------------------------------

      subroutine reflect_rind(nvar,f,e,facew,wbc)
      include 'SIZE'
      include 'CMTBCDATA'
      include 'CMTDATA'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'DG'
      include 'MASS'
      include 'TSTEP'
      integer  f,e
! JH091614 facew now has intent(inout)...
! JH031315 not anymore. nobody changes qminus here. that's dumb
      real facew(nx1*nz1,2*ndim,nelt,nvar)
      real wbc(nx1*nz1,2*ndim,nelt,nvar)
      integer i, nxz, fdim
      real nx,ny,nz,rl,ul,vl,wl,pl,fs
      parameter (lfd1=lxd*lzd,lfc1=lx1*lz1)

      nxz=nx1*nz1
      nxzd=nxd*nzd
      fdim=ndim-1
      ieg=lglel(e)
      ifield=1

! I know this says slipwall, but to the inviscid terms all walls are
! slip. or something.
      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      do l=1,nxz
         nx = unx(l,1,f,e)
         ny = uny(l,1,f,e)
         nz = unz(l,1,f,e)
         rl = facew(l,f,e,irho)
         rr = rl
         ul = facew(l,f,e,iux)
         vl = facew(l,f,e,iuy)
         wl = facew(l,f,e,iuz)
         fs = 0.0 ! no moving grid for awhile, and it will not look anything
                  ! like RocFlu

! JH111516 Mirror a la' Dolejsi & Feistauer (2015) section 8.3.1.2
! JH021717 This is for inviscid fluxes, which are produced by the Riemann
!          solver. Presently, this is AUSM, which acts on primitive variables
!          as-coded. Still, it always makes sense to form UBC, so we do so
!          here even though it is different for viscous BC
         udotn = ul*nx+vl*ny+wl*nz
         ur = ul-2.0*udotn*nx
         vr = vl-2.0*udotn*ny
         wr = wl-2.0*udotn*nz
         wbc(l,f,e,irho)= rr
         wbc(l,f,e,iux) = ur
         wbc(l,f,e,iuy) = vr
         wbc(l,f,e,iuz) = wr
         wbc(l,f,e,ipr) = facew(l,f,e,ipr)
         wbc(l,f,e,ithm)= facew(l,f,e,ithm)
         wbc(l,f,e,isnd)= facew(l,f,e,isnd)
         wbc(l,f,e,iph) = facew(l,f,e,iph)
         wbc(l,f,e,icvf)= facew(l,f,e,icvf)
         wbc(l,f,e,icpf)= facew(l,f,e,icpf)
         wbc(l,f,e,iu1)= facew(l,f,e,iu1)
         wbc(l,f,e,iu2)= rr*ur
         wbc(l,f,e,iu3)= rr*vr
         wbc(l,f,e,iu4)= rr*wr
         wbc(l,f,e,iu5)= facew(l,f,e,iu5)
      enddo

      return
      end

!--------------------------------------------------------------------
! NOT LONG FOR THIS WORLD

      subroutine slipwall_rflu(nvar,f,e,facew,wbc,fluxw)
      include 'SIZE'
      include 'CMTBCDATA'
      include 'CMTDATA'
      include 'GEOM'
      include 'NEKUSE'
      include 'INPUT'
      include 'PARALLEL'
      include 'DG'
      include 'MASS'
      include 'TSTEP'
      integer  f,e
! JH091614 facew now has intent(inout)...
! JH031315 not anymore. nobody changes qminus here. that's dumb
      real facew(nx1*nz1,2*ndim,nelt,nvar)
      real wbc(nx1*nz1,2*ndim,nelt,nvar)
      real fluxw(nx1*nz1,2*ndim,nelt,*)
      integer i, nxz, fdim
      real nx,ny,nz,rl,ul,vl,wl,pl,fs
      parameter (lfd1=lxd*lzd,lfc1=lx1*lz1)
      common /SCRNS/ nxf(lfd1),nyf(lfd1),nzf(lfd1),fs2(lfd1),
     >               ufacel(lfd1,5),plc(lfc1),ufacer(lfd1,5),prc(lfd1),
     >               flx(lfd1,5),plf(lfd1),jaco_c(lfc1),
     >               jaco_f(lfd1),dumminus(lfd1,5)
      real nxf,nyf,nzf,ufacel,ufacer,plc,prc,plf,jaco_c,jaco_f,dumminus

      nxz=nx1*nz1
      nxzd=nxd*nzd
      fdim=ndim-1
      ieg=lglel(e)
      ifield=1

! I know this says slipwall, but to the inviscid terms all walls are
! slip. or something.
      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! get molarmass, phi
c                                     ! ux,uy,uz someday
         l=l+1
         nx = unx(l,1,f,e)
         ny = uny(l,1,f,e)
         nz = unz(l,1,f,e)
         rl = facew(l,f,e,irho)
         ul = facew(l,f,e,iux)
         vl = facew(l,f,e,iuy)
         wl = facew(l,f,e,iuz)
         plc(l)= facew(l,f,e,ipr)
         fs = 0.0 ! no moving grid for awhile, and it will not look anything
                  ! like RocFlu
         call RFLU_SetRindStateSlipWallPerf(cp,molarmass,nx,ny,nz,
     >                                      rl,ul,vl,wl,fs,plc(l))
         wbc(l,f,e,irho)=rl

!-----------------------------------------------------------------
! JH111516 INVISCID HARDCODING SLIP WALL. DO THIS SMARTER SOON
! JH111516 Mirror a la' Dolejsi & Feistauer (2015) section 8.3.1.2
!-----------------------------------------------------------------
         udotn=ul*nx+vl*ny+wl*nz
         ur=ul-2.0*udotn*nx
         vr=vl-2.0*udotn*ny
         wr=wl-2.0*udotn*nz
         wbc(l,f,e,iux)=ur
         wbc(l,f,e,iuy)=vr
         wbc(l,f,e,iuz)=wr
! JH111516 SHOULD BE SET TO WALL SPEED i.e. 0 FOR NO-SLIP WALLS!!!
!-----------------------------------------------------------------
         wbc(l,f,e,ipr)=plc(l)! from RFLU_SetRindStateSlipWallPerf
         wbc(l,f,e,iph)=phi
         wbc(l,f,e,iu1)=facew(l,f,e,iu1)
         wbc(l,f,e,iu2)=wbc(l,f,e,iu1)*ur
         wbc(l,f,e,iu3)=wbc(l,f,e,iu1)*vr
         wbc(l,f,e,iu4)=wbc(l,f,e,iu1)*wr
!-------------------------------------------------------------
! JH111516 INVISCID HARDCODING ADIABATIC WALL. DO SMARTER SOON
         wbc(l,f,e,iu5)=facew(l,f,e,iu5)
! JH111516 INVISCID HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
! need a different place to set dirichlet BC for viscous fluxes
!           wbc(l,f,e,iux)=ux ! better b
!           wbc(l,f,e,iuy)=uy
!           wbc(l,f,e,iuz)=uz
!        if (cbc(f,e,ifield) .eq. 'W  ') wbc(l,f,e,ithm)=temp
         plc(l)=plc(l)*phi
      enddo
      enddo
      enddo

! Inviscid flux at walls is due to pressure only. should probably just
! hardcode that instead of calling CentralInviscid so trivially
      if (nxd.gt.nx1) then
         call map_faced(nxf,unx(1,1,f,e),nx1,nxd,fdim,0)
         call map_faced(nyf,uny(1,1,f,e),nx1,nxd,fdim,0)
         call map_faced(nzf,unz(1,1,f,e),nx1,nxd,fdim,0)
         call map_faced(plf,plc,nx1,nxd,fdim,0)

         call invcol3(jaco_c,area(1,1,f,e),wghtc,nxz)
         call map_faced(jaco_f,jaco_c,nx1,nxd,fdim,0)
         call col2(jaco_f,wghtf,nxzd)
      else
         call copy(nxf,unx(1,1,f,e),nxz)
         call copy(nyf,uny(1,1,f,e),nxz)
         call copy(nzf,unz(1,1,f,e),nxz)
         call copy(plf,plc,nxz)

         call copy(jaco_f,area(1,1,f,e),nxz)
      endif
      call rzero(dumminus,toteq*nxzd)
      call map_faced(dumminus(1,1),facew(1,f,e,iu1),nx1,nxd,fdim,0)
      call rzero(fs2,nxzd)
! START BY GETTING RID OF THESE TRIVIAL CENTRAL CALLS AND CENTRAL ALTOGETHER
      call CentralInviscid_FluxFunction(nxzd,nxf,nyf,nzf,fs2,dumminus,
     >                                    plf,dumminus,plf,flx)

      do ieq=1,toteq-1
         call col2(flx(1,ieq),jaco_f,nxzd)
      enddo

      if (nxd.gt.nx1) then
         do j=1,toteq-1
            call map_faced(fluxw(1,f,e,j),flx(1,j),nx1,nxd,fdim,1)
         enddo
         if(cbc(f,e,ifield).ne.'I  ') call map_faced(fluxw(1,f,e,toteq),
     >                              flx(1,toteq),nx1,nxd,fdim,1)
      else
         do j=1,toteq-1
            call copy(fluxw(1,f,e,j),flx(1,j),nxz)
         enddo
         if (cbc(f,e,ifield).ne.'I  ') call copy(fluxw(1,f,e,toteq),
     >                              flx(1,toteq),nxz)
      endif

      return
      end

!-----------------------------------------------------------------------
! ******************************************************************************
!
! Purpose: Set rind state for slip-wall boundaries and perfect gas.
!
! Description: Torn bleeding from RocFlu. I think "rind" means the same thing as
!              "ghost," but I gotta admit that it's a better way of putting it.
!              Not sure if low-order reconstruction lurks here.
! Input:
!   cpGas       Specific heat at constant pressure
!   mmGas       Molecular mass
!   nx,ny,nz    Components of unit normal vector
!   rl          Density
!   ul         x-velocity component
!   vl         y-velocity component
!   wl         z-velocity component
!   fs          Grid speed
!   pl          Pressure
!
! Output: 
!   pl          Pressure
!
! Notes: 
!   1. Valid only for thermally and calorically perfect gas.
!
! ******************************************************************************

      SUBROUTINE RFLU_SetRindStateSlipWallPerf(cpGas,mmGas,nx,ny,nz,rl,
     >                                         ul,vl,wl,fs,pl)

! ******************************************************************************
!   Declarations and definitions
! ******************************************************************************

      real MixtPerf_R_M, MixtPerf_G_CpR, MixtPerf_C_DGP
      external MixtPerf_R_M, MixtPerf_G_CpR, MixtPerf_C_DGP

! ==============================================================================  
!   Arguments 
! ==============================================================================  

      REAL cpGas,fs,mmGas,nx,ny,nz,rl,ul,vl,wl
      REAL pl

! ==============================================================================  
!   Locals 
! ==============================================================================  

      REAL al,gGas,irl,ql,rGas,term
          
! ******************************************************************************
!   Compute wall pressure
! ******************************************************************************

      rGas = MixtPerf_R_M(mmGas)
      gGas = MixtPerf_G_CpR(cpGas,rGas)
 
      irl = 1.0/rl
      ql  = ul*nx + vl*ny + wl*nz - fs
 
      al  = MixtPerf_C_DGP(rl,gGas,pl)

      IF ( ql .lt. 0.0 ) THEN
         term = 1.0 + 0.5*(gGas-1.0)*ql/al
         pl   = pl*term**(2.0*gGas/(gGas-1.0))
      ELSE
         term = (gGas+1.0)/4.0
         pl   = pl + term*rl*ql*(ql + SQRT(al*al+term*term*ql*ql)/term)
      END IF ! ql
 
! ******************************************************************************
!   End
! ******************************************************************************

      end
! NOT LONG FOR THIS WORLD
!--------------------------------------------------------------------

C> @file wall_bc.f Dirichlet states for wall boundary conditions
! FUN FACT: Did you know that bdry.f has a subroutine called
!           BCNEUSC
      subroutine wallbc2(nstate,f,e,facew,wbc)
! DIRICHLET WALL CONDITIONS BECAUSE I DONT KNOW HOW TO INDEX
! UNX in userbc with volume indices instead of face indices
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      include 'SOLN'
      include 'NEKUSE'
      include 'PARALLEL'
      include 'CMTBCDATA'
      include 'CMTDATA'

      integer nstate,f,e
      real    facew(nx1*nz1,2*ndim,nelt,nstate)
      real    wbc(nx1*nz1,2*ndim,nelt,nstate) 
      common /nekcb/ cb
      character*3 cb

      tol=1.0e-10
! JH112116
! rind state for inviscid fluxes is different from viscous fluxes? not
! sure what the right thing to do is.
! JH031617 Collis (CTR 2002-ish), Hartmann & Houston (2006-8) probably BR
!          and Dolejsi and Feistatuer (2015) (check that)
!          all say YES, inviscid rind and viscous rind are different.
      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      ieg=lglel(e)
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg)
         l=l+1

! bring this outside of the face point loop you moron
         if (abs(vdiff(ix,iy,iz,e,ilam)) .gt. tol) then ! physical not artvisc

            wbc(l,f,e,iux)=ux
            wbc(l,f,e,iuy)=uy
            wbc(l,f,e,iuz)=uz
!-----------------------------------------------------------------
! JH112116 I need to check wbc(:,ipr) to make sure it is unchanced
!          from inviscid computation (assuming that's the right
!          answer for general viscous BC).
!-----------------------------------------------------------------
            wbc(l,f,e,iph) =phi
            wbc(l,f,e,ithm)=temp
            wbc(l,f,e,iu1) =facew(l,f,e,iu1)
            wbc(l,f,e,iu2) =wbc(l,f,e,iu1)*ux
            wbc(l,f,e,iu3) =wbc(l,f,e,iu1)*uy
            wbc(l,f,e,iu4) =wbc(l,f,e,iu1)*uz
            if (cb .eq. 'W  ') then ! consider taking properties from userbc too
!           wbc(l,f,e,iu5)=wbc(l,f,e,iu1)*facew(l,f,e,icvf)
               wbc(l,f,e,iu5)=phi*facew(l,f,e,icvf)*temp+
     >         0.5/wbc(l,f,e,iu1)*(wbc(l,f,e,iu2)**2+wbc(l,f,e,iu3)**2+
     >                          wbc(l,f,e,iu4)**2)
            else ! BETTA JUST BE 'I  '
!-------------------------------------------------------------
! JH111516 HARDCODING ADIABATIC WALL. DO SMARTER SOON
!          METHOD "B"
!           wbc(l,f,e,iu5)=facew(l,f,e,iu5)-0.5/facew(l,f,e,iu1)*
!    >     (facew(l,f,e,iu2)**2+facew(l,f,e,iu3)**2+facew(l,f,e,iu4)**2)
!          METHOD "A"
               wbc(l,f,e,iu5)=facew(l,f,e,iu5)
            endif
! JH111516 INVISCID HARDCODING ADIABATIC WALL. DO SMARTER SOON
!-------------------------------------------------------------
         else ! artvisc only

! JH031617 For now, artificial viscosity does not directly contribute to
!          boundary fluxes at all. This means dU=0 for IGTU and gradU is
!          strictly interior for IGU
            do m=1,nqq ! TEST FOR vDIFF OUTSIDE WHAT THE EHLLO IS WRONG WITH YOU
               wbc(l,f,e,m)=facew(l,f,e,m)
            enddo

         endif ! physical viscosity or artvisc
      enddo
      enddo
      enddo

      return
      end

!--------------------------------------------------------------------

      subroutine inflow_rflu(nvar,f,e,facew,wbc)
      include 'SIZE'
      include 'INPUT'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'GEOM'
      include 'PARALLEL'
      include 'DG'
      include 'PERFECTGAS'

      integer f,e,fdim ! intent(in)
      integer i,bcOptType
      real facew(nx1*nz1,2*ldim,nelt,nvar) ! intent(in)
      real wbc  (nx1*nz1,2*ldim,nelt,nvar)   ! intent(out)
      real snx,sny,snz,rhou,rhov,rhow,pl,rhob,rhoub,rhovb,rhowb
     >     ,rhoeb, mach

      nxz=nx1*nz1
      nxzd=nxd*nzd
      fdim=ndim-1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! get molarmass asnd phi t0in p0in cp cv
c                                     !     ux,uy,uz
         l=l+1

         bcOptType=0
         snx  = unx(l,1,f,e)
         sny  = uny(l,1,f,e)

         rho  = facew(l,f,e,iu1)/facew(l,f,e,iph)
         rhou = facew(l,f,e,iu2)/facew(l,f,e,iph)
         rhov = facew(l,f,e,iu3)/facew(l,f,e,iph)
         rhow = facew(l,f,e,iu4)/facew(l,f,e,iph)
         rhoe = facew(l,f,e,iu5)/facew(l,f,e,iph)

         if (if3d) then
            mach = sqrt(ux**2+uy**2+uz**2)/asnd
            snz  = unz(l,1,f,e)
         else
            mach = sqrt(ux**2+uy**2)/asnd
            snz=0.0
         endif
         if (mach.lt.1.0) bcOptType=1

         call BcondInflowPerf(bcOptType,0,p0in,t0in,ux,uy,uz
     >                       ,mach,snx,sny,snz,cp
     >                       ,molarmass,rho,rhou,rhov,rhow,rhob,rhoub
     >                       ,rhovb,rhowb,rhoeb,pres,asnd,temp)
         
         wbc(l,f,e,irho) = rhob
         wbc(l,f,e,iux)  = ux
         wbc(l,f,e,iuy)  = uy
         wbc(l,f,e,iuz)  = uz
         wbc(l,f,e,isnd) = asnd ! overwritten by Bcond
         wbc(l,f,e,ipr)  = pres ! overwritten by Bcond
         wbc(l,f,e,ithm) = temp ! overwritten by Bcond
         wbc(l,f,e,icpf) = rho*cp
         wbc(l,f,e,icvf) = rho*cv
         wbc(l,f,e,iu1)  = rhob*phi
         wbc(l,f,e,iu2)  = rhoub*phi
         wbc(l,f,e,iu3)  = rhovb*phi
         wbc(l,f,e,iu4)  = rhowb*phi
         wbc(l,f,e,iu5)  = rhoeb*phi
      enddo
      enddo
      enddo

      return
      end

!--------------------------------------------------------------------

      subroutine inflow_inviscid(nvar,f,e,facew,wbc)
! JH021717 more conventional Dolejsi & Feistauer (2015),
!          Hartmann & Houston (2006) type boundary conditions
!          Emergency fallback if Holmes just doesn't play nice with DG
      include 'SIZE'
!     include 'TSTEP' ! diagnostics
      include 'INPUT'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'GEOM'
      include 'PARALLEL'
      include 'DG'
      include 'PERFECTGAS'

      integer f,e,fdim ! intent(in)
      integer i
      real facew(nx1*nz1,2*ldim,nelt,nvar) ! intent(in)
      real wbc  (nx1*nz1,2*ldim,nelt,nvar)   ! intent(out)
      real snx,sny,snz,rhou,rhov,rhow,pl,rhob,rhoub,rhovb,rhowb
     >     ,rhoeb, mach

      nxz=nx1*nz1
      nxzd=nxd*nzd
      fdim=ndim-1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! get molarmass asnd phi t0in p0in cp cv
c                                     !     ux,uy,uz
         l=l+1
         wbc(l,f,e,irho) = rho  ! Dirichlet, userbc
         wbc(l,f,e,iux)  = ux   ! Dirichlet, userbc
         wbc(l,f,e,iuy)  = uy   ! Dirichlet, userbc
         wbc(l,f,e,iuz)  = uz   ! Dirichlet, userbc
         wbc(l,f,e,iph)  = phi  ! Dirichlet, userbc
         rhob   = rho*phi
         rhoub  = rho*ux*phi
         rhovb  = rho*uy*phi
         rhowb  = rho*uz*phi
         wbc(l,f,e,iu1)  = rhob
         wbc(l,f,e,iu2)  = rhoub
         wbc(l,f,e,iu3)  = rhovb
         wbc(l,f,e,iu4)  = rhowb

         if (if3d) then ! shouldn't this be normal Mach number?
            mach = sqrt(ux**2+uy**2+uz**2)/asnd
            snz  = unz(l,1,f,e)
         else
            mach = sqrt(ux**2+uy**2)/asnd
            snz=0.0
         endif

         snx  = unx(l,1,f,e)
         sny  = uny(l,1,f,e)

         if (mach.lt.1.0) then

            pres  = facew(l,f,e,ipr) ! extrapolated, overwritten
            temp = pres/rho/(cp-cv) ! definitely too perfect!
            wbc(l,f,e,ipr)  = pres
            wbc(l,f,e,isnd) = sqrt(cp/cv*pres/rho) ! too perfect?
            wbc(l,f,e,ithm) = temp      ! definitely too perfect!
            wbc(l,f,e,icpf) = rho*cp ! NEED EOS WITH TEMP Dirichlet, userbc
            wbc(l,f,e,icvf) = rho*cv ! NEED EOS WITH TEMP Dirichlet, userbc

         else ! supersonic inflow

            wbc(l,f,e,ipr)  = pres
            wbc(l,f,e,isnd) = asnd
            wbc(l,f,e,ithm) = temp
            wbc(l,f,e,icpf) = rho*cp
            wbc(l,f,e,icvf) = rho*cv

         endif

! find a smarter way of doing this. fold it into usr file if you must
         wbc(l,f,e,iu5)  = phi*rho*cv*temp+0.5/rhob*(rhoub**2+rhovb**2+
     >                                               rhowb**2)

      enddo
      enddo
      enddo

      return
      end

!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
!******************************************************************************
!
! Purpose: set inflow boundary condition for one cell.
!
! Description: The subsonic boundary condition is based on the extrapolation of 
!   the Riemann invariant from the interior field (Holmes, D.G.: Inviscid 2D 
!   Solutions on Unstructured, Adaptive Grids. VKI Lecture Series 1989-06, 1989). 
!   The supersonic inflow boundary condition computes the conservative variables
!   from given velocity components, density and pressure.
!
! Input: bcOptType  = boundary treatment: subsonic, supersonic, or mixed
!        bcOptFixed = whether _computed_ inflow angle should be fixed or not
!        ptot       = given total pressure
!        ttot       = given total temperature
!        sx/y/zn    = components of normalized face vector (outward facing)
!        cpgas      = specific heat at constant pressure (boundary cell)
!        mm         = molecular mass at boundary cell
!        rl         = given density
!        ru/v/wl    = given velocity components
!
! Output: rr      = density at boundary | velocity compts inout from userbc
!         ru/v/wr = density * velocity components at boundary
!         rer     = density * total internal energy at boundary
!         pr      = pressure at boundary
!
! Notes: 
!   1. This condition is valid only for thermally and calorically perfect  
!      gas.
!   2. Important to avoid division by MakeNonZero(sl) when computing eta 
!      because that computation can become undefined for reservoir inflow
!      conditions, i.e., for a vanishing velocity vector. 
!
!******************************************************************************
!
! $Id: BcondInflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002-2004 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondInflowPerf(bcOptType,bcOptFixed,ptot,ttot,ur,vr,wr
     >                          ,mach,sxn,syn,szn,cpgas,mm,rl,rul
     >                          ,rvl,rwl,rr,rur,rvr,rwr,rer,pr,ar,tr)
      IMPLICIT NONE

      integer BCOPT_SUBSONIC, BCOPT_MIXED, BCOPT_FIXED_NO
! faked flags
      parameter (BCOPT_SUBSONIC = 1)
      parameter (BCOPT_MIXED = 2)
      parameter (BCOPT_FIXED_NO = 2)
      real  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G
      external  MixtPerf_C_Co2GUVW, MixtPerf_C_DGP, MixtPerf_C_GRT,
     >  MixtPerf_Co2_CGUVW, MixtPerf_C2_GRT,MixtPerf_D_PRT,
     >  MixtPerf_Eo_DGPUVW, MixtPerf_Eo_DGPVm, MixtPerf_G_CpR,
     >  MixtPerf_P_GMaPo, MixtPerf_P_GPoTTo, MixtPerf_Po_GPTTo,
     >  MixtPerf_Po_CGPUVW, MixtPerf_R_M, MixtPerf_T_CGR,
     >  MixtPerf_T_GMaTo,  MixtPerf_Vm_C2Co2G


! ... parameters
      INTEGER bcOptFixed,bcOptType!, INTENT(IN) ::
!          .      .     .                  .   .   .   from userbc
      REAL cpgas, mach, mm, sxn, syn, szn, ur, vr, wr,!, INTENT(IN) ::
     >                        ptot, rl, rul, rvl, rwl, ttot

      REAL rer, rr, rur, rvr, rwr, pr!, INTENT(OUT) ::

      REAL al, ar, a02, cp, disc, eta, g, gm1, igm1, ql, rgas, Rm,
     >            sl, sr, tr, ul,  vl,  wl

!******************************************************************************
! gas properties

      rgas = MixtPerf_R_M(mm)
      g    = MixtPerf_G_CpR(cpgas,rgas)

! subsonic or mixed -----------------------------------------------------------   

      IF ( bcOptType.eq.BCOPT_SUBSONIC.OR.bcOptType.eq.BCOPT_MIXED) THEN
         gm1  = g - 1.0
         igm1 = 1.0/gm1

         ul = rul/rl
         vl = rvl/rl
         wl = rwl/rl

         a02 = MixtPerf_C2_GRT(g,rgas,ttot)

         al = MixtPerf_C_Co2GUVW(a02,g,ul,vl,wl) ! make al consistent with a02
         ql = ul*sxn + vl*syn + wl*szn

! - subsonic

         IF ( ABS(ql) .lt. al ) THEN
            sl = SQRT(ul*ul + vl*vl + wl*wl)
         
            IF ( bcOptFixed .eq. BCOPT_FIXED_NO ) THEN 
               IF ( sl .gt. 1.0E-6 ) THEN ! Avoid ill-defined angle computation
                  eta = ql/sl        
               ELSE 
                  eta = -1.0
               END IF ! sl
            ELSE 
               eta = -1.0
            END IF ! bcOptFixed

            Rm   = al - 0.5*gm1*ql
            disc = 0.5*gm1*eta**2*
     >       (a02/(Rm*Rm)*(1.0 + 0.5*gm1*eta**2) - 1.0)     
        
            IF ( disc .lt. 0.0 ) THEN ! discriminant cannot be negative
               ar = SQRT(a02)
               tr = ttot
               pr = ptot
               sr = 0.0
            ELSE
               ar = Rm/(1.0 + 0.5*gm1*eta*eta)*(1.0+SQRT(disc))
               tr = MixtPerf_T_CGR(ar,g,rgas)
               pr = MixtPerf_P_GPoTTo(g,ptot,tr,ttot)
               sr = MixtPerf_Vm_C2Co2G(ar*ar,a02,g)
            END IF ! disc    
                     
            rr = MixtPerf_D_PRT( pr,rgas,tr )

            rer = rr*MixtPerf_Eo_DGPVm(rr,g,pr,sr)
            rur = rr*ur
            rvr = rr*vr
            rwr = rr*wr

! - supersonic

         ELSE
            IF ( bcOptType .eq. BCOPT_MIXED ) THEN
               pr = mixtPerf_P_GMaPo(g,mach,ptot)
               tr = mixtPerf_T_GMaTo(g,mach,ttot)
               rr = mixtPerf_D_PRT(pr,rgas,tr)
               ar = mixtPerf_C_GRT(g,rgas,tr)

               rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
               rur = rr*ur
               rvr = rr*vr
               rwr = rr*wr
            END IF ! bcOptType
         END IF ! ql < al

! supersonic ------------------------------------------------------------------

      ELSE ! bcOptType .eq. BCOPT_SUPERSONIC
         pr = mixtPerf_P_GMaPo(g,mach,ptot)
         tr = mixtPerf_T_GMaTo(g,mach,ttot)
         rr = mixtPerf_D_PRT(pr,rgas,tr)
         ar = mixtPerf_C_GRT(g,rgas,tr)

         rer = rr*MixtPerf_Eo_DGPUVW(rr,g,pr,ur,vr,wr)
         rur = rr*ur
         rvr = rr*vr
         rwr = rr*wr
      END IF ! bcOptType

      return
      end

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondInflowPerf.F90,v $
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2004/12/01 16:47:56  haselbac
! Initial revision after changing case
!
! Revision 1.7  2004/01/29 22:52:36  haselbac
! Added bcOptFixed, fixed bug, clean-up
!
! Revision 1.6  2003/12/04 03:22:56  haselbac
! Fixed bug in formulation, added partial fix for eta
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************

      subroutine outflow_rflu(nvar,f,e,facew,wbc)
      include 'SIZE'
      include 'NEKUSE'
      include 'CMTDATA'
      include 'CMTBCDATA'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'
      include 'DG'

      integer i,bcOpt
      integer  f,e,fdim
      real facew(nx1*nz1,2*ldim,nelt,nvar)
      real wbc(nx1*nz1,2*ldim,nelt,nvar)
      real sxn,syn,szn,rhou,rhov,rhow,pl,rhob,rhoub,rhovb,rhowb,rhoeb

      nxz=nx1*nz1
      nxzd=nxd*nzd
      fdim=ndim-1
      ieg=lglel(e)

      call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)    
      l=0
      do iz=k0,k1
      do iy=j0,j1
      do ix=i0,i1
         call nekasgn(ix,iy,iz,e)     ! gives us phi- and rho-
         call cmtasgn(ix,iy,iz,e)
         call userbc (ix,iy,iz,f,ieg) ! just for molarmass, and
                                      ! pres
         l=l+1
         sxn = unx(l,1,f,e)
         syn = uny(l,1,f,e)
         szn = unz(l,1,f,e)
         rhou= facew(l,f,e,iu2)/phi
         rhov= facew(l,f,e,iu3)/phi
         rhow= facew(l,f,e,iu4)/phi
         rhoe= facew(l,f,e,iu5)/phi
         pl= facew(l,f,e,ipr) ! P- here
         wbc(l,f,e,icpf)=facew(l,f,e,icpf)
         wbc(l,f,e,icvf)=facew(l,f,e,icvf)
         cp=facew(l,f,e,icpf)/rho
         cv=facew(l,f,e,icvf)/rho
c        fs = 0.0
         if(outflsub)then
            pres= pinfty
         else
            pres= facew(l,f,e,ipr)
         endif
         call BcondOutflowPerf(1,pres,sxn,syn,szn,cp,molarmass,
     >                         rho,rhou,rhov,rhow,rhoe,pl,
     >                         rhob,rhoub,rhovb,rhowb,rhoeb )
         wbc(l,f,e,irho)=rhob
         wbc(l,f,e,iux)=rhoub/rhob
         wbc(l,f,e,iuy)=rhovb/rhob
         wbc(l,f,e,iuz)=rhowb/rhob
! dammit fix this. tdstate to the rescue?
         wbc(l,f,e,ithm)=(rhoeb-0.5*(rhoub**2+rhovb**2+rhowb**2)/rhob)/
     >                   cv
! dammit fix that
         wbc(l,f,e,iu1)=rhob*phi
         wbc(l,f,e,iu2)=rhoub*phi
         wbc(l,f,e,iu3)=rhovb*phi
         wbc(l,f,e,iu4)=rhowb*phi
         wbc(l,f,e,iu5)=rhoeb*phi
         wbc(l,f,e,iph)=phi
         wbc(l,f,e,ipr)=pres
! dammit fix this. tdstate to the rescue?
         wbc(l,f,e,isnd)=sqrt(cp/cv*pres/rho)
! dammit fix that
      enddo
      enddo
      enddo

      return
      end

!******************************************************************************
!
! Purpose: set outflow boundary condition for one cell.
!
! Description: the subsonic boundary condition is based on non-reflecting,
!              characteristics method of Whitfield and Janus: Three-Dimensional
!              Unsteady Euler Equations Solution Using Flux Vector Splitting.
!              AIAA Paper 84-1552, 1984. The supersonic boundary condition
!              consists of simple extrapolation.
!
! Input: bcOpt    = boundary treatment: subsonic, supersonic, or mixed
!        pout     = given static outlet pressure
!        sx/y/zn  = components of ortho-normalized face vector (outward facing)
!        cpgas    = specific heat at constant pressure (boundary cell)
!        mol      = molecular mass at boundary cell
!        rho      = density at boundary cell
!        rhou/v/w = density * velocity components at boundary cell
!        rhoe     = density * total energy at boundary cell
!        press    = static pressure at boundary cell
!
! Output: rhob      = density at boundary
!         rhou/v/wb = density * velocity components at boundary
!         rhoeb     = density * total energy at boundary
!
! Notes: this condition is valid only for thermally and calorically
!        perfect gas (supersonic outflow valid for all gases).
!
!******************************************************************************
!
! $Id: BcondOutflowPerf.F90,v 1.1.1.1 2014/05/05 21:47:47 tmish Exp $
!
! Copyright: (c) 2002 by the University of Illinois
!
!******************************************************************************

      SUBROUTINE BcondOutflowPerf( bcOpt,pout,sxn,syn,szn,cpgas,mol,
     >                       rho,rhou,rhov,rhow,rhoe,press,
     >                       rhob,rhoub,rhovb,rhowb,rhoeb )

      IMPLICIT NONE
      integer bcopt_subsonic,bcopt_mixed
      parameter (BCOPT_SUBSONIC=1)
      parameter (BCOPT_MIXED=0)
      real MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,MixtPerf_R_M
     >   , MixtPerf_P_DEoGVm2
      external MixtPerf_C_DGP,MixtPerf_Eo_DGPUVW,MixtPerf_G_CpR,
     >     MixtPerf_R_M,MixtPerf_P_DEoGVm2

! ... parameters
      INTEGER bcOpt

      REAL pout
      REAL rho, rhou, rhov, rhow, rhoe, press
      REAL sxn, syn, szn, cpgas, mol
      REAL rhob, rhoub, rhovb, rhowb, rhoeb

! ... local variables
      REAL csound, rgas, gamma, gam1, u, v, w, mach, rrhoc, deltp,
     >            ub, vb, wb, vnd

!******************************************************************************
! gas properties; velocity components; Mach number

      rgas  = MixtPerf_R_M( mol )
      gamma = MixtPerf_G_CpR( cpgas,rgas )
      gam1  = gamma - 1.0

      u      = rhou/rho
      v      = rhov/rho
      w      = rhow/rho
      csound = MixtPerf_C_DGP( rho,gamma,press )
      mach   = SQRT(u*u+v*v+w*w)/csound

! subsonic flow ---------------------------------------------------------------

      IF (mach .lt. 1.0 .AND.
     >(bcOpt .eq. BCOPT_SUBSONIC .OR. bcOpt .eq. BCOPT_MIXED)) THEN
         rrhoc = 1.0/(rho*csound)
         deltp = press - pout
         rhob  = rho - deltp/(csound*csound)
         ub    = u + sxn*deltp*rrhoc
         vb    = v + syn*deltp*rrhoc
         wb    = w + szn*deltp*rrhoc

! - special treatment to prevent "deltp" from changing the sign
!   of velocity components. This may happen for very small u, v, w.

         vnd = ub*sxn + vb*syn + wb*szn
         IF ( vnd .lt. 0.0 ) THEN ! inflow at outflow boundary
            ub = SIGN(1.0,u)*MAX(ABS(ub),ABS(u))
            vb = SIGN(1.0,v)*MAX(ABS(vb),ABS(v))
            wb = SIGN(1.0,w)*MAX(ABS(wb),ABS(w))
         END IF ! vnd

         rhoub = rhob*ub
         rhovb = rhob*vb
         rhowb = rhob*wb
         rhoeb = rhob*MixtPerf_Eo_DGPUVW( rhob,gamma,pout,ub,vb,wb )

! supersonic flow -------------------------------------------------------------

      ELSE
         rhob  = rho
         rhoub = rhou
         rhovb = rhov
         rhowb = rhow
         rhoeb = rhoe
      END IF ! mach

      end

!******************************************************************************
!
! RCS Revision history:
!
! $Log: BcondOutflowPerf.F90,v $
! Revision 1.1.1.1  2014/05/05 21:47:47  tmish
! Initial checkin for rocflu macro.
!
! Revision 1.3  2008/12/06 08:43:31  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:46  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:24  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.2  2006/03/26 20:21:09  haselbac
! Fix mistake in declarations
!
! Revision 1.1  2004/12/01 16:48:04  haselbac
! Initial revision after changing case
!
! Revision 1.5  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.2  2002/06/22 00:49:50  jblazek
! Modified interfaces to BC routines.
!
! Revision 1.1  2002/06/10 21:19:34  haselbac
! Initial revision
!
!******************************************************************************

      subroutine InviscidFlux(wminus,wplus,flux,nstate,nflux)
!-------------------------------------------------------------------------------
! JH091514 A fading copy of RFLU_ModAUSM.F90 from RocFlu
!-------------------------------------------------------------------------------

!#ifdef SPEC
!      USE ModSpecies, ONLY: t_spec_type
!#endif
      include 'SIZE'
      include 'INPUT' ! do we need this?
      include 'GEOM' ! for unx
      include 'CMTDATA' ! do we need this without outflsub?
      include 'TSTEP' ! for ifield?
      include 'DG'

! ==============================================================================
! Arguments
! ==============================================================================
      integer nstate,nflux
      real wminus(nx1*nz1,2*ndim,nelt,nstate),
     >     wplus(nx1*nz1,2*ndim,nelt,nstate),
     >     flux(nx1*nz1,2*ndim,nelt,nflux)

! ==============================================================================
! Locals
! ==============================================================================

      integer e,f,fdim,i,k,nxz,nface
      parameter (lfd=lxd*lzd)
! JH111815 legacy rocflu names.
!
! nx,ny,nz : outward facing unit normal components
! fs       : face speed. zero until we have moving grid
! jaco_c   : fdim-D GLL grid Jacobian
! nm       : jaco_c, fine grid
!
! State on the interior (-, "left") side of the face
! rl       : density
! ul,vl,wl : velocity
! tl       : temperature
! al       : sound speed
! pl       : pressure, then phi
! cpl      : rho*cp
! State on the exterior (+, "right") side of the face
! rr       : density
! ur,vr,wr : velocity
! tr       : temperature
! ar       : sound speed
! pr       : pressure
! cpr      : rho*cp

      COMMON /SCRNS/ nx(lfd), ny(lfd), nz(lfd), rl(lfd), ul(lfd),
     >               vl(lfd), wl(lfd), pl(lfd), tl(lfd), al(lfd),
     >               cpl(lfd),rr(lfd), ur(lfd), vr(lfd), wr(lfd),
     >               pr(lfd),tr(lfd), ar(lfd),cpr(lfd),phl(lfd),fs(lfd),
     >               jaco_f(lfd),flx(lfd,toteq),jaco_c(lx1*lz1)
      real nx, ny, nz, rl, ul, vl, wl, pl, tl, al, cpl, rr, ur, vr, wr,
     >                pr,tr, ar,cpr,phl,fs,jaco_f,flx,jaco_c

!     REAL vf(3)
      real nTol
      character*132 deathmessage
      common /nekcb/ cb
      character*3 cb

      nTol = 1.0E-14

      fdim=ndim-1
      nface = 2*ndim
      nxz   = nx1*nz1
      nxzd  = nxd*nzd
      ifield= 1 ! You need to figure out the best way of dealing with
                ! this variable

!     if (outflsub)then
!        call maxMachnumber
!     endif
      do e=1,nelt
      do f=1,nface

! JH021717 Finally corrected BC wrongheadedness. Riemann solver handles
!          all fluxes with even the slightest Dirichlet interpretation,
!          and BC fill wplus sanely for the Riemann solver to provide
!          a good flux for weak BC.
! JH111715 now with dealiased surface integrals. I am too lazy to write
!          something better

! diagnostic
!        if (cbc(f,e,ifield).eq.'v  '.or.cbc(f,e,ifield).eq.'V  ')then
!        if (istep .eq. 1000) then
!           do i=1,nxz
!              write(10+istep,'(2i3,a3,16e15.7)') e,f,cbc(f,e,ifield),
!    .         wminus(i,f,e,irho),wplus(i,f,e,irho),
!    .      wminus(i,f,e,iux), wplus(i,f,e,iux), wminus(i,f,e,iuy),
!    .      wplus(i,f,e,iuy), wminus(i,f,e,ipr), wplus(i,f,e,ipr),
!    .      wminus(i,f,e,ithm), wplus(i,f,e,ithm), wminus(i,f,e,isnd),
!    .      wplus(i,f,e,isnd), wminus(i,f,e,icpf), wplus(i,f,e,icpf),
!    .      wminus(i,f,e,iph), wplus(i,f,e,iph)
!           enddo
!        endif
! diagnostic

         if (nxd.gt.nx1) then
            call map_faced(nx,unx(1,1,f,e),nx1,nxd,fdim,0)
            call map_faced(ny,uny(1,1,f,e),nx1,nxd,fdim,0)
            call map_faced(nz,unz(1,1,f,e),nx1,nxd,fdim,0)

            call map_faced(rl,wminus(1,f,e,irho),nx1,nxd,fdim,0)
            call map_faced(ul,wminus(1,f,e,iux),nx1,nxd,fdim,0)
            call map_faced(vl,wminus(1,f,e,iuy),nx1,nxd,fdim,0)
            call map_faced(wl,wminus(1,f,e,iuz),nx1,nxd,fdim,0)
            call map_faced(pl,wminus(1,f,e,ipr),nx1,nxd,fdim,0)
            call map_faced(tl,wminus(1,f,e,ithm),nx1,nxd,fdim,0)
            call map_faced(al,wminus(1,f,e,isnd),nx1,nxd,fdim,0)
            call map_faced(cpl,wminus(1,f,e,icpf),nx1,nxd,fdim,0)

            call map_faced(rr,wplus(1,f,e,irho),nx1,nxd,fdim,0)
            call map_faced(ur,wplus(1,f,e,iux),nx1,nxd,fdim,0)
            call map_faced(vr,wplus(1,f,e,iuy),nx1,nxd,fdim,0)
            call map_faced(wr,wplus(1,f,e,iuz),nx1,nxd,fdim,0)
            call map_faced(pr,wplus(1,f,e,ipr),nx1,nxd,fdim,0)
            call map_faced(tr,wplus(1,f,e,ithm),nx1,nxd,fdim,0)
            call map_faced(ar,wplus(1,f,e,isnd),nx1,nxd,fdim,0)
            call map_faced(cpr,wplus(1,f,e,icpf),nx1,nxd,fdim,0)

            call map_faced(phl,wminus(1,f,e,iph),nx1,nxd,fdim,0)

            call invcol3(jaco_c,area(1,1,f,e),wghtc,nxz)
            call map_faced(jaco_f,jaco_c,nx1,nxd,fdim,0) 
            call col2(jaco_f,wghtf,nxzd)
         else

            call copy(nx,unx(1,1,f,e),nxz)
            call copy(ny,uny(1,1,f,e),nxz)
            call copy(nz,unz(1,1,f,e),nxz)

            call copy(rl,wminus(1,f,e,irho),nxz)
            call copy(ul,wminus(1,f,e,iux),nxz)
            call copy(vl,wminus(1,f,e,iuy),nxz)
            call copy(wl,wminus(1,f,e,iuz),nxz)
            call copy(pl,wminus(1,f,e,ipr),nxz)
            call copy(tl,wminus(1,f,e,ithm),nxz)
            call copy(al,wminus(1,f,e,isnd),nxz)
            call copy(cpl,wminus(1,f,e,icpf),nxz)

            call copy(rr,wplus(1,f,e,irho),nxz)
            call copy(ur,wplus(1,f,e,iux),nxz)
            call copy(vr,wplus(1,f,e,iuy),nxz)
            call copy(wr,wplus(1,f,e,iuz),nxz)
            call copy(pr,wplus(1,f,e,ipr),nxz)
            call copy(tr,wplus(1,f,e,ithm),nxz)
            call copy(ar,wplus(1,f,e,isnd),nxz)
            call copy(cpr,wplus(1,f,e,icpf),nxz)

            call copy(phl,wminus(1,f,e,iph),nxz)

            call copy(jaco_f,area(1,1,f,e),nxz) 
         endif
         call rzero(fs,nxzd) ! moving grid stuff later

         call AUSM_FluxFunction(nxzd,nx,ny,nz,jaco_f,fs,rl,ul,vl,wl,pl,
     >                          al,tl,rr,ur,vr,wr,pr,ar,tr,flx,cpl,cpr)

         do j=1,toteq
            call col2(flx(1,j),phl,nxzd)
         enddo

         if (nxd.gt.nx1) then
            do j=1,toteq
               call map_faced(flux(1,f,e,j),flx(1,j),nx1,nxd,fdim,1)
            enddo
         else
            do j=1,toteq
               call copy(flux(1,f,e,j),flx(1,j),nxz)
            enddo
         endif

      enddo
      enddo

      end

!-------------------------------------------------------------------------------

      subroutine fill_All_q(fatface)
      include 'SIZE'
      include 'CMTDATA'
      real fatface(*)
      integer eq
      nfq=nx1*nz1*2*ndim*nelt
      nstate = nqq
! where different things live
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq

      call fillq(irho,vtrans,fatface(iwm),fatface(iwp))
      call fillq(iux, vx,    fatface(iwm),fatface(iwp))
      call fillq(iuy, vy,    fatface(iwm),fatface(iwp))
      call fillq(iuz, vz,    fatface(iwm),fatface(iwp))
      call fillq(ipr, pr,    fatface(iwm),fatface(iwp))
      call fillq(ithm,t,     fatface(iwm),fatface(iwp))
      call fillq(isnd,csound,fatface(iwm),fatface(iwp))
      call fillq(iph, phig,  fatface(iwm),fatface(iwp))
      call fillq(icvf,vtrans(1,1,1,1,icv),fatface(iwm),fatface(iwp))
      call fillq(icpf,vtrans(1,1,1,1,icp),fatface(iwm),fatface(iwp))
      call fillq(imuf, vdiff(1,1,1,1,imu), fatface(iwm),fatface(iwp))
      call fillq(ikndf,vdiff(1,1,1,1,iknd),fatface(iwm),fatface(iwp))
      call fillq(ilamf,vdiff(1,1,1,1,ilam),fatface(iwm),fatface(iwp))

      return
      end
