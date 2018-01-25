C> @file drive1_cmt.f high-level driver for CMT-nek
C> \defgroup convhvol Volume integral for inviscid fluxes
C> \defgroup bcond Surface integrals due to boundary conditions
C> \defgroup diffhvol Volume integral for viscous fluxes
C> \defgroup vfjac Jacobians for viscous fluxes
C> \defgroup isurf Inviscid surface terms
C> \defgroup vsurf Viscous surface terms
C> \defgroup faceops utility functions for manipulating face data
C> Branch from subroutine nek_advance in core/drive1.f
C> Advance CMT-nek one time step within nek5000 time loop
      subroutine cmt_nek_advance
c     Solve a conservation law whose flux functions, riemann solvers,
c     jacobians and state relations all live in the usr file

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'SOLN'
      include 'GEOM'
      include 'CTIMER'
      include 'CMTDATA'
      include 'CMTTIMERS'
      
      integer e,eq
      character*32 dumchars

      ftime_dum = dnekclock()
      nxyz1=lx1*ly1*lz1
      n = nxyz1*lelt*toteq
      nfldpart = ndim*npart

      if(istep.eq.1) then
         call cmt_ics
         if (ifrestart) then
            time_cmt=time
         else
            time_cmt=0.0 !time !0.0 ! until we can get settime to behave
         endif
         call cmt_flow_ics
         call init_cmt_timers
         call userchk ! need more ifdefs
         call compute_mesh_h(meshh,xm1,ym1,zm1)
         call compute_grid_h(gridh,xm1,ym1,zm1)
         call compute_primitive_vars ! get good mu
         call entropy_viscosity      ! for high diffno
         call compute_transport_props! at t=0

      endif

      nstage = 3
      do stage=1,nstage
         if (stage.eq.1) call copy(res3(1,1,1,1,1),U(1,1,1,1,1),n)

         rhst_dum = dnekclock()
         call compute_rhs_and_dt
         rhst = rhst + dnekclock() - rhst_dum
c particle equations of motion are solved (also includes forcing)
c In future this subroutine may compute the back effect of particles
c on the fluid and suitably modify the residue computed by 
c compute_rhs_dt for the 5 conserved variables
         call usr_particles_solver

! JH111815 soon....
! JH082316 someday...maybe?
!        do eq=1,toteq
!           call fbinvert(res1(1,1,1,1,eq))
!        enddo

         do e=1,nelt
            do eq=1,toteq
            do i=1,nxyz1
c multiply u with bm1 as res has been multiplied by bm1 in compute_rhs
               u(i,1,1,eq,e) = bm1(i,1,1,e)*tcoef(1,stage)
     >                     *res3(i,1,1,eq,e)+bm1(i,1,1,e)*
     >                     tcoef(2,stage)*u(i,1,1,eq,e)-
     >                     tcoef(3,stage)*res1(i,1,1,e,eq)
c              u(i,1,1,eq,e) = bm1(i,1,1,e)*u(i,1,1,eq,e) - DT *
c    >                        (c1*res1(i,1,1,e,eq) + c2*res2(i,1,1,e,eq)
c    >                       + c3*res3(i,1,1,e,eq))
c-----------------------------------------------------------------------
! JH111815 in fact, I'd like to redo the time marching stuff above and
!          have an fbinvert call for res1
               u(i,1,1,eq,e) = u(i,1,1,eq,e)/bm1(i,1,1,e)
c-----------------------------------------------------------------------
            enddo
            enddo
         enddo
      enddo

      call compute_primitive_vars ! for next time step? Not sure anymore
      call copy(t(1,1,1,1,2),vtrans(1,1,1,1,irho),nxyz1*nelt)
      ftime = ftime + dnekclock() - ftime_dum

      if (mod(istep,iostep).eq.0.or.istep.eq.1)then
         call out_fld_nek
         call mass_balance(if3d)
c dump out particle information. 
         call usr_particles_io(istep)
      end if

!     call print_cmt_timers ! NOT NOW

 101  format(4(2x,e18.9))
      return
      end

c-----------------------------------------------------------------------

C> Compute right-hand-side of the semidiscrete conservation law
C> Store it in res1
      subroutine compute_rhs_and_dt
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'CMTDATA'
      include 'CTIMER'

      integer lfq,heresize,hdsize
      parameter (lfq=lx1*lz1*2*ldim*lelt,
     >                   heresize=nqq*3*lfq,! guarantees transpose of Q+ fits
     >                   hdsize=toteq*3*lfq) ! might not need ldim
! not sure if viscous surface fluxes can live here yet
      common /CMTSURFLX/ flux(heresize),graduf(hdsize)
      real graduf

      integer e,eq
      real wkj(lx1+lxd)
      character*32  dumchars

      call compute_mesh_h(meshh,xm1,ym1,zm1)
      call compute_grid_h(gridh,xm1,ym1,zm1)

      if (nxd.gt.nx1) then
         call set_dealias_face
      else
         call set_alias_rx(istep)
      endif

!     call set_dealias_rx ! done in set_convect_cons,
! JH113015                ! now called from compute_primitive_variables

!     filter the conservative variables before start of each
!     time step
!     if(IFFLTR)  call filter_cmtvar(IFCNTFILT)
!        primitive vars = rho, u, v, w, p, T, phi_g

      call compute_primitive_vars

!-----------------------------------------------------------------------
! JH072914 We can really only proceed with dt once we have current
!          primitive variables. Only then can we compute CFL and/or dt.
!-----------------------------------------------------------------------
      if(stage.eq.1) then
         call setdtcmt
         call set_tstep_coef
      endif

      call entropy_viscosity ! accessed through uservp. computes
                             ! entropy residual and max wave speed
      call compute_transport_props ! everything inside rk stage
!     call smoothing(vdiff(1,1,1,1,imu)) ! still done in usr file
! you have GOT to figure out where phig goes!!!!

      ntot = lx1*ly1*lz1*lelt*toteq
      call rzero(res1,ntot)
      call rzero(flux,heresize)
      call rzero(graduf,hdsize)

!     !Total_eqs = 5 (we will set this up so that it can be a user 
!     !defined value. 5 will be its default value)
!     !eq = 1 -------- Mass balance
!     !eq = 2 -------- x  momentum 
!     !eq = 3 -------- y  momentum 
!     !eq = 4 -------- z  momentum 
!     !eq = 5 -------- Energy Equation 

C> Restrict via \f$\mathbf{E}\f$ to get primitive and conserved variables
C> on interior faces \f$\mathbf{U}^-\f$ and neighbor faces
C> \f$\mathbf{U}^+\f$; store in CMTSURFLX
      call fluxes_full_field

C> res1+=\f$\oint \mathbf{H}^{c\ast}\cdot\mathbf{n}dA\f$ on face points
      nstate=nqq
      nfq=nx1*nz1*2*ndim*nelt
      iwm =1
      iwp =iwm+nstate*nfq
      iflx=iwp+nstate*nfq
      do eq=1,toteq
         ieq=(eq-1)*ndg_face+iflx
         call surface_integral_full(res1(1,1,1,1,eq),flux(ieq))
      enddo
!     dumchars='after_inviscid'
!     call dumpresidue(dumchars,999)

               !                   -
      iuj=iflx ! overwritten with U -{{U}}
!-----------------------------------------------------------------------
!                          /     1  T \
! JH082316 imqqtu computes | I - -QQ  | U for all 5 conserved variables
!                          \     2    /
! which I now make the following be'neon-billboarded assumption:
!***********************************************************************
! ASSUME CONSERVED VARS U1 THROUGH U5 ARE CONTIGUOUSLY STORED
! SEQUENTIALLY IN /CMTSURFLX/ i.e. that iu2=iu1+1, etc.
! CMTDATA BETTA REFLECT THIS!!!
!***********************************************************************
C> res1+=\f$\int_{\Gamma} \{\{\mathbf{A}^{\intercal}\nabla v\}\} \cdot \left[\mathbf{U}\right] dA\f$
      ium=(iu1-1)*nfq+iwm
      iup=(iu1-1)*nfq+iwp
      call   imqqtu(flux(iuj),flux(ium),flux(iup))
      call   imqqtu_dirichlet(flux(iuj),flux(iwm),flux(iwp))
      call igtu_cmt(flux(iwm),flux(iuj),graduf) ! [[u]].{{gradv}}
!     dumchars='after_igtu'
!     call dumpresidue(dumchars,999)

C> res1+=\f$\int \left(\nabla v\right) \cdot \left(\mathbf{H}^c+\mathbf{H}^d\right)dV\f$ 
C> for each equation (inner), one element at a time (outer)
      do e=1,nelt
!-----------------------------------------------------------------------
! JH082216 Since the dawn of CMT-nek we have called this particular loop
!***********************************************************************
!*         "THE" ELEMENT LOOP                                          *
!***********************************************************************
!          since it does several operations, mostly for volume integrals,
!          for all equations, one element at a time. If we get memory
!          under control and GPUs really need to act on gigabytes all
!          at once, then this and its dependents can still have their
!          loop order flipped and things like totalh declared for
!          15 full fields or more.
!-----------------------------------------------------------------------
! Get user defined forcing from userf defined in usr file
         call cmtusrf(e)
         call compute_gradients(e) ! gradU
         do eq=1,toteq
            call convective_cmt(e,eq)        ! convh & totalh -> res1
            call    viscous_cmt(e,eq) ! diffh -> half_iku_cmt -> res1
                                             !       |
                                             !       -> diffh2graduf
! Compute the forcing term in each of the 5 eqs
            call compute_forcing(e,eq)
         enddo
      enddo
      dumchars='after_elm'
      call dumpresidue(dumchars,999)
      call exitt

C> res1+=\f$\int_{\Gamma} \{\{\mathbf{A}\nabla \mathbf{U}\}\} \cdot \left[v\right] dA\f$
      call igu_cmt(flux(iwp),graduf,flux(iwm))
      do eq=1,toteq
         ieq=(eq-1)*ndg_face+iwp
!Finally add viscous surface flux functions of derivatives to res1.
         call surface_integral_full(res1(1,1,1,1,eq),flux(ieq))
      enddo
      dumchars='end_of_rhs'
      call dumpresidue(dumchars,999)
      call exitt

      return
      end
!-----------------------------------------------------------------------
C> Compute coefficients for Runge-Kutta stages \cite{TVDRK}
      subroutine set_tstep_coef

      real tcoef(3,3),dt_cmt,time_cmt
      COMMON /TIMESTEPCOEF/ tcoef,dt_cmt,time_cmt

      tcoef(1,1) = 0.0
      tcoef(2,1) = 1.0 
      tcoef(3,1) = dt_cmt
      tcoef(1,2) = 3.0/4.0
      tcoef(2,2) = 1.0/4.0 
      tcoef(3,2) = dt_cmt/4.0 
      tcoef(1,3) = 1.0/3.0
      tcoef(2,3) = 2.0/3.0 
      tcoef(3,3) = dt_cmt*2.0/3.0 

      return
      end

!-----------------------------------------------------------------------

      subroutine print_cmt_timers
      include 'SIZE'
      include 'CMTTIMERS'
      include 'TSTEP'
      include 'PARALLEL'

c we need our own IO features. Until then we use the default nek routines
      if ((mod(istep,flio_freq).eq.0.and.istep.gt.0)
     $                               .or.istep.eq.nstep)then
         dmtime1 = ftime/istep
         dtime_ = glsum(dmtime1,1)
         if(nio.eq.0) write(6,*) 'fluid rhs compute time(Avg)  '
     $               ,dtime_/np
      endif
      return 
      end
!-----------------------------------------------------------------------

      subroutine init_cmt_timers
      include 'CMTTIMERS'

      rhst    = 0.00
      ftime   = 0.00

      return
      end
!-----------------------------------------------------------------------
