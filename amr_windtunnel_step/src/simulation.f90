!> Wind tunnel with step simulation (Woodward and Colella, 1984) with AMR code
module simulation
   use precision,           only: WP
   use amrgrid_class,       only: amrgrid
   use amrcomp_class,       only: amrcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use ideal_gas_class,     only: ideal_gas
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible solver
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> IBs
   type(amrdata), target :: VF

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: drho_tag=huge(1.0_WP)

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile

   !> Materials
   type(ideal_gas), target :: fluid

   !> Flow parameters
   real(WP) :: P_in,rho,Mach_in,U

contains

   !> Levelset function for step
   function step_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G = max(0.6_WP-xyz(1),xyz(2)-0.2_WP) ! step starts at x=0.6 and is 0.2 tall in y
   end function step_levelset

   !> Set viscosities to zero
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
      call fs%diff%setval(val=0.0_WP)
   end subroutine get_viscosities
    
   !> User init callback - set initial Q
   subroutine user_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,x,y,IE
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointer to Q
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get physical coordinates
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
            ! Set initial conditions
            IE=fluid%get_e_from_p_rho(p=P_in,rho=rho,y=[1.0_WP])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=0.0_WP
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine user_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      select case (face)
       case (1)  ! x-low
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=U
            end do; end do; end do
          case ('V','W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=rho
               p(i,j,k,2)=rho*U
               p(i,j,k,3)=0.0_WP
               p(i,j,k,4)=0.0_WP
               p(i,j,k,5)=rho*fluid%get_e_from_p_rho(p=P_in,rho=rho,y=[1.0_WP])
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   !> Initialize fluid volume fraction
   subroutine init_VF(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      integer, parameter :: nref=3
      real(WP), parameter :: VFlo=1.0e-12_WP
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVF=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=step_levelset,time=time,level=nref,VFlo=VFlo,VF=pVF(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VF

   !> Tagger based on density discontinuities
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: rho, drho_ind_x, drho_ind_y, drho_ind
      integer :: i,j,k
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Local density
            rho=pQ(i,j,k,1)
            ! Centered density jumps
            drho_ind_x=0.5_WP*(pQ(i+1,j,k,1)-pQ(i-1,j,k,1))/max(pQ(i,j,k,1),solver%rho_floor)
            drho_ind_y=0.5_WP*(pQ(i,j+1,k,1)-pQ(i,j-1,k,1))/max(pQ(i,j,k,1),solver%rho_floor)
            ! Tag
            drho_ind=sqrt(drho_ind_x**2+drho_ind_y**2)
            if (drho_ind.gt.drho_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Post-regrid dispatcher for automatic VF filling
   subroutine postregrid_VF(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrdata), pointer :: this
      call c_f_pointer(ctx,this)
      call this%fill(time=time,lbase=lbase)
   end subroutine postregrid_VF

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         real(WP) :: Gamma,Cv
         call param_read('Mach',Mach_in)
         Gamma=1.4_WP
         Cv=1.0_WP/(Gamma-1.0_WP) ! such that gas constant R = 1, i.e. T = p/rho
         rho=1.4_WP; P_in=1.0_WP  ! such that rho/P=1.4=Gamma, and sound speed c=1
         call fluid%initialize(gamma=Gamma,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         ! do this after initializing the fluid!
         U=Mach_in*fluid%get_c_from_p_rho(p=P_in,rho=rho,y=[1.0_WP])
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='tunnel_step'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP; amr%xhi=3.0_WP
         amr%ylo=0.0_WP; amr%yhi=1.0_WP
         amr%zlo=0.0_WP; amr%zhi=1.0_WP
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         call param_read('Max levels',amr%maxlvl)
         ! Enable quasi-2D: set z extent to match one cell at finest level
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax,default=4.0_WP)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         call param_read('Subiterations',time%itmax,default=2)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrdata_class,    only: interp_face_lin
         ! Assign material and create flow solver
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set initial conditions
         fs%user_init=>user_init
         ! Set inflow at X_lo and outflow at X_hi
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir; fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir; fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir; fs%W%hi_bc(1,1)=amrex_bc_foextrap
         fs%user_bc=>shock_dirichlet
         ! Set reflective walls at Y_lo and Y_hi
         fs%Q%lo_bc(2,:)=amrex_bc_foextrap;    fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%Q%lo_bc(2,3)=amrex_bc_reflect_odd; fs%Q%hi_bc(2,3)=amrex_bc_reflect_odd
         fs%U%lo_bc(2,1)=amrex_bc_foextrap;    fs%U%hi_bc(2,1)=amrex_bc_foextrap
         fs%V%lo_bc(2,1)=amrex_bc_reflect_odd; fs%V%hi_bc(2,1)=amrex_bc_reflect_odd
         fs%W%lo_bc(2,1)=amrex_bc_foextrap;    fs%W%hi_bc(2,1)=amrex_bc_foextrap
      end block create_solver

      ! Create VF for IB
      create_VF: block
         use amrex_amr_module, only: amrex_bc_foextrap
         use amrdata_class, only: interp_const
         use iso_c_binding, only: c_loc
         ! Create VF field with constant interpolation
         call VF%initialize(amr,name='VF',ncomp=1,ng=fs%nover,interp=interp_const)
         call VF%register()
         call amr%add_postregrid(postregrid_VF,c_loc(VF))
         VF%user_init=>init_VF
         VF%lo_bc(1,1)=amrex_bc_foextrap; VF%hi_bc(1,1)=amrex_bc_foextrap
         VF%lo_bc(2,1)=amrex_bc_foextrap; VF%hi_bc(2,1)=amrex_bc_foextrap
      end block create_VF

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging Criterion',drho_tag)
         ! Create initial grid
         call amr%init_from_scratch(time=time%t)
         ! Initialize primitive variables and face velocity
         call fs%get_primitive(Q=fs%Q); call fs%get_face_velocity()
         ! Compute viscosities
         call get_viscosities()
         ! Add artificial bulk and shear viscosities
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'TunnelStep',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%T,1,'T')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(VF,1,'VF')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmin,'Pmin')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%Qmin(1),'RHOmin')
         call mfile%add_column(fs%Qmax(1),'RHOmax')
         call mfile%add_column(fs%Imin,'Imin')
         call mfile%add_column(fs%Imax,'Imax')
         call mfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Mass')
         call consfile%add_column(fs%Qint(2),'U Momentum')
         call consfile%add_column(fs%Qint(3),'V Momentum')
         call consfile%add_column(fs%Qint(4),'W Momentum')
         call consfile%add_column(fs%Qint(5),'Internal energy')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old conserved variables and face velocities
         call fs%Qold%copy(src=fs%Q)
         call fs%Uold%copy(src=fs%U)
         call fs%Vold%copy(src=fs%V)
         call fs%Wold%copy(src=fs%W)

         ! Sub-iterations (same scheme as amrcomp_devol/amrcomp_pd)
         do while (time%it.le.time%itmax)

            ! Build midpoint state: Q^{mid}=0.5*(Q+Qold), U^{mid}=0.5*(U+Uold)
            call fs%Q%lincomb(a=0.5_WP,src1=fs%Qold,b=0.5_WP,src2=fs%Q)
            call fs%U%lincomb(a=0.5_WP,src1=fs%Uold,b=0.5_WP,src2=fs%U)
            call fs%V%lincomb(a=0.5_WP,src1=fs%Vold,b=0.5_WP,src2=fs%V)
            call fs%W%lincomb(a=0.5_WP,src1=fs%Wold,b=0.5_WP,src2=fs%W)

            ! Get primitive variables at midpoint
            call fs%get_primitive(Q=fs%Q)

            ! Advance Q using Q^{mid}
            call fs%get_dQdt(dQdt=dQdt)
            call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
            call fs%Q%average_down(); call fs%Q%fill(time%t)

            ! Compute face velocities and ensure C/F consistency
            call fs%get_face_velocity(); call fs%average_down_velocity()

            ! Increment both velocities with current pressure term
            call fs%get_primitive(Q=fs%Q)
            call fs%add_pressure(scale=time%dt,phi=fs%P,mask=VF)

            ! Apply IB direct forcing
            call apply_ib_forcing()

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Pressure correction (only if 'Use projection' is enabled; off by default)
            if (fs%use_projection) then
               call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
               call fs%prepare_psolver(dt=time%dt)
               call fs%psolver%solve(rhs=fs%div)
               call fs%add_pressure(scale=time%dt,mask=VF)
               call apply_ib_forcing()
               call fs%Q%average_down(); call fs%Q%fill(time=time%t)
               call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            end if

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute primitive variables
         call fs%get_primitive(Q=fs%Q)
         
         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add artificial bulk and shear viscosities
         call fs%add_viscartif(dt=time%dt)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         
         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do

   contains

      !> Apply IB forcing - zero Q inside solid and apply quasi-Neumann
      subroutine apply_ib_forcing()
         use amrex_amr_module, only: amrex_mfiter,amrex_box
         type(amrex_mfiter) :: mfi
         type(amrex_box) :: bx
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVF
         real(WP), dimension(:,:,:,:), allocatable :: pQold
         real(WP) :: sum_VF,sum_VFQ1,sum_VFQ5
         integer :: i,j,k,lvl,ii,jj,kk
         ! Compressible IB scheme requires updated ghosts for Q
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Apply IB scheme in solid region
         do lvl=0,amr%clvl()
            call amr%mfiter_build(lvl,mfi)
            do while (mfi%next())
               ! Get pointers to data
               pQ=>fs%Q%mf(lvl)%dataptr(mfi)
               pU=>fs%U%mf(lvl)%dataptr(mfi)
               pV=>fs%V%mf(lvl)%dataptr(mfi)
               pW=>fs%W%mf(lvl)%dataptr(mfi)
               pVF=>VF%mf(lvl)%dataptr(mfi)
               ! Get interior tilebox
               bx=mfi%tilebox()
               ! Create backup of Q
               allocate(pQold,source=pQ)
               ! Loop over tile interior
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Skip pure fluid cells
                  if (pVF(i,j,k,1).eq.1.0_WP) cycle
                  ! Scale Q(2-4) by VF
                  pQ(i,j,k,2)=pVF(i,j,k,1)*pQ(i,j,k,2)
                  pQ(i,j,k,3)=pVF(i,j,k,1)*pQ(i,j,k,3)
                  pQ(i,j,k,4)=pVF(i,j,k,1)*pQ(i,j,k,4)
                  ! VF-weighted neighbor average for Q(1) and Q(5)
                  sum_VF=0.0_WP; sum_VFQ1=0.0_WP; sum_VFQ5=0.0_WP
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                     if (ii.eq.0.and.jj.eq.0.and.kk.eq.0) cycle
                     sum_VF  =sum_VF  +pVF(i+ii,j+jj,k+kk,1)
                     sum_VFQ1=sum_VFQ1+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,1)
                     sum_VFQ5=sum_VFQ5+pVF(i+ii,j+jj,k+kk,1)*pQold(i+ii,j+jj,k+kk,5)
                  end do; end do; end do
                  if (sum_VF.gt.0.0_WP) then
                     pQ(i,j,k,1)=pVF(i,j,k,1)*pQold(i,j,k,1)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ1/sum_VF
                     pQ(i,j,k,5)=pVF(i,j,k,1)*pQold(i,j,k,5)+(1.0_WP-pVF(i,j,k,1))*sum_VFQ5/sum_VF
                  end if
               end do; end do; end do
               ! Deallocate pQold
               deallocate(pQold)
               ! Force face velocities
               bx=mfi%nodaltilebox(1)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pU(i,j,k,1)=0.5_WP*sum(pVF(i-1:i,j,k,1))*pU(i,j,k,1)
               end do; end do; end do
               bx=mfi%nodaltilebox(2)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pV(i,j,k,1)=0.5_WP*sum(pVF(i,j-1:j,k,1))*pV(i,j,k,1)
               end do; end do; end do
               bx=mfi%nodaltilebox(3)
               do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  pW(i,j,k,1)=0.5_WP*sum(pVF(i,j,k-1:k,1))*pW(i,j,k,1)
               end do; end do; end do
            end do
            call amr%mfiter_destroy(mfi)
         end do
      end subroutine apply_ib_forcing
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call VF%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize material
      call fluid%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
   end subroutine simulation_final
   
end module simulation
