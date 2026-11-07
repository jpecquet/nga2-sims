!> 2D Riemann problem simulation with AMR code
module simulation
   use precision,           only: WP
   use amrgrid_class,       only: amrgrid
   use amrcomp_class,       only: amrcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use messager,            only: log,die
   use ideal_gas_class,     only: ideal_gas
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible solver
   type(timetracker) :: time
   type(amrcomp), target :: fs
   type(amrdata) :: dQdt

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Regrid parameters
   type(event) :: regrid_evt
   real(WP) :: drho_tag=huge(1.0_WP)

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile

   !> Flow parameters
   integer :: config    !< Lax & Liu (1998) configuration number (1-19)
   real(WP) :: tend_ref !< End time for selected configuration
   real(WP) :: P1,rho1,U1,V1
   real(WP) :: P2,rho2,U2,V2
   real(WP) :: P3,rho3,U3,V3
   real(WP) :: P4,rho4,U4,V4

   !> Materials
   type(ideal_gas),     target :: fluid

contains

   !> Set viscosities to zero
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
      call fs%diff%setval(val=0.0_WP)
   end subroutine get_viscosities
    
   !> User init callback - set initial Q for 2D Riemann problem
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
      real(WP) :: dx,dy,x,y,rho,P,U,V,W,IE
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
            ! Find quadrant conditions
            if (x > 0.5_WP .and. y > 0.5_WP) then ! quadrant 1
                P=P1; rho=rho1; U=U1; V=V1
            else if (x < 0.5_WP .and. y > 0.5_WP) then ! quadrant 2
                P=P2; rho=rho2; U=U2; V=V2
            else if (x < 0.5_WP .and. y < 0.5_WP) then ! quadrant 3
                P=P3; rho=rho3; U=U3; V=V3
            else ! quadrant 4
                P=P4; rho=rho4; U=U4; V=V4
            end if
            ! Zero out-of-plane velocity
            W=0.0_WP
            ! Internal energy for polytropic gas
            IE=fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=rho*V
            pQ(i,j,k,4)=rho*W
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine user_init

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

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         real(WP) :: Gamma,Cv
         call param_read('Gamma',Gamma)
         call param_read('Configuration',config)
         select case (config)
         case (1) ! t_end=0.2
            P2=0.4_WP; rho2=0.5197_WP; P1=1.0_WP; rho1=1.0_WP
            U2=-0.7259_WP; V2=0.0_WP; U1=0.0_WP; V1=0.0_WP
            P3=0.0439_WP; rho3=0.1072_WP; P4=0.15_WP; rho4=0.2579_WP
            U3=-0.7259_WP; V3=-1.4045_WP; U4=0.0_WP; V4=-1.4045_WP
            tend_ref=0.2_WP
         case (2) ! t_end=0.2
            P2=0.4_WP; rho2=0.5197_WP; P1=1.0_WP; rho1=1.0_WP
            U2=-0.7259_WP; V2=0.0_WP; U1=0.0_WP; V1=0.0_WP
            P3=1.0_WP; rho3=1.0_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=-0.7259_WP; V3=-0.7259_WP; U4=0.0_WP; V4=-0.7259_WP
            tend_ref=0.2_WP
         case (3) ! t_end=0.3
            P2=0.3_WP; rho2=0.5323_WP; P1=1.5_WP; rho1=1.5_WP
            U2=1.206_WP; V2=0.0_WP; U1=0.0_WP; V1=0.0_WP
            P3=0.029_WP; rho3=0.138_WP; P4=0.3_WP; rho4=0.5323_WP
            U3=1.206_WP; V3=1.206_WP; U4=0.0_WP; V4=1.206_WP
            tend_ref=0.3_WP
         case (4) ! t_end=0.25
            P2=0.35_WP; rho2=0.5065_WP; P1=1.1_WP; rho1=1.1_WP
            U2=0.8939_WP; V2=0.0_WP; U1=0.0_WP; V1=0.0_WP
            P3=1.1_WP; rho3=1.1_WP; P4=0.35_WP; rho4=0.5065_WP
            U3=0.8939_WP; V3=0.8939_WP; U4=0.0_WP; V4=0.8939_WP
            tend_ref=0.25_WP
         case (5) ! t_end=0.23
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=-0.75_WP; V2=0.5_WP; U1=-0.75_WP; V1=-0.5_WP
            P3=1.0_WP; rho3=1.0_WP; P4=1.0_WP; rho4=3.0_WP
            U3=0.75_WP; V3=0.5_WP; U4=0.75_WP; V4=-0.5_WP
            tend_ref=0.23_WP
         case (6) ! t_end=0.3
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.75_WP; V2=0.5_WP; U1=0.75_WP; V1=-0.5_WP
            P3=1.0_WP; rho3=1.0_WP; P4=1.0_WP; rho4=3.0_WP
            U3=-0.75_WP; V3=0.5_WP; U4=-0.75_WP; V4=-0.5_WP
            tend_ref=0.3_WP
         case (7) ! t_end=0.25
            P2=0.4_WP; rho2=0.5197_WP; P1=1.0_WP; rho1=1.0_WP
            U2=-0.6259_WP; V2=0.1_WP; U1=0.1_WP; V1=0.1_WP
            P3=0.4_WP; rho3=0.8_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=0.1_WP; V3=0.1_WP; U4=0.1_WP; V4=-0.6259_WP
            tend_ref=0.25_WP
         case (8) ! t_end=0.25
            P2=1.0_WP; rho2=1.0_WP; P1=0.4_WP; rho1=0.5197_WP
            U2=-0.6259_WP; V2=0.1_WP; U1=0.1_WP; V1=0.1_WP
            P3=1.0_WP; rho3=0.8_WP; P4=1.0_WP; rho4=1.0_WP
            U3=0.1_WP; V3=0.1_WP; U4=0.1_WP; V4=-0.6259_WP
            tend_ref=0.25_WP
         case (9) ! t_end=0.3
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=-0.3_WP; U1=0.0_WP; V1=0.3_WP
            P3=0.4_WP; rho3=1.039_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=0.0_WP; V3=-0.8133_WP; U4=0.0_WP; V4=-0.4259_WP
            tend_ref=0.3_WP
         case (10) ! t_end=0.15
            P2=1.0_WP; rho2=0.5_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=0.6076_WP; U1=0.0_WP; V1=0.4297_WP
            P3=0.3333_WP; rho3=0.2281_WP; P4=0.3333_WP; rho4=0.4562_WP
            U3=0.0_WP; V3=-0.6076_WP; U4=0.0_WP; V4=-0.4297_WP
            tend_ref=0.15_WP
         case (11) ! t_end=0.3
            P2=0.4_WP; rho2=0.5313_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.8276_WP; V2=0.0_WP; U1=0.1_WP; V1=0.0_WP
            P3=0.4_WP; rho3=0.8_WP; P4=0.4_WP; rho4=0.5313_WP
            U3=0.1_WP; V3=0.0_WP; U4=0.1_WP; V4=0.7276_WP
            tend_ref=0.3_WP
         case (12) ! t_end=0.25
            P2=1.0_WP; rho2=1.0_WP; P1=0.4_WP; rho1=0.5313_WP
            U2=0.7276_WP; V2=0.0_WP; U1=0.0_WP; V1=0.0_WP
            P3=1.0_WP; rho3=0.8_WP; P4=1.0_WP; rho4=1.0_WP
            U3=0.0_WP; V3=0.0_WP; U4=0.0_WP; V4=0.7276_WP
            tend_ref=0.25_WP
         case (13) ! t_end=0.3
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=0.3_WP; U1=0.0_WP; V1=-0.3_WP
            P3=0.4_WP; rho3=1.0625_WP; P4=0.4_WP; rho4=0.5313_WP
            U3=0.0_WP; V3=0.8145_WP; U4=0.0_WP; V4=0.4276_WP
            tend_ref=0.3_WP
         case (14) ! t_end=0.1
            P2=8.0_WP; rho2=1.0_WP; P1=8.0_WP; rho1=2.0_WP
            U2=0.0_WP; V2=-1.2172_WP; U1=0.0_WP; V1=-0.5606_WP
            P3=2.6667_WP; rho3=0.4736_WP; P4=2.6667_WP; rho4=0.9474_WP
            U3=0.0_WP; V3=1.2172_WP; U4=0.0_WP; V4=1.1606_WP
            tend_ref=0.1_WP
         case (15) ! t_end=0.2
            P2=0.4_WP; rho2=0.5197_WP; P1=1.0_WP; rho1=1.0_WP
            U2=-0.6259_WP; V2=-0.3_WP; U1=0.1_WP; V1=-0.3_WP
            P3=0.4_WP; rho3=0.8_WP; P4=0.4_WP; rho4=0.5313_WP
            U3=0.1_WP; V3=-0.3_WP; U4=0.1_WP; V4=0.4276_WP
            tend_ref=0.2_WP
         case (16) ! t_end=0.2
            P2=1.0_WP; rho2=1.0222_WP; P1=0.4_WP; rho1=0.5313_WP
            U2=-0.6179_WP; V2=0.1_WP; U1=0.1_WP; V1=0.1_WP
            P3=1.0_WP; rho3=0.8_WP; P4=1.0_WP; rho4=1.0_WP
            U3=0.1_WP; V3=0.1_WP; U4=0.1_WP; V4=0.8276_WP
            tend_ref=0.2_WP
         case (17) ! t_end=0.3
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=-0.3_WP; U1=0.0_WP; V1=-0.4_WP
            P3=0.4_WP; rho3=1.0625_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=0.0_WP; V3=0.2145_WP; U4=0.0_WP; V4=-1.1259_WP
            tend_ref=0.3_WP
         case (18) ! t_end=0.2
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=-0.3_WP; U1=0.0_WP; V1=1.0_WP
            P3=0.4_WP; rho3=1.0625_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=0.0_WP; V3=0.2145_WP; U4=0.0_WP; V4=0.2741_WP
            tend_ref=0.2_WP
         case (19) ! t_end=0.3
            P2=1.0_WP; rho2=2.0_WP; P1=1.0_WP; rho1=1.0_WP
            U2=0.0_WP; V2=-0.3_WP; U1=0.0_WP; V1=0.3_WP
            P3=0.4_WP; rho3=1.0625_WP; P4=0.4_WP; rho4=0.5197_WP
            U3=0.0_WP; V3=0.2145_WP; U4=0.0_WP; V4=-0.4259_WP
            tend_ref=0.3_WP
         case default
            call die('[simulation_init] Configuration must be between 1 and 19')
         end select
         ! Log the chosen configuration
         log_config: block
            use string, only: str_long
            character(len=str_long) :: message
            write(message,'("[2D Riemann] Lax-Liu config. ",i0,", T=",f5.3)') config,tend_ref
            call log(message)
         end block log_config
         Cv=1.0_WP/(Gamma-1.0_WP) ! Such that gas constant R = 1, i.e. T = p/rho
         call fluid%initialize(gamma=Gamma,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='riemann2d'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP; amr%xhi=1.0_WP
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
         call param_read('Max time',time%tmax,default=tend_ref)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         call param_read('Subiterations',time%itmax,default=2)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Initialize compressible solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap
         use amrdata_class,    only: interp_face_lin
         ! Assign material and create flow solver
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set initial conditions
         fs%user_init=>user_init
         ! Set BCs: outflow (foextrap) on X and Y boundaries
         if (.not.amr%xper) then
            fs%Q%lo_bc(1,:)=amrex_bc_foextrap
            fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%U%lo_bc(1,:)=amrex_bc_foextrap; fs%U%hi_bc(1,:)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)=amrex_bc_foextrap; fs%V%hi_bc(1,:)=amrex_bc_foextrap
            fs%W%lo_bc(1,:)=amrex_bc_foextrap; fs%W%hi_bc(1,:)=amrex_bc_foextrap
         end if
         if (.not.amr%yper) then
            fs%Q%lo_bc(2,:)=amrex_bc_foextrap
            fs%Q%hi_bc(2,:)=amrex_bc_foextrap
            fs%U%lo_bc(2,:)=amrex_bc_foextrap; fs%U%hi_bc(2,:)=amrex_bc_foextrap
            fs%V%lo_bc(2,:)=amrex_bc_foextrap; fs%V%hi_bc(2,:)=amrex_bc_foextrap
            fs%W%lo_bc(2,:)=amrex_bc_foextrap; fs%W%hi_bc(2,:)=amrex_bc_foextrap
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=interp_none); call dQdt%register()
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
      end block init_regridding

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'Riemann2D',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%T,1,'T')
         call viz%add_scalar(fs%beta,1,'beta')
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

         ! Sub-iterations (same scheme as amrcomp_devol/amrcomp_pd, without IB/solid coupling)
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
            call fs%add_pressure(scale=time%dt,phi=fs%P)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

            ! Pressure correction (only if 'Use projection' is enabled; off by default)
            if (fs%use_projection) then
               call fs%get_div(); call fs%div%mult(val=1.0_WP/time%dt)
               call fs%prepare_psolver(dt=time%dt)
               call fs%psolver%solve(rhs=fs%div)
               call fs%add_pressure(scale=time%dt)
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
         
         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         
      end do
      
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
