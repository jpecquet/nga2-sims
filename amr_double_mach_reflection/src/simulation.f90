!> Simulation of DOUBLE MACH REFLECTION (Woodward & Colella, 1984) with AMR
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
   real(WP) :: Ms                               !< Shock Mach number
   real(WP) :: rho_pre ,P_pre ,U_pre ,V_pre     !< Pre-shock state
   real(WP) :: rho_post,P_post,U_post,V_post    !< Post-shock state
   real(WP) :: x0_shock                         !< Shock x-intercept with y=0 at t=0
   real(WP) :: theta_shock                      !< Shock angle with respect to the x-axis (rad)
   real(WP) :: Ushock_x                         !< Speed of the shock x-intercept along any line y=const

contains

   !> Heaviside function for shock: 0 pre-shock, 1 post-shock
   real(WP) function shock_heaviside(x,y,theta,x0)
      real(WP), intent(in) :: x,y,theta,x0
      real(WP) :: n
      ! Shock-normal coordinate, x0 is the x-intercept of the shock at y=0
      n=-(x-x0)*sin(theta)+y*cos(theta)
      if (n.lt.0.0_WP) then
         shock_heaviside=0.0_WP
      else
         shock_heaviside=1.0_WP
      end if
   end function shock_heaviside

   !> Set physical viscosities to zero (inviscid)
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=0.0_WP)
      call fs%beta%setval(val=0.0_WP)
      call fs%diff%setval(val=0.0_WP)
   end subroutine get_viscosities

   !> Update viscosities and derived output fields from current primitive variables
   subroutine update_derived()
      implicit none
      ! Compute viscosities and add artificial bulk viscosity
      call get_viscosities()
      call fs%add_viscartif(dt=time%dt)
      ! Compute Umag and Mach number
      call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
      call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
   end subroutine update_derived

   !> Solver init callback - set initial Q
   subroutine fs_init(solver,lvl,t,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: t
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,x,y,H,rho,U,V,P,IE
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointer to Q
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get physical coordinates
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
            ! Set initial conditions
            H=shock_heaviside(x,y,theta_shock,x0_shock)
            rho=rho_pre+(rho_post-rho_pre)*H
            U  =U_pre  +(U_post  -U_pre  )*H
            V  =V_pre  +(V_post  -V_pre  )*H
            P  =P_pre  +(P_post  -P_pre  )*H
            IE =fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=rho*V
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine fs_init

   !> Solver BC callback - post-shock inflow at x-low, half-inflow/half-reflective wall
   !> at y-low, and time-dependent moving-shock inflow at y-high
   subroutine fs_bc(solver,lvl,t,face,bx,comp,pBC)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: t
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pBC
      real(WP) :: dx,dy,x,y,H,rho,U,V,P,IE
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl)
      select case (face)
       case (1) ! x-low: post-shock inflow
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=U_post
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=V_post
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=rho_post
               pBC(i,j,k,2)=rho_post*U_post
               pBC(i,j,k,3)=rho_post*V_post
               pBC(i,j,k,4)=0.0_WP
               pBC(i,j,k,5)=rho_post*fluid%get_e_from_p_rho(p=P_post,rho=rho_post,y=[1.0_WP])
            end do; end do; end do
         end select
       case (3) ! y-low: post-shock inflow for x<x0_shock, reflective wall for x>x0_shock
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+real(i,WP)*dx
               if (x.lt.x0_shock) then
                  pBC(i,j,k,1)=U_post
               else
                  pBC(i,j,k,1)=pBC(i,-j-1,k,1)
               end if
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x.lt.x0_shock) then
                  pBC(i,j,k,1)=V_post
               else
                  if (j.eq.0) then
                     pBC(i,j,k,1)=0.0_WP
                  else
                     pBC(i,j,k,1)=-pBC(i,-j,k,1)
                  end if
               end if
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x.lt.x0_shock) then
                  pBC(i,j,k,1)=rho_post
                  pBC(i,j,k,2)=rho_post*U_post
                  pBC(i,j,k,3)=rho_post*V_post
                  pBC(i,j,k,4)=0.0_WP
                  pBC(i,j,k,5)=rho_post*fluid%get_e_from_p_rho(p=P_post,rho=rho_post,y=[1.0_WP])
               else
                  pBC(i,j,k,1)= pBC(i,-j-1,k,1)
                  pBC(i,j,k,2)= pBC(i,-j-1,k,2)
                  pBC(i,j,k,3)=-pBC(i,-j-1,k,3)
                  pBC(i,j,k,4)= 0.0_WP
                  pBC(i,j,k,5)= pBC(i,-j-1,k,5)
               end if
            end do; end do; end do
         end select
       case (4) ! y-high: time-dependent inflow following the moving shock
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+real(i,WP)*dx
               y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
               H=shock_heaviside(x,y,theta_shock,x0_shock+Ushock_x*t)
               pBC(i,j,k,1)=U_pre+(U_post-U_pre)*H
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               y=solver%amr%ylo+real(j,WP)*dy
               H=shock_heaviside(x,y,theta_shock,x0_shock+Ushock_x*t)
               pBC(i,j,k,1)=V_pre+(V_post-V_pre)*H
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
               H=shock_heaviside(x,y,theta_shock,x0_shock+Ushock_x*t)
               rho=rho_pre+(rho_post-rho_pre)*H
               U  =U_pre  +(U_post  -U_pre  )*H
               V  =V_pre  +(V_post  -V_pre  )*H
               P  =P_pre  +(P_post  -P_pre  )*H
               IE =fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
               pBC(i,j,k,1)=rho
               pBC(i,j,k,2)=rho*U
               pBC(i,j,k,3)=rho*V
               pBC(i,j,k,4)=0.0_WP
               pBC(i,j,k,5)=rho*IE
            end do; end do; end do
         end select
      end select
   end subroutine fs_bc

   !> Solver tagging callback - tag cells with large normalized density gradient
   subroutine fs_tag(solver,lvl,t,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: t
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: drho_ind_x,drho_ind_y,drho_ind
      integer :: i,j,k
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         bx=mfi%tilebox()
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Centered density jumps
            drho_ind_x=0.5_WP*(pQ(i+1,j,k,1)-pQ(i-1,j,k,1))/max(pQ(i,j,k,1),solver%rho_floor)
            drho_ind_y=0.5_WP*(pQ(i,j+1,k,1)-pQ(i,j-1,k,1))/max(pQ(i,j,k,1),solver%rho_floor)
            ! Tag
            drho_ind=sqrt(drho_ind_x**2+drho_ind_y**2)
            if (drho_ind.gt.drho_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine fs_tag

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Initialize EoS and flow parameters
      init_flow: block
         real(WP) :: Gamma,Cv,Wshock,Mn_post,Un_post
         ! Read parameters
         call param_read('Shock Mach number',Ms)
         call param_read('Shock angle',theta_shock)
         theta_shock=4.0_WP*atan(1.0_WP)/180.0_WP*theta_shock ! deg to rad
         x0_shock=1.0_WP/6.0_WP
         ! Create fluid
         call param_read('Gamma',Gamma,default=1.4_WP)
         Cv=1.0_WP/(Gamma-1.0_WP) ! such that gas constant R=1, i.e., T=P/rho
         call fluid%initialize(gamma=Gamma,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         ! Pre-shock state
         rho_pre=1.4_WP; P_pre=1.0_WP ! such that c=1 for Gamma=1.4
         U_pre=0.0_WP; V_pre=0.0_WP
         ! Shock speed (normal to the shock, and of its x-intercept)
         Wshock=Ms*fluid%get_c_from_p_rho(p=P_pre,rho=rho_pre,y=[1.0_WP])
         Ushock_x=Wshock/sin(theta_shock)
         ! Post-shock state from Rankine-Hugoniot
         rho_post=rho_pre*(Gamma+1.0_WP)*Ms**2/(2.0_WP+(Gamma-1.0_WP)*Ms**2)
         P_post=P_pre*(2.0_WP*Gamma*Ms**2-(Gamma-1.0_WP))/(Gamma+1.0_WP)
         Mn_post=sqrt(((Gamma-1.0_WP)*Ms**2+2.0_WP)/(2.0_WP*Gamma*Ms**2-(Gamma-1.0_WP)))
         ! Post-shock normal velocity in the lab frame
         Un_post=Wshock-Mn_post*fluid%get_c_from_p_rho(p=P_post,rho=rho_post,y=[1.0_WP])
         U_post= Un_post*sin(theta_shock)
         V_post=-Un_post*cos(theta_shock)
      end block init_flow

      ! Initialize AMR grid
      init_amrgrid: block
         amr%name='grid'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Max levels',amr%maxlvl)
         amr%xlo=0.0_WP; amr%xhi=4.0_WP
         amr%ylo=0.0_WP; amr%yhi=1.0_WP
         ! Enable quasi-2D: set z extent to match one cell at finest level
         amr%nz=1
         amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         ! Set periodicity
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         call amr%initialize()
      end block init_amrgrid

      ! Initialize time tracker
      init_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax,default=0.2_WP)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         call param_read('Subiterations',time%itmax,default=2)
         time%dt=time%dtmax
      end block init_timetracker

      ! Initialize compressible solver
      init_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrdata_class,    only: interp_face_lin
         ! Assign material and create flow solver
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set initial conditions
         fs%user_init=>fs_init
         ! Set user-defined Dirichlet inflow at x-low and zero-gradient outflow at x-high
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir; fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir; fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir; fs%W%hi_bc(1,1)=amrex_bc_foextrap
         ! Set user-defined BCs at y-low and y-high
         fs%Q%lo_bc(2,:)=amrex_bc_ext_dir; fs%Q%hi_bc(2,:)=amrex_bc_ext_dir
         fs%U%lo_bc(2,1)=amrex_bc_ext_dir; fs%U%hi_bc(2,1)=amrex_bc_ext_dir
         fs%V%lo_bc(2,1)=amrex_bc_ext_dir; fs%V%hi_bc(2,1)=amrex_bc_ext_dir
         fs%W%lo_bc(2,1)=amrex_bc_ext_dir; fs%W%hi_bc(2,1)=amrex_bc_ext_dir
         ! Point to user-defined boundary conditions
         fs%user_bc=>fs_bc
      end block init_solver

      ! Initialize workspaces
      init_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=5,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block init_workspace

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>fs_tag
         call param_read('Tagging Criterion',drho_tag)
      end block init_regridding

      ! Initialize flow state
      init_state: block
         ! Create initial grid and fill Q
         call amr%init_from_scratch(time=time%t)
         ! Initialize primitive variables and face velocity
         call fs%get_primitive(Q=fs%Q); call fs%get_face_velocity()
         ! Update viscosities and derived fields
         call update_derived()
      end block init_state

      ! Initialize visualization
      init_viz: block
         use string, only: str_medium
         character(len=str_medium) :: viz_name
         ! Create visualization object (output in amrviz/<viz_name>/)
         call param_read('Viz name',viz_name,default='flow')
         call viz%initialize(amr,trim(viz_name),use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'rho')
         call viz%add_scalar(fs%I,1,'I')
         call viz%add_scalar(fs%P,1,'P')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(fs%T,1,'T')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block init_viz

      ! Initialize monitors
      init_monitors: block
         ! Get solver info and CFL
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
         call mfile%add_column(fs%Qmin(1),'rhomin')
         call mfile%add_column(fs%Qmax(1),'rhomax')
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
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'Timestep size')
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
         call gridfile%add_column(time%n,'Timestep number')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
      end block init_monitors

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

         ! Iterated midpoint scheme: each sub-iteration re-evaluates dQdt at 0.5*(Qold+Q)
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
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)

            ! Compute face velocities and ensure C/F consistency
            call fs%get_face_velocity(); call fs%average_down_velocity()

            ! Increment both velocities with current pressure term
            call fs%get_primitive(Q=fs%Q)
            call fs%add_pressure(scale=time%dt,phi=fs%P)

            ! Average down and fill ghosts
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)

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

         ! Update viscosities and derived fields
         call update_derived()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

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
