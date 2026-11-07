!> Double Mach reflection simulation (Woodward and Colella, 1984) with AMR code
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
   real(WP) :: M1,P1,rho1,U1,V1
   real(WP) :: M2,P2,rho2,U2,V2
   real(WP) :: x0_shock,theta_shock,Wshock,Ushock

contains

   !> Heaviside function for shock
   real(WP) function Hshock(x,y,theta,x0,delta)
      real(WP), intent(in) :: x,y,theta,x0,delta
      real(WP) :: n ! shock-normal coordinate
      n = -(x-x0)*sin(theta) + y*cos(theta) ! x0 is x-intercept of shock (y=0)
      if (n.lt.0.0_WP) then
          Hshock=0.0_WP ! pre-shock
      else
          Hshock=1.0_WP ! post-shock
      end if
   end function Hshock

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
      real(WP) :: dx,dy,x,y,P,rho,U,V,IE
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
            U=U1+(U2-U1)*Hshock(x,y,theta_shock,x0_shock,dx)
            V=V1+(V2-V1)*Hshock(x,y,theta_shock,x0_shock,dx)
            P=P1+(P2-P1)*Hshock(x,y,theta_shock,x0_shock,dx)
            rho=rho1+(rho2-rho1)*Hshock(x,y,theta_shock,x0_shock,dx)
            IE=fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=rho*V
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine user_init

   !> Apply Dirichlet BCs for the inflow, moving-shock top boundary,
   !> and half-inflow, half-reflective wall bottom boundary
   subroutine dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: dx,dy,x,y,pres,rho,U,V,IE,Hshock_t
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl)
      select case (face)
       case (1)  ! Xlo, post-shock conditions inflow
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=U2
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=V2
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=rho2
               p(i,j,k,2)=rho2*U2
               p(i,j,k,3)=rho2*V2
               p(i,j,k,4)=0.0_WP
               p(i,j,k,5)=rho2*fluid%get_e_from_p_rho(p=P2,rho=rho2,y=[1.0_WP])
            end do; end do; end do
         end select
     case (3) ! Ylo, post-shock conditions inflow for x<1/6, reflective wall for x>1/6
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+real(i,WP)*dx
               if (x.lt.x0_shock) then ! post-shock conditions
                 p(i,j,k,1)=U2
               else ! reflective wall
                 p(i,j,k,1)=p(i,-j-1,k,1)
               end if
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x.lt.x0_shock) then ! post-shock conditions
                 p(i,j,k,1)=V2
               else ! reflective wall 
                 p(i,0,k,1)=0
                 p(i,j,k,1)=-p(i,-j,k,1)
               end if
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               if (x.lt.x0_shock) then ! post-shock conditions
                 p(i,j,k,1)=rho2
                 p(i,j,k,2)=rho2*U2
                 p(i,j,k,3)=rho2*V2
                 p(i,j,k,4)=0.0_WP
                 p(i,j,k,5)=rho2*fluid%get_e_from_p_rho(p=P2,rho=rho2,y=[1.0_WP])
               else ! reflective wall
                 p(i,j,k,1)=p(i,-j-1,k,1)
                 p(i,j,k,2)=p(i,-j-1,k,2)
                 p(i,j,k,3)=-p(i,-j-1,k,3)
                 p(i,j,k,4)=0.0_WP
                 p(i,j,k,5)=p(i,-j-1,k,5)
               end if
            end do; end do; end do
         end select
     case (4) ! Yhi, time-dependent inflow condition (post- and pre-shock)
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+real(i,WP)*dx
               y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
               Hshock_t=Hshock(x,y,theta_shock,x0_shock+Ushock*time,dx)
               U=U1+(U2-U1)*Hshock_t
               p(i,j,k,1)=U
            end do; end do; end do
          case ('V')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               y=solver%amr%ylo+real(j,WP)*dy
               Hshock_t=Hshock(x,y,theta_shock,x0_shock+Ushock*time,dx)
               V=V1+(V2-V1)*Hshock_t
               p(i,j,k,1)=V
            end do; end do; end do
          case ('W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
               y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
               Hshock_t=Hshock(x,y,theta_shock,x0_shock+Ushock*time,dx)
               pres=P1+(P2-P1)*Hshock_t
               U=U1+(U2-U1)*Hshock_t
               V=V1+(V2-V1)*Hshock_t
               rho=rho1+(rho2-rho1)*Hshock_t
               IE=fluid%get_e_from_p_rho(p=pres,rho=rho,y=[1.0_WP])
               p(i,j,k,1)=rho
               p(i,j,k,2)=rho*U
               p(i,j,k,3)=rho*V
               p(i,j,k,4)=0.0_WP
               p(i,j,k,5)=rho*IE
            end do; end do; end do
         end select
      end select
   end subroutine dirichlet

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
         real(WP) :: Cv,Gamma,Un2
         call param_read('Mach number',M1)
         call param_read('Shock angle',theta_shock)
         theta_shock=4.0_WP*atan(1.0_WP)/180.0_WP*theta_shock ! deg to rad
         x0_shock=1.0_WP/6.0_WP
         Gamma=1.4_WP
         Cv=1.0_WP/(Gamma-1.0_WP) ! such that gas constant R = 1, i.e. T = p/rho
         call fluid%initialize(gamma=Gamma,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         ! Pre-shock conditions
         U1=0.0_WP; V1=0.0_WP
         rho1=1.4_WP; P1=1.0_WP ! such that rho/P=1.4, c=1
         Wshock=M1*fluid%get_c_from_p_rho(p=P1,rho=rho1,y=[1.0_WP])
         Ushock=Wshock/sin(theta_shock)
         ! Post-shock conditions from Rankine-Hugoniot
         rho2=rho1*(Gamma+1)*M1**2/(2+(Gamma-1)*M1**2)
         P2=P1*(2*Gamma*M1**2-(Gamma-1))/(Gamma+1)
         M2=sqrt(((Gamma-1)*M1**2+2)/(2*Gamma*M1**2-(Gamma-1)))
         ! Post-shock normal velocity relative to fixed frame
         Un2=Wshock-M2*fluid%get_c_from_p_rho(p=P2,rho=rho2,y=[1.0_WP])
         U2=Un2*sin(theta_shock)
         V2=-Un2*cos(theta_shock)
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
          amr%name='channel'
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
      end block create_amrgrid

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax,default=0.2_WP)
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
         ! Set user-defined Dirichlet inflow at X_lo and zero-gradient outflow at X_hi
         fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,1)=amrex_bc_ext_dir; fs%U%hi_bc(1,1)=amrex_bc_foextrap
         fs%V%lo_bc(1,1)=amrex_bc_ext_dir; fs%V%hi_bc(1,1)=amrex_bc_foextrap
         fs%W%lo_bc(1,1)=amrex_bc_ext_dir; fs%W%hi_bc(1,1)=amrex_bc_foextrap
         ! Set user-defined Dirichlet BCs at Y_lo and Y_hi
         fs%Q%lo_bc(2,:)=amrex_bc_ext_dir; fs%Q%hi_bc(2,:)=amrex_bc_ext_dir
         fs%U%lo_bc(2,1)=amrex_bc_ext_dir; fs%U%hi_bc(2,1)=amrex_bc_ext_dir
         fs%V%lo_bc(2,1)=amrex_bc_ext_dir; fs%V%hi_bc(2,1)=amrex_bc_ext_dir
         fs%W%lo_bc(2,1)=amrex_bc_ext_dir; fs%W%hi_bc(2,1)=amrex_bc_ext_dir
         ! Point to user-defined boundary condition
         fs%user_bc=>dirichlet
      end block create_solver

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
         call viz%initialize(amr,'DoubleMach',use_hdf5=.false.)
         call viz%add_scalar(fs%Q,1,'RHO')
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
