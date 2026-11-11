!> Simulation of SHOCK-BUBBLE INTERACTION (Quirk & Karni, 1996) with AMR
module simulation
   use precision,          only: WP
   use amrgrid_class,      only: amrgrid
   use amrcomp_class,      only: amrcomp
   use amrviz_class,       only: amrviz
   use amrdata_class,      only: amrdata
   use timetracker_class,  only: timetracker
   use event_class,        only: event
   use monitor_class,      only: monitor
   use igmix_class,        only: igmix
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
   type(igmix), target :: fluid

   !> Flow parameters
   real(WP) :: rho_pre ,P_pre          !< Pre-shock state
   real(WP) :: rho_post,P_post,U_post  !< Post-shock state
   real(WP) :: x0_shock                !< Initial shock position
   real(WP) :: xb,Db,rhob              !< Bubble properties

contains

   !> Heaviside function for shock initialization
   !> 0 pre-shock, 1 post-shock
   real(WP) function shock_heaviside(x)
      real(WP), intent(in) :: x
      if (x.lt.x0_shock) then
         shock_heaviside=0.0_WP
      else
         shock_heaviside=1.0_WP
      end if
   end function shock_heaviside

   !> Levelset function for 2D bubble
   real(WP) function bubble_levelset(xyz)
      real(WP), dimension(3), intent(in) :: xyz
      bubble_levelset=sqrt((xyz(1)-xb)**2+xyz(2)**2)-0.5_WP*Db
   end function bubble_levelset

   !> Set zero physical viscosities (inviscid)
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
      real(WP) :: dx,dy,x,y,H,rho,P,U,IE,Yb
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.true.)
      do while (mfi%next())
         ! Get pointers to Q
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         ! Get tilebox with overlap
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get physical coordinates
            x=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            y=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
            ! Find species fraction
            Yb=0.0_WP
            if (bubble_levelset([x,y,0.0_WP]).lt.0.0_WP) then
               Yb=1.0_WP
            end if
            ! Set initial conditions
            H=shock_heaviside(x)
            U=U_post*H
            rho=rho_pre+(rho_post-rho_pre)*H
            P  =P_pre  +(P_post  -P_pre  )*H
            if (Yb.gt.0.0_WP) then
               rho=rhob; P=P_pre; U=0.0_WP
            end if
            IE =fluid%get_e_from_p_rho(p=P,rho=rho,y=[Yb,1.0_WP-Yb])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=rho*U
            pQ(i,j,k,3)=0.0_WP
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho*IE
            pQ(i,j,k,6)=rho*Yb
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine fs_init

   !> Solver BC callback - post-shock inflow at x-high
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
       case (2) ! x-high: post-shock inflow
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=U_post
            end do; end do; end do
          case ('V','W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pBC(i,j,k,1)=rho_post
               pBC(i,j,k,2)=rho_post*U_post
               pBC(i,j,k,3)=0.0_WP
               pBC(i,j,k,4)=0.0_WP
               pBC(i,j,k,5)=rho_post*fluid%get_e_from_p_rho(p=P_post,rho=rho_post,y=[0.0_WP,1.0_WP])
               pBC(i,j,k,6)=0.0_WP
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
         real(WP) :: gamma,R,Cv,Ms,Ws,M_post,gammab,Rb,Cvb
         call param_read('Shock location',x0_shock)
         call param_read('Shock Mach number',Ms)
         call param_read('Density pre-shock',rho_pre)
         call param_read('Pressure pre-shock',P_pre)
         call param_read('Specific heat ratio',gamma)
         call param_read('Gas constant',R)
         call param_read('Bubble specific heat ratio',gammab)
         call param_read('Bubble gas constant',Rb)
         call param_read('Bubble diameter',Db)
         ! Initialize shock tube and bubble fluids
         Cv=R/(gamma-1.0_WP)
         Cvb=Rb/(gammab-1.0_WP)
         call fluid%initialize(gamma=[gammab,gamma],cv=[Cvb,Cv], &
                               q=[0.0_WP,0.0_WP],qp=[0.0_WP,0.0_WP], &
                               species_names=[character(len=6)::'bubble','tube'], name='mix')
         ! Shock speed in lab frame
         Ws=Ms*fluid%get_c_from_p_rho(p=P_pre,rho=rho_pre,y=[0.0_WP,1.0_WP])
         ! Post-shock state from Rankine-Hugoniot
         rho_post=rho_pre*(gamma+1.0_WP)*Ms**2/(2.0_WP+(gamma-1.0_WP)*Ms**2)
         P_post=P_pre*(2.0_WP*gamma*Ms**2-(gamma-1.0_WP))/(gamma+1.0_WP)
         M_post=sqrt(((gamma-1.0_WP)*Ms**2+2.0_WP)/(2.0_WP*gamma*Ms**2-(gamma-1.0_WP)))
         U_post=-Ws+M_post*fluid%get_c_from_p_rho(p=P_post,rho=rho_post,y=[0.0_WP,1.0_WP])
         ! Initial bubble density
         rhob=rho_pre*R/Rb ! for thermal equilibrium
      end block init_flow

      ! Initialize AMR grid
      init_amrgrid: block
         amr%name='grid'
         call param_read('Length',amr%xhi)
         call param_read('Height',amr%yhi)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Max levels',amr%maxlvl)
         amr%xlo=0.0_WP; amr%ylo=0.0_WP
         ! Enable quasi-2D: set z extent to match one cell at finest level
         amr%nz=1
         amr%zlo=-0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         amr%zhi=+0.5_WP*(amr%xhi-amr%xlo)/real(amr%nx*2**amr%maxlvl,WP)
         ! Set periodicity
         amr%xper=.false.; amr%yper=.false.; amr%zper=.true.
         ! Place bubble at center of domain
         xb=amr%xhi*0.5_WP
         call amr%initialize()
      end block init_amrgrid

      ! Initialize time tracker
      init_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         call param_read('Subiterations',time%itmax,default=2)
         time%dt=time%dtmax
      end block init_timetracker

      ! Initialize compressible solver
      init_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrex_amr_module, only: amrex_bc_reflect_odd,amrex_bc_reflect_even
         use amrdata_class,    only: interp_face_lin
         ! Assign material and create flow solver
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set initial conditions
         fs%user_init=>fs_init
         ! Set reflective wall at x-low
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap
         fs%Q%lo_bc(1,2)=amrex_bc_reflect_odd
         fs%Q%lo_bc(1,3)=amrex_bc_reflect_even
         fs%U%lo_bc(1,1)=amrex_bc_reflect_odd
         fs%V%lo_bc(1,1)=amrex_bc_reflect_even
         fs%W%lo_bc(1,1)=amrex_bc_foextrap
         ! Set post-shock inflow at x-high
         fs%Q%hi_bc(1,:)=amrex_bc_ext_dir
         fs%U%hi_bc(1,1)=amrex_bc_ext_dir
         fs%V%hi_bc(1,1)=amrex_bc_ext_dir
         fs%W%hi_bc(1,1)=amrex_bc_ext_dir
         ! Set symmetry at y-low
         fs%Q%lo_bc(2,:)=amrex_bc_foextrap
         fs%Q%lo_bc(2,2)=amrex_bc_reflect_even
         fs%Q%lo_bc(2,3)=amrex_bc_reflect_odd
         fs%U%lo_bc(2,1)=amrex_bc_reflect_even
         fs%V%lo_bc(2,1)=amrex_bc_reflect_odd
         fs%W%lo_bc(2,1)=amrex_bc_foextrap
         ! Set reflective wall at y-high
         fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%Q%hi_bc(2,2)=amrex_bc_reflect_even
         fs%Q%hi_bc(2,3)=amrex_bc_reflect_odd
         fs%U%hi_bc(2,1)=amrex_bc_reflect_even
         fs%V%hi_bc(2,1)=amrex_bc_reflect_odd
         fs%W%hi_bc(2,1)=amrex_bc_foextrap
         ! Point to user-defined boundary conditions
         fs%user_bc=>fs_bc
      end block init_solver

      ! Initialize workspaces
      init_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none)
         call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none)
         call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none)
         call Mach%register()
      end block init_workspace

      ! Initialize regridding
      init_regridding: block
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>fs_tag
         call param_read('Tagging criterion',drho_tag)
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
         call viz%add_scalar(fs%Y,1,'Yb')
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
