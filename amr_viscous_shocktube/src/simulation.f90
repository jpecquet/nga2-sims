!> Simulation of VISCOUS SHOCKTUBE (Daru & Tenaud, 2009) with AMR
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
   real(WP) :: x0_diaphragm
   real(WP) :: rho_left,rho_right,P_left,P_right
   real(WP) :: Cp,Pr,mu

contains

   !> Heaviside function for shock tube initialization
   !> 0 left of diaphragm, 1 right of diaphragm
   real(WP) function diaphragm_heaviside(x)
      real(WP), intent(in) :: x
      if (x.lt.x0_diaphragm) then
         diaphragm_heaviside=0.0_WP
      else
         diaphragm_heaviside=1.0_WP
      end if
   end function diaphragm_heaviside

   !> Set constant physical viscosities (Stokes hypothesis)
   subroutine get_viscosities()
      implicit none
      call fs%visc%setval(val=mu)
      call fs%beta%setval(val=0.0_WP)
      call fs%diff%setval(val=mu*Cp/Pr)
   end subroutine get_viscosities

   !> Update viscosities and derived output fields from current primitive variables
   subroutine update_derived()
      implicit none
      ! Compute viscosities
      call get_viscosities()
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
      real(WP) :: dx,x,H,rho,P,IE
      integer :: i,j,k
      ! Get mesh size
      dx=solver%amr%dx(lvl)
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
            ! Set initial conditions
            H=diaphragm_heaviside(x)
            rho=rho_left+(rho_right-rho_left)*H
            P  =P_left  +(P_right  -P_left  )*H
            IE =fluid%get_e_from_p_rho(p=P,rho=rho,y=[1.0_WP])
            ! Set conserved variable Q
            pQ(i,j,k,1)=rho
            pQ(i,j,k,2)=0.0_WP
            pQ(i,j,k,3)=0.0_WP
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=rho*IE
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine fs_init

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
         real(WP) :: gamma,Cv,Re
         call param_read('Diaphragm location',x0_diaphragm)
         call param_read('Density left',rho_left)
         call param_read('Density right',rho_right)
         call param_read('Specific heat ratio',gamma)
         call param_read('Prandtl number',Pr)
         call param_read('Reynolds number',Re)
         Cv=1.0_WP/(gamma-1.0_WP)
         call fluid%initialize(gamma=gamma,cv=Cv,q=0.0_WP,qp=0.0_WP,name='fluid')
         mu=1.0_WP/Re
         Cp=gamma*Cv
         P_left =rho_left/gamma
         P_right=rho_right/gamma
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
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_odd,amrex_bc_reflect_even
         use amrdata_class,    only: interp_face_lin
         ! Assign material and create flow solver
         fs%mat=>fluid; call fs%initialize(amr=amr)
         ! Use face-linear interp (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Set initial conditions
         fs%user_init=>fs_init
         ! Set no-slip wall at x-low
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap
         fs%Q%lo_bc(1,2)=amrex_bc_reflect_odd
         fs%Q%lo_bc(1,3)=amrex_bc_reflect_odd
         fs%U%lo_bc(1,1)=amrex_bc_reflect_odd
         fs%V%lo_bc(1,1)=amrex_bc_reflect_odd
         fs%W%lo_bc(1,1)=amrex_bc_foextrap
         ! Set no-slip wall at x-high
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%Q%hi_bc(1,2)=amrex_bc_reflect_odd
         fs%Q%hi_bc(1,3)=amrex_bc_reflect_odd
         fs%U%hi_bc(1,1)=amrex_bc_reflect_odd
         fs%V%hi_bc(1,1)=amrex_bc_reflect_odd
         fs%W%hi_bc(1,1)=amrex_bc_foextrap
         ! Set no-slip wall at y-low
         fs%Q%lo_bc(2,:)=amrex_bc_foextrap
         fs%Q%lo_bc(2,2)=amrex_bc_reflect_odd
         fs%Q%lo_bc(2,3)=amrex_bc_reflect_odd
         fs%U%lo_bc(2,1)=amrex_bc_reflect_odd
         fs%V%lo_bc(2,1)=amrex_bc_reflect_odd
         fs%W%lo_bc(2,1)=amrex_bc_foextrap
         ! Set symmetry at y-hi
         fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%Q%hi_bc(2,2)=amrex_bc_reflect_even
         fs%Q%hi_bc(2,3)=amrex_bc_reflect_odd
         fs%U%hi_bc(2,1)=amrex_bc_reflect_even
         fs%V%hi_bc(2,1)=amrex_bc_reflect_odd
         fs%W%hi_bc(2,1)=amrex_bc_foextrap
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
