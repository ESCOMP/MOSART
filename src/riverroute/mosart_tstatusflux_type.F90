module mosart_tstatusflux_type

   ! status and flux variables

   use shr_kind_mod, only : r8 => shr_kind_r8, CL => SHR_KIND_CL

   implicit none
   private

   public :: TstatusFlux_type
   type TstatusFlux_type
      ! hillsloope
      !! states
      real(r8), pointer :: wh(:,:)      ! storage of surface water, [m]
      real(r8), pointer :: dwh(:,:)     ! change of water storage, [m/s]
      real(r8), pointer :: yh(:,:)      ! depth of surface water, [m]
      real(r8), pointer :: wsat(:,:)    ! storage of surface water within saturated area at hillslope [m]
      real(r8), pointer :: wunsat(:,:)  ! storage of surface water within unsaturated area at hillslope [m]
      real(r8), pointer :: qhorton(:,:) ! Infiltration excess runoff generated from hillslope, [m/s] NOT_USED
      real(r8), pointer :: qdunne(:,:)  ! Saturation excess runoff generated from hillslope, [m/s] NOT_USED
      real(r8), pointer :: qsur(:,:)    ! Surface runoff generated from hillslope, [m/s]
      real(r8), pointer :: qsub(:,:)    ! Subsurface runoff generated from hillslope, [m/s]
      real(r8), pointer :: qgwl(:,:)    ! gwl runoff term from glacier, wetlands and lakes, [m/s]
      !! fluxes
      real(r8), pointer :: ehout(:,:)   ! overland flow from hillslope into the sub-channel, [m/s]
      real(r8), pointer :: asat(:,:)    ! saturated area fraction from hillslope, [-]
      real(r8), pointer :: esat(:,:)    ! evaporation from saturated area fraction at hillslope, [m/s]

      ! subnetwork channel
      !! states
      real(r8), pointer :: tarea(:,:)   ! area of channel water surface, [m2]
      real(r8), pointer :: wt(:,:)      ! storage of surface water, [m3]
      real(r8), pointer :: dwt(:,:)     ! change of water storage, [m3]
      real(r8), pointer :: yt(:,:)      ! water depth, [m]
      real(r8), pointer :: mt(:,:)      ! cross section area, [m2]
      real(r8), pointer :: rt(:,:)      ! hydraulic radii, [m]
      real(r8), pointer :: pt(:,:)      ! wetness perimeter, [m]
      real(r8), pointer :: vt(:,:)      ! flow velocity, [m/s]
      real(r8), pointer :: tt(:,:)      ! mean travel time of the water within the channel, [s] NOT_USED
      !! fluxes
      real(r8), pointer :: etin(:,:)    ! lateral inflow from hillslope, including surface and subsurface runoff generation components, [m3/s]
      real(r8), pointer :: etout(:,:)   ! discharge from sub-network into the main reach, [m3/s]

      ! main channel
      !! states
      real(r8), pointer :: rarea(:,:)   ! area of channel water surface, [m2]
      real(r8), pointer :: wr(:,:)      ! storage of surface water, [m3]
      real(r8), pointer :: dwr(:,:)     ! change of water storage, [m3]
      real(r8), pointer :: yr(:,:)      ! water depth. [m]
      real(r8), pointer :: mr(:,:)      ! cross section area, [m2]
      real(r8), pointer :: rr(:,:)      ! hydraulic radius, [m]
      real(r8), pointer :: pr(:,:)      ! wetness perimeter, [m]
      real(r8), pointer :: vr(:,:)      ! flow velocity, [m/s]
      real(r8), pointer :: tr(:,:)      ! mean travel time of the water within the channel, [s] NOT_USED
      !! exchange fluxes
      real(r8), pointer :: erlateral(:,:)   ! lateral flow from hillslope, including surface and subsurface runoff generation components, [m3/s]
      real(r8), pointer :: erin(:,:)        ! inflow from upstream links, [m3/s]
      real(r8), pointer :: erout(:,:)       ! outflow into downstream links, [m3/s]
      real(r8), pointer :: erout_prev(:,:)  ! outflow into downstream links from previous timestep, [m3/s]
      real(r8), pointer :: eroutUp(:,:)     ! outflow sum of upstream gridcells, instantaneous (m3/s)
      real(r8), pointer :: eroutUp_avg(:,:) ! outflow sum of upstream gridcells, average [m3/s]
      real(r8), pointer :: erlat_avg(:,:)   ! erlateral average [m3/s]
      real(r8), pointer :: flow(:,:)        ! streamflow from the outlet of the reach, [m3/s]
      real(r8), pointer :: erin1(:,:)       ! inflow from upstream links during previous step, used for Muskingum method, [m3/s] NOT_USED
      real(r8), pointer :: erin2(:,:)       ! inflow from upstream links during current step, used for Muskingum method, [m3/s] NOT_USED
      real(r8), pointer :: ergwl(:,:)       ! flux item for the adjustment of water balance residual in glacie, wetlands and lakes dynamics [m3/s] NOT_USED

      !! for Runge-Kutta algorithm NOT_USED
      real(r8), pointer :: wrtemp(:,:)  ! temporary storage item, for 4th order Runge-Kutta  algorithm;
      real(r8), pointer :: erintemp(:,:)
      real(r8), pointer :: erouttemp(:,:)
      real(r8), pointer :: k1(:,:)
      real(r8), pointer :: k2(:,:)
      real(r8), pointer :: k3(:,:)
      real(r8), pointer :: k4(:,:)
   contains
      procedure, public :: Init
   end type TstatusFlux_type

contains

   subroutine Init(this, begr, endr, ntracers)
      class(TstatusFlux_type) :: this
      integer, intent(in) :: begr, endr, ntracers

      ! Initialize water states and fluxes
      allocate (this%wh(begr:endr,ntracers))
      this%wh = 0._r8
      allocate (this%dwh(begr:endr,ntracers))
      this%dwh = 0._r8
      allocate (this%yh(begr:endr,ntracers))
      this%yh = 0._r8
      allocate (this%qsur(begr:endr,ntracers))
      this%qsur = 0._r8
      allocate (this%qsub(begr:endr,ntracers))
      this%qsub = 0._r8
      allocate (this%qgwl(begr:endr,ntracers))
      this%qgwl = 0._r8
      allocate (this%ehout(begr:endr,ntracers))
      this%ehout = 0._r8
      allocate (this%tarea(begr:endr,ntracers))
      this%tarea = 0._r8
      allocate (this%wt(begr:endr,ntracers))
      this%wt= 0._r8
      allocate (this%dwt(begr:endr,ntracers))
      this%dwt = 0._r8
      allocate (this%yt(begr:endr,ntracers))
      this%yt = 0._r8
      allocate (this%mt(begr:endr,ntracers))
      this%mt = 0._r8
      allocate (this%rt(begr:endr,ntracers))
      this%rt = 0._r8
      allocate (this%pt(begr:endr,ntracers))
      this%pt = 0._r8
      allocate (this%vt(begr:endr,ntracers))
      this%vt = 0._r8
      allocate (this%tt(begr:endr,ntracers))
      this%tt = 0._r8
      allocate (this%etin(begr:endr,ntracers))
      this%etin = 0._r8
      allocate (this%etout(begr:endr,ntracers))
      this%etout = 0._r8
      allocate (this%rarea(begr:endr,ntracers))
      this%rarea = 0._r8
      allocate (this%wr(begr:endr,ntracers))
      this%wr = 0._r8
      allocate (this%dwr(begr:endr,ntracers))
      this%dwr = 0._r8
      allocate (this%yr(begr:endr,ntracers))
      this%yr = 0._r8
      allocate (this%mr(begr:endr,ntracers))
      this%mr = 0._r8
      allocate (this%rr(begr:endr,ntracers))
      this%rr = 0._r8
      allocate (this%pr(begr:endr,ntracers))
      this%pr = 0._r8
      allocate (this%vr(begr:endr,ntracers))
      this%vr = 0._r8
      allocate (this%tr(begr:endr,ntracers))
      this%tr = 0._r8
      allocate (this%erlateral(begr:endr,ntracers))
      this%erlateral = 0._r8
      allocate (this%erin(begr:endr,ntracers))
      this%erin = 0._r8
      allocate (this%erout(begr:endr,ntracers))
      this%erout = 0._r8
      allocate (this%erout_prev(begr:endr,ntracers))
      this%erout_prev = 0._r8
      allocate (this%eroutUp(begr:endr,ntracers))
      this%eroutUp = 0._r8
      allocate (this%eroutUp_avg(begr:endr,ntracers))
      this%eroutUp_avg = 0._r8
      allocate (this%erlat_avg(begr:endr,ntracers))
      this%erlat_avg = 0._r8
      allocate (this%ergwl(begr:endr,ntracers))
      this%ergwl = 0._r8
      allocate (this%flow(begr:endr,ntracers))
      this%flow = 0._r8

   end subroutine Init

end module mosart_tstatusflux_type
