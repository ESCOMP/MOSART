module mosart_histflds

   ! Module containing initialization of history fields and files
   ! This is the module that the user must modify in order to add new
   ! history fields or modify defaults associated with existing history
   ! fields.

   use shr_kind_mod    , only : r8 => shr_kind_r8
   use mosart_histfile , only : mosart_hist_addfld, mosart_hist_printflds
   use mosart_data     , only : ctl, Trunoff

   implicit none
   private

   public :: mosart_histflds_init
   public :: mosart_histflds_set

  type, public ::  hist_pointer_type
     real(r8), pointer :: data(:) => null()
  end type hist_pointer_type

  type(hist_pointer_type), allocatable :: h_runofflnd(:)
  type(hist_pointer_type), allocatable :: h_runoffocn(:)
  type(hist_pointer_type), allocatable :: h_runofftot(:)
  type(hist_pointer_type), allocatable :: h_direct(:)
  type(hist_pointer_type), allocatable :: h_direct_glc(:)
  type(hist_pointer_type), allocatable :: h_dvolrdtlnd(:)
  type(hist_pointer_type), allocatable :: h_dvolrdtocn(:)
  type(hist_pointer_type), allocatable :: h_volr(:)
  type(hist_pointer_type), allocatable :: h_qsur(:)
  type(hist_pointer_type), allocatable :: h_qsub(:)
  type(hist_pointer_type), allocatable :: h_qgwl(:)

  real(r8), pointer :: h_volr_mch(:)
  real(r8), pointer :: h_qglc_liq_input(:)
  real(r8), pointer :: h_qglc_ice_input(:)

!------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine mosart_histflds_init(begr, endr, ntracers)

      ! Arguments
      integer, intent(in) :: begr
      integer, intent(in) :: endr
      integer, intent(in) :: ntracers

      ! Local variables
      integer :: nt

      !-------------------------------------------------------
      ! Allocate memory for module variables
      !-------------------------------------------------------

      allocate(h_runofflnd(ntracers))
      allocate(h_runoffocn(ntracers))
      allocate(h_runofftot(ntracers))
      allocate(h_direct(ntracers))
      allocate(h_dvolrdtlnd(ntracers))
      allocate(h_dvolrdtocn(ntracers))
      allocate(h_volr(ntracers))
      allocate(h_qsur(ntracers))
      allocate(h_qsub(ntracers))
      allocate(h_qgwl(ntracers))
      allocate(h_direct_glc(2))

      do nt = 1,ntracers
         allocate(h_runofflnd(nt)%data(begr:endr))
         allocate(h_runoffocn(nt)%data(begr:endr))
         allocate(h_runofftot(nt)%data(begr:endr))
         allocate(h_direct(nt)%data(begr:endr))
         allocate(h_dvolrdtlnd(nt)%data(begr:endr))
         allocate(h_dvolrdtocn(nt)%data(begr:endr))
         allocate(h_volr(nt)%data(begr:endr))
         allocate(h_qsur(nt)%data(begr:endr))
         allocate(h_qsub(nt)%data(begr:endr))
         allocate(h_qgwl(nt)%data(begr:endr))
      end do
      allocate(h_direct_glc(ctl%nt_liq)%data(begr:endr))
      allocate(h_direct_glc(ctl%nt_ice)%data(begr:endr))

      allocate(h_volr_mch(begr:endr))
      allocate(h_qglc_liq_input(begr:endr))
      allocate(h_qglc_ice_input(begr:endr))

      !-------------------------------------------------------
      ! Build master field list of all possible fields in a history file.
      ! Each field has associated with it a ``long\_name'' netcdf attribute that
      ! describes what the field is, and a ``units'' attribute. A subroutine is
      ! called to add each field to the masterlist.
      !-------------------------------------------------------

      do nt = 1,ctl%ntracers

         call mosart_hist_addfld (fname='RIVER_DISCHARGE_OVER_LAND'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART river basin flow: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runofflnd(nt)%data, default='active')

         call mosart_hist_addfld (fname='RIVER_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART river discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runoffocn(nt)%data, default='active')

         call mosart_hist_addfld (fname='TOTAL_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s', &
              avgflag='A', long_name='MOSART total discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_runofftot(nt)%data, default='active')

         call mosart_hist_addfld (fname='DIRECT_DISCHARGE_TO_OCEAN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s', &
              avgflag='A', long_name='MOSART direct discharge into ocean: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_direct(nt)%data, default='active')

         call mosart_hist_addfld (fname='DIRECT_DISCHARGE_TO_OCEAN_GLC'//'_'//trim(ctl%tracer_names(nt)), units='m3/s', &
              avgflag='A', long_name='MOSART direct discharge into ocean from glc: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_direct_glc(nt)%data, default='active')

         call mosart_hist_addfld (fname='STORAGE'//'_'//trim(ctl%tracer_names(nt)), units='m3',  &
              avgflag='A', long_name='MOSART storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_volr(nt)%data, default='inactive')

         call mosart_hist_addfld (fname='DVOLRDT_LND'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART land change in storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_dvolrdtlnd(nt)%data, default='inactive')

         call mosart_hist_addfld (fname='DVOLRDT_OCN'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART ocean change of storage: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_dvolrdtocn(nt)%data, default='inactive')

         call mosart_hist_addfld (fname='QSUR'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART input surface runoff: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_qsur(nt)%data, default='inactive')

         call mosart_hist_addfld (fname='QSUB'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART input subsurface runoff: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_qsub(nt)%data, default='inactive')

         call mosart_hist_addfld (fname='QGWL'//'_'//trim(ctl%tracer_names(nt)), units='m3/s',  &
              avgflag='A', long_name='MOSART input GWL runoff: '//trim(ctl%tracer_names(nt)), &
              ptr_rof=h_qgwl(nt)%data, default='inactive')
      end do

      call mosart_hist_addfld (fname='STORAGE_MCH', units='m3',  &
           avgflag='A', long_name='MOSART main channelstorage', &
           ptr_rof=h_volr_mch, default='inactive')

      call mosart_hist_addfld (fname='QIRRIG_FROM_COUPLER', units='m3/s',  &
           avgflag='A', long_name='Amount of water used for irrigation (total flux received from coupler)', &
           ptr_rof=ctl%qirrig, default='inactive')

      call mosart_hist_addfld (fname='QIRRIG_ACTUAL', units='m3/s',  &
           avgflag='A', long_name='Actual irrigation (if limited by river storage)', &
           ptr_rof=ctl%qirrig_actual, default='inactive')

      call mosart_hist_addfld (fname='QGLC_LIQ_INPUT', units='m3',  &
           avgflag='A', long_name='liquid runoff from glc input', &
           ptr_rof=h_qglc_liq_input, default='active')

      call mosart_hist_addfld (fname='QGLC_ICE_INPUT', units='m3',  &
           avgflag='A', long_name='ice runoff from glc input', &
           ptr_rof=h_qglc_ice_input, default='active')

      ! print masterlist of history fields
      call mosart_hist_printflds()

   end subroutine mosart_histflds_init

   !-----------------------------------------------------------------------

   subroutine mosart_histflds_set(ntracers)

      !-----------------------------------------------------------------------
      ! Set mosart history fields as 1d pointer arrays
      !-----------------------------------------------------------------------

      ! Arguments
      integer, intent(in) :: ntracers

      ! Local variables
      integer :: nt
      integer :: nt_liq, nt_ice

      nt_liq = ctl%nt_liq
      nt_ice = ctl%nt_ice

      do nt = 1,ntracers
         h_runofflnd(nt)%data(:)  = ctl%runofflnd(:,nt)
         h_runoffocn(nt)%data(:)  = ctl%runoffocn(:,nt)
         h_runofftot(nt)%data(:)  = ctl%runofftot(:,nt)
         h_direct(nt)%data(:)     = ctl%direct(:,nt)
         h_dvolrdtlnd(nt)%data(:) = ctl%dvolrdtlnd(:,nt)
         h_dvolrdtocn(nt)%data(:) = ctl%dvolrdtocn(:,nt)
         h_qsub(nt)%data(:)       = ctl%qsub(:,nt)
         h_qsur(nt)%data(:)       = ctl%qsur(:,nt)
         h_qgwl(nt)%data(:)       = ctl%qgwl(:,nt)
      end do
      h_volr_mch(:) = Trunoff%wr(:,1)
      h_qglc_liq_input(:) = ctl%qglc_liq(:)
      h_qglc_ice_input(:) = ctl%qglc_ice(:)
      h_direct_glc(nt_liq)%data(:) = ctl%direct_glc(:,nt_liq)
      h_direct_glc(nt_ice)%data(:) = ctl%direct_glc(:,nt_ice)

   end subroutine mosart_histflds_set

end module mosart_histflds
