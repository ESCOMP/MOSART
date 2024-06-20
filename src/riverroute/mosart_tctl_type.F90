module mosart_tctl_type

   use shr_kind_mod, only : r8 => shr_kind_r8, CL => SHR_KIND_CL

   implicit none
   private

   type Tctl_type
      real(r8) :: DeltaT        ! Time step in seconds
      integer  :: DLevelH2R     ! The base number of channel routing sub-time-steps within one hillslope routing step.
                                ! Usually channel routing requires small time steps than hillslope routing.
      integer  :: DLevelR       ! The number of channel routing sub-time-steps at a higher level within one channel routing step at a lower level.
      integer  :: RoutingMethod ! Flag for routing methods. 1 --> variable storage method from SWAT model
   contains
      procedure :: Init
   end type Tctl_type
   public :: Tctl_type

contains

   subroutine Init(this)
      class(Tctl_type) :: this

      this%RoutingMethod = 1
      this%DLevelH2R = 5
      this%DLevelR = 3

   end subroutine Init

end module mosart_tctl_type
