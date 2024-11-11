module mosart_tparameter_type

   ! parameters to be calibrated. Ideally, these parameters are supposed to be uniform for one region

   use shr_kind_mod, only : r8 => shr_kind_r8, CL => SHR_KIND_CL

   implicit none
   private

   public :: Tparameter_type
   type Tparameter_type
      real(r8), pointer :: c_nr(:)       ! coefficient to adjust the manning's roughness of channels NOT_USED
      real(r8), pointer :: c_nh(:)       ! coefficient to adjust the manning's roughness of overland flow across hillslopes NOT_USED
      real(r8), pointer :: c_twid(:)     ! coefficient to adjust the width of sub-reach channel
   contains
      procedure, public :: Init
   end type Tparameter_type

contains

   subroutine Init(this, begr, endr)

      ! Arguments
      class(tparameter_type) :: this
      integer, intent(in) :: begr, endr

      ! Initialize TPara
      allocate (this%c_twid(begr:endr))
      this%c_twid = 1.0_r8

   end subroutine Init

end module mosart_tparameter_type
