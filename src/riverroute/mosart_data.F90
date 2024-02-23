module mosart_data

   use mosart_control_type,      only : control_type
   use mosart_tctl_type,         only : tctl_type
   use mosart_tspatialunit_type, only : tspatialunit_type
   use mosart_tstatusflux_type,  only : tstatusflux_type
   use mosart_tparameter_type,   only : tparameter_type

   implicit none
   private

   ! Derived types
   type(Tctl_type),         public :: Tctl
   type(Tspatialunit_type), public :: TUnit
   type(TstatusFlux_type),  public :: TRunoff
   type(Tparameter_type),   public :: TPara
   type(control_type),      public :: ctl

end module mosart_data
