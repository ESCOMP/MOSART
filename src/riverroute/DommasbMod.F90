MODULE DommasbMod
  !Description: core code of Dissolved Organic Matter mass balance utilizing river routing models
  !Developed by Dev Narayanappa 05/03/2021
  !This module is currently made interact with MOSART routing model
  ! USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_abort
  use RunoffMod     , only : TRunoff, Tdom

  implicit none

  public hillslopeRoutingDOM
  public subnetworkRoutingDOM
  public mainchannelRoutingDOM
  !--------------------------------------------------------------------

  ! ! PUBLIC MEMBER FUNCTIONS:
  contains

  !----------------------------------------------------------------------
    subroutine hillslopeRoutingDOM(iunit,nt,theDeltaT)
      ! ! DESCRIPTION: solve the ODEs with Euler algorithm for hillslope routing
      implicit none
      integer,  intent(in) :: iunit, nt
      real(r8), intent(in) :: theDeltaT
      ! assume no chemical reaction in the water hence sink term is zero implies domH ~= domHout
      Tdom%domH(iunit,nt) = Tdom%domH(iunit,nt) + (-TRunoff%ehout(iunit,nt) * Tdom%domH(iunit,nt) + Tdom%domSource(iunit,nt)) * theDeltaT/TRunoff%wh(iunit,nt)
    end subroutine hillslopeRoutingDOM

    subroutine subnetworkRoutingDOM(iunit,nt,theDeltaT)
      ! solve the ODEs with Euler algorithm for subnetwork routing
      implicit none
      integer, intent(in) :: iunit,nt
      real(r8), intent(in) :: theDeltaT
      Tdom%domT(iunit,nt) = Tdom%domT(iunit,nt) + (TRunoff%etin(iunit,nt) * Tdom%domH(iunit,nt) - TRunoff%etout(iunit,nt) * Tdom%domT(iunit,nt)) * theDeltaT/TRunoff%wt(iunit,nt)
    end subroutine subnetworkRoutingDOM

    subroutine mainchannelRoutingDOM(iunit,nt,theDeltaT)
      ! solve the ODE with Euler algorithm for main-channel routing
      implicit none
      integer, intent(in) :: iunit, nt
      real(r8), intent(in) :: theDeltaT
      real(r8)  :: mainchinT, mainchinUp
      real(r8)  :: mainchout
      mainchinT  = TRunoff%etout(iunit,nt) - TRunoff%erlateral(iunit,nt) !input to main channel from Tributaries
      mainchinUp = TRunoff%eroutUp(iunit,nt)   !inflow to main channel from Upstream grid cells of main channel
      mainchout = TRunoff%flow(iunit,nt)     ! flow to the outlet of the reach in the main channel
      Tdom%domR(iunit,nt) = Tdom%domR(iunit,nt) + ( (mainchinT*Tdom%domT(iunit,nt) + mainchinUp*Tdom%domRout(iunit,nt)) - TRunoff%flow(iunit,nt) * Tdom%domR(iunit,nt))*theDeltaT/TRunoff%wr(iunit,nt)
    end subroutine mainchannelRoutingDOM
!-------------------------------------------------------------------------
end MODULE DommasbMod
