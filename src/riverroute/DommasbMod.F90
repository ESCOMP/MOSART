MODULE DommasbMod
  !Description: core code of Dissolved Organic Matter mass balance utilizing river routing models
  !Developed by Dev Narayanappa 05/03/2021
  !This module is currently made interact with MOSART routing model
  ! USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_abort
  use RunoffMod     , only : TRunoff, Tdom
  use RtmVar        , only : iulog

  implicit none

  public hillslopeRoutingDOM
  public subnetworkRoutingDOM
  public mainchannelRoutingDOM
  !--------------------------------------------------------------------

  ! ! PUBLIC MEMBER FUNCTIONS:
  contains

  !----------------------------------------------------------------------
    subroutine hillslopeRoutingDOM(iunit,nt,ntdom,theDeltaT,Darea,Dfrac)
      ! ! DESCRIPTION: solve the ODEs with Euler algorithm for hillslope routing
      implicit none
      integer,  intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT, Darea, Dfrac
      ! assume no chemical reaction in the water hence sink term is zero implies domsur = domR*flow
      ! ehout is negative
      Tdom%domH(iunit,ntdom) = Tdom%domH(iunit,ntdom) + (TRunoff%ehout(iunit,nt) * Darea * Dfrac * Tdom%domH(iunit,ntdom) + Tdom%domsur(iunit,ntdom)) * theDeltaT/(TRunoff%wh(iunit,nt)*Darea*Dfrac)
    end subroutine hillslopeRoutingDOM

    subroutine subnetworkRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! solve the ODEs with Euler algorithm for subnetwork routing
      ! etin is positive and etout is negative
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      Tdom%domT(iunit,ntdom) = Tdom%domT(iunit,ntdom) + (TRunoff%etin(iunit,nt) * Tdom%domH(iunit,ntdom) + TRunoff%etout(iunit,nt) * Tdom%domT(iunit,ntdom)) * theDeltaT/TRunoff%wt(iunit,nt)
    end subroutine subnetworkRoutingDOM

    subroutine mainchannelRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! solve the ODE with Euler algorithm for main-channel routing
      ! erout is negative, while erlateral and erin are positive
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      Tdom%domR(iunit,ntdom) = Tdom%domR(iunit,ntdom) + (TRunoff%erlateral(iunit,nt)*Tdom%domT(iunit,ntdom) + TRunoff%erin(iunit,nt)*Tdom%domRUp(iunit,ntdom) + TRunoff%erout(iunit,nt)*Tdom%domR(iunit,ntdom))*theDeltaT/TRunoff%wr(iunit,nt)
    end subroutine mainchannelRoutingDOM
!-------------------------------------------------------------------------
end MODULE DommasbMod
