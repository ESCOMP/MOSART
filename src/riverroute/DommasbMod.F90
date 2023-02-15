MODULE DommasbMod
  !Description: core code of Dissolved Organic Matter mass balance utilizing river routing models
  !Developed by Marius Lambert 02-02-2023
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
    subroutine hillslopeRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! ! DESCRIPTION: solve the ODEs with Euler algorithm for hillslope routing
      implicit none
      integer,  intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      ! assume no chemical reaction in the water hence sink term is zero implies domsur = domR*flow
      ! ehout is negative
      Tdom%domH(iunit,ntdom) = (Tdom%domH(iunit,ntdom)*max(0._r8,TRunoff%wh(iunit,nt) - TRunoff%dwh(iunit,nt) * theDeltaT) + (TRunoff%ehout(iunit,nt) * Tdom%domH(iunit,ntdom) + Tdom%domsur(iunit,ntdom)) * theDeltaT)/TRunoff%wh(iunit,nt)
    end subroutine hillslopeRoutingDOM

    subroutine subnetworkRoutingDOM(iunit,nt,ntdom,theDeltaT,temp_ehout)
      ! solve the ODEs with Euler algorithm for subnetwork routing
      ! etin is positive and etout is negative
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT,temp_ehout
      Tdom%domT(iunit,ntdom) = (Tdom%domT(iunit,ntdom)*max(0._r8,TRunoff%wt(iunit,nt) - TRunoff%dwt(iunit,nt) * theDeltaT) + (Tdom%domsub(iunit,ntdom) + temp_ehout * Tdom%domH(iunit,ntdom) + TRunoff%etout(iunit,nt) * Tdom%domT(iunit,ntdom)) * theDeltaT)/TRunoff%wt(iunit,nt)
    end subroutine subnetworkRoutingDOM

    subroutine mainchannelRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! solve the ODE with Euler algorithm for main-channel routing
      ! erout is negative, while erlateral and erin are positive
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      Tdom%domR(iunit,ntdom) = (Tdom%domR(iunit,ntdom)*max(0._r8,TRunoff%wr(iunit,nt) - TRunoff%dwr(iunit,nt) * theDeltaT) + (TRunoff%erlateral(iunit,nt)*Tdom%domT(iunit,ntdom) + Tdom%domRUp(iunit,ntdom) + TRunoff%erout(iunit,nt)*Tdom%domR(iunit,ntdom))*theDeltaT)/TRunoff%wr(iunit,nt)
    end subroutine mainchannelRoutingDOM
!-------------------------------------------------------------------------
end MODULE DommasbMod
