MODULE DommasbMod
  !Description: core code of Dissolved Organic Matter mass balance utilizing river routing models
  !Developed by Marius Lambert 02-02-2023
  !This module is currently made interact with MOSART routing model
  ! USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
  use shr_sys_mod   , only : shr_sys_abort
  use RunoffMod     , only : TRunoff, Tdom, TUnit
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
      !domsur (kg/m2*s) ,domH  (kg/m2), ehout (m/s), domHout (kg/m2*s), qsur (m/s), wh (m)
      Tdom%domHout(iunit,ntdom) = -TRunoff%ehout(iunit,nt) * (Tdom%domH(iunit,ntdom) + Tdom%domsur(iunit,ntdom) * theDeltaT)/(TRunoff%wh(iunit,nt)-TRunoff%dwh(iunit,nt)*theDeltaT+TRunoff%qsur(iunit,nt)*theDeltaT)
      !we dont want a too high out
      Tdom%domHout(iunit,ntdom) = min(-TRunoff%ehout(iunit,nt) * 0.3_r8, Tdom%domHout(iunit,ntdom))
      !cannot be more be less than 0, lower boundary
      Tdom%domHout(iunit,ntdom) = max(0._r8,Tdom%domHout(iunit,ntdom))
      !cannot be more than available carbon, upper boundary
      Tdom%domHout(iunit,ntdom) = min((Tdom%domH(iunit,ntdom)+Tdom%domsur(iunit,ntdom)*theDeltaT)/theDeltaT,Tdom%domHout(iunit,ntdom))

      Tdom%domH(iunit,ntdom) = max(0._r8,Tdom%domH(iunit,ntdom) + (Tdom%domsur(iunit,ntdom) - Tdom%domHout(iunit,ntdom))* theDeltaT)
    end subroutine hillslopeRoutingDOM

    subroutine subnetworkRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! solve the ODEs with Euler algorithm for subnetwork routing
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      Tdom%domTout(iunit,ntdom) = -TRunoff%etout(iunit,nt) * (Tdom%domT(iunit,ntdom) + (Tdom%domsub(iunit,ntdom)+Tdom%domHout(iunit,ntdom)) * theDeltaT)/(TRunoff%wt(iunit,nt)-TRunoff%dwt(iunit,nt)*theDeltaT+TRunoff%etin(iunit,nt)*theDeltaT)
      Tdom%domTout(iunit,ntdom) = min(-TRunoff%etout(iunit,nt) *0.3_r8,Tdom%domTout(iunit,ntdom))
      Tdom%domTout(iunit,ntdom) = max(0._r8,Tdom%domTout(iunit,ntdom))
      Tdom%domTout(iunit,ntdom) = min((Tdom%domT(iunit,ntdom)+(Tdom%domsub(iunit,ntdom)+Tdom%domHout(iunit,ntdom))* theDeltaT)/theDeltaT,Tdom%domTout(iunit,ntdom))

      
      Tdom%domT(iunit,ntdom) = max(0._r8,Tdom%domT(iunit,ntdom) + ( Tdom%domsub(iunit,ntdom) + Tdom%domHout(iunit,ntdom) - Tdom%domTout(iunit,ntdom) ) * theDeltaT)
    end subroutine subnetworkRoutingDOM

    subroutine mainchannelRoutingDOM(iunit,nt,ntdom,theDeltaT)
      ! solve the ODE with Euler algorithm for main-channel routing
      implicit none
      integer, intent(in) :: iunit, nt, ntdom
      real(r8), intent(in) :: theDeltaT
      real(r8) :: temp_gwl
      temp_gwl = TRunoff%qgwl(iunit,nt) * TUnit%area(iunit) * TUnit%frac(iunit)

      Tdom%domRout(iunit,ntdom) = -TRunoff%erout(iunit,nt) * (Tdom%domR(iunit,ntdom) + (Tdom%domRUp(iunit,ntdom) + Tdom%domToutLat(iunit,ntdom)) * theDeltaT)/(TRunoff%wr(iunit,nt)-TRunoff%dwr(iunit,nt)*theDeltaT+(TRunoff%erlateral(iunit,nt)+TRunoff%erin(iunit,nt)+temp_gwl)*theDeltaT)
      Tdom%domRout(iunit,ntdom) = min(-TRunoff%erout(iunit,nt) *0.3_r8,Tdom%domRout(iunit,ntdom))
      Tdom%domRout(iunit,ntdom) = max(0._r8,Tdom%domRout(iunit,ntdom))
      Tdom%domRout(iunit,ntdom) = min((Tdom%domR(iunit,ntdom)+(Tdom%domRUp(iunit,ntdom) + Tdom%domToutLat(iunit,ntdom))* theDeltaT)/theDeltaT,Tdom%domRout(iunit,ntdom))

      
      Tdom%domR(iunit,ntdom) = max(0._r8,Tdom%domR(iunit,ntdom) + (Tdom%domRUp(iunit,ntdom) + Tdom%domToutLat(iunit,ntdom) - Tdom%domRout(iunit,ntdom)) * theDeltaT)
    end subroutine mainchannelRoutingDOM
!-------------------------------------------------------------------------
end MODULE DommasbMod
