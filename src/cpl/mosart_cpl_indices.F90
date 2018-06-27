module mosart_cpl_indices

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  !    Module containing the indices for the fields passed between MOSART and
  !    the driver.
  !-----------------------------------------------------------------------

  ! !USES:
  use shr_sys_mod,    only : shr_sys_abort
  implicit none

  private                              ! By default make data private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: mosart_cpl_indices_set        ! Set the coupler indices

  ! !PUBLIC DATA MEMBERS:
  integer, public :: index_x2r_Flrl_rofsur = 0 ! lnd->rof liquid surface runoff forcing from land
  integer, public :: index_x2r_Flrl_rofgwl = 0 ! lnd->rof liquid gwl runoff from land
  integer, public :: index_x2r_Flrl_rofsub = 0 ! lnd->rof liquid subsurface runoff from land
  integer, public :: index_x2r_Flrl_rofdto = 0 ! lnd->rof liquid direct to ocean runoff
  integer, public :: index_x2r_Flrl_rofi   = 0 ! lnd->rof ice runoff forcing from land
  integer, public :: index_x2r_Flrl_irrig  = 0 ! lnd->rof fraction of volr to be removed for irrigation
  integer, public :: nflds_x2r             = 0

  !TODO - nt_rtm and rtm_tracers need to be removed and set by access to the index array
  integer, parameter, public :: nt_rtm = 2    ! number of tracers
  character(len=3), parameter, public :: rtm_tracers(nt_rtm) =  (/'LIQ','ICE'/)

  ! roff to driver (part of land for now) (optional if ROF is off)

  integer, public :: index_r2x_Forr_rofl    = 0 ! rof->ocn liquid runoff to ocean
  integer, public :: index_r2x_Forr_rofi    = 0 ! rof->ocn ice runoff to ocean
  integer, public :: index_r2x_Flrr_flood   = 0 ! rof->lnd flood runoff (>fthresh) back to land
  integer, public :: index_r2x_Flrr_volr    = 0 ! rof->lnd volr total volume back to land
  integer, public :: index_r2x_Flrr_volrmch = 0 ! rof->lnd volr main channel back to land
  integer, public :: nflds_r2x              = 0

!=======================================================================
contains
!=======================================================================

  subroutine mosart_cpl_indices_set(flds_x2r, flds_r2x )

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set the coupler indices needed by the mosart model coupler interface.
    ! runoff - (mosart -> ocn) and (mosart->lnd)
    !
    use mct_mod, only: mct_aVect, mct_aVect_init, mct_avect_indexra
    use mct_mod, only: mct_aVect_clean, mct_avect_nRattr
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: flds_x2r
    character(len=*), intent(in) :: flds_r2x
    !
    ! !LOCAL VARIABLES:
    type(mct_aVect)   :: avtmp      ! temporary av
    character(len=32) :: subname = 'mosart_cpl_indices_set'  ! subroutine name
    !-----------------------------------------------------------------------

    !-------------------------------------------------------------
    ! driver -> mosart
    !-------------------------------------------------------------

    call mct_aVect_init(avtmp, rList=flds_x2r, lsize=1)

    index_x2r_Flrl_rofsur = mct_avect_indexra(avtmp,'Flrl_rofsur',perrwith='quiet')
    if (index_x2r_Flrl_rofsur == 0) then
       index_x2r_Flrl_rofsur = mct_avect_indexra(avtmp,'Frxx_rofsur',perrwith='quiet')
    end if
    if (index_x2r_Flrl_rofsur == 0) then
       call shr_sys_abort('ERROR: subname: index_x2r_Flrl_rofsur is 0')
    end if

    index_x2r_Flrl_rofgwl = mct_avect_indexra(avtmp,'Flrl_rofgwl', perrwith='quiet')
    if (index_x2r_Flrl_rofgwl == 0) then
       index_x2r_Flrl_rofgwl = mct_avect_indexra(avtmp,'Frxx_rofgwl',perrwith='quiet')
    end if
    if (index_x2r_Flrl_rofgwl == 0) then
       call shr_sys_abort('ERROR: subname: index_x2r_Flrl_rofgwl is 0')
    end if

    index_x2r_Flrl_rofsub = mct_avect_indexra(avtmp,'Flrl_rofsub', perrwith='quiet')
    if (index_x2r_Flrl_rofsub == 0) then
       index_x2r_Flrl_rofsub = mct_avect_indexra(avtmp,'Frxx_rofsub',perrwith='quiet')
    end if
    if (index_x2r_Flrl_rofsub == 0) then
       call shr_sys_abort('ERROR: subname: index_x2r_Flrl_rofsub is 0')
    end if

    index_x2r_Flrl_rofdto = mct_avect_indexra(avtmp,'Flrl_rofdto',perrwith='quiet')
    if (index_x2r_Flrl_rofdto == 0) then
       index_x2r_Flrl_rofdto = mct_avect_indexra(avtmp,'Frxx_rofdto',perrwith='quiet')
    end if

    index_x2r_Flrl_rofi = mct_avect_indexra(avtmp,'Flrl_rofi', perrwith='quiet')
    if (index_x2r_Flrl_rofi == 0) then
       index_x2r_Flrl_rofi = mct_avect_indexra(avtmp,'Frxx_rofi',perrwith='quiet')
    end if
    if (index_x2r_Flrl_rofi == 0) then
       call shr_sys_abort('ERROR: subname: index_x2r_Flrl_rofi is 0')
    end if

    index_x2r_Flrl_irrig = mct_avect_indexra(avtmp,'Flrl_irrig', perrwith='quiet')
    if (index_x2r_Flrl_irrig == 0) then
       index_x2r_Flrl_irrig = mct_avect_indexra(avtmp,'Frxx_irrig',perrwith='quiet')
    end if
    if (index_x2r_Flrl_irrig == 0) then
       call shr_sys_abort('ERROR: subname: index_x2r_Flrl_irrig is 0')
    end if

    nflds_x2r = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

    !-------------------------------------------------------------
    ! mosart -> driver
    !-------------------------------------------------------------

    call mct_aVect_init(avtmp, rList=flds_r2x, lsize=1)

    index_r2x_Forr_rofl    = mct_avect_indexra(avtmp,'Forr_rofl')
    index_r2x_Forr_rofi    = mct_avect_indexra(avtmp,'Forr_rofi')
    index_r2x_Flrr_flood   = mct_avect_indexra(avtmp,'Flrr_flood')
    index_r2x_Flrr_volr    = mct_avect_indexra(avtmp,'Flrr_volr')
    index_r2x_Flrr_volrmch = mct_avect_indexra(avtmp,'Flrr_volrmch')

    nflds_r2x = mct_avect_nRattr(avtmp)

    call mct_aVect_clean(avtmp)

  end subroutine mosart_cpl_indices_set

end module mosart_cpl_indices
