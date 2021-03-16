module rof_import_export

  use ESMF            , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF            , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF            , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF            , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF            , only : operator(/=), operator(==)
  use NUOPC           , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model     , only : NUOPC_ModelGet
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_sys_mod     , only : shr_sys_abort
  use RunoffMod       , only : rtmCTL, TRunoff
  use RtmVar          , only : iulog, nt_rtm, rtm_tracers
  use RtmSpmd         , only : masterproc, mpicom_rof
  use RtmTimeManager  , only : get_nstep
  use nuopc_shr_methods , only : chkerr

  implicit none
  private ! except

  public  :: advertise_fields
  public  :: realize_fields
  public  :: import_fields
  public  :: export_fields

  private :: fldlist_add
  private :: fldlist_realize
  private :: state_getimport
  private :: state_setexport
  private :: check_for_nans

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToRof_num = 0
  integer                :: fldsFrRof_num = 0
  type (fld_list_type)   :: fldsToRof(fldsMax)
  type (fld_list_type)   :: fldsFrRof(fldsMax)

  ! area correction factors for fluxes send and received from mediator
  real(r8), allocatable :: mod2med_areacor(:)
  real(r8), allocatable :: med2mod_areacor(:)

  character(*),parameter :: F01 = "('(mosart_import_export) ',a,i5,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine advertise_fields(gcomp, flds_scalar_name, do_rtm, do_rtmflood, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    logical          , intent(in)  :: do_rtm
    logical          , intent(in)  :: do_rtmflood
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: stdname
    character(ESMF_MAXSTR) :: cvalue
    integer                :: n, num
    character(len=128)     :: fldname
    character(len=*), parameter :: subname='(rof_import_export:advertise_fields)'
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Advertise export fields
    !--------------------------------

    if (do_rtm) then
       call fldlist_add(fldsFrRof_num, fldsFrRof, trim(flds_scalar_name))
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_flood')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volr')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volrmch')
    end if

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    if (do_rtm) then
       call fldlist_add(fldsToRof_num, fldsToRof, trim(flds_scalar_name))
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur')
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofgwl')
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsub')
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofdto')
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofi')
       call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_irrig')
    end if

    do n = 1,fldsToRof_num
       call NUOPC_Advertise(importState, standardName=fldsToRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

  end subroutine advertise_fields

  !===============================================================================
  subroutine realize_fields(gcomp, Emesh, flds_scalar_name, flds_scalar_num, rc)

    use ESMF          , only : ESMF_GridComp, ESMF_StateGet
    use ESMF          , only : ESMF_Mesh, ESMF_MeshGet
    use ESMF          , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldRegridGetArea
    use shr_const_mod , only : shr_const_rearth
    use shr_mpi_mod   , only : shr_mpi_min, shr_mpi_max

    ! input/output variables
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_Mesh)     , intent(in)    :: Emesh
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_State)      :: importState
    type(ESMF_State)      :: exportState
    type(ESMF_Field)      :: lfield
    integer               :: numOwnedElements
    integer               :: n,g
    real(r8), allocatable :: mesh_areas(:)
    real(r8), allocatable :: model_areas(:)
    real(r8), pointer     :: dataptr(:)
    real(r8)              :: re = shr_const_rearth*0.001_r8 ! radius of earth (km)
    real(r8)              :: max_mod2med_areacor
    real(r8)              :: max_med2mod_areacor
    real(r8)              :: min_mod2med_areacor
    real(r8)              :: min_med2mod_areacor
    real(r8)              :: max_mod2med_areacor_glob
    real(r8)              :: max_med2mod_areacor_glob
    real(r8)              :: min_mod2med_areacor_glob
    real(r8)              :: min_med2mod_areacor_glob
    character(len=*), parameter :: subname='(rof_import_export:realize_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call NUOPC_ModelGet(gcomp, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=ExportState, &
         fldList=fldsFrRof, &
         numflds=fldsFrRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MosartExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldlist_realize( &
         state=importState, &
         fldList=fldsToRof, &
         numflds=fldsToRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':MosartImport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine areas for regridding
    call ESMF_MeshGet(Emesh, numOwnedElements=numOwnedElements, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_StateGet(exportState, itemName=trim(fldsFrRof(2)%stdname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldRegridGetArea(lfield, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=dataptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    allocate(mesh_areas(numOwnedElements))
    mesh_areas(:) = dataptr(:)

    ! Determine model areas
    allocate(model_areas(numOwnedElements))
    allocate(mod2med_areacor(numOwnedElements))
    allocate(med2mod_areacor(numOwnedElements))
    n = 0
    do g = rtmCTL%begr,rtmCTL%endr
       n = n + 1
       model_areas(n) = rtmCTL%area(g)*1.0e-6_r8/(re*re)
       mod2med_areacor(n) = model_areas(n) / mesh_areas(n)
       med2mod_areacor(n) = mesh_areas(n) / model_areas(n)
    end do
    deallocate(model_areas)
    deallocate(mesh_areas)

    min_mod2med_areacor = minval(mod2med_areacor)
    max_mod2med_areacor = maxval(mod2med_areacor)
    min_med2mod_areacor = minval(med2mod_areacor)
    max_med2mod_areacor = maxval(med2mod_areacor)
    call shr_mpi_max(max_mod2med_areacor, max_mod2med_areacor_glob, mpicom_rof)
    call shr_mpi_min(min_mod2med_areacor, min_mod2med_areacor_glob, mpicom_rof)
    call shr_mpi_max(max_med2mod_areacor, max_med2mod_areacor_glob, mpicom_rof)
    call shr_mpi_min(min_med2mod_areacor, min_med2mod_areacor_glob, mpicom_rof)

    if (masterproc) then
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'MOSART'
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'MOSART'
    end if

  end subroutine realize_fields

  !===============================================================================
  subroutine import_fields( gcomp, rc )

    !---------------------------------------------------------------------------
    ! Obtain the runoff input from the mediator and convert from kg/m2s to m3/s
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: importState
    integer          :: n,nt
    integer          :: begr, endr
    integer          :: nliq, nfrz
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will convert from input kg/m2s to m3/s

    call state_getimport(importState, 'Flrl_rofsur', begr, endr, rtmCTL%area, output=rtmCTL%qsur(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofsub', begr, endr, rtmCTL%area, output=rtmCTL%qsub(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofgwl', begr, endr, rtmCTL%area, output=rtmCTL%qgwl(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofi', begr, endr, rtmCTL%area, output=rtmCTL%qsur(:,nfrz), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_irrig', begr, endr, rtmCTL%area, output=rtmCTL%qirrig(:), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    rtmCTL%qsub(begr:endr, nfrz) = 0.0_r8
    rtmCTL%qgwl(begr:endr, nfrz) = 0.0_r8

  end subroutine import_fields

  !====================================================================================
  subroutine export_fields (gcomp, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! uses
    use RtmVar, only : ice_runoff

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)  :: exportState
    integer           :: n,nt
    integer           :: begr,endr
    integer           :: nliq, nfrz
    real(r8), pointer :: rofl(:)
    real(r8), pointer :: rofi(:)
    real(r8), pointer :: flood(:)
    real(r8), pointer :: volr(:)
    real(r8), pointer :: volrmch(:)
    logical, save     :: first_time = .true.
    character(len=*), parameter :: subname='(rof_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    if (first_time) then
       if (masterproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
       endif
       first_time = .false.
    end if

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    allocate(rofl(begr:endr))
    allocate(rofi(begr:endr))
    allocate(flood(begr:endr))
    allocate(volr(begr:endr))
    allocate(volrmch(begr:endr))

    if ( ice_runoff )then
       ! separate liquid and ice runoff
       do n = begr,endr
          rofl(n) =  rtmCTL%direct(n,nliq) / (rtmCTL%area(n)*0.001_r8)
          rofi(n) =  rtmCTL%direct(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             ! liquid and ice runoff are treated separately - this is what goes to the ocean
             rofl(n) = rofl(n) + rtmCTL%runoff(n,nliq) / (rtmCTL%area(n)*0.001_r8)
             rofi(n) = rofi(n) + rtmCTL%runoff(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          end if
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = begr,endr
          rofl(n) = (rtmCTL%direct(n,nfrz) + rtmCTL%direct(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             rofl(n) = rofl(n) + (rtmCTL%runoff(n,nfrz) + rtmCTL%runoff(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          endif
          rofi(n) = 0._r8
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ! scs: is there a reason for the wr+wt rather than volr (wr+wt+wh)?
    ! volr(n) = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / rtmCTL%area(n)

    do n = begr, endr
       flood(n)   = -rtmCTL%flood(n)    / (rtmCTL%area(n)*0.001_r8)
       volr(n)    =  rtmCTL%volr(n,nliq)/ rtmCTL%area(n)
       volrmch(n) =  Trunoff%wr(n,nliq) / rtmCTL%area(n)
    end do

    call state_setexport(exportState, 'Forr_rofl', begr, endr, input=rofl, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofi', begr, endr, input=rofi, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_flood', begr, endr, input=flood, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volr', begr, endr, input=volr, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(rofl, rofi, flood, volr, volrmch)

  end subroutine export_fields

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname

    ! local variables
    integer :: rc
    character(len=*), parameter :: subname='(rof_import_export:fldlist_add)'
    !-------------------------------------------------------------------------------

    ! Set up a list of field information

    num = num + 1
    if (num > fldsMax) then
       call ESMF_LogWrite(trim(subname)//": ERROR num > fldsMax "//trim(stdname), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__)
       return
    endif
    fldlist(num)%stdname = trim(stdname)

  end subroutine fldlist_add

  !===============================================================================
  subroutine fldlist_realize(state, fldList, numflds, flds_scalar_name, flds_scalar_num, mesh, tag, rc)

    use NUOPC , only : NUOPC_IsConnected, NUOPC_Realize
    use ESMF  , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF  , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Mesh, ESMF_StateRemove
    use ESMF  , only : ESMF_LogFoundError, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LOGERR_PASSTHRU

    type(ESMF_State)    , intent(inout) :: state
    type(fld_list_type) , intent(in)    :: fldList(:)
    integer             , intent(in)    :: numflds
    character(len=*)    , intent(in)    :: flds_scalar_name
    integer             , intent(in)    :: flds_scalar_num
    character(len=*)    , intent(in)    :: tag
    type(ESMF_Mesh)     , intent(in)    :: mesh
    integer             , intent(inout) :: rc

    ! local variables
    integer                :: dbrc
    integer                :: n
    type(ESMF_Field)       :: field
    character(len=80)      :: stdname
    character(len=*),parameter  :: subname='(rof_import_export:fldlist_realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    do n = 1, numflds
       stdname = fldList(n)%stdname
       if (NUOPC_IsConnected(state, fieldName=stdname)) then
          if (stdname == trim(flds_scalar_name)) then
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected on root pe", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the scalar field
             call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          else
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             ! Create the field
             field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          endif

          ! NOW call NUOPC_Realize
          call NUOPC_Realize(state, field=field, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       else
          if (stdname /= trim(flds_scalar_name)) then
             call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(stdname) // " is not connected.", &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             call ESMF_StateRemove(state, (/stdname/), rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          end if
       end if
    end do

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8

      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(rof_import_export:SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), grid=grid, typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), ungriddedUBound=(/flds_scalar_num/), gridToFieldMap=(/2/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine fldlist_realize

  !===============================================================================
  subroutine state_getimport(state, fldname, begr, endr, area, output, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_Field

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: area(begr:endr)
    logical             , intent(in)    :: do_area_correction
    real(r8)            , intent(out)   :: output(begr:endr)
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    character(len=*), parameter :: subname='(rof_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine output array and scale by unit convertsion
    if (do_area_correction) then
       fldptr(:) = fldptr(:) * med2mod_areacor(:)
    end if
    do g = begr,endr
       output(g) = fldptr(g-begr+1) * area(g)*0.001_r8
    end do

    ! check for nans
    call check_for_nans(fldptr, trim(fldname), begr)

  end subroutine state_getimport

  !===============================================================================
  subroutine state_setexport(state, fldname, begr, endr, input, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map input array to export state field
    ! ----------------------------------------------

    use ESMF         , only : ESMF_StateGet, ESMF_FieldGet, ESMF_Field
    use shr_const_mod, only : fillvalue=>shr_const_spval

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: input(begr:endr)
    logical             , intent(in)    :: do_area_correction
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    character(len=*), parameter :: subname='(rof_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set fldptr values to input array
    fldptr(:) = 0._r8
    do g = begr,endr
       fldptr(g-begr+1) = input(g)
    end do
    if (do_area_correction) then
       fldptr(:) = fldptr(:) * mod2med_areacor(:)
    end if

    ! check for nans
    call check_for_nans(fldptr, trim(fldname), begr)

  end subroutine state_setexport

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)

    ! uses
    use shr_infnan_mod, only : isnan => shr_infnan_isnan

    ! input/output variables
    real(r8), pointer             :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg

    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    ! Check if any input from mediator or output to mediator is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from MOSART to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module rof_import_export
