module rof_import_export

  use ESMF                , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF                , only : ESMF_LogWrite, ESMF_LOGMSG_ERROR, ESMF_LogFoundError
  use ESMF                , only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF                , only : operator(/=), operator(==)
  use NUOPC               , only : NUOPC_CompAttributeGet, NUOPC_Advertise, NUOPC_IsConnected
  use NUOPC_Model         , only : NUOPC_ModelGet
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_abort
  use mosart_vars         , only : iulog, mainproc, mpicom_rof, ice_runoff
  use mosart_data         , only : ctl, TRunoff, TUnit
  use mosart_timemanager  , only : get_nstep
  use nuopc_shr_methods   , only : chkerr

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
  private :: fldchk

  type fld_list_type
     character(len=128) :: stdname
  end type fld_list_type

  integer, parameter     :: fldsMax = 100
  integer                :: fldsToRof_num = 0
  integer                :: fldsFrRof_num = 0
  logical                :: flds_r2l_stream_channel_depths = .false.   ! If should pass the channel depth fields needed for the hillslope model
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

  subroutine advertise_fields(gcomp, flds_scalar_name, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(out) :: rc

    ! local variables
    type(ESMF_State)       :: importState
    type(ESMF_State)       :: exportState
    character(ESMF_MAXSTR) :: cvalue          ! Character string read from driver attribute
    logical                :: isPresent       ! Atribute is present
    logical                :: isSet           ! Atribute is set
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

    call NUOPC_CompAttributeGet(gcomp, name="flds_r2l_stream_channel_depths", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) read(cvalue,*) flds_r2l_stream_channel_depths

    call fldlist_add(fldsFrRof_num, fldsFrRof, trim(flds_scalar_name))
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl_glc')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi_glc')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_flood')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volr')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volrmch')
    if ( flds_r2l_stream_channel_depths )then
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Sr_tdepth')
       call fldlist_add(fldsFrRof_num, fldsFrRof, 'Sr_tdepth_max')
    end if

    do n = 1,fldsFrRof_num
       call NUOPC_Advertise(exportState, standardName=fldsFrRof(n)%stdname, &
            TransferOfferGeomObject='will provide', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    enddo

    !--------------------------------
    ! Advertise import fields
    !--------------------------------

    call fldlist_add(fldsToRof_num, fldsToRof, trim(flds_scalar_name))
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofgwl')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsub')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofi')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_irrig')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Fgrg_rofl') ! liq runoff from glc
    call fldlist_add(fldsToRof_num, fldsToRof, 'Fgrg_rofi') ! ice runoff from glc

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
    do g = ctl%begr,ctl%endr
       n = n + 1
       model_areas(n) = ctl%area(g)*1.0e-6_r8/(re*re)
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

    if (mainproc) then
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_mod2med_areacor, max_mod2med_areacor ',&
            min_mod2med_areacor_glob, max_mod2med_areacor_glob, 'MOSART'
       write(iulog,'(2A,2g23.15,A )') trim(subname),' :  min_med2mod_areacor, max_med2mod_areacor ',&
            min_med2mod_areacor_glob, max_med2mod_areacor_glob, 'MOSART'
    end if

    if (fldchk(importState, 'Fgrg_rofl') .and. fldchk(importState, 'Fgrg_rofl')) then
      ctl%rof_from_glc = .true.
    else
      ctl%rof_from_glc = .false.
    end if
    if (mainproc) then
      write(iulog,'(A,l1)') trim(subname) //' rof from glc is ',ctl%rof_from_glc
    end if

  end subroutine realize_fields

  !===============================================================================
  subroutine import_fields( gcomp, begr, endr, rc )

    !---------------------------------------------------------------------------
    ! Obtain the runoff input from the mediator and convert from kg/m2s to m3/s
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: begr, endr
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: importState
    integer          :: n,nt
    integer          :: nliq, nice
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    nliq = ctl%nt_liq
    nice = ctl%nt_ice

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will convert from input kg/m2s to m3/s

    call state_getimport(importState, 'Flrl_rofsur', begr, endr, ctl%area, output=ctl%qsur(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofsub', begr, endr, ctl%area, output=ctl%qsub(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofgwl', begr, endr, ctl%area, output=ctl%qgwl(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofi', begr, endr, ctl%area, output=ctl%qsur(:,nice), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_irrig', begr, endr, ctl%area, output=ctl%qirrig(:), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ctl%qsub(begr:endr, nice) = 0.0_r8
    ctl%qgwl(begr:endr, nice) = 0.0_r8

    if (ctl%rof_from_glc) then
      call state_getimport(importState, 'Fgrg_rofl', begr, endr, ctl%area, output=ctl%qglc_liq(:), &
           do_area_correction=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call state_getimport(importState, 'Fgrg_rofi', begr, endr, ctl%area, output=ctl%qglc_ice(:), &
           do_area_correction=.true., rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
      ctl%qglc_liq(:) = 0._r8
      ctl%qglc_ice(:) = 0._r8
    end if

  end subroutine import_fields

  !====================================================================================
  subroutine export_fields (gcomp, begr, endr, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: begr, endr
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: exportState
    integer          :: n,nt
    integer          :: nliq, nice
    real(r8)         :: rofl(begr:endr)
    real(r8)         :: rofi(begr:endr)
    real(r8)         :: rofl_glc(begr:endr)
    real(r8)         :: rofi_glc(begr:endr)
    real(r8)         :: flood(begr:endr)
    real(r8)         :: volr(begr:endr)
    real(r8)         :: volrmch(begr:endr)
    real(r8)         :: tdepth(begr:endr)
    real(r8)         :: tdepth_max(begr:endr)
    logical, save    :: first_time = .true.
    character(len=*), parameter :: subname='(rof_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call NUOPC_ModelGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = ctl%nt_liq
    nice = ctl%nt_ice

    if (first_time) then
       if (mainproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
       endif
       first_time = .false.
    end if

    if ( ice_runoff )then
       ! separate liquid and ice runoff
       do n = begr,endr
          rofl(n) =  ctl%direct(n,nliq) / (ctl%area(n)*0.001_r8)
          rofi(n) =  ctl%direct(n,nice) / (ctl%area(n)*0.001_r8)
          if (ctl%mask(n) >= 2) then
             ! liquid and ice runoff are treated separately - this is what goes to the ocean
             rofl(n) = rofl(n) + ctl%runoff(n,nliq) / (ctl%area(n)*0.001_r8)
             rofi(n) = rofi(n) + ctl%runoff(n,nice) / (ctl%area(n)*0.001_r8)
          end if
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = begr,endr
          rofl(n) = (ctl%direct(n,nice) + ctl%direct(n,nliq)) / (ctl%area(n)*0.001_r8)
          if (ctl%mask(n) >= 2) then
             rofl(n) = rofl(n) + (ctl%runoff(n,nice) + ctl%runoff(n,nliq)) / (ctl%area(n)*0.001_r8)
          endif
          rofi(n) = 0._r8
       end do
    end if

    do n = begr,endr
      rofl_glc(n) = ctl%direct_glc(n,nliq) / (ctl%area(n)*0.001_r8)
      rofi_glc(n) = ctl%direct_glc(n,nice) / (ctl%area(n)*0.001_r8)
    end do

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ! scs: is there a reason for the wr+wt rather than volr (wr+wt+wh)?
    ! volr(n) = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / ctl%area(n)

    do n = begr, endr
       flood(n)   = -ctl%flood(n)    / (ctl%area(n)*0.001_r8)
       volr(n)    =  ctl%volr(n,nliq)/ ctl%area(n)
       volrmch(n) =  Trunoff%wr(n,nliq) / ctl%area(n)
       if ( flds_r2l_stream_channel_depths )then
          tdepth(n)  = Trunoff%yt(n,nliq)
          ! assume height to width ratio is the same for tributaries and main channel
          tdepth_max(n) = max(TUnit%twidth0(n),0._r8)*(TUnit%rdepth(n)/TUnit%rwidth(n))
        end if
    end do

    call state_setexport(exportState, 'Forr_rofl', begr, endr, input=rofl, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofi', begr, endr, input=rofi, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofl_glc', begr, endr, input=rofl_glc, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofi_glc', begr, endr, input=rofi_glc, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_flood', begr, endr, input=flood, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volr', begr, endr, input=volr, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( flds_r2l_stream_channel_depths ) then
       call state_setexport(exportState, 'Sr_tdepth', begr, endr, input=tdepth, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_setexport(exportState, 'Sr_tdepth_max', begr, endr, input=tdepth_max, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

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
    real(r8)         , pointer    :: array(:)
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

  !===============================================================================
  logical function fldchk(state, fldname)
    ! ----------------------------------------------
    ! Determine if field with fldname is in the input state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------
    call ESMF_StateGet(state, trim(fldname), itemFlag)
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       fldchk = .true.
    else
       fldchk = .false.
    endif
  end function fldchk

end module rof_import_export
