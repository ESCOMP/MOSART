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
  private :: state_getimport1d
  private :: state_getimport2d
  private :: state_setexport1d
  private :: state_setexport2d

  type fld_list_type
     character(len=128) :: stdname
     integer :: ungridded_lbound = 0
     integer :: ungridded_ubound = 0
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

  subroutine advertise_fields(gcomp, flds_scalar_name, ntracers, rc)

    ! input/output variables
    type(ESMF_GridComp)            :: gcomp
    character(len=*) , intent(in)  :: flds_scalar_name
    integer          , intent(in)  :: ntracers
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

    call fldlist_add(fldsFrRof_num, fldsFrRof, trim(flds_scalar_name))
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl', ungridded_lbound=1, ungridded_ubound=ntracers-1)
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_flood')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volr')
    call fldlist_add(fldsFrRof_num, fldsFrRof, 'Flrr_volrmch')

    call NUOPC_CompAttributeGet(gcomp, name="flds_r2l_stream_channel_depths", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) read(cvalue,*) flds_r2l_stream_channel_depths
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
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsur', ungridded_lbound=1, ungridded_ubound=ntracers-1)
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofgwl')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofsub')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_rofi')
    call fldlist_add(fldsToRof_num, fldsToRof, 'Flrl_irrig')

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
    !call ESMF_StateGet(exportState, itemName=trim(fldsFrRof(2)%stdname), field=lfield, rc=rc)
    call ESMF_StateGet(exportState, itemName=trim(fldsFrRof(3)%stdname), field=lfield, rc=rc)
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
    integer          :: nliq, nfrz, nno3, nh2o
    integer          :: begg, endg
    integer          :: ntracers_liq, ntracers_tot
    real(r8), pointer:: fld2d(:,:)
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call NUOPC_ModelGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    begg = ctl%begr
    endg = ctl%endr
    ntracers_tot = ctl%ntracers
    ntracers_liq = ctl%ntracers - 1

    ! Set tracers
    nliq = 1
    nfrz = ntracers_tot

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will convert from input kg/m2s to m3/s
    ! NOTE: all liquid tracers need to come BEFORE the ice tracers -

    allocate(fld2d(ntracers_liq, begg:endg))
    call state_getimport2d(importState, 'Flrl_rofsur', begr, endr, ntracers_liq, ctl%area, output=fld2d, &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,ntracers_liq
       ctl%qsur(begg:endg,n) = fld2d(n,begg:endg)
       if (n > 1) then
          ctl%qsub(begg:endg,n) = 0._r8
          ctl%qgwl(begg:endg,n) = 0._r8
       end if
    end do
    deallocate(fld2d)

    call state_getimport1d(importState, 'Flrl_rofsub', begr, endr, ctl%area, output=ctl%qsub(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_rofgwl', begr, endr, ctl%area, output=ctl%qgwl(:,nliq), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_rofi', begr, endr, ctl%area, output=ctl%qsur(:,nfrz), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport1d(importState, 'Flrl_irrig', begr, endr, ctl%area, output=ctl%qirrig(:), &
         do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ctl%qsub(begr:endr, nfrz) = 0.0_r8
    ctl%qgwl(begr:endr, nfrz) = 0.0_r8

  end subroutine import_fields

  !====================================================================================
  subroutine export_fields (gcomp, begr, endr, ntracers_tot, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: begr, endr, ntracers_tot
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: exportState
    integer          :: n,nt
    integer          :: nliq, nfrz
    integer          :: ntracers_liq
    real(r8)         :: rofl(begr:endr,ntracers_tot-1)
    real(r8)         :: rofi(begr:endr)
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
    nliq = 1
    nfrz = ntracers_tot
    ntracers_liq = ntracers_tot - 1

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
       ! liquid and ice runoff are treated separately - this is what goes to the ocean
       do nt = 1,ntracers_liq
          do n = begr,endr
             rofl(n,nt) =  ctl%direct(n,nt) / (ctl%area(n)*0.001_r8)
             if (ctl%mask(n) >= 2) then
                rofl(n,nt) = rofl(n,nt) + ctl%runoff(n,nt) / (ctl%area(n)*0.001_r8)
             end if
          end do
       end do
       do n = begr,endr
          rofi(n) =  ctl%direct(n,nfrz) / (ctl%area(n)*0.001_r8)
          if (ctl%mask(n) >= 2) then
             rofi(n) = rofi(n) + ctl%runoff(n,nfrz) / (ctl%area(n)*0.001_r8)
          end if
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = begr,endr
          rofl(n,nliq) = (ctl%direct(n,nfrz) + ctl%direct(n,nliq)) / (ctl%area(n)*0.001_r8)
          if (ctl%mask(n) >= 2) then
             rofl(n,nliq) = rofl(n,nliq) + (ctl%runoff(n,nfrz) + ctl%runoff(n,nliq)) / (ctl%area(n)*0.001_r8)
          endif
          rofi(n) = 0._r8
       end do
    end if

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

    ! minus 1 below to only include liquid tracers
    call state_setexport2d(exportState, 'Forr_rofl', begr, endr, ntracers_liq, input=rofl, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Forr_rofi', begr, endr, input=rofi, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_flood', begr, endr, input=flood, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_volr', begr, endr, input=volr, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport1d(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, do_area_correction=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( flds_r2l_stream_channel_depths ) then
       call state_setexport1d(exportState, 'Sr_tdepth', begr, endr, input=tdepth, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call state_setexport1d(exportState, 'Sr_tdepth_max', begr, endr, input=tdepth_max, do_area_correction=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine export_fields

  !===============================================================================
  subroutine fldlist_add(num, fldlist, stdname, ungridded_lbound, ungridded_ubound)
    integer,                    intent(inout) :: num
    type(fld_list_type),        intent(inout) :: fldlist(:)
    character(len=*),           intent(in)    :: stdname
    integer,          optional, intent(in)    :: ungridded_lbound
    integer,          optional, intent(in)    :: ungridded_ubound

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


    if (present(ungridded_lbound) .and. present(ungridded_ubound)) then
       fldlist(num)%ungridded_lbound = ungridded_lbound
       fldlist(num)%ungridded_ubound = ungridded_ubound
    end if

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
             ! Create the field
             if (fldlist(n)%ungridded_lbound > 0 .and. fldlist(n)%ungridded_ubound > 0) then
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, &
                     ungriddedLbound=(/fldlist(n)%ungridded_lbound/), &
                     ungriddedUbound=(/fldlist(n)%ungridded_ubound/), &
                     gridToFieldMap=(/2/), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=stdname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(stdname)//" is connected using mesh", &
                  ESMF_LOGMSG_INFO)
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
  subroutine state_getimport1d(state, fldname, begr, endr, area, output, do_area_correction, rc)

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
    type(ESMF_Field)  :: lfield
    integer           :: g, i
    real(R8), pointer :: fldptr1d(:)
    character(len=*), parameter :: subname='(rof_import_export:state_getimport1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! check for nans
    if (any(isnan(fldptr1d))) then
       write(iulog,*) '# of NaNs = ', count(isnan(fldptr1d))
       write(iulog,*) 'Which are NaNs = ', isnan(fldptr1d)
       call shr_sys_abort(trim(subname //' One or more of the output from MOSART to the coupler are NaN'))
    end if

    if (do_area_correction) then
       fldptr1d(:) = fldptr1d(:) * med2mod_areacor(:)
    end if
    do g = begr,endr
       output(g) = fldptr1d(g-begr+1) * area(g)*0.001_r8
    end do

 end subroutine state_getimport1d

  !===============================================================================
  subroutine state_getimport2d(state, fldname, begr, endr, ntracers, area, output, do_area_correction, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_Field

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    integer             , intent(in)    :: ntracers
    real(r8)            , intent(in)    :: area(begr:endr)
    logical             , intent(in)    :: do_area_correction
    real(r8)            , intent(out)   :: output(begr:endr,ntracers)
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)  :: lfield
    integer           :: g, i, nt
    real(R8), pointer :: fldptr2d(:,:)
    character(len=*), parameter :: subname='(rof_import_export:state_getimport2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (any(isnan(fldptr2d))) then
       write(iulog,*) '# of NaNs = ', count(isnan(fldptr2d))
       write(iulog,*) 'Which are NaNs = ', isnan(fldptr2d)
       call shr_sys_abort(trim(subname) //' One or more of the output from MOSART to the coupler are NaN ' )
    end if

    if (do_area_correction) then
        fldptr2d(ntracers,:) = fldptr2d(ntracers,:) * med2mod_areacor(:)
    end if

    do nt = 1,ntracers
       do g = begr,endr
          output(g,nt) = fldptr2d(nt,g-begr+1) * area(g)*0.001_r8
       end do
    end do

 end subroutine state_getimport2d

  !===============================================================================

  subroutine state_setexport1d(state, fldname, begr, endr, input, do_area_correction, rc)

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
    real(R8), pointer           :: fldptr1d(:)
    character(len=*), parameter :: subname='(rof_import_export:state_setexport1d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Check for nans in input
    if (any(isnan(input))) then
       write(iulog,*) '# of NaNs = ', count(isnan(input))
       write(iulog,*) 'Which are NaNs = ', isnan(input)
       call shr_sys_abort(trim(subname) //' One or more of the output from MOSART to the coupler are NaN ' )
    end if

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set fldptr1d values to input array
    fldptr1d(:) = 0._r8
    do g = begr,endr
       fldptr1d(g-begr+1) = input(g)
    end do
    if (do_area_correction) then
       fldptr1d(:) = fldptr1d(:) * mod2med_areacor(:)
    end if

 end subroutine state_setexport1d

  !===============================================================================

  subroutine state_setexport2d(state, fldname, begr, endr, ntracers, input, do_area_correction, rc)

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
    integer             , intent(in)    :: ntracers
    real(r8)            , intent(in)    :: input(begr:endr,ntracers)
    logical             , intent(in)    :: do_area_correction
    integer             , intent(out)   :: rc

    ! local variables
    type(ESMF_Field)            :: lfield
    integer                     :: g, i, nt
    real(R8), pointer           :: fldptr2d(:,:)
    character(len=*), parameter :: subname='(rof_import_export:state_setexport2d)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Check for nans in input
    if (any(isnan(input))) then
       write(iulog,*) '# of NaNs = ', count(isnan(input))
       write(iulog,*) 'Which are NaNs = ', isnan(input)
       call shr_sys_abort(trim(subname) //' One or more of the output from MOSART to the coupler are NaN ' )
    end if

    ! get field pointer
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set fldptr2d values to input array
    fldptr2d(:,:) = 0._r8
    do nt = 1,ntracers
       do g = begr,endr
          fldptr2d(nt,g-begr+1) = input(g,nt)
       end do
       if (do_area_correction) then
          fldptr2d(nt,:) = fldptr2d(nt,:) * mod2med_areacor(:)
       end if
    end do

 end subroutine state_setexport2d

end module rof_import_export
