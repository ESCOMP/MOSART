module rof_comp_nuopc
  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for MOSART
  !----------------------------------------------------------------------------

  use shr_kind_mod          , only : R8=>SHR_KIND_R8
  use shr_kind_mod          , only : CL=>SHR_KIND_CL, CXX => shr_kind_CXX
  use shr_sys_mod           , only : shr_sys_abort
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use shr_string_mod        , only : shr_string_listGetNum
  use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use esmFlds               , only : fldListFr, fldListTo, comprof, compname
  use esmFlds               , only : flds_scalar_name, flds_scalar_num
  use esmFlds               , only : flds_scalar_index_nx, flds_scalar_index_ny
  use esmFlds               , only : flds_scalar_index_rofice_present
  use esmFlds               , only : flds_scalar_index_flood_present
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_chkerr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS        => SetServices,       &
    model_label_Advance     => label_Advance,     &
    model_label_SetRunClock => label_SetRunClock, &
    model_label_Finalize    => label_Finalize

  use mosart_import_export , only : mosart_import, mosart_export
  use mosart_cpl_indices   , only : mosart_cpl_indices_set, nt_rtm, rtm_tracers
  use mosart_cpl_indices   , only : index_x2r_Flrl_rofsur, index_x2r_Flrl_rofi
  use mosart_cpl_indices   , only : index_x2r_Flrl_rofgwl, index_x2r_Flrl_rofsub
  use mosart_cpl_indices   , only : index_x2r_Flrl_irrig
  use mosart_cpl_indices   , only : index_r2x_Forr_rofl, index_r2x_Forr_rofi, index_r2x_Flrr_flood
  use mosart_cpl_indices   , only : index_r2x_Flrr_volr, index_r2x_Flrr_volrmch
  use RunoffMod            , only : rtmCTL, TRunoff
  use RtmVar               , only : rtmlon, rtmlat, ice_runoff, iulog
  use RtmVar               , only : nsrStartup, nsrContinue, nsrBranch
  use RtmVar               , only : inst_index, inst_suffix, inst_name, RtmVarSet
  use RtmSpmd              , only : RtmSpmdInit, masterproc, mpicom_rof, ROFID, iam, npes
  use RtmMod               , only : Rtmini, Rtmrun
  use RtmTimeManager       , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep
  use perf_mod             , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  public :: SetServices

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(CXX)             :: flds_r2x = ''
  character(CXX)             :: flds_x2r = ''
  real(r8), allocatable      :: x2r(:,:)
  real(r8), allocatable      :: r2x(:,:)
  integer                    :: nflds_r2x
  integer                    :: nflds_x2r
  character(len=*),parameter :: grid_option = "mesh" ! grid_de, grid_arb, grid_reg, mesh
  integer, parameter         :: dbug = 10
  integer                    :: dbrc
  logical                    :: flood_present        ! flag
  logical                    :: rof_prognostic       ! flag

  !----- formats -----
  character(*),parameter :: modName =  "(mosart_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Time)         :: currTime              ! Current time
    type(ESMF_Time)         :: startTime             ! Start time
    type(ESMF_Time)         :: stopTime              ! Stop time
    type(ESMF_Time)         :: refTime               ! Ref time
    type(ESMF_TimeInterval) :: timeStep              ! Model timestep
    type(ESMF_Calendar)     :: esmf_calendar         ! esmf calendar     
    type(ESMF_CalKind_Flag) :: esmf_caltype          ! esmf calendar type
    integer                 :: ref_ymd               ! reference date (YYYYMMDD)
    integer                 :: ref_tod               ! reference time of day (sec)
    integer                 :: yy,mm,dd              ! Temporaries for time query 
    integer                 :: start_ymd             ! start date (YYYYMMDD)
    integer                 :: start_tod             ! start time of day (sec)
    integer                 :: stop_ymd              ! stop date (YYYYMMDD)
    integer                 :: stop_tod              ! stop time of day (sec)
    integer                 :: curr_ymd              ! Start date (YYYYMMDD)
    integer                 :: curr_tod              ! Start time of day (sec)
    type(ESMF_VM)           :: vm
    integer                 :: mpicom
    character(CL)           :: cvalue
    logical                 :: exists
    integer                 :: lsize                 ! local array size
    integer                 :: ierr                  ! error code
    integer                 :: shrlogunit            ! original log unit
    integer                 :: shrloglev             ! original log level
    integer                 :: nsrest                ! restart type
    character(CL)           :: calendar              ! calendar type name
    character(CL)           :: username              ! user name
    character(CL)           :: caseid                ! case identifier name
    character(CL)           :: ctitle                ! case description title
    character(CL)           :: hostname              ! hostname of machine running on
    character(CL)           :: model_version         ! Model version
    character(CL)           :: starttype             ! start-type (startup, continue, branch, hybrid)
    logical                 :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    logical                 :: isPresent
    integer                 :: n
    character(len=512)      :: diro
    character(len=512)      :: logfile
    logical                 :: activefld
    character(CL)           :: stdname, shortname
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! initialize MOSART MPI communicator (mpicom_rof in RtmSpmd)
    !----------------------------------------------------------------------------

    call RtmSpmdInit(mpicom)

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ROFID  ! convert from string to integer

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) inst_index 

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ''
    end if

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    if (masterproc) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (len_trim(logfile) > 0) then
          iulog = shr_file_getUnit()
          open(iulog,file=trim(diro)//"/"//trim(logfile))
       else
          iulog = shrlogunit
       endif
    endif

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,format) "MOSART river model initialization"
       write(iulog,*) ' mosart npes = ',npes
       write(iulog,*) ' mosart iam  = ',iam
       write(iulog,*) ' inst_name = ',trim(inst_name)
    endif

    !----------------------
    ! Obtain attribute values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle=trim(caseid)

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) model_version

    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname

    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) username

    !----------------------
    ! Get properties from clock
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    !----------------------
    ! Initialize time info
    !----------------------

    call timemgr_setup(&
         calendar_in=calendar, &
         start_ymd_in=start_ymd, &
         start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, &
         ref_tod_in=ref_tod, &
         stop_ymd_in=stop_ymd, &
         stop_tod_in=stop_tod)

    !----------------------
    ! Initialize RtmVar module variables
    !----------------------

    !TODO: the following strings must not be hard-wired - must have module variables
    ! like seq_infodata_start_type_type - maybe another entry in seq_flds_mod?
    if (     trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim('branch')) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    call RtmVarSet(&
         caseid_in=caseid, &
         ctitle_in=ctitle,   &
         brnch_retain_casename_in=brnch_retain_casename, &
         nsrest_in=nsrest, &
         version_in=model_version,     &
         hostname_in=hostname, &
         username_in=username)

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    ! Note - all of the above needs to be done here to obtain rof_prognostic and flood_present
    ! This is a limitatin of the current mosart initialization

    call Rtmini(rtm_active=rof_prognostic, flood_active=flood_present)

    !--------------------------------
    ! create import and export field list and determine field indices of import/export arrays
    !--------------------------------

    call shr_nuopc_fldList_Concat(fldListFr(comprof), fldListTo(comprof), flds_r2x, flds_x2r, flds_scalar_name)

    call mosart_cpl_indices_set(flds_x2r, flds_r2x)

    !--------------------------------
    ! Advertise import and export fields
    !--------------------------------

    do n = 1,shr_nuopc_fldList_Getnumflds(fldListFr(comprof))
       call shr_nuopc_fldList_Getfldinfo(fldListFr(comprof), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(exportState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':Fr_'//trim(compname(comprof))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    do n = 1,shr_nuopc_fldList_Getnumflds(fldListTo(comprof))
       call shr_nuopc_fldList_Getfldinfo(fldListTo(comprof), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(importState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite(subname//':To_'//trim(compname(comprof))//': '//trim(shortname), ESMF_LOGMSG_INFO)
    end do

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer , allocatable  :: gindex(:)
    real(r8), pointer      :: elemCoords(:,:)
    real(r8), pointer      :: elemCornerCoords(:,:,:)
    real(r8)               :: dx,dy
    character(CL)          :: cvalue
    character(ESMF_MAXSTR) :: convCIM, purpComp
    type(ESMF_Mesh)        :: Emesh
    integer                :: shrlogunit            ! original log unit
    integer                :: shrloglev             ! original log level
    type(ESMF_VM)          :: vm
    logical                :: connected             ! is field connected?
    integer                :: lsize                 ! local size ofarrays
    integer                :: g,i,j,n,m,ni          ! indices
    integer                :: lbnum                 ! input to memory diagnostic
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)


#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_init_mct:start::',lbnum)
    endif
#endif

    !----------------------
    ! Determine field arays
    !----------------------

    nflds_r2x = shr_string_listGetNum(flds_r2x)
    nflds_x2r = shr_string_listGetNum(flds_x2r)
    lsize = rtmCTL%lnumr
    allocate(r2x(nflds_r2x,lsize))
    allocate(x2r(nflds_x2r,lsize))
    r2x(:,:)  = 0._r8
    x2r(:,:)  = 0._r8

    !--------------------------------
    ! Determine bounds numbering consistency
    !--------------------------------

    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       if (ni > rtmCTL%lnumr) then
          write(iulog,*) subname, ' : ERROR runoff count',n,ni,rtmCTL%lnumr
          call shr_sys_abort( subname //' ERROR: runoff > expected' )
       end if
    end do
    if (ni /= rtmCTL%lnumr) then
       write(iulog,*) subname, ' : ERROR runoff total count',ni,rtmCTL%lnumr
       call shr_sys_abort( subname//' ERROR: runoff not equal to expected' )
    endif

    !--------------------------------
    ! generate the mesh
    !--------------------------------

    allocate(gindex(lsize))
    allocate(elemCoords(2,lsize))          ! (lon+lat) * n_gridcells
    allocate(elemCornerCoords(2,4,lsize))  ! (lon+lat) * n_corners * n_gridcells

    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       gindex(ni) = rtmCTL%gindex(n)
       elemCoords(1,ni) = rtmCTL%lonc(n)
       elemCoords(2,ni) = rtmCTL%latc(n)
       ! TBD
       ! tcraig, MOSART does not define corner values and there is no info about grid sizes (ie. dx, dy)
       ! anywhere so make something up for now.  this has to be fixed if weights are generated on the fly!
       ! someone from MOSART has to define the corner lon and lat for use here.
       ! corners are defined counterclockwise
       do m = 1,4
          if (m == 1 .or. m == 4) dx = -0.05
          if (m == 2 .or. m == 3) dx =  0.05
          if (m == 1 .or. m == 2) dy = -0.05
          if (m == 3 .or. m == 4) dy =  0.05
          elemCornerCoords(1,m,ni) = rtmCTL%lonc(n) + dx
          elemCornerCoords(2,m,ni) = rtmCTL%latc(n) + dy
       enddo
    end do

    Emesh = ESMF_MeshCreate(parametricDim=2, &
       coordSys=ESMF_COORDSYS_SPH_DEG, &
       elementIds=gindex, &
       elementType=ESMF_MESHELEMTYPE_QUAD, &
       elementCoords=elemCoords, &
       elementCornerCoords=elemCornerCoords, &
       rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    deallocate(gindex)
    deallocate(elemCoords)
    deallocate(elemCornerCoords)

    !--------------------------------
    ! realize the actively coupled fields
    !--------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(comprof), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':mosartImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(comprof), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':mosartExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Create MOSART export state
    !--------------------------------

    call mosart_export( r2x )

    !--------------------------------
    ! Pack export state
    ! Copy from r2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(r2x, flds_r2x, exportState, grid_option, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(rtmlon), flds_scalar_index_nx, exportState, mpicom_rof, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(rtmlat), flds_scalar_index_ny, exportState, mpicom_rof, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flood_present) then
       call shr_nuopc_methods_State_SetScalar(1.0_r8, flds_scalar_index_flood_present, exportState, mpicom_rof, &
            flds_scalar_name, flds_scalar_num, rc)
    else
       call shr_nuopc_methods_State_SetScalar(0.0_r8, flds_scalar_index_flood_present, exportState, mpicom_rof, &
            flds_scalar_name, flds_scalar_num, rc)
    end if
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    ! For now hard-wire rofice_present to be false
    call shr_nuopc_methods_State_SetScalar(0.0_r8, flds_scalar_index_rofice_present, exportState, mpicom_rof, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "MOSART", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "MOSART River Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", "MOSART River Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "River", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", TBD, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(subname) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_int_mct:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! Run MOSART

    ! arguments:
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables:
    type(ESMF_Clock)  :: clock
    type(ESMF_Time)   :: currTime      ! Current time
    type(ESMF_Alarm)  :: alarm
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(CL)     :: cvalue
    integer           :: shrlogunit    ! original log unit
    integer           :: shrloglev     ! original log level
    integer           :: dtime         ! time step size
    integer           :: ymd_sync, ymd ! current date (YYYYMMDD)
    integer           :: yr_sync, yr   ! current year
    integer           :: mon_sync, mon ! current month
    integer           :: day_sync, day ! current day
    integer           :: tod_sync, tod ! current time of day (sec)
    logical           :: rstwr         ! .true. ==> write restart file before returning
    logical           :: nlend         ! .true. ==> signaling last time-step
    integer           :: lbnum         ! input to memory diagnostic
    integer           :: g,i           ! indices
    character(len=32) :: rdate         ! date char string for restart file names
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mosart_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1) then
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    endif

    !--------------------------------
    ! Unpack import state from mediator
    !--------------------------------

    call t_startf ('lc_mosart_import')
    call shr_nuopc_grid_StateToArray(importState, x2r, flds_x2r, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call mosart_import( x2r )
    call t_stopf ('lc_mosart_import')

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_restart', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='seq_timemgr_alarm_stop', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       nlend = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       nlend = .false.
    endif

    !--------------------------------
    ! Run MOSART
    !--------------------------------

    ! First advance mosart time step
    call advance_timestep()

    call ESMF_ClockGet( clock, currTime=currTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( currTime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync

    ! Export data is in rtmCTL and Trunoff data types
    call Rtmrun(rstwr, nlend, rdate)

    !--------------------------------
    ! Pack export state to mediator
    !--------------------------------

    ! (input is rtmCTL%runoff, output is r2x)
    call t_startf ('lc_rof_export')
    call mosart_export( r2x )

    call shr_nuopc_grid_ArrayToState(r2x, flds_r2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_rof_export')

    !--------------------------------
    ! Check that internal clock is in sync with master clock
    !--------------------------------

    ! Check that internal clock is in sync with master clock
    ! Note that the driver clock has not been updated yet - so at this point
    ! MOSART is actually 1 coupling intervals ahead of the driver clock

    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       write(iulog,*)' mosart ymd=',ymd     ,'  mosart tod= ',tod
       write(iulog,*)'   sync ymd=',ymd_sync,'    sync tod= ',tod_sync
       rc = ESMF_FAILURE
       call ESMF_LogWrite(subname//" MOSART clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR, rc=dbrc)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ROF from: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    endif

    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mosart_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=128)       :: mtimestring, dtimestring
    type(ESMF_Alarm),pointer :: alarmList(:)
    type(ESMF_Alarm)         :: dalarm
    integer                  :: alarmcount, n
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return


    !--------------------------------
    ! check that the current time in the model and driver are the same
    !--------------------------------

    if (mcurrtime /= dcurrtime) then
      call ESMF_TimeGet(dcurrtime, timeString=dtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_TimeGet(mcurrtime, timeString=mtimestring, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      call ESMF_LogWrite(subname//" ERROR in time consistency; "//trim(dtimestring)//" ne "//trim(mtimestring),  &
           ESMF_LOGMSG_ERROR, rc=dbrc)
      rc=ESMF_Failure
      return
    endif

    !--------------------------------                                                                                 
    ! force model clock currtime and timestep to match driver and set stoptime                                        
    !--------------------------------                                                                                 

    mstoptime = mcurrtime + dtimestep

    call ESMF_ClockSet(mclock, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    !--------------------------------
    ! copy alarms from driver to model clock if model clock has no alarms (do this only once!)
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    if (alarmCount == 0) then
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      allocate(alarmList(alarmCount))
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmList=alarmList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      do n = 1, alarmCount
         ! call ESMF_AlarmPrint(alarmList(n), rc=rc)
         ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
         dalarm = ESMF_AlarmCreate(alarmList(n), rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
         call ESMF_AlarmSet(dalarm, clock=mclock, rc=rc)
         if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
      enddo

      deallocate(alarmList)
    endif

    !--------------------------------                                                                                 
    ! Advance model clock to trigger alarms then reset model clock back to currtime                                   
    !--------------------------------                                                                                 

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(mosart_comp_nuopc) ',8a)"
    character(*), parameter :: F91   = "('(mosart_comp_nuopc) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'MOSART: end of main integration loop'
       write(iulog,F91)
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

end module rof_comp_nuopc
