module rof_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for MOSART
  !----------------------------------------------------------------------------

  use ESMF               , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSetEntryPoint, ESMF_State, &
                                  ESMF_Clock, ESMF_ClockGet, ESMF_ClockSet, ESMF_ClockAdvance, &
                                  ESMF_ClockGetAlarm, ESMF_ClockGetNextTime, ESMF_ClockGet, ESMF_ClockGetAlarm, &
                                  ESMF_ClockGetAlarmList, ESMF_Alarm, ESMF_AlarmSet, ESMF_AlarmIsRinging, &
                                  ESMF_AlarmRingerOff, ESMF_Time, ESMF_TimeGet, ESMF_TimeInterval, &
                                  ESMF_CalKind_Flag, ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN, &
                                  ESMF_Mesh, ESMF_MeshCreate, ESMF_FILEFORMAT_ESMFMESH, &
                                  ESMF_DistGrid, ESMF_DistGridCreate, &
                                  ESMF_MethodRemove,ESMF_VM, ESMF_VMGet, ESMF_LogFoundError, &
                                  ESMF_SUCCESS, ESMF_LogWrite, ESMF_FAILURE, ESMF_LogFoundError, &
                                  ESMF_LOGERR_PASSTHRU, ESMF_LOGMSG_ERROR, ESMF_LOGMSG_INFO, &
                                  ESMF_METHOD_INITIALIZE, ESMF_ALARMLIST_ALL, &
                                  operator(==), operator(/=), operator(<), operator(<=), &
                                  operator(>), operator(>=), operator(-), operator(+)
  use NUOPC              , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC              , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model        , only : model_routine_SS           => SetServices, &
                                  model_label_Advance        => label_Advance, &
                                  model_label_DataInitialize => label_DataInitialize, &
                                  model_label_SetRunClock    => label_SetRunClock, &
                                  model_label_Finalize       => label_Finalize, &
                                  SetVM, NUOPC_ModelGet
  use shr_kind_mod       , only : R8=>SHR_KIND_R8, CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use shr_sys_mod        , only : shr_sys_abort
  use shr_log_mod        , only : shr_log_getlogunit, shr_log_setlogunit
  use shr_cal_mod        , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use mosart_vars        , only : nsrStartup, nsrContinue, nsrBranch, &
                                  inst_index, inst_suffix, inst_name, &
                                  mainproc, mpicom_rof, iam, npes, iulog, &
                                  nsrest, caseid, ctitle, version, hostname, username
  use mosart_data        , only : ctl
  use mosart_driver      , only : mosart_read_namelist, mosart_init1, mosart_init2, mosart_run
  use mosart_timemanager , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep
  use mosart_io          , only : ncd_pio_init
  use mosart_restfile    , only : brnch_retain_casename
  use rof_import_export  , only : import_fields, export_fields, advertise_fields, realize_fields
  use nuopc_shr_methods  , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit, &
                                  set_component_logging, get_component_instance, log_clock_advance
  use perf_mod           , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  ! Module routines
  public  :: SetServices
  public  :: SetVM
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelSetRunClock
  private :: ModelAdvance
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=CL)       :: flds_scalar_name = ''
  integer                 :: flds_scalar_num = 0
  integer                 :: flds_scalar_index_nx = 0
  integer                 :: flds_scalar_index_ny = 0
  integer                 :: flds_scalar_index_nextsw_cday = 0._r8
  integer                 :: nthrds

  integer     , parameter :: debug = 1
  character(*), parameter :: modName =  "(rof_comp_nuopc)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output arguments
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
    type(ESMF_CalKind_Flag) :: esmf_caltype          ! esmf calendar type
    type(ESMF_VM)           :: vm                    ! esmf virtual machine
    integer                 :: ref_ymd               ! reference date (YYYYMMDD)
    integer                 :: ref_tod               ! reference time of day (sec)
    integer                 :: yy,mm,dd              ! Temporaries for time query
    integer                 :: start_ymd             ! start date (YYYYMMDD)
    integer                 :: start_tod             ! start time of day (sec)
    integer                 :: stop_ymd              ! stop date (YYYYMMDD)
    integer                 :: stop_tod              ! stop time of day (sec)
    integer                 :: curr_ymd              ! Start date (YYYYMMDD)
    integer                 :: curr_tod              ! Start time of day (sec)
    integer                 :: shrlogunit            ! original log unit
    integer                 :: n                     ! indices
    character(CL)           :: calendar              ! calendar type name
    character(CL)           :: starttype             ! start-type (startup, continue, branch, hybrid)
    logical                 :: isPresent, isSet
    character(CL)           :: cvalue
    character(len=CL)       :: logmsg
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_rof, peCount=npes, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    mainproc = (iam == 0)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    inst_name = "ROF"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! reset shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, mainproc, iulog, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_log_setLogUnit (iulog)

    !----------------------------------------------------------------------------
    ! advertise fields
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNX')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxGridNY')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nextsw_cday
       write(logmsg,*) flds_scalar_index_nextsw_cday
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nextsw_cday = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    ! Need to run the initial phase of mosart here in order to
    ! get the advertise phase correct

    !----------------------
    ! Obtain attribute values
    !----------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    read(cvalue,*) caseid
    ctitle=trim(caseid)

    call NUOPC_CompAttributeGet(gcomp, name='brnch_retain_casename', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) brnch_retain_casename

    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) starttype

    if (     trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim('branch')) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    call NUOPC_CompAttributeGet(gcomp, name='model_version', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) version

    call NUOPC_CompAttributeGet(gcomp, name='hostname', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) hostname

    call NUOPC_CompAttributeGet(gcomp, name='username', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) username

    !----------------------
    ! Get properties from clock
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, s=curr_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,curr_ymd)

    call ESMF_TimeGet( startTime, yy=yy, mm=mm, dd=dd, s=start_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,start_ymd)

    call ESMF_TimeGet( stopTime, yy=yy, mm=mm, dd=dd, s=stop_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,stop_ymd)

    call ESMF_TimeGet( refTime, yy=yy, mm=mm, dd=dd, s=ref_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yy,mm,dd,ref_ymd)

    call ESMF_TimeGet( currTime, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort( subname//'ERROR:: bad calendar for ESMF' )
    end if

    !----------------------
    ! Set time manager module variables
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
    ! Read namelist, grid and surface data
    !----------------------

    if (mainproc) then
       write(iulog,*) "MOSART river model initialization"
       write(iulog,*) ' mosart npes = ',npes
       write(iulog,*) ' mosart iam  = ',iam
       write(iulog,*) ' inst_name = ',trim(inst_name)
    endif

    !----------------------
    ! Read in mosart namelist
    !----------------------

    call mosart_read_namelist()

    !----------------------------------------------------------------------------
    ! Now advertise fields
    !----------------------------------------------------------------------------

    call advertise_fields(gcomp, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)       :: Emesh
    type(ESMF_VM)         :: vm
    integer , allocatable :: gindex(:)             ! global index space on my processor
    integer               :: lbnum                 ! input to memory diagnostic
    character(CL)         :: cvalue                ! temporary
    integer               :: shrlogunit            ! original log unit
    integer               :: lsize
    integer               :: n,ni
    integer               :: localPet
    integer               :: localPeCount
    integer               :: rofid                 ! component id for pio
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call shr_log_setLogUnit (iulog)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, localPet=localPet, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, pet=localPet, peCount=localPeCount, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize threading.  If ESMF_AWARE_THREADING is used localPeCount will be
    ! the thread count, otherwise the nthreads attribute is used.
    !----------------------------------------------------------------------------


    if(localPeCount == 1) then
       call NUOPC_CompAttributeGet(gcomp, "nthreads", value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) nthrds
    else
       nthrds = localPeCount
    endif

#if (defined _MEMTRACE)
    if (mainproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_nuopc_InitializeRealize:start::',lbnum)
    endif
#endif

    !-------------------------------------------------------
    ! Initialize mosart pio
    !-------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) rofid  ! convert from string to integer

    call ncd_pio_init(rofid)

    !-------------------------------------------------------
    ! Call first phase of mosart initialization (set decomp, grid)
    !-------------------------------------------------------

    ! - Initialize mosart time manager
    ! - Initialize number of mosart tracers
    ! - Read input data (river direction file) (global)
    ! - Deriver gridbox edges (global)
    ! - Determine mosart ocn/land mask (global)
    ! - Compute total number of basins and runoff ponts
    ! - Compute river basins, actually compute ocean outlet gridcell
    ! - Allocate basins to pes
    ! - Count and distribute cells to rglo2gdc (determine ctl%begr, ctl%endr)
    ! - Adjust area estimation from DRT algorithm for those outlet grids
    !     - useful for grid-based representation only
    !     - need to compute areas where they are not defined in input file

    call mosart_init1(rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! generate the mesh and realize fields
    !--------------------------------

    ! determine global index array
    lsize = ctl%endr - ctl%begr + 1
    allocate(gindex(lsize))
    ni = 0
    do n = ctl%begr,ctl%endr
       ni = ni + 1
       gindex(ni) = ctl%gindex(n)
    end do

    ! create distGrid from global index array
    ctl%DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gindex)

    ! read in the mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_rof', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (mainproc) then
       write(iulog,*)'mesh file for domain is ',trim(cvalue)
    end if

    EMesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, elementDistgrid=ctl%Distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize actively coupled fields
    !--------------------------------

    call realize_fields(gcomp,  Emesh, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------------------------------
    ! Initialize mosart maps and restart
    ! This must be called after the ESMF mesh is read in
    !-------------------------------------------------------

    call t_startf('mosarti_mosart_init')
    call mosart_init2(Emesh, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('mosarti_mosart_init')

    !--------------------------------
    ! Create MOSART export state
    !--------------------------------

    call export_fields(gcomp, ctl%begr, ctl%endr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set global grid size scalars in export state
    call State_SetScalar(dble(ctl%nlon), flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(ctl%nlat), flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Reset shr logging
    !----------------------------------------------------------------------------

    call shr_log_setLogUnit (shrlogunit)
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#if (defined _MEMTRACE)
    if(mainproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_nuopc_InitializeRealize:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    !------------------------
    ! Run MOSART
    !------------------------

    ! arguments:
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables:
    type(ESMF_Clock)  :: clock
    type(ESMF_Alarm)  :: alarm
    type(ESMF_Time)   :: currTime
    type(ESMF_Time)   :: nextTime
    type(ESMF_State)  :: importState
    type(ESMF_State)  :: exportState
    character(CL)     :: cvalue
    integer           :: shrlogunit    ! original log unit
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
    character(len=CS) :: rdate         ! date char string for restart file names
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call shr_log_getLogUnit (shrlogunit)
    call shr_log_setLogUnit (iulog)

#if (defined _MEMTRACE)
    if(mainproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mosart_comp_nuopc_ModelAdvance:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       rstwr = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       rstwr = .false.
    endif

    !--------------------------------
    ! Determine if time to stop
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_stop', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       nlend = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       nlend = .false.
    endif

    !--------------------------------
    ! Unpack import state from mediator
    !--------------------------------

    call t_startf ('lc_mosart_import')

    call import_fields(gcomp, ctl%begr, ctl%endr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf ('lc_mosart_import')

    !--------------------------------
    ! Run MOSART
    !--------------------------------

    ! Restart File - use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for restart file names

    call ESMF_ClockGetNextTime(clock, nextTime=nextTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(nexttime, yy=yr_sync, mm=mon_sync, dd=day_sync, s=tod_sync, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_cal_ymd2date(yr_sync, mon_sync, day_sync, ymd_sync)
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync, mon_sync, day_sync, tod_sync

    ! Advance mosart time step then run MOSART (export data is in ctl and Trunoff data types)
    call advance_timestep()
    call mosart_run(ctl%begr, ctl%endr, ctl%ntracers, rstwr, nlend, rdate, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state to mediator
    !--------------------------------

    call t_startf ('lc_rof_export')
    call export_fields(gcomp, ctl%begr, ctl%endr, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_rof_export')

    !--------------------------------
    ! Check that internal clock is in sync with sync clock
    !--------------------------------

    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    tod = tod

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       write(iulog,*)' mosart ymd=',ymd     ,'  mosart tod= ',tod
       write(iulog,*)'   sync ymd=',ymd_sync,'    sync tod= ',tod_sync
       rc = ESMF_FAILURE
       call ESMF_LogWrite(subname//" MOSART clock not in sync with sync clock",ESMF_LOGMSG_ERROR)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    end if

    if (debug > 1) then
       call log_clock_advance(clock, 'MOSART', iulog, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif


    !--------------------------------
    ! Reset shr logging to my original values
    !--------------------------------

    call shr_log_setLogUnit (shrlogunit)

#if (defined _MEMTRACE)
    if(mainproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','mosart_comp_nuopc_ModelAdvance:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option ! Restart option units
    integer                  :: restart_n      ! Number until restart interval
    integer                  :: restart_ymd    ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    character(len=256)       :: stop_option    ! Stop option units
    integer                  :: stop_n         ! Number until stop interval
    integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
    type(ESMF_ALARM)         :: stop_alarm
    character(len=128)       :: name
    integer                  :: alarmcount
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! query the Component for its clocks
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! set restart and stop alarms
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then

       call ESMF_GridCompGet(gcomp, name=name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite(subname//'setting alarms for' // trim(name), ESMF_LOGMSG_INFO)

       !----------------
       ! Restart alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="restart_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !----------------
       ! Stop alarm
       !----------------
       call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_n

       call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) stop_ymd

       call alarmInit(mclock, stop_alarm, stop_option, &
            opt_n   = stop_n,           &
            opt_ymd = stop_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_stop', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (mainproc) then
       write(iulog,F91)
       write(iulog,F00) 'MOSART: end of main integration loop'
       write(iulog,F91)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

end module rof_comp_nuopc
