module rof_comp_esmf

  !----------------------------------------------------------------------------
  ! This is the LILAC cap for MOSART
  !----------------------------------------------------------------------------

  ! external libraries
  use ESMF
  use mpi               , only : MPI_BCAST, MPI_CHARACTER
  use mct_mod           , only : mct_world_init
  use perf_mod          , only : t_startf, t_stopf, t_barrierf
  use lilac_utils       , only : lilac_field_bundle_to_land, lilac_field_bundle_fr_land

  ! cime share code
  use shr_kind_mod      , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_sys_mod       , only : shr_sys_abort
  use shr_file_mod      , only : shr_file_setLogUnit, shr_file_getLogUnit
  use shr_cal_mod       , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use shr_nl_mod        , only : shr_nl_find_group_name

  ! mosart code
  use RtmVar            , only : rtmlon, rtmlat, iulog
  use RtmVar            , only : nsrStartup, nsrContinue, nsrBranch
  use RtmVar            , only : inst_index, inst_suffix, inst_name, RtmVarSet
  use RtmSpmd           , only : RtmSpmdInit, masterproc, mpicom_rof, rofid, iam, npes
  use RunoffMod         , only : rtmCTL
  use RtmMod            , only : Rtmini, Rtmrun
  use RtmTimeManager    , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep
  use rof_import_export , only : import_fields, export_fields
  use rof_shr_methods   , only : chkerr, state_diagnose

  implicit none
  private ! except

  ! Module routines
  public :: rof_register ! register mosart initial, run, final methods
  public :: rof_init     ! mosart initialization
  public :: rof_run      ! mosart run phase
  public :: rof_final    ! mosart finalization/cleanup

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer     , parameter :: debug = 1
  character(*), parameter :: modName =  "(rof_comp_esmf)"
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine rof_register(comp, rc)

    ! Register the mosart initial, run, and final phase methods with ESMF.

    ! input/output argumenents
    type(ESMF_GridComp)  :: comp  ! MOSART grid component
    integer, intent(out) :: rc    ! return status
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogSet ( flush =.true.)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, rof_init, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, rof_run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, rof_final, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

    call ESMF_LogWrite("rof gridcompset entry points finished!", ESMF_LOGMSG_INFO)

  end subroutine rof_register

  !===============================================================================

  subroutine rof_init(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    integer                 :: mpicom_vm
    logical                 :: flood_present     ! flag
    logical                 :: rof_prognostic    ! flag
    integer                 :: shrlogunit        ! original log unit
    integer                 :: lsize             ! local size ofarrays
    integer                 :: n,ni              ! indices
    integer                 :: lbnum             ! input to memory diagnostic
    integer                 :: nsrest            ! restart type
    integer                 :: ierr              ! error code
    character(CL)           :: cvalue            ! temporary 
    character(CL)           :: caseid            ! case identifier name
    character(CL)           :: starttype         ! start-type (startup, continue, branch, hybrid)

    ! mesh generation
    type(ESMF_Mesh)         :: rof_mesh
    character(ESMF_MAXSTR)  :: rof_mesh_filename ! full filepath of river mesh file
    type(ESMF_DistGrid)     :: distgrid          ! esmf global index space descriptor
    integer                 :: fileunit          ! input fileunit for reading mesh
    integer , allocatable   :: gindex(:)         ! global index space on my processor

    ! generation of field bundles
    type(ESMF_FieldBundle)  :: c2r_fb            ! field bundle in import state from river
    type(ESMF_FieldBundle)  :: r2c_fb            ! field bundle in export state to river

    ! clock info
    type(ESMF_Time)         :: currTime          ! Current time
    type(ESMF_Time)         :: startTime         ! Start time
    type(ESMF_Time)         :: stopTime          ! Stop time
    type(ESMF_Time)         :: refTime           ! Ref time
    type(ESMF_TimeInterval) :: timeStep          ! Model timestep
    type(ESMF_Calendar)     :: esmf_calendar     ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype      ! esmf calendar type
    character(CL)           :: calendar          ! calendar type name
    integer                 :: dtime_lilac
    integer                 :: ref_ymd           ! reference date (YYYYMMDD)
    integer                 :: ref_tod           ! reference time of day (sec)
    integer                 :: yy,mm,dd          ! Temporaries for time query
    integer                 :: start_ymd         ! start date (YYYYMMDD)
    integer                 :: start_tod         ! start time of day (sec)
    integer                 :: stop_ymd          ! stop date (YYYYMMDD)
    integer                 :: stop_tod          ! stop time of day (sec)
    integer                 :: curr_ymd          ! Start date (YYYYMMDD)
    integer                 :: curr_tod          ! Start time of day (sec)
    
    ! input namelist read for mosart mesh and run info
    namelist /lilac_rof_input/ rof_mesh_filename
    namelist /lilac_run_input/ caseid, starttype

    character(len=*), parameter :: subname=trim(modName)//':(rof_init) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !------------------------------------------------------------------------
    ! Query VM for local PET and mpi communicator
    !------------------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    call ESMF_LogWrite(subname//"ESMF_VMGet", ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! initialize MOSART MPI communicator
    !----------------------------------------------------------------------------

    ! The following call initializees the module variable mpicom_rof in RtmSpmd
    call RtmSpmdInit(mpicom_vm)

    ! Set ROFID - needed for the mosart code that requires MCT
    rofid = 1

    !------------------------------------------------------------------------
    !--- Log File ---
    !------------------------------------------------------------------------

    ! TODO: by default iulog = 6 in mosart_varctl - this should be generalized so that we
    ! can control the output log file for mosart running with a lilac driver

    inst_name  = 'LND'; inst_index  = 1; inst_suffix = ""

    ! Initialize io log unit
    call shr_file_getLogUnit (shrlogunit)
    if (.not. masterproc) then
       iulog = shrlogunit  ! All shr code output will go to iulog for masterproc
    end if
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,*) "========================================="
       write(iulog,*) " starting (rof_comp_esmf): rof_comp_init "
       write(iulog,*) " MOSART river model initialization"
    end if

#if (defined _MEMTRACE)
    if (masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_esmf::',lbnum)
    endif
#endif

    !----------------------
    ! Set time manager module variables
    !----------------------

    call ESMF_ClockGet( clock, &
         currTime=currTime, startTime=startTime, stopTime=stopTime, refTime=RefTime, timeStep=timeStep, rc=rc)
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

    call timemgr_setup(&
         calendar_in=calendar, &
         start_ymd_in=start_ymd, &
         start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, &
         ref_tod_in=ref_tod, &
         stop_ymd_in=stop_ymd, &
         stop_tod_in=stop_tod)

    !--------------------------------
    ! read in lilac_in namelists 
    !--------------------------------

    if (masterproc) then
       open(newunit=fileunit, status="old", file="lilac_in")
       call shr_nl_find_group_name(fileunit, 'lilac_run_input', ierr)
       if (ierr == 0) then
          read(fileunit, lilac_run_input, iostat=ierr)
          if (ierr > 0) then
             call shr_sys_abort( 'problem on read of lilac_run_input')
          end if
       end if
       call shr_nl_find_group_name(fileunit, 'lilac_rof_input', ierr)
       if (ierr == 0) then
          read(fileunit, lilac_rof_input, iostat=ierr)
          if (ierr > 0) then
             call shr_sys_abort( 'problem on read of lilac_rof_input')
          end if
       end if
       close(fileunit)
    end if
    call mpi_bcast(rof_mesh_filename, len(rof_mesh_filename), MPI_CHARACTER, 0, mpicom_rof, ierr)
    call mpi_bcast(starttype, len(starttype), MPI_CHARACTER, 0, mpicom_rof, ierr)
    call mpi_bcast(caseid, len(caseid), MPI_CHARACTER, 0, mpicom_rof, ierr)

    !--------------------------------
    ! Initialize RtmVar module variables
    !--------------------------------

    if (trim(starttype) == trim('startup')) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim('continue') ) then
       nsrest = nsrContinue
    else
       call shr_sys_abort( subname//' ERROR: unknown starttype' )
    end if

    !----------------------
    ! Read namelist, grid and surface data
    !----------------------

    if (masterproc) then
       write(iulog,*) "MOSART river model initialization"
       write(iulog,*) ' mosart npes   = ',npes
       write(iulog,*) ' mosart caseid = ',trim(caseid)
       write(iulog,*) ' mosart nsrest = ',nsrest
    endif

    call RtmVarSet(caseid_in=trim(caseid), ctitle_in=trim(caseid), nsrest_in=nsrest)

    !----------------------
    ! Initialize Mosart
    !----------------------

    ! - Read in mosart namelist
    ! - Initialize mosart time manager
    ! - Initialize number of mosart tracers
    ! - Read input data (river direction file) (global)
    ! - Deriver gridbox edges (global)
    ! - Determine mosart ocn/land mask (global)
    ! - Compute total number of basins and runoff ponts
    ! - Compute river basins, actually compute ocean outlet gridcell
    ! - Allocate basins to pes
    ! - Count and distribute cells to rglo2gdc (determine rtmCTL%begr, rtmCTL%endr)
    ! - Adjust area estimation from DRT algorithm for those outlet grids
    !     - useful for grid-based representation only
    !     - need to compute areas where they are not defined in input file
    ! - Initialize runoff datatype (rtmCTL)

    ! TODO: are not handling rof_prognostic = .false. for now

    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet(timeStep, s=dtime_lilac, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call Rtmini(rtm_active=rof_prognostic, flood_active=flood_present, dtime_driver=dtime_lilac)

    !--------------------------------
    ! generate the mesh and realize fields
    !--------------------------------

    ! determine global index array
    lsize = rtmCTL%endr - rtmCTL%begr + 1
    allocate(gindex(lsize))
    ni = 0
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       gindex(ni) = rtmCTL%gindex(n)
    end do

    ! create distGrid from global index array
    DistGrid = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(gindex)

    ! create esmf mesh using distgrid and rof_mesh_filename
    rof_mesh = ESMF_MeshCreate(filename=trim(rof_mesh_filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, elementDistgrid=Distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    if (masterproc) then
       write(iulog,*)'mesh file for domain is ',trim(rof_mesh_filename)
    end if
    call ESMF_LogWrite(subname//" Create Mesh using file ...."//trim(rof_mesh_filename), ESMF_LOGMSG_INFO)

    !--------------------------------
    ! Create mosart import state
    !--------------------------------

    ! create an empty field bundle of fields received from land
    c2r_fb =  ESMF_FieldBundleCreate (name='c2r_fb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! now add fields on lnd_mesh to this field bundle
    call fldbundle_add('Flrl_rofsur'    , c2r_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofgwl'    , c2r_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofsub'    , c2r_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rofi'      , c2r_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrl_rof_irrig' , c2r_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add the field bundle to the import
    call ESMF_StateAdd(importState, fieldbundleList = (/c2r_fb/))

    !--------------------------------
    ! Create mosart export state
    !--------------------------------

    ! create an empty field bundle of field sent to alnd
    r2c_fb = ESMF_FieldBundleCreate(name='r2c_fb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! create an empty field bundle for the import of rof fields
    r2c_fb = ESMF_FieldBundleCreate (name='r2c_fb', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call fldbundle_add('Flrr_flood', r2c_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrr_volr', r2c_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call fldbundle_add('Flrr_volrmch', r2c_fb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add the field bundle to the state
    call ESMF_StateAdd(exportState, fieldbundleList = (/r2c_fb/))

    !--------------------------------
    ! fill in mosart export state
    !--------------------------------

    call export_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add global size attributes to rof gridded component
    call ESMF_AttributeAdd(gcomp, convention="custom", purpose="global grid sizes", rc=rc)
    write(cvalue,*) dble(rtmlon)
    call ESMF_AttributeSet(gcomp, "global_nx", cvalue, rc=rc)
    write(cvalue,*) dble(rtmlat)
    call ESMF_AttributeSet(gcomp, "global_ny", cvalue, rc=rc)

    !--------------------------------
    ! diagnostics
    !--------------------------------

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    if (debug > 1) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_comp_esmf::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  !---------------------------
  contains
  !---------------------------

    subroutine fldbundle_add(stdname, fieldbundle, rc)
      !---------------------------
      ! Create an empty input field with name 'stdname' to add to fieldbundle
      !---------------------------

      ! input/output variables
      character(len=*)        , intent(in)    :: stdname
      type (ESMF_FieldBundle) , intent(inout) :: fieldbundle
      integer                 , intent(out)   :: rc
      ! local variables
      type(ESMF_Field) :: field
      !-------------------------------------------------------------------------------
      rc = ESMF_SUCCESS
      field = ESMF_FieldCreate(rof_mesh, ESMF_TYPEKIND_R8 , meshloc=ESMF_MESHLOC_ELEMENT , name=trim(stdname), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldBundleAdd(fieldbundle, (/field/), rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end subroutine fldbundle_add

  end subroutine rof_init

  !===============================================================================

  subroutine rof_run(gcomp, import_state, export_state, clock, rc)

    !------------------------
    ! Run MOSART
    !------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp           ! MOSART gridded component
    type(ESMF_State)     :: import_state    ! MOSART import state
    type(ESMF_State)     :: export_state    ! MOSART export state
    type(ESMF_Clock)     :: clock           ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code

    ! local variables:
    type(ESMF_Alarm)  :: alarm
    type(ESMF_Time)   :: currTime
    type(ESMF_Time)   :: nextTime
    character(CL)     :: cvalue
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
    character(len=*),parameter  :: subname=trim(modName)//':(rof_run) '
    !-------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','runf_run:start::',lbnum)
    endif
#endif

    !--------------------------------
    ! Unpack import state 
    !--------------------------------

    call t_startf ('lc_mosart_import')
    call import_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_mosart_import')

    !--------------------------------
    ! Determine if time to write restart
    !--------------------------------

    call ESMF_ClockGetAlarm(clock, alarmname='lilac_restart_alarm', alarm=alarm, rc=rc)
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

    call ESMF_ClockGetAlarm(clock, alarmname='lilac_stop_alarm', alarm=alarm, rc=rc)
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

    ! Advance mosart time
    call advance_timestep()

    ! Run MOSART (export data is in rtmCTL and Trunoff data types)
    call Rtmrun(rstwr, nlend, rdate)

    !--------------------------------
    ! Pack export state to mediator
    !--------------------------------

    ! (input is rtmCTL%runoff, output is r2x)
    call t_startf ('lc_rof_export')
    call export_fields(gcomp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf ('lc_rof_export')

    !--------------------------------
    ! Check that internal clock is in sync with master clock
    !--------------------------------

    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    tod = tod

    if ( (ymd /= ymd_sync) .and. (tod /= tod_sync) ) then
       write(iulog,*)' mosart ymd=',ymd     ,'  mosart tod= ',tod
       write(iulog,*)'   sync ymd=',ymd_sync,'    sync tod= ',tod_sync
       rc = ESMF_FAILURE
       call ESMF_LogWrite(subname//" MOSART clock not in sync with Master Sync clock",ESMF_LOGMSG_ERROR)
    end if

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug > 1) then
       call State_diagnose(export_state,subname//':ES',rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    end if

    if (masterproc) then
       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ROF from: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    endif

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_run:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine rof_run

  !===============================================================================

  subroutine rof_final(gcomp, import_state, export_state, clock, rc)

    !---------------------------------
    ! Finalize MOSART
    !---------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp        ! MOSART gridded component
    type(ESMF_State)     :: import_state ! MOSART import state
    type(ESMF_State)     :: export_state ! MOSART export state
    type(ESMF_Clock)     :: clock        ! ESMF synchronization clock
    integer, intent(out) :: rc           ! Return code

    ! local variables
    character(*), parameter :: F00   = "('(rof_final) ',8a)"
    character(*), parameter :: F91   = "('(rof_final) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    if (masterproc) then
       write(iulog,F91)
       write(iulog,F00) 'MOSART: end of main integration loop'
       write(iulog,F91)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine rof_final

end module rof_comp_esmf
