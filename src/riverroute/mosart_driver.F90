module mosart_driver

   !-----------------------------------------------------------------------
   ! Mosart Routing Model
   !-----------------------------------------------------------------------

   use shr_kind_mod       , only : r8 => shr_kind_r8, CS => shr_kind_cs, CL => shr_kind_CL
   use shr_sys_mod        , only : shr_sys_abort
   use shr_const_mod      , only : SHR_CONST_PI, SHR_CONST_CDAY
   use mosart_vars        , only : re, spval, iulog, ice_runoff, &
                                   frivinp, nsrContinue, nsrBranch, nsrStartup, nsrest, &
                                   inst_index, inst_suffix, inst_name, decomp_option, &
                                   bypass_routing_option, qgwl_runoff_option, barrier_timers, &
                                   mainproc, npes, iam, mpicom_rof, budget_frq, isecspday
   use mosart_data        , only : ctl, Tctl, Tunit, TRunoff, Tpara
   use mosart_budget_type , only : budget_type
   use mosart_fileutils   , only : getfil
   use mosart_timemanager , only : timemgr_init, get_nstep, get_curr_date
   use mosart_histflds    , only : mosart_histflds_init, mosart_histflds_set
   use mosart_histfile    , only : mosart_hist_updatehbuf, mosart_hist_htapeswrapup, mosart_hist_htapesbuild, &
                                   ndens, mfilt, nhtfrq, avgflag_pertape, avgflag_pertape, &
                                   fincl1, fincl2, fincl3, fexcl1, fexcl2, fexcl3, max_tapes, max_namlen
   use mosart_restfile    , only : mosart_rest_timemanager, mosart_rest_getfile, mosart_rest_fileread, &
                                   mosart_rest_filewrite, mosart_rest_filename, finidat, nrevsn
   use mosart_physics     , only : updatestate_hillslope, updatestate_subnetwork, updatestate_mainchannel, Euler
   use perf_mod           , only : t_startf, t_stopf
   use nuopc_shr_methods  , only : chkerr
   use ESMF               , only : ESMF_SUCCESS, ESMF_FieldGet, ESMF_FieldSMMStore, ESMF_FieldSMM, &
                                   ESMF_TERMORDER_SRCSEQ, ESMF_Mesh, ESMF_Time
   use mosart_io          , only : ncd_pio_openfile, ncd_inqdid, ncd_inqdlen, ncd_pio_closefile, ncd_decomp_init, &
                                   pio_subsystem
   use pio                , only : file_desc_t
   use mpi

   implicit none
   private

   ! public member functions:
   public :: mosart_read_namelist ! Read in mosart namelist
   public :: mosart_init1         ! Initialize mosart grid
   public :: mosart_init2         ! Initialize mosart maps
   public :: mosart_run           ! River routing model

   ! mosart namelists
   integer           :: coupling_period         ! mosart coupling period
   integer           :: delt_mosart             ! mosart internal timestep (->nsub)
   logical           :: use_halo_option         ! enable halo capability using ESMF
   character(len=CS) :: mosart_tracers          ! colon delimited string of tracer names
   character(len=CS) :: mosart_euler_calc       ! colon delimited string of logicals for using Euler  algorithm

   ! subcycling
   integer   :: nsub_save ! previous nsub
   real(r8)  :: delt_save ! previous delt

   ! global (glo)
   integer , allocatable :: IDkey(:) ! translation key from ID to gindex

   ! budget
   type(budget_type), public :: budget  ! type containing vars and routines for budget checking

   character(len=CL) :: nlfilename_rof = 'mosart_in'
   character(len=CL) :: fnamer              ! name of netcdf restart file

   integer :: nt_liq, nt_ice

   character(*), parameter :: u_FILE_u = &
        __FILE__
   !-----------------------------------------------------------------------

contains

   !-----------------------------------------------------------------------
   subroutine mosart_read_namelist()
      !
      ! Read and distribute mosart namelist
      !
      ! local variables
      integer           :: i
      integer           :: ier       ! error code
      integer           :: unitn     ! unit for namelist file
      logical           :: lexist    ! File exists
      character(len=CS) :: runtyp(4) ! run type
      logical, allocatable :: do_euler_calc(:) ! turn on euler algorithm
      character(len=*),parameter :: subname = '(mosart_read_namelist) '
      !-----------------------------------------------------------------------

      !-------------------------------------------------------
      ! Read in mosart namelist
      !-------------------------------------------------------

      namelist /mosart_inparm / frivinp, finidat, nrevsn, coupling_period, ice_runoff, &
           ndens, mfilt, nhtfrq, fincl1,  fincl2, fincl3, fexcl1,  fexcl2, fexcl3, &
           avgflag_pertape, decomp_option, bypass_routing_option, qgwl_runoff_option, &
           use_halo_option, delt_mosart, mosart_tracers, mosart_euler_calc, budget_frq

      ! Preset values
      ice_runoff  = .true.
      finidat = ' '
      nrevsn  = ' '
      coupling_period   = -1
      delt_mosart = 3600
      decomp_option = 'basin'
      bypass_routing_option = 'direct_in_place'
      qgwl_runoff_option = 'threshold'
      use_halo_option = .false.
      mosart_tracers = 'LIQ:ICE'
      mosart_euler_calc = 'T:F'

      nlfilename_rof = "mosart_in" // trim(inst_suffix)
      inquire (file = trim(nlfilename_rof), exist = lexist)
      if ( .not. lexist ) then
         write(iulog,*) subname // ' ERROR: nlfilename_rof does NOT exist: '//trim(nlfilename_rof)
         call shr_sys_abort(trim(subname)//' ERROR nlfilename_rof does not exist')
      end if
      if (mainproc) then
         write(iulog,*) 'Reading mosart_inparm namelist from: ', trim(nlfilename_rof)
         open( newunit=unitn, file=trim(nlfilename_rof), status='old' )
         ier = 1
         do while ( ier /= 0 )
            read(unitn, mosart_inparm, iostat=ier)
            if (ier < 0) then
               call shr_sys_abort( subname//' encountered end-of-file on mosart_inparm read' )
            endif
         end do
         close(unitn)
      end if

      call mpi_bcast (finidat, len(finidat), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (frivinp, len(frivinp), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (nrevsn, len(nrevsn), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (decomp_option, len(decomp_option), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (use_halo_option, 1, MPI_LOGICAL, 0, mpicom_rof, ier)
      call mpi_bcast (coupling_period, 1, MPI_INTEGER, 0, mpicom_rof, ier)
      call mpi_bcast (delt_mosart, 1, MPI_INTEGER, 0, mpicom_rof, ier)
      call mpi_bcast (bypass_routing_option, len(bypass_routing_option), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (qgwl_runoff_option, len(qgwl_runoff_option), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (ice_runoff, 1, MPI_LOGICAL, 0, mpicom_rof, ier)
      call mpi_bcast (nhtfrq, size(nhtfrq), MPI_INTEGER, 0, mpicom_rof, ier)
      call mpi_bcast (mfilt, size(mfilt), MPI_INTEGER, 0, mpicom_rof, ier)
      call mpi_bcast (ndens, size(ndens), MPI_INTEGER, 0, mpicom_rof, ier)
      call mpi_bcast (fexcl1, (max_namlen+2)*size(fexcl1), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (fexcl2, (max_namlen+2)*size(fexcl2), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (fexcl3, (max_namlen+2)*size(fexcl3), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (fincl1, (max_namlen+2)*size(fincl1), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (fincl2, (max_namlen+2)*size(fincl2), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (fincl3, (max_namlen+2)*size(fincl3), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (avgflag_pertape, size(avgflag_pertape), MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (mosart_tracers, CS, MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (mosart_euler_calc, CS, MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (budget_frq,1,MPI_INTEGER,0,mpicom_rof,ier)

      ! Determine number of tracers and array of tracer names and initialize module variables
      call ctl%init_tracer_names(mosart_tracers)
      nt_liq = ctl%nt_liq
      nt_ice = ctl%nt_ice

      runtyp(:)               = 'missing'
      runtyp(nsrStartup  + 1) = 'initial'
      runtyp(nsrContinue + 1) = 'restart'
      runtyp(nsrBranch   + 1) = 'branch '

      if (mainproc) then
         write(iulog,*) 'define run:'
         write(iulog,'(a)'   ) '   run type                = '//trim(runtyp(nsrest+1))
         write(iulog,'(a,i8)') '   coupling_period         = ',coupling_period
         write(iulog,'(a,i8)') '   delt_mosart             = ',delt_mosart
         write(iulog,'(a)'   ) '   decomp option           = '//trim(decomp_option)
         write(iulog,'(a,l1)') '   use_halo_option         = ',use_halo_option
         write(iulog,'(a)'   ) '   bypass_routing option   = '//trim(bypass_routing_option)
         write(iulog,'(a)'   ) '   qgwl runoff option      = '//trim(qgwl_runoff_option)
         write(iulog,'(a)'   ) '   mosart tracers          = '//trim(mosart_tracers)
         write(iulog,'(a)'   ) '   mosart euler calc       = '//trim(mosart_euler_calc)
         if (nsrest == nsrStartup .and. finidat /= ' ') then
           write(iulog,'(a)') '   mosart initial data     = '//trim(finidat)
         end if
      endif

      if (frivinp == ' ') then
         call shr_sys_abort( subname//' ERROR: frivinp NOT set' )
      else
         if (mainproc) then
            write(iulog,*) '   mosart river data       = ',trim(frivinp)
         endif
      end if

      if (trim(bypass_routing_option) == 'direct_to_outlet') then
         if (trim(qgwl_runoff_option) == 'threshold') then
            call shr_sys_abort( subname//' ERROR: qgwl_runoff_option &
                 CANNOT be threshold if bypass_routing_option==direct_to_outlet' )
         end if
      else if (trim(bypass_routing_option) == 'none') then
         if (trim(qgwl_runoff_option) /= 'all') then
            call shr_sys_abort( subname//' ERROR: qgwl_runoff_option &
                 can only be all if bypass_routing_option==none' )
         end if
      end if

      if (coupling_period <= 0) then
         write(iulog,*) subname,' ERROR mosart coupling_period invalid',coupling_period
         call shr_sys_abort( subname//' ERROR: coupling_period invalid' )
      endif

      if (delt_mosart <= 0) then
         write(iulog,*) subname,' ERROR mosart delt_mosart invalid',delt_mosart
         call shr_sys_abort( subname//' ERROR: delt_mosart invalid' )
      endif

      do i = 1, max_tapes
         if (nhtfrq(i) == 0) then
            mfilt(i) = 1
         else if (nhtfrq(i) < 0) then
            nhtfrq(i) = nint(-nhtfrq(i)*SHR_CONST_CDAY/(24._r8*coupling_period))
         endif
      end do

   end subroutine mosart_read_namelist

   !-----------------------------------------------------------------------

   subroutine mosart_init1(currTime, rc)

      !-------------------------------------------------
      ! Initialize mosart grid, mask, decomp
      !
      ! Arguments
      type(ESMF_Time), intent(in) :: currTime
      integer, intent(out) :: rc
      !
      ! Local variables
      integer            :: n, nr, nt         ! indices
      type(file_desc_t)  :: ncid              ! netcdf file id
      character(len=CL)  :: trstr             ! tracer string
      character(len=CL)  :: locfn             ! local file
      integer            :: dimid             ! netcdf dimension identifier
      character(len=*), parameter :: subname = '(mosart_init1) '
      !-------------------------------------------------

      rc = ESMF_SUCCESS
      call timemgr_init(dtime_in=coupling_period, curr_date=currTime)

      !-------------------------------------------------------
      ! Obtain restart file if appropriate
      !-------------------------------------------------------
      if ((nsrest == nsrStartup .and. finidat /= ' ') .or. &
          (nsrest == nsrContinue) .or. (nsrest == nsrBranch  )) then
         call mosart_rest_getfile( file=fnamer )
      endif

      !-------------------------------------------------------
      ! Write out tracers to stdout
      !-------------------------------------------------------
      if (mainproc) then
         trstr = trim(ctl%tracer_names(1))
         do n = 2,ctl%ntracers
            trstr = trim(trstr)//':'//trim(ctl%tracer_names(n))
         enddo
         write(iulog,*)'mosart tracers = ',ctl%ntracers,trim(trstr)
      end if

      !-------------------------------------------------------
      ! Obtain global sizes of grid from river direction file
      !-------------------------------------------------------
      call getfil(frivinp, locfn, 0 )
      call ncd_pio_openfile(ncid, trim(locfn), 0)
      call ncd_inqdid(ncid,'lon',dimid)
      call ncd_inqdlen(ncid,dimid,ctl%nlon)
      call ncd_inqdid(ncid,'lat',dimid)
      call ncd_inqdlen(ncid,dimid,ctl%nlat)
      call ncd_pio_closefile(ncid)
      if (mainproc) then
         write(iulog,'(a)') 'MOSART river data file name: ',trim(frivinp)
         write(iulog,'(a)') 'Successfully read mosart dimensions'
         write(iulog,'(a,i8,2x,i8)') 'Values for global nlon/nlat: ',ctl%nlon,ctl%nlat
      endif

      !-------------------------------------------------------
      ! Initialize ctl derived type allocatable variables
      !-------------------------------------------------------
      allocate(IDkey(ctl%nlon*ctl%nlat))
      call ctl%Init(locfn, decomp_option, use_halo_option, IDkey, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      !-------------------------------------------------------
      ! Initialize pio compDOF (module variable in mosart_io)
      !-------------------------------------------------------
      call ncd_decomp_init(ctl%begr, ctl%endr, ctl%numr, ctl%gindex)

   end subroutine mosart_init1

   !-----------------------------------------------------------------------

   subroutine mosart_init2(Emesh, rc)

      ! Second phyas of mosart initialization
      !
      ! Arguments
      type(ESMF_Mesh), intent(in)  :: Emesh
      integer        , intent(out) :: rc
      !
      ! Local variables
      integer :: nr, nt
      integer :: begr, endr
      integer :: ntracers
      character(len=*),parameter :: subname = '(mosart_init2)'
      !-----------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! Set up local variables to be used below
      begr  = ctl%begr
      endr  = ctl%endr
      ntracers = ctl%ntracers

      !-------------------------------------------------------
      ! Initialize MOSART types TCtl, Tpara, TUnit and Trunoff
      !-------------------------------------------------------

      call Tctl%Init()

      call Tpara%Init(begr, endr)

      call TRunoff%Init(begr, endr, ntracers)

      call Tunit%Init(begr, endr, ntracers, &
           mosart_euler_calc, ctl%nlon, ctl%nlat, Emesh, trim(frivinp), IDKey, &
           Tpara%c_twid, Tctl%DLevelR, ctl%area, ctl%gindex, ctl%outletg, pio_subsystem, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      !-------------------------------------------------------
      ! Read restart/initial info
      !-------------------------------------------------------

      call t_startf('mosarti_restart')
      if ((nsrest == nsrStartup .and. finidat  /= ' ') .or. &
          (nsrest == nsrContinue) .or. &
          (nsrest == nsrBranch  )) then
         call mosart_rest_fileread( file=fnamer )
      endif

      do nt = 1,ntracers
         do nr = begr,endr
            call UpdateState_hillslope(nr,nt)
            call UpdateState_subnetwork(nr,nt)
            call UpdateState_mainchannel(nr,nt)
            ctl%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + TRunoff%wh(nr,nt)*ctl%area(nr))
         enddo
      enddo
      call t_stopf('mosarti_restart')

      !-------------------------------------------------------
      ! Initialize mosart history handler and fields
      !-------------------------------------------------------

      call t_startf('mosarti_histinit')
      call mosart_histflds_init(begr, endr, ntracers)
      if (nsrest==nsrStartup .or. nsrest==nsrBranch) then
         call mosart_hist_HtapesBuild()
      end if
      call mosart_histflds_set(ntracers)
      if (mainproc) write(iulog,*) subname,' done'
      call t_stopf('mosarti_histinit')

      !-------------------------------------------------------
      ! Initialize mosart budget
      !-------------------------------------------------------

      call t_startf('mosarti_budgetinit')
      call budget%Init(begr, endr, ntracers)
      call t_stopf('mosarti_budgetinit')

   end subroutine mosart_init2

   !-----------------------------------------------------------------------

   subroutine mosart_run(begr, endr, ntracers, rstwr, nlend, rdate, rc)

      ! Run mosart river routing model
      !
      ! Arguments
      integer          , intent(in)  :: begr, endr, ntracers
      logical          , intent(in)  :: rstwr ! true => write restart file this step)
      logical          , intent(in)  :: nlend ! true => end of run on this step
      character(len=*) , intent(in)  :: rdate ! restart file time stamp for name
      integer          , intent(out) :: rc
      !
      ! Local variables
      integer            :: i, j, n, nr, ns, nt, n2, nf ! indices
      logical            :: budget_check                ! if budget check needs to be performed
      real(r8)           :: volr_init                   ! temporary storage to compute dvolrdt
      integer            :: yr, mon, day, ymd, tod      ! time information
      integer            :: nsub                        ! subcyling for cfl
      real(r8)           :: delt                        ! delt associated with subcycling
      real(r8)           :: delt_coupling               ! real value of coupling_period
      character(len=CL)  :: filer                       ! restart file name
      integer            :: cnt                         ! counter for gridcells
      integer            :: ier                         ! error code
      real(r8), pointer  :: src_direct(:,:)
      real(r8), pointer  :: dst_direct(:,:)

      ! parameters used in negative runoff partitioning algorithm
      real(r8) :: river_depth_minimum = 1.e-4 ! gridcell average minimum river depth [m]
      real(r8) :: river_volume_minimum        ! gridcell area multiplied by average river_depth_minimum [m3]
      real(r8) :: qgwl_volume                 ! volume of runoff during time step [m3]
      real(r8) :: irrig_volume                ! volume of irrigation demand during time step [m3]
      logical, save :: first_call = .true.    ! first time flag (for backwards compatibility)
      character(len=*),parameter :: subname = ' (mosart_run) '
      !-----------------------------------------------------------------------

      call t_startf('mosartr_tot')

      rc = ESMF_SUCCESS

      !-----------------------------------------------------
      ! Get date info
      !-----------------------------------------------------

      call get_curr_date(yr, mon, day, tod)
      ymd = yr*10000 + mon*100 + day
      if (tod == 0) then
         if (mainproc) then
            write(iulog,*) ' '
            write(iulog,'(2a,i10,i6)') trim(subname),' model date is',ymd,tod
         end if
      endif

      delt_coupling = coupling_period*1.0_r8

      if (first_call) then
         delt_save = delt_mosart
         if (mainproc) then
            write(iulog,'(2a,g20.12)') trim(subname),' mosart coupling period ',delt_coupling
         end if
      end if


      ! BUDGET

      budget_check = .false.
      if (budget_frq == 0) then
        if (day == 1 .and. tod == 0) then
          budget_check = .true.
        endif
      else if (budget_frq < 0) then
        if (mod(get_nstep() * coupling_period, abs(budget_frq) * 3600) == 0) then
          budget_check = .true.
        endif
      else
        if (mod(get_nstep() , budget_frq) == 0) then
          budget_check = .true.
        endif
      endif
      if (first_call) then ! ignore budget during the first timestep
        budget_check = .false.
      endif
      if (budget_check) then
        call t_startf('mosartr_budgetset')
        call  budget%set_budget(begr,endr,ntracers, delt_coupling)
        call t_stopf('mosartr_budgetset')
      endif

      ! initialize data for euler solver, in m3/s here
      do nr = begr,endr
         do nt = 1,ntracers
            TRunoff%qsur(nr,nt) = ctl%qsur(nr,nt)
            TRunoff%qsub(nr,nt) = ctl%qsub(nr,nt)
            TRunoff%qgwl(nr,nt) = ctl%qgwl(nr,nt)
         enddo
      enddo

      !-----------------------------------
      ! Compute irrigation flux based on demand from clm
      ! Must be calculated before volr is updated to be consistent with lnd
      ! Just consider land points and only remove liquid water
      !-----------------------------------

      call t_startf('mosartr_irrig')
      ctl%qirrig_actual = 0._r8
      do nr = begr,endr

         ! calculate volume of irrigation flux during timestep
         irrig_volume = -ctl%qirrig(nr) * coupling_period

         ! compare irrig_volume to main channel storage;
         ! add overage to subsurface runoff
         if(irrig_volume > TRunoff%wr(nr,nt_liq)) then
            ctl%qsub(nr,nt_liq) = ctl%qsub(nr,nt_liq) + (TRunoff%wr(nr,nt_liq) - irrig_volume) / coupling_period
            TRunoff%qsub(nr,nt_liq) = ctl%qsub(nr,nt_liq)
            irrig_volume = TRunoff%wr(nr,nt_liq)
         endif

         ! actual irrigation rate [m3/s]
         ! i.e. the rate actually removed from the main channel
         ! if irrig_volume is greater than TRunoff%wr
         ctl%qirrig_actual(nr) = - irrig_volume / coupling_period

         ! remove irrigation from wr (main channel)
         TRunoff%wr(nr,nt_liq) = TRunoff%wr(nr,nt_liq) - irrig_volume

      enddo
      call t_stopf('mosartr_irrig')

      !-----------------------------------
      ! Compute flood
      ! Remove water from mosart and send back to clm
      ! Just consider land points and only remove liquid water
      ! ctl%flood is m3/s here
      !-----------------------------------

      call t_startf('mosartr_flood')
      ctl%flood = 0._r8
      do nr = begr,endr
         ! initialize ctl%flood to zero
         if (ctl%mask(nr) == 1) then
            if (ctl%volr(nr,nt_liq) > ctl%fthresh(nr)) then
               ! determine flux that is sent back to the land this is in m3/s
               ctl%flood(nr) = (ctl%volr(nr,nt_liq)-ctl%fthresh(nr)) / (delt_coupling)

               ! ctl%flood will be sent back to land - so must subtract this
               ! from the input runoff from land
               ! tcraig, comment - this seems like an odd approach, you
               !   might create negative forcing.  why not take it out of
               !   the volr directly?  it's also odd to compute this
               !   at the initial time of the time loop.  why not do
               !   it at the end or even during the run loop as the
               !   new volume is computed.  fluxout depends on volr, so
               !   how this is implemented does impact the solution.
               TRunoff%qsur(nr,nt_liq) = TRunoff%qsur(nr,nt_liq) - ctl%flood(nr)
            endif
         endif
      enddo
      call t_stopf('mosartr_flood')

      !-----------------------------------------------------
      ! DIRECT transfer to outlet point
      ! Remember to subtract water from TRunoff forcing
      !-----------------------------------------------------

      if (barrier_timers) then
         call t_startf('mosartr_SMdirect_barrier')
         call mpi_barrier(mpicom_rof,ier)
         call t_stopf ('mosartr_SMdirect_barrier')
      endif

      call t_startf('mosartr_SMdirect')

      !-----------------------------------------------------
      ! Set up pointer arrays into srcfield and dstfield
      !-----------------------------------------------------

      call ESMF_FieldGet(Tunit%srcfield, farrayPtr=src_direct, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(Tunit%dstfield, farrayPtr=dst_direct, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      !-----------------------------------------------------
      !--- initialize ctl%direct
      !-----------------------------------------------------

      ctl%direct(:,:) = 0._r8

      !-----------------------------------------------------
      !--- direct to outlet: all liquid and frozen runoff from glc
      !-----------------------------------------------------

      if (ctl%rof_from_glc) then
        src_direct(:,:) = 0._r8
        dst_direct(:,:) = 0._r8

        cnt = 0
        do nr = begr,endr
          cnt = cnt + 1
          src_direct(nt_liq,cnt) = ctl%qglc_liq(nr)
          src_direct(nt_ice,cnt) = ctl%qglc_ice(nr)
        enddo

        call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
        if (chkerr(rc,__LINE__,u_FILE_u)) return

        ! copy direct transfer water to output field
        cnt = 0
        do nr = begr,endr
          cnt = cnt + 1
          ctl%direct_glc(nr,nt_liq) = dst_direct(nt_liq,cnt)
          ctl%direct_glc(nr,nt_ice) = dst_direct(nt_ice,cnt)
        enddo
      else
        ctl%direct_glc(:,:) = 0._r8
        ctl%direct_glc(:,:) = 0._r8
      end if

      !-----------------------------------------------------
      !--- direct to outlet: all frozen runoff from lnd
      !-----------------------------------------------------

      src_direct(:,:) = 0._r8
      dst_direct(:,:) = 0._r8

      cnt = 0
      do nr = begr,endr
         cnt = cnt + 1
         src_direct(nt_ice,cnt) = TRunoff%qsur(nr,nt_ice) + TRunoff%qsub(nr,nt_ice) + TRunoff%qgwl(nr,nt_ice)
      enddo

      call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! copy direct transfer water to output field
      cnt = 0
      do nr = begr,endr
         cnt = cnt + 1
         ctl%direct(nr,nt_ice) = ctl%direct(nr,nt_ice) + dst_direct(nt_ice,cnt)
      enddo

      ! set euler_calc = false for frozen runoff
      ! TODO: will be reworked after addition of multiple tracers
      Tunit%euler_calc(nt_ice) = .false.

      ! Set Trunoff%qsur, TRunoff%qsub and Trunoff%qgwl to zero for nt_ice
      TRunoff%qsur(:,nt_ice) = 0._r8
      TRunoff%qsub(:,nt_ice) = 0._r8
      TRunoff%qgwl(:,nt_ice) = 0._r8

      !-----------------------------------------------------
      !--- direct to outlet: qgwl
      !-----------------------------------------------------

      !-- liquid runoff components
      if (trim(bypass_routing_option) == 'direct_to_outlet') then

         src_direct(:,:) = 0._r8
         dst_direct(:,:) = 0._r8

         !--- copy direct transfer fields, convert kg/m2s to m3/s
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            if (trim(qgwl_runoff_option) == 'all') then
               src_direct(nt_liq,cnt) = TRunoff%qgwl(nr,nt_liq)
               TRunoff%qgwl(nr,nt_liq) = 0._r8
            else if (trim(qgwl_runoff_option) == 'negative') then
               if(TRunoff%qgwl(nr,nt_liq) < 0._r8) then
                  src_direct(nt_liq,cnt) = TRunoff%qgwl(nr,nt_liq)
                  TRunoff%qgwl(nr,nt_liq) = 0._r8
               endif
            endif
         enddo

         call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         !--- copy direct transfer water to output field ---
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            ctl%direct(nr,nt_liq) = ctl%direct(nr,nt_liq) + dst_direct(nt_liq,cnt)
         enddo
      endif

      !-----------------------------------------------------
      !--- direct in place qgwl, qgwl
      !-----------------------------------------------------

      if (trim(bypass_routing_option) == 'direct_in_place') then
         do nr = begr,endr
            if (trim(qgwl_runoff_option) == 'all') then
               ctl%direct(nr,nt_liq) = TRunoff%qgwl(nr,nt_liq)
               TRunoff%qgwl(nr,nt_liq) = 0._r8
            else if (trim(qgwl_runoff_option) == 'negative') then
               if(TRunoff%qgwl(nr,nt_liq) < 0._r8) then
                  ctl%direct(nr,nt_liq) = TRunoff%qgwl(nr,nt_liq)
                  TRunoff%qgwl(nr,nt_liq) = 0._r8
               endif
            else if (trim(qgwl_runoff_option) == 'threshold') then
               ! --- calculate volume of qgwl flux during timestep
               qgwl_volume = TRunoff%qgwl(nr,nt_liq) * ctl%area(nr) * coupling_period
               river_volume_minimum = river_depth_minimum * ctl%area(nr)

               ! if qgwl is negative, and adding it to the main channel
               ! would bring main channel storage below a threshold,
               ! send qgwl directly to ocean
               if (((qgwl_volume + TRunoff%wr(nr,nt_liq)) < river_volume_minimum) .and. (TRunoff%qgwl(nr,nt_liq) < 0._r8)) then
                  ctl%direct(nr,nt_liq) = TRunoff%qgwl(nr,nt_liq)
                  TRunoff%qgwl(nr,nt_liq) = 0._r8
               endif
            endif
         enddo
      endif

      !-------------------------------------------------------
      !--- direct in place: add other direct terms, e.g. inputs outside of mosart mask, negative qsur
      !-------------------------------------------------------

      if (trim(bypass_routing_option) == 'direct_in_place') then
         do nt = 1,ntracers
            do nr = begr,endr
               if (TRunoff%qsub(nr,nt) < 0._r8) then
                  ctl%direct(nr,nt) = ctl%direct(nr,nt) + TRunoff%qsub(nr,nt)
                  TRunoff%qsub(nr,nt) = 0._r8
               endif
               if (TRunoff%qsur(nr,nt) < 0._r8) then
                  ctl%direct(nr,nt) = ctl%direct(nr,nt) + TRunoff%qsur(nr,nt)
                  TRunoff%qsur(nr,nt) = 0._r8
               endif
               ! Note Tunit%mask is set in Tunit%init and is obtained from reading in fdir
               ! if fdir<0 then mask=0 (ocean), if fdir=0 then mask=2 (outlet) and if fdir>0 then mask=1 (land)
               if (Tunit%mask(nr) > 0) then
                  ! mosart euler
               else
                  ctl%direct(nr,nt) = ctl%direct(nr,nt) + TRunoff%qsub(nr,nt) + TRunoff%qsur(nr,nt) + TRunoff%qgwl(nr,nt)
                  TRunoff%qsub(nr,nt) = 0._r8
                  TRunoff%qsur(nr,nt) = 0._r8
                  TRunoff%qgwl(nr,nt) = 0._r8
               endif
            enddo
         enddo
      endif

      !-------------------------------------------------------
      !--- direct to outlet: add other direct terms, e.g. inputs outside of mosart mask, negative qsur
      !-------------------------------------------------------

      if (trim(bypass_routing_option) == 'direct_to_outlet') then
         src_direct(:,:) = 0._r8
         dst_direct(:,:) = 0._r8
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            do nt = 1,ntracers
               !---- negative qsub water, remove from TRunoff ---
               if (TRunoff%qsub(nr,nt) < 0._r8) then
                  src_direct(nt,cnt) = src_direct(nt,cnt) + TRunoff%qsub(nr,nt)
                  TRunoff%qsub(nr,nt) = 0._r8
               endif

               !---- negative qsur water, remove from TRunoff ---
               if (TRunoff%qsur(nr,nt) < 0._r8) then
                  src_direct(nt,cnt) = src_direct(nt,cnt) + TRunoff%qsur(nr,nt)
                  TRunoff%qsur(nr,nt) = 0._r8
               endif

               !---- water outside the basin ---
               !---- *** DO NOT TURN THIS ONE OFF, conservation will fail *** ---

               ! Note Tunit%mask is set in Tunit%init and is obtained from reading in fdir
               ! if fdir<0 then mask=0 (ocean), if fdir=0 then mask=2 (outlet) and if fdir>0 then mask=1 (land)
               if (Tunit%mask(nr) > 0) then
                  ! mosart euler
               else
                  ! NOTE: that when nt = nt_ice, the TRunoff terms
                  ! below have already been set to zero in the frozen
                  ! runoff calculation above - where frozen runoff is always set to the outlet
                  src_direct(nt,cnt) = src_direct(nt,cnt) + TRunoff%qsub(nr,nt) + TRunoff%qsur(nr,nt) + TRunoff%qgwl(nr,nt)
                  TRunoff%qsub(nr,nt) = 0._r8
                  TRunoff%qsur(nr,nt) = 0._r8
                  TRunoff%qgwl(nr,nt) = 0._r8
               end if
            enddo
         enddo

         call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         !--- copy direct transfer water to output field ---
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            do nt = 1,ntracers
               ctl%direct(nr,nt) = ctl%direct(nr,nt) + dst_direct(nt,cnt)
            enddo
         enddo
      endif
      call t_stopf('mosartr_SMdirect')

      !-----------------------------------
      ! mosart Subcycling
      !-----------------------------------

      call t_startf('mosartr_subcycling')

      if (first_call .and. mainproc) then
         do nt = 1,ntracers
            write(iulog,'(2a,i6,l4)') trim(subname),' euler_calc for nt = ',nt,Tunit%euler_calc(nt)
         enddo
      endif

      nsub = coupling_period/delt_mosart
      if (nsub*delt_mosart < coupling_period) then
         nsub = nsub + 1
      end if
      delt = delt_coupling/float(nsub)
      if (delt /= delt_save) then
         if (mainproc) then
            write(iulog,'(2a,2g20.12,2i12)') trim(subname),' mosart delt update from/to',&
                 delt_save,delt,nsub_save,nsub
         end if
      endif

      nsub_save = nsub
      delt_save = delt
      Tctl%DeltaT = delt

      !-----------------------------------
      ! mosart euler solver
      !-----------------------------------

      ! convert TRunoff fields from m3/s to m/s before calling Euler
      do nt = 1,ntracers
         do nr = begr,endr
            TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) / ctl%area(nr)
            TRunoff%qsub(nr,nt) = TRunoff%qsub(nr,nt) / ctl%area(nr)
            TRunoff%qgwl(nr,nt) = TRunoff%qgwl(nr,nt) / ctl%area(nr)
         enddo
      enddo

      ! Subcycle the call to Euler
      call t_startf('mosartr_euler')
      ctl%flow = 0._r8
      ctl%erout_prev = 0._r8
      ctl%eroutup_avg = 0._r8
      ctl%erlat_avg = 0._r8
      do ns = 1,nsub
         ! solve the ODEs with Euler algorithm
         call Euler(rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         ! accumulate local flow field
         do nt = 1,ntracers
            do nr = begr,endr
               ctl%flow(nr,nt)        = ctl%flow(nr,nt)        + TRunoff%flow(nr,nt)
               ctl%erout_prev(nr,nt)  = ctl%erout_prev(nr,nt)  + TRunoff%erout_prev(nr,nt)
               ctl%eroutup_avg(nr,nt) = ctl%eroutup_avg(nr,nt) + TRunoff%eroutup_avg(nr,nt)
               ctl%erlat_avg(nr,nt)   = ctl%erlat_avg(nr,nt)   + TRunoff%erlat_avg(nr,nt)
            enddo
         enddo
      enddo ! nsub
      call t_stopf('mosartr_euler')

      ! average flow over subcycling
      ctl%flow        = ctl%flow        / float(nsub)
      ctl%erout_prev  = ctl%erout_prev  / float(nsub)
      ctl%eroutup_avg = ctl%eroutup_avg / float(nsub)
      ctl%erlat_avg   = ctl%erlat_avg   / float(nsub)

      ! update states when subsycling completed
      ctl%runoff = 0._r8
      ctl%runofflnd = spval
      ctl%runoffocn = spval
      ctl%dvolrdt = 0._r8
      ctl%dvolrdtlnd = spval
      ctl%dvolrdtocn = spval
      do nt = 1,ntracers
         do nr = begr,endr
            volr_init = ctl%volr(nr,nt)
            ctl%volr(nr,nt) = (TRunoff%wt(nr,nt) + TRunoff%wr(nr,nt) + TRunoff%wh(nr,nt)*ctl%area(nr))
            ctl%dvolrdt(nr,nt) = (ctl%volr(nr,nt) - volr_init) / delt_coupling
            ctl%runoff(nr,nt) = ctl%flow(nr,nt)
            ctl%runofftot(nr,nt) = ctl%direct(nr,nt)
            if (ctl%mask(nr) == 1) then
               ctl%runofflnd(nr,nt) = ctl%runoff(nr,nt)
               ctl%dvolrdtlnd(nr,nt)= ctl%dvolrdt(nr,nt)
            elseif (ctl%mask(nr) >= 2) then
               ctl%runoffocn(nr,nt) = ctl%runoff(nr,nt)
               ctl%runofftot(nr,nt) = ctl%runofftot(nr,nt) + ctl%runoff(nr,nt)
               ctl%dvolrdtocn(nr,nt)= ctl%dvolrdt(nr,nt)
            endif
         enddo
      enddo

      ! final update from glc input
      do nr = begr,endr
        ctl%runofftot(nr,nt_liq) = ctl%runofftot(nr,nt_liq) + ctl%direct_glc(nr,nt_liq)
        ctl%runofftot(nr,nt_ice) = ctl%runofftot(nr,nt_ice) + ctl%direct_glc(nr,nt_ice)
      end do

      call t_stopf('mosartr_subcycling')

      !-----------------------------------
      ! BUDGET
      !-----------------------------------
      if (budget_check) then
        call t_startf('mosartr_budgetcheck')
        call budget%check_budget(begr,endr,ntracers,delt_coupling)
        call t_stopf('mosartr_budgetcheck')
      endif

      !-----------------------------------
      ! Write out mosart history file
      !-----------------------------------

      call t_startf('mosartr_hbuf')
      call mosart_histflds_set(ntracers)
      call mosart_hist_updatehbuf()
      call t_stopf('mosartr_hbuf')

      call t_startf('mosartr_htapes')
      call mosart_hist_htapeswrapup( rstwr, nlend )
      call t_stopf('mosartr_htapes')

      !-----------------------------------
      ! Write out mosart restart file
      !-----------------------------------

      if (rstwr) then
         call t_startf('mosartr_rest')
         filer = mosart_rest_filename(rdate=rdate)
         call mosart_rest_filewrite( filer, rdate=rdate )
         call t_stopf('mosartr_rest')
      end if

      !-----------------------------------
      ! Done
      !-----------------------------------

      first_call = .false.

      call t_stopf('mosartr_tot')

   end subroutine mosart_run

end module mosart_driver
