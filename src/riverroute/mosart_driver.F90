module mosart_driver

   !-----------------------------------------------------------------------
   ! Mosart Routing Model
   !-----------------------------------------------------------------------

   use shr_kind_mod       , only : r8 => shr_kind_r8, CS => shr_kind_cs, CL => shr_kind_CL
   use shr_sys_mod        , only : shr_sys_abort
   use shr_mpi_mod        , only : shr_mpi_sum, shr_mpi_max
   use shr_const_mod      , only : SHR_CONST_PI, SHR_CONST_CDAY
   use shr_string_mod     , only : shr_string_listGetNum, shr_string_listGetName
   use mosart_vars        , only : re, spval, iulog, ice_runoff, &
                                   frivinp, nsrContinue, nsrBranch, nsrStartup, nsrest, &
                                   inst_index, inst_suffix, inst_name, decomp_option, &
                                   bypass_routing_option, qgwl_runoff_option, barrier_timers, &
                                   mainproc, npes, iam, mpicom_rof
   use mosart_data        , only : ctl, Tctl, Tunit, TRunoff, Tpara
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
                                   ESMF_TERMORDER_SRCSEQ, ESMF_Mesh
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
   integer :: coupling_period ! mosart coupling period
   integer :: delt_mosart     ! mosart internal timestep (->nsub)
   logical :: use_halo_option ! enable halo capability using ESMF
   character(len=CS) :: mosart_tracers    ! colon delimited string of tracer names
   character(len=CS) :: mosart_euler_calc ! colon delimited string of logicals for using Euler  algorithm

   ! subcycling
   integer   :: nsub_save ! previous nsub
   real(r8)  :: delt_save ! previous delt

   ! global (glo)
   integer , allocatable :: IDkey(:) ! translation key from ID to gindex

   ! budget accumulation
   real(r8), allocatable :: budget_accum(:)  ! BUDGET accumulator over run
   integer               :: budget_accum_cnt ! counter for budget_accum

   character(len=CL) :: nlfilename_rof = 'mosart_in'
   character(len=CL) :: fnamer              ! name of netcdf restart file

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
           use_halo_option, delt_mosart, mosart_tracers, mosart_euler_calc

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

      ! Determine number of tracers and array of tracer names
      ctl%ntracers = shr_string_listGetNum(mosart_tracers)
      allocate(ctl%tracer_names(ctl%ntracers))
      do i = 1,ctl%ntracers
         call shr_string_listGetName(mosart_tracers, i, ctl%tracer_names(i))
      end do

      runtyp(:)               = 'missing'
      runtyp(nsrStartup  + 1) = 'initial'
      runtyp(nsrContinue + 1) = 'restart'
      runtyp(nsrBranch   + 1) = 'branch '

      if (mainproc) then
         write(iulog,*) 'define run:'
         write(iulog,'(a)'   ) '   run type              = '//trim(runtyp(nsrest+1))
         write(iulog,'(a,i8)') '   coupling_period       = ',coupling_period
         write(iulog,'(a,i8)') '   delt_mosart           = ',delt_mosart
         write(iulog,'(a)'   ) '   decomp option         = '//trim(decomp_option)
         write(iulog,'(a,l)' ) '   use_halo_option       = ',use_halo_option
         write(iulog,'(a)'   ) '   bypass_routing option = '//trim(bypass_routing_option)
         write(iulog,'(a)'   ) '   qgwl runoff option    = '//trim(qgwl_runoff_option)
         write(iulog,'(a)'   ) '   mosart tracers        = '//trim(mosart_tracers)
         write(iulog,'(a)'   ) '   mosart euler calc     = '//trim(mosart_euler_calc)
         if (nsrest == nsrStartup .and. finidat /= ' ') then
            write(iulog,'(a)') '   mosart initial data   = '//trim(finidat)
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

   subroutine mosart_init1(rc)

      !-------------------------------------------------
      ! Initialize mosart grid, mask, decomp
      !
      ! Arguments
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

      !-------------------------------------------------------
      ! Obtain restart file if appropriate
      !-------------------------------------------------------
      if ((nsrest == nsrStartup .and. finidat /= ' ') .or. &
          (nsrest == nsrContinue) .or. (nsrest == nsrBranch  )) then
         call mosart_rest_getfile( file=fnamer )
      endif

      !-------------------------------------------------------
      ! Initialize time manager
      !-------------------------------------------------------
      if (nsrest == nsrStartup) then
         call timemgr_init(dtime_in=coupling_period)
      else
         call mosart_rest_timemanager(file=fnamer)
      end if

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
      ! BUDGET terms 1-10 are for volumes (m3)
      ! BUDGET terms 11-30 are for flows (m3/s)
      ! even (2n) budget terms refer to current state odd terms (2n-1) rever to previous state.
      integer            :: i, j, n, nr, ns, nt, n2, nf ! indices
      real(r8)           :: budget_terms(30,ntracers)   ! BUDGET terms
      real(r8)           :: budget_input
      real(r8)           :: budget_output
      real(r8)           :: budget_volume
      real(r8)           :: budget_total
      real(r8)           :: budget_euler
      real(r8)           :: budget_eroutlag
      real(r8)           :: budget_global(30,ntracers)  ! global budget sum
      logical            :: budget_check                ! do global budget check
      real(r8),parameter :: budget_tolerance = 1.0e-6   ! budget tolerance, m3/day
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
         budget_accum = 0._r8
         budget_accum_cnt = 0
         delt_save = delt_mosart
         allocate(budget_accum(ntracers))
         if (mainproc) then
            write(iulog,'(2a,g20.12)') trim(subname),' mosart coupling period ',delt_coupling
         end if
      end if

      budget_check = .false.
      !TODO make budget check frequency adjustable
      if (day == 1 .and. mon == 1) budget_check = .true.
      if (tod == 0) budget_check = .true.
      budget_terms = 0._r8

      ! BUDGET
      ! BUDGET terms 1-10 are for volumes (m3)
      ! BUDGET terms 11-30 are for flows (m3/s)
      call t_startf('mosartr_budget')
      do nt = 1,ntracers
         do nr = begr,endr
            budget_terms( 1,nt) = budget_terms( 1,nt) + ctl%volr(nr,nt)
            budget_terms( 3,nt) = budget_terms( 3,nt) + TRunoff%wt(nr,nt)
            budget_terms( 5,nt) = budget_terms( 5,nt) + TRunoff%wr(nr,nt)
            budget_terms( 7,nt) = budget_terms( 7,nt) + TRunoff%wh(nr,nt)*ctl%area(nr)
            budget_terms(13,nt) = budget_terms(13,nt) + ctl%qsur(nr,nt)
            budget_terms(14,nt) = budget_terms(14,nt) + ctl%qsub(nr,nt)
            budget_terms(15,nt) = budget_terms(15,nt) + ctl%qgwl(nr,nt)
            budget_terms(17,nt) = budget_terms(17,nt) + ctl%qsur(nr,nt) + ctl%qsub(nr,nt)+ ctl%qgwl(nr,nt)
            if (nt==1) then
               budget_terms(16,nt) = budget_terms(16,nt) + ctl%qirrig(nr)
               budget_terms(17,nt) = budget_terms(17,nt) + ctl%qirrig(nr)
            endif
         enddo
      enddo
      call t_stopf('mosartr_budget')

      ! data for euler solver, in m3/s here
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
      nt = 1
      ctl%qirrig_actual = 0._r8
      do nr = begr,endr

         ! calculate volume of irrigation flux during timestep
         irrig_volume = -ctl%qirrig(nr) * coupling_period

         ! compare irrig_volume to main channel storage;
         ! add overage to subsurface runoff
         if(irrig_volume > TRunoff%wr(nr,nt)) then
            ctl%qsub(nr,nt) = ctl%qsub(nr,nt) + (TRunoff%wr(nr,nt) - irrig_volume) / coupling_period
            TRunoff%qsub(nr,nt) = ctl%qsub(nr,nt)
            irrig_volume = TRunoff%wr(nr,nt)
         endif

         ! actual irrigation rate [m3/s]
         ! i.e. the rate actually removed from the main channel
         ! if irrig_volume is greater than TRunoff%wr
         ctl%qirrig_actual(nr) = - irrig_volume / coupling_period

         ! remove irrigation from wr (main channel)
         TRunoff%wr(nr,nt) = TRunoff%wr(nr,nt) - irrig_volume

      enddo
      call t_stopf('mosartr_irrig')

      !-----------------------------------
      ! Compute flood
      ! Remove water from mosart and send back to clm
      ! Just consider land points and only remove liquid water
      ! ctl%flood is m3/s here
      !-----------------------------------

      call t_startf('mosartr_flood')
      nt = 1
      ctl%flood = 0._r8
      do nr = begr,endr
         ! initialize ctl%flood to zero
         if (ctl%mask(nr) == 1) then
            if (ctl%volr(nr,nt) > ctl%fthresh(nr)) then
               ! determine flux that is sent back to the land this is in m3/s
               ctl%flood(nr) = (ctl%volr(nr,nt)-ctl%fthresh(nr)) / (delt_coupling)

               ! ctl%flood will be sent back to land - so must subtract this
               ! from the input runoff from land
               ! tcraig, comment - this seems like an odd approach, you
               !   might create negative forcing.  why not take it out of
               !   the volr directly?  it's also odd to compute this
               !   at the initial time of the time loop.  why not do
               !   it at the end or even during the run loop as the
               !   new volume is computed.  fluxout depends on volr, so
               !   how this is implemented does impact the solution.
               TRunoff%qsur(nr,nt) = TRunoff%qsur(nr,nt) - ctl%flood(nr)
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
      !--- all frozen runoff passed direct to outlet
      !-----------------------------------------------------

      nt = 2
      src_direct(:,:) = 0._r8
      dst_direct(:,:) = 0._r8

      ! set euler_calc = false for frozen runoff
      ! TODO: will be reworked after addition of multiple tracers 
      Tunit%euler_calc(nt) = .false.

      cnt = 0
      do nr = begr,endr
         cnt = cnt + 1
         src_direct(nt,cnt) = TRunoff%qsur(nr,nt) + TRunoff%qsub(nr,nt) + TRunoff%qgwl(nr,nt)
         TRunoff%qsur(nr,nt) = 0._r8
         TRunoff%qsub(nr,nt) = 0._r8
         TRunoff%qgwl(nr,nt) = 0._r8
      enddo

      call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! copy direct transfer water to output field
      ctl%direct = 0._r8
      cnt = 0
      do nr = begr,endr
         cnt = cnt + 1
         ctl%direct(nr,nt) = ctl%direct(nr,nt) + dst_direct(nt,cnt)
      enddo

      !-----------------------------------------------------
      !--- direct to outlet qgwl
      !-----------------------------------------------------

      !-- liquid runoff components
      if (trim(bypass_routing_option) == 'direct_to_outlet') then

         nt = 1
         src_direct(:,:) = 0._r8
         dst_direct(:,:) = 0._r8

         !--- copy direct transfer fields, convert kg/m2s to m3/s
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            if (trim(qgwl_runoff_option) == 'all') then
               src_direct(nt,cnt) = TRunoff%qgwl(nr,nt)
               TRunoff%qgwl(nr,nt) = 0._r8
            else if (trim(qgwl_runoff_option) == 'negative') then
               if(TRunoff%qgwl(nr,nt) < 0._r8) then
                  src_direct(nt,cnt) = TRunoff%qgwl(nr,nt)
                  TRunoff%qgwl(nr,nt) = 0._r8
               endif
            endif
         enddo

         call ESMF_FieldSMM(Tunit%srcfield, Tunit%dstfield, Tunit%rh_direct, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         !--- copy direct transfer water to output field ---
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            ctl%direct(nr,nt) = ctl%direct(nr,nt) + dst_direct(nt,cnt)
         enddo
      endif

      !-----------------------------------------------------
      !--- direct in place qgwl
      !-----------------------------------------------------

      if (trim(bypass_routing_option) == 'direct_in_place') then

         nt = 1
         do nr = begr,endr

            if (trim(qgwl_runoff_option) == 'all') then
               ctl%direct(nr,nt) = TRunoff%qgwl(nr,nt)
               TRunoff%qgwl(nr,nt) = 0._r8
            else if (trim(qgwl_runoff_option) == 'negative') then
               if(TRunoff%qgwl(nr,nt) < 0._r8) then
                  ctl%direct(nr,nt) = TRunoff%qgwl(nr,nt)
                  TRunoff%qgwl(nr,nt) = 0._r8
               endif
            else if (trim(qgwl_runoff_option) == 'threshold') then
               ! --- calculate volume of qgwl flux during timestep
               qgwl_volume = TRunoff%qgwl(nr,nt) * ctl%area(nr) * coupling_period
               river_volume_minimum = river_depth_minimum * ctl%area(nr)

               ! if qgwl is negative, and adding it to the main channel
               ! would bring main channel storage below a threshold,
               ! send qgwl directly to ocean
               if (((qgwl_volume + TRunoff%wr(nr,nt)) < river_volume_minimum) .and. (TRunoff%qgwl(nr,nt) < 0._r8)) then
                  ctl%direct(nr,nt) = TRunoff%qgwl(nr,nt)
                  TRunoff%qgwl(nr,nt) = 0._r8
               endif
            endif
         enddo

      endif

      !-------------------------------------------------------
      !--- add other direct terms, e.g. inputs outside of
      !--- mosart mask, negative qsur
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
               if (Tunit%mask(nr) > 0) then
                  ! mosart euler
               else
                  src_direct(nt,cnt) = src_direct(nt,cnt) + TRunoff%qsub(nr,nt) + TRunoff%qsur(nr,nt) &
                       + TRunoff%qgwl(nr,nt)
                  TRunoff%qsub(nr,nt) = 0._r8
                  TRunoff%qsur(nr,nt) = 0._r8
                  TRunoff%qgwl(nr,nt) = 0._r8
               endif
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

      call t_startf('mosartr_budget')
      do nt = 1,ntracers
         do nr = begr,endr
            budget_terms(20,nt) = budget_terms(20,nt) + TRunoff%qsur(nr,nt) + TRunoff%qsub(nr,nt) + TRunoff%qgwl(nr,nt)
            budget_terms(29,nt) = budget_terms(29,nt) + TRunoff%qgwl(nr,nt)
         enddo
      enddo
      call t_stopf('mosartr_budget')

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
      ! TODO: move of this to hist_set_flds
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
      call t_stopf('mosartr_subcycling')

      !-----------------------------------
      ! BUDGET
      !-----------------------------------

      ! BUDGET terms 1-10 are for volumes (m3)
      ! BUDGET terms 11-30 are for flows (m3/s)
      ! BUDGET only ocean runoff and direct gets out of the system

      call t_startf('mosartr_budget')
      do nt = 1,ntracers
         do nr = begr,endr
            budget_terms( 2,nt) = budget_terms( 2,nt) + ctl%volr(nr,nt)
            budget_terms( 4,nt) = budget_terms( 4,nt) + TRunoff%wt(nr,nt)
            budget_terms( 6,nt) = budget_terms( 6,nt) + TRunoff%wr(nr,nt)
            budget_terms( 8,nt) = budget_terms( 8,nt) + TRunoff%wh(nr,nt)*ctl%area(nr)
            budget_terms(21,nt) = budget_terms(21,nt) + ctl%direct(nr,nt)
            if (ctl%mask(nr) >= 2) then
               budget_terms(18,nt) = budget_terms(18,nt) + ctl%runoff(nr,nt)
               budget_terms(26,nt) = budget_terms(26,nt) - ctl%erout_prev(nr,nt)
               budget_terms(27,nt) = budget_terms(27,nt) + ctl%flow(nr,nt)
            else
               budget_terms(23,nt) = budget_terms(23,nt) - ctl%erout_prev(nr,nt)
               budget_terms(24,nt) = budget_terms(24,nt) + ctl%flow(nr,nt)
            endif
            budget_terms(25,nt) = budget_terms(25,nt) - ctl%eroutup_avg(nr,nt)
            budget_terms(28,nt) = budget_terms(28,nt) - ctl%erlat_avg(nr,nt)
            budget_terms(22,nt) = budget_terms(22,nt) + ctl%runoff(nr,nt) + ctl%direct(nr,nt) + ctl%eroutup_avg(nr,nt)
         enddo
      enddo
      nt = 1
      do nr = begr,endr
         budget_terms(19,nt) = budget_terms(19,nt) + ctl%flood(nr)
         budget_terms(22,nt) = budget_terms(22,nt) + ctl%flood(nr)
      enddo

      ! accumulate the budget total over the run to make sure it's decreasing on avg
      budget_accum_cnt = budget_accum_cnt + 1
      do nt = 1,ntracers
         budget_volume = (budget_terms( 2,nt) - budget_terms( 1,nt)) / delt_coupling
         budget_input  = (budget_terms(13,nt) + budget_terms(14,nt) + budget_terms(15,nt) + budget_terms(16,nt))
         budget_output = (budget_terms(18,nt) + budget_terms(19,nt) + budget_terms(21,nt))
         budget_total  = budget_volume - budget_input + budget_output
         budget_accum(nt) = budget_accum(nt) + budget_total
         budget_terms(30,nt) = budget_accum(nt)/budget_accum_cnt
      enddo
      call t_stopf('mosartr_budget')

      if (budget_check) then
         call t_startf('mosartr_budget')
         !--- check budget

         ! convert fluxes from m3/s to m3 by mult by coupling_period
         budget_terms(11:30,:) = budget_terms(11:30,:) * delt_coupling

         ! convert terms from m3 to million m3
         budget_terms(:,:) = budget_terms(:,:) * 1.0e-6_r8

         ! global sum
         call shr_mpi_sum(budget_terms,budget_global,mpicom_rof,'mosart global budget',all=.false.)

         ! write budget
         if (mainproc) then
            write(iulog,'(2a,i10,i6)') trim(subname),' mosart BUDGET diagnostics (million m3) for ',ymd,tod
            do nt = 1,ntracers
               budget_volume = (budget_global( 2,nt) - budget_global( 1,nt))
               budget_input  = (budget_global(13,nt) + budget_global(14,nt) + budget_global(15,nt))
               budget_output = (budget_global(18,nt) + budget_global(19,nt) + budget_global(21,nt))
               budget_total  = budget_volume - budget_input + budget_output
               budget_euler  = budget_volume - budget_global(20,nt) + budget_global(18,nt)
               budget_eroutlag = budget_global(23,nt) - budget_global(24,nt)
               write(iulog,'(2a,i4)')       trim(subname),'  tracer = ',nt
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   volume   init = ',nt,budget_global(1,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   volume  final = ',nt,budget_global(2,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   input surface = ',nt,budget_global(13,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   input subsurf = ',nt,budget_global(14,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   input gwl     = ',nt,budget_global(15,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   input irrig   = ',nt,budget_global(16,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   input total   = ',nt,budget_global(17,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   output flow   = ',nt,budget_global(18,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   output direct = ',nt,budget_global(21,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   output flood  = ',nt,budget_global(19,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   output total  = ',nt,budget_global(22,nt)
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   sum input     = ',nt,budget_input
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   sum dvolume   = ',nt,budget_volume
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   sum output    = ',nt,budget_output
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   net (dv-i+o)  = ',nt,budget_total
               write(iulog,'(2a,i4,f22.6)') trim(subname),'   eul erout lag = ',nt,budget_eroutlag
               if ((budget_total-budget_eroutlag) > 1.0e-6) then
                  write(iulog,'(2a,i4)') trim(subname),' ***** BUDGET WARNING error gt 1. m3 for nt = ',nt
               endif
               if ((budget_total+budget_eroutlag) >= 1.0e-6) then
                  if ((budget_total-budget_eroutlag)/(budget_total+budget_eroutlag) > 0.001_r8) then
                     write(iulog,'(2a,i4)') trim(subname),' ***** BUDGET WARNING out of balance for nt = ',nt
                  endif
               endif
            enddo
            write(iulog,'(a)') '----------------------------------- '
         endif

         call t_stopf('mosartr_budget')
      endif  ! budget_check

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
