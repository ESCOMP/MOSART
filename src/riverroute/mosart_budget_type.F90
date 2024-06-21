module mosart_budget_type

   ! Variables and routines used for
   ! calculating and checking tracer global and local budgets

   use shr_kind_mod,       only: r8 => shr_kind_r8, CL => SHR_KIND_CL
   use shr_sys_mod,        only: shr_sys_abort
   use shr_mpi_mod,        only: shr_mpi_sum, shr_mpi_max
   use mosart_vars,        only: re, spval, barrier_timers, iulog, mainproc, npes, iam, mpicom_rof
   use mosart_data,        only: ctl, Tctl, Tunit, TRunoff, Tpara
   use mosart_timemanager, only: get_nstep, get_curr_date

   implicit none
   private

   type budget_type
      ! accumulated budget over run (not used for now)
      real(r8), pointer :: accum_grc(:, :)            ! Gridcell level budget accumulator per tracer over the run (m3)
      real(r8), pointer :: accum_glob(:)              ! Global budget accumulator (1e6 m3)

      ! budget terms per gridcell
      real(r8), pointer :: beg_vol_grc(:, :)          ! volume begining of the timestep (m3)
      real(r8), pointer :: end_vol_grc(:, :)          ! volume end of the timestep (m3)
      real(r8), pointer :: in_grc(:, :)               ! budget in terms (m3)
      real(r8), pointer :: out_grc(:, :)              ! budget out terms (m3)
      real(r8), pointer :: net_grc(:, :)              ! net budget (dvolume + inputs - outputs) (m3)
      real(r8), pointer :: lag_grc(:, :)              ! euler erout lagged (m3)

      ! budget global terms
      real(r8), pointer :: beg_vol_glob(:)            ! volume begining of the timestep (1e6 m3)
      real(r8), pointer :: end_vol_glob(:)            ! volume end of the timestep (1e6 m3)
      real(r8), pointer :: in_glob(:)                 ! budget in terms (1e6 m3)
      real(r8), pointer :: out_glob(:)                ! budget out terms (1e6 m3)
      real(r8), pointer :: net_glob(:)                ! net budget (dvolume + inputs - outputs) (1e6 m3)
      real(r8), pointer :: lag_glob(:)                ! euler erout lagged (1e6 m3)

      ! budget parameters
      real(r8)             :: tolerance = 1e-6_r8     ! budget absolute tolerance
      real(r8)             :: rel_tolerance = 1e-6_r8 ! budget relative tolerance
      logical(r8), pointer :: do_budget(:)            ! if budget should be checked (per tracer)
   contains
      procedure, public :: Init
      procedure, public :: set_budget
      procedure, public :: check_budget
   end type budget_type
   public :: budget_type

   integer, parameter :: index_beg_vol_grc = 1
   integer, parameter :: index_end_vol_grc = 2
   integer, parameter :: index_in_grc      = 3
   integer, parameter :: index_out_grc     = 4
   integer, parameter :: index_net_grc     = 5
   integer, parameter :: index_lag_grc     = 6

   character(*), parameter :: u_FILE_u = &
                              __FILE__

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine Init(this, begr, endr, ntracers)

      ! Initialize budget type

      ! Arguments
      class(budget_type) :: this
      integer, intent(in) :: begr, endr, ntracers
      !-------------------------------------------------

      ! gridcell level:
      allocate (this%accum_grc(begr:endr, ntracers))
      this%accum_grc = 0._r8

      allocate (this%beg_vol_grc(begr:endr, ntracers))
      this%beg_vol_grc = 0._r8

      allocate (this%end_vol_grc(begr:endr, ntracers))
      this%end_vol_grc = 0._r8

      allocate (this%in_grc(begr:endr, ntracers))
      this%in_grc = 0._r8

      allocate (this%out_grc(begr:endr, ntracers))
      this%out_grc = 0._r8

      allocate (this%net_grc(begr:endr, ntracers))
      this%net_grc = 0._r8

      allocate (this%lag_grc(begr:endr, ntracers))
      this%lag_grc = 0._r8

      ! global level:
      allocate (this%accum_glob(ntracers))
      this%accum_glob = 0._r8

      allocate (this%beg_vol_glob(ntracers))
      this%beg_vol_glob = 0._r8

      allocate (this%end_vol_glob(ntracers))
      this%end_vol_glob = 0._r8

      allocate (this%in_glob(ntracers))
      this%in_glob = 0._r8

      allocate (this%out_glob(ntracers))
      this%out_glob = 0._r8

      allocate (this%net_glob(ntracers))
      this%net_glob = 0._r8

      allocate (this%lag_glob(ntracers))
      this%lag_glob = 0._r8

      allocate (this%do_budget(ntracers))
      this%do_budget = .true.

   end subroutine Init

   !-----------------------------------------------------------------------

   subroutine set_budget(this, begr, endr, ntracers, dt)

      ! Arguments
      class(budget_type)   :: this
      integer, intent(in)  :: begr, endr, ntracers
      real(r8), intent(in) :: dt

      ! local variables
      integer :: nr, nt   !indices
      integer :: nt_liq, nt_ice
      !-------------------------------------------------

      nt_liq = ctl%nt_liq
      nt_ice = ctl%nt_ice
      do nr = begr, endr
         do nt = 1, ntracers
            this%beg_vol_grc(nr, nt) = ctl%volr(nr, nt)
            if (nt == nt_ice) then
               this%in_grc(nr, nt) = (ctl%qsur(nr, nt) + ctl%qsub(nr, nt) + ctl%qgwl(nr, nt) + ctl%qglc_ice(nr)) * dt
            else if (nt == nt_liq) then
               this%in_grc(nr, nt) = (ctl%qsur(nr, nt) + ctl%qsub(nr, nt) + ctl%qgwl(nr, nt) + ctl%qglc_liq(nr)) * dt
            end if
            ! this was for budget_terms(17)
            !if (nt==1) then
            !  this%in_grc(nr,nt)=this%in_grc(nr,nt) +ctl%qirrig(nr)
            !endif
         end do
      end do

      this%end_vol_grc(:,:) = 0.0_r8
      this%out_grc(:,:) = 0.0_r8
      this%net_grc(:,:) = 0.0_r8
      this%lag_grc(:,:) = 0.0_r8

      this%beg_vol_glob(:) = 0.0_r8
      this%end_vol_glob(:) = 0.0_r8
      this%in_glob(:) = 0.0_r8
      this%out_glob(:) = 0.0_r8
      this%net_glob(:) = 0.0_r8
      this%lag_glob(:) = 0.0_r8

   end subroutine set_budget

   !-----------------------------------------------------------------------

   subroutine check_budget(this, begr, endr, ntracers, dt)

      ! Arguments
      class(budget_type)   :: this
      integer, intent(in)  :: begr, endr, ntracers
      real(r8), intent(in) :: dt

      ! Local variables
      integer  :: nr, nt                !indecies
      integer  :: nt_liq                ! liquid index
      integer  :: yr,mon,day,ymd,tod    !time vars
      real(r8) :: tmp_in(6, ntracers)   ! array to pass to mpi_sum
      real(r8) :: tmp_glob(6, ntracers) ! array from mpi_sum
      logical  :: error_budget          ! flag for an error
      real(r8) :: abserr, relerr
      !-------------------------------------------------

      call get_curr_date(yr, mon, day, tod)
      ymd = yr*10000 + mon*100 + day
      tmp_in = 0.0_r8
      tmp_glob = 0.0_r8

      nt_liq = ctl%nt_liq
      do nr = begr, endr
         do nt = 1, ntracers
            this%end_vol_grc(nr, nt) = ctl%volr(nr, nt)
            this%out_grc(nr, nt) = this%out_grc(nr, nt) + ctl%direct(nr, nt) + ctl%direct_glc(nr, nt)
            if (nt == nt_liq) then
               this%out_grc(nr, nt) = this%out_grc(nr, nt) + ctl%flood(nr)
            end if
            if (ctl%mask(nr) >= 2) then
               this%out_grc(nr, nt) = this%out_grc(nr, nt) + ctl%runoff(nr, nt)
            else
               this%lag_grc(nr, nt) = this%lag_grc(nr, nt) - ctl%erout_prev(nr, nt) - ctl%flow(nr, nt)
            end if
            this%out_grc(nr,nt) = this%out_grc(nr,nt) * dt
            this%lag_grc(nr,nt) = this%lag_grc(nr,nt) * dt
            this%net_grc(nr,nt) = this%end_vol_grc(nr,nt) - this%beg_vol_grc(nr,nt) - (this%in_grc(nr,nt)-this%out_grc(nr,nt))
            this%accum_grc(nr,nt) = this%accum_grc(nr,nt) + this%net_grc(nr,nt)
         end do
      end do

      do nt = 1, ntracers
         tmp_in(index_beg_vol_grc, nt) = sum(this%beg_vol_grc(:, nt))
         tmp_in(index_end_vol_grc, nt) = sum(this%end_vol_grc(:, nt))
         tmp_in(index_in_grc, nt)      = sum(this%in_grc(:, nt))
         tmp_in(index_out_grc, nt)     = sum(this%out_grc(:, nt))
         tmp_in(index_net_grc, nt)     = sum(this%net_grc(:, nt))
         tmp_in(index_lag_grc, nt)     = sum(this%lag_grc(:, nt))
      end do

      tmp_in = tmp_in*1e-6_r8 !convert to million m3
      call shr_mpi_sum(tmp_in, tmp_glob, mpicom_rof, 'mosart global budget', all=.false.)

      do nt = 1, ntracers
         error_budget = .false.
         abserr = 0.0_r8
         relerr = 0.0_r8
         this%beg_vol_glob(nt) = tmp_glob(index_beg_vol_grc, nt)
         this%end_vol_glob(nt) = tmp_glob(index_end_vol_grc, nt)
         this%in_glob(nt)      = tmp_glob(index_in_grc, nt)
         this%out_glob(nt)     = tmp_glob(index_out_grc, nt)
         this%net_glob(nt)     = tmp_glob(index_net_grc, nt)
         this%lag_glob(nt)     = tmp_glob(index_lag_grc, nt)
         if (this%do_budget(nt)) then
            if (abs(this%net_glob(nt) - this%lag_glob(nt)*dt) > this%tolerance) then
               error_budget = .true.
               abserr = abs(this%net_glob(nt) - this%lag_glob(nt))
            end if
            if (abs(this%net_glob(nt) + this%lag_glob(nt)) > 1e-6) then
               if (  abs(this%net_glob(nt) - this%lag_glob(nt)) &
                    /abs(this%net_glob(nt) + this%lag_glob(nt)) > this%rel_tolerance) then
                  error_budget = .true.
                  relerr = abs(this%net_glob(nt) - this%lag_glob(nt)) /abs(this%net_glob(nt) + this%lag_glob(nt))
               end if
            end if
            if (mainproc) then
              write (iulog, '(a)')         '-----------------------------------'
              write (iulog, '(a)')         '*****MOSART BUDGET DIAGNOSTICS*****'
              write (iulog,'(a,i10,i6)')   '   diagnostics for ', ymd, tod
              write (iulog, '(a,i4,2a)')   '   tracer                      = ', nt, '  ', ctl%tracer_names(nt)
              write (iulog, '(a,f22.6,a)') '   time step size              = ', dt, ' sec'
              write (iulog, '(a,f22.6,a)') '   volume begining of the step = ', this%beg_vol_glob(nt), ' (mil m3)'
              write (iulog, '(a,f22.6,a)') '   volume end of the step      = ', this%end_vol_glob(nt), ' (mil m3)'
              write (iulog, '(a,f22.6,a)') '   inputs                      = ', this%in_glob(nt), ' (mil m3)'
              write (iulog, '(a,f22.6,a)') '   outputs                     = ', this%out_glob(nt), ' (mil m3)'
              write (iulog, '(a,f22.6,a)') '   net budget (dv -i + o)      = ', this%net_glob(nt), ' (mil m3)'
              write (iulog, '(a,f22.6,a)') '   eul erout lag               = ', this%lag_glob(nt), '(mil m3)'
              write (iulog, '(a,f22.6)')   '   absolute budget error       = ', abserr
              write (iulog, '(a,f22.6)')   '   relative budget error       = ', relerr
              if (error_budget) then
                  write(iulog,'(a)')       '   BUDGET OUT OF BALANCE WARNING   '
              endif
              write (iulog, '(a)')         '-----------------------------------'
            end if
         end if
      end do

   end subroutine check_budget

end module mosart_budget_type
