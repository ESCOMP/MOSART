MODULE MOSART_physics_mod

   !-----------------------------------------------------------------------
   ! Description: core code of MOSART.
   ! Contains routines for solving diffusion wave and update the state of 
   ! hillslope, subnetwork and main channel variables
   ! Developed by Hongyi Li, 12/29/2011.
   !-----------------------------------------------------------------------

   use shr_kind_mod      , only : r8 => shr_kind_r8
   use shr_const_mod     , only : SHR_CONST_REARTH, SHR_CONST_PI
   use shr_sys_mod       , only : shr_sys_abort
   use mosart_vars       , only : iulog, barrier_timers, mpicom_rof
   use mosart_data       , only : Tctl, TUnit, TRunoff, TPara, ctl
   use perf_mod          , only : t_startf, t_stopf
   use nuopc_shr_methods , only : chkerr
   use ESMF              , only : ESMF_FieldGet, ESMF_FieldSMM, ESMF_Finalize, &
                                  ESMF_SUCCESS, ESMF_END_ABORT, ESMF_TERMORDER_SRCSEQ

   implicit none
   private

   public :: Euler
   public :: updatestate_hillslope
   public :: updatestate_subnetwork
   public :: updatestate_mainchannel
   public :: hillsloperouting
   public :: subnetworkrouting
   public :: mainchannelrouting

   real(r8), parameter :: TINYVALUE = 1.0e-14_r8 ! double precision variable has a significance of about 16 decimal digits
   real(r8), parameter :: SLOPE1def = 0.1_r8     ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.

   character(*), parameter :: u_FILE_u = &
        __FILE__

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine Euler(rc)

      ! solve the ODEs with Euler algorithm

      ! Arguments
      integer, intent(out) :: rc

      ! Local variables
      integer           :: nt, nr, m, k, unitUp, cnt, ier   !local index
      real(r8)          :: temp_erout, localDeltaT
      real(r8)          :: negchan
      real(r8), pointer :: src_eroutUp(:,:)
      real(r8), pointer :: dst_eroutUp(:,:)

      !------------------
      ! hillslope
      !------------------

      rc = ESMF_SUCCESS

      call t_startf('mosartr_hillslope')
      do nt=1,ctl%ntracers
         if (TUnit%euler_calc(nt)) then
            do nr=ctl%begr,ctl%endr
               if(TUnit%mask(nr) > 0) then
                  call hillslopeRouting(nr,nt,Tctl%DeltaT)
                  TRunoff%wh(nr,nt) = TRunoff%wh(nr,nt) + TRunoff%dwh(nr,nt) * Tctl%DeltaT

                  call UpdateState_hillslope(nr,nt)
                  TRunoff%etin(nr,nt) = (-TRunoff%ehout(nr,nt) + TRunoff%qsub(nr,nt)) * TUnit%area(nr) * TUnit%frac(nr)
               endif
            end do
         endif
      end do
      call t_stopf('mosartr_hillslope')

      call ESMF_FieldGet(Tunit%srcfield, farrayPtr=src_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(Tunit%dstfield, farrayPtr=dst_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      src_eroutUp(:,:) = 0._r8
      dst_eroutUp(:,:) = 0._r8

      TRunoff%flow = 0._r8
      TRunoff%erout_prev = 0._r8
      TRunoff%eroutup_avg = 0._r8
      TRunoff%erlat_avg = 0._r8
      negchan = 9999.0_r8

      do m=1,Tctl%DLevelH2R

         !--- accumulate/average erout at prior timestep (used in eroutUp calc) for budget analysis
         do nt=1,ctl%ntracers
            if (TUnit%euler_calc(nt)) then
               do nr=ctl%begr,ctl%endr
                  TRunoff%erout_prev(nr,nt) = TRunoff%erout_prev(nr,nt) + TRunoff%erout(nr,nt)
               end do
            end if
         end do

         !------------------
         ! subnetwork
         !------------------

         call t_startf('mosartr_subnetwork')
         TRunoff%erlateral(:,:) = 0._r8
         do nt=1,ctl%ntracers
            if (TUnit%euler_calc(nt)) then
               do nr=ctl%begr,ctl%endr
                  if(TUnit%mask(nr) > 0) then
                     localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(nr)
                     do k=1,TUnit%numDT_t(nr)
                        call subnetworkRouting(nr,nt,localDeltaT)
                        TRunoff%wt(nr,nt) = TRunoff%wt(nr,nt) + TRunoff%dwt(nr,nt) * localDeltaT
                        call UpdateState_subnetwork(nr,nt)
                        TRunoff%erlateral(nr,nt) = TRunoff%erlateral(nr,nt)-TRunoff%etout(nr,nt)
                     end do ! numDT_t
                     TRunoff%erlateral(nr,nt) = TRunoff%erlateral(nr,nt) / TUnit%numDT_t(nr)
                  endif
               end do ! nr
            endif  ! euler_calc
         end do ! nt
         call t_stopf('mosartr_subnetwork')

         !------------------
         ! upstream interactions
         !------------------

         if (barrier_timers) then
            call t_startf('mosartr_SMeroutUp_barrier')
            call mpi_barrier(mpicom_rof,ier)
            call t_stopf('mosartr_SMeroutUp_barrier')
         endif

         call t_startf('mosartr_SMeroutUp')

         !--- copy erout into src_eroutUp ---
         TRunoff%eroutUp = 0._r8
         src_eroutUp(:,:) = 0._r8
         cnt = 0
         do nr = ctl%begr,ctl%endr
            cnt = cnt + 1
            do nt = 1,ctl%ntracers
               src_eroutUp(nt,cnt) = TRunoff%erout(nr,nt)
            enddo
         enddo

         ! --- map src_eroutUp to dst_eroutUp
         call ESMF_FieldSMM(TUnit%srcfield, TUnit%dstField, TUnit%rh_eroutUp, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         !--- copy mapped eroutUp to TRunoff ---
         cnt = 0
         do nr = ctl%begr,ctl%endr
            cnt = cnt + 1
            do nt = 1,ctl%ntracers
               TRunoff%eroutUp(nr,nt) = dst_eroutUp(nt,cnt)
            enddo
         enddo

         call t_stopf('mosartr_SMeroutUp')

         TRunoff%eroutup_avg = TRunoff%eroutup_avg + TRunoff%eroutUp
         TRunoff%erlat_avg   = TRunoff%erlat_avg   + TRunoff%erlateral

         !------------------
         ! channel routing
         !------------------

         call t_startf('mosartr_chanroute')
         do nt=1,ctl%ntracers
            if (TUnit%euler_calc(nt)) then
               do nr=ctl%begr,ctl%endr
                  if(TUnit%mask(nr) > 0) then
                     localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(nr)
                     temp_erout = 0._r8
                     do k=1,TUnit%numDT_r(nr)
                        call mainchannelRouting(nr,nt,localDeltaT)
                        TRunoff%wr(nr,nt) = TRunoff%wr(nr,nt) + TRunoff%dwr(nr,nt) * localDeltaT
                        ! check for negative channel storage
                        ! if(TRunoff%wr(nr,1) < -1.e-10) then
                        !    write(iulog,*) 'Negative channel storage! ', nr, TRunoff%wr(nr,1)
                        !    call shr_sys_abort('mosart: negative channel storage')
                        ! end if
                        call UpdateState_mainchannel(nr,nt)
                        ! erout here might be inflow to some downstream subbasin, so treat it differently than erlateral
                        temp_erout = temp_erout + TRunoff%erout(nr,nt)
                     end do
                     temp_erout = temp_erout / TUnit%numDT_r(nr)
                     TRunoff%erout(nr,nt) = temp_erout
                     TRunoff%flow(nr,nt) = TRunoff%flow(nr,nt) - TRunoff%erout(nr,nt)
                  endif
               end do ! nr
            endif  ! euler_calc
         end do ! nt
         negchan = min(negchan, minval(TRunoff%wr(:,:)))

         call t_stopf('mosartr_chanroute')
      end do

      ! check for negative channel storage
      if (negchan < -1.e-10) then
         write(iulog,*) 'Warning: Negative channel storage found! ',negchan
         ! call shr_sys_abort('mosart: negative channel storage')
      endif
      TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
      TRunoff%erout_prev = TRunoff%erout_prev / Tctl%DLevelH2R
      TRunoff%eroutup_avg = TRunoff%eroutup_avg / Tctl%DLevelH2R
      TRunoff%erlat_avg = TRunoff%erlat_avg / Tctl%DLevelH2R

   end subroutine Euler

   !-----------------------------------------------------------------------

   subroutine hillslopeRouting(nr, nt, theDeltaT)
      !  Hillslope routing considering uniform runoff generation across hillslope

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

      TRunoff%ehout(nr,nt) = -CREHT_nosqrt(TUnit%hslpsqrt(nr), TUnit%nh(nr), TUnit%Gxr(nr), TRunoff%yh(nr,nt))
      if(TRunoff%ehout(nr,nt) < 0._r8 .and. &
           TRunoff%wh(nr,nt) + (TRunoff%qsur(nr,nt) + TRunoff%ehout(nr,nt)) * theDeltaT < TINYVALUE) then
         TRunoff%ehout(nr,nt) = -(TRunoff%qsur(nr,nt) + TRunoff%wh(nr,nt) / theDeltaT)
      end if
      TRunoff%dwh(nr,nt) = (TRunoff%qsur(nr,nt) + TRunoff%ehout(nr,nt))

   end subroutine hillslopeRouting

   !-----------------------------------------------------------------------

   subroutine subnetworkRouting(nr,nt,theDeltaT)
      !  subnetwork channel routing

      ! Arguments
      integer, intent(in) :: nr,nt
      real(r8), intent(in) :: theDeltaT

      if(TUnit%tlen(nr) <= TUnit%hlen(nr)) then ! if no tributaries, not subnetwork channel routing
         TRunoff%etout(nr,nt) = -TRunoff%etin(nr,nt)
      else
         TRunoff%vt(nr,nt) = CRVRMAN_nosqrt(TUnit%tslpsqrt(nr), TUnit%nt(nr), TRunoff%rt(nr,nt))
         TRunoff%etout(nr,nt) = -TRunoff%vt(nr,nt) * TRunoff%mt(nr,nt)
         if(TRunoff%wt(nr,nt) + (TRunoff%etin(nr,nt) + TRunoff%etout(nr,nt)) * theDeltaT < TINYVALUE) then
            TRunoff%etout(nr,nt) = -(TRunoff%etin(nr,nt) + TRunoff%wt(nr,nt)/theDeltaT)
            if(TRunoff%mt(nr,nt) > 0._r8) then
               TRunoff%vt(nr,nt) = -TRunoff%etout(nr,nt)/TRunoff%mt(nr,nt)
            end if
         end if
      end if
      TRunoff%dwt(nr,nt) = TRunoff%etin(nr,nt) + TRunoff%etout(nr,nt)

      ! check stability
      !    if(TRunoff%vt(nr,nt) < -TINYVALUE .or. TRunoff%vt(nr,nt) > 30) then
      !       write(iulog,*) "Numerical error in subnetworkRouting, ", nr,nt,TRunoff%vt(nr,nt)
      !    end if

   end subroutine subnetworkRouting

   !-----------------------------------------------------------------------

   subroutine mainchannelRouting(nr, nt, theDeltaT)
      !  main channel routing

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

      if(Tctl%RoutingMethod == 1) then
         call Routing_KW(nr, nt, theDeltaT)
      else if(Tctl%RoutingMethod == 2) then
         call Routing_MC(nr, nt, theDeltaT)
      else if(Tctl%RoutingMethod == 3) then
         call Routing_THREW(nr, nt, theDeltaT)
      else if(Tctl%RoutingMethod == 4) then
         call Routing_DW(nr, nt, theDeltaT)
      else
         call shr_sys_abort( "mosart: Please check the routing method! There are only 4 methods available." )
      end if

   end subroutine mainchannelRouting

   !-----------------------------------------------------------------------

   subroutine Routing_KW(nr, nt, theDeltaT)
      !  classic kinematic wave routing method

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

      ! Local variables
      integer  :: k
      real(r8) :: temp_gwl, temp_dwr, temp_gwl0

      ! estimate the inflow from upstream units
      TRunoff%erin(nr,nt) = 0._r8
      TRunoff%erin(nr,nt) = TRunoff%erin(nr,nt) - TRunoff%eroutUp(nr,nt)

      ! estimate the outflow
      if(TUnit%rlen(nr) <= 0._r8) then ! no river network, no channel routing
         TRunoff%vr(nr,nt) = 0._r8
         TRunoff%erout(nr,nt) = -TRunoff%erin(nr,nt)-TRunoff%erlateral(nr,nt)
      else
         if(TUnit%areaTotal2(nr)/TUnit%rwidth(nr)/TUnit%rlen(nr) > 1e6_r8) then
            TRunoff%erout(nr,nt) = -TRunoff%erin(nr,nt)-TRunoff%erlateral(nr,nt)
         else
            TRunoff%vr(nr,nt) = CRVRMAN_nosqrt(TUnit%rslpsqrt(nr), TUnit%nr(nr), TRunoff%rr(nr,nt))
            TRunoff%erout(nr,nt) = -TRunoff%vr(nr,nt) * TRunoff%mr(nr,nt)
            if(-TRunoff%erout(nr,nt) > TINYVALUE .and. TRunoff%wr(nr,nt) + &
              (TRunoff%erlateral(nr,nt) + TRunoff%erin(nr,nt) + TRunoff%erout(nr,nt)) * theDeltaT < TINYVALUE) then
               TRunoff%erout(nr,nt) = &
                    -(TRunoff%erlateral(nr,nt) + TRunoff%erin(nr,nt) + TRunoff%wr(nr,nt) / theDeltaT)
               if(TRunoff%mr(nr,nt) > 0._r8) then
                  TRunoff%vr(nr,nt) = -TRunoff%erout(nr,nt) / TRunoff%mr(nr,nt)
               end if
            end if
         end if
      end if

      temp_gwl = TRunoff%qgwl(nr,nt) * TUnit%area(nr) * TUnit%frac(nr)

      TRunoff%dwr(nr,nt) = TRunoff%erlateral(nr,nt) + TRunoff%erin(nr,nt) + TRunoff%erout(nr,nt) + temp_gwl

      if((TRunoff%wr(nr,nt)/theDeltaT &
           + TRunoff%dwr(nr,nt)) < -TINYVALUE) then
         write(iulog,*) 'mosart: ERROR main channel going negative: ', nr, nt
         write(iulog,*) theDeltaT, TRunoff%wr(nr,nt), &
              TRunoff%wr(nr,nt)/theDeltaT, TRunoff%dwr(nr,nt), temp_gwl
         write(iulog,*) ' '
      endif

      ! check for stability
      !    if(TRunoff%vr(nr,nt) < -TINYVALUE .or. TRunoff%vr(nr,nt) > 30) then
      !       write(iulog,*) "Numerical error inRouting_KW, ", nr,nt,TRunoff%vr(nr,nt)
      !    end if

      ! check for negative wr
      !    if(TRunoff%wr(nr,nt) > 1._r8 .and. &
      !      (TRunoff%wr(nr,nt)/theDeltaT + TRunoff%dwr(nr,nt))/TRunoff%wr(nr,nt) < -TINYVALUE) then
      !       write(iulog,*) 'negative wr!', TRunoff%wr(nr,nt), TRunoff%dwr(nr,nt), temp_dwr, temp_gwl, temp_gwl0, theDeltaT
      !       stop
      !    end if

   end subroutine Routing_KW

   !-----------------------------------------------------------------------

   subroutine Routing_MC(nr, nt, theDeltaT)
      !  Muskingum-Cunge routing method

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

   end subroutine Routing_MC

   !-----------------------------------------------------------------------

   subroutine Routing_THREW(nr, nt, theDeltaT)
      !  kinematic wave routing method from THREW model

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

   end subroutine Routing_THREW

   !-----------------------------------------------------------------------

   subroutine Routing_DW(nr, nt, theDeltaT)
      !  classic diffusion wave routing method

      ! Arguments
      integer, intent(in) :: nr, nt
      real(r8), intent(in) :: theDeltaT

   end subroutine Routing_DW

   !-----------------------------------------------------------------------

   subroutine updateState_hillslope(nr,nt)
      !  update the state variables at hillslope

      ! Arguments
      integer, intent(in) :: nr, nt

      TRunoff%yh(nr,nt) = TRunoff%wh(nr,nt) !/ TUnit%area(nr) / TUnit%frac(nr)

   end subroutine updateState_hillslope

   !-----------------------------------------------------------------------

   subroutine updateState_subnetwork(nr,nt)
      !  update the state variables in subnetwork channel

      ! Arguments
      integer, intent(in) :: nr,nt

      if(TUnit%tlen(nr) > 0._r8 .and. TRunoff%wt(nr,nt) > 0._r8) then
         TRunoff%mt(nr,nt) = GRMR(TRunoff%wt(nr,nt), TUnit%tlen(nr))
         TRunoff%yt(nr,nt) = GRHT(TRunoff%mt(nr,nt), TUnit%twidth(nr))
         TRunoff%pt(nr,nt) = GRPT(TRunoff%yt(nr,nt), TUnit%twidth(nr))
         TRunoff%rt(nr,nt) = GRRR(TRunoff%mt(nr,nt), TRunoff%pt(nr,nt))
      else
         TRunoff%mt(nr,nt) = 0._r8
         TRunoff%yt(nr,nt) = 0._r8
         TRunoff%pt(nr,nt) = 0._r8
         TRunoff%rt(nr,nt) = 0._r8
      end if
   end subroutine updateState_subnetwork

   !-----------------------------------------------------------------------

   subroutine updateState_mainchannel(nr, nt)
      !  update the state variables in main channel

      ! Arguments
      integer, intent(in) :: nr, nt

      if(TUnit%rlen(nr) > 0._r8 .and. TRunoff%wr(nr,nt) > 0._r8) then
         TRunoff%mr(nr,nt) = GRMR(TRunoff%wr(nr,nt), TUnit%rlen(nr))
         TRunoff%yr(nr,nt) = GRHR(TRunoff%mr(nr,nt), TUnit%rwidth(nr), TUnit%rwidth0(nr), TUnit%rdepth(nr))
         TRunoff%pr(nr,nt) = GRPR(TRunoff%yr(nr,nt), TUnit%rwidth(nr), TUnit%rwidth0(nr), TUnit%rdepth(nr))
         TRunoff%rr(nr,nt) = GRRR(TRunoff%mr(nr,nt), TRunoff%pr(nr,nt))
      else
         TRunoff%mr(nr,nt) = 0._r8
         TRunoff%yr(nr,nt) = 0._r8
         TRunoff%pr(nr,nt) = 0._r8
         TRunoff%rr(nr,nt) = 0._r8
      end if
   end subroutine updateState_mainchannel

   !-----------------------------------------------------------------------

   function CRVRMAN_nosqrt(sqrtslp_, n_, rr_) result(v_)
      ! Function for calculating channel velocity according to Manning's equation.

      ! Arguments
      real(r8), intent(in) :: sqrtslp_, n_, rr_ ! sqrt(slope), manning's roughness coeff., hydraulic radius
      real(r8)             :: v_            ! v_ is  discharge

      ! Local varaibles
      real(r8) :: ftemp, vtemp

      if(rr_ <= 0._r8) then
         v_ = 0._r8
      else
         v_ = ((rr_*rr_)**(1._r8/3._r8)) * sqrtslp_ / n_
      end if

   end function CRVRMAN_nosqrt

   !-----------------------------------------------------------------------

   function CREHT_nosqrt(sqrthslp_, nh_, Gxr_, yh_) result(eht_)
      ! Function for overland from hillslope into the sub-network channels

      ! Arguments
      real(r8), intent(in) :: sqrthslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
      real(r8)                   :: eht_            ! velocity, specific discharge

      real(r8) :: vh_
      vh_ = CRVRMAN_nosqrt(sqrthslp_,nh_,yh_)
      eht_ = Gxr_*yh_*vh_

   end function CREHT_nosqrt

   !-----------------------------------------------------------------------

   function GRMR(wr_, rlen_) result(mr_)
      ! Function for estimate wetted channel area

      ! Arguments
      real(r8), intent(in) :: wr_, rlen_      ! storage of water, channel length
      real(r8)             :: mr_             ! wetted channel area

      mr_ = wr_ / rlen_
   end function GRMR

   !-----------------------------------------------------------------------

   function GRHT(mt_, twid_) result(ht_)
      ! Function for estimating water depth assuming rectangular channel

      ! Arguments
      real(r8), intent(in) :: mt_, twid_      ! wetted channel area, channel width
      real(r8)             :: ht_             ! water depth

      if(mt_ <= TINYVALUE) then
         ht_ = 0._r8
      else
         ht_ = mt_ / twid_
      end if
   end function GRHT

   !-----------------------------------------------------------------------

   function GRPT(ht_, twid_) result(pt_)
      ! Function for estimating wetted perimeter assuming rectangular channel

      ! Arguments
      real(r8), intent(in) :: ht_, twid_      ! water depth, channel width
      real(r8)             :: pt_             ! wetted perimeter

      if(ht_ <= TINYVALUE) then
         pt_ = 0._r8
      else
         pt_ = twid_ + 2._r8 * ht_
      end if
   end function GRPT

   !-----------------------------------------------------------------------

   function GRRR(mr_, pr_) result(rr_)
      ! Function for estimating hydraulic radius

      ! Arguments
      real(r8), intent(in) :: mr_, pr_        ! wetted area and perimeter
      real(r8)             :: rr_             ! hydraulic radius

      if(pr_ <= TINYVALUE) then
         rr_ = 0._r8
      else
         rr_ = mr_ / pr_
      end if
   end function GRRR

   !-----------------------------------------------------------------------

   function GRHR(mr_, rwidth_, rwidth0_, rdepth_) result(hr_)
      ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
      ! here assuming the channel cross-section consists of three parts, from bottom to up,
      ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
      ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
      ! part 3 is a rectagular with the width rwid0

      ! Arguments
      real(r8), intent(in) :: mr_, rwidth_, rwidth0_, rdepth_ ! wetted channel area, channel width, flood plain wid, water depth
      real(r8)             :: hr_                             ! water depth

      ! Local variables
      real(r8) :: SLOPE1  ! slope of flood plain, TO DO
      real(r8) :: deltamr_

      SLOPE1 = SLOPE1def
      if(mr_ <= TINYVALUE) then
         hr_ = 0._r8
      else
         if(mr_ - rdepth_*rwidth_ <= TINYVALUE) then ! not flooded
            hr_ = mr_/rwidth_
         else ! if flooded, the find out the equivalent depth
            if(mr_ > rdepth_*rwidth_ + (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_-rwidth_)/2._r8)/2._r8 + TINYVALUE) then
               deltamr_ = mr_ - rdepth_*rwidth_ - (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_ - rwidth_)/2._r8)/2._r8;
               hr_ = rdepth_ + SLOPE1*((rwidth0_ - rwidth_)/2._r8) + deltamr_/(rwidth0_);
            else
               deltamr_ = mr_ - rdepth_*rwidth_;
               hr_ = rdepth_ + (-rwidth_+sqrt((rwidth_*rwidth_)+4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
            end if
         end if
      end if
   end function GRHR

   !-----------------------------------------------------------------------

   function GRPR(hr_, rwidth_, rwidth0_,rdepth_) result(pr_)
      ! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
      ! here assuming the channel cross-section consists of three parts, from bottom to up,
      ! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
      ! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
      ! part 3 is a rectagular with the width rwid0

      ! Arguments
      real(r8), intent(in) :: hr_, rwidth_, rwidth0_, rdepth_ ! wwater depth, channel width, flood plain wid, water depth
      real(r8)             :: pr_                             ! water depth

      ! Local variables
      real(r8) :: SLOPE1  ! slope of flood plain, TO DO
      real(r8) :: deltahr_
      real(r8) :: sinatanSLOPE1defr      ! 1.0/sin(atan(slope1))
      logical, save :: first_call = .true.

      SLOPE1 = SLOPE1def
!scs      if (first_call) then
         sinatanSLOPE1defr = 1.0_r8/(sin(atan(SLOPE1def)))
!scs      endif
      first_call = .false.

      if(hr_ < TINYVALUE) then
         pr_ = 0._r8
      else
         if(hr_ <= rdepth_ + TINYVALUE) then ! not flooded
            pr_ = rwidth_ + 2._r8*hr_
         else
            if(hr_ > rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1 + TINYVALUE) then
               deltahr_ = hr_ - rdepth_ - ((rwidth0_-rwidth_)/2._r8)*SLOPE1
               pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1*sinatanSLOPE1defr + deltahr_)
            else
               pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)*sinatanSLOPE1defr)
            end if
         end if
      end if
   end function GRPR

end MODULE MOSART_physics_mod
