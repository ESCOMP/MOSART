module mosart_tspatialunit_type

  ! Topographic and geometric properties, applicable for both grid- and subbasin-based representations

   use shr_kind_mod,      only : r8=>shr_kind_r8, CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
   use shr_sys_mod,       only : shr_sys_abort
   use shr_mpi_mod,       only : shr_mpi_sum, shr_mpi_max
   use shr_string_mod,    only : shr_string_listGetName
   use mosart_io,         only : ncd_pio_openfile, compDOF
   use mosart_vars,       only : mainproc, mpicom_rof, iulog
   use nuopc_shr_methods, only : chkerr
   use ESMF,              only : ESMF_Field, ESMF_RouteHandle, ESMF_Mesh, ESMF_FieldCreate, &
                                 ESMF_FieldSMMStore, ESMF_FieldGet, ESMF_FieldSMM, &
                                 ESMF_SUCCESS, ESMF_TYPEKIND_R8, ESMF_MESHLOC_ELEMENT, ESMF_TERMORDER_SRCSEQ
   use pio,               only : iosystem_desc_t, var_desc_t, io_desc_t, file_desc_t, pio_seterrorhandling, &
                                 pio_inq_varid, pio_inq_vardimid, pio_inq_dimlen, pio_initdecomp, pio_closefile, &
                                 pio_int, pio_double, PIO_INTERNAL_ERROR, pio_read_darray, pio_freedecomp

   implicit none
   private

  type Tspatialunit_type

     ! euler computation
     logical , pointer :: euler_calc(:) ! flag for calculating tracers in euler

     ! frozen water runoff index
     integer :: nice

     ! grid properties
     integer , pointer :: mask(:)      ! mosart mask of mosart cell, 0=null, 1=land with dnID, 2=outlet
     integer , pointer :: ID0(:)
     real(r8), pointer :: lat(:)       ! latitude of the centroid of the cell
     real(r8), pointer :: lon(:)       ! longitude of the centroid of the cell
     real(r8), pointer :: area(:)      ! area of local cell, [m2]
     real(r8), pointer :: areaTotal(:) ! total upstream drainage area, [m2]
     real(r8), pointer :: areaTotal2(:)! computed total upstream drainage area, [m2]
     real(r8), pointer :: rlenTotal(:) ! length of all reaches, [m]
     real(r8), pointer :: Gxr(:)       ! drainage density within the cell, [1/m]
     real(r8), pointer :: frac(:)      ! fraction of cell included in the study area, [-]

     ! hillslope properties
     real(r8), pointer :: nh(:)        ! manning's roughness of the hillslope (channel network excluded)
     real(r8), pointer :: hslp(:)      ! slope of hillslope, [-]
     real(r8), pointer :: hslpsqrt(:)  ! sqrt of slope of hillslope, [-]
     real(r8), pointer :: hlen(:)      ! length of hillslope within the cell, [m]

     ! subnetwork channel properties
     real(r8), pointer :: nt(:)        ! manning's roughness of the subnetwork at hillslope
     real(r8), pointer :: tslp(:)      ! average slope of tributaries, [-]
     real(r8), pointer :: tslpsqrt(:)  ! sqrt of average slope of tributaries, [-]
     real(r8), pointer :: tlen(:)      ! length of all sub-network reach within the cell, [m]
     real(r8), pointer :: twidth(:)    ! bankfull width of the sub-reach, [m]
     real(r8), pointer :: twidth0(:)   ! unadjusted twidth

     ! main channel properties
     real(r8), pointer :: nr(:)        ! manning's roughness of the main reach
     real(r8), pointer :: rlen(:)      ! length of main river reach, [m]
     real(r8), pointer :: rslp(:)      ! slope of main river reach, [-]
     real(r8), pointer :: rslpsqrt(:)  ! sqrt of slope of main river reach, [-]
     real(r8), pointer :: rwidth(:)    ! bankfull width of main reach, [m]
     real(r8), pointer :: rwidth0(:)   ! total width of the flood plain, [m]
     real(r8), pointer :: rdepth(:)    ! bankfull depth of river cross section, [m]
     !
     integer , pointer :: dnID(:)      ! IDs of the downstream units, corresponding to the subbasin ID in the input table
     integer , pointer :: iUp(:,:)     ! IDs of upstream units, corresponding to the subbasin ID in the input table
     integer , pointer :: nUp(:)       ! number of upstream units, maximum 8
     integer , pointer :: indexDown(:) ! indices of the downstream units in the ID array. sometimes subbasins IDs may not be continuous
     integer , pointer :: numDT_r(:)   ! for a main reach, the number of sub-time-steps needed for numerical stability
     integer , pointer :: numDT_t(:)   ! for a subnetwork reach, the number of sub-time-steps needed for numerical stability
     real(r8), pointer :: phi_r(:)     ! the indicator used to define numDT_r
     real(r8), pointer :: phi_t(:)     ! the indicator used to define numDT_t

     ! mapping
     type(ESMF_Field)        :: srcField
     type(ESMF_Field)        :: dstField
     type(ESMF_RouteHandle)  :: rh_direct
     type(ESMF_RouteHandle)  :: rh_eroutUp

  contains

     procedure, public  :: Init
     procedure, private :: set_routehandles
     procedure, private :: set_subtimesteps
     procedure, private :: set_areatotal2

  end type Tspatialunit_type
  public :: Tspatialunit_type

   character(*), parameter :: u_FILE_u = &
        __FILE__
   !-----------------------------------------------------------------------

contains

   !-----------------------------------------------------------------------
   subroutine Init(this, begr, endr, ntracers, nlon, nlat, EMesh, &
        frivinp, IDkey, c_twid, DLevelR, area, gindex, outletg, pio_subsystem, rc)

      ! Arguments
      class(Tspatialunit_type)            :: this
      integer               , intent(in)  :: begr, endr
      integer               , intent(in)  :: ntracers
      real(r8)              , intent(in)  :: area(begr:endr)
      integer               , intent(in)  :: nlon, nlat
      character(len=*)      , intent(in)  :: frivinp
      integer               , intent(in)  :: IDkey(:)
      real(r8)              , intent(in)  :: c_twid(begr:endr)
      integer               , intent(in)  :: DLevelR
      type(iosystem_desc_t) , pointer     :: pio_subsystem
      type(ESMF_Mesh)       , intent(in)  :: Emesh
      integer               , intent(in)  :: gindex(begr:endr)
      integer               , intent(in)  :: outletg(begr:endr)
      integer               , intent(out) :: rc

      ! Local variables
      integer           :: n
      integer           :: ier
      type(file_desc_t) :: ncid          ! pio file desc
      type(var_desc_t)  :: vardesc       ! pio variable desc
      type(io_desc_t)   :: iodesc_dbl    ! pio io desc
      type(io_desc_t)   :: iodesc_int    ! pio io desc
      integer           :: dids(2)       ! variable dimension ids
      integer           :: dsizes(2)     ! variable dimension lengths
      real(r8)          :: hlen_max, rlen_min
      character(len=CS) :: ctemp
      character(len=*),parameter :: FORMI = '(2A,2i10)'
      character(len=*),parameter :: FORMR = '(2A,2g15.7)'
      character(len=*),parameter :: subname = '(mosart_tspatialunit_type_init) '
      !--------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! Read in routing parameters
      call ncd_pio_openfile (ncid, trim(frivinp), 0)
      call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

      ! Setup iodesc based on frac dids
      ier = pio_inq_varid(ncid, name='frac', vardesc=vardesc)
      ier = pio_inq_vardimid(ncid, vardesc, dids)
      ier = pio_inq_dimlen(ncid, dids(1),dsizes(1))
      ier = pio_inq_dimlen(ncid, dids(2),dsizes(2))
      call pio_initdecomp(pio_subsystem, pio_double, dsizes, compDOF, iodesc_dbl)
      call pio_initdecomp(pio_subsystem, pio_int   , dsizes, compDOF, iodesc_int)

      ! Set ice index
      this%nice = ntracers

      ! For now assume that frozen runoff is the last tracer
      ! set euler_calc = false for frozen runoff
      allocate(this%euler_calc(ntracers))
      do n = 1,ntracers
         if (n < ntracers) then
            this%euler_calc(n) = .true.
         else
            this%euler_calc(n) = .false.
         end if
      end do

      allocate(this%frac(begr:endr))
      ier = pio_inq_varid(ncid, name='frac', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%frac, ier)
      if (mainproc) then
         write(iulog,FORMR) trim(subname),' read frac ',minval(this%frac),maxval(this%frac)
      end if

      ! read fdir, convert to mask
      ! fdir <0 ocean, 0=outlet, >0 land
      ! tunit mask is 0=ocean, 1=land, 2=outlet for mosart calcs

      allocate(this%mask(begr:endr))
      ier = pio_inq_varid(ncid, name='fdir', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_int, this%mask, ier)
      if (mainproc) then
         write(iulog,'(2A,2i10)') trim(subname),' read fdir mask ',minval(this%mask),maxval(this%mask)
      end if

      do n = begr, endr
         if (this%mask(n) < 0) then
            this%mask(n) = 0
         elseif (this%mask(n) == 0) then
            this%mask(n) = 2
            if (abs(this%frac(n)-1.0_r8)>1.0e-9) then
               write(iulog,*) subname,' ERROR frac ne 1.0',n,this%frac(n)
               call shr_sys_abort(subname//' ERROR frac ne 1.0')
            endif
         elseif (this%mask(n) > 0) then
            this%mask(n) = 1
            if (abs(this%frac(n)-1.0_r8)>1.0e-9) then
               write(iulog,*) subname,' ERROR frac ne 1.0',n,this%frac(n)
               call shr_sys_abort(subname//' ERROR frac ne 1.0')
            endif
         else
            call shr_sys_abort(subname//' this mask error')
         endif
      enddo

      allocate(this%ID0(begr:endr))
      ier = pio_inq_varid(ncid, name='ID', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_int, this%ID0, ier)
      if (mainproc) write(iulog,'(2A,2i10)') trim(subname),' read ID0 ',minval(this%ID0),maxval(this%ID0)

      allocate(this%dnID(begr:endr))
      ier = pio_inq_varid(ncid, name='dnID', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_int, this%dnID, ier)
      if (mainproc) write(iulog,'(2A,2i10)') trim(subname),' read dnID ',minval(this%dnID),maxval(this%dnID)

      ! RESET ID0 and dnID indices using the IDkey to be consistent with standard gindex order
      do n=begr, endr
         this%ID0(n)  = IDkey(this%ID0(n))
         if (this%dnID(n) > 0 .and. this%dnID(n) <= nlon*nlat) then
            if (IDkey(this%dnID(n)) > 0 .and. IDkey(this%dnID(n)) <= nlon*nlat) then
               this%dnID(n) = IDkey(this%dnID(n))
            else
               write(iulog,*) subname,' ERROR bad IDkey for this%dnID',n,this%dnID(n),IDkey(this%dnID(n))
               call shr_sys_abort(subname//' ERROR bad IDkey for this%dnID')
            endif
         endif
      enddo

      allocate(this%area(begr:endr))
      ier = pio_inq_varid(ncid, name='area', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%area, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read area ',minval(this%area),maxval(this%area)

      do n=begr, endr
         if (this%area(n) < 0._r8) this%area(n) = area(n)
         if (this%area(n) /= area(n)) then
            write(iulog,*) subname,' ERROR area mismatch',this%area(n),area(n)
            call shr_sys_abort(subname//' ERROR area mismatch')
         endif
      enddo

      allocate(this%areaTotal(begr:endr))
      ier = pio_inq_varid(ncid, name='areaTotal', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%areaTotal, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read areaTotal ',minval(this%areaTotal),maxval(this%areaTotal)

      allocate(this%rlenTotal(begr:endr))
      this%rlenTotal = 0._r8

      allocate(this%nh(begr:endr))
      ier = pio_inq_varid(ncid, name='nh', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%nh, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read nh ',minval(this%nh),maxval(this%nh)

      allocate(this%hslp(begr:endr))
      ier = pio_inq_varid(ncid, name='hslp', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%hslp, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read hslp ',minval(this%hslp),maxval(this%hslp)

      allocate(this%hslpsqrt(begr:endr))
      this%hslpsqrt = 0._r8

      allocate(this%gxr(begr:endr))
      ier = pio_inq_varid(ncid, name='gxr', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%gxr, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read gxr ',minval(this%gxr),maxval(this%gxr)

      allocate(this%hlen(begr:endr))
      this%hlen = 0._r8

      allocate(this%tslp(begr:endr))
      ier = pio_inq_varid(ncid, name='tslp', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%tslp, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read tslp ',minval(this%tslp),maxval(this%tslp)

      allocate(this%tslpsqrt(begr:endr))
      this%tslpsqrt = 0._r8

      allocate(this%tlen(begr:endr))
      this%tlen = 0._r8

      allocate(this%twidth(begr:endr))
      ier = pio_inq_varid(ncid, name='twid', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%twidth, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read twidth ',minval(this%twidth),maxval(this%twidth)

      ! save twidth before adjusted below
      allocate(this%twidth0(begr:endr))
      this%twidth0(begr:endr)=this%twidth(begr:endr)

      allocate(this%nt(begr:endr))
      ier = pio_inq_varid(ncid, name='nt', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%nt, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read nt ',minval(this%nt),maxval(this%nt)

      allocate(this%rlen(begr:endr))
      ier = pio_inq_varid(ncid, name='rlen', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%rlen, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read rlen ',minval(this%rlen),maxval(this%rlen)

      allocate(this%rslp(begr:endr))
      ier = pio_inq_varid(ncid, name='rslp', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%rslp, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read rslp ',minval(this%rslp),maxval(this%rslp)

      allocate(this%rslpsqrt(begr:endr))
      this%rslpsqrt = 0._r8

      allocate(this%rwidth(begr:endr))
      ier = pio_inq_varid(ncid, name='rwid', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%rwidth, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read rwidth ',minval(this%rwidth),maxval(this%rwidth)

      allocate(this%rwidth0(begr:endr))
      ier = pio_inq_varid(ncid, name='rwid0', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%rwidth0, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read rwidth0 ',minval(this%rwidth0),maxval(this%rwidth0)

      allocate(this%rdepth(begr:endr))
      ier = pio_inq_varid(ncid, name='rdep', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%rdepth, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read rdepth ',minval(this%rdepth),maxval(this%rdepth)

      allocate(this%nr(begr:endr))
      ier = pio_inq_varid(ncid, name='nr', vardesc=vardesc)
      call pio_read_darray(ncid, vardesc, iodesc_dbl, this%nr, ier)
      if (mainproc) write(iulog,FORMR) trim(subname),' read nr ',minval(this%nr),maxval(this%nr)

      allocate(this%nUp(begr:endr))
      this%nUp = 0
      allocate(this%iUp(begr:endr,8))
      this%iUp = 0
      allocate(this%indexDown(begr:endr))
      this%indexDown = 0

      ! control parameters and some other derived parameters
      ! estimate derived input variables

      ! add minimum value to rlen (length of main channel); rlen values can
      ! be too small, leading to tlen values that are too large

      do n=begr,endr
         rlen_min = sqrt(this%area(n))
         if(this%rlen(n) < rlen_min) then
            this%rlen(n) = rlen_min
         end if
      end do

      do n=begr,endr
         if(this%Gxr(n) > 0._r8) then
            this%rlenTotal(n) = this%area(n)*this%Gxr(n)
         end if
      end do

      do n=begr,endr
         if(this%rlen(n) > this%rlenTotal(n)) then
            this%rlenTotal(n) = this%rlen(n)
         end if
      end do

      do n=begr,endr

         if(this%rlen(n) > 0._r8) then
            this%hlen(n) = this%area(n) / this%rlenTotal(n) / 2._r8

            ! constrain hlen (hillslope length) values based on cell area
            hlen_max = max(1000.0_r8, sqrt(this%area(n)))
            if(this%hlen(n) > hlen_max) then
               this%hlen(n) = hlen_max   ! allievate the outlier in drainag\e density estimation. TO DO
            end if

            this%tlen(n) = this%area(n) / this%rlen(n) / 2._r8 - this%hlen(n)

            if (this%twidth(n) < 0._r8) then
               this%twidth(n) = 0._r8
            end if
            if ( this%tlen(n) > 0._r8 .and. &
                 (this%rlenTotal(n)-this%rlen(n))/this%tlen(n) > 1._r8 ) then
               this%twidth(n) = c_twid(n)*this%twidth(n) * &
                    ((this%rlenTotal(n)-this%rlen(n))/this%tlen(n))
            end if
            if (this%tlen(n) > 0._r8 .and. this%twidth(n) <= 0._r8) then
               this%twidth(n) = 0._r8
            end if
         else
            this%hlen(n) = 0._r8
            this%tlen(n) = 0._r8
            this%twidth(n) = 0._r8
         end if
         if(this%rslp(n) <= 0._r8) then
            this%rslp(n) = 0.0001_r8
         end if
         if(this%tslp(n) <= 0._r8) then
            this%tslp(n) = 0.0001_r8
         end if
         if(this%hslp(n) <= 0._r8) then
            this%hslp(n) = 0.005_r8
         end if

         this%rslpsqrt(n) = sqrt(this%rslp(n))
         this%tslpsqrt(n) = sqrt(this%tslp(n))
         this%hslpsqrt(n) = sqrt(this%hslp(n))
      end do

      call pio_freedecomp(ncid, iodesc_dbl)
      call pio_freedecomp(ncid, iodesc_int)
      call pio_closefile(ncid)

      ! Create srcfield and dstfield - needed for mapping
      this%srcfield = ESMF_FieldCreate(EMesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
           ungriddedLBound=(/1/), ungriddedUBound=(/ntracers/), gridToFieldMap=(/2/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      this%dstfield = ESMF_FieldCreate(EMesh, ESMF_TYPEKIND_R8, meshloc=ESMF_MESHLOC_ELEMENT, &
           ungriddedLBound=(/1/), ungriddedUBound=(/ntracers/), gridToFieldMap=(/2/), rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Create route handles
      call this%set_routehandles(begr, endr, gindex, outletg, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Compute areatotal2
      ! this basically advects upstream areas downstream and
      ! adds them up as it goes until all upstream areas are accounted for
      allocate(this%areatotal2(begr:endr))
      call this%set_areatotal2(begr, endr, nlon, nlat, area, rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      ! Determine subcycling time steps
      allocate(this%numDT_r(begr:endr))
      allocate(this%numDT_t(begr:endr))
      allocate(this%phi_r(begr:endr))
      allocate(this%phi_t(begr:endr))
      call this%set_subtimesteps(begr, endr, DLevelR)

   end subroutine Init

   !-----------------------------------------------------------------------

   subroutine set_routehandles(this, begr, endr, gindex, outletg, rc)

      ! Arguments
      class(Tspatialunit_type) :: this
      integer , intent(in)     :: begr, endr
      integer , intent(in)     :: gindex(begr:endr)
      integer , intent(in)     :: outletg(begr:endr)
      integer , intent(out)    :: rc

      ! Local variables
      integer               :: nn, n, cnt, nr, nt
      real(r8), pointer     :: src_direct(:,:)
      real(r8), pointer     :: dst_direct(:,:)
      real(r8), pointer     :: src_eroutUp(:,:)
      real(r8), pointer     :: dst_eroutUp(:,:)
      real(r8), allocatable :: factorList(:)
      integer , allocatable :: factorIndexList(:,:)
      integer               :: srcTermProcessing_Value = 0
      !--------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! ---------------------------------------
      ! Calculate map for direct to outlet mapping
      ! ---------------------------------------

      ! Set up pointer arrays into srcfield and dstfield
      call ESMF_FieldGet(this%srcfield, farrayPtr=src_direct, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(this%dstfield, farrayPtr=dst_direct, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      src_direct(:,:) = 0._r8
      dst_direct(:,:) = 0._r8

      ! The route handle rh_direct will then be used in mosart_run
      cnt = endr - begr + 1
      allocate(factorList(cnt))
      allocate(factorIndexList(2,cnt))
      cnt = 0
      do nr = begr,endr
         cnt = cnt + 1
         if (outletg(nr) > 0) then
            factorList(cnt) = 1.0_r8
            factorIndexList(1,cnt) = gindex(nr)
            factorIndexList(2,cnt) = outletg(nr)
         else
            factorList(cnt) = 1.0_r8
            factorIndexList(1,cnt) = gindex(nr)
            factorIndexList(2,cnt) = gindex(nr)
         endif
      enddo

      call ESMF_FieldSMMStore(this%srcField, this%dstField, this%rh_direct, factorList, factorIndexList, &
           ignoreUnmatchedIndices=.true., srcTermProcessing=srcTermProcessing_Value, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      deallocate(factorList)
      deallocate(factorIndexList)

      if (mainproc) write(iulog,*) " Done initializing rh_direct "

      ! ---------------------------------------
      ! Compute map rh_eroutUp
      ! ---------------------------------------

      ! Set up pointer arrays into srcfield and dstfield
      call ESMF_FieldGet(this%srcfield, farrayPtr=src_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(this%dstfield, farrayPtr=dst_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      src_eroutUp(:,:) = 0._r8
      dst_eroutUp(:,:) = 0._r8

      cnt = 0
      do nr = begr,endr
         if (this%dnID(nr) > 0) then
            cnt = cnt + 1
         end if
      end do
      allocate(factorList(cnt))
      allocate(factorIndexList(2,cnt))
      cnt = 0
      do nr = begr,endr
         if (this%dnID(nr) > 0) then
            cnt = cnt + 1
            factorList(cnt) = 1.0_r8
            factorIndexList(1,cnt) = this%ID0(nr)
            factorIndexList(2,cnt) = this%dnID(nr)
         endif
      enddo
      if (mainproc) write(iulog,*) " Done initializing rh_eroutUp"

      call ESMF_FieldSMMStore(this%srcfield, this%dstfield, this%rh_eroutUp, factorList, factorIndexList, &
           ignoreUnmatchedIndices=.true., srcTermProcessing=srcTermProcessing_Value, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

      deallocate(factorList)
      deallocate(factorIndexList)

   end subroutine set_routehandles

   !-----------------------------------------------------------------------

   subroutine set_areatotal2(this, begr, endr, nlon, nlat, area, rc)

      ! Arguments
      class(Tspatialunit_type) :: this
      integer  , intent(in)    :: begr, endr
      integer  , intent(in)    :: nlon,nlat
      real(r8) , intent(in)    :: area(begr:endr)
      integer  , intent(out)   :: rc

      ! Local variables
      integer           :: nr, cnt, tcnt ! indices
      real(r8)          :: areatot_prev, areatot_tmp, areatot_new
      real(r8), pointer :: src_direct(:,:)
      real(r8), pointer :: dst_direct(:,:)
      real(r8), pointer :: src_eroutUp(:,:)
      real(r8), pointer :: dst_eroutUp(:,:)
      character(len=*),parameter :: subname = '(mosart_tspatialunit_type_set_areatotal2) '
      ! --------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! ---------------------------------------
      ! compute areatot from area using dnID
      ! ---------------------------------------

      ! Set up pointer arrays into srcfield and dstfield
      call ESMF_FieldGet(this%srcfield, farrayPtr=src_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      call ESMF_FieldGet(this%dstfield, farrayPtr=dst_eroutUp, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      src_eroutUp(:,:) = 0._r8
      dst_eroutUp(:,:) = 0._r8

      ! this basically advects upstream areas downstream and
      ! adds them up as it goes until all upstream areas are accounted for

      this%areatotal2(:) = 0._r8

      ! initialize dst_eroutUp to local area and add that to areatotal2
      cnt = 0
      dst_eroutUp(:,:) = 0._r8
      do nr = begr,endr
         cnt = cnt + 1
         dst_eroutUp(1,cnt) = area(nr)
         this%areatotal2(nr) = area(nr)
      enddo

      tcnt = 0
      areatot_prev = -99._r8
      areatot_new = -50._r8
      do while (areatot_new /= areatot_prev .and. tcnt < nlon*nlat)

         tcnt = tcnt + 1

         ! copy dst_eroutUp to src_eroutUp for next downstream step
         src_eroutUp(:,:) = 0._r8
         cnt = 0
         do nr = begr,endr
            cnt = cnt + 1
            src_eroutUp(1,cnt) = dst_eroutUp(1,cnt)
         enddo

         dst_eroutUp(:,:) = 0._r8
         call ESMF_FieldSMM(this%srcfield, this%dstField, this%rh_eroutUp, termorderflag=ESMF_TERMORDER_SRCSEQ, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return

         ! add dst_eroutUp to areatot and compute new global sum
         cnt = 0
         areatot_prev = areatot_new
         areatot_tmp = 0._r8
         do nr = begr,endr
            cnt = cnt + 1
            this%areatotal2(nr) = this%areatotal2(nr) + dst_eroutUp(1,cnt)
            areatot_tmp = areatot_tmp + this%areatotal2(nr)
         enddo
         call shr_mpi_sum(areatot_tmp, areatot_new, mpicom_rof, 'areatot_new', all=.true.)

         if (mainproc) then
            write(iulog,*) trim(subname),' areatot calc ',tcnt,areatot_new
         endif
      enddo

      if (areatot_new /= areatot_prev) then
         write(iulog,*) trim(subname),' MOSART ERROR: areatot incorrect ',areatot_new, areatot_prev
         call shr_sys_abort(trim(subname)//' MOSART ERROR areatot incorrect')
      endif

   end subroutine set_areatotal2

   !-----------------------------------------------------------------------

   subroutine set_subtimesteps(this, begr, endr, DLevelR)

      ! Set the sub-time-steps for channel routing

      ! Arguments
      class(Tspatialunit_type) :: this
      integer, intent(in) :: begr, endr
      integer, intent(in) :: DLevelR

      ! Local variables
      integer :: nr   !local index
      integer :: numDT_r, numDT_t
      character(len=*),parameter :: subname = '(mosart_tspatialunit_type_subtimestep) '
      ! --------------------------------------------------------------

      this%numDT_r(:) = 1
      this%numDT_t(:) = 1
      this%phi_r(:) = 0._r8
      this%phi_t(:) = 0._r8

      do nr = begr,endr
         if (this%mask(nr) > 0 .and. this%rlen(nr) > 0._r8) then
            this%phi_r(nr) = this%areaTotal2(nr)*sqrt(this%rslp(nr))/(this%rlen(nr)*this%rwidth(nr))
            if (this%phi_r(nr) >= 10._r8) then
               this%numDT_r(nr) = (this%numDT_r(nr)*log10(this%phi_r(nr))*DLevelR) + 1
            else
               this%numDT_r(nr) = this%numDT_r(nr)*1.0_r8*DLevelR + 1
            end if
         end if
         if (this%numDT_r(nr) < 1) this%numDT_r(nr) = 1

         if (this%tlen(nr) > 0._r8) then
            this%phi_t(nr) = this%area(nr)*sqrt(this%tslp(nr))/(this%tlen(nr)*this%twidth(nr))
            if (this%phi_t(nr) >= 10._r8) then
               this%numDT_t(nr) = (this%numDT_t(nr)*log10(this%phi_t(nr))*DLevelR) + 1
            else
               this%numDT_t(nr) = (this%numDT_t(nr)*1.0*DLevelR) + 1
            end if
         end if
         if (this%numDT_t(nr) < 1) this%numDT_t(nr) = 1
      end do

      call shr_mpi_max(maxval(this%numDT_r),numDT_r,mpicom_rof,'numDT_r',all=.false.)
      call shr_mpi_max(maxval(this%numDT_t),numDT_t,mpicom_rof,'numDT_t',all=.false.)
      if (mainproc) then
         write(iulog,*) subname,' DLevelR = ',DlevelR
         write(iulog,*) subname,' numDT_r   = ',minval(this%numDT_r),maxval(this%numDT_r)
         write(iulog,*) subname,' numDT_r max  = ',numDT_r
         write(iulog,*) subname,' numDT_t   = ',minval(this%numDT_t),maxval(this%numDT_t)
         write(iulog,*) subname,' numDT_t max  = ',numDT_t
      endif

   end subroutine set_subtimesteps

end module mosart_tspatialunit_type
