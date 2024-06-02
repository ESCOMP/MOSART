module mosart_control_type

  use shr_kind_mod,      only : r8 => shr_kind_r8, CL => SHR_KIND_CL, CS => SHR_KIND_CS
  use shr_sys_mod,       only : shr_sys_abort
  use shr_const_mod,     only : shr_const_pi, shr_const_rearth
  use shr_string_mod,    only : shr_string_listGetNum, shr_string_listGetName
  use shr_mpi_mod,       only : shr_mpi_sum, shr_mpi_max
  use mosart_io,         only : ncd_io, ncd_pio_openfile, ncd_pio_closefile
  use mosart_vars,       only : mainproc, iam, npes, mpicom_rof, iulog, spval, re
  use pio,               only : file_desc_t, PIO_BCAST_ERROR, pio_seterrorhandling
  use ESMF,              only : ESMF_DistGrid, ESMF_Array, ESMF_RouteHandle, ESMF_SUCCESS, &
                                ESMF_DistGridCreate, ESMF_ArrayCreate, ESMF_ArrayHaloStore, &
                                ESMF_ArrayHalo, ESMF_ArrayGet
  use perf_mod,          only : t_startf, t_stopf
  use nuopc_shr_methods, only : chkerr

  implicit none
  private

  type control_type

     ! grid sizes
     integer :: lnumr                                ! local  number of cells
     integer :: numr                                 ! global number of cells
     integer :: nlon = -999                          ! number of longitudes
     integer :: nlat = -999                          ! number of latitudes

     ! tracers
     integer :: ntracers = -999                      ! number of tracers
     character(len=3), allocatable :: tracer_names(:)! tracer names
     integer :: nt_liq                               ! index of liquid tracer in tracer_names
     integer :: nt_ice                               ! index of ice tracer in tracer_names

     ! decomp info
     integer :: begr                                 ! local start index
     integer :: endr                                 ! local stop indices
     integer , pointer :: gindex(:) => null()        ! global index consistent with map file
     type(ESMF_DistGrid) :: distgrid                 ! esmf global index space descriptor

     ! grid
     real(r8), pointer :: rlon(:) => null()          ! longitude list, 1d
     real(r8), pointer :: rlat(:) => null()          ! latitude list, 1d
     real(r8), pointer :: lonc(:) => null()          ! lon of cell
     real(r8), pointer :: latc(:) => null()          ! lat of cell
     integer , pointer :: dsig(:) => null()          ! downstream index, global index
     integer , pointer :: outletg(:) => null()       ! outlet index, global index
     real(r8), pointer :: area(:) => null()          ! area of cell
     integer , pointer :: mask(:) => null()          ! general mask of cell 1=land, 2=ocean, 3=outlet
     real(r8)          :: totarea                    ! global area

     ! inputs to MOSART
     real(r8), pointer :: qsur(:,:) => null()        ! surface runoff from coupler [m3/s] (lnd)
     real(r8), pointer :: qsub(:,:) => null()        ! subsurfacer runoff from coupler [m3/s] (lnd)
     real(r8), pointer :: qgwl(:,:) => null()        ! glacier/wetland/lake runoff from coupler [m3/s] (lnd)
     real(r8), pointer :: qirrig(:) => null()        ! irrigation flow from coupler [m3/s]
     real(r8), pointer :: qglc_liq(:) => null()      ! glacier liquid runoff from coupler [m3/s] (glc)
     real(r8), pointer :: qglc_ice(:) => null()      ! glacier ice runoff from coupler [m3/s] (glc)

     ! outputs from MOSART
     real(r8), pointer :: flood(:) => null()         ! flood water to coupler [m3/s] (lnd)
     real(r8), pointer :: runoff(:,:) => null()      ! runoff (from outlet to reach) to coupler [m3/s]
     real(r8), pointer :: direct(:,:) => null()      ! direct flow to outlet from land input [m3/s]
     real(r8), pointer :: qirrig_actual(:) => null() ! minimum of irrigation and available main channel storage [m3/s]
     real(r8), pointer :: direct_glc(:,:) =>null()   ! direct flow to outlet from glc input [m3/s]

     ! storage, runoff
     real(r8), pointer :: runofflnd(:,:) => null()   ! runoff masked for land [m3/s]
     real(r8), pointer :: runoffocn(:,:) => null()   ! runoff masked for ocn  [m3/s]
     real(r8), pointer :: runofftot(:,:) => null()   ! total runoff masked for ocn  [m3/s]
     real(r8), pointer :: dvolrdt(:,:) => null()     ! change in storage (mm/s)
     real(r8), pointer :: dvolrdtlnd(:,:) => null()  ! dvolrdt masked for land (mm/s)
     real(r8), pointer :: dvolrdtocn(:,:) => null()  ! dvolrdt masked for ocn  (mm/s)
     real(r8), pointer :: volr(:,:) => null()        ! storage (m3)
     real(r8), pointer :: fthresh(:) => null()       ! water flood threshold

     ! flux variables
     real(r8), pointer :: flow(:,:) => null()        ! stream flow out of gridcell (m3/s)
     real(r8), pointer :: evel(:,:) => null()        ! effective tracer velocity (m/s) NOT_USED
     real(r8), pointer :: erout_prev(:,:) => null()  ! erout previous timestep (m3/s)
     real(r8), pointer :: eroutup_avg(:,:) => null() ! eroutup average over coupling period (m3/s)
     real(r8), pointer :: erlat_avg(:,:) => null()   ! erlateral average over coupling period (m3/s)
     real(r8), pointer :: effvel(:) => null()        ! effective velocity for a tracer NOT_USED

     ! halo operations
     type(ESMF_RouteHandle) :: haloHandle
     type(ESMF_Array)       :: fld_halo_array
     type(ESMF_Array)       :: lon_halo_array
     type(ESMF_Array)       :: lat_halo_array
     integer , pointer      :: halo_arrayptr_index(:,:) => null() ! index into halo_arrayptr that corresponds to a halo point
     real(r8), pointer      :: fld_halo_arrayptr(:) => null()         ! preallocated memory for exclusive region plus halo
     real(r8), pointer      :: lon_halo_arrayptr(:) => null()     ! preallocated memory for exclusive region plus halo
     real(r8), pointer      :: lat_halo_arrayptr(:) => null()     ! preallocated memory for exclusive region plus halo

   contains

     procedure, public  :: Init
     procedure, public  :: init_tracer_names
     procedure, private :: init_decomp
     procedure, private :: test_halo
     procedure, public  :: calc_gradient

  end type control_type
  public :: control_type

  private :: init_decomp
  public  :: calc_gradient

#ifdef NDEBUG
  integer,parameter :: dbug = 0 ! 0 = none, 1=normal, 2=much, 3=max
#else
  integer,parameter :: dbug = 3 ! 0 = none, 1=normal, 2=much, 3=max
#endif

  integer, public :: max_num_halo = 8
  ! eight surrounding indices ordered as [N,NE,E,SE,S,SW,W,NW]
  integer, public :: halo_n  = 1
  integer, public :: halo_ne = 2
  integer, public :: halo_e  = 3
  integer, public :: halo_se = 4
  integer, public :: halo_s  = 5
  integer, public :: halo_sw = 6
  integer, public :: halo_w  = 7
  integer, public :: halo_nw = 8

  ! The following are set from

  character(*), parameter :: u_FILE_u = &
       __FILE__

!========================================================================
contains
!========================================================================

  subroutine init_tracer_names(this, mosart_tracers)

    ! Arguments
    class(control_type) :: this
    character(len=CS)   :: mosart_tracers    ! colon delimited string of tracer names

    ! Local variables
    integer :: nt
    character(len=*),parameter :: subname = '(mosart_control_type: init_tracer_names)'
    !-----------------------------------------------------------------------

    ! Determine number of tracers and array of tracer names
    this%ntracers = shr_string_listGetNum(mosart_tracers)
    allocate(this%tracer_names(this%ntracers))
    do nt = 1,this%ntracers
      call shr_string_listGetName(mosart_tracers, nt, this%tracer_names(nt))
    end do

    ! Set tracers
    this%nt_liq = 0
    this%nt_ice = 0
    do nt = 1,this%ntracers
       if (trim(this%tracer_names(nt)) == 'LIQ') this%nt_liq = nt
       if (trim(this%tracer_names(nt)) == 'ICE') this%nt_ice = nt
    enddo
    if (this%nt_liq == 0 .or. this%nt_ice == 0) then
       write(iulog,*) trim(subname),': ERROR in tracers LIQ ICE ',this%nt_liq,this%nt_ice,this%tracer_names(:)
       call shr_sys_abort()
    endif

  end subroutine init_tracer_names


  !========================================================================
  subroutine Init(this, locfn, decomp_option, use_halo_option, IDkey, rc)

    ! Arguments
    class(control_type)            :: this
    character(len=*) , intent(in)  :: locfn
    character(len=*) , intent(in)  :: decomp_option   ! decomposition option
    logical          , intent(in)  :: use_halo_option ! create ESMF array and route handle for halos
    integer          , intent(out) :: IDkey(:)        ! translation key from ID to gindex
    integer          , intent(out) :: rc

    ! Local variables
    real(r8)          :: area_global(this%nlon*this%nlat) ! area
    real(r8)          :: tempr(this%nlon,this%nlat)       ! temporary buffer
    real(r8)          :: rlats(this%nlat)                 ! latitude of 1d south grid cell edge (deg)
    real(r8)          :: rlatn(this%nlat)                 ! latitude of 1d north grid cell edge (deg)
    real(r8)          :: rlonw(this%nlon)                 ! longitude of 1d west grid cell edge (deg)
    real(r8)          :: rlone(this%nlon)                 ! longitude of 1d east grid cell edge (deg)
    real(r8)          :: larea                            ! tmp local sum of area
    real(r8)          :: deg2rad                          ! pi/180
    integer           :: g, n, i, j, nr, nt               ! iterators
    real(r8)          :: edgen                            ! North edge of the direction file
    real(r8)          :: edgee                            ! East edge of the direction file
    real(r8)          :: edges                            ! South edge of the direction file
    real(r8)          :: edgew                            ! West edge of the direction file
    real(r8)          :: dx                               ! lon dist. betn grid cells (m)
    real(r8)          :: dy                               ! lat dist. betn grid cells (m)
    type(file_desc_t) :: ncid                             ! pio file desc
    logical           :: found                            ! flag
    integer           :: ntracers                         ! used to simplify code
    integer           :: ier                              ! error status
    integer           :: begr, endr                       ! used to simplify code
    integer           :: nlon,nlat
    real(r8)          :: effvel0 = 10.0_r8                ! default velocity (m/s)
    character(len=*),parameter :: subname = '(mosart_control_type: Init)'
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    nlon = this%nlon
    nlat = this%nlat

    !---------------------------------------
    ! Read the routing parameters
    !---------------------------------------

    call ncd_pio_openfile (ncid, trim(locfn), 0)
    call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)

    call ncd_io(ncid=ncid, varname='longxy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read mosart longitudes')
    if (mainproc) write(iulog,*) 'Read longxy ',minval(tempr),maxval(tempr)
    allocate(this%rlon(this%nlon))
    do i=1,nlon
       this%rlon(i) = tempr(i,1)
    enddo
    if (mainproc) write(iulog,*) 'rlon center ',minval(this%rlon),maxval(this%rlon)

    call ncd_io(ncid=ncid, varname='latixy', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read mosart latitudes')
    if (mainproc) write(iulog,*) 'Read latixy ',minval(tempr),maxval(tempr)
    allocate(this%rlat(this%nlat))
    do j=1,this%nlat
       this%rlat(j) = tempr(1,j)
    end do
    if (mainproc) write(iulog,*) 'rlat center ',minval(this%rlat),maxval(this%rlat)

    call ncd_io(ncid=ncid, varname='area', flag='read', data=tempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read mosart area')
    if (mainproc) write(iulog,*) 'Read area ',minval(tempr),maxval(tempr)
    do j=1,this%nlat
       do i=1,nlon
          n = (j-1)*nlon + i
          area_global(n) = tempr(i,j)
       end do
    end do
    if (mainproc) write(iulog,*) 'area ',minval(area_global),maxval(area_global)
    call ncd_pio_closefile(ncid)

    !-------------------------------------------------------
    ! adjust area estimation from DRT algorithm for those outlet grids
    ! useful for grid-based representation only
    ! need to compute areas where they are not defined in input file
    !-------------------------------------------------------

    ! Derive gridbox edges
    ! assuming equispaced grid, calculate edges from nlat/nlon
    ! w/o assuming a global grid
    edgen = maxval(this%rlat) + 0.5*abs(this%rlat(1) - this%rlat(2))
    edges = minval(this%rlat) - 0.5*abs(this%rlat(1) - this%rlat(2))
    edgee = maxval(this%rlon) + 0.5*abs(this%rlon(1) - this%rlon(2))
    edgew = minval(this%rlon) - 0.5*abs(this%rlon(1) - this%rlon(2))
    if (edgen .ne.  90._r8)then
       if (mainproc ) write(iulog,*) 'Regional grid: edgen = ', edgen
    end if
    if (edges .ne. -90._r8)then
       if (mainproc ) write(iulog,*) 'Regional grid: edges = ', edges
    end if
    if (edgee .ne. 180._r8)then
       if (mainproc ) write(iulog,*) 'Regional grid: edgee = ', edgee
    end if
    if (edgew .ne.-180._r8)then
       if ( mainproc ) write(iulog,*) 'Regional grid: edgew = ', edgew
    end if

    ! Set edge latitudes (assumes latitudes are constant for a given longitude)
    rlats(:) = edges
    rlatn(:) = edgen
    do j = 2, nlat
       if (this%rlat(2) > this%rlat(1)) then ! South to North grid
          rlats(j)   = (this%rlat(j-1) + this%rlat(j)) / 2._r8
          rlatn(j-1) = rlats(j)
       else  ! North to South grid
          rlatn(j)   = (this%rlat(j-1) + this%rlat(j)) / 2._r8
          rlats(j-1) = rlatn(j)
       end if
    end do

    ! Set edge longitudes
    rlonw(:) = edgew
    rlone(:) = edgee
    dx = (edgee - edgew) / nlon
    do i = 2, nlon
       rlonw(i)   = rlonw(i) + (i-1)*dx
       rlone(i-1) = rlonw(i)
    end do

    ! adjust area estimation from DRT algorithm for those outlet grids
    deg2rad = shr_const_pi / 180._r8
    do n=1,nlon*nlat
       if (area_global(n) <= 0._r8) then
          i = mod(n-1,nlon) + 1
          j = (n-1)/nlon + 1
          dx = (rlone(i) - rlonw(i)) * deg2rad
          dy = sin(rlatn(j)*deg2rad) - sin(rlats(j)*deg2rad)
          area_global(n) = abs(1.e6_r8 * dx*dy*re*re)
          if (mainproc .and. area_global(n) <= 0) then
             write(iulog,*) 'Warning! Zero area for unit ', n, area_global(n),dx,dy,re
          end if
       end if
    end do

    ! ---------------------------------------------
    ! Determine decomposition
    ! ---------------------------------------------

    ! memory for this%gindex, this%mask and this%dsig is allocated in init_decomp
    call t_startf('mosarti_decomp')
    call this%init_decomp(locfn, decomp_option, use_halo_option, &
         nlon, nlat, this%begr, this%endr, this%lnumr, this%numr, IDkey, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('mosarti_decomp')

    ! ---------------------------------------------
    ! Allocate and initialize remaining variables
    ! ---------------------------------------------

    begr = this%begr
    endr = this%endr
    ntracers = this%ntracers

    allocate(this%area(begr:endr),            &
         !
         this%volr(begr:endr,ntracers),       &
         this%dvolrdt(begr:endr,ntracers),    &
         this%dvolrdtlnd(begr:endr,ntracers), &
         this%dvolrdtocn(begr:endr,ntracers), &
         !
         this%runoff(begr:endr,ntracers),     &
         this%runofflnd(begr:endr,ntracers),  &
         this%runoffocn(begr:endr,ntracers),  &
         this%runofftot(begr:endr,ntracers),  &
         !
         this%fthresh(begr:endr),             &
         this%flood(begr:endr),               &
         !
         this%direct(begr:endr,ntracers),     &
         this%qsur(begr:endr,ntracers),       &
         this%qsub(begr:endr,ntracers),       &
         this%qgwl(begr:endr,ntracers),       &
         this%qirrig(begr:endr),              &
         this%qirrig_actual(begr:endr),       &
         this%qglc_liq(begr:endr),            &
         this%qglc_ice(begr:endr),            &
         !
         this%evel(begr:endr,ntracers),       &
         this%flow(begr:endr,ntracers),       &
         this%erout_prev(begr:endr,ntracers), &
         this%eroutup_avg(begr:endr,ntracers),&
         this%erlat_avg(begr:endr,ntracers),  &
         !
         this%effvel(ntracers),               &
         this%direct_glc(begr:endr,2), &
         stat=ier)
    if (ier /= 0) then
       write(iulog,*)'mosarart_control_type allocation error'
       call shr_sys_abort
    end if

    this%runoff(:,:)      = 0._r8
    this%runofflnd(:,:)   = spval
    this%runoffocn(:,:)   = spval
    this%runofftot(:,:)   = spval
    this%dvolrdt(:,:)     = 0._r8
    this%dvolrdtlnd(:,:)  = spval
    this%dvolrdtocn(:,:)  = spval
    this%volr(:,:)        = 0._r8
    this%flood(:)         = 0._r8
    this%direct(:,:)      = 0._r8
    this%qirrig(:)        = 0._r8
    this%qirrig_actual(:) = 0._r8
    this%qsur(:,:)        = 0._r8
    this%qsub(:,:)        = 0._r8
    this%qgwl(:,:)        = 0._r8
    this%qglc_liq(:)      = 0._r8
    this%qglc_ice(:)      = 0._r8
    this%fthresh(:)       = abs(spval)
    this%flow(:,:)        = 0._r8
    this%erout_prev(:,:)  = 0._r8
    this%eroutup_avg(:,:) = 0._r8
    this%erlat_avg(:,:)   = 0._r8
    this%direct_glc(:,:)  = 0._r8

    this%effvel(:) = effvel0  ! downstream velocity (m/s)
    do nt = 1,ntracers
       do nr = begr,endr
          this%evel(nr,nt) = this%effvel(nt)
       enddo
    enddo

    do nr = begr,endr
       n = this%gindex(nr)
       i = mod(n-1,nlon) + 1
       j = (n-1)/nlon + 1
       this%lonc(nr) = this%rlon(i)
       this%latc(nr) = this%rlat(j)
       this%area(nr) = area_global(n)
    enddo

    larea = 0.0_r8
    do nr = begr,endr
       larea = larea + this%area(nr)
    end do
    if (minval(this%mask) < 1) then
       write(iulog,*) subname,'ERROR this mask lt 1 ',minval(this%mask),maxval(this%mask)
       call shr_sys_abort(subname//' ERROR this mask')
    endif
    call shr_mpi_sum(larea, this%totarea, mpicom_rof, 'mosart totarea', all=.true.)
    if (mainproc) then
       write(iulog,*) subname,'  earth area ',4.0_r8*shr_const_pi*1.0e6_r8*re*re
       write(iulog,*) subname,' mosart area ',this%totarea
    end if

  end subroutine Init

  !========================================================================
  subroutine init_decomp(this, locfn, decomp_option, use_halo_option, &
       nlon, nlat, begr, endr, lnumr, numr, IDkey, rc)

    ! Arguments
    class(control_type)             :: this
    character(len=*)  , intent(in)  :: locfn      ! local routing filename
    character(len=*)  , intent(in)  :: decomp_option
    logical           , intent(in)  :: use_halo_option
    integer           , intent(in)  :: nlon
    integer           , intent(in)  :: nlat
    integer           , intent(out) :: begr
    integer           , intent(out) :: endr
    integer           , intent(out) :: lnumr
    integer           , intent(out) :: numr
    integer           , intent(out) :: IDkey(:)   ! translation key from ID to gindex
    integer           , intent(out) :: rc

    ! Local variables
    integer                    :: n, nr, i, j, g           ! indices
    integer                    :: nl,nloops                ! used for decomp search
    real(r8),allocatable       :: rtempr(:,:)              ! global temporary buffer - real
    integer, allocatable       :: gmask(:)                 ! global mask
    integer, allocatable       :: glo2loc(:)               ! global global->local mapping
    integer, allocatable       :: dnID_global(:)           ! global downstream ID based on ID0
    integer, allocatable       :: idxocn(:)                ! global downstream ocean outlet cell
    integer, allocatable       :: nupstrm(:)               ! number of upstream cells including own cell
    integer, allocatable       :: pocn(:)                  ! pe number assigned to basin
    integer                    :: ID0_global               ! global (local) ID index
    integer                    :: nop(0:npes-1)            ! number of gridcells on a pe
    integer                    :: nba(0:npes-1)            ! number of basins on each pe
    integer                    :: nrs(0:npes-1)            ! begr on each pe
    integer                    :: maxgcells_per_pe         ! max num of points per pe for decomp
    integer                    :: minbas,maxbas            ! used for decomp search
    integer                    :: pid,np,npmin,npmax,npint ! log loop control
    integer                    :: nmos                     ! number of mosart points
    integer                    :: nout                     ! number of basin with outlets
    integer                    :: nbas                     ! number of basin/ocean points
    integer                    :: nrof                     ! num of active mosart points
    integer                    :: baspe                    ! pe with min number of mosart cells
    logical                    :: found                    ! flag
    integer                    :: ier                      ! error status
    type(file_desc_t)          :: ncid                     ! pio file desc
    integer                    :: procid
    integer                    :: im1,ip1
    integer                    :: jm1,jp1
    integer                    :: n_sw, n_s, n_se
    integer                    :: n_nw, n_n, n_ne
    integer                    :: n_e, n_w
    integer                    :: num_halo
    integer, pointer           :: halo_list(:)
    integer, pointer           :: seqlist(:)
    integer, allocatable       :: store_halo_index(:)
    integer                    :: nglob
    character(len=*),parameter :: subname = '(mosart_control_type: init_decomp) '
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------------------------------------------
    ! Read ID and DnID from routing file
    !-------------------------------------------------------

    ! RESET dnID indices based on ID0
    ! rename the dnID values to be consistent with global grid indexing.
    ! where 1 = lower left of grid and nlon*nlat is upper right.
    ! ID0 is the "key", modify dnID based on that.  keep the IDkey around
    ! for as long as needed.  This is a key that translates the ID0 value
    ! to the gindex value.  compute the key, then apply the key to dnID_global.
    ! As part of this, check that each value of ID0 is unique and within
    ! the range of 1 to nlon*nlat.

    call ncd_pio_openfile(ncid, trim(locfn), 0)

    allocate(rtempr(nlon,nlat))
    allocate(dnID_global(nlon*nlat))

    call ncd_io(ncid=ncid, varname='ID', flag='read', data=rtempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read mosart ID')
    if (mainproc) write(iulog,*) 'Read ID ',minval(rtempr),maxval(rtempr)

    IDkey(:) = 0
    do j=1,nlat
       do i=1,nlon
          n = (j-1)*nlon + i
          ID0_global = int(rtempr(i,j))
          if (ID0_global < 0 .or. ID0_global > nlon*nlat) then
             write(iulog,*) subname,' ERROR ID0 out of range',n,ID0_global
             call shr_sys_abort(subname//' ERROR error ID0 out of range')
          endif
          if (IDkey(ID0_global) /= 0) then
             write(iulog,*) subname,' ERROR ID0 value occurs twice',n,ID0_global
             call shr_sys_abort(subname//' ERROR ID0 value occurs twice')
          endif
          IDkey(ID0_global) = n
       end do
    end do
    if (minval(IDkey) < 1) then
       write(iulog,*) subname,' ERROR IDkey incomplete'
       call shr_sys_abort(subname//' ERROR IDkey incomplete')
    endif

    call ncd_io(ncid=ncid, varname='dnID', flag='read', data=rtempr, readvar=found)
    if ( .not. found ) call shr_sys_abort( trim(subname)//' ERROR: read mosart dnID')
    if (mainproc) write(iulog,*) 'Read dnID ',minval(rtempr),maxval(rtempr)
    do j=1,nlat
       do i=1,nlon
          n = (j-1)*nlon + i
          dnID_global(n) = int(rtempr(i,j))
          if (dnID_global(n) > 0 .and. dnID_global(n) <= nlon*nlat) then
             if (IDkey(dnID_global(n)) > 0 .and. IDkey(dnID_global(n)) <= nlon*nlat) then
                dnID_global(n) = IDkey(dnID_global(n))
             else
                write(iulog,*) subname,' ERROR bad IDkey',n,dnID_global(n),IDkey(dnID_global(n))
                call shr_sys_abort(subname//' ERROR bad IDkey')
             endif
          endif
       end do
    end do
    if (mainproc) write(iulog,*) 'dnID ',minval(rtempr),maxval(rtempr)
    deallocate(rtempr)

    call ncd_pio_closefile(ncid)

    !-------------------------------------------------------
    ! Determine mosart ocn/land mask (global, all procs)
    !-------------------------------------------------------

    ! 1=land, 2=ocean, 3=ocean outlet from land
    allocate(gmask(nlon*nlat))
    gmask(:) = 2     ! assume ocean point
    do n=1,nlon*nlat ! mark all downstream points as outlet
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= nlon*nlat)) then
          gmask(nr) = 3          ! <- nr
       end if
    enddo
    do n=1,nlon*nlat         ! now mark all points with downstream points as land
       nr = dnID_global(n)
       if ((nr > 0) .and. (nr <= nlon*nlat)) then
          gmask(n) = 1           ! <- n
       end if
    enddo

    !-------------------------------------------------------
    ! Compute total number of basins and runoff points
    !-------------------------------------------------------

    nbas = 0
    nrof = 0
    nout = 0
    nmos = 0
    do nr=1,nlon*nlat
       if (gmask(nr) == 3) then
          nout = nout + 1
          nbas = nbas + 1
          nmos = nmos + 1
          nrof = nrof + 1
       elseif (gmask(nr) == 2) then
          nbas = nbas + 1
          nrof = nrof + 1
       elseif (gmask(nr) == 1) then
          nmos = nmos + 1
          nrof = nrof + 1
       endif
    enddo
    if (mainproc) then
       write(iulog,*) 'Number of outlet basins = ',nout
       write(iulog,*) 'Number of total  basins = ',nbas
       write(iulog,*) 'Number of mosart points = ',nmos
       write(iulog,*) 'Number of runoff points = ',nrof
    endif

    !-------------------------------------------------------
    ! Compute river basins, actually compute ocean outlet gridcell
    !-------------------------------------------------------

    ! idxocn = final downstream cell, index is global 1d ocean gridcell
    ! nupstrm = number of source gridcells upstream including self
    allocate(idxocn(nlon*nlat))
    allocate(nupstrm(nlon*nlat))
    idxocn(:)  = 0
    nupstrm(:) = 0
    do nr=1,nlon*nlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < nlon*nlat)  ! follow downstream
             nupstrm(n) = nupstrm(n) + 1
             n = dnID_global(n)
             g = g + 1
          end do
          if (gmask(n) == 3) then           ! found ocean outlet
             nupstrm(n) = nupstrm(n) + 1    ! one more land cell for n
             idxocn(nr) = n                 ! set ocean outlet or nr to n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) subname,' ERROR closed basin found', &
                  g,nr,gmask(nr),dnID_global(nr),n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR closed basin found')
          elseif (gmask(n) == 2) then
             write(iulog,*) subname,' ERROR found invalid ocean cell ',nr
             call shr_sys_abort(subname//' ERROR found invalid ocean cell')
          else
             write(iulog,*) subname,' ERROR downstream cell is unknown', &
                  g,nr,gmask(nr),dnID_global(nr),n,gmask(n),dnID_global(n)
             call shr_sys_abort(subname//' ERROR downstream cell is unknown')
          endif
       elseif (gmask(n) >= 2) then  ! ocean, give to self
          nupstrm(n) = nupstrm(n) + 1
          idxocn(nr) = n
       endif
    enddo

    !-------------------------------------------------------
    !--- Now allocate those basins to pes
    !-------------------------------------------------------

    ! this is the heart of the decomp, need to set pocn and nop by the end of this
    ! pocn is the pe that gets the basin associated with ocean outlet nr
    ! nop is a running count of the number of mosart cells/pe
    allocate(pocn(nlon*nlat))
    pocn(:) = -99
    nop(0:npes-1) = 0
    if (trim(decomp_option) == 'basin') then

       baspe = 0
       maxgcells_per_pe = int(float(nrof)/float(npes)*0.445) + 1
       nloops = 3
       minbas = nrof
       do nl=1,nloops
          maxbas = minbas - 1
          minbas = maxval(nupstrm)/(2**nl)
          if (nl == nloops) minbas = min(minbas,1)
          do nr=1,nlon*nlat
             if (gmask(nr) >= 2 .and. nupstrm(nr) > 0 .and. nupstrm(nr) >= minbas .and. nupstrm(nr) <= maxbas) then
                ! Decomp options
                !   find min pe (implemented but scales poorly)
                !   use increasing thresholds (implemented, ok load balance for l2r or calc)
                !   distribute basins using above methods but work from max to min basin size
                ! find next pe below maxgcells_per_pe threshhold and increment
                do while (nop(baspe) > maxgcells_per_pe)
                   baspe = baspe + 1
                   if (baspe > npes-1) then
                      baspe = 0
                      ! 3 loop, .445 and 1.5 chosen carefully
                      maxgcells_per_pe = max(maxgcells_per_pe*1.5, maxgcells_per_pe+1.0)
                   endif
                enddo
                if (baspe > npes-1 .or. baspe < 0) then
                   write(iulog,*) 'ERROR in decomp for mosart ',nr,npes,baspe
                   call shr_sys_abort('ERROR mosart decomp')
                endif
                nop(baspe) = nop(baspe) + nupstrm(nr)
                pocn(nr) = baspe
             endif
          enddo ! nr
       enddo ! nl

       ! set pocn for land cells, was set for ocean above
       do nr=1,nlon*nlat
          if (idxocn(nr) > 0) then
             pocn(nr) = pocn(idxocn(nr))
             if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
                write(iulog,*) subname,' ERROR pocn lnd setting ',&
                     nr,idxocn(nr),idxocn(idxocn(nr)),pocn(idxocn(nr)),pocn(nr),npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    elseif (trim(decomp_option) == '1d') then

       ! distribute active points in 1d fashion to pes
       ! baspe is the pe assignment
       ! maxgcells_per_pe is the maximum number of points to assign to each pe
       baspe = 0
       maxgcells_per_pe = (nrof-1)/npes + 1
       do nr=1,nlon*nlat
          pocn(nr) = baspe
          nop(baspe) = nop(baspe) + 1
          if (nop(baspe) >= maxgcells_per_pe) then
             baspe = (mod(baspe+1,npes))
             if (baspe < 0 .or. baspe > npes-1) then
                write(iulog,*) subname,' ERROR basepe ',baspe,npes
                call shr_sys_abort(subname//' ERROR pocn lnd')
             endif
          endif
       enddo

    elseif (trim(decomp_option) == 'roundrobin') then

       ! distribute active points in roundrobin fashion to pes
       ! baspe is the pe assignment
       ! maxgcells_per_pe is the maximum number of points to assign to each pe
       baspe = 0
       do nr=1,nlon*nlat
          pocn(nr) = baspe
          nop(baspe) = nop(baspe) + 1
          baspe = (mod(baspe+1,npes))
          if (baspe < 0 .or. baspe > npes-1) then
             write(iulog,*) subname,' ERROR basepe ',baspe,npes
             call shr_sys_abort(subname//' ERROR pocn lnd')
          endif
       enddo
       do nr = 1,nlon*nlat
          if (pocn(nr) < 0) then
             write(6,*)'WARNING: nr,pocn(nr) is < 0',nr,pocn(nr)
          end if
       end do

    else
       write(iulog,*) subname,' ERROR decomp option unknown ',trim(decomp_option)
       call shr_sys_abort(subname//' ERROR pocn lnd')
    endif  ! decomp_option

    if (mainproc) then
       write(iulog,*) 'mosart cells and basins total  = ',nrof,nbas
       write(iulog,*) 'mosart cells per basin avg/max = ',nrof/nbas,maxval(nupstrm)
       write(iulog,*) 'mosart cells per pe    min/max = ',minval(nop),maxval(nop)
       write(iulog,*) 'mosart basins per pe   min/max = ',minval(nba),maxval(nba)
    endif
    deallocate(nupstrm)

    !-------------------------------------------------------
    ! Determine begr, endr, numr and lnumr
    !-------------------------------------------------------

    numr = 0
    do n = 0,npes-1
       if (iam == n) then
          begr  = numr + 1
          endr  = begr + nop(n)  - 1
       endif
       numr = numr  + nop(n)
    enddo
    lnumr = endr - begr + 1

    !-------------------------------------------------------
    ! Determine glo2loc (global to local)
    !-------------------------------------------------------

    ! pocn(nlon*nlat) pe number assigned to basin
    ! nop(0:npes-1)   number of gridcells on a pe
    ! nba(0:npes-1)   number of basins on each pe
    ! nrs(0:npes-1)   begr on each pe

    ! Determine glo2loc
    ! nrs is begr on each pe
    ! reuse nba for nop-like counter here, pocn -99 is unused cell

    nrs(:) = 0
    nrs(0) = 1
    do n = 1,npes-1
       ! nop is number of cells per pe
       ! so loop through the pes and determine begr on each pe
       nrs(n) = nrs(n-1) + nop(n-1)
    enddo

    allocate(glo2loc(nlon*nlat))
    glo2loc(:) = 0
    nba(:) = 0
    do nr = 1,nlon*nlat
       procid = pocn(nr)
       if (procid >= 0) then
          glo2loc(nr) = nrs(procid) + nba(procid)
          nba(procid) = nba(procid) + 1
       endif
    enddo
    do n = 0,npes-1
       if (nba(n) /= nop(n)) then
          write(iulog,*) subname,' ERROR mosart cell count ',n,nba(n),nop(n)
          call shr_sys_abort(subname//' ERROR mosart cell count')
       endif
    enddo

    ! Determine gindex
    allocate(this%gindex(begr:endr))
    do j = 1,nlat
       do i = 1,nlon
          n = (j-1)*nlon + i
          if (dnID_global(n) > 0) then
             if (glo2loc(dnID_global(n)) == 0) then
                write(iulog,*) subname,' ERROR glo2loc dnID_global ',&
                     nr,n,dnID_global(n),glo2loc(dnID_global(n))
                call shr_sys_abort(subname//' ERROT glo2loc dnID_global')
             end if
          end if
          nr = glo2loc(n)
          if (nr >= begr .and. nr <= endr) then
             this%gindex(nr) = n
          endif
       end do
    end do

    !-------------------------------------------------------
    ! Create distGrid from global index array
    !-------------------------------------------------------

    allocate(seqlist(endr-begr+1))
    n = 0
    do nr = begr,endr
       n = n + 1
       seqlist(n) = this%gindex(nr)
    end do
    this%DistGrid = ESMF_DistGridCreate(arbSeqIndexList=seqlist, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    deallocate(seqlist)

    !-------------------------------------------------------
    ! Determine local lonc and latc
    !-------------------------------------------------------

    allocate(this%lonc(begr:endr), this%latc(begr:endr))
    do nr = begr,endr
       n = this%gindex(nr)
       i = mod(n-1,nlon) + 1
       j = (n-1)/nlon + 1
       this%lonc(nr) = this%rlon(i)
       this%latc(nr) = this%rlat(j)
    end do

    !-------------------------------------------------------
    ! Determine halo points and create halo route handle
    !-------------------------------------------------------
    if( use_halo_option ) then
       ! note that for each gridcell below there are nhalo extra elements that need to be allocated
       ! Need to keep track of the global index of each halo point
       ! temporary allocatable array store_halo_index = size((endr-begr+1)*nhalo) (nhalo is the number of halo points)
       !
       ! Allocate halo_arrayptr_index - local index (starting at 1) into this%halo_arrayptr on my pe
       allocate(this%halo_arrayptr_index(endr-begr+1,max_num_halo))
       this%halo_arrayptr_index(:,:) = -999

       allocate(store_halo_index((endr-begr+1)*max_num_halo))
       store_halo_index(:) = 0

       do nr = begr,endr
          n = this%gindex(nr)
          i = mod(n-1,nlon) + 1
          j = (n-1)/nlon + 1
          jm1 = j-1
          jp1 = j+1
          im1 = i-1
          ip1 = i+1
          if (i == 1) im1 = 1
          if (j == 1) jm1 = 1
          if (i == nlon) ip1 = nlon
          if (j == nlat) jp1 = nlat
          n_sw = (jm1-1)*nlon + im1
          n_s  = (jm1-1)*nlon + i
          n_se = (jm1-1)*nlon + ip1
          n_e  = (  j-1)*nlon + ip1
          n_ne = (jp1-1)*nlon + ip1
          n_n  = (jp1-1)*nlon + i
          n_nw = (jp1-1)*nlon + im1
          n_w  = (  j-1)*nlon + im1
          call set_halo_index(n_sw, halo_sw, glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_s , halo_s , glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_se, halo_se, glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_e , halo_e , glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_ne, halo_ne, glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_n , halo_n , glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_nw, halo_nw, glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
          call set_halo_index(n_w , halo_w , glo2loc, nr, begr, endr, pocn, store_halo_index, this%halo_arrayptr_index)
       end do

       ! Allocate halo_list - global indices of the halo points on my pe
       num_halo = count(store_halo_index /= 0)
       allocate(halo_list(num_halo))
       halo_list(1:num_halo) = store_halo_index(1:num_halo)

       ! Create halo route handle using predefined allocatable memory
       allocate(this%fld_halo_arrayptr(endr-begr+1+num_halo))
       this%fld_halo_arrayptr(:) = 0.
       this%fld_halo_array = ESMF_ArrayCreate(this%distgrid, this%fld_halo_arrayptr, haloSeqIndexList=halo_list, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create a halo route handle - only need one
       call ESMF_ArrayHaloStore(this%fld_halo_array, routehandle=this%haloHandle, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Create ESMF arrays for lon, lat and fld
       allocate(this%lon_halo_arrayptr(endr-begr+1+num_halo))
       this%lon_halo_arrayptr(:) = 0.
       this%lon_halo_array = ESMF_ArrayCreate(this%distgrid, this%lon_halo_arrayptr, haloSeqIndexList=halo_list, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       allocate(this%lat_halo_arrayptr(endr-begr+1+num_halo))
       this%lat_halo_arrayptr(:) = 0.
       this%lat_halo_array = ESMF_ArrayCreate(this%distgrid, this%lat_halo_arrayptr, haloSeqIndexList=halo_list, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Set halo array for lon and lat - these do not change with time
       n = 0
       do nr = this%begr,this%endr
          n = n + 1
          this%lon_halo_arrayptr(n) = this%lonc(nr)
          this%lat_halo_arrayptr(n) = this%latc(nr)
       end do
       call ESMF_ArrayHalo(this%lon_halo_array, routehandle=this%haloHandle, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ArrayHalo(this%lat_halo_array, routehandle=this%haloHandle, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Deallocate memory
       deallocate(halo_list)
       deallocate(store_halo_index)

       ! Now do a test of the halo operation
       call this%test_halo(rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    endif
    deallocate(glo2loc)
    deallocate(pocn)

    !-------------------------------------------------------
    ! Determine mask, outletg and dsig
    !-------------------------------------------------------

    allocate(this%mask(begr:endr), this%outletg(begr:endr), this%dsig(begr:endr))
    do nr = begr,endr
       n = this%gindex(nr)
       this%mask(nr) = gmask(n)
       this%outletg(nr) = idxocn(n)
       if (dnID_global(n) <= 0) then
          this%dsig(nr) = 0
       else
          this%dsig(nr) = dnID_global(n)
       endif
    end do
    deallocate(gmask)
    deallocate(dnID_global)
    deallocate(idxocn)

    !-------------------------------------------------------
    ! Write per-processor runoff bounds depending on dbug level
    !-------------------------------------------------------

    if (mainproc) then
       write(iulog,*) 'total runoff cells numr  = ',numr
    endif
    call mpi_barrier(mpicom_rof,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    elseif (dbug == 3) then
       npint = 1
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)
       if (iam == pid) then
          write(iulog,'(2a,i9,a,i9,a,i9,a,i9)')' mosart decomp info',&
               ' proc = ',iam,' begr = ',begr,' endr = ',endr,' numr = ',lnumr
       endif
       call mpi_barrier(mpicom_rof,ier)
    enddo

  end subroutine init_decomp

  !========================================================================

  subroutine set_halo_index(global_index, halo_index, glo2loc, nr, begr, endr, pocn, store_halo_index, halo_arrayptr_index)

    ! Arguments
    integer, intent(in)    :: global_index
    integer, intent(in)    :: halo_index
    integer, intent(in)    :: glo2loc(:)
    integer, intent(in)    :: nr
    integer, intent(in)    :: begr, endr
    integer, intent(in)    :: pocn(:)
    integer, intent(inout) :: store_halo_index(:)
    integer, intent(inout) :: halo_arrayptr_index(:,:)

    ! Local variables
    integer :: n
    logical :: found_index
    integer :: nsize
    integer :: num_halo
    !-----------------------------------------------------------------------

    nsize = endr-begr+1
    if (pocn(global_index) /= iam) then
       found_index = .false.
       do n = 1,size(store_halo_index)
          if (store_halo_index(n) == global_index) then
             num_halo = n
             found_index = .true.
             exit
          else if (store_halo_index(n) == 0) then
             store_halo_index(n) = global_index
             num_halo = n
             found_index = .true.
             exit
          end if
       end do
       if (.not. found_index) then
          call shr_sys_abort('ERROR: global halo index not found')
       end if
       halo_arrayptr_index(nr-begr+1,halo_index) = nsize + num_halo
    else
       halo_arrayptr_index(nr-begr+1,halo_index) = glo2loc(global_index) - begr + 1
    end if

  end subroutine set_halo_index

  !========================================================================
  subroutine test_halo(this, rc)

    ! Arguments
    class(control_type)  :: this
    integer, intent(out) :: rc

    ! Local variables
    integer  :: i,j
    integer  :: n, nr
    integer  :: nglob
    integer  :: halo_value
    integer  :: valid_value
    real(r8) :: lon, lon_p1, lon_m1
    real(r8) :: lat, lat_p1, lat_m1
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    n = 0
    do nr = this%begr,this%endr
       n = n + 1
       this%fld_halo_arrayptr(n) = this%latc(nr)*10. + this%lonc(nr)/100.
    end do

    call ESMF_ArrayHalo(this%fld_halo_array, routehandle=this%haloHandle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    n = 0
    do nr = this%begr,this%endr
       n = n+1
       nglob = this%gindex(nr)
       i = mod(nglob-1,this%nlon) + 1
       j = (nglob-1)/this%nlon + 1
       if (j== 1) then
          lat_m1 = this%rlat(1)
       else
          lat_m1 = this%rlat(j-1)
       end if
       if (j == this%nlat) then
          lat_p1 = this%rlat(this%nlat)
       else
          lat_p1 = this%rlat(j+1)
       end if
       lat = this%rlat(j)
       if (i == 1) then
          lon_m1 = this%rlon(1)
       else
          lon_m1 = this%rlon(i-1)
       end if
       if (i == this%nlon) then
          lon_p1 = this%rlon(this%nlon)
       else
          lon_p1 = this%rlon(i+1)
       end if
       lon = this%rlon(i)
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_sw))
       valid_value = lat_m1*10 + lon_m1/100.
       if (halo_value /= valid_value) then
          write(6,'(a,2f20.10)')'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_s))
       valid_value = lat_m1*10 + lon/100.
       if (halo_value /= valid_value) then
          write(6,'(a,2f20.10)')'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_se))
       valid_value = lat_m1*10 + lon_p1/100.
       if (halo_value /= valid_value) then
          write(6,'(a,2f20.10)')'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_e))
       valid_value = lat*10 + lon_p1/100.
       if (halo_value /= valid_value) then
          write(6,'(a,2f20.10)')'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_ne))
       valid_value = lat_p1*10 + lon_p1/100.
       if (halo_value /= valid_value) then
          write(6,'(a,2f20.10)')'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
       !
       halo_value = this%fld_halo_arrayptr(this%halo_arrayptr_index(n,halo_nw))
       valid_value = lat_p1*10 + lon_m1/100.
       if (halo_value /= valid_value) then
          write(6,*)'ERROR: halo, valid not the same = ',halo_value, valid_value
          call shr_sys_abort('ERROR: invalid halo')
       end if
    end do

  end subroutine test_halo

  !========================================================================

  subroutine calc_gradient(this, fld, fld_halo_array, dfld_dx, dfld_dy, rc)

    ! Calculate gradient from nine gridcells (center and surrounding)

    ! Uses

    ! Arguments:
    class(control_type)   :: this
    real(r8), intent(in)  :: fld(this%begr:this%endr)
    type(ESMF_Array)      :: fld_halo_array
    real(r8), intent(out) :: dfld_dx(:)      ! gradient x component
    real(r8), intent(out) :: dfld_dy(:)      ! gradient y component
    integer , intent(out) :: rc

    ! Local variables
    integer  :: i, n, nr               ! local indices
    real(r8) :: deg2rad
    real(r8) :: mean_dx, mean_dy, dlon, dlat
    real(r8) :: ax_indices(4)                 ! x indices to add
    real(r8) :: sx_indices(4)                 ! x indices to subtract
    real(r8) :: ay_indices(4)                 ! y indices to add
    real(r8) :: sy_indices(4)                 ! y indices to subtract
    real(r8) :: fld_surrounding(max_num_halo)
    real(r8) :: dx(max_num_halo)
    real(r8) :: dy(max_num_halo)
    integer  :: index
    real(r8), pointer :: fld_halo_arrayptr(:)
    !-----------------------------------------------------------------------

    call t_startf('gradient')

    rc = ESMF_SUCCESS

    ! Define indices for addition/subtraction
    ax_indices(:) = (/halo_ne,halo_e,halo_e,halo_se/) ! x indices to add
    sx_indices(:) = (/halo_nw,halo_w,halo_w,halo_sw/) ! x indices to subtract
    ay_indices(:) = (/halo_ne,halo_n,halo_n,halo_nw/) ! y indices to add
    sy_indices(:) = (/halo_se,halo_s,halo_s,halo_sw/) ! y indices to subtract

    ! degrees to radians
    deg2rad = shr_const_pi / 180._r8

    ! Get pointer to data in ESMF array
    call ESMF_ArrayGet(fld_halo_array, farrayPtr=fld_halo_arrayptr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! update halo array for fld
    n = 0
    do nr = this%begr,this%endr
       n = n + 1
       fld_halo_arrayptr(n) = fld(nr)
    end do
    call ESMF_ArrayHalo(fld_halo_array, routehandle=this%haloHandle, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialize gradient components
    dfld_dx(:) = 0._r8
    dfld_dy(:) = 0._r8

    n = 0
    do nr = this%begr,this%endr
       n = n+1

       ! extract neighbors from halo array
       do i = 1,max_num_halo
          index = this%halo_arrayptr_index(n,i)
          fld_surrounding(i) = fld_halo_arrayptr(index)
          dlon = (this%lon_halo_arrayptr(n) - this%lon_halo_arrayptr(index))
          dlat = (this%lat_halo_arrayptr(n) - this%lat_halo_arrayptr(index))
          dx(i) = shr_const_rearth * abs(deg2rad*dlon) * cos(deg2rad*this%latc(nr))
          dy(i) = shr_const_rearth * abs(deg2rad*dlat)
       enddo

       ! calculate mean spacing
       mean_dx = 0.5_r8 * (dx(halo_w)+dx(halo_e)) ! average dx west and east
       mean_dy = 0.5_r8 * (dy(halo_s)+dy(halo_n)) ! average dy south and north

       ! compute gradient values
       ! for x gradient sum [NE,2xE,SE,-NW,-2xW,-SW]
       ! for y gradient sum [NE,2xN,NW,-SE,-2xS,-SW]
       do i = 1,4
          dfld_dx(n) = dfld_dx(n) + (fld_surrounding(ax_indices(i)) - fld_surrounding(sx_indices(i)))
          dfld_dy(n) = dfld_dy(n) + (fld_surrounding(ay_indices(i)) - fld_surrounding(sy_indices(i)))
       enddo

       dfld_dx(n) = dfld_dx(n) / (8._r8*mean_dx)
       dfld_dy(n) = dfld_dy(n) / (8._r8*mean_dy)

    enddo ! end of nr loop

    call t_stopf('gradient')

  end subroutine calc_gradient

end module mosart_control_type
