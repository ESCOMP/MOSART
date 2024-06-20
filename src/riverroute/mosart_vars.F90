module mosart_vars

   use shr_kind_mod             , only : r8 => shr_kind_r8, CL => SHR_KIND_CL, CS => shr_kind_CS
   use shr_const_mod            , only : SHR_CONST_CDAY,SHR_CONST_REARTH
   use shr_sys_mod              , only : shr_sys_abort
   use ESMF                     , only : ESMF_VM

   implicit none
   public

   ! MPI
   logical :: mainproc                 ! proc 0 logical for printing msgs
   integer :: iam                      ! processor number
   integer :: npes                     ! number of processors for mosart
   integer :: mpicom_rof               ! communicator group for mosart
   logical :: barrier_timers = .false. ! barrier timers
   type(ESMF_VM) :: vm                 ! ESMF VM

   ! Constants
   integer  , parameter :: iundef    = -9999999
   integer  , parameter :: rundef    = -9999999._r8
   real(r8) , parameter :: secspday  = SHR_CONST_CDAY     ! Seconds per day
   integer  , parameter :: isecspday = secspday           ! Integer seconds per day
   real(r8) , parameter :: spval     = 1.e36_r8           ! special value for real data
   integer  , parameter :: ispval    = -9999              ! special value for int data

   real(r8)             :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)

   ! Run startup
   integer  , parameter :: nsrStartup  = 0 ! Startup from initial conditions
   integer  , parameter :: nsrContinue = 1 ! Continue from restart files
   integer  , parameter :: nsrBranch   = 2 ! Branch from restart files
   integer              :: nsrest = iundef ! Type of run

   ! Namelist variables
   character(len=CL)  :: frivinp                ! MOSART input data file name
   logical            :: ice_runoff             ! true => runoff is split into liquid and ice, otherwise just liquid
   character(len=CS)  :: decomp_option          ! decomp option
   character(len=CS)  :: bypass_routing_option  ! bypass routing model method
   character(len=CS)  :: qgwl_runoff_option     ! method for handling qgwl runoff
   integer            :: budget_frq = -24       ! budget check frequency

   ! Metadata variables used in history and restart generation
   character(len=CL)  :: caseid  = ' '          ! case id
   character(len=CL)  :: ctitle  = ' '          ! case title
   character(len=CL)  :: hostname = ' '         ! Hostname of machine running on
   character(len=CL)  :: username = ' '         ! username of user running program
   character(len=CL)  :: version  = " "         ! version of program
   character(len=CL)  :: conventions = "CF-1.0" ! dataset conventions
   character(len=CL)  :: model_doi_url          ! Web address of the Digital Object Identifier (DOI) for this model version
   character(len=CL)  :: source   = "Model for Scale Adaptive River Transport MOSART1.0" ! description of this source

   ! Stdout
   integer :: iulog = 6        ! "stdout" log file unit number, default is 6

   ! Instance control
   integer           :: inst_index
   character(len=CS) :: inst_name
   character(len=CS) :: inst_suffix

end module mosart_vars
