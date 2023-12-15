module RtmSpmd

   ! SPMD initialization

   implicit none
   private

   ! Default settings valid even if there is no mpi

   logical, public :: mainproc              ! proc 0 logical for printing msgs
   integer, public :: iam                   ! processor number
   integer, public :: npes                  ! number of processors for rtm
   integer, public :: mpicom_rof            ! communicator group for rtm
   integer, public :: ROFID                 ! component id needed for PIO

   ! Public methods
   public :: RtmSpmdInit                ! Initialization

contains

   !-----------------------------------------------------------------------

   subroutine RtmSpmdInit(mpicom)

      !-----------------------------------------------------------------------
      ! MPI initialization (number of processes, etc)
      !
      ! Arguments
      integer, intent(in) :: mpicom
      !
      ! Local variables
      integer :: ier  ! return error status
      integer :: maintask
      !-----------------------------------------------------------------------

      ! Initialize mpi communicator group
      mpicom_rof = mpicom

      ! Get my processor id
      call mpi_comm_rank(mpicom_rof, iam, ier)
      maintask = 0
      if (iam == maintask) then
         mainproc = .true.
      else
         mainproc = .false.
      end if

      ! Get number of processors
      call mpi_comm_size(mpicom_rof, npes, ier)

   end subroutine RtmSpmdInit

end module RtmSpmd
