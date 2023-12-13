module RtmSpmd

   ! SPMD initialization

   implicit none
   private

#include <mpif.h>

   ! Default settings valid even if there is no mpi

   logical, public :: mainproc              ! proc 0 logical for printing msgs
   integer, public :: iam                   ! processor number
   integer, public :: npes                  ! number of processors for rtm
   integer, public :: mpicom_rof            ! communicator group for rtm
   integer, public :: ROFID                 ! component id needed for PIO
   integer, public, parameter :: MAINTASK=0 ! the value of iam which is assigned
                                            ! the mainproc duties

   ! Public methods
   public :: RtmSpmdInit                ! Initialization

   ! Values from mpif.h that can be used
   public :: MPI_INTEGER
   public :: MPI_REAL8
   public :: MPI_LOGICAL
   public :: MPI_CHARACTER
   public :: MPI_SUM
   public :: MPI_MIN
   public :: MPI_MAX

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
      !-----------------------------------------------------------------------

      ! Initialize mpi communicator group
      mpicom_rof = mpicom

      ! Get my processor id
      call mpi_comm_rank(mpicom_rof, iam, ier)
      if (iam == MAINTASK) then
         mainproc = .true.
      else
         mainproc = .false.
      end if

      ! Get number of processors
      call mpi_comm_size(mpicom_rof, npes, ier)

   end subroutine RtmSpmdInit

end module RtmSpmd
