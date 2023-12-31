module RtmDateTime

   implicit none
   public

contains

   !-----------------------------------------------------------------------
   subroutine getdatetime (cdate, ctime)
      !
      ! A generic Date and Time routine
      !
      use RtmSpmd, only : mpicom_rof, mainproc
      use mpi
      !
      ! Arguments
      character(len=8), intent(out) :: cdate  !current date
      character(len=8), intent(out) :: ctime  !current time
      !
      ! Local variables
      character(len=8)      :: date   !current date
      character(len=10)     :: time   !current time
      character(len=5)      :: zone   !zone
      integer, dimension(8) :: values !temporary
      integer               :: ier    !MPI error code
      !-----------------------------------------------------------------------

      if (mainproc) then
         call date_and_time (date, time, zone, values)

         cdate(1:2) = date(5:6)
         cdate(3:3) = '/'
         cdate(4:5) = date(7:8)
         cdate(6:6) = '/'
         cdate(7:8) = date(3:4)

         ctime(1:2) = time(1:2)
         ctime(3:3) = ':'
         ctime(4:5) = time(3:4)
         ctime(6:6) = ':'
         ctime(7:8) = time(5:6)
      endif

      call mpi_bcast (cdate,len(cdate),MPI_CHARACTER, 0, mpicom_rof, ier)
      call mpi_bcast (ctime,len(ctime),MPI_CHARACTER, 0, mpicom_rof, ier)

   end subroutine getdatetime

end module RtmDateTime
