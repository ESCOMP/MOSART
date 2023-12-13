module RtmFileUtils

   !-----------------------------------------------------------------------
   ! Module containing file I/O utilities
   !
   ! !USES:
   use shr_sys_mod , only : shr_sys_abort
   use RtmSpmd     , only : mainproc
   use RtmVar      , only : iulog
   !
   ! !PUBLIC TYPES:
   implicit none
   private
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   public :: get_filename  !Returns filename given full pathname
   public :: getfil        !Obtain local copy of file
   !
   !-----------------------------------------------------------------------

contains

   !-----------------------------------------------------------------------
   character(len=256) function get_filename (fulpath)

      ! !DESCRIPTION:
      ! Returns filename given full pathname
      !
      ! !ARGUMENTS:
      character(len=*), intent(in)  :: fulpath !full pathname
      !
      ! !LOCAL VARIABLES:
      integer i     !loop index
      integer klen  !length of fulpath character string
      !----------------------------------------------------------

      klen = len_trim(fulpath)
      do i = klen, 1, -1
         if (fulpath(i:i) == '/') go to 10
      end do
      i = 0
10    get_filename = fulpath(i+1:klen)

   end function get_filename

   !------------------------------------------------------------------------

   subroutine getfil (fulpath, locfn, iflag)

      ! !DESCRIPTION:
      ! Obtain local copy of file. First check current working directory,
      ! Next check full pathname[fulpath] on disk
      !
      ! !ARGUMENTS:
      implicit none
      character(len=*), intent(in)  :: fulpath !Archival or permanent disk full pathname
      character(len=*), intent(out) :: locfn   !output local file name
      integer,          intent(in)  :: iflag   !0=>abort if file not found 1=>do not abort

      ! !LOCAL VARIABLES:
      integer i               !loop index
      logical lexist          !true if local file exists
      !--------------------------------------------------

      ! get local file name from full name
      locfn = get_filename( fulpath )
      if (len_trim(locfn) == 0) then
         if (mainproc) write(iulog,*)'(GETFIL): local filename has zero length'
         call shr_sys_abort()
      else
         if (mainproc) write(iulog,*)'(GETFIL): attempting to find local file ',trim(locfn)
      endif

      ! first check if file is in current working directory.
      inquire (file=locfn,exist=lexist)
      if (lexist) then
         if (mainproc) write(iulog,*) '(GETFIL): using ',trim(locfn),' in current working directory'
         RETURN
      endif

      ! second check for full pathname on disk
      locfn = fulpath

      inquire (file=fulpath,exist=lexist)
      if (lexist) then
         if (mainproc) write(iulog,*) '(GETFIL): using ',trim(fulpath)
         RETURN
      else
         if (mainproc) write(iulog,*)'(GETFIL): failed getting file from full path: ', fulpath
         if (iflag==0) then
            call shr_sys_abort ('GETFIL: FAILED to get '//trim(fulpath))
         else
            RETURN
         endif
      endif

   end subroutine getfil

end module RtmFileUtils
