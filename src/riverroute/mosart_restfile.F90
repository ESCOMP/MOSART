module mosart_restfile

   ! Read from and write to the MOSART restart file.

   use shr_kind_mod,       only : r8 => shr_kind_r8, CL => shr_kind_cl, CS => shr_kind_cs
   use shr_sys_mod,        only : shr_sys_abort
   use mosart_vars,        only : iulog, inst_suffix, caseid, nsrest, &
                                  spval, mainproc, nsrContinue, nsrBranch, nsrStartup, &
                                  ctitle, version, username, hostname, conventions, source
   use mosart_data,        only : ctl, Trunoff
   use mosart_histfile,    only : mosart_hist_restart
   use mosart_fileutils,   only : getfil
   use mosart_timemanager, only : timemgr_restart, get_nstep, get_curr_date
   use mosart_io,          only : ncd_pio_createfile, ncd_enddef, ncd_pio_openfile, ncd_pio_closefile, &
                                  ncd_defdim, ncd_putatt, ncd_defvar, ncd_io, ncd_global, ncd_double, &
                                  ncd_getdatetime
   use pio,                only : file_desc_t

   implicit none
   private

   ! public member functions:
   public :: mosart_rest_FileName
   public :: mosart_rest_FileRead
   public :: mosart_rest_FileWrite
   public :: mosart_rest_Getfile
   public :: mosart_rest_TimeManager
   public :: mosart_rest_restart
   !
   ! private member functions:
   private :: restFile_read_pfile
   private :: restFile_write_pfile    ! Writes restart pointer file
   private :: restFile_dimset

   ! true => allow case name to remain the same for branch run
   ! by default this is not allowed
   logical, public :: brnch_retain_casename = .false.

   ! file name for local restart pointer file
   character(len=CL) :: rpntfil = 'rpointer.rof'

   ! initial conditions file name
   character(len=CL), public :: finidat

   ! restart data file name for branch run
   character(len=CL), public :: nrevsn

!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------

   subroutine mosart_rest_FileWrite( file, rdate )

      !-------------------------------------
      ! Read/write MOSART restart file.

      ! Arguments:
      character(len=*) , intent(in) :: file            ! output netcdf restart file
      character(len=*) , intent(in) :: rdate           ! restart file time stamp for name

      ! Local variables
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i       ! index
      logical :: ptrfile ! write out the restart pointer file
      !-------------------------------------

      ! Define dimensions and variables

      if (mainproc) then
         write(iulog,*)
         write(iulog,*)'restFile_open: writing MOSART restart dataset '
         write(iulog,*)
      end if
      call ncd_pio_createfile(ncid, trim(file))
      call restFile_dimset( ncid )
      call mosart_rest_restart ( ncid, flag='define' )
      call mosart_hist_restart ( ncid, flag='define', rdate=rdate )
      call timemgr_restart( ncid, flag='define' )
      call ncd_enddef(ncid)

      ! Write restart file variables
      call mosart_rest_restart( ncid, flag='write' )
      call mosart_hist_restart ( ncid, flag='write' )
      call timemgr_restart( ncid, flag='write' )
      call ncd_pio_closefile(ncid)

      if (mainproc) then
         write(iulog,*) 'Successfully wrote local restart file ',trim(file)
         write(iulog,'(72a1)') ("-",i=1,60)
         write(iulog,*)
      end if

      ! Write restart pointer file
      call restFile_write_pfile( file )

      ! Write out diagnostic info

      if (mainproc) then
         write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
         write(iulog,'(72a1)') ("-",i=1,60)
      end if

   end subroutine mosart_rest_FileWrite

   !-----------------------------------------------------------------------

   subroutine mosart_rest_FileRead( file )

      !-------------------------------------
      ! Read a MOSART restart file.
      !
      ! Arguments
      character(len=*), intent(in) :: file  ! output netcdf restart file
      !
      ! Local variables
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i              ! index
      !-------------------------------------

      ! Read file
      if (mainproc) write(iulog,*) 'Reading restart dataset'
      call ncd_pio_openfile (ncid, trim(file), 0)
      call mosart_rest_restart(ncid, flag='read')
      call mosart_hist_restart(ncid, flag='read')
      call ncd_pio_closefile(ncid)

      ! Write out diagnostic info
      if (mainproc) then
         write(iulog,'(72a1)') ("-",i=1,60)
         write(iulog,*) 'Successfully read restart data for restart run'
         write(iulog,*)
      end if

   end subroutine mosart_rest_FileRead

   !-----------------------------------------------------------------------

   subroutine mosart_rest_TimeManager( file )

      !-------------------------------------
      ! Read a MOSART restart file.
      !
      ! Arguments
      character(len=*), intent(in) :: file  ! output netcdf restart file
      !
      ! Local Variables:
      type(file_desc_t) :: ncid ! netcdf id
      integer :: i              ! index
      !-------------------------------------

      ! Read file
      if (mainproc) write(iulog,*) 'Reading restart Timemanger'
      call ncd_pio_openfile (ncid, trim(file), 0)
      call timemgr_restart(ncid, flag='read')
      call ncd_pio_closefile(ncid)

      ! Write out diagnostic info
      if (mainproc) then
         write(iulog,'(72a1)') ("-",i=1,60)
         write(iulog,*) 'Successfully read restart data for restart run'
         write(iulog,*)
      end if

   end subroutine mosart_rest_TimeManager

   !-----------------------------------------------------------------------

   subroutine mosart_rest_Getfile( file )

      !-------------------------------------
      ! Determine and obtain netcdf restart file

      ! Arguments:
      character(len=*), intent(out) :: file  ! name of netcdf restart file

      ! Local variables:
      integer :: status                 ! return status
      integer :: length                 ! temporary
      character(len=CL) :: ftest,ctest ! temporaries
      character(len=CL) :: path        ! full pathname of netcdf restart file
      !-------------------------------------

      ! Continue run:
      ! Restart file pathname is read restart pointer file
      if (nsrest==nsrContinue) then
         call restFile_read_pfile( path )
         call getfil( path, file, 0 )
      end if

      ! Branch run:
      ! Restart file pathname is obtained from namelist "nrevsn"
      if (nsrest==nsrBranch) then
         length = len_trim(nrevsn)
         if (nrevsn(length-2:length) == '.nc') then
            path = trim(nrevsn)
         else
            path = trim(nrevsn) // '.nc'
         end if
         call getfil( path, file, 0 )

         ! Check case name consistency (case name must be different
         ! for branch run, unless brnch_retain_casename is set)
         ctest = 'xx.'//trim(caseid)//'.mosart'
         ftest = 'xx.'//trim(file)
         status = index(trim(ftest),trim(ctest))
         if (status /= 0 .and. .not.(brnch_retain_casename)) then
            write(iulog,*) 'Must change case name on branch run if ',&
                 'brnch_retain_casename namelist is not set'
            write(iulog,*) 'previous case filename= ',trim(file),&
                 ' current case = ',trim(caseid), ' ctest = ',trim(ctest), &
                 ' ftest = ',trim(ftest)
            call shr_sys_abort()
         end if
      end if

      ! Initial run
      if (nsrest==nsrStartup) then
         call getfil( finidat, file, 0 )
      end if

   end subroutine mosart_rest_Getfile

   !-----------------------------------------------------------------------

   subroutine restFile_read_pfile( pnamer )

      !-------------------------------------
      ! Setup restart file and perform necessary consistency checks

      ! Arguments
      character(len=*), intent(out) :: pnamer ! full path of restart file

      ! Local variables
      integer :: nio             ! restart unit
      integer :: ier             ! error return from fortran open
      integer :: i               ! index
      character(len=CL) :: locfn ! Restart pointer file name
      !-------------------------------------

      ! Obtain the restart file from the restart pointer file.
      ! For restart runs, the restart pointer file contains the full pathname
      ! of the restart file. For branch runs, the namelist variable
      ! [nrevsn] contains the full pathname of the restart file.
      ! New history files are always created for branch runs.

      if (mainproc) then
         write(iulog,*) 'Reading restart pointer file....'
      endif
      locfn = './'// trim(rpntfil)//trim(inst_suffix)
      open (newunit=nio, file=trim(locfn), status='unknown', form='formatted', iostat=ier)
      if (ier /= 0) then
         write(iulog,'(a,i8)')'(restFile_read_pfile): failed to open file '//trim(locfn)//' ierr=',ier
         call shr_sys_abort()
      end if
      read (nio,'(a256)') pnamer
      close(nio)
      if (mainproc) then
         write(iulog,'(a)') 'Reading restart data.....'
         write(iulog,'(72a1)') ("-",i=1,60)
      end if

   end subroutine restFile_read_pfile

   !-----------------------------------------------------------------------

   subroutine restFile_write_pfile( fnamer )

      !-------------------------------------
      ! Open restart pointer file. Write names of current netcdf restart file.
      !
      ! Arguments
      character(len=*), intent(in) :: fnamer
      !
      ! Local variables
      integer :: nio ! restart pointer file unit number
      integer :: ier ! error return from fortran open
      character(len=CL) :: filename  ! local file name
      !-------------------------------------

      if (mainproc) then
         filename= './'// trim(rpntfil)//trim(inst_suffix)
         open (newunit=nio, file=trim(filename), status='unknown', form='formatted', iostat=ier)
         if (ier /= 0) then
            write(iulog,'(a,i8)')'(restFile_write_pfile): failed to open file '//trim(filename)//' ierr=',ier
            call shr_sys_abort()
         end if
         write(nio,'(a)') fnamer
         close(nio)
         write(iulog,*)'Successfully wrote local restart pointer file'
      end if

   end subroutine restFile_write_pfile

   !-----------------------------------------------------------------------

   character(len=CL) function mosart_rest_FileName( rdate )

      ! Arguments
      character(len=*), intent(in) :: rdate   ! input date for restart file name

      mosart_rest_FileName = "./"//trim(caseid)//".mosart"//trim(inst_suffix)//".r."//trim(rdate)//".nc"
      if (mainproc) then
         write(iulog,*)'writing restart file ',trim(mosart_rest_FileName),' for model date = ',rdate
      end if

   end function mosart_rest_FileName

   !------------------------------------------------------------------------

   subroutine restFile_dimset( ncid )

      !-------------------------------------
      ! Read/Write initial data from/to netCDF instantaneous initial data file

      ! Arguments
      type(file_desc_t), intent(inout) :: ncid

      ! Local Variables:
      integer :: dimid              ! netCDF dimension id
      integer :: ier                ! error status
      character(len= 8) :: curdate  ! current date
      character(len= 8) :: curtime  ! current time
      character(len=CL) :: str
      character(len=*),parameter :: subname='restFile_dimset'
      !-------------------------------------

      ! Define dimensions

      call ncd_defdim(ncid, 'nlon'  , ctl%nlon  , dimid)
      call ncd_defdim(ncid, 'nlat'  , ctl%nlat  , dimid)
      call ncd_defdim(ncid, 'string_length', CS , dimid)

      ! Define global attributes

      call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
      call ncd_getdatetime(curdate, curtime)
      str = 'created on ' // curdate // ' ' // curtime
      call ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str))
      call ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username))
      call ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname))
      call ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version))
      call ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source))
      call ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle))
      call ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(caseid))
      call ncd_putatt(ncid, NCD_GLOBAL, 'title', &
           'MOSART Restart information, required to continue a simulation' )

   end subroutine restFile_dimset

   !-----------------------------------------------------------------------

   subroutine mosart_rest_restart(ncid, flag)

      !-------------------------------------
      ! Read/write MOSART restart data.
      !
      ! Arguments:
      type(file_desc_t), intent(inout) :: ncid ! netcdf id
      character(len=*) , intent(in)    :: flag ! 'read' or 'write'

      ! Local variables
      logical            :: readvar ! determine if variable is on initial file
      integer            :: n,nt,nv ! indices
      integer            :: nvariables
      real(r8) , pointer :: dfld(:) ! temporary array
      character(len=CS)  :: vname,uname
      character(len=CL)  :: lname
      !-------------------------------------

      nvariables = 7
      do nv = 1,nvariables
         do nt = 1,ctl%ntracers

            if (nv == 1) then
               vname = 'VOLR_'//trim(ctl%tracer_names(nt))
               lname = 'water volume in cell (volr)'
               uname = 'm3'
               dfld  => ctl%volr(:,nt)
            elseif (nv == 2) then
               vname = 'RUNOFF_'//trim(ctl%tracer_names(nt))
               lname = 'runoff (runoff)'
               uname = 'm3/s'
               dfld  => ctl%runoff(:,nt)
            elseif (nv == 3) then
               vname = 'DVOLRDT_'//trim(ctl%tracer_names(nt))
               lname = 'water volume change in cell (dvolrdt)'
               uname = 'mm/s'
               dfld  => ctl%dvolrdt(:,nt)
            elseif (nv == 4) then
               vname = 'WH_'//trim(ctl%tracer_names(nt))
               lname = 'surface water storage at hillslopes in cell'
               uname = 'm'
               dfld  => Trunoff%wh(:,nt)
            elseif (nv == 5) then
               vname = 'WT_'//trim(ctl%tracer_names(nt))
               lname = 'water storage in tributary channels in cell'
               uname = 'm3'
               dfld  => Trunoff%wt(:,nt)
            elseif (nv == 6) then
               vname = 'WR_'//trim(ctl%tracer_names(nt))
               lname = 'water storage in main channel in cell'
               uname = 'm3'
               dfld  => Trunoff%wr(:,nt)
            elseif (nv == 7) then
               vname = 'EROUT_'//trim(ctl%tracer_names(nt))
               lname = 'instataneous flow out of main channel in cell'
               uname = 'm3/s'
               dfld  => Trunoff%erout(:,nt)
            else
               write(iulog,*) 'ERROR: illegal nv value a ',nv
               call shr_sys_abort()
            endif

            if (flag == 'define') then
               call ncd_defvar(ncid=ncid, varname=trim(vname), &
                    xtype=ncd_double,  dim1name='nlon', dim2name='nlat', &
                    long_name=trim(lname), units=trim(uname), fill_value=spval)
            else if (flag == 'read' .or. flag == 'write') then
               call ncd_io(varname=trim(vname), data=dfld, dim1name='allrof', &
                    ncid=ncid, flag=flag, readvar=readvar)
               if (flag=='read' .and. .not. readvar) then
                  if (nsrest == nsrContinue) then
                     call shr_sys_abort()
                  else
                     dfld = 0._r8
                  end if
               end if
            end if

         enddo
      enddo

      if (flag == 'read') then
         do n = ctl%begr,ctl%endr
            do nt = 1,ctl%ntracers
               if (abs(ctl%volr(n,nt))      > 1.e30) ctl%volr(n,nt) = 0.
               if (abs(ctl%runoff(n,nt))    > 1.e30) ctl%runoff(n,nt) = 0.
               if (abs(ctl%dvolrdt(n,nt))   > 1.e30) ctl%dvolrdt(n,nt) = 0.
               if (abs(Trunoff%wh(n,nt))    > 1.e30) Trunoff%wh(n,nt) = 0.
               if (abs(Trunoff%wt(n,nt))    > 1.e30) Trunoff%wt(n,nt) = 0.
               if (abs(Trunoff%wr(n,nt))    > 1.e30) Trunoff%wr(n,nt) = 0.
               if (abs(Trunoff%erout(n,nt)) > 1.e30) Trunoff%erout(n,nt) = 0.
            end do
            if (ctl%mask(n) == 1) then
               do nt = 1,ctl%ntracers
                  ctl%runofflnd(n,nt) = ctl%runoff(n,nt)
                  ctl%dvolrdtlnd(n,nt)= ctl%dvolrdt(n,nt)
               end do
            elseif (ctl%mask(n) >= 2) then
               do nt = 1,ctl%ntracers
                  ctl%runoffocn(n,nt) = ctl%runoff(n,nt)
                  ctl%dvolrdtocn(n,nt)= ctl%dvolrdt(n,nt)
               enddo
            endif
         enddo
      endif

   end subroutine mosart_rest_restart

end module mosart_restfile
