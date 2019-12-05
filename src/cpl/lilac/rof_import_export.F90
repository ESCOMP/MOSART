module rof_import_export

  use ESMF
  use shr_kind_mod    , only : r8 => shr_kind_r8, cx=>shr_kind_cx, cxx=>shr_kind_cxx, cs=>shr_kind_cs
  use shr_sys_mod     , only : shr_sys_abort
  use rof_shr_methods , only : chkerr
  use RunoffMod       , only : rtmCTL, TRunoff
  use RtmVar          , only : iulog, nt_rtm, rtm_tracers
  use RtmSpmd         , only : masterproc
  use RtmTimeManager  , only : get_nstep

  implicit none
  private ! except

  public  :: import_fields
  public  :: export_fields

  private :: state_getimport
  private :: state_setexport
  private :: state_getfldptr
  private :: check_for_nans

  integer     ,parameter :: debug = 0 ! internal debug level
  character(*),parameter :: F01 = "('(mosart_import_export) ',a,i5,2x,i8,2x,d21.14)"
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine import_fields( gcomp, rc )

    !---------------------------------------------------------------------------
    ! Obtain the runoff input from the mediator and convert from kg/m2s to m3/s
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State) :: importState
    integer          :: n,nt
    integer          :: begr, endr
    integer          :: nliq, nfrz 
    character(len=*), parameter :: subname='(rof_import_export:import_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get import state
    call ESMF_GridCompGet(gcomp, importState=importState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    ! determine output array and scale by unit convertsion
    ! NOTE: the call to state_getimport will convert from input kg/m2s to m3/s
    
    call state_getimport(importState, 'Flrl_rofsur', begr, endr, rtmCTL%area, output=rtmCTL%qsur(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofsub', begr, endr, rtmCTL%area, output=rtmCTL%qsub(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofgwl', begr, endr, rtmCTL%area, output=rtmCTL%qgwl(:,nliq), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_rofi', begr, endr, rtmCTL%area, output=rtmCTL%qsur(:,nfrz), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(importState, 'Flrl_irrig', begr, endr, rtmCTL%area, output=rtmCTL%qirrig(:), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    rtmCTL%qsub(begr:endr, nfrz) = 0.0_r8
    rtmCTL%qgwl(begr:endr, nfrz) = 0.0_r8

    if (debug > 0 .and. masterproc .and. get_nstep() < 5) then
       do n = begr,endr
          write(iulog,F01)'import: nstep, n, Flrl_rofsur = ',get_nstep(),n,rtmCTL%qsur(n,nliq)
          write(iulog,F01)'import: nstep, n, Flrl_rofsub = ',get_nstep(),n,rtmCTL%qsub(n,nliq)
          write(iulog,F01)'import: nstep, n, Flrl_rofgwl = ',get_nstep(),n,rtmCTL%qgwl(n,nliq)
          write(iulog,F01)'import: nstep, n, Flrl_rofi   = ',get_nstep(),n,rtmCTL%qsur(n,nfrz)
          write(iulog,F01)'import: nstep, n, Flrl_irrig  = ',get_nstep(),n,rtmCTL%qirrig(n)
       end do
    end if

  end subroutine import_fields

  !====================================================================================

  subroutine export_fields (gcomp, rc)

    !---------------------------------------------------------------------------
    ! Send the runoff model export state to the mediator and convert from m3/s to kg/m2s
    !---------------------------------------------------------------------------

    ! uses
    use RtmVar, only : ice_runoff

    ! input/output/variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_State)  :: exportState
    integer           :: n,nt
    integer           :: begr,endr
    integer           :: nliq, nfrz
    real(r8), pointer :: rofl(:)
    real(r8), pointer :: rofi(:)
    real(r8), pointer :: flood(:)
    real(r8), pointer :: volr(:)
    real(r8), pointer :: volrmch(:)
    logical, save     :: first_time = .true.
    character(len=*), parameter :: subname='(rof_import_export:export_fields)'
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Get export state
    call ESMF_GridCompGet(gcomp, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set tracers
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') nliq = nt
       if (trim(rtm_tracers(nt)) == 'ICE') nfrz = nt
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(subname),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    if (first_time) then
       if (masterproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
       endif
       first_time = .false.
    end if

    begr = rtmCTL%begr
    endr = rtmCTL%endr

    allocate(rofl(begr:endr))
    allocate(rofi(begr:endr))
    allocate(flood(begr:endr))
    allocate(volr(begr:endr))
    allocate(volrmch(begr:endr))

    if ( ice_runoff )then
       ! separate liquid and ice runoff
       do n = begr,endr
          rofl(n) =  rtmCTL%direct(n,nliq) / (rtmCTL%area(n)*0.001_r8)
          rofi(n) =  rtmCTL%direct(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             ! liquid and ice runoff are treated separately - this is what goes to the ocean
             rofl(n) = rofl(n) + rtmCTL%runoff(n,nliq) / (rtmCTL%area(n)*0.001_r8)
             rofi(n) = rofi(n) + rtmCTL%runoff(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          end if
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = begr,endr
          rofl(n) = (rtmCTL%direct(n,nfrz) + rtmCTL%direct(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             rofl(n) = rofl(n) + (rtmCTL%runoff(n,nfrz) + rtmCTL%runoff(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          endif
          rofi(n) = 0._r8
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ! scs: is there a reason for the wr+wt rather than volr (wr+wt+wh)?
    ! volr(n) = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / rtmCTL%area(n)

    do n = begr, endr
       flood(n)   = -rtmCTL%flood(n)    / (rtmCTL%area(n)*0.001_r8)
       volr(n)    =  rtmCTL%volr(n,nliq)/ rtmCTL%area(n)
       volrmch(n) =  Trunoff%wr(n,nliq) / rtmCTL%area(n)
    end do

    call state_setexport(exportState, 'Forr_rofl', begr, endr, input=rofl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Forr_rofi', begr, endr, input=rofi, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_flood', begr, endr, input=flood, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volr', begr, endr, input=volr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_setexport(exportState, 'Flrr_volrmch', begr, endr, input=volrmch, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (debug > 0 .and. masterproc .and. get_nstep() <  5) then
       do n = begr,endr
          write(iulog,F01)'export: nstep, n, Flrr_flood   = ',get_nstep(), n, flood(n)
          write(iulog,F01)'export: nstep, n, Flrr_volr    = ',get_nstep(), n, volr(n)
          write(iulog,F01)'export: nstep, n, Flrr_volrmch = ',get_nstep(), n, volrmch(n)
          write(iulog,F01)'export: nstep, n, Forr_rofl    = ',get_nstep() ,n, rofl(n)
          write(iulog,F01)'export: nstep, n, Forr_rofi    = ',get_nstep() ,n, rofi(n)
       end do
    end if

    deallocate(rofl, rofi, flood, volr, volrmch)

  end subroutine export_fields

  !===============================================================================

  subroutine state_getimport(state, fldname, begr, endr, area, output, rc)

    ! ----------------------------------------------
    ! Map import state field to output array
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(in)    :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr 
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: area(begr:endr)
    real(r8)            , intent(out)   :: output(begr:endr)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=*), parameter :: subname='(rof_import_export:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr,  rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! determine output array and scale by unit convertsion
       do g = begr,endr
          output(g) = fldptr(g-begr+1) * area(g)*0.001_r8
       end do

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
          do g = begr,endr
             i = 1 + g - begr
             if (output(g) /= 0._r8) then
                ! write(iulog,F01)'import: nstep, n, '//trim(fldname)//' = ',get_nstep(),g,output(g)
             end if
          end do
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), begr)
    end if

  end subroutine state_getimport

  !===============================================================================

  subroutine state_setexport(state, fldname, begr, endr, input, rc)

    use shr_const_mod, only : fillvalue=>SHR_CONST_SPVAL

    ! ----------------------------------------------
    ! Map input array to export state field 
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State)    , intent(inout) :: state
    character(len=*)    , intent(in)    :: fldname
    integer             , intent(in)    :: begr
    integer             , intent(in)    :: endr
    real(r8)            , intent(in)    :: input(begr:endr)
    integer             , intent(out)   :: rc

    ! local variables
    integer                     :: g, i
    real(R8), pointer           :: fldptr(:)
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=*), parameter :: subname='(rof_import_export:state_setexport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine if field with name fldname exists in state
    call ESMF_StateGet(state, trim(fldname), itemFlag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! if field exists then create output array - else do nothing
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then

       ! get field pointer
       call state_getfldptr(state, trim(fldname), fldptr, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       fldptr(:) = 0._r8

       ! set fldptr values to input array
       do g = begr,endr
          fldptr(g-begr+1) = input(g)
       end do

       ! write debug output if appropriate
       if (masterproc .and. debug > 0 .and. get_nstep() < 5) then
          do g = begr,endr
             i = 1 + g - begr
             if (input(g) /= 0._r8) then
!                write(iulog,F01)'export: nstep, n, '//trim(fldname)//' = ',get_nstep(),i,input(g)
             end if
          end do
       end if

       ! check for nans
       call check_for_nans(fldptr, trim(fldname), begr)
    end if

  end subroutine state_setexport

  !===============================================================================

  subroutine state_getfldptr(State, fldname, fldptr, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    type(ESMF_State),  intent(in)    :: State
    character(len=*),  intent(in)    :: fldname
    real(R8), pointer, intent(out)   :: fldptr(:)
    integer,           intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    type(ESMF_Mesh)             :: lmesh
    integer                     :: nnodes, nelements
    character(len=*), parameter :: subname='(rof_import_export:state_getfldptr)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_FieldGet(lfield, status=status, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (status /= ESMF_FIELDSTATUS_COMPLETE) then
       call ESMF_LogWrite(trim(subname)//": ERROR data not allocated ", ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    else
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_MeshGet(lmesh, numOwnedNodes=nnodes, numOwnedElements=nelements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (nnodes == 0 .and. nelements == 0) then
          call ESMF_LogWrite(trim(subname)//": no local nodes or elements ", ESMF_LOGMSG_INFO)
          rc = ESMF_FAILURE
          return
       end if

       call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif  ! status

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)

  end subroutine state_getfldptr

  !===============================================================================

  subroutine check_for_nans(array, fname, begg)

    ! uses
    use shr_infnan_mod, only : isnan => shr_infnan_isnan

    ! input/output variables
    real(r8), pointer             :: array(:)
    character(len=*) , intent(in) :: fname
    integer          , intent(in) :: begg

    ! local variables
    integer :: i
    !-------------------------------------------------------------------------------

    ! Check if any input from mediator or output to mediator is NaN

    if (any(isnan(array))) then
       write(iulog,*) '# of NaNs = ', count(isnan(array))
       write(iulog,*) 'Which are NaNs = ', isnan(array)
       do i = 1, size(array)
          if (isnan(array(i))) then
             write(iulog,*) "NaN found in field ", trim(fname), ' at gridcell index ',begg+i-1
          end if
       end do
       call shr_sys_abort(' ERROR: One or more of the output from MOSART to the coupler are NaN ' )
    end if
  end subroutine check_for_nans

end module rof_import_export
