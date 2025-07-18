<hr>
# Tag name:  mosart1.1.10
### Originator(s): slevis
### Date: Jul 14, 2025
### One-line Summary: Fix MOSART h0i file metadata

Fixes ESCOMP/MOSART#120 MOSART h0i files end up with the same metadata as h0a files

Testing: standard testing
  izumi ---- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/121 for more details

<hr>
# Tag name:  mosart1.1.09
### Originator(s): slevis
### Date: Jul 03, 2025
### One-line Summary: Separate instantaneous and non-inst. history files

This is the mosart equivalent of ESCOMP/CTSM#2445.
Also includes the merge of #118, though this DART-related one-line update seems to have been present in master already.

Contributors: Erik Kluzek, Kevin Raeder

Fixes ESCOMP/MOSART#52 Separate instantaneous from non-inst. history tapes
Fixes ESCOMP/MOSART#116 Make st_archive handle output files from future MOSART+DART experiments

Testing: standard testing
  izumi ---- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/117 for more details
See https://github.com/ESCOMP/MOSART/pull/118 for more details
Contributes to https://github.com/ESCOMP/CTSM/pull/2445

<hr>
# Tag name:  mosart1.1.08
### Originator(s): samrabin
### Date: Jan 14, 2025
### One-line Summary: Standardize time metadata

Standardizes a dimension name of output variable time_bounds, as well as attributes for that plus mcdate, mcsec, mdcur, and mscur.

Contributors: Adam Phillips, Erik Kluzek

Fixes ESCOMP/MOSART#53
Contributes to https://github.com/ESCOMP/CTSM/issues/1693

Testing: standard testing (ekluzek)
  izumi ---- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/66 for more details
Contributes to https://github.com/ESCOMP/CTSM/pull/2052

<hr>
# Tag name:  mosart1.1.07
### Originator(s): erik
### Date: Jan 11, 2025
### One-line Summary: Fix for Nag compiler

Fix nag compiler, change ChangeLog to markdown format.
Remove complexity of having curr_date_in in timemgr_init
be optional and make it required. As was done in RTM.

Fixes ESCOMP/MOSART#110
Fixes ESCOMP/MOSART#111

Testing: standard testing (ekluzek)
  izumi ---- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/112 for more details

<hr>
# Tag name:  mosart1.1.06
### Originator(s): jedwards4b
### Date: Dec 24, 2024
### One-line Summary: Add simulation timestamp to rpointer filenames

Add timestamps to rpointer files, initialize curr date from driver and compare 
to what is in the restart files (instead of initializing from restart files).

Testing: standard testing (ekluzek)
  izumi ---- OK (but problems with baseline compare)
  derecho -- OK

<hr>
# Tag name:  mosart1.1.05
### Originator(s): slevis
### Date: Nov 12, 2024
### One-line Summary: Stop running 0th timestep

For consistency with CAM and CTSM.

Issues addressed:
Relates to issue ESCOMP/CTSM#925 and PR ESCOMP/CTSM#2084.

Testing: standard testing
  izumi -- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/67 for more details

<hr>
# Tag name:  mosart1.1.04
### Originator(s): slevis
### Date: Nov 11, 2024
### One-line Summary: time now equals the middle of the time_bounds

For consistency with CAM and CTSM and for more intuitive history output.

Issues addressed:
Relates to issue #52 but does not separate instantaneous fields into separate history tapes.

Testing: standard testing
  izumi -- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/106 for more details

<hr>
# Tag name:  mosart1.1.03
### Originator(s): ekluzek
### Date: Nov 6, 2024
### One-line Summary: Add more tests and merge in some other PRs


Issues addressed:
Resolves #68 add more gnu and izumi tests
Resolves #79 nobypass test fails on derecho
Resolves #61 incorrect compset in config_compsets.xml
Resolves #91 change clm51 tests to clm60
Resolves #103 rof_to_glc issue (we had to revert later)
Resolves #104 reduce logging noise
Resolves #107 change all mosart tests to clm60?

Testing: standard testing
  izumi -- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/70 for more details

<hr>
# Tag name:  mosart1.1.02
### Originator(s): mvertens
### Date: Jun 21, 2024
### One-line Summary: cism runoff will be now routed to ocn via mosart

Enables CISM runoff to be routed to the ocean via mosart.

All runoff from CISM will be routed directly to the outlet points
New fields will be advertised in the mosart cap
call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofl_glc')
call fldlist_add(fldsFrRof_num, fldsFrRof, 'Forr_rofi_glc')

Issues addressed:
Fixes #92
Fixes #102

Testing: standard testing
  izumi -- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/94 for more details

<hr>
# Tag name:  mosart1.1.01
### Originator(s): mvertens
### Date: Jun 06, 2024
### One-line Summary: major mosart refactor including addition of new halo capability

Removed all references to rtm

files have been renamed and namelists no longer contain rtm in the name

New modularity:

introduced new modules with new derived types and methods
mosart_control_type.F90
mosart_tctl_type.F90
mosart_tparameter_type.F90
mosart_tspatialunit_type.F90
mosart_tstatusflux_type.F90

the new modules modularize a lot of the complexity and variables that were previously found in RunOffMod.F90 and permit decomposition initialization to be more flexible and transparent.

New halo capability

Ability to have halo regions and communication using ESMF. This is needed for computing derivatives in upcoming new additions to MOSART.
New halo namelist - use_halo_option. When this is set to true halos can be activated. See the test_halo subroutine in mosart_control_type.F90 module.
Verified that the results for the halos are bfb identical regardless of the number of processors that are used.
-To set the values for the exclusive region that will be used in halo operations - you need to access the pointer as is done in the test_halo routine in mosart_control_type.F90:
     n = 0
      do nr = this%begr,this%endr
         n = n + 1
         this%halo_arrayptr(n) = this%latc(nr)*10. + this%lonc(nr)/100.
      end do

      call ESMF_ArrayHalo(this%haloArray, routehandle=this%haloHandle, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return

Issues addressed:
Fixes #93 "if (ierr /= PIO_NOERR)" will not be invoked unless PIO_BCAST_ERROR is explicitly set
Fixes #98 Change namelist items to remove "rtm" in the names of namelist variables
Fixes #97 Remove RTM in more of the MOSART code (filenames, subroutines, variables etc.)
Fixes #99 Add a new mosart_noresm testlist

Testing: standard testing
  izumi -- OK
  derecho -- OK

See https://github.com/ESCOMP/MOSART/pull/76 for more details

<hr>
# Tag name:  mosart1_0_49
### Originator(s): mvertens
### Date: Feb 02, 2024
### One-line Summary: Remove MCT, some cleanup and high level refactoring

Removes all MCT references from the code and replaces them with ESMF routehandles and mapping calls
major changes to RtmMod.F90 along with other code cleanup described below

RtmVar
Now contains new ESMF data types needed for the MOSART mapping
       type(ESMF_Field)       , public :: srcField
       type(ESMF_Field)       , public :: dstField
       type(ESMF_RouteHandle) , public :: rh_dnstream
       type(ESMF_RouteHandle) , public :: rh_direct
       type(ESMF_RouteHandle) , public :: rh_eroutUp

RtmMod:
now have two new init phases for mosart. The first init phase is now called MOSART_init1 and replaces Rtmini. This has mostly what was there before but moves the creation of all routehandles to the second init phase - MOSART_init2 which must be called after the mesh has been read in. Also - moved the section of code for MOSART_init2 to be right below the section for MOSART_init1.
removed the mapping for Smatp_dnstrm since it was not used and there is no reason to create a map that is not needed. The associated code that was commented out for this has also been removed.
renamed RtmRun to MOSART_run
new indentation
MOSART_physics.F90
now using the computed routehandle rh_eroutUp
new indentation
Removed namelist variable do_rtmflood and xml variable MOSART_FLOOD_MODE. Also removed subroutine MOSART_FloodInit in RtmMod.F90 which was never activated and in fact the model aborted if you tried to invoke it.
Verified that this was no longer needed in consult with @swensosc.
masterproc -> mainproc
updated the MOSART testlist for derecho and betzy (betzy is a NorESM platform) and added a PFS test

Issues resolved:
  Resolves #65 -- Remove MCT
  Resolves #75 -- masterproc to mainproc
  Resolves #73 -- testlist to Derecho
  Resolved #85 -- Remove RtmFileUtils

Testing: standard testing
  izumi -- PASS
  derecho -- PASS  (following change answers but determined to be OK)
ERP_D.f10_f10_mg37.I1850Clm50Bgc.derecho_intel.mosart-qgrwlOpts
PEM_D.f10_f10_mg37.I1850Clm50Sp.derecho_intel.mosart-inplacethreshold
SMS_D.f10_f10_mg37.I1850Clm50Bgc.derecho_intel.mosart-decompOpts

(first two due to baseline not having history output, so rerunning shows b4b)
(Last one shows roundoff level answer changes)

See https://github.com/ESCOMP/MOSART/pull/74 for more details

<hr>
# Tag name:  mosart1_0_48
### Originator(s): erik
### Date: Nov 13, 2022
### One-line Summary: Change bypass_routing option to direct_to_outlet

Change the bypass_routing option to negative flow is sent to
the river outlets.

Fixes #58 -- change bypass_routing default

<hr>
# Tag name:  mosart1_0_47
### Originator(s): erik
### Date: Nov 12, 2022
### One-line Summary: Some fixes for the direct_to_outlet option

This fixes the balance errors evident in the mosart log file when 
the bypass_routing_option='direct_to_outlet' method is used.
Note that these balance errors are warnings and do not stop the model.

Code from @swensosc, PR by @olyson

Fixes #54 -- update cime directory in buildlib/buildnml
Fixes #56 -- Some issues with the direct_to_outlet option

Also don't allow the threshold option to be set with
direct_to_outlet is choosen, as well as any option other
than "all" for "none".

Answers change when direct_to_outlet option is used

<hr>
# Tag name:  mosart1_0_46
### Originator(s): jedwards
### Date: Jun 12, 2022
### One-line Summary: PIO ascynchronus changes

Currently mosart initializes and try's to use IO in the advertise phase, but asyncio requires that this is not done until the
realize phase. This PR splits the mosart initialization so that the namelist is still read in the advertise phase but IO waits until
the realize phase. Tested with cesm prealpha tests on cheyenne, a modification was also required in the mct driver to handle the new
initialize subroutine.

<hr>
# Tag name:  mosart1_0_45
### Originator(s): erik
### Date: Nov 12, 2021
### One-line Summary: Handle CLM_ACCELERATED_SPINUP option differently

Handle CLM_ACCELERATED_SPINUP option differently so it doesn't change 
MOSART_MODE, but dies (unless the new xml variable MOSART_IGNORE_WARNINGS 
is set to TRUE). Also only output information about decomposition from all 
PE's only if running with DEBUG compiler options set.

Fixes #48 -- Remove setting of MOSART_MODE from buildnml
Fixes #49 -- Reduce MOSART output to cesm.log

<hr>
# Tag name:  mosart1_0_44
### Originator(s): erik
### Date: Nov 02, 2021
### One-line Summary: Nuopc now the default driver, pass channel depth fields, remove do_rtm namelist

Assume nuopc is the default driver now, so switch testlist for special driver tests
to be mct rather than nuopc. Pass channel depth fields if flds_r2l_stream_channel_depths
is set in the driver. Remove do_rtm option and trigger MOSART_MODE==NULL instead.

Fixes #47 -- do_rtm not working with NUOPC
Fixes #45 -- frivinp_rtm is still required even if MOSART is turned off
Fixes #44 -- Remove reference to Flrl_rofdto 
Fixes #36 -- pass water depths to LND
Fixes #20 -- MOSART_SIMYR not setup correctly

This test fails in the build unless there is an appropriate verison of CMEPS:

SMS_D.f10_f10_mg37.I1850Clm50Bgc.cheyenne_intel.mosart-passChannelDepths
   SHARED_BUILD -- FAIL

<hr>
# Tag name:  mosart1_0_43
### Originator(s): erik
### Date: July 16, 2021
### One-line Summary: Change shebang lines to python3, abort if direction file not set, add compset

Change the shebang lines in buildnml and buildlib for python3
Abort if a direction file is not set
Add a compset for running MOSART on it's own (not tested yet)

Fixes #42 -- Change shebang lines for python
Fixes #35 -- point to directory
Fixes #28 -- add R2000MOSART compset

See https://github.com/ESCOMP/MOSART/pull/43

<hr>
# Tag name:  mosart1_0_42
### Originator(s): jedwards
### Date: April 5, 2021
### One-line Summary: Support for threading with NUOPC/CMEPS

See https://github.com/ESCOMP/MOSART/pull/41

Testing:
Mariana Vertenstein ran
ERS_Vnuopc_Ld5.f09_g17.I1850Clm50Sp.cheyenne_intel.clm-default; this was
bit-for-bit against baselines.

<hr>
# Tag name:  mosart1_0_41
### Originator(s): jedwards
### Date: March 17, 2021
### One-line Summary: Minor changes for pio2 compatibility

Testing: mosart test suite on cheyenne & izumi, in the context of
ctsm5.1.dev027. All tests passed (once I removed the two non-debug nag
tests, which now fail) and were bit-for-bit with baselines (baselines
were generated using mosart1_0_40 in the context of ctsm5.1.dev027).

<hr>
# Tag name:  mosart1_0_40
### Originator(s): mvertens, jedwards
### Date: March 16, 2021
### One-line Summary: Addition of flux area correction factors

Adds area flux correction factors as is done in CPL7. Remove some debug
writes of coupler fields. And remove state_getfldptr as it was unused.

Also sets the number of OpenMP threads according to the coupler namelist
rather than env variables.

Testing:
SMS_D_Ld5_Vnuopc.f10_f10_mg37.I2000Clm50BgcCrop.cheyenne_intel.mosart-default
in the context of CESM's nuopc_dev branch. (No baseline comparisons done
for full set of changes.)

<hr>
# Tag name:  mosart1_0_39
### Originator(s): jedwards
### Date: December 22, 2020
### One-line Summary: Remove nuopc_cap_share from build

nuopc_cap_share is now included in the share library.

Testing: mosart test suite on cheyenne & izumi, in the context of
ctsm5.1.dev019. All tests passed and were bit-for-bit with baselines.


<hr>
# Tag name:  mosart1_0_38
### Originator(s): sacks
### Date: November 01, 2020
### One-line Summary: Change compset in test

In an upcoming CTSM tag, the compset I2000Clm50BgcCropGs will be renamed
to I2000Clm50BgcCrop (while keeping the same meaning as before). This
tag changes a test to use the new alias, to keep MOSART compatible with
this upcoming CTSM tag.

Note that this change will require the upcoming CTSM tag (probably
ctsm5.1.dev011) in order to run this test from the mosart test list.

<hr>
# Tag name:  mosart1_0_37
### Originator(s): mvertens, erik, jedwards
### Date: August 06, 2020
### One-line Summary: Bring updates needed for NUOPC

Some more changes from @mvertens with updates for the NUOPC cap. 
There's a new nuopc_cap_share directory added. SetVM is made public
in rof_comp_nuopc.F90 and some of the logging is slightly changed. 
RtmIO is now using ROFID rather than instance name.

Also convert the nldas2 file into NetCDF5 format from Jim Edwards.

<hr>
# Tag name:  mosart1_0_36
### Originator(s): mvertens, sacks
### Date: March 16, 2020
### One-line Summary: Updates to NUOPC cap

Main purpose is updates to NUOPC cap. Also a fix for PIO2.

Most changes are from Mariana Vertenstein, brought to master by Bill
Sacks.

These changes are documented in https://github.com/ESCOMP/MOSART/pull/30

Testing: mosart test suite (cheyenne & izumi), in the context of a CTSM
checkout (https://github.com/ESCOMP/CTSM/pull/939). Baseline comparisons
done against ctsm1.0.dev085.

<hr>
# Tag name:  mosart1_0_35
### Originator(s): erik
### Date: Nov 08, 2019
### One-line Summary: Fix cold-start, 8th degree file, wallclock, add COC, change SHR_KIND

    Fix for cold start #24 and add cold start test

    Fix 8th degree file, #26

    Change wallclock time in tests and change hobart for izumi

    Add code of conduct

    Change all instances of SHR_KIND_CL to use CL that is renamed from
    shr_kind_cl in shr_kind_mod. This is important for a recent cime
    update.

<hr>
# Tag name:  mosart1_0_34
### Originator(s): mvertens/erik
### Date: Jul 25, 2019
### One-line Summary: Add nuopc cap for NUOPC coupler option

Add another driver cap for NUOPC in addition to the MCT cap.
As part of this rof_cpl_indices.F90 which was shared under
"src/cpl" was moved inside of the "src/cpl/mct" with "rof"
renamed to "mosart" since this isn't shared in between the
driver caps. The NUOPC cap doesn't use the coupler integer
indexes, only the character names.

A helper module for NUOPC ESMF functionality was added called
rof_comp_nuopc.

This also meant that nt_rtm and rtm_tracers was moved from
the rof_cpl_indices file to RtmVar.F90 so it could be shared
between the two caps.

Also to match the nuopc cap the import_export part of the
MCT cap was pulled out in it's own file.

Testing: Run mosart test suite (with ctsm tag ctsm1.0.dev049 and
         branch_tags/cime5.8.3_chint17-03
   izumi ---- PASS (compare to hobart)

Pull Request: #25
  #25 -- nuopc cap

<hr>
<hr>
# Tag name:  mosart1_0_33
### Originator(s): erik
### Date: Jun 11, 2019
### One-line Summary: buildlib updates for cime5.8 and nldas grid

buildlib changes needed for cime5.8 from @jedwards.

Also adds a nldas routing grid from @billsacks.

Testing: Run mosart test suite (with cime tag branch_tags/cime5.8.3_chint17-02)
   hobart ---- PASS
   cheyenne -- PASS

<hr>
# Tag name:  mosart1_0_32
### Originator(s): erik
### Date: May 07, 2019
### One-line Summary:  Move release-cesm2.0.03 to mosart master

Don't allow the namelist option rtmhist_ndens to be set to 2, because this
option doesn't currently function. The simple fix we put into place for it
is not robust.

Fix for python3, using the floor operator for integer division.

Add in 8th degree routine file. Run pylint and check for py3 compatability.
There's explict setting of fill type. And also explicit use of shr_kind_r4 for 
kind rather than real(4), which is a better mechanism. Also ncd_getiodesc 
will read in PIO_DOUBLE for input xtype= PIO_DOUBLE or PIO_REAL. Most of this
is direct from jedwards4b (other than r8th addition).

Issues Fixed: #18 #10
  #18 -- rtmhist_ndens=2 does NOT work (so don't allow it as an option)
  #10 -- With python3, coupling_period is real rather than int

Science  changes since: mosart1_0_31
   Added in 8th degree routing file (r8th)
Software changes since: mosart1_0_31
   Run pylint on python buildlib and buildnml scripts, check for py3 compatibility.
   Corrects the integer fill value. Needed for pio2.

Testing: Run mosart test suite
   hobart ---- PASS
   cheyenne -- PASS

<hr>
# Tag name:  mosart1_0_31
### Originator(s): erik
### Date: May 15, 2018
### One-line Summary:  Add model_doi_url

Add config_archive for mosart, delete rof_comp_esmf, add model_doi_url, 
change some instances of RTM/CLM in documentation to MOSART

Pull request #11

https://github.com/ESCOMP/mosart/pull/11

<hr>
# Tag name:  mosart1_0_30
### Originator(s): erik
### Date: Jan 24, 2018
### One-line Summary:  Fix testlist so options are under test not machines

Fix testlist so options are under <test> rather than <machines>.
Also remove aux_clm category (any of those tests should be in CLM's testlist).

 M cime_config/testdefs/testlist_mosart.xml

<hr>
# Tag name:  mosart1_0_29
### Originator(s): erik
### Date: Jan 17, 2018
### One-line Summary:  Update testlist to version 2 format, remove ys tests

 Add buildnmlc to .gitignore

Convert the test list to version 2 format, add comments and wallclock to it.
Remove ys tests. Remove aux_clm tests (put those in CLM itself).  Add
two more test mods for currently untested options.

 M cime_config/testdefs/testlist_mosart.xml
 A cime_config/testdefs/testmods_dirs/mosart/decompOpts/user_nl_mosart
 A cime_config/testdefs/testmods_dirs/mosart/qgrwlOpts/user_nl_mosart

<hr>
# Tag name:  mosart1_0_28
### Originator(s): erik
### Date: Oct 05, 2017
### One-line Summary:  Explicitly truncate the river length to a mininum value, rather than
                   just change the calculation of tlen for rlen_min.

Testing: Sean and Keith both ran long simulations with this change (in a branch off of 1_0_26)

Previous simulations showed that short rivers would accumulate water in them slowly over 
time. This doesn't happen in longer rivers. But it affects the total water balance of the
system. Simulations by Sean Swenson and Keith Oleson show that this no longer happens
with this change.

Changes from Sean Swenson to add minimum value to rlen (length of main channel); rlen values 
can be too small, leading to tlen values that are too large

M    src/riverroute/RtmIO.F90 - Explicitly set rlen to rlen_min early and remove code that was
        calculating tlen based on rlen_min when rlen<rlen_min. Since rlen is changed this affects code
        that is based on rlen, and can change the if statements that are triggered.

<hr>
# Tag name:  mosart1_0_27
### Originator(s): erik
### Date: Oct 03, 2017
### One-line Summary: Upgrade config_component to version 3, allow output file format
                  to change, fix a couple bugs

Bugs fixed: 2477, 2494

Testing: aux_clm testlist with cfgcompsetv3_n01_clm4_5_16_r251 on yellowstone/cheyenne/hobart

M    src/riverroute/RtmIO.F90 - Changes from Jim Edwards to change format of output NetCDF
           files, and remove public attribute from io_type

           Remove unneeded setting of namelist group name in test namelists (bug 2494)
M    cime_config/testdefs/testmods_dirs/mosart/default/user_nl_mosart
M    cime_config/testdefs/testmods_dirs/mosart/iceOff/user_nl_mosart
M    cime_config/testdefs/testmods_dirs/mosart/mosartOff/user_nl_mosart

M    cime_config/buildnml ------------- Fix bug 2477 for MOSART so that
        if NINST_RTM > 1, will check if REFCASE has instance name
        in it and then use it, but if not use it without instance name
M    cime_config/config_component.xml - Upgrade to version 3 format

<hr>
# Tag name:  mosart1_0_26
### Originator(s): erik
### Date: July 07, 2017
### One-line Summary: Update areas on routing file and add some comments

Testing: aux_clm testlist with clm4_5_16_r250 on yellowstone

Point to a new routing file that has area's updated by Hongyi. Sean
Swenson showed that they are within single precision roundoff to what he 
expects them to be. The dataset also updates lon to be -180-180
rather than 0-360. Straighten out formatting of lines about hlen
hillslope length in RtmMod and add some comments to it.

M       cime_config/namelist_definition_mosart.xml -- new frivinp_rtm
            converted to NetCDF3 format
M       src/riverroute/RtmMod.F90 -- straighten out formatting, remove
            tabs, add comments about hlen

<hr>
# Tag name:  mosart1_0_25
### Originator(s): andre
### Date: July 05, 2017
### One-line Summary: bugfix and update testlist

Update the testlist to use clm5 compset alias naming conventions.

Bugfix for bugzilla 2481 from Erik Kluzek.

Testing: clm test list run on branch before merge.

Modified files:

M    cime_config/testdefs/testlist_mosart.xml
M    src/riverroute/RtmTimeManager.F90


<hr>
# Tag name:  mosart1_0_24
### Originator(s): erik
### Date: Jun 01, 2017
### One-line Summary: Update routing file with metadata and in NetCDF3 format

Update routing file with extra metadata and in NetCDF3 format rather than
NetCDF4.

Testing:

Modified files: 

M       cime_config/namelist_definition_mosart.xml

<hr>
# Tag name:  mosart1_0_23
### Originator(s): andre
### Date: March 17, 2017
### One-line Summary: answer changing improvements and bugfix

Answer changing improvements to channel storage from HongYi Li for
faster spinup. Bugfix from Tony Craig to use correct delta time for
qgwl flux.

Changes reviewed by Sean Swenson.

Testing:

    Test suite: mosart - yellowstone gnu, intel, pgi
    Test baseline mosart1_0_22_clm4_5_14_r229
    Test status: pass, answer changing as expected when mosart on

    Test suite: aux_clm4{0,5} - yellowstone gnu, intel, pgi
    Test baseline: clm4_5_14_r229
    Test status: ok - only expected failures. Bit for bit except for
                 tests with mosart as expected.

    Dave Lawrence and Sean Swenson ran branch in a long cesm2 coupled
    simulation and verified the results were ok.

Modified files: 

    components/mosart/src/riverroute/RtmMod.F90


<hr>
# Tag name:  mosart1_0_22
### Originator(s): andre
### Date: March 17, 2017
### One-line Summary: histfile pointer bugfix bugz-2184

Bugfix from bugzilla issue 2184. The histfile infrastructure mosart
inherited from clm via rtm needs to use pointers to refer to
pio/netcdf file handles. Bug 2184 has never occured in mosart, but the
same problem existed.

Testing:

    NOTE: testing done with a franken-branch of clm4_5_14_r229 and
    mosart1_0_21.

    mosart - pass - yellowstone gnu, intel, pgi. baseline mosart1_0_21_clm4_5_14_r229.

    aux_clm4{0,5}: ok vs clm4_5_14_r229 baselines yellowstone gnu, intel, pgi
        Only clm expected fails


Modified files: fix for bugzilla #2184

    components/mosart/src/riverroute/RtmHistFile.F90


<hr>
# Tag name:  mosart1_0_21
### Originator(s): mvertens, andre
### Date: March 14, 2017
### One-line Summary: cime5 python namelist generation

Changes from Mariana Vertenstein to convert namelist generation to use
the cime5 python namelist infrastructure.

Testing:

  aux_clm45 vs clm4_5_14_r227 - passed except for expected failures.

  mosart vs clm4_5_14_r227 generated baselines - pass

  mvertens
    ERI_Ld9.f10_f10.IMCRUCLM50BGC.yellowstone_intel.clm-default
    ERS_Ld7.f19_g16.B1850.yellowstone_intel.allactive-defaultio -
        was bfb with cesm2_0_alpha06f

Deleted files: perl based namelist generation and xml files.
   bld
   bld/build-namelist
   bld/namelist_files
   bld/namelist_files/namelist_defaults_mosart.xml
   bld/namelist_files/namelist_definition_mosart.xml

Added files: python based namelist generation and xml files
   cime_config/buildlib
   cime_config/buildnml
   cime_config/config_component.xml
   cime_config/namelist_definition_mosart.xml

<hr>
# Tag name:  mosart1_0_20
### Originator(s): erik
### Date: Feb 21 2017
### One-line Summary: Fix an issue with nag, lower amount of log output

Fix bugs: 2393, 3424

Add some if(masterproc) for write statements, change some writes to aborts, bug 2393
Fix bug for nag compiler on hobart moving global save after mpif bug 2424


M       src/riverroute/RtmMod.F90 ---- Add if(masterproc)
M       src/riverroute/MOSART_physics_mod.F90 ---- Remove print statements before abort
M       src/riverroute/RtmSpmd.F90 -- Move the global save until after the mpif.h include
           so it can work on the nag compiler on hobart.
M       src/riverroute/RtmIO.F90 ---- Add if(masterproc)

<hr>
# Tag name:  mosart1_0_19
### Originator(s): swenson, sacks
### Date: Oct 17 2016
### One-line Summary: Treat irrigation specially, fix volr

(1) Together with corresponding cime and clm changes, treats irrigation as a
    separate flux. The point of this is to map irrigation withdrawals normalized
    by volr, to help prevent river channels from going to negative volumes.

(2) Fixes the volr field sent to the coupler

Requires the cime changes from https://github.com/ESMCI/cime/pull/681

In order for the new irrigation changes to have an effect, also requires the CLM
changes from the branch
https://svn-ccsm-models.cgd.ucar.edu/clm2/branches/limitirrig - although this
code will work correctly without the CLM changes (in this case, irrigation will
simply be 0, and the irrigation flux will be folded in to a different runoff
flux).

Changes are from Sean Swenson; reviewed, tested and brought to the trunk by Bill
Sacks.

Testing: mosart test suite (on yellowstone), from
https://svn-ccsm-models.cgd.ucar.edu/clm2/branch_tags/limitirrig_tags/limitirrig_n10_clm4_5_12_r196

All tests passed, but with answer changes in r2x_Flrr_volr and x2l_Flrr_volr, as
expected.

Also ran aux_clm45 test suite; most tests passed; I'm looking into remaining
failures and will update mosart if it turns out it was the source of any
problems.

M       src/riverroute/RunoffMod.F90
M       src/riverroute/RtmMod.F90
M       src/riverroute/RtmHistFlds.F90
M       src/cpl/rof_comp_mct.F90
M       src/cpl/rof_cpl_indices.F90

<hr>
# Tag name:  mosart1_0_18
### Originator(s): erik
### Date: Sep 13 2016
### One-line Summary: Add output frequency to history files


M   src/riverroute/RtmHistFile.Fi90 -- add output history frequency "time_period_freq"

<hr>
# Tag name:  mosart1_0_17
### Originator(s): erik
### Date: Apr 14 2016
### One-line Summary: Turn off for CLM_ACCELERATED_SPINUP="on" and fix a few bugs

Have MOSART react to CLM_ACCELERATED_SPINUP setting from CLM and turn itself off
by default if it's "on".

Also fix following bugs:

2307 -- assumes history files are no_leap calendar
2299 -- missing timer call
2230 -- for MOSART delt_save not initialized

M       bld/build-namelist ---- Check CLM_ACCELERATED_SPINUP, read do_rtm from
           defaults file and check based on CLM_ACCELERATED_SPINUP and MOSART_MODE
M       bld/namelist_files/namelist_defaults_mosart.xml - add do_rtm settings
M       src/riverroute/RtmRestFile.F90 -- Increase string length to 255
M       src/riverroute/RtmMod.F90 ------- Increase string length to 255
M       src/riverroute/RtmHistFile.F90 -- Increase filename length to 255
           check calendar type before writing to history file
M       src/riverroute/RtmIO.F90 -------- Remove commented out lines, add timer call (2299), 
           initialize delt_save (2230)
M       src/riverroute/RtmTimeManager.F90 Add NO_LEAP_C and GREGORIAN_C as public constants (2307)

<hr>
# Tag name:  mosart1_0_16
### Originator(s): swenson
### Date: Feb 8 2016
### One-line Summary: bugfix for budget diagnostic output

  bugfix: budget diagnostic output


M       src/riverroute/RtmMod.F90

<hr>
# Tag name:  mosart1_0_14
### Originator(s): swenson
### Date: Dec 22 2015
### One-line Summary: bugfix for negative river runoff

  Note(bja, 20151223) this was made directly as a trunk tag, but not put on the trunk.

  bugfix: minimum channel length

M       src/riverroute/RtmMod.F90

<hr>
# Tag name:  mosart1_0_13
### Originator(s): swenson, andre
### Date: Dec 21 2015
### One-line Summary: bugfix for negative river runoff

  bugfix: calculation of qgwl_volume must be multiplied by area.

M       src/riverroute/RtmMod.F90

<hr>
# Tag name:  mosart1_0_12
### Originator(s): swenson, andre
### Date: Dec 18 2015
### One-line Summary: negative river runoff changes

- Update default routing file.
- Add new namelist options bypass_routing_option,
  qgwl_runoff_option
- New method of handeling runoff terms to avoid negative runoff.
- Error checking on max length of history filenames.
- Add testdefs dir, testmods and testlist for integration to cime
  test system.
- remove rofdto from coupler interface fields.

M       bld/build-namelist
M       bld/namelist_files/namelist_defaults_mosart.xml
M       bld/namelist_files/namelist_definition_mosart.xml
A  +    cime_config/testdefs
A  +    cime_config/testdefs/testlist_mosart.xml
A  +    cime_config/testdefs/testmods_dirs
A  +    cime_config/testdefs/testmods_dirs/mosart
A  +    cime_config/testdefs/testmods_dirs/mosart/mosartOff
A  +    cime_config/testdefs/testmods_dirs/mosart/mosartOff/include_user_mods
A  +    cime_config/testdefs/testmods_dirs/mosart/mosartOff/user_nl_mosart
A  +    cime_config/testdefs/testmods_dirs/mosart/default
A  +    cime_config/testdefs/testmods_dirs/mosart/default/user_nl_mosart
A  +    cime_config/testdefs/testmods_dirs/mosart/iceOff
A  +    cime_config/testdefs/testmods_dirs/mosart/iceOff/include_user_mods
A  +    cime_config/testdefs/testmods_dirs/mosart/iceOff/user_nl_mosart
M       doc/ChangeLog
M       src/riverroute/RtmMod.F90
M       src/riverroute/RtmHistFlds.F90
M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RtmHistFile.F90
M       src/riverroute/RtmVar.F90
M       src/riverroute/RunoffMod.F90
M       src/cpl/rof_comp_mct.F90
M       src/cpl/rof_comp_esmf.F90

<hr>
# Tag name:  mosart1_0_11
### Originator(s): katec
### Date: Dec 7 2015
### One-line Summary: Changes to get Mosart to build with Nag on Hobart

- Shortened a few lines of code that were too long
- Added '.' for some literals that Nag thought were integers
- Removed call to 'system' and used open status='replace' instead

M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RtmMod.F90

<hr>
# Tag name:  mosart1_0_10
### Originator(s): tcraig
### Date: Dec 2 2015
### One-line Summary: New input file, update direct terms, update history file

- switch to input dataset MOSART_Global_half_20151130a.nc
- update direct sparse matrix to include non basin points in order
  to pass data from any grid cell directly to the ocean.
- modify the direct term and push all direct water to outlet points
- set all tracer 2 water (frozen water) to be a direct term
- add ability to skip some tracers in the euler solver via euler_calc flag
- add a budget accumulator term
- update history file fields and fieldnames, add new history fields
- compare total upstream basin areas from the input area field, use the
  computed total areas instead of the input total areas

M       bld/namelist_files/namelist_defaults_mosart.xml
M       src/riverroute/RtmMod.F90
M       src/riverroute/RtmHistFlds.F90
M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RtmHistFile.F90
M       src/riverroute/RunoffMod.F90

<hr>
# Tag name:  mosart1_0_09
### Originator(s): tcraig
### Date: Nov 29 2015
### One-line Summary: Code cleanup, add budget diagnostics, history files

M       src/riverroute/RtmMod.F90
M       src/riverroute/RtmHistFlds.F90
M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RtmHistFile.F90
M       src/riverroute/RunoffMod.F90
M       src/riverroute/RtmRestFile.F90
M       src/cpl/rof_cpl_indices.F90
M       src/cpl/rof_comp_mct.F90
M       src/cpl/rof_comp_esmf.F90
<hr>
# Tag name:  mosart1_0_08
### Originator(s): tcraig
### Date: Nov 24 2015
### One-line Summary: Fix exact restart in atan slope calc

M       src/riverroute/MOSART_physics_mod.F90

<hr>
# Tag name:  mosart1_0_07
### Originator(s): tcraig
### Date: Nov 22 2015
### One-line Summary: update rtmini and rtmrun routine and add budget
- code cleanup of rtmini and rtmrun
- works with mosart input files with scrambled IDs
- moved dto term into rtmrun
- added direct-to-outlet tranfer capability
- removed a bunch of old rtm code
- fixed esmf interfaces and tested in DEBUG mode
- added budget calculation (still being validated)
- has a known exact restart error that introduces a roundoff
  difference at the first timestep at a handful of gridcells.
  This is probably not going to impact science, will be fixed
  next.

M       src/riverroute/RtmMod.F90
M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RunoffMod.F90
M       src/riverroute/RtmRestFile.F90
M       src/cpl/rof_comp_mct.F90
M       src/cpl/rof_comp_esmf.F90
<hr>
# Tag name:  mosart1_0_06
### Originator(s): tcraig
### Date: Nov 19 2015
### One-line Summary: merge ACME fixes to decomp and performance

This works with MOSART_Global_half_20130604a.nc, NOT
MOSART_Global_half_20151015a.nc.  This will be fixed in the next
commit.  Probably shouldn't use this tag for now.

Not bit for big with previous tag.

M       bld/build-namelist
M       bld/namelist_files/namelist_defaults_mosart.xml
M       bld/namelist_files/namelist_definition_mosart.xml
M       src/riverroute/RtmMod.F90
M       src/riverroute/MOSART_physics_mod.F90
M       src/riverroute/RtmSpmd.F90
M       src/riverroute/RtmVar.F90
M       src/riverroute/RunoffMod.F90

<hr>
# Tag name:  mosart1_0_05
### Originator(s): andre
### Date: Oct 15 2015
### One-line Summary: swenson bugfix for mosart direction file

The old mosart direction file, rof/mosart/MOSART_Global_half_20130604a.nc,
has antarctica shifted by 180 degrees.

Tested with pre-clm4_5_3_r140:
  ERS_D_Ld5.f10_f10.IMCRUCLM50BGC.yellowstone_intel - runs to completion

Not bit for bit with previous tag.

<hr>
# Tag name:  mosart1_0_04
### Originator(s): andre
### Date: Oct 15 2015
### One-line Summary: swenson river volume normalization bugfix

Tested by Sean Swenson. Verified to compile and run SMS_D.f10_f10.IMCRUCLM50BGC.yellowstone_intel.clm-default.

Not bit for bit with previous tag.

<hr>
# Tag name:  mosart1_0_03
### Originator(s): andre
### Date: Oct 13, 2015
### One-line Summary: update mosart

Updates to mosart:

* cime compatible infrastructure from Mariana Vertenstein

* Add direct to ocean runoff flux from Sean Swenson.

* PIO2 updates from Jim Edwards

Tested against version of clm4_5_3_r135

  ERS_D_Ld5.f10_f10.IMCRUCLM50BGC.yellowstone_intel.clm-default
  SMS_D_Ld3.f10_f10.IMCRUCLM50BGC.yellowstone_intel.clm-default


Not expected to be bit for bit with previous tag.

M       bld/build-namelist
D       bld/mosart.buildlib
D       bld/mosart.buildnml
D       bld/user_nl_mosart
A  +    cime_config
A  +    cime_config/buildlib
A  +    cime_config/buildnml
A  +    cime_config/config_component.xml
A  +    cime_config/user_nl_mosart
M       src/riverroute/RtmMod.F90
M       src/riverroute/RtmHistFlds.F90
M       src/riverroute/RtmIO.F90
M       src/riverroute/RunoffMod.F90
M       src/cpl/rof_cpl_indices.F90
M       src/cpl/rof_comp_mct.F90
M       src/cpl/rof_comp_esmf.F90

<hr>
# Tag name:  mosart1_0_00
### Originator(s): tcraig
### Date: May 1, 2015
### One-line Summary: add mosart to CESM repository

This is based on the following version from PNL but has been updated
to fit into cesm1_4.

URL: https://svn.pnl.gov/svn/iRESM/cesm1/trunk/models/rof/mosart
Repository Root: https://svn.pnl.gov/svn/iRESM
Repository UUID: 97a048bb-0f8f-0410-8848-820bd1cc90bf
Revision: 1311
Node Kind: directory
Schedule: normal
Last Changed Author: tcraig
Last Changed Rev: 1310
Last Changed Date: 2014-10-03 10:38:55 -0600 (Fri, 03 Oct 2014)

A       bld
A       bld/build-namelist
A       bld/mosart.buildlib
A       bld/mosart.buildnml
A       bld/user_nl_mosart
A       bld/namelist_files
A       bld/namelist_files/namelist_defaults_mosart.xml
A       bld/namelist_files/namelist_definition_mosart.xml
A       doc
A       doc/ChangeLog
A       src
A       src/riverroute
A       src/riverroute/RtmMod.F90
A       src/riverroute/RtmFileUtils.F90
A       src/riverroute/RtmHistFlds.F90
A       src/riverroute/MOSART_physics_mod.F90
A       src/riverroute/RtmSpmd.F90
A       src/riverroute/RtmHistFile.F90
A       src/riverroute/RtmIO.F90
A       src/riverroute/RtmVar.F90
A       src/riverroute/RtmTimeManager.F90
A       src/riverroute/RtmDateTime.F90
A       src/riverroute/RunoffMod.F90
A       src/riverroute/RtmRestFile.F90
A       src/cpl
A       src/cpl/rof_cpl_indices.F90
A       src/cpl/rof_comp_mct.F90
A       src/cpl/rof_comp_esmf.F90

<hr>
