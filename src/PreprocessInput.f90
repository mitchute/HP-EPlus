! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
!This module pre-processes the input for the simulation
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This does not represent a physical component of the system

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module removes previously backed-up input files and backs up the current files.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! Does not open any files or output any files
! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! NONE

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 methods:
!    PUBLIC PreProcessInput
!      ORNLSolver

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! no apparent issues

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! YEAR-12-29 | JEH | Header completion

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! No apparent to do items, besides removing extraneous lines
    
    MODULE InputPreProcessor

IMPLICIT NONE !require explicit declaration of everything

PRIVATE !unless otherwise declared public

PUBLIC PreProcessInput

CONTAINS

    SUBROUTINE PreProcessInput
    
        !this routine will rename hpdata.idf, run epmacro to get the fluid properties, and rename that to in.idf for the E+ input processor
    
!DEC$ IF DEFINED(_WIN32)
        ! delete the previously backed up version of the input file
        call system('if exist inBackup.idf del inBackup.idf')
        ! backup the last input file
        call system('if exist in.idf rename in.idf inBackup.idf')
        ! rename the heat pump input file in preparation for epmacro
        !call system('copy Minimal_HPSim.idf in.imf')
        call system('copy ZoneWSHP_wDOAS_HPSim.idf in.imf') !RS: Debugging: Temporarily commenting it out
        !call system('copy HeatPump_Unitary.idf in.imf') !RS: Debugging: Temporarily testing this case
        ! call epmacro on it
        call system('EPMacro.exe')
        ! now rename the file to be read by the E+ input processor
        call system('rename out.idf in.idf')
!DEC$ ELSEIF DEFINED(__linux)
        ! delete the previously backed up version of the input file
        call system('rm -f inBackup.idf > /dev/null')
        ! backup the last input file
        call system('if [ -f in.idf ]; then mv in.idf inBackup.idf > /dev/null; fi')
        ! rename the heat pump input file in preparation for epmacro
        call system('cp HPdataUnix.idf in.imf > /dev/null')
        ! call epmacro on it
        call system('./epmacro')
        ! now rename the file to be read by the E+ input processor
        call system('mv out.idf in.idf > /dev/null')
!DEC$ ELSE
        !write(*,*) 'Get off your mac!! :)'
!DEC$ ENDIF
    
    END SUBROUTINE

END MODULE
