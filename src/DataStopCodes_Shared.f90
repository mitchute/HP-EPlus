MODULE DataStopCodes      ! EnergyPlus Data-Only Module

  USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength  !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)

  !globally declare program exit codes to alert calling routine to exit status
  !will become useful during parametrics, test suites, or otherwise
  !exit codes should be 0 for success (not necessary), or 1-255 for specific error codes
  !note that sending the exit status back to the OS is not mandatory according to FORTRAN standards, but Intel 11 should do it properly

  IMPLICIT NONE

  PUBLIC

  TYPE StopCodeInformation
    INTEGER :: ExitCode !0 to 255 for POSIX standards with 0 being a successful exit
    CHARACTER(LEN=MaxNameLength) :: Message
  END TYPE

  ! normal exit status
  INTEGER, PARAMETER :: exit_Normal = 0

  ! File I/O type exits
  INTEGER, PARAMETER :: exit_FileIO_Missing_HPData = 2

  ! Code diagnostics
  INTEGER, PARAMETER :: exit_Diagnostic_RefrigerantName = 52

  ! Other exit types
  INTEGER, PARAMETER :: exit_SimProblem_BadInitialization = 102
  INTEGER, PARAMETER :: exit_SimProblem_EnergyPlusProblem = 108


  ! Main stop code array
  TYPE(StopCodeInformation), DIMENSION(5), PARAMETER :: StopCodes = (/ &
  StopCodeInformation(exit_Normal,                       "Normal Exit"),                              &
  StopCodeInformation(exit_FileIO_Missing_HPData,        "Could not find heat pump data input file"), &
  StopCodeInformation(exit_Diagnostic_RefrigerantName,   "Bad refrigerant name input"),               &
  StopCodeInformation(exit_SimProblem_BadInitialization, "Bad initialization for simulation"),        &
  StopCodeInformation(exit_SimProblem_EnergyPlusProblem, "A problem was encountered by the input processor") /)

END MODULE
