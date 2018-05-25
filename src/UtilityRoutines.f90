
! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! Provide a 1 or 2 sentence overview of this module.
! In most cases, it is probably not a useful entry and can be inferred from the name of the module anyway.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component represents something...or nothing...in a heat pump system.
! A description of the component is found at:
! some website
! From that website:
!  - It does something

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! Here's a one line summary of what this does for the simulation itself.
! This module takes inputs such as...and modifies them like so...and outputs these things

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! Check for any OPEN statements in the code
! Check for any WRITE statements in the code
!  Note that writing to unit "*" or "6" means just write to the terminal, not to a file

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! What vars and structures are defined at the *module* level...are units defined?  Any other notes?

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains X methods:
!    PUBLIC InitSomething -- What does this do (in one line)?
!      Called by what other modules?

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! Are there any interesting issues with this, unfuddle ticket numbers?

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! YEAR-MM-DD | ABC | Some other log message?

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Put some notes for what needs to happen here
! Silly things are fine
! Somethings these small silly things are great to grab on to when starting up with the code after being off for a while


!SUBROUTINE IssueHPFatalError(exitCode)
!
!! the fortran keyword STOP cannot accept a variable, only a literal or a parameter
!! thus we need a ridiculous case statement for all possibilities found in DataStopCodes.f90
!
!    USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength  !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
!    USE DataStopCodes
!    implicit none
!
!    INTEGER, INTENT(IN) :: exitCode
!
!    INTEGER :: Counter
!    CHARACTER(LEN=MaxNameLength) :: CodeMessage
!
!    DO Counter = 1, SIZE(StopCodes)
!        IF (exitCode == StopCodes(Counter)%ExitCode) THEN
!            CodeMessage = StopCodes(Counter)%Message
!            EXIT
!        END IF
!    END DO
!
!    WRITE(*,*) '-+-+-+-+-+-+-+-'
!    WRITE(*,*) 'Heat pump simulation fatal error!'
!    WRITE(*,*) 'Error explanation: '//TRIM(CodeMessage)
!    WRITE(*,*) 'Exit code follows:'
!    SELECT CASE (exitCode)
!    CASE (exit_FileIO_Missing_HPData)
!        STOP exit_FileIO_Missing_HPData
!    CASE (exit_Diagnostic_RefrigerantName)
!        STOP exit_Diagnostic_RefrigerantName
!    CASE (exit_SimProblem_BadInitialization)
!        STOP exit_SimProblem_BadInitialization
!    CASE (exit_SimProblem_EnergyPlusProblem)
!        STOP exit_SimProblem_EnergyPlusProblem
!    CASE DEFAULT
!        WRITE(*,*) '-+-Diagnostic-+- Unimplemented stop code in UtilityRoutines::IssueHPFatalError'
!        STOP 1
!    END SELECT
!
!END SUBROUTINE
!
!LOGICAL FUNCTION IssueRefPropError(RefPropErrValue, CallingRoutine, ValueIfErrorFound, VariableToSet1, VariableToSet2) RESULT (ErrorFound)
!    implicit none
!
!    INTEGER(2), INTENT(IN) :: RefPropErrValue ! the value that was returned from the RefProp call
!    CHARACTER(len=*), INTENT(IN) :: CallingRoutine ! an identifier to the routine calling me, for reporting
!    INTEGER, INTENT(IN), OPTIONAL :: ValueIfErrorFound ! if RefProp was erroneous, this is the signaling value to be used
!    INTEGER, INTENT(INOUT), OPTIONAL :: VariableToSet1 ! if RefProp was erroneous, this will be set to the signal value
!    REAL, INTENT(INOUT), OPTIONAL :: VariableToSet2 ! another variable to set...optionally
!
!    IF ( (PRESENT(VariableToSet1) .OR. PRESENT(VariableToSet2)) .AND.  .NOT. PRESENT(ValueIfErrorFound) ) THEN
!        !malformed, how are we going to assign variables if we don't have a value to assign with
!        WRITE(*,*) '-+-Diagnostic-+- Improper call to IssueRefPropError, callingroutine = '//CallingRoutine
!    END IF
!
!    ErrorFound = .FALSE.
!    IF (RefPropErrValue .GT. 0) THEN
!        CALL ShowWarningError(CallingRoutine//': RefProp lookup error')
!        IF ( PRESENT ( VariableToSet1 ) ) THEN
!            VariableToSet1 = ValueIfErrorFound
!        END IF
!        IF ( PRESENT ( VariableToSet2 ) ) THEN
!            VariableToSet2 = REAL(ValueIfErrorFound)
!        END IF
!        ErrorFound = .TRUE.
!    END IF
!
!    RETURN
!
!END FUNCTION
!
!SUBROUTINE IssueOutputMessage(Message)
!    implicit none
!
!    CHARACTER(LEN=*), INTENT(IN) :: Message
!
!    WRITE(6,*) Message
!    WRITE(*,*) Message
!
!END SUBROUTINE
SUBROUTINE IssueHPFatalError(exitCode)



  ! the fortran keyword STOP cannot accept a variable, only a literal or a parameter
  ! thus we need a ridiculous case statement for all possibilities found in DataStopCodes.f90

  USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength  !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  USE DataStopCodes
  implicit none

  INTEGER, INTENT(IN) :: exitCode

  INTEGER :: Counter
  CHARACTER(LEN=MaxNameLength) :: CodeMessage

  logical :: exist
  inquire(file='Log.audit', exist=exist)
  if (exist) then
  else
    open(122, file="Log.audit", status="new", action="write")
  end if
  close(122)

  DO Counter = 1, SIZE(StopCodes)
    IF (exitCode == StopCodes(Counter)%ExitCode) THEN
      CodeMessage = StopCodes(Counter)%Message
      EXIT
    END IF
  END DO

  WRITE(*,*) '-+-+-+-+-+-+-+-'
  WRITE(*,*) 'Heat pump simulation fatal error!'
  WRITE(*,*) 'Error explanation: '//TRIM(CodeMessage)
  WRITE(*,*) 'Exit code follows:'
  SELECT CASE (exitCode)
  CASE (exit_FileIO_Missing_HPData)
    OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
    WRITE(19,*) 'Initializing "Not Converged" file'
    WRITE(*,*)'Initializing "Not Converged" file'
    CLOSE(19)
    OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
    CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
    STOP exit_FileIO_Missing_HPData
  CASE (exit_Diagnostic_RefrigerantName)
    OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
    WRITE(19,*) 'Initializing "Not Converged" file'
    WRITE(*,*) 'Initializing "Not Converged" file'
    CLOSE(19)
    OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
    CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
    STOP exit_Diagnostic_RefrigerantName
  CASE (exit_SimProblem_BadInitialization)
    OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
    WRITE(19,*) 'Initializing "Not Converged" file'
    WRITE(*,*) 'Initializing "Not Converged" file'
    CLOSE(19)
    OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
    CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
    STOP exit_SimProblem_BadInitialization
  CASE (exit_SimProblem_EnergyPlusProblem)
    OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
    WRITE(19,*) 'Initializing "Not Converged" file'
    WRITE(*,*) 'Initializing "Not Converged" file'
    CLOSE(19)
    OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
    CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
    STOP exit_SimProblem_EnergyPlusProblem
  CASE DEFAULT
    OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
    WRITE(19,*) 'Initializing "Not Converged" file'
    WRITE(*,*) 'Initializing "Not Converged" file'
    CLOSE(19)
    WRITE(*,*) '-+-Diagnostic-+- Unimplemented stop code in UtilityRoutines::IssueHPFatalError'
    OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
    CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
    STOP 1
  END SELECT
  !READ(*,*)
END SUBROUTINE

LOGICAL FUNCTION IssueRefPropError(RefPropErrValue, CallingRoutine, ValueIfErrorFound, VariableToSet1, VariableToSet2) RESULT (ErrorFound)
  implicit none

  INTEGER(2), INTENT(IN) :: RefPropErrValue ! the value that was returned from the RefProp call
  CHARACTER(len=*), INTENT(IN) :: CallingRoutine ! an identifier to the routine calling me, for reporting
  INTEGER, INTENT(IN), OPTIONAL :: ValueIfErrorFound ! if RefProp was erroneous, this is the signaling value to be used
  INTEGER, INTENT(INOUT), OPTIONAL :: VariableToSet1 ! if RefProp was erroneous, this will be set to the signal value
  REAL, INTENT(INOUT), OPTIONAL :: VariableToSet2 ! another variable to set...optionally
  logical :: exist
  inquire(file='Log.audit', exist=exist)
  if (exist) then
  else
    open(122, file="Log.audit", status="new", action="write")
  end if
  close(122)
  IF ( (PRESENT(VariableToSet1) .OR. PRESENT(VariableToSet2)) .AND.  .NOT. PRESENT(ValueIfErrorFound) ) THEN
    !malformed, how are we going to assign variables if we don't have a value to assign with
    WRITE(*,*) '-+-Diagnostic-+- Improper call to IssueRefPropError, callingroutine = '//CallingRoutine
  END IF

  ErrorFound = .FALSE.
  IF (RefPropErrValue .GT. 0) THEN    !RS: Debugging: Having an issue with CallingRoutine being an undefined address (11/4/14)
    CALL ShowWarningError(': RefProp lookup error') !CallingRoutine//': RefProp lookup error')
    !IF ( PRESENT ( VariableToSet1 ) ) THEN
    !    VariableToSet1 = ValueIfErrorFound
    !END IF
    !IF ( PRESENT ( VariableToSet2 ) ) THEN
    !    VariableToSet2 = REAL(ValueIfErrorFound)
    !END IF
    ErrorFound = .TRUE.
  END IF

  RETURN
  !READ(*,*)
END FUNCTION

SUBROUTINE IssueOutputMessage(Message)
  implicit none
  CHARACTER(LEN=*), INTENT(IN) :: Message
  logical :: exist
  inquire(file='Log.audit', exist=exist)
  if (exist) then
  else
    open(122, file="Log.audit", status="new", action="write")
  end if
  close(122)
  open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
  WRITE(111,*) Message
  WRITE(*,*) Message
  close(111)
  !READ(*,*)
END SUBROUTINE
!
!SUBROUTINE AbortEnergyPlus
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   December 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine causes the program to halt due to a fatal error.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Puts a message on output files.
!          ! Closes files.
!          ! Stops the program.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
!  USE DataStopCodes
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!          ! na
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!    SUBROUTINE ShowMessage(Message)
!        CHARACTER(len=*) Message
!    END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  INTEGER tempfl
!  INTEGER, EXTERNAL :: GetNewUnitNumber
!  CHARACTER(len=20) NumWarnings
!  CHARACTER(len=20) NumSevere
!  logical :: exist
!  inquire(file='Log.audit', exist=exist)
!  if (exist) then
!  else
!    open(122, file="Log.audit", status="new", action="write")
!  end if
!  close(122)
!  open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
!  WRITE(111,*) TotalWarningErrors
!  WRITE(*,*) TotalWarningErrors
!  NumWarnings=ADJUSTL(NumWarnings)
!  WRITE(111,*) TotalSevereErrors
!  WRITE(*,*) TotalSevereErrors
!  NumSevere=ADJUSTL(NumSevere)
!  close(111)
!  CALL ShowMessage('EnergyPlus Terminated--Fatal Error Detected. '//TRIM(NumWarnings)//' Warning; '//  &
!                           TRIM(NumSevere)//' Severe Errors')
!  tempfl=GetNewUnitNumber()
!  open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
!
!  write(111,*) 'EnergyPlus Terminated--Fatal Error Detected. '//TRIM(NumWarnings)//' Warning; '//  &
!                           TRIM(NumSevere)//' Severe Errors'
!  write(*,*) 'EnergyPlus Terminated--Fatal Error Detected. '//TRIM(NumWarnings)//' Warning; '//  &
!                           TRIM(NumSevere)//' Severe Errors'
!  close(111)
!  CALL CloseMiscOpenFiles
!  CALL IssueHPFatalError(exit_SimProblem_EnergyPlusProblem)
!  !READ(*,*)
!END SUBROUTINE AbortEnergyPlus
!
!SUBROUTINE CloseMiscOpenFiles
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   December 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine scans potential unit numbers and closes
!          ! any that are still open.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Use INQUIRE to determine if file is open.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
! ! USE DaylightingManager, ONLY: CloseReportIllumMaps
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!          ! na
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!   INTEGER, PARAMETER :: MaxUnitNumber = 1000
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!          ! na
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!
!      LOGICAL :: exists, opened
!      INTEGER :: UnitNumber
!      INTEGER :: ios
!      logical :: exist
!  inquire(file='Log.audit', exist=exist)
!  if (exist) then
!  else
!    open(122, file="Log.audit", status="new", action="write")
!  end if
!  close(122)
!      !CALL CloseReportIllumMaps
!      DO UnitNumber = 1, MaxUnitNumber
!         INQUIRE (UNIT = UnitNumber, EXIST = exists,  OPENED = opened, IOSTAT = ios)
!         IF (exists .and. opened .and. ios == 0) THEN
!             CLOSE(UnitNumber)
!         END IF
!      END DO
!  RETURN
!END SUBROUTINE CloseMiscOpenFiles
!
!SUBROUTINE EndEnergyPlus
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   December 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine causes the program to terminate when complete (no errors).
!
!          ! METHODOLOGY EMPLOYED:
!          ! Puts a message on output files.
!          ! Closes files.
!          ! Stops the program.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
!  USE InputProcessor_HPSim
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!          ! na
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!    SUBROUTINE ShowMessage(Message)
!        CHARACTER(len=*) Message
!    END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  INTEGER tempfl
!  INTEGER, EXTERNAL :: GetNewUnitNumber
!  CHARACTER(len=20) NumWarnings
!  CHARACTER(len=20) NumSevere
!  CHARACTER(len=25) Elapsed
!  INTEGER Hours   ! Elapsed Time Hour Reporting
!  INTEGER Minutes ! Elapsed Time Minute Reporting
!  INTEGER Seconds ! Elapsed Time Second Reporting
!logical :: exist
!  inquire(file='Log.audit', exist=exist)
!  if (exist) then
!  else
!    open(122, file="Log.audit", status="new", action="write")
!  end if
!  close(122)
!  open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
!  NumWarnings=ADJUSTL(NumWarnings)
!  WRITE(111,*) TotalSevereErrors
!  WRITE(*,*) TotalSevereErrors
!  close(111)
!  NumSevere=ADJUSTL(NumSevere)
!  Hours=Elapsed_Time/3600.
!  Elapsed_Time=Elapsed_Time-Hours*3600
!  Minutes=Elapsed_Time/60.
!  Elapsed_Time=Elapsed_Time-Minutes*60
!  Seconds=Elapsed_Time
!  WRITE(Elapsed,"(I2.2,'hr ',I2.2,'min ',I2.2,'sec')") Hours,Minutes,Seconds
!
!  CALL ShowMessage('EnergyPlus Completed Successfully-- '//TRIM(NumWarnings)//' Warning; '//TRIM(NumSevere)//' Severe Errors;'// &
!                   ' Elapsed Time='//TRIM(Elapsed))
!  tempfl=GetNewUnitNumber()
!  open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
!  write(111,'(A)') 'EnergyPlus Completed Successfully-- '//TRIM(NumWarnings)//' Warning; '//TRIM(NumSevere)//' Severe Errors'
!  write(*,'(A)') 'EnergyPlus Completed Successfully-- '//TRIM(NumWarnings)//' Warning; '//TRIM(NumSevere)//' Severe Errors'
!  close(111)
!  CALL CloseMiscOpenFiles
!  CALL DeallocateArrays !-ISI 02/23/04
!!  READ(*,*)
!  RETURN
!
!END SUBROUTINE EndEnergyPlus
!
!FUNCTION GetNewUnitNumber ()  RESULT (UnitNumber)
!
!          ! FUNCTION INFORMATION:
!          !       AUTHOR         Linda K. Lawrie, adapted from reference
!          !       DATE WRITTEN   September 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS FUNCTION:
!          ! Returns a unit number of a unit that can exist and is not connected.  Note
!          ! this routine does not magically mark that unit number in use.  In order to
!          ! have the unit "used", the source code must OPEN the file.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Use Inquire function to find out if proposed unit: exists or is opened.
!          ! If not, can be used for a new unit number.
!
!          ! REFERENCES:
!          ! Copyright (c) 1994 Unicomp, Inc.  All rights reserved.
!          !
!          ! Developed at Unicomp, Inc.
!          !
!          ! Permission to use, copy, modify, and distribute this
!          ! software is freely granted, provided that this notice
!          ! is preserved.
!
!          ! USE STATEMENTS:
!          ! na
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! FUNCTION ARGUMENT DEFINITIONS:
!  INTEGER UnitNumber  ! Result from scanning currently open files
!
!          ! FUNCTION PARAMETER DEFINITIONS:
!!  IO Status Values:
!
!  INTEGER, PARAMETER :: END_OF_RECORD = -2
!  INTEGER, PARAMETER :: END_OF_FILE = -1
!
!!  Indicate default input and output units:
!
!  INTEGER, PARAMETER :: DEFAULT_INPUT_UNIT = 5
!  INTEGER, PARAMETER :: DEFAULT_OUTPUT_UNIT = 6
!
!!  Indicate number and value of preconnected units
!
!  INTEGER, PARAMETER :: NUMBER_OF_PRECONNECTED_UNITS = 2
!  INTEGER, PARAMETER :: PRECONNECTED_UNITS (NUMBER_OF_PRECONNECTED_UNITS) = (/ 5, 6 /)
!
!!  Largest allowed unit number (or a large number, if none)
!  INTEGER, PARAMETER :: MaxUnitNumber = 1000
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!          ! na
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
!  LOGICAL :: exists  ! File exists
!  LOGICAL :: opened  ! Unit is open
!  INTEGER :: ios     ! return value from Inquire intrinsic
!  logical :: exist
!  inquire(file='Log.audit', exist=exist)
!  if (exist) then
!  else
!    open(122, file="Log.audit", status="new", action="write")
!  end if
!  close(122)
!  DO UnitNumber = 1, MaxUnitNumber
!    IF (UnitNumber == DEFAULT_INPUT_UNIT .or. &
!        UnitNumber == DEFAULT_OUTPUT_UNIT) THEN
!        CYCLE
!    END IF
!    IF (ANY (UnitNumber == PRECONNECTED_UNITS)) THEN
!        CYCLE
!    END IF
!    INQUIRE (UNIT = UnitNumber, EXIST = exists,  OPENED = opened, IOSTAT = ios)
!    IF (exists .and. .not. opened .and. ios == 0) THEN
!        RETURN      ! result is set in UnitNumber
!    END IF
!  END DO
!  UnitNumber = -1
!END FUNCTION GetNewUnitNumber

SUBROUTINE ShowFatalError(ErrorMessage)
  IMPLICIT NONE
  CHARACTER(len=*) ErrorMessage
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message)
      CHARACTER(len=*) Message
    END SUBROUTINE
  END INTERFACE
  CALL ShowErrorMessage(' **  Fatal  ** '//ErrorMessage)
  CALL AbortEnergyPlus
  RETURN
END SUBROUTINE ShowFatalError

SUBROUTINE ShowSevereError(ErrorMessage)
  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  IMPLICIT NONE
  CHARACTER(len=*) ErrorMessage
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message)
      CHARACTER(len=*) Message
    END SUBROUTINE
  END INTERFACE
  TotalSevereErrors=TotalSevereErrors+1
  CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage)
  RETURN
END SUBROUTINE ShowSevereError

SUBROUTINE ShowContinueError(Message)
  IMPLICIT NONE
  CHARACTER(len=*) Message
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message)
      CHARACTER(len=*) Message
    END SUBROUTINE
  END INTERFACE
  CALL ShowErrorMessage(' **   ~~~   ** '//Message)
  RETURN
END SUBROUTINE ShowContinueError

SUBROUTINE ShowMessage(MessageVal)
  CHARACTER(len=*) MessageVal
  INTERFACE
    SUBROUTINE ShowErrorMessage(MessageVal)
      CHARACTER(len=*) MessageVal
    END SUBROUTINE
  END INTERFACE
  CALL ShowErrorMessage(' ************* '//MessageVal)
  RETURN
END SUBROUTINE ShowMessage

SUBROUTINE ShowWarningError(ErrorMessage)
  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  IMPLICIT NONE
  CHARACTER(len=*) ErrorMessage
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message)
      CHARACTER(len=*) Message
    END SUBROUTINE
  END INTERFACE
  TotalWarningErrors=TotalWarningErrors+1
  CALL ShowErrorMessage(' ** Warning ** '//ErrorMessage)
  RETURN
END SUBROUTINE ShowWarningError
!
SUBROUTINE ShowErrorMessage(ErrorMessage)
  USE DataGlobals_HPSimIntegrated
  IMPLICIT NONE
  CHARACTER(len=*) ErrorMessage
  CHARACTER(len=120) :: VerString='EnergyPlus, Version 1.1.1'
  CHARACTER(len=*), PARAMETER :: ErrorFormat='(2X,A)'
  INTEGER  :: TotalErrors=0        ! used to determine when to open standard error output file.
  INTEGER  :: StandardErrorOutput
  INTEGER,EXTERNAL  :: GetNewUnitNumber
  logical :: exist
  SAVE     TotalErrors,StandardErrorOutput
  inquire(file='Log.audit', exist=exist)
  if (exist) then
  else
    open(122, file="Log.audit", status="new", action="write")
  end if
  close(122)
  IF (TotalErrors .eq. 0) THEN
    StandardErrorOutput=GetNewUnitNumber()
    open(111,file='Log.audit',status='old',action='write',form='formatted',position="append")
    WRITE(111,'(A)') 'Program Version,'//TRIM(VerString)
    WRITE(*,'(A)') 'Program Version,'//TRIM(VerString)
  ENDIF
  TotalErrors=TotalErrors+1
  WRITE(111,ErrorFormat) TRIM(ErrorMessage)
  WRITE(*,ErrorFormat) TRIM(ErrorMessage)
  close(111)
  !READ(*,*)
  RETURN
END SUBROUTINE ShowErrorMessage


SUBROUTINE AbortEnergyPlus

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   December 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine causes the program to halt due to a fatal error.

  ! METHODOLOGY EMPLOYED:
  ! Puts a message on output files.
  ! Closes files.
  ! Stops the program.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  USE DataStopCodes

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE

    SUBROUTINE ShowMessage(Message,Unit1,Unit2)
      CHARACTER(len=*) Message
      INTEGER, OPTIONAL :: Unit1
      INTEGER, OPTIONAL :: Unit2
    END SUBROUTINE

  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER tempfl
  INTEGER, EXTERNAL :: GetNewUnitNumber
  CHARACTER(len=20) NumWarnings
  CHARACTER(len=20) NumSevere

  WRITE(NumWarnings,*) TotalWarningErrors
  NumWarnings=ADJUSTL(NumWarnings)
  WRITE(NumSevere,*) TotalSevereErrors
  NumSevere=ADJUSTL(NumSevere)

  CALL ShowMessage('EnergyPlus Terminated--Fatal Error Detected. '//TRIM(NumWarnings)//' Warning; '//  &
  TRIM(NumSevere)//' Severe Errors')
  tempfl=GetNewUnitNumber()
  open(tempfl,file='eplusout.end')
  write(tempfl,*) 'EnergyPlus Terminated--Fatal Error Detected. '//TRIM(NumWarnings)//' Warning; '//  &
  TRIM(NumSevere)//' Severe Errors'
  close(tempfl)
  CALL CloseMiscOpenFiles

  CALL IssueHPFatalError(exit_SimProblem_EnergyPlusProblem)

END SUBROUTINE AbortEnergyPlus

SUBROUTINE CloseMiscOpenFiles

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   December 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine scans potential unit numbers and closes
  ! any that are still open.

  ! METHODOLOGY EMPLOYED:
  ! Use INQUIRE to determine if file is open.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! USE DaylightingManager, ONLY: CloseReportIllumMaps

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  INTEGER, PARAMETER :: MaxUnitNumber = 1000

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  LOGICAL :: exists, opened
  INTEGER :: UnitNumber
  INTEGER :: ios

  !CALL CloseReportIllumMaps

  DO UnitNumber = 1, MaxUnitNumber
    INQUIRE (UNIT = UnitNumber, EXIST = exists,  OPENED = opened, IOSTAT = ios)
    IF (exists .and. opened .and. ios == 0) THEN
      CLOSE(UnitNumber)
    END IF
  END DO

  RETURN

END SUBROUTINE CloseMiscOpenFiles

SUBROUTINE EndEnergyPlus

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   December 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine causes the program to terminate when complete (no errors).

  ! METHODOLOGY EMPLOYED:
  ! Puts a message on output files.
  ! Closes files.
  ! Stops the program.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataGlobals_HPSimIntegrated   !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  USE InputProcessor

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE

    SUBROUTINE ShowMessage(Message,Unit1,Unit2)
      CHARACTER(len=*) Message
      INTEGER, OPTIONAL :: Unit1
      INTEGER, OPTIONAL :: Unit2
    END SUBROUTINE

  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER tempfl
  INTEGER, EXTERNAL :: GetNewUnitNumber
  CHARACTER(len=20) NumWarnings
  CHARACTER(len=20) NumSevere
  CHARACTER(len=25) Elapsed
  INTEGER Hours   ! Elapsed Time Hour Reporting
  INTEGER Minutes ! Elapsed Time Minute Reporting
  INTEGER Seconds ! Elapsed Time Second Reporting

  WRITE(NumWarnings,*) TotalWarningErrors
  NumWarnings=ADJUSTL(NumWarnings)
  WRITE(NumSevere,*) TotalSevereErrors
  NumSevere=ADJUSTL(NumSevere)
  Hours=Elapsed_Time/3600.
  Elapsed_Time=Elapsed_Time-Hours*3600
  Minutes=Elapsed_Time/60.
  Elapsed_Time=Elapsed_Time-Minutes*60
  Seconds=Elapsed_Time
  WRITE(Elapsed,"(I2.2,'hr ',I2.2,'min ',I2.2,'sec')") Hours,Minutes,Seconds

  CALL ShowMessage('EnergyPlus Completed Successfully-- '//TRIM(NumWarnings)//' Warning; '//TRIM(NumSevere)//' Severe Errors;'// &
  ' Elapsed Time='//TRIM(Elapsed))
  tempfl=GetNewUnitNumber()
  open(tempfl,file='eplusout.end')
  write(tempfl,'(A)') 'EnergyPlus Completed Successfully-- '//TRIM(NumWarnings)//' Warning; '//TRIM(NumSevere)//' Severe Errors'
  close(tempfl)
  CALL CloseMiscOpenFiles
  CALL DeallocateArrays !-ISI 02/23/04

  RETURN

END SUBROUTINE EndEnergyPlus

FUNCTION GetNewUnitNumber ()  RESULT (UnitNumber)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie, adapted from reference
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! Returns a unit number of a unit that can exist and is not connected.  Note
  ! this routine does not magically mark that unit number in use.  In order to
  ! have the unit "used", the source code must OPEN the file.

  ! METHODOLOGY EMPLOYED:
  ! Use Inquire function to find out if proposed unit: exists or is opened.
  ! If not, can be used for a new unit number.

  ! REFERENCES:
  ! Copyright (c) 1994 Unicomp, Inc.  All rights reserved.
  !
  ! Developed at Unicomp, Inc.
  !
  ! Permission to use, copy, modify, and distribute this
  ! software is freely granted, provided that this notice
  ! is preserved.

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  INTEGER UnitNumber  ! Result from scanning currently open files

  ! FUNCTION PARAMETER DEFINITIONS:
  !  IO Status Values:

  INTEGER, PARAMETER :: END_OF_RECORD = -2
  INTEGER, PARAMETER :: END_OF_FILE = -1

  !  Indicate default input and output units:

  INTEGER, PARAMETER :: DEFAULT_INPUT_UNIT = 5
  INTEGER, PARAMETER :: DEFAULT_OUTPUT_UNIT = 6

  !  Indicate number and value of preconnected units

  INTEGER, PARAMETER :: NUMBER_OF_PRECONNECTED_UNITS = 2
  INTEGER, PARAMETER :: PRECONNECTED_UNITS (NUMBER_OF_PRECONNECTED_UNITS) = (/ 5, 6 /)

  !  Largest allowed unit number (or a large number, if none)
  INTEGER, PARAMETER :: MaxUnitNumber = 1000

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  LOGICAL :: exists  ! File exists
  LOGICAL :: opened  ! Unit is open
  INTEGER :: ios     ! return value from Inquire intrinsic

  DO UnitNumber = 1, MaxUnitNumber
    IF (UnitNumber == DEFAULT_INPUT_UNIT .or. &
    UnitNumber == DEFAULT_OUTPUT_UNIT) THEN
    CYCLE
  END IF
  IF (ANY (UnitNumber == PRECONNECTED_UNITS)) THEN
    CYCLE
  END IF
  INQUIRE (UNIT = UnitNumber, EXIST = exists,  OPENED = opened, IOSTAT = ios)
  IF (exists .and. .not. opened .and. ios == 0) THEN
    RETURN      ! result is set in UnitNumber
  END IF
END DO

UnitNumber = -1

END FUNCTION GetNewUnitNumber

!SUBROUTINE ShowFatalError(ErrorMessage,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   September 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine puts ErrorMessage with a Fatal designation on
!          ! designated output files.  Then, the program is aborted.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Calls ShowErrorMessage utility routine.
!          ! Calls AbortEnergyPlus
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!          ! na
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) ErrorMessage
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!    SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)
!        CHARACTER(len=*) Message
!        INTEGER, OPTIONAL :: Unit1
!        INTEGER, OPTIONAL :: Unit2
!    END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!          ! na
!
!  CALL ShowErrorMessage(' **  Fatal  ** '//ErrorMessage,OutUnit1,OutUnit2)
!  CALL AbortEnergyPlus
!
!  RETURN
!
!END SUBROUTINE ShowFatalError
!
!SUBROUTINE ShowSevereError(ErrorMessage,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   September 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine puts ErrorMessage with a Severe designation on
!          ! designated output files.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Calls ShowErrorMessage utility routine.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) ErrorMessage
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!    SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)
!        CHARACTER(len=*) Message
!        INTEGER, OPTIONAL :: Unit1
!        INTEGER, OPTIONAL :: Unit2
!    END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!
!  TotalSevereErrors=TotalSevereErrors+1
!  CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage,OutUnit1,OutUnit2)
!
!  !  Could set a variable here that gets checked at some point?
!
!  RETURN
!
!END SUBROUTINE ShowSevereError

!SUBROUTINE ShowContinueError(Message,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   October 2001
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine displays a 'continued error' message on designated output files.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Calls ShowErrorMessage utility routine.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!          ! na
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) Message
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!    SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)
!        CHARACTER(len=*) Message
!        INTEGER, OPTIONAL :: Unit1
!        INTEGER, OPTIONAL :: Unit2
!    END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!          ! na
!
!  CALL ShowErrorMessage(' **   ~~~   ** '//Message,OutUnit1,OutUnit2)
!
!  RETURN
!
!END SUBROUTINE ShowContinueError
!
!SUBROUTINE ShowMessage(Message,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   September 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine displays a simple message on designated output files.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Calls ShowErrorMessage utility routine.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!          ! na
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) Message
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!  SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)
!    CHARACTER(len=*) Message
!    INTEGER, OPTIONAL :: Unit1
!    INTEGER, OPTIONAL :: Unit2
!  END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!          ! na
!
!  CALL ShowErrorMessage(' ************* '//Message,OutUnit1,OutUnit2)
!
!  RETURN
!
!END SUBROUTINE ShowMessage
!
!SUBROUTINE ShowWarningError(ErrorMessage,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   September 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine puts ErrorMessage with a Warning designation on
!          ! designated output files.
!
!          ! METHODOLOGY EMPLOYED:
!          ! Calls ShowErrorMessage utility routine.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!  USE DataGlobals_HPSimIntegrated !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) ErrorMessage
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!  INTERFACE
!
!  SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)
!    CHARACTER(len=*) Message
!    INTEGER, OPTIONAL :: Unit1
!    INTEGER, OPTIONAL :: Unit2
!  END SUBROUTINE
!
!  END INTERFACE
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!
!  TotalWarningErrors=TotalWarningErrors+1
!  CALL ShowErrorMessage(' ** Warning ** '//ErrorMessage,OutUnit1,OutUnit2)
!
!  RETURN
!
!END SUBROUTINE ShowWarningError

!SUBROUTINE ShowErrorMessage(ErrorMessage,OutUnit1,OutUnit2)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   December 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine displays the error messages on the indicated
!          ! file unit numbers, in addition to the "standard error output"
!          ! unit.
!
!          ! METHODOLOGY EMPLOYED:
!          ! If arguments OutUnit1 and/or OutUnit2 are present the
!          ! error message is written to these as well and the standard one.
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!  USE DataStringGlobals !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12), ONLY: VerString
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS:
!  CHARACTER(len=*) ErrorMessage
!  INTEGER, OPTIONAL :: OutUnit1
!  INTEGER, OPTIONAL :: OutUnit2
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!  CHARACTER(len=*), PARAMETER :: ErrorFormat='(2X,A)'
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!          ! na
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  INTEGER  :: TotalErrors=0        ! used to determine when to open standard error output file.
!  INTEGER  :: StandardErrorOutput
!  INTEGER,EXTERNAL  :: GetNewUnitNumber
!  SAVE     TotalErrors,StandardErrorOutput
!
!  IF (TotalErrors .eq. 0) THEN
!    StandardErrorOutput=GetNewUnitNumber()
!    OPEN(StandardErrorOutput,FILE='eplusout.err')
!    WRITE(StandardErrorOutput,'(A)') 'Program Version,'//TRIM(VerString)
!  ENDIF
!
!  TotalErrors=TotalErrors+1
!  WRITE(StandardErrorOutput,ErrorFormat) TRIM(ErrorMessage)
!  IF (PRESENT(OutUnit1)) THEN
!    WRITE(OutUnit1,ErrorFormat) TRIM(ErrorMessage)
!  ENDIF
!  IF (PRESENT(OutUnit2)) THEN
!    WRITE(OutUnit2,ErrorFormat) TRIM(ErrorMessage)
!  ENDIF
!
!  RETURN
!
!END SUBROUTINE ShowErrorMessage

SUBROUTINE ShowSevereError1(ErrorMessage,OutUnit1)  !RS: Trying this to debug since the Optional values don't seem to be working.

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine puts ErrorMessage with a Severe designation on
  ! designated output files.

  ! METHODOLOGY EMPLOYED:
  ! Calls ShowErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataGlobals_HPSimIntegrated

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) ErrorMessage
  INTEGER:: OutUnit1

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message,Unit1)   !RS: Just commenting out to see what will happen.
      CHARACTER(len=*) Message
      INTEGER:: Unit1
    END SUBROUTINE
  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  TotalSevereErrors=TotalSevereErrors+1

  CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage,OutUnit1) !RS:The optional integers weren't being defined properly, so I took them out for now.
  !CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage)

  !  Could set a variable here that gets checked at some point?

  RETURN

END SUBROUTINE ShowSevereError1

SUBROUTINE ShowSevereError2(ErrorMessage,OutUnit1,OutUnit2)  !RS: Trying this to debug since the Optional values don't seem to be working.

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine puts ErrorMessage with a Severe designation on
  ! designated output files.

  ! METHODOLOGY EMPLOYED:
  ! Calls ShowErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataGlobals_HPSimIntegrated

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) ErrorMessage
  INTEGER :: OutUnit1
  INTEGER :: OutUnit2

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE
    SUBROUTINE ShowErrorMessage(Message,Unit1,Unit2)   !RS: Just commenting out to see what will happen.
      CHARACTER(len=*) Message
      INTEGER :: Unit1
      INTEGER :: Unit2
    END SUBROUTINE
  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  TotalSevereErrors=TotalSevereErrors+1

  CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage,OutUnit1,OutUnit2) !RS:The optional integers weren't being defined properly, so I took them out for now.
  !CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage)

  !  Could set a variable here that gets checked at some point?

  RETURN

END SUBROUTINE ShowSevereError2

!RS: The following subroutines are being added in for the Energy+ side

SUBROUTINE ConvertCasetoLower(InputString,OutputString)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Convert a string to lower case

  ! METHODOLOGY EMPLOYED:
  ! This routine is not dependant upon the ASCII
  ! code.  It works by storing the upper and lower case alphabet.  It
  ! scans the whole input string.  If it finds a character in the lower
  ! case alphabet, it makes an appropriate substitution.


  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataStringGLobals

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: InputString    ! Input string
  CHARACTER(len=*), INTENT(OUT) :: OutputString  ! Output string (in LowerCase)

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER A,B

  OutputString=' '

  DO A=1,LEN_TRIM(InputString)
    B=INDEX(UpperCase,InputString(A:A))
    IF (B .NE. 0) THEN
      OutputString(A:A)=LowerCase(B:B)
    ELSE
      OutputString(A:A)=InputString(A:A)
    ENDIF
  END DO

  RETURN

END SUBROUTINE ConvertCasetoLower

SUBROUTINE ShowContinueErrorTimeStamp(Message,OutUnit1,OutUnit2)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   February 2004
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine displays a 'continued error' timestamp message on designated output files.

  ! METHODOLOGY EMPLOYED:
  ! Calls ShowErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE General,                         ONLY : CreateSysTimeIntervalString
  USE DataEnvironment,                 ONLY : EnvironmentName,CurMnDy
  USE DataGlobals_HPSimIntegrated,                     ONLY : WarmupFlag,DoingSizing
  USE DataInterfaces, ONLY: ShowErrorMessage
  USE SQLiteProcedures, ONLY: UpdateSQLiteErrorRecord, WriteOutputToSQLite

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) Message
  INTEGER, OPTIONAL :: OutUnit1
  INTEGER, OPTIONAL :: OutUnit2

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! see DataInterfaces

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=100) :: cEnvHeader

  IF (WarmupFlag) THEN
    IF (.not. DoingSizing) THEN
      cEnvHeader=' During Warmup, Environment='
    ELSE
      cEnvHeader=' During Warmup & Sizing, Environment='
    ENDIF
  ELSE
    IF (.not. DoingSizing) THEN
      cEnvHeader=' Environment='
    ELSE
      cEnvHeader=' During Sizing, Environment='
    ENDIF
  ENDIF

  IF (Len_Trim(Message) < 50) THEN
    CALL ShowErrorMessage(' **   ~~~   ** '//TRIM(Message)//TRIM(cEnvHeader)//TRIM(EnvironmentName)//', at Simulation time='//  &
    TRIM(CurMnDy)//' '//TRIM(CreateSysTimeIntervalString()),  &
    OutUnit1,OutUnit2)
    IF(WriteOutputToSQLite) THEN
      CALL UpdateSQLiteErrorRecord(TRIM(Message)//TRIM(cEnvHeader)//TRIM(EnvironmentName)//', at Simulation time='//  &
      TRIM(CurMnDy)//' '//TRIM(CreateSysTimeIntervalString()))
    ENDIF

  ELSE
    CALL ShowErrorMessage(' **   ~~~   ** '//TRIM(Message))
    CALL ShowErrorMessage(' **   ~~~   ** '//TRIM(cEnvHeader)//TRIM(EnvironmentName)//', at Simulation time='//  &
    TRIM(CurMnDy)//' '//TRIM(CreateSysTimeIntervalString()),  &
    OutUnit1,OutUnit2)
    IF(WriteOutputToSQLite) THEN
      CALL UpdateSQLiteErrorRecord(TRIM(Message)// &
      TRIM(cEnvHeader)//TRIM(EnvironmentName)//', at Simulation time='//  &
      TRIM(CurMnDy)//' '//TRIM(CreateSysTimeIntervalString()))
    ENDIF
  ENDIF

  RETURN

END SUBROUTINE ShowContinueErrorTimeStamp

SUBROUTINE ShowRecurringSevereErrorAtEnd(Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)


  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Michael J. Witte
  !       DATE WRITTEN   August 2004
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine stores a recurring ErrorMessage with a Severe designation
  ! for output at the end of the simulation with automatic tracking of number
  ! of occurences and optional tracking of associated min, max, and sum values

  ! METHODOLOGY EMPLOYED:
  ! Calls StoreRecurringErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataStringGlobals
  USE DataErrorTracking

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*)              :: Message     ! Message automatically written to "error file" at end of simulation
  INTEGER, INTENT(INOUT)        :: MsgIndex    ! Recurring message index, if zero, next available index is assigned
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMaxOf ! Track and report the max of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMinOf ! Track and report the min of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportSumOf ! Track and report the sum of the values passed to this argument
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMaxUnits ! optional char string (<=15 length) of units for max value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMinUnits ! optional char string (<=15 length) of units for min value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportSumUnits ! optional char string (<=15 length) of units for sum value

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE
    SUBROUTINE StoreRecurringErrorMessage(ErrorMessage,ErrorMsgIndex,ErrorReportMaxOf,ErrorReportMinOf,ErrorReportSumOf,  &
      ErrorReportMaxUnits,ErrorReportMinUnits,ErrorReportSumUnits)
      USE DataPrecisionGlobals
      !  Use for recurring "warning" error messages shown once at end of simulation
      !  with count of occurences and optional max, min, sum
      CHARACTER(len=*) ErrorMessage    ! Message automatically written to "error file" at end of simulation
      INTEGER, INTENT(INOUT)        :: ErrorMsgIndex    ! Recurring message index, if zero, next available index is assigned
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMaxOf ! Track and report the max of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMinOf ! Track and report the min of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportSumOf ! Track and report the sum of the values passed to this argument
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMaxUnits ! Units for "max" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMinUnits ! Units for "min" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportSumUnits ! Units for "sum" reporting
    END SUBROUTINE
  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Loop

  DO Loop=1,SearchCounts
    IF (INDEX(Message,TRIM(MessageSearch(Loop))) > 0) MatchCounts(Loop)=MatchCounts(Loop)+1
  ENDDO

  TotalSevereErrors=TotalSevereErrors+1
  CALL StoreRecurringErrorMessage(' ** Severe  ** '//Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)

  RETURN

END SUBROUTINE ShowRecurringSevereErrorAtEnd

SUBROUTINE ShowRecurringWarningErrorAtEnd(Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Michael J. Witte
  !       DATE WRITTEN   August 2004
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine stores a recurring ErrorMessage with a Warning designation
  ! for output at the end of the simulation with automatic tracking of number
  ! of occurences and optional tracking of associated min, max, and sum values

  ! METHODOLOGY EMPLOYED:
  ! Calls StoreRecurringErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataStringGlobals
  USE DataErrorTracking

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*)              :: Message     ! Message automatically written to "error file" at end of simulation
  INTEGER, INTENT(INOUT)        :: MsgIndex    ! Recurring message index, if zero, next available index is assigned
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMaxOf ! Track and report the max of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMinOf ! Track and report the min of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportSumOf ! Track and report the sum of the values passed to this argument
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMaxUnits ! optional char string (<=15 length) of units for max value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMinUnits ! optional char string (<=15 length) of units for min value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportSumUnits ! optional char string (<=15 length) of units for sum value

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE
    SUBROUTINE StoreRecurringErrorMessage(ErrorMessage,ErrorMsgIndex,ErrorReportMaxOf,ErrorReportMinOf,ErrorReportSumOf,  &
      ErrorReportMaxUnits,ErrorReportMinUnits,ErrorReportSumUnits)
      USE DataPrecisionGlobals
      !  Use for recurring "warning" error messages shown once at end of simulation
      !  with count of occurences and optional max, min, sum
      CHARACTER(len=*) ErrorMessage    ! Message automatically written to "error file" at end of simulation
      INTEGER, INTENT(INOUT)        :: ErrorMsgIndex    ! Recurring message index, if zero, next available index is assigned
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMaxOf ! Track and report the max of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMinOf ! Track and report the min of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportSumOf ! Track and report the sum of the values passed to this argument
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMaxUnits ! Units for "max" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMinUnits ! Units for "min" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportSumUnits ! Units for "sum" reporting
    END SUBROUTINE
  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Loop

  DO Loop=1,SearchCounts
    IF (INDEX(Message,TRIM(MessageSearch(Loop))) > 0) MatchCounts(Loop)=MatchCounts(Loop)+1
  ENDDO

  TotalWarningErrors=TotalWarningErrors+1
  CALL StoreRecurringErrorMessage(' ** Warning ** '//Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)

  RETURN

END SUBROUTINE ShowRecurringWarningErrorAtEnd

SUBROUTINE ShowSevereMessage(ErrorMessage,OutUnit1,OutUnit2)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 2009
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine puts ErrorMessage with a Severe designation on
  ! designated output files.
  ! But does not bump the error count so can be used in conjunction with recurring
  ! error calls.

  ! METHODOLOGY EMPLOYED:
  ! Calls ShowErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataStringGlobals
  USE DataErrorTracking
  USE DataInterfaces, ONLY: ShowErrorMessage
  USE SQLiteProcedures, ONLY: CreateSQLiteErrorRecord, WriteOutputToSQLite

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) ErrorMessage
  INTEGER, OPTIONAL :: OutUnit1
  INTEGER, OPTIONAL :: OutUnit2

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! see DataInterfaces

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Loop

  DO Loop=1,SearchCounts
    IF (INDEX(ErrorMessage,TRIM(MessageSearch(Loop))) > 0) MatchCounts(Loop)=MatchCounts(Loop)+1
  ENDDO

  CALL ShowErrorMessage(' ** Severe  ** '//ErrorMessage,OutUnit1,OutUnit2)
  LastSevereError=ErrorMessage

  !  Could set a variable here that gets checked at some point?

  IF(WriteOutputToSQLite) THEN
    CALL CreateSQLiteErrorRecord(1,1,ErrorMessage,0)
  ENDIF
  RETURN

END SUBROUTINE ShowSevereMessage

SUBROUTINE ShowWarningMessage(ErrorMessage,OutUnit1,OutUnit2)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 2009
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine puts ErrorMessage with a Warning designation on
  ! designated output files.
  ! But does not bump the error count so can be used in conjunction with recurring
  ! error calls.

  ! METHODOLOGY EMPLOYED:
  ! Calls ShowErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataStringGlobals
  USE DataErrorTracking
  USE DataInterfaces, ONLY: ShowErrorMessage
  USE SQLiteProcedures, ONLY: CreateSQLiteErrorRecord, WriteOutputToSQLite

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) ErrorMessage
  INTEGER, OPTIONAL :: OutUnit1
  INTEGER, OPTIONAL :: OutUnit2

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! see DataInterfaces

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Loop

  DO Loop=1,SearchCounts
    IF (INDEX(ErrorMessage,TRIM(MessageSearch(Loop))) > 0) MatchCounts(Loop)=MatchCounts(Loop)+1
  ENDDO

  CALL ShowErrorMessage(' ** Warning ** '//ErrorMessage,OutUnit1,OutUnit2)
  IF(WriteOutputToSQLite) THEN
    CALL CreateSQLiteErrorRecord(1,0,ErrorMessage,0)
  ENDIF

  RETURN

END SUBROUTINE ShowWarningMessage

SUBROUTINE ShowRecurringContinueErrorAtEnd(Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Michael J. Witte
  !       DATE WRITTEN   August 2004
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine stores a recurring ErrorMessage with a continue designation
  ! for output at the end of the simulation with automatic tracking of number
  ! of occurences and optional tracking of associated min, max, and sum values

  ! METHODOLOGY EMPLOYED:
  ! Calls StoreRecurringErrorMessage utility routine.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataStringGlobals
  USE DataErrorTracking

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*)              :: Message     ! Message automatically written to "error file" at end of simulation
  INTEGER, INTENT(INOUT)        :: MsgIndex    ! Recurring message index, if zero, next available index is assigned
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMaxOf ! Track and report the max of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportMinOf ! Track and report the min of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ReportSumOf ! Track and report the sum of the values passed to this argument
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMaxUnits ! optional char string (<=15 length) of units for max value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportMinUnits ! optional char string (<=15 length) of units for min value
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ReportSumUnits ! optional char string (<=15 length) of units for sum value

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  INTERFACE
    SUBROUTINE StoreRecurringErrorMessage(ErrorMessage,ErrorMsgIndex,ErrorReportMaxOf,ErrorReportMinOf,ErrorReportSumOf,  &
      ErrorReportMaxUnits,ErrorReportMinUnits,ErrorReportSumUnits)
      USE DataPrecisionGlobals
      !  Use for recurring "warning" error messages shown once at end of simulation
      !  with count of occurences and optional max, min, sum
      CHARACTER(len=*) ErrorMessage    ! Message automatically written to "error file" at end of simulation
      INTEGER, INTENT(INOUT)        :: ErrorMsgIndex    ! Recurring message index, if zero, next available index is assigned
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMaxOf ! Track and report the max of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMinOf ! Track and report the min of the values passed to this argument
      REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportSumOf ! Track and report the sum of the values passed to this argument
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMaxUnits ! Units for "max" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMinUnits ! Units for "min" reporting
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportSumUnits ! Units for "sum" reporting
    END SUBROUTINE
  END INTERFACE

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Loop

  DO Loop=1,SearchCounts
    IF (INDEX(Message,TRIM(MessageSearch(Loop))) > 0) MatchCounts(Loop)=MatchCounts(Loop)+1
  ENDDO

  CALL StoreRecurringErrorMessage(' **   ~~~   ** '//Message,MsgIndex,ReportMaxOf,ReportMinOf,ReportSumOf,  &
  ReportMaxUnits,ReportMinUnits,ReportSumUnits)

  RETURN

END SUBROUTINE ShowRecurringContinueErrorAtEnd

SUBROUTINE StoreRecurringErrorMessage(ErrorMessage,ErrorMsgIndex,ErrorReportMaxOf,ErrorReportMinOf,ErrorReportSumOf,  &
  ErrorReportMaxUnits,ErrorReportMinUnits,ErrorReportSumUnits)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Michael J. Witte
  !       DATE WRITTEN   August 2004
  !       MODIFIED       September 2005;LKL;Added Units
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine stores a recurring ErrorMessage with
  ! for output at the end of the simulation with automatic tracking of number
  ! of occurences and optional tracking of associated min, max, and sum values

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataPrecisionGlobals
  USE DataStringGlobals
  USE DataErrorTracking
  USE DataGlobals_HPSimIntegrated,         ONLY : WarmupFlag,DoingSizing

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) ErrorMessage    ! Message automatically written to "error file" at end of simulation
  INTEGER, INTENT(INOUT)        :: ErrorMsgIndex    ! Recurring message index, if zero, next available index is assigned
  REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMaxOf ! Track and report the max of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportMinOf ! Track and report the min of the values passed to this argument
  REAL(r64),    INTENT(IN), OPTIONAL :: ErrorReportSumOf ! Track and report the sum of the values passed to this argument
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMaxUnits ! Units for "max" reporting
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportMinUnits ! Units for "min" reporting
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ErrorReportSumUnits ! Units for "sum" reporting

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  TYPE (RecurringErrorData),   ALLOCATABLE, DIMENSION(:) :: TempRecurringErrors

  ! If Index is zero, then assign next available index and reallocate array
  IF (ErrorMsgIndex == 0) THEN
    NumRecurringErrors = NumRecurringErrors + 1
    ErrorMsgIndex = NumRecurringErrors
    IF (NumRecurringErrors == 1) THEN
      ALLOCATE(RecurringErrors(NumRecurringErrors))
    ELSEIF (NumRecurringErrors > 1) THEN
      ALLOCATE(TempRecurringErrors(NumRecurringErrors))
      TempRecurringErrors(1:NumRecurringErrors-1)=RecurringErrors(1:NumRecurringErrors-1)
      DEALLOCATE(RecurringErrors)
      ALLOCATE(RecurringErrors(NumRecurringErrors))
      RecurringErrors = TempRecurringErrors
      DEALLOCATE(TempRecurringErrors)
    ENDIF
    ! The message string only needs to be stored once when a new recurring message is created
    RecurringErrors(ErrorMsgIndex)%Message = TRIM(ErrorMessage)
    RecurringErrors(ErrorMsgIndex)%Count   = 1
    IF (WarmupFlag) RecurringErrors(ErrorMsgIndex)%WarmupCount   = 1
    IF (DoingSizing) RecurringErrors(ErrorMsgIndex)%SizingCount   = 1


    ! For max, min, and sum values, store the current value when a new recurring message is created
    IF (PRESENT(ErrorReportMaxOf)) THEN
      RecurringErrors(ErrorMsgIndex)%MaxValue = ErrorReportMaxOf
      RecurringErrors(ErrorMsgIndex)%ReportMax = .TRUE.
      IF (PRESENT(ErrorReportMaxUnits)) THEN
        RecurringErrors(ErrorMsgIndex)%MaxUnits=ErrorReportMaxUnits
      ENDIF
    ENDIF
    IF (PRESENT(ErrorReportMinOf)) THEN
      RecurringErrors(ErrorMsgIndex)%MinValue = ErrorReportMinOf
      RecurringErrors(ErrorMsgIndex)%ReportMin = .TRUE.
      IF (PRESENT(ErrorReportMinUnits)) THEN
        RecurringErrors(ErrorMsgIndex)%MinUnits=ErrorReportMinUnits
      ENDIF
    ENDIF
    IF (PRESENT(ErrorReportSumOf)) THEN
      RecurringErrors(ErrorMsgIndex)%SumValue = ErrorReportSumOf
      RecurringErrors(ErrorMsgIndex)%ReportSum = .TRUE.
      IF (PRESENT(ErrorReportSumUnits)) THEN
        RecurringErrors(ErrorMsgIndex)%SumUnits=ErrorReportSumUnits
      ENDIF
    ENDIF

  ELSEIF (ErrorMsgIndex > 0) THEN
    ! Do stats and store
    RecurringErrors(ErrorMsgIndex)%Count = RecurringErrors(ErrorMsgIndex)%Count + 1
    IF (WarmupFlag) RecurringErrors(ErrorMsgIndex)%WarmupCount = RecurringErrors(ErrorMsgIndex)%WarmupCount + 1
    IF (DoingSizing) RecurringErrors(ErrorMsgIndex)%SizingCount = RecurringErrors(ErrorMsgIndex)%SizingCount + 1

    IF (PRESENT(ErrorReportMaxOf)) THEN
      RecurringErrors(ErrorMsgIndex)%MaxValue = MAX(ErrorReportMaxOf,RecurringErrors(ErrorMsgIndex)%MaxValue)
      RecurringErrors(ErrorMsgIndex)%ReportMax = .TRUE.
    ENDIF
    IF (PRESENT(ErrorReportMinOf)) THEN
      RecurringErrors(ErrorMsgIndex)%MinValue = MIN(ErrorReportMinOf,RecurringErrors(ErrorMsgIndex)%MinValue)
      RecurringErrors(ErrorMsgIndex)%ReportMin = .TRUE.
    ENDIF
    IF (PRESENT(ErrorReportSumOf)) THEN
      RecurringErrors(ErrorMsgIndex)%SumValue = ErrorReportSumOf + RecurringErrors(ErrorMsgIndex)%SumValue
      RecurringErrors(ErrorMsgIndex)%ReportSum = .TRUE.
    ENDIF
  ELSE
    ! If ErrorMsgIndex < 0, then do nothing
  ENDIF

  RETURN

END SUBROUTINE StoreRecurringErrorMessage

SUBROUTINE ConvertCasetoUpper(InputString,OutputString)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Convert a string to upper case

  ! METHODOLOGY EMPLOYED:
  ! This routine is not dependant upon the ASCII
  ! code.  It works by storing the upper and lower case alphabet.  It
  ! scans the whole input string.  If it finds a character in the lower
  ! case alphabet, it makes an appropriate substitution.


  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataStringGlobals

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: InputString    ! Input string
  CHARACTER(len=*), INTENT(OUT) :: OutputString  ! Output string (in UpperCase)

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER A,B

  OutputString=' '

  DO A=1,LEN_TRIM(InputString)
    B=INDEX(LowerCase,InputString(A:A))
    IF (B .NE. 0) THEN
      OutputString(A:A)=UpperCase(B:B)
    ELSE
      OutputString(A:A)=InputString(A:A)
    ENDIF
  END DO

  RETURN

END SUBROUTINE ConvertCasetoUpper

INTEGER FUNCTION FindNonSpace(String)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function finds the first non-space character in the passed string
  ! and returns that position as the result to the calling program.

  ! METHODOLOGY EMPLOYED:
  ! Scan string for character not equal to blank.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: String  ! String to be scanned

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER I,ILEN

  FindNonSpace=0
  ILEN=LEN_TRIM(String)
  DO I=1,ILEN
    IF (String(I:I) .NE. ' ') THEN
      FindNonSpace=I
      EXIT
    END IF
  END DO

  RETURN

END FUNCTION FindNonSpace

FUNCTION FindUnitNumber (FileName) RESULT (UnitNumber)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997, adapted from reference
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! Returns a unit number for the file name that is either opened or exists.

  ! METHODOLOGY EMPLOYED:
  ! Use Inquire function to find out if proposed unit: exists or is opened.
  ! If not, can be used for a new unit number.

  ! REFERENCES:
  ! Copyright (c) 1994 Unicomp, Inc.  All rights reserved.
  !
  ! Developed at Unicomp, Inc.
  !
  ! Permission to use, copy, modify, and distribute this
  ! software is freely granted, provided that this notice
  ! is preserved.

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*) FileName  ! File name to be searched.
  INTEGER UnitNumber         ! Unit number that should be used

  ! FUNCTION PARAMETER DEFINITIONS:
  !  Largest allowed unit number (or a large number, if none)
  INTEGER, PARAMETER :: MaxUnitNumber = 1000

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  CHARACTER(Len=255) TestFileName       ! File name returned from opened file
  INTEGER TestFileLength                ! Length from INQUIRE intrinsic
  INTEGER,EXTERNAL :: GetNewUnitNumber  ! Function to call if file not opened
  LOGICAL :: exists                     ! True if file already exists
  LOGICAL :: opened                     ! True if file is open
  INTEGER Pos                           ! Position pointer
  INTEGER FileNameLength                ! Length of requested file
  INTEGER :: ios                        ! Status indicator from INQUIRE intrinsic

  INQUIRE (FILE=FileName, EXIST = exists,  OPENED = opened, IOSTAT = ios)
  IF (.not. OPENED) THEN
    UnitNumber=GetNewUnitNumber()
    OPEN(UNIT=UnitNumber,FILE=FileName,POSITION='APPEND',iostat=ios)
    IF (ios /= 0) THEN
      CALL DisplayString('FindUnitNumber: Could not open file "'//trim(FileName)//'" for append.')
    ENDIF
  ELSE
    FileNameLength=LEN_TRIM(FileName)
    DO UnitNumber=1,MaxUnitNumber
      INQUIRE(UNIT=UnitNumber,NAME=TestFileName,OPENED=opened)
      !  Powerstation returns just file name
      !  DVF (Digital Fortran) returns whole path
      TestFileLength=LEN_TRIM(TestFileName)
      Pos=INDEX(TestFileName,FileName)
      IF (Pos .ne. 0) THEN
        !  Must be the last part of the file
        IF (Pos+FileNameLength-1 .eq. TestFileLength) EXIT
      ENDIF
    END DO
  ENDIF

  RETURN

END FUNCTION FindUnitNumber


!     NOTICE
!
!     Copyright © 1996-2003 The Board of Trustees of the University of Illinois
!     and The Regents of the University of California through Ernest Orlando Lawrence
!     Berkeley National Laboratory.  All rights reserved.
!
!     Portions of the EnergyPlus software package have been developed and copyrighted
!     by other individuals, companies and institutions.  These portions have been
!     incorporated into the EnergyPlus software package under license.   For a complete
!     list of contributors, see "Notice" located in EnergyPlus.f90.
!
!     NOTICE: The U.S. Government is granted for itself and others acting on its
!     behalf a paid-up, nonexclusive, irrevocable, worldwide license in this data to
!     reproduce, prepare derivative works, and perform publicly and display publicly.
!     Beginning five (5) years after permission to assert copyright is granted,
!     subject to two possible five year renewals, the U.S. Government is granted for
!     itself and others acting on its behalf a paid-up, non-exclusive, irrevocable
!     worldwide license in this data to reproduce, prepare derivative works,
!     distribute copies to the public, perform publicly and display publicly, and to
!     permit others to do so.
!
!     TRADEMARKS: EnergyPlus is a trademark of the US Department of Energy.
!
