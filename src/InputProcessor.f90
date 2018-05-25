MODULE InputProcessor
  ! Module containing the input processor routines

  ! MODULE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   August 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS MODULE:
  ! To provide the capabilities of reading the input data dictionary
  ! and input file and supplying the simulation routines with the data
  ! contained therein.

  ! METHODOLOGY EMPLOYED:
  !

  ! REFERENCES:
  ! The input syntax is designed to allow for future flexibility without
  ! necessitating massive (or any) changes to this code.  Two files are
  ! used as key elements: (1) the input data dictionary will specify the
  ! sections and objects that will be allowed in the actual simulation
  ! input file and (2) the simulation input data file will be processed
  ! with the data therein being supplied to the actual simulation routines.



  ! OTHER NOTES:
  !
  !

  ! USE STATEMENTS:
  ! Use statements for data only modules
  USE DataPrecisionGlobals
  USE DataStringGlobals
  USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength,AutoCalculate,rTinyValue, DisplayAllWarnings,DisplayUnusedObjects,  &
  CacheIPErrorFile,DoingInputProcessing
  USE DataSizing, ONLY: AutoSize
  USE DataIPShortCuts
  USE DataSystemVariables, ONLY: SortedIDD, iASCII_CR, iUnicode_end

  ! Use statements for access to subroutines in other modules

  IMPLICIT NONE         ! Enforce explicit typing of all variables

  PRIVATE

  !MODULE PARAMETER DEFINITIONS
  INTEGER, PARAMETER         :: ObjectDefAllocInc=100     ! Starting number of Objects allowed in IDD as well as the increment
  ! when max is reached
  INTEGER, PARAMETER         :: ANArgsDefAllocInc=500     ! The increment when max total args is reached
  INTEGER, PARAMETER         :: SectionDefAllocInc=20     ! Starting number of Sections allowed in IDD as well as the increment
  ! when max is reached
  INTEGER, PARAMETER         :: SectionsIDFAllocInc=20    ! Initial number of Sections allowed in IDF as well as the increment
  ! when max is reached
  INTEGER, PARAMETER         :: ObjectsIDFAllocInc=500    ! Initial number of Objects allowed in IDF as well as the increment
  ! when max is reached
  INTEGER, PARAMETER         :: MaxObjectNameLength=MaxNameLength    ! Maximum number of characters in an Object Name
  INTEGER, PARAMETER         :: MaxSectionNameLength=MaxNameLength   ! Maximum number of characters in a Section Name
  INTEGER, PARAMETER         :: MaxAlphaArgLength=MaxNameLength  ! Maximum number of characters in an Alpha Argument
  INTEGER, PARAMETER         :: MaxInputLineLength=500    ! Maximum number of characters in an input line (in.idf, energy+.idd)
  INTEGER, PARAMETER         :: MaxFieldNameLength=140    ! Maximum number of characters in a field name string
  CHARACTER(len=1), PARAMETER :: Blank=' '
  CHARACTER(len=*), PARAMETER :: AlphaNum='ANan'     ! Valid indicators for Alpha or Numeric fields (A or N)
  CHARACTER(len=*), PARAMETER :: fmta='(A)'
  REAL(r64), PARAMETER :: DefAutoSizeValue=AutoSize
  REAL(r64), PARAMETER :: DefAutoCalculateValue=AutoCalculate

  ! DERIVED TYPE DEFINITIONS
  TYPE RangeCheckDef
    LOGICAL :: MinMaxChk                            =.false.   ! true when Min/Max has been added
    INTEGER :: FieldNumber                          =0         ! which field number this is
    CHARACTER(len=MaxFieldNameLength) :: FieldName =Blank       ! Name of the field
    CHARACTER(len=32), DIMENSION(2) :: MinMaxString =Blank       ! appropriate Min/Max Strings
    REAL(r64), DIMENSION(2) :: MinMaxValue          =0.0       ! appropriate Min/Max Values
    INTEGER, DIMENSION(2) :: WhichMinMax            =0         !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
    LOGICAL :: DefaultChk                           =.false.   ! true when default has been entered
    REAL(r64)  :: Default                           =0.0       ! Default value
    LOGICAL :: DefAutoSize                          =.false.   ! Default value is "autosize"
    LOGICAL :: AutoSizable                          =.false.   ! True if this field can be autosized
    REAL(r64)  :: AutoSizeValue                     =0.0       ! Value to return for autosize field
    LOGICAL :: DefAutoCalculate                     =.false.   ! Default value is "autocalculate"
    LOGICAL :: AutoCalculatable                     =.false.   ! True if this field can be autocalculated
    REAL(r64)  :: AutoCalculateValue                =0.0       ! Value to return for autocalculate field
  END TYPE

  TYPE ObjectsDefinition
    CHARACTER(len=MaxObjectNameLength) :: Name =Blank ! Name of the Object
    INTEGER :: NumParams                       =0   ! Number of parameters to be processed for each object
    INTEGER :: NumAlpha                        =0   ! Number of Alpha elements in the object
    INTEGER :: NumNumeric                      =0   ! Number of Numeric elements in the object
    INTEGER :: MinNumFields                    =0   ! Minimum number of fields to be passed to the Get routines
    LOGICAL :: NameAlpha1                  =.false. ! True if the first alpha appears to "name" the object for error messages
    LOGICAL :: UniqueObject                =.false. ! True if this object has been designated \unique-object
    LOGICAL :: RequiredObject              =.false. ! True if this object has been designated \required-object
    LOGICAL :: ExtensibleObject            =.false. ! True if this object has been designated \extensible
    INTEGER :: ExtensibleNum                   =0   ! how many fields to extend
    INTEGER :: LastExtendAlpha                 =0   ! Count for extended alpha fields
    INTEGER :: LastExtendNum                   =0   ! Count for extended numeric fields
    INTEGER :: ObsPtr                          =0   ! If > 0, object is obsolete and this is the
    ! Pointer to ObsoleteObjectRepNames Array for replacement object
    LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: AlphaorNumeric ! Positionally, whether the argument
    ! is alpha (true) or numeric (false)
    LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: ReqField ! True for required fields
    LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: AlphRetainCase ! true if retaincase is set for this field (alpha fields only)
    CHARACTER(len=MaxFieldNameLength),  &
    ALLOCATABLE, DIMENSION(:) :: AlphFieldChks ! Field names for alphas
    CHARACTER(len=MaxNameLength),  &
    ALLOCATABLE, DIMENSION(:) :: AlphFieldDefs ! Defaults for alphas
    TYPE(RangeCheckDef), ALLOCATABLE, DIMENSION(:) :: NumRangeChks  ! Used to range check and default numeric fields
    INTEGER :: NumFound                        =0   ! Number of this object found in IDF
  END TYPE

  TYPE SectionsDefinition
    CHARACTER(len=MaxSectionNameLength) :: Name =Blank ! Name of the Section
    INTEGER :: NumFound                         =0   ! Number of this object found in IDF
  END TYPE

  TYPE FileSectionsDefinition
    CHARACTER(len=MaxSectionNameLength) :: Name =Blank ! Name of this section
    INTEGER :: FirstRecord                      =0   ! Record number of first object in section
    INTEGER :: FirstLineNo                      =0   ! Record number of first object in section
    INTEGER :: LastRecord                       =0   ! Record number of last object in section
  END TYPE

  TYPE LineDefinition      ! Will be saved for each "object" input
    ! The arrays (Alphas, Numbers) will be dimensioned to be
    ! the size expected from the definition.
    CHARACTER(len=MaxObjectNameLength) :: Name  =Blank ! Object name for this record
    INTEGER :: NumAlphas                        =0   ! Number of alphas on this record
    INTEGER :: NumNumbers                       =0   ! Number of numbers on this record
    INTEGER :: ObjectDefPtr                     =0   ! Which Object Def is this
    CHARACTER(len=MaxAlphaArgLength), ALLOCATABLE, DIMENSION(:) :: Alphas ! Storage for the alphas
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: AlphBlank  ! Set to true if this field was blank on input
    REAL(r64), ALLOCATABLE, DIMENSION(:) :: Numbers       ! Storage for the numbers
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: NumBlank   ! Set to true if this field was blank on input
  END TYPE

  TYPE SecretObjects
    CHARACTER(len=MaxObjectNameLength) :: OldName = Blank    ! Old Object Name
    CHARACTER(len=MaxObjectNameLength) :: NewName = Blank    ! New Object Name if applicable
    LOGICAL                            :: Deleted =.false. ! true if this (old name) was deleted
    LOGICAL                            :: Used    =.false. ! true when used (and reported) in this input file
    LOGICAL                            :: Transitioned =.false. ! true if old name will be transitioned to new object within IP
    LOGICAL                            :: TransitionDefer =.false. ! true if old name will be transitioned to new object within IP
  END TYPE

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! MODULE VARIABLE DECLARATIONS:

  !Integer Variables for the Module
  INTEGER :: NumObjectDefs       =0 ! Count of number of object definitions found in the IDD
  INTEGER :: NumSectionDefs      =0 ! Count of number of section defintions found in the IDD
  INTEGER :: MaxObjectDefs       =0 ! Current "max" object defs (IDD), when reached will be reallocated and new Max set
  INTEGER :: MaxSectionDefs      =0 ! Current "max" section defs (IDD), when reached will be reallocated and new Max set
  INTEGER :: IDDFile             =0 ! Unit number for reading IDD (Energy+.idd)
  INTEGER :: IDFFile             =0 ! Unit number for reading IDF (in.idf)
  INTEGER :: NumLines            =0 ! Count of number of lines in IDF
  INTEGER :: MaxIDFRecords       =0 ! Current "max" IDF records (lines), when reached will be reallocated and new Max set
  INTEGER :: NumIDFRecords       =0 ! Count of number of IDF records
  INTEGER :: MaxIDFSections      =0 ! Current "max" IDF sections (lines), when reached will be reallocated and new Max set
  INTEGER :: NumIDFSections      =0 ! Count of number of IDF records
  INTEGER, EXTERNAL :: GetNewUnitNumber  ! External  function to "get" a unit number
  INTEGER :: EchoInputFile       =0 ! Unit number of the file echoing the IDD and input records (eplusout.audit)
  INTEGER :: InputLineLength     =0 ! Actual input line length or position of comment character
  INTEGER :: MaxAlphaArgsFound   =0 ! Count of max alpha args found in the IDD
  INTEGER :: MaxNumericArgsFound =0 ! Count of max numeric args found in the IDD
  INTEGER :: MaxAlphaIDFArgsFound   =0 ! Count of max alpha args found in the IDF
  INTEGER :: MaxNumericIDFArgsFound =0 ! Count of max numeric args found in the IDF
  INTEGER :: MaxAlphaIDFDefArgsFound   =0 ! Count of max alpha args found in the IDF
  INTEGER :: MaxNumericIDFDefArgsFound =0 ! Count of max numeric args found in the IDF
  INTEGER,PUBLIC :: NumOutOfRangeErrorsFound=0  ! Count of number of "out of range" errors found
  INTEGER,PUBLIC :: NumBlankReqFieldFound=0 ! Count of number of blank required field errors found
  INTEGER,PUBLIC :: NumMiscErrorsFound  =0  ! Count of other errors found
  INTEGER :: MinimumNumberOfFields=0 ! When ReadLine discovers a "minimum" number of fields for an object, this variable is set
  INTEGER :: NumObsoleteObjects=0    ! Number of \obsolete objects
  INTEGER :: TotalAuditErrors=0      ! Counting some warnings that go onto only the audit file
  INTEGER :: NumSecretObjects=0      ! Number of objects in "Secret Mode"
  LOGICAL :: ProcessingIDD=.false.   ! True when processing IDD, false when processing IDF

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  !Real Variables for Module
  !na

  !Character Variables for Module
  CHARACTER(len=MaxInputLineLength+50) :: InputLine=Blank        ! Each line can be up to MaxInputLineLength characters long
  CHARACTER(len=MaxSectionNameLength), ALLOCATABLE, DIMENSION(:) :: ListofSections
  CHARACTER(len=MaxObjectNameLength),  ALLOCATABLE, DIMENSION(:) :: ListofObjects
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iListOfObjects
  INTEGER,  ALLOCATABLE, DIMENSION(:) :: ObjectGotCount
  INTEGER,  ALLOCATABLE, DIMENSION(:) :: ObjectStartRecord
  CHARACTER(len=MaxObjectNameLength) :: CurrentFieldName=Blank   ! Current Field Name (IDD)
  CHARACTER(len=MaxObjectNameLength), ALLOCATABLE, DIMENSION(:) ::   &
  ObsoleteObjectsRepNames  ! Array of Replacement names for Obsolete objects
  CHARACTER(len=MaxObjectNameLength) :: ReplacementName=Blank

  !Logical Variables for Module
  LOGICAL,PUBLIC :: OverallErrorFlag =.false.     ! If errors found during parse of IDF, will fatal at end
  LOGICAL :: EchoInputLine=.true.          ! Usually True, if the IDD is backspaced, then is set to false, then back to true
  LOGICAL :: ReportRangeCheckErrors=.true. ! Module level reporting logical, can be turned off from outside the module (and then
  ! must be turned back on.
  LOGICAL :: FieldSet=.false.              ! Set to true when ReadInputLine has just scanned a "field"
  LOGICAL :: RequiredField=.false.         ! Set to true when ReadInputLine has determined that this field is required
  LOGICAL :: RetainCaseFlag=.false.        ! Set to true when ReadInputLine has determined that this field should retain case
  LOGICAL :: ObsoleteObject=.false.        ! Set to true when ReadInputLine has an obsolete object
  LOGICAL :: RequiredObject=.false.        ! Set to true when ReadInputLine has a required object
  LOGICAL :: UniqueObject=.false.          ! Set to true when ReadInputLine has a unique object
  LOGICAL :: ExtensibleObject=.false.      ! Set to true when ReadInputLine has an extensible object
  LOGICAL :: StripCR=.false.               ! If true, strip last character (<cr> off each schedule:file line)
  INTEGER :: ExtensibleNumFields=0         ! set to number when ReadInputLine has an extensible object
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: IDFRecordsGotten  ! Denotes that this record has been "gotten" from the IDF

  !Derived Types Variables

  TYPE (ObjectsDefinition), ALLOCATABLE, DIMENSION(:)      :: ObjectDef   ! Contains all the Valid Objects on the IDD
  TYPE (SectionsDefinition), ALLOCATABLE, DIMENSION(:)     :: SectionDef ! Contains all the Valid Sections on the IDD
  TYPE (FileSectionsDefinition), ALLOCATABLE, DIMENSION(:) :: SectionsonFile  ! lists the sections on file (IDF)
  TYPE (LineDefinition), SAVE :: LineItem                                          ! Description of current record
  TYPE (LineDefinition), ALLOCATABLE, DIMENSION(:)         :: IDFRecords     ! All the objects read from the IDF
  TYPE (SecretObjects), ALLOCATABLE, DIMENSION(:)          :: RepObjects         ! Secret Objects that could replace old ones

  !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)
  TYPE (ObjectsDefinition), ALLOCATABLE, DIMENSION(:)      :: ObjectDef2   ! Contains all the Valid Objects on the IDD
  TYPE (SectionsDefinition), ALLOCATABLE, DIMENSION(:)     :: SectionDef2 ! Contains all the Valid Sections on the IDD
  TYPE (FileSectionsDefinition), ALLOCATABLE, DIMENSION(:) :: SectionsonFile2  ! lists the sections on file (IDF)
  TYPE (LineDefinition):: LineItem2                        ! Description of current record
  TYPE (LineDefinition), ALLOCATABLE, DIMENSION(:)         :: IDFRecords2     ! All the objects read from the IDF
  TYPE (SecretObjects), ALLOCATABLE, DIMENSION(:)          :: RepObjects2         ! Secret Objects that could replace old ones

  CHARACTER(len=MaxSectionNameLength), ALLOCATABLE, DIMENSION(:) :: ListofSections2
  CHARACTER(len=MaxObjectNameLength),  ALLOCATABLE, DIMENSION(:) :: ListofObjects2
  CHARACTER(len=MaxObjectNameLength) :: CurrentFieldName2   ! Current Field Name (IDD)
  CHARACTER(len=MaxObjectNameLength), ALLOCATABLE, DIMENSION(:) ::   &
  ObsoleteObjectsRepNames2  ! Array of Replacement names for Obsolete objects

  INTEGER NumObjectDefs2       ! Count of number of object definitions found in the IDD
  INTEGER NumSectionDefs2      ! Count of number of section defintions found in the IDD
  INTEGER NumIDFRecords2       ! Current "max" IDF records (lines), when reached will be reallocated and new Max set

  PUBLIC  ProcessInput

  PUBLIC  GetNumSectionsFound
  PRIVATE GetNumSectionsinInput
  PUBLIC  FindIteminList
  PUBLIC  FindIteminSortedList
  PUBLIC  FindItem
  PUBLIC  SameString
  PUBLIC  MakeUPPERCase
  PUBLIC  ProcessNumber
  PUBLIC  RangeCheck
  PUBLIC  VerifyName

  PUBLIC  GetNumObjectsFound
  PUBLIC  GetObjectItem
  PUBLIC  GetObjectItemNum

  PUBLIC DeallocateArrays !RS: Debugging: Added in from InputProcessor_HPSim

  PUBLIC GetObjectItem2   !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)
  PUBLIC ProcessInput2    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (10/6/14)


  PRIVATE GetObjectItemfromFile
  PRIVATE GetRecordLocations
  PRIVATE TellMeHowManyObjectItemArgs
  PUBLIC  GetNumRangeCheckErrorsFound

  PUBLIC  GetObjectDefMaxArgs
  PRIVATE GetIDFRecordsStats
  PUBLIC  ReportOrphanRecordObjects
  PRIVATE InitSecretObjects
  PRIVATE MakeTransition
  PRIVATE AddRecordFromSection
  PUBLIC  PreProcessorCheck
  PUBLIC  CompactObjectsCheck
  PUBLIC  ParametricObjectsCheck
  PRIVATE DumpCurrentLineBuffer
  PRIVATE ShowAuditErrorMessage
  PUBLIC  PreScanReportingVariables
  PRIVATE IPTrimSigDigits

CONTAINS

  ! MODULE SUBROUTINES:
  !*************************************************************************

  SUBROUTINE ProcessInput

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine processes the input for EnergyPlus.  First, the
    ! input data dictionary is read and interpreted.  Using the structure
    ! from the data dictionary, the actual simulation input file is read.
    ! This file is processed according to the "rules" in the data dictionary
    ! and stored in a local data structure which will be used during the simulation.

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    USE SortAndStringUtilities, ONLY: SetupAndSort
    USE DataOutputs,            ONLY: iNumberOfRecords,iNumberOfDefaultedFields,iTotalFieldsWithDefaults,  &
    iNumberOfAutosizedFields,iTotalAutoSizableFields,iNumberOfAutoCalcedFields,iTotalAutoCalculatableFields

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    ! na

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL FileExists ! Check variable for .idd/.idf files
    LOGICAL :: ErrorsInIDD=.false.   ! to check for any errors flagged during data dictionary processing
    INTEGER :: Loop
    INTEGER :: CountErr
    INTEGER :: Num1
    INTEGER :: Which
    INTEGER :: endcol
    INTEGER :: write_stat
    INTEGER :: read_stat

    CALL InitSecretObjects

    IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10 .OR. DebugFile .EQ. 12) THEN
      WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
    END IF

    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

    EchoInputFile=GetNewUnitNumber()
    OPEN(unit=EchoInputFile,file='eplusout.audit',action='write',iostat=write_stat)
    IF (write_stat /= 0) THEN
      CALL DisplayString('Could not open (write) eplusout.audit.')
      CALL ShowFatalError('ProcessInput: Could not open file "eplusout.audit" for output (write).')
    ENDIF

    INQUIRE(FILE='eplusout.iperr',EXIST=FileExists)
    IF (FileExists) THEN
      CacheIPErrorFile=GetNewUnitNumber()
      open(unit=CacheIPErrorFile,file='eplusout.iperr',action='read', iostat=read_stat)
      IF (read_stat /= 0) THEN
        CALL ShowFatalError('EnergyPlus: Could not open file "eplusout.iperr" for input (read).')
      ENDIF
      close(unit=CacheIPErrorFile,status='delete')
    ENDIF
    CacheIPErrorFile=GetNewUnitNumber()
    OPEN(unit=CacheIPErrorFile,file='eplusout.iperr',action='write',iostat=write_stat)
    IF (write_stat /= 0) THEN
      CALL DisplayString('Could not open (write) eplusout.iperr.')
      CALL ShowFatalError('ProcessInput: Could not open file "eplusout.audit" for output (write).')
    ENDIF

    !               FullName from StringGlobals is used to build file name with Path
    IF (LEN_TRIM(ProgramPath) == 0) THEN     !RS: Line 76244 of the file starts the HPSim part of the IDD
      !FullName='Energy+.idd'
      FullName='Energy+ HPSim.idd'
      !FullName='Energy+ base.idd'    !RS: Reading in the baseline IDD
    ELSE
      !FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+.idd'
      Fullname=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+ HPSim.idd'
      !FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+ base.idd' !RS: Reading in the baseline IDD
    ENDIF
    INQUIRE(file=FullName,EXIST=FileExists)
    IF (.not. FileExists) THEN
      CALL DisplayString('Missing '//TRIM(FullName))
      CALL ShowFatalError('ProcessInput: Energy+.idd missing. Program terminates. Fullname='//TRIM(FullName))
    ENDIF
    IDDFile=GetNewUnitNumber()
    StripCR=.false.
    Open (unit=IDDFile, file=FullName, action='read', iostat=read_stat)
    IF (read_stat /= 0) THEN
      CALL DisplayString('Could not open (read) Energy+.idd.')
      CALL ShowFatalError('ProcessInput: Could not open file "Energy+.idd" for input (read).')
    ENDIF
    IF(IDDFile .EQ. 9 .OR. IDDFile .EQ. 10 .OR. IDDFile .EQ. 12) THEN
      WRITE(*,*) 'Error with IDDFile'    !RS: Debugging: Searching for a mis-set file number
    END IF
    READ(Unit=IDDFile, FMT=fmta) InputLine
    endcol=LEN_TRIM(InputLine)
    IF (endcol > 0) THEN
      IF (ICHAR(InputLine(endcol:endcol)) == iASCII_CR) THEN
        StripCR=.true.
      ENDIF
      IF (ICHAR(InputLine(endcol:endcol)) == iUnicode_end) THEN
        CALL ShowSevereError('ProcessInput: "Energy+.idd" appears to be a Unicode file.')
        CALL ShowContinueError('...This file cannot be read by this program. Please save as PC or Unix file and try again')
        CALL ShowFatalError('Program terminates due to previous condition.')
      ENDIF
    ENDIF
    BACKSPACE(Unit=IDDFile)
    NumLines=0

    DoingInputProcessing=.true.
    IF(EchoInputFile .EQ. 9 .OR. EchoInputFile .EQ. 10 .OR. EchoInputFile .EQ. 12) THEN
      WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
    END IF
    WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Start'
    CALL DisplayString('Processing Data Dictionary')
    ProcessingIDD=.true.

    Call ProcessDataDicFile(ErrorsInIDD)

    ALLOCATE (ListofObjects(NumObjectDefs))
    ListofObjects=ObjectDef(1:NumObjectDefs)%Name
    IF (SortedIDD) THEN
      ALLOCATE (iListofObjects(NumObjectDefs))
      iListOfObjects=0
      CALL SetupAndSort(ListOfObjects,iListOfObjects)
    ENDIF
    ALLOCATE (ObjectStartRecord(NumObjectDefs))
    ObjectStartRecord=0
    ALLOCATE (ObjectGotCount(NumObjectDefs))
    ObjectGotCount=0

    Close (unit=IDDFile)

    IF (NumObjectDefs == 0) THEN
      CALL ShowFatalError('ProcessInput: No objects found in IDD.  Program will terminate.')
      ErrorsInIDD=.true.
    ENDIF
    !  If no fatal to here, rewind EchoInputFile -- only keep processing data...
    IF (.not. ErrorsInIDD) THEN
      REWIND(Unit=EchoInputFile)
    ENDIF

    ProcessingIDD=.false.
    WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Complete'

    WRITE(EchoInputFile,*) ' Maximum number of Alpha Args=',MaxAlphaArgsFound
    WRITE(EchoInputFile,*) ' Maximum number of Numeric Args=',MaxNumericArgsFound
    WRITE(EchoInputFile,*) ' Number of Object Definitions=',NumObjectDefs
    WRITE(EchoInputFile,*) ' Number of Section Definitions=',NumSectionDefs


    WRITE(EchoInputFile,*) ' Processing Input Data File (in.idf) -- Start'

    INQUIRE(file='in.idf',EXIST=FileExists)
    IF (.not. FileExists) THEN
      CALL DisplayString('Missing '//TRIM(CurrentWorkingFolder)//'in.idf')
      CALL ShowFatalError('ProcessInput: in.idf missing. Program terminates.')
    ENDIF

    StripCR=.false.
    IDFFile=GetNewUnitNumber()
    Open (unit=IDFFile, file='in.idf', action='READ', iostat=read_stat)
    IF (read_stat /= 0) THEN
      CALL DisplayString('Could not open (read) in.idf.')
      CALL ShowFatalError('ProcessInput: Could not open file "in.idf" for input (read).')
    ENDIF
    IF(IDFFile .EQ. 9 .OR. IDFFile .EQ. 10 .OR. IDFFile .EQ. 12) THEN
      WRITE(*,*) 'Error with IDFFile'    !RS: Debugging: Searching for a mis-set file number
    END IF
    READ(Unit=IDFFile, FMT=fmta) InputLine
    endcol=LEN_TRIM(InputLine)
    IF (endcol > 0) THEN
      IF (ICHAR(InputLine(endcol:endcol)) == iASCII_CR) THEN
        StripCR=.true.
      ENDIF
      IF (ICHAR(InputLine(endcol:endcol)) == iUnicode_end) THEN
        CALL ShowSevereError('ProcessInput: "in.idf" appears to be a Unicode file.')
        CALL ShowContinueError('...This file cannot be read by this program. Please save as PC or Unix file and try again')
        CALL ShowFatalError('Program terminates due to previous condition.')
      ENDIF
    ENDIF
    BACKSPACE(Unit=IDFFile)
    NumLines=0
    EchoInputLine=.true.
    CALL DisplayString('Processing Input File')

    Call ProcessInputDataFile

    ALLOCATE (ListofSections(NumSectionDefs))
    ListofSections=SectionDef(1:NumSectionDefs)%Name

    Close (unit=IDFFile)


    !!RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)
    !
    !EchoInputFile=GetNewUnitNumber()
    !OPEN(unit=EchoInputFile,file='HPSimVar.audit')
    !!               FullName from StringGlobals is used to build file name with Path
    !!IF (LEN_TRIM(ProgramPath) == 0) THEN
    !!  FullName='HPSim_Variables.idd'
    !!ELSE
    !!  FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'HPSim_Variables.idd'
    !!ENDIF
    !
    !FullName='C:/Users/lab303user/Desktop/GenOpt/HPSim/HPSim_Variables.idd'
    !
    !INQUIRE(file=FullName,EXIST=FileExists)
    !IF (.not. FileExists) THEN
    !  CALL ShowFatalError('Energy+.idd missing. Program terminates. Fullname='//TRIM(FullName))
    !ENDIF
    !IDDFile=GetNewUnitNumber()
    !Open (unit=IDDFile, file=FullName, action='READ')
    !NumLines=0
    !
    !WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Start'
    !
    !Call ProcessDataDicFile2(ErrorsInIDD)
    !
    !ALLOCATE (ListofSections2(NumSectionDefs2), ListofObjects2(NumObjectDefs2))
    !ListofSections2=SectionDef2(1:NumSectionDefs2)%Name
    !ListofObjects2=ObjectDef2(1:NumObjectDefs2)%Name
    !
    !Close (unit=IDDFile)
    !
    !WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Complete'
    !
    !WRITE(EchoInputFile,*) ' Maximum number of Alpha Args=',MaxAlphaArgsFound
    !WRITE(EchoInputFile,*) ' Maximum number of Numeric Args=',MaxNumericArgsFound
    !WRITE(EchoInputFile,*) ' Number of Object Definitions=',NumObjectDefs2
    !WRITE(EchoInputFile,*) ' Number of Section Definitions=',NumSectionDefs2
    !
    !!If no fatal to here, rewind EchoInputFile -- only keep processing data...
    !IF (.not. ErrorsInIDD) THEN
    !  REWIND(Unit=EchoInputFile)
    !ENDIF
    !
    !!IF (LEN_TRIM(ProgramPath) == 0) THEN
    !!  FullName='HPSim_Variables.idf'
    !!ELSE
    !!  FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'HPSim_Variables.idf'
    !!END IF
    !
    !FullName='C:/Users/lab303user/Desktop/GenOpt/HPSim/HPSim_Variables.idf'
    !
    !!FileName = "in.idf"
    !!FileName = "in_longtubes.idf"
    !
    !INQUIRE(file=FullName,EXIST=FileExists)
    !IF (.not. FileExists) THEN
    !   CALL ShowFatalError('Input file missing. Program terminates.')
    !ENDIF
    !
    !IDFFile=GetNewUnitNumber()
    !Open (unit=IDFFile, file = FullName, action='READ')
    !NumLines=0
    !
    !Call ProcessInputDataFile2
    !
    !Close (unit=IDFFile)
    !
    !!RS: Debugging: End of duplication of IDF and IDD-reading code (9/22/14)

    ALLOCATE(cAlphaFieldNames(MaxAlphaIDFDefArgsFound))
    cAlphaFieldNames=Blank
    ALLOCATE(cAlphaArgs(MaxAlphaIDFDefArgsFound))
    cAlphaArgs=Blank
    ALLOCATE(lAlphaFieldBlanks(MaxAlphaIDFDefArgsFound))
    lAlphaFieldBlanks=.false.
    ALLOCATE(cNumericFieldNames(MaxNumericIDFDefArgsFound))
    cNumericFieldNames=Blank
    ALLOCATE(rNumericArgs(MaxNumericIDFDefArgsFound))
    rNumericArgs=0.0d0
    ALLOCATE(lNumericFieldBlanks(MaxNumericIDFDefArgsFound))
    lNumericFieldBlanks=.false.

    ALLOCATE(IDFRecordsGotten(NumIDFRecords))
    IDFRecordsGotten=.false.


    WRITE(EchoInputFile,*) ' Processing Input Data File (in.idf) -- Complete'
    WRITE(EchoInputFile,*) ' Number of IDF "Lines"=',NumIDFRecords
    WRITE(EchoInputFile,*) ' Maximum number of Alpha IDF Args=',MaxAlphaIDFArgsFound
    WRITE(EchoInputFile,*) ' Maximum number of Numeric IDF Args=',MaxNumericIDFArgsFound
    CALL GetIDFRecordsStats(iNumberOfRecords,iNumberOfDefaultedFields,iTotalFieldsWithDefaults,  &
    iNumberOfAutosizedFields,iTotalAutoSizableFields,  &
    iNumberOfAutoCalcedFields,iTotalAutoCalculatableFields)
    WRITE(EchoInputFile,*) ' Number of IDF "Lines"=',iNumberOfRecords
    WRITE(EchoInputFile,*) ' Number of Defaulted Fields=',iNumberOfDefaultedFields
    WRITE(EchoInputFile,*) ' Number of Fields with Defaults=',iTotalFieldsWithDefaults
    WRITE(EchoInputFile,*) ' Number of Autosized Fields=',iNumberOfAutosizedFields
    WRITE(EchoInputFile,*) ' Number of Autosizable Fields =',iTotalAutoSizableFields
    WRITE(EchoInputFile,*) ' Number of Autocalculated Fields=',iNumberOfAutoCalcedFields
    WRITE(EchoInputFile,*) ' Number of Autocalculatable Fields =',iTotalAutoCalculatableFields

    CountErr=0
    DO Loop=1,NumIDFSections
      IF (SectionsonFile(Loop)%LastRecord /= 0) CYCLE
      IF (MakeUPPERCase(SectionsonFile(Loop)%Name) == 'REPORT VARIABLE DICTIONARY') CYCLE
      IF (CountErr == 0) THEN
        !       CALL ShowSevereError('IP: Potential errors in IDF processing -- see .audit file for details.')  !RS: Secret Search String
        WRITE(EchoInputFile,fmta) ' Potential errors in IDF processing:'
        WRITE(DebugFile,*) CountErr  !RS: Debugging
      ENDIF
      CountErr=CountErr+1
      Which=SectionsOnFile(Loop)%FirstRecord
      IF (Which > 0) THEN
        IF (SortedIDD) THEN
          Num1=FindItemInSortedList(IDFRecords(Which)%Name,ListOfObjects,NumObjectDefs)
          IF (Num1 /= 0) Num1=iListOfObjects(Num1)
        ELSE
          Num1=FindItemInList(IDFRecords(Which)%Name,ListOfObjects,NumObjectDefs)
        ENDIF
        IF (ObjectDef(Num1)%NameAlpha1 .and. IDFRecords(Which)%NumAlphas > 0) THEN
          WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
          TRIM(SectionsonFile(Loop)%Name)//  &
          ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
          '], Object Type Preceding='//TRIM(IDFRecords(Which)%Name)//   &
          ', Object Name='//TRIM(IDFRecords(Which)%Alphas(1))
        ELSE
          WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
          TRIM(SectionsonFile(Loop)%Name)//  &
          ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
          '], Object Type Preceding='//TRIM(IDFRecords(Which)%Name)//   &
          ', Name field not recorded for Object.'
        ENDIF
      ELSE
        WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
        TRIM(SectionsonFile(Loop)%Name)//  &
        ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
        '], No prior Objects.'
      ENDIF
    ENDDO

    IF (NumIDFRecords == 0) THEN
      CALL ShowSevereError('IP: The IDF file has no records.')
      NumMiscErrorsFound=NumMiscErrorsFound+1
    ENDIF

    ! Check for required objects
    DO Loop=1,NumObjectDefs
      IF (.not. ObjectDef(Loop)%RequiredObject) CYCLE
      IF (ObjectDef(Loop)%NumFound > 0) CYCLE
      !     CALL ShowSevereError('IP: Required Object="'//trim(ObjectDef(Loop)%Name)//'" not found in IDF.')  !RS: Secret Search String
      IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
        WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
      END IF
      WRITE(DebugFile,*) 'Required Object="'//TRIM(ObjectDef(Loop)%Name)//'" not found in IDF.'
      NumMiscErrorsFound=NumMiscErrorsFound+1
    ENDDO

    IF (TotalAuditErrors > 0) THEN
      !CALL ShowWarningError('IP: Note -- Some missing fields have been filled with defaults.'//  &  !RS: Secret Search String
      !   ' See the audit output file for details.')
      WRITE(DebugFile,*) 'IP: Note -- Some missing fields have been filled with defaults.'//  &
      ' See the audit output file for details.'
    ENDIF

    IF (NumOutOfRangeErrorsFound > 0) THEN
      CALL ShowSevereError('IP: Out of "range" values found in input')
    ENDIF

    IF (NumBlankReqFieldFound > 0) THEN
      CALL ShowSevereError('IP: Blank "required" fields found in input')
    ENDIF

    IF (NumMiscErrorsFound > 0) THEN
      !CALL ShowSevereError('IP: Other miscellaneous errors found in input') !RS: Secret Search String
      WRITE(DebugFile,*) 'Other miscellaneous errors found in input'
    ENDIF

    IF (OverallErrorFlag) THEN
      !CALL ShowSevereError('IP: Possible incorrect IDD File')
      !CALL ShowContinueError('IDD Version:"'//TRIM(IDDVerString)//'"')
      WRITE(DebugFile,*) 'IP: Possible incorrect IDD File'   !RS: Secret Search String
      WRITE(DebugFile,*) 'IDD Version: "'//TRIM(IDDVerString)//'"'
      DO Loop=1,NumIDFRecords
        IF (SameString(IDFRecords(Loop)%Name,'Version')) THEN
          Num1=LEN_TRIM(MatchVersion)
          IF (MatchVersion(Num1:Num1) == '0') THEN
            Which=INDEX(IDFRecords(Loop)%Alphas(1)(1:Num1-2),MatchVersion(1:Num1-2))
          ELSE
            Which=INDEX(IDFRecords(Loop)%Alphas(1),MatchVersion)
          ENDIF
          IF (Which /= 1) THEN
            CALL ShowContinueError('Version in IDF="'//TRIM(IDFRecords(Loop)%Alphas(1))//  &
            '" not the same as expected="'//TRIM(MatchVersion)//'"')
          ENDIF
          EXIT
        ENDIF
      ENDDO
      CALL ShowContinueError('Possible Invalid Numerics or other problems')
      ! Fatal error will now occur during post IP processing check in Simulation manager.
    ENDIF
    RETURN

  END SUBROUTINE ProcessInput

  SUBROUTINE ProcessInput2

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine processes the input for EnergyPlus.  First, the
    ! input data dictionary is read and interpreted.  Using the structure
    ! from the data dictionary, the actual simulation input file is read.
    ! This file is processed according to the "rules" in the data dictionary
    ! and stored in a local data structure which will be used during the simulation.

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    USE SortAndStringUtilities, ONLY: SetupAndSort
    USE DataOutputs,            ONLY: iNumberOfRecords,iNumberOfDefaultedFields,iTotalFieldsWithDefaults,  &
    iNumberOfAutosizedFields,iTotalAutoSizableFields,iNumberOfAutoCalcedFields,iTotalAutoCalculatableFields

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    ! na

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL FileExists ! Check variable for .idd/.idf files
    LOGICAL :: ErrorsInIDD=.false.   ! to check for any errors flagged during data dictionary processing
    INTEGER :: Loop
    INTEGER :: CountErr
    INTEGER :: Num1
    INTEGER :: Which
    INTEGER :: endcol
    INTEGER :: write_stat
    INTEGER :: read_stat
    LOGICAL :: file_exists
    CHARACTER(LEN=500) :: FolderPath
    INTEGER :: ios = 0
    INTEGER :: pos
    !CALL InitSecretObjects

    IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10 .OR. DebugFile .EQ. 12) THEN
      WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
    END IF

    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

    INQUIRE(FILE='HPSim_Vars.audit',EXIST=FileExists)
    IF (FileExists) THEN
      !EchoInputFile=GetNewUnitNumber()
      open(unit=EchoInputFile,file='HPSim_Vars.audit',action='write', iostat=write_stat)
      !IF (read_stat /= 0) THEN
      !  CALL ShowFatalError('EnergyPlus: Could not open file "eplusout.iperr" for input (read).')
      !ENDIF
      !     close(unit=EchoInputFile,status='delete')
      !EchoInputFile=GetNewUnitNumber()
      !OPEN(unit=EchoInputFile,file='HPSim_Vars.audit') !,action='write',iostat=write_stat)
    ELSE
      EchoInputFile=GetNewUnitNumber()
      OPEN(unit=EchoInputFile,file='HPSim_Vars.audit',action='write',iostat=write_stat)
    ENDIF
    !IF (write_stat /= 0) THEN
    !  CALL DisplayString('Could not open (write) eplusout.audit.')
    !  CALL ShowFatalError('ProcessInput: Could not open file "eplusout.audit" for output (write).')
    !ENDIF

    INQUIRE(FILE='HPSim_Vars.iperr',EXIST=FileExists)
    IF (FileExists) THEN
      CacheIPErrorFile=GetNewUnitNumber()
      open(unit=CacheIPErrorFile,file='HPSim_Vars.iperr',action='write',iostat=write_stat)
      !IF (read_stat /= 0) THEN
      !  CALL ShowFatalError('EnergyPlus: Could not open file "eplusout.iperr" for input (read).')
      !ENDIF
      !close(unit=CacheIPErrorFile,status='delete')
      !CacheIPErrorFile=GetNewUnitNumber()
      !OPEN(unit=CacheIPErrorFile,file='HPSim_Vars.iperr',action='write',iostat=write_stat)
    ELSE
      CacheIPErrorFile=GetNewUnitNumber()
      OPEN(unit=CacheIPErrorFile,file='HPSim_Vars.iperr',action='write',iostat=write_stat)
    ENDIF

    !IF (write_stat /= 0) THEN
    !  CALL DisplayString('Could not open (write) eplusout.iperr.')
    !  CALL ShowFatalError('ProcessInput: Could not open file "eplusout.audit" for output (write).')
    !ENDIF

    !!               FullName from StringGlobals is used to build file name with Path
    !IF (LEN_TRIM(ProgramPath) == 0) THEN     !RS: Line 76244 of the file starts the HPSim part of the IDD
    !  !FullName='Energy+.idd'
    !  FullName='Energy+ HPSim.idd'
    !  !FullName='Energy+ base.idd'    !RS: Reading in the baseline IDD
    !ELSE
    !  !FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+.idd'
    !  Fullname=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+ HPSim.idd'
    !  !FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'Energy+ base.idd' !RS: Reading in the baseline IDD
    !ENDIF
    !INQUIRE(file=FullName,EXIST=FileExists)
    !IF (.not. FileExists) THEN
    !  CALL DisplayString('Missing '//TRIM(FullName))
    !  CALL ShowFatalError('ProcessInput: Energy+.idd missing. Program terminates. Fullname='//TRIM(FullName))
    !ENDIF
    !IDDFile=GetNewUnitNumber()
    !StripCR=.false.
    !Open (unit=IDDFile, file=FullName, action='read', iostat=read_stat)
    !IF (read_stat /= 0) THEN
    !  CALL DisplayString('Could not open (read) Energy+.idd.')
    !  CALL ShowFatalError('ProcessInput: Could not open file "Energy+.idd" for input (read).')
    !ENDIF
    !IF(IDDFile .EQ. 9 .OR. IDDFile .EQ. 10 .OR. IDDFile .EQ. 12) THEN
    ! WRITE(*,*) 'Error with IDDFile'    !RS: Debugging: Searching for a mis-set file number
    !END IF
    !READ(Unit=IDDFile, FMT=fmta) InputLine
    !endcol=LEN_TRIM(InputLine)
    !IF (endcol > 0) THEN
    !  IF (ICHAR(InputLine(endcol:endcol)) == iASCII_CR) THEN
    !    StripCR=.true.
    !  ENDIF
    !  IF (ICHAR(InputLine(endcol:endcol)) == iUnicode_end) THEN
    !    CALL ShowSevereError('ProcessInput: "Energy+.idd" appears to be a Unicode file.')
    !    CALL ShowContinueError('...This file cannot be read by this program. Please save as PC or Unix file and try again')
    !    CALL ShowFatalError('Program terminates due to previous condition.')
    !  ENDIF
    !ENDIF
    !BACKSPACE(Unit=IDDFile)
    !NumLines=0
    !
    !DoingInputProcessing=.true.
    !IF(EchoInputFile .EQ. 9 .OR. EchoInputFile .EQ. 10 .OR. EchoInputFile .EQ. 12) THEN
    ! WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
    !END IF
    !WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Start'
    !CALL DisplayString('Processing Data Dictionary')
    !ProcessingIDD=.true.
    !
    !Call ProcessDataDicFile(ErrorsInIDD)
    !
    !ALLOCATE (ListofObjects(NumObjectDefs))
    !ListofObjects=ObjectDef(1:NumObjectDefs)%Name
    !IF (SortedIDD) THEN
    !  ALLOCATE (iListofObjects(NumObjectDefs))
    !  iListOfObjects=0
    !  CALL SetupAndSort(ListOfObjects,iListOfObjects)
    !ENDIF
    !ALLOCATE (ObjectStartRecord(NumObjectDefs))
    !ObjectStartRecord=0
    !ALLOCATE (ObjectGotCount(NumObjectDefs))
    !ObjectGotCount=0
    !
    !Close (unit=IDDFile)
    !
    !IF (NumObjectDefs == 0) THEN
    !  CALL ShowFatalError('ProcessInput: No objects found in IDD.  Program will terminate.')
    !  ErrorsInIDD=.true.
    !ENDIF
    !!  If no fatal to here, rewind EchoInputFile -- only keep processing data...
    !IF (.not. ErrorsInIDD) THEN
    !  REWIND(Unit=EchoInputFile)
    !ENDIF
    !
    !ProcessingIDD=.false.
    !WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Complete'
    !
    !WRITE(EchoInputFile,*) ' Maximum number of Alpha Args=',MaxAlphaArgsFound
    !WRITE(EchoInputFile,*) ' Maximum number of Numeric Args=',MaxNumericArgsFound
    !WRITE(EchoInputFile,*) ' Number of Object Definitions=',NumObjectDefs
    !WRITE(EchoInputFile,*) ' Number of Section Definitions=',NumSectionDefs
    !
    !
    !WRITE(EchoInputFile,*) ' Processing Input Data File (in.idf) -- Start'
    !
    !INQUIRE(file='in.idf',EXIST=FileExists)
    !IF (.not. FileExists) THEN
    !  CALL DisplayString('Missing '//TRIM(CurrentWorkingFolder)//'in.idf')
    !  CALL ShowFatalError('ProcessInput: in.idf missing. Program terminates.')
    !ENDIF
    !
    !StripCR=.false.
    !IDFFile=GetNewUnitNumber()
    !Open (unit=IDFFile, file='in.idf', action='READ', iostat=read_stat)
    !IF (read_stat /= 0) THEN
    !  CALL DisplayString('Could not open (read) in.idf.')
    !  CALL ShowFatalError('ProcessInput: Could not open file "in.idf" for input (read).')
    !ENDIF
    !IF(IDFFile .EQ. 9 .OR. IDFFile .EQ. 10 .OR. IDFFile .EQ. 12) THEN
    ! WRITE(*,*) 'Error with IDFFile'    !RS: Debugging: Searching for a mis-set file number
    !END IF
    !READ(Unit=IDFFile, FMT=fmta) InputLine
    !endcol=LEN_TRIM(InputLine)
    !IF (endcol > 0) THEN
    !  IF (ICHAR(InputLine(endcol:endcol)) == iASCII_CR) THEN
    !    StripCR=.true.
    !  ENDIF
    !  IF (ICHAR(InputLine(endcol:endcol)) == iUnicode_end) THEN
    !    CALL ShowSevereError('ProcessInput: "in.idf" appears to be a Unicode file.')
    !    CALL ShowContinueError('...This file cannot be read by this program. Please save as PC or Unix file and try again')
    !    CALL ShowFatalError('Program terminates due to previous condition.')
    !  ENDIF
    !ENDIF
    !BACKSPACE(Unit=IDFFile)
    !NumLines=0
    !EchoInputLine=.true.
    !CALL DisplayString('Processing Input File')
    !
    !Call ProcessInputDataFile
    !
    !ALLOCATE (ListofSections(NumSectionDefs))
    !ListofSections=SectionDef(1:NumSectionDefs)%Name
    !
    !Close (unit=IDFFile)


    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

    !INQUIRE(FILE='HPSimVar.audit',EXIST=FileExists)
    !IF (FileExists) THEN
    !  EchoInputFile=GetNewUnitNumber()
    !  open(unit=EchoInputFile,file='HPSimVar.audit',action='read', iostat=read_stat)
    !  IF (read_stat /= 0) THEN
    !    CALL ShowFatalError('EnergyPlus: Could not open file "HPSimVar.audit" for input (read).')
    !  ENDIF
    !  close(unit=CacheIPErrorFile,status='delete')
    !ENDIF
    !EchoInputFile=GetNewUnitNumber()
    !OPEN(unit=EchoInputFile,file='HPSimVar.audit',action='write',iostat=write_stat)
    !IF (write_stat /= 0) THEN
    !  CALL DisplayString('Could not open (write) HPSimVar.audit.')
    !  CALL ShowFatalError('ProcessInput: Could not open file "eplusout.audit" for output (write).')
    !ENDIF

    !EchoInputFile=GetNewUnitNumber()
    !OPEN(unit=EchoInputFile,file='eplusout.audit')

    !EchoInputFile=GetNewUnitNumber()
    !OPEN(unit=EchoInputFile,file='HPSimVar.audit')
    !               FullName from StringGlobals is used to build file name with Path
    IF (LEN_TRIM(ProgramPath) == 0) THEN
      FullName='HPSim_Variables.idd'
    ELSE
      FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'HPSim_Variables.idd'
    ENDIF

    !FullName='C:/Users/lab303user/Desktop/GenOpt/HPSim/HPSim_Variables.idd'

    INQUIRE(file=FullName,EXIST=FileExists)
    IF (.not. FileExists) THEN
      CALL ShowFatalError('Energy+.idd missing. Program terminates. Fullname='//TRIM(FullName))
    ENDIF
    IDDFile=GetNewUnitNumber()
    Open (unit=IDDFile, file=FullName, action='READ')
    NumLines=0

    WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Start'

    Call ProcessDataDicFile2(ErrorsInIDD)

    IF ( .not. ALLOCATED(ListofSections2) .AND. .not. ALLOCATED(ListofObjects2)) THEN
      ALLOCATE (ListofSections2(NumSectionDefs2), ListofObjects2(NumObjectDefs2))
    END IF
    ListofSections2=SectionDef2(1:NumSectionDefs2)%Name
    ListofObjects2=ObjectDef2(1:NumObjectDefs2)%Name

    Close (unit=IDDFile)

    WRITE(EchoInputFile,*) ' Processing Data Dictionary (Energy+.idd) File -- Complete'

    WRITE(EchoInputFile,*) ' Maximum number of Alpha Args=',MaxAlphaArgsFound
    WRITE(EchoInputFile,*) ' Maximum number of Numeric Args=',MaxNumericArgsFound
    WRITE(EchoInputFile,*) ' Number of Object Definitions=',NumObjectDefs2
    WRITE(EchoInputFile,*) ' Number of Section Definitions=',NumSectionDefs2

    !If no fatal to here, rewind EchoInputFile -- only keep processing data...
    IF (.not. ErrorsInIDD) THEN
      REWIND(Unit=EchoInputFile)
    ENDIF

    !IF (LEN_TRIM(ProgramPath) == 0) THEN
    !  FullName='HPSim_Variables.idf'
    !ELSE
    !  FullName=ProgramPath(1:LEN_TRIM(ProgramPath))//'HPSim_Variables.idf'
    !END IF

    INQUIRE(FILE="FilePathBufferProgram.txt", EXIST=file_exists)
    if (file_exists) THEN
      OPEN (UNIT=580, FILE="FilePathBufferProgram.txt", STATUS="OLD")   ! Current directory
      read(580,'(A)', iostat=ios) FolderPath
      read(580,'(A)', iostat=ios) FolderPath
      CLOSE (UNIT=580)
      FolderPath=TRIM(FolderPath)
      pos = scan(FolderPath, '\')
      FolderPath = TRIM(ADJUSTL(FolderPath(1:pos-1)))
    end if

    FullName=TRIM(ADJUSTL(FolderPath(1:pos-1))) // '\\HPSim_Variables.idf'
    !FullName='C:/GenOptFiles/HPSim_Variables.idf'

    !FileName = "in.idf"
    !FileName = "in_longtubes.idf"

    INQUIRE(file=FullName,EXIST=FileExists)
    IF (.not. FileExists) THEN
      CALL ShowFatalError('Input file missing. Program terminates.')
    ENDIF

    IDFFile=GetNewUnitNumber()
    Open (unit=IDFFile, file = FullName, action='READ')
    NumLines=0

    Call ProcessInputDataFile2

    Close (unit=IDFFile)

    !RS: Debugging: End of duplication of IDF and IDD-reading code (9/22/14)

    !   ALLOCATE(cAlphaFieldNames(MaxAlphaIDFDefArgsFound))
    !   cAlphaFieldNames=Blank
    !   ALLOCATE(cAlphaArgs(MaxAlphaIDFDefArgsFound))
    !   cAlphaArgs=Blank
    !   ALLOCATE(lAlphaFieldBlanks(MaxAlphaIDFDefArgsFound))
    !   lAlphaFieldBlanks=.false.
    !   ALLOCATE(cNumericFieldNames(MaxNumericIDFDefArgsFound))
    !   cNumericFieldNames=Blank
    !   ALLOCATE(rNumericArgs(MaxNumericIDFDefArgsFound))
    !   rNumericArgs=0.0d0
    !   ALLOCATE(lNumericFieldBlanks(MaxNumericIDFDefArgsFound))
    !   lNumericFieldBlanks=.false.
    !
    !   ALLOCATE(IDFRecordsGotten(NumIDFRecords))
    !   IDFRecordsGotten=.false.
    !
    !
    !   WRITE(EchoInputFile,*) ' Processing Input Data File (in.idf) -- Complete'
    !   WRITE(EchoInputFile,*) ' Number of IDF "Lines"=',NumIDFRecords
    !   WRITE(EchoInputFile,*) ' Maximum number of Alpha IDF Args=',MaxAlphaIDFArgsFound
    !   WRITE(EchoInputFile,*) ' Maximum number of Numeric IDF Args=',MaxNumericIDFArgsFound
    !   CALL GetIDFRecordsStats(iNumberOfRecords,iNumberOfDefaultedFields,iTotalFieldsWithDefaults,  &
    !                              iNumberOfAutosizedFields,iTotalAutoSizableFields,  &
    !                              iNumberOfAutoCalcedFields,iTotalAutoCalculatableFields)
    !   WRITE(EchoInputFile,*) ' Number of IDF "Lines"=',iNumberOfRecords
    !   WRITE(EchoInputFile,*) ' Number of Defaulted Fields=',iNumberOfDefaultedFields
    !   WRITE(EchoInputFile,*) ' Number of Fields with Defaults=',iTotalFieldsWithDefaults
    !   WRITE(EchoInputFile,*) ' Number of Autosized Fields=',iNumberOfAutosizedFields
    !   WRITE(EchoInputFile,*) ' Number of Autosizable Fields =',iTotalAutoSizableFields
    !   WRITE(EchoInputFile,*) ' Number of Autocalculated Fields=',iNumberOfAutoCalcedFields
    !   WRITE(EchoInputFile,*) ' Number of Autocalculatable Fields =',iTotalAutoCalculatableFields
    !
    !   CountErr=0
    !   DO Loop=1,NumIDFSections
    !     IF (SectionsonFile(Loop)%LastRecord /= 0) CYCLE
    !     IF (MakeUPPERCase(SectionsonFile(Loop)%Name) == 'REPORT VARIABLE DICTIONARY') CYCLE
    !     IF (CountErr == 0) THEN
    !!       CALL ShowSevereError('IP: Potential errors in IDF processing -- see .audit file for details.')  !RS: Secret Search String
    !       WRITE(EchoInputFile,fmta) ' Potential errors in IDF processing:'
    !       WRITE(DebugFile,*) CountErr  !RS: Debugging
    !     ENDIF
    !     CountErr=CountErr+1
    !     Which=SectionsOnFile(Loop)%FirstRecord
    !     IF (Which > 0) THEN
    !       IF (SortedIDD) THEN
    !         Num1=FindItemInSortedList(IDFRecords(Which)%Name,ListOfObjects,NumObjectDefs)
    !         IF (Num1 /= 0) Num1=iListOfObjects(Num1)
    !       ELSE
    !         Num1=FindItemInList(IDFRecords(Which)%Name,ListOfObjects,NumObjectDefs)
    !       ENDIF
    !       IF (ObjectDef(Num1)%NameAlpha1 .and. IDFRecords(Which)%NumAlphas > 0) THEN
    !         WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
    !               TRIM(SectionsonFile(Loop)%Name)//  &
    !               ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
    !               '], Object Type Preceding='//TRIM(IDFRecords(Which)%Name)//   &
    !               ', Object Name='//TRIM(IDFRecords(Which)%Alphas(1))
    !       ELSE
    !         WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
    !               TRIM(SectionsonFile(Loop)%Name)//  &
    !               ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
    !               '], Object Type Preceding='//TRIM(IDFRecords(Which)%Name)//   &
    !               ', Name field not recorded for Object.'
    !       ENDIF
    !     ELSE
    !       WRITE(EchoInputFile,fmta) ' Potential "semi-colon" misplacement='//  &
    !             TRIM(SectionsonFile(Loop)%Name)//  &
    !             ', at about line number=['//TRIM(IPTrimSigDigits(SectionsonFile(Loop)%FirstLineNo))//  &
    !             '], No prior Objects.'
    !     ENDIF
    !   ENDDO

    IF (NumIDFRecords == 0) THEN
      CALL ShowSevereError('IP: The IDF file has no records.')
      NumMiscErrorsFound=NumMiscErrorsFound+1
    ENDIF

    ! Check for required objects
    DO Loop=1,NumObjectDefs
      IF (.not. ObjectDef(Loop)%RequiredObject) CYCLE
      IF (ObjectDef(Loop)%NumFound > 0) CYCLE
      !     CALL ShowSevereError('IP: Required Object="'//trim(ObjectDef(Loop)%Name)//'" not found in IDF.')  !RS: Secret Search String
      IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
        WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
      END IF
      WRITE(DebugFile,*) 'Required Object="'//TRIM(ObjectDef(Loop)%Name)//'" not found in IDF.'
      NumMiscErrorsFound=NumMiscErrorsFound+1
    ENDDO

    IF (TotalAuditErrors > 0) THEN
      !CALL ShowWarningError('IP: Note -- Some missing fields have been filled with defaults.'//  &  !RS: Secret Search String
      !   ' See the audit output file for details.')
      WRITE(DebugFile,*) 'IP: Note -- Some missing fields have been filled with defaults.'//  &
      ' See the audit output file for details.'
    ENDIF

    IF (NumOutOfRangeErrorsFound > 0) THEN
      CALL ShowSevereError('IP: Out of "range" values found in input')
    ENDIF

    IF (NumBlankReqFieldFound > 0) THEN
      CALL ShowSevereError('IP: Blank "required" fields found in input')
    ENDIF

    IF (NumMiscErrorsFound > 0) THEN
      !CALL ShowSevereError('IP: Other miscellaneous errors found in input') !RS: Secret Search String
      WRITE(DebugFile,*) 'Other miscellaneous errors found in input'
    ENDIF

    IF (OverallErrorFlag) THEN
      !CALL ShowSevereError('IP: Possible incorrect IDD File')
      !CALL ShowContinueError('IDD Version:"'//TRIM(IDDVerString)//'"')
      WRITE(DebugFile,*) 'IP: Possible incorrect IDD File'   !RS: Secret Search String
      WRITE(DebugFile,*) 'IDD Version: "'//TRIM(IDDVerString)//'"'
      DO Loop=1,NumIDFRecords
        IF (SameString(IDFRecords(Loop)%Name,'Version')) THEN
          Num1=LEN_TRIM(MatchVersion)
          IF (MatchVersion(Num1:Num1) == '0') THEN
            Which=INDEX(IDFRecords(Loop)%Alphas(1)(1:Num1-2),MatchVersion(1:Num1-2))
          ELSE
            Which=INDEX(IDFRecords(Loop)%Alphas(1),MatchVersion)
          ENDIF
          IF (Which /= 1) THEN
            CALL ShowContinueError('Version in IDF="'//TRIM(IDFRecords(Loop)%Alphas(1))//  &
            '" not the same as expected="'//TRIM(MatchVersion)//'"')
          ENDIF
          EXIT
        ENDIF
      ENDDO
      CALL ShowContinueError('Possible Invalid Numerics or other problems')
      ! Fatal error will now occur during post IP processing check in Simulation manager.
    ENDIF
    RETURN

  END SUBROUTINE ProcessInput2

  SUBROUTINE ProcessDataDicFile(ErrorsFound)

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine processes data dictionary file for EnergyPlus.
    ! The structure of the sections and objects are stored in derived
    ! types (SectionDefs and ObjectDefs)

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if any errors flagged during IDD processing

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL  :: EndofFile = .false.        ! True when End of File has been reached (IDD or IDF)
    INTEGER Pos                            ! Test of scanning position on the current input line
    TYPE (SectionsDefinition), ALLOCATABLE :: TempSectionDef(:)  ! Like SectionDef, used during Re-allocation
    TYPE (ObjectsDefinition), ALLOCATABLE :: TempObjectDef(:)    ! Like ObjectDef, used during Re-allocation
    LOGICAL BlankLine


    MaxSectionDefs=SectionDefAllocInc
    MaxObjectDefs=ObjectDefAllocInc

    ALLOCATE (SectionDef(MaxSectionDefs))

    ALLOCATE(ObjectDef(MaxObjectDefs))

    NumObjectDefs=0
    NumSectionDefs=0
    EndofFile=.false.

    ! Read/process first line (Version info)
    DO WHILE (.not. EndofFile)
      CALL ReadInputLine(IDDFile,Pos,BlankLine,InputLineLength,EndofFile)
      IF (EndofFile) CYCLE
      Pos=INDEX(InputLine,'!IDD_Version')
      IF (Pos /= 0) THEN
        IDDVerString=InputLine(2:LEN_TRIM(InputLine))
      ENDIF
      EXIT
    ENDDO

    DO WHILE (.not. EndofFile)
      CALL ReadInputLine(IDDFile,Pos,BlankLine,InputLineLength,EndofFile)
      IF (BlankLine .or. EndofFile) CYCLE
      Pos=SCAN(InputLine(1:InputLineLength),',;')
      If (Pos /= 0) then

        If (InputLine(Pos:Pos) == ';') then
          CALL AddSectionDef(InputLine(1:Pos-1),ErrorsFound)
          IF (NumSectionDefs == MaxSectionDefs) THEN
            ALLOCATE (TempSectionDef(MaxSectionDefs+SectionDefAllocInc))
            TempSectionDef(1:MaxSectionDefs)=SectionDef
            DEALLOCATE (SectionDef)
            ALLOCATE (SectionDef(MaxSectionDefs+SectionDefAllocInc))
            SectionDef=TempSectionDef
            DEALLOCATE (TempSectionDef)
            MaxSectionDefs=MaxSectionDefs+SectionDefAllocInc
          ENDIF
        else
          CALL AddObjectDefandParse(InputLine(1:Pos-1),Pos,EndofFile,ErrorsFound)
          IF (NumObjectDefs == MaxObjectDefs) THEN
            ALLOCATE (TempObjectDef(MaxObjectDefs+ObjectDefAllocInc))
            TempObjectDef(1:MaxObjectDefs)=ObjectDef
            DEALLOCATE (ObjectDef)
            ALLOCATE (ObjectDef(MaxObjectDefs+ObjectDefAllocInc))
            ObjectDef=TempObjectDef
            DEALLOCATE (TempObjectDef)
            MaxObjectDefs=MaxObjectDefs+ObjectDefAllocInc
          ENDIF
        endif

      else
        CALL ShowSevereError('IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' , or ; expected on this line',EchoInputFile)
        ErrorsFound=.true.
      endif

    END DO

    RETURN

  END SUBROUTINE ProcessDataDicFile

  SUBROUTINE AddSectionDef(ProposedSection,ErrorsFound)

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine adds a new section to SectionDefs.

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    CHARACTER(len=*), INTENT(IN) :: ProposedSection  ! Proposed Section to be added
    LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if errors found here

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    CHARACTER(len=MaxSectionNameLength) SqueezedSection  ! Input Argument, Left-Justified and Uppercase
    LOGICAL ErrFlag  ! Local error flag.  When True, Proposed Section is not added to global list

    SqueezedSection=MakeUPPERCase(ADJUSTL(ProposedSection))
    IF (LEN_TRIM(ADJUSTL(ProposedSection)) > MaxSectionNameLength) THEN
      CALL ShowWarningError('IP: Section length exceeds maximum, will be truncated='//TRIM(ProposedSection),EchoInputFile)
      CALL ShowContinueError('Will be processed as Section='//TRIM(SqueezedSection),EchoInputFile)
      ErrorsFound=.true.
    ENDIF
    ErrFlag=.false.

    IF (SqueezedSection /= Blank) THEN
      IF (FindItemInList(SqueezedSection,SectionDef%Name,NumSectionDefs) > 0) THEN
        CALL ShowSevereError('IP: Already a Section called '//TRIM(SqueezedSection)//'. This definition ignored.',EchoInputFile)
        ! Error Condition
        ErrFlag=.true.
        ErrorsFound=.true.
      ENDIF
    ELSE
      CALL ShowSevereError('IP: Blank Sections not allowed.  Review eplusout.audit file.',EchoInputFile)
      ErrFlag=.true.
      ErrorsFound=.true.
    ENDIF

    IF (.not. ErrFlag) THEN
      NumSectionDefs=NumSectionDefs+1
      SectionDef(NumSectionDefs)%Name=SqueezedSection
      SectionDef(NumSectionDefs)%NumFound=0
    ENDIF

    RETURN

  END SUBROUTINE AddSectionDef

  SUBROUTINE AddObjectDefandParse(ProposedObject,CurPos,EndofFile,ErrorsFound)

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine processes data dictionary file for EnergyPlus.
    ! The structure of the sections and objects are stored in derived
    ! types (SectionDefs and ObjectDefs)

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS
    CHARACTER(len=*), INTENT(IN) :: ProposedObject  ! Proposed Object to Add
    INTEGER, INTENT(INOUT) :: CurPos ! Current position (initially at first ',') of InputLine
    LOGICAL, INTENT(INOUT) :: EndofFile ! End of File marker
    LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if errors found here

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    CHARACTER(len=MaxObjectNameLength) SqueezedObject  ! Input Object, Left Justified, UpperCase
    INTEGER Count  ! Count on arguments, loop
    INTEGER Pos    ! Position scanning variable
    LOGICAL EndofObjectDef   ! Set to true when ; has been found
    LOGICAL ErrFlag   ! Local Error condition flag, when true, object not added to Global list
    CHARACTER(len=1) TargetChar   ! Single character scanned to test for current field type (A or N)
    LOGICAL BlankLine ! True when this line is "blank" (may have comment characters as first character on line)
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: AlphaorNumeric    ! Array of argument designations, True is Alpha,
    ! False is numeric, saved in ObjectDef when done
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempAN            ! Array (ref: AlphaOrNumeric) for re-allocation procedure
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: RequiredFields    ! Array of argument required fields
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempRqF           ! Array (ref: RequiredFields) for re-allocation procedure
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: AlphRetainCase    ! Array of argument for retain case
    LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempRtC           ! Array (ref: AlphRetainCase) for re-allocation procedure
    CHARACTER(len=MaxFieldNameLength),   &
    ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldChecks   ! Array with alpha field names
    CHARACTER(len=MaxFieldNameLength),   &
    ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFC           ! Array (ref: AlphFieldChecks) for re-allocation procedure
    CHARACTER(len=MaxObjectNameLength),   &
    ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldDefaults ! Array with alpha field defaults
    CHARACTER(len=MaxObjectNameLength),   &
    ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFD           ! Array (ref: AlphFieldDefaults) for re-allocation procedure
    TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: NumRangeChecks  ! Structure for Range Check, Defaults of numeric fields
    TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: TempChecks ! Structure (ref: NumRangeChecks) for re-allocation procedure
    LOGICAL MinMax   ! Set to true when MinMax field has been found by ReadInputLine
    LOGICAL Default  ! Set to true when Default field has been found by ReadInputLine
    LOGICAL AutoSize ! Set to true when Autosizable field has been found by ReadInputLine
    LOGICAL AutoCalculate ! Set to true when Autocalculatable field has been found by ReadInputLine
    CHARACTER(len=32) MinMaxString ! Set from ReadInputLine
    CHARACTER(len=MaxObjectNameLength) AlphDefaultString
    INTEGER WhichMinMax   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
    REAL(r64) Value  ! Value returned by ReadInputLine (either min, max, default, autosize or autocalculate)
    LOGICAL MinMaxError  ! Used to see if min, max, defaults have been set appropriately (True if error)
    INTEGER,SAVE   :: MaxANArgs=7700  ! Current count of Max args to object  (9/2010)
    LOGICAL ErrorsFoundFlag
    INTEGER,SAVE :: PrevSizeNumNumeric = -1
    INTEGER,SAVE :: PrevCount  = -1
    INTEGER,SAVE :: PrevSizeNumAlpha = -1
    INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

    IF (.not. ALLOCATED(AlphaorNumeric)) THEN
      ALLOCATE (AlphaorNumeric(0:MaxANArgs))
      ALLOCATE (RequiredFields(0:MaxANArgs))
      ALLOCATE (AlphRetainCase(0:MaxANArgs))
      ALLOCATE (NumRangeChecks(MaxANArgs))
      ALLOCATE (AlphFieldChecks(MaxANArgs))
      ALLOCATE (AlphFieldDefaults(MaxANArgs))
      ALLOCATE (ObsoleteObjectsRepNames(0))
    ENDIF

    SqueezedObject=MakeUPPERCase(ADJUSTL(ProposedObject))
    IF (LEN_TRIM(ADJUSTL(ProposedObject)) > MaxObjectNameLength) THEN
      CALL ShowWarningError('IP: Object length exceeds maximum, will be truncated='//TRIM(ProposedObject),EchoInputFile)
      CALL ShowContinueError('Will be processed as Object='//TRIM(SqueezedObject),EchoInputFile)
      ErrorsFound=.true.
    ENDIF

    ! Start of Object parse, set object level items
    ErrFlag=.false.
    ErrorsFoundFlag=.false.
    MinimumNumberOfFields=0
    ObsoleteObject=.false.
    UniqueObject=.false.
    RequiredObject=.false.
    ExtensibleObject=.false.
    ExtensibleNumFields=0
    MinMax=.false.
    Default=.false.
    AutoSize=.false.
    AutoCalculate=.false.
    WhichMinMax=0


    IF (SqueezedObject /= Blank) THEN
      IF (FindItemInList(SqueezedObject,ObjectDef%Name,NumObjectDefs) > 0) THEN
        CALL ShowSevereError('IP: Already an Object called '//TRIM(SqueezedObject)//'. This definition ignored.',EchoInputFile)
        ! Error Condition
        ErrFlag=.true.
        ! Rest of Object has to be processed. Error condition will be caught
        ! at end
        ErrorsFound=.true.
      ENDIF
    ELSE
      ErrFlag=.true.
      ErrorsFound=.true.
    ENDIF

    NumObjectDefs=NumObjectDefs+1
    ObjectDef(NumObjectDefs)%Name=SqueezedObject
    ObjectDef(NumObjectDefs)%NumParams=0
    ObjectDef(NumObjectDefs)%NumAlpha=0
    ObjectDef(NumObjectDefs)%NumNumeric=0
    ObjectDef(NumObjectDefs)%NumFound=0
    ObjectDef(NumObjectDefs)%MinNumFields=0
    ObjectDef(NumObjectDefs)%NameAlpha1=.false.
    ObjectDef(NumObjectDefs)%ObsPtr=0
    ObjectDef(NumObjectDefs)%UniqueObject=.false.
    ObjectDef(NumObjectDefs)%RequiredObject=.false.
    ObjectDef(NumObjectDefs)%ExtensibleObject=.false.
    ObjectDef(NumObjectDefs)%ExtensibleNum=0

    IF (PrevCount .EQ. -1) THEN
      PrevCount = MaxANArgs
    END IF

    AlphaorNumeric(1:PrevCount)=.true.
    RequiredFields(1:PrevCount)=.false.
    AlphRetainCase(1:PrevCount)=.false.

    IF (PrevSizeNumAlpha .EQ. -1) THEN
      PrevSizeNumAlpha = MaxANArgs
    END IF

    AlphFieldChecks(1:PrevSizeNumAlpha)=Blank
    AlphFieldDefaults(1:PrevSizeNumAlpha)=Blank

    IF (PrevSizeNumNumeric .EQ. -1) THEN
      PrevSizeNumNumeric = MaxANArgs
    END IF

    !clear only portion of NumRangeChecks array used in the previous
    !call to reduce computation time to clear this large array.
    NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxChk=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%WhichMinMax(1)=0
    NumRangeChecks(1:PrevSizeNumNumeric)%WhichMinMax(2)=0
    NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxString(1)=Blank
    NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxString(2)=Blank
    NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxValue(1)=0.0
    NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxValue(2)=0.0
    NumRangeChecks(1:PrevSizeNumNumeric)%Default=0.0
    NumRangeChecks(1:PrevSizeNumNumeric)%DefaultChk=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%DefAutoSize=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%DefAutoCalculate=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%FieldNumber=0
    NumRangeChecks(1:PrevSizeNumNumeric)%FieldName=Blank
    NumRangeChecks(1:PrevSizeNumNumeric)%AutoSizable=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%AutoSizeValue=DefAutoSizeValue
    NumRangeChecks(1:PrevSizeNumNumeric)%AutoCalculatable=.false.
    NumRangeChecks(1:PrevSizeNumNumeric)%AutoCalculateValue=DefAutoCalculateValue

    Count=0
    EndofObjectDef=.false.
    ! Parse rest of Object Definition

    DO WHILE (.not. EndofFile .and. .not. EndofObjectDef)

      IF (CurPos <= InputLineLength) THEN
        Pos=SCAN(InputLine(CurPos:InputLineLength),AlphaNum)
        IF (Pos > 0) then

          Count=Count+1
          RequiredField=.false.
          RetainCaseFlag=.false.

          IF (Count > MaxANArgs) THEN   ! Reallocation
            ALLOCATE(TempAN(0:MaxANArgs+ANArgsDefAllocInc))
            TempAN=.false.
            TempAN(0:MaxANArgs)=AlphaorNumeric
            DEALLOCATE(AlphaorNumeric)
            ALLOCATE(TempRqF(0:MaxANArgs+ANArgsDefAllocInc))
            TempRqF=.false.
            TempRqF(0:MaxANArgs)=RequiredFields
            DEALLOCATE(RequiredFields)
            ALLOCATE(TempRtC(0:MaxANArgs+ANArgsDefAllocInc))
            TempRtC=.false.
            TempRtC(0:MaxANArgs)=AlphRetainCase
            DEALLOCATE(AlphRetainCase)
            ALLOCATE(TempChecks(MaxANArgs+ANArgsDefAllocInc))
            TempChecks(1:MaxANArgs)=NumRangeChecks(1:MaxANArgs)
            DEALLOCATE(NumRangeChecks)
            ALLOCATE(TempAFC(MaxANArgs+ANArgsDefAllocInc))
            TempAFC=Blank
            TempAFC(1:MaxANArgs)=AlphFieldChecks
            DEALLOCATE(AlphFieldChecks)
            ALLOCATE(TempAFD(MaxANArgs+ANArgsDefAllocInc))
            TempAFD=Blank
            TempAFD(1:MaxANArgs)=AlphFieldDefaults
            DEALLOCATE(AlphFieldDefaults)
            ALLOCATE(AlphaorNumeric(0:MaxANArgs+ANArgsDefAllocInc))
            AlphaorNumeric=TempAN
            DEALLOCATE(TempAN)
            ALLOCATE(RequiredFields(0:MaxANArgs+ANArgsDefAllocInc))
            RequiredFields=TempRqF
            DEALLOCATE(TempRqF)
            ALLOCATE(AlphRetainCase(0:MaxANArgs+ANArgsDefAllocInc))
            AlphRetainCase=TempRtC
            DEALLOCATE(TempRtC)
            ALLOCATE(NumRangeChecks(MaxANArgs+ANArgsDefAllocInc))
            NumRangeChecks=TempChecks
            DEALLOCATE(TempChecks)
            ALLOCATE(AlphFieldChecks(MaxANArgs+ANArgsDefAllocInc))
            AlphFieldChecks=TempAFC
            DEALLOCATE(TempAFC)
            ALLOCATE(AlphFieldDefaults(MaxANArgs+ANArgsDefAllocInc))
            AlphFieldDefaults=TempAFD
            DEALLOCATE(TempAFD)
            MaxANArgs=MaxANArgs+ANArgsDefAllocInc
          ENDIF

          TargetChar=InputLine(CurPos+Pos-1:CurPos+Pos-1)

          IF (TargetChar == 'A' .or. TargetChar == 'a') THEN
            AlphaorNumeric(Count)=.true.
            ObjectDef(NumObjectDefs)%NumAlpha=ObjectDef(NumObjectDefs)%NumAlpha+1
            IF (FieldSet) AlphFieldChecks(ObjectDef(NumObjectDefs)%NumAlpha)=CurrentFieldName
            IF (ObjectDef(NumObjectDefs)%NumAlpha == 1) THEN
              IF (INDEX(MakeUpperCase(CurrentFieldName),'NAME') /= 0) ObjectDef(NumObjectDefs)%NameAlpha1=.true.
            ENDIF
          ELSE
            AlphaorNumeric(Count)=.false.
            ObjectDef(NumObjectDefs)%NumNumeric=ObjectDef(NumObjectDefs)%NumNumeric+1
            IF (FieldSet) NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%FieldName=CurrentFieldName
          ENDIF

        ELSE
          CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
          MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
          Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
          AutoCalculatable=AutoCalculate,RetainCase=RetainCaseFlag,ErrorsFound=ErrorsFoundFlag)
          IF (.not. AlphaorNumeric(Count)) THEN
            ! only record for numeric fields
            IF (MinMax) THEN
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxChk=.true.
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%FieldNumber=Count
              IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%WhichMinMax(1)=WhichMinMax
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxString(1)=MinMaxString
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxValue(1)=Value
              ELSE
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%WhichMinMax(2)=WhichMinMax
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxString(2)=MinMaxString
                NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxValue(2)=Value
              ENDIF
            ENDIF   ! End Min/Max
            IF (Default) THEN
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefaultChk=.true.
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%Default=Value
              IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefAutoSize=.true.
              IF (AlphDefaultString == 'AUTOCALCULATE')  NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefAutoCalculate=.true.
            ENDIF
            IF (AutoSize) THEN
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoSizable=.true.
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoSizeValue=Value
            ENDIF
            IF (AutoCalculate) THEN
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoCalculatable=.true.
              NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoCalculateValue=Value
            ENDIF
          ELSE  ! Alpha Field
            IF (Default) THEN
              AlphFieldDefaults(ObjectDef(NumObjectDefs)%NumAlpha)=AlphDefaultString
            ENDIF
          ENDIF
          IF (ErrorsFoundFlag) THEN
            ErrFlag=.true.
            ErrorsFoundFlag=.false.
          ENDIF
          IF (RequiredField) THEN
            RequiredFields(Count)=.true.
            MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
          ENDIF
          IF (RetainCaseFlag) THEN
            AlphRetainCase(Count)=.true.
          ENDIF
          CYCLE
        ENDIF

        !  For the moment dont care about descriptions on each object
        IF (CurPos <= InputLineLength) THEN
          CurPos=CurPos+Pos
          Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
          IF (Pos == 0) THEN
            CALL ShowSevereError('IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' , or ; expected on this line'//  &
            ',position="'//InputLine(CurPos:InputLineLength)//'"',EchoInputFile)
            ErrFlag=.true.
            ErrorsFound=.true.
          ENDIF
          IF (InputLine(InputLineLength:InputLineLength) /= '\') THEN
            !CALL ShowWarningError('IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' \ expected on this line',EchoInputFile)    !RS: Secret Search String
            IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 13) THEN
              WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
            END IF
            WRITE(DebugFile,*) 'IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' \ expected on this line',EchoInputFile
          ENDIF
        ELSE
          CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
          IF (BlankLine .or. EndofFile) CYCLE
          Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
        ENDIF
      ELSE
        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
        CYCLE
      ENDIF

      IF (Pos <= 0) THEN
        ! must be time to read another line
        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
        IF (BlankLine .or. EndofFile) CYCLE
      ELSE
        IF (InputLine(CurPos+Pos-1:CurPos+Pos-1) == ';') THEN
          EndofObjectDef=.true.
        ENDIF
        CurPos=CurPos+Pos
      ENDIF

    END DO

    ! Reached end of object def but there may still be more \ lines to parse....
    ! Goes until next object is encountered ("not blankline") or end of IDDFile
    ! If last object is not numeric, then exit immediately....
    BlankLine=.true.
    DO WHILE (BlankLine .and. .not.EndofFile)
      ! It's a numeric object as last one...
      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
      MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
      Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
      AutoCalculatable=AutoCalculate,RetainCase=RetainCaseFlag,ErrorsFound=ErrorsFoundFlag)
      IF (MinMax) THEN
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxChk=.true.
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%FieldNumber=Count
        IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%WhichMinMax(1)=WhichMinMax
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxString(1)=MinMaxString
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxValue(1)=Value
        ELSE
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%WhichMinMax(2)=WhichMinMax
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxString(2)=MinMaxString
          NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%MinMaxValue(2)=Value
        ENDIF
      ENDIF
      IF (Default .and. .not. AlphaorNumeric(Count)) THEN
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefaultChk=.true.
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%Default=Value
        IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefAutoSize=.true.
        IF (AlphDefaultString == 'AUTOCALCULATE') NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%DefAutoCalculate=.true.
      ELSEIF (Default .and. AlphaorNumeric(Count)) THEN
        AlphFieldDefaults(ObjectDef(NumObjectDefs)%NumAlpha)=AlphDefaultString
      ENDIF
      IF (AutoSize) THEN
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoSizable=.true.
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoSizeValue=Value
      ENDIF
      IF (AutoCalculate) THEN
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoCalculatable=.true.
        NumRangeChecks(ObjectDef(NumObjectDefs)%NumNumeric)%AutoCalculateValue=Value
      ENDIF
      IF (ErrorsFoundFlag) THEN
        ErrFlag=.true.
        ErrorsFoundFlag=.false.
      ENDIF
    ENDDO

    IF (.not. BlankLine) THEN
      BACKSPACE(Unit=IDDFile)
      EchoInputLine=.false.
    ENDIF

    IF (RequiredField) THEN
      RequiredFields(Count)=.true.
      MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
    ENDIF
    IF (RetainCaseFlag) THEN
      AlphRetainCase(Count)=.true.
    ENDIF

    ObjectDef(NumObjectDefs)%NumParams=Count  ! Also the total of ObjectDef(..)%NumAlpha+ObjectDef(..)%NumNumeric
    ObjectDef(NumObjectDefs)%MinNumFields=MinimumNumberOfFields
    IF (ObsoleteObject) THEN
      ALLOCATE(TempAFD(NumObsoleteObjects+1))
      IF (NumObsoleteObjects > 0) THEN
        TempAFD(1:NumObsoleteObjects)=ObsoleteObjectsRepNames
      ENDIF
      TempAFD(NumObsoleteObjects+1)=ReplacementName
      DEALLOCATE(ObsoleteObjectsRepNames)
      NumObsoleteObjects=NumObsoleteObjects+1
      ALLOCATE(ObsoleteObjectsRepNames(NumObsoleteObjects))
      ObsoleteObjectsRepNames=TempAFD
      ObjectDef(NumObjectDefs)%ObsPtr=NumObsoleteObjects
      DEALLOCATE(TempAFD)
    ENDIF
    IF (RequiredObject) THEN
      ObjectDef(NumObjectDefs)%RequiredObject=.true.
    ENDIF
    IF (UniqueObject) THEN
      ObjectDef(NumObjectDefs)%UniqueObject=.true.
    ENDIF
    IF (ExtensibleObject) THEN
      ObjectDef(NumObjectDefs)%ExtensibleObject=.true.
      ObjectDef(NumObjectDefs)%ExtensibleNum=ExtensibleNumFields
    ENDIF

    MaxAlphaArgsFound=MAX(MaxAlphaArgsFound,ObjectDef(NumObjectDefs)%NumAlpha)
    MaxNumericArgsFound=MAX(MaxNumericArgsFound,ObjectDef(NumObjectDefs)%NumNumeric)
    ALLOCATE(ObjectDef(NumObjectDefs)%AlphaorNumeric(Count))
    ObjectDef(NumObjectDefs)%AlphaorNumeric=AlphaorNumeric(1:Count)
    ALLOCATE(ObjectDef(NumObjectDefs)%AlphRetainCase(Count))
    ObjectDef(NumObjectDefs)%AlphRetainCase=AlphRetainCase(1:Count)
    PrevCount = Count
    ALLOCATE(ObjectDef(NumObjectDefs)%NumRangeChks(ObjectDef(NumObjectDefs)%NumNumeric))
    IF (ObjectDef(NumObjectDefs)%NumNumeric > 0) THEN
      ObjectDef(NumObjectDefs)%NumRangeChks=NumRangeChecks(1:ObjectDef(NumObjectDefs)%NumNumeric)
    ENDIF
    PrevSizeNumNumeric = ObjectDef(NumObjectDefs)%NumNumeric !used to clear only portion of NumRangeChecks array
    ALLOCATE(ObjectDef(NumObjectDefs)%AlphFieldChks(ObjectDef(NumObjectDefs)%NumAlpha))
    IF (ObjectDef(NumObjectDefs)%NumAlpha > 0) THEN
      ObjectDef(NumObjectDefs)%AlphFieldChks=AlphFieldChecks(1:ObjectDef(NumObjectDefs)%NumAlpha)
    ENDIF
    ALLOCATE(ObjectDef(NumObjectDefs)%AlphFieldDefs(ObjectDef(NumObjectDefs)%NumAlpha))
    IF (ObjectDef(NumObjectDefs)%NumAlpha > 0) THEN
      ObjectDef(NumObjectDefs)%AlphFieldDefs=AlphFieldDefaults(1:ObjectDef(NumObjectDefs)%NumAlpha)
    ENDIF
    PrevSizeNumAlpha = ObjectDef(NumObjectDefs)%NumAlpha
    ALLOCATE(ObjectDef(NumObjectDefs)%ReqField(Count))
    ObjectDef(NumObjectDefs)%ReqField=RequiredFields(1:Count)
    DO Count=1,ObjectDef(NumObjectDefs)%NumNumeric
      IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxChk) THEN
        ! Checking MinMax Range (min vs. max and vice versa)
        MinMaxError=.false.
        ! check min against max
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
          ! min
          Value=ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)
          IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
            IF (Value > ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
          ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
            IF (Value == ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
          ENDIF
        ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
          ! min>
          Value=ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1) + rTinyValue  ! infintesimally bigger than min
          IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
            IF (Value > ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
          ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
            IF (Value == ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
          ENDIF
        ENDIF
        ! check max against min
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
          ! max
          Value=ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)
          ! Check max value against min
          IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
            IF (Value < ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
          ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
            IF (Value == ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
          ENDIF
        ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
          ! max<
          Value=ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2) - rTinyValue  ! infintesimally bigger than min
          IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
            IF (Value < ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
          ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
            IF (Value == ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
          ENDIF
        ENDIF
        ! check if error condition
        IF (MinMaxError) THEN
          !  Error stated min is not in range with stated max
          MinMaxString=IPTrimSigDigits(ObjectDef(NumObjectDefs)%NumRangeChks(Count)%FieldNumber)
          CALL ShowSevereError('IP: IDD: Field #'//TRIM(MinMaxString)//' conflict in Min/Max specifications/values, in class='//  &
          TRIM(ObjectDef(NumObjectDefs)%Name),EchoInputFile)
          ErrFlag=.true.
        ENDIF
      ENDIF
      IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%DefaultChk) THEN
        ! Check Default against MinMaxRange
        !  Don't check when default is autosize...
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%Autosizable .and.   &
        ObjectDef(NumObjectDefs)%NumRangeChks(Count)%DefAutoSize) CYCLE
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%Autocalculatable .and.   &
        ObjectDef(NumObjectDefs)%NumRangeChks(Count)%DefAutoCalculate) CYCLE
        MinMaxError=.false.
        Value=ObjectDef(NumObjectDefs)%NumRangeChks(Count)%Default
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
          IF (Value < ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
          IF (Value <= ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ENDIF
        IF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
          IF (Value > ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ELSEIF (ObjectDef(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
          IF (Value >= ObjectDef(NumObjectDefs)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ENDIF
        IF (MinMaxError) THEN
          !  Error stated default is not in min/max range
          MinMaxString=IPTrimSigDigits(ObjectDef(NumObjectDefs)%NumRangeChks(Count)%FieldNumber)
          CALL ShowSevereError('IP: IDD: Field #'//TRIM(MinMaxString)//' default is invalid for Min/Max values, in class='//  &
          TRIM(ObjectDef(NumObjectDefs)%Name),EchoInputFile)
          ErrFlag=.true.
        ENDIF
      ENDIF
    ENDDO

    IF (ErrFlag) THEN
      CALL ShowContinueError('IP: Errors occured in ObjectDefinition for Class='//TRIM(ObjectDef(NumObjectDefs)%Name)// &
      ', Object not available for IDF processing.',EchoInputFile)
      DEALLOCATE(ObjectDef(NumObjectDefs)%AlphaorNumeric)
      DEALLOCATE(ObjectDef(NumObjectDefs)%NumRangeChks)
      DEALLOCATE(ObjectDef(NumObjectDefs)%AlphFieldChks)
      DEALLOCATE(ObjectDef(NumObjectDefs)%AlphFieldDefs)
      DEALLOCATE(ObjectDef(NumObjectDefs)%ReqField)
      DEALLOCATE(ObjectDef(NumObjectDefs)%AlphRetainCase)
      NumObjectDefs=NumObjectDefs-1
      ErrorsFound=.true.
    ENDIF

    RETURN

  END SUBROUTINE AddObjectDefandParse

  SUBROUTINE ProcessInputDataFile

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   August 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine processes input data file for EnergyPlus.  Each "record" is
    ! parsed into the LineItem data structure and, if okay, put into the
    ! IDFRecords data structure.

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    TYPE (FileSectionsDefinition), ALLOCATABLE :: TempSectionsonFile(:)   ! Used during reallocation procedure
    TYPE (LineDefinition), ALLOCATABLE :: TempIDFRecords(:)   ! Used during reallocation procedure

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

    LOGICAL :: EndofFile = .false.
    LOGICAL BlankLine
    INTEGER Pos

    INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

    OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

    MaxIDFRecords=ObjectsIDFAllocInc
    NumIDFRecords=0
    MaxIDFSections=SectionsIDFAllocInc
    NumIDFSections=0

    ALLOCATE (SectionsonFile(MaxIDFSections))
    ALLOCATE (IDFRecords(MaxIDFRecords))
    ALLOCATE (LineItem%Numbers(MaxNumericArgsFound))
    ALLOCATE (LineItem%NumBlank(MaxNumericArgsFound))
    ALLOCATE (LineItem%Alphas(MaxAlphaArgsFound))
    ALLOCATE (LineItem%AlphBlank(MaxAlphaArgsFound))
    EndofFile=.false.

    DO WHILE (.not. EndofFile)
      CALL ReadInputLine(IDFFile,Pos,BlankLine,InputLineLength,EndofFile)
      IF (BlankLine .or. EndofFile) CYCLE
      Pos=SCAN(InputLine,',;')
      If (Pos /= 0) then
        If (InputLine(Pos:Pos) == ';') then
          CALL ValidateSection(InputLine(1:Pos-1),NumLines)
          IF (NumIDFSections == MaxIDFSections) THEN
            ALLOCATE (TempSectionsonFile(MaxIDFSections+SectionsIDFAllocInc))
            TempSectionsonFile(1:MaxIDFSections)=SectionsonFile
            DEALLOCATE (SectionsonFile)
            ALLOCATE (SectionsonFile(MaxIDFSections+SectionsIDFAllocInc))
            SectionsonFile=TempSectionsonFile
            DEALLOCATE (TempSectionsonFile)
            MaxIDFSections=MaxIDFSections+SectionsIDFAllocInc
          ENDIF
        else
          CALL ValidateObjectandParse(InputLine(1:Pos-1),Pos,EndofFile)
          IF (NumIDFRecords == MaxIDFRecords) THEN
            ALLOCATE(TempIDFRecords(MaxIDFRecords+ObjectsIDFAllocInc))
            TempIDFRecords(1:MaxIDFRecords)=IDFRecords
            DEALLOCATE(IDFRecords)
            ALLOCATE(IDFRecords(MaxIDFRecords+ObjectsIDFAllocInc))
            IDFRecords=TempIDFRecords
            DEALLOCATE(TempIDFRecords)
            MaxIDFRecords=MaxIDFRecords+ObjectsIDFAllocInc
          ENDIF
        endif
      else
        !Error condition, no , or ; on first line
        CALL ShowMessage('IP: IDF Line~'//TRIM(IPTrimSigDigits(NumLines))//' '//TRIM(InputLine))
        CALL ShowSevereError(', or ; expected on this line',EchoInputFile)
      endif

    END DO

    !IF (NumIDFSections > 0) THEN
    !  SectionsonFile(NumIDFSections)%LastRecord=NumIDFRecords
    !ENDIF

    IF (NumIDFRecords > 0) THEN
      DO Pos=1,NumObjectDefs
        IF (ObjectDef(Pos)%RequiredObject .and. ObjectDef(Pos)%NumFound == 0) THEN
          !         CALL ShowSevereError('IP: No items found for Required Object='//TRIM(ObjectDef(Pos)%Name)) !RS: Debugging: Removing error msg. call so it won't crash
          WRITE(DebugFile,*) 'IP: No items found for Required Object=' //TRIM(ObjectDef(Pos)%Name)   !RS: Secret Search String
          NumMiscErrorsFound=NumMiscErrorsFound+1
        ENDIF
      ENDDO
    ENDIF

    RETURN

  END SUBROUTINE ProcessInputDataFile

  SUBROUTINE ValidateSection(ProposedSection,LineNo)

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   September 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine validates the section from the input data file
    ! with the list of objects from the data dictionary file.

    ! METHODOLOGY EMPLOYED:
    ! A "squeezed" string is formed and checked against the list of
    ! sections.

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    CHARACTER(len=*), INTENT(IN) :: ProposedSection
    INTEGER, INTENT(IN)          :: LineNo

    ! SUBROUTINE PARAMETER DEFINITIONS:
    ! na

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    CHARACTER(len=MaxSectionNameLength) SqueezedSection
    INTEGER Found
    TYPE (SectionsDefinition), ALLOCATABLE :: TempSectionDef(:)  ! Like SectionDef, used during Re-allocation
    INTEGER OFound

    SqueezedSection=MakeUPPERCase(ADJUSTL(ProposedSection))
    IF (LEN_TRIM(ADJUSTL(ProposedSection)) > MaxSectionNameLength) THEN
      CALL ShowWarningError('IP: Section length exceeds maximum, will be truncated='//TRIM(ProposedSection),EchoInputFile)
      CALL ShowContinueError('Will be processed as Section='//TRIM(SqueezedSection),EchoInputFile)
    ENDIF
    IF (SqueezedSection(1:3) /= 'END') THEN
      Found=FindIteminList(SqueezedSection,SectionDef%Name,NumSectionDefs)
      IF (Found == 0) THEN
        ! Make sure this Section not an object name
        IF (SortedIDD) THEN
          OFound=FindItemInSortedList(SqueezedSection,ListOfObjects,NumObjectDefs)
          IF (OFound /= 0) OFound=iListOfObjects(OFound)
        ELSE
          OFound=FindItemInList(SqueezedSection,ListOfObjects,NumObjectDefs)
        ENDIF
        IF (OFound /= 0) THEN
          CALL AddRecordFromSection(OFound)
        ELSEIF (NumSectionDefs == MaxSectionDefs) THEN
          ALLOCATE (TempSectionDef(MaxSectionDefs+SectionDefAllocInc))
          TempSectionDef(1:MaxSectionDefs)=SectionDef
          DEALLOCATE (SectionDef)
          ALLOCATE (SectionDef(MaxSectionDefs+SectionDefAllocInc))
          SectionDef=TempSectionDef
          DEALLOCATE (TempSectionDef)
          MaxSectionDefs=MaxSectionDefs+SectionDefAllocInc
        ENDIF
        NumSectionDefs=NumSectionDefs+1
        SectionDef(NumSectionDefs)%Name=SqueezedSection
        SectionDef(NumSectionDefs)%NumFound=1
        ! Add to "Sections on file" if appropriate
        IF (.not. ProcessingIDD) THEN
          NumIDFSections=NumIDFSections+1
          SectionsonFile(NumIDFSections)%Name=SqueezedSection
          SectionsonFile(NumIDFSections)%FirstRecord=NumIDFRecords
          SectionsonFile(NumIDFSections)%FirstLineNo=LineNo
        ENDIF
      ELSE
        !      IF (NumIDFSections > 0) THEN
        !        SectionsonFile(NumIDFSections)%LastRecord=NumIDFRecords
        !      ENDIF
        SectionDef(Found)%NumFound=SectionDef(Found)%NumFound+1
        IF (.not. ProcessingIDD) THEN
          NumIDFSections=NumIDFSections+1
          SectionsonFile(NumIDFSections)%Name=SqueezedSection
          SectionsonFile(NumIDFSections)%FirstRecord=NumIDFRecords
          SectionsonFile(NumIDFSections)%FirstLineNo=LineNo
        ENDIF
      ENDIF
    ELSE  ! End ...
      IF (.not. ProcessingIDD) THEN
        SqueezedSection=SqueezedSection(4:)
        SqueezedSection=ADJUSTL(SqueezedSection)
        DO Found=NumIDFSections,1,-1
          IF (.not. SameString(SectionsonFile(Found)%Name,SqueezedSection)) CYCLE
          SectionsonFile(Found)%LastRecord=NumIDFRecords
        ENDDO
      ENDIF
    ENDIF

    RETURN

  END SUBROUTINE ValidateSection

  SUBROUTINE ValidateObjectandParse(ProposedObject,CurPos,EndofFile)
    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Linda K. Lawrie
    !       DATE WRITTEN   September 1997
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine validates the proposed object from the IDF and then
    ! parses it, putting it into the internal InputProcessor Data structure.

    ! METHODOLOGY EMPLOYED:
    ! na

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    ! na

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    CHARACTER(len=*), INTENT(IN) :: ProposedObject
    INTEGER, INTENT(INOUT) :: CurPos
    LOGICAL, INTENT(INOUT) :: EndofFile

    ! SUBROUTINE PARAMETER DEFINITIONS:
    INTEGER, PARAMETER :: dimLineBuf=10

    ! INTERFACE BLOCK SPECIFICATIONS
    ! na

    ! DERIVED TYPE DEFINITIONS
    ! na

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    CHARACTER(len=MaxObjectNameLength) SqueezedObject
    CHARACTER(len=MaxAlphaArgLength) SqueezedArg
    INTEGER Found
    INTEGER NumArg
    INTEGER NumArgExpected
    INTEGER NumAlpha
    INTEGER NumNumeric
    INTEGER Pos
    LOGICAL EndofObject
    LOGICAL BlankLine
    LOGICAL,SAVE  :: ErrFlag=.false.
    INTEGER LenLeft
    INTEGER Count
    CHARACTER(len=32) FieldString
    CHARACTER(len=MaxFieldNameLength) FieldNameString
    CHARACTER(len=300) Message
    CHARACTER(len=300) cStartLine
    CHARACTER(len=300) cStartName
    CHARACTER(len=300), DIMENSION(dimLineBuf), SAVE :: LineBuf
    INTEGER, SAVE :: StartLine
    INTEGER, SAVE :: NumConxLines
    INTEGER, SAVE :: CurLines
    INTEGER, SAVE :: CurQPtr

    CHARACTER(len=52) :: String
    LOGICAL IDidntMeanIt
    LOGICAL TestingObject
    LOGICAL TransitionDefer
    INTEGER TFound
    INTEGER, EXTERNAL :: FindNonSpace
    INTEGER NextChr
    CHARACTER(len=32) :: String1

    INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

    !OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

    SqueezedObject=MakeUPPERCase(ADJUSTL(ProposedObject))
    IF (LEN_TRIM(ADJUSTL(ProposedObject)) > MaxObjectNameLength) THEN
      CALL ShowWarningError('IP: Object name length exceeds maximum, will be truncated='//TRIM(ProposedObject),EchoInputFile)
      CALL ShowContinueError('Will be processed as Object='//TRIM(SqueezedObject),EchoInputFile)
    ENDIF
    IDidntMeanIt=.false.

    TestingObject=.true.
    TransitionDefer=.false.
    DO WHILE (TestingObject)
      ErrFlag=.false.
      IDidntMeanIt=.false.
      IF (SortedIDD) THEN
        Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
        IF (Found /= 0) Found=iListofObjects(Found)
      ELSE
        Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
      ENDIF
      IF (Found /= 0) THEN
        IF (ObjectDef(Found)%ObsPtr > 0) THEN
          TFound=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
          IF (TFound /= 0) THEN
            IF (RepObjects(TFound)%Transitioned) THEN
              IF (.not. RepObjects(TFound)%Used)  &
              CALL ShowWarningError('IP: Objects="'//TRIM(ADJUSTL(ProposedObject))//  &
              '" are being transitioned to this object="'//  &
              TRIM(RepObjects(TFound)%NewName)//'"')
              RepObjects(TFound)%Used=.true.
              IF (SortedIDD) THEN
                Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
                IF (Found /= 0) Found=iListofObjects(Found)
              ELSE
                Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
              ENDIF
            ELSEIF (RepObjects(TFound)%TransitionDefer) THEN
              IF (.not. RepObjects(TFound)%Used)  &
              CALL ShowWarningError('IP: Objects="'//TRIM(ADJUSTL(ProposedObject))//  &
              '" are being transitioned to this object="'//  &
              TRIM(RepObjects(TFound)%NewName)//'"')
              RepObjects(TFound)%Used=.true.
              IF (SortedIDD) THEN
                Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
                IF (Found /= 0) Found=iListofObjects(Found)
              ELSE
                Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
              ENDIF
              TransitionDefer=.true.
            ELSE
              Found=0    ! being handled differently for this obsolete object
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      TestingObject=.false.
      IF (Found == 0) THEN
        ! Check to see if it's a "secret" object
        Found=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
        IF (Found == 0) THEN
          !CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          !   ' Did not find "'//TRIM(ADJUSTL(ProposedObject))//'" in list of Objects',EchoInputFile) !RS: Secret Search String
          IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
            WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
          END IF
          WRITE(DebugFile,*) 'IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//' Did not find "'&
          //TRIM(ADJUSTL(ProposedObject))//'" in list of Objects'
          ! Will need to parse to next ;
          ErrFlag=.true.
        ELSEIF (RepObjects(Found)%Deleted) THEN
          IF (.not. RepObjects(Found)%Used) THEN
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" have been deleted from the IDD.  Will be ignored.')
            RepObjects(Found)%Used=.true.
          ENDIF
          IDidntMeanIt=.true.
          ErrFlag=.true.
          Found=0
        ELSEIF (RepObjects(Found)%TransitionDefer) THEN

        ELSE ! This name is replaced with something else
          IF (.not. RepObjects(Found)%Used) THEN
            IF (.not. RepObjects(Found)%Transitioned) THEN
              CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being replaced with this object="'//  &
              TRIM(RepObjects(Found)%NewName)//'"')
              RepObjects(Found)%Used=.true.
              SqueezedObject=RepObjects(Found)%NewName
              TestingObject=.true.
            ELSE
              CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being transitioned to this object="'//  &
              TRIM(RepObjects(Found)%NewName)//'"')
              RepObjects(Found)%Used=.true.
              IF (SortedIDD) THEN
                Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
                IF (Found /= 0) Found=iListofObjects(Found)
              ELSE
                Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
              ENDIF
            ENDIF
          ELSEIF (.not. RepObjects(Found)%Transitioned) THEN
            SqueezedObject=RepObjects(Found)%NewName
            TestingObject=.true.
          ELSE
            IF (SortedIDD) THEN
              Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
              IF (Found /= 0) Found=iListofObjects(Found)
            ELSE
              Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
            ENDIF
          ENDIF
        ENDIF
      ELSE

        ! Start Parsing the Object according to definition

        ErrFlag=.false.
        LineItem%Name=SqueezedObject
        LineItem%Alphas=Blank
        LineItem%AlphBlank=.false.
        LineItem%NumAlphas=0
        LineItem%Numbers=0.0
        LineItem%NumNumbers=0
        LineItem%NumBlank=.false.
        LineItem%ObjectDefPtr=Found
        NumArgExpected=ObjectDef(Found)%NumParams
        ObjectDef(Found)%NumFound=ObjectDef(Found)%NumFound+1
        IF (ObjectDef(Found)%UniqueObject .and. ObjectDef(Found)%NumFound > 1) THEN
          CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' Multiple occurrences of Unique Object='//TRIM(ADJUSTL(ProposedObject)))
          NumMiscErrorsFound=NumMiscErrorsFound+1
        ENDIF
        IF (ObjectDef(Found)%ObsPtr > 0) THEN
          TFound=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
          IF (TFound == 0) THEN
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Obsolete object='//TRIM(ADJUSTL(ProposedObject))//  &
            ', encountered.  Should be replaced with new object='//  &
            TRIM(ObsoleteObjectsRepNames(ObjectDef(Found)%ObsPtr)))
          ELSEIF (.not. RepObjects(TFound)%Used .and. RepObjects(TFound)%Transitioned) THEN
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being transitioned to this object="'//  &
            TRIM(RepObjects(TFound)%NewName)//'"')
            RepObjects(TFound)%Used=.true.
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    NumArg=0
    NumAlpha=0
    NumNumeric=0
    EndofObject=.false.
    CurPos=CurPos+1

    !  Keep context buffer in case of errors
    LineBuf=Blank
    NumConxLines=0
    StartLine=NumLines
    cStartLine=InputLine(1:300)
    cStartName=Blank
    NumConxLines=0
    CurLines=NumLines
    CurQPtr=0

    DO WHILE (.not. EndofFile .and. .not. EndofObject)
      IF (CurLines /= NumLines) THEN
        NumConxLines=MIN(NumConxLines+1,dimLineBuf)
        CurQPtr=CurQPtr+1
        IF (CurQPtr == 1 .and. cStartName == Blank .and. InputLine /= Blank) THEN
          IF (Found > 0) THEN
            IF (ObjectDef(Found)%NameAlpha1) THEN
              Pos=INDEX(InputLine,',')
              cStartName=InputLine(1:Pos-1)
              cStartName=ADJUSTL(cStartName)
            ENDIF
          ENDIF
        ENDIF
        IF (CurQPtr > dimLineBuf) CurQPtr=1
        LineBuf(CurQPtr)=InputLine(1:300)
        CurLines=NumLines
      ENDIF
      IF (CurPos <= InputLineLength) THEN
        Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
        IF (Pos == 0) THEN
          IF (InputLine(InputLineLength:InputLineLength) == '!') THEN
            LenLeft=LEN_TRIM(InputLine(CurPos:InputLineLength-1))
          ELSE
            LenLeft=LEN_TRIM(InputLine(CurPos:InputLineLength))
          ENDIF
          IF (LenLeft == 0) THEN
            CurPos=InputLineLength+1
            CYCLE
          ELSE
            IF (InputLine(InputLineLength:InputLineLength) == '!') THEN
              Pos=InputLineLength-CurPos+1
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Comma being inserted after:"'//InputLine(CurPos:InputLineLength-1)//   &
              '" in Object='//TRIM(SqueezedObject),EchoInputFile)
            ELSE
              Pos=InputLineLength-CurPos+2
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Comma being inserted after:"'//InputLine(CurPos:InputLineLength)// &
              '" in Object='//TRIM(SqueezedObject),EchoInputFile)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        CALL ReadInputLine(IDFFile,CurPos,BlankLine,InputLineLength,EndofFile)
        CYCLE
      ENDIF
      IF (Pos > 0) THEN
        IF (.not. ErrFlag) THEN
          IF (CurPos <= CurPos+Pos-2) THEN
            SqueezedArg=MakeUPPERCase(ADJUSTL(InputLine(CurPos:CurPos+Pos-2)))
            IF (LEN_TRIM(ADJUSTL(InputLine(CurPos:CurPos+Pos-2))) > MaxAlphaArgLength) THEN
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Alpha Argument length exceeds maximum, will be truncated='// &
              TRIM(InputLine(CurPos:CurPos+Pos-2)), EchoInputFile)
              CALL ShowContinueError('Will be processed as Alpha='//TRIM(SqueezedArg),EchoInputFile)
            ENDIF
          ELSE
            SqueezedArg=Blank
          ENDIF
          IF (NumArg == NumArgExpected .and. .not. ObjectDef(Found)%ExtensibleObject) THEN
            CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Error detected for Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
            CALL ShowContinueError(' Maximum arguments reached for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
            EchoInputFile)
            ErrFlag=.true.
          ELSE
            IF (NumArg == NumArgExpected .and. ObjectDef(Found)%ExtensibleObject) THEN
              CALL ExtendObjectDefinition(Found,NumArgExpected)
            ENDIF
            NumArg=NumArg+1
            IF (ObjectDef(Found)%AlphaorNumeric(NumArg)) THEN
              IF (NumAlpha == ObjectDef(Found)%NumAlpha) THEN
                CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                ' Error detected for Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
                CALL ShowContinueError(' Too many Alphas for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
                EchoInputFile)
                ErrFlag=.true.
              ELSE
                NumAlpha=NumAlpha+1
                LineItem%NumAlphas=NumAlpha
                IF (ObjectDef(Found)%AlphRetainCase(NumArg)) THEN
                  SqueezedArg=InputLine(CurPos:CurPos+Pos-2)
                  SqueezedArg=ADJUSTL(SqueezedArg)
                ENDIF
                IF (SqueezedArg /= Blank) THEN
                  LineItem%Alphas(NumAlpha)=SqueezedArg
                ELSEIF (ObjectDef(Found)%ReqField(NumArg)) THEN  ! Blank Argument
                  IF (ObjectDef(Found)%AlphFieldDefs(NumAlpha) /= Blank) THEN
                    LineItem%Alphas(NumAlpha)=ObjectDef(Found)%AlphFieldDefs(NumAlpha)
                  ELSE
                    IF (ObjectDef(Found)%NameAlpha1 .and. NumAlpha /= 1) THEN
                      CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                      CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                      ' Error detected in Object='//TRIM(ObjectDef(Found)%Name)//', name='//  &
                      TRIM(LineItem%Alphas(1)),EchoInputFile)
                    ELSE
                      CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                      CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                      ' Error detected in Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
                    ENDIF
                    CALL ShowContinueError('Field ['//TRIM(ObjectDef(Found)%AlphFieldChks(NumAlpha))//  &
                    '] is required but was blank',EchoInputFile)
                    NumBlankReqFieldFound=NumBlankReqFieldFound+1
                  ENDIF
                ELSE
                  LineItem%AlphBlank(NumAlpha)=.true.
                  IF (ObjectDef(Found)%AlphFieldDefs(NumAlpha) /= Blank) THEN
                    LineItem%Alphas(NumAlpha)=ObjectDef(Found)%AlphFieldDefs(NumAlpha)
                  ENDIF
                ENDIF
              ENDIF
            ELSE
              IF (NumNumeric == ObjectDef(Found)%NumNumeric) THEN
                CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                ' Error detected for Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
                CALL ShowContinueError(' Too many Numbers for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
                EchoInputFile)
                ErrFlag=.true.
              ELSE
                NumNumeric=NumNumeric+1
                LineItem%NumNumbers=NumNumeric
                IF (SqueezedArg /= Blank) THEN
                  IF (.not. ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoSizable .and.   &
                  .not. ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoCalculatable) THEN
                  LineItem%Numbers(NumNumeric)=ProcessNumber(SqueezedArg,Errflag)
                ELSEIF (SqueezedArg == 'AUTOSIZE') THEN
                  LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
                ELSEIF (SqueezedArg == 'AUTOCALCULATE') THEN
                  LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
                ELSE
                  LineItem%Numbers(NumNumeric)=ProcessNumber(SqueezedArg,Errflag)
                ENDIF
              ELSE  ! numeric arg is blank.
                IF (ObjectDef(Found)%NumRangeChks(NumNumeric)%DefaultChk) THEN  ! blank arg has default
                  IF (.not. ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoSize .and.   &
                  .not. ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoCalculatable) THEN
                  LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%Default
                  LineItem%NumBlank(NumNumeric)=.true.
                ELSEIF (ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoSize) THEN
                  LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
                  LineItem%NumBlank(NumNumeric)=.true.
                ELSEIF (ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
                  LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
                  LineItem%NumBlank(NumNumeric)=.true.
                ENDIF
                ErrFlag=.false.
              ELSE ! blank arg does not have default
                IF (ObjectDef(Found)%ReqField(NumArg)) THEN  ! arg is required
                  IF (ObjectDef(Found)%NameAlpha1) THEN  ! object has name field - more context for error
                    CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                    CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                    ' Error detected in Object='//TRIM(ObjectDef(Found)%Name)// &
                    ', name='//TRIM(LineItem%Alphas(1)),EchoInputFile)
                    ErrFlag=.true.
                  ELSE  ! object does not have name field
                    CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                    CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                    ' Error detected in Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
                    ErrFlag=.true.
                  ENDIF
                  CALL ShowContinueError('Field ['//TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName)//  &
                  '] is required but was blank',EchoInputFile)
                  NumBlankReqFieldFound=NumBlankReqFieldFound+1
                ENDIF
                LineItem%Numbers(NumNumeric)=0.0
                LineItem%NumBlank(NumNumeric)=.true.
                !LineItem%Numbers(NumNumeric)=-999999.  !0.0
                !CALL ShowWarningError('Default number in Input, in object='//TRIM(ObjectDef(Found)%Name))
              ENDIF
            ENDIF
            IF (ErrFlag) THEN
              IF (SqueezedArg(1:1) /= '=') THEN  ! argument does not start with "=" (parametric)
                FieldString=IPTrimSigDigits(NumNumeric)
                FieldNameString=ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName
                IF (FieldNameString /= Blank) THEN
                  Message='Invalid Number in Numeric Field#'//TRIM(FieldString)//' ('//TRIM(FieldNameString)//  &
                  '), value='//TRIM(SqueezedArg)
                ELSE ! Field Name not recorded
                  Message='Invalid Number in Numeric Field#'//TRIM(FieldString)//', value='//TRIM(SqueezedArg)
                ENDIF
                Message=TRIM(Message)//', in '//TRIM(ObjectDef(Found)%Name)
                IF (ObjectDef(Found)%NameAlpha1) THEN
                  Message=TRIM(Message)//'='//TRIM(LineItem%Alphas(1))
                ENDIF
                CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                ' '//TRIM(Message),EchoInputFile)
              ELSE  ! parametric in Numeric field
                ErrFlag=.false.
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF

    IF (InputLine(CurPos+Pos-1:CurPos+Pos-1) == ';') THEN
      EndofObject=.true.
      ! Check if more characters on line -- and if first is a comment character
      IF (InputLine(CurPos+Pos:) /= Blank) THEN
        NextChr=FindNonSpace(InputLine(CurPos+Pos:))
        IF (InputLine(CurPos+Pos+NextChr-1:CurPos+Pos+NextChr-1) /= '!') THEN
          CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
          CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' End of Object="'//TRIM(ObjectDef(Found)%Name)//  &
          '" reached, but next part of line not comment.',EchoInputFile)
          CALL ShowContinueError('Final line above shows line that contains error.')
        ENDIF
      ENDIF
    ENDIF
    CurPos=CurPos+Pos
  ENDIF

END DO

! Store to IDFRecord Data Structure, ErrFlag is true if there was an error
! Check out MinimumNumberOfFields
IF (.not. ErrFlag .and. .not. IDidntMeanIt) THEN
  IF (NumArg < ObjectDef(Found)%MinNumFields) THEN
    IF (ObjectDef(Found)%NameAlpha1) THEN
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef(Found)%Name)//  &
      ', name='//TRIM(LineItem%Alphas(1))//       &
      ', entered with less than minimum number of fields.')
    ELSE
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef(Found)%Name)//  &
      ', entered with less than minimum number of fields.')
    ENDIF
    CALL ShowAuditErrorMessage(' **   ~~~   ** ','Attempting fill to minimum.')
    NumAlpha=0
    NumNumeric=0
    IF (ObjectDef(Found)%MinNumFields > ObjectDef(Found)%NumParams) THEN
      String=IPTrimSigDigits(ObjectDef(Found)%MinNumFields)
      String1=IPTrimSigDigits(ObjectDef(Found)%NumParams)
      CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object \min-fields > number of fields specified, Object='//TRIM(ObjectDef(Found)%Name))
      CALL ShowContinueError('..\min-fields='//TRIM(String)//  &
      ', total number of fields in object definition='//TRIM(String1))
      ErrFlag=.true.
    ELSE
      DO Count=1,ObjectDef(Found)%MinNumFields
        IF (ObjectDef(Found)%AlphaOrNumeric(Count)) THEN
          NumAlpha=NumAlpha+1
          IF (NumAlpha <= LineItem%NumAlphas) CYCLE
          LineItem%NumAlphas=LineItem%NumAlphas+1
          IF (ObjectDef(Found)%AlphFieldDefs(LineItem%NumAlphas) /= Blank) THEN
            LineItem%Alphas(LineItem%NumAlphas)=ObjectDef(Found)%AlphFieldDefs(LineItem%NumAlphas)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(ObjectDef(Found)%AlphFieldDefs(LineItem%NumAlphas))//   &
            '   ! field=>'//TRIM(ObjectDef(Found)%AlphFieldChks(NumAlpha)))
          ELSEIF (ObjectDef(Found)%ReqField(Count)) THEN
            IF (ObjectDef(Found)%NameAlpha1) THEN
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef(Found)%Name)//  &
              ', name='//TRIM(LineItem%Alphas(1))// &
              ', Required Field=['//  &
              TRIM(ObjectDef(Found)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ELSE
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef(Found)%Name)//  &
              ', Required Field=['//  &
              TRIM(ObjectDef(Found)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ENDIF
            ErrFlag=.true.
          ELSE
            LineItem%Alphas(LineItem%NumAlphas)=Blank
            LineItem%AlphBlank(LineItem%NumAlphas)=.true.
            CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
            TRIM(ObjectDef(Found)%AlphFieldChks(NumAlpha)))
          ENDIF
        ELSE
          NumNumeric=NumNumeric+1
          IF (NumNumeric <= LineItem%NumNumbers) CYCLE
          LineItem%NumNumbers=LineItem%NumNumbers+1
          LineItem%NumBlank(NumNumeric)=.true.
          IF (ObjectDef(Found)%NumRangeChks(NumNumeric)%Defaultchk) THEN
            IF (.not. ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoSize .and.   &
            .not. ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%Default
            WRITE(String,*) ObjectDef(Found)%NumRangeChks(NumNumeric)%Default
            String=ADJUSTL(String)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(String)//  &
            '   ! field=>'//TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoSize) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autosize    ! field=>'//  &
            TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autocalculate    ! field=>'//  &
            TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName))
          ENDIF
        ELSEIF (ObjectDef(Found)%ReqField(Count)) THEN
          IF (ObjectDef(Found)%NameAlpha1) THEN
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef(Found)%Name)//  &
            ', name='//TRIM(LineItem%Alphas(1))// &
            ', Required Field=['//  &
            TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ELSE
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef(Found)%Name)//  &
            ', Required Field=['//  &
            TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ENDIF
          ErrFlag=.true.
        ELSE
          LineItem%Numbers(NumNumeric)=0.0
          LineItem%NumBlank(NumNumeric)=.true.
          CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
          TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName))
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF
ENDIF

IF (.not. ErrFlag .and. .not. IDidntMeanIt) THEN
  IF (TransitionDefer) THEN
    CALL MakeTransition(Found)
  ENDIF
  NumIDFRecords=NumIDFRecords+1
  IF (ObjectStartRecord(Found) == 0) ObjectStartRecord(Found)=NumIDFRecords
  MaxAlphaIDFArgsFound=MAX(MaxAlphaIDFArgsFound,LineItem%NumAlphas)
  MaxNumericIDFArgsFound=MAX(MaxNumericIDFArgsFound,LineItem%NumNumbers)
  MaxAlphaIDFDefArgsFound=MAX(MaxAlphaIDFDefArgsFound,ObjectDef(Found)%NumAlpha)
  MaxNumericIDFDefArgsFound=MAX(MaxNumericIDFDefArgsFound,ObjectDef(Found)%NumNumeric)
  IDFRecords(NumIDFRecords)%Name=LineItem%Name
  IDFRecords(NumIDFRecords)%NumNumbers=LineItem%NumNumbers
  IDFRecords(NumIDFRecords)%NumAlphas=LineItem%NumAlphas
  IDFRecords(NumIDFRecords)%ObjectDefPtr=LineItem%ObjectDefPtr
  ALLOCATE(IDFRecords(NumIDFRecords)%Alphas(LineItem%NumAlphas))
  ALLOCATE(IDFRecords(NumIDFRecords)%AlphBlank(LineItem%NumAlphas))
  ALLOCATE(IDFRecords(NumIDFRecords)%Numbers(LineItem%NumNumbers))
  ALLOCATE(IDFRecords(NumIDFRecords)%NumBlank(LineItem%NumNumbers))
  IDFRecords(NumIDFRecords)%Alphas(1:LineItem%NumAlphas)=LineItem%Alphas(1:LineItem%NumAlphas)
  IDFRecords(NumIDFRecords)%AlphBlank(1:LineItem%NumAlphas)=LineItem%AlphBlank(1:LineItem%NumAlphas)
  IDFRecords(NumIDFRecords)%Numbers(1:LineItem%NumNumbers)=LineItem%Numbers(1:LineItem%NumNumbers)
  IDFRecords(NumIDFRecords)%NumBlank(1:LineItem%NumNumbers)=LineItem%NumBlank(1:LineItem%NumNumbers)
  IF (LineItem%NumNumbers > 0) THEN
    DO Count=1,LineItem%NumNumbers
      IF (ObjectDef(Found)%NumRangeChks(Count)%MinMaxChk .and. .not. LineItem%NumBlank(Count)) THEN
        CALL InternalRangeCheck(LineItem%Numbers(Count),Count,Found,LineItem%Alphas(1),  &
        ObjectDef(Found)%NumRangeChks(Count)%AutoSizable,        &
        ObjectDef(Found)%NumRangeChks(Count)%AutoCalculatable)
      ENDIF
    ENDDO
  ENDIF
ELSEIF (.not. IDidntMeanIt) THEN
  OverallErrorFlag=.true.
ENDIF

RETURN

END SUBROUTINE ValidateObjectandParse

SUBROUTINE ValidateSectionsInput

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine uses the data structure that is set up during
  ! IDF processing and makes sure that record pointers are accurate.
  ! They could be inaccurate if a 'section' is input without any
  ! 'objects' following.  The invalidity will show itself in the
  ! values of the FirstRecord and Last Record pointer.
  ! If FirstRecord>LastRecord, then no records (Objects) have been
  ! written to the SIDF file for that Section.

  ! METHODOLOGY EMPLOYED:
  ! Scan the SectionsonFile data structure and look for invalid
  ! FirstRecord,LastRecord items.  Reset those items to -1.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count

  IF(EchoInputFile .EQ. 9 .OR. EchoInputFile .EQ. 10 .OR. EchoInputFile .EQ. 12) THEN
    WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
  END IF

  DO Count=1,NumIDFSections
    IF (SectionsonFile(Count)%FirstRecord > SectionsonFile(Count)%LastRecord) THEN
      WRITE(EchoInputFile,*) ' Section ',Count,' ',TRIM(SectionsonFile(Count)%Name),' had no object records'
      SectionsonFile(Count)%FirstRecord=-1
      SectionsonFile(Count)%LastRecord=-1
    ENDIF
  END DO

  RETURN

END SUBROUTINE ValidateSectionsInput

INTEGER FUNCTION GetNumSectionsFound(SectionWord)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This function returns the number of a particular section (in input data file)
  ! found in the current run.  If it can't find the section in list
  ! of sections, a -1 will be returned.

  ! METHODOLOGY EMPLOYED:
  ! Look up section in list of sections.  If there, return the
  ! number of sections of that kind found in the current input.  If not, return
  ! -1.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: SectionWord

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Found

  Found=FindIteminList(MakeUPPERCase(SectionWord),ListofSections,NumSectionDefs)
  IF (Found == 0) THEN
    !    CALL ShowFatalError('Requested Section not found in Definitions: '//TRIM(SectionWord))
    GetNumSectionsFound=0
  ELSE
    GetNumSectionsFound=SectionDef(Found)%NumFound
  ENDIF

  RETURN

END FUNCTION GetNumSectionsFound

INTEGER FUNCTION GetNumSectionsinInput()

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This function returns the number of sections in the entire input data file
  ! of the current run.

  ! METHODOLOGY EMPLOYED:
  ! Return value of NumIDFSections.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  GetNumSectionsinInput=NumIDFSections

  RETURN

END FUNCTION GetNumSectionsinInput

INTEGER FUNCTION GetNumObjectsFound(ObjectWord)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This function returns the number of objects (in input data file)
  ! found in the current run.  If it can't find the object in list
  ! of objects, a 0 will be returned.

  ! METHODOLOGY EMPLOYED:
  ! Look up object in list of objects.  If there, return the
  ! number of objects found in the current input.  If not, return 0.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: ObjectWord

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Found

  IF (SortedIDD) THEN
    Found=FindIteminSortedList(MakeUPPERCase(ObjectWord),ListofObjects,NumObjectDefs)
    IF (Found /= 0) Found=iListofObjects(Found)
  ELSE
    Found=FindIteminList(MakeUPPERCase(ObjectWord),ListofObjects,NumObjectDefs)
  ENDIF

  IF (Found /= 0) THEN
    GetNumObjectsFound=ObjectDef(Found)%NumFound
  ELSE
    GetNumObjectsFound=0
    CALL ShowWarningError('Requested Object not found in Definitions: '//TRIM(ObjectWord))
  ENDIF

  RETURN

END FUNCTION GetNumObjectsFound

SUBROUTINE GetRecordLocations(Which,FirstRecord,LastRecord)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine returns the record location values (which will be
  ! passed to 'GetObjectItem') for a section from the list of inputted
  ! sections (sequential).

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: Which
  INTEGER, INTENT(OUT) :: FirstRecord
  INTEGER, INTENT(OUT) :: LastRecord

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  FirstRecord=SectionsonFile(Which)%FirstRecord
  LastRecord=SectionsonFile(Which)%LastRecord

  RETURN

END SUBROUTINE GetRecordLocations

SUBROUTINE GetObjectItem(Object,Number,Alphas,NumAlphas,Numbers,NumNumbers,Status,NumBlank,AlphaBlank,   &
  AlphaFieldNames,NumericFieldNames)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine gets the 'number' 'object' from the IDFRecord data structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: Object
  INTEGER, INTENT(IN) :: Number
  CHARACTER(len=*), INTENT(OUT), DIMENSION(:) :: Alphas
  INTEGER, INTENT(OUT) :: NumAlphas
  REAL(r64), INTENT(OUT), DIMENSION(:) :: Numbers
  !REAL, INTENT(OUT), DIMENSION(:) :: Numbers
  INTEGER, INTENT(OUT) :: NumNumbers
  INTEGER, INTENT(OUT) :: Status
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaBlank
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumBlank
  CHARACTER(len=*), DIMENSION(:), OPTIONAL :: AlphaFieldNames
  CHARACTER(len=*), DIMENSION(:), OPTIONAL :: NumericFieldNames

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count
  INTEGER LoopIndex
  CHARACTER(len=MaxObjectNameLength) ObjectWord
  CHARACTER(len=MaxObjectNameLength) UCObject
  CHARACTER(len=MaxObjectNameLength), SAVE, ALLOCATABLE, DIMENSION(:) :: AlphaArgs
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: NumberArgs
  LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: AlphaArgsBlank
  LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: NumberArgsBlank
  INTEGER MaxAlphas,MaxNumbers
  INTEGER Found
  INTEGER StartRecord
  CHARACTER(len=32) :: cfld1
  CHARACTER(len=32) :: cfld2

  MaxAlphas=SIZE(Alphas,1)
  MaxNumbers=SIZE(Numbers,1)

  IF (.not. ALLOCATED(AlphaArgs)) THEN
    IF (NumObjectDefs == 0) THEN
      CALL ProcessInput
    ENDIF
    ALLOCATE(AlphaArgs(MaxAlphaArgsFound))
    ALLOCATE(NumberArgs(MaxNumericArgsFound))
    ALLOCATE(NumberArgsBlank(MaxNumericArgsFound))
    ALLOCATE(AlphaArgsBlank(MaxAlphaArgsFound))
  ENDIF

  Count=0
  Status=-1
  UCOBject=MakeUPPERCase(Object)
  IF (SortedIDD) THEN
    Found=FindIteminSortedList(UCOBject,ListofObjects,NumObjectDefs)
    IF (Found /= 0) Found=iListofObjects(Found)
  ELSE
    Found=FindIteminList(UCOBject,ListofObjects,NumObjectDefs)
  ENDIF
  IF (Found == 0) THEN   !  This is more of a developer problem
    CALL ShowFatalError('Requested object='//TRIM(UCObject)//', not found in Object Definitions -- incorrect IDD attached.')
  ENDIF

  IF (ObjectDef(Found)%NumAlpha > 0) THEN
    IF (ObjectDef(Found)%NumAlpha > MaxAlphas) THEN
      cfld1=IPTrimSigDigits(ObjectDef(Found)%NumAlpha)
      cfld2=IPTrimSigDigits(MaxAlphas)
      CALL ShowFatalError('GetObjectItem: '//TRIM(Object)//', Number of ObjectDef Alpha Args ['//TRIM(cfld1)//  &
      '] > Size of AlphaArg array ['//TRIM(cfld2)//'].')
    ENDIF
    Alphas(1:ObjectDef(Found)%NumAlpha)=Blank
  ENDIF
  IF (ObjectDef(Found)%NumNumeric > 0) THEN
    IF (ObjectDef(Found)%NumNumeric > MaxNumbers) THEN
      cfld1=IPTrimSigDigits(ObjectDef(Found)%NumNumeric)
      cfld2=IPTrimSigDigits(MaxNumbers)
      CALL ShowFatalError('GetObjectItem: '//TRIM(Object)//', Number of ObjectDef Numeric Args ['//TRIM(cfld1)//  &
      '] > Size of NumericArg array ['//TRIM(cfld2)//'].')
    ENDIF
    Numbers(1:ObjectDef(Found)%NumNumeric)=0.0
  ENDIF

  StartRecord=ObjectStartRecord(Found)
  IF (StartRecord == 0) THEN
    CALL ShowWarningError('Requested object='//TRIM(UCObject)//', not found in IDF.')
    Status=-1
    StartRecord=NumIDFRecords+1
  ENDIF

  IF (ObjectGotCount(Found) == 0) THEN
    WRITE(EchoInputFile,*) 'Getting object=',TRIM(UCObject)
  ENDIF
  ObjectGotCount(Found)=ObjectGotCount(Found)+1

  DO LoopIndex=StartRecord,NumIDFRecords
    IF (IDFRecords(LoopIndex)%Name == UCObject) THEN
      Count=Count+1
      IF (Count == Number) THEN
        IDFRecordsGotten(LoopIndex)=.true.  ! only object level "gets" recorded
        ! Read this one
        CALL GetObjectItemfromFile(LoopIndex,ObjectWord,AlphaArgs,NumAlphas,NumberArgs,NumNumbers,AlphaArgsBlank,NumberArgsBlank)
        IF (NumAlphas > MaxAlphas .or. NumNumbers > MaxNumbers) THEN
          CALL ShowFatalError('Too many actual arguments for those expected on Object: '//TRIM(ObjectWord)//     &
          ' (GetObjectItem)',EchoInputFile)
        ENDIF
        NumAlphas=MIN(MaxAlphas,NumAlphas)
        NumNumbers=MIN(MaxNumbers,NumNumbers)
        IF (NumAlphas > 0) THEN
          Alphas(1:NumAlphas)=AlphaArgs(1:NumAlphas)
        ENDIF
        IF (NumNumbers > 0) THEN
          Numbers(1:NumNumbers)=NumberArgs(1:NumNumbers)
        ENDIF
        IF (PRESENT(NumBlank)) THEN
          NumBlank=.true.
          IF (NumNumbers > 0) &
          NumBlank(1:NumNumbers)=NumberArgsBlank(1:NumNumbers)
        ENDIF
        IF (PRESENT(AlphaBlank)) THEN
          AlphaBlank=.true.
          IF (NumAlphas > 0) &
          AlphaBlank(1:NumAlphas)=AlphaArgsBlank(1:NumAlphas)
        ENDIF
        IF (PRESENT(AlphaFieldNames)) THEN
          AlphaFieldNames(1:ObjectDef(Found)%NumAlpha)=ObjectDef(Found)%AlphFieldChks(1:ObjectDef(Found)%NumAlpha)
        ENDIF
        IF (PRESENT(NumericFieldNames)) THEN
          NumericFieldNames(1:ObjectDef(Found)%NumNumeric)=ObjectDef(Found)%NumRangeChks(1:ObjectDef(Found)%NumNumeric)%FieldName
        ENDIF
        Status=1
        EXIT
      ENDIF
    ENDIF
  END DO


  RETURN

END SUBROUTINE GetObjectItem

INTEGER FUNCTION GetObjectItemNum(ObjType,ObjName)

  ! SUBROUTINE INFORMATION
  !             AUTHOR:  Fred Buhl
  !       DATE WRITTEN:  Jan 1998
  !           MODIFIED:  Lawrie, September 1999. Take advantage of internal
  !                      InputProcessor structures to speed search.
  !      RE-ENGINEERED:  This is new code, not reengineered

  ! PURPOSE OF THIS SUBROUTINE:
  ! Get the occurrence number of an object of type ObjType and name ObjName

  ! METHODOLOGY EMPLOYED:
  ! Use internal IDF record structure for each object occurrence
  ! and compare the name with ObjName.

  ! REFERENCES:
  ! na

  IMPLICIT NONE

  ! SUBROUTINE ARGUMENTS:
  CHARACTER(len=*), INTENT(IN) :: ObjType   ! Object Type (ref: IDD Objects)
  CHARACTER(len=*), INTENT(IN) :: ObjName   ! Name of the object type

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK DEFINITIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DEFINITIONS
  INTEGER                                 :: NumObjOfType ! Total number of Object Type in IDF
  INTEGER                                 :: ObjNum       ! Loop index variable
  INTEGER                                 :: ItemNum      ! Item number for Object Name
  INTEGER                                 :: Found        ! Indicator for Object Type in list of Valid Objects
  CHARACTER(len=MaxObjectNameLength)      :: UCObjType    ! Upper Case for ObjType
  LOGICAL                                 :: ItemFound    ! Set to true if item found
  LOGICAL                                 :: ObjectFound  ! Set to true if object found
  INTEGER                                 :: StartRecord  ! Start record for objects

  ItemNum = 0
  ItemFound=.false.
  ObjectFound=.false.
  UCObjType=MakeUPPERCase(ObjType)
  IF (SortedIDD) THEN
    Found=FindIteminSortedList(UCObjType,ListofObjects,NumObjectDefs)
    IF (Found /= 0) Found=iListofObjects(Found)
  ELSE
    Found=FindIteminList(UCObjType,ListofObjects,NumObjectDefs)
  ENDIF

  IF (Found /= 0) THEN

    ObjectFound=.true.
    NumObjOfType=ObjectDef(Found)%NumFound
    ItemNum=0
    StartRecord=ObjectStartRecord(Found)

    IF (StartRecord > 0) THEN
      DO ObjNum=StartRecord,NumIDFRecords
        IF (IDFRecords(ObjNum)%Name /= UCObjType) CYCLE
        ItemNum=ItemNum+1
        IF (ItemNum > NumObjOfType) EXIT
        IF (IDFRecords(ObjNum)%Alphas(1) == ObjName) THEN
          ItemFound=.true.
          EXIT
        ENDIF
      END DO
    ENDIF
  ENDIF

  IF (ObjectFound) THEN
    IF (.not. ItemFound) ItemNum=0
  ELSE
    ItemNum=-1  ! if object not found, then flag it
  ENDIF

  GetObjectItemNum = ItemNum

  RETURN

END FUNCTION GetObjectItemNum


SUBROUTINE TellMeHowManyObjectItemArgs(Object,Number,NumAlpha,NumNumbers,Status)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine returns the number of arguments (alpha and numeric) for
  ! the referenced 'number' Object.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: Object
  INTEGER, INTENT(IN) :: Number
  INTEGER, INTENT(OUT) :: NumAlpha
  INTEGER, INTENT(OUT) :: NumNumbers
  INTEGER, INTENT(OUT) :: Status


  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count
  INTEGER LoopIndex
  CHARACTER(len=MaxObjectNameLength) ObjectWord

  Count=0
  Status=-1
  DO LoopIndex=1,NumIDFRecords
    IF (SameString(IDFRecords(LoopIndex)%Name,Object)) THEN
      Count=Count+1
      IF (Count == Number) THEN
        ! Read this one
        CALL GetObjectItemfromFile(LoopIndex,ObjectWord,NumAlpha=NumAlpha,NumNumeric=NumNumbers)
        Status=1
        EXIT
      ENDIF
    ENDIF
  END DO


  RETURN

END SUBROUTINE TellMeHowManyObjectItemArgs

SUBROUTINE GetObjectItemfromFile(Which,ObjectWord,AlphaArgs,NumAlpha,NumericArgs,NumNumeric,AlphaBlanks,NumericBlanks)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine "gets" the object instance from the data structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: Which
  CHARACTER(len=*), INTENT(OUT) :: ObjectWord
  CHARACTER(len=*), INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaArgs
  INTEGER, INTENT(OUT) :: NumAlpha
  !REAL(r64), INTENT(OUT), DIMENSION(:), OPTIONAL :: NumericArgs
  REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumericArgs
  INTEGER, INTENT(OUT) :: NumNumeric
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaBlanks
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumericBlanks

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  TYPE (LineDefinition):: xLineItem                        ! Description of current record

  IF (Which > 0 .and. Which <= NumIDFRecords) THEN
    xLineItem=IDFRecords(Which)
    ObjectWord=xLineItem%Name
    NumAlpha=xLineItem%NumAlphas
    NumNumeric=xLineItem%NumNumbers
    IF (PRESENT(AlphaArgs)) THEN
      IF (NumAlpha >=1) THEN
        AlphaArgs(1:NumAlpha)=xLineItem%Alphas(1:NumAlpha)
      ENDIF
    ENDIF
    IF (PRESENT(AlphaBlanks)) THEN
      IF (NumAlpha >=1) THEN
        AlphaBlanks(1:NumAlpha)=xLineItem%AlphBlank(1:NumAlpha)
      ENDIF
    ENDIF
    IF (PRESENT(NumericArgs)) THEN
      IF (NumNumeric >= 1) THEN
        NumericArgs(1:NumNumeric)=xLineItem%Numbers(1:NumNumeric)
      ENDIF
    ENDIF
    IF (PRESENT(NumericBlanks)) THEN
      IF (NumNumeric >= 1) THEN
        NumericBlanks(1:NumNumeric)=xLineItem%NumBlank(1:NumNumeric)
      ENDIF
    ENDIF
  ELSE
    WRITE(EchoInputFile,*) ' Requested Record',Which,' not in range, 1 -- ',NumIDFRecords
  ENDIF

  RETURN

END SUBROUTINE GetObjectItemfromFile

! Utility Functions/Routines for Module

SUBROUTINE ReadInputLine(UnitNumber,CurPos,BlankLine,InputLineLength,EndofFile,  &
  MinMax,WhichMinMax,MinMaxString,Value,Default,DefString,AutoSizable,  &
  AutoCalculatable,RetainCase,ErrorsFound)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine reads a line in the specified file and checks for end of file

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: UnitNumber
  INTEGER, INTENT(INOUT) :: CurPos
  LOGICAL, INTENT(INOUT) :: EndofFile
  LOGICAL, INTENT(INOUT) :: BlankLine
  INTEGER, INTENT(INOUT) :: InputLineLength
  LOGICAL, INTENT(INOUT), OPTIONAL :: MinMax
  INTEGER, INTENT(INOUT), OPTIONAL :: WhichMinMax   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
  CHARACTER(len=*), INTENT(INOUT), OPTIONAL :: MinMaxString
  REAL(r64), INTENT(INOUT), OPTIONAL :: Value
  LOGICAL, INTENT(INOUT), OPTIONAL :: Default
  CHARACTER(len=*), INTENT(INOUT), OPTIONAL :: DefString
  LOGICAL, INTENT(INOUT), OPTIONAL :: AutoSizable
  LOGICAL, INTENT(INOUT), OPTIONAL :: AutoCalculatable
  LOGICAL, INTENT(INOUT), OPTIONAL :: RetainCase
  LOGICAL, INTENT(INOUT), OPTIONAL :: ErrorsFound

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=1), PARAMETER :: TabChar=CHAR(9)

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER ReadStat
  INTEGER Pos
  INTEGER Slash
  INTEGER P1
  CHARACTER(len=MaxInputLineLength) UCInputLine        ! Each line can be up to MaxInputLineLength characters long
  LOGICAL TabsInLine
  INTEGER NSpace
  LOGICAL ErrFlag
  INTEGER, EXTERNAL :: FindNonSpace
  INTEGER ErrLevel
  INTEGER endcol
  CHARACTER(len=52) cNumLines
  LOGICAL LineTooLong

  ErrFlag=.false.
  LineTooLong=.false.

  !IF(UnitNumber .EQ. 9 .OR. UnitNumber .EQ. 10 .OR. UnitNumber .EQ. 12) THEN
  !  WRITE(*,*) 'Error with UnitNumber'    !RS: Debugging: Searching for a mis-set file number
  !END IF

  READ(UnitNumber,fmta,IOSTAT=ReadStat) InputLine

  IF (ReadStat /= 0) InputLine=Blank

  ! Following section of code allows same software to read Win or Unix files without translating
  IF (StripCR) THEN
    endcol=LEN_TRIM(InputLine)
    IF (ICHAR(InputLine(endcol:endcol)) == iASCII_CR) InputLine(endcol:endcol)=Blank
  ENDIF

  IF (InputLine(MaxInputLineLength+1:) /= Blank) THEN
    LineTooLong=.true.
    InputLine=InputLine(1:MaxInputLineLength)
  ENDIF

  P1=SCAN(InputLine,TabChar)
  TabsInLine=.false.
  DO WHILE (P1>0)
    TabsInLine=.true.
    InputLine(P1:P1)=Blank
    P1=SCAN(InputLine,TabChar)
  ENDDO
  BlankLine=.false.
  CurPos=1
  IF (ReadStat < 0) THEN
    EndofFile=.true.
  ELSE
    IF (EchoInputLine) THEN
      NumLines=NumLines+1
      IF (NumLines < 100000) THEN
        WRITE(EchoInputFile,'(2X,I5,1X,A)') NumLines,TRIM(InputLine)
      ELSE
        cNumLines=IPTrimSigDigits(NumLines)
        WRITE(EchoInputFile,'(1X,A,1X,A)') TRIM(cNumLines),TRIM(InputLine)
      ENDIF
      IF (TabsInLine) WRITE(EchoInputFile,"(6X,'***** Tabs eliminated from above line')")
      IF (LineTooLong) THEN
        CALL ShowSevereError('Input line longer than maximum length allowed='//TRIM(IPTrimSigDigits(MaxInputLineLength))//  &
        ' characters. Other errors may follow.')
        CALL ShowContinueError('.. at line='//TRIM(IPTrimSigDigits(NumLines))//', first 50 characters='//  &
        TRIM(InputLine(1:50)))
        WRITE(EchoInputFile,"(6X,'***** Previous line is longer than allowed length for input line')")
      ENDIF
    ENDIF
    EchoInputLine=.true.
    InputLineLength=LEN_TRIM(InputLine)
    IF (InputLineLength == 0) THEN
      BlankLine=.true.
    ENDIF
    IF (ProcessingIDD) THEN
      Pos=SCAN(InputLine,'!\')  ! 4/30/09 remove ~
      Slash=INDEX(InputLine,'\')
    ELSE
      Pos=SCAN(InputLine,'!')  ! 4/30/09 remove ~
      Slash=0
    ENDIF
    IF (Pos /= 0) THEN
      InputLineLength=Pos
      IF (Pos-1 > 0) THEN
        IF (LEN_TRIM(InputLine(1:Pos-1)) == 0) THEN
          BlankLine=.true.
        ENDIF
      ELSE
        BlankLine=.true.
      ENDIF
      IF (Slash /= 0 .and. Pos == Slash) THEN
        UCInputLine=MakeUPPERCase(InputLine)
        IF (UCInputLine(Slash:Slash+5) == '\FIELD') THEN
          ! Capture Field Name
          CurrentFieldName=InputLine(Slash+6:)
          CurrentFieldName=ADJUSTL(CurrentFieldName)
          P1=SCAN(CurrentFieldName,'!')
          IF (P1 /= 0) CurrentFieldName(P1:)=Blank
          FieldSet=.true.
        ELSE
          FieldSet=.false.
        ENDIF
        IF (UCInputLine(Slash:Slash+14) == '\REQUIRED-FIELD') THEN
          RequiredField=.true.
        ENDIF  ! Required-field arg
        IF (UCInputLine(Slash:Slash+15) == '\REQUIRED-OBJECT') THEN
          RequiredObject=.true.
        ENDIF  ! Required-object arg
        IF (UCInputLine(Slash:Slash+13) == '\UNIQUE-OBJECT') THEN
          UniqueObject=.true.
        ENDIF  ! Unique-object arg
        IF (UCInputLine(Slash:Slash+10) == '\EXTENSIBLE') THEN
          ExtensibleObject=.true.
          IF (UCInputLine(Slash+11:Slash+11) /= ':') THEN
            CALL ShowFatalError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Illegal definition for extensible object, should be "\extensible:<num>"',EchoInputFile)
          ELSE
            ! process number
            NSpace=SCAN(UCInputLine(Slash+12:),' !')
            ExtensibleNumFields=INT(ProcessNumber(UCInputLine(Slash+12:Slash+12+NSpace-1),ErrFlag))
            IF (ErrFlag) THEN
              CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Illegal Number for \extensible:<num>',EchoInputFile)
            ENDIF
          ENDIF
        ENDIF  ! Extensible arg
        IF (UCInputLine(Slash:Slash+10) == '\RETAINCASE') THEN
          RetainCase=.true.
        ENDIF  ! Unique-object arg
        IF (UCInputLine(Slash:Slash+10) == '\MIN-FIELDS') THEN
          !              RequiredField=.true.
          NSpace=FindNonSpace(UCInputLine(Slash+11:))
          IF (NSpace == 0) THEN
            CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Need number for \Min-Fields',EchoInputFile)
            ErrFlag=.true.
            MinimumNumberOfFields=0
          ELSE
            Slash=Slash+11+NSpace-1
            NSpace=SCAN(UCInputLine(Slash:),' !')
            MinimumNumberOfFields=INT(ProcessNumber(UCInputLine(Slash:Slash+NSpace-1),ErrFlag))
            IF (ErrFlag) THEN
              CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Illegal Number for \Min-Fields',EchoInputFile)
            ENDIF
          ENDIF
        ENDIF  ! Min-Fields Arg
        IF (UCInputLine(Slash:Slash+9) == '\OBSOLETE') THEN
          NSpace=INDEX(UCInputLine(Slash+9:),'=>')
          IF (NSpace == 0) THEN
            CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Need replacement object for \Obsolete objects',EchoInputFile)
            ErrFlag=.true.
          ELSE
            Slash=Slash+9+NSpace+1
            NSpace=SCAN(UCInputLine(Slash:),'!')
            IF (NSpace == 0) THEN
              ReplacementName=InputLine(Slash:)
            ELSE
              ReplacementName=InputLine(Slash:Slash+NSpace-2)
            ENDIF
            ObsoleteObject=.true.
          ENDIF
        ENDIF  ! Obsolete Arg
        IF (PRESENT(MinMax)) THEN
          IF (UCInputLine(Pos:Pos+7)=='\MINIMUM' .or.  &
          UCInputLine(Pos:Pos+7)=='\MAXIMUM') THEN
          MinMax=.true.
          CALL ProcessMinMaxDefLine(UCInputLine(Pos:),WhichMinMax,MinMaxString,Value,DefString,ErrLevel)
          IF (ErrLevel > 0) THEN
            CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Error in Minimum/Maximum designation -- invalid number='//TRIM(UCInputLine(Pos:)),  &
            EchoInputFile)
            ErrFlag=.true.
          ENDIF
        ELSE
          MinMax=.false.
        ENDIF
      ENDIF  ! Min/Max Args
      IF (PRESENT(Default)) THEN
        IF (UCInputLine(Pos:Pos+7)=='\DEFAULT') THEN
          ! WhichMinMax, MinMaxString not filled here
          Default=.true.
          CALL ProcessMinMaxDefLine(InputLine(Pos:),WhichMinMax,MinMaxString,Value,DefString,ErrLevel)
          IF (.not. RetainCase .and. DefString /= Blank) DefString=MakeUPPERCase(DefString)
          IF (ErrLevel > 1) THEN
            CALL ShowContinueError('Blank Default Field Encountered',EchoInputFile)
            ErrFlag=.true.
          ENDIF
        ELSE
          Default=.false.
        ENDIF
      ENDIF  ! Default Arg
      IF (PRESENT(AutoSizable)) THEN
        IF (UCInputLine(Pos:Pos+5)=='\AUTOS') THEN
          AutoSizable=.true.
          CALL ProcessMinMaxDefLine(UCInputLine(Pos:),WhichMinMax,MinMaxString,Value,DefString,ErrLevel)
          IF (ErrLevel > 0) THEN
            CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Error in Autosize designation -- invalid number='//TRIM(UCInputLine(Pos:)),EchoInputFile)
            ErrFlag=.true.
          ENDIF
        ELSE
          AutoSizable=.false.
        ENDIF
      ENDIF  ! AutoSizable Arg
      IF (PRESENT(AutoCalculatable)) THEN
        IF (UCInputLine(Pos:Pos+5)=='\AUTOC') THEN
          AutoCalculatable=.true.
          CALL ProcessMinMaxDefLine(UCInputLine(Pos:),WhichMinMax,MinMaxString,Value,DefString,ErrLevel)
          IF (ErrLevel > 0) THEN
            CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Error in Autocalculate designation -- invalid number='//  &
            TRIM(UCInputLine(Pos:)),EchoInputFile)
            ErrFlag=.true.
          ENDIF
        ELSE
          AutoCalculatable=.false.
        ENDIF
      ENDIF  ! AutoCalculatable Arg
    ENDIF
  ENDIF
ENDIF
IF (ErrFlag) THEN
  IF (PRESENT(ErrorsFound)) THEN
    ErrorsFound=.true.
  ENDIF
ENDIF

RETURN

END SUBROUTINE ReadInputLine

SUBROUTINE ExtendObjectDefinition(ObjectNum,NumNewArgsLimit)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   Sep 2008
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine expands the object definition according to the extensible "rules" entered
  ! by the developer.  The developer should enter the number of fields to be duplicated.
  ! See References section for examples.

  ! METHODOLOGY EMPLOYED:
  ! The routine determines the type of the fields to be added (A or N) and reallocates the
  ! appropriate arrays in the object definition structure.

  ! REFERENCES:
  ! Extensible objects have a \extensible:<num> specification
  ! \extensible:3 -- the last 3 fields are "extended"
  ! Works on this part of the definition:
  !   INTEGER :: NumParams                       =0   ! Number of parameters to be processed for each object
  !   INTEGER :: NumAlpha                        =0   ! Number of Alpha elements in the object
  !   INTEGER :: NumNumeric                      =0   ! Number of Numeric elements in the object
  !   LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: AlphaorNumeric ! Positionally, whether the argument
  !                                                           ! is alpha (true) or numeric (false)
  !   LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: ReqField ! True for required fields
  !   LOGICAL(1), ALLOCATABLE, DIMENSION(:) :: AlphRetainCase ! true if retaincase is set for this field (alpha fields only)
  !   CHARACTER(len=MaxNameLength+40),  &
  !               ALLOCATABLE, DIMENSION(:) :: AlphFieldChks ! Field names for alphas
  !   CHARACTER(len=MaxNameLength),  &
  !               ALLOCATABLE, DIMENSION(:) :: AlphFieldDefs ! Defaults for alphas
  !   TYPE(RangeCheckDef), ALLOCATABLE, DIMENSION(:) :: NumRangeChks  ! Used to range check and default numeric fields
  !   INTEGER :: LastExtendAlpha                 =0   ! Count for extended alpha fields
  !   INTEGER :: LastExtendNum                   =0   ! Count for extended numeric fields


  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)    :: ObjectNum        ! Number of the object definition to be extended.
  INTEGER, INTENT(INOUT) :: NumNewArgsLimit  ! Number of the parameters after extension


  ! SUBROUTINE PARAMETER DEFINITIONS:
  INTEGER, PARAMETER :: NewAlloc=1000  ! number of new items to allocate (* number of fields)

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumAlphaField
  INTEGER :: NumNumericField
  INTEGER :: NumNewAlphas
  INTEGER :: NumNewNumerics
  INTEGER :: NumNewParams
  INTEGER :: NumExtendFields
  INTEGER :: NumParams
  INTEGER :: Loop
  INTEGER :: Count
  INTEGER :: Item
  !  LOGICAL :: MaxArgsChanged
  LOGICAL, DIMENSION(:), ALLOCATABLE :: AorN
  LOGICAL, DIMENSION(:), ALLOCATABLE :: TempLogical
  REAL(r64), DIMENSION(:), ALLOCATABLE :: TempReals
  CHARACTER(len=MaxFieldNameLength), DIMENSION(:), ALLOCATABLE :: TempFieldCharacter
  CHARACTER(len=MaxNameLength), DIMENSION(:), ALLOCATABLE :: TempCharacter
  CHARACTER(len=32) :: charout
  TYPE(RangeCheckDef), ALLOCATABLE, DIMENSION(:) :: TempChecks
  CHARACTER(len=MaxNameLength), SAVE :: CurObject

  write(EchoInputFile,'(A)') 'Attempting to auto-extend object='//TRIM(ObjectDef(ObjectNum)%Name)
  IF (CurObject /= ObjectDef(ObjectNum)%Name) THEN
    CALL DisplayString('Auto-extending object="'//trim(ObjectDef(ObjectNum)%Name)//'", input processing may be slow.')
    CurObject=ObjectDef(ObjectNum)%Name
  ENDIF

  NumAlphaField=0
  NumNumericField=0
  NumParams=ObjectDef(ObjectNum)%NumParams
  Count=NumParams-ObjectDef(ObjectNum)%ExtensibleNum+1
  !  MaxArgsChanged=.false.

  ALLOCATE(AorN(ObjectDef(ObjectNum)%ExtensibleNum))
  AorN=.false.
  do Loop=NumParams,Count,-1
    if (ObjectDef(ObjectNum)%AlphaOrNumeric(Loop)) then
      NumAlphaField=NumAlphaField+1
    else
      NumNumericField=NumNumericField+1
    endif
  enddo
  Item=0
  do Loop=Count,NumParams
    Item=Item+1
    AorN(Item)=ObjectDef(ObjectNum)%AlphaOrNumeric(Loop)
  enddo
  NumNewAlphas=NumAlphaField*NewAlloc
  NumNewNumerics=NumNumericField*NewAlloc
  NumNewParams=NumParams+NumNewAlphas+NumNewNumerics
  NumExtendFields=NumAlphaField+NumNumericField
  ALLOCATE(TempLogical(NumNewParams))
  TempLogical(1:NumParams)=ObjectDef(ObjectNum)%AlphaOrNumeric
  TempLogical(NumParams+1:NumNewParams)=.false.
  DEALLOCATE(ObjectDef(ObjectNum)%AlphaOrNumeric)
  ALLOCATE(ObjectDef(ObjectNum)%AlphaOrNumeric(NumNewParams))
  ObjectDef(ObjectNum)%AlphaOrNumeric=TempLogical
  DEALLOCATE(TempLogical)
  do Loop=NumParams+1,NumNewParams,NumExtendFields
    ObjectDef(ObjectNum)%AlphaOrNumeric(Loop:Loop+NumExtendFields-1)=AorN
  enddo
  DEALLOCATE(AorN)  ! done with this object AorN array.

  ! required fields -- can't be extended and required.
  ALLOCATE(TempLogical(NumNewParams))
  TempLogical(1:NumParams)=ObjectDef(ObjectNum)%ReqField
  TempLogical(NumParams+1:NumNewParams)=.false.
  DEALLOCATE(ObjectDef(ObjectNum)%ReqField)
  ALLOCATE(ObjectDef(ObjectNum)%ReqField(NumNewParams))
  ObjectDef(ObjectNum)%ReqField=TempLogical
  DEALLOCATE(TempLogical)

  ALLOCATE(TempLogical(NumNewParams))
  TempLogical(1:NumParams)=ObjectDef(ObjectNum)%AlphRetainCase
  TempLogical(NumParams+1:NumNewParams)=.false.
  DEALLOCATE(ObjectDef(ObjectNum)%AlphRetainCase)
  ALLOCATE(ObjectDef(ObjectNum)%AlphRetainCase(NumNewParams))
  ObjectDef(ObjectNum)%AlphRetainCase=TempLogical
  DEALLOCATE(TempLogical)


  if (NumAlphaField > 0) then
    ALLOCATE(TempFieldCharacter(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
    TempFieldCharacter(1:ObjectDef(ObjectNum)%NumAlpha)=ObjectDef(ObjectNum)%AlphFieldChks
    TempFieldCharacter(ObjectDef(ObjectNum)%NumAlpha+1:ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas)=Blank
    DEALLOCATE(ObjectDef(ObjectNum)%AlphFieldChks)
    ALLOCATE(ObjectDef(ObjectNum)%AlphFieldChks(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
    ObjectDef(ObjectNum)%AlphFieldChks=TempFieldCharacter
    DEALLOCATE(TempFieldCharacter)
    do Loop=ObjectDef(ObjectNum)%NumAlpha+1,ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas
      ObjectDef(ObjectNum)%LastExtendAlpha=ObjectDef(ObjectNum)%LastExtendAlpha+1
      charout=IPTrimSigDigits(ObjectDef(ObjectNum)%LastExtendAlpha)
      ObjectDef(ObjectNum)%AlphFieldChks(Loop)='Extended Alpha Field '//TRIM(charout)
    enddo

    ALLOCATE(TempCharacter(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
    TempCharacter(1:ObjectDef(ObjectNum)%NumAlpha)=ObjectDef(ObjectNum)%AlphFieldDefs
    TempCharacter(ObjectDef(ObjectNum)%NumAlpha+1:ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas)=Blank
    DEALLOCATE(ObjectDef(ObjectNum)%AlphFieldDefs)
    ALLOCATE(ObjectDef(ObjectNum)%AlphFieldDefs(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
    ObjectDef(ObjectNum)%AlphFieldDefs=TempCharacter
    DEALLOCATE(TempCharacter)

    if (ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas > MaxAlphaArgsFound) then
      ! must redimension LineItem args
      ALLOCATE(TempCharacter(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
      TempCharacter(1:ObjectDef(ObjectNum)%NumAlpha)=LineItem%Alphas
      TempCharacter(ObjectDef(ObjectNum)%NumAlpha+1:ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas)=Blank
      DEALLOCATE(LineItem%Alphas)
      ALLOCATE(LineItem%Alphas(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
      LineItem%Alphas=TempCharacter
      DEALLOCATE(TempCharacter)

      ALLOCATE(TempLogical(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
      TempLogical(1:ObjectDef(ObjectNum)%NumAlpha)=LineItem%AlphBlank
      TempLogical(ObjectDef(ObjectNum)%NumAlpha+1:ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas)=.true.
      DEALLOCATE(LineItem%AlphBlank)
      ALLOCATE(LineItem%AlphBlank(ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas))
      LineItem%AlphBlank=TempLogical
      DEALLOCATE(TempLogical)

      MaxAlphaArgsFound=ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas
      !      MaxArgsChanged=.true.
    endif

  endif

  if (NumNumericField > 0) then
    ALLOCATE(TempChecks(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
    TempChecks(1:ObjectDef(ObjectNum)%NumNumeric)=ObjectDef(ObjectNum)%NumRangeChks
    DEALLOCATE(ObjectDef(ObjectNum)%NumRangeChks)
    ALLOCATE(ObjectDef(ObjectNum)%NumRangeChks(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
    ObjectDef(ObjectNum)%NumRangeChks=TempChecks
    DEALLOCATE(TempChecks)
    do Loop=ObjectDef(ObjectNum)%NumNumeric+1,ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics
      ObjectDef(ObjectNum)%NumRangeChks(Loop)%FieldNumber=Loop
      ObjectDef(ObjectNum)%LastExtendNum=ObjectDef(ObjectNum)%LastExtendNum+1
      charout=IPTrimSigDigits(ObjectDef(ObjectNum)%LastExtendNum)
      ObjectDef(ObjectNum)%NumRangeChks(Loop)%FieldName='Extended Numeric Field '//TRIM(charout)
    enddo

    if (ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics > MaxNumericArgsFound) then
      ! must redimension LineItem args
      ALLOCATE(TempReals(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
      TempReals(1:ObjectDef(ObjectNum)%NumNumeric)=LineItem%Numbers
      TempReals(ObjectDef(ObjectNum)%NumNumeric+1:ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics)=0.0
      DEALLOCATE(LineItem%Numbers)
      ALLOCATE(LineItem%Numbers(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
      LineItem%Numbers=TempReals
      DEALLOCATE(TempReals)

      ALLOCATE(TempLogical(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
      TempLogical(1:ObjectDef(ObjectNum)%NumNumeric)=LineItem%NumBlank
      TempLogical(ObjectDef(ObjectNum)%NumNumeric+1:ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics)=.true.
      DEALLOCATE(LineItem%NumBlank)
      ALLOCATE(LineItem%NumBlank(ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics))
      LineItem%NumBlank=TempLogical
      DEALLOCATE(TempLogical)

      MaxNumericArgsFound=ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics
      !      MaxArgsChanged=.true.
    endif

  endif

  ObjectDef(ObjectNum)%NumParams=NumNewParams
  NumNewArgsLimit=NumNewParams
  ObjectDef(ObjectNum)%NumAlpha=ObjectDef(ObjectNum)%NumAlpha+NumNewAlphas
  ObjectDef(ObjectNum)%NumNumeric=ObjectDef(ObjectNum)%NumNumeric+NumNewNumerics


  RETURN

END SUBROUTINE ExtendObjectDefinition

FUNCTION ProcessNumber(String,ErrorFlag) RESULT(rProcessNumber)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function processes a string that should be numeric and
  ! returns the real value of the string.

  ! METHODOLOGY EMPLOYED:
  ! FUNCTION ProcessNumber translates the argument (a string)
  ! into a real number.  The string should consist of all
  ! numeric characters (except a decimal point).  Numerics
  ! with exponentiation (i.e. 1.2345E+03) are allowed but if
  ! it is not a valid number an error message along with the
  ! string causing the error is printed out and 0.0 is returned
  ! as the value.

  ! The Fortran input processor is used to make the conversion.

  ! REFERENCES:
  ! List directed Fortran input/output.

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: String
  LOGICAL, INTENT(OUT)         :: ErrorFlag
  REAL(r64) :: rProcessNumber

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER  :: ValidNumerics='0123456789.+-EeDd'//CHAR(9)

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  REAL(r64) Temp
  INTEGER IoStatus
  INTEGER VerNumber
  INTEGER StringLen
  CHARACTER(len=MaxNameLength) :: PString


  rProcessNumber=0.0
  !  Make sure the string has all what we think numerics should have
  PString=ADJUSTL(String)
  StringLen=LEN_TRIM(PString)
  ErrorFlag=.false.
  IoStatus=0
  IF (StringLen == 0) RETURN
  VerNumber=VERIFY(PString(1:StringLen),ValidNumerics)
  IF (VerNumber == 0) THEN
    Read(PString,*,IOSTAT=IoStatus) Temp
    rProcessNumber=Temp
    ErrorFlag=.false.
  ELSE
    rProcessNumber=0.0
    ErrorFlag=.true.
  ENDIF
  IF (IoStatus /= 0) THEN
    rProcessNumber=0.0
    ErrorFlag=.true.
  ENDIF

  RETURN

END FUNCTION ProcessNumber

SUBROUTINE ProcessMinMaxDefLine(UCInputLine,WhichMinMax,MinMaxString,Value,DefaultString,ErrLevel)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine processes the IDD lines that start with
  ! \minimum or \maximum and set up the parameters so that it can
  ! be automatically checked.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! IDD Statements.
  !  \minimum         Minimum that includes the following value
  !  i.e. min >=
  !  \minimum>        Minimum that must be > than the following value
  !
  !  \maximum         Maximum that includes the following value
  !  i.e. max <=
  !  \maximum<        Maximum that must be < than the following value
  !
  !  \default         Default for field (when field is blank)

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: UCInputLine ! part of input line starting \min or \max
  INTEGER, INTENT(OUT)          :: WhichMinMax  !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
  CHARACTER(len=*), INTENT(OUT) :: MinMaxString
  REAL(r64), INTENT(OUT)             :: Value
  CHARACTER(len=*), INTENT(OUT) :: DefaultString
  INTEGER, INTENT(OUT)          :: ErrLevel

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Pos
  INTEGER NSpace
  INTEGER, EXTERNAL :: FindNonSpace
  LOGICAL ErrFlag

  ErrLevel=0
  Pos=SCAN(UCInputLine,' ')

  SELECT CASE (MakeUPPERCase(UCInputLine(1:4)))

  CASE('\MIN')
    WhichMinMax=1
    IF (SCAN(UCInputLine,'>') /= 0) THEN
      Pos=SCAN(UCInputLine,'>')+1
      WhichMinMax=2
    ENDIF
    IF (WhichMinMax == 1) THEN
      MinMaxString='>='
    ELSE
      MinMaxString='>'
    ENDIF

  CASE('\MAX')
    WhichMinMax=3
    IF (SCAN(UCInputLine,'<') /= 0) THEN
      POS=SCAN(UCInputLine,'<')+1
      WhichMinMax=4
    ENDIF
    IF (WhichMinMax == 3) THEN
      MinMaxString='<='
    ELSE
      MinMaxString='<'
    ENDIF

  CASE('\DEF')
    WhichMinMax=5
    MinMaxString=Blank

  CASE('\AUT')
    WhichMinMax=6
    MinMaxString=Blank

  CASE DEFAULT
    WhichMinMax=0  ! invalid field
    MinMaxString=Blank
    Value=-999999.d0

  END SELECT

  IF (WhichMinMax /= 0) THEN
    NSpace=FindNonSpace(UCInputLine(Pos:))
    IF (NSpace == 0) THEN
      IF (WhichMinMax /= 6) THEN  ! Only autosize/autocalculate can't have argument
        CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
        'Min/Max/Default field cannot be blank -- must have value',EchoInputFile)
        ErrLevel=2
      ELSEIF (UCINPUTLINE(1:6) == '\AUTOS') THEN
        Value=DefAutosizeValue
      ELSEIF (UCINPUTLINE(1:6) == '\AUTOC') THEN
        Value=DefAutocalculateValue
      ENDIF
    ELSE
      Pos=Pos+NSpace-1
      NSpace=SCAN(UCInputLine(Pos:),' !')
      MinMaxString=TRIM(MinMaxString)//TRIM(UCInputLine(Pos:Pos+NSpace-1))
      Value=ProcessNumber(UCInputLine(Pos:Pos+NSpace-1),ErrFlag)
      IF (ErrFlag) ErrLevel=1
      NSpace=Scan(UCInputLine(Pos:),'!')
      IF (NSpace > 0) THEN
        DefaultString=UCInputLine(Pos:Pos+NSpace-2)
      ELSE
        DefaultString=UCInputLine(Pos:)
      ENDIF
      DefaultString=ADJUSTL(DefaultString)
      IF (DefaultString == Blank) THEN
        IF (WhichMinMax == 6) THEN
          IF (UCINPUTLINE(1:6) == '\AUTOS') THEN
            Value=DefAutosizeValue
          ELSE
            Value=DefAutoCalculateValue
          ENDIF
        ELSE
          CALL ShowSevereError('IP: IDD Line='//TRIM(IPTrimSigDigits(NumLines))//  &
          'Min/Max/Default field cannot be blank -- must have value',EchoInputFile)
          ErrLevel=2
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  RETURN

END SUBROUTINE ProcessMinMaxDefLine

INTEGER FUNCTION FindIteminList(String,ListofItems,NumItems)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function looks up a string in a similar list of
  ! items and returns the index of the item in the list, if
  ! found.  This routine is not case insensitive and doesn't need
  ! for most inputs -- they are automatically turned to UPPERCASE.
  ! If you need case insensitivity use FindItem.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: String
  CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: ListofItems
  INTEGER, INTENT(IN) :: NumItems

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count

  FindIteminList=0

  DO Count=1,NumItems
    IF (String == ListofItems(Count)) THEN
      FindIteminList=Count
      EXIT
    ENDIF
  END DO

  RETURN

END FUNCTION FindIteminList

INTEGER FUNCTION FindIteminSortedList(String,ListofItems,NumItems)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function looks up a string in a similar list of
  ! items and returns the index of the item in the list, if
  ! found.  This routine is not case insensitive and doesn't need
  ! for most inputs -- they are automatically turned to UPPERCASE.
  ! If you need case insensitivity use FindItem.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: String
  CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: ListofItems
  INTEGER, INTENT(IN) :: NumItems

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: LBnd
  INTEGER :: UBnd
  INTEGER :: Probe
  LOGICAL :: Found

  LBnd=0
  UBnd=NumItems+1
  Found=.false.

  DO WHILE (.not. found .or. Probe /= 0)
    Probe=(UBnd-LBnd)/2
    IF (Probe == 0) EXIT
    Probe=LBnd+Probe
    IF (SameString(String,ListOfItems(Probe))) THEN
      Found=.true.
      EXIT
    ELSEIF (MakeUPPERCase(String) < MakeUPPERCase(ListOfItems(Probe))) THEN
      UBnd=Probe
    ELSE
      LBnd=Probe
    ENDIF
  ENDDO

  FindIteminSortedList=Probe

  RETURN

END FUNCTION FindIteminSortedList

INTEGER FUNCTION FindItem(String,ListofItems,NumItems)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   April 1999
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function looks up a string in a similar list of
  ! items and returns the index of the item in the list, if
  ! found.  This routine is case insensitive -- it uses the
  ! SameString function to assure that both strings are in
  ! all upper case.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: String
  CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: ListofItems
  INTEGER, INTENT(IN) :: NumItems

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count
  CHARACTER(len=MaxInputLineLength) StringUC
  CHARACTER(len=MaxInputLineLength) ListUC

  FindItem=0
  FindItem=FindItemInList(String,ListofItems,NumItems)
  IF (FindItem /= 0) RETURN
  !
  StringUC=MakeUPPERCase(String)

  DO Count=1,NumItems
    ListUC=MakeUPPERCase(ListofItems(Count))
    IF (StringUC == ListUC) THEN
      FindItem=Count
      EXIT
    ENDIF
  END DO

  RETURN

END FUNCTION FindItem

FUNCTION MakeUPPERCase(InputString) RESULT (ResultString)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This function returns the Upper Case representation of the InputString.

  ! METHODOLOGY EMPLOYED:
  ! Uses the Intrinsic SCAN function to scan the lowercase representation of
  ! characters (DataStringGlobals) for each character in the given string.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine


  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: InputString    ! Input String
  CHARACTER(len=len(InputString)) ResultString ! Result String, string is limited to
  ! MaxInputLineLength because of PowerStation Compiler
  ! otherwise could say (CHARACTER(len=LEN(InputString))


  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER i                  ! Loop Counter
  INTEGER :: CurCharVal

  ResultString=InputString

  do i = 1, LEN_TRIM(InputString)
    curCharVal = ICHAR(InputString(i:i))
    SELECT CASE (curCharVal)
    CASE (97:122,224:255) !lowercase ASCII and accented characters
      ResultString(i:i) = CHAR(curCharVal-32)
    CASE DEFAULT
    END SELECT
  end do

  RETURN

END FUNCTION MakeUPPERCase

LOGICAL FUNCTION SameString(TestString1,TestString2)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   November 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This function returns true if the two strings are equal (case insensitively)

  ! METHODOLOGY EMPLOYED:
  ! Make both strings uppercase.  Do internal compare.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine


  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: TestString1  ! First String to Test
  CHARACTER(len=*), INTENT(IN) :: TestString2  ! Second String to Test


  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  ! na

  IF (LEN_TRIM(TestString1) /= LEN_TRIM(TestString2)) THEN
    SameString=.false.
  ELSEIF (TestString1 == TestString2) THEN
    SameString=.true.
  ELSE
    SameString=MakeUPPERCase(TestString1) == MakeUPPERCase(TestString2)
  ENDIF

  RETURN

END FUNCTION SameString

SUBROUTINE VerifyName(NameToVerify,NamesList,NumOfNames,ErrorFound,IsBlank,StringToDisplay)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   February 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine verifys that a new name can be added to the
  ! list of names for this item (i.e., that there isn't one of that
  ! name already and that this name is not blank).

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)               :: NameToVerify
  CHARACTER(len=*), DIMENSION(:), INTENT(IN) :: NamesList
  INTEGER, INTENT(IN)                        :: NumOfNames
  LOGICAL, INTENT(OUT)                       :: ErrorFound
  LOGICAL, INTENT(OUT)                       :: IsBlank
  CHARACTER(len=*), INTENT(IN)               :: StringToDisplay

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Found

  ErrorFound=.false.
  IF (NumOfNames > 0) THEN
    Found=FindItem(NameToVerify,NamesList,NumOfNames)
    IF (Found /= 0) THEN
      !CALL ShowSevereError(TRIM(StringToDisplay)//', duplicate name='//TRIM(NameToVerify))  !RS: Secret Search String
      IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
        WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
      END IF
      WRITE(DebugFile,*) TRIM(StringToDisplay)//', duplicate name='//TRIM(NameToVerify)
      ErrorFound=.true.
    ENDIF
  ENDIF

  IF (NameToVerify == '     ') THEN
    CALL ShowSevereError(TRIM(StringToDisplay)//', cannot be blank')
    ErrorFound=.true.
    IsBlank=.true.
  ELSE
    IsBlank=.false.
  ENDIF

  RETURN

END SUBROUTINE VerifyName

SUBROUTINE RangeCheck(ErrorsFound,WhatFieldString,WhatObjectString,ErrorLevel,  &
  LowerBoundString,LowerBoundCondition,UpperBoundString,UpperBoundCondition,  &
  ValueString)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine is a general purpose "range check" routine for GetInput routines.
  ! Using the standard "ErrorsFound" logical, this routine can produce a reasonable
  ! error message to describe the situation in addition to setting the ErrorsFound variable
  ! to true.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT)                 :: ErrorsFound          ! Set to true if error detected
  CHARACTER(len=*), INTENT(IN)           :: WhatFieldString      ! Descriptive field for string
  CHARACTER(len=*), INTENT(IN)           :: WhatObjectString     ! Descriptive field for object, Zone Name, etc.
  CHARACTER(len=*), INTENT(IN)           :: ErrorLevel           ! 'Warning','Severe','Fatal')
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: LowerBoundString     ! String for error message, if applicable
  LOGICAL, INTENT(IN), OPTIONAL          :: LowerBoundCondition  ! Condition for error condition, if applicable
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: UpperBoundString     ! String for error message, if applicable
  LOGICAL, INTENT(IN), OPTIONAL          :: UpperBoundCondition  ! Condition for error condition, if applicable
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: ValueString          ! Value with digits if to be displayed with error

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=7) ErrorString  ! Uppercase representation of ErrorLevel
  LOGICAL Error
  CHARACTER(len=300) Message

  Error=.false.
  IF (PRESENT(UpperBoundCondition)) THEN
    IF (.not. UpperBoundCondition) Error=.true.
  ENDIF
  IF (PRESENT(LowerBoundCondition)) THEN
    IF (.not. LowerBoundCondition) Error=.true.
  ENDIF

  IF (Error) THEN
    CALL ConvertCasetoUPPER(ErrorLevel,ErrorString)
    Message='Out of range value field='//TRIM(WhatFieldString)
    IF (PRESENT(ValueString)) Message=trim(Message)//', Value=['//trim(ValueString)//']'
    Message=trim(Message)//', range={'
    IF (PRESENT(LowerBoundString)) Message=TRIM(Message)//TRIM(LowerBoundString)
    IF (PRESENT(LowerBoundString) .and. PRESENT(UpperBoundString)) THEN
      Message=TRIM(Message)//' and '//TRIM(UpperBoundString)
    ELSEIF (PRESENT(UpperBoundString)) THEN
      Message=TRIM(Message)//TRIM(UpperBoundString)
    ENDIF
    Message=TRIM(Message)//'}, for item='//TRIM(WhatObjectString)

    SELECT CASE(ErrorString(1:1))

    CASE('W','w')
      CALL ShowWarningError(TRIM(Message))

    CASE('S','s')
      CALL ShowSevereError(TRIM(Message))
      ErrorsFound=.true.

    CASE('F','f')
      CALL ShowFatalError(TRIM(Message))

    CASE DEFAULT
      CALL ShowSevereError(TRIM(Message))
      ErrorsFound=.true.

    END SELECT

  ENDIF

  RETURN

END SUBROUTINE RangeCheck

SUBROUTINE InternalRangeCheck(Value,FieldNumber,WhichObject,PossibleAlpha,AutoSizable,AutoCalculatable)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine is an internal range check that checks fields which have
  ! the \min and/or \max values set for appropriate values.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  REAL(r64), INTENT(IN)             :: Value
  INTEGER, INTENT(IN)          :: FieldNumber
  INTEGER, INTENT(IN)          :: WhichObject
  CHARACTER(len=*), INTENT(IN) :: PossibleAlpha
  LOGICAL, INTENT(IN)          :: AutoSizable
  LOGICAL, INTENT(IN)          :: AutoCalculatable

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL Error
  CHARACTER(len=32) FieldString
  CHARACTER(len=MaxFieldNameLength) FieldNameString
  CHARACTER(len=25) ValueString
  CHARACTER(len=300) Message

  Error=.false.
  IF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(1) == 1) THEN
    IF (Value < ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxValue(1)) Error=.true.
  ELSEIF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(1) == 2) THEN
    IF (Value <= ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxValue(1)) Error=.true.
  ENDIF
  IF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(2) == 3) THEN
    IF (Value > ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxValue(2)) Error=.true.
  ELSEIF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(2) == 4) THEN
    IF (Value >= ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxValue(2)) Error=.true.
  ENDIF

  IF (Error) THEN
    IF (.not. (AutoSizable .and. Value == ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%AutoSizeValue) .and.   &
    .not. (AutoCalculatable .and. Value == ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%AutoCalculateValue)) THEN
    NumOutOfRangeErrorsFound=NumOutOfRangeErrorsFound+1
    IF (ReportRangeCheckErrors) THEN
      FieldString=IPTrimSigDigits(FieldNumber)
      FieldNameString=ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%FieldName
      WRITE(ValueString,'(F20.5)') Value
      ValueString=ADJUSTL(ValueString)
      IF (FieldNameString /= Blank) THEN
        Message='Out of range value Numeric Field#'//TRIM(FieldString)//' ('//TRIM(FieldNameString)//  &
        '), value='//TRIM(ValueString)//', range={'
      ELSE ! Field Name not recorded
        Message='Out of range value Numeric Field#'//TRIM(FieldString)//', value='//TRIM(ValueString)//', range={'
      ENDIF
      IF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(1) /= 0) &
      Message=TRIM(Message)//ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxString(1)
      IF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(1) /= 0 .and. &
      ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(2) /= 0) THEN
      Message=TRIM(Message)//' and '//ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxString(2)
    ELSEIF (ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%WhichMinMax(2) /= 0) THEN
      Message=TRIM(Message)//ObjectDef(WhichObject)%NumRangeChks(FieldNumber)%MinMaxString(2)
    ENDIF
    Message=TRIM(Message)//'}, in '//TRIM(ObjectDef(WhichObject)%Name)
    IF (ObjectDef(WhichObject)%NameAlpha1) THEN
      Message=TRIM(Message)//'='//PossibleAlpha
    ENDIF
    CALL ShowSevereError(TRIM(Message),EchoInputFile)
  ENDIF
ENDIF
ENDIF


RETURN

END SUBROUTINE InternalRangeCheck

INTEGER FUNCTION GetNumRangeCheckErrorsFound()

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   July 2000
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function returns the number of OutOfRange errors found during
  ! input processing.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  ! na

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  ! na

  GetNumRangeCheckErrorsFound=NumOutOfRangeErrorsFound

  RETURN

END FUNCTION GetNumRangeCheckErrorsFound

!==============================================================================
! The following routines allow access to the definition lines of the IDD and
! thus can be used to "report" on expected arguments for the Input Processor.

SUBROUTINE GetObjectDefMaxArgs(ObjectWord,NumArgs,NumAlpha,NumNumeric)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   October 2001
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine returns maximum argument limits (total, alphas, numerics) of an Object from the IDD.
  ! These dimensions (not sure what one can use the total for) can be used to dynamically dimension the
  ! arrays in the GetInput routines.

  ! METHODOLOGY EMPLOYED:
  ! Essentially allows outside access to internal variables of the InputProcessor.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: ObjectWord ! Object for definition
  INTEGER, INTENT(OUT) :: NumArgs                              ! How many arguments (max) this Object can have
  INTEGER, INTENT(OUT) :: NumAlpha                             ! How many Alpha arguments (max) this Object can have
  INTEGER, INTENT(OUT) :: NumNumeric                           ! How many Numeric arguments (max) this Object can have

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Which  ! to determine which object definition to use

  IF (SortedIDD) THEN
    Which=FindItemInSortedList(MakeUPPERCase(ObjectWord),ListOfObjects,NumObjectDefs)
    IF (Which/= 0) Which=iListofObjects(Which)
  ELSE
    Which=FindItemInList(MakeUPPERCase(ObjectWord),ListOfObjects,NumObjectDefs)
  ENDIF

  IF (Which > 0) THEN
    NumArgs=ObjectDef(Which)%NumParams
    NumAlpha=ObjectDef(Which)%NumAlpha
    NumNumeric=ObjectDef(Which)%NumNumeric
  ELSE
    NumArgs=0
    NumAlpha=0
    NumNumeric=0
    CALL ShowSevereError('GetObjectDefMaxArgs: Did not find object="'//TRIM(ObjectWord)//'" in list of objects.')
  END IF

  RETURN

END SUBROUTINE GetObjectDefMaxArgs

SUBROUTINE GetIDFRecordsStats(iNumberOfRecords,iNumberOfDefaultedFields,iTotalFieldsWithDefaults,  &
  iNumberOfAutosizedFields,iTotalAutoSizableFields,                    &
  iNumberOfAutoCalcedFields,iTotalAutoCalculatableFields)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   February 2009
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine provides some statistics on the current IDF, such as number of records, total fields with defaults,
  ! number of fields that overrode the default (even if it was default value), and similarly for Autosize.

  ! METHODOLOGY EMPLOYED:
  ! Traverses the IDF Records looking at each field vs object definition for defaults and autosize.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(INOUT) :: iNumberOfRecords             ! Number of IDF Records
  INTEGER, INTENT(INOUT) :: iNumberOfDefaultedFields     ! Number of defaulted fields in IDF
  INTEGER, INTENT(INOUT) :: iTotalFieldsWithDefaults     ! Total number of fields that could be defaulted
  INTEGER, INTENT(INOUT) :: iNumberOfAutosizedFields     ! Number of autosized fields in IDF
  INTEGER, INTENT(INOUT) :: iTotalAutoSizableFields      ! Total number of autosizeable fields
  INTEGER, INTENT(INOUT) :: iNumberOfAutoCalcedFields    ! Total number of autocalculate fields
  INTEGER, INTENT(INOUT) :: iTotalAutoCalculatableFields ! Total number of autocalculatable fields

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: iRecord
  INTEGER :: iField
  INTEGER :: iObjectDef

  iNumberOfRecords=NumIDFRecords
  iNumberOfDefaultedFields     =0
  iTotalFieldsWithDefaults     =0
  iNumberOfAutosizedFields     =0
  iTotalAutoSizableFields      =0
  iNumberOfAutoCalcedFields    =0
  iTotalAutoCalculatableFields =0

  DO iRecord=1,NumIDFRecords
    IF (IDFRecords(iRecord)%ObjectDefPtr <= 0 .or. IDFRecords(iRecord)%ObjectDefPtr > NumObjectDefs) CYCLE
    iObjectDef=IDFRecords(iRecord)%ObjectDefPtr
    DO iField=1,IDFRecords(iRecord)%NumAlphas
      IF (ObjectDef(iObjectDef)%AlphFieldDefs(iField) /= Blank) iTotalFieldsWithDefaults=iTotalFieldsWithDefaults+1
      IF (ObjectDef(iObjectDef)%AlphFieldDefs(iField) /= Blank .and. IDFRecords(iRecord)%AlphBlank(iField)) &
      iNumberOfDefaultedFields=iNumberOfDefaultedFields+1
    ENDDO
    DO iField=1,IDFRecords(iRecord)%NumNumbers
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%DefaultChk) iTotalFieldsWithDefaults=iTotalFieldsWithDefaults+1
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%DefaultChk .and. IDFRecords(iRecord)%NumBlank(iField)) &
      iNumberOfDefaultedFields=iNumberOfDefaultedFields+1
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoSizable) iTotalAutoSizableFields=iTotalAutoSizableFields+1
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoSizable .and.   &
      IDFRecords(iRecord)%Numbers(iField) == ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoSizeValue) &
      iNumberOfAutosizedFields=iNumberOfAutosizedFields+1
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoCalculatable) iTotalAutoCalculatableFields=iTotalAutoCalculatableFields+1
      IF (ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoCalculatable .and.   &
      IDFRecords(iRecord)%Numbers(iField) == ObjectDef(iObjectDef)%NumRangeChks(iField)%AutoCalculateValue) &
      iNumberOfAutoCalcedFields=iNumberOfAutoCalcedFields+1
    ENDDO
  ENDDO

  RETURN

END SUBROUTINE GetIDFRecordsStats

SUBROUTINE ReportOrphanRecordObjects

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   August 2002
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine reports "orphan" objects that are in the IDF but were
  ! not "gotten" during the simulation.

  ! METHODOLOGY EMPLOYED:
  ! Uses internal (to InputProcessor) IDFRecordsGotten array, cross-matched with Object
  ! names -- puts those into array to be printed (not adding dups).

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na
  CHARACTER(len=MaxAlphaArgLength), ALLOCATABLE, DIMENSION(:) :: OrphanObjectNames
  CHARACTER(len=MaxNameLength), ALLOCATABLE, DIMENSION(:) :: OrphanNames
  INTEGER Count
  INTEGER Found
  INTEGER ObjFound
  INTEGER NumOrphObjNames

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.

  !OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging


  ALLOCATE(OrphanObjectNames(NumIDFRecords),OrphanNames(NumIDFRecords))
  OrphanObjectNames=Blank
  OrphanNames=Blank
  NumOrphObjNames=0

  DO Count=1,NumIDFRecords
    IF (IDFRecordsGotten(Count)) CYCLE
    !  This one not gotten
    Found=FindIteminList(IDFRecords(Count)%Name,OrphanObjectNames,NumOrphObjNames)
    IF (Found == 0) THEN
      IF (SortedIDD) THEN
        ObjFound=FindItemInSortedList(IDFRecords(Count)%Name,ListOfObjects,NumObjectDefs)
        IF (ObjFound /= 0) ObjFound=iListofObjects(ObjFound)
      ELSE
        ObjFound=FindItemInList(IDFRecords(Count)%Name,ListOfObjects,NumObjectDefs)
      ENDIF
      IF (ObjFound > 0) THEN
        IF (ObjectDef(ObjFound)%ObsPtr > 0) CYCLE   ! Obsolete object, don't report "orphan"
        NumOrphObjNames=NumOrphObjNames+1
        OrphanObjectNames(NumOrphObjNames)=IDFRecords(Count)%Name
        IF (ObjectDef(ObjFound)%NameAlpha1) THEN
          OrphanNames(NumOrphObjNames)=IDFRecords(Count)%Alphas(1)
        ENDIF
      ELSE
        CALL ShowWarningError('object not found='//trim(IDFRecords(Count)%Name))
      ENDIF
    ELSEIF (DisplayAllWarnings) THEN
      IF (SortedIDD) THEN
        ObjFound=FindItemInSortedList(IDFRecords(Count)%Name,ListOfObjects,NumObjectDefs)
        IF (ObjFound /= 0) ObjFound=iListofObjects(ObjFound)
      ELSE
        ObjFound=FindItemInList(IDFRecords(Count)%Name,ListOfObjects,NumObjectDefs)
      ENDIF
      IF (ObjFound > 0) THEN
        IF (ObjectDef(ObjFound)%ObsPtr > 0) CYCLE   ! Obsolete object, don't report "orphan"
        NumOrphObjNames=NumOrphObjNames+1
        OrphanObjectNames(NumOrphObjNames)=IDFRecords(Count)%Name
        IF (ObjectDef(ObjFound)%NameAlpha1) THEN
          OrphanNames(NumOrphObjNames)=IDFRecords(Count)%Alphas(1)
        ENDIF
      ELSE
        CALL ShowWarningError('ReportOrphanRecordObjects: object not found='//trim(IDFRecords(Count)%Name))
      ENDIF
    ENDIF
  ENDDO

  IF (NumOrphObjNames > 0 .and. DisplayUnusedObjects) THEN
    WRITE(EchoInputFile,*) 'Unused Objects -- Objects in IDF that were never "gotten"'
    DO Count=1,NumOrphObjNames
      IF (OrphanNames(Count) /= Blank) THEN
        WRITE(EchoInputFile,fmta) ' '//TRIM(OrphanObjectNames(Count))//'='//TRIM(OrphanNames(Count))
      ELSE
        WRITE(EchoInputFile,*) TRIM(OrphanObjectNames(Count))
      ENDIF
    ENDDO
    CALL ShowWarningError('The following lines are "Unused Objects".  These objects are in the idf')
    CALL ShowContinueError(' file but are never obtained by the simulation and therefore are NOT used.')
    IF (.not. DisplayAllWarnings) THEN
      CALL ShowContinueError(' Only the first unused named object of an object class is shown.  '//  &
      'Use Output:Diagnostics,DisplayAllWarnings to see all.')
    ELSE
      CALL ShowContinueError(' Each unused object is shown.')
    ENDIF
    CALL ShowContinueError(' See InputOutputReference document for more details.')
    IF (OrphanNames(1) /= Blank) THEN
      CALL ShowMessage('Object='//TRIM(OrphanObjectNames(1))//'='//TRIM(OrphanNames(1)))
    ELSE
      CALL ShowMessage('Object='//TRIM(OrphanObjectNames(1)))
    ENDIF
    DO Count=2,NumOrphObjNames
      IF (OrphanNames(Count) /= Blank) THEN
        CALL ShowContinueError('Object='//TRIM(OrphanObjectNames(Count))//'='//TRIM(OrphanNames(Count)))
      ELSE
        CALL ShowContinueError('Object='//TRIM(OrphanObjectNames(Count)))
      ENDIF
    ENDDO
  ELSEIF (NumOrphObjNames > 0) THEN
    !CALL ShowMessage('There are '//trim(IPTrimSigDigits(NumOrphObjNames))//' unused objects in input.')
    !CALL ShowMessage('Use Output:Diagnostics,DisplayUnusedObjects; to see them.')  !RS: Secret Search String
    WRITE(DebugFile,*) 'There are '//TRIM(IPTrimSigDigits(NumOrphObjNames))//' unused objects in input.'
    WRITE(DebugFile,*) 'Use Output:Diagnostics,DisplayUnusedObjects; to see them.'
  ENDIF

  IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
    WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
  END IF

  WRITE(DebugFile,*) 'EchoInputFile=',EchoInputFile    !RS: Debugging: Trying to find error in WeatherDataFileNumber

  DEALLOCATE(OrphanObjectNames)
  DEALLOCATE(OrphanNames)

  RETURN

END SUBROUTINE ReportOrphanRecordObjects

SUBROUTINE InitSecretObjects

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   March 2003
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine holds a set of objects that are either exact replacements for existing
  ! objects or objects which are deleted.  If these are encountered in a user input file, they
  ! will be flagged with a warning message but will not cause termination.  This routine allocates
  ! and builds an internal structure used by the InputProcessor.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  NumSecretObjects=5
  ALLOCATE(RepObjects(NumSecretObjects))

  RepObjects(1)%OldName='SKY RADIANCE DISTRIBUTION'
  RepObjects(1)%Deleted=.true.

  RepObjects(2)%OldName='SURFACE:SHADING:DETACHED'
  RepObjects(2)%NewName='Shading:Site:Detailed'

  RepObjects(3)%OldName='AIRFLOW MODEL'
  RepObjects(3)%Deleted=.true.

  RepObjects(4)%OldName='AIRFLOWNETWORK:MULTIZONE:SITEWINDCONDITIONS'
  RepObjects(4)%Deleted=.true.

  RepObjects(5)%OldName='OUTPUT:REPORTS'
  RepObjects(5)%NewName='various - depends on fields'
  RepObjects(5)%Deleted=.true.
  RepObjects(5)%TransitionDefer=.true.  ! defer transition until ready to write IDF Record

  RETURN

END SUBROUTINE InitSecretObjects

SUBROUTINE MakeTransition(ObjPtr)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   March 2009
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! For those who keep Output:Reports in their input files, this will make a
  ! transition before storing in IDF Records

  ! METHODOLOGY EMPLOYED:
  ! Manipulates LineItem structure

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(INOUT) :: ObjPtr   ! Pointer to Object Definition

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  !
  IF (MakeUPPERCase(LineItem%Name) /= 'OUTPUT:REPORTS')   &
  CALL ShowFatalError('Invalid object for deferred transition='//trim(LineItem%Name))
  IF (LineItem%NumAlphas < 1) &
  CALL ShowFatalError('Invalid object for deferred transition='//trim(LineItem%Name))

  SELECT CASE (MakeUPPERCase(LineItem%Alphas(1)))

  CASE ('VARIABLEDICTIONARY')
    LineItem%Name='OUTPUT:VARIABLEDICTIONARY'
    IF (SameString(LineItem%Alphas(2),'IDF')) THEN
      LineItem%Alphas(1)='IDF'
    ELSE
      LineItem%Alphas(1)='REGULAR'
    ENDIF
    LineItem%NumAlphas=1
    IF (SameString(LineItem%Alphas(3),'Name')) THEN
      LineItem%Alphas(2)='NAME'
      LineItem%NumAlphas=2
    ELSE
      LineItem%Alphas(2)='NONE'
      LineItem%NumAlphas=2
    ENDIF

  CASE ('SURFACES')
    ! Depends on first Alpha
    SELECT CASE(MakeUPPERCase(LineItem%Alphas(2)))

    CASE ('DXF', 'DXF:WIREFRAME', 'VRML')
      LineItem%Name='OUTPUT:SURFACES:DRAWING'
      LineItem%Alphas(1)=LineItem%Alphas(2)
      LineItem%NumAlphas=1
      IF (LineItem%Alphas(3) /= Blank) THEN
        LineItem%NumAlphas=LineItem%NumAlphas+1
        LineItem%Alphas(2)=LineItem%Alphas(3)
      ENDIF
      IF (LineItem%Alphas(4) /= Blank) THEN
        LineItem%NumAlphas=LineItem%NumAlphas+1
        LineItem%Alphas(3)=LineItem%Alphas(4)
      ENDIF

    CASE ('LINES', 'DETAILS', 'VERTICES', 'DETAILSWITHVERTICES', 'VIEWFACTORINFO', 'COSTINFO')
      LineItem%Name='OUTPUT:SURFACES:LIST'
      LineItem%Alphas(1)=LineItem%Alphas(2)
      LineItem%NumAlphas=1
      IF (LineItem%Alphas(3) /= Blank) THEN
        LineItem%NumAlphas=LineItem%NumAlphas+1
        LineItem%Alphas(2)=LineItem%Alphas(3)
      ENDIF

    CASE DEFAULT
      CALL ShowSevereError('MakeTransition: Cannot transition='//trim(LineItem%Name)//  &
      ', first field='//trim(LineItem%Alphas(1))//', second field='//trim(LineItem%Alphas(2)))

    END SELECT

  CASE ('CONSTRUCTIONS', 'CONSTRUCTION')
    LineItem%Name='OUTPUT:CONSTRUCTIONS'
    LineItem%Alphas(1)='CONSTRUCTIONS'
    LineItem%NumAlphas=1

  CASE ('MATERIALS', 'MATERIAL')
    LineItem%Name='OUTPUT:CONSTRUCTIONS'
    LineItem%Alphas(1)='MATERIALS'
    LineItem%NumAlphas=1

  CASE ('SCHEDULES')
    LineItem%Name='OUTPUT:SCHEDULES'
    LineItem%Alphas(1)=LineItem%Alphas(2)
    LineItem%NumAlphas=1

  CASE DEFAULT
    CALL ShowSevereError('MakeTransition: Cannot transition='//trim(LineItem%Name)//  &
    ', first field='//trim(LineItem%Alphas(1)))

  END SELECT

  ObjectDef(ObjPtr)%NumFound=ObjectDef(ObjPtr)%NumFound-1
  ObjPtr=FindItemInList(LineItem%Name,ListOfObjects,NumObjectDefs)
  ObjPtr=iListofObjects(ObjPtr)

  IF (ObjPtr == 0) CALL ShowFatalError('No Object Def for '//trim(LineItem%Name))
  ObjectDef(ObjPtr)%NumFound=ObjectDef(ObjPtr)%NumFound+1

  RETURN

END SUBROUTINE MakeTransition

SUBROUTINE AddRecordFromSection(Which)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   March 2009
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! When an object is entered like a section (i.e., <objectname>;), try to add a record
  ! of the object using minfields, etc.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: Which ! Which object was matched

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumArg
  INTEGER :: NumAlpha
  INTEGER :: NumNumeric
  INTEGER :: Count
  CHARACTER(len=52) :: String

  NumArg=0
  LineItem%Name=ObjectDef(Which)%Name
  LineItem%Alphas=Blank
  LineItem%AlphBlank=.false.
  LineItem%NumAlphas=0
  LineItem%Numbers=0.0
  LineItem%NumNumbers=0
  LineItem%NumBlank=.false.
  LineItem%ObjectDefPtr=Which

  ObjectDef(Which)%NumFound=ObjectDef(Which)%NumFound+1

  ! Check out MinimumNumberOfFields
  IF (NumArg < ObjectDef(Which)%MinNumFields) THEN
    IF (ObjectDef(Which)%NameAlpha1) THEN
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef(Which)%Name)//  &
      ', name='//TRIM(LineItem%Alphas(1))//       &
      ', entered with less than minimum number of fields.')
    ELSE
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef(Which)%Name)//  &
      ', entered with less than minimum number of fields.')
    ENDIF
    CALL ShowAuditErrorMessage(' **   ~~~   ** ','Attempting fill to minimum.')
    NumAlpha=0
    NumNumeric=0
    IF (ObjectDef(Which)%MinNumFields > ObjectDef(Which)%NumParams) THEN
      CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object \min-fields > number of fields specified, Object='//TRIM(ObjectDef(Which)%Name))
      CALL ShowContinueError('..\min-fields='//TRIM(IPTrimSigDigits(ObjectDef(Which)%MinNumFields))//  &
      ', total number of fields in object definition='//TRIM(IPTrimSigDigits(ObjectDef(Which)%NumParams)))
      !      ErrFlag=.true.
    ELSE
      DO Count=1,ObjectDef(Which)%MinNumFields
        IF (ObjectDef(Which)%AlphaOrNumeric(Count)) THEN
          NumAlpha=NumAlpha+1
          IF (NumAlpha <= LineItem%NumAlphas) CYCLE
          LineItem%NumAlphas=LineItem%NumAlphas+1
          IF (ObjectDef(Which)%AlphFieldDefs(LineItem%NumAlphas) /= Blank) THEN
            LineItem%Alphas(LineItem%NumAlphas)=ObjectDef(Which)%AlphFieldDefs(LineItem%NumAlphas)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(ObjectDef(Which)%AlphFieldDefs(LineItem%NumAlphas))//   &
            '   ! field=>'//TRIM(ObjectDef(Which)%AlphFieldChks(NumAlpha)))
          ELSEIF (ObjectDef(Which)%ReqField(Count)) THEN
            IF (ObjectDef(Which)%NameAlpha1) THEN
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef(Which)%Name)//  &
              ', name='//TRIM(LineItem%Alphas(1))// &
              ', Required Field=['//  &
              TRIM(ObjectDef(Which)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ELSE
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef(Which)%Name)//  &
              ', Required Field=['//  &
              TRIM(ObjectDef(Which)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ENDIF
            !            ErrFlag=.true.
          ELSE
            LineItem%Alphas(LineItem%NumAlphas)=Blank
            LineItem%AlphBlank(LineItem%NumAlphas)=.true.
            CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
            TRIM(ObjectDef(Which)%AlphFieldChks(NumAlpha)))
          ENDIF
        ELSE
          NumNumeric=NumNumeric+1
          IF (NumNumeric <= LineItem%NumNumbers) CYCLE
          LineItem%NumNumbers=LineItem%NumNumbers+1
          LineItem%NumBlank(NumNumeric)=.true.
          IF (ObjectDef(Which)%NumRangeChks(NumNumeric)%Defaultchk) THEN
            IF (.not. ObjectDef(Which)%NumRangeChks(NumNumeric)%DefAutoSize .and.   &
            .not. ObjectDef(Which)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Which)%NumRangeChks(NumNumeric)%Default
            WRITE(String,*) ObjectDef(Which)%NumRangeChks(NumNumeric)%Default
            String=ADJUSTL(String)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(String)//  &
            '   ! field=>'//TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef(Which)%NumRangeChks(NumNumeric)%DefAutoSize) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Which)%NumRangeChks(NumNumeric)%AutoSizeValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autosize    ! field=>'//  &
            TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef(Which)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem%Numbers(NumNumeric)=ObjectDef(Which)%NumRangeChks(NumNumeric)%AutoCalculateValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autocalculate    ! field=>'//  &
            TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName))
          ENDIF
        ELSEIF (ObjectDef(Which)%ReqField(Count)) THEN
          IF (ObjectDef(Which)%NameAlpha1) THEN
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef(Which)%Name)//  &
            ', name='//TRIM(LineItem%Alphas(1))// &
            ', Required Field=['//  &
            TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ELSE
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef(Which)%Name)//  &
            ', Required Field=['//  &
            TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ENDIF
          !            ErrFlag=.true.
        ELSE
          LineItem%Numbers(NumNumeric)=0.0
          LineItem%NumBlank(NumNumeric)=.true.
          CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
          TRIM(ObjectDef(Which)%NumRangeChks(NumNumeric)%FieldName))
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF

!  IF (TransitionDefer) THEN
!    CALL MakeTransition(Which)
!  ENDIF
NumIDFRecords=NumIDFRecords+1
IF (ObjectStartRecord(Which) == 0) ObjectStartRecord(Which)=NumIDFRecords
MaxAlphaIDFArgsFound=MAX(MaxAlphaIDFArgsFound,LineItem%NumAlphas)
MaxNumericIDFArgsFound=MAX(MaxNumericIDFArgsFound,LineItem%NumNumbers)
MaxAlphaIDFDefArgsFound=MAX(MaxAlphaIDFDefArgsFound,ObjectDef(Which)%NumAlpha)
MaxNumericIDFDefArgsFound=MAX(MaxNumericIDFDefArgsFound,ObjectDef(Which)%NumNumeric)
IDFRecords(NumIDFRecords)%Name=LineItem%Name
IDFRecords(NumIDFRecords)%NumNumbers=LineItem%NumNumbers
IDFRecords(NumIDFRecords)%NumAlphas=LineItem%NumAlphas
IDFRecords(NumIDFRecords)%ObjectDefPtr=LineItem%ObjectDefPtr
ALLOCATE(IDFRecords(NumIDFRecords)%Alphas(LineItem%NumAlphas))
ALLOCATE(IDFRecords(NumIDFRecords)%AlphBlank(LineItem%NumAlphas))
ALLOCATE(IDFRecords(NumIDFRecords)%Numbers(LineItem%NumNumbers))
ALLOCATE(IDFRecords(NumIDFRecords)%NumBlank(LineItem%NumNumbers))
IDFRecords(NumIDFRecords)%Alphas(1:LineItem%NumAlphas)=LineItem%Alphas(1:LineItem%NumAlphas)
IDFRecords(NumIDFRecords)%AlphBlank(1:LineItem%NumAlphas)=LineItem%AlphBlank(1:LineItem%NumAlphas)
IDFRecords(NumIDFRecords)%Numbers(1:LineItem%NumNumbers)=LineItem%Numbers(1:LineItem%NumNumbers)
IDFRecords(NumIDFRecords)%NumBlank(1:LineItem%NumNumbers)=LineItem%NumBlank(1:LineItem%NumNumbers)
IF (LineItem%NumNumbers > 0) THEN
  DO Count=1,LineItem%NumNumbers
    IF (ObjectDef(Which)%NumRangeChks(Count)%MinMaxChk .and. .not. LineItem%NumBlank(Count)) THEN
      CALL InternalRangeCheck(LineItem%Numbers(Count),Count,Which,LineItem%Alphas(1),  &
      ObjectDef(Which)%NumRangeChks(Count)%AutoSizable,        &
      ObjectDef(Which)%NumRangeChks(Count)%AutoCalculatable)
    ENDIF
  ENDDO
ENDIF

RETURN

END SUBROUTINE AddRecordFromSection

SUBROUTINE PreProcessorCheck(PreP_Fatal)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   August 2005
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine checks for existance of "Preprocessor Message" object and
  ! performs appropriate action.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! Preprocessor Message,
  !    \memo This object does not come from a user input.  This is generated by a pre-processor
  !    \memo so that various conditions can be gracefully passed on by the InputProcessor.
  !    A1,        \field preprocessor name
  !    A2,        \field error severity
  !               \note Depending on type, InputProcessor may terminate the program.
  !               \type choice
  !               \key warning
  !               \key severe
  !               \key fatal
  !    A3,        \field message line 1
  !    A4,        \field message line 2
  !    A5,        \field message line 3
  !    A6,        \field message line 4
  !    A7,        \field message line 5
  !    A8,        \field message line 6
  !    A9,        \field message line 7
  !    A10,       \field message line 8
  !    A11,       \field message line 9
  !    A12;       \field message line 10


  ! USE STATEMENTS:
  USE DataIPShortCuts

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: PreP_Fatal  ! True if a preprocessor flags a fatal error

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: NumAlphas     ! Used to retrieve names from IDF
  INTEGER :: NumNumbers    ! Used to retrieve rNumericArgs from IDF
  INTEGER :: IOStat        ! Could be used in the Get Routines, not currently checked
  INTEGER :: NumParams     ! Total Number of Parameters in 'Output:PreprocessorMessage' Object
  INTEGER :: NumPrePM      ! Number of Preprocessor Message objects in IDF
  INTEGER :: CountP
  INTEGER :: CountM
  CHARACTER(len=1) :: Multiples
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.

  !OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  cCurrentModuleObject='Output:PreprocessorMessage'
  NumPrePM=GetNumObjectsFound(TRIM(cCurrentModuleObject))
  IF (NumPrePM > 0) THEN
    CALL GetObjectDefMaxArgs(TRIM(cCurrentModuleObject),NumParams,NumAlphas,NumNumbers)
    cAlphaArgs(1:NumAlphas)=Blank
    DO CountP=1,NumPrePM
      CALL GetObjectItem(TRIM(cCurrentModuleObject),CountP,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStat,  &
      AlphaBlank=lAlphaFieldBlanks,NumBlank=lNumericFieldBlanks,  &
      AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)
      IF (cAlphaArgs(1) == Blank) cAlphaArgs(1)='Unknown'
      IF (NumAlphas > 3) THEN
        Multiples='s'
      ELSE
        Multiples=Blank
      ENDIF
      IF (cAlphaArgs(2) == Blank) cAlphaArgs(2)='Unknown'
      SELECT CASE (MakeUPPERCase(cAlphaArgs(2)))
      CASE('INFORMATION')
        CALL ShowMessage(TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &
        '" has the following Information message'//TRIM(Multiples)//':')
      CASE('WARNING')
        CALL ShowWarningError(TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &
        '" has the following Warning condition'//TRIM(Multiples)//':')
      CASE('SEVERE')
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &
        '" has the following Severe condition'//TRIM(Multiples)//':')
      CASE('FATAL')
        !CALL ShowSevereError(TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &  !RS: Secret Search String
        !   '" has the following Fatal condition'//TRIM(Multiples)//':')
        IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
          WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
        END IF
        WRITE(DebugFile,*) TRIM(cCurrentModuleObject), ', ', TRIM(cALphaArgs(1)), ', ', TRIM(Multiples)
        PreP_Fatal=.true.
      CASE DEFAULT
        CALL ShowSevereError(TRIM(cCurrentModuleObject)//'="'//TRIM(cAlphaArgs(1))//  &
        '" has the following '//TRIM(cAlphaArgs(2))//' condition'//TRIM(Multiples)//':')
      END SELECT
      CountM=3
      IF (CountM > NumAlphas) THEN
        CALL ShowContinueError(TRIM(cCurrentModuleObject)//' was blank.  Check '//TRIM(cAlphaArgs(1))//  &
        ' audit trail or error file for possible reasons.')
      ENDIF
      DO WHILE (CountM <= NumAlphas)
        IF (LEN_TRIM(cAlphaArgs(CountM)) == MaxNameLength) THEN
          CALL ShowContinueError(TRIM(cAlphaArgs(CountM))//TRIM(cAlphaArgs(CountM+1)))
          CountM=CountM+2
        ELSE
          IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
            WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
          END IF
          !CALL ShowContinueError(TRIM(cAlphaArgs(CountM))) !RS: Secret Search String
          WRITE(DebugFile,*) TRIM(cAlphaArgs(CountM))
          CountM=CountM+1
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  RETURN

END SUBROUTINE PreProcessorCheck

SUBROUTINE CompactObjectsCheck

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   December 2005
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Check to see if Compact Objects (i.e. CompactHVAC and its ilk) exist in the
  ! input file.  If so, expandobjects was not run and there's a possible problem.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL :: CompactObjectsFound
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopefully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  CompactObjectsFound=.false.
  IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
    WRITE(*,*) 'Error with DebugFile'    !RS: Debugging: Searching for a mis-set file number
  END IF

  IF ( ANY(IDFRecords%Name(1:13) == 'HVACTEMPLATE:') .or. ANY(IDFRecords%Name(1:13) == 'HVACTemplate:') ) THEN
    !    CALL ShowSevereError('HVACTemplate objects are found in the IDF File.')    !RS: Secret Search String
    WRITE(DebugFile, *) 'HVACTemplate objects are found in the IDF File.'
    CompactObjectsFound=.true.
  ENDIF

  IF (CompactObjectsFound) THEN
    !CALL ShowFatalError('Program Terminates: The ExpandObjects program has'// &
    !  ' not been run or is not in your EnergyPlus.exe folder.')    !RS: Secret Search String
    WRITE (DebugFile, *) 'They wanted the program to terminate, but we are forcing it through anyhow'
  ENDIF

  RETURN

END SUBROUTINE CompactObjectsCheck

SUBROUTINE ParametricObjectsCheck

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Jason Glazer (based on CompactObjectsCheck by Linda Lawrie)
  !       DATE WRITTEN   September 2009
  !       MODIFIED
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! Check to see if Parametric Objects exist in the input file.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging


  IF ( ANY(IDFRecords%Name(1:11) == 'PARAMETRIC:') .or. ANY(IDFRecords%Name(1:11) == 'Parametric:') .or.   &
  ANY(IDFRecords%Name(1:11) == 'parametric:') ) THEN
  !CALL ShowSevereError('Parametric objects are found in the IDF File.')
  !CALL ShowFatalError('Program Terminates: The ParametricPreprocessor program has'// &
  !  ' not been run.')    !RS: Secret Search String
  WRITE(DebugFile,*) 'Parametric objects are found in the IDF file. The ParametricPreprocessor program has not been run.'
ENDIF

RETURN

END SUBROUTINE ParametricObjectsCheck

SUBROUTINE PreScanReportingVariables

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine scans the input records and determines which output variables
  ! are actually being requested for the run so that the OutputProcessor will only
  ! consider those variables for output.  (At this time, all metered variables are
  ! allowed to pass through).

  ! METHODOLOGY EMPLOYED:
  ! Uses internal records and structures.
  ! Looks at:
  ! Output:Variable
  ! Meter:Custom
  ! Meter:CustomDecrement
  ! Meter:CustomDifference
  ! Output:Table:Monthly
  ! Output:Table:TimeBins
  ! Output:Table:SummaryReports
  ! EnergyManagementSystem:Sensor
  ! EnergyManagementSystem:OutputVariable

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataOutputs

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: OutputVariable='OUTPUT:VARIABLE'
  CHARACTER(len=*), PARAMETER :: MeterCustom='METER:CUSTOM'
  CHARACTER(len=*), PARAMETER :: MeterCustomDecrement='METER:CUSTOMDECREMENT'
  CHARACTER(len=*), PARAMETER :: MeterCustomDifference='METER:CUSTOMDIFFERENCE'
  CHARACTER(len=*), PARAMETER :: OutputTableMonthly='OUTPUT:TABLE:MONTHLY'
  CHARACTER(len=*), PARAMETER :: OutputTableTimeBins='OUTPUT:TABLE:TIMEBINS'
  CHARACTER(len=*), PARAMETER :: OutputTableSummaries='OUTPUT:TABLE:SUMMARYREPORTS'
  CHARACTER(len=*), PARAMETER :: EMSSensor='ENERGYMANAGEMENTSYSTEM:SENSOR'
  CHARACTER(len=*), PARAMETER :: EMSOutputVariable='ENERGYMANAGEMENTSYSTEM:OUTPUTVARIABLE'

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: CurrentRecord
  INTEGER :: Loop
  INTEGER :: Loop1

  ALLOCATE(OutputVariablesForSimulation(10000))
  MaxConsideredOutputVariables=10000

  ! Output Variable
  CurrentRecord=FindFirstRecord(OutputVariable)
  DO WHILE (CurrentRecord /= 0)
    IF (IDFRecords(CurrentRecord)%NumAlphas < 2) CYCLE  ! signals error condition for later on
    IF (.not. IDFRecords(CurrentRecord)%AlphBlank(1)) THEN
      CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(1),IDFRecords(CurrentRecord)%Alphas(2))
    ELSE
      CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(2))
    ENDIF
    CurrentRecord=FindNextRecord(OutputVariable,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(MeterCustom)
  DO WHILE (CurrentRecord /= 0)
    DO Loop=3,IDFRecords(CurrentRecord)%NumAlphas,2
      IF (Loop > IDFRecords(CurrentRecord)%NumAlphas .or. Loop+1 > IDFRecords(CurrentRecord)%NumAlphas) CYCLE  ! error condition
      IF (.not. IDFRecords(CurrentRecord)%AlphBlank(Loop)) THEN
        CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(Loop),IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ELSE
        CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ENDIF
    ENDDO
    CurrentRecord=FindNextRecord(MeterCustom,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(MeterCustomDecrement)
  DO WHILE (CurrentRecord /= 0)
    DO Loop=4,IDFRecords(CurrentRecord)%NumAlphas,2
      IF (Loop > IDFRecords(CurrentRecord)%NumAlphas .or. Loop+1 > IDFRecords(CurrentRecord)%NumAlphas) CYCLE  ! error condition
      IF (.not. IDFRecords(CurrentRecord)%AlphBlank(Loop)) THEN
        CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(Loop),IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ELSE
        CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ENDIF
    ENDDO
    CurrentRecord=FindNextRecord(MeterCustomDecrement,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(MeterCustomDifference)
  DO WHILE (CurrentRecord /= 0)
    DO Loop=4,IDFRecords(CurrentRecord)%NumAlphas,2
      IF (Loop > IDFRecords(CurrentRecord)%NumAlphas .or. Loop+1 > IDFRecords(CurrentRecord)%NumAlphas) CYCLE  ! error condition
      IF (.not. IDFRecords(CurrentRecord)%AlphBlank(Loop)) THEN
        CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(Loop),IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ELSE
        CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(Loop+1))
      ENDIF
    ENDDO
    CurrentRecord=FindNextRecord(MeterCustomDifference,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(EMSSensor)
  DO WHILE (CurrentRecord /= 0)
    IF (IDFRecords(CurrentRecord)%NumAlphas < 2) CurrentRecord=FindNextRecord(EMSSensor,CurrentRecord)
    IF (IDFRecords(CurrentRecord)%Alphas(2) /= blank) THEN
      CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(2),IDFRecords(CurrentRecord)%Alphas(3))
    ELSE
      CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(3))
    ENDIF
    CurrentRecord=FindNextRecord(EMSSensor,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(EMSOutputVariable)
  DO WHILE (CurrentRecord /= 0)
    IF (IDFRecords(CurrentRecord)%NumAlphas < 2) CurrentRecord=FindNextRecord(EMSOutputVariable,CurrentRecord)
    CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(1))
    CurrentRecord=FindNextRecord(EMSOutputVariable,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(OutputTableTimeBins)
  DO WHILE (CurrentRecord /= 0)
    IF (IDFRecords(CurrentRecord)%NumAlphas < 2) CurrentRecord=FindNextRecord(OutputTableTimeBins,CurrentRecord)
    IF (.not. IDFRecords(CurrentRecord)%AlphBlank(1)) THEN
      CALL AddRecordToOutputVariableStructure(IDFRecords(CurrentRecord)%Alphas(1),IDFRecords(CurrentRecord)%Alphas(2))
    ELSE
      CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(2))
    ENDIF
    CurrentRecord=FindNextRecord(OutputTableTimeBins,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(OutputTableMonthly)
  DO WHILE (CurrentRecord /= 0)
    DO Loop=2,IDFRecords(CurrentRecord)%NumAlphas,2
      IF (IDFRecords(CurrentRecord)%NumAlphas < 2) CYCLE
      CALL AddRecordToOutputVariableStructure('*',IDFRecords(CurrentRecord)%Alphas(Loop))
    ENDDO
    CurrentRecord=FindNextRecord(OutputTableMonthly,CurrentRecord)
  ENDDO

  CurrentRecord=FindFirstRecord(OutputTableSummaries)  ! summary tables, not all add to variable structure
  DO WHILE (CurrentRecord /= 0)
    DO Loop=1,IDFRecords(CurrentRecord)%NumAlphas
      IF (IDFRecords(CurrentRecord)%Alphas(Loop) == 'ALLMONTHLY' .or.  &
      IDFRecords(CurrentRecord)%Alphas(Loop) == 'ALLSUMMARYANDMONTHLY') THEN
      DO Loop1=1,NumMonthlyReports
        CALL AddVariablesForMonthlyReport(MonthlyNamedReports(Loop1))
      ENDDO
    ELSE
      CALL AddVariablesForMonthlyReport(IDFRecords(CurrentRecord)%Alphas(Loop))
    ENDIF

  ENDDO
  CurrentRecord=FindNextRecord(OutputTableSummaries,CurrentRecord)
ENDDO

IF (NumConsideredOutputVariables > 0) THEN
  ALLOCATE(TempOutputVariablesForSimulation(NumConsideredOutputVariables))
  TempOutputVariablesForSimulation(1:NumConsideredOutputVariables)=OutputVariablesForSimulation(1:NumConsideredOutputVariables)
  DEALLOCATE(OutputVariablesForSimulation)
  ALLOCATE(OutputVariablesForSimulation(NumConsideredOutputVariables))
  OutputVariablesForSimulation(1:NumConsideredOutputVariables)=TempOutputVariablesForSimulation(1:NumConsideredOutputVariables)
  DEALLOCATE(TempOutputVariablesForSimulation)
  MaxConsideredOutputVariables=NumConsideredOutputVariables
ENDIF
RETURN

END SUBROUTINE PreScanReportingVariables

SUBROUTINE AddVariablesForMonthlyReport(ReportName)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine adds specific variables to the Output Variables for Simulation
  ! Structure. Note that only non-metered variables need to be added here.  Metered
  ! variables are automatically included in the minimized output variable structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: ReportName

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  SELECT CASE(ReportName)

  CASE ('ZONECOOLINGSUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE/SYS SENSIBLE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','ZONE TOTAL INTERNAL LATENT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE TOTAL INTERNAL LATENT GAIN RATE')

  CASE ('ZONEHEATINGSUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE/SYS SENSIBLE HEATING ENERGY')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','ZONE/SYS SENSIBLE HEATING RATE')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')

  CASE ('ZONEELECTRICSUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE LIGHTS ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','ZONE ELECTRIC EQUIPMENT ELECTRIC CONSUMPTION')

  CASE ('SPACEGAINSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE PEOPLE TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE LIGHTS TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE ELECTRIC EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE GAS EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE HOT WATER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE STEAM EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE OTHER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT LOSS')

  CASE ('PEAKSPACEGAINSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE PEOPLE TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE LIGHTS TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE ELECTRIC EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE GAS EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE HOT WATER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE STEAM EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE OTHER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT LOSS')

  CASE ('SPACEGAINCOMPONENTSATCOOLINGPEAKMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE/SYS SENSIBLE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','ZONE PEOPLE TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE LIGHTS TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE ELECTRIC EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE GAS EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE HOT WATER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE STEAM EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE OTHER EQUIPMENT TOTAL HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE INFILTRATION SENSIBLE HEAT LOSS')

  CASE ('SETPOINTSNOTMETWITHTEMPERATURESMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','TIME HEATING SETPOINT NOT MET')
    CALL AddRecordToOutputVariableStructure('*','ZONE MEAN AIR TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','TIME HEATING SETPOINT NOT MET WHILE OCCUPIED')
    CALL AddRecordToOutputVariableStructure('*','TIME COOLING SETPOINT NOT MET')
    CALL AddRecordToOutputVariableStructure('*','TIME COOLING SETPOINT NOT MET WHILE OCCUPIED')

  CASE ('COMFORTREPORTSIMPLE55MONTHLY')
    CALL AddRecordToOutputVariableStructure('*','TIME NOT COMFORTABLE SUMMER CLOTHES')
    CALL AddRecordToOutputVariableStructure('*','ZONE MEAN AIR TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','TIME NOT COMFORTABLE WINTER CLOTHES')
    CALL AddRecordToOutputVariableStructure('*','TIME NOT COMFORTABLE SUMMER OR WINTER CLOTHES')

  CASE ('UNGLAZEDTRANSPIREDSOLARCOLLECTORSUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','UTSC OVERALL EFFICIENCY')
    CALL AddRecordToOutputVariableStructure('*','UTSC AVERAGE SUCTION FACE VELOCITY')
    CALL AddRecordToOutputVariableStructure('*','UTSC SENSIBLE HEATING RATE')

  CASE ('OCCUPANTCOMFORTDATASUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','PEOPLE NUMBER OF OCCUPANTS')
    CALL AddRecordToOutputVariableStructure('*','PEOPLE AIR TEMPERATURES')
    CALL AddRecordToOutputVariableStructure('*','PEOPLE AIR RELATIVE HUMIDITY')
    CALL AddRecordToOutputVariableStructure('*','FANGERPMV')
    CALL AddRecordToOutputVariableStructure('*','FANGERPPD')

  CASE ('CHILLERREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','CHILLER ELECTRIC CONSUMPTION')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','CHILLER ELECTRIC POWER')
    CALL AddRecordToOutputVariableStructure('*','CHILLER EVAP HEAT TRANS') ! on meter
    CALL AddRecordToOutputVariableStructure('*','CHILLER COND HEAT TRANS') ! on meter
    CALL AddRecordToOutputVariableStructure('*','CHILLER COP')

  CASE ('TOWERREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','TOWER FAN ELECTRIC CONSUMPTION') ! on meter
    CALL AddRecordToOutputVariableStructure('*','TOWER FAN ELECTRIC POWER')
    CALL AddRecordToOutputVariableStructure('*','TOWER HEAT TRANSFER')
    CALL AddRecordToOutputVariableStructure('*','TOWER WATER INLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','TOWER WATER OUTLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','TOWER WATER MASS FLOW RATE')

  CASE ('BOILERREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','BOILER HEATING OUTPUT ENERGY')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','BOILER GAS CONSUMPTION')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','BOILER HEATING OUTPUT ENERGY') ! on meter
    CALL AddRecordToOutputVariableStructure('*','BOILER HEATING OUTPUT RATE')
    CALL AddRecordToOutputVariableStructure('*','BOILER GAS CONSUMPTION RATE')
    CALL AddRecordToOutputVariableStructure('*','BOILER WATER INLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','BOILER WATER OUTLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','BOILER WATER MASS FLOW RATE')
    CALL AddRecordToOutputVariableStructure('*','BOILER PARASITIC ELECTRIC CONSUMPTION RATE')

  CASE ('DXREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','DX COIL TOTAL COOLING ENERGY')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','DX COOLING COIL ELECTRIC CONSUMPTION')  ! on meter
    CALL AddRecordToOutputVariableStructure('*','DX COIL SENSIBLE COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','DX COIL LATENT COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','DX COOLING COIL CRANKCASE HEATER CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','DX COOLING COIL RUNTIME FRACTION')
    CALL AddRecordToOutputVariableStructure('*','DX COIL TOTAL COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','DX COIL SENSIBLE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','DX COIL LATENT COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','DX COOLING COIL ELECTRIC POWER')
    CALL AddRecordToOutputVariableStructure('*','DX COOLING COIL CRANKCASE HEATER POWER')

  CASE ('WINDOWREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED SOLAR')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED BEAM SOLAR')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','WINDOW HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','WINDOW HEAT LOSS')
    CALL AddRecordToOutputVariableStructure('*','INSIDE GLASS CONDENSATION FLAG')
    CALL AddRecordToOutputVariableStructure('*','FRACTION OF TIME SHADING DEVICE IS ON')
    CALL AddRecordToOutputVariableStructure('*','STORM WINDOW ON/OFF FLAG')

  CASE ('WINDOWENERGYREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED SOLAR ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED BEAM SOLAR ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW TRANSMITTED DIFFUSE SOLAR ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW HEAT GAIN ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW HEAT LOSS ENERGY')

  CASE ('WINDOWZONESUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE WINDOW HEAT GAIN')
    CALL AddRecordToOutputVariableStructure('*','ZONE WINDOW HEAT LOSS')
    CALL AddRecordToOutputVariableStructure('*','ZONE TRANSMITTED SOLAR')
    CALL AddRecordToOutputVariableStructure('*','ZONE BEAM SOLAR FROM EXTERIOR WINDOWS')
    CALL AddRecordToOutputVariableStructure('*','ZONE DIFF SOLAR FROM EXTERIOR WINDOWS')
    CALL AddRecordToOutputVariableStructure('*','ZONE DIFF SOLAR FROM INTERIOR WINDOWS')
    CALL AddRecordToOutputVariableStructure('*','ZONE BEAM SOLAR FROM INTERIOR WINDOWS')

  CASE ('WINDOWENERGYZONESUMMARYMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE WINDOW HEAT GAIN ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE WINDOW HEAT LOSS ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE TRANSMITTED SOLAR ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE BEAM SOLAR FROM EXTERIOR WINDOWS ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE DIFF SOLAR FROM EXTERIOR WINDOWS ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE DIFF SOLAR FROM INTERIOR WINDOWS ENERGY')
    CALL AddRecordToOutputVariableStructure('*','ZONE BEAM SOLAR FROM INTERIOR WINDOWS ENERGY')

  CASE ('AVERAGEOUTDOORCONDITIONSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DEW POINT')
    CALL AddRecordToOutputVariableStructure('*','WIND SPEED')
    CALL AddRecordToOutputVariableStructure('*','SKY TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','DIRECT SOLAR')
    CALL AddRecordToOutputVariableStructure('*','RAINING')

  CASE ('OUTDOORCONDITIONSMAXIMUMDRYBULBMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DEW POINT')
    CALL AddRecordToOutputVariableStructure('*','WIND SPEED')
    CALL AddRecordToOutputVariableStructure('*','SKY TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','DIRECT SOLAR')

  CASE ('OUTDOORCONDITIONSMINIMUMDRYBULBMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DEW POINT')
    CALL AddRecordToOutputVariableStructure('*','WIND SPEED')
    CALL AddRecordToOutputVariableStructure('*','SKY TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','DIRECT SOLAR')

  CASE ('OUTDOORCONDITIONSMAXIMUMWETBULBMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DEW POINT')
    CALL AddRecordToOutputVariableStructure('*','WIND SPEED')
    CALL AddRecordToOutputVariableStructure('*','SKY TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','DIRECT SOLAR')

  CASE ('OUTDOORCONDITIONSMAXIMUMDEWPOINTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DEW POINT')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR DRY BULB')
    CALL AddRecordToOutputVariableStructure('*','OUTDOOR WET BULB')
    CALL AddRecordToOutputVariableStructure('*','WIND SPEED')
    CALL AddRecordToOutputVariableStructure('*','SKY TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DIFFUSE SOLAR')
    CALL AddRecordToOutputVariableStructure('*','DIRECT SOLAR')

  CASE ('OUTDOORGROUNDCONDITIONSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','GROUND TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','SURFACE GROUND TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','DEEP GROUND TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','WATER MAINS TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','GROUND REFLECTED SOLAR')
    CALL AddRecordToOutputVariableStructure('*','SNOW ON GROUND')

  CASE ('WINDOWACREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC TOTAL ZONE COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC TOTAL ZONE COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC SENSIBLE ZONE COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC LATENT ZONE COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC TOTAL ZONE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC SENSIBLE ZONE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC LATENT ZONE COOLING RATE')
    CALL AddRecordToOutputVariableStructure('*','WINDOW AC ELECTRIC POWER')

  CASE ('WATERHEATERREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER TOTAL DEMAND ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER USE ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER BURNER HEATING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER GAS CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER TOTAL DEMAND ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER LOSS DEMAND ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER LOSS ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER TANK TEMPERATURE')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER HEAT RECOVERY SUPPLY ENERGY')
    CALL AddRecordToOutputVariableStructure('*','WATER HEATER SOURCE ENERGY')

  CASE ('GENERATORREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR ELECTRIC ENERGY PRODUCED')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR DIESEL CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR GAS CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR ELECTRIC ENERGY PRODUCED')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR TOTAL HEAT RECOVERY')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR JACKET HEAT RECOVERY')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR LUBE HEAT RECOVERY')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR EXHAUST HEAT RECOVERY')
    CALL AddRecordToOutputVariableStructure('*','GENERATOR EXHAUST STACK TEMP')

  CASE ('DAYLIGHTINGREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','EXTERIOR BEAM NORMAL ILLUMINANCE')
    CALL AddRecordToOutputVariableStructure('*','LTG POWER MULTIPLIER FROM DAYLIGHTING')
    CALL AddRecordToOutputVariableStructure('*','LTG POWER MULTIPLIER FROM DAYLIGHTING')
    CALL AddRecordToOutputVariableStructure('*','DAYLIGHT ILLUM AT REF POINT 1')
    CALL AddRecordToOutputVariableStructure('*','GLARE INDEX AT REF POINT 1')
    CALL AddRecordToOutputVariableStructure('*','TIME EXCEEDING GLARE INDEX SETPOINT AT REF POINT 1')
    CALL AddRecordToOutputVariableStructure('*','TIME EXCEEDING DAYLIGHT ILLUMINANCE SETPOINT AT REF POINT 1')
    CALL AddRecordToOutputVariableStructure('*','DAYLIGHT ILLUM AT REF POINT 2')
    CALL AddRecordToOutputVariableStructure('*','GLARE INDEX AT REF POINT 2')
    CALL AddRecordToOutputVariableStructure('*','TIME EXCEEDING GLARE INDEX SETPOINT AT REF POINT 2')
    CALL AddRecordToOutputVariableStructure('*','TIME EXCEEDING DAYLIGHT ILLUMINANCE SETPOINT AT REF POINT 2')

  CASE ('COILREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','TOTAL WATER HEATING COIL ENERGY')
    CALL AddRecordToOutputVariableStructure('*','TOTAL WATER HEATING COIL RATE')
    CALL AddRecordToOutputVariableStructure('*','SENSIBLE WATER COOLING COIL ENERGY')
    CALL AddRecordToOutputVariableStructure('*','TOTAL WATER COOLING COIL ENERGY')
    CALL AddRecordToOutputVariableStructure('*','TOTAL WATER COOLING COIL RATE')
    CALL AddRecordToOutputVariableStructure('*','SENSIBLE WATER COOLING COIL RATE')
    CALL AddRecordToOutputVariableStructure('*','COOLING COIL AREA WET FRACTION')

  CASE ('PLANTLOOPDEMANDREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','PLANT LOOP COOLING DEMAND')
    CALL AddRecordToOutputVariableStructure('*','PLANT LOOP HEATING DEMAND')

  CASE ('FANREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','FAN ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','FAN DELTA TEMP')
    CALL AddRecordToOutputVariableStructure('*','FAN ELECTRIC POWER')

  CASE ('PUMPREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','PUMP ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','PUMP HEAT TO FLUID ENERGY')
    CALL AddRecordToOutputVariableStructure('*','PUMP ELECTRIC POWER')
    CALL AddRecordToOutputVariableStructure('*','PUMP SHAFT POWER')
    CALL AddRecordToOutputVariableStructure('*','PUMP HEAT TO FLUID')
    CALL AddRecordToOutputVariableStructure('*','PUMP OUTLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','PUMP MASS FLOW RATE')

  CASE ('CONDLOOPDEMANDREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','COND LOOP COOLING DEMAND')
    CALL AddRecordToOutputVariableStructure('*','COND LOOP INLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','COND LOOP OUTLET TEMP')
    CALL AddRecordToOutputVariableStructure('*','COND LOOP HEATING DEMAND')

  CASE ('ZONETEMPERATUREOSCILLATIONREPORTMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','TIME ZONE TEMPERATURE OSCILLATING')
    CALL AddRecordToOutputVariableStructure('*','ZONE PEOPLE NUMBER OF OCCUPANTS')

  CASE ('AIRLOOPSYSTEMENERGYANDWATERUSEMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HOT WATER CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP STEAM CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP CHILLED WATER CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP GAS CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP WATER CONSUMPTION')

  CASE ('AIRLOOPSYSTEMCOMPONENTLOADSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP FAN HEATING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL COOLING COIL ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL HEATING COIL ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL HEAT EXCHANGER HEATING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL HEAT EXCHANGER COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL HUMIDIFIER HEATING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL EVAP COOLER COOLING ENERGY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP TOTAL DESICCANT DEHUMIDIFIER COOLING ENERGY')

  CASE ('AIRLOOPSYSTEMCOMPONENTENERGYUSEMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP FAN ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HEATING COIL HOT WATER CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP COOLING COIL CHILLED WATER CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP DX HEATING COIL ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP DX COOLING COIL ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HEATING COIL ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HEATING COIL GAS CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HEATING COIL STEAM CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP HUMIDIFIER ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP EVAP COOLER ELECTRIC CONSUMPTION')
    CALL AddRecordToOutputVariableStructure('*','AIR LOOP DESICCANT DEHUMIDIFIER ELECTRIC CONSUMPTION')

  CASE ('MECHANICALVENTILATIONLOADSMONTHLY')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION NO LOAD HEAT REMOVAL')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION COOLING LOAD INCREASE')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECH VENTILATION COOLING LOAD INCREASE: OVERHEATING')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION COOLING LOAD DECREASE')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION NO LOAD HEAT ADDITION')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION HEATING LOAD INCREASE')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECH VENTILATION HEATING LOAD INCREASE: OVERCOOLING')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION HEATING LOAD DECREASE')
    CALL AddRecordToOutputVariableStructure('*','ZONE MECHANICAL VENTILATION AIR CHANGE RATE')

  CASE DEFAULT

  END SELECT

  RETURN

END SUBROUTINE AddVariablesForMonthlyReport

FUNCTION FindFirstRecord(UCObjType) RESULT (StartPointer)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! Finds next record of Object Name

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: UCObjType
  INTEGER                      :: StartPointer

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: Found

  IF (SortedIDD) THEN
    Found=FindIteminSortedList(UCObjType,ListofObjects,NumObjectDefs)
    IF (Found /= 0) Found=iListofObjects(Found)
  ELSE
    Found=FindIteminList(UCObjType,ListofObjects,NumObjectDefs)
  ENDIF

  IF (Found /= 0) THEN
    StartPointer=ObjectStartRecord(Found)
  ELSE
    StartPointer=0
  ENDIF

  RETURN

END FUNCTION FindFirstRecord

FUNCTION FindNextRecord(UCObjType,StartPointer) RESULT (NextPointer)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! Finds next record of Object Name

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: UCObjType
  INTEGER, INTENT(IN)          :: StartPointer
  INTEGER                      :: NextPointer

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  INTEGER :: ObjNum

  NextPointer=0
  DO ObjNum=StartPointer+1,NumIDFRecords
    IF (IDFRecords(ObjNum)%Name /= UCObjType) CYCLE
    NextPointer=ObjNum
    EXIT
  END DO

  RETURN

END FUNCTION FindNextRecord

SUBROUTINE AddRecordToOutputVariableStructure(KeyValue,VariableName)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   July 2010
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine adds a new record (if necessary) to the Output Variable
  ! reporting structure.  DataOutputs, OutputVariablesForSimulation

  ! METHODOLOGY EMPLOYED:
  ! OutputVariablesForSimulation is a linked list structure for later
  ! semi-easy perusal.

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataOutputs

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: KeyValue
  CHARACTER(len=*), INTENT(IN) :: VariableName

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER :: CurNum
  INTEGER :: NextNum
  LOGICAL :: FoundOne
  INTEGER :: vnameLen  ! if < length, there were units on the line/name

  vnameLen=INDEX(VariableName,'[')
  IF (vnameLen == 0) THEN
    vnameLen=LEN_TRIM(VariableName)
  ELSE
    vnameLen=vnameLen-1
    vnameLen=LEN_TRIM(VariableName(1:vnameLen))
  ENDIF

  FoundOne=.false.
  DO CurNum=1,NumConsideredOutputVariables
    IF (VariableName(1:vnameLen) == OutputVariablesForSimulation(CurNum)%VarName) THEN
      FoundOne=.true.
      EXIT
    ENDIF
  ENDDO

  IF (.not. FoundOne) THEN
    IF (NumConsideredOutputVariables == MaxConsideredOutputVariables) THEN
      CALL ReAllocateAndPreserveOutputVariablesForSimulation
    ENDIF
    NumConsideredOutputVariables=NumConsideredOutputVariables+1
    OutputVariablesForSimulation(NumConsideredOutputVariables)%Key=KeyValue
    OutputVariablesForSimulation(NumConsideredOutputVariables)%VarName=VariableName(1:vnameLen)
    OutputVariablesForSimulation(NumConsideredOutputVariables)%Previous=0
    OutputVariablesForSimulation(NumConsideredOutputVariables)%Next=0
  ELSE
    IF (KeyValue /= OutputVariablesForSimulation(CurNum)%Key) THEN
      NextNum=CurNum
      IF (OutputVariablesForSimulation(NextNum)%Next /= 0) THEN
        DO WHILE (OutputVariablesForSimulation(NextNum)%Next /= 0)
          CurNum=NextNum
          NextNum=OutputVariablesForSimulation(NextNum)%Next
        ENDDO
        IF (NumConsideredOutputVariables == MaxConsideredOutputVariables) THEN
          CALL ReAllocateAndPreserveOutputVariablesForSimulation
        ENDIF
        NumConsideredOutputVariables=NumConsideredOutputVariables+1
        OutputVariablesForSimulation(NumConsideredOutputVariables)%Key=KeyValue
        OutputVariablesForSimulation(NumConsideredOutputVariables)%VarName=VariableName(1:vnameLen)
        OutputVariablesForSimulation(NumConsideredOutputVariables)%Previous=NextNum
        OutputVariablesForSimulation(NextNum)%Next=NumConsideredOutputVariables
      ELSE
        IF (NumConsideredOutputVariables == MaxConsideredOutputVariables) THEN
          CALL ReAllocateAndPreserveOutputVariablesForSimulation
        ENDIF
        NumConsideredOutputVariables=NumConsideredOutputVariables+1
        OutputVariablesForSimulation(NumConsideredOutputVariables)%Key=KeyValue
        OutputVariablesForSimulation(NumConsideredOutputVariables)%VarName=VariableName(1:vnameLen)
        OutputVariablesForSimulation(NumConsideredOutputVariables)%Previous=CurNum
        OutputVariablesForSimulation(CurNum)%Next=NumConsideredOutputVariables
      ENDIF
    ENDIF
  ENDIF

  RETURN

END SUBROUTINE AddRecordToOutputVariableStructure

SUBROUTINE ReAllocateAndPreserveOutputVariablesForSimulation

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   April 2011
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This routine does a simple reallocate for the OutputVariablesForSimulation structure, preserving
  ! the data that is already in the structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  USE DataOutputs

  IMPLICIT NONE ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  ! na

  ! SUBROUTINE PARAMETER DEFINITIONS:
  INTEGER, PARAMETER :: OutputVarAllocInc=ObjectsIDFAllocInc

  ! INTERFACE BLOCK SPECIFICATIONS:
  ! na

  ! DERIVED TYPE DEFINITIONS:
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  ! up allocation by OutputVarAllocInc
  ALLOCATE(TempOutputVariablesForSimulation(MaxConsideredOutputVariables+OutputVarAllocInc))
  TempOutputVariablesForSimulation(1:NumConsideredOutputVariables)=OutputVariablesForSimulation(1:NumConsideredOutputVariables)
  DEALLOCATE(OutputVariablesForSimulation)
  MaxConsideredOutputVariables=MaxConsideredOutputVariables+OutputVarAllocInc
  ALLOCATE(OutputVariablesForSimulation(MaxConsideredOutputVariables))
  OutputVariablesForSimulation(1:MaxConsideredOutputVariables)=TempOutputVariablesForSimulation(1:MaxConsideredOutputVariables)
  DEALLOCATE(TempOutputVariablesForSimulation)

  RETURN

END SUBROUTINE ReAllocateAndPreserveOutputVariablesForSimulation

SUBROUTINE DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,CurLine,NumConxLines,LineBuf,CurQPtr)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda Lawrie
  !       DATE WRITTEN   February 2003
  !       MODIFIED       March 2012 - Que lines instead of holding all.
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine dumps the "context" lines for error messages detected by
  ! the input processor.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)                        :: StartLine
  CHARACTER(len=*), INTENT(IN)               :: cStartLine
  CHARACTER(len=*), INTENT(IN)               :: cStartName
  INTEGER, INTENT(IN)                        :: CurLine
  INTEGER, INTENT(IN)                        :: NumConxLines
  CHARACTER(len=*), INTENT(IN), DIMENSION(:) :: LineBuf
  INTEGER, INTENT(IN)                        :: CurQPtr

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Line
  !  INTEGER PLine
  INTEGER SLine
  INTEGER CurPos
  CHARACTER(len=32) :: cLineNo
  CHARACTER(len=300) TextLine

  CALL ShowMessage('IDF Context for following error/warning message:')
  CALL ShowMessage('Note -- lines truncated at 300 characters, if necessary...')
  IF (StartLine <= 99999) THEN
    WRITE(TextLine,'(1X,I5,1X,A)') StartLine,trim(cStartLine)
  ELSE
    WRITE(cLineNo,*) StartLine
    cLineNo=ADJUSTL(cLineNo)
    WRITE(TextLine,'(1X,A,1X,A)') trim(cLineNo),trim(cStartLine)
  ENDIF
  CALL ShowMessage(TRIM(TextLine))
  IF (cStartName /= Blank) THEN
    CALL ShowMessage('indicated Name='//trim(cStartName))
  ENDIF
  CALL ShowMessage('Only last '//trim(IPTrimSigDigits(NumConxLines))//' lines before error line shown.....')
  SLine=CurLine-NumConxLines+1
  IF (NumConxLines == SIZE(LineBuf)) THEN
    CurPos=CurQPtr+1
    IF (CurQPtr+1 > SIZE(LineBuf)) CurPos=1
  ELSE
    CurPos=1
  ENDIF
  DO Line=1,NumConxLines
    IF (SLine <= 99999) THEN
      WRITE(TextLine,'(1X,I5,1X,A)') SLine,trim(LineBuf(CurPos))
    ELSE
      WRITE(cLineNo,*) SLine
      cLineNo=ADJUSTL(cLineNo)
      WRITE(TextLine,'(1X,A,1X,A)') trim(cLineNo),trim(LineBuf(CurPos))
    ENDIF
    CALL ShowMessage(TRIM(TextLine))
    CurPos=CurPos+1
    IF (CurPos > SIZE(LineBuf)) CurPos=1
    SLine=SLine+1
  ENDDO

  RETURN

END SUBROUTINE DumpCurrentLineBuffer

SUBROUTINE ShowAuditErrorMessage(Severity,ErrorMessage)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   March 2003
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine is just for messages that will be displayed on the audit trail
  ! (echo of the input file).  Errors are counted and a summary is displayed after
  ! finishing the scan of the input file.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*) Severity     ! if blank, does not add to sum
  CHARACTER(len=*) ErrorMessage

  ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: ErrorFormat='(2X,A)'

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  ! na

  IF (Severity /= Blank) THEN
    TotalAuditErrors=TotalAuditErrors+1
    WRITE(EchoInputFile,ErrorFormat) Severity//TRIM(ErrorMessage)
  ELSE
    WRITE(EchoInputFile,ErrorFormat) ' ************* '//TRIM(ErrorMessage)
  ENDIF


  RETURN

END SUBROUTINE ShowAuditErrorMessage

FUNCTION IPTrimSigDigits(IntegerValue) RESULT(OutputString)

  ! FUNCTION INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   March 2002
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function accepts a number as parameter as well as the number of
  ! significant digits after the decimal point to report and returns a string
  ! that is appropriate.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: IntegerValue
  CHARACTER(len=32) OutputString

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=32) String     ! Working string

  WRITE(String,*) IntegerValue
  OutputString=ADJUSTL(String)

  RETURN

END FUNCTION IPTrimSigDigits

SUBROUTINE DeallocateArrays !RS: Debugging: Carried over from HPSim InputProcessor

  IF (ALLOCATED(ObjectDef)) THEN
    DEALLOCATE(ObjectDef)
  END IF
  IF (ALLOCATED(SectionDef)) THEN
    DEALLOCATE(SectionDef)
  END IF
  IF (ALLOCATED(SectionsonFile )) THEN
    DEALLOCATE(SectionsonFile )
  END IF
  IF (ALLOCATED(IDFRecords)) THEN
    DEALLOCATE(IDFRecords)
  END IF
  IF (ALLOCATED(RepObjects)) THEN
    DEALLOCATE(RepObjects)
  END IF
  IF (ALLOCATED(ListofSections)) THEN
    DEALLOCATE(ListofSections)
  END IF
  IF (ALLOCATED(ListofObjects)) THEN
    DEALLOCATE(ListofObjects)
  END IF
  IF (ALLOCATED(ObsoleteObjectsRepNames)) THEN
    DEALLOCATE(ObsoleteObjectsRepNames)
  END IF
  IF (ALLOCATED(IDFRecordsGotten)) THEN
    DEALLOCATE(IDFRecordsGotten)
  END IF

  RETURN

END SUBROUTINE DeallocateArrays

!SUBROUTINE AddObjectDefandParse2(ProposedObject,CurPos,EndofFile,ErrorsFound)   !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)
!
!          ! SUBROUTINE INFORMATION:
!          !       AUTHOR         Linda K. Lawrie
!          !       DATE WRITTEN   August 1997
!          !       MODIFIED       na
!          !       RE-ENGINEERED  na
!
!          ! PURPOSE OF THIS SUBROUTINE:
!          ! This subroutine processes data dictionary file for EnergyPlus.
!          ! The structure of the sections and objects are stored in derived
!          ! types (SectionDefs and ObjectDefs)
!
!          ! METHODOLOGY EMPLOYED:
!          ! na
!
!          ! REFERENCES:
!          ! na
!
!          ! USE STATEMENTS:
!          ! na
!
!  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine
!
!          ! SUBROUTINE ARGUMENT DEFINITIONS
!  CHARACTER(len=*), INTENT(IN) :: ProposedObject  ! Proposed Object to Add
!  INTEGER, INTENT(INOUT) :: CurPos ! Current position (initially at first ',') of InputLine
!  LOGICAL, INTENT(INOUT) :: EndofFile ! End of File marker
!  LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if errors found here
!
!          ! SUBROUTINE PARAMETER DEFINITIONS:
!          ! na
!
!          ! INTERFACE BLOCK SPECIFICATIONS
!          ! na
!
!          ! DERIVED TYPE DEFINITIONS
!          ! na
!
!          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
!  CHARACTER(len=MaxObjectNameLength) SqueezedObject  ! Input Object, Left Justified, UpperCase
!  INTEGER Count  ! Count on arguments, loop
!  INTEGER Pos    ! Position scanning variable
!  LOGICAL EndofObjectDef   ! Set to true when ; has been found
!  LOGICAL ErrFlag   ! Local Error condition flag, when true, object not added to Global list
!  CHARACTER(len=1) TargetChar   ! Single character scanned to test for current field type (A or N)
!  LOGICAL BlankLine ! True when this line is "blank" (may have comment characters as first character on line)
!  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: AlphaorNumeric    ! Array of argument designations, True is Alpha,
!                                                                   ! False is numeric, saved in ObjectDef when done
!  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempAN            ! Array (ref: AlphaOrNumeric) for re-allocation procedure
!  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: RequiredFields    ! Array of argument required fields
!  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempRqF           ! Array (ref: RequiredFields) for re-allocation procedure
!  CHARACTER(len=MaxObjectNameLength),   &
!              ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldChecks   ! Array with alpha field names
!  CHARACTER(len=MaxObjectNameLength),   &
!              ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFC           ! Array (ref: AlphFieldChecks) for re-allocation procedure
!  CHARACTER(len=MaxObjectNameLength),   &
!              ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldDefaults ! Array with alpha field defaults
!  CHARACTER(len=MaxObjectNameLength),   &
!              ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFD           ! Array (ref: AlphFieldDefaults) for re-allocation procedure
!  TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: NumRangeChecks  ! Structure for Range Check, Defaults of numeric fields
!  TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: TempChecks ! Structure (ref: NumRangeChecks) for re-allocation procedure
!  LOGICAL MinMax   ! Set to true when MinMax field has been found by ReadInputLine
!  LOGICAL Default  ! Set to true when Default field has been found by ReadInputLine
!  LOGICAL AutoSize ! Set to true when Autosizable field has been found by ReadInputLine
!  CHARACTER(len=20) MinMaxString ! Set from ReadInputLine
!  CHARACTER(len=MaxObjectNameLength) AlphDefaultString
!  INTEGER WhichMinMax   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
!  REAL Value  ! Value returned by ReadInputLine (either min, max, default or autosize)
!  LOGICAL MinMaxError  ! Used to see if min, max, defaults have been set appropriately (True if error)
!  INTEGER,SAVE   :: MaxANArgs=100  ! Current count of Max args to object
!  LOGICAL ErrorsFoundFlag
!
!  IF (.not. ALLOCATED(AlphaorNumeric)) THEN
!    ALLOCATE (AlphaorNumeric(0:MaxANArgs))
!    ALLOCATE (RequiredFields(0:MaxANArgs))
!    ALLOCATE (NumRangeChecks(MaxANArgs))
!    ALLOCATE (AlphFieldChecks(MaxANArgs))
!    ALLOCATE (AlphFieldDefaults(MaxANArgs))
!    ALLOCATE (ObsoleteObjectsRepNames2(0))
!  ENDIF
!
!  SqueezedObject=MakeUPPERCase(ADJUSTL(ProposedObject))
!  IF (LEN_TRIM(ADJUSTL(ProposedObject)) > MaxObjectNameLength) THEN
!    CALL ShowWarningError('Object length exceeds maximum, will be truncated='//TRIM(ProposedObject),EchoInputFile)
!    CALL ShowContinueError('Will be processed as Object='//TRIM(SqueezedObject),EchoInputFile)
!    ErrorsFound=.true.
!  ENDIF
!
!  ! Start of Object parse, set object level items
!  ErrFlag=.false.
!  MinimumNumberOfFields=0
!  ObsoleteObject=.false.
!  UniqueObject=.false.
!  RequiredObject=.false.
!
!  IF (SqueezedObject /= Blank) THEN
!    IF (FindItemInList(SqueezedObject,ObjectDef2%Name,NumObjectDefs2) > 0) THEN
!      CALL ShowSevereError('Already an Object called '//TRIM(SqueezedObject)//'. This definition ignored.',EchoInputFile)
!      ! Error Condition
!      ErrFlag=.true.
!      ! Rest of Object has to be processed. Error condition will be caught
!      ! at end
!      ErrorsFound=.true.
!    ENDIF
!  ELSE
!    ErrFlag=.true.
!    ErrorsFound=.true.
!  ENDIF
!
!  NumObjectDefs2=NumObjectDefs2+1
!  ObjectDef2(NumObjectDefs2)%Name=SqueezedObject
!  ObjectDef2(NumObjectDefs2)%NumParams=0
!  ObjectDef2(NumObjectDefs2)%NumAlpha=0
!  ObjectDef2(NumObjectDefs2)%NumNumeric=0
!  ObjectDef2(NumObjectDefs2)%NumFound=0
!  ObjectDef2(NumObjectDefs2)%MinNumFields=0
!  ObjectDef2(NumObjectDefs2)%NameAlpha1=.false.
!  ObjectDef2(NumObjectDefs2)%ObsPtr=0
!  ObjectDef2(NumObjectDefs2)%UniqueObject=.false.
!  ObjectDef2(NumObjectDefs2)%RequiredObject=.false.
!
!  AlphaorNumeric=.true.
!  RequiredFields=.false.
!  AlphFieldChecks=Blank
!  AlphFieldDefaults=Blank
!
!  NumRangeChecks%MinMaxChk=.false.
!  NumRangeChecks%WhichMinMax(1)=0
!  NumRangeChecks%WhichMinMax(2)=0
!  NumRangeChecks%MinMaxString(1)=Blank
!  NumRangeChecks%MinMaxString(2)=Blank
!  NumRangeChecks%MinMaxValue(1)=0.0
!  NumRangeChecks%MinMaxValue(2)=0.0
!  NumRangeChecks%Default=0.0
!  NumRangeChecks%DefaultChk=.false.
!  NumRangeChecks%DefAutoSize=.false.
!  NumRangeChecks%FieldName=Blank
!  NumRangeChecks%AutoSizable=.false.
!  NumRangeChecks%AutoSizeValue=DefAutoSizeValue
!
!  Count=0
!  EndofObjectDef=.false.
!  ! Parse rest of Object Definition
!
!  DO WHILE (.not. EndofFile .and. .not. EndofObjectDef)
!
!    IF (CurPos <= InputLineLength) THEN
!      Pos=SCAN(InputLine(CurPos:InputLineLength),AlphaNum)
!      IF (Pos > 0) then
!
!        Count=Count+1
!        RequiredField=.false.
!
!        IF (Count > MaxANArgs) THEN   ! Reallocation
!          ALLOCATE(TempAN(0:MaxANArgs+ObjectDefAllocInc))
!          TempAN=.false.
!          TempAN(0:MaxANArgs)=AlphaorNumeric
!          DEALLOCATE(AlphaorNumeric)
!          ALLOCATE(TempRqF(0:MaxANArgs+ObjectDefAllocInc))
!          TempRqF=.false.
!          TempRqF(1:MaxANArgs)=RequiredFields
!          DEALLOCATE(RequiredFields)
!          ALLOCATE(TempChecks(MaxANArgs+ObjectDefAllocInc))
!          TempChecks%MinMaxChk=.false.
!          TempChecks%WhichMinMax(1)=0
!          TempChecks%WhichMinMax(2)=0
!          TempChecks%MinMaxString(1)=Blank
!          TempChecks%MinMaxString(2)=Blank
!          TempChecks%MinMaxValue(1)=0.0
!          TempChecks%MinMaxValue(2)=0.0
!          TempChecks%Default=0.0
!          TempChecks%DefaultChk=.false.
!          TempChecks%DefAutoSize=.false.
!          TempChecks%FieldName=Blank
!          TempChecks(1:MaxANArgs)=NumRangeChecks(1:MaxANArgs)
!          DEALLOCATE(NumRangeChecks)
!          ALLOCATE(TempAFC(MaxANArgs+ObjectDefAllocInc))
!          TempAFC=Blank
!          TempAFC(1:MaxANArgs)=AlphFieldChecks
!          DEALLOCATE(AlphFieldChecks)
!          ALLOCATE(TempAFD(MaxANArgs+ObjectDefAllocInc))
!          TempAFD=Blank
!          TempAFD(1:MaxANArgs)=AlphFieldDefaults
!          DEALLOCATE(AlphFieldDefaults)
!          ALLOCATE(AlphaorNumeric(0:MaxANArgs+ObjectDefAllocInc))
!          AlphaorNumeric=TempAN
!          DEALLOCATE(TempAN)
!          ALLOCATE(RequiredFields(0:MaxANArgs+ObjectDefAllocInc))
!          RequiredFields=TempRqF
!          DEALLOCATE(TempRqF)
!          ALLOCATE(NumRangeChecks(MaxANArgs+ObjectDefAllocInc))
!          NumRangeChecks=TempChecks
!          DEALLOCATE(TempChecks)
!          ALLOCATE(AlphFieldChecks(MaxANArgs+ObjectDefAllocInc))
!          AlphFieldChecks=TempAFC
!          DEALLOCATE(TempAFC)
!          ALLOCATE(AlphFieldDefaults(MaxANArgs+ObjectDefAllocInc))
!          AlphFieldDefaults=TempAFD
!          DEALLOCATE(TempAFD)
!          MaxANArgs=MaxANArgs+ObjectDefAllocInc
!        ENDIF
!
!        TargetChar=InputLine(CurPos+Pos-1:CurPos+Pos-1)
!
!        IF (TargetChar == 'A' .or. TargetChar == 'a') THEN
!          AlphaorNumeric(Count)=.true.
!          ObjectDef2(NumObjectDefs2)%NumAlpha=ObjectDef2(NumObjectDefs2)%NumAlpha+1
!          IF (FieldSet) THEN
!              AlphFieldChecks(ObjectDef2(NumObjectDefs2)%NumAlpha)=CurrentFieldName
!          END IF
!          IF (ObjectDef2(NumObjectDefs2)%NumAlpha == 1) THEN
!            IF (INDEX(MakeUpperCase(CurrentFieldName),'NAME') /= 0) THEN
!                ObjectDef2(NumObjectDefs2)%NameAlpha1=.true.
!            END IF
!          ENDIF
!        ELSE
!          AlphaorNumeric(Count)=.false.
!          ObjectDef2(NumObjectDefs2)%NumNumeric=ObjectDef2(NumObjectDefs2)%NumNumeric+1
!          IF (FieldSet) THEN
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldName=CurrentFieldName
!          END IF
!        ENDIF
!
!      ELSE
!        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
!                           MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
!                           Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
!                           ErrorsFound=ErrorsFoundFlag)
!        IF (.not. AlphaorNumeric(Count)) THEN
!          ! only record for numeric fields
!          IF (MinMax) THEN
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxChk=.true.
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldNumber=Count
!            IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(1)=WhichMinMax
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(1)=MinMaxString
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(1)=Value
!            ELSE
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(2)=WhichMinMax
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(2)=MinMaxString
!              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(2)=Value
!            ENDIF
!          ENDIF   ! End Min/Max
!          IF (Default) THEN
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefaultChk=.true.
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%Default=Value
!            IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoSize=.true.
!          ENDIF
!          IF (AutoSize) THEN
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizable=.true.
!            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizeValue=Value
!          ENDIF
!        ELSE  ! Alpha Field
!          IF (Default) THEN
!            AlphFieldDefaults(ObjectDef2(NumObjectDefs2)%NumAlpha)=AlphDefaultString
!          ENDIF
!        ENDIF
!        IF (ErrorsFoundFlag) THEN
!          ErrFlag=.true.
!          ErrorsFoundFlag=.false.
!        ENDIF
!        IF (RequiredField) THEN
!          RequiredFields(Count)=.true.
!          MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
!        ENDIF
!        CYCLE
!      ENDIF
!
!      !  For the moment dont care about descriptions on each object
!      IF (CurPos <= InputLineLength) THEN
!        CurPos=CurPos+Pos
!        Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
!      ELSE
!        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
!        IF (BlankLine .or. EndofFile) THEN
!            CYCLE
!        END IF
!        Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
!      ENDIF
!    ELSE
!      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
!      CYCLE
!    ENDIF
!
!    IF (Pos <= 0) THEN
!                   ! must be time to read another line
!      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
!      IF (BlankLine .or. EndofFile) THEN
!          CYCLE
!      END IF
!    ELSE
!      IF (InputLine(CurPos+Pos-1:CurPos+Pos-1) == ';') THEN
!        EndofObjectDef=.true.
!      ENDIF
!      CurPos=CurPos+Pos
!    ENDIF
!
!  END DO
!
!  ! Reached end of object def but there may still be more \ lines to parse....
!  ! Goes until next object is encountered ("not blankline") or end of IDDFile
!  ! If last object is not numeric, then exit immediately....
!    BlankLine=.true.
!    DO WHILE (BlankLine .and. .not.EndofFile)
!    ! It's a numeric object as last one...
!      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
!                         MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
!                         Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
!                         ErrorsFound=ErrorsFoundFlag)
!      IF (MinMax) THEN
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxChk=.true.
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldNumber=Count
!        IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(1)=WhichMinMax
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(1)=MinMaxString
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(1)=Value
!        ELSE
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(2)=WhichMinMax
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(2)=MinMaxString
!          NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(2)=Value
!        ENDIF
!      ENDIF
!      IF (Default .and. .not. AlphaorNumeric(Count)) THEN
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefaultChk=.true.
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%Default=Value
!        IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoSize=.true.
!      ELSEIF (Default .and. AlphaorNumeric(Count)) THEN
!        AlphFieldDefaults(ObjectDef2(NumObjectDefs2)%NumAlpha)=AlphDefaultString
!      ENDIF
!      IF (AutoSize) THEN
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizable=.true.
!        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizeValue=Value
!      ENDIF
!      IF (ErrorsFoundFlag) THEN
!        ErrFlag=.true.
!        ErrorsFoundFlag=.false.
!      ENDIF
!    ENDDO
!    IF (.not. BlankLine) THEN
!      BACKSPACE(Unit=IDDFile)
!      EchoInputLine=.false.
!    ENDIF
!  IF (RequiredField) THEN
!    RequiredFields(Count)=.true.
!    MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
!  ENDIF
!
!  ObjectDef2(NumObjectDefs2)%NumParams=Count  ! Also the total of ObjectDef(..)%NumAlpha+ObjectDef(..)%NumNumeric
!  ObjectDef2(NumObjectDefs2)%MinNumFields=MinimumNumberOfFields
!  IF (ObsoleteObject) THEN
!    ALLOCATE(TempAFD(NumObsoleteObjects+1))
!    IF (NumObsoleteObjects > 0) THEN
!      TempAFD(1:NumObsoleteObjects)=ObsoleteObjectsRepNames2
!    ENDIF
!    TempAFD(NumObsoleteObjects+1)=ReplacementName
!    DEALLOCATE(ObsoleteObjectsRepNames2)
!    NumObsoleteObjects=NumObsoleteObjects+1
!    ALLOCATE(ObsoleteObjectsRepNames2(NumObsoleteObjects))
!    ObsoleteObjectsRepNames2=TempAFD
!    ObjectDef2(NumObjectDefs)%ObsPtr=NumObsoleteObjects
!    DEALLOCATE(TempAFD)
!  ENDIF
!  IF (RequiredObject) THEN
!    ObjectDef2(NumObjectDefs)%RequiredObject=.true.
!  ENDIF
!  IF (UniqueObject) THEN
!    ObjectDef2(NumObjectDefs)%UniqueObject=.true.
!  ENDIF
!
!  MaxAlphaArgsFound=MAX(MaxAlphaArgsFound,ObjectDef2(NumObjectDefs2)%NumAlpha)
!  MaxNumericArgsFound=MAX(MaxNumericArgsFound,ObjectDef2(NumObjectDefs2)%NumNumeric)
!  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphaorNumeric(Count))
!  ObjectDef2(NumObjectDefs2)%AlphaorNumeric=AlphaorNumeric(1:Count)
!  ALLOCATE(ObjectDef2(NumObjectDefs2)%NumRangeChks(ObjectDef2(NumObjectDefs2)%NumNumeric))
!  ObjectDef2(NumObjectDefs2)%NumRangeChks=NumRangeChecks(1:ObjectDef2(NumObjectDefs2)%NumNumeric)
!  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldChks(ObjectDef2(NumObjectDefs2)%NumAlpha))
!  ObjectDef2(NumObjectDefs2)%AlphFieldChks=AlphFieldChecks(1:ObjectDef2(NumObjectDefs2)%NumAlpha)
!  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldDefs(ObjectDef2(NumObjectDefs2)%NumAlpha))
!  ObjectDef2(NumObjectDefs2)%AlphFieldDefs=AlphFieldDefaults(1:ObjectDef2(NumObjectDefs2)%NumAlpha)
!  ALLOCATE(ObjectDef2(NumObjectDefs2)%ReqField(Count))
!  ObjectDef2(NumObjectDefs2)%ReqField=RequiredFields(1:Count)
!  DO Count=1,ObjectDef2(NumObjectDefs2)%NumNumeric
!    IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxChk) THEN
!    ! Checking MinMax Range (min vs. max and vice versa)
!      MinMaxError=.false.
!      ! check min against max
!      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
!        ! min
!        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)
!        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
!          IF (Value > ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!              MinMaxError=.true.
!          END IF
!        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
!          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!              MinMaxError=.true.
!          END IF
!        ENDIF
!      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
!        ! min>
!        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1) + TINY(Value)  ! infintesimally bigger than min
!        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
!          IF (Value > ObjectDef(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!              MinMaxError=.true.
!          END IF
!        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
!          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!              MinMaxError=.true.
!          END IF
!        ENDIF
!      ENDIF
!      ! check max against min
!      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
!        ! max
!        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)
!        ! Check max value against min
!        IF (ObjectDef2(NumObjectDefs)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
!          IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!              MinMaxError=.true.
!          END IF
!        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
!          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!              MinMaxError=.true.
!          END IF
!        ENDIF
!      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
!        ! max<
!        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2) - TINY(Value)  ! infintesimally bigger than min
!        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
!          IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!              MinMaxError=.true.
!          END IF
!        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
!          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!              MinMaxError=.true.
!          END IF
!        ENDIF
!      ENDIF
!      ! check if error condition
!      IF (MinMaxError) THEN
!        !  Error stated min is not in range with stated max
!        WRITE(MinMaxString,*) ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%FieldNumber
!        MinMaxString=ADJUSTL(MinMaxString)
!        CALL ShowSevereError('Field #'//TRIM(MinMaxString)//' conflict in Min/Max specifications/values, in class='//  &
!                             TRIM(ObjectDef2(NumObjectDefs2)%Name),EchoInputFile)
!        ErrFlag=.true.
!      ENDIF
!    ENDIF
!    IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%DefaultChk) THEN
!    ! Check Default against MinMaxRange
!      MinMaxError=.false.
!      Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%Default
!      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
!        IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!            MinMaxError=.true.
!        END IF
!      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
!        IF (Value <= ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) THEN
!            MinMaxError=.true.
!        END IF
!      ENDIF
!      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
!        IF (Value > ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!            MinMaxError=.true.
!        END IF
!      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
!        IF (Value >= ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) THEN
!            MinMaxError=.true.
!        END IF
!      ENDIF
!      IF (MinMaxError) THEN
!        !  Error stated default is not in min/max range
!        WRITE(MinMaxString,*) ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%FieldNumber
!        MinMaxString=ADJUSTL(MinMaxString)
!        CALL ShowSevereError('Field #'//TRIM(MinMaxString)//' default is invalid for Min/Max values, in class='//  &
!                             TRIM(ObjectDef2(NumObjectDefs)%Name),EchoInputFile)
!        ErrFlag=.true.
!      ENDIF
!    ENDIF
!  ENDDO
!
!  IF (ErrFlag) THEN
!    CALL ShowContinueError('Errors occured in ObjectDefinition for Class='//TRIM(ObjectDef2(NumObjectDefs)%Name)// &
!                           ', Object not available for IDF processing.',EchoInputFile)
!    DEALLOCATE(ObjectDef2(NumObjectDefs2)%AlphaorNumeric)
!    NumObjectDefs=NumObjectDefs-1
!    ErrorsFound=.true.
!  ENDIF
!
!  RETURN
!
!END SUBROUTINE AddObjectDefandParse2    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

SUBROUTINE AddObjectDefandParse2(ProposedObject,CurPos,EndofFile,ErrorsFound)    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   August 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine processes data dictionary file for EnergyPlus.
  ! The structure of the sections and objects are stored in derived
  ! types (SectionDefs and ObjectDefs)

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS
  CHARACTER(len=*), INTENT(IN) :: ProposedObject  ! Proposed Object to Add
  INTEGER, INTENT(INOUT) :: CurPos ! Current position (initially at first ',') of InputLine
  LOGICAL, INTENT(INOUT) :: EndofFile ! End of File marker
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if errors found here

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=MaxObjectNameLength) SqueezedObject  ! Input Object, Left Justified, UpperCase
  INTEGER Count  ! Count on arguments, loop
  INTEGER Pos    ! Position scanning variable
  LOGICAL EndofObjectDef   ! Set to true when ; has been found
  LOGICAL ErrFlag   ! Local Error condition flag, when true, object not added to Global list
  CHARACTER(len=1) TargetChar   ! Single character scanned to test for current field type (A or N)
  LOGICAL BlankLine ! True when this line is "blank" (may have comment characters as first character on line)
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: AlphaorNumeric    ! Array of argument designations, True is Alpha,
  ! False is numeric, saved in ObjectDef when done
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempAN            ! Array (ref: AlphaOrNumeric) for re-allocation procedure
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: RequiredFields    ! Array of argument required fields
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempRqF           ! Array (ref: RequiredFields) for re-allocation procedure
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: AlphRetainCase    ! Array of argument for retain case
  LOGICAL(1), ALLOCATABLE, SAVE, DIMENSION(:) :: TempRtC           ! Array (ref: AlphRetainCase) for re-allocation procedure
  CHARACTER(len=MaxFieldNameLength),   &
  ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldChecks   ! Array with alpha field names
  CHARACTER(len=MaxFieldNameLength),   &
  ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFC           ! Array (ref: AlphFieldChecks) for re-allocation procedure
  CHARACTER(len=MaxObjectNameLength),   &
  ALLOCATABLE, SAVE, DIMENSION(:) :: AlphFieldDefaults ! Array with alpha field defaults
  CHARACTER(len=MaxObjectNameLength),   &
  ALLOCATABLE, SAVE, DIMENSION(:) :: TempAFD           ! Array (ref: AlphFieldDefaults) for re-allocation procedure
  TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: NumRangeChecks  ! Structure for Range Check, Defaults of numeric fields
  TYPE(RangeCheckDef), ALLOCATABLE, SAVE, DIMENSION(:) :: TempChecks ! Structure (ref: NumRangeChecks) for re-allocation procedure
  LOGICAL MinMax   ! Set to true when MinMax field has been found by ReadInputLine
  LOGICAL Default  ! Set to true when Default field has been found by ReadInputLine
  LOGICAL AutoSize ! Set to true when Autosizable field has been found by ReadInputLine
  LOGICAL AutoCalculate ! Set to true when Autocalculatable field has been found by ReadInputLine
  CHARACTER(len=32) MinMaxString ! Set from ReadInputLine
  CHARACTER(len=MaxObjectNameLength) AlphDefaultString
  INTEGER WhichMinMax   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
  REAL(r64) Value  ! Value returned by ReadInputLine (either min, max, default, autosize or autocalculate)
  LOGICAL MinMaxError  ! Used to see if min, max, defaults have been set appropriately (True if error)
  INTEGER,SAVE   :: MaxANArgs=7700  ! Current count of Max args to object  (9/2010)
  LOGICAL ErrorsFoundFlag
  INTEGER,SAVE :: PrevSizeNumNumeric = -1
  INTEGER,SAVE :: PrevCount  = -1
  INTEGER,SAVE :: PrevSizeNumAlpha = -1
  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  IF (.not. ALLOCATED(AlphaorNumeric)) THEN
    ALLOCATE (AlphaorNumeric(0:MaxANArgs))
    ALLOCATE (RequiredFields(0:MaxANArgs))
    ALLOCATE (AlphRetainCase(0:MaxANArgs))
    ALLOCATE (NumRangeChecks(MaxANArgs))
    ALLOCATE (AlphFieldChecks(MaxANArgs))
    ALLOCATE (AlphFieldDefaults(MaxANArgs))
    ALLOCATE (ObsoleteObjectsRepNames2(0))
  ENDIF

  SqueezedObject=MakeUPPERCase(ADJUSTL(ProposedObject))
  IF (LEN_TRIM(ADJUSTL(ProposedObject)) > MaxObjectNameLength) THEN
    CALL ShowWarningError('IP: Object length exceeds maximum, will be truncated='//TRIM(ProposedObject),EchoInputFile)
    CALL ShowContinueError('Will be processed as Object='//TRIM(SqueezedObject),EchoInputFile)
    ErrorsFound=.true.
  ENDIF

  ! Start of Object parse, set object level items
  ErrFlag=.false.
  ErrorsFoundFlag=.false.
  MinimumNumberOfFields=0
  ObsoleteObject=.false.
  UniqueObject=.false.
  RequiredObject=.false.
  ExtensibleObject=.false.
  ExtensibleNumFields=0
  MinMax=.false.
  Default=.false.
  AutoSize=.false.
  AutoCalculate=.false.
  WhichMinMax=0


  IF (SqueezedObject /= Blank) THEN
    IF (FindItemInList(SqueezedObject,ObjectDef2%Name,NumObjectDefs2) > 0) THEN
      CALL ShowSevereError('IP: Already an Object called '//TRIM(SqueezedObject)//'. This definition ignored.',EchoInputFile)
      ! Error Condition
      ErrFlag=.true.
      ! Rest of Object has to be processed. Error condition will be caught
      ! at end
      ErrorsFound=.true.
    ENDIF
  ELSE
    ErrFlag=.true.
    ErrorsFound=.true.
  ENDIF

  NumObjectDefs2=NumObjectDefs2+1
  ObjectDef2(NumObjectDefs2)%Name=SqueezedObject
  ObjectDef2(NumObjectDefs2)%NumParams=0
  ObjectDef2(NumObjectDefs2)%NumAlpha=0
  ObjectDef2(NumObjectDefs2)%NumNumeric=0
  ObjectDef2(NumObjectDefs2)%NumFound=0
  ObjectDef2(NumObjectDefs2)%MinNumFields=0
  ObjectDef2(NumObjectDefs2)%NameAlpha1=.false.
  ObjectDef2(NumObjectDefs2)%ObsPtr=0
  ObjectDef2(NumObjectDefs2)%UniqueObject=.false.
  ObjectDef2(NumObjectDefs2)%RequiredObject=.false.
  ObjectDef2(NumObjectDefs2)%ExtensibleObject=.false.
  ObjectDef2(NumObjectDefs2)%ExtensibleNum=0

  IF (PrevCount .EQ. -1) THEN
    PrevCount = MaxANArgs
  END IF

  AlphaorNumeric(1:PrevCount)=.true.
  RequiredFields(1:PrevCount)=.false.
  AlphRetainCase(1:PrevCount)=.false.

  IF (PrevSizeNumAlpha .EQ. -1) THEN
    PrevSizeNumAlpha = MaxANArgs
  END IF

  AlphFieldChecks(1:PrevSizeNumAlpha)=Blank
  AlphFieldDefaults(1:PrevSizeNumAlpha)=Blank

  IF (PrevSizeNumNumeric .EQ. -1) THEN
    PrevSizeNumNumeric = MaxANArgs
  END IF

  !clear only portion of NumRangeChecks array used in the previous
  !call to reduce computation time to clear this large array.
  NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxChk=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%WhichMinMax(1)=0
  NumRangeChecks(1:PrevSizeNumNumeric)%WhichMinMax(2)=0
  NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxString(1)=Blank
  NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxString(2)=Blank
  NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxValue(1)=0.0
  NumRangeChecks(1:PrevSizeNumNumeric)%MinMaxValue(2)=0.0
  NumRangeChecks(1:PrevSizeNumNumeric)%Default=0.0
  NumRangeChecks(1:PrevSizeNumNumeric)%DefaultChk=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%DefAutoSize=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%DefAutoCalculate=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%FieldNumber=0
  NumRangeChecks(1:PrevSizeNumNumeric)%FieldName=Blank
  NumRangeChecks(1:PrevSizeNumNumeric)%AutoSizable=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%AutoSizeValue=DefAutoSizeValue
  NumRangeChecks(1:PrevSizeNumNumeric)%AutoCalculatable=.false.
  NumRangeChecks(1:PrevSizeNumNumeric)%AutoCalculateValue=DefAutoCalculateValue

  Count=0
  EndofObjectDef=.false.
  ! Parse rest of Object Definition

  DO WHILE (.not. EndofFile .and. .not. EndofObjectDef)

    IF (CurPos <= InputLineLength) THEN
      Pos=SCAN(InputLine(CurPos:InputLineLength),AlphaNum)
      IF (Pos > 0) then

        Count=Count+1
        RequiredField=.false.
        RetainCaseFlag=.false.

        IF (Count > MaxANArgs) THEN   ! Reallocation
          ALLOCATE(TempAN(0:MaxANArgs+ANArgsDefAllocInc))
          TempAN=.false.
          TempAN(0:MaxANArgs)=AlphaorNumeric
          DEALLOCATE(AlphaorNumeric)
          ALLOCATE(TempRqF(0:MaxANArgs+ANArgsDefAllocInc))
          TempRqF=.false.
          TempRqF(0:MaxANArgs)=RequiredFields
          DEALLOCATE(RequiredFields)
          ALLOCATE(TempRtC(0:MaxANArgs+ANArgsDefAllocInc))
          TempRtC=.false.
          TempRtC(0:MaxANArgs)=AlphRetainCase
          DEALLOCATE(AlphRetainCase)
          ALLOCATE(TempChecks(MaxANArgs+ANArgsDefAllocInc))
          TempChecks(1:MaxANArgs)=NumRangeChecks(1:MaxANArgs)
          DEALLOCATE(NumRangeChecks)
          ALLOCATE(TempAFC(MaxANArgs+ANArgsDefAllocInc))
          TempAFC=Blank
          TempAFC(1:MaxANArgs)=AlphFieldChecks
          DEALLOCATE(AlphFieldChecks)
          ALLOCATE(TempAFD(MaxANArgs+ANArgsDefAllocInc))
          TempAFD=Blank
          TempAFD(1:MaxANArgs)=AlphFieldDefaults
          DEALLOCATE(AlphFieldDefaults)
          ALLOCATE(AlphaorNumeric(0:MaxANArgs+ANArgsDefAllocInc))
          AlphaorNumeric=TempAN
          DEALLOCATE(TempAN)
          ALLOCATE(RequiredFields(0:MaxANArgs+ANArgsDefAllocInc))
          RequiredFields=TempRqF
          DEALLOCATE(TempRqF)
          ALLOCATE(AlphRetainCase(0:MaxANArgs+ANArgsDefAllocInc))
          AlphRetainCase=TempRtC
          DEALLOCATE(TempRtC)
          ALLOCATE(NumRangeChecks(MaxANArgs+ANArgsDefAllocInc))
          NumRangeChecks=TempChecks
          DEALLOCATE(TempChecks)
          ALLOCATE(AlphFieldChecks(MaxANArgs+ANArgsDefAllocInc))
          AlphFieldChecks=TempAFC
          DEALLOCATE(TempAFC)
          ALLOCATE(AlphFieldDefaults(MaxANArgs+ANArgsDefAllocInc))
          AlphFieldDefaults=TempAFD
          DEALLOCATE(TempAFD)
          MaxANArgs=MaxANArgs+ANArgsDefAllocInc
        ENDIF

        TargetChar=InputLine(CurPos+Pos-1:CurPos+Pos-1)

        IF (TargetChar == 'A' .or. TargetChar == 'a') THEN
          AlphaorNumeric(Count)=.true.
          ObjectDef2(NumObjectDefs2)%NumAlpha=ObjectDef2(NumObjectDefs2)%NumAlpha+1
          IF (FieldSet) AlphFieldChecks(ObjectDef2(NumObjectDefs2)%NumAlpha)=CurrentFieldName
          IF (ObjectDef2(NumObjectDefs2)%NumAlpha == 1) THEN
            IF (INDEX(MakeUpperCase(CurrentFieldName),'NAME') /= 0) ObjectDef2(NumObjectDefs2)%NameAlpha1=.true.
          ENDIF
        ELSE
          AlphaorNumeric(Count)=.false.
          ObjectDef2(NumObjectDefs2)%NumNumeric=ObjectDef2(NumObjectDefs2)%NumNumeric+1
          IF (FieldSet) NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldName=CurrentFieldName
        ENDIF

      ELSE
        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
        MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
        Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
        AutoCalculatable=AutoCalculate,RetainCase=RetainCaseFlag,ErrorsFound=ErrorsFoundFlag)
        IF (.not. AlphaorNumeric(Count)) THEN
          ! only record for numeric fields
          IF (MinMax) THEN
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxChk=.true.
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldNumber=Count
            IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(1)=WhichMinMax
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(1)=MinMaxString
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(1)=Value
            ELSE
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(2)=WhichMinMax
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(2)=MinMaxString
              NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(2)=Value
            ENDIF
          ENDIF   ! End Min/Max
          IF (Default) THEN
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefaultChk=.true.
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%Default=Value
            IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoSize=.true.
            IF (AlphDefaultString == 'AUTOCALCULATE')  NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoCalculate=.true.
          ENDIF
          IF (AutoSize) THEN
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizable=.true.
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizeValue=Value
          ENDIF
          IF (AutoCalculate) THEN
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoCalculatable=.true.
            NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoCalculateValue=Value
          ENDIF
        ELSE  ! Alpha Field
          IF (Default) THEN
            AlphFieldDefaults(ObjectDef2(NumObjectDefs2)%NumAlpha)=AlphDefaultString
          ENDIF
        ENDIF
        IF (ErrorsFoundFlag) THEN
          ErrFlag=.true.
          ErrorsFoundFlag=.false.
        ENDIF
        IF (RequiredField) THEN
          RequiredFields(Count)=.true.
          MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
        ENDIF
        IF (RetainCaseFlag) THEN
          AlphRetainCase(Count)=.true.
        ENDIF
        CYCLE
      ENDIF

      !  For the moment dont care about descriptions on each object
      IF (CurPos <= InputLineLength) THEN
        CurPos=CurPos+Pos
        Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
        IF (Pos == 0) THEN
          CALL ShowSevereError('IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' , or ; expected on this line'//  &
          ',position="'//InputLine(CurPos:InputLineLength)//'"',EchoInputFile)
          ErrFlag=.true.
          ErrorsFound=.true.
        ENDIF
        IF (InputLine(InputLineLength:InputLineLength) /= '\') THEN
          !CALL ShowWarningError('IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' \ expected on this line',EchoInputFile)    !RS: Secret Search String
          IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 13) THEN
            WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
          END IF
          WRITE(DebugFile,*) 'IP: IDD line~'//TRIM(IPTrimSigDigits(NumLines))//' \ expected on this line',EchoInputFile
        ENDIF
      ELSE
        CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
        IF (BlankLine .or. EndofFile) CYCLE
        Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
      ENDIF
    ELSE
      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
      CYCLE
    ENDIF

    IF (Pos <= 0) THEN
      ! must be time to read another line
      CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile)
      IF (BlankLine .or. EndofFile) CYCLE
    ELSE
      IF (InputLine(CurPos+Pos-1:CurPos+Pos-1) == ';') THEN
        EndofObjectDef=.true.
      ENDIF
      CurPos=CurPos+Pos
    ENDIF

  END DO

  ! Reached end of object def but there may still be more \ lines to parse....
  ! Goes until next object is encountered ("not blankline") or end of IDDFile
  ! If last object is not numeric, then exit immediately....
  BlankLine=.true.
  DO WHILE (BlankLine .and. .not.EndofFile)
    ! It's a numeric object as last one...
    CALL ReadInputLine(IDDFile,CurPos,BlankLine,InputLineLength,EndofFile,  &
    MinMax=MinMax,WhichMinMax=WhichMinMax,MinMaxString=MinMaxString,  &
    Value=Value,Default=Default,DefString=AlphDefaultString,AutoSizable=AutoSize, &
    AutoCalculatable=AutoCalculate,RetainCase=RetainCaseFlag,ErrorsFound=ErrorsFoundFlag)
    IF (MinMax) THEN
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxChk=.true.
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%FieldNumber=Count
      IF (WhichMinMax <= 2) THEN   !=0 (none/invalid), =1 \min, =2 \min>, =3 \max, =4 \max<
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(1)=WhichMinMax
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(1)=MinMaxString
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(1)=Value
      ELSE
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%WhichMinMax(2)=WhichMinMax
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxString(2)=MinMaxString
        NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%MinMaxValue(2)=Value
      ENDIF
    ENDIF
    IF (Default .and. .not. AlphaorNumeric(Count)) THEN
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefaultChk=.true.
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%Default=Value
      IF (AlphDefaultString == 'AUTOSIZE') NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoSize=.true.
      IF (AlphDefaultString == 'AUTOCALCULATE') NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%DefAutoCalculate=.true.
    ELSEIF (Default .and. AlphaorNumeric(Count)) THEN
      AlphFieldDefaults(ObjectDef2(NumObjectDefs2)%NumAlpha)=AlphDefaultString
    ENDIF
    IF (AutoSize) THEN
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizable=.true.
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoSizeValue=Value
    ENDIF
    IF (AutoCalculate) THEN
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoCalculatable=.true.
      NumRangeChecks(ObjectDef2(NumObjectDefs2)%NumNumeric)%AutoCalculateValue=Value
    ENDIF
    IF (ErrorsFoundFlag) THEN
      ErrFlag=.true.
      ErrorsFoundFlag=.false.
    ENDIF
  ENDDO

  IF (.not. BlankLine) THEN
    BACKSPACE(Unit=IDDFile)
    EchoInputLine=.false.
  ENDIF

  IF (RequiredField) THEN
    RequiredFields(Count)=.true.
    MinimumNumberOfFields=MAX(Count,MinimumNumberOfFields)
  ENDIF
  IF (RetainCaseFlag) THEN
    AlphRetainCase(Count)=.true.
  ENDIF

  ObjectDef2(NumObjectDefs2)%NumParams=Count  ! Also the total of ObjectDef(..)%NumAlpha+ObjectDef(..)%NumNumeric
  ObjectDef2(NumObjectDefs2)%MinNumFields=MinimumNumberOfFields
  IF (ObsoleteObject) THEN
    ALLOCATE(TempAFD(NumObsoleteObjects+1))
    IF (NumObsoleteObjects > 0) THEN
      TempAFD(1:NumObsoleteObjects)=ObsoleteObjectsRepNames2
    ENDIF
    TempAFD(NumObsoleteObjects+1)=ReplacementName
    DEALLOCATE(ObsoleteObjectsRepNames2)
    NumObsoleteObjects=NumObsoleteObjects+1
    ALLOCATE(ObsoleteObjectsRepNames2(NumObsoleteObjects))
    ObsoleteObjectsRepNames=TempAFD
    ObjectDef2(NumObjectDefs2)%ObsPtr=NumObsoleteObjects
    DEALLOCATE(TempAFD)
  ENDIF
  IF (RequiredObject) THEN
    ObjectDef2(NumObjectDefs2)%RequiredObject=.true.
  ENDIF
  IF (UniqueObject) THEN
    ObjectDef2(NumObjectDefs2)%UniqueObject=.true.
  ENDIF
  IF (ExtensibleObject) THEN
    ObjectDef2(NumObjectDefs2)%ExtensibleObject=.true.
    ObjectDef2(NumObjectDefs2)%ExtensibleNum=ExtensibleNumFields
  ENDIF

  MaxAlphaArgsFound=MAX(MaxAlphaArgsFound,ObjectDef2(NumObjectDefs2)%NumAlpha)
  MaxNumericArgsFound=MAX(MaxNumericArgsFound,ObjectDef2(NumObjectDefs2)%NumNumeric)
  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphaorNumeric(Count))
  ObjectDef2(NumObjectDefs2)%AlphaorNumeric=AlphaorNumeric(1:Count)
  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphRetainCase(Count))
  ObjectDef2(NumObjectDefs2)%AlphRetainCase=AlphRetainCase(1:Count)
  PrevCount = Count
  ALLOCATE(ObjectDef2(NumObjectDefs2)%NumRangeChks(ObjectDef2(NumObjectDefs2)%NumNumeric))
  IF (ObjectDef2(NumObjectDefs2)%NumNumeric > 0) THEN
    ObjectDef2(NumObjectDefs2)%NumRangeChks=NumRangeChecks(1:ObjectDef2(NumObjectDefs2)%NumNumeric)
  ENDIF
  PrevSizeNumNumeric = ObjectDef2(NumObjectDefs2)%NumNumeric !used to clear only portion of NumRangeChecks array
  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldChks(ObjectDef2(NumObjectDefs2)%NumAlpha))
  IF (ObjectDef2(NumObjectDefs2)%NumAlpha > 0) THEN
    ObjectDef2(NumObjectDefs2)%AlphFieldChks=AlphFieldChecks(1:ObjectDef2(NumObjectDefs2)%NumAlpha)
  ENDIF
  ALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldDefs(ObjectDef2(NumObjectDefs2)%NumAlpha))
  IF (ObjectDef2(NumObjectDefs2)%NumAlpha > 0) THEN
    ObjectDef2(NumObjectDefs2)%AlphFieldDefs=AlphFieldDefaults(1:ObjectDef2(NumObjectDefs2)%NumAlpha)
  ENDIF
  PrevSizeNumAlpha = ObjectDef2(NumObjectDefs2)%NumAlpha
  ALLOCATE(ObjectDef2(NumObjectDefs2)%ReqField(Count))
  ObjectDef2(NumObjectDefs2)%ReqField=RequiredFields(1:Count)
  DO Count=1,ObjectDef2(NumObjectDefs2)%NumNumeric
    IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxChk) THEN
      ! Checking MinMax Range (min vs. max and vice versa)
      MinMaxError=.false.
      ! check min against max
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
        ! min
        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)
        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
          IF (Value > ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ENDIF
      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
        ! min>
        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1) + rTinyValue  ! infintesimally bigger than min
        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
          IF (Value > ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
        ENDIF
      ENDIF
      ! check max against min
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
        ! max
        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)
        ! Check max value against min
        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
          IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ENDIF
      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
        ! max<
        Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2) - rTinyValue  ! infintesimally bigger than min
        IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
          IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
          IF (Value == ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
        ENDIF
      ENDIF
      ! check if error condition
      IF (MinMaxError) THEN
        !  Error stated min is not in range with stated max
        MinMaxString=IPTrimSigDigits(ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%FieldNumber)
        CALL ShowSevereError('IP: IDD: Field #'//TRIM(MinMaxString)//' conflict in Min/Max specifications/values, in class='//  &
        TRIM(ObjectDef2(NumObjectDefs2)%Name),EchoInputFile)
        ErrFlag=.true.
      ENDIF
    ENDIF
    IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%DefaultChk) THEN
      ! Check Default against MinMaxRange
      !  Don't check when default is autosize...
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%Autosizable .and.   &
      ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%DefAutoSize) CYCLE
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%Autocalculatable .and.   &
      ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%DefAutoCalculate) CYCLE
      MinMaxError=.false.
      Value=ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%Default
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 1) THEN
        IF (Value < ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(1) == 2) THEN
        IF (Value <= ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(1)) MinMaxError=.true.
      ENDIF
      IF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 3) THEN
        IF (Value > ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
      ELSEIF (ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%WhichMinMax(2) == 4) THEN
        IF (Value >= ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%MinMaxValue(2)) MinMaxError=.true.
      ENDIF
      IF (MinMaxError) THEN
        !  Error stated default is not in min/max range
        MinMaxString=IPTrimSigDigits(ObjectDef2(NumObjectDefs2)%NumRangeChks(Count)%FieldNumber)
        CALL ShowSevereError('IP: IDD: Field #'//TRIM(MinMaxString)//' default is invalid for Min/Max values, in class='//  &
        TRIM(ObjectDef(NumObjectDefs)%Name),EchoInputFile)
        ErrFlag=.true.
      ENDIF
    ENDIF
  ENDDO

  IF (ErrFlag) THEN
    CALL ShowContinueError('IP: Errors occured in ObjectDefinition for Class='//TRIM(ObjectDef2(NumObjectDefs2)%Name)// &
    ', Object not available for IDF processing.',EchoInputFile)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%AlphaorNumeric)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%NumRangeChks)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldChks)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%AlphFieldDefs)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%ReqField)
    DEALLOCATE(ObjectDef2(NumObjectDefs2)%AlphRetainCase)
    NumObjectDefs=NumObjectDefs-1
    ErrorsFound=.true.
  ENDIF

  RETURN

END SUBROUTINE AddObjectDefandParse2

SUBROUTINE GetObjectItem2(Object,Number,Alphas,NumAlphas,Numbers,NumNumbers,Status,NumBlank,AlphaBlank,   &
  AlphaFieldNames,NumericFieldNames) !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine gets the 'number' 'object' from the IDFRecord data structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na


  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: Object
  INTEGER, INTENT(IN) :: Number
  CHARACTER(len=*), INTENT(OUT), DIMENSION(:) :: Alphas
  INTEGER, INTENT(OUT) :: NumAlphas
  REAL(r64), INTENT(OUT), DIMENSION(:) :: Numbers
  !REAL, INTENT(OUT), DIMENSION(:) :: Numbers
  INTEGER, INTENT(OUT) :: NumNumbers
  INTEGER, INTENT(OUT) :: Status
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaBlank
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumBlank
  CHARACTER(len=*), DIMENSION(:), OPTIONAL :: AlphaFieldNames
  CHARACTER(len=*), DIMENSION(:), OPTIONAL :: NumericFieldNames

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER Count
  INTEGER LoopIndex
  CHARACTER(len=MaxObjectNameLength) ObjectWord
  CHARACTER(len=MaxObjectNameLength) UCObject
  CHARACTER(len=MaxObjectNameLength), SAVE, ALLOCATABLE, DIMENSION(:) :: AlphaArgs
  REAL, SAVE, ALLOCATABLE, DIMENSION(:) :: NumberArgs
  LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: AlphaArgsBlank
  LOGICAL, SAVE, ALLOCATABLE, DIMENSION(:) :: NumberArgsBlank
  INTEGER MaxAlphas,MaxNumbers
  INTEGER Found
  INTEGER StartRecord
  CHARACTER(len=32) :: cfld1
  CHARACTER(len=32) :: cfld2

  MaxAlphas=SIZE(Alphas,1)
  MaxNumbers=SIZE(Numbers,1)

  IF (.not. ALLOCATED(AlphaArgs)) THEN
    IF (NumObjectDefs2 == 0) THEN
      CALL ProcessInput
    ENDIF
    ALLOCATE(AlphaArgs(MaxAlphaArgsFound))
    ALLOCATE(NumberArgs(MaxNumericArgsFound))
    ALLOCATE(NumberArgsBlank(MaxNumericArgsFound))
    ALLOCATE(AlphaArgsBlank(MaxAlphaArgsFound))
  ENDIF

  Alphas(1:MaxAlphas)=' '
  Numbers(1:MaxNumbers)=0.0
  Count=0
  Status=-1
  UCOBject=MakeUPPERCase(Object)
  Found=FindIteminList(UCOBject,ListofObjects2,NumObjectDefs2)
  !IF (SortedIDD) THEN
  !  Found=FindIteminSortedList(UCOBject,ListofObjects2,NumObjectDefs2)
  !  IF (Found /= 0) Found=iListofObjects(Found)
  !ELSE
  !  Found=FindIteminList(UCOBject,ListofObjects2,NumObjectDefs2)
  !ENDIF
  IF (Found == 0) THEN   !  This is more of a developer problem
    CALL ShowFatalError('Requested object='//TRIM(UCObject)//', not found in Object Definitions -- incorrect IDD attached.')
  ENDIF

  IF (ObjectDef2(Found)%NumAlpha > 0) THEN
    IF (ObjectDef2(Found)%NumAlpha > MaxAlphas) THEN
      cfld1=IPTrimSigDigits(ObjectDef2(Found)%NumAlpha)
      cfld2=IPTrimSigDigits(MaxAlphas)
      CALL ShowFatalError('GetObjectItem: '//TRIM(Object)//', Number of ObjectDef Alpha Args ['//TRIM(cfld1)//  &
      '] > Size of AlphaArg array ['//TRIM(cfld2)//'].')
    ENDIF
    Alphas(1:ObjectDef2(Found)%NumAlpha)=Blank
  ENDIF
  IF (ObjectDef2(Found)%NumNumeric > 0) THEN
    IF (ObjectDef2(Found)%NumNumeric > MaxNumbers) THEN
      cfld1=IPTrimSigDigits(ObjectDef2(Found)%NumNumeric)
      cfld2=IPTrimSigDigits(MaxNumbers)
      CALL ShowFatalError('GetObjectItem: '//TRIM(Object)//', Number of ObjectDef Numeric Args ['//TRIM(cfld1)//  &
      '] > Size of NumericArg array ['//TRIM(cfld2)//'].')
    ENDIF
    Numbers(1:ObjectDef2(Found)%NumNumeric)=0.0
  ENDIF

  StartRecord=ObjectStartRecord(Found)
  IF (StartRecord == 0) THEN
    CALL ShowWarningError('Requested object='//TRIM(UCObject)//', not found in IDF.')
    Status=-1
    StartRecord=NumIDFRecords+1
  ENDIF

  IF (ObjectGotCount(Found) == 0) THEN
    WRITE(EchoInputFile,*) 'Getting object=',TRIM(UCObject)
  ENDIF
  ObjectGotCount(Found)=ObjectGotCount(Found)+1

  DO LoopIndex=StartRecord,NumIDFRecords2
    IF (IDFRecords2(LoopIndex)%Name == UCObject) THEN
      Count=Count+1
      IF (Count == Number) THEN
        ! Read this one
        CALL GetObjectItemfromFile2(LoopIndex,ObjectWord,AlphaArgs,NumAlphas,NumberArgs,NumNumbers)
        IF (NumAlphas > MaxAlphas .or. NumNumbers > MaxNumbers) THEN
          CALL ShowWarningError('Too many actual arguments for those expected on Object: '//TRIM(ObjectWord)//     &
          ' (GetObjectItem)',EchoInputFile)
        ENDIF
        NumAlphas=MIN(MaxAlphas,NumAlphas)
        NumNumbers=MIN(MaxNumbers,NumNumbers)
        IF (NumAlphas > 0) THEN
          Alphas(1:NumAlphas)=AlphaArgs(1:NumAlphas)
        ENDIF
        IF (NumNumbers > 0) THEN
          Numbers(1:NumNumbers)=NumberArgs(1:NumNumbers)
        ENDIF
        Status=1
        EXIT
      ENDIF
      !IF (Count == Number) THEN
      !  IDFRecordsGotten(LoopIndex)=.true.  ! only object level "gets" recorded
      !  ! Read this one
      !  CALL GetObjectItemfromFile(LoopIndex,ObjectWord,AlphaArgs,NumAlphas,NumberArgs,NumNumbers,AlphaArgsBlank,NumberArgsBlank)
      !  IF (NumAlphas > MaxAlphas .or. NumNumbers > MaxNumbers) THEN
      !    CALL ShowFatalError('Too many actual arguments for those expected on Object: '//TRIM(ObjectWord)//     &
      !                           ' (GetObjectItem)',EchoInputFile)
      !  ENDIF
      !  NumAlphas=MIN(MaxAlphas,NumAlphas)
      !  NumNumbers=MIN(MaxNumbers,NumNumbers)
      !  IF (NumAlphas > 0) THEN
      !    Alphas(1:NumAlphas)=AlphaArgs(1:NumAlphas)
      !  ENDIF
      !  IF (NumNumbers > 0) THEN
      !    Numbers(1:NumNumbers)=NumberArgs(1:NumNumbers)
      !  ENDIF
      !  IF (PRESENT(NumBlank)) THEN
      !    NumBlank=.true.
      !    IF (NumNumbers > 0) &
      !      NumBlank(1:NumNumbers)=NumberArgsBlank(1:NumNumbers)
      !  ENDIF
      !  IF (PRESENT(AlphaBlank)) THEN
      !    AlphaBlank=.true.
      !    IF (NumAlphas > 0) &
      !      AlphaBlank(1:NumAlphas)=AlphaArgsBlank(1:NumAlphas)
      !  ENDIF
      !  IF (PRESENT(AlphaFieldNames)) THEN
      !    AlphaFieldNames(1:ObjectDef2(Found)%NumAlpha)=ObjectDef2(Found)%AlphFieldChks(1:ObjectDef2(Found)%NumAlpha)
      !  ENDIF
      !  IF (PRESENT(NumericFieldNames)) THEN
      !    NumericFieldNames(1:ObjectDef2(Found)%NumNumeric)=ObjectDef2(Found)%NumRangeChks(1:ObjectDef2(Found)%NumNumeric)%FieldName
      !  ENDIF
      !  Status=1
      !  EXIT
      !ENDIF
    ENDIF
  END DO

  RETURN

END SUBROUTINE GetObjectItem2   !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

SUBROUTINE GetObjectItemfromFile2(Which,ObjectWord,AlphaArgs,NumAlpha,NumericArgs,NumNumeric,AlphaBlanks,NumericBlanks)    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine "gets" the object instance from the data structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: Which
  CHARACTER(len=*), INTENT(OUT) :: ObjectWord
  CHARACTER(len=*), INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaArgs
  INTEGER, INTENT(OUT) :: NumAlpha
  REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumericArgs
  INTEGER, INTENT(OUT) :: NumNumeric
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: AlphaBlanks
  LOGICAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: NumericBlanks

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  TYPE (LineDefinition):: xLineItem                        ! Description of current record

  IF (Which > 0 .and. Which <= NumIDFRecords) THEN
    xLineItem=IDFRecords2(Which)
    ObjectWord=xLineItem%Name
    NumAlpha=xLineItem%NumAlphas
    NumNumeric=xLineItem%NumNumbers
    IF (PRESENT(AlphaArgs)) THEN
      IF (NumAlpha >=1) THEN
        AlphaArgs(1:NumAlpha)=xLineItem%Alphas(1:NumAlpha)
      ENDIF
    ENDIF
    IF (PRESENT(AlphaBlanks)) THEN
      IF (NumAlpha >=1) THEN
        AlphaBlanks(1:NumAlpha)=xLineItem%AlphBlank(1:NumAlpha)
      ENDIF
    ENDIF
    IF (PRESENT(NumericArgs)) THEN
      IF (NumNumeric >= 1) THEN
        NumericArgs(1:NumNumeric)=xLineItem%Numbers(1:NumNumeric)
      ENDIF
    ENDIF
    IF (PRESENT(NumericBlanks)) THEN
      IF (NumNumeric >= 1) THEN
        NumericBlanks(1:NumNumeric)=xLineItem%NumBlank(1:NumNumeric)
      ENDIF
    ENDIF
  ELSE
    WRITE(EchoInputFile,*) ' Requested Record',Which,' not in range, 1 -- ',NumIDFRecords2
  ENDIF

  RETURN

END SUBROUTINE GetObjectItemfromFile2    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

SUBROUTINE ProcessDataDicFile2(ErrorsFound) !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   August 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine processes data dictionary file for EnergyPlus.
  ! The structure of the sections and objects are stored in derived
  ! types (SectionDefs and ObjectDefs)

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  LOGICAL, INTENT(INOUT) :: ErrorsFound ! set to true if any errors flagged during IDD processing

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL  :: EndofFile = .false.        ! True when End of File has been reached (IDD or IDF)
  INTEGER Pos                            ! Test of scanning position on the current input line
  TYPE (SectionsDefinition), ALLOCATABLE :: TempSectionDef2(:)  ! Like SectionDef, used during Re-allocation
  TYPE (ObjectsDefinition), ALLOCATABLE :: TempObjectDef2(:)    ! Like ObjectDef, used during Re-allocation
  LOGICAL BlankLine

  MaxSectionDefs=SectionDefAllocInc
  MaxObjectDefs=ObjectDefAllocInc

  IF (ALLOCATED(SectionDef2)) THEN !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(SectionDef2)
  END IF

  IF (ALLOCATED(ObjectDef2)) THEN  !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(ObjectDef2)
  END IF

  ALLOCATE (SectionDef2(MaxSectionDefs))
  SectionDef2%Name=' '        ! Name of the section
  SectionDef2%NumFound=0      ! Number of this section found in IDF

  ALLOCATE(ObjectDef2(MaxObjectDefs))
  ObjectDef2%Name=' '                ! Name of the object
  ObjectDef2%NumParams=0             ! Number of parameters to be processed for each object
  ObjectDef2%NumAlpha=0              ! Number of Alpha elements in the object
  ObjectDef2%NumNumeric=0            ! Number of Numeric elements in the object
  ObjectDef2%MinNumFields=0          ! Minimum number of fields
  ObjectDef2%NameAlpha1=.false.      ! by default, not the "name"
  ObjectDef2%ObsPtr=0.               ! by default, not obsolete
  ObjectDef2%NumFound=0              ! Number of this object found in IDF
  ObjectDef2%UniqueObject=.false.    ! by default, not Unique
  ObjectDef2%RequiredObject=.false.  ! by default, not Required

  NumObjectDefs2=0
  NumSectionDefs=0
  EndofFile=.false.

  DO WHILE (.not. EndofFile)
    CALL ReadInputLine(IDDFile,Pos,BlankLine,InputLineLength,EndofFile)
    IF (BlankLine .or. EndofFile) THEN
      CYCLE
    END IF
    Pos=SCAN(InputLine(1:InputLineLength),',;')
    If (Pos /= 0) then

      If (InputLine(Pos:Pos) == ';') then
        CALL AddSectionDef(InputLine(1:Pos-1),ErrorsFound)
        IF (NumSectionDefs == MaxSectionDefs) THEN
          ALLOCATE (TempSectionDef2(MaxSectionDefs+SectionDefAllocInc))
          TempSectionDef2%Name=' '
          TempSectionDef2%NumFound=0
          TempSectionDef2(1:MaxSectionDefs)=SectionDef2
          DEALLOCATE (SectionDef2)
          ALLOCATE (SectionDef2(MaxSectionDefs+SectionDefAllocInc))
          SectionDef=TempSectionDef2
          DEALLOCATE (TempSectionDef2)
          MaxSectionDefs=MaxSectionDefs+SectionDefAllocInc
        ENDIF
      else
        CALL AddObjectDefandParse2(InputLine(1:Pos-1),Pos,EndofFile,ErrorsFound)
        IF (NumObjectDefs2 == MaxObjectDefs) THEN
          ALLOCATE (TempObjectDef2(MaxObjectDefs+ObjectDefAllocInc))
          TempObjectDef2%Name=' '         ! Name of the object
          TempObjectDef2%NumParams=0      ! Number of parameters to be processed for each object
          TempObjectDef2%NumAlpha=0       ! Number of Alpha elements in the object
          TempObjectDef2%NumNumeric=0     ! Number of Numeric elements in the object
          TempObjectDef2%MinNumFields=0   ! Minimum number of fields
          TempObjectDef2%NameAlpha1=.false.  ! by default, not the "name"
          TempObjectDef2%ObsPtr=0.        ! by default, not obsolete
          TempObjectDef2%NumFound=0       ! Number of this object found in IDF
          TempObjectDef2%UniqueObject=.false.    ! by default, not Unique
          TempObjectDef2%RequiredObject=.false.  ! by default, not Required
          TempObjectDef2(1:MaxObjectDefs)=ObjectDef2
          DEALLOCATE (ObjectDef2)
          ALLOCATE (ObjectDef2(MaxObjectDefs+ObjectDefAllocInc))
          ObjectDef2=TempObjectDef2
          DEALLOCATE (TempObjectDef2)
          MaxObjectDefs=MaxObjectDefs+ObjectDefAllocInc
        ENDIF
      endif

    else
      CALL ShowSevereError(', or ; expected on this line',EchoInputFile)
      ErrorsFound=.true.
    endif

  END DO

  RETURN

END SUBROUTINE ProcessDataDicFile2

SUBROUTINE ProcessInputDataFile2    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   August 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine processes data dictionary file for EnergyPlus.
  ! The structure of the sections and objects are stored in derived
  ! types (SectionDefs and ObjectDefs)

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  TYPE (FileSectionsDefinition), ALLOCATABLE :: TempSectionsonFile(:)   ! Used during reallocation procedure
  TYPE (LineDefinition), ALLOCATABLE :: TempIDFRecords(:)   ! Used during reallocation procedure

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  LOGICAL :: EndofFile = .false.
  LOGICAL BlankLine
  INTEGER Pos
  CHARACTER(len=25) LineNum

  MaxIDFRecords=ObjectsIDFAllocInc
  NumIDFRecords2=0
  MaxIDFSections=SectionsIDFAllocInc
  NumIDFSections=0

  IF (ALLOCATED(SectionsonFile2)) THEN !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(SectionsonFile2)
  END IF

  IF (ALLOCATED(IDFRecords2)) THEN  !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(IDFRecords2)
  END IF

  ALLOCATE (SectionsonFile2(MaxIDFSections))
  SectionsonFile2%Name=' '        ! Name of this section
  SectionsonFile2%FirstRecord=0   ! Record number of first object in section
  SectionsonFile2%LastRecord=0    ! Record number of last object in section
  ALLOCATE (IDFRecords2(MaxIDFRecords))
  IDFRecords2%Name=' '          ! Object name for this record
  IDFRecords2%NumAlphas=0       ! Number of alphas on this record
  IDFRecords2%NumNumbers=0      ! Number of numbers on this record

  IF (ALLOCATED(LineItem2%Numbers)) THEN !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(LineItem2%Numbers)
  END IF

  IF (ALLOCATED(LineItem2%NumBlank)) THEN  !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(LineItem2%NumBlank)
  END IF

  IF (ALLOCATED(LineItem2%Alphas)) THEN !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(LineItem2%Alphas)
  END IF

  IF (ALLOCATED(LineItem2%AlphBlank)) THEN  !RS: Debugging: To handle multiple calls of this routine (10/7/14)
    DEALLOCATE(LineItem2%AlphBlank)
  END IF

  ALLOCATE (LineItem2%Numbers(MaxNumericArgsFound))
  ALLOCATE (LineItem2%NumBlank(MaxNumericArgsFound))
  ALLOCATE (LineItem2%Alphas(MaxAlphaArgsFound))
  ALLOCATE (LineItem2%AlphBlank(MaxAlphaArgsFound))
  EndofFile=.false.

  DO WHILE (.not. EndofFile)
    CALL ReadInputLine(IDFFile,Pos,BlankLine,InputLineLength,EndofFile)
    IF (BlankLine .or. EndofFile) THEN
      CYCLE
    END IF
    Pos=SCAN(InputLine,',;')
    If (Pos /= 0) then
      If (InputLine(Pos:Pos) == ';') then
        CALL ValidateSection(InputLine(1:Pos-1),NumLines)
        IF (NumIDFSections == MaxIDFSections) THEN
          ALLOCATE (TempSectionsonFile(MaxIDFSections+SectionsIDFAllocInc))
          TempSectionsonFile%Name=' '        ! Name of this section
          TempSectionsonFile%FirstRecord=0   ! Record number of first object in section
          TempSectionsonFile%LastRecord=0    ! Record number of last object in section
          TempSectionsonFile(1:MaxIDFSections)=SectionsonFile
          DEALLOCATE (SectionsonFile)
          ALLOCATE (SectionsonFile(MaxIDFSections+SectionsIDFAllocInc))
          SectionsonFile2=TempSectionsonFile
          DEALLOCATE (TempSectionsonFile)
          MaxIDFSections=MaxIDFSections+SectionsIDFAllocInc
        ENDIF
      else
        CALL ValidateObjectandParse2(InputLine(1:Pos-1),Pos,EndofFile)
        IF (NumIDFRecords2 == MaxIDFRecords) THEN
          ALLOCATE(TempIDFRecords(MaxIDFRecords+ObjectsIDFAllocInc))
          TempIDFRecords%Name=' '          ! Object name for this record
          TempIDFRecords%NumAlphas=0       ! Number of alphas on this record
          TempIDFRecords%NumNumbers=0      ! Number of numbers on this record
          TempIDFRecords(1:MaxIDFRecords)=IDFRecords
          DEALLOCATE(IDFRecords)
          ALLOCATE(IDFRecords(MaxIDFRecords+ObjectsIDFAllocInc))
          IDFRecords2=TempIDFRecords
          DEALLOCATE(TempIDFRecords)
          MaxIDFRecords=MaxIDFRecords+ObjectsIDFAllocInc
        ENDIF
      endif
    else
      !Error condition, no , or ; on first line
      WRITE(LineNum,*) NumLines
      LineNum=ADJUSTL(LineNum)
      CALL ShowMessage('IDF Line='//TRIM(LineNum)//' '//TRIM(InputLine))
      CALL ShowSevereError(', or ; expected on this line',EchoInputFile)
    endif

  END DO

  IF (NumIDFSections > 0) THEN
    SectionsonFile2(NumIDFSections)%LastRecord=NumIDFRecords
  ENDIF

  IF (OverallErrorFlag) THEN
    CALL ShowSevereError('Possible incorrect IDD File')
    CALL ShowContinueError('Possible Invalid Numerics or other problems')
    CALL ShowFatalError('Errors occurred on processing IDF file. Preceding condition(s) cause termination.')
  ENDIF

  IF (NumIDFRecords2 > 0) THEN
    DO Pos=1,NumObjectDefs2
      IF (ObjectDef(Pos)%RequiredObject .and. ObjectDef(Pos)%NumFound == 0) THEN
        CALL ShowSevereError('No items found for Required Object='//TRIM(ObjectDef(Pos)%Name))
        NumMiscErrorsFound=NumMiscErrorsFound+1
      ENDIF
    ENDDO
  ENDIF

  RETURN

END SUBROUTINE ProcessInputDataFile2    !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

SUBROUTINE ValidateObjectandParse2(ProposedObject,CurPos,EndofFile) !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Linda K. Lawrie
  !       DATE WRITTEN   September 1997
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine validates the proposed object from the IDF and then
  ! parses it, putting it into the internal InputProcessor Data structure.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! na

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN) :: ProposedObject
  INTEGER, INTENT(INOUT) :: CurPos
  LOGICAL, INTENT(INOUT) :: EndofFile

  ! SUBROUTINE PARAMETER DEFINITIONS:
  INTEGER, PARAMETER :: dimLineBuf=10

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  CHARACTER(len=MaxObjectNameLength) SqueezedObject
  CHARACTER(len=MaxAlphaArgLength) SqueezedArg
  INTEGER, SAVE :: Found
  INTEGER NumArg
  INTEGER NumArgExpected
  INTEGER NumAlpha
  INTEGER NumNumeric
  INTEGER Pos
  LOGICAL EndofObject
  LOGICAL BlankLine
  LOGICAL,SAVE  :: ErrFlag=.false.
  INTEGER LenLeft
  INTEGER Count
  CHARACTER(len=32) FieldString
  CHARACTER(len=MaxFieldNameLength) FieldNameString
  CHARACTER(len=300) Message
  CHARACTER(len=300) cStartLine
  CHARACTER(len=300) cStartName
  CHARACTER(len=300), DIMENSION(dimLineBuf), SAVE :: LineBuf
  INTEGER, SAVE :: StartLine
  INTEGER, SAVE :: NumConxLines
  INTEGER, SAVE :: CurLines
  INTEGER, SAVE :: CurQPtr

  CHARACTER(len=52) :: String
  LOGICAL IDidntMeanIt
  LOGICAL TestingObject
  LOGICAL TransitionDefer
  INTEGER TFound
  INTEGER, EXTERNAL :: FindNonSpace
  INTEGER NextChr
  CHARACTER(len=32) :: String1

  INTEGER :: DebugFile       =150 !RS: Debugging file denotion, hopfully this works.

  !OPEN(unit=DebugFile,file='Debug.txt')    !RS: Debugging

  SqueezedObject=MakeUPPERCase(ADJUSTL(ProposedObject))
  IF (LEN_TRIM(ADJUSTL(ProposedObject)) > MaxObjectNameLength) THEN
    CALL ShowWarningError('IP: Object name length exceeds maximum, will be truncated='//TRIM(ProposedObject),EchoInputFile)
    CALL ShowContinueError('Will be processed as Object='//TRIM(SqueezedObject),EchoInputFile)
  ENDIF
  IDidntMeanIt=.false.

  TestingObject=.true.
  TransitionDefer=.false.
  DO WHILE (TestingObject)
    ErrFlag=.false.
    IDidntMeanIt=.false.
    Found=FindIteminList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    IF (Found /= 0 .and. ObjectDef2(Found)%ObsPtr > 0) THEN
      TFound=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
      IF (TFound /= 0) THEN
        Found=0    ! being handled differently for this obsolete object
      END IF
    ENDIF
    !ErrFlag=.false.
    !IDidntMeanIt=.false.
    !IF (SortedIDD) THEN
    !  Found=FindIteminSortedList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !  IF (Found /= 0) Found=iListofObjects(Found)
    !ELSE
    !  Found=FindIteminList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !ENDIF
    !IF (Found /= 0) THEN
    !  IF (ObjectDef(Found)%ObsPtr > 0) THEN
    !    TFound=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
    !    IF (TFound /= 0) THEN
    !      IF (RepObjects(TFound)%Transitioned) THEN
    !        IF (.not. RepObjects(TFound)%Used)  &
    !           CALL ShowWarningError('IP: Objects="'//TRIM(ADJUSTL(ProposedObject))//  &
    !              '" are being transitioned to this object="'//  &
    !              TRIM(RepObjects(TFound)%NewName)//'"')
    !        RepObjects(TFound)%Used=.true.
    !        IF (SortedIDD) THEN
    !          Found=FindIteminSortedList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !          IF (Found /= 0) Found=iListofObjects(Found)
    !        ELSE
    !          Found=FindIteminList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !        ENDIF
    !      ELSEIF (RepObjects(TFound)%TransitionDefer) THEN
    !        IF (.not. RepObjects(TFound)%Used)  &
    !           CALL ShowWarningError('IP: Objects="'//TRIM(ADJUSTL(ProposedObject))//  &
    !              '" are being transitioned to this object="'//  &
    !              TRIM(RepObjects(TFound)%NewName)//'"')
    !        RepObjects(TFound)%Used=.true.
    !        IF (SortedIDD) THEN
    !          Found=FindIteminSortedList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !          IF (Found /= 0) Found=iListofObjects(Found)
    !        ELSE
    !          Found=FindIteminList(SqueezedObject,ListofObjects2,NumObjectDefs2)
    !        ENDIF
    !        TransitionDefer=.true.
    !      ELSE
    !        Found=0    ! being handled differently for this obsolete object
    !      ENDIF
    !    ENDIF
    !  ENDIF
    !ENDIF

    TestingObject=.false.
    IF (Found == 0) THEN
      ! Check to see if it's a "secret" object
      Found=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
      IF (Found == 0) THEN
        !CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
        !   ' Did not find "'//TRIM(ADJUSTL(ProposedObject))//'" in list of Objects',EchoInputFile) !RS: Secret Search String
        IF(DebugFile .EQ. 9 .OR. DebugFile .EQ. 10) THEN
          WRITE(*,*) 'Error with OutputFileDebug'    !RS: Debugging: Searching for a mis-set file number
        END IF
        WRITE(DebugFile,*) 'IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//' Did not find "'&
        //TRIM(ADJUSTL(ProposedObject))//'" in list of Objects'
        ! Will need to parse to next ;
        ErrFlag=.true.
      ELSEIF (RepObjects(Found)%Deleted) THEN
        IF (.not. RepObjects(Found)%Used) THEN
          CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" have been deleted from the IDD.  Will be ignored.')
          RepObjects(Found)%Used=.true.
        ENDIF
        IDidntMeanIt=.true.
        ErrFlag=.true.
        Found=0
      ELSEIF (RepObjects(Found)%TransitionDefer) THEN

      ELSE ! This name is replaced with something else
        IF (.not. RepObjects(Found)%Used) THEN
          IF (.not. RepObjects(Found)%Transitioned) THEN
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being replaced with this object="'//  &
            TRIM(RepObjects(Found)%NewName)//'"')
            RepObjects(Found)%Used=.true.
            SqueezedObject=RepObjects(Found)%NewName
            TestingObject=.true.
          ELSE
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being transitioned to this object="'//  &
            TRIM(RepObjects(Found)%NewName)//'"')
            RepObjects(Found)%Used=.true.
            IF (SortedIDD) THEN
              Found=FindIteminSortedList(SqueezedObject,ListofObjects,NumObjectDefs)
              IF (Found /= 0) Found=iListofObjects(Found)
            ELSE
              Found=FindIteminList(SqueezedObject,ListofObjects,NumObjectDefs)
            ENDIF
          ENDIF
        ELSEIF (.not. RepObjects(Found)%Transitioned) THEN
          SqueezedObject=RepObjects(Found)%NewName
          TestingObject=.true.
        ELSE
          IF (SortedIDD) THEN
            Found=FindIteminSortedList(SqueezedObject,ListofObjects2,NumObjectDefs2)
            IF (Found /= 0) Found=iListofObjects(Found)
          ELSE
            Found=FindIteminList(SqueezedObject,ListofObjects2,NumObjectDefs2)
          ENDIF
        ENDIF
      ENDIF
    ELSE

      ! Start Parsing the Object according to definition

      ErrFlag=.false.
      LineItem2%Name=SqueezedObject
      LineItem2%Alphas=Blank
      LineItem2%AlphBlank=.false.
      LineItem2%NumAlphas=0
      LineItem2%Numbers=0.0
      LineItem2%NumNumbers=0
      LineItem2%NumBlank=.false.
      LineItem2%ObjectDefPtr=Found
      NumArgExpected=ObjectDef2(Found)%NumParams
      ObjectDef2(Found)%NumFound=ObjectDef2(Found)%NumFound+1
      IF (ObjectDef2(Found)%UniqueObject .and. ObjectDef2(Found)%NumFound > 1) THEN
        CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
        ' Multiple occurrences of Unique Object='//TRIM(ADJUSTL(ProposedObject)))
        NumMiscErrorsFound=NumMiscErrorsFound+1
      ENDIF
      IF (ObjectDef2(Found)%ObsPtr > 0) THEN
        TFound=FindItemInList(SqueezedObject,RepObjects%OldName,NumSecretObjects)
        IF (TFound == 0) THEN
          CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' Obsolete object='//TRIM(ADJUSTL(ProposedObject))//  &
          ', encountered.  Should be replaced with new object='//  &
          TRIM(ObsoleteObjectsRepNames(ObjectDef(Found)%ObsPtr)))
        ELSEIF (.not. RepObjects(TFound)%Used .and. RepObjects(TFound)%Transitioned) THEN
          CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' Objects="'//TRIM(ADJUSTL(ProposedObject))//'" are being transitioned to this object="'//  &
          TRIM(RepObjects(TFound)%NewName)//'"')
          RepObjects(TFound)%Used=.true.
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  NumArg=0
  NumAlpha=0
  NumNumeric=0
  EndofObject=.false.
  CurPos=CurPos+1

  !  Keep context buffer in case of errors
  LineBuf=Blank
  NumConxLines=0
  StartLine=NumLines
  cStartLine=InputLine(1:300)
  cStartName=Blank
  NumConxLines=0
  CurLines=NumLines
  CurQPtr=0

  DO WHILE (.not. EndofFile .and. .not. EndofObject)
    IF (CurLines /= NumLines) THEN
      NumConxLines=MIN(NumConxLines+1,dimLineBuf)
      CurQPtr=CurQPtr+1
      IF (CurQPtr == 1 .and. cStartName == Blank .and. InputLine /= Blank) THEN
        IF (Found > 0) THEN
          IF (ObjectDef(Found)%NameAlpha1) THEN
            Pos=INDEX(InputLine,',')
            cStartName=InputLine(1:Pos-1)
            cStartName=ADJUSTL(cStartName)
          ENDIF
        ENDIF
      ENDIF
      IF (CurQPtr > dimLineBuf) CurQPtr=1
      LineBuf(CurQPtr)=InputLine(1:300)
      CurLines=NumLines
    ENDIF
    IF (CurPos <= InputLineLength) THEN
      Pos=SCAN(InputLine(CurPos:InputLineLength),',;')
      IF (Pos == 0) THEN
        IF (InputLine(InputLineLength:InputLineLength) == '!') THEN
          LenLeft=LEN_TRIM(InputLine(CurPos:InputLineLength-1))
        ELSE
          LenLeft=LEN_TRIM(InputLine(CurPos:InputLineLength))
        ENDIF
        IF (LenLeft == 0) THEN
          CurPos=InputLineLength+1
          CYCLE
        ELSE
          IF (InputLine(InputLineLength:InputLineLength) == '!') THEN
            Pos=InputLineLength-CurPos+1
            CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Comma being inserted after:"'//InputLine(CurPos:InputLineLength-1)//   &
            '" in Object='//TRIM(SqueezedObject),EchoInputFile)
          ELSE
            Pos=InputLineLength-CurPos+2
            CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Comma being inserted after:"'//InputLine(CurPos:InputLineLength)// &
            '" in Object='//TRIM(SqueezedObject),EchoInputFile)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL ReadInputLine(IDFFile,CurPos,BlankLine,InputLineLength,EndofFile)
      CYCLE
    ENDIF
    IF (Pos > 0) THEN
      IF (.not. ErrFlag) THEN
        IF (CurPos <= CurPos+Pos-2) THEN
          SqueezedArg=MakeUPPERCase(ADJUSTL(InputLine(CurPos:CurPos+Pos-2)))
          IF (LEN_TRIM(ADJUSTL(InputLine(CurPos:CurPos+Pos-2))) > MaxAlphaArgLength) THEN
            CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
            CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Alpha Argument length exceeds maximum, will be truncated='// &
            TRIM(InputLine(CurPos:CurPos+Pos-2)), EchoInputFile)
            CALL ShowContinueError('Will be processed as Alpha='//TRIM(SqueezedArg),EchoInputFile)
          ENDIF
        ELSE
          SqueezedArg=Blank
        ENDIF
        IF (NumArg == NumArgExpected .and. .not. ObjectDef2(Found)%ExtensibleObject) THEN
          CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
          CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
          ' Error detected for Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
          CALL ShowContinueError(' Maximum arguments reached for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
          EchoInputFile)
          ErrFlag=.true.
        ELSE
          IF (NumArg == NumArgExpected .and. ObjectDef2(Found)%ExtensibleObject) THEN
            CALL ExtendObjectDefinition(Found,NumArgExpected)
          ENDIF
          NumArg=NumArg+1
          IF (ObjectDef2(Found)%AlphaorNumeric(NumArg)) THEN
            IF (NumAlpha == ObjectDef2(Found)%NumAlpha) THEN
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Error detected for Object='//TRIM(ObjectDef2(Found)%Name),EchoInputFile)
              CALL ShowContinueError(' Too many Alphas for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
              EchoInputFile)
              ErrFlag=.true.
            ELSE
              NumAlpha=NumAlpha+1
              LineItem2%NumAlphas=NumAlpha
              IF (ObjectDef2(Found)%AlphRetainCase(NumArg)) THEN
                SqueezedArg=InputLine(CurPos:CurPos+Pos-2)
                SqueezedArg=ADJUSTL(SqueezedArg)
              ENDIF
              IF (SqueezedArg /= Blank) THEN
                LineItem2%Alphas(NumAlpha)=SqueezedArg
              ELSEIF (ObjectDef2(Found)%ReqField(NumArg)) THEN  ! Blank Argument
                IF (ObjectDef2(Found)%AlphFieldDefs(NumAlpha) /= Blank) THEN
                  LineItem2%Alphas(NumAlpha)=ObjectDef2(Found)%AlphFieldDefs(NumAlpha)
                ELSE
                  IF (ObjectDef2(Found)%NameAlpha1 .and. NumAlpha /= 1) THEN
                    CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                    CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                    ' Error detected in Object='//TRIM(ObjectDef2(Found)%Name)//', name='//  &
                    TRIM(LineItem2%Alphas(1)),EchoInputFile)
                  ELSE
                    CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                    CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                    ' Error detected in Object='//TRIM(ObjectDef2(Found)%Name),EchoInputFile)
                  ENDIF
                  CALL ShowContinueError('Field ['//TRIM(ObjectDef2(Found)%AlphFieldChks(NumAlpha))//  &
                  '] is required but was blank',EchoInputFile)
                  NumBlankReqFieldFound=NumBlankReqFieldFound+1
                ENDIF
              ELSE
                LineItem2%AlphBlank(NumAlpha)=.true.
                IF (ObjectDef2(Found)%AlphFieldDefs(NumAlpha) /= Blank) THEN
                  LineItem2%Alphas(NumAlpha)=ObjectDef2(Found)%AlphFieldDefs(NumAlpha)
                ENDIF
              ENDIF
            ENDIF
          ELSE
            IF (NumNumeric == ObjectDef2(Found)%NumNumeric) THEN
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Error detected for Object='//TRIM(ObjectDef2(Found)%Name),EchoInputFile)
              CALL ShowContinueError(' Too many Numbers for this object, trying to process ->'//TRIM(SqueezedArg)//'<-',  &
              EchoInputFile)
              ErrFlag=.true.
            ELSE
              NumNumeric=NumNumeric+1
              LineItem2%NumNumbers=NumNumeric
              IF (SqueezedArg /= Blank) THEN
                IF (.not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoSizable .and.   &
                .not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoCalculatable) THEN
                LineItem2%Numbers(NumNumeric)=ProcessNumber(SqueezedArg,Errflag)
              ELSEIF (SqueezedArg == 'AUTOSIZE') THEN
                LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
              ELSEIF (SqueezedArg == 'AUTOCALCULATE') THEN
                LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
              ELSE
                LineItem2%Numbers(NumNumeric)=ProcessNumber(SqueezedArg,Errflag)
              ENDIF
            ELSE  ! numeric arg is blank.
              IF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefaultChk) THEN  ! blank arg has default
                IF (.not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoSize .and.   &
                .not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoCalculatable) THEN
                LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%Default
                LineItem2%NumBlank(NumNumeric)=.true.
              ELSEIF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoSize) THEN
                LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
                LineItem2%NumBlank(NumNumeric)=.true.
              ELSEIF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
                LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
                LineItem2%NumBlank(NumNumeric)=.true.
              ENDIF
              ErrFlag=.false.
            ELSE ! blank arg does not have default
              IF (ObjectDef2(Found)%ReqField(NumArg)) THEN  ! arg is required
                IF (ObjectDef2(Found)%NameAlpha1) THEN  ! object has name field - more context for error
                  CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                  CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                  ' Error detected in Object='//TRIM(ObjectDef(Found)%Name)// &
                  ', name='//TRIM(LineItem%Alphas(1)),EchoInputFile)
                  ErrFlag=.true.
                ELSE  ! object does not have name field
                  CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
                  CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
                  ' Error detected in Object='//TRIM(ObjectDef(Found)%Name),EchoInputFile)
                  ErrFlag=.true.
                ENDIF
                CALL ShowContinueError('Field ['//TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName)//  &
                '] is required but was blank',EchoInputFile)
                NumBlankReqFieldFound=NumBlankReqFieldFound+1
              ENDIF
              LineItem2%Numbers(NumNumeric)=0.0
              LineItem2%NumBlank(NumNumeric)=.true.
              !LineItem%Numbers(NumNumeric)=-999999.  !0.0
              !CALL ShowWarningError('Default number in Input, in object='//TRIM(ObjectDef(Found)%Name))
            ENDIF
          ENDIF
          IF (ErrFlag) THEN
            IF (SqueezedArg(1:1) /= '=') THEN  ! argument does not start with "=" (parametric)
              FieldString=IPTrimSigDigits(NumNumeric)
              FieldNameString=ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName
              IF (FieldNameString /= Blank) THEN
                Message='Invalid Number in Numeric Field#'//TRIM(FieldString)//' ('//TRIM(FieldNameString)//  &
                '), value='//TRIM(SqueezedArg)
              ELSE ! Field Name not recorded
                Message='Invalid Number in Numeric Field#'//TRIM(FieldString)//', value='//TRIM(SqueezedArg)
              ENDIF
              Message=TRIM(Message)//', in '//TRIM(ObjectDef2(Found)%Name)
              IF (ObjectDef2(Found)%NameAlpha1) THEN
                Message=TRIM(Message)//'='//TRIM(LineItem2%Alphas(1))
              ENDIF
              CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' '//TRIM(Message),EchoInputFile)
            ELSE  ! parametric in Numeric field
              ErrFlag=.false.
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF (InputLine(CurPos+Pos-1:CurPos+Pos-1) == ';') THEN
    EndofObject=.true.
    ! Check if more characters on line -- and if first is a comment character
    IF (InputLine(CurPos+Pos:) /= Blank) THEN
      NextChr=FindNonSpace(InputLine(CurPos+Pos:))
      IF (InputLine(CurPos+Pos+NextChr-1:CurPos+Pos+NextChr-1) /= '!') THEN
        CALL DumpCurrentLineBuffer(StartLine,cStartLine,cStartName,NumLines,NumConxLines,LineBuf,CurQPtr)
        CALL ShowWarningError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
        ' End of Object="'//TRIM(ObjectDef(Found)%Name)//  &
        '" reached, but next part of line not comment.',EchoInputFile)
        CALL ShowContinueError('Final line above shows line that contains error.')
      ENDIF
    ENDIF
  ENDIF
  CurPos=CurPos+Pos
ENDIF

END DO

! Store to IDFRecord Data Structure, ErrFlag is true if there was an error
! Check out MinimumNumberOfFields
IF (.not. ErrFlag .and. .not. IDidntMeanIt) THEN
  IF (NumArg < ObjectDef2(Found)%MinNumFields) THEN
    IF (ObjectDef2(Found)%NameAlpha1) THEN
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
      ', name='//TRIM(LineItem%Alphas(1))//       &
      ', entered with less than minimum number of fields.')
    ELSE
      CALL ShowAuditErrorMessage(' ** Warning ** ','IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
      ', entered with less than minimum number of fields.')
    ENDIF
    CALL ShowAuditErrorMessage(' **   ~~~   ** ','Attempting fill to minimum.')
    NumAlpha=0
    NumNumeric=0
    IF (ObjectDef2(Found)%MinNumFields > ObjectDef2(Found)%NumParams) THEN
      String=IPTrimSigDigits(ObjectDef2(Found)%MinNumFields)
      String1=IPTrimSigDigits(ObjectDef2(Found)%NumParams)
      CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
      ' Object \min-fields > number of fields specified, Object='//TRIM(ObjectDef(Found)%Name))
      CALL ShowContinueError('..\min-fields='//TRIM(String)//  &
      ', total number of fields in object definition='//TRIM(String1))
      ErrFlag=.true.
    ELSE
      DO Count=1,ObjectDef2(Found)%MinNumFields
        IF (ObjectDef2(Found)%AlphaOrNumeric(Count)) THEN
          NumAlpha=NumAlpha+1
          IF (NumAlpha <= LineItem2%NumAlphas) CYCLE
          LineItem2%NumAlphas=LineItem2%NumAlphas+1
          IF (ObjectDef2(Found)%AlphFieldDefs(LineItem2%NumAlphas) /= Blank) THEN
            LineItem2%Alphas(LineItem2%NumAlphas)=ObjectDef2(Found)%AlphFieldDefs(LineItem2%NumAlphas)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(ObjectDef2(Found)%AlphFieldDefs(LineItem2%NumAlphas))//   &
            '   ! field=>'//TRIM(ObjectDef2(Found)%AlphFieldChks(NumAlpha)))
          ELSEIF (ObjectDef2(Found)%ReqField(Count)) THEN
            IF (ObjectDef2(Found)%NameAlpha1) THEN
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
              ', name='//TRIM(LineItem2%Alphas(1))// &
              ', Required Field=['//  &
              TRIM(ObjectDef2(Found)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ELSE
              CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
              ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
              ', Required Field=['//  &
              TRIM(ObjectDef2(Found)%AlphFieldChks(NumAlpha))//   &
              '] was blank.',EchoInputFile)
            ENDIF
            ErrFlag=.true.
          ELSE
            LineItem2%Alphas(LineItem2%NumAlphas)=Blank
            LineItem2%AlphBlank(LineItem2%NumAlphas)=.true.
            CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
            TRIM(ObjectDef2(Found)%AlphFieldChks(NumAlpha)))
          ENDIF
        ELSE
          NumNumeric=NumNumeric+1
          IF (NumNumeric <= LineItem2%NumNumbers) CYCLE
          LineItem2%NumNumbers=LineItem2%NumNumbers+1
          LineItem2%NumBlank(NumNumeric)=.true.
          IF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%Defaultchk) THEN
            IF (.not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoSize .and.   &
            .not. ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%Default
            WRITE(String,*) ObjectDef2(Found)%NumRangeChks(NumNumeric)%Default
            String=ADJUSTL(String)
            CALL ShowAuditErrorMessage(' **   Add   ** ',TRIM(String)//  &
            '   ! field=>'//TRIM(ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoSize) THEN
            LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoSizeValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autosize    ! field=>'//  &
            TRIM(ObjectDef(Found)%NumRangeChks(NumNumeric)%FieldName))
          ELSEIF (ObjectDef2(Found)%NumRangeChks(NumNumeric)%DefAutoCalculate) THEN
            LineItem2%Numbers(NumNumeric)=ObjectDef2(Found)%NumRangeChks(NumNumeric)%AutoCalculateValue
            CALL ShowAuditErrorMessage(' **   Add   ** ','autocalculate    ! field=>'//  &
            TRIM(ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName))
          ENDIF
        ELSEIF (ObjectDef2(Found)%ReqField(Count)) THEN
          IF (ObjectDef2(Found)%NameAlpha1) THEN
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
            ', name='//TRIM(LineItem2%Alphas(1))// &
            ', Required Field=['//  &
            TRIM(ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ELSE
            CALL ShowSevereError('IP: IDF line~'//TRIM(IPTrimSigDigits(NumLines))//  &
            ' Object='//TRIM(ObjectDef2(Found)%Name)//  &
            ', Required Field=['//  &
            TRIM(ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName)//   &
            '] was blank.',EchoInputFile)
          ENDIF
          ErrFlag=.true.
        ELSE
          LineItem2%Numbers(NumNumeric)=0.0
          LineItem2%NumBlank(NumNumeric)=.true.
          CALL ShowAuditErrorMessage(' **   Add   ** ','<blank field>   ! field=>'//  &
          TRIM(ObjectDef2(Found)%NumRangeChks(NumNumeric)%FieldName))
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF
ENDIF

IF (.not. ErrFlag .and. .not. IDidntMeanIt) THEN
  IF (TransitionDefer) THEN
    CALL MakeTransition(Found)
  ENDIF
  NumIDFRecords2=NumIDFRecords2+1
  IF (ObjectStartRecord(Found) == 0) ObjectStartRecord(Found)=NumIDFRecords2
  MaxAlphaIDFArgsFound=MAX(MaxAlphaIDFArgsFound,LineItem2%NumAlphas)
  MaxNumericIDFArgsFound=MAX(MaxNumericIDFArgsFound,LineItem2%NumNumbers)
  MaxAlphaIDFDefArgsFound=MAX(MaxAlphaIDFDefArgsFound,ObjectDef2(Found)%NumAlpha)
  MaxNumericIDFDefArgsFound=MAX(MaxNumericIDFDefArgsFound,ObjectDef2(Found)%NumNumeric)
  IDFRecords2(NumIDFRecords2)%Name=LineItem2%Name
  IDFRecords2(NumIDFRecords2)%NumNumbers=LineItem2%NumNumbers
  IDFRecords2(NumIDFRecords2)%NumAlphas=LineItem2%NumAlphas
  IDFRecords2(NumIDFRecords2)%ObjectDefPtr=LineItem2%ObjectDefPtr
  ALLOCATE(IDFRecords2(NumIDFRecords2)%Alphas(LineItem2%NumAlphas))
  ALLOCATE(IDFRecords2(NumIDFRecords2)%AlphBlank(LineItem2%NumAlphas))
  ALLOCATE(IDFRecords2(NumIDFRecords2)%Numbers(LineItem2%NumNumbers))
  ALLOCATE(IDFRecords2(NumIDFRecords2)%NumBlank(LineItem2%NumNumbers))
  IDFRecords2(NumIDFRecords2)%Alphas(1:LineItem2%NumAlphas)=LineItem2%Alphas(1:LineItem2%NumAlphas)
  IDFRecords2(NumIDFRecords2)%AlphBlank(1:LineItem2%NumAlphas)=LineItem2%AlphBlank(1:LineItem2%NumAlphas)
  IDFRecords2(NumIDFRecords2)%Numbers(1:LineItem2%NumNumbers)=LineItem2%Numbers(1:LineItem2%NumNumbers)
  IDFRecords2(NumIDFRecords2)%NumBlank(1:LineItem2%NumNumbers)=LineItem2%NumBlank(1:LineItem2%NumNumbers)
  IF (LineItem2%NumNumbers > 0) THEN
    DO Count=1,LineItem2%NumNumbers
      IF (ObjectDef2(Found)%NumRangeChks(Count)%MinMaxChk .and. .not. LineItem2%NumBlank(Count)) THEN
        CALL InternalRangeCheck(LineItem2%Numbers(Count),Count,Found,LineItem2%Alphas(1),  &
        ObjectDef2(Found)%NumRangeChks(Count)%AutoSizable,        &
        ObjectDef2(Found)%NumRangeChks(Count)%AutoCalculatable)
      ENDIF
    ENDDO
  ENDIF
ELSEIF (.not. IDidntMeanIt) THEN
  OverallErrorFlag=.true.
ENDIF

RETURN

END SUBROUTINE ValidateObjectandParse2  !RS: Debugging: Testing to see if we can use more than one IDD and IDF here (9/22/14)


!     NOTICE
!
!     Copyright © 1996-2012 The Board of Trustees of the University of Illinois
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

END MODULE InputProcessor
