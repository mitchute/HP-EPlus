MODULE FreeCoolingHeatExchanger

  ! Module containing the routines dealing with free cooling hydronic heat
  ! exchanger components.

  ! MODULE INFORMATION:
  !       AUTHOR         Simon Rees
  !       DATE WRITTEN   August 2002
  !       MODIFIED       May 2003
  !       RE-ENGINEERED  March 2012, B.Griffith, reimplemented for new plant

  ! PURPOSE OF THIS MODULE:
  ! The purpose of this module is to simulate a free cooling hydronic heat
  ! exchanger.This component is a simple hydronic heat exchanger that can be
  ! used to couple plant side and condenser side hydronic loops. The purpose
  ! is primarily to enable implementation of free cooling. In other words,
  ! coupling the zone plant directly to the heat rejection equipment. For example,
  ! coupling a radiant hydronic system to a cooling tower and bypassing the
  ! heat pump or chiller.

  ! METHODOLOGY EMPLOYED:
  ! The component is concieved as a simple hydronic heat exchanger but with some
  ! control abilities. Another component (chiller or heat pump) on the
  ! loops is switched on or off by this component. When the other component is
  ! operating the loops are 'uncoupled' and no heat exchange takes place. When
  ! the control conditions are satisfied this heat exchanger switches the other
  ! component off and heat exchange between the loops is calculated. This is to
  ! emulate a bypass arrangement. The heat exchanger can be operated in 'Ideal'
  ! mode, that is the max possible heat exchanger occurs. If the loop fluids are
  ! the same, and the flow rates are the same, complete coupling will occur i.e.
  ! it will be as if the condenser fluid flowed through the plant loop.
  
  ! March 2012, reengineered so that chiller control is moved up to plant loop operation scheme.
  ! Change to use lagged node temperature for stability. 

  ! REFERENCES:

  ! OTHER NOTES: none

  ! USE STATEMENTS:
  ! Use statements for data only modules
USE DataPrecisionGlobals
USE DataGlobals_HPSimIntegrated,       ONLY : MaxNameLength, InitConvTemp,BeginTimeStepFlag, BeginEnvrnFlag, SecInHour !, &
                              !ShowWarningError, ShowSevereError, ShowFatalError, ShowContinueError
USE DataInterfaces,    ONLY : SetupOutputVariable

  ! Use statements for access to subroutines in other modules

IMPLICIT NONE         ! Enforce explicit typing of all variables

PRIVATE ! Everything private unless explicitly made public

  ! MODULE PARAMETER DEFINITIONS
! zero load tolerance used in equipment operation decisions (Watts)
REAL(r64), PARAMETER                         :: zeroloadtol = 1.0d0
!zero capacity tolerance used in heat exchanger cal (kg/s * J/kg/K)
REAL(r64), PARAMETER                         :: zerocaptol  = .00001d0
! free cooling HX type string
CHARACTER(len=*), PARAMETER :: ThisType = 'HeatExchanger:Hydronic'
CHARACTER(len=*), PARAMETER :: ThisTypeUC = 'HEATEXCHANGER:HYDRONIC'
CHARACTER(len=*), PARAMETER :: Blank = ' '
! Control Mode Types
INTEGER, PARAMETER :: ControlMode_WetBulb = 1
INTEGER, PARAMETER :: ControlMode_DryBulb = 2
INTEGER, PARAMETER :: ControlMode_Loop    = 3
INTEGER, PARAMETER :: ControlMode_Schedule= 4
! HX Modes
INTEGER, PARAMETER :: HXMode_UA = 1
INTEGER, PARAMETER :: HXMode_Ideal = 2

  ! DERIVED TYPE DEFINITIONS
TYPE Locator
  INTEGER :: LoopNum    ! Loop number for location
  INTEGER :: LoopSideNum ! Loop side number for location
  INTEGER :: BranchNum  ! Branch number for location (on LoopNum)
  INTEGER :: CompNum    ! Component number number for location (on branch)
END TYPE Locator

TYPE FreeCoolHXData
  ! Input data
  CHARACTER(len=MaxNameLength) :: Name                  = Blank  ! name of the free cooling HX component
  CHARACTER(len=MaxNameLength) :: ComponentName         = Blank  ! name of the component controlled by HX
  CHARACTER(len=MaxNameLength) :: ComponentType         = Blank  ! type of the component controlled by HX
  INTEGER                      :: ComponentIndex        =0       ! Index of component controlled by HX
  CHARACTER(len=MaxNameLength) :: ScheduleName          = Blank  ! Loop schedule name
  CHARACTER(len=MaxNameLength) :: CondInletNode         = Blank  ! condenser side inlet node name
  CHARACTER(len=MaxNameLength) :: CondOutletNode        = Blank  ! condenser side outlet node name
  CHARACTER(len=MaxNameLength) :: PlantInletNode        = Blank  ! plant side inlet node name
  CHARACTER(len=MaxNameLength) :: PlantOutletNode       = Blank  ! plant side outlet node name
  CHARACTER(len=MaxNameLength) :: OtherCondInletNode    = Blank  ! controlled equip. condenser inlet node name
  CHARACTER(len=MaxNameLength) :: OtherCondOutletNode   = Blank  ! controlled equip. condenser outlet node name
  CHARACTER(len=MaxNameLength) :: OtherPlantInletNode   = Blank  ! controlled equip. plant inlet node name
  CHARACTER(len=MaxNameLength) :: OtherPlantOutletNode  = Blank  ! controlled equip. plant outlet node name
  CHARACTER(len=MaxNameLength) :: ControlMode           = Blank  ! use specified control mode (Wetbulb, drybulb etc.)
  INTEGER                      :: ControlMode_Num       =0       ! number equivalend for control node
  INTEGER                      :: HeatXMode             =0       ! heat exchange mode (ideal, NTU-effectiveness)
  REAL(r64)                    :: UA                    =0.0     ! UA for heat exchanger (ignored in ideal mode)
  REAL(r64)                    :: CondenserSideFlowRate =0.0     ! volumetric flow rate through condenser side of unit
  REAL(r64)                    :: CondDesignMassFlowRate =0.d0   ! kg/s
  REAL(r64)                    :: PlantSideFlowRate     =0.0     ! volumetric flow rate through plant side of unit
  REAL(r64)                    :: PlantDesignMassFlowRate =0.d0  ! kg/s
  ! records
  INTEGER                      :: ScheduleNum           =0       ! Loop schedule number
  INTEGER                      :: ComponentSchedNum     =0       ! HX schedule number
  INTEGER                      :: CondInletNodeNum      =0       ! condenser side inlet node number
  INTEGER                      :: CondOutletNodeNum     =0       ! condenser side outlet node number
  INTEGER                      :: PlantInletNodeNum     =0       ! plant side inlet node number
  INTEGER                      :: PlantOutletNodeNum    =0       ! plant side outlet node number
  LOGICAL                      :: CouplingOn             =.false.! coupling status logical flag
  INTEGER                      :: CouplingIntFlag        =0      ! coupling status integer flag - for reporting
  INTEGER                      :: OtherCondOutletNodeNum =0      ! controlled condenser side outlet node number
  INTEGER                      :: OtherPlantOutletNodeNum=0      ! controlled plant side outlet node number
  INTEGER                      :: OtherCondInletNodeNum  =0      ! controlled condenser side inlet node number
  INTEGER                      :: OtherPlantInletNodeNum =0      ! controlled plant side inlet node number
  TYPE (Locator)               :: PlantLocCtrlEquip      =Locator(0,0,0,0)  ! Location on plant side (supply) of controlled equip
  TYPE (Locator)               :: CondLocCtrlEquip       =Locator(0,0,0,0)  ! Location on cond side of controlled equip
!  TYPE (Locator)               :: PlantLocThisEquip      =Locator(0,0,0,0)  ! Location on plant side - this heat exchanger
!  TYPE (Locator)               :: CondLocThisEquip       =Locator(0,0,0,0)  ! Location on cond side - this heat exchanger
  ! Report data
  REAL(r64)                    :: CondInletTemp          =0.0    ! condenser inlet temperature
  REAL(r64)                    :: CondOutletTemp         =0.0    ! condenser outlet temperature
  REAL(r64)                    :: PlantInletTemp         =0.0    ! plant inlet temperature
  REAL(r64)                    :: PlantOutletTemp        =0.0    ! plant outlet temperature
  REAL(r64)                    :: CondMassFlowRate       =0.0    ! condenser mass flow rate
  REAL(r64)                    :: PlantMassFlowRate      =0.0    ! plant mass flow rate
  REAL(r64)                    :: HeatTransRate          =0.0    ! total heat transfer rate, Watts
  REAL(r64)                    :: HeatTransEnergy        =0.0    ! total heat transfer energy, Joules
  !loop topology variables
  INTEGER                      :: CondLoopNum            =0 ! condenser side plant loop number
  INTEGER                      :: CondLoopSideNum        =0 ! condenser side plant loop side number
  INTEGER                      :: CondBranchNum          =0 ! condenser side plant loop branch number
  INTEGER                      :: CondCompNum            =0 ! condenser side plant component number
  INTEGER                      :: PlantLoopNum            =0 ! plant side plant loop number
  INTEGER                      :: PlantLoopSideNum        =0 ! plant side plant loop side number
  INTEGER                      :: PlantBranchNum          =0 ! plant side plant loop branch number
  INTEGER                      :: PlantCompNum            =0 ! plant side plant component number
END TYPE FreeCoolHXData

TYPE(FreeCoolHXData), DIMENSION(:), ALLOCATABLE :: FreeCool
LOGICAL, ALLOCATABLE, DIMENSION(:) :: CheckEquipName

  ! MODULE VARIABLE DECLARATIONS:

INTEGER :: NumOfFreeCools         =0    ! Number of free cooling heat exchangers
INTEGER :: CondInletNodeNum       =0    ! module variable for condenser side inlet node number
INTEGER :: CondOutletNodeNum      =0    ! module variable for condenser side outlet node number
INTEGER :: PlantInletNodeNum      =0    ! module variable for plant side inlet node number
INTEGER :: PlantOutletNodeNum     =0    ! module variable for plant side outlet node number
INTEGER :: OtherCondOutletNodeNum =0    ! module variable for controlled condenser side outlet node number
INTEGER :: OtherPlantOutletNodeNum=0    ! module variable for controlled plant side outlet node number
INTEGER :: OtherCondInletNodeNum  =0    ! module variable for controlled condenser side outlet node number
INTEGER :: OtherPlantInletNodeNum =0    ! module variable for controlled plant side outlet node number
INTEGER :: PlantBranchNum         =0    ! module variable for plant branch number
INTEGER :: PlantCompNum           =0    ! module variable for plant component number
INTEGER :: CondLoopNum            =0    ! module variable for condenser loop number
INTEGER :: OtherCondBranchNum     =0    ! module variable for condenser branch number
INTEGER :: OtherCondCompNum       =0    ! module variable for condenser component number
INTEGER :: PlantLoopNum           =0    ! module variable for plant loop number
INTEGER :: OtherPlantBranchNum    =0    ! module variable for plant branch number
INTEGER :: OtherPlantCompNum      =0    ! module variable for plant component number

  ! SUBROUTINE SPECIFICATIONS FOR MODULE PlantOutsideHX

PUBLIC  SimFreeCoolingHeatExchanger
PRIVATE GetFreeCoolingHeatExchanger
PRIVATE InitFreeCoolingHeatExchanger
PRIVATE CalcFreeCoolingHeatExchanger
PRIVATE UpdateFreeCoolingHeatExchanger
PRIVATE ReportFreeCoolingHeatExchanger


CONTAINS

!==============================================================================

SUBROUTINE SimFreeCoolingHeatExchanger(CompName,CompIndex,FirstHVACIteration,InitLoopEquip)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine is the public interface to this free cooling
          ! heat exchanger component. Control operation and requirement
          ! for coupling the plant side and condenser side are checked here.
          ! Other calcs are made by calling private routines.

          ! METHODOLOGY EMPLOYED:
          ! na

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,  ONLY : FindItemInList
  USE DataEnvironment, ONLY : OutDryBulbTemp, OutWetBulbTemp
  USE ScheduleManager, ONLY : GetCurrentScheduleValue
  USE DataLoopNode,    ONLY : Node
  USE DataHVACGlobals, ONLY : NumCondLoops
  USE General,         ONLY : TrimSigDigits

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  CHARACTER(len=*), INTENT(IN)  :: CompName            ! name of the free cooling heat exchanger.
  INTEGER,        INTENT(INOUT) :: CompIndex           ! index in local derived types
  LOGICAL,          INTENT(IN)  :: FirstHVACIteration  ! TRUE if 1st HVAC simulation of system timestep !unused1208
  LOGICAL,          INTENT(IN)  :: InitLoopEquip       ! for init

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  LOGICAL, SAVE :: GetInputFlag = .TRUE.    ! First time, input is "gotten"
  REAL(r64)     :: DryBulb                  ! Outside drybulb temp
  REAL(r64)     :: WetBulb                  ! Outside wetbulb temp
  REAL(r64)     :: SchedSetpoint            ! scheduled setpoint temp
  REAL(r64)     :: CondTemp                 ! condenser inlet temperature
!  LOGICAL       :: SetCoupled               ! status flag for checking purposes
  INTEGER       :: Item
  INTEGER       :: FreeCoolNum           ! index in local derived types

          ! check for input
  IF (GetInputFlag) THEN
    CALL GetFreeCoolingHeatExchanger
    GetInputFlag=.FALSE.
  ENDIF

  ! Find the correct Equipment
  IF (CompIndex == 0) THEN
    FreeCoolNum = FindItemInList(CompName,FreeCool%Name,NumOfFreeCools)
    IF (FreeCoolNum == 0) THEN
      CALL ShowFatalError('SimFreeCoolingHeatExchanger: Unit not found='//TRIM(CompName))
    ENDIF
    CompIndex=FreeCoolNum
  ELSE
    FreeCoolNum=CompIndex
    IF (FreeCoolNum > NumOfFreeCools .or. FreeCoolNum < 1) THEN
      CALL ShowFatalError('SimFreeCoolingHeatExchanger:  Invalid CompIndex passed='//  &
                          TRIM(TrimSigDigits(FreeCoolNum))// &
                          ', Number of Units='//TRIM(TrimSigDigits(NumOfFreeCools))//  &
                          ', Entered Unit name='//TRIM(CompName))
    ENDIF
    IF (CheckEquipName(FreeCoolNum)) THEN
      IF (CompName /= FreeCool(FreeCoolNum)%Name) THEN
        CALL ShowFatalError('SimFreeCoolingHeatExchanger: Invalid CompIndex passed='//  &
                            TRIM(TrimSigDigits(FreeCoolNum))// &
                            ', Unit name='//TRIM(CompName)//', stored Unit Name for that index='//  &
                            TRIM(FreeCool(FreeCoolNum)%Name))
      ENDIF
      CheckEquipName(FreeCoolNum)=.false.
    ENDIF
  ENDIF

  IF (InitLoopEquip) THEN
    RETURN
  ENDIF

  ! check for condenser loops - may not have been initialized yet
  IF(NumCondLoops<=0)THEN
    DO Item = 1, NumOfFreeCools
      Node(FreeCool(Item)%PlantOutletNodeNum) = Node(FreeCool(Item)%PlantInletNodeNum)
    END DO
  RETURN
  END IF

  ! initialize - set node numbers etc.
  CALL InitFreeCoolingHeatExchanger(FreeCoolNum)

  ! get reference temperatures
  DryBulb = OutDryBulbTemp
  WetBulb = OutWetBulbTemp
  ! get setpoint data
  SchedSetpoint = GetCurrentScheduleValue(FreeCool(FreeCoolNum)%ScheduleNum)

  ! check control mode and sim accordingly
!  IF(Flowlock == 0)THEN
  FreeCool(FreeCoolNum)%CouplingOn = .FALSE. ! local status flag
    ! check control conditions according to control mode
  SELECT CASE(FreeCool(FreeCoolNum)%ControlMode_Num)

  CASE(ControlMode_WetBulb)  ! 'WetBulbTemperature'
    IF(WetBulb < SchedSetpoint)THEN
      FreeCool(FreeCoolNum)%CouplingOn = .TRUE.
    END IF

  CASE(ControlMode_DryBulb)  ! 'DryBulbTemperature'
    IF(DryBulb < SchedSetpoint)THEN
      FreeCool(FreeCoolNum)%CouplingOn = .TRUE.
    END IF

  CASE(ControlMode_Loop) ! 'Loop'

    CondTemp = Node(OtherCondInletNodeNum)%TempLastTimestep
    IF(CondTemp <= SchedSetpoint)THEN
      FreeCool(FreeCoolNum)%CouplingOn = .TRUE.
    END IF

  CASE DEFAULT
        ! invalid mode
    CALL ShowFatalError('SimFreeCoolingHeatExchanger: Invalid free cooling control mode: '//TRIM(CompName))

  END SELECT

  ! calc the heat transfer rate
  CALL CalcFreeCoolingHeatExchanger(FreeCoolNum)
  ! update nodes
  CALL UpdateFreeCoolingHeatExchanger(FreeCoolNum)
  ! update report variables
  CALL ReportFreeCoolingHeatExchanger(FreeCoolNum)

  RETURN

END SUBROUTINE SimFreeCoolingHeatExchanger

!==============================================================================

SUBROUTINE GetFreeCoolingHeatExchanger

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine reads the input for free cooling heat exchangers
          ! from the user input file.  This will contain all of the information
          ! needed to define and simulate the heat exchanger. Some input data
          ! checking is done, and report variables set up.

          ! METHODOLOGY EMPLOYED:
          ! Standard EnergyPlus methodology.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE InputProcessor,     ONLY : GetNumObjectsFound, GetObjectItem, &
                                 FindItemInList, SameString, MakeUPPERcase
  USE DataIPShortCuts
  USE NodeInputManager,   ONLY : GetOnlySingleNode
  USE ScheduleManager,    ONLY : GetScheduleIndex
  USE BranchNodeConnections, ONLY : TestCompSet
  USE DataLoopNode
  USE General,            ONLY: RoundSigDigits
  USE PlantUtilities,     ONLY: RegisterPlantCompDesignFlow

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
          ! na

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='GetFreeCoolingHeatExchanger:'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER                        :: IOStatus   ! Used in GetObjectItem
  INTEGER                        :: Item       ! Item to be "gotten"
  INTEGER                        :: NumAlphas  ! Number of Alphas for each GetObjectItem call
  INTEGER                        :: NumNumbers ! Number of Numbers for each GetObjectItem call
  !INTEGER                        :: NumFluids  ! number of fluids in sim.  !RS: Debugging: Extraneous
  LOGICAL                        :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine
  LOGICAL                        :: IsNotOK=.false.  ! Used to validate component

  cCurrentModuleObject = ThisType

          ! Initializations and allocations
  NumOfFreeCools = GetNumObjectsFound(TRIM(cCurrentModuleObject)) ! 'HeatExchanger:Hydronic'


  ALLOCATE(FreeCool(NumOfFreeCools))
  ALLOCATE(CheckEquipName(NumOfFreeCools))
  CheckEquipName=.true.


          ! Obtain all of the user data
  DO Item = 1, NumOfFreeCools

    CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus, &
                           NumBlank=lNumericFieldBlanks,AlphaBlank=lAlphaFieldBlanks, &
                           AlphaFieldNames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

          ! General user input data
    FreeCool(Item)%Name = cAlphaArgs(1)
    FreeCool(Item)%ComponentName = cAlphaArgs(2)

    FreeCool(Item)%ComponentType = cAlphaArgs(3)
    CALL ValidateComponent(FreeCool(Item)%ComponentType,FreeCool(Item)%ComponentName,IsNotOK,TRIM(cCurrentModuleObject))
    IF (IsNotOK) THEN
      CALL ShowContinueError('Occurs on '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)))
      ErrorsFound=.true.
    ENDIF

    FreeCool(Item)%ScheduleName = cAlphaArgs(4)
          ! Get schedule number
    FreeCool(Item)%ScheduleNum = GetScheduleIndex(FreeCool(Item)%ScheduleName)
    IF (FreeCool(Item)%ScheduleNum == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(4))//'="'//trim(cAlphaArgs(4))//'" not found.')
      ErrorsFound=.true.
    END IF

        ! get condense side inlet node data
    FreeCool(Item)%CondInletNode = cAlphaArgs(5)
    FreeCool(Item)%CondInletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(5),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Inlet, 1, ObjectIsNotParent)
    IF (FreeCool(Item)%CondInletNodeNum == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(5))//'="'//trim(cAlphaArgs(5))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

        ! get condenser side outlet node data
    FreeCool(Item)%CondOutletNode = cAlphaArgs(6)
    FreeCool(Item)%CondOutletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(6),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Outlet, 1, ObjectIsNotParent)
    IF (FreeCool(Item)%CondOutletNodeNum == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(6))//'="'//trim(cAlphaArgs(6))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

    CALL TestCompSet(TRIM(cCurrentModuleObject),cAlphaArgs(1),cAlphaArgs(5),cAlphaArgs(6),'Condenser Water Nodes')


       ! get Plant side inlet node data
    FreeCool(Item)%PlantInletNode = cAlphaArgs(7)
    FreeCool(Item)%PlantInletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(7),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Inlet, 2, ObjectIsNotParent)
    IF (FreeCool(Item)%PlantInletNodeNum == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(7))//'="'//trim(cAlphaArgs(7))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

        ! get Plant side outlet node data
    FreeCool(Item)%PlantOutletNode = cAlphaArgs(8)
    FreeCool(Item)%PlantOutletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(8),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Outlet, 2, ObjectIsNotParent)
    IF (FreeCool(Item)%PlantOutletNodeNum == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(8))//'="'//trim(cAlphaArgs(8))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

    CALL TestCompSet(TRIM(cCurrentModuleObject),cAlphaArgs(1),cAlphaArgs(7),cAlphaArgs(8),'Water Nodes')

        ! Get controlled equipment node data

    FreeCool(Item)%OtherCondInletNode     = cAlphaArgs(12)
    FreeCool(Item)%OtherCondInletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(12),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Actuator, 3, ObjectIsNotParent)
        ! This node is both a sensor and an actuator node, call GetOnlySingleNode to register second connection type
    FreeCool(Item)%OtherCondInletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(12),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Sensor, 3, ObjectIsNotParent)
    IF (FreeCool(Item)%OtherCondInletNodeNum  == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(12))//'="'//trim(cAlphaArgs(12))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

    FreeCool(Item)%OtherCondOutletNode     = cAlphaArgs(13)
    FreeCool(Item)%OtherCondOutletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(13),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Actuator, 3, ObjectIsNotParent)
    IF (FreeCool(Item)%OtherCondOutletNodeNum  == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(13))//'="'//trim(cAlphaArgs(13))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

    FreeCool(Item)%OtherPlantInletNode     = cAlphaArgs(14)
    FreeCool(Item)%OtherPlantInletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(14),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Actuator, 4, ObjectIsNotParent)
    IF (FreeCool(Item)%OtherPlantInletNodeNum  == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(14))//'="'//trim(cAlphaArgs(14))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF

    FreeCool(Item)%OtherPlantOutletNode     = cAlphaArgs(15)
    FreeCool(Item)%OtherPlantOutletNodeNum  =   &
       GetOnlySingleNode(cAlphaArgs(15),ErrorsFound,'Free cooling heat exchanger',cAlphaArgs(1), &
               NodeType_Water,NodeConnectionType_Actuator, 4, ObjectIsNotParent)
    IF (FreeCool(Item)%OtherPlantOutletNodeNum  == 0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(15))//'="'//trim(cAlphaArgs(15))//'" not assigned correctly.')
      ErrorsFound=.true.
    END IF


       ! get control mode
    FreeCool(Item)%ControlMode = MakeUPPERCase(cAlphaArgs(9))
    FreeCool(Item)%ControlMode_Num = 0
    IF (SameString(cAlphaArgs(9),'WetBulbTemperature')) THEN
      FreeCool(Item)%ControlMode_Num=ControlMode_Wetbulb
    ELSEIF (SameString(cAlphaArgs(9),'DryBulbTemperature')) THEN
      FreeCool(Item)%ControlMode_Num=ControlMode_Drybulb
    ELSEIF (SameString(cAlphaArgs(9),'Schedule')) THEN
      FreeCool(Item)%ControlMode_Num=ControlMode_Schedule
    ELSEIF (SameString(cAlphaArgs(9),'Loop')) THEN
      FreeCool(Item)%ControlMode_Num=ControlMode_Loop
    ELSE
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(9))//'="'//trim(cAlphaArgs(9))//'" invalid.')
      CALL ShowContinueError('...value must be one of "WetBulbTemperature", "DryBulbTemperature"'// &
          ' or "Loop"')
      ErrorsFound=.true.
    END IF
    ! overide schedule control - set to loop
    IF(SameString(cAlphaArgs(9),'Schedule')) THEN
      CALL ShowWarningError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(9))//' is no longer used.  Setting to LOOP.')
      FreeCool(Item)%ControlMode = 'LOOP'
      FreeCool(Item)%ControlMode_Num=ControlMode_Loop
    ENDIF

    IF (.not. lAlphaFieldBlanks(10)) THEN
      CALL ShowWarningError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(10))//' is no longer used.')
      CALL ShowContinueError('...ignoring entry="'//trim(cAlphaArgs(10))//'" and using '//  &
         trim(cAlphaFieldNames(4))//'="'//trim(cAlphaArgs(4))//'".')
    ENDIF

       ! get heat exchange mode
    IF (SameString(cAlphaArgs(11),'UFactorTimesAreaEffectiveness')) THEN
      FreeCool(Item)%HeatXMode = HXMode_UA
    ELSEIF (SameString(cAlphaArgs(11),'Ideal')) THEN
      FreeCool(Item)%HeatXMode = HXMode_Ideal
    ENDIF

    IF (.NOT. SameString(cAlphaArgs(11),'Ideal') .AND. &
        .NOT. SameString(cAlphaArgs(11),'UFactorTimesAreaEffectiveness')) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cAlphaFieldNames(11))//'="'//trim(cAlphaArgs(11))//'" invalid.')
      CALL ShowContinueError('...value must be "Ideal" or "UFactorTimesAreaEffectiveness".')
      ErrorsFound=.true.
    END IF

        ! UA effectiveness data
    FreeCool(Item)%UA = rNumericArgs(1)

        ! Flow Rate data
    FreeCool(Item)%CondenserSideFlowRate = rNumericArgs(2)
    FreeCool(Item)%PlantSideFlowRate = rNumericArgs(3)

    IF (rNumericArgs(1) <= 0.0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cNumericFieldNames(1))//'=['//trim(RoundSigDigits(rNumericArgs(1),2))//'] invalid.')
      CALL ShowContinueError('...must be > 0.0.')
      ErrorsFound=.true.
    END IF

    IF (rNumericArgs(2) <= 0.0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cNumericFieldNames(1))//'=['//trim(RoundSigDigits(rNumericArgs(2),2))//'] invalid.')
      CALL ShowContinueError('...must be > 0.0.')
      ErrorsFound=.true.
    END IF

    IF (rNumericArgs(3) <= 0.0) THEN
      CALL ShowSevereError(RoutineName//ThisType//'="'//trim(cAlphaArgs(1))// &
           '", '//TRIM(cNumericFieldNames(3))//'=['//trim(RoundSigDigits(rNumericArgs(3),2))//'] invalid.')
      CALL ShowContinueError('...must be > 0.0.')
      ErrorsFound=.true.
    END IF


  END DO  ! end of input loop



  ! Set up the output variables, CurrentModuleObject='HeatExchanger:Hydronic'
  DO Item = 1, NumOfFreeCools

    ! heat transfer
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Heat Transfer Rate[W]',    &
                              FreeCool(Item)%HeatTransRate,'System','Average', &
                              FreeCool(Item)%Name)
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Heat Transfer Energy[J]',    &
                              FreeCool(Item)%HeatTransEnergy,'System','Sum', FreeCool(Item)%Name, &
                              ResourceTypeKey='ENERGYTRANSFER',EndUseKey='FREECOOLING',GroupKey='Plant')

    ! flow rates
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Condenser Side Mass Flow rate[kg/s]',      &
                              FreeCool(Item)%CondMassFlowRate,'System','Average', &
                              FreeCool(Item)%Name)
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Plant Side Mass Flow rate[kg/s]',      &
                              FreeCool(Item)%PlantMassFlowRate,'System','Average', &
                              FreeCool(Item)%Name)
    ! condenser side temps
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Condenser Side Inlet Temp[C]',     &
                              FreeCool(Item)%CondInletTemp,'System','Average', &
                              FreeCool(Item)%Name)
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Condenser Side Outlet Temp[C]',     &
                              FreeCool(Item)%CondOutletTemp,'System','Average', &
                              FreeCool(Item)%Name)
    ! plant side temps
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Plant Side Inlet Temp[C]',     &
                              FreeCool(Item)%PlantInletTemp,'System','Average', &
                              FreeCool(Item)%Name)
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Plant Side Outlet Temp[C]',     &
                              FreeCool(Item)%PlantOutletTemp,'System','Average', &
                              FreeCool(Item)%Name)
    ! coupling state
    CALL SetupOutputVariable('Free Cooling Heat Exchanger Operation Status[0=Off 1=On]',     &
                              FreeCool(Item)%CouplingIntFlag,'System','Average', &
                              FreeCool(Item)%Name)
  END DO
  ! save the design demand side flow rate for use by plant loop sizing algorithms
  DO Item = 1, NumOfFreeCools
    CALL RegisterPlantCompDesignFlow(FreeCool(Item)%CondInletNodeNum,FreeCool(Item)%CondenserSideFlowRate)
    CALL RegisterPlantCompDesignFlow(FreeCool(Item)%PlantInletNodeNum,FreeCool(Item)%PlantSideFlowRate)
  END DO

  IF (ErrorsFound) THEN
    CALL ShowFatalError('Free Cooling Heat Exchanger: Errors found in input. Preceding conditions cause termination.')
  END IF

  RETURN

END SUBROUTINE GetFreeCoolingHeatExchanger

!==============================================================================

SUBROUTINE InitFreeCoolingHeatExchanger(FreeCoolNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine just sets some local variables used to store
          ! node numbers for current component.

          ! METHODOLOGY EMPLOYED:
          ! Reads free cooling HX component data structure

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:

!  USE DataHVACGlobals, ONLY : NumPlantLoops, NumCondLoops
  USE DataLoopNode,    ONLY : Node
  USE FluidProperties, ONLY : GetDensityGlycol
  USE DataPlant
  USE PlantUtilities,  ONLY : InterConnectTwoPlantLoopSides, InitComponentNodes
  USE InputProcessor,  ONLY : FindItem
  USE ScheduleManager, ONLY : GetCurrentScheduleValue
  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:

  INTEGER, INTENT(IN) :: FreeCoolNum             ! Index for the free cooling heat exchanger.

          ! SUBROUTINE PARAMETER DEFINITIONS:
  CHARACTER(len=*), PARAMETER :: RoutineName='SetLoopCoupled:'

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  INTEGER             :: Num                     ! counter
  !REAL(r64)           :: DesignCondMassFlowRate  ! Mass flow rate used to initialize condenser nodes   !RS: Debugging: Extraneous

  LOGICAL, SAVE       :: OneTimeFlag=.TRUE.
  LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: MyPlantScanFlag
  LOGICAL, ALLOCATABLE, DIMENSION(:), SAVE :: MyEnvrnFlag
  REAL(r64)   :: rho
  LOGICAL :: errFlag
  INTEGER :: ThisTypeOfNum


  ! Do the one time initializations
  IF (OneTimeFlag) THEN
    ALLOCATE(MyPlantScanFlag(NumOfFreeCools))
    ALLOCATE(MyEnvrnFlag(NumOfFreeCools))
    MyPlantScanFlag = .TRUE.
    MyEnvrnFlag = .TRUE.
    OneTimeFlag = .false.
  END IF

  IF (MyPlantScanFlag(FreeCoolNum)) THEN
    ! Locate the hx on the plant loops for later usage
    errFlag=.false.
    CALL ScanPlantLoopsForObject(FreeCool(FreeCoolNum)%Name, &
                                 TypeOf_FreeCoolingHtExchg, &
                                 FreeCool(FreeCoolNum)%CondLoopNum, &
                                 FreeCool(FreeCoolNum)%CondLoopSideNum, &
                                 FreeCool(FreeCoolNum)%CondBranchNum, &
                                 FreeCool(FreeCoolNum)%CondCompNum, &
                                 InletNodeNumber = FreeCool(FreeCoolNum)%CondInletNodeNum,  &
                                 errFlag=errFlag)
    CALL ScanPlantLoopsForObject(FreeCool(FreeCoolNum)%Name, &
                                 TypeOf_FreeCoolingHtExchg, &
                                 FreeCool(FreeCoolNum)%PlantLoopNum, &
                                 FreeCool(FreeCoolNum)%PlantLoopSideNum, &
                                 FreeCool(FreeCoolNum)%PlantBranchNum, &
                                 FreeCool(FreeCoolNum)%PlantCompNum, &
                                 InletNodeNumber = FreeCool(FreeCoolNum)%PlantInletNodeNum,  &
                                 errFlag=errFlag)
    !now test if HX is nominally connected to the proper loop side CR 8276
    IF (FreeCool(FreeCoolNum)%CondLoopSideNum /= DemandSide) THEN ! throw error
      CALL ShowSevereError('Invalid connections for '// ThisType //' = ' //FreeCool(FreeCoolNum)%name )
      CALL ShowContinueError(' The condenser side of component is not connected to the demand side of the loop')
      errFlag = .TRUE.
    ENDIF

    IF (FreeCool(FreeCoolNum)%PlantLoopSideNum /= SupplySide) THEN !throw error
      CALL ShowSevereError('Invalid connections for '// ThisType //' = ' //FreeCool(FreeCoolNum)%name )
      CALL ShowContinueError(' The plant side of component is not connected to the supply side of the loop')
      errFlag = .TRUE.
    ENDIF
    CALL InterConnectTwoPlantLoopSides( FreeCool(FreeCoolNum)%CondLoopNum,      &
                                        FreeCool(FreeCoolNum)%CondLoopSideNum,  &
                                        FreeCool(FreeCoolNum)%PlantLoopNum,     &
                                        FreeCool(FreeCoolNum)%PlantLoopSideNum, &
                                        TypeOf_FreeCoolingHtExchg , .FALSE. )

    ! now get topology locations for the other controlled equipment
    ThisTypeOfNum = FindItem(FreeCool(FreeCoolNum)%ComponentType, SimPlantEquipTypes, NumSimPlantEquipTypes)
    CALL ScanPlantLoopsForObject (FreeCool(FreeCoolNum)%ComponentName, ThisTypeOfNum, &
                                  FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum, &
                                  FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum, &
                                  FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum,  &
                                  FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum, &
                                  InletNodeNumber = FreeCool(FreeCoolNum)%OtherPlantInletNodeNum,  &
                                  errFlag=errFlag)

    CALL ScanPlantLoopsForObject (FreeCool(FreeCoolNum)%ComponentName, ThisTypeOfNum, &
                                  FreeCool(FreeCoolNum)%CondLocCtrlEquip%LoopNum, &
                                  FreeCool(FreeCoolNum)%CondLocCtrlEquip%LoopSideNum, &
                                  FreeCool(FreeCoolNum)%CondLocCtrlEquip%BranchNum,  &
                                  FreeCool(FreeCoolNum)%CondLocCtrlEquip%CompNum, &
                                  InletNodeNumber =  FreeCool(FreeCoolNum)%OtherCondInletNodeNum,  &
                                  errFlag=errFlag)


    IF (errFlag) THEN
      CALL ShowFatalError('InitFreeCoolingHeatExchanger: Program terminated due to previous condition(s).')
    ENDIF
    ! revise how loads served category for other controlled equipment
    SELECT CASE (PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
       Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)%HowLoadServed)

    CASE (HowMet_ByNominalCap)
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)%HowLoadServed &
          = HowMet_ByNominalCapFreeCoolCntrl
    CASE (HowMet_ByNominalCapLowOutLimit)
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)%HowLoadServed &
          = HowMet_ByNominalCapLowOutLimitFreeCoolCntrl
    END SELECT
    
    SELECT CASE(FreeCool(FreeCoolNum)%ControlMode_Num)
    CASE(ControlMode_WetBulb)  
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)% &
          FreeCoolCntrlMode = FreeCoolControlMode_WetBulb
    CASE(ControlMode_DryBulb)  
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)% &
          FreeCoolCntrlMode = FreeCoolControlMode_DryBulb
    CASE(ControlMode_Loop)
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)% &
         FreeCoolCntrlMode = FreeCoolControlMode_Loop
      PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
        Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)% &
          FreeCoolCntrlNodeNum =  FreeCool(FreeCoolNum)%OtherCondInletNodeNum
    END SELECT

    MyPlantScanFlag(FreeCoolNum)=.FALSE.
  ENDIF

  IF (BeginEnvrnFlag .and. MyEnvrnFlag(FreeCoolNum) .AND. .NOT. MyPlantScanFlag(FreeCoolNum) ) THEN  ! .OR. OneTimeInitFlag) THEN
    DO Num = 1, NumOfFreeCools

      rho = GetDensityGlycol(PlantLoop( FreeCool(FreeCoolNum)%CondLoopNum )%FluidName, &
                            InitConvTemp, &
                            PlantLoop( FreeCool(FreeCoolNum)%CondLoopNum )%FluidIndex , &
                            'InitFreeCoolingHeatExchanger')

      FreeCool(Num)%CondDesignMassFlowRate = FreeCool(Num)%CondenserSideFlowRate * rho

      CALL InitComponentNodes(0.d0, FreeCool(Num)%CondDesignMassFlowRate, &
                             FreeCool(Num)%CondInletNodeNum, &
                             FreeCool(Num)%CondOutletNodeNum, &
                             FreeCool(Num)%CondLoopNum, &
                             FreeCool(Num)%CondLoopSideNum, &
                             FreeCool(Num)%CondBranchNum, &
                             FreeCool(Num)%CondCompNum)
      rho = GetDensityGlycol(PlantLoop( FreeCool(FreeCoolNum)%PlantLoopNum )%FluidName, &
                            InitConvTemp, &
                            PlantLoop( FreeCool(FreeCoolNum)%PlantLoopNum )%FluidIndex , &
                            'InitFreeCoolingHeatExchanger')

      FreeCool(Num)%PlantDesignMassFlowRate = FreeCool(Num)%PlantSideFlowRate * rho

      CALL InitComponentNodes(0.d0, FreeCool(Num)%PlantDesignMassFlowRate, &
                             FreeCool(Num)%PlantInletNodeNum, &
                             FreeCool(Num)%PlantOutletNodeNum, &
                             FreeCool(Num)%PlantLoopNum, &
                             FreeCool(Num)%PlantLoopSideNum, &
                             FreeCool(Num)%PlantBranchNum, &
                             FreeCool(Num)%PlantCompNum)
    END DO
    MyEnvrnFlag(FreeCoolNum) = .FALSE.
  END IF
  IF (.not. BeginEnvrnFlag) THEN
    MyEnvrnFlag(FreeCoolNum)=.true.
  ENDIF

  ! Init more variables

  !  set module variables for node numbers for this heat exchanger
  CondInletNodeNum   = FreeCool(FreeCoolNum)%CondInletNodeNum
  CondOutletNodeNum  = FreeCool(FreeCoolNum)%CondOutletNodeNum
  PlantInletNodeNum  = FreeCool(FreeCoolNum)%PlantInletNodeNum
  PlantOutletNodeNum = FreeCool(FreeCoolNum)%PlantOutletNodeNum
  ! variables for controlled equipment
  OtherCondInletNodeNum   = FreeCool(FreeCoolNum)%OtherCondInletNodeNum
  OtherCondOutletNodeNum  = FreeCool(FreeCoolNum)%OtherCondOutletNodeNum
  OtherPlantInletNodeNum  = FreeCool(FreeCoolNum)%OtherPlantInletNodeNum
  OtherPlantOutletNodeNum = FreeCool(FreeCoolNum)%OtherPlantOutletNodeNum

  ! store current value for control schedule in central plant loop data structure
  PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)%LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)% &
    Branch(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)% &
      FreeCoolCntrlMinCntrlTemp =  GetCurrentScheduleValue(FreeCool(FreeCoolNum)%ScheduleNum)
      
END SUBROUTINE InitFreeCoolingHeatExchanger

!==============================================================================

SUBROUTINE CalcFreeCoolingHeatExchanger(FreeCoolNum)

          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  March 2012, B.Griffith, revised for new plant, new approach

          ! PURPOSE OF THIS SUBROUTINE:
          ! This calculates the total heat transfer rate between the
          ! two loop fluids streams. This heat transfer rate is used
          ! in the update routine to calc the node temps.

          ! METHODOLOGY EMPLOYED:
          ! NTU-effectiveness heat exchanger model. Effectiveness is
          ! calculated from the user supplied UA value. If 'Ideal' mode
          ! has been set, effectiveness is set to 1.0.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataLoopNode,    ONLY  : Node
  USE FluidProperties, ONLY  : GetSpecificHeatGlycol
  USE DataPlant,       ONLY : PlantLoop
  USE PlantUtilities,  ONLY : SetComponentFlowRate

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN)  :: FreeCoolNum          ! Index for the free cooling heat exchanger.

          ! SUBROUTINE PARAMETER DEFINITIONS:

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  REAL(r64)    :: PlantFluidCp        ! Specific heat of Plant side fluid
  REAL(r64)    :: PlantMassFlowRate   ! Flow rate of Plant side fluid
  REAL(r64)    :: PlantCapRate        ! Capacity rate (mdot*Cp) of Plant side fluid
  REAL(r64)    :: PlantInletTemp      ! Plant side inlet temperature
  REAL(r64)    :: CondFluidCp         ! Specific heat of condenser side fluid

  REAL(r64)    :: CondMassFlowRate    ! Flow rate of condenser side fluid
  REAL(r64)    :: CondCapRate         ! Capacity rate (mdot*Cp) of condenser side fluid
  REAL(r64)    :: CondInletTemp       ! condenser side inlet temperature
  REAL(r64)    :: MinCapRate          ! minimum capacity rate
  REAL(r64)    :: CapRatio            ! capacity ratio (min/max)
  REAL(r64)    :: Effectiveness       ! heat exchanger effectiveness
  REAL(r64)    :: NTU                 ! dimensionless NTU calculated from UA
  !REAL(r64)    :: ChillerLoad         ! current load on chiller (Myload)   !RS: Debugging: Extraneous
  INTEGER :: LoopNum
  INTEGER :: LoopSideNum
  REAL(r64)  :: mdot
  LOGICAL :: ChillerShutDown = .FALSE. 


  LoopNum = FreeCool(FreeCoolNum)%PlantLoopNum
  LoopSideNum = FreeCool(FreeCoolNum)%PlantLoopSideNum
  ChillerShutDown = PlantLoop(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopNum)% &
     LoopSide(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%LoopSideNum)%Branch(FreeCool(FreeCoolNum)% &
     PlantLocCtrlEquip%BranchNum)%Comp(FreeCool(FreeCoolNum)%PlantLocCtrlEquip%CompNum)%FreeCoolCntrlShutDown

  IF(ChillerShutDown .AND. FreeCool(FreeCoolNum)%CouplingOn )THEN

    ! set the heat exchanger plant inlet  flow rate to heat exchanger design flow
    ! this will allow the plant loop flow resolvers to set flow rate properly
    mdot = FreeCool(FreeCoolNum)%PlantDesignMassFlowRate
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%PlantInletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantOutletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantLoopNum, &
                                FreeCool(FreeCoolNum)%PlantLoopSideNum, &
                                FreeCool(FreeCoolNum)%PlantBranchNum,  &
                                FreeCool(FreeCoolNum)%PlantCompNum)

    ! set the heat exchanger condenser inlet side flow rate to heat exchanger design flow
    mdot = FreeCool(FreeCoolNum)%CondDesignMassFlowRate
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%CondInletNodeNum, &
                                FreeCool(FreeCoolNum)%CondOutletNodeNum, &
                                FreeCool(FreeCoolNum)%CondLoopNum, &
                                FreeCool(FreeCoolNum)%CondLoopSideNum, &
                                FreeCool(FreeCoolNum)%CondBranchNum,  &
                                FreeCool(FreeCoolNum)%CondCompNum)


    ! If there is a NO load on the chiller, if coupling is 'On', and if flowlock = 0 (request mode),
    ! then set all flow rates to zero
  ELSE IF( .NOT. ChillerShutDown .AND. FreeCool(FreeCoolNum)%CouplingOn )THEN
    ! set the heat exchanger plant inlet  flow rate to 0.0
    ! this will allow the plant loop flow resolvers to set flow rate properly
    mdot = 0.d0
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%PlantInletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantOutletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantLoopNum, &
                                FreeCool(FreeCoolNum)%PlantLoopSideNum, &
                                FreeCool(FreeCoolNum)%PlantBranchNum,  &
                                FreeCool(FreeCoolNum)%PlantCompNum)

    ! set the heat exchanger condenser inlet side flow rate to 0.0
    mdot = 0.d0
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%CondInletNodeNum, &
                                FreeCool(FreeCoolNum)%CondOutletNodeNum, &
                                FreeCool(FreeCoolNum)%CondLoopNum, &
                                FreeCool(FreeCoolNum)%CondLoopSideNum, &
                                FreeCool(FreeCoolNum)%CondBranchNum,  &
                                FreeCool(FreeCoolNum)%CondCompNum)



    ! If coupling is 'Off', and if flowlock = 0 (request mode),
    ! then set all heat exchanger flow rates to zero, but don't change chiller flow rates.
    ! The chiller, which is now set to 'Active' mode will reset it's own flow rates.
  ELSE IF(.NOT. FreeCool(FreeCoolNum)%CouplingOn )THEN

    ! set the heat exchanger plant inlet  flow rate to 0.0
    mdot = 0.d0
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%PlantInletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantOutletNodeNum, &
                                FreeCool(FreeCoolNum)%PlantLoopNum, &
                                FreeCool(FreeCoolNum)%PlantLoopSideNum, &
                                FreeCool(FreeCoolNum)%PlantBranchNum,  &
                                FreeCool(FreeCoolNum)%PlantCompNum)

    ! set the heat exchanger condenser inlet side flow rate to 0.0
    mdot = 0.d0
    CALL SetComponentFlowRate( mdot , &
                                FreeCool(FreeCoolNum)%CondInletNodeNum, &
                                FreeCool(FreeCoolNum)%CondOutletNodeNum, &
                                FreeCool(FreeCoolNum)%CondLoopNum, &
                                FreeCool(FreeCoolNum)%CondLoopSideNum, &
                                FreeCool(FreeCoolNum)%CondBranchNum,  &
                                FreeCool(FreeCoolNum)%CondCompNum)

  ENDIF
    !set local variables for heat exchanger calculation
  CondInletTemp     = Node(CondInletNodeNum)%Temp
  PlantInletTemp    = Node(PlantInletNodeNum)%Temp
  CondMassFlowRate  = Node(CondInletNodeNum)%MassFlowRate
  PlantMassFlowRate = Node(PlantInletNodeNum)%MassFlowRate
  CondFluidCp    = 0.0
  PlantFluidCp   = 0.0

  IF (CondMassFlowRate /= 0.0) &  ! When Mass flow rate = 0, end result will be 0.0
    CondFluidCp    = GetSpecificHeatGlycol(PlantLoop(FreeCool(FreeCoolNum)%CondLoopNum)%FluidName,CondInletTemp, &
                                           PlantLoop(FreeCool(FreeCoolNum)%CondLoopNum)%FluidIndex,  &
                                           'FreeCoolHeatEx='//TRIM(FreeCool(FreeCoolNum)%Name))
  IF (PlantMassFlowRate /= 0.0) &  ! When Mass flow rate = 0, end result will be 0.0
    PlantFluidCp   = GetSpecificHeatGlycol(PlantLoop( FreeCool(FreeCoolNum)%PlantLoopNum)%FluidName, PlantInletTemp, &
                                           PlantLoop( FreeCool(FreeCoolNum)%PlantLoopNum)%FluidIndex,  &
                                           'FreeCoolHeatEx='//TRIM(FreeCool(FreeCoolNum)%Name))

  CondCapRate  = CondFluidCp * CondMassFlowRate
  PlantCapRate = PlantFluidCp * PlantMassFlowRate
  MinCapRate = MIN(CondCapRate, PlantCapRate)

    !If there is no flow rate on either the condenser or plant side, turn off heat exchanger and return
  IF (CondCapRate <= zerocaptol .OR. PlantCapRate <= zerocaptol)THEN
    FreeCool(FreeCoolNum)%HeatTransRate = 0.0
    RETURN
  END IF

  ! calc effectiveness - 1.0 if in ideal mode
  IF(FreeCool(FreeCoolNum)%CouplingOn)THEN
    IF(FreeCool(FreeCoolNum)%HeatXMode == HXMode_UA)THEN
      ! assume cross flow, both mixed
      CapRatio = MinCapRate/MAX(CondCapRate, PlantCapRate)
      NTU = FreeCool(FreeCoolNum)%UA/MinCapRate

      Effectiveness = 1.0 - EXP( (NTU**0.22d0/CapRatio) * &
                  (EXP(-CapRatio*NTU**0.78d0) - 1.0) )
    ELSE
      ! must be in ideal mode
      Effectiveness = 1.0
    END IF
    ! overall heat transfer rate
    ! convention is +ve rate is rejected from plant side to condenser side
    FreeCool(FreeCoolNum)%HeatTransRate = Effectiveness*MinCapRate*(PlantInletTemp-CondInletTemp)
  ELSE
    ! explicitly set heat transfer rate to zero if not coupled
    FreeCool(FreeCoolNum)%HeatTransRate = 0.0
  END IF

  RETURN

END SUBROUTINE CalcFreeCoolingHeatExchanger

!==============================================================================

SUBROUTINE UpdateFreeCoolingHeatExchanger(FreeCoolNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! This subroutine takes the inlet conditions and the previously
          ! calculated heat transfer rate to calculate the outlet temperatures.
          ! All flow rates are passed through. If the loops are not coupled
          ! All node info is passed through from inlet to outlet

          ! METHODOLOGY EMPLOYED:
          ! use previously calcultated heat transfer rate and update node data.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataLoopNode,    ONLY : Node
  USE FluidProperties, ONLY : GetSpecificHeatGlycol
  USE PlantUtilities,  ONLY : SafeCopyPlantNode
  USE DataPlant,       ONLY: PlantLoop

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: FreeCoolNum  ! Index for the free cooling heat exchanger.

          ! SUBROUTINE PARAMETER DEFINITIONS:
          ! na

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  REAL(r64)    :: CondFluidCp         ! Specific heat of condenser side fluid
  REAL(r64)    :: CondMassFlowRate    ! Flow rate of condenser side fluid
  REAL(r64)    :: PlantFluidCp        ! Specific heat of Plant side fluid
  REAL(r64)    :: PlantMassFlowRate   ! Flow rate of Plant side fluid
  REAL(r64)    :: CondInletTemp       ! condenser side inlet temperature
  REAL(r64)    :: PlantInletTemp      ! Plant side inlet temperature


  ! flow rates
  CondMassFlowRate  = Node(CondInletNodeNum)%MassFlowRate
  PlantMassFlowRate = Node(PlantInletNodeNum)%MassFlowRate

  ! inlet temps and specific heats
  CondInletTemp  = Node(CondInletNodeNum)%Temp
  PlantInletTemp = Node(PlantInletNodeNum)%Temp
  CondFluidCp    = 0.0
  PlantFluidCp   = 0.0

!  IF (CondMassFlowRate /= 0.0) &
  CondFluidCp    = GetSpecificHeatGlycol(PlantLoop(FreeCool(FreeCoolNum)%CondLoopNum)%FluidName,CondInletTemp, &
                                         PlantLoop(FreeCool(FreeCoolNum)%CondLoopNum)%FluidIndex,  &
                                         'FreeCoolHeatEx='//TRIM(FreeCool(FreeCoolNum)%Name))
!  IF (PlantMassFlowRate /= 0.0) &
  PlantFluidCp   = GetSpecificHeatGlycol(PlantLoop(FreeCool(FreeCoolNum)%PlantLoopNum)%FluidName, PlantInletTemp, &
                                         PlantLoop(FreeCool(FreeCoolNum)%PlantLoopNum)%FluidIndex,  &
                                         'FreeCoolHeatEx='//TRIM(FreeCool(FreeCoolNum)%Name))

  ! always mass through flow rates
  CALL SafeCopyPlantNode(CondInletNodeNum, CondOutletNodeNum)

  CALL SafeCopyPlantNode(PlantInletNodeNum, PlantOutletNodeNum)

  ! check if coupled or zero heat transfer rate
  IF(FreeCool(FreeCoolNum)%CouplingOn .AND. FreeCool(FreeCoolNum)%HeatTransRate /= 0.0) THEN
    ! calc outlet temps from heat transfer rate
    Node(CondOutletNodeNum)%Temp  = Node(CondInletNodeNum)%Temp + FreeCool(FreeCoolNum)%HeatTransRate/ &
                                                           (CondMassFlowRate * CondFluidCp)

    Node(PlantOutletNodeNum)%Temp = Node(PlantInletNodeNum)%Temp - FreeCool(FreeCoolNum)%HeatTransRate/ &
                                                            (PlantMassFlowRate * PlantFluidCp)
  ELSE
    ! just pass through
    Node(CondOutletNodeNum)%Temp  = Node(CondInletNodeNum)%Temp
    Node(PlantOutletNodeNum)%Temp = Node(PlantInletNodeNum)%Temp
  END IF

  RETURN


END SUBROUTINE UpdateFreeCoolingHeatExchanger

!==============================================================================

SUBROUTINE ReportFreeCoolingHeatExchanger(FreeCoolNum)

          ! SUBROUTINE INFORMATION:
          !       AUTHOR         Simon Rees
          !       DATE WRITTEN   August 2002
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS SUBROUTINE:
          ! Updates the free cooling heat exchanger variables used for reporting.

          ! METHODOLOGY EMPLOYED:
          ! Update variables from node data.

          ! REFERENCES:
          ! na

          ! USE STATEMENTS:
  USE DataLoopNode,    ONLY : Node
  USE DataHVACGlobals, ONLY : TimeStepSys

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

          ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: FreeCoolNum  ! Index for the free cooling heat exchanger.

          ! SUBROUTINE PARAMETER DEFINITIONS:

          ! INTERFACE BLOCK SPECIFICATIONS
          ! na

          ! DERIVED TYPE DEFINITIONS
          ! na

          ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

  ! update condenser side variables
  FreeCool(FreeCoolNum)%CondInletTemp     = Node(CondInletNodeNum)%Temp
  FreeCool(FreeCoolNum)%CondOutletTemp    = Node(CondOutletNodeNum)%Temp
  FreeCool(FreeCoolNum)%CondMassFlowRate  = Node(CondInletNodeNum)%MassFlowRate

  ! update plant side variables
  FreeCool(FreeCoolNum)%PlantInletTemp    = Node(PlantInletNodeNum)%Temp
  FreeCool(FreeCoolNum)%PlantOutletTemp   = Node(PlantOutletNodeNum)%Temp
  FreeCool(FreeCoolNum)%PlantMassFlowRate = Node(PlantInletNodeNum)%MassFlowRate

  ! update the energy reporting variable
  FreeCool(FreeCoolNum)%HeatTransEnergy   = FreeCool(FreeCoolNum)%HeatTransRate*TimeStepSys*SecInHour


  ! update status record
  IF(FreeCool(FreeCoolNum)%CouplingOn)THEN
    FreeCool(FreeCoolNum)%CouplingIntFlag = 1
  ELSE
    FreeCool(FreeCoolNum)%CouplingIntFlag = 0
  END IF

  RETURN

END SUBROUTINE ReportFreeCoolingHeatExchanger

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

END MODULE FreeCoolingHeatExchanger
