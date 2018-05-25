MODULE PlateCoolingHeatExchanger

  ! Module containing the routines dealing with Plate heat
  ! exchanger components required for coupling two condenser
  ! loops (reqd. for Hybrid systems)

  ! MODULE INFORMATION:
  !       AUTHOR         Sankaranarayanan K P
  !       DATE WRITTEN   NOV 2003
  !       MODIFIED       JAN 2004
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS MODULE:
  ! The purpose of this module is to simulate a Plate heat exchanger.
  ! This component is a simple hydronic heat exchanger that can be
  ! used to couple two condenser hydronic loops as in hybrid systems. The purpose
  ! is primarily to enable implementation of free cooling.

  ! METHODOLOGY EMPLOYED:
  ! The component is concieved as a simple counter flow heat exchanger.
  ! The coupling is based on whether the component on second loop is on or off.
  ! When the other component is operating the loops are 'coupled' and heat exchange
  ! takes place. The heat exchanger can be operated in 'Ideal' mode, that is the max
  ! possible heat exchanger occurs. If the loop fluids are the same, and the flow rates
  ! are the same, complete coupling will occur i.e. it will be as if the condenser fluid
  ! flowed through the other condenser loop.

  ! REFERENCES:

  ! OTHER NOTES: none

  ! USE STATEMENTS:
  ! Use statements for data only modules
  USE DataPrecisionGlobals
  USE DataGlobals_HPSimIntegrated,       ONLY : MaxNameLength, InitConvTemp, BeginTimeStepFlag, BeginEnvrnFlag !, &
  !ShowWarningError, ShowSevereError, ShowFatalError, ShowContinueError
  USE DataInterfaces,    ONLY : SetupOutputVariable
  !USE DataPlant,         ONLY : CondSupplySide
  USE DataHVACGlobals,    ONLY : NumCondLoops
  USE Psychrometrics
  USE DataLoopNode

  ! Use statements for access to subroutines in other modules

  IMPLICIT NONE         ! Enforce explicit typing of all variables

  PRIVATE ! Everything private unless explicitly made public

  ! MODULE PARAMETER DEFINITIONS
  ! zero load tolerance used in equipment operation decisions (Watts)
  REAL(r64), PARAMETER                         :: zeroloadtol = 1.d0
  !zero capacity tolerance used in heat exchanger cal (kg/s * J/kg/K)
  REAL(r64), PARAMETER                         :: zerocaptol  = .00001d0

  !  CHARACTER(len=MaxNameLength), PARAMETER :: ThisType = 'HeatExchanger:Plate'

  ! DERIVED TYPE DEFINITIONS
  TYPE PlateCoolHXData
    ! Input data
    CHARACTER(len=MaxNameLength) :: Name                 =' '    ! name of the free cooling HX component
    CHARACTER(len=MaxNameLength) :: ComponentName        =' '    ! name of the component controlled by HX
    CHARACTER(len=MaxNameLength) :: ComponentType        =' '    ! type of the component controlled by HX
    CHARACTER(len=MaxNameLength) :: DemandSideInletNode  =' '    ! condenser side inlet node name
    CHARACTER(len=MaxNameLength) :: DemandSideOutletNode =' '    ! condenser side outlet node name
    CHARACTER(len=MaxNameLength) :: SupplySideInletNode  =' '    ! plant side inlet node name
    CHARACTER(len=MaxNameLength) :: SupplySideOutletNode =' '    ! plant side outlet node name
    CHARACTER(len=MaxNameLength) :: HeatXMode            =' '    ! heat exchange mode (ideal, NTU-effectiveness)
    REAL(r64)                    :: UA                   =0.0    ! UA for heat exchanger (ignored in ideal mode)
    REAL(r64)                    :: DemandSideFlowRate   =0.0    ! volumetric flow rate through condenser side of unit
    REAL(r64)                    :: SupplySideFlowRate   =0.0    ! volumetric flow rate through plant side of unit
    INTEGER                      :: DemandSideInletNodeNum  =0   ! condenser side inlet node number
    INTEGER                      :: DemandSideOutletNodeNum =0   ! condenser side outlet node number
    INTEGER                      :: SupplySideInletNodeNum  =0   ! plant side inlet node number
    INTEGER                      :: SupplySideOutletNodeNum =0   ! plant side outlet node number
    ! Report data
    REAL(r64)                    :: DemandSideInletTemp     =0.0 ! condenser inlet temperature
    REAL(r64)                    :: DemandSideOutletTemp    =0.0 ! condenser outlet temperature
    REAL(r64)                    :: SupplySideInletTemp     =0.0 ! plant inlet temperature
    REAL(r64)                    :: SupplySideOutletTemp    =0.0 ! plant outlet temperature
    REAL(r64)                    :: DemandSideMassFlowRate  =0.0 ! condenser mass flow rate
    REAL(r64)                    :: SupplySideMassFlowRate  =0.0 ! plant mass flow rate
    REAL(r64)                    :: HeatTransRate           =0.0 ! total heat transfer rate, Watts
    !loop topology variables
    INTEGER                      :: CondLoopNum            =0 ! condenser side plant loop number
    INTEGER                      :: CondLoopSideNum        =0 ! condenser side plant loop side number
    INTEGER                      :: CondBranchNum          =0 ! condenser side plant loop branch number
    INTEGER                      :: CondCompNum            =0 ! condenser side plant component number
    INTEGER                      :: PlantLoopNum            =0 ! plant side plant loop number
    INTEGER                      :: PlantLoopSideNum        =0 ! plant side plant loop side number
    INTEGER                      :: PlantBranchNum          =0 ! plant side plant loop branch number
    INTEGER                      :: PlantCompNum            =0 ! plant side plant component number
    LOGICAL                      :: MyFlag                 =.TRUE.
  END TYPE PlateCoolHXData

  TYPE(PlateCoolHXData), DIMENSION(:), ALLOCATABLE :: PlateCool

  ! MODULE VARIABLE DECLARATIONS:

  INTEGER   :: NumOfFreeCools           =0   ! Number of free cooling heat exchangers
  INTEGER   :: DemandSideInletNodeNum   =0   ! module variable for condenser side inlet node number
  INTEGER   :: DemandSideOutletNodeNum  =0   ! module variable for condenser side outlet node number
  INTEGER   :: SupplySideInletNodeNum   =0   ! module variable for plant side inlet node number
  INTEGER   :: SupplySideOutletNodeNum  =0   ! module variable for plant side outlet node number
  INTEGER   :: SupplySideBranchNum      =0   ! module variable for plant branch number
  INTEGER   :: SupplySideCompNum        =0   ! module variable for plant component number
  INTEGER   :: DemandSideLoopNum        =0   ! module variable for condenser loop number
  INTEGER   :: SupplySideLoopNum        =0   ! module variable for plant loop number
  REAL(r64) :: SupplySideMassFlowRate   =0.0 ! Flow rate of Plant side fluid
  REAL(r64) :: DemandSideMassFlowRate   =0.0 ! Flow rate of condenser side fluid
  REAL(r64) :: SupplySideInletTemp      =0.0 ! SupplySide side inlet temperature
  REAL(r64) :: DemandSideInletTemp      =0.0 ! condenser side inlet temperature
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: CheckEquipName


  ! SUBROUTINE SPECIFICATIONS FOR MODULE PlantOutsideHX

  PUBLIC  SimPlateCoolingHeatExchanger
  PRIVATE GetPlateCoolingHeatExchanger
  PRIVATE InitPlateCoolingHeatExchanger
  PRIVATE CalcPlateCoolingHeatExchanger
  PRIVATE TurnOffPHXIfNotNeeded
  PRIVATE UpdatePlateCoolingHeatExchanger
  PRIVATE ReportPlateCoolingHeatExchanger

CONTAINS

  !==============================================================================

  SUBROUTINE SimPlateCoolingHeatExchanger(CompName,CompIndex,FirstHVACIteration,RunFlag,InitLoopEquip,DemandSideCheck)

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Sankaranarayanan K P
    !       DATE WRITTEN   Nov 2003
    !       MODIFIED       na
    !       RE-ENGINEERED  na

    ! PURPOSE OF THIS SUBROUTINE:
    ! This subroutine is the public interface to this Plate cooling
    ! heat exchanger component. Control operation and requirement
    ! for coupling the plant side and condenser side are checked here.
    ! Other calcs are made by calling private routines.

    ! METHODOLOGY EMPLOYED:
    ! Needs description, as appropriate.

    ! REFERENCES:
    ! na

    ! USE STATEMENTS:
    USE InputProcessor,  ONLY : FindItemInList
    USE DataEnvironment, ONLY : OutDryBulbTemp, OutWetBulbTemp
    USE ScheduleManager, ONLY : GetCurrentScheduleValue
    USE DataLoopNode,    ONLY : Node
    USE DataHVACGlobals, ONLY : NumCondLoops
    USE General,         ONLY : TrimSigDigits
    USE PlantUtilities,  ONLY : SafeCopyPlantNode

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE ARGUMENT DEFINITIONS:
    CHARACTER(len=*), INTENT(IN)    :: CompName            ! name of the Plate cooling heat exchanger.
    INTEGER,          INTENT(INOUT) :: CompIndex           ! index in local derived types
    LOGICAL,          INTENT(IN)    :: FirstHVACIteration  ! TRUE if 1st HVAC simulation of system timestep
    LOGICAL,          INTENT(IN)    :: RunFlag             ! Run condition flag
    LOGICAL                         :: InitLoopEquip
    LOGICAL, OPTIONAL, INTENT(IN)   :: DemandSideCheck

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
    LOGICAL, SAVE :: GetInputFlag = .TRUE.    ! First time, input is "gotten"
    INTEGER       :: PlateCoolNum             ! index in local derived types
    INTEGER       :: Item

    ! check for input
    IF (GetInputFlag) THEN
      CALL GetPlateCoolingHeatExchanger
      GetInputFlag=.FALSE.
    ENDIF

    IF (CompIndex == 0) THEN
      PlateCoolNum = FindItemInList(CompName,PlateCool%Name,NumOfFreeCools)
      IF (PlateCoolNum == 0) THEN
        CALL ShowFatalError('SimPlateCoolingHeatExchanger: Unit not found='//TRIM(CompName))
      ENDIF
      CompIndex=PlateCoolNum
    ELSE
      PlateCoolNum=CompIndex
      IF (PlateCoolNum > NumOfFreeCools .or. PlateCoolNum < 1) THEN
        CALL ShowFatalError('SimPlateCoolingHeatExchanger:  Invalid CompIndex passed='//  &
        TRIM(TrimSigDigits(PlateCoolNum))// &
        ', Number of Units='//TRIM(TrimSigDigits(NumOfFreeCools))//  &
        ', Entered Unit name='//TRIM(CompName))
      ENDIF
      IF (CheckEquipName(PlateCoolNum)) THEN
        IF (CompName /= PlateCool(PlateCoolNum)%Name) THEN
          CALL ShowFatalError('SimPlateCoolingHeatExchanger: Invalid CompIndex passed='//  &
          TRIM(TrimSigDigits(PlateCoolNum))// &
          ', Unit name='//TRIM(CompName)//', stored Unit Name for that index='//  &
          TRIM(PlateCool(PlateCoolNum)%Name))
        ENDIF
        CheckEquipName(PlateCoolNum)=.false.
      ENDIF
    ENDIF

    IF (InitLoopEquip) THEN
      RETURN
    ENDIF

    ! check for condenser loops - may not have been initialized yet so just pass information across
    IF(NumCondLoops<=0)THEN
      DO Item = 1, NumOfFreeCools
        CALL SafeCopyPlantNode(PlateCool(Item)%SupplySideInletNodeNum, PlateCool(Item)%SupplySideOutletNodeNum)
      END DO
      RETURN
    END IF

    ! If we are being called from the demand side, just check to see if we need to run and turn off if not
    IF (PRESENT(DemandSideCheck)) THEN
      CALL TurnOffPHXIfNotNeeded(PlateCoolNum)
      RETURN
    END IF

    ! initialize - set node numbers etc.
    CALL InitPlateCoolingHeatExchanger(PlateCoolNum)
    ! calc the heat transfer rate
    CALL CalcPlateCoolingHeatExchanger(PlateCoolNum)
    ! update nodes
    CALL UpdatePlateCoolingHeatExchanger(PlateCoolNum)
    ! update report variables
    CALL ReportPlateCoolingHeatExchanger(PlateCoolNum)

    RETURN

  END SUBROUTINE SimPlateCoolingHeatExchanger

  !==============================================================================

  SUBROUTINE GetPlateCoolingHeatExchanger

    ! SUBROUTINE INFORMATION:
    !       AUTHOR         Sankaranarayanan K P
    !       DATE WRITTEN   Nov 2003
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
    USE InputProcessor,        ONLY : GetNumObjectsFound, GetObjectItem, FindItemInList, FindItem, SameString
    USE DataIPShortCuts
    USE NodeInputManager,      ONLY : GetOnlySingleNode
    USE ScheduleManager,       ONLY : GetScheduleIndex
    USE BranchNodeConnections, ONLY : TestCompSet
    USE General,               ONLY : RoundSigDigits
    USE PlantUtilities,        ONLY : RegisterPlantCompDesignFlow

    IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

    ! SUBROUTINE PARAMETER DEFINITIONS:
    CHARACTER(len=*), PARAMETER :: Blank = ' '

    ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

    INTEGER                        :: IOStatus   ! Used in GetObjectItem
    INTEGER                        :: Item       ! Item to be "gotten"
    INTEGER                        :: NumAlphas  ! Number of Alphas for each GetObjectItem call
    INTEGER                        :: NumNumbers ! Number of Numbers for each GetObjectItem call
    LOGICAL                        :: ErrorsFound=.false.  ! Set to true if errors in input, fatal at end of routine

    ! Initializations and allocations
    cCurrentModuleObject = 'HeatExchanger:Plate'
    NumOfFreeCools = GetNumObjectsFound(TRIM(cCurrentModuleObject))


    ALLOCATE(PlateCool(NumOfFreeCools))
    ALLOCATE(CheckEquipName(NumOfFreeCools))
    CheckEquipName=.true.


    ! Obtain all of the user data
    DO Item = 1, NumOfFreeCools

      CALL GetObjectItem(TRIM(cCurrentModuleObject),Item,cAlphaArgs,NumAlphas,rNumericArgs,NumNumbers,IOStatus, &
      AlphaFieldnames=cAlphaFieldNames,NumericFieldNames=cNumericFieldNames)

      ! General user input data
      PlateCool(Item)%Name = cAlphaArgs(1)
      PlateCool(Item)%ComponentName = cAlphaArgs(2)
      PlateCool(Item)%ComponentType = cAlphaArgs(3)

      ! get Demand side inlet node data
      PlateCool(Item)%DemandSideInletNode = cAlphaArgs(4)
      PlateCool(Item)%DemandSideInletNodeNum  =   &
      GetOnlySingleNode(cAlphaArgs(4),ErrorsFound,TRIM(cCurrentModuleObject),cAlphaArgs(1),  &
      NodeType_Water,NodeConnectionType_Inlet,1,ObjectIsNotParent)
      IF (PlateCool(Item)%DemandSideInletNodeNum == 0) THEN
        CALL ShowSevereError('Invalid '//TRIM(cAlphaFieldNames(4))//'='//TRIM(cAlphaArgs(4)) )
        CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)) )
        ErrorsFound=.true.
      END IF

      ! get Demand side outlet node data
      PlateCool(Item)%DemandSideOutletNode = cAlphaArgs(5)
      PlateCool(Item)%DemandSideOutletNodeNum  =   &
      GetOnlySingleNode(cAlphaArgs(5),ErrorsFound,TRIM(cCurrentModuleObject),cAlphaArgs(1),  &
      NodeType_Water,NodeConnectionType_Outlet,1,ObjectIsNotParent)
      IF (PlateCool(Item)%DemandSideOutletNodeNum == 0) THEN
        CALL ShowSevereError('Invalid '//TRIM(cAlphaFieldNames(5))//'='//TRIM(cAlphaArgs(5)) )
        CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)) )
        ErrorsFound=.true.
      END IF

      CALL TestCompSet(TRIM(cCurrentModuleObject),cAlphaArgs(1),cAlphaArgs(4),cAlphaArgs(5),'Condenser Water Nodes')


      ! get Supply side inlet node data
      PlateCool(Item)%SupplySideInletNode = cAlphaArgs(6)
      PlateCool(Item)%SupplySideInletNodeNum  =   &
      GetOnlySingleNode(cAlphaArgs(6),ErrorsFound,TRIM(cCurrentModuleObject),cAlphaArgs(1),  &
      NodeType_Water,NodeConnectionType_Inlet,2,ObjectIsNotParent)
      IF (PlateCool(Item)%SupplySideInletNodeNum == 0) THEN
        CALL ShowSevereError('Invalid '//TRIM(cAlphaFieldNames(6))//'='//TRIM(cAlphaArgs(6)) )
        CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)) )
        ErrorsFound=.true.
      END IF

      ! get Supply side outlet node data
      PlateCool(Item)%SupplySideOutletNode = cAlphaArgs(7)
      PlateCool(Item)%SupplySideOutletNodeNum  =   &
      GetOnlySingleNode(cAlphaArgs(7),ErrorsFound,TRIM(cCurrentModuleObject),cAlphaArgs(1),  &
      NodeType_Water,NodeConnectionType_Outlet,2,ObjectIsNotParent)
      IF (PlateCool(Item)%SupplySideOutletNodeNum == 0) THEN
        CALL ShowSevereError('Invalid '//TRIM(cAlphaFieldNames(7))//'='//TRIM(cAlphaArgs(7)) )
        CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)) )
        ErrorsFound=.true.
      END IF

      CALL TestCompSet(TRIM(cCurrentModuleObject),cAlphaArgs(1),cAlphaArgs(6),cAlphaArgs(7),'Water Nodes')

      ! get heat exchange mode
      PlateCool(Item)%HeatXMode = cAlphaArgs(8)
      IF (.NOT. SameString(cAlphaArgs(8),'Ideal') .AND. &
      .NOT. SameString(cAlphaArgs(8),'UFactorTimesAreaEffectiveness')) THEN
      CALL ShowSevereError('Invalid '//TRIM(cAlphaFieldNames(8))//'='//TRIM(cAlphaArgs(8)) )
      CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)) )
      ErrorsFound=.true.
    END IF

    ! UA effectiveness data
    PlateCool(Item)%UA = rNumericArgs(1)

    ! Flow Rate data
    PlateCool(Item)%DemandSideFlowRate = rNumericArgs(2)
    PlateCool(Item)%SupplySideFlowRate = rNumericArgs(3)

    IF (rNumericArgs(1) == 0.0) THEN
      CALL ShowSevereError('Invalid '//TRIM(cNumericFieldNames(1))//'='//TRIM(RoundSigDigits(rNumericArgs(1),2)))
      CALL ShowContinueError('Entered in '//TRIM(cCurrentModuleObject)//'='//TRIM(cAlphaArgs(1)))
      CALL ShowContinueError('Value must be greater than 0.0')
      ErrorsFound=.true.
    END IF

  END DO  ! end of input loop

  ! Set up the output variables
  DO Item = 1, NumOfFreeCools

    ! heat ransfer rates
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Heat Transfer Rate[W]',    &
    PlateCool(Item)%HeatTransRate,'System','Average', &
    PlateCool(Item)%Name)
    ! flow rates
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Demand Side Mass Flow rate[kg/s]',      &
    PlateCool(Item)%DemandSideMassFlowRate,'System','Average', &
    PlateCool(Item)%Name)
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Supply Side Mass Flow rate[kg/s]',      &
    PlateCool(Item)%SupplySideMassFlowRate,'System','Average', &
    PlateCool(Item)%Name)
    ! Demand side temps
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Demand Side Inlet Temp[C]',     &
    PlateCool(Item)%DemandSideInletTemp,'System','Average', &
    PlateCool(Item)%Name)
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Demand Side Outlet Temp[C]',     &
    PlateCool(Item)%DemandSideOutletTemp,'System','Average', &
    PlateCool(Item)%Name)
    ! Supply side temps
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Supply Side Inlet Temp[C]',     &
    PlateCool(Item)%SupplySideInletTemp,'System','Average', &
    PlateCool(Item)%Name)
    CALL SetupOutputVariable('Plate Cooling Heat Exchanger Supply Side Outlet Temp[C]',     &
    PlateCool(Item)%SupplySideOutletTemp,'System','Average', &
    PlateCool(Item)%Name)
  END DO

  ! save the design demand side flow rate for use by plant loop sizing algorithms
  DO Item = 1, NumOfFreeCools
    CALL RegisterPlantCompDesignFlow(PlateCool(Item)%DemandSideInletNodeNum,PlateCool(Item)%DemandSideFlowRate)
  END DO

  ! Issue fatal if errors were encountered
  IF (ErrorsFound) THEN
    CALL ShowFatalError('Errors found in processing input for '//TRIM(cCurrentModuleObject) )
  END IF

  RETURN

END SUBROUTINE GetPlateCoolingHeatExchanger

!==============================================================================

SUBROUTINE InitPlateCoolingHeatExchanger(PlateCoolNum)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Sankaranarayanan K P
  !       DATE WRITTEN   Nov 2003
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine just sets some local variables used to store
  ! node numbers for current component.

  ! METHODOLOGY EMPLOYED:
  ! Reads free cooling HX component data structure

  ! USE STATEMENTS:
  USE DataLoopNode,    ONLY : Node
  USE DataPlant,       ONLY : TypeOf_HtExchgPlateFreeClng, ScanPlantLoopsForObject, &
  PlantLoop
  USE PlantUtilities,  ONLY : InterConnectTwoPlantLoopSides, InitComponentNodes
  USE FluidProperties, ONLY : GetDensityGlycol


  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  INTEGER, INTENT(IN) :: PlateCoolNum             ! Index for the free cooling heat exchanger.

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  INTEGER      :: Num                     ! counter
  REAL(r64)    :: DesignDemandSideMassFlowRate  ! Mass flow rate used to initialize condenser nodes
  REAL(r64)    :: DesignSupplySideMassFlowRate
  LOGICAL,SAVE :: MyEnvrnFlag = .TRUE.
  LOGICAL      :: errFlag
  REAL(r64)    :: rho ! local fluid density


  PlateCool(PlateCoolNum)%HeatTransRate = 0.0

  ! Init variables
  IF (PlateCool(PlateCoolNum)%MyFlag) THEN
    ! Locate the hx on the plant loops for later usage
    errFlag=.false.
    CALL ScanPlantLoopsForObject(PlateCool(PlateCoolNum)%Name, &
    TypeOf_HtExchgPlateFreeClng, &
      PlateCool(PlateCoolNum)%CondLoopNum, &
      PlateCool(PlateCoolNum)%CondLoopSideNum, &
      PlateCool(PlateCoolNum)%CondBranchNum, &
      PlateCool(PlateCoolNum)%CondCompNum, &
      InletNodeNumber = PlateCool(PlateCoolNum)%DemandSideInletNodeNum,  &
      errFlag=errFlag)
      CALL ScanPlantLoopsForObject(PlateCool(PlateCoolNum)%Name, &
      TypeOf_HtExchgPlateFreeClng, &
        PlateCool(PlateCoolNum)%PlantLoopNum, &
        PlateCool(PlateCoolNum)%PlantLoopSideNum, &
        PlateCool(PlateCoolNum)%PlantBranchNum, &
        PlateCool(PlateCoolNum)%PlantCompNum, &
        InletNodeNumber = PlateCool(PlateCoolNum)%SupplySideInletNodeNum,  &
        errFlag=errFlag)
        CALL InterConnectTwoPlantLoopSides( PlateCool(PlateCoolNum)%CondLoopNum,      &
        PlateCool(PlateCoolNum)%CondLoopSideNum,  &
        PlateCool(PlateCoolNum)%PlantLoopNum,     &
        PlateCool(PlateCoolNum)%PlantLoopSideNum, &
        TypeOf_HtExchgPlateFreeClng , .FALSE.)

          IF (errFlag) THEN
            CALL ShowFatalError('InitPlateCoolingHeatExchanger: Program terminated due to previous condition(s).')
          ENDIF
          PlateCool(PlateCoolNum)%MyFlag=.FALSE.
        ENDIF

        IF (BeginEnvrnFlag .AND. MyEnvrnFlag) THEN

          PlateCool(PlateCoolNum)%DemandSideInletTemp      =0.0
          PlateCool(PlateCoolNum)%DemandSideOutletTemp     =0.0
          PlateCool(PlateCoolNum)%SupplySideInletTemp      =0.0
          PlateCool(PlateCoolNum)%SupplySideOutletTemp     =0.0
          PlateCool(PlateCoolNum)%DemandSideMassFlowRate   =0.0
          PlateCool(PlateCoolNum)%SupplySideMassFlowRate   =0.0

          SupplySideMassFlowRate   =0.0
          DemandSideMassFlowRate   =0.0
          SupplySideInletTemp      =0.0
          DemandSideInletTemp      =0.0

          DO Num = 1, NumOfFreeCools

            rho = GetDensityGlycol(PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidName, &
            InitConvTemp, &
            PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidIndex, &
            'InitPlateCoolingHeatExchanger')

            DesignDemandSideMassFlowRate = PlateCool(Num)%DemandSideFlowRate * rho
            CALL InitComponentNodes( &
            0.0d0, &
            DesignDemandSideMassFlowRate, &
            PlateCool(Num)%DemandSideInletNodeNum, &
            PlateCool(Num)%DemandSideOutletNodeNum, &
            PlateCool(Num)%CondLoopNum, &
            PlateCool(Num)%CondLoopSideNum, &
            PlateCool(Num)%CondBranchNum, &
            PlateCool(Num)%CondCompNum &
            )

            rho = GetDensityGlycol(PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidName, &
            InitConvTemp, &
            PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidIndex, &
            'InitPlateCoolingHeatExchanger')

            DesignSupplySideMassFlowRate = PlateCool(Num)%SupplySideFlowRate * rho
            CALL InitComponentNodes( &
            0.0d0, &
            DesignSupplySideMassFlowRate, &
            PlateCool(Num)%SupplySideInletNodeNum, &
            PlateCool(Num)%SupplySideOutletNodeNum, &
            PlateCool(Num)%PlantLoopNum, &
            PlateCool(Num)%PlantLoopSideNum, &
            PlateCool(Num)%PlantBranchNum, &
            PlateCool(Num)%PlantCompNum &
            )

          END DO

          MyEnvrnFlag = .FALSE.

        END IF

        IF (.not. BeginEnvrnFlag) THEN
          MyEnvrnFlag=.true.
        ENDIF

        !  set module variables for node numbers for this heat exchanger
        DemandSideInletNodeNum  = PlateCool(PlateCoolNum)%DemandSideInletNodeNum
        DemandSideOutletNodeNum = PlateCool(PlateCoolNum)%DemandSideOutletNodeNum
        SupplySideInletNodeNum  = PlateCool(PlateCoolNum)%SupplySideInletNodeNum
        SupplySideOutletNodeNum = PlateCool(PlateCoolNum)%SupplySideOutletNodeNum

      END SUBROUTINE InitPlateCoolingHeatExchanger

      !==============================================================================

      SUBROUTINE CalcPlateCoolingHeatExchanger(PlateCoolNum)

        !       AUTHOR         Sankaranarayanan K P
        !       DATE WRITTEN   Nov 2003
        !       MODIFIED       na
        !       RE-ENGINEERED  na

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
        USE FluidProperties, ONLY  : GetSpecificHeatGlycol, GetDensityGlycol
        USE DataPlant,       ONLY  : PlantLoop
        USE DataBranchAirLoopPlant, ONLY  : MassFlowTolerance
        USE PlantUtilities,  ONLY  : SetComponentFlowRate

        IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

        ! SUBROUTINE ARGUMENT DEFINITIONS:
        INTEGER, INTENT(IN)  :: PlateCoolNum          ! Index for the free cooling heat exchanger.

        ! SUBROUTINE PARAMETER DEFINITIONS:
        CHARACTER(len=*), PARAMETER :: CalledFrom='CalcPlateCoolingHeatExchanger'

        ! INTERFACE BLOCK SPECIFICATIONS
        ! na

        ! DERIVED TYPE DEFINITIONS
        ! na

        ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        REAL(r64)    :: SupplySideFluidCp        ! Specific heat of Plant side fluid
        REAL(r64)    :: SupplySideCapRate        ! Capacity rate (mdot*Cp) of SupplySide side fluid
        REAL(r64)    :: DemandSideFluidCp         ! Specific heat of condenser side fluid
        REAL(r64)    :: SupplySideInletdensity   ! density on SupplySide side
        REAL(r64)    :: DemandSideInletdensity    ! density on condenser side
        REAL(r64)    :: DemandSideCapRate         ! Capacity rate (mdot*Cp) of condenser side fluid
        REAL(r64)    :: MinCapRate          ! minimum capacity rate
        REAL(r64)    :: CapRatio            ! capacity ratio (min/max)
        REAL(r64)    :: Effectiveness       ! heat exchanger effectiveness
        REAL(r64)    :: NTU                 ! dimensionless NTU calculated from UA
        LOGICAL      :: ItemNotFound        ! error flag
        INTEGER      :: LoopNum
        INTEGER      :: LoopSideNum

        ItemNotFound = .FALSE.
        SupplySideInletdensity = GetDensityGlycol(PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidName, &
        Node(SupplySideInletNodeNum)%Temp, &
        PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidIndex, &
        'CalcPlateCoolingHeatExchanger')

        DemandSideInletdensity  = GetDensityGlycol(PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidName, &
        Node(DemandSideInletNodeNum)%Temp, &
        PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidIndex, &
        'CalcPlateCoolingHeatExchanger')
        LoopNum = PlateCool(PlateCoolNum)%PlantLoopNum
        LoopSideNum = PlateCool(PlateCoolNum)%PlantLoopSideNum

        ! Calculate design supply (plant) flow rate, and ask for this
        SupplySideMassFlowRate = PlateCool(PlateCoolNum)%SupplySideFlowRate*SupplySideInletdensity
        CALL SetComponentFlowRate( &
        SupplySideMassFlowRate, &
        PlateCool(PlateCoolNum)%SupplySideInletNodeNum, &
        PlateCool(PlateCoolNum)%SupplySideOutletNodeNum, &
        PlateCool(PlateCoolNum)%PlantLoopNum, &
        PlateCool(PlateCoolNum)%PlantLoopSideNum, &
        PlateCool(PlateCoolNum)%PlantBranchNum, &
        PlateCool(PlateCoolNum)%PlantCompNum &
        )

        ! Calculate the design demand (condenser) flow rate, but try to turn off if we don't need it to run
        DemandSideMassFlowRate = PlateCool(PlateCoolNum)%DemandSideFlowRate*DemandSideInletdensity
        IF(SupplySideMassFlowRate < MassFlowTolerance) DemandSideMassFlowRate = 0.0
        CALL SetComponentFlowRate( &
        DemandSideMassFlowRate, &
        PlateCool(PlateCoolNum)%DemandSideInletNodeNum, &
        PlateCool(PlateCoolNum)%DemandSideOutletNodeNum, &
        PlateCool(PlateCoolNum)%CondLoopNum, &
        PlateCool(PlateCoolNum)%CondLoopSideNum, &
        PlateCool(PlateCoolNum)%CondBranchNum, &
        PlateCool(PlateCoolNum)%CondCompNum &
        )

        !set local variables for heat exchanger calculation
        DemandSideInletTemp = Node(DemandSideInletNodeNum)%Temp
        SupplySideInletTemp = Node(SupplySideInletNodeNum)%Temp

        DemandSideFluidCp   = GetSpecificHeatGlycol(PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidName,DemandSideInletTemp, &
        PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidIndex,CalledFrom)
        SupplySideFluidCp   = GetSpecificHeatGlycol(PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidName, SupplySideInletTemp, &
        PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidIndex,CalledFrom)

        DemandSideCapRate  = DemandSideFluidCp * DemandSideMassFlowRate
        SupplySideCapRate = SupplySideFluidCp * SupplySideMassFlowRate
        MinCapRate = MIN(DemandSideCapRate, SupplySideCapRate)

        !If there is no flow rate on either the condenser or plant side, turn off heat exchanger and return
        IF (DemandSideCapRate <= zerocaptol .OR. SupplySideCapRate <= zerocaptol)THEN
          !    PlateCool(PlateCoolNum)%HeatTransRate = 0.0
          RETURN
        END IF

        ! calc effectiveness - 1.0 if in ideal mode
        IF(PlateCool(PlateCoolNum)%HeatXMode == 'UFACTORTIMESAREAEFFECTIVENESS')THEN
          ! assume cross flow, both unmixed
          CapRatio = MinCapRate/MAX(DemandSideCapRate, SupplySideCapRate)
          NTU = PlateCool(PlateCoolNum)%UA/MinCapRate

          Effectiveness = 1.0d0 - EXP( (NTU**0.22d0/CapRatio) * &
          (EXP(-CapRatio*NTU**0.78d0) - 1.0d0) )
        ELSE
          ! must be in ideal mode
          Effectiveness = 1.0
        END IF
        ! overall heat transfer rate
        ! convention is +ve rate is rejected from plant side to condenser side
        PlateCool(PlateCoolNum)%HeatTransRate = Effectiveness*MinCapRate*(SupplySideInletTemp-DemandSideInletTemp)

        RETURN

      END SUBROUTINE CalcPlateCoolingHeatExchanger

      !==============================================================================

      SUBROUTINE TurnOffPHXIfNotNeeded(PlateCoolNum)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Edwin Lee
        !       DATE WRITTEN   September 2009
        !       MODIFIED       na
        !       RE-ENGINEERED  na

        ! PURPOSE OF THIS SUBROUTINE:
        ! This subroutine is triggered via the demand side solver
        ! The PHX would tend to keep a flow request on the demand side loop, even
        !  if the supply side didn't want any flow, which would keep the pump running
        ! This routine will tell the demand side loop that it doesn't want flow
        !  if the supply side of the PHX doesn't have any flow

        ! METHODOLOGY EMPLOYED:
        ! Check supply node flow rates, turn off demand side if no flow

        ! NOTE:
        ! This fix does work, but it is not ideal.  The complete rewrite of the demand
        !  side solver, which is underway now, will fix this problem without these little
        !  band-aids.  However, this was placed in for the v4.0 release as a short term fix,
        !  by recommendation of BG.

        ! REFERENCES:
        ! na

        ! USE STATEMENTS:
        USE DataLoopNode,    ONLY : Node
        USE DataBranchAirLoopPlant, ONLY : MassFlowTolerance
        USE PlantUtilities,  ONLY : SetComponentFlowRate

        IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

        ! SUBROUTINE ARGUMENT DEFINITIONS:
        INTEGER, INTENT(IN) :: PlateCoolNum  ! Index for the free cooling heat exchanger.

        ! SUBROUTINE PARAMETER DEFINITIONS:
        ! na

        ! INTERFACE BLOCK SPECIFICATIONS
        ! na

        ! DERIVED TYPE DEFINITIONS
        ! na

        ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        REAL(r64) :: DemandSideMassFlowRate

        !check supply side for flow, turn off if not needed
        IF (Node(PlateCool(PlateCoolNum)%SupplySideInletNodeNum)%MassFlowRate < MassFlowTolerance) THEN
          DemandSideMassFlowRate = 0.0
          CALL SetComponentFlowRate( &
          DemandSideMassFlowRate, &
          PlateCool(PlateCoolNum)%DemandSideInletNodeNum, &
          PlateCool(PlateCoolNum)%DemandSideOutletNodeNum, &
          PlateCool(PlateCoolNum)%CondLoopNum, &
          PlateCool(PlateCoolNum)%CondLoopSideNum, &
          PlateCool(PlateCoolNum)%CondBranchNum, &
          PlateCool(PlateCoolNum)%CondCompNum &
          )
        END IF

        RETURN

      END SUBROUTINE


      SUBROUTINE UpdatePlateCoolingHeatExchanger(PlateCoolNum)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Sankaranarayanan K P
        !       DATE WRITTEN   Nov 2003
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
        USE DataPlant,       ONLY : PlantLoop

        IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

        ! SUBROUTINE ARGUMENT DEFINITIONS:
        INTEGER, INTENT(IN) :: PlateCoolNum  ! Index for the free cooling heat exchanger.

        ! SUBROUTINE PARAMETER DEFINITIONS:
        CHARACTER(len=*), PARAMETER :: CalledFrom='UpdatePlateCoolingHeatExchanger'

        ! INTERFACE BLOCK SPECIFICATIONS
        ! na

        ! DERIVED TYPE DEFINITIONS
        ! na

        ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
        REAL(r64)    :: DemandSideFluidCp         ! Specific heat of condenser side fluid
        REAL(r64)    :: SupplySideFluidCp        ! Specific heat of SupplySide side fluid


        ! flow rates
        !   Node(DemandSideInletNodeNum)%MassFlowRate = DemandSideMassFlowRate
        !   Node(SupplySideInletNodeNum)%MassFlowRate = SupplySideMassFlowRate

        ! specific heats
        DemandSideFluidCp   = GetSpecificHeatGlycol(PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidName,DemandSideInletTemp, &
        PlantLoop(PlateCool(PlateCoolNum)%CondLoopNum)%FluidIndex,CalledFrom)
        SupplySideFluidCp   = GetSpecificHeatGlycol(PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidName, SupplySideInletTemp, &
        PlantLoop(PlateCool(PlateCoolNum)%PlantLoopNum)%FluidIndex,CalledFrom)

        ! always mass through flow rates
        !  Node(DemandSideOutletNodeNum)%MassFlowRate          = Node(DemandSideInletNodeNum)%MassFlowRate
        !  Node(DemandSideOutletNodeNum)%MassFlowRateMaxAvail  = Node(DemandSideInletNodeNum)%MassFlowRateMaxAvail
        !  Node(DemandSideOutletNodeNum)%MassFlowRateMinAvail  = Node(DemandSideInletNodeNum)%MassFlowRateMinAvail

        !  Node(SupplySideOutletNodeNum)%MassFlowRate         = Node(SupplySideInletNodeNum)%MassFlowRate
        !  Node(SupplySideOutletNodeNum)%MassFlowRateMaxAvail = Node(SupplySideInletNodeNum)%MassFlowRateMaxAvail
        !  Node(SupplySideOutletNodeNum)%MassFlowRateMinAvail = Node(SupplySideInletNodeNum)%MassFlowRateMinAvail

        ! check if coupled or zero heat transfer rate
        IF(PlateCool(PlateCoolNum)%HeatTransRate /= 0.0) THEN
          ! calc outlet temps from heat transfer rate
          Node(DemandSideOutletNodeNum)%Temp  = Node(DemandSideInletNodeNum)%Temp + PlateCool(PlateCoolNum)%HeatTransRate/ &
          (DemandSideMassFlowRate * DemandSideFluidCp)

          Node(SupplySideOutletNodeNum)%Temp = Node(SupplySideInletNodeNum)%Temp - PlateCool(PlateCoolNum)%HeatTransRate/ &
          (SupplySideMassFlowRate * SupplySideFluidCp)
        ELSE
          ! just pass through
          Node(DemandSideOutletNodeNum)%Temp  = Node(DemandSideInletNodeNum)%Temp
          Node(SupplySideOutletNodeNum)%Temp = Node(SupplySideInletNodeNum)%Temp
        END IF

        RETURN


      END SUBROUTINE UpdatePlateCoolingHeatExchanger

      !==============================================================================

      SUBROUTINE ReportPlateCoolingHeatExchanger(PlateCoolNum)

        ! SUBROUTINE INFORMATION:
        !       AUTHOR         Sankaranarayanan K P
        !       DATE WRITTEN   Nov 2003
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

        IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

        ! SUBROUTINE ARGUMENT DEFINITIONS:
        INTEGER, INTENT(IN) :: PlateCoolNum  ! Index for the free cooling heat exchanger.

        ! SUBROUTINE PARAMETER DEFINITIONS:

        ! INTERFACE BLOCK SPECIFICATIONS
        ! na

        ! DERIVED TYPE DEFINITIONS
        ! na

        ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:

        ! update condenser side variables
        PlateCool(PlateCoolNum)%DemandSideInletTemp     = Node(DemandSideInletNodeNum)%Temp
        PlateCool(PlateCoolNum)%DemandSideOutletTemp    = Node(DemandSideOutletNodeNum)%Temp
        PlateCool(PlateCoolNum)%DemandSideMassFlowRate  = Node(DemandSideInletNodeNum)%MassFlowRate

        ! update plant side variables
        PlateCool(PlateCoolNum)%SupplySideInletTemp    = Node(SupplySideInletNodeNum)%Temp
        PlateCool(PlateCoolNum)%SupplySideOutletTemp   = Node(SupplySideOutletNodeNum)%Temp
        PlateCool(PlateCoolNum)%SupplySideMassFlowRate = Node(SupplySideInletNodeNum)%MassFlowRate

        RETURN

      END SUBROUTINE ReportPlateCoolingHeatExchanger

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

    END MODULE PlateCoolingHeatExchanger
