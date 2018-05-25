! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module contains a subroutine to bring in the necessary inputs for a heat pump simulation.

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component has no physical representation.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module merely reads in the inputs, and puts them into property arrays as needed.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! This reads in from an already open input file---see "InputProcessor" for file information.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! Error and Manufacturer flags are the only things defined at the module level.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method:
!    PUBLIC GetInputs -- Reads in all the input data for the HP Simulation
!      Called by ORNLsolver

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! NA

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-12 | RAS | Updated header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! The arrays, particularly EvaporatorPAR and CondenserPAR, should be cleaned up.
! The local variables should probably also be cleaned up.

!***********************************************************************************

MODULE HeatPumpInput

  USE DataSimulation
  USE DataGlobals_HPSimIntegrated, ONLY: RefName    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  implicit none

  PRIVATE

  !Error Flags
  INTEGER,PARAMETER :: NOERROR       = 0
  INTEGER,PARAMETER :: CONVERGEERROR = 1
  INTEGER,PARAMETER :: REFPROPERROR  = 2
  INTEGER,PARAMETER :: CKTFILEERROR  = 3
  INTEGER,PARAMETER :: COILTUBEERROR = 4
  INTEGER,PARAMETER :: COILFINERROR  = 5
  INTEGER,PARAMETER :: AIRSIDEERROR  = 6
  INTEGER,PARAMETER :: ZEROLENCOILERROR  = 7
  INTEGER,PARAMETER :: DPERROR       = 8

  !Compressor Manufacturer
  INTEGER,PARAMETER :: COPELAND  = 1 !ISI - 10/05/06
  INTEGER,PARAMETER :: BRISTOL   = 2
  INTEGER,PARAMETER :: DANFOSS   = 3
  INTEGER,PARAMETER :: PANASONIC = 4

  !Unit flags !ISI - 07/14/06
  INTEGER,PARAMETER :: SI=1
  INTEGER,PARAMETER :: IP=2

  PUBLIC GetInputs
  PUBLIC GetGenOptInputs

CONTAINS

  SUBROUTINE GetInputs

    ! ----------------------------------------------------------------------
    !
    !   Description: To collect all the input data for
    !                steady-state heat pump simulation
    !
    !   Author:
    !   Cavlin Iu
    !   Mechanical and Aerospace Engineering
    !   Oklahoma State University, Stillwater
    !
    !   Date:
    !   November, 2002
    !
    ! ----------------------------------------------------------------------

    USE DataStopCodes
    USE InputProcessor
    USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength, RefName !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
    USE DataSimulation, ONLY: IsCoolingMode !RS: Debugging: Global variable now

    IMPLICIT NONE

    !Local variables
    REAL CoolingShTbPAR(5)		!Cooling mode short tube model input data
    REAL HeatingShTbPAR(5)		!Heating mode short tube model input data
    REAL CoolingCapTubePAR(5)	!Cooling mode cap. tube model input data
    REAL HeatingCapTubePAR(5)	!Heating mode cap. tube model input data
    REAL CoolingDistubeLength	!Distributor tube length for cooling, mm or ft
    REAL HeatingDistubeLength	!Distributor tube length for heating, mm or ft
    REAL SucLnPAR(7)			!Suction line parameter
    REAL DisLnPAR(7)			!Discharge line parameter
    REAL LiqLnPAR(7)			!Liquid line parameter
    REAL ValveIDCLnPAR(7)		!Valve to indoor coil line parameter
    REAL ValveODCLnPAR(7)		!Valve to outdoor coil line parameter
    REAL :: EqLength=0.0				!Equivalent length, m or ft
    REAL :: EqDiameter=0.0				!Equivalent diameter, mm or in
    REAL EqThickness			!Equivalent thickenss, mm or mil
    REAL TotElevation			!Total elevation, m or ft
    REAL TotHeatGain			!Total heat gain, w or Btu/hr
    REAL TotTempChange			!Total temperature change, C or F
    REAL TotAddDP				!Total additional pressure drop, kPa or psia
    REAL VolSucLn				!Suction line volume, m^3 or ft^3
    REAL VolDisLn				!Discharge line volume, m^3 or ft^3
    REAL VolValveIDCLn			!Valve to IDC line volume, m^3 or ft^3
    REAL VolValveODCLn			!Valve to ODC line volume, m^3 or ft^3
    REAL TotVolume				!Total line volume, m^3 or ft^3
    REAL,PARAMETER :: PI=3.14159265 !Pi
    REAL,PARAMETER :: CopperDensity=8920 !Density of copper, kg/m3

    CHARACTER*150 LineData

    REAL :: CopperVol !Copper volume, m3
    INTEGER(2) :: CoolingExpDevice !Cooling Expansion device: 1=short tube; 2=TXV; 3=Cap. tube
    INTEGER(2) :: HeatingExpDevice !Heating Expansion device: 1=short tube; 2=TXV; 3=Cap. tube
    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"

    REAL RefChg    !Design Refrigerant Charge Mass

    CHARACTER(len=MaxNameLength)SucLn_RefrigerantLine
    CHARACTER(len=MaxNameLength)SucLn_TubeType
    REAL SucLn_KTube    !Conductivity of Suction Line Tube
    REAL SucLn_TubeID   !Inner Diameter of Suction Line Tube
    REAL SucLn_Charge   !Suction Line Charge
    CHARACTER(len=MaxNameLength)DisLn_RefrigerantLine
    CHARACTER(len=MaxNameLength)DisLn_TubeType
    REAL DisLn_KTube    !Conductivity of Discharge Line Tube
    REAL DisLn_TubeID   !Inner Diameter of Discharge Line Tube
    REAL DisLn_Charge   !Discharge Line Charge
    CHARACTER(len=MaxNameLength)LiqLn_RefrigerantLine
    CHARACTER(len=MaxNameLength)LiqLn_TubeType
    REAL LiqLn_KTube    !Conductivity of Liquid Line Tube
    REAL LiqLn_TubeID   !Inner Diameter of Liquid Line Tube
    REAL LiqLn_Charge   !Liquid Line Charge
    CHARACTER(len=MaxNameLength)ValveIDCLn_RefrigerantLine
    CHARACTER(len=MaxNameLength)ValveIDCLn_TubeType
    REAL ValveIDCLn_KTube    !Conductivity of Valve to IDC Line Tube
    REAL ValveIDCLn_TubeID   !Inner Diameter of Valve to IDC Line Tube
    REAL ValveIDCLn_Charge   !Charge of Valve to IDC Line Tube
    CHARACTER(len=MaxNameLength)ValveODCLn_RefrigerantLine
    CHARACTER(len=MaxNameLength)ValveODCLn_TubeType
    REAL ValveODCLn_KTube    !Conductivity of Valve to ODC Line Tube
    REAL ValveODCLn_TubeID   !Inner Diameter of Valve to ODC Line Tube
    INTEGER CompressorManufacturer

    !Compressor Manufacturer
    INTEGER,PARAMETER :: COPELAND  = 1 !ISI - 10/05/06
    INTEGER,PARAMETER :: BRISTOL   = 2
    INTEGER,PARAMETER :: DANFOSS   = 3
    INTEGER,PARAMETER :: PANASONIC = 4

    !Flow:
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)


    !***************** System data *****************  !RS: Debugging: Moving: Stay here

    CALL GetObjectItem('MainDesignData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    SELECT CASE (Alphas(1)(1:1))
    CASE ('F','f')
      Unit = IP
    CASE ('T','t')
      Unit = SI
    CASE DEFAULT
      !FAIL
    END SELECT

    !Calculation mode
    SELECT CASE(TRIM(MakeUPPERcase(Alphas(2))))
    CASE('TXVSIMULATION')
      MODE=TXVSIMULATION
    CASE('FIXEDORIFICESIMULATION')
      MODE=FIXEDORIFICESIM
    CASE('FIXEDSUPERHEATSIMULATION')
      MODE=FIXEDSUPERHEATSIM
    CASE DEFAULT
      !FAIL
    END SELECT

    SystemType = Numbers(1)

    SELECT CASE (Alphas(3)(1:1))  !RS: Debugging: This cooling/heating mode should depend on what's being brought in from EPlus
    CASE ('F','f')
      IsCoolingMode=0
    CASE DEFAULT
      IsCoolingMode=1
    END SELECT

    Ref$=Alphas(4)    !Refrigerant Name
    RefChg = Numbers(2)    !Design Refrigerant Charge Mass

    !***************** Compressor data *****************  !RS: Debugging: Moving: Stay here

    CALL GetObjectItem('CompressorData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    SELECT CASE (Alphas(2)(1:1))
    CASE ('C','c')
      CompressorManufacturer=COPELAND
    CASE ('B','b')
      CompressorManufacturer=BRISTOL
    CASE ('D','d')
      CompressorManufacturer=DANFOSS
    CASE ('P','p')
      CompressorManufacturer=PANASONIC
    CASE DEFAULT
      CompressorManufacturer=BRISTOL
    END SELECT

    !RS: Debugging: Moving array data up since useless data has been removed
    EvapPAR%EvapCompMan=CompressorManufacturer !EvapPAR(52)=CompressorManufacturer !ISI - 10/05/06
    CondPAR%CondCompMan=CompressorManufacturer    !RS: Debugging: Formerly CONDPAR(60)

    TsiCmp = Numbers(26) !UserSpecifiedRatingEvapTemperature
    TsoCmp = Numbers(27) !UserSpecifiedRatingCondTemperature
    Subcool = Numbers(28) !UserSpecifiedRatingSubcooling
    Super = Numbers(29) !UserSpecifiedRatingSuperheat


    !***************** Outdoor fan data ***************** !RS: Debugging: Moving: Evaporator & Condenser

    CALL GetObjectItem('OutdoorFanData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    !PwrODfan = Numbers(1) !Fan Power
    IF (IsCoolingMode .GT. 0) THEN    !Populating arrays
      CFMcnd = Numbers(2)    !Fan Air Flow Rate
    ELSE
      CFMevp = Numbers(2)
    END IF
    !ODdrawBlow = Numbers(3)   !Draw Through (1) or Blow Through (2)

    FanOut%MotorEff=Numbers(4)
    FanOut%FanEff=Numbers(5)
    FanOut%DeltaPress=Numbers(6)
    FanOut%MotInAirFrac=Numbers(7)


    !***************** Indoor fan data *****************  !RS: Debugging: Moving: Evaporator & Condenser

    CALL GetObjectItem('IndoorFanData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    !PwrIDfan = Numbers(1) !Fan Power
    IF (IsCoolingMode .GT. 0) THEN    !Populating arrays
      CFMevp = Numbers(2)    !Fan Air Flow Rate
    ELSE
      CFMcnd = Numbers(2)
    END IF
    !IDdrawBlow = Numbers(3)   !Draw Through or Blow Through

    FanOutE%MotorEff=Numbers(4)
    FanOutE%FanEff=Numbers(5)
    FanOutE%DeltaPress=Numbers(6)
    FanOutE%MotInAirFrac=Numbers(7)

    !***************** Expansion device data *****************    !RS: Debugging: Moving: Stay here

    CALL GetObjectItem('ExpansionDeviceData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    SELECT CASE (Alphas(1)(1:1))
    CASE ('C','c')
      CoolingExpDevice=3 !Cap. tube
    CASE ('T','t')
      CoolingExpDevice=2 !TXV
    CASE DEFAULT
      CoolingExpDevice=1 !Short tube orifice
    END SELECT

    SELECT CASE (Alphas(2)(1:1))
    CASE ('C','c')
      HeatingExpDevice=3 !Cap. tube
    CASE ('T','t')
      HeatingExpDevice=2 !TXV
    CASE DEFAULT
      HeatingExpDevice=1 !Short tube orifice
    END SELECT

    !Short tube orifice

    !Cooling mode

    CoolingShTbPAR(1) = Numbers(1)    !Length
    CoolingShTbPAR(2) = Numbers(2)    !Diameter
    CoolingShTbPAR(3) = Numbers(3)    !Chamfer Depth

    !Heating mode

    HeatingShTbPAR(1) = Numbers(4)    !Length
    HeatingShTbPAR(2) = Numbers(5)    !Diameter
    HeatingShTbPAR(3) = Numbers(6)    !Chamfer Depth

    !Capillary Tube

    !Cooling Mode

    CoolingCapTubePAR(2) = Numbers(7) !Length
    CoolingCapTubePAR(1) = Numbers(8) !Diameter
    CoolingCapTubePAR(3) = Numbers(9) !Coil Diameter

    !Heating Mode

    HeatingCapTubePAR(2) = Numbers(10)    !Length
    HeatingCapTubePAR(1) = Numbers(11)    !Diameter
    HeatingCapTubePAR(3) = Numbers(12)    !Coil Diameter

    !Distributor tubes

    CoolingDistubeLength = Numbers(13)

    IF (Unit .EQ. SI) THEN
      CoolingDistubeLength=CoolingDistubeLength/1000  !RS Comment: Unit Conversion
    ELSEIF (Unit .EQ. IP) THEN
      CoolingDistubeLength=CoolingDistubeLength/12    !RS Comment: Unit Conversion, from in to ft?
    END IF

    HeatingDistubeLength = Numbers(14)

    IF (Unit .EQ. SI) THEN
      HeatingDistubeLength=HeatingDistubeLength/1000  !RS Comment: Unit Conversion
    ELSEIF (Unit .EQ. IP) THEN
      HeatingDistubeLength=HeatingDistubeLength/12    !RS Comment: Unit Conversion, from in to ft?
    END IF

    !*****************Refrigerant line data****************** !RS: Debugging: Moving: Split into Evaporator & Condenser?

    CALL GetObjectItem('RefrigerantLineData',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    !Suction Line

    !SucLn_RefrigerantLine = Alphas(1)
    !SucLn_TubeType = Alphas(2)

    SucLnPAR(1) = Numbers(1)  !Refrigerant Line Length
    SucLnPAR(4) = Numbers(2)  !Refrigerant Line Elevation
    SucLnPAR(5) = Numbers(3)  !Refrigerant Line Heat Loss
    SucLnPAR(6) = Numbers(4)  !Refrigerant Line Temperature Change
    SucLnPAR(2) = Numbers(6) !5)  !Tube Outside Diameter
    !SucLn_KTube = Numbers(5)
    SucLn_TubeID = Numbers(5) !6)
    SucLnPAR(7) = Numbers(7)  !Additional Pressure Drop
    !SucLn_Charge = Numbers(9) !Charge in Line

    !Suction line tube wall thickness, mm or mil
    SucLnPAR(3)=(SucLnPAR(2)-SucLn_TubeID)/2
    IF (Unit .EQ. IP) THEN
      SucLnPAR(3)=SucLnPAR(3)*1000  !RS Comment: Unit Conversion
    END IF

    !*************** Accumulator **************** !RS: Debugging: Moving: ORNLSolver?
    AccumPAR%AccDTube=(SucLnPAR(2)-SucLnPAR(3)/1000*2) !J-tube diameter, mm or in   !RS: Debugging: Formerly AccumPAR(6)
    !********************************************

    !Discharge Line

    !DisLn_RefrigerantLine = Alphas(3)
    !DisLn_TubeType = Alphas(4)

    DisLnPAR(1) = Numbers(8) !(10) !Refrigerant Line Length
    DisLnPAR(4) = Numbers(9) !(11) !Refrigerant Line Elevation
    DisLnPAR(5) = Numbers(10) !2) !Refrigerant Line Heat Loss
    DisLnPAR(6) = Numbers(11) !3) !Refrigerant Line Temperature Change
    !DisLn_Ktube = Numbers(14)
    DisLn_TubeID = Numbers(12) !5)
    DisLnPAR(2) = Numbers(13) !6) !Tube Outside Diameter
    DisLnPAR(7) = Numbers(14) !7) !Additional Pressure Drop
    !DisLn_Charge = Numbers(18)    !Charge in Line

    !Discharge line tube wall thickness, mm or mil
    DisLnPAR(3)=(DisLnPAR(2)-DisLn_TubeID)/2
    IF (Unit .EQ. IP) THEN
      DisLnPAR(3)=DisLnPAR(3)*1000  !RS Comment: Unit Conversion
    END IF

    !Liquid Line

    !LiqLn_RefrigerantLine = Alphas(5)
    !LiqLn_TubeType = Alphas(6)

    LiqLnPAR(1) = Numbers(15) !9) !Refrigerant Line Length
    LiqLnPAR(4) = Numbers(16) !0) !Refrigerant Line Elevation
    LiqLnPAR(5) = Numbers(17) !21) !Refrigerant Line Heat Loss
    LiqLnPAR(6) = Numbers(18) !22) !Refrigerant Line Temperature Change
    !LiqLn_Ktube = Numbers(23) !Tube Conductivity
    LiqLn_TubeID = Numbers(19) !24)
    LiqLnPAR(2) = Numbers(20) !5) !Tube Outside Diameter
    LiqLnPAR(7) = Numbers(21) !6) !Additional Pressure Drop
    !LiqLn_Charge = Numbers(27)    !Charge in Line

    !Liquid line tube wall thickness, mm or mil
    LiqLnPAR(3)=(LiqLnPAR(2)-LiqLn_TubeID)/2
    IF (Unit .EQ. IP) THEN
      LiqLnPAR(3)=LiqLnPAR(3)*1000  !RS Comment: Unit Conversion
    END IF

    !Reversing Valve to IDC

    !ValveIDCLn_RefrigerantLine = Alphas(7)
    !ValveIDCLn_TubeType = Alphas(8)

    ValveIDCLnPAR(1) = Numbers(22) !8)    !Refrigerant Line Length
    ValveIDCLnPAR(4) = Numbers(23) !9)    !Refrigerant Line Elevation
    ValveIDCLnPAR(5) = Numbers(24) !30)    !Refrigerant Line Heat Loss
    ValveIDCLnPAR(6) = Numbers(25)! 31)    !Refrigerant Line Temperature Change
    !ValveIDCLn_Ktube = Numbers(32)
    ValveIDCLn_TubeID = Numbers(26) !33)
    ValveIDCLnPAR(2) = Numbers(27) !34)    !Tube Outside Diameter
    ValveIDCLnPAR(7) = Numbers(28) !35)    !Additional Pressure Drop
    !ValveIDCLn_Charge = Numbers(36)   !Charge in Line

    !Valve to IDC line tube wall thickness, mm or mil
    ValveIDCLnPAR(3)=(ValveIDCLnPAR(2)-ValveIDCLn_TubeID)/2
    IF (Unit .EQ. IP) THEN
      ValveIDCLnPAR(3)=ValveIDCLnPAR(3)*1000    !RS Comment: Unit Conversion
    END IF

    !Valve to ODC Line

    !ValveODCLn_RefrigerantLine = Alphas(9)
    !ValveODCLn_TubeType = Alphas(10)

    ValveODCLnPAR(1) = Numbers(29) !37)    !Refrigerant Line Length
    ValveODCLnPAR(4) = Numbers(30) !8)    !Refrigerant Line Elevation
    ValveODCLnPAR(5) = Numbers(31) !9)    !Refrigerant Line Heat Loss
    ValveODCLnPAR(6) = Numbers(32) !40)    !Refrigerant Line Temperature Change
    !ValveODCLn_Ktube = Numbers(41)
    ValveODCLn_TubeID = Numbers(33) !42)
    ValveODCLnPAR(2) = Numbers(34) !43)    !Tube Outside Diameter
    ValveODCLnPAR(7) = Numbers(35) !44)    !Additional Pressure Drop

    !Valve to ODC line tube wall thickness, mm or mil
    ValveODCLnPAR(3)=(ValveODCLnPAR(2)-ValveODCLn_TubeID)/2
    IF (Unit .EQ. IP) THEN
      ValveODCLnPAR(3)=ValveODCLnPAR(3)*1000    !RS Comment: Unit Conversion
    END IF

    !********************Refrigerant Cycle Data (Heating)***********************  !RS: Debugging: Moving: Stay here? Compressor? ORNLSolver?

    CALL GetObjectItem('RefrigerantCycleData(Heating)',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    !Indoor Coil Outlet

    BaroPressure = Numbers(1)  !Barometric Pressure
    !IsCmpInAirStream = Numbers(2) !Is Compressor in Air Stream

    !***** Calculate weight of interconnecting pipes **** !RS: Debugging: Moving: These are only ever used to report the weights out

    IF (Unit .EQ. SI) THEN
      CopperVol=PI*(DisLnPAR(2)**2-(DisLnPAR(2)-2*DisLnPAR(3))**2)/4/1000*DisLnPAR(1)
      WeightDisLn=CopperVol*CopperDensity

      CopperVol=PI*(SucLnPAR(2)**2-(SucLnPAR(2)-2*SucLnPAR(3))**2)/4/1000*SucLnPAR(1)
      WeightSucLn=CopperVol*CopperDensity

      CopperVol=PI*(LiqLnPAR(2)**2-(LiqLnPAR(2)-2*LiqLnPAR(3))**2)/4/1000*LiqLnPAR(1)
      WeightLiqLn=CopperVol*CopperDensity

      CopperVol=PI*(ValveIDCLnPAR(2)**2-(ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3))**2)/4/1000*ValveIDCLnPAR(1)
      WeightValveIDCLn=CopperVol*CopperDensity

      CopperVol=PI*(ValveODCLnPAR(2)**2-(ValveODCLnPAR(2)-2*ValveODCLnPAR(3))**2)/4/1000*ValveODCLnPAR(1)
      WeightValveODCLn=CopperVol*CopperDensity
    ELSE
      CopperVol=PI*(DisLnPAR(2)**2-(DisLnPAR(2)-2*DisLnPAR(3)/1000)**2)/4/144*DisLnPAR(1)
      CopperVol=CopperVol/35.31467 !Convert from ft3 to m3
      WeightDisLn=CopperVol*CopperDensity

      CopperVol=PI*(SucLnPAR(2)**2-(SucLnPAR(2)-2*SucLnPAR(3)/1000)**2)/4/144*SucLnPAR(1)
      CopperVol=CopperVol/35.31467 !Convert from ft3 to m3
      WeightSucLn=CopperVol*CopperDensity

      CopperVol=PI*(LiqLnPAR(2)**2-(LiqLnPAR(2)-2*LiqLnPAR(3)/1000)**2)/4/144*LiqLnPAR(1)
      CopperVol=CopperVol/35.31467 !Convert from ft3 to m3
      WeightLiqLn=CopperVol*CopperDensity

      CopperVol=PI*(ValveIDCLnPAR(2)**2-(ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3)/1000)**2)/4/144*ValveIDCLnPAR(1)
      CopperVol=CopperVol/35.31467 !Convert from ft3 to m3
      WeightValveIDCLn=CopperVol*CopperDensity

      CopperVol=PI*(ValveODCLnPAR(2)**2-(ValveODCLnPAR(2)-2*ValveODCLnPAR(3)/1000)**2)/4/144*ValveODCLnPAR(1)
      CopperVol=CopperVol/35.31467 !Convert from ft3 to m3
      WeightValveODCLn=CopperVol*CopperDensity

    END IF

    IF (SystemType .EQ. HEATPUMP) THEN

      !CondPAR(58)=HeatingDistubeLength   !RS: Debugging: Not really used

      IF (IsCoolingMode .GT. 0) THEN    !Cooling Mode

        IF (Unit .EQ. SI) THEN !SI unit   !Equilibrium Discharge line, combines compressor discharge line and valve to ODC line
          VolDisLn=PI*((DisLnPAR(2)-2*DisLnPAR(3))*1e-3)**2/4*DisLnPAR(1)
          VolValveODCLn=PI*((ValveODCLnPAR(2)-2*ValveODCLnPAR(3))*1e-3)**2/4*ValveODCLnPAR(1)
        ELSE
          VolDisLn=PI*((DisLnPAR(2)-2*DisLnPAR(3)*1e-3)/12)**2/4*DisLnPAR(1)
          VolValveODCLn=PI*((ValveODCLnPAR(2)-2*ValveODCLnPAR(3)*1e-3)/12)**2/4*ValveODCLnPAR(1)
        END IF

        TotVolume=VolDisLn+VolValveODCLn

        IF (VolValveODCLn .LE. 0) THEN
          EqDiameter=DisLnPAR(2)
          EqThickness=DisLnPAR(3)
          TotElevation=DisLnPAR(4)
          TotHeatGain=DisLnPAR(5)
          TotTempChange=DisLnPAR(6)
          TotAddDP=DisLnPAR(7)
        ELSE
          EqLength=DisLnPAR(1)+ValveODCLnPAR(1) !ISI - 08/03/06
          EqThickness=(DisLnPAR(3)+ValveODCLnPAR(3))/2
          TotElevation=DisLnPAR(4)+ValveODCLnPAR(4)
          TotHeatGain=DisLnPAR(5)+ValveODCLnPAR(5)
          TotTempChange=DisLnPAR(6)+ValveODCLnPAR(6)
          TotAddDP=DisLnPAR(7)+ValveODCLnPAR(7)
        END IF

        IF (UNIT .EQ. SI) THEN
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*1000+2*EqThickness !ISI - 08/03/06
        ELSE
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*12+2*EqThickness*1e-3 !ISI - 08/03/06
        END IF

        DisLnPAR(1)=EqLength
        DisLnPAR(2)=EqDiameter
        DisLnPAR(3)=EqThickness
        DisLnPAR(4)=TotElevation
        DisLnPAR(5)=TotHeatGain
        DisLnPAR(6)=TotTempChange
        DisLnPAR(7)=TotAddDP

        IF (UNIT .EQ. SI) THEN        !Equilibrium suction line, combines compressor suction line and valve to IDC line
          VolSucLn=PI*((SucLnPAR(2)-2*SucLnPAR(3))*1e-3)**2/4*SucLnPAR(1)
          VolValveIDCLn=PI*((ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3))*1e-3)**2/4*ValveIDCLnPAR(1)
        ELSE
          VolSucLn=PI*((SucLnPAR(2)-2*SucLnPAR(3)*1e-3)/12)**2/4*SucLnPAR(1)
          VolValveIDCLn=PI*((ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3)*1e-3)/12)**2/4*ValveIDCLnPAR(1)
        END IF

        TotVolume=VolSucLn+VolValveIDCLn

        IF (VolValveIDCLn .LE. 0) THEN
          EqDiameter=SucLnPAR(2)
          EqThickness=SucLnPAR(3)
          TotElevation=SucLnPAR(4)
          TotHeatGain=SucLnPAR(5)
          TotTempChange=SucLnPAR(6)
          TotAddDP=SucLnPAR(7)
        ELSE
          EqLength=SucLnPAR(1)+ValveIDCLnPAR(1) !ISI - 08/03/06
          EqThickness=(SucLnPAR(3)+ValveIDCLnPAR(3))/2
          TotElevation=SucLnPAR(4)+ValveIDCLnPAR(4)
          TotHeatGain=SucLnPAR(5)+ValveIDCLnPAR(5)
          TotTempChange=SucLnPAR(6)+ValveIDCLnPAR(6)
          TotAddDP=SucLnPAR(7)+ValveIDCLnPAR(7)
        END IF

        IF (UNIT .EQ. SI) THEN
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*1000+2*EqThickness !ISI - 08/03/06
        ELSE
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*12+2*EqThickness*1e-3 !ISI - 08/03/06
        END IF
        SucLnPAR(1)=EqLength
        SucLnPAR(2)=EqDiameter
        SucLnPAR(3)=EqThickness
        SucLnPAR(4)=TotElevation
        SucLnPAR(5)=TotHeatGain
        SucLnPAR(6)=TotTempChange
        SucLnPAR(7)=TotAddDP

      ELSE !Heating mode

        IF (Unit .EQ. SI) THEN !SI unit   !Equilibrium Discharge line, combines compressor discharge line and valve to IDC line
          VolDisLn=PI*((DisLnPAR(2)-2*DisLnPAR(3))*1e-3)**2/4*DisLnPAR(1)
          VolValveIDCLn=PI*((ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3))*1e-3)**2/4*ValveIDCLnPAR(1)
        ELSE
          VolDisLn=PI*((DisLnPAR(2)-2*DisLnPAR(3)*1e-3)/12)**2/4*DisLnPAR(1)
          VolValveIDCLn=PI*((ValveIDCLnPAR(2)-2*ValveIDCLnPAR(3)*1e-3)/12)**2/4*ValveIDCLnPAR(1)
        END IF

        TotVolume=VolDisLn+VolValveIDCLn

        IF (VolValveIDCLn .LE. 0) THEN
          EqDiameter=DisLnPAR(2)
          EqThickness=DisLnPAR(3)
          TotElevation=DisLnPAR(4)
          TotHeatGain=DisLnPAR(5)
          TotTempChange=DisLnPAR(6)
          TotAddDP=DisLnPAR(7)
        ELSE
          EqLength=DisLnPAR(1)+ValveIDCLnPAR(1) !ISI - 08/03/06
          EqThickness=(DisLnPAR(3)+ValveIDCLnPAR(3))/2
          TotElevation=DisLnPAR(4)+ValveIDCLnPAR(4)
          TotHeatGain=DisLnPAR(5)+ValveIDCLnPAR(5)
          TotTempChange=DisLnPAR(6)+ValveIDCLnPAR(6)
          TotAddDP=DisLnPAR(7)+ValveIDCLnPAR(7)
        END IF

        IF (UNIT .EQ. SI) THEN
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*1000+2*EqThickness !ISI - 08/03/06
        ELSE
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*12+2*EqThickness*1e-3 !ISI - 08/03/06
        END IF

        DisLnPAR(1)=EqLength
        DisLnPAR(2)=EqDiameter
        DisLnPAR(3)=EqThickness
        DisLnPAR(4)=TotElevation
        DisLnPAR(5)=TotHeatGain
        DisLnPAR(6)=TotTempChange
        DisLnPAR(7)=TotAddDP

        IF (UNIT .EQ. SI) THEN !Equilibrium suction line, combines compressor suction line and valve to ODC line
          VolSucLn=PI*((SucLnPAR(2)-2*SucLnPAR(3))*1e-3)**2/4*SucLnPAR(1)
          VolValveODCLn=PI*((ValveODCLnPAR(2)-2*ValveODCLnPAR(3))*1e-3)**2/4*ValveODCLnPAR(1)
        ELSE
          VolSucLn=PI*((SucLnPAR(2)-2*SucLnPAR(3)*1e-3)/12)**2/4*SucLnPAR(1)
          VolValveODCLn=PI*((ValveODCLnPAR(2)-2*ValveODCLnPAR(3)*1e-3)/12)**2/4*ValveODCLnPAR(1)
        END IF

        TotVolume=VolSucLn+VolValveODCLn

        IF (VolValveODCLn .LE. 0) THEN
          EqDiameter=SucLnPAR(2)
          EqThickness=SucLnPAR(3)
          TotElevation=SucLnPAR(4)
          TotHeatGain=SucLnPAR(5)
          TotTempChange=SucLnPAR(6)
          TotAddDP=SucLnPAR(7)
        ELSE
          EqLength=SucLnPAR(1)+ValveODCLnPAR(1) !ISI - 08/03/06
          EqThickness=(SucLnPAR(3)+ValveODCLnPAR(3))/2
          TotElevation=SucLnPAR(4)+ValveODCLnPAR(4)
          TotHeatGain=SucLnPAR(5)+ValveODCLnPAR(5)
          TotTempChange=SucLnPAR(6)+ValveODCLnPAR(6)
          TotAddDP=SucLnPAR(7)+ValveODCLnPAR(7)
        END IF

        IF (UNIT .EQ. SI) THEN
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*1000+2*EqThickness !ISI - 08/03/06
        ELSE
          EqDiameter=SQRT(4*TotVolume/(PI*EqLength))*12+2*EqThickness*1e-3 !ISI - 08/03/06
        END IF

        SucLnPAR(1)=EqLength
        SucLnPAR(2)=EqDiameter
        SucLnPAR(3)=EqThickness
        SucLnPAR(4)=TotElevation
        SucLnPAR(5)=TotHeatGain
        SucLnPAR(6)=TotTempChange
        SucLnPAR(7)=TotAddDP

      END IF

    END IF

    CondPAR%CondDisLnLen=DisLnPAR(1) !Discharge line length, m or ft !RS: Debugging: Formerly CondPAR(1)
    CondPAR%CondDisLnOD=DisLnPAR(2) !Discharge line outside diameter, mm or in !RS: Debugging: Formerly CondPAR(2)
    CondPAR%CondDisLnTWThick=DisLnPAR(3) !Discharge line tube wall thickness, mm or mil !RS: Debugging: Formerly CondPAR(3)
    CondPAR%CondDisLnElev=DisLnPAR(4) !Discharge line elevation, m or ft !RS: Debugging: Formerly CondPAR(4)
    CondPAR%CondDisLnQLoss=DisLnPAR(5) !Discharge line heat loss, W or Btu/hr !RS: Debugging: Formerly CondPAR(5)
    CondPAR%CondDisLnTempChg=DisLnPAR(6) !Discharge line temperature change, C or F !RS: Debugging: Formerly CondPAR(6)
    CondPAR%CondDisLnAddPD=DisLnPAR(7) !Discharge line additional pressure drop   !RS: Debugging: Formerly CondPAR(7)

    CondPAR%CondLiqLnLen=LiqLnPAR(1)  !Liquid line length, m or ft  !RS: Debugging: Formerly CondPAR(8)
    CondPAR%CondLiqLnOD=LiqLnPAR(2)  !Liquid line outside diameter, mm or in   !RS: Debugging: Formerly CondPAR(9)
    CondPAR%CondLiqLnTWThick=LiqLnPAR(3) !Liquid line tube wall thickness, mm or mil   !RS: Debugging: Formerly CondPAR(10)
    CondPAR%CondLiqLnElev=LiqLnPAR(4) !Liquid line elevation, m or ft   !RS: Debugging: Formerly CondPAR(11)
    CondPAR%CondLiqLnQLoss=LiqLnPAR(5) !Liquid line heat loss, W or Btu/hr   !RS: Debugging: Formerly CondPAR(12)
    CondPAR%CondLiqLnTempChg=LiqLnPAR(6) !Liquid line temperature change, C or F   !RS: Debugging: Formerly CondPAR(13)
    CondPAR%CondLiqLnAddPD=LiqLnPAR(7) !Liquid line additional pressure drop !RS: Debugging: Formerly CondPAR(14)

    EvapPAR%EvapSucLnLen=SucLnPAR(1) !Suction line length, m or ft  !RS: Debugging: Formerly EvapPAR(1)
    EvapPAR%EvapSucLnOD=SucLnPAR(2) !Suction line outside diameter, mm or in   !RS: Debugging: Formerly EvapPAR(2)
    EvapPAR%EvapSucLnTWThick=SucLnPAR(3) !Suction line tube wall thickness, mm or mil   !RS: Debugging: Formerly EvapPAR(3)
    EvapPAR%EvapSucLnElev=SucLnPAR(4) !Suction line elevation, m or ft   !RS: Debugging: Formerly EvapPAR(4)
    EvapPAR%EvapSucLnQLoss=SucLnPAR(5) !Suction line heat loss, W or Btu/hr   !RS: Debugging: Formerly EvapPAR(5)
    EvapPAR%EvapSucLnTempChg=SucLnPAR(6) !Suction line temperature change, C or F   !RS: Debugging: Formerly EvapPAR(6)
    EvapPAR%EvapSucLnAddPD=SucLnPAR(7) !Suction line additional pressure drop !RS: Debugging: Formerly EvapPAR(7)

    IF (IsCoolingMode .GT. 0) THEN    !Populating arrays
      !***************** Indoor coil data ***************** !RS: Debugging: Evaporator & Condenser

      CALL GetObjectItem('IndoorCoilData',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      EvapPAR%EvapNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly EvapPAR(19)

      !----------

      CoilParams(1)%AirFlowRate=CFMevp
      CoilParams(2)%AirFlowRate=CFMcnd

      ShTbPAR%ShTbTLen=CoolingShTbPAR(1)
      ShTbPAR%ShTbTID=CoolingShTbPAR(2)
      ShTbPAR%ShTbChamDep=CoolingShTbPAR(3)
      ShTbPAR%ShTbECktNum=EvapPAR%EvapNumCkt !Number of circuits in evaporator    !RS: Debugging: Formerly EvapPAR(19), ShTbPAR(4)
      ShTbPAR%ShTbDTubeLen=CoolingDistubeLength !RS: Debugging: Formerly ShTbPAR(5)

      CapTubePAR%CTTubeID=CoolingCapTubePAR(1)
      CapTubePAR%CTTubeLen=CoolingCapTubePAR(2)
      CapTubePAR%CTTubeCoilD=CoolingCapTubePAR(3)
      CapTubePAR%CTEvapCktNum=EvapPAR%EvapNumCkt !Number of circuits in evaporator !RS: Debugging: Formerly EvapPAR(19), CapTubePAR(4)
      CapTubePAR%CTDisTubeLen=CoolingDistubeLength  !RS: Debugging: Formerly CapTubePAR(5)

      ExpDevice=CoolingExpDevice
    ELSE
      CALL GetObjectItem('OutdoorCoilData',1,Alphas,NumAlphas, &  !RS: Debugging:
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      EvapPAR%EvapNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly EvapPAR(19)

      !-------

      CoilParams(1)%AirFlowRate=CFMcnd
      CoilParams(2)%AirFlowRate=CFMevp

      ShTbPAR%ShTbTLen=HeatingShTbPAR(1)
      ShTbPAR%ShTbTID=HeatingShTbPAR(2)
      ShTbPAR%ShTbChamDep=HeatingShTbPAR(3)
      ShTbPAR%ShTbECktNum=EvapPAR%EvapNumCkt !Number of circuits in evaporator    !RS: Debugging: Formerly EvapPAR(19), ShTbPAR(4)
      ShTbPAR%ShTbDTubeLen=HeatingDistubeLength !RS: Debugging: Formerly ShTbPAR(5)

      CapTubePAR%CTTubeID=HeatingCapTubePAR(1)
      CapTubePAR%CTTubeLen=HeatingCapTubePAR(2)
      CapTubePAR%CTTubeCoilD=HeatingCapTubePAR(3)
      CapTubePAR%CTEvapCktNum=EvapPAR%EvapNumCkt !Number of circuits in evaporator !RS: Debugging: Formerly EvapPAR(19), CapTubePAR(4)
      CapTubePAR%CTDisTubeLen=HeatingDistubeLength  !RS: Debugging: CapTubePAR(5)

      ExpDevice=HeatingExpDevice
    END IF

    EvapPAR%EvapBarPress=BaroPressure    !RS: Debugging: Formerly EvapPAR(31)
    CondPAR%CondBarPress=BaroPressure    !RS: Debugging: Formerly CondPAR(38)

    !IF (UNIT .EQ. IP) THEN
    !    CondPAR%CondBarPress=CondPAR%CondBarPress*UnitP
    !    EvapPAR%EvapBarPress=EvapPAR%EvapBarPress*UnitP
    !END IF

    EvapPAR%EvapSysType=SystemType !RS: Debugging: Formerly EvapPAR(34)
    CondPAR%CondSysType=SystemType !ISI - 07/14/06    !RS: Debugging: Formerly CONDPAR(57)

    IF (LineData(1:17) .EQ. 'Microchannel Coil') THEN
      IF (IsCoolingMode .GT. 0) THEN
        ODCcoilType = MCCONDENSER
      ELSE
        ODCcoilType = MCEVAPORATOR
      END IF
    ELSE
      IF (IsCoolingMode .GT. 0) THEN
        ODCcoilType = CONDENSERCOIL
      ELSE
        ODCcoilType = EVAPORATORCOIL
      END IF
    END IF

    !CLOSE(11)    !RS: Debugging: Removing because the file number isn't actually attached to the file (5/21/14)

    RETURN

  END SUBROUTINE

  !***********************************************************************************

  SUBROUTINE GetGenOptInputs()    !RS: Debugging: Trying to read from the modified GenOpt IDF (9/22/14)

    USE InputProcessor
    USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
    USE DataSimulation, ONLY: inputratio    !RS: The compressor capacity ratio for use in the new compressor model (10/3/14)

    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    !CHARACTER(len=MaxNameLength)InputLine   !RS: Debugging: Trying to read from the modified GenOpt IDF (9/22/14)
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    CALL ProcessInput2

    !RS: Debugging: Trying to read from the modified GenOpt IDF (9/22/14)
    CALL GetObjectItem2('HPSimVariables',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    !Numbers(1) is the compressor capacity input ratio
    inputratio=Numbers(1)

    IF (IsCoolingMode .GT. 0) THEN    !Populating arrays
      CFMcnd = Numbers(2)    !Fan Air Flow Rate
      CFMevp = Numbers(3)
    ELSE
      CFMevp = Numbers(2)
      CFMcnd = Numbers(3)
    END IF

  END SUBROUTINE GetGenOptInputs

END MODULE HeatPumpInput

!***********************************************************************************
