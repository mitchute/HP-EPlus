! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module models a compressor.

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component represents an ARI compressor by use of a 2nd order polynomial equation fit.
! A description of the component is found at:
! http://www.e-refrigeration.com/learn-refrigeration/system-components/compressors
! From that website: 
!  - A compressor compresses and cycles refrigerant, while also increasing the pressure.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module reads in entering refrigerant and compressor properties, calculates, and returns exiting refrigerant properties and power.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no input or output files directly connected to this module.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! DataSimulation is called in; otherwise, nothing is defined at the module level.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 2 methods:
!    PUBLIC Compressor -- Models a 2nd order polynomail equation fit compressor model
!      Called by ORNLSolver.f90
!      Called by FlowRateLoop.f90
!    PRIVATE X -- A 10 coefficient polynomial equation calculation
!      Called by Compressor

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
! Some more documentation or narrative might be useful.

    MODULE CompressorMod

    USE DataSimulation
    implicit none

    PUBLIC  Compressor
    PRIVATE X
    PUBLIC CompressorSimple

    CONTAINS

    SUBROUTINE CompressorSimple(Ref$) !Compressor(Ref$) !,XIN,PAR,OUT) !(Ref$,PureRef,XIN,PAR,OUT) !RS: Debugging: Extraneous PureRef
        !RS: Implementation: Commenting out Compressor call since adding in a different comp model (10/1/14)
    ! ----------------------------------------------------------------------
    !
    !   Description: ARI Compressor model
    !
    !   Method: 2nd order polynomial equation fit 
    !
    !   Inputs:
    !       Ref$=Refrigerant name
    !       PureRef=Refrigerant flag: 1=pure refrigerant
    !                                 0=refrigerant mixture
    !       XIN(1) = Suction refrigerant pressure, kPa
    !       XIN(2) = Discharge refrigerant pressure, kPa
    !       XIN(3) = Suction enthalpy, kJ/kg
    !
    !       PAR(1..10)  = Coefficients for power calc.
    !       PAR(11..20) = Coefficients for mass flow rate calc.
    !       PAR(21) = Shell heat loss fraction of power consumption
    !       PAR(22) = Shell heat loss W
    !       PAR(23) = Internal void volume, m^3
    !       PAR(24) = Power Correction, 1 - 230 VAC; 2 - 208 VAC
    !       PAR(25) = Power multiplier
    !       PAR(26) = Mass flow rate multiplier
    !
    !   Outputs:
    !       OUT(1) = Power, KW
    !       OUT(2) = Mass flow rate, kg/s
    !       OUT(3) = Discharge enthalpy, kJ/kg
    !       OUT(4) = Discharge quality
    !       OUT(5) = Discharge temperature, C
    !       OUT(6) = Mass in compressor, kg
    !       OUT(7) = Error flag: 0-No error
    !                            1-Compressor solution error
    !                            2-Refprop error
    !
    !   Author:
    !   Ipseng Iu
    !   Mechanical and Aerospace Engineering
    !   Oklahoma State University, Stillwater
    !
    !   Date: July 2003
    !
    ! ----------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex, MaxNameLength,RefName   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
    USE InputProcessor    !RS: Debugging: Brought over from GetHPSimInputs
    USE UnitConvertMod

    IMPLICIT NONE

    !Subroutine argument declarations
    CHARACTER*80,     INTENT(IN) :: Ref$    !Refrigerant name
    !0-refrigerant mixture

    !REAL, INTENT(IN) :: XIN(3)
    !REAL :: PAR(26) !, INTENT(IN)
    !REAL, INTENT(OUT) :: OUT(7)
    INTEGER,PARAMETER :: SI=1

    !Subroutine local variables
    REAL Temperature,Quality,Pressure,Enthalpy,Entropy

    REAL A(10),B(10)
    REAL TDPsuc      !Suction dew point temperature, C
    REAL TDPdis      !Discharge dew point temperature, C
    REAL TDPsucF     !Suction dew point temperature, F
    REAL TDPdisF     !Discharge dew point temperature, F
    REAL Power       !Compressor power consumption, KW
    REAL PowerMap    !Map based Compressor power consumption, KW
    REAL mdot        !Refrigerant mass flow rate, kg/s
    REAL mdotMap     !Map based Refrigerant mass flow rate, kg/s
    REAL Tsuc        !Suction temp, C
    REAL TsucMap     !Map based suction temp, C
    REAL Tdis        !Discharge temperature, C
    REAL Psuc        !Suction pressure, kPa
    REAL Pdis        !Discharge pressure, kPa
    REAL Xdis        !Discharge quality
    REAL Hsuc        !Suction enthalpy, kJ/kg
    REAL Hdis        !Discharge enthalpy, kJ/kg
    REAL Qshellfrac  !Compressor shell heat loss fraction of power consumption
    REAL Qshell      !Compressor shell heat loss, kW
    REAL VolCmp      !Compressor internal volume, m3
    REAL MassCmp     !Refrigerant mass in compressor, kg
    REAL rhoDis      !Discharge density, kg/m3
    REAL rhoSuc      !Suction density, kg/m3
    REAL rhoMap      !Map based density value, kg/m3
    REAL Ssuc        !Suction entropy
    REAL SsucMap     !Map based entropy value
    REAL HdisIsen    !Isentropic discharge enthalpy, kJ/kg
    REAL HdisIsenMap !Map based insentropic discharge enthalpy, kJ/kg
    REAL HsucMap     !Map based suction enthalpy, kJ/kg
    REAL Wcorrect    !Correction factor for power calc. with different input voltage
!    REAL Mcorrect    !Correction factor for mass flow rate !RS: Debugging: Removed because it's hardcoded to 1
    REAL PwrMultiplier    !Power multiplier
    REAL mdotMultiplier   !Mass flow rate multiplier
    INTEGER I !Loop control
    INTEGER ErrorFlag          !0-No error
    !1-Compressor solution error
    !2-Refprop error
    REAL, SAVE:: FirstTime=1 !Setting a first time variable

    INTEGER(2) RefPropErr  !Error flag:1-error; 0-no error
    LOGICAL, EXTERNAL :: IssueRefPropError
    
    !RS: Debugging: Bringing this over from GetHPSimInputs
    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    INTEGER CompressorManufacturer
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    
    !Compressor Manufacturer
INTEGER,PARAMETER :: COPELAND  = 1 !ISI - 10/05/06
INTEGER,PARAMETER :: BRISTOL   = 2
INTEGER,PARAMETER :: DANFOSS   = 3
INTEGER,PARAMETER :: PANASONIC = 4
    
    !RS: Debugging: Bringing this over from GetHPSimInputs
      !***************** Compressor data *****************  !RS: Debugging: Moving: Compressor
IF (FirstTime .EQ. 1) THEN
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

  CompPAR%CompQLossFrac = Numbers(1) !CompressorHeatLossFraction  !RS: Debugging: Formerly PAR(21)
  CompPAR%CompQLoss = Numbers(2) !CompressorHeatLoss  !RS: Debugging: Formerly PAR(22)
  CompPAR%CompIntVol = Numbers(3) !CompressorVolume    !RS: Debugging: Formerly PAR(23)
  CompPAR%CompCoeffM1 = Numbers(4) !CompressorMassCoefficient1  !RS: Debugging: Formerly PAR(11)
  CompPAR%CompCoeffM2 = Numbers(5) !CompressorMassCoefficient2  !RS: Debugging: Formerly PAR(12)
  CompPAR%CompCoeffM3 = Numbers(6) !CompressorMassCoefficient3  !RS: Debugging: Formerly PAR(13)
  CompPAR%CompCoeffM4 = Numbers(7) !CompressorMassCoefficient4  !RS: Debugging: Formerly PAR(14)
  CompPAR%CompCoeffM5 = Numbers(8) !CompressorMassCoefficient5  !RS: Debugging: Formerly PAR(15)
  CompPAR%CompCoeffM6 = Numbers(9) !CompressorMassCoefficient6  !RS: Debugging: Formerly PAR(16)
  CompPAR%CompCoeffM7 = Numbers(10) !CompressorMassCoefficient7 !RS: Debugging: Formerly PAR(17)
  CompPAR%CompCoeffM8 = Numbers(11) !CompressorMassCoefficient8 !RS: Debugging: Formerly PAR(18)
  CompPAR%CompCoeffM9 = Numbers(12) !CompressorMassCoefficient9 !RS: Debugging: Formerly PAR(19)
  CompPAR%CompCoeffM10 = Numbers(13) !CompressorMassCoefficient10    !RS: Debugging: Formerly PAR(20)
  CompPAR%CompCoeffP1 = Numbers(14) !CompressorPowerCoefficient1 !RS: Debugging: Formerly PAR(1)
  CompPAR%CompCoeffP2 = Numbers(15) !CompressorPowerCoefficient2 !RS: Debugging: Formerly PAR(2)
  CompPAR%CompCoeffP3 = Numbers(16) !CompressorPowerCoefficient3 !RS: Debugging: Formerly PAR(3)
  CompPAR%CompCoeffP4 = Numbers(17) !CompressorPowerCoefficient4 !RS: Debugging: Formerly PAR(4)
  CompPAR%CompCoeffP5 = Numbers(18) !CompressorPowerCoefficient5 !RS: Debugging: Formerly PAR(5)
  CompPAR%CompCoeffP6 = Numbers(19) !CompressorPowerCoefficient6 !RS: Debugging: Formerly PAR(6)
  CompPAR%CompCoeffP7 = Numbers(20) !CompressorPowerCoefficient7 !RS: Debugging: Formerly PAR(7)
  CompPAR%CompCoeffP8 = Numbers(21) !CompressorPowerCoefficient8 !RS: Debugging: Formerly PAR(8)
  CompPAR%CompCoeffP9 = Numbers(22) !CompressorPowerCoefficient9 !RS: Debugging: Formerly PAR(9)
  CompPAR%CompCoeffP10 = Numbers(23) !CompressorPowerCoefficient10   !RS: Debugging: Formerly PAR(10)
  
  CompPAR%CompPwrMult = Numbers(24) !PowerMultiplier    !RS: Debugging: Formerly PAR(25)
  CompPAR%CompMFRMult = Numbers(25) !MassFlowRateMultiplier !RS: Debugging: Formerly PAR(26)
  
  FirstTime=2
END IF
  !TsiCmp = Numbers(26) !UserSpecifiedRatingEvapTemperature
  !TsoCmp = Numbers(27) !UserSpecifiedRatingCondTemperature
  !Subcool = Numbers(28) !UserSpecifiedRatingSubcooling
  !Super = Numbers(29) !UserSpecifiedRatingSuperheat

    !Flow:

    Psuc = CompIN%CompInPsuc   !RS: Debugging: Formerly XIN(1)
    Pdis = CompIN%CompInPdis   !RS: Debugging: Formerly XIN(2)
    Hsuc = CompIN%CompInHsuc   !RS: Debugging: Formerly XIN(3)

    DO I=1,10
        A(I)= Numbers(13+I) !CompPAR%(I)
        B(I)= Numbers(3+I) !CompPAR%(I+10) 
    END DO
    
    IF (Unit .EQ. SI)THEN !SI unit inputs   !RS: Debugging: 
    	CompPAR%CompIntVol=CompPAR%CompIntVol/(100**3) !Compressor internal volume, m^3   !RS: Formerly CompPAR(23)
    ELSE
        CompPAR%CompQLoss=CompPAR%CompQLoss*UnitPwr*1000 !Compressor shell heat loss W  !RS: Debugging: Formerly CompPAR(22)
        CompPAR%CompIntVol=CompPAR%CompIntVol/(12**3)*(UnitL**3) !Compressor internal volume, m^3 !RS: Debugging: Formerly CompPAR(23)
    END IF
    
    Qshellfrac = CompPAR%CompQLossFrac    !RS: Debugging: Formerly PAR(21)
    Qshell = CompPAR%CompQLoss    !RS: Debugging: Formerly PAR(22)
    VolCmp = CompPAR%CompIntVol    !RS: Debugging: Formerly PAR(23)
    !Wcorrect = CompPAR%CompPwrCor  !RS: Debugging: Formerly PAR(24)    !RS: Debugging: CompPwrCor is never set
    PwrMultiplier=CompPAR%CompPwrMult   !RS: Debugging: Formerly PAR(25)
    mdotMultiplier=CompPAR%CompMFRMult  !RS: Debugging: Formerly PAR(26)
    
    Wcorrect = 1 !1.21 !1.25

    ErrorFlag=0 !Initialize

    Pressure=Psuc*1000  !RS Comment: Unit Conversion
    Quality=1
    TDPsuc=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Suction Dew Point Temperature, C
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    Pressure=Pdis*1000  !RS Comment: Unit Conversion
    Quality=1
    TDPdis=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Discharge Dew Point Temperature, C
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    TDPsucF=TDPsuc*1.8+32   !Suction Dew Point Temperature, F
    TDPdisF=TDPdis*1.8+32   !Discharge Dew Point Temperature, F

    PowerMap=X(A,TDPdisF,TDPsucF)/1000      !RS Comment: Unit Conversion
    mdotMap=X(B,TDPdisF,TDPsucF)*0.454/3600 !RS Comment: Unit Conversion?

    TsucMap=((TDPsucF+20)-32)*5/9 !20-rated superheat

    Temperature=TsucMap
    Pressure=Psuc*1000  !RS Comment: Unit Conversion
    HsucMap=TP(Ref$,Temperature,Pressure,'enthalpy',RefrigIndex,RefPropErr) !Map-Based Suction Enthalpy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    HsucMap=HsucMap/1000    !RS Comment: Unit Conversion
    rhoMap=TP(Ref$,Temperature,Pressure,'density',RefrigIndex,RefPropErr)   !Map-Based Density
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    SsucMap=TP(Ref$,Temperature,Pressure,'entropy',RefrigIndex,RefPropErr)  !Map-Based Suction Entropy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    SsucMap=SsucMap/1000    !RS Comment: Unit Conversion

    Pressure=Pdis*1000      !RS Comment: Unit Conversion
    Entropy=SsucMap*1000    !RS Comment: Unit Conversion
    HdisIsenMap=PS(Ref$,Pressure,Entropy,'enthalpy',RefrigIndex,RefPropErr) !Map-Based Isentropic Discharge Enthalpy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    HdisIsenMap=HdisIsenMap/1000    !RS Comment: Unit Conversion

    Pressure=Psuc*1000  !RS Comment: Unit Conversion
    Enthalpy=Hsuc*1000  !RS Comment: Unit Conversion
    Tsuc=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)    !Suction Temperature
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    rhosuc=PH(Ref$,Pressure,Enthalpy,'density',RefrigIndex,RefPropErr)  !Suction Density
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    Ssuc=PH(Ref$,Pressure,Enthalpy,'entropy',RefrigIndex,RefPropErr)    !Suction Entropy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    Ssuc=Ssuc/1000  !RS Comment: Unit Conversion

    Pressure=Pdis*1000  !RS Comment: Unit Conversion
    Entropy=Ssuc*1000   !RS Comment: Unit Conversion
    HdisIsen=PS(Ref$,Pressure,Entropy,'enthalpy',RefrigIndex,RefPropErr)    !Isentropic Discharge Enthalpy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    HdisIsen=HdisIsen/1000  !RS Comment: Unit Conversion

    mdot=mdotMap*(rhosuc/rhoMap)

    !RS: If Mcorrect is hardcoded to 1, then it is useless (12/16/13)
    !Mcorrect=1 !0.95
    !mdot=mdot*Mcorrect

    Power=PowerMap*(mdot/mdotMap)*(HdisIsen-Hsuc)/(HdisIsenMap-HsucMap)
    
    IF (HdisIsenMap .EQ. HsucMap .AND. HdisIsen .EQ. HSuc) THEN !RS: Debugging: If both are the same, the ratio should be 1
        Power=PowerMap*(mdot/mdotMap)
    END IF

    Power=Power/Wcorrect

    mdot=mdot*mdotMultiplier
    Power=Power*PwrMultiplier

    IF (Qshellfrac .NE. 0) THEN
        Qshell = Qshellfrac * Power !Fraction of power input
    ELSE
        Qshell = Qshell * 1E-3  !Convert to kW
    END IF
    Hdis=(Power-Qshell)/mdot+Hsuc

    Pressure=Pdis*1000  !RS Comment: Unit Conversion
    Enthalpy=Hdis*1000  !RS Comment: Unit Conversion
    Tdis=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)    !Discharge Temperature
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
            RETURN
    END IF

    Xdis=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)    !Discharge Quality
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    rhoDis=PH(Ref$,Pressure,Enthalpy,'density',RefrigIndex,RefPropErr)  !Discharge Density
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    MassCmp=VolCmp*(rhoDis+rhoSuc)/2

    IF (Power .LT. 0 .OR. MassCmp .LT. 0) THEN
        ErrorFlag=1
        CompOUT%CmpOErrFlag=ErrorFlag    !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF

    CompOUT%CmpOPwr=Power    !RS: Debugging: Formerly OUT(1)
    CompOUT%CmpOMdot=mdot !RS: Debugging: Formerly OUT(2)
    CompOUT%CmpOHdis=Hdis !RS: Debugging: Formerly OUT(3)
    !OUT(4)=Xdis    !RS: Debugging: Never used
    CompOUT%CmpOTdis=Tdis !RS: Debugging: Formerly OUT(5)
    CompOUT%CmpOMCmp=MassCmp  !RS: Debugging: Formerly OUT(6)
    CompOUT%CmpOErrFlag=ErrorFlag    !RS: Debugging: Formerly OUT(7)

    RETURN

        END SUBROUTINE CompressorSimple !Compressor 
        !RS: Implementation: Commenting out Compressor call since adding in a different comp model (10/1/14)

    !***********************************************************************

    REAL FUNCTION X(C,D,S)
    implicit none

    REAL C(10) !Coefficients
    REAL D !Discharge dew point temperature
    REAL S !Suction dew point temperature

    X=C(1)+C(2)*S+C(3)*D+C(4)*S**2+C(5)*(S*D)+C(6)*D**2+ &
    C(7)*S**3+C(8)*D*S**2+C(9)*S*D**2+C(10)*D**3

    END FUNCTION

    !***********************************************************************
    
    SUBROUTINE Compressor(Ref$) !CompressorSimple(Ref$) 
        !RS: Implementation: Commenting out Compressor call since adding in a different comp model (10/1/14)
    !Algorithm and structure by John Gall
    !Program written by Rachel Spitler
    
    USE FluidProperties_HPSim
    USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex, MaxNameLength,RefName   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
    USE InputProcessor    !RS: Debugging: Brought over from GetHPSimInputs
    USE UnitConvertMod
    USE DataSimulation, ONLY: inputratio    !RS: The compressor capacity ratio for use in the new compressor model (10/3/14)
    
    CHARACTER*80,     INTENT(IN) :: Ref$    !Refrigerant name
    
    REAL Temperature,Quality,Pressure,Enthalpy,Entropy
    INTEGER(2) RefPropErr  !Error flag:1-error; 0-no error
    LOGICAL, EXTERNAL :: IssueRefPropError
    
    REAL volcoeff1
    REAL volcoeff2
    REAL volcoeff3
    REAL volcoeff4
    REAL volcoeff5
    REAL isocoeff1
    REAL isocoeff2
    REAL isocoeff3
    REAL isocoeff4
    REAL isocoeff5
    REAL isocoeff6
    REAL isocoeff7
    
    REAL full_load_power    !Power at full load
    !REAL, SAVE:: inputratio    !0 to 1, percentage of capacity required to handle load
    REAL suction_density    !kg/m^3
    REAL psuction   !MPa
    REAL pratio     !Ratio of discharge to suction pressure
    REAL isEff
    REAL m_dot
    REAL suction_volume !m^3
    REAL volumetric_efficiency
    REAL n  !Cycles per Second (60 Hz)
    REAL hsd
    REAL SATS   !Saturation suction temperature
    REAL SATD   !Saturated discharge temperature
    REAL SATSC
    REAL SATDC
    REAL Ssuc
    
    REAL, SAVE:: FirstTime=1 !Setting a first time variable

!    INTEGER(2) RefPropErr  !Error flag:1-error; 0-no error
!    LOGICAL:: IssueRefPropError
    INTEGER ErrorFlag          !0-No error
    
    !RS: Debugging: Bringing this over from GetHPSimInputs
    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    INTEGER CompressorManufacturer
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    
    !Psuc = CompIN%CompInPsuc   !RS: Debugging: Formerly XIN(1)
    !Pdis = CompIN%CompInPdis   !RS: Debugging: Formerly XIN(2)
    !Hsuc = CompIN%CompInHsuc   !RS: Debugging: Formerly XIN(3)
    
    volcoeff1 = 0.09605984
    volcoeff2 = 0.347651
    volcoeff3 = 0.4530326
    volcoeff4 = -0.001568792
    volcoeff5 = -6.44515E-06

    isocoeff1 = -0.454209
    isocoeff2 = 0.1517346
    isocoeff3 = 0.09022652
    isocoeff4 = 0.2497634
    isocoeff5 = -0.2206927
    isocoeff6 = 0.508598
    isocoeff7 = -0.08659438
    
IF (FirstTime .EQ. 1) THEN
  CALL GetObjectItem('CompressorData',1,Alphas,NumAlphas, &
                      TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    
    CompPAR%CompCoeffP1 = Numbers(14) !CompressorPowerCoefficient1 !RS: Debugging: Formerly PAR(1)
    CompPAR%CompCoeffP2 = Numbers(15) !CompressorPowerCoefficient2 !RS: Debugging: Formerly PAR(2)
    CompPAR%CompCoeffP3 = Numbers(16) !CompressorPowerCoefficient3 !RS: Debugging: Formerly PAR(3)
    CompPAR%CompCoeffP4 = Numbers(17) !CompressorPowerCoefficient4 !RS: Debugging: Formerly PAR(4)
    CompPAR%CompCoeffP5 = Numbers(18) !CompressorPowerCoefficient5 !RS: Debugging: Formerly PAR(5)
    CompPAR%CompCoeffP6 = Numbers(19) !CompressorPowerCoefficient6 !RS: Debugging: Formerly PAR(6)
    CompPAR%CompCoeffP7 = Numbers(20) !CompressorPowerCoefficient7 !RS: Debugging: Formerly PAR(7)
    CompPAR%CompCoeffP8 = Numbers(21) !CompressorPowerCoefficient8 !RS: Debugging: Formerly PAR(8)
    CompPAR%CompCoeffP9 = Numbers(22) !CompressorPowerCoefficient9 !RS: Debugging: Formerly PAR(9)
    CompPAR%CompCoeffP10 = Numbers(23) !CompressorPowerCoefficient10   !RS: Debugging: Formerly PAR(10)
    
    !inputratio=Numbers(30)    !To start with; this will need to vary

    FirstTime=2
END IF

    Pressure=CompIN%CompInPsuc*1000  !RS Comment: Unit Conversion
    Quality=1
    SATSC=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Suction Dew Point Temperature, C
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    
    SATS=SATSC*1.8+32 

    Enthalpy=CompIN%CompInHSuc*1000
    suction_density=PH(Ref$,Pressure,Enthalpy,'density',RefrigIndex,RefPropErr)  !Suction Density
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF 
    
    Enthalpy=CompIN%CompInHsuc*1000
    Ssuc=PH(Ref$,Pressure,Enthalpy,'entropy',RefrigIndex,RefPropErr)    !Suction Entropy
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    
    Pressure=CompIN%CompInPdis*1000  !RS Comment: Unit Conversion
    Quality=1
    SATDC=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Discharge Dew Point Temperature, C
    IF (IssueRefPropError(RefPropErr, 'Compressor', 2, ErrorFlag, CompOUT%CmpOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
        RETURN
    END IF
    
    SATD=SATDC*1.8+32 

    psuction=CompIN%CompInPsuc/(1000)  !From kPa to MPa
    pratio=CompIN%CompInPdis/CompIN%CompInPsuc

    full_load_power = CompPAR%CompCoeffP1+CompPAR%CompCoeffP2*SATS+CompPAR%CompCoeffP3*SATD+CompPAR%CompCoeffP4*SATS**2+ &
            CompPAR%CompCoeffP5*SATS*SATD+CompPAR%CompCoeffP6*SATD**2+CompPAR%CompCoeffP7*SATS**3+ &
            CompPAR%CompCoeffP8*SATD*SATS**2+CompPAR%CompCoeffP9*SATS*SATD**2+CompPAR%CompCoeffP10*SATD**3
    
    CompOUT%CmpOPwr = (full_load_power*inputratio + (1-inputratio)*0.1*full_load_power)/1000    !From W to kW
    
    
    
    volumetric_efficiency = volcoeff1+ volcoeff2*inputratio + volcoeff3*inputratio**2 + volcoeff4*suction_density+ &
            volcoeff5*suction_density**2
    
    suction_volume = 3.2984e-5 !(m^3)
    
    n=60
    
    m_dot=suction_density*n*suction_volume*volumetric_efficiency
    
    isEff = isocoeff1 + isocoeff2*inputratio + isocoeff3*inputratio**2 + isocoeff4*psuction + &
            isocoeff5*(psuction)**2 + isocoeff6*pratio + isocoeff7*pratio**2
    
    Pressure=CompIN%CompInPdis*1000
    Entropy=Ssuc
    hsd=PS(Ref$,Pressure,Entropy,'enthalpy',RefrigIndex,RefPropErr)
    CompOUT%CmpOHdis = CompIN%CompInHsuc + ((hsd/1000) - CompIN%CompInHsuc)/isEff !Discharge enthalpy
     
    Pressure=CompIN%CompInPdis*1000  !RS Comment: Unit Conversion
    Enthalpy=CompOUT%CmpOHdis*1000  !RS Comment: Unit Conversion
    CompOUT%CmpOTdis=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)    !Discharge Temperature
     
    CompOUT%CmpOMdot=m_dot !RS: Debugging: Formerly OUT(2)
    
END SUBROUTINE Compressor !CompressorSimple
        !RS: Implementation: Commenting out Compressor call since adding in a different comp model (10/1/14)
        
    !***********************************************************************

    END MODULE CompressorMod
