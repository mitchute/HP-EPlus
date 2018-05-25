MODULE DataSimulation
implicit none

!Conversion factors
REAL, PARAMETER :: UnitPwr = 0.2927     !(Btu/hr X Upower = KW)
REAL, PARAMETER :: UrefFlow = 0.4536    !(Lbm/hr X UrefFlow = kg/hr)
REAL, PARAMETER :: UairFlow = 0.472E-3  !(CFM X UairFlow = kg/s)
REAL, PARAMETER :: UnitArFlw = 0.0004719474 !(CFM X UnitArFlw = m^3/s)
REAL, PARAMETER :: Upressure = 6.895    !(psi X Upressure = kPa)
REAL, PARAMETER :: UairPres = 249.1   !(in-H2o X UairPres = Pa)
REAL, PARAMETER :: Umass = 0.4536       !(lbm X Umass = kg)
REAL, PARAMETER :: Ulength = 0.3048     !(ft X Ulength = m)
REAL, PARAMETER :: Uenthalpy = 2.326     !(Btu/lbm X Uenthalpy = kJ/kg)
REAL, PARAMETER :: Ucndct = 0.1442E-3   !(Btu-in/hr-ft2-F X Ucndct = kW/m-C)
REAL, PARAMETER :: UnitP=6.895 !Pressure unit conversion (psi * UnitP = kPa)
REAL, PARAMETER :: UnitH=2.326 !Enthalpy unit conversion (Btu/lbm * UnitH = kJ/kg)
REAL, PARAMETER :: UnitM=0.4536 !Mass unit conversion (lbm * UnitM = kg)
REAL, PARAMETER :: UnitL=0.3048 !Length unit conversion (ft * UnitL = m) 

!Coil types
INTEGER,PARAMETER :: CONDENSERCOIL  = 1
INTEGER,PARAMETER :: EVAPORATORCOIL = 2
INTEGER,PARAMETER :: HIGHSIDETUBE   = 3
INTEGER,PARAMETER :: LOWSIDETUBE    = 4
INTEGER,PARAMETER :: MCCONDENSER    = 5
INTEGER,PARAMETER :: MCEVAPORATOR   = 6

!System types
INTEGER,PARAMETER :: AIRCONDITIONER = 1
INTEGER,PARAMETER :: HEATPUMP       = 2
INTEGER,PARAMETER :: CONDENSERUNIT  = 3
INTEGER,PARAMETER :: REHEAT         = 4
INTEGER,PARAMETER :: EVAPORATORONLY = 5

!Calculation Mode
INTEGER,PARAMETER :: FIXEDORIFICESIM     = 1
INTEGER,PARAMETER :: ORIFICEANDTXVDESIGN = 2
INTEGER,PARAMETER :: FIXEDSUPERHEATSIM   = 3
INTEGER,PARAMETER :: TXVSIMULATION       = 4
INTEGER,PARAMETER :: CONDENSERUNITSIM    = 5
INTEGER,PARAMETER :: COILONLYSIM       = 6

INTEGER,PARAMETER :: MaxIter = 30 !Maximum number of iterations
!INTEGER,PARAMETER :: NumTimeSteps = 5 !Sankar added transient
!Frosting Simulation Parameters
LOGICAL :: FrostingPeriod=.TRUE.
LOGICAL :: DeFrostingPeriod=.FALSE.
LOGICAL :: DefrostInitiate
REAL :: CurSimTime=0.0
REAL :: PrevSimTime=0.0
REAL :: TimeInterval=0.0
!Defrost Parameters
REAL :: DefrostControlTemp
REAL :: DefrostSetPoint
REAL :: TimeInhibit
REAL :: DefrostInhibitTime = 2500
REAL :: LastDefrostInitTime
REAL,DIMENSION(2) :: AmbTemp
REAL,DIMENSION(2) :: InitiateTemp

INTEGER :: TimeStep=0

INTEGER(2), SAVE:: IsCoolingMode !1=yes; 0=no   !RS: Debugging: Saving this throughout
INTEGER(2), SAVE:: CoilMode !0=condenser, 1=evaporator  !RS: Debugging: This is for a test in ZeroConvergence (11/8/13)


!Expansion device
INTEGER(2) ExpDevice !1=Orifice; 2=TXV; 3=Cap. Tube

INTEGER IDCcoilType !Indoor coil coil type
INTEGER ODCcoilType !Outdoor coil coil type

REAL REFCHG								!Specified refrigerant charge, lbm
REAL CALCHG								!Calculated refrigerant charge, lbm
REAL CondLiqTubeLength		!Condenser liquid tube length, m
REAL CondVapTubeLength		!Condenser vapor tube length, m
REAL CondTwoPhaseTubeLength !Condenser two-phase tube length, m
REAL CondNumLiqTubes		!Number of liquid tubes in condenser
REAL EvapLiqTubeLength		!Evaporator liquid tube length, m
REAL EvapVapTubeLength		!Evaporaotr vapor tube length, m
REAL EvapTwoPhaseTubeLength !Evaporator two-phase tube length, m
REAL EvapNumLiqTubes		!Number of liquid tube length in evaporator
INTEGER IsChargeTuning    !Flag to indicate if charge tuning is performed
REAL ChargeCurveSlope     !Charge curve slope, kg/m
REAL ChargeCurveIntercept !Charge curve intercept, kg
REAL RefLiquidLength      !Liquid length at reference point, m

REAL BaroPressure   !Barometric pressure, kPa

REAL PwrODfan					!Outdoor fan power, W
REAL PwrIDfan					!Indoor fan power, W
REAL XMaC						!Condenser inlet air flow rate, kg/s
REAL TAIC						!Condenser inlet DB temp. F
REAL RHIC						!Condenser inlet relative humidity
REAL XMaE						!Evaporator inlet air flow rate, kg/s
REAL TAIE						!Evaporator inlet DB temp. F
REAL RHIE						!Evaporator inlet relative humidity
REAL DTROC
REAL CFMcnd                  !Standard condenser CFM, m3/s
REAL CFMevp                  !Standard evaperator CFM, m3/s

CHARACTER*80 Ref$	!Refrigerant name 
CHARACTER*80 Rref   !Referance Refrigerant name

INTEGER(2) PureRef	!Pure refrigerant flag: 1=Pure; 0=mixture

!Compressor model passing parameters
REAL TSICMP								!Compressor inlet saturation temperature, F
REAL TSOCMP								!Compressor outlet saturation temperature, F 
REAL SUPER								!Superheat (F) or quality
REAL SUPERE
REAL SUBCOOL							!Subcooling, F

REAL AMBCON						!Convergence criterion for ambient temp.
REAL CNDCON						!Convergence criterion for condenser subcooling
REAL EVPCON						!Convergence criterion for evaporator superheat
REAL FLOCON						!Convergence criterion for flow rate
Logical CoarseConvergenceCriteriaMet  !Flag to refine convergence criteria for final solution

CHARACTER*15  Punit
CHARACTER*15  Hunit
CHARACTER*15  Tunit
CHARACTER*15  DTunit
CHARACTER*15  MdotUnit
CHARACTER*15  MassUnit
CHARACTER*15  PwrUnit
CHARACTER*15  CapUnit
CHARACTER*15  EERunit
CHARACTER*15  SysUnit
CHARACTER*15  NoUnit
CHARACTER*15  Xunit
CHARACTER*15  Lunit
CHARACTER*15  MiniLunit

REAL TICMP,PICMP,HICMP,XICMP	!Compressor inlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TOCMP,POCMP,HOCMP,XOCMP	!Compressor outlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TICND,PICND,HICND,XICND	!Condenser inlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TOCND,POCND,HOCND,XOCND	!Condenser outlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TIEXP,PIEXP,HIEXP,XIEXP	!Exp.device inlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TOEXP,POEXP,HOEXP,XOEXP	!Exp.device outlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TIEVP,PIEVP,HIEVP,XIEVP	!Evaporator inlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality
REAL TOEVP,POEVP,HOEVP,XOEVP	!Evaporator outlet temp. (F); pressure (psi); enthalpy (Btu/lbm); Quality

INTEGER(2)	:: Unit !Unit flag: 1=SI; 2=IP

INTEGER		:: MODE !1-Design mode: Given superheat, subcooling; Compute charge, orifice size   
					!2-Simulation mode: Given charge, orifice size; Compute superheat, subcooling
					!3-Design mode: Given expansion device, superheat; Compute charge, subcooling
					!4-Simulation mode: Given superheat, charge; Compute subcooling, expansion deivce

LOGICAL FirstTimeAirTempLoop	!First time to run AirTempLoop?
LOGICAL FirstTimeFlowRateLoop   !First time to run FlowRateLoop?
LOGICAL FirstTimeHPdesignMode	!First time to run HPdesignMode?
LOGICAL FirstTimeChargeLoop     !First time to run ChargeLoop?
INTEGER(2) :: SystemType !1=A/C; 2=Heat Pump; 3=Condenser Unit; 4=Reheat; 5=Evaporator Only

REAL RhoAiE  
REAL RhoAiC
REAL RhoAoC
REAL RhoAoE

REAL MdotR

REAL :: WeightSucLn      !Weight of suction line, kg or lbm
REAL :: WeightDisLn      !Weight of discharge line, kg or lbm  
REAL :: WeightLiqLn      !Weight of liquid line, kg or lbm
REAL :: WeightValveIDCLn !Weight of Valve to IDC line, kg or lbm
REAL :: WeightValveODCLn !Weight of Valve to ODC line, kg or lbm

REAL TimeStepCorrection         ! Added Sankar Transient Code

REAL Tdis !Discharge temperature, C
REAL Tliq !Liquid temperature, C

REAL DTDISLN					!Discharge line temperature change, F
REAL DTSUCLN					!Suction line temperature change, F
REAL DTLIQLN					!Liquid line temperature change, F
REAL QDISLN						!Discharge line heat loss, Btu/hr
REAL QSUCLN						!Suction line heat gain, Btu/hr
REAL QLIQLN						!Liquid line heat loss, Btu/hr

INTEGER(2) :: CoolHeatModeFlag = -1

TYPE FrostParameters
 REAL :: Density = 0.0
 REAL :: Conductivity = 0.0
 REAL :: Thickness = 0.0
END TYPE FrostParameters

TYPE CoilParameters
 REAL CoilFaceArea
 REAL CoilFinArea
 REAL CoilTubeArea
 REAL CoilFreeFlowArea
 REAL CoilAirSideArea
 REAL CoilInletRefTemp
 REAL CoilOutletRefTemp
 REAL CoilInletAirTemp
 REAL CoilOutletAirTemp
 REAL DPair
 REAL DPRef
 REAL AirFlowRate
 REAL FinPitch
 REAL FinThickness
 REAL TSurfCoil
END TYPE CoilParameters

TYPE(FrostParameters) :: FrostParam
TYPE(CoilParameters), DIMENSION(2):: CoilParams

REAL EvapTubeArea
REAL EvapFinArea
REAL EvapTotArea
REAL EvapBareArea
REAL EvapMinArea
REAL CondTubeArea
REAL CondFinArea

!CondPAR Variables
TYPE CondParameters
    REAL CondDisLnLen
    REAL CondDisLnOD
    REAL CondDisLnTWThick
    REAL CondDisLnElev
    REAL CondDisLnQLoss
    REAL CondDisLnTempChg
    REAL CondDisLnAddPD
    REAL CondLiqLnLen
    REAL CondLiqLnOD
    REAL CondLiqLnTWThick
    REAL CondLiqLnElev
    REAL CondLiqLnQLoss
    REAL CondLiqLnTempChg
    REAL CondLiqLnAddPD
    REAL CondCoilTOD
    REAL CondCoilTWThick
    REAL CondCoilSTLen
    REAL CondCoilTThermCon
    REAL CondTspc
    REAL CondRspc
    REAL CondFinThick
    REAL CondFinPitch
    REAL CondFinThermCon
    INTEGER CondNt
    INTEGER CondNl
    INTEGER CondNumCkt
    INTEGER CondCoolMode
    INTEGER CondNumMod
    INTEGER CondFinType
    REAL CondMultRefQT
    REAL CondMultRefPD
    REAL CondMultAirQT
    REAL CondMultAirPD
    REAL CondFanPwr
    REAL CondFanLoc
    REAL CondSurfAbs
    REAL CondTube
    REAL CondBarPress
    REAL CondCompQLoss
    REAL CondPressTolConv
    INTEGER CondSysType
    REAL CondOilMassFrac
    INTEGER CondCompMan
    INTEGER CondSimpCoil
    INTEGER CondFirstTime
END TYPE CondParameters

TYPE(CondParameters) :: CondPAR
!INTEGER, SAVE:: CondDisLnLen=1
!INTEGER, SAVE:: CondDisLnOD=2
!INTEGER, SAVE:: CondDisLnTWThick=3
!INTEGER, SAVE:: CondDisLnElev=4
!INTEGER, SAVE:: CondDisLnQLoss=5
!INTEGER, SAVE:: CondDisLnTempChg=6
!INTEGER, SAVE:: CondDisLnAddPD=7
!INTEGER, SAVE:: CondLiqLnLen=8
!INTEGER, SAVE:: CondLiqLnOD=9
!INTEGER, SAVE:: CondLiqLnTWThick=10
!INTEGER, SAVE:: CondLiqLnElev=11
!INTEGER, SAVE:: CondLiqLnQLoss=12
!INTEGER, SAVE:: CondLiqLnTempChg=13
!INTEGER, SAVE:: CondLiqLnAddPD=14
!INTEGER, SAVE:: CondCoilTOD=15
!INTEGER, SAVE:: CondCoilTWThick=16
!INTEGER, SAVE:: CondCoilSTLen=17
!INTEGER, SAVE:: CondCoilTThermCon=18
!INTEGER, SAVE:: CondTspc=19
!INTEGER, SAVE:: CondRspc=20
!INTEGER, SAVE:: CondFinThick=21
!INTEGER, SAVE:: CondFinPitch=22
!INTEGER, SAVE:: CondFinThermCon=23
!INTEGER, SAVE:: CondNt=24
!INTEGER, SAVE:: CondNl=25
!INTEGER, SAVE:: CondNumCkt=26
!INTEGER, SAVE:: CondCoolMode=27
!INTEGER, SAVE:: CondNumMod=28
!INTEGER, SAVE:: CondFinType=29
!INTEGER, SAVE:: CondMultRefQT=30
!INTEGER, SAVE:: CondMultRefPD=31
!INTEGER, SAVE:: CondMultAirQT=32
!INTEGER, SAVE:: CondMultAirPD=33
!INTEGER, SAVE:: CondFanPwr=34
!INTEGER, SAVE:: CondFanLoc=35
!INTEGER, SAVE:: CondSurfAbs=36
!INTEGER, SAVE:: CondTube=37
!INTEGER, SAVE:: CondBarPress=38
!INTEGER, SAVE:: CondCompQLoss=39
!INTEGER, SAVE:: CondPressTolConv=40
!INTEGER, SAVE:: CondSysType=41
!INTEGER, SAVE:: CondOilMassFrac=42
!INTEGER, SAVE:: CondCompMan=43
!INTEGER, SAVE:: CondSimpCoil=44
!INTEGER, SAVE:: CondFirstTime=45

!CondIN variables
TYPE CondInlet
    REAL CInmRef
    REAL CInpRo
    REAL CInhRo
    REAL CInmAi
    REAL CIntAi
    REAL CInrhAi
    REAL CInSolFlux
END TYPE CondInlet

TYPE(CondInlet) :: CondIN
!INTEGER, SAVE:: CInmRef=1
!INTEGER, SAVE:: CInpRo=2
!INTEGER, SAVE:: CInhRo=3
!INTEGER, SAVE:: CInmAi=4
!INTEGER, SAVE:: CIntAi=5
!INTEGER, SAVE:: CInrhAi=6
!INTEGER, SAVE:: CInSolFlux=7

!CondOUT variables
TYPE CondOutlet
    REAL COutpRiC
    REAL COuthRiC
    REAL COuttAoC
    REAL COutrhAoC
    REAL COutpRoC
    REAL COuthRoC
    REAL COuttRoC
    REAL COutWtAl
    REAL COutWtCu
    REAL COutpRiE
    REAL COuthRiE
    REAL COuttRiE
    REAL COutxRiE
    REAL COuttSCiE
    REAL COutQC
    REAL COutMDisLn
    REAL COutMLiqLn
    REAL COutMC
    REAL COutDPAir
    REAL COutErrFlag
END TYPE CondOutlet

TYPE(CondOutlet) :: CondOUT
!INTEGER, SAVE:: COutpRiC=1
!INTEGER, SAVE:: COuthRiC=2
!INTEGER, SAVE:: COuttAoC=3
!INTEGER, SAVE:: COutrhAoC=4
!INTEGER, SAVE:: COutpRoC=5
!INTEGER, SAVE:: COuthRoC=6
!INTEGER, SAVE:: COuttRoC=7
!INTEGER, SAVE:: COutWtAl=8
!INTEGER, SAVE:: COutWtCu=9
!INTEGER, SAVE:: COutpRiE=10
!INTEGER, SAVE:: COuthRiE=11
!INTEGER, SAVE:: COuttRiE=12
!INTEGER, SAVE:: COutxRiE=13
!INTEGER, SAVE:: COuttSCiE=14
!INTEGER, SAVE:: COutQC=15
!INTEGER, SAVE:: COutMDisLn=16
!INTEGER, SAVE:: COutMLiqLn=17
!INTEGER, SAVE:: COutMC=18
!INTEGER, SAVE:: COutDPAir=19
!INTEGER, SAVE:: COutErrFlag=20

!EvapPAR variables
TYPE EvapParameters
    REAL EvapSucLnLen
    REAL EvapSucLnOD
    REAL EvapSucLnElev
    REAL EvapSucLnQLoss
    REAL EvapSucLnTempChg
    REAL EvapSucLnAddPD
    REAL EvapSucLnTWThick
    REAL EvapCoilTOD
    REAL EvapCoilTWThick
    REAL EvapCoilSTLen
    REAL EvapCoilTThermCon
    REAL EvapTspc
    REAL EvapRspc
    REAL EvapFinThick
    REAL EvapFinPitch
    REAL EvapFinThermCon
    INTEGER EvapNt
    INTEGER EvapNl
    INTEGER EvapNumCkt
    INTEGER EvapCoolMode
    INTEGER EvapNumMod
    INTEGER EvapFinType
    REAL EvapMultRefQT
    REAL EvapMultRefPD
    REAL EvapMultAirQT
    REAL EvapMultAirPD
    REAL EvapFanPwr
    REAL EvapFanLoc
    REAL EvapSurfAbs
    REAL EvapTube
    REAL EvapBarPress
    REAL EvapCompQLoss
    REAL EvapPressTolConv
    INTEGER EvapSysType
    REAL EvapOilMassFrac
    INTEGER EvapCompMan
    INTEGER EvapSimpCoil
    INTEGER EvapFirstTime
END TYPE EvapParameters

TYPE(EvapParameters) :: EvapPAR
!INTEGER, SAVE:: EvapSucLnLen=1
!INTEGER, SAVE:: EvapSucLnOD=2
!INTEGER, SAVE:: EvapSucLnTWThick=3
!INTEGER, SAVE:: EvapSucLnElev=4
!INTEGER, SAVE:: EvapSucLnQLoss=5
!INTEGER, SAVE:: EvapSucLnTempChg=6
!INTEGER, SAVE:: EvapSucLnAddPD=7
!INTEGER, SAVE:: EvapCoilTOD=8
!INTEGER, SAVE:: EvapCoilTWThick=9
!INTEGER, SAVE:: EvapCoilSTLen=10
!INTEGER, SAVE:: EvapCoilTThermCon=11
!INTEGER, SAVE:: EvapTspc=12
!INTEGER, SAVE:: EvapRspc=13
!INTEGER, SAVE:: EvapFinThick=14
!INTEGER, SAVE:: EvapFinPitch=15
!INTEGER, SAVE:: EvapFinThermCon=16
!INTEGER, SAVE:: EvapNt=17
!INTEGER, SAVE:: EvapNl=18
!INTEGER, SAVE:: EvapNumCkt=19
!INTEGER, SAVE:: EvapCoolMode=20
!INTEGER, SAVE:: EvapNumMod=21
!INTEGER, SAVE:: EvapFinType=22
!INTEGER, SAVE:: EvapMultRefQT=23
!INTEGER, SAVE:: EvapMultRefPD=24
!INTEGER, SAVE:: EvapMultAirQT=25
!INTEGER, SAVE:: EvapMultAirPD=26
!INTEGER, SAVE:: EvapFanPwr=27
!INTEGER, SAVE:: EvapFanLoc=28
!INTEGER, SAVE:: EvapSurfAbs=29
!INTEGER, SAVE:: EvapTube=30
!INTEGER, SAVE:: EvapBarPress=31
!INTEGER, SAVE:: EvapCompQLoss=32
!INTEGER, SAVE:: EvapSysType=33
!INTEGER, SAVE:: EvapPressTolConv=34
!INTEGER, SAVE:: EvapOilMassFrac=35
!INTEGER, SAVE:: EvapCompMan=36
!INTEGER, SAVE:: EvapSimpCoil=37
!INTEGER, SAVE:: EvapFirstTime=38

!EvapIN variables
TYPE EvapInlet
    REAL EInmRef
    REAL EInpRi
    REAL EInhRi
    REAL EInmAi
    REAL EIntAi
    REAL EInrhAi
    REAL EInSolFlux
    REAL EIntRdis
END TYPE EvapInlet

TYPE(EvapInlet) :: EvapIN
!INTEGER, SAVE:: EInmRef=1
!INTEGER, SAVE:: EInpRi=2
!INTEGER, SAVE:: EInhRi=3
!INTEGER, SAVE:: EInmAi=4
!INTEGER, SAVE:: EIntAi=5
!INTEGER, SAVE:: EInrhAi=6
!INTEGER, SAVE:: EInSolFlux=8
!INTEGER, SAVE:: EIntRdis=9
!
!EvapOUT variables
TYPE EvapOutlet
    REAL EOutpRoC
    REAL EOuthRoC
    REAL EOuttAoC
    REAL EOutrhAoC
    REAL EOutDPAir
    REAL EOutpRiC
    REAL EOuthRiC
    REAL EOuttRiC
    REAL EOutxRiC
    REAL EOuttSHiC
    REAL EOutQC
    REAL EOutQCSens
    REAL EOutMSucLn
    REAL EOutMC
    REAL EOutWtAl
    REAL EOutWtCu
    REAL EOutErrFlag
    REAL AirEnth    !RS: Debugging: Setting to pass out the exiting air enthalpy in the output file (10/29/19)
    REAL AirHumRat  !RS: Debugging: Setting to pass out the exiting air humidity ratio in the output file (10/29/19)
END TYPE EvapOutlet

TYPE(EvapOutlet) :: EvapOUT
!INTEGER, SAVE:: EOutpRoC=1
!INTEGER, SAVE:: EOuthRoC=2
!INTEGER, SAVE:: EOuttAoC=3
!INTEGER, SAVE:: EOutrhAoC=4
!INTEGER, SAVE:: EOutDPAir=5
!INTEGER, SAVE:: EOutpRiC=6
!INTEGER, SAVE:: EOuthRiC=7
!INTEGER, SAVE:: EOuttRiC=8
!INTEGER, SAVE:: EOutxRiC=9
!INTEGER, SAVE:: EOuttSHiC=10
!INTEGER, SAVE:: EOutQC=11
!INTEGER, SAVE:: EOutQCSens=12
!INTEGER, SAVE:: EOutMSucLn=13
!INTEGER, SAVE:: EOutMC=14
!INTEGER, SAVE:: EOutWtAl=15
!INTEGER, SAVE:: EOutWtCu=16
!INTEGER, SAVE:: EOutErrFlag=17

!CompPAR variables
TYPE CompParameters
    REAL CompCoeffP1
    REAL CompCoeffP2
    REAL CompCoeffP3
    REAL CompCoeffP4
    REAL CompCoeffP5
    REAL CompCoeffP6
    REAL CompCoeffP7
    REAL CompCoeffP8
    REAL CompCoeffP9
    REAL CompCoeffP10
    REAL CompCoeffM1
    REAL CompCoeffM2
    REAL CompCoeffM3
    REAL CompCoeffM4
    REAL CompCoeffM5
    REAL CompCoeffM6
    REAL CompCoeffM7
    REAL CompCoeffM8
    REAL CompCoeffM9
    REAL CompCoeffM10
    REAL CompQLossFrac
    REAL CompQLoss
    REAL CompIntVol
    !REAL CompPwrCor !RS: Debugging: This is never set; used once in Compressor 262
    REAL CompPWrMult
    REAL CompMFRMult
END TYPE CompParameters

TYPE(CompParameters) :: CompPAR
!INTEGER, SAVE:: CompCoeffP1=1
!INTEGER, SAVE:: CompCoeffP2=2
!INTEGER, SAVE:: CompCoeffP3=3
!INTEGER, SAVE:: CompCoeffP4=4
!INTEGER, SAVE:: CompCoeffP5=5
!INTEGER, SAVE:: CompCoeffP6=6
!INTEGER, SAVE:: CompCoeffP7=7
!INTEGER, SAVE:: CompCoeffP8=8
!INTEGER, SAVE:: CompCoeffP9=9
!INTEGER, SAVE:: CompCoeffP10=10
!INTEGER, SAVE:: CompCoeffM1=11
!INTEGER, SAVE:: CompCoeffM2=12
!INTEGER, SAVE:: CompCoeffM3=13
!INTEGER, SAVE:: CompCoeffM4=14
!INTEGER, SAVE:: CompCoeffM5=15
!INTEGER, SAVE:: CompCoeffM6=16
!INTEGER, SAVE:: CompCoeffM7=17
!INTEGER, SAVE:: CompCoeffM8=18
!INTEGER, SAVE:: CompCoeffM9=19
!INTEGER, SAVE:: CompCoeffM10=20
!INTEGER, SAVE:: CompQLossFrac=21
!INTEGER, SAVE:: CompQLoss=22
!INTEGER, SAVE:: CompIntVol=23
!INTEGER, SAVE:: CompPwrCor=24
!INTEGER, SAVE:: CompPwrMult=25
!INTEGER, SAVE:: CompMFRMult=26

!CompIn variables
TYPE CompInlet
    REAL CompInPSuc
    REAL CompInPdis
    REAL CompInHsuc
END TYPE CompInlet

TYPE(CompInlet) :: CompIn
!INTEGER, SAVE:: CompInPsuc=1
!INTEGER, SAVE:: CompInPdis=2
!INTEGER, SAVE:: CompInHsuc=3

!CompOUT variables
TYPE CompOutlet
    REAL CmpOPwr
    REAL CmpOMdot
    REAL CmpOHdis
    REAL CmpOTdis
    REAL CmpOMcmp
    REAL CmpOErrFlag
END TYPE CompOutlet

TYPE(CompOutlet) :: CompOUT
!INTEGER, SAVE:: CmpOPwr=1
!INTEGER, SAVE:: CmpOMdot=2
!INTEGER, SAVE:: CmpOHdis=3
!INTEGER, SAVE:: CmpOTdis=5
!INTEGER, SAVE:: CmpOMCmp=6
!INTEGER, SAVE:: CmpOErrFlag=7

!CapTube PAR variables
TYPE CapTubeParameters
    REAL CTTubeID
    REAL CTTubeLen
    REAL CTTubeCoilD
    INTEGER CTEvapCktNum
    REAL CTDisTubeLen
END TYPE CapTubeParameters

TYPE(CapTubeParameters) :: CapTubePAR
!INTEGER, SAVE:: CTTubeID=1
!INTEGER, SAVE:: CTTubeLen=2
!INTEGER, SAVE:: CTTubeCoilD=3
!INTEGER, SAVE:: CTEvapCktNum=4
!INTEGER, SAVE:: CTDisTubeLen=5

!CapTube IN variables
TYPE CapTubeInlet
    REAL CTIMdot
    REAL CTIPiEx
    REAL CTIHiEx
    REAL CTIPiEv
    REAL CTIPoEv
END TYPE CapTubeInlet

TYPE(CapTubeInlet) :: CapTubeIN
!INTEGER, SAVE:: CTIMdot=1
!INTEGER, SAVE:: CTIPiEx=2
!INTEGER, SAVE:: CTIHiEx=3
!INTEGER, SAVE:: CTIPiEv=4
!INTEGER, SAVE:: CTIPoEv=5

!CapTube OUT variables
TYPE CapTubeOutlet
    REAL CTOMdot
    INTEGER CTOErrFlag
    REAL CTOToE
    REAL CTOXoE
    REAL CTOMDT
    REAL CTOPoE
END TYPE CapTubeOutlet

TYPE(CapTubeOutlet) :: CapTubeOUT
!INTEGER, SAVE:: CTOMdot=1
!INTEGER, SAVE:: CTOErrFlag=2
!INTEGER, SAVE:: CTOToE=3
!INTEGER, SAVE:: CTOXoE=4
!INTEGER, SAVE:: CTOMDT=5
!INTEGER, SAVE:: CTOPoE=6

!ShortTube PAR variables
TYPE ShTbParameters
    REAL ShTbTLen
    REAL ShTbTID
    REAL ShTbChamDep
    INTEGER ShTbECktNum
    REAL ShTbDTubeLen
END TYPE ShTbParameters

TYPE(ShTbParameters) :: ShTbPAR
!INTEGER, SAVE:: ShTbTLen=1
!INTEGER, SAVE:: ShTbTID=2
!INTEGER, SAVE:: ShTbChamDep=3
!INTEGER, SAVE:: ShTbECktNum=4
!INTEGER, SAVE:: ShTbDTubeLen=5

!ShortTube IN variables
TYPE ShTbInlet
    REAL ShTbINMdotC
    REAL ShTbINPiE
    REAL ShTbINHiE
    REAL ShTbINPiEv
    REAL ShTbINPoEv
END TYPE ShTbInlet

TYPE(ShTbInlet) :: ShTbIN
!INTEGER, SAVE:: ShTbINMdotC=1
!INTEGER, SAVE:: ShTbINPiE=2
!INTEGER, SAVE:: ShTbINHiE=3
!INTEGER, SAVE:: ShTbINPiEv=4
!INTEGER, SAVE:: ShTbINPoEv=5

!ShortTube OUT variables
TYPE ShTbOutlet
    REAL ShTbOMdotE
    REAL ShTbOPoE
    REAL ShTbOToE
    REAL ShTbOXoE
    REAL ShTbOMDT
    INTEGER ShTbOErrFlag
END TYPE ShTbOutlet

TYPE(ShTbOutlet) :: ShTbOUT
!INTEGER, SAVE:: ShTbOMdotE=1
!INTEGER, SAVE:: ShTbOPoE=2
!INTEGER, SAVE:: ShTbOToE=3
!INTEGER, SAVE:: ShTbOXoE=4
!INTEGER, SAVE:: ShTbOMDT=5
!INTEGER, SAVE:: ShTbOErrFlag=6

!AirProp variables
TYPE AirProperties
    REAL APTDB
    REAL APHumRat
    REAL APRelHum
    REAL APEnth
    REAL APTWB
    REAL APTDP
    REAL APDryDens
    REAL APWetDens
END TYPE AirProperties

TYPE(AirProperties) :: AirProp
!INTEGER, SAVE:: APTDB=1
!INTEGER, SAVE:: APHumRat=2
!INTEGER, SAVE:: APRelHum=3
!INTEGER, SAVE:: APEnth=4
!INTEGER, SAVE:: APTWB=5
!INTEGER, SAVE:: APTDP=6
!INTEGER, SAVE:: APDryDens=7
!INTEGER, SAVE:: APWetDens=8

!Accumulator PAR variables
TYPE AccumParameters
    REAL AccD
    REAL AccH
    REAL AccD1
    REAL AccD2
    REAL AccHDis
    REAL AccDTube
    REAL AccDP
    REAL AccDT
    REAL AccCM
    REAL AccCB
END TYPE AccumParameters

TYPE(AccumParameters) :: AccumPAR
!INTEGER, SAVE:: AccD=1
!INTEGER, SAVE:: AccH=2
!INTEGER, SAVE:: AccD1=3
!INTEGER, SAVE:: AccD2=4
!INTEGER, SAVE:: AccHDis=5
!INTEGER, SAVE:: AccDTube=6
!INTEGER, SAVE:: AccDP=7
!INTEGER, SAVE:: AccDT=8
!INTEGER, SAVE:: AccCM=9
!INTEGER, SAVE:: AccCB=10

!Accumulator IN variables
TYPE AccumInlet
    REAL AccImdot
    REAL AccIpRo
    REAL AccIhRo
END TYPE AccumInlet

TYPE(AccumInlet) :: AccumIN
!INTEGER, SAVE:: AccImdot=1
!INTEGER, SAVE:: AccIpRo=2
!INTEGER, SAVE:: AccIhRo=3

!Accumulator OUT variables
TYPE AccumOutlet
    REAL AccOMass
    REAL AccODP
END TYPE AccumOutlet

TYPE(AccumOutlet) :: AccumOUT
!INTEGER, SAVE:: AccOMass=1
!INTEGER, SAVE:: AccODP=2

!Filter PAR variables
TYPE FilterParameters
    REAL FilFlowCap
    REAL FilRatDP
END TYPE FilterParameters

TYPE(FilterParameters) :: FilterPAR
!INTEGER, SAVE:: FilFlowCap=1
!INTEGER, SAVE:: FilRatDP=2
!
!Filter IN variable
TYPE FilterInlet
    REAL FIDP
END TYPE FilterInlet

TYPE(FilterInlet) :: FilterIN
!INTEGER, SAVE:: FIDP=1

!Filter OUT variable
TYPE FilterOutlet
    REAL FODP
END TYPE FilterOutlet

TYPE(FilterOutlet) :: FilterOUT
!INTEGER, SAVE:: FODP=1

!TXV PAR variable
TYPE TXVOutput !Parameters  !RS: Changing from parameters to output because that's what TXVQ is
    REAL TXVQ
END TYPE TXVOutput !Parameters

TYPE(TXVOutput) :: TXVOUT !RS: Changing this from TXVPar to TVXOUT
INTEGER, SAVE:: TXVQ=1

TYPE FanData
    REAL HumRat
    REAL RhoAir
    REAL MotorEff
    REAL FanEff
    REAL DeltaPress
    REAL MotInAirFrac
    REAL Power
END TYPE FanData

TYPE(FanData) :: FanOut
TYPE(FanData) :: FanOutE

REAL, SAVE :: inputratio    !RS: The compressor capacity ratio for use in the new compressor model (10/3/14)

END MODULE DataSimulation
