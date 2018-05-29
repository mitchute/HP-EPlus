! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module simulates the condenser in the heat pump cycle.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component is the condenser in the Heat Pump Vapor-Compression cycle.
! A description of the component is found at:
! http://en.wikipedia.org/wiki/Condenser_%28heat_transfer%29
! From that website:
!  - This component takes refrigerant from a gas to a liquid state

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module takes the properties of the air entering the condenser and the evaporator properties
! and ouputs the leaving air properties and heat transfer rates to the fluid.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! Condenser.csv is an associated output file.
! The IDF file is used to read in some condenser properties.
!
! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There are a number of variables defined at the module level.
! These include convergence parameters, error and location flags,
! circuiting variables, refrigerant property variables,
! coil and segment variables, and geometry properties.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 23 methods:
!   PUBLIC Condenser-- Models the condenser in the heat pump cycle
!       Called by FlowRateLoop.f90
!       Called by ORNLsolver.f90
!   PUBLIC MicrochannelCondenser --  A segment-by-segment microchannel condenser model
!       Called internally only
!   PUBLIC CalcCondenserInventory -- Calculates the refrigerant inventory of the condenser
!       Called by HPdesignMod.f90
!       Called by ORNLsolver.f90
!   PUBLIC PrintCondenserResult - Prints simulation result to output file "Condenser.csv"
!       Called by ORNLsolver.f90
!   PUBLIC EndCondenserCoil -- Cleans and frees allocated arrays
!       Called by ORNLsolver.f90
!   PRIVATE InitCondenserCoil -- Initializes condenser geometry and circuiting
!       Called internally only
!   PRIVATE RefrigerantParameters -- Converts three refrigerant parameters
!       Called internally only
!   PRIVATE DischargeLine -- Calculates discharge line outlet conditions
!       Called internally only
!   PRIVATE LiquidLine -- Calculates liquid line outlet conditions
!       Called internally only
!   PRIVATE LoadMicrochannelInputs -- Brings in the input data to the microchannel model
!       Called internally only
!   PRIVATE LoadMicrochannelOutputs -- Sends out the output data from the microchannel model
!       Called internally only
!   PRIVATE InitBoundaryConditions -- Initializes the segment boundary conditions
!       Called internally only
!   PRIVATE CalcCircuitRefInletConditions -- Calculates the circuit refrigerant inlet conditions
!       Called internally only
!   PRIVATE CalcSegmentRefInletConditions -- Calculates the segment inlet refrigerant pressure and enthalpy
!       Called internally only
!   PRIVATE CalcSegmentAirInletConditions -- Calculates the inlet air temp and relative humidity
!       Called internally only
!   PRIVATE CalcCoilSegment -- Calculates heat exchanger equations for a segment
!       Called internally only
!   PRIVATE CalcSegmentOutletConditions -- Calculates the segment outlet conditions
!       Called internally only
!   PRIVATE CalcTransitionSegment -- Calculates the heat transfer in the transition segment
!       Called internally only
!   PRIVATE FindTransitionBoundary -- Finds the transition boundary in the transition segment
!       Called internally only
!   PRIVATE CalcRefProperty -- Calculates the refrigerant properties
!       Called internally only
!   PRIVATE CalcSegmentRefOutletPressure -- Calculates the refrigerant outlet pressure for a segment
!       Called internally only
!   PRIVATE UpdateTubeDataFromCircuitData -- Updates tube data from circuit data
!       Called internally only
!   PRIVATE InitCondenserStructures -- Allocates and initializes various condenser arrays
!       Called internally only

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! Unit conversions may need to be changed.

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2013-12-18 | RAS | Filled out the header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Some more documentation would be helpful. The ability of the code to switch from
! outdoor circuit to indoor and vice-versa needs to be explored; it currently seems
! to require a manual changing of certain inputs and outputs for the switch.

MODULE CondenserMod

  USE DataGlobals_HPSimIntegrated, ONLY: RefName, RefrigIndex    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  USE DataSimulation, ONLY: IsCoolingMode
  USE CoilCalcMod

  IMPLICIT NONE

  PRIVATE

  !Parameters
  INTEGER,PARAMETER  :: MdotMaxIter=10    !Max. number of iterations for mass flow rate calc.
  INTEGER,PARAMETER  :: RefBCmaxIter=20   !Max. number of iterations for refrigerant boundary condition calc.
  INTEGER,PARAMETER  :: AirBCmaxIter=20   !Max. number of iterations for air boundary condition calc.
  INTEGER,PARAMETER  :: PressureMaxIter=20        !Max. number of iterations
  REAL,PARAMETER :: SMALL=1.0E-4  !Small number
  REAL,PARAMETER :: BIG=1.0E20    !Big number
  REAL, PARAMETER :: Hout = 0.009937536 !Bare tube outside film heat transfer coefficient, kW/(m2-K)
  INTEGER,PARAMETER :: SI=1   !RS: Debugging: For unit conversion

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

  !System Types ISI - 07/13/06
  !INTEGER,PARAMETER :: ACUNIT        = 1
  !INTEGER,PARAMETER :: HEATPUMPUNIT  = 2
  !INTEGER,PARAMETER :: CONDENSERUNIT = 3
  !INTEGER,PARAMETER :: REHEATUNIT    = 4

  !Fan locations
  INTEGER,PARAMETER :: DRAWTHROUGH = 1
  INTEGER,PARAMETER :: BLOWTHROUGH = 2

  !Coil orientations
  INTEGER,PARAMETER :: HORIZONTAL = 1
  INTEGER,PARAMETER :: VERTICAL   = 2

  REAL,ALLOCATABLE,DIMENSION(:),SAVE :: mRefIter !Circuit flow rate for iteration check, kg/s
  INTEGER,ALLOCATABLE,DIMENSION(:),SAVE :: JoinTubes         !Joined tube numbers

  !Subcooling circuits variables, ISI - 06/05/07
  INTEGER :: SubcoolingTube !Subcooling tube number

  !Circuitry variables
  INTEGER I,J,K,II,III,IV !Loop control
  INTEGER NumOfTubes !Total number of tubes
  INTEGER TubeNum    !Tube number

  INTEGER ErrorFlag   !0-No error
  !1-Condenser solution not converge
  !2-Refprop error
  !3-Circuit file error
  !4,5-Coil size misdefined

  INTEGER FirstTube           !First simulation tube
  INTEGER LastTube            !Last simulation tube
  INTEGER EqCircuits          !1=Equivalent circuits; otherwise=no
  LOGICAL,SAVE :: IsUniformVelProfile !Is velocity profile uniform?

  !Refprop Table variable
  INTEGER :: RefID !1-R22; 2-R410A; 3-R407C; 4-R134a; 5-Propane; 6-R417A; 7-R509A
  REAL Temperature,Quality,Pressure,Enthalpy

  !Nomenclature:
  !R-Refrigerant;          A-Air;
  !i-Inlet;                o-Outlet;
  !f-Liquid phase;         g-Vapor phase;
  !t-Temperature(C);       p-Pressure(kPa);      x-Quality; h-Enthalpy(kJ/kg);
  !m-Mass flow rate(kg/s); rh-Relative humidity; v-Specific volume, m^3/kg;
  !sat-saturation;		 prev-previous iteration value

  !Variables for module
  REAL tRiMod,tRoMod,tRmod
  REAL pRiMod,pRoMod
  REAL xRmod,xRiMod,xRoMod
  REAL hRiMod,hRoMod
  REAL hfRiMod,hfRoMod
  REAL hgRiMod,hgRoMod
  REAL hfgRmod,hfgRiMod,hfgRoMod
  REAL vgRmod,vgRiMod,vgRoMod
  REAL vfRmod,vfRiMod,vfRoMod
  REAL vRiMod,vRoMod
  REAL muRmod,muRiMod,muRoMod
  REAL mugRmod,mugRiMod,mugRoMod
  REAL mufRmod,mufRiMod,mufRoMod
  REAL kRmod,kRiMod,kRoMod
  REAL kfRmod,kfRiMod,kfRoMod
  REAL kgRmod,kgRiMod,kgRoMod
  REAL cpRmod,cpRiMod,cpRoMod
  REAL cpfRmod,cpfRiMod,cpfRoMod
  REAL cpgRmod,cpgRiMod,cpgRoMod
  REAL mAiMod
  REAL tAiMod,tAoMod,tAmod
  REAL rhAiMod,rhAoMod
  REAL wbAoMod
  REAL hAiMod,hAoMod

  REAL DPmod !Pressure drop in module, kPa
  REAL SigmaMod !Surface tension, N/m
  REAL DTmod   !Temperature difference between saturated vapor and wall, C
  REAL tRoCkt  !Circuit outlet temp. C
  REAL pRiCkt  !Circuit Inlet pressure, kPa
  REAL pRoCkt  !Circuit outlet pressure, kPa
  REAL SumpRoCkt !Sum of outlet circuit outlet pressure, kPa
  REAL hRoCkt  !Circuit outlet enthalpy, kJ/kg
  REAL xRoCkt  !Circuit outlet quality
  REAL tRoAVG  !Average of ref outlet temperature, C
  REAL tSiSUM  !Sum of inlet surface temperature, C
  REAL tSoSUM  !Sum of outlet surface temperature, C

  !Compressor variables
  REAL tRoCmp
  REAL pRoCmp
  REAL hRoCmp
  REAL CpgRoCmp

  !Expansion device variables
  REAL tRiExp
  REAL pRiExp
  REAL hRiExp
  REAL xRiExp

  !Heat transfer calc. variables
  REAL mRefTot  !Refrigerant mass flow rate, kg/s
  REAL mRefMod  !Module refrigerant mass flow rate, kg/s
  REAL mRefJoin !Join tube total mass flow rate, kg/s
  REAL Cmin     !Min. capacity rate, kW/C
  REAL DT       !Temperature difference, C
  REAL NTU      !Number of transfer unit
  REAL Cratio   !Ratio of min to max capacity rate
  REAL EPS      !Heat exchanger effectiveness
  REAL Qcoil    !Total coil heat transfer, kW
  REAL PrevQcoil !Previous value of total coil heat transfer, kW
  REAL DiffQcoil !Difference of Qcoil in iteration, kW
  REAL Qckt     !Circuit heat transfer, kW
  REAL Qmod     !Module heat transfer, kW
  REAL QmodPrev !Previous module heat transfer, kW
  REAL QmodTP   !Heat transfer in two-phase region, kW
  REAL QmodSH   !Heat transfer in superheated region, kW
  REAL cAir     !Capacity rate of air, kW/C
  REAL cRef     !Capacity rate of refrigerant, kW/C
  REAL UA       !Overall heat transfer coefficient, kW/C
  REAL hcRef    !Refrigerant film coefficent, W/m^2-C
  REAL Rtube    !Thermal resistance of tube, K/W
  REAL Rair     !Module air film resistance, K/W
  REAL Rrefrig  !Module refrigerant film resistance, K/W
  REAL hco      !Air side heat tranfer coefficient, kW/m2-K
  REAL hci      !Refrigerant side heat tranfer coefficient, kW/m2-K
  REAL hcoMod   !Air side heat tranfer coefficient, kW/m2-K
  REAL hciMod   !Refrigerant side heat tranfer coefficient, kW/m2-K
  REAL EFref    !Refrigerant side heat tranfer enhancement factor
  REAL Velavg   !Average face velocity, m/s
  REAL ReVap    !Module Reynolds number vapor
  REAL ReLiq    !Module Reynolds number liquid
  REAL Const    !A constant
  REAL MolWeight !Molecular weight, kg/kmol
  REAL tSat      !Saturation temp., C
  REAL hciMultiplier   !Multiplier for hci
  REAL hcoMultiplier   !Multiplier for hco
  REAL DPrefMultiplier !Multiplier for DPref
  REAL DPairMultiplier !Multiplier for DPair
  REAL DPfric     !Frictional pressure drop, kPa
  REAL DPgrav     !Gravitational pressure drop, kPa
  REAL DPmom      !Momentum pressure drop, kPa
  REAL FaceVel    !Face velocity, m/s
  REAL DPair      !Air side pressure drop, kPa
  REAL SurfAbsorptivity !Surface absorptivity
  REAL SolarFlux  !Solar heat flux, kW/m2
  REAL QlossCmp   !Compressor heat loss, kW
  REAL IsCmpInAirStream !Is compressor in air stream, 1=yes, 0=no
  !INTEGER(2) SystemType !1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat !ISI - 07/14/06
  INTEGER,SAVE :: CompManufacturer !Compressor manufacturer: 1=Copeland
  !2=Bristol
  !3=Danfoss
  !4=Panasonic

  REAL AddDPDisLn !Discharge line additional pressure drop, kPa

  REAL AddDPLiqLn !Liquid line additional pressure drop, kPa

  !Properties
  REAL mu       !Bulk viscosity, Pa-s
  REAL muf      !Liquid viscosity, Pa-s
  REAL mug      !Vapor viscosity, Pa-s
  REAL kRef     !Refrigerant bulk conductivity, kW/m-K
  REAL cpRef    !Ref. specific heat, kJ/(kg-K)
  REAL rhoRef   !Ref. density, kg/m3
  REAL CPAir    !Specific heat of air, kJ/kg-K
  REAL DensityIn   !Inlet air density, kg/m3
  REAL DensityOut  !Outlet air density, kg/m3
  REAL Wabsolute   !Asolute oil mass fraction
  INTEGER(2) AirPropOpt     !Air prop calc. option
  INTEGER(2) AirPropErr     !Error flag:1-error; 0-no error
  !REAL AirProp(8)

  !Variables for coil
  REAL mRiCoil
  REAL mAiCoil
  REAL tRiCoil,tRoCoil
  REAL pRiCoil,pRoCoil
  REAL pRoCoilPrev !previous iteration value
  REAL hRiCoil,hRoCoil
  REAL xRiCoil,xRoCoil
  REAL tAiCoil
  REAL,SAVE :: tAoCoil
  REAL tSiCoil,tSoCoil
  REAL rhAiCoil,rhAoCoil
  REAL wbAiCoil
  REAL hAiCoil,hAoCoil
  REAL tSCoCoil     !Coil outlet subcooling, C
  REAL tSCiExp      !Exp.device inlet subcooling, C

  !Geometry variables
  REAL Aface       !Coil face area
  REAL,SAVE :: AiCoil !Inside coil surface area, m^2
  REAL,SAVE :: AoCoil !Outside coil surface area, m^2
  REAL,SAVE :: AfCoil !Coil fin surface area, m^2
  REAL,SAVE :: AmCoil !Coil tube mean surface area, m^2
  REAL AoMod       !Module outside surface area
  REAL AfMod       !Module fin surface area
  REAL AiMod       !Module inside surface area
  REAL AiModLiq    !Module inside surface area for liquid line
  REAL AiModDis    !Module inside surface area for discharge line
  REAL AmMod       !Module tube mean surface area
  REAL Lcoil       !Total tube length, m
  REAL LmodTube    !Module length of tube, m
  REAL,SAVE :: LmodTP !Two-phase module length, m
  REAL LmodTPmin   !Minimum two-phase module length, m
  REAL LmodTPmax   !Maximum two-phase module length, m
  REAL LmodTPratio !Ratio of two-phase module length to total module length
  REAL,SAVE :: LmodSH !Superheated module length, m
  REAL LmodSHmin   !Minimum Superheated module length, m
  REAL LmodSHmax   !Maximum Superheated module length, m
  REAL LmodSHratio !Ratio of Superheated module length to total module length

  REAL LmodDis     !Module length of discharge line, m
  REAL LdisLn      !Reversing valve to coil line length, m
  REAL ElevDisLn   !Discharge line elevation, m
  REAL IDdisLn     !Inside diameter of discharge line, m
  REAL ODdisLn     !Outside diameter of discharge line, m
  REAL DisLnThk    !Discharge line tube wall thickness, m

  REAL LmodLiq     !Module length of liquid line, m
  REAL LliqLn      !Liquid line length, m
  REAL ElevLiqLn   !Liquid line elevation, m
  REAL IDliqLn     !Inside diameter of liquid line, m
  REAL ODliqLn     !Outside diameter of liquid line, m
  REAL LiqLnThk    !Liquid line tube wall thickness, m

  REAL DreturnBend !Return bend diameter, m
  REAL LreturnBend !Return bend length, m

  REAL HtCoil      !Coil height, m
  REAL FinSpg      !Fin spacing, m
  REAL SurfEff     !Surface effecitiveness
  REAL FinEff      !Fin effecitiveness

  INTEGER,SAVE :: FinType       !1=Plain; 2=Wavy; 3=Louver; 4-11-element
  REAL,SAVE    :: FinPitch      !Fin pitch, fins/m
  REAL,SAVE    :: Kfin          !Fin thermal conductivity, kW/m-K
  REAL,SAVE    :: FinThk        !Fin thickness, m
  REAL,SAVE    :: FinHeight     !Fin height, m
  INTEGER,SAVE :: TubeType      !1=Plain; 2=General Micro Fin; 3=Herringbone; 4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
  REAL,SAVE    :: TubeHeight    !Tube height, m
  REAL,SAVE    :: TubeDepth     !Tube depth, m
  REAL,SAVE    :: TubeThk       !Coil tube wall thickness, m
  REAL,SAVE    :: Ktube         !Tube thermal conductivity, kW/m-K
  REAL,SAVE    :: Pl            !Tube spacing in longitudinal direction, m
  REAL,SAVE    :: Pt            !Tube spacing in lateral direction, m
  INTEGER,SAVE :: Nl            !Number of tubes in longitudinal direction
  INTEGER,SAVE :: Nt            !Number of tubes in traverse direction
  REAL,SAVE    :: ODtube        !Outside diameter of coil tube, m
  REAL,SAVE    :: IDtube        !Inside diameter of coil tube, m
  REAL,SAVE    :: Ltube         !Tube length, m
  INTEGER,SAVE :: NumOfMods	  !Number of modules per tube
  INTEGER,SAVE :: NumOfChannels !Number cf channels
  REAL,SAVE    :: Dchannel      !Channel diameter, m
  INTEGER,SAVE :: NumOfCkts     !Number of circuits
  INTEGER,SAVE :: ShiftTube     !1= last row lower than 2nd last row
  !0= last row higher than 2nd last row
  INTEGER NmodLast      !Total number of modules in the last row
  INTEGER IsParallelSlabs !Parallel microchannel slabs (1=yes; 0=no)
  INTEGER RowNum        !Coil row number
  INTEGER Ntube         !Tube number !Loop counter
  INTEGER Nckt          !Circuit number !Loop counter
  INTEGER Nmod          !Module number !Loop counter
  INTEGER NcktLast      !Total number of outlet circuits
  INTEGER NcktFirst     !Total number of inlet circuits
  INTEGER Nnode         !Number of split and joint nodes
  LOGICAL,SAVE :: IsSameNumOfTubes !Flag to check if same number of tubes
  !in all circuit branches
  REAL DrawBlow  !Fan location, 1=draw through; 2=blow through
  REAL PwrFan	   !Fan power, kW
  INTEGER WetFlag !1=Wet; 0=dry

  INTEGER(2)       :: RefPropErr  !Error flag:1-error; 0-no error
  REAL Psat,Pcr,Tcr

  !Mass inventory
  REAL MassDisLn   !Total refrigerant inventory in discharge line, kg
  REAL MassLiqLn   !Total refrigerant inventory in liquid line, kg
  REAL MassMod     !Refrigerant inventory in a module, kg
  REAL MassLiqMod  !Mass in liquid phase, kg
  REAL MassVapMod  !Mass in vapor phase, kg

  REAL, SAVE :: WeightAluminum !Weight of aluminum, kg
  REAL, SAVE :: WeightCopper   !Weight of copper, kg

  REAL TestH  !RS: Debugging: Air Entering enthalpy, hopefully.

  INTEGER FirstTime !Flag to indicate the first time of execution
  !1=yes, otherwise=no
  INTEGER Counter   !Iteration loop counter

  INTEGER IsSimpleCoil !Flag to indicate if it is simple coil, i.e. ignoring circuiting
  !1=Simple coil
  !otherwise=detailed
  INTEGER NumOfSections !Number of sections, ISI - 09/10/07

  INTEGER,SAVE :: CoilType  !1=Condenser; 2=Evaporator;           !Indicates what kind of set-up it is
  !3=High side interconnecting pipes;
  !4=Low side interconnecting pipes
  !5=Microchannel condenser
  !6=Microchannel evaporator

  TYPE (SlabInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Slab      !Coil slab pointer
  TYPE (CktInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Ckt        !Circuit pointer
  TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Tube      !Tube pointer
  TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:,:),SAVE :: Tube2D  !2-dimensional Tube pointer
  TYPE (ModInfo),ALLOCATABLE,DIMENSION(:),SAVE :: DisLnSeg   !Discharge line pointer
  TYPE (ModInfo),ALLOCATABLE,DIMENSION(:),SAVE :: LiqLnSeg   !Liquid line pointer
  TYPE (NodeInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Node      !Split or joint node
  TYPE (SectionInfo),ALLOCATABLE,DIMENSION(:) :: CoilSection !Coil section, ISI - 09/10/07

  PUBLIC  Condenser
  PUBLIC  MicrochannelCondenser
  PUBLIC  CalcCondenserInventory
  PUBLIC  PrintCondenserResult
  PUBLIC  EndCondenserCoil
  PRIVATE InitCondenserCoil
  PRIVATE RefrigerantParameters
  PRIVATE LoadMicrochannelInputs
  PRIVATE LoadMicrochannelOutputs
  PRIVATE DischargeLine
  PRIVATE LiquidLine
  PRIVATE InitBoundaryConditions
  PRIVATE CalcCircuitRefInletConditions
  PRIVATE CalcSegmentRefInletConditions
  PRIVATE CalcSegmentAirInletConditions
  PRIVATE CalcCoilSegment
  PRIVATE CalcSegmentOutletConditions
  PRIVATE CalcTransitionSegment
  PRIVATE FindTransitionBoundary
  PRIVATE CalcRefProperty
  PRIVATE CalcSegmentRefOutletPressure
  PRIVATE UpdateTubeDataFromCircuitData
  PUBLIC GetCondProp  !RS: Integration: Trying to carry over the properties to output
  PRIVATE InitCondenserStructures !RS: Debugging:
  PRIVATE PsyTdbFnHWLocal !RS: Debugging: Adding in condenser fan model (9/4/14)
  PRIVATE SimSimpleFan  !RS: Debugging: Adding in condenser fan model (9/4/14)
  PRIVATE PsyRhoAirFnPbTdbWLocal !RS: Debugging: Adding in condenser fan model (9/4/14)
  PRIVATE PsyWFnTdbHLocal !RS: Debugging: Adding in condenser fan model (9/4/14)

CONTAINS

  !***********************************************************************************

  SUBROUTINE Condenser(Ref$) !(Ref$,PureRef,XIN,PAR,OUT)  !RS: Debugging: Extraneous PureRef

    !-----------------------------------------------------------------------------------
    !
    !  Description:
    !  Ragazzi's modular coil model (Fixed length version)
    !  To predict coil air side and refrigerant side properties, heat transfer,
    !  and prssure drop
    !
    !  Inputs:
    !  Ref$=Refrigerant name
    !  PureRef=Refrigerant flag: 1=pure refrigerant
    !                            0=refrigerant mixture
    !  XIN(1)=Refrigerant side mass flow rate, kg/s
    !  XIN(2)=Refrigerant side inlet (compressor outlet) pressure, kPa
    !  XIN(3)=Refrigerant side inlet (compressor outlet) enthalpy, kJ/kg
    !  XIN(4)=Air side mass flow rate, kg/s
    !  XIN(5)=Air side inlet temp. C
    !  XIN(6)=Air side inlet relative humidity
    !  XIN(7)=Solar heat flux, kW/m^2
    !  XIN(8)=Evaporator outlet temperature, C
    !  XIN(9)=Evaporator entering air temperature, C
    !
    !  Parameters:
    !  PAR(1)=Discharge line length, m
    !  PAR(2)=Discharge line outside diameter, m
    !  PAR(3)=Discharge line tube wall thickness, m
    !  PAR(4)=Discharge line elevation, m
    !  PAR(5)=Discharge line heat loss, kW
    !  PAR(6)=Discharge line temperature change, C
    !  PAR(7)=Discharge line additional pressure drop, kPa
    !  PAR(8)=Liquid line length, m
    !  PAR(9)=Liquid line outside diameter, m
    !  PAR(10)=Liquid line tube wall thickness, m
    !  PAR(11)=Liquid line elevation, m
    !  PAR(12)=Liquid line heat loss, kW
    !  PAR(13)=Liquid line temperature change, C
    !  PAR(14)=Liquid line additional pressure drop, kPa
    !  PAR(15)=Coil tube outside diameter, m
    !  PAR(16)=Coil tube wall thickness, m
    !  PAR(17)=Coil single tube length, m
    !  PAR(18)=Coil tube thermal conductivity, kW/m-C
    !  PAR(19)=Tube spacing in transverse direction, m (normal to air flow)
    !  PAR(20)=Row spacing in longitudinal direction, m (parallel to air flow)
    !  PAR(21)=Fin thickness, m
    !  PAR(22)=Fin pitch, fin/m
    !  PAR(23)=Fin thermal conductivity, kW/m-C
    !  PAR(24)=Number of tubes in transverse direction (normal to air flow)
    !  PAR(25)=Number of rows in longitudinal direction (parallel to air flow)
    !  PAR(26)=Number of circuits
    !  PAR(27)=Cooling mode? 1=yes; 0=no
    !  PAR(28)=Number of modules per tube
    !  PAR(29)=Fin type: 1=Plain; 2=Wavy; 3=Louver
    !  PAR(30)=Multiplier for ref. side heat transfer correlation
    !  PAR(31)=Multiplier for ref. side pressure drop correlation
    !  PAR(32)=Multiplier for air side heat transfer correlation
    !  PAR(33)=Multiplier for air side pressure drop correlation
    !  PAR(34)=Fan power, kW
    !  PAR(35)=Fan location, 1=draw through; 2=blow through
    !  PAR(36)=Surface absorptivity
    !  PAR(37)=Tube tube: 1-Smooth; 2-Microfin; 3=Herringbone; 4=Crosshatch;
    !                     5=Herringbone w/crosshatch; 6=Rifle; 7=Helical;
    !                     8=42F HXH
    !  PAR(38)=Barometric pressure, kPa
    !  PAR(39)=Compressor heat loss, kW
    !  PAR(40)=Pressure tolerance convergence Criteria, kPa
    !  PAR(41)=System type !1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat
    !  PAR(42)=Oil mass fraction
    !  PAR(43)=Compressor manufacturer: 1=Copeland; 2=Bristol;
    !                                   3=Danfoss;  4=Panasonic
    !  PAR(44)=Simple coil flag: 1=Simple coil; otherwise=Detailed coil
    !  PAR(45)=First time to run this model flag: 1=yes, otherwise=no
    !          for component validation, set it to 1
    !          for system validation, set it to 1 first, then zero
    !
    !  Outputs:
    !  OUT(1)=Coil inlet pressure, kPa
    !  OUT(2)=Coil inlet enthalpy, kJ/kg
    !  OUT(3)=Coil inlet temperature, C
    !  OUT(4)=Coil inlet quality
    !  OUT(5)=Coil outlet pressure, kPa
    !  OUT(6)=Coil outlet enthalpy, kJ/kg
    !  OUT(7)=Coil outlet temperature, C
    !  OUT(8)=Coil outlet quality
    !  OUT(9)=Coil outlet subcooling, C
    !  OUT(10)=Liquid line outlet pressure, kPa
    !  OUT(11)=Liquid line outlet enthalpy, kJ/kg
    !  OUT(12)=Liquid line outlet temperature, C
    !  OUT(13)=Liquid line outlet quality
    !  OUT(14)=Liquid line outlet subcooling, C
    !  OUT(15)=Coil capacity, kW
    !  OUT(16)=Mass in discharge line, kg
    !  OUT(17)=Mass in liquid line, kg
    !  OUT(18)=Mass in coil, kg
    !  OUT(19)=Aluminum weight, kg
    !  OUT(20)=Copper weight, kg
    !  OUT(21)=Air side outlet temperature, C
    !  OUT(22)=Air side outlet relative humidity
    !  OUT(23)=Air side pressure drop, kPa
    !  OUT(24)=Error flag: 0-No error
    !                      1-Condenser solution not converge
    !                      2-Refprop error
    !                      3-Circuit file error
    !
    !
    !  Reference:
    !  Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
    !    of an air-cooled condenser. ACRC technical report 07.
    !    University of Illinois, Urbana,IL.
    !
    !  Author:
    !  Ipseng Iu
    !  Mechanical and Aerospace Engineering
    !  Oklahoma State University, Stillwater
    !
    !  Date: June 2002
    !
    !-----------------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod
    USE AirPropMod
    USE OilMixtureMod
    USE InputProcessor    !RS: Debugging: GetObjectItem

    IMPLICIT NONE

    !Subroutine argument declarations
    CHARACTER*80,     INTENT(IN)  :: Ref$
    !INTEGER(2),       INTENT(IN)  :: PureRef   !RS: Debugging: Extraneous
    !REAL, INTENT(IN)  :: XIN(7) !RS: Debugging: Formerly XIN(9)
    !REAL, INTENT(IN)  :: PAR(45) !ISI - 12/21/06 !RS: Debugging: Formerly PAR(62)
    !REAL, INTENT(OUT) :: OUT(20)    !RS: Debugging: Formerly OUT(29), OUT(24)

    !Subroutine lcoal variables
    REAL :: MCXIN(7)  !Microchannel coil input data
    REAL :: MCPAR(39) !Microchannel coil input parameters
    REAL :: MCOUT(22) !Microchannel coil output data

    !INTEGER,SAVE :: CoilType  !1=Condenser; 2=Evaporator;
    !                          !3=High side interconnecting pipes;
    !                          !4=Low side interconnecting pipes
    !                          !5=Microchannel condenser
    !                          !6=Microchannel evaporator

    INTEGER Iter                  !Iteration loop counter
    INTEGER AirBCiter             !Iteration loop counter
    LOGICAL Converged             !Solution convergence flag
    REAL MaxResidual !Maximum residual in iteration
    REAL PTol !Condenser Outlet Pressure Convergence criteria, kPa

    !ISI - 07/14/06

    INTEGER, PARAMETER :: MaxNameLength = 200
    REAL, PARAMETER :: UnitP     = 6.8947453 !(psi X UnitP = kPa)
    REAL, PARAMETER :: UnitM     = 0.4536    !(lbm X UnitM = kg)
    REAL, PARAMETER :: UnitL     = 0.3048    !(ft X UnitL = m)
    REAL, PARAMETER :: UnitK     = 0.1442E-3 !(Btu-in/hr-ft2-F X UnitK = kW/m-C)

    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL :: LmodTubeKeep    !RS: Debugging: Placeholder for LmodTube that keeps the previous one

    CHARACTER(LEN=11),PARAMETER :: FMT_107 = "(A66,F10.3)"

    mRefTot =CondIN%CInmRef   !RS: Debugging: Formerly XIN(1)
    pRoCmp  =CondIN%CInpRo !RS: Debugging: Formerly XIN(2)
    hRoCmp  =CondIN%CInhRo !RS: Debugging: Formerly XIN(3)
    mAiCoil =CondIN%CInmAi !RS: Debugging: Formerly XIN(4)
    tAiCoil =CondIN%CIntAi !RS: Debugging: Formerly XIN(5)
    rhAiCoil=CondIN%CInrhAi !RS: Debugging: Formerly XIN(6)
    SolarFlux=CondIN%CInSolFlux    !RS: Debugging: Used once but was set equal to 0    !RS: Debugging: Formerly XIN(7)

    LdisLn    = CondPAR%CondDisLnLen  !RS: Debugging: Formerly PAR(1)
    ODdisLn   = CondPAR%CondDisLnOD  !RS: Debugging: Formerly PAR(2)
    DisLnThk  = CondPAR%CondDisLnTWThick  !RS: Debugging: Formerly PAR(3)
    ElevDisLn = CondPAR%CondDisLnElev  !RS: Debugging: Formerly PAR(4)
    QdisLn    = CondPAR%CondDisLnQLoss  !RS: Debugging: Formerly PAR(5)
    DTdisLn   = CondPAR%CondDisLnTempChg  !RS: Debugging: Formerly PAR(6)
    AddDPdisLn = CondPAR%CondDisLnAddPD !RS: Debugging: Formerly PAR(7)

    LliqLn    = CondPAR%CondLiqLnLen  !RS: Debugging: Formerly PAR(8)
    ODliqLn   = CondPAR%CondLiqLnOD  !RS: Debugging: Formerly PAR(9)
    LiqLnThk  = CondPAR%CondLiqLnTWThick !RS: Debugging: Formerly PAR(10)
    ElevLiqLn = CondPAR%CondLiqLnElev !RS: Debugging: Formerly PAR(11)
    QliqLn    = CondPAR%CondLiqLnQLoss !RS: Debugging: Formerly PAR(12)
    DTliqLn   = CondPAR%CondLiqLnTempChg !RS: Debugging: Formerly PAR(13)
    AddDPLiqLn = CondPAR%CondLiqLnAddPD    !RS: Debugging: Formerly PAR(14)

    IsSimpleCoil=CondPAR%CondSimpCoil !ISI - 12/22/06    !RS: Debugging: Formerly PAR(44) !RS: Debugging: IsSimple
    FirstTime=CondPAR%CondFirstTime    !ISI - 12/22/06    !RS: Debugging: Formerly PAR(45)

    !Initialize circuiting and refrigerant parameters
    IF (FirstTime .EQ. 1) THEN
      ODtube      = CondPAR%CondCoilTOD   !RS: Debugging: Formerly PAR(15)
      TubeThk     = CondPAR%CondCoilTWThick   !RS: Debugging: Formerly PAR(16)
      Ltube       = CondPAR%CondCoilSTLen   !RS: Debugging: Formerly PAR(17)
      Ktube       = CondPAR%CondCoilTThermCon   !RS: Debugging: Formerly PAR(18)
      Pt          = CondPAR%CondTspc   !RS: Debugging: Formerly PAR(19)
      Pl          = CondPAR%CondRspc   !RS: Debugging: Formerly PAR(20)
      FinThk      = CondPAR%CondFinThick   !RS: Debugging: Formerly PAR(21)
      FinPitch    = CondPAR%CondFinPitch   !RS: Debugging: Formerly PAR(22)
      Kfin        = CondPAR%CondFinThermCon   !RS: Debugging: Formerly PAR(23)
      Nt          = CondPAR%CondNt   !RS: Debugging: Formerly PAR(24)
      Nl          = CondPAR%CondNl   !RS: Debugging: Formerly PAR(25)
      !NumOfCkts   = CondPAR%CondNumCkt   !RS: Debugging: Formerly PAR(26) !These are set in InitCondenserCoil
      !NumOfMods   = CondPAR%CondNumMod   !RS: Debugging: Formerly PAR(28)
      FinType     = CondPAR%CondFinType   !RS: Debugging: Formerly PAR(29)
      TubeType    = CondPAR%CondTube   !RS: Debugging: Formerly PAR(37)
      CALL InitCondenserCoil(CoilType)
      CALL CalcMaterialWeight(CoilType,Ltube,IDtube,ODtube,TubeHeight,TubeDepth, &
      Dchannel,NumOfChannels,Pt,Pl,Nt,Nl,NumOfCkts, &
      FinThk,FinPitch,WeightAluminum,WeightCopper)
      IF (ErrorFlag .NE. NOERROR) THEN
        CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
        CALL Condenser_Helper_1
        RETURN
      END IF
      CALL RefrigerantParameters(Ref$)
      CALL GetRefID(Ref$,RefID)
      tAoCoil=tAiCoil !ISI - 05/27/2008


      !********************Refrigerant Cycle Data (Heating)***********************  !RS: Debugging: Moving: Stay here? Compressor? ORNLSolver?

      CALL GetObjectItem('RefrigerantCycleData(Heating)',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IsCmpInAirStream = Numbers(2) !Is Compressor in Air Stream

    END IF

    hciMultiplier   = CondPAR%CondMultRefQT   !RS: Debugging: Formerly PAR(30)
    DPrefMultiplier = CondPAR%CondMultRefPD   !RS: Debugging: Formerly PAR(31)
    hcoMultiplier   = CondPAR%CondMultAirQT   !RS: Debugging: Formerly PAR(32)
    DPairMultiplier = CondPAR%CondMultAirPD   !RS: Debugging: Formerly PAR(33)

    PwrFan           = CondPAR%CondFanPwr/1000  !RS: Debugging: Formerly PAR(34)    !RS: Debugging: Converting
    DrawBlow         = CondPAR%CondFanLoc  !RS: Debugging: Formerly PAR(35)
    SurfAbsorptivity = CondPAR%CondSurfAbs  !RS: Debugging: Formerly PAR(36)

    BaroPressure     = CondPAR%CondBarPress  !RS: Debugging: Formerly PAR(38)
    QlossCmp         = CondPAR%CondCompQLoss  !RS: Debugging: Formerly PAR(39)
    PTol             = CondPAR%CondPressTolConv !(56)  !RS: Debugging: Formerly PAR(40)
    SystemType       = CondPAR%CondSysType !ISI - 07/14/06  !RS: Debugging: Formerly PAR(41)
    Wabsolute        = CondPAR%CondOilMassFrac  !RS: Debugging: Formerly PAR(42)
    CompManufacturer = CondPAR%CondCompMan  !RS: Debugging: Formerly PAR(43)

    IF (CoilType .EQ. MCCONDENSER) THEN
      CALL LoadMicrochannelInputs(MCXIN,MCPAR) !(XIN,PAR,MCXIN,MCPAR)
      CALL MicrochannelCondenser(MCXIN,MCPAR,MCOUT) !(Ref$,MCXIN,MCPAR,MCOUT)    !RS: Debugging: Extraneous Ref$
      CALL LoadMicrochannelOutputs(MCOUT) !,OUT)
      RETURN
    END IF

    EqCircuits=0 !Equivalent circuit flag

    ErrorFlag=NOERROR !Initialize

    !Tube inside diameter
    IDdisLn=ODdisLn-DisLnThk*2
    IDliqLn=ODliqLn-LiqLnThk*2

    !Return bend length
    Dreturnbend=Pt
    Lreturnbend=Dreturnbend*PI/2

    !Total number of modules in the last row
    NmodLast=Nt*NumOfMods

    !Coil height
    HtCoil=Nt*Pt

    !Coil length
    Lcoil=Nl*Nt*Ltube

    !Face area
    Aface=Ltube*Nt*Pt

    !Tube information
    LmodTube=Ltube/NumOfMods

    !Discharge line info
    LmodDis=Ldisln
    AiModDis=PI*IDdisLn*LmodDis

    !liquid line info
    LmodLiq=Lliqln
    AiModLiq=PI*IDliqLn*LmodLiq

    CALL InitBoundaryConditions(CoilType)
    IF (ErrorFlag .NE. NOERROR) THEN
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    !****** Start coil calculation ******
    Converged=.TRUE.
    pRoCoilPrev=pRiCoil
    MaxResidual=0

    DO Iter=1, MdotMaxIter !RS: Debugging: Temporarily setting in an Epsilon-NTU method

      !Initialize
      mRefJoin=0
      PrevQcoil=BIG
      QmodPrev=0

      DO AirBCiter=1, AirBCmaxIter

        Qcoil=0.0;  tSiSUM=0.0; tSoSUM=0.0;

        DO I=1,NumOfCkts

          Qckt=0.0

          CALL CalcCircuitRefInletConditions(I,I,CoilType)

          !Find first and last simulation tubes
          IF (IsSimpleCoil .EQ. 1) THEN
            FirstTube=1
            LastTube=1
          ELSE
            FirstTube=1
            LastTube=Ckt(I)%Ntube
            IF (Ckt(I)%InSplit .GT. 1) THEN
              FirstTube=2 !Skip first tube
            END IF
            IF (Ckt(I)%OutJoin .GT. 1) THEN
              LastTube=Ckt(I)%Ntube-1 !Ignore last tube
            END IF
          END IF

          mRefMod=Ckt(I)%mRef

          DO J=FirstTube,LastTube

            IF (IsSimpleCoil .EQ. 1) THEN
              TubeNum=1
            ELSE
              TubeNum=Ckt(I)%TubeSequence(J)
            END IF

            DO K=1,NumOfMods

              IF (IsSimpleCoil .EQ. 1) THEN
                !SELECT CASE(K)
                !CASE (1)
                !    LmodTube=Lcoil/NumOfCkts !Start with guessing the whole length
                !CASE (2)
                !    LmodTube=Lcoil/NumOfCkts-Ckt(I)%Tube(J)%Seg(1)%Len
                !CASE (3)
                !    LmodTube=Lcoil/NumOfCkts-(Ckt(I)%Tube(J)%Seg(1)%Len+Ckt(I)%Tube(J)%Seg(2)%Len)
                !END SELECT
                IF (K .EQ. 1) THEN  !RS: Debugging: Handling every module
                  LmodTube=Lcoil/NumofCkts
                  LmodTubeKeep=LmodTube
                ELSE
                  LmodTube=LmodTubeKeep-Ckt(I)%Tube(J)%Seg(1)%Len
                  LmodTubeKeep=LmodTube
                END IF
              END IF

              CALL CalcCoilSegment(I,I,J,K,CoilType) !RS: Debugging: Temporarily setting in an Epsilon-NTU method

              IF (ErrorFlag .GT. CONVERGEERROR) THEN
                CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
                CALL Condenser_Helper_1
                RETURN
              END IF

              !Calc. circuit heat transfer
              Qckt=Qckt+Ckt(I)%Tube(J)%Seg(K)%Qmod

              !Calc. sum of surface air temperature
              IF (Ckt(I)%Tube(J)%Back .EQ. 0) THEN
                tSoSUM=tSoSUM+Ckt(I)%Tube(J)%Seg(K)%tSo
              END IF

              IF (Ckt(I)%Tube(J)%Fup .EQ. 0 .AND. Ckt(I)%Tube(J)%Fdown .EQ. 0) THEN
                tSiSUM=tSiSUM+Ckt(I)%Tube(J)%Seg(K)%tSi
              END IF

            END DO !End mod

          END DO !End tube

          pRoCkt=Ckt(I)%Tube(LastTube)%Seg(NumOfMods)%pRo !Circuit outlet pressure
          hRoCkt=Ckt(I)%Tube(LastTube)%Seg(NumOfMods)%hRo !Circuit outlet enthalpy

          Pressure=pRoCkt*1000    !RS Comment: Unit Conversion
          Enthalpy=hRoCkt*1000    !RS Comment: Unit Conversion
          tRoCkt=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !RS Comment: Circuit outlet temperature
          IF (RefPropErr .GT. 0) THEN
            WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
            ErrorFlag=REFPROPERROR
            CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
            CALL Condenser_Helper_1
            RETURN
          END IF
          xRoCkt=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !RS Comment: Circuit outlet quality
          IF (RefPropErr .GT. 0) THEN
            WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
            ErrorFlag=REFPROPERROR
            CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
            CALL Condenser_Helper_1
            RETURN
          END IF

          Pressure=pRoCkt*1000    !RS Comment: Unit Conversion
          Quality=0
          tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
          IF (RefPropErr .GT. 0) THEN
            WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
            ErrorFlag=REFPROPERROR
            CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
            CALL Condenser_Helper_1
            RETURN
          END IF

          IF (xRoCkt .LE. 0.0) THEN
            Ckt(I)%tSC=tSat-tRoCkt !Subcooling
          ELSE
            Ckt(I)%tSC=0.0
          END IF

          Ckt(I)%tRo=tRoCkt
          Ckt(I)%pRo=pRoCkt
          Ckt(I)%hRo=hRoCkt
          Ckt(I)%Qckt=Qckt        !Circuit capacity
          Qcoil=Qcoil+Ckt(I)%Qckt !Total coil capacity

          IF (EqCircuits .EQ. 1 .AND. IsUniformVelProfile .OR. IsSimpleCoil .EQ. 1) THEN  !Equivalent circuit and Uniform velocity profile
            Qcoil=Qcoil*NumOfCkts
            pRoCoil=pRoCkt
            hRoCoil=hRoCkt
            EXIT
          END IF

          IF (Ckt(I)%OutSplit .LE. 1 .AND. Ckt(I)%OutJoin .LE. 1) THEN
            mRefJoin=mRefJoin+Ckt(I)%mRef
          END IF

        END DO !End circuit

        IF (IsSimpleCoil .EQ. 1) THEN
          EXIT
        END IF

        CALL CalcMeanProp(tAiCoil,tAoCoil,tAmod)    !RS Comment: Mean Air Coil Temperature

        !CPair=CPA(REAL(tAmod))  !RS Comment: Finding the specific heat of air   !RS: Replace: CPA (2/19/14)
        CPair=CPAirFunction(tAMod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
        Cair=mAicoil*CPAir      !RS Comment: Finding the capacity rate of air

        tAoCoil=tAiCoil+Qcoil/Cair  !RS Comment: Air Coil Outlet Temperature

        DiffQcoil=ABS((Qcoil-PrevQcoil)/PrevQcoil)
        IF (DiffQcoil .GT. 1E-3) THEN
          PrevQcoil=Qcoil
        ELSE
          EXIT
        END IF

      END DO !end AirBCiter

      IF (IsSimpleCoil .EQ. 1) THEN
        EXIT
      END IF

      IF (AirBCiter .GT. AirBCmaxIter) THEN
        !AirBCiter not converged.
        ErrorFlag=CONVERGEERROR
      END IF

      IF (EqCircuits .EQ. 1 .AND. IsUniformVelProfile) THEN
        EXIT !for equivalent circuit and uniform velocity profile no need to update mdot ref.
      END IF

      !Synchronize from circuit array to tube array
      DO I=1, NumOfCkts

        FirstTube=1
        LastTube=Ckt(I)%Ntube
        IF (Ckt(I)%InSplit .GT. 1) THEN
          FirstTube=2 !Skip first tube
        END IF
        IF (Ckt(I)%OutJoin .GT. 1) THEN
          LastTube=Ckt(I)%Ntube-1 !Ignore last tube
        END IF

        DO J=FirstTube, LastTube
          TubeNum=Ckt(I)%TubeSequence(J)
          Tube(TubeNum)=Ckt(I)%Tube(J)    !RS Comment: Implementing Tube Numbers into the Circuit
        END DO

        IF (Ckt(I)%InSplit .GT. 1) THEN
          Ckt(I)%Tube(1)=Tube(Ckt(I)%TubeSequence(1))
        END IF
        IF (Ckt(I)%OutJoin .GT. 1) THEN
          Ckt(I)%Tube(Ckt(I)%Ntube)=Tube(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
        END IF

      END DO !End circuit

      !Calculate maximum residual
      pRoCkt=0
      SumpRoCkt=0
      DO I=1, NumOfCkts
        IF (Ckt(I)%OutSplit .EQ. 0) THEN !outlet circuit
          SumpRoCkt=SumpRoCkt+Ckt(I)%pRo  !RS Comment: Summing the outlet circuit refrigerant pressure
        END IF
      END DO !End Circuit
      IF (SumpRoCkt .EQ. 0) THEN
        SumpRoCkt=Ckt(1)%pRo !At least 1 circuit, ISI - 07/28/06
      END IF
      pRoCoil=SumpRoCkt/NcktLast

      IF (ABS(pRoCoil-pRoCoilPrev) .GT. Ptol) THEN
        MaxResidual=ABS(pRoCoilPrev-pRoCoil)
        pRoCoilPrev=pRoCoil
        Converged=.FALSE.   !No convergence
      END IF

      !RS: Debugging: The only reason I can think why this would be useful is if you calculated and converged the
      ! pRoCoil for one circuit and then saved and set all the rest of the circuits to the first's converged values.
      ! But it doesn't look like that's what's happening, so I've removed this for now to try to get a better set of
      ! outputs to match Iu's validation files. (11/4/13)
      IF (IsSameNumOfTubes .AND. IsUniformVelProfile) THEN   !RS: Debugging: I'm not sure why this would be useful
        EXIT
      END IF

      IF (NOT(Converged) .OR. Iter .LE. 2) THEN
        Converged=.TRUE. ! Reinitialize
        !Moved this subroutine to CoilCalc and share with evaporator ISI - 06/05/07
        CALL UpdateRefMassFlowRate(Iter,Ckt,NumOfCkts,pRiCoil,mRefTot,Nnode,Node)
      ELSE
        EXIT
      END IF

    END DO !End iter

    IF (Iter .GT. MdotMaxIter) THEN
      WRITE(*,FMT_107)'-- WARNING -- Condenser: Solution not converged. Max. Residual = ',MaxResidual
      ErrorFlag=CONVERGEERROR
    END IF

    !Surface temperature
    tSiCoil=tSiSUM/NmodLast
    tSoCoil=tSoSUM/NmodLast

    !Coil air side outlet conditions
    !CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAicoil*CPAir

    tAoCoil=tAiCoil+Qcoil/Cair
    hAoCoil=hAiCoil+Qcoil/mAiCoil

    !Fan air side inlet conditions
    !CPair=CPA(REAL(tAoCoil))    !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAiCoil*CPAir

    IF (DrawBlow .EQ. DRAWTHROUGH) THEN !Draw through
      !tAoCoil=tAoCoil+PwrFan/Cair
      !hAoCoil=hAoCoil+PwrFan/mAiCoil
      CALL SimSimpleFan(tAoCoil,hAoCoil,mAiCoil,tAoCoil,hAoCoil)  !RS: Debugging: Adding in condenser fan model (9/4/14)
      CondPAR%CondFanPwr=FanOut%Power
    END IF

    !RS Comment: Inlet and Outlet Air Properties
    AirPropOpt=1
    AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APEnth=hAiCoil  !RS: Debugging: Formerly AirProp(4)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
    rhAiCoil=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
    DensityIn=AirProp%APDryDens    !RS: Debugging: Formerly AirProp(7)

    AirPropOpt=1
    AirProp%APTDB=tAoCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APEnth=hAoCoil  !RS: Debugging: Formerly AirProp(4)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
    rhAoCoil=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
    DensityOut=AirProp%APDryDens   !RS: Debugging: Formerly AirProp(7)

    WetFlag=0
    RowNum=0
    CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
    IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair, &
    AirProp%APEnth)

    !CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
    !IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC, &
    !Poly1HTC,Poly2HTC,Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP, &
    !Poly1DP,Poly2DP,Poly3DP,Poly4DP,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair)

    DPair=DPair*DPairMultiplier

    hRoCoil=hRiCoil-Qcoil/mRefTot

    Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
    Enthalpy=hRoCoil*1000   !RS Comment: Unit Conversion
    tRoCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Circuit Outlet Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    xRoCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Circuit Outlet Quality
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
    Quality=0
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    IF (xRoCoil .LE. 0.0) THEN
      tSCoCoil=tSat-tRoCoil !Subcooling
    ELSE
      tSCoCoil=0.0
    END IF

    !Distributor pressure drop - ISI - 07/14/06
    !RS: The "IF" and "ELSE" sections are identical, and thus the statement is useless (12/16/13)
    !IF (SystemType .EQ. HEATPUMP) THEN !Heat Pump

    !****** Liquid line calculation ******
    IF (LliqLn .GT. 0) THEN
      CALL LiquidLine
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        WRITE(*,*)'LiquidLine: Refprop error.'
        CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
        CALL Condenser_Helper_1
        RETURN
      END IF
    ELSE
      pRiExp=pRoCoil !pRoCoilTemp
      hRiExp=hRoCoil
    END IF

    !ELSE
    !
    !    !****** Liquid line calculation ******
    !    IF (LliqLn .GT. 0) THEN
    !        CALL LiquidLine
    !        IF (ErrorFlag .GT. CONVERGEERROR) THEN
    !            WRITE(*,*)'LiquidLine: Refprop error.'
    !            CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
    !            CALL Condenser_Helper_1
    !            RETURN
    !        END IF
    !    ELSE
    !        pRiExp=pRoCoil
    !        hRiExp=hRoCoil
    !    END IF
    !
    !END IF

    Pressure=pRiExp*1000    !RS Comment: Unit Conversion
    Enthalpy=hRiExp*1000    !RS Comment: Unit Conversion
    tRiExp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Refrigerant Inlet Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF
    xRiExp=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Expansion Device Refrigerant Inlet Quality
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    Pressure=pRiExp*1000    !RS Comment: Unit Conversion
    Quality=0
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error.'
      ErrorFlag=REFPROPERROR
      CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(24)
      CALL Condenser_Helper_1
      RETURN
    END IF

    !IF (Wabsolute .GT. 0) THEN
    !  Wlocal=LocalOilMassFraction(Wabsolute,xRiExp)
    !  tSat=OilMixtureTsat(RefName,Wlocal,Psat/1000)
    !END IF

    IF (xRiExp .LE. 0) THEN
      tSCiExp=tSat-tRiExp
    ELSE
      tSCiExp=0
    END IF

    CondOUT%COutpRiC=pRiCoil  !RS: Debugging: Formerly OUT(1)
    CondOUT%COuthRiC=hRiCoil  !RS: Debugging: Used only for output    !RS: Debugging: Formerly OUT(1)
    CondOUT%COuttAoC=tAoCoil !RS: Debugging: Only used to be output  !RS: Debugging: Formerly OUT(3)
    CondOUT%COutrhAoC=rhAoCoil    !RS: Debugging: Only used to be output  !RS: Debugging: Formerly OUT(4)
    CondOUT%COutpRoC=pRoCoil  !RS: Debugging: Used only for output    !RS: Debugging: Formerly OUT(5)
    CondOUT%COuthRoC=hRoCoil  !RS: Debugging: Formerly OUT(6)
    CondOUT%COuttRoC=tRoCoil  !RS: Debugging: Used only for output    !RS: Debugging: Formerly OUT(7)
    CondOUT%COutWtAl=WeightAluminum  !RS: Debugging: Only used to be output  !RS: Debugging: Formerly OUT(8)
    CondOUT%COutWtCu=WeightCopper    !RS: Debugging: Only used to be output  !RS: Debugging: Formerly OUT(9)
    CondOUT%COutpRiE=pRiExp  !RS: Debugging: Formerly OUT(10)
    CondOUT%COuthRiE=hRiExp  !RS: Debugging: Formerly OUT(11)
    CondOUT%COuttRiE=tRiExp  !RS: Debugging: Formerly OUT(12)
    CondOUT%COutxRiE=xRiExp  !RS: Debugging: Formerly OUT(13)
    CondOUT%COuttSCiE=tSCiExp !RS: Debugging: Formerly OUT(14)
    CondOUT%COutQC=Qcoil   !RS: Debugging: Formerly OUT(15)
    CondOUT%COutMDisLn=MassDisLn   !RS: Debugging: Formerly OUT(16)
    CondOUT%COutMLiqLn=MassLiqLn   !RS: Debugging: Formerly OUT(17)
    !OUT(18)=0  !RS: Debugging: Set elsewhere
    CondOUT%COutDPAir=DPair   !RS: Debugging: Only used to be output  !RS: Debugging: Formerly OUT(19)

    CondOUT%COutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(20)

    CALL Condenser_Helper_1

    RETURN

  END SUBROUTINE Condenser


  SUBROUTINE Condenser_Helper_1

    IF(CoolHeatModeFlag == 1) THEN
      CoilParams(2)%CoilInletRefTemp=tRiCoil
      CoilParams(2)%CoilOutletRefTemp=tRoCoil
      CoilParams(2)%CoilInletAirTemp=tAiCoil
      CoilParams(2)%CoilOutletAirTemp=tAoCoil
      CoilParams(2)%DPAir=DPair
      CoilParams(2)%DPRef=pRicoil-pRocoil
    ELSE IF(CoolHeatModeFlag == 0) THEN
      CoilParams(1)%CoilInletRefTemp=tRiCoil
      CoilParams(1)%CoilOutletRefTemp=tRoCoil
      CoilParams(1)%CoilInletAirTemp=tAiCoil
      CoilParams(1)%CoilOutletAirTemp=tAoCoil
      CoilParams(1)%DPAir=DPair
      CoilParams(1)%DPRef=pRicoil-pRocoil
    END IF

  END SUBROUTINE Condenser_Helper_1

  !************************************************************************

  SUBROUTINE CalcCondenserInventory(MassCoil,MassLiqCoil,MassVapCoil, &
    LiqTubeLength,VapTubeLength,TwoPhaseTubeLength,NumLiqTubes)

    !------------------------------------------------------------------------
    !Purpose:
    !To get calculate refrigerant inventory in condenser
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod

    IMPLICIT NONE

    REAL,INTENT(OUT) :: MassCoil           !Total refrigerant inventory in coil, kg
    REAL,INTENT(OUT) :: MassLiqCoil        !Total liquid refrigerant inventory in coil, kg
    REAL,INTENT(OUT) :: MassVapCoil        !Total vapor refrigerant inventory in coil, kg
    REAL,INTENT(OUT) :: NumLiqTubes        !Number of Liquid tubes
    REAL,INTENT(OUT) :: LiqTubeLength      !Liquid tube length, m
    REAL,INTENT(OUT) :: VapTubeLength      !Vapor tube length, m
    REAL,INTENT(OUT) :: TwoPhaseTubeLength !Two-phase tube length, m

    INTEGER :: CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel evaporator

    INTEGER I,J,K,II,III,IV !Loop Counter
    REAL Lregion !Region length, m

    LiqTubeLength=0.0
    VapTubeLength=0.0
    TwoPhaseTubeLength=0.0

    MassCoil=0
    MassLiqCoil=0
    MassVapCoil=0

    IF (NumOfChannels .GT. 1) THEN
      CoilType = MCCONDENSER !Microchannel coil
    ELSE
      CoilType = CONDENSERCOIL !Fin-tube coil
    END IF

    IF (CoilType .NE. MCCONDENSER) THEN

      DO I=1, NumOfCkts
        Ckt(I)%Qckt=0.0
        Ckt(I)%tSat=0.0

        !Find first and last simulation tubes
        IF (IsSimpleCoil .EQ. 1) THEN
          FirstTube=1
          LastTube=1
        ELSE
          FirstTube=1
          LastTube=Ckt(I)%Ntube
          IF (Ckt(I)%InSplit .GT. 1) THEN
            FirstTube=2
          END IF

          IF (Ckt(I)%OutJoin .GT. 1) THEN
            LastTube=Ckt(I)%Ntube-1
          END IF
        END IF

        DO J=1,LastTube !Ckt(I)%Ntube !ISI - 10/30/06
          DO K=1,NumOfMods

            !Defining the module properties
            pRiMod=Ckt(I)%Tube(J)%Seg(K)%pRi
            hRiMod=Ckt(I)%Tube(J)%Seg(K)%hRi
            pRoMod=Ckt(I)%Tube(J)%Seg(K)%pRo
            hRoMod=Ckt(I)%Tube(J)%Seg(K)%hRo

            LmodTube=Ckt(I)%Tube(J)%Seg(K)%Len
            AiMod=AiCoil*LmodTube/Lcoil

            Pressure=pRiMod*1000    !RS Comment: Unit Conversion
            Enthalpy=hRiMod*1000    !RS Comment: Unit Conversion
            tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Temperature
            xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Quality
            vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
            vRiMod=1/vRiMod !Module Refrigerant Inlet Specific Volume
            cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Specific Heat
            cpRiMod=cpRiMod/1000    !RS Comment: Unit Conversion
            muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Dynamic Viscosity
            kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Thermal Conductivity
            kRiMod=kRiMod/1000  !RS Comment: Unit Conversion

            Quality=1
            vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vgRiMod=1/vgRiMod   !Module Refrigerant Inlet Vapor Specific Volume

            Quality=0
            vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vfRiMod=1/vfRiMod   !Module Refrigerant Inlet Liquid Specific Volume

            IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
              cpRiMod=0
              muRiMod=0
              kRiMod=0
            END IF

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion
            tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Temperature
            xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Quality
            vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
            vRoMod=1/vRoMod !Module Refrigerant Outlet Specific Volume
            cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Specific Heat
            cpRoMod=cpRoMod/1000    !RS Comment: Unit Conversion
            muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
            kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Thermal Conductivity
            kRoMod=kRoMod/1000  !RS Comment: Unit Conversion

            Quality=1
            vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vgRoMod=1/vgRoMod   !Module Refrigerant Outlet Vapor Specific Volume

            Quality=0
            vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vfRoMod=1/vfRoMod   !Module Refrigerant Outlet Liquid Specific Volume

            Temperature=tRoMod
            Quality=1
            IF (tRoMod+273.15 .GT. Tcr) THEN
              Psat=pRoMod
            ELSE
              Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !Saturation Pressure
              Psat=Psat/1000  !RS Comment: Unit Conversion
            END IF

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Quality=0
            hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Enthalpy
            hfRoMod=hfRoMod/1000    !RS Comment: Unit Conversion
            CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Liquid Specific Heat
            CpfRoMod=CpfRoMod/1000  !RS Comment: Unit Conversion
            mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Dynamic Viscosity
            kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Thermal Conductivity
            kfRoMod=kfRoMod/1000    !RS Comment: Unit Conversion

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Quality=1
            hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Enthalpy
            hgRoMod=hgRoMod/1000    !RS Comment: Unit Conversion
            CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Vapor Specific Heat
            CpgRoMod=CpgRoMod/1000  !RS Comment: Unit Conversion
            mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Vapor Dynamic Viscosity
            kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Thermal Conductivity
            kgRoMod=kgRoMod/1000    !RS Comment: Unit Conversion

            IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
              cpRoMod=0
              muRoMod=0
              kRoMod=0
            END IF

            !Average values
            mu=(muRiMod+muRoMod)/2
            kRef=(kRiMod+kRoMod)/2
            cpRef=(cpRiMod+cpRoMod)/2
            rhoRef=(1/vRiMod+1/vRoMod)/2

            mRefMod=Ckt(I)%mRef

            tAiMod=Ckt(I)%Tube(J)%Seg(K)%tAi
            rhAiMod=Ckt(I)%Tube(J)%Seg(K)%rhAi

            tAoMod=Ckt(I)%Tube(J)%Seg(K)%tAo
            rhAoMod=Ckt(I)%Tube(J)%Seg(K)%rhAo

            Qmod=Ckt(I)%Tube(J)%Seg(K)%Qmod
            hciMod=Ckt(I)%Tube(J)%Seg(K)%hci
            hcoMod=Ckt(I)%Tube(J)%Seg(K)%hco
            ReVap=Ckt(I)%Tube(J)%Seg(K)%ReVap
            ReLiq=Ckt(I)%Tube(J)%Seg(K)%ReLiq
            Cair=Ckt(I)%Tube(J)%Seg(K)%cAir
            Rair=Ckt(I)%Tube(J)%Seg(K)%Rair
            Rtube=Ckt(I)%Tube(J)%Seg(K)%Rtube

            IF (K .EQ. NumOfMods .OR. (J .EQ. LastTube .AND. (Ckt(I)%OutSplit .GT. 1 .OR. Ckt(I)%OutJoin .GT. 1))) THEN
              !Include return bend length
              CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, &
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube+Lreturnbend,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
              !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
              !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
              !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

            ELSE
              CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
              !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
              !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
              !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

            END IF

            Ckt(I)%Tube(J)%Seg(K)%Mass=MassMod

            IF (Ckt(I)%InSplit .GT. 1 .AND. J .EQ. 1) THEN
              LmodTube=0 !RS: Debugging: Keeping skipped tubes from having large default/unassigned values
            ENDIF

            IF ((K .EQ. NumOfMods .AND. J .NE. LastTube) .OR. (J .EQ. LastTube .AND. (Ckt(I)%OutSplit .GT. 1 .OR. Ckt(I)%OutJoin .GT. 1))) THEN !ISI - 02/05/07
              Lregion=LmodTube+Lreturnbend
            ELSE
              Lregion=LmodTube
            END IF

            IF (IsCoolingMode .GT. 0) THEN !Condenser
              IF (xRoMod .GE. 1) THEN
                VapTubeLength=VapTubeLength+Lregion !Superheated region
              ELSEIF (xRiMod .LE. 0) THEN
                LiqTubeLength=LiqTubeLength+Lregion !Subcooled region
              ELSEIF (xRiMod .LT. 1 .AND. xRoMod .GT. 0) THEN
                TwoPhaseTubeLength=TwoPhaseTubeLength+Lregion !Two-phase region
              ELSEIF (xRiMod .GT. 0 .AND. xRoMod .LT. 0) THEN !Condenser outlet
                LiqTubeLength=LiqTubeLength+(hfRoMod-hRoMod)/(hRiMod-hRoMod)*Lregion
                TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hfRoMod-hRoMod)/(hRiMod-hRoMod))*Lregion
              ELSEIF (xRiMod .GT. 1 .AND. xRoMod .LT. 1) THEN !Condenser inlet
                VapTubeLength=VapTubeLength+(hRiMod-hgRoMod)/(hRiMod-hRoMod)*Lregion
                TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hRiMod-hgRoMod)/(hRiMod-hRoMod))*Lregion
              END IF
            ELSE !Evaporator
              IF (xRiMod .GE. 1) THEN
                VapTubeLength=VapTubeLength+Lregion !Superheated region
              ELSEIF (xRoMod .LE. 0) THEN
                LiqTubeLength=LiqTubeLength+Lregion !Subcooled region
              ELSEIF (xRiMod .GT. 0 .AND. xRoMod .LT. 1) THEN
                TwoPhaseTubeLength=TwoPhaseTubeLength+Lregion !Two-phase region
              ELSEIF (xRiMod .LT. 1 .AND. xRoMod .GT. 1) THEN !Evaporator outlet
                VapTubeLength=VapTubeLength+(hRoMod-hgRoMod)/(hRoMod-hRiMod)*Lregion
                TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hRoMod-hgRoMod)/(hRoMod-hRiMod))*Lregion
              ELSEIF (xRiMod .LT. 0 .AND. xRoMod .GT. 0) THEN !Evaporator inlet
                LiqTubeLength=LiqTubeLength+(hfRoMod-hRiMod)/(hRoMod-hRiMod)*Lregion
                TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hfRoMod-hRiMod)/(hRoMod-hRiMod))*Lregion
              END IF
            END IF

            !Total mass inventory
            IF (J .GE. FirstTube .AND. J .LE. LastTube) THEN
              MassCoil=MassCoil+MassMod
              MassLiqCoil=MassLiqCoil+MassLiqMod
              MassVapCoil=MassVapCoil+MassVapMod
            END IF

          END DO !end Nmod

        END DO !end Ntube

        IF (EqCircuits .EQ. 1 .AND. IsUniformVelProfile .OR. IsSimpleCoil .EQ. 1) THEN
          MassCoil=MassCoil*NumOfCkts
          MassLiqCoil=MassLiqCoil*NumOfCkts
          MassVapCoil=MassVapCoil*NumOfCkts
          LiqTubeLength=LiqTubeLength*NumOfCkts
          TwoPhaseTubeLength=TwoPhaseTubeLength*NumOfCkts !ISI - 02/05/07
          VapTubeLength=VapTubeLength*NumOfCkts !ISI - 02/05/07
          EXIT
        END IF

      END DO !end circuit

      NumLiqTubes=LiqTubeLength/(LiqTubeLength+TwoPhaseTubeLength+VapTubeLength) !ISI - 02/05/07

    ELSE !Microchannel coil

      DO I=1, Nl

        DO II=1,Slab(I)%Npass

          DO III=1,1 !NumOfTubes

            DO IV=1, NumOfMods

              pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
              hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi
              pRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo
              hRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo

              LmodTube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len
              AiMod=AiCoil*LmodTube/Lcoil

              Pressure=pRiMod*1000    !RS Comment: Unit Conversion
              Enthalpy=hRiMod*1000    !RS Comment: Unit Conversion
              tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Temperature
              xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Quality
              vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
              vRiMod=1/vRiMod !Module Refrigerant Inlet Specific Volume
              cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Specific Heat
              cpRiMod=cpRiMod/1000    !RS Comment: Unit Conversion
              muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Dynamic Viscosity
              kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Thermal Conductivity
              kRiMod=kRiMod/1000  !RS Comment: Unit Conversion

              Quality=1
              vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vgRiMod=1/vgRiMod   !Module Refrigerant Inlet Vapor Specific Volume

              Quality=0
              vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vfRiMod=1/vfRiMod   !Module Refrigerant Inlet Liquid Specific Volume

              IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
                cpRiMod=0
                muRiMod=0
                kRiMod=0
              END IF

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion
              tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Temperature
              xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Quality
              vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
              vRoMod=1/vRoMod !Module Refrigerant Outlet Specific Volume
              cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Specific Heat
              cpRoMod=cpRoMod/1000    !RS Comment: Unit Conversion
              muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
              kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Thermal Conductivity
              kRoMod=kRoMod/1000  !RS Comment: Unit Conversion

              Quality=1
              vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vgRoMod=1/vgRoMod   !Module Refrigerant Outlet Vapor Specific Volume

              Quality=0
              vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vfRoMod=1/vfRoMod   !Module Refrigerant Outlet Liquid Specific Volume

              Temperature=tRoMod
              Quality=1
              IF (tRoMod+273.15 .GT. Tcr) THEN
                Psat=pRoMod
              ELSE
                Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !Saturation Pressure
                Psat=Psat/1000  !RS Comment: Unit Conversion
              END IF

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Quality=0
              hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Enthalpy
              hfRoMod=hfRoMod/1000    !RS Comment: Unit Conversion
              CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Liquid Specific Heat
              CpfRoMod=CpfRoMod/1000  !RS Comment: Unit Conversion
              mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Dynamic Viscosity
              kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Thermal Conductivity
              kfRoMod=kfRoMod/1000    !RS Comment: Unit Conversion

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Quality=1
              hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Enthalpy
              hgRoMod=hgRoMod/1000    !RS Comment: Unit Conversion
              CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Vapor Specific Heat
              CpgRoMod=CpgRoMod/1000  !RS Comment: Unit Conversion
              mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
              kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Thermal Conductivity
              kgRoMod=kgRoMod/1000    !RS Comment: Unit Conversion

              IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
                cpRoMod=0
                muRoMod=0
                kRoMod=0
              END IF

              !Average values
              mu=(muRiMod+muRoMod)/2
              kRef=(kRiMod+kRoMod)/2
              cpRef=(cpRiMod+cpRoMod)/2
              rhoRef=(1/vRiMod+1/vRoMod)/2

              !Defining more of the module values
              mRefMod=mRefTot/Slab(I)%Pass(II)%Ntube !Ckt(I)%mRef

              tAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi
              rhAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi

              tAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo
              rhAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAo

              Qmod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod
              hciMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hci
              hcoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco
              ReVap=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReVap
              ReLiq=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReLiq
              Cair=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir
              Rair=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair
              Rtube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube

              CALL Inventory(CoilType,Dchannel,mRefMod/NumOfChannels,hgRoMod,hfRoMod, & !hRiMod,hRoMod, & !RS: Debugging: Extraneous
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
              !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
              !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
              !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

              MassMod=MassMod*NumOfChannels
              Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Mass=MassMod

              Lregion=LmodTube

              IF (IsCoolingMode .GT. 0) THEN !Condenser
                IF (xRoMod .GE. 1) THEN
                  VapTubeLength=VapTubeLength+Lregion !Superheated region
                ELSEIF (xRiMod .LE. 0) THEN
                  LiqTubeLength=LiqTubeLength+Lregion !Subcooled region
                ELSEIF (xRiMod .LT. 1 .AND. xRoMod .GT. 0) THEN
                  TwoPhaseTubeLength=TwoPhaseTubeLength+Lregion !Two-phase region
                ELSEIF (xRiMod .GT. 0 .AND. xRoMod .LT. 0) THEN !Condenser outlet
                  LiqTubeLength=LiqTubeLength+(hfRoMod-hRoMod)/(hRiMod-hRoMod)*Lregion
                  TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hfRoMod-hRoMod)/(hRiMod-hRoMod))*Lregion
                ELSEIF (xRiMod .GT. 1 .AND. xRoMod .LT. 1) THEN !Condenser inlet
                  VapTubeLength=VapTubeLength+(hRiMod-hgRoMod)/(hRiMod-hRoMod)*Lregion
                  TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hRiMod-hgRoMod)/(hRiMod-hRoMod))*Lregion
                END IF
              ELSE !Evaporator
                IF (xRiMod .GE. 1) THEN
                  VapTubeLength=VapTubeLength+Lregion !Superheated region
                ELSEIF (xRoMod .LE. 0) THEN
                  LiqTubeLength=LiqTubeLength+Lregion !Subcooled region
                ELSEIF (xRiMod .GT. 0 .AND. xRoMod .LT. 1) THEN
                  TwoPhaseTubeLength=TwoPhaseTubeLength+Lregion !Two-phase region
                ELSEIF (xRiMod .LT. 1 .AND. xRoMod .GT. 1) THEN !Evaporator outlet
                  VapTubeLength=VapTubeLength+(hRoMod-hgRoMod)/(hRoMod-hRiMod)*Lregion
                  TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hRoMod-hgRoMod)/(hRoMod-hRiMod))*Lregion
                ELSEIF (xRiMod .LT. 0 .AND. xRoMod .GT. 0) THEN !Evaporator inlet
                  LiqTubeLength=LiqTubeLength+(hfRoMod-hRiMod)/(hRoMod-hRiMod)*Lregion
                  TwoPhaseTubeLength=TwoPhaseTubeLength+(1-(hfRoMod-hRiMod)/(hRoMod-hRiMod))*Lregion
                END IF
              END IF

              !Total mass inventory

              MassCoil=MassCoil+MassMod*Slab(I)%Pass(II)%Ntube
              MassLiqCoil=MassLiqCoil+MassLiqMod*Slab(I)%Pass(II)%Ntube
              MassVapCoil=MassVapCoil+MassVapMod*Slab(I)%Pass(II)%Ntube

            END DO !end Nmod

          END DO !end Ntube

        END DO !end pass

      END DO !end slab

      NumLiqTubes=LiqTubeLength/(LiqTubeLength+TwoPhaseTubeLength+VapTubeLength) !ISI - 02/05/07

    END IF

    RETURN

  END SUBROUTINE CalcCondenserInventory

  !************************************************************************

  SUBROUTINE PrintCondenserResult

    !------------------------------------------------------------------------
    !Purpose:
    !To print simulation result to output file "condenser.csv"
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod

    IMPLICIT NONE

    INTEGER CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel evaporator

    REAL :: MassCoil    !Total refrigerant inventory in coil, kg
    REAL :: MassLiqCoil !Total liquid refrigerant inventory in coil, kg
    REAL :: MassVapCoil !Total vapor refrigerant inventory in coil, kg

    INTEGER I,J,K,II,III,IV !Loop Counter

    CHARACTER(LEN=13),PARAMETER :: FMT_100 = "(50(A12,','))"
    CHARACTER(LEN=25),PARAMETER :: FMT_104 = "(3(I3,','),50(F10.3,','))"

    OPEN (16,FILE='Condenser.csv')
    !OPEN (16,FILE='Condenser_PlainFin.csv')   !RS: Test case output file

    MassCoil=0
    MassLiqCoil=0
    MassVapCoil=0

    IF (NumOfChannels .GT. 1) THEN
      CoilType = MCCONDENSER
    ELSE
      CoilType = CONDENSERCOIL
    END IF

    IF (CoilType .NE. MCCONDENSER) THEN

      WRITE(16,FMT_100)'Nckt','Ntube','Nmod','tRi(C)','tRo(C)','pRi(kPa)','pRo(kPa)', &
      'hRi(kJ/kg)','hRo(kJ/kg)','xRi','xRo','tAi(C)','tAo(C)', &
      'rhAi','rhAo','hci(W/m2K)','EF','hco(W/m2K)', &
      'mu(uPa-s)','k(W/mK)','cp(kJ/kgK)','rho(kg/m3)','ReVap','ReLiq', &
      'Qmod(W)','MassLiq(g)','MassVap(g)','MassTot(g)','mdot(kg/h)' !, &
      !'cpAir', 'mAiMod (kg/s)', 'hAiMod', 'hAoMod' !RS: Debugging: Adding cp, air m_dot and h's for a heat balance

      DO I=1, NumOfCkts
        Ckt(I)%Qckt=0.0
        Ckt(I)%tSat=0.0

        !Find first and last simulation tubes
        FirstTube=1
        LastTube=Ckt(I)%Ntube
        IF (Ckt(I)%InSplit .GT. 1) THEN
          FirstTube=2
        END IF
        IF (Ckt(I)%OutJoin .GT. 1) THEN
          LastTube=Ckt(I)%Ntube-1
        END IF

        DO J=1,Ckt(I)%Ntube
          DO K=1,NumOfMods

            !Defining the module properties
            pRiMod=Ckt(I)%Tube(J)%Seg(K)%pRi
            hRiMod=Ckt(I)%Tube(J)%Seg(K)%hRi
            pRoMod=Ckt(I)%Tube(J)%Seg(K)%pRo
            hRoMod=Ckt(I)%Tube(J)%Seg(K)%hRo

            LmodTube=Ckt(I)%Tube(J)%Seg(K)%Len
            AiMod=AiCoil*LmodTube/Lcoil

            Pressure=pRiMod*1000    !RS Comment: Unit Conversion
            Enthalpy=hRiMod*1000    !RS Comment: Unit Conversion
            tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Temperature
            xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Quality
            vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
            vRiMod=1/vRiMod !Module Refrigerant Inlet Specific Volume
            cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Specific Heat
            cpRiMod=cpRiMod/1000    !RS Comment: Unit Conversion
            muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Dynamic Viscosity
            kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Thermal Conductivity
            kRiMod=kRiMod/1000  !RS Comment: Unit Conversion

            Quality=1
            vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vgRiMod=1/vgRiMod   !Module Refrigerant Inlet Vapor Specific Volume

            Quality=0
            vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vfRiMod=1/vfRiMod   !Module Refrigerant Inlet Fluid Specific Volume

            IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
              cpRiMod=0
              muRiMod=0
              kRiMod=0
            END IF

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion
            tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Temperature
            xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Quality
            vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
            vRoMod=1/vRoMod !Module Refrigerant Outlet Specific Volume
            cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Specific Heat
            cpRoMod=cpRoMod/1000    !RS Comment: Unit Conversion
            muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
            kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Conductivity
            kRoMod=kRoMod/1000  !RS Comment: Unit Conversion

            Quality=1
            vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vgRoMod=1/vgRoMod   !Module Refrigerant Outlet Vapor Specific Heat

            Quality=0
            vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
            vfRoMod=1/vfRoMod   !Module Refrigerant Outlet Liquid Specific Heat

            Temperature=tRoMod
            Quality=1
            IF (tRoMod+273.15 .GT. Tcr) THEN
              Psat=pRoMod
            ELSE
              Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !Saturation Pressure
              Psat=Psat/1000  !RS Comment: Unit Conversion
            END IF

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Quality=0
            hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)      !Module Refrigerant Outlet Liquid Enthalpy
            hfRoMod=hfRoMod/1000    !RS Comment: Unit Conversion
            CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Liquid Specific Heat
            CpfRoMod=CpfRoMod/1000  !RS Comment: Unit Conversion
            mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Dynamic Viscosity
            kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Thermal Conductivity
            kfRoMod=kfRoMod/1000    !RS Comment: Unit Conversion

            Pressure=pRoMod*1000    !RS Comment: Unit Conversion
            Quality=1
            hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Enthalpy
            hgRoMod=hgRoMod/1000    !RS Comment: Unit Conversion
            CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Vapor Specific Heat
            CpgRoMod=CpgRoMod/1000  !RS Comment: Unit Conversion
            mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Vapor Dynamic Viscosity
            kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Thermal Conductivity
            kgRoMod=kgRoMod/1000    !RS Comment: Unit Conversion

            IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
              cpRoMod=0
              muRoMod=0
              kRoMod=0
            END IF

            !Averag values
            mu=(muRiMod+muRoMod)/2
            kRef=(kRiMod+kRoMod)/2
            cpRef=(cpRiMod+cpRoMod)/2
            rhoRef=(1/vRiMod+1/vRoMod)/2

            !Defining more of the module values
            mRefMod=Ckt(I)%mRef

            tAiMod=Ckt(I)%Tube(J)%Seg(K)%tAi
            rhAiMod=Ckt(I)%Tube(J)%Seg(K)%rhAi

            tAoMod=Ckt(I)%Tube(J)%Seg(K)%tAo
            rhAoMod=Ckt(I)%Tube(J)%Seg(K)%rhAo

            Qmod=Ckt(I)%Tube(J)%Seg(K)%Qmod
            hciMod=Ckt(I)%Tube(J)%Seg(K)%hci
            EFref=Ckt(I)%Tube(J)%Seg(K)%EFref
            hcoMod=Ckt(I)%Tube(J)%Seg(K)%hco
            ReVap=Ckt(I)%Tube(J)%Seg(K)%ReVap
            ReLiq=Ckt(I)%Tube(J)%Seg(K)%ReLiq
            Cair=Ckt(I)%Tube(J)%Seg(K)%cAir
            Rair=Ckt(I)%Tube(J)%Seg(K)%Rair
            Rtube=Ckt(I)%Tube(J)%Seg(K)%Rtube

            IF (K .EQ. NumOfMods .OR. (J .EQ. LastTube .AND. (Ckt(I)%OutSplit .GT. 1 .OR. Ckt(I)%OutJoin .GT. 1))) THEN
              !Include return bend length
              CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, & !RS: Debugging: Extraneous
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube+Lreturnbend,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
              !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
              !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
              !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)
            ELSE
              CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, & !RS: Debugging: Extraneous
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
            END IF

            Ckt(I)%Tube(J)%Seg(K)%Mass=MassMod

            !Total mass inventory
            IF (J .GE. FirstTube .AND. J .LE. LastTube) THEN
              MassCoil=MassCoil+MassMod
              MassLiqCoil=MassLiqCoil+MassLiqMod
              MassVapCoil=MassVapCoil+MassVapMod
            END IF

            IF (xRiMod .EQ. -100) THEN
              xRiMod=0.0
            END IF
            IF (xRiMod .EQ. 100) THEN
              xRiMod=1.0
            END IF
            IF (xRoMod .EQ. -100) THEN
              xRoMod=0.0
            END IF
            IF (xRoMod .EQ. 100) THEN
              xRoMod=1.0
            END IF

            MassMod=Ckt(I)%Tube(J)%Seg(K)%Mass
            WRITE(16,FMT_104)I,J,K,tRiMod,tRoMod,pRiMod,pRoMod,hRiMod,hRoMod, &
            xRiMod,xRoMod,tAiMod,tAoMod,rhAiMod,rhAoMod, &
            hciMod*1000,EFref,hcoMod*1000,mu*1e6,kRef*1e3,cpRef,rhoRef,ReVap,ReLiq, &
            Qmod*1000,MassLiqMod*1000,MassVapMod*1000,MassMod*1000, &
            mRefMod*3600 !, CPAir, mAiMod, hAiMod, hAoMod, AirProp(4), TestH !RS: Debugging: Adding the cp, air m_dot and h's for a heat balance

          END DO !end Nmod

        END DO !end Ntube

      END DO !end circuit

    ELSE

      WRITE(16,FMT_100)'Nslab','Npass','Nmod','tRi(C)','tRo(C)','pRi(kPa)','pRo(kPa)', &
      'hRi(kJ/kg)','hRo(kJ/kg)','xRi','xRo','tAi(C)','tAo(C)', &
      'rhAi','rhAo','hci(W/m2K)','hco(W/m2K)', &
      'mu(uPa-s)','k(W/mK)','cp(kJ/kgK)','rho(kg/m3)','ReVap','ReLiq', &
      'Qmod(W)','MassLiq(g)','MassVap(g)','MassTot(g)','mdot(kg/h)'

      DO I=1, Nl

        DO II=1,Slab(I)%Npass

          DO III=1,1 !Num Of tubes

            DO IV=1, NumOfMods

              !Defininf module properties
              pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
              hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi
              pRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo
              hRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo

              LmodTube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len
              AiMod=AiCoil*LmodTube/Lcoil

              Pressure=pRiMod*1000    !RS Comment: Unit Conversion
              Enthalpy=hRiMod*1000    !RS Comment: Unit Conversion
              tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Temperature
              xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Quality
              vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
              vRiMod=1/vRiMod !Module Refrigerant Inlet Specific Volume
              cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Specific Heat
              cpRiMod=cpRiMod/1000    !RS Comment: Unit Conversion
              muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Dynamic Viscosity
              kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Thermal Conductivity
              kRiMod=kRiMod/1000  !RS Comment: Unit Conversion

              Quality=1
              vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vgRiMod=1/vgRiMod   !Module Refrigerant Inlet Vapor Specific Volume

              Quality=0
              vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vfRiMod=1/vfRiMod   !Module Refrigerant Inlet Liquid Specific Volume

              IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
                cpRiMod=0
                muRiMod=0
                kRiMod=0
              END IF

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion
              tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Temperature
              xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Quality
              vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
              vRoMod=1/vRoMod !Module Refrigerant Outlet Specific Volume
              cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Specific Heat
              cpRoMod=cpRoMod/1000    !RS Comment: Unit Conversion
              muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
              kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Thermal Conductivity
              kRoMod=kRoMod/1000  !RS Comment: Unit Conversion

              Quality=1
              vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vgRoMod=1/vgRoMod   !Module Refrigerant Outlet Vapor Dynamic Viscosity

              Quality=0
              vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
              vfRoMod=1/vfRoMod   !Module Refrigerant Outlet Liquid Dynamic Viscosity

              Temperature=tRoMod
              Quality=1
              IF (tRoMod+273.15 .GT. Tcr) THEN
                Psat=pRoMod
              ELSE
                Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !Saturation Pressure
                Psat=Psat/1000  !RS Comment: Unit Conversion
              END IF

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Quality=0
              hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Enthalpy
              hfRoMod=hfRoMod/1000    !RS Comment: Unit Conversion
              CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Liquid Specific Heat
              CpfRoMod=CpfRoMod/1000  !RS Comment: Unit Conversion
              mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Dynamic Viscosity
              kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Thermal Conductivity
              kfRoMod=kfRoMod/1000    !RS Comment: Unit Conversion

              Pressure=pRoMod*1000    !RS Comment: Unit Conversion
              Quality=1
              hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Enthalpy
              hgRoMod=hgRoMod/1000    !RS Comment: Unit Conversion
              CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Vapor Specific Heat
              CpgRoMod=CpgRoMod/1000  !RS Comment: Unit Conversion
              mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Dynamic Viscosity
              kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Thermal Conductivity
              kgRoMod=kgRoMod/1000    !RS Comment: Unit Conversion

              IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
                cpRoMod=0
                muRoMod=0
                kRoMod=0
              END IF

              !Average values
              mu=(muRiMod+muRoMod)/2
              kRef=(kRiMod+kRoMod)/2
              cpRef=(cpRiMod+cpRoMod)/2
              rhoRef=(1/vRiMod+1/vRoMod)/2

              !Defining more of the module values
              mRefMod=mRefTot/Slab(I)%Pass(II)%Ntube !Ckt(I)%mRef

              tAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi
              rhAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi

              tAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo
              rhAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAo

              Qmod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod
              hciMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hci
              hcoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco
              ReVap=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReVap
              ReLiq=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReLiq
              Cair=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir
              Rair=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair
              Rtube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube

              CALL Inventory(CoilType,Dchannel,mRefMod/NumOfChannels,hgRoMod,hfRoMod, & !hRiMod,hRoMod, & !RS: Debugging: Extraneous
              xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
              LmodTube,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
              !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
              !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
              !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

              MassMod=MassMod*NumOfChannels
              Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Mass=MassMod

              !Total mass inventory
              IF (J .GE. FirstTube .AND. J .LE. LastTube) THEN
                MassCoil=MassCoil+MassMod*Slab(I)%Pass(II)%Ntube
                MassLiqCoil=MassLiqCoil+MassLiqMod*Slab(I)%Pass(II)%Ntube
                MassVapCoil=MassVapCoil+MassVapMod*Slab(I)%Pass(II)%Ntube
              END IF

              !Keeping the quality as a decimal
              IF (xRiMod .EQ. -100) THEN
                xRiMod=0.0
              END IF
              IF (xRiMod .EQ. 100) THEN
                xRiMod=1.0
              END IF
              IF (xRoMod .EQ. -100) THEN
                xRoMod=0.0
              END IF
              IF (xRoMod .EQ. 100) THEN
                xRoMod=1.0
              END IF

              MassMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Mass
              WRITE(16,FMT_104)I,II,IV,tRiMod,tRoMod,pRiMod,pRoMod,hRiMod,hRoMod, &
              xRiMod,xRoMod,tAiMod,tAoMod,rhAiMod,rhAoMod, &
              hciMod*1000,hcoMod*1000,mu*1e6,kRef*1e3,cpRef,rhoRef,ReVap,ReLiq, &
              Qmod*1000,MassLiqMod*1000,MassVapMod*1000,MassMod*1000, &
              mRefMod*3600

            END DO !end Nmod

          END DO !end Ntube

        END DO !end circuit

      END DO !end Slab

    END IF

    CLOSE(16)
    RETURN

  END SUBROUTINE PrintCondenserResult

  !************************************************************************

  SUBROUTINE InitCondenserCoil(CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To initialize condenser geometry and circuiting
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE UnitConvertMod
    USE InputProcessor

    IMPLICIT NONE

    INTEGER,INTENT(OUT) :: CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel Evaporator

    INTEGER I,J,II,III,IV !Loop counter
    INTEGER NumOfPasses !Number of passes
    INTEGER Ntubes !Number of tubes
    INTEGER NumOfInlets !Number of inlets
    LOGICAL IsSIunit !SI unit input flag
    LOGICAL IsShift !Is shift tube flag (for staggered tubes)

    INTEGER, PARAMETER :: MaxNameLength = 200

    REAL, PARAMETER :: UnitP     = 6.8947453 !(psi X UnitP = kPa)
    REAL, PARAMETER :: UnitM     = 0.4536    !(lbm X UnitM = kg)
    REAL, PARAMETER :: UnitL     = 0.3048    !(ft X UnitL = m)
    REAL, PARAMETER :: UnitK     = 0.1442E-3 !(Btu-in/hr-ft2-F X UnitK = kW/m-C)

    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    CHARACTER(len=MaxNameLength) :: ModelName !Model Name tells how to address Fin-Tube Coil or MicroChannel, etc.

    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    INTEGER :: TempNumofMods  !RS: Debugging: Temporary variable to hold real number of modules

    !RS: Debugging: Bringing over from GetHPSimInputs
    CHARACTER(len=MaxNameLength)ODC_FinName
    CHARACTER(len=MaxNameLength)IDC_FinName
    REAL :: ODC_TubeID
    REAL :: IDC_TubeID

    INTEGER,PARAMETER :: SI=1
    INTEGER,PARAMETER :: IP=2

    !FLOW:

    NumOfSections=1 !ISI - 09/10/07

    If (errorflag .EQ. 1) Then
      errorflag = 0 !JG set a zero error
    End IF

    !***** Get circuit info *****
    IF (ErrorFlag .NE. NOERROR) THEN
      ErrorFlag=CKTFILEERROR !RS: Debugging: Could also be (and often is) a convergence error
      CALL InitCondenserCoil_Helper_1
      RETURN
    END IF

    !**************************** Model *************************************

    CALL GetObjectItem('ODCcktModel',1,Alphas,NumAlphas, &
    TmpNumbers,NumNumbers,Status)
    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    ModelName = Alphas(1)

    IF (ModelName .EQ. 'MICROCHANNEL COIL') THEN
      CoilType = MCCONDENSER
    ELSE
      CoilType = CONDENSERCOIL
    END IF

    IF (CoilType .EQ. CONDENSERCOIL) THEN !Fin-tube coil

      IF (IsCoolingMode .GT. 0) THEN  !ODC or IDC circuits?
        !ODC circuit here

        !**************************** Geometry *************************************

        !Reading in the variable values
        CALL GetObjectItem('ODCcktGeometry',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        SELECT CASE (Alphas(1)(1:1))
        CASE ('F','f')
          IsSIunit=.FALSE.
        CASE DEFAULT
          IsSIunit=.TRUE.
        END SELECT

        !Defining the variables
        FinType = Numbers(1)
        FinPitch = Numbers(2)
        Kfin = Numbers(3)
        FinThk = Numbers(4)
        TubeType = Numbers(5)
        ODtube = Numbers(6)
        IDtube = Numbers(7)
        Ktube = Numbers(8)
        Pt = Numbers(9)
        Pl = Numbers(10)
        Nl = Numbers(11)
        Nt = Numbers(12)
        Ltube = Numbers(13)

        IF (Ltube .LE. 1e-3) THEN
          ErrorFlag=ZEROLENCOILERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        NumOfMods = Numbers(14)
        TempNumofMods=NumOfMods
        NumOfCkts = Numbers(15)

        SELECT CASE (Alphas(5)(1:1))    !Tube Shift Flag
        CASE ('F','f')
          IsShift=.FALSE.
        CASE DEFAULT
          IsShift=.TRUE.
        END SELECT

        !************************** Inputs ****************************************

        !***************** Outdoor coil data *****************    !RS: Debugging: Moving: Evaporator & Condenser

        CALL GetObjectItem('OutdoorCoilData',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        !Fin type (1-smooth; 2-Wavy; 3-louvered)

        CondPAR%CondFinType = Numbers(1)  !RS: Debugging: Formerly CondPAR(29)

        ODC_FinName = Alphas(1)

        CondPAR%CondFinPitch = Numbers(2)  !RS: Debugging: Formerly CondPAR(22)
        CondPAR%CondFinThermCon = Numbers(3) !Conductivity of Fin !RS: Debugging: Formerly CondPAR(23)
        CondPAR%CondFinThick = Numbers(4)   !Fin Thickness !RS: Debugging: Formerly CondPAR(21)
        CondPAR%CondTube = Numbers(5) !Numerical Denotion of Tube Type    !RS: Debugging: Formerly CondPAR(37)
        ODC_TubeID = Numbers(6)   !Tube Inner Diameter
        CondPAR%CondCoilTOD = Numbers(7)   !Tube Outer Diameter   !RS: Debugging: Formerly CondPAR(15)
        CondPAR%CondCoilTThermCon = Numbers(8)    !Tube Conductivity    !RS: Debugging: Formerly CondPAR(18)
        CondPAR%CondRspc = Numbers(9)   !Tube Lateral Spacing  !RS: Debugging: Formerly CondPAR(20)
        CondPAR%CondTspc = Numbers(10)   !Tube Vertical Spacing    !RS: Debugging: Formerly CondPAR(19)
        CondPAR%CondNl= Numbers(11)  !Number of Rows !RS: Debugging: Formerly CondPAR(25)
        CondPAR%CondNt = Numbers(12)  !Number of Tubes per Row   !RS: Debugging: Formerly CondPAR(24)
        CondPAR%CondNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly CondPAR(26)
        CondPAR%CondNumMod = Numbers(14)    !Number of Segments  !RS: Debugging: Formerly CondPAR(28)
        CondPAR%CondCoilSTLen = Numbers(15)   !Single Tube Length   !RS: Debugging: Formerly CondPAR(17)
        CondPAR%CondMultRefQT = Numbers(16)   !Ref Side Heat Transfer Multiplier    !RS: Debugging: Formerly CondPAR(30)
        CondPAR%CondMultRefPD = Numbers(17) !Ref Side Pressure Drop Multiplier  !RS: Debugging: Formerly CondPAR(31)
        CondPAR%CondMultAirQT = Numbers(18)   !Air Side Heat Transfer Multiplier    !RS: Debugging: Formerly CondPAR(32)
        CondPAR%CondMultAirPD = Numbers(19) !Air Side Pressure Drop Multiplier  !RS: Debugging: Formerly CondPAR(33)

        !  IF (UNIT .EQ. SI) THEN  !RS: Debugging: Unit conversions !RS: Debugging: The following code doesn't seem to make a difference
        ! CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/1000   !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
        !    CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick/1000   !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
        !    CondPAR%CondTspc=CondPAR%CondTspc/1000   !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
        !    CondPAR%CondRspc=CondPAR%CondRspc/1000   !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
        !    CondPAR%CondFinThick=CondPAR%CondFinThick/1000   !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
        ! CondPAR%CondFanPwr=CondPAR%CondFanPwr/1000   !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
        !ELSE
        ! CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/12*UnitL       !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
        !    CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*0.001/12*UnitL !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
        !    CondPAR%CondCoilSTLen=CondPAR%CondCoilSTLen/12*UnitL       !Tube length, m  !RS: Debugging: Formerly CondPAR(17)
        !    CondPAR%CondCoilTThermCon=CondPAR%CondCoilTThermCon*UnitK          !Tube thermal conductivity, kW/m-C   !RS: Debugging: Formerly CondPAR(18)
        !    CondPAR%CondTspc=CondPAR%CondTspc/12*UnitL       !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
        !    CondPAR%CondRspc=CondPAR%CondRspc/12*UnitL       !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
        !    CondPAR%CondFinThick=CondPAR%CondFinThick*0.001/12*UnitL !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
        !    CondPAR%CondFinPitch=CondPAR%CondFinPitch*12/UnitL       !Fin pitch, fin/m    !RS: Debugging: Formerly CondPAR(22)
        !    CondPAR%CondFinThermCon=CondPAR%CondFinThermCon*UnitK          !Fin thermal conductivity, kW/m-C    !RS: Debugging: Formerly CondPAR(23)
        ! CondPAR%CondFanPwr=CondPAR%CondFanPwr*1E-3           !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
        ! !CondPAR%CondBarPress=CondPAR%CondBarPress*UnitP          !Barometric pressure, kPa    !RS: Debugging: Formerly CondPAR(38)
        !END IF
        !
        !Tube wall thickness, mm or mil
        CondPAR%CondCoilTWThick=(CondPAR%CondCoilTOD-ODC_TubeID)/2    !RS: Debugging: Formerly CondPAR(16)
        IF (Unit .EQ. IP) THEN
          CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*1000
        END IF

        CondPAR%CondCoolMode=IsCoolingMode   !RS: Debugging: Formerly CondPAR(27)
        !CondPAR(36)=ODC_SurfAbsorptivity   !RS: Debugging: Extraneous

        !***************** Outdoor fan data ***************** !RS: Debugging: Moving: Evaporator & Condenser

        CALL GetObjectItem('OutdoorFanData',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        CondPAR%CondFanPwr = Numbers(1) !Fan Power   !RS: Debugging: Formerly CondPAR(34)
        !VdotODfan = Numbers(2)    !Fan Air Flow Rate
        CondPAR%CondFanLoc = Numbers(3)   !Draw Through (1) or Blow Through (2)  !RS: Debugging: Formerly CondPAR(35)


        !*************************** Circuiting ************************************

        CALL FinTubeCoilUnitConvert(IsSIUnit,FinPitch,Kfin,FinThk, &
        ODtube,IDtube,Ktube,Pt,Pl,Ltube)

        TubeThk=(ODtube-IDtube)/2

        IF (IsShift) THEN
          ShiftTube = 1
        ELSE
          ShiftTube = 0
        END IF

        IF (Pl .EQ. 0) THEN
          Pl = Pt
        END IF

        !Branch#,#Tubes
        IF (ErrorFlag .NE. NOERROR) THEN
          ErrorFlag=CKTFILEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        NumOfTubes=Nl*Nt

        !Fin spacing
        FinSpg = 1/FinPitch-FinThk

        !For plate finned tube
        FinHeight=0 !No use for plate finned tube
        TubeDepth=ODtube
        TubeHeight=ODtube
        NumOfChannels=1
        Dchannel=IDtube

        IF (FinSpg .LT. FinThk) THEN
          ErrorFlag=COILFINERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        IF (Pt .LT. ODtube+2*FinThk) THEN
          ErrorFlag=COILTUBEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        IF (IsSimpleCoil .EQ. 1) THEN   !This is an open block currently; it will need fixing. !RS Comment: No longer an open block!
          !NumOfMods=3
          IF (.NOT. ALLOCATED(Ckt)) THEN
            ALLOCATE(Ckt(NumOfCkts))
            ALLOCATE(DisLnSeg(NumOfMods))
            ALLOCATE(LiqLnSeg(NumOfMods))
            DO I=1, NumOfCkts
              Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
              ALLOCATE(Ckt(I)%Tube(1))
              ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
            END DO
          END IF
        END IF

        IF (.NOT. ALLOCATED(Ckt)) THEN
          !RS: Debugging: The following comes from the evaporator module
          CALL InitCondenserStructures(TempNumofMods)
          !IF (TempNumofMods .NE. 3) THEN
          !   NumOfMods=TempNumofMods
          !   DEALLOCATE(DisLnSeg)
          !   ALLOCATE(DisLnSeg(NumOfMods))
          !   DEALLOCATE(LiqLnSeg)
          !   ALLOCATE(LiqLnSeg(NumOfMods))
          !   DO I=1, NumOfCkts
          !       DO J=1,Ckt(I)%Ntube
          !           ALLOCATE(Ckt(I)%Tube(J)%Seg(NumOfMods))
          !       END DO
          !   END DO
          !END IF
        END IF

        IF (IsSimpleCoil .NE. 1) THEN

          !IF (IsSimpleCoil .EQ. 1) THEN   !This is an open block currently; it will need fixing. !RS Comment: No longer an open block!
          !        NumOfMods=3
          !        ALLOCATE(Ckt(NumOfCkts))
          !        DO I=1, NumOfCkts
          !            Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
          !            ALLOCATE(Ckt(I)%Tube(1))
          !            ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
          !        END DO
          !ELSE
          !        !ISI - 09/10/07
          !        ALLOCATE(CoilSection(NumOfSections))
          !        ALLOCATE(Tube(NumOfTubes))
          !        ALLOCATE(Tube2D(Nl,Nt))
          !        ALLOCATE(JoinTubes(NumOfTubes))
          !
          !        DO I=1, NumOfTubes
          !            ALLOCATE(Tube(I)%Seg(NumOfMods))
          !        END DO
          !
          !        DO I=1, Nl
          !            DO J=1, Nt
          !                ALLOCATE(Tube2D(I,J)%Seg(NumOfMods))
          !            END DO
          !        END DO
          !
          !        ALLOCATE(Ckt(NumOfCkts))
          !        ALLOCATE(mRefIter(NumOfCkts))
          !        CALL GetObjectItem('ODCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
          !                            TmpNumbers,NumNumbers,Status)
          !        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
          !
          !        DO I=1, NumOfCkts
          !            Ckt(I)%Ntube=Numbers(I)
          !            IF (ErrorFlag .NE. NOERROR) THEN
          !                ErrorFlag=CKTFILEERROR
          !                CALL InitCondenserCoil_Helper_1
          !                RETURN
          !            END IF
          !            ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
          !            ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
          !        END DO

          !Check if all circuit have the same number of tubes !ISI - 09/12/06
          IsSameNumOfTubes=.TRUE.
          Ntubes=Ckt(1)%Ntube
          DO I=2, NumOfCkts
            IF (Ntubes .NE. Ckt(I)%Ntube) THEN
              IsSameNumOfTubes=.FALSE.
              EXIT
            END IF
          END DO

          IF (ErrorFlag .NE. NOERROR) THEN
            ErrorFlag=CKTFILEERROR
            CALL InitCondenserCoil_Helper_1
            RETURN
          END IF

          DO I=1, NumOfCkts
            IF (I .EQ. 1) THEN
              CALL GetObjectItem('ODCcktCircuit1_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 2) THEN
              CALL GetObjectItem('ODCcktCircuit2_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 3) THEN
              CALL GetObjectItem('ODCcktCircuit3_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 4) THEN
              CALL GetObjectItem('ODCcktCircuit4_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 5) THEN
              CALL GetObjectItem('ODCcktCircuit5_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 6) THEN
              CALL GetObjectItem('ODCcktCircuit6_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 7) THEN
              CALL GetObjectItem('ODCcktCircuit7_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 8) THEN
              CALL GetObjectItem('ODCcktCircuit8_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSE
              CALL GetObjectItem('ODCcktCircuit9_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            END IF
            Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            DO J=1, Ckt(I)%Ntube
              Ckt(I)%TubeSequence(J)=Numbers(J)   !RS Comment: Populating the tube sequence arrays
            END DO
            IF (ErrorFlag .NE. NOERROR) THEN
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
          END DO

          !************************* Velocity Profile ********************************

          CoilSection(NumOfSections)%NumOfCkts=NumOfCkts

          DO I=1,2
            IF (ErrorFlag .NE. NOERROR) THEN  !Tube#, velocity Deviation from mean value
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
          END DO

          IsUniformVelProfile=.TRUE.
          CALL GetObjectItem('ODCcktVelocityProfile',1,Alphas,NumAlphas, &
          TmpNumbers,NumNumbers,Status)
          Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

          DO I=Nt*(Nl-1)+1,Nt*Nl !last row faces air inlet (Cross flow HX)
            DO J=1, NumOfMods
              Tube(I)%Seg(J)%VelDev = Numbers(J)  !RS Comment: Bringing in the velocity deviation values
            END DO
            IF (ErrorFlag .NE. NOERROR) THEN
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
            IF (IsUniformVelProfile) THEN
              DO J=1,NumOfMods
                IF (Tube(I)%Seg(J)%VelDev .NE. 1) THEN
                  IsUniformVelProfile=.FALSE.
                  EXIT
                END IF
              END DO
            END IF
          END DO

          !Synchronize 1-D and 2-D arrays
          DO I=1, Nl
            DO J=1, Nt
              Tube2D(I,J)=Tube((I-1)*Nt+J)
            END DO
          END DO

          !Propagate velocity profile to suceeding rows
          DO I=Nl-1,1,-1
            DO J=1, Nt
              DO k=1, NumOfMods
                Tube2D(I,J)%Seg(k)%VelDev=Tube2D(I+1,J)%Seg(k)%VelDev
              END DO
            END DO
          END DO

          IF (ErrorFlag .NE. NOERROR) THEN
            ErrorFlag=CKTFILEERROR
            CALL InitCondenserCoil_Helper_1
            RETURN
          END IF

          !Determine inlet and outlet flags, split, joint or nothing
          !SubcoolingTube=0 !Subcooling node number
          Nnode=1
          DO I=1, NumOfCkts
            !Initialize
            Ckt(I)%InJoin=0
            Ckt(I)%InSplit=0
            Ckt(I)%OutJoin=0
            Ckt(I)%OutSplit=0
            DO J=1, NumOfCkts
              IF (Ckt(I)%TubeSequence(1) .EQ. Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
                Ckt(I)%InJoin=Ckt(I)%InJoin+1
              END IF
              IF (Ckt(I)%TubeSequence(1) .EQ. Ckt(J)%TubeSequence(1)) THEN
                Ckt(I)%InSplit=Ckt(I)%InSplit+1
              END IF
              IF (Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
                Ckt(I)%OutJoin=Ckt(I)%OutJoin+1
              END IF
              IF (Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Ckt(J)%TubeSequence(1)) THEN
                Ckt(I)%OutSplit=Ckt(I)%OutSplit+1
              END IF
            END DO

            IF (Ckt(I)%InJoin .GT. 1 .OR. Ckt(I)%OutSplit .GT. 1) THEN
              Nnode=Nnode+1
            END IF
          END DO

          IF (.NOT. ALLOCATED(Node)) THEN
            ALLOCATE(Node(Nnode))
          END IF

          !Find split and joint tube numbers
          J=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%InJoin .GT. 1 .AND. Ckt(I)%TubeSequence(1) .NE. SubcoolingTube) THEN
              J=J+1
              Node(J)%Num=Ckt(I)%TubeSequence(1)
            ELSEIF (Ckt(I)%OutSplit .GT. 1 .AND. Ckt(I)%TubeSequence(Ckt(I)%Ntube) .NE. SubcoolingTube) THEN
              J=J+1
              Node(J)%Num=Ckt(I)%TubeSequence(Ckt(I)%Ntube)
            END IF
            IF (J .GT. Nnode) THEN
              EXIT
            END IF
          END DO

          Node(Nnode)%Num=0 !coil outlet

          !Find surrounding tubes
          DO I=1, Nl
            DO J=1, Nt
              IF (FinType .EQ. 4 .OR. FinType .EQ. 7 .OR. FinType .EQ. 6) THEN
                Tube2D(I,J)%RowNum=Nl+1-I !Corrected by ISI 07/11/06
              ELSE
                Tube2D(I,J)%RowNum=0
              END IF
              IF (ShiftTube .EQ. 0) THEN
                IF (MOD(I,2) .EQ. 1) THEN
                  Tube2D(I,J)%Fup=I*Nt+(J-1)
                  Tube2D(I,J)%Fdown=I*Nt+J
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 1 .AND. J .EQ. 1) THEN
                    Tube2D(I,J)%Fup=0 !Odd row, first tube
                  END IF
                ELSE
                  Tube2D(I,J)%Fup=I*Nt+J
                  Tube2D(I,J)%Fdown=I*Nt+(J+1)
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 0 .AND. J .EQ. Nt) THEN
                    Tube2D(I,J)%Fdown=0 !even row, first tube
                  END IF
                END IF
              ELSE
                IF (MOD(I,2) .EQ. 1) THEN
                  Tube2D(I,J)%Fup=I*Nt+J
                  Tube2D(I,J)%Fdown=I*Nt+(J+1)
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 1 .AND. J .EQ. Nt) THEN
                    Tube2D(I,J)%Fdown=0 !Odd row, first tube
                  END IF
                ELSE
                  Tube2D(I,J)%Fup=I*Nt+(J-1)
                  Tube2D(I,J)%Fdown=I*Nt+J
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 0 .AND. J .EQ. 1) THEN
                    Tube2D(I,J)%Fup=0 !even row, first tube
                  END IF
                END IF
              END IF

              IF (I .EQ. Nl) THEN
                Tube2D(I,J)%Fup=0
                Tube2D(I,J)%Fdown=0
              END IF
              IF (I .EQ. 1) THEN
                Tube2D(I,J)%Back=0
              END IF
            END DO !End of J
          END DO !End of I

          !Synchronize 1-D and 2-D arrays
          DO I=1, Nl
            DO J=1, Nt
              Tube((I-1)*Nt+J)=Tube2D(I,J)
            END DO
          END DO

          !Find even tubes
          DO I=1, NumOfCkts

            !Find first and last simulation tubes
            FirstTube=1
            LastTube=Ckt(I)%Ntube
            IF (Ckt(I)%InSplit .GT. 1) THEN
              FirstTube=2 !Skip first tube
            END IF
            IF (Ckt(I)%OutJoin .GT. 1) THEN
              LastTube=Ckt(I)%Ntube-1 !Ignore last tube
            END IF

            DO J=FirstTube, LastTube
              TubeNum=Ckt(I)%TubeSequence(J)
              Tube(TubeNum)%Even=0
              IF (FirstTube .EQ. 2 ) THEN
                IF (MOD(J,2) .EQ. 1) THEN
                  Tube(TubeNum)%Even=1
                END IF
              ELSE
                IF (MOD(J,2) .EQ. 0) THEN
                  Tube(TubeNum)%Even=1
                END IF
              END IF
            END DO !End of J
          END DO !End of I

          !Find empty tubes
          Tube%Empty = 1 !Initialize
          DO I=1, NumOfCkts
            DO J=1, Ckt(I)%Ntube
              TubeNum=Ckt(I)%TubeSequence(J)
              Tube(TubeNum)%Empty=0
            END DO
          END DO

          !Number of inlet circuits
          NcktFirst=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%InJoin .LT. 1) THEN
              NcktFirst=NcktFirst+1
            END IF
          END DO
          IF (NcktFirst .EQ. 0) THEN
            NcktFirst = 1 !At least one circuit, ISI - 07/28/06
          END IF

          !Number of outlet circuits
          NcktLast=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%OutSplit .EQ. 0) THEN
              NcktLast=NcktLast+1
            END IF
          END DO
          IF (NcktLast .EQ. 0) THEN
            NcktLast = 1 !At least one circuit, ISI - 07/28/06
          END IF

        END IF !RS Comment: Adding in an END IF to close the above open block (Simple Condenser or not)

      ELSE !IDC ckt

        !**************************** Geometry *************************************

        !Reading in the values for the variables
        CALL GetObjectItem('IDCcktGeometry',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        SELECT CASE (Alphas(1)(1:1))
        CASE ('F','f')
          IsSIunit=.FALSE.
        CASE DEFAULT
          IsSIunit=.TRUE.
        END SELECT

        !Defining variables
        FinType = Numbers(1)
        FinPitch = Numbers(2)
        Kfin = Numbers(3)
        FinThk = Numbers(4)
        TubeType = Numbers(5)
        ODtube = Numbers(6)
        IDtube = Numbers(7)
        Ktube = Numbers(8)
        Pt = Numbers(9)
        Pl = Numbers(10)
        Nl = Numbers(11)
        Nt = Numbers(12)
        Ltube = Numbers(13)

        IF (Ltube .LE. 1e-3) THEN
          ErrorFlag=ZEROLENCOILERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        NumOfMods = Numbers(14)
        NumOfCkts = Numbers(15)

        SELECT CASE (Alphas(5)(1:1))    !Tube Shift Flag
        CASE ('F','f')
          IsShift=.FALSE.
        CASE DEFAULT
          IsShift=.TRUE.
        END SELECT

        !************************** Inputs ****************************************

        !***************** Indoor coil data ***************** !RS: Debugging: Evaporator & Condenser

        CALL GetObjectItem('IndoorCoilData',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        CondPAR%CondFinType = Numbers(1)  !IDC_FinType    !RS: Debugging: Formerly CondPAR(29)

        IDC_FinName = Alphas(1)

        CondPAR%CondFinPitch = Numbers(2)  !RS: Debugging: Formerly CondPAR(22)
        CondPAR%CondFinThermCon = Numbers(3) !Fin Conductivity    !RS: Debugging: Formerly CondPAR(23)
        CondPAR%CondFinThick = Numbers(4)   !Fin Thickness !RS: Debugging: Formerly CondPAR(21)
        CondPAR%CondTube = Numbers(5) !Numerical Denotion of the tube type    !RS: Debugging: Formerly CondPAR(37)
        IDC_TubeID = Numbers(6)   !Tube Inner Diameter
        CondPAR%CondCoilTOD = Numbers(7)   !Tube Outer Diameter   !RS: Debugging: Formerly CondPAR(15)
        CondPAR%CondCoilTThermCon = Numbers(8)    !Tube Conductivity    !RS: Debugging: Formerly CondPAR(18)
        CondPAR%CondRspc = Numbers(9)   !Tube Lateral Spacing  !RS: Debugging: Formerly CondPAR(20)
        CondPAR%CondTspc = Numbers(10)   !Tube Vertical Spacing    !RS: Debugging: Formerly CondPAR(19)
        CondPAR%CondNl = Numbers(11)  !Number of Rows    !RS: Debugging: Formerly CondPAR(25)
        CondPAR%CondNt = Numbers(12)  !Number of Tubes Per Row   !RS: Debugging: Formerly CondPAR(24)
        CondPAR%CondNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly CondPAR(26)
        CondPAR%CondNumMod = Numbers(14)    !Number of Segments  !RS: Debugging: Formerly CondPAR(28)
        CondPAR%CondCoilSTLen = Numbers(15)   !Length of Tube   !RS: Debugging: Formerly CondPAR(17)
        CondPAR%CondMultRefQT = Numbers(16)   !Ref Side Heat Transfer Multiplier    !RS: Debugging: Formerly CondPAR(30)
        CondPAR%CondMultRefPD = Numbers(17) !Ref Side Pressure Drop Multiplier  !RS: Debugging: Formerly CondPAR(31)
        CondPAR%CondMultAirQT = Numbers(18)   !Air Side Heat Transfer Multiplier    !RS: Debugging: Formerly CondPAR(32)
        CondPAR%CondMultAirPD = Numbers(19) !Air Side Pressure Drop Multiplier  !RS: Debugging: Formerly CondPAR(33)

        IF (UNIT .EQ. SI) THEN  !RS: Debugging: Unit conversions
          CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/1000   !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
          CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick/1000   !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
          CondPAR%CondTspc=CondPAR%CondTspc/1000   !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
          CondPAR%CondRspc=CondPAR%CondRspc/1000   !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
          CondPAR%CondFinThick=CondPAR%CondFinThick/1000   !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
          CondPAR%CondFanPwr=CondPAR%CondFanPwr/1000   !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
        ELSE
          CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/12*UnitL       !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
          CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*0.001/12*UnitL !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
          CondPAR%CondCoilSTLen=CondPAR%CondCoilSTLen/12*UnitL       !Tube length, m  !RS: Debugging: Formerly CondPAR(17)
          CondPAR%CondCoilTThermCon=CondPAR%CondCoilTThermCon*UnitK          !Tube thermal conductivity, kW/m-C   !RS: Debugging: Formerly CondPAR(18)
          CondPAR%CondTspc=CondPAR%CondTspc/12*UnitL       !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
          CondPAR%CondRspc=CondPAR%CondRspc/12*UnitL       !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
          CondPAR%CondFinThick=CondPAR%CondFinThick*0.001/12*UnitL !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
          CondPAR%CondFinPitch=CondPAR%CondFinPitch*12/UnitL       !Fin pitch, fin/m    !RS: Debugging: Formerly CondPAR(22)
          CondPAR%CondFinThermCon=CondPAR%CondFinThermCon*UnitK          !Fin thermal conductivity, kW/m-C    !RS: Debugging: Formerly CondPAR(23)
          CondPAR%CondFanPwr=CondPAR%CondFanPwr*1E-3           !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
          !CondPAR%CondBarPress=CondPAR%CondBarPress*UnitP          !Barometric pressure, kPa    !RS: Debugging: Formerly CondPAR(38)
        END IF

        !Tube wall thickness, mm or mil
        CondPAR%CondCoilTWThick=(CondPAR%CondCoilTOD-IDC_TubeID)/2    !RS: Debugging: Formerly CondPAR(16)
        IF (Unit .EQ. IP) THEN
          CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*1000
        END IF

        CondPAR%CondCoolMode=IsCoolingMode   !RS: Debugging: Formerly CondPAR(27)

        !***************** Indoor fan data *****************  !RS: Debugging: Moving: Evaporator & Condenser

        CALL GetObjectItem('IndoorFanData',1,Alphas,NumAlphas, &
        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        CondPAR%CondFanPwr = Numbers(1) !Fan Power   !RS: Debugging: Formerly CondPAR(34)
        !VdotIDfan = Numbers(2)    !Fan Air Flow Rate
        CondPAR%CondFanLoc = Numbers(3)   !Draw Through or Blow Through  !RS: Debugging: Formerly CondPAR(35)

        !IF (UNIT .EQ. SI) THEN  !RS: Debugging: Unit conversions
        !   CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/1000   !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
        !      CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick/1000   !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
        !      CondPAR%CondTspc=CondPAR%CondTspc/1000   !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
        !      CondPAR%CondRspc=CondPAR%CondRspc/1000   !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
        !      CondPAR%CondFinThick=CondPAR%CondFinThick/1000   !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
        !   CondPAR%CondFanPwr=CondPAR%CondFanPwr/1000   !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
        !  ELSE
        !   CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/12*UnitL       !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
        !      CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*0.001/12*UnitL !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
        !      CondPAR%CondCoilSTLen=CondPAR%CondCoilSTLen/12*UnitL       !Tube length, m  !RS: Debugging: Formerly CondPAR(17)
        !      CondPAR%CondCoilTThermCon=CondPAR%CondCoilTThermCon*UnitK          !Tube thermal conductivity, kW/m-C   !RS: Debugging: Formerly CondPAR(18)
        !      CondPAR%CondTspc=CondPAR%CondTspc/12*UnitL       !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
        !      CondPAR%CondRspc=CondPAR%CondRspc/12*UnitL       !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
        !      CondPAR%CondFinThick=CondPAR%CondFinThick*0.001/12*UnitL !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
        !      CondPAR%CondFinPitch=CondPAR%CondFinPitch*12/UnitL       !Fin pitch, fin/m    !RS: Debugging: Formerly CondPAR(22)
        !      CondPAR%CondFinThermCon=CondPAR%CondFinThermCon*UnitK          !Fin thermal conductivity, kW/m-C    !RS: Debugging: Formerly CondPAR(23)
        !   CondPAR%CondFanPwr=CondPAR%CondFanPwr*1E-3           !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
        !   !CondPAR%CondBarPress=CondPAR%CondBarPress*UnitP          !Barometric pressure, kPa    !RS: Debugging: Formerly CondPAR(38)
        !  END IF



        !*************************** Circuiting ************************************

        CALL FinTubeCoilUnitConvert(IsSIUnit,FinPitch,Kfin,FinThk, &
        ODtube,IDtube,Ktube,Pt,Pl,Ltube)

        TubeThk=(ODtube-IDtube)/2

        IF (IsShift) THEN
          ShiftTube = 1
        ELSE
          ShiftTube = 0
        END IF

        IF (Pl .EQ. 0) THEN
          Pl = Pt
        END IF

        !Branch#,#Tubes
        IF (ErrorFlag .NE. NOERROR) THEN
          ErrorFlag=CKTFILEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        NumOfTubes=Nl*Nt    !Number of rows * Number of tubes per row

        !Fin spacing
        FinSpg = 1/FinPitch-FinThk

        !For plate finned tube
        FinHeight=0 !No use for plate finned tube
        TubeDepth=ODtube
        TubeHeight=ODtube
        NumOfChannels=1
        Dchannel=IDtube

        IF (FinSpg .LT. FinThk) THEN
          ErrorFlag=COILFINERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        IF (Pt .LT. ODtube+2*FinThk) THEN
          ErrorFlag=COILTUBEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF

        IF (.NOT. ALLOCATED(Ckt)) THEN
          CALL InitCondenserStructures(TempNumofMods)
        END IF

        IF (IsSimpleCoil .NE. 1) THEN

          !IF (IsSimpleCoil .EQ. 1) THEN
          !    NumOfMods=3
          !    ALLOCATE(Ckt(NumOfCkts))
          !    DO I=1, NumOfCkts
          !        Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
          !        ALLOCATE(Ckt(I)%Tube(1))
          !        ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
          !    END DO
          !ELSE
          !
          !    !ISI - 09/10/07
          !    ALLOCATE(CoilSection(NumOfSections))
          !
          !    ALLOCATE(Tube(NumOfTubes))
          !    ALLOCATE(Tube2D(Nl,Nt))
          !    ALLOCATE(JoinTubes(NumOfTubes))
          !
          !    DO I=1, NumOfTubes
          !        ALLOCATE(Tube(I)%Seg(NumOfMods))
          !    END DO
          !
          !    DO I=1, Nl
          !        DO J=1, Nt
          !            ALLOCATE(Tube2D(I,J)%Seg(NumOfMods))
          !        END DO
          !    END DO
          !
          !    ALLOCATE(Ckt(NumOfCkts))
          !    ALLOCATE(mRefIter(NumOfCkts))
          !    CALL GetObjectItem('IDCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
          !                        TmpNumbers,NumNumbers,Status)
          !    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
          !
          !    DO I=1, NumOfCkts
          !    Ckt(I)%Ntube=Numbers(I)
          !    IF (ErrorFlag .NE. NOERROR) THEN
          !        ErrorFlag=CKTFILEERROR
          !        CALL InitCondenserCoil_Helper_1
          !        RETURN
          !    END IF
          !    ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
          !    ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
          !    END DO

          !Check if all circuit have the same number of tubes !ISI - 09/12/06
          IsSameNumOfTubes=.TRUE.
          Ntubes=Ckt(1)%Ntube
          DO I=2, NumOfCkts
            IF (Ntubes .NE. Ckt(I)%Ntube) THEN
              IsSameNumOfTubes=.FALSE.
              EXIT
            END IF
          END DO

          IF (ErrorFlag .NE. NOERROR) THEN
            ErrorFlag=CKTFILEERROR
            CALL InitCondenserCoil_Helper_1
            RETURN
          END IF

          DO I=1, NumOfCkts
            IF (I .EQ. 1) THEN
              CALL GetObjectItem('IDCcktCircuit1_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 2) THEN
              CALL GetObjectItem('IDCcktCircuit2_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 3) THEN
              CALL GetObjectItem('IDCcktCircuit3_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 4) THEN
              CALL GetObjectItem('IDCcktCircuit4_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 5) THEN
              CALL GetObjectItem('IDCcktCircuit5_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 6) THEN
              CALL GetObjectItem('IDCcktCircuit6_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 7) THEN
              CALL GetObjectItem('IDCcktCircuit7_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSEIF (I .EQ. 8) THEN
              CALL GetObjectItem('IDCcktCircuit8_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            ELSE
              CALL GetObjectItem('IDCcktCircuit9_TubeSequence',1,Alphas,NumAlphas, &
              TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            END IF
            Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            DO J=1, Ckt(I)%Ntube
              Ckt(I)%TubeSequence(J)=Numbers(J)   !RS Comment: Populating the tube sequence arrays
            END DO
            IF (ErrorFlag .NE. NOERROR) THEN
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
          END DO

          !************************* Velocity Profile ********************************

          CoilSection(NumOfSections)%NumOfCkts=NumOfCkts

          DO I=1,2
            IF (ErrorFlag .NE. NOERROR) THEN  !Tube#, velocity Deviation from mean value
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
          END DO

          IsUniformVelProfile=.TRUE.
          DO I=Nt*(Nl-1)+1,Nt*Nl !last row faces air inlet (Cross flow HX)
            CALL GetObjectItem('IDCcktVelocityProfile',1,Alphas,NumAlphas, &
            TmpNumbers,NumNumbers,Status)
            Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            DO J=1,NumOfMods    !RS: Debugging: Adding in a J-value
              Tube(I)%Seg(J)%VelDev = Numbers(J)  !RS Comment: Bringing in the velocity deviation values
            END DO
            IF (ErrorFlag .NE. NOERROR) THEN
              ErrorFlag=CKTFILEERROR
              CALL InitCondenserCoil_Helper_1
              RETURN
            END IF
            IF (IsUniformVelProfile) THEN
              DO J=1,NumOfMods
                IF (Tube(I)%Seg(J)%VelDev .NE. 1) THEN
                  IsUniformVelProfile=.FALSE.
                  EXIT
                END IF
              END DO
            END IF
          END DO

          !Synchronize 1-D and 2-D arrays
          DO I=1, Nl
            DO J=1, Nt
              Tube2D(I,J)=Tube((I-1)*Nt+J)
            END DO
          END DO

          !Propagate velocity profile to suceeding rows
          DO I=Nl-1,1,-1
            DO J=1, Nt
              DO k=1, NumOfMods
                Tube2D(I,J)%Seg(k)%VelDev=Tube2D(I+1,J)%Seg(k)%VelDev
              END DO
            END DO
          END DO

          IF (ErrorFlag .NE. NOERROR) THEN
            ErrorFlag=CKTFILEERROR
            CALL InitCondenserCoil_Helper_1
            RETURN
          END IF

          !Determine inlet and outlet flags, split, joint or nothing
          Nnode=1
          DO I=1, NumOfCkts
            !Initialize
            Ckt(I)%InJoin=0
            Ckt(I)%InSplit=0
            Ckt(I)%OutJoin=0
            Ckt(I)%OutSplit=0
            DO J=1, NumOfCkts
              IF (Ckt(I)%TubeSequence(1) .EQ. Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
                Ckt(I)%InJoin=Ckt(I)%InJoin+1
              END IF
              IF (Ckt(I)%TubeSequence(1) .EQ. Ckt(J)%TubeSequence(1)) THEN
                Ckt(I)%InSplit=Ckt(I)%InSplit+1
              END IF
              IF (Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
                Ckt(I)%OutJoin=Ckt(I)%OutJoin+1
              END IF
              IF (Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Ckt(J)%TubeSequence(1)) THEN
                Ckt(I)%OutSplit=Ckt(I)%OutSplit+1
              END IF
            END DO

            IF (Ckt(I)%InJoin .GT. 1 .OR. Ckt(I)%OutSplit .GT. 1) THEN
              Nnode=Nnode+1
            END IF
          END DO

          IF (.NOT. ALLOCATED(Node)) THEN
            ALLOCATE(Node(Nnode))
          END IF

          !Find split and joint tube numbers
          J=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%InJoin .GT. 1 .AND. Ckt(I)%TubeSequence(1) .NE. SubcoolingTube) THEN
              J=J+1
              Node(J)%Num=Ckt(I)%TubeSequence(1)
            ELSEIF (Ckt(I)%OutSplit .GT. 1 .AND. Ckt(I)%TubeSequence(Ckt(I)%Ntube) .NE. SubcoolingTube) THEN
              J=J+1
              Node(J)%Num=Ckt(I)%TubeSequence(Ckt(I)%Ntube)
            END IF
            IF (J .GT. Nnode) THEN
              EXIT
            END IF
          END DO

          Node(Nnode)%Num=0 !coil outlet

          !Find surrounding tubes
          DO I=1, Nl
            DO J=1, Nt
              IF (FinType .EQ. 4 .OR. FinType .EQ. 7 .OR. FinType .EQ. 6) THEN
                Tube2D(I,J)%RowNum=Nl+1-I !Corrected by ISI 07/11/06
              ELSE
                Tube2D(I,J)%RowNum=0
              END IF
              IF (ShiftTube .EQ. 0) THEN
                IF (MOD(I,2) .EQ. 1) THEN
                  Tube2D(I,J)%Fup=I*Nt+(J-1)
                  Tube2D(I,J)%Fdown=I*Nt+J
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 1 .AND. J .EQ. 1) THEN
                    Tube2D(I,J)%Fup=0 !Odd row, first tube
                  END IF
                ELSE
                  Tube2D(I,J)%Fup=I*Nt+J
                  Tube2D(I,J)%Fdown=I*Nt+(J+1)
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 0 .AND. J .EQ. Nt) THEN
                    Tube2D(I,J)%Fdown=0 !even row, first tube
                  END IF
                END IF
              ELSE
                IF (MOD(I,2) .EQ. 1) THEN
                  Tube2D(I,J)%Fup=I*Nt+J
                  Tube2D(I,J)%Fdown=I*Nt+(J+1)
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 1 .AND. J .EQ. Nt) THEN
                    Tube2D(I,J)%Fdown=0 !Odd row, first tube
                  END IF
                ELSE
                  Tube2D(I,J)%Fup=I*Nt+(J-1)
                  Tube2D(I,J)%Fdown=I*Nt+J
                  Tube2D(I,J)%Back=1
                  IF (MOD(I,2) .EQ. 0 .AND. J .EQ. 1) THEN
                    Tube2D(I,J)%Fup=0 !even row, first tube
                  END IF
                END IF
              END IF

              IF (I .EQ. Nl) THEN
                Tube2D(I,J)%Fup=0
                Tube2D(I,J)%Fdown=0
              END IF
              IF (I .EQ. 1) THEN
                Tube2D(I,J)%Back=0
              END IF
            END DO !End of J
          END DO !End of I

          !Synchronize 1-D and 2-D arrays
          DO I=1, Nl
            DO J=1, Nt
              Tube((I-1)*Nt+J)=Tube2D(I,J)
            END DO
          END DO

          !Find even tubes
          DO I=1, NumOfCkts

            !Find first and last simulation tubes
            FirstTube=1
            LastTube=Ckt(I)%Ntube
            IF (Ckt(I)%InSplit .GT. 1) THEN
              FirstTube=2 !Skip first tube
            END IF
            IF (Ckt(I)%OutJoin .GT. 1) THEN
              LastTube=Ckt(I)%Ntube-1 !Ignore last tube
            END IF

            DO J=FirstTube, LastTube
              TubeNum=Ckt(I)%TubeSequence(J)
              Tube(TubeNum)%Even=0
              IF (FirstTube .EQ. 2 ) THEN
                IF (MOD(J,2) .EQ. 1) THEN
                  Tube(TubeNum)%Even=1
                END IF
              ELSE
                IF (MOD(J,2) .EQ. 0) THEN
                  Tube(TubeNum)%Even=1
                END IF
              END IF
            END DO !End of J
          END DO !End of I

          !Find empty tubes
          Tube%Empty = 1 !Initialize
          DO I=1, NumOfCkts
            DO J=1, Ckt(I)%Ntube
              TubeNum=Ckt(I)%TubeSequence(J)
              Tube(TubeNum)%Empty=0
            END DO
          END DO

          !Number of inlet circuits
          NcktFirst=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%InJoin .LT. 1) THEN
              NcktFirst=NcktFirst+1
            END IF
          END DO
          IF (NcktFirst .EQ. 0) THEN
            NcktFirst = 1 !At least one circuit, ISI - 07/28/06
          END IF

          !Number of outlet circuits
          NcktLast=0
          DO I=1, NumOfCkts
            IF (Ckt(I)%OutSplit .EQ. 0) THEN
              NcktLast=NcktLast+1
            END IF
          END DO
          IF (NcktLast .EQ. 0) THEN
            NcktLast = 1 !At least one circuit, ISI - 07/28/06
          END IF

        END IF !RS Comment: Adding in an END IF to close the above open block (Simple Condenser or not)
      END IF  !End of the IDC or ODC if-statement

    ELSE !Microchannel coil

      !**************************** Geometry *************************************

      !RS Comment: Updating input data for the microchannel option
      !Reading in the values for the variables
      CALL GetObjectItem('ODCcktGeometry',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      SELECT CASE (Alphas(1)(1:1))
      CASE ('F','f')
        IsSIunit=.FALSE.
      CASE DEFAULT
        IsSIunit=.TRUE.
      END SELECT

      !Defining variables
      FinType = Numbers(1)
      FinPitch = Numbers(2)
      Kfin = Numbers(3)
      FinThk = Numbers(4)
      TubeType = Numbers(5)
      TubeHeight = Numbers(6)
      TubeDepth = Numbers(7)
      Tubethk = Numbers(8)
      Ktube = Numbers(9)
      Pt = Numbers(10)
      Nl = Numbers(11)
      Nt = Numbers(12)
      Ltube = Numbers(13)

      IF (Ltube .LE. 1e-3) THEN
        ErrorFlag=ZEROLENCOILERROR
        CALL InitCondenserCoil_Helper_1
        RETURN
      END IF

      !SELECT CASE (Alphas(5)(1:1))   !RS: Debugging: Set once but never used
      !CASE ('V','v')
      !    TubeOrientation=VERTICAL
      !CASE ('H','h')
      !    TubeOrientation=HORIZONTAL
      !CASE DEFAULT
      !    TubeOrientation=HORIZONTAL
      !END SELECT

      NumOfMods = Numbers(14) !Number of segments or modules
      NumOfChannels = Numbers(15) !Number of channels
      Dchannel = Numbers(16)

      !*************************** Circuiting ************************************

      CALL GetObjectItem('ODCcktCircuiting_Slab1',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers=DBLE(TmpNumbers)    !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      ALLOCATE(Slab(Nl))
      !Slab#,#Passes,Tubes per Pass
      DO I=1, Nl

        NumOfPasses = Numbers(2)
        ALLOCATE(Slab(I)%Pass(NumOfPasses))
        Slab(I)%Npass=NumOfPasses

        DO II=1, NumOfPasses

          J=2+II      !+2 because the first two Numbers() are Slab# and # of Passes
          Ntubes=Numbers(J)   !Allows for Ntubes to vary for the different passes

          !cooling and heating are different flow direction - ISI 02/06/2009
          IF (IsCoolingMode .GT. 0) THEN
            ALLOCATE(Slab(I)%Pass(II)%Tube(Ntubes))
            Slab(I)%Pass(II)%Ntube=Ntubes

            DO III=1, Ntubes
              ALLOCATE(Slab(I)%Pass(II)%Tube(III)%Seg(NumOfMods))
            END DO

          ELSE
            ALLOCATE(Slab(I)%Pass(NumOfPasses-II+1)%Tube(Ntubes))
            Slab(I)%Pass(NumOfPasses-II+1)%Ntube=Ntubes

            DO III=1, Ntubes
              ALLOCATE(Slab(I)%Pass(NumOfPasses-II+1)%Tube(III)%Seg(NumOfMods))
            END DO

          END IF

        END DO

      END DO

      !Inlet pass
      !Slab#,#Inlets,Tubes per Inlet

      DO I=1, Nl

        NumOfInlets = Numbers(5)
        ALLOCATE(Slab(I)%InletPass(NumOfInlets))
        Slab(I)%Ninlet=NumOfInlets

        DO II=1, NumOfInlets

          Ntubes = Numbers(6)

          IF (IsCoolingMode .GT. 0) THEN
            Slab(I)%InletPass(II)%Ntube=Ntubes
          ELSE !For heating, set it to equal number of inlet tubes, at least for now - ISI - 02/06/2009
            Slab(I)%InletPass(II)%Ntube=Slab(I)%Pass(1)%Ntube
          END IF

        END DO

      END DO

      !************************* Velocity Profile ********************************

      !Tube# ,velocity Deviation from mean value

      CALL GetObjectItem('ODCcktVelocityProfile',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers=DBLE(TmpNumbers)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      IsUniformVelProfile=.TRUE.
      DO II=1,Slab(Nl)%Npass
        DO III=1,Slab(Nl)%Pass(II)%Ntube
          DO IV=1, NumOfMods
            Slab(Nl)%Pass(II)%Tube(III)%Seg(IV)%VelDev=Numbers(IV)  !Will need to be rewritten for the passes---is reading the same five numbers every time
          END DO
          IF (IsUniformVelProfile) THEN
            DO IV=1,NumOfMods
              IF (Slab(Nl)%Pass(II)%Tube(III)%Seg(IV)%VelDev .NE. 1) THEN
                IsUniformVelProfile=.FALSE.
                EXIT
              END IF
            END DO
          END IF
        END DO
      END DO

      CALL MicroChannelCoilUnitConvert(IsSIUnit,FinPitch,Kfin,FinThk, &
      TubeHeight,TubeDepth,TubeThk,Ktube, &
      Pt,Ltube,Dchannel)

      ODtube=TubeHeight
      IDtube=Dchannel

      FinHeight=Pt-TubeHeight
      FinSpg=1/FinPitch-FinThk

    END IF

    !CLOSE(11) !Circuit file    !RS: Debugging: Removing because the file number isn't actually attached to the file (5/21/14)

    !Discharge line info
    IDdisLn=ODdisLn-DisLnThk*2
    LmodDis=Ldisln
    AiModDis=PI*IDdisLn*LmodDis

    !liquid line info
    IDliqLn=ODliqLn-LiqLnThk*2
    LmodLiq=Lliqln
    AiModLiq=PI*IDliqLn*LmodLiq

    !!***** Allocate pointer for discharge lines *****
    !ALLOCATE(DisLnSeg(NumOfMods))
    !
    !!***** Allocate pointer for liquid lines *****
    !ALLOCATE(LiqLnSeg(NumOfMods))
    !
    CALL InitCondenserCoil_Helper_1

    RETURN

  END SUBROUTINE InitCondenserCoil

  SUBROUTINE InitCondenserStructures(TempNumofMods)
    !RS: Debugging: Initialize and allocates the circuit, coil, segment, and tube structures

    USE InputProcessor

    INTEGER, PARAMETER :: MaxNameLength = 200

    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    INTEGER :: TempNumofMods

    !NumOfMods=TempNumofMods

    ALLOCATE(Ckt(NumOfCkts))
    ALLOCATE(CoilSection(NumOfSections))
    ALLOCATE(Tube(NumOfTubes))
    ALLOCATE(Tube2D(Nl,Nt))
    ALLOCATE(JoinTubes(NumOfTubes))
    ALLOCATE(mRefIter(NumOfCkts))
    ALLOCATE(DisLnSeg(NumOfMods))
    ALLOCATE(LiqLnSeg(NumOfMods))

    !NumOfMods=3
    !NumOfMods=TempNumofMods

    DO I=1, NumOfTubes
      ALLOCATE(Tube(I)%Seg(NumOfMods))
    END DO

    DO I=1, Nl
      DO J=1, Nt
        ALLOCATE(Tube2D(I,J)%Seg(NumOfMods))
      END DO
    END DO

    IF (IsCoolingMode .GT. 0) THEN

      CALL GetObjectItem('ODCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      DO I=1, NumOfCkts
        Ckt(I)%Ntube=Numbers(I)
        IF (ErrorFlag .NE. NOERROR) THEN
          ErrorFlag=CKTFILEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF
        ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
        ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
        DO J=1,Ckt(I)%Ntube
          ALLOCATE(Ckt(I)%Tube(J)%Seg(NumOfMods))
        END DO
      END DO

    ELSE
      CALL GetObjectItem('IDCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
      TmpNumbers,NumNumbers,Status)
      Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

      DO I=1, NumOfCkts
        Ckt(I)%Ntube=Numbers(I)
        IF (ErrorFlag .NE. NOERROR) THEN
          ErrorFlag=CKTFILEERROR
          CALL InitCondenserCoil_Helper_1
          RETURN
        END IF
        ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
        ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
        DO J=1,Ckt(I)%Ntube
          ALLOCATE(Ckt(I)%Tube(J)%Seg(NumOfMods))
        END DO
      END DO

    END IF

  END SUBROUTINE InitCondenserStructures

  !************************************************************************

  SUBROUTINE InitCondenserCoil_Helper_1

    IMPLICIT NONE

    IF (ErrorFlag .NE. NOERROR) THEN
      IF (ErrorFlag .EQ. CKTFILEERROR) THEN
        WRITE(*,*)'## ERROR ## Condenser: Circuit file error.'
      ELSEIF (ErrorFlag .EQ. COILTUBEERROR) THEN
        WRITE(*,*)'## ERROR ## Condenser: Coil geometry misdefined.'
        WRITE(*,*)'Tube spacing is less than tube diameter.'
      ELSEIF (ErrorFlag .EQ. COILFINERROR) THEN
        WRITE(*,*)'## ERROR ## Condenser: Coil geometry misdefined.'
        WRITE(*,*)'Fin spacing is less than fin thickness.'
      END IF
    END IF

  END SUBROUTINE InitCondenserCoil_Helper_1

  !************************************************************************

  SUBROUTINE EndCondenserCoil

    !------------------------------------------------------------------------
    !Purpose:
    !To free allocated arrays
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER I,II,III,J !Loop counters

    IF (CoilType .EQ. MCCONDENSER) THEN
      IF (IsSimpleCoil .EQ. 1) THEN   !IsSimpleCoil doesn't seem to really simplify for the microchannel case so much as cause errors
        IsSimpleCoil=2
      END IF
    END IF

    IF (IsSimpleCoil .EQ. 1) THEN
      DO I=1, NumOfCkts
        DEALLOCATE(Ckt(I)%Tube(1)%Seg)
        DEALLOCATE(Ckt(I)%Tube)
      END DO
      DEALLOCATE(Ckt)
      IF (ALLOCATED(DisLnSeg)) THEN
        DEALLOCATE(DisLnSeg) !Discharge line
      END IF
      IF (ALLOCATED(LiqLnSeg)) THEN
        DEALLOCATE(LiqLnSeg) !Line line
      END IF

      RETURN
    END IF

    IF (NumOfChannels .GT. 1) THEN

      IF (ALLOCATED(Slab)) THEN
        DO I=1, Nl
          IF (ErrorFlag .NE. CKTFILEERROR) THEN
            DO II=1,Slab(I)%Npass
              DO III=1,Slab(I)%Pass(II)%Ntube
                DEALLOCATE(Slab(I)%Pass(II)%Tube(III)%Seg)
              END DO !end tube
              DEALLOCATE(Slab(I)%Pass(II)%Tube)
            END DO !end pass
            DEALLOCATE(Slab(I)%Pass)
            DEALLOCATE(Slab(I)%InletPass)
          END IF
        END DO !end slab
        DEALLOCATE(Slab)
      END IF

    ELSE

      IF (ErrorFlag .EQ. CKTFILEERROR) THEN
        DO I=1, NumOfTubes
          DEALLOCATE(Tube(I)%Seg)
        END DO
      END IF

      IF (ALLOCATED(Ckt)) THEN
        DO I=1, NumOfCkts
          IF (ErrorFlag .NE. CKTFILEERROR) THEN
            DO J=1,Ckt(I)%Ntube
              IF (Ckt(I)%OutJoin .LE. 1) THEN
                DEALLOCATE(Ckt(I)%Tube(J)%Seg)
              END IF
            END DO
          END IF
          DEALLOCATE(Ckt(I)%Tube)
          DEALLOCATE(Ckt(I)%TubeSequence)
        END DO
        DEALLOCATE(Ckt)
      END IF

      !ISI - 09/10/07
      IF (ALLOCATED(CoilSection)) THEN
        DO I=1, NumOfSections
          IF (ALLOCATED(CoilSection(I)%CktNum)) THEN
            DEALLOCATE(CoilSection(I)%CktNum)
          END IF
        END DO
        DEALLOCATE(CoilSection)
      END IF

      IF (ALLOCATED(Tube)) THEN
        DEALLOCATE(Tube)
      END IF

      IF (ALLOCATED(Tube2D)) THEN
        DEALLOCATE(Tube2D)
      END IF

      IF (ALLOCATED(JoinTubes)) THEN
        DEALLOCATE(JoinTubes)
      END IF

      IF (ALLOCATED(mRefIter)) THEN
        DEALLOCATE(mRefIter)
      END IF

      IF (ALLOCATED(Node)) THEN
        DEALLOCATE(Node)
      END IF

    END IF

    IF (ALLOCATED(DisLnSeg)) THEN
      DEALLOCATE(DisLnSeg) !Discharge line
    END IF

    IF (ALLOCATED(LiqLnSeg)) THEN
      DEALLOCATE(LiqLnSeg) !Line line
    END IF

    IF (CoilType .EQ. MCCONDENSER) THEN
      IF (IsSimpleCoil .EQ. 2) THEN   !IsSimpleCoil doesn't seem to really simplify for the microchannel case so much as cause errors
        IsSimpleCoil=1
      END IF
    END IF

    RETURN

  END SUBROUTINE EndCondenserCoil

  !************************************************************************

  SUBROUTINE DischargeLine

    !------------------------------------------------------------------------
    !Purpose:
    !To calculation discharge line outlet conditions
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !Andrews, J.W. 2001. �Impact of refrigerant line length on system
    !efficiency in residential heating and cooling system using refrigerant
    !distributor.� BNL-68550, Brookhaven National Laboratory, Upton, NY.
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod

    IMPLICIT NONE

    !INTEGER,PARAMETER :: CoilType = 3  !1=Condenser; 2=Evaporator;
    !                                   !3=High side interconnecting pipes;
    !								    !4=Low side interconnecting pipes
    !								    !5=Microchannel condenser

    INTEGER CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes; 4=Low side interconnecting pipes
    INTEGER TubeType !1=Plain; 2=General Micro Fin; 3=Herringbone; 4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A

    REAL DPman !Manifold pressure drop, kPa

    !FLOW:

    CoilType=HIGHSIDETUBE
    TubeType=SMOOTH

    MassDisLn=0

    IF (DTdisLn .NE. 0) THEN !Given discharge line temperature change
      Pressure=pRoCmp*1000    !RS Comment: Unit Conversion
      Enthalpy=hRoCmp*1000    !RS Comment: Unit Conversion
      CpgRoCmp=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)
      CpgRoCmp=CpgRoCmp/1000  !RS Comment: Unit Conversion
      tRoCmp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)
      IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 1869'
        ErrorFlag=REFPROPERROR
        RETURN
      END IF

      QdisLn=mRefTot*CpgRoCmp*(-DTdisLn)
    END IF

    DO K=1,1

      IF (K .EQ. 1) THEN
        !Equal to compressor outlet condition
        DisLnSeg(K)%pRi=pRoCmp
        DisLnSeg(K)%hRi=hRoCmp
      ELSE !Equal to outlet of previous module(section)
        DisLnSeg(K)%pRi=DisLnSeg(K-1)%pRo
        DisLnSeg(K)%hRi=DisLnSeg(K-1)%hRo
      END IF

      pRiMod=DisLnSeg(K)%pRi
      hRiMod=DisLnSeg(K)%hRi

      CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
      vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
      muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      !Discharge line heat loss
      !        Tloss = PI*LmodDis*ODdisLn*Hout*(tRiMod-tAoCoil)/(mRefTot*cpgRiMod)
      !        Qloss = mRefTot*cpgRiMod*Tloss
      !        QdisLn = QdisLn+Qloss

      hRoMod=-QdisLn/mRefTot+hRiMod

      IF (xRiMod .GT. 1) THEN
        LmodTPratio=0
      ELSE
        LmodTPratio=1
      END IF

      !Find outlet ref. pressure
      CALL CalcSegmentRefOutletPressure(CoilType,TubeType,pRiMod,hgRiMod,hfRiMod, &    !CoilType,TubeType,tRiMod,pRiMod,hgRiMod,hfRiMod, &
      hRiMod,hRoMod,xRiMod,vRiMod,vgRiMod,vfRiMod,mRefTot, &
      muRiMod,mugRiMod,mufRiMod,LmodDis,LmodTPratio, & !muRiMod,mugRiMod,mufRiMod,SigmaMod,LmodDis,LmodTPratio, & !RS: Debugging: Extraneous SigmaMod
      IDdisLn,ElevDisLn,LmodDis,DPrefMultiplier,pRoMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      pRoMod=pRoMod-AddDPdisLn

      CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
      vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
      muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      mu=(muRiMod+muRoMod)/2  !Average viscosity

      CALL Inventory(CoilType,IDdisLn,mRefTot,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &  !RS: Debugging: Extraneous
      xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
      LmodDis,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
      !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
      !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
      !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

      DisLnSeg(K)%Mass=MassMod
      DisLnSeg(K)%pRo=pRoMod
      DisLnSeg(K)%hRo=hRoMod

      MassDisLn=MassDisLn+MassMod

    END DO !End Nmod

    CALL manifold(CoilType,IDdisLn,mRefTot,xRiMod,vgRimod,vfRimod,mugRiMod,mufRiMod,dPman)
    !(CoilType,ID,mRef,xRef,vRef,vgi,vfi,mug,muf,dPman)

    pRiCoil=DisLnSeg(1)%pRo-DPman
    hRiCoil=DisLnSeg(1)%hRo

    RETURN

  END SUBROUTINE DischargeLine

  !************************************************************************

  SUBROUTINE LiquidLine

    !------------------------------------------------------------------------
    !Purpose:
    !To calculation liquid line outlet conditions
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !Andrews, J.W. 2001. �Impact of refrigerant line length on system
    !efficiency in residential heating and cooling system using refrigerant
    !distributor.� BNL-68550, Brookhaven National Laboratory, Upton, NY.
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod

    IMPLICIT NONE

    !INTEGER,PARAMETER :: CoilType = 3  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser

    INTEGER CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes; 4=Low side interconnecting pipes
    INTEGER TubeType !1=Plain; 2=General Micro Fin; 3=Herringbone; 4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
    INTEGER K !Loop counter !ISI - 08/25/06

    CoilType=HIGHSIDETUBE
    TubeType=SMOOTH

    MassLiqLn=0
    IF (DTliqLn .NE. 0) THEN !Given liquid line temperature change
      QliqLn=mRefTot*CpfRoMod*(-DTliqLn)
    END IF

    DO K=1,1 !ISI - 12/17/2009

      IF (K .EQ. 1) THEN
        !Equal to coil outlet condition
        LiqLnSeg(K)%pRi=pRoCoil
        LiqLnSeg(K)%hRi=hRoCoil
      ELSE !Equal to outlet of previous module(section)
        LiqLnSeg(K)%pRi=LiqLnSeg(K-1)%pRo
        LiqLnSeg(K)%hRi=LiqLnSeg(K-1)%hRo
      END IF

      pRiMod=LiqLnSeg(K)%pRi
      hRiMod=LiqLnSeg(K)%hRi

      !Calculate inlet ref. property
      CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
      vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
      muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      hRoMod=-QliqLn/mRefTot+hRiMod   !Determining Module Refrigerant Outlet Enthalpy

      IF (xRiMod .GT. 0) THEN
        LmodTPratio=1
      ELSE
        LmodTPratio=0
      END IF

      !Find outlet ref. pressure
      CALL CalcSegmentRefOutletPressure(CoilType,TubeType,pRiMod,hgRiMod,hfRiMod, & !CoilType,TubeType,tRiMod,pRiMod,hgRiMod,hfRiMod, &
      hRiMod,hRoMod,xRiMod,vRiMod,vgRiMod,vfRiMod,mRefTot, &
      muRiMod,mugRiMod,mufRiMod,LmodLiq,LmodTPratio, & !muRiMod,mugRiMod,mufRiMod,SigmaMod,LmodLiq,LmodTPratio, & !RS: Debugging: Extraneous SigmaMod
      IDliqLn,ElevLiqLn,LmodLiq,DPrefMultiplier,pRoMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      pRoMod=pRoMod-AddDPLiqLn    !Determining Module Refrigerant Outlet Pressure

      !Calculate outlet ref. property
      CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
      vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
      muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      !Average viscosities
      mu=(muRiMod+muRoMod)/2
      muf=(mufRiMod+mufRoMod)/2
      mug=(mugRiMod+mugRoMod)/2

      CALL Inventory(CoilType,IDliqLn,mRefTot,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &  !RS: Debugging: Extraneous
      xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
      Lmodliq,LmodTP,LmodSH,MassLiqMod,MassVapMod,MassMod)
      !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
      !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
      !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)
      LiqLnSeg(K)%Mass=MassMod
      LiqLnSeg(K)%pRo=pRoMod
      LiqLnSeg(K)%hRo=hRoMod

      MassLiqLn=MassLiqLn+MassMod

    END DO !End Nmod

    pRiExp=LiqLnSeg(1)%pRo !ISI - 08/25/06
    hRiExp=LiqLnSeg(1)%hRo !ISI - 08/25/06

    RETURN

  END SUBROUTINE LiquidLine

  !************************************************************************

  SUBROUTINE RefrigerantParameters(Ref$)

    !------------------------------------------------------------------------
    !Purpose:
    !To set refrigerant parameters
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    IMPLICIT NONE

    CHARACTER*80, INTENT(IN)  :: Ref$

    RefName=Ref$

    MolWeight=MW(RefName)*1000  !RS Comment: Unit Conversion

    Tcr=Tcrit(RefName)+273.15   !RS Comment: Unit Conversion
    Pcr=Pcrit(RefName)/1000     !RS Comment: Unit Conversion

  END SUBROUTINE RefrigerantParameters

  !************************************************************************

  SUBROUTINE InitBoundaryConditions(CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To initialize segment boundary conditions
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE AirPropMod
    USE CoilCalcMod

    IMPLICIT NONE


    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser

    !FLOW:

    AirPropOpt=2
    AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APRelHum=rhAiCoil !RS: Debugging: Formerly AirProp(3)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
    hAiCoil=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)
    wbAiCoil=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)

    !air side inlet conditions
    !RS: Replace: Moving after above PsyChart call to update humidity ratio
    !CPair=CPA(REAL(tAiCoil))    !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAiCoil*CPAir

    IF (DrawBlow .EQ. BLOWTHROUGH) THEN !Blow through
      !tAiCoil=tAiCoil+PwrFan/Cair
      !hAiCoil=hAiCoil+PwrFan/mAiCoil
      CALL SimSimpleFan(tAiCoil,hAiCoil,mAiCoil,tAiCoil,hAiCoil)  !RS: Debugging: Adding in condenser fan model (9/4/14)
      CondPAR%CondFanPwr=FanOut%Power
    END IF

    IF (IsCmpInAirStream .NE. 0) THEN !Compressor in air stream
      tAiCoil=tAiCoil+QlossCmp/Cair
      hAiCoil=hAiCoil+QlossCmp/mAiCoil
    END IF

    AirPropOpt=1
    AirProp%APRelHum=tAiCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APEnth=hAiCoil  !RS: Debugging: Formerly AirProp(4)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
    DensityIn=AirProp%APDryDens    !RS: Debugging: Formerly AirProp(7)

    !Area calculations
    CALL CalcCoilarea(CoilType,Nl,Nt,Pt,Pl,TubeDepth, &
    Ltube,IDtube,TubeHeight,Dchannel,NumOfChannels, &
    FinThk,FinSpg,Lcoil,AfCoil, &
    AoCoil,AiCoil,AmCoil)

    !Surface areas
    AoMod=AoCoil*LmodTube/Lcoil
    AfMod=AfCoil*LmodTube/Lcoil
    AiMod=AiCoil*LmodTube/Lcoil
    AmMod=AmCoil*LmodTube/Lcoil

    !****** Discharge line calculation ******
    IF (LdisLn .GT. 0) THEN
      CALL DischargeLine
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        WRITE(*,*)'-- WARNING -- DischargeLine: Refprop error. Line 2216'
        RETURN
      END IF
    ELSE
      pRiCoil=pRoCmp
      hRiCoil=hRoCmp
    END IF

    Pressure=pRiCoil*1000   !RS Comment: Unit Conversion
    Enthalpy=hRiCoil*1000   !RS Comment: Unit Conversion
    tRiCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Coil Inlet Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3904'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF
    xRiCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Coil Inlet Quality
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3910'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF

    IF (IsSimpleCoil .EQ. 1) THEN
      !Initialize
      DO I=1, NumOfCkts
        Ckt(I)%Tube(1)%Seg%Qmod=0
        Ckt(I)%Tube(1)%Seg%tAi=tAiCoil
        Ckt(I)%Tube(1)%Seg%wbAi=wbAiCoil
        Ckt(I)%Tube(1)%Seg%rhAi=rhAiCoil
        Ckt(I)%Tube(1)%Seg%tAo=tAiCoil
        Ckt(I)%Tube(1)%Seg%wbAo=wbAiCoil
        Ckt(I)%Tube(1)%Seg%rhAo=rhAiCoil
        Ckt(I)%Tube(1)%Seg%pRo=pRiCoil
        Ckt(I)%Tube(1)%Seg%hRo=hRiCoil
        Ckt(I)%Tube(1)%Seg%pRi=pRiCoil
        Ckt(I)%Tube(1)%Seg%hRi=hRiCoil

        Ckt(I)%pRi=pRiCoil
        Ckt(I)%pRo=pRiCoil
        Ckt(I)%hRi=hRiCoil
        Ckt(I)%hRo=hRiCoil
        Ckt(I)%mRef=mRefTot/NumOfCkts
      END DO

    ELSE
      !Calc. air side mass flow rate at the front row
      VelAvg=(mAiCoil/DensityIn)/Aface
      DO I=1,Nt*Nl
        DO J=1,NumOfMods
          Tube(I)%Seg(J)%Aface=LmodTube/(Ltube*Nt)*Aface
          Tube(I)%Seg(J)%mAi=mAiCoil*LmodTube/(Ltube*Nt)*Tube(I)%Seg(J)%VelDev
        END DO
      END DO

      !Initialize boundary conditions
      DO I=1, NumOfTubes
        Tube(I)%Seg%Qmod=0
        Tube(I)%Seg%tAi=tAiCoil
        Tube(I)%Seg%wbAi=wbAiCoil
        Tube(I)%Seg%rhAi=rhAiCoil
        Tube(I)%Seg%tAo=tAiCoil
        Tube(I)%Seg%wbAo=wbAiCoil
        Tube(I)%Seg%rhAo=rhAiCoil
        Tube(I)%Seg%pRo=pRiCoil
        Tube(I)%Seg%hRo=hRiCoil
        Tube(I)%Seg%pRi=pRiCoil
        Tube(I)%Seg%hRi=hRiCoil
        Tube(I)%ID=IDtube !ISI - 06/05/07
        Tube(I)%NumOfMods=NumOfMods !ISI - 06/05/07
      END DO

      !Synchronize boundary conditions from tube to circuits
      DO I=1, NumOfCkts
        DO J=1, Ckt(I)%Ntube
          TubeNum=Ckt(I)%TubeSequence(J)
          Ckt(I)%Tube(J)=Tube(TubeNum)
        END DO
        !Initialize
        Ckt(I)%pRi=pRiCoil
        Ckt(I)%pRo=pRiCoil
        Ckt(I)%hRi=hRiCoil
        Ckt(I)%hRo=hRiCoil
      END DO

      mRefIter=Ckt%mRef !Ref. flow rate iteration value

      DO I=1, NumOfCkts

        IF (Ckt(I)%OutSplit .EQ. 0) THEN !Outlet circuit
          Ckt(I)%mRef=mRefTot/NcktLast
        ELSEIF (Ckt(I)%InJoin .EQ. 0) THEN !Inlet circuit
          Ckt(I)%mRef=mRefTot/NcktFirst
        ELSEIF (Ckt(I)%InSplit .GT. 1) THEN !Split inlet
          DO J=1, NumOfCkts
            IF (Ckt(I)%TubeSequence(1) .EQ. Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
              Ckt(I)%mRef=Ckt(J)%mRef/Ckt(I)%InSplit
              EXIT !Found the split tube
            END IF
          END DO
        ELSE IF (Ckt(I)%InJoin .GT. 1) THEN !Joint inlet
          DO J=1, NumOfCkts
            IF (Ckt(J)%TubeSequence(Ckt(J)%Ntube) .EQ. Ckt(I)%TubeSequence(1)) THEN
              Ckt(I)%mRef=Ckt(J)%mRef*Ckt(I)%InJoin
              EXIT !Found the joined tube
            END IF
          END DO
        ELSE
          Ckt(I)%mRef=mRefTot/NumOfCkts !to take care of one tube case, ISI - 07/28/06
        END IF

      END DO

      !mRefTotPrev = MRefTot  !RS: Debugging: Set once but never used

    END IF

    RETURN

  END SUBROUTINE InitBoundaryConditions

  !************************************************************************

  SUBROUTINE CalcCircuitRefInletConditions(I,II,CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To calculate circuit refrigerant inlet condition according to circuitry
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: I   !Slab number
    INTEGER,INTENT(IN) :: II  !Circuit,pass number

    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel Evaporator

    REAL SumMref    !Sum of mdot ref.
    REAL SumMrefHri !Sum of mdot ref x hri
    INTEGER J,K !Loop counters

    !FLOW:

    IF (IsSimpleCoil .EQ. 1) THEN
      Ckt(II)%pRi=pRiCoil
      Ckt(II)%hRi=hRiCoil
      RETURN
    END IF

    IF (CoilType .NE. MCCONDENSER) THEN

      IF (Ckt(II)%InSplit .GT. 1) THEN !Split inlet
        DO J=1, NumOfCkts
          IF (Ckt(J)%TubeSequence(Ckt(J)%Ntube) .EQ. Ckt(II)%TubeSequence(1)) THEN
            IF (Ckt(J)%OutJoin .LE. 1) THEN !No joint at outlet
              !Inlet conditions
              Ckt(II)%pRi=Ckt(J)%pRo
              Ckt(II)%hRi=Ckt(J)%hRo
            ELSE !Outlet has joint
              !Inlet conditions
              Ckt(II)%pRi=Ckt(J)%Tube(Ckt(J)%Ntube-1)%Seg(NumOfMods)%pRo
              Ckt(II)%hRi=Ckt(J)%Tube(Ckt(J)%Ntube-1)%Seg(NumOfMods)%hRo
            END IF
            EXIT !Found the split tube
          END IF
        END DO
      ELSE IF (Ckt(II)%InJoin .GT. 1) THEN !Joint inlet
        K=0
        Ckt(II)%pRi=0
        Ckt(II)%hRi=0
        SumMref=0
        SumMrefHri=0
        DO J=1, NumOfCkts
          IF (Ckt(J)%TubeSequence(Ckt(J)%Ntube) .EQ. Ckt(II)%TubeSequence(1)) THEN
            K=K+1
            !Inlet conditions
            Ckt(II)%pRi=Ckt(II)%pRi+Ckt(J)%pRo
            Ckt(II)%hRi=Ckt(II)%hRi+Ckt(J)%hRo
            JoinTubes(K)=J
            SumMref=SumMref+Ckt(J)%mRef
            SumMrefHri=SumMrefHri+Ckt(J)%mRef*Ckt(J)%hRo
            IF (K .EQ. Ckt(II)%InJoin) THEN
              EXIT !Found all joined tubes
            END IF
          END IF
        END DO
        !Calculate according to energy balance
        Ckt(II)%pRi=Ckt(II)%pRi/Ckt(II)%InJoin !Calc inlet pressure
        Ckt(II)%hRi=SumMrefHri/SumMref       !Calc inlet enthalpy

      ELSE !Coil inlet
        Ckt(II)%pRi=pRiCoil
        Ckt(II)%hRi=hRiCoil
      END IF

    ELSE

      IF (IsParallelSlabs) THEN
        !ISI - 07/13/07
        IF (II .EQ. 1) THEN !1st pass
          Slab(I)%Pass(II)%pRi=pRiCoil
          Slab(I)%Pass(II)%hRi=hRiCoil
        ELSE
          Slab(I)%Pass(II)%pRi=Slab(I)%Pass(II-1)%pRo
          Slab(I)%Pass(II)%hRi=Slab(I)%Pass(II-1)%hRo
        END IF

      ELSE !Series

        IF (I .EQ. 1) THEN !1st slab
          IF (II .EQ. 1) THEN !1st pass
            Slab(I)%Pass(II)%pRi=pRiCoil
            Slab(I)%Pass(II)%hRi=hRiCoil
          ELSE
            Slab(I)%Pass(II)%pRi=Slab(I)%Pass(II-1)%pRo
            Slab(I)%Pass(II)%hRi=Slab(I)%Pass(II-1)%hRo
          END IF
        ELSE
          IF (II .EQ. 1) THEN !1st pass
            Slab(I)%Pass(II)%pRi=Slab(I-1)%Pass(Slab(I)%Npass)%pRo
            Slab(I)%Pass(II)%hRi=Slab(I-1)%Pass(Slab(I)%Npass)%hRo
          ELSE
            Slab(I)%Pass(II)%pRi=Slab(I)%Pass(II-1)%pRo
            Slab(I)%Pass(II)%hRi=Slab(I)%Pass(II-1)%hRo
          END IF
        END IF

      END IF

    END IF

    RETURN

  END SUBROUTINE CalcCircuitRefInletConditions

  !************************************************************************

  SUBROUTINE CalcCoilSegment(I,II,III,IV,CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To perform heat exchanger calculation for a segment
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod
    USE AirPropMod

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: I   !Slab number
    INTEGER,INTENT(IN) :: II  !Circuit,pass number
    INTEGER,INTENT(IN) :: III !Tube number
    INTEGER,INTENT(IN) :: IV  !Segment number

    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel Evaporator

    !FLOW:

    IF (CoilType .NE. MCCONDENSER) THEN

      Ckt(II)%Tube(III)%Seg(IV)%Len=LmodTube

      IF (IsSimpleCoil .EQ. 1) THEN

        IF (LmodTube .LT. SMALL) THEN !For zero length, ISI - 02/08/08

          Ckt(II)%Tube(III)%Seg(IV)%Len=0
          Ckt(II)%Tube(III)%Seg(IV)%Qmod=0
          Ckt(II)%Tube(III)%Seg(IV)%pRo=pRoMod    !1
          Ckt(II)%Tube(III)%Seg(IV)%hRo=hRoMod
          Ckt(II)%Tube(III)%Seg(IV)%tAo=tAoMod
          Ckt(II)%Tube(III)%Seg(IV)%rhAo=rhAoMod
          Ckt(II)%Tube(III)%Seg(IV)%wbAo=wbAoMod

          Ckt(II)%Tube(III)%Seg(IV)%hci=0
          Ckt(II)%Tube(III)%Seg(IV)%EFref=0
          Ckt(II)%Tube(III)%Seg(IV)%hco=0

          Ckt(II)%Tube(III)%Seg(IV)%ReVap=0
          Ckt(II)%Tube(III)%Seg(IV)%ReLiq=0

          Ckt(II)%Tube(III)%Seg(IV)%cAir=0
          Ckt(II)%Tube(III)%Seg(IV)%Rair=0
          Ckt(II)%Tube(III)%Seg(IV)%Rtube=0

          !Surface temperature
          Ckt(II)%Tube(III)%Seg(IV)%tSi=0
          Ckt(II)%Tube(III)%Seg(IV)%tSo=0

          RETURN
        END IF
      END IF

      CALL CalcSegmentRefInletConditions(II,II,III,IV,CoilType)

      CALL CalcSegmentAirInletConditions(II,II,III,IV,CoilType)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      !Defining module values
      mAiMod=Ckt(II)%Tube(III)%Seg(IV)%mAi
      tAiMod=Ckt(II)%Tube(III)%Seg(IV)%tAi
      rhAiMod=Ckt(II)%Tube(III)%Seg(IV)%rhAi

      WetFlag=0
      IF (IsSimpleCoil .EQ. 1) THEN
        RowNum=0
        IF (FinType .EQ. 4) FinType=3 !Use regular louver fin correlation, ISI - 02/12/08
      ELSE
        RowNum=Ckt(II)%Tube(III)%RowNum
      END IF

      !RS: Replace: Copying this up here so that the enthalpy will be correct for the AirSideCalc CALL
      AirPropOpt=2
      AirProp%APTDB=Ckt(II)%Tube(III)%Seg(IV)%tAi    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=Ckt(II)%Tube(III)%Seg(IV)%rhAi   !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      hAiMod=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

      IF (RowNum .EQ. 0) THEN
        CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiCoil,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
        IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair, &
        hAiMod)

        !CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
        !IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC, &
        !Poly1HTC,Poly2HTC,Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP, &
        !Poly1DP,Poly2DP,Poly3DP,Poly4DP,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair)

      ELSE
        CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiMod,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
        IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair, &
        hAiMod)
      END IF
      !Surface areas
      AoMod=AoCoil*LmodTube/Lcoil
      AfMod=AfCoil*LmodTube/Lcoil
      AiMod=AiCoil*LmodTube/Lcoil
      AmMod=AmCoil*LmodTube/Lcoil

      hco=hco*hcoMultiplier
      DPair=DPair*DPairMultiplier

      hcoMod=Ckt(II)%Tube(III)%Seg(IV)%VelDev*hco !*LmodTube/Lcoil

      Ckt(II)%Tube(III)%Seg(IV)%hco=hcoMod

      AirPropOpt=2
      AirProp%APTDB=Ckt(II)%Tube(III)%Seg(IV)%tAi    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=Ckt(II)%Tube(III)%Seg(IV)%rhAi   !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      hAiMod=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

      mRefMod=Ckt(II)%mRef
      pRiMod=Ckt(II)%Tube(III)%Seg(IV)%pRi
      hRiMod=Ckt(II)%Tube(III)%Seg(IV)%hRi

      CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
      vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
      muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      IF (IsSimpleCoil .EQ. 1) THEN
        IF (IV .EQ. 2) THEN
          hRiMod=hgRiMod*0.999 !Perturb a little to make sure it is in single phase region
        ELSEIF (IV .EQ. 3) THEN
          hRiMod=hfRiMod*0.999 !Perturb a little to make sure it is in single phase region
        END IF

        IF (IV .EQ. 2 .OR. IV .EQ. 3) THEN
          CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
          vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
          muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
          IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
          END IF
        END IF
      END IF

      tAoMod=Ckt(II)%Tube(III)%Seg(IV)%tAo
      CALL CalcMeanProp(tAiMod,tAoMod,tAmod)

      CALL CalcSegmentOutletConditions(II,III,IV,CoilType) !(II,II,III,IV,CoilType)   !RS: Debugging: Extraneous I or II
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      QmodPrev=Qmod
      Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiMod !ISI - 12/05/06
      Ckt(II)%Tube(III)%Seg(IV)%Len=LmodTube
      Ckt(II)%Tube(III)%Seg(IV)%Qmod=Qmod
      Ckt(II)%Tube(III)%Seg(IV)%pRo=pRoMod    !2
      Ckt(II)%Tube(III)%Seg(IV)%hRo=hRoMod
      Ckt(II)%Tube(III)%Seg(IV)%tAo=tAoMod
      Ckt(II)%Tube(III)%Seg(IV)%rhAo=rhAoMod
      Ckt(II)%Tube(III)%Seg(IV)%wbAo=wbAoMod

      Ckt(II)%Tube(III)%Seg(IV)%hci=hciMod
      Ckt(II)%Tube(III)%Seg(IV)%EFref=EFref
      Ckt(II)%Tube(III)%Seg(IV)%hco=hcoMod

      IF (xRmod .GE. 1) THEN
        Ckt(II)%Tube(III)%Seg(IV)%ReVap=ReVap
        Ckt(II)%Tube(III)%Seg(IV)%ReLiq=0
      ELSE IF (xRmod .LE. 0) THEN
        Ckt(II)%Tube(III)%Seg(IV)%ReVap=0
        Ckt(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
      ELSE
        Ckt(II)%Tube(III)%Seg(IV)%ReVap=ReVap
        Ckt(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
      END IF

      Ckt(II)%Tube(III)%Seg(IV)%cAir=cAir
      Ckt(II)%Tube(III)%Seg(IV)%Rair=Rair
      Ckt(II)%Tube(III)%Seg(IV)%Rtube=Rtube

      !Surface temperature
      Ckt(II)%Tube(III)%Seg(IV)%tSi=tAiMod+ABS(Qmod)*Rair
      Ckt(II)%Tube(III)%Seg(IV)%tSo=tAoMod+ABS(Qmod)*Rair

      IF (IsSimpleCoil .NE. 1) THEN
        CALL UpdateTubeDataFromCircuitData(II,III)
      END IF

    ELSE    !For the case of the microchannel condenser

      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len=LmodTube

      IF (IsSimpleCoil .EQ. 1) THEN
        IsSimpleCoil=2  !RS Comment: Prevents tripping of the following calls due to the simple case not being applicable for the microchannel case
        CALL CalcSegmentRefInletConditions(I,II,III,IV,CoilType)
        CALL CalcSegmentAirInletConditions(I,II,III,IV,CoilType)
        IsSimpleCoil=1
      ELSE
        CALL CalcSegmentRefInletConditions(I,II,III,IV,CoilType)
        CALL CalcSegmentAirInletConditions(I,II,III,IV,CoilType)
      END IF

      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      mAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%mAi !*Slab(I)%Pass(II)%Ntube
      tAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi
      rhAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi

      !RS: Replace: Copying this up here so that the enthalpy will be correct for the AirSideCalc CALL
      AirPropOpt=2
      AirProp%APTDB=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi   !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi  !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      hAiMod=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

      WetFlag=0
      RowNum=0 !Ckt(I)%Tube(J)%RowNum
      CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiMod,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
      IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair, &
      hAiMod)

      !CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
      !IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC, &
      !Poly1HTC,Poly2HTC,Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP, &
      !Poly1DP,Poly2DP,Poly3DP,Poly4DP,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair)

      hco=hco*hcoMultiplier
      DPair=DPair*DPairMultiplier

      hcoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%VelDev*hco

      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco=hcoMod

      AirPropOpt=2
      AirProp%APTDB=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi   !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi  !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      hAiMod=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

      pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
      hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi

      CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
      vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
      muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      tAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo
      CALL CalcMeanProp(tAiMod,tAoMod,tAmod)

      CALL CalcSegmentOutletConditions(II,III,IV,CoilType) !(I,II,III,IV,CoilType)    !RS: Debugging: Extraneous I
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

      QmodPrev=Qmod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod=Qmod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo=pRoMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo=hRoMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo=tAoMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAo=rhAoMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%wbAo=wbAoMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir=cAir
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=Rair
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=Rtube
      IF (xRmod .GE. 1) THEN
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReVap=ReVap
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReLiq=0
      ELSE IF (xRmod .LE. 0) THEN
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReVap=0
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
      ELSE
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReVap=ReVap
        Slab(I)%Pass(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
      END IF
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hci=hciMod
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%EFref=EFref
      Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco=hcoMod

    END IF

    RETURN

  END SUBROUTINE CalcCoilSegment

  !************************************************************************

  SUBROUTINE CalcSegmentAirInletConditions(I,II,III,IV,CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To calculate inlet air temp and relative humidity
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: I   !Slab number
    INTEGER,INTENT(IN) :: II  !Circuit,pass number
    INTEGER,INTENT(IN) :: III !Tube number
    INTEGER,INTENT(IN) :: IV  !Segment number

    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel evaporator

    REAL tAiFavg   !Average front tube inlet air temp. C
    REAL tAoFavg   !Average front tube outlet air temp. C
    REAL rhAiFavg  !Average front tube inlet RH
    REAL rhAoFavg  !Average front tube outlet RH
    REAL tAiFup    !Upper front tube inlet air temp. C
    REAL tAiFdown  !Lower front tube inlet air temp. C
    REAL rhAiFup   !Upper front tube inlet air humidity
    REAL rhAiFdown !Lower front tube inlet air humidity
    REAL tAoFup    !Upper front tube outlet air temp. C
    REAL tAoFdown  !Lower front tube outlet air temp. C
    REAL rhAoFup   !Upper front tube outlet air humidity
    REAL rhAoFdown !Lower front tube outlet air humidity
    REAL mAiFup    !Upper front tube inlet air mass flow rate, kg/s
    REAL mAiFdown  !Lower front tube inlet air mass flow rate, kg/s
    REAL VelDevFup    !Upper front tube Velocity deviation
    REAL VelDevFdown  !Lower front tube Velocity deviation

    !FLOW:

    IF (IsSimpleCoil .EQ. 1) THEN
      Ckt(II)%Tube(III)%Seg(IV)%tAi=tAiCoil
      Ckt(II)%Tube(III)%Seg(IV)%rhAi=rhAiCoil
      Ckt(II)%Tube(III)%Seg(IV)%VelDev=1
      Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiCoil*Lmodtube/Lcoil !ISI - 12/05/06
      mAiMod=Ckt(II)%Tube(III)%Seg(IV)%mAi
      RETURN
    END IF

    IF (CoilType .NE. MCCONDENSER) THEN
      IF (Tube(TubeNum)%Fup .NE. 0) THEN !Upper front tubes
        IF (Tube(TubeNum)%Even .EQ. 0) THEN !Odd tubes
          IF (Tube(Tube(TubeNum)%Fup)%Empty .NE. 0) THEN !Empty tubes

            IF (Tube(TubeNum)%Fdown .EQ. 0) THEN !Coil bottom, ISI - 07/29/07

              IF (Tube(Tube(TubeNum)%Fup)%Fdown .EQ. 0) THEN
                tAoFdown=tAiCoil
                rhAoFdown=rhAiCoil
                tAiFdown=tAiCoil
                rhAiFdown=rhAiCoil
                mAiFdown=Ckt(II)%Tube(III)%Seg(IV)%mAi
                VelDevFdown=Ckt(II)%Tube(III)%Seg(IV)%VelDev
              ELSE
                tAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAo
                rhAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAo
                tAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAi
                rhAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAi
                mAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%mAi
                VelDevFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%VelDev !Velocity deviation
              END IF

            ELSE IF (Tube(Tube(TubeNum)%Fdown)%Fup .NE. 0) THEN
              tAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAo   !Temperature
              rhAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAo !Relative humidity
              !wbAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(K)%wbAo !Wet bulb temp.
              tAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAi   !Temperature
              rhAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAi !Relative humidity
              !wbAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(K)%wbAi !Wet bulb temp.
              mAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%mAi   !mass flow rate
              VelDevFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%VelDev !Velocity deviation
            ELSE !Frontal tubes
              tAoFup=tAiCoil
              rhAoFup=rhAiCoil
              !wbAoFup=wbAiCoil
              tAiFup=tAiCoil
              rhAiFup=rhAiCoil
              !wbAiFup=wbAiCoil
              mAiFup=Ckt(II)%Tube(III)%Seg(IV)%mAi
              VelDevFup=Ckt(II)%Tube(III)%Seg(IV)%VelDev
            END IF
          ELSE
            tAoFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%tAo   !Temperature
            rhAoFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%rhAo !Relative humidity
            !wbAoFup=Tube(Tube(TubeNum)%Fup)%Seg(K)%wbAo !wet bulb temp.
            tAiFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%tAi   !Temperature
            rhAiFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%rhAi !Relative humidity
            !wbAiFup=Tube(Tube(TubeNum)%Fup)%Seg(K)%wbAi !wet bulb temp.
            mAiFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%mAi   !mass flow rate
            VelDevFup=Tube(Tube(TubeNum)%Fup)%Seg(IV)%VelDev !Velocity deviation
          END IF
        ELSE !Even tubes
          IF (Tube(Tube(TubeNum)%Fup)%Empty .NE. 0) THEN !Empty tubes

            IF (Tube(TubeNum)%Fdown .EQ. 0) THEN !Coil bottom, ISI - 07/29/07

              IF (Tube(Tube(TubeNum)%Fup)%Fdown .EQ. 0) THEN
                tAoFdown=tAiCoil
                rhAoFdown=rhAiCoil
                tAiFdown=tAiCoil
                rhAiFdown=rhAiCoil
                mAiFdown=Ckt(II)%Tube(III)%Seg(IV)%mAi
                VelDevFdown=Ckt(II)%Tube(III)%Seg(IV)%VelDev
              ELSE
                tAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAo
                rhAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAo
                tAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAi
                rhAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAi
                mAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%mAi
                VelDevFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%VelDev !Velocity deviation
              END IF

            ELSE IF (Tube(Tube(TubeNum)%Fdown)%Fup .NE. 0) THEN
              tAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%tAo   !Temperature
              rhAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%rhAo !Relative humidity
              !wbAoFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-K)%wbAo !wet bulb temp.
              tAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%tAi   !Temperature
              rhAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%rhAi !Relative humidity
              !wbAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-K)%wbAi !wet bulb temp.
              mAiFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%mAi   !mass flow rate
              VelDevFup=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(NumOfMods+1-IV)%VelDev !Velocity deviation
            ELSE !Frontal tubes
              tAoFup=tAiCoil
              rhAoFup=rhAiCoil
              !wbAoFup=wbAiCoil
              tAiFup=tAiCoil
              rhAiFup=rhAiCoil
              !wbAiFup=wbAiCoil
              mAiFup=Ckt(II)%Tube(III)%Seg(IV)%mAi
              VelDevFup=Ckt(II)%Tube(III)%Seg(IV)%VelDev
            END IF
          ELSE
            tAoFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%tAo   !Temperature
            rhAoFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%rhAo !Relative humidity
            !wbAoFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-K)%wbAo !wet bulb temp.
            tAiFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%tAi   !Temperature
            rhAiFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%rhAi !Relative humidity
            !wbAiFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-K)%wbAi !Relative humidity
            mAiFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%mAi   !mass flow rate
            VelDevFup=Tube(Tube(TubeNum)%Fup)%Seg(NumOfMods+1-IV)%VelDev !Velocity deviation
          END IF
        END IF
      ELSE !Front row tubes
        tAoFup=tAiCoil
        rhAoFup=rhAiCoil
        !wbAoFup=wbAiCoil
        tAiFup=tAiCoil
        rhAiFup=rhAiCoil
        !wbAiFup=wbAiCoil
        mAiFup=Ckt(II)%Tube(III)%Seg(IV)%mAi
        VelDevFup=Ckt(II)%Tube(III)%Seg(IV)%VelDev
      END IF

      IF (Tube(TubeNum)%Fdown .NE. 0) THEN !Lower front tube
        IF (Tube(TubeNum)%Even .EQ. 0) THEN !Odd tubes
          IF (Tube(Tube(TubeNum)%Fdown)%Empty .NE. 0) THEN !Empty tubes

            IF (Tube(TubeNum)%Fup .EQ. 0) THEN !Coil top, ISI - 07/29/07

              tAoFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAo
              rhAoFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAo
              tAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAi
              rhAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAi
              mAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%mAi
              VelDevFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%VelDev !Velocity deviation

            ELSE IF (Tube(Tube(TubeNum)%Fup)%Fdown .NE. 0) THEN
              tAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAo
              rhAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAo
              !wbAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(K)%wbAo
              tAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAi
              rhAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAi
              !wbAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(K)%wbAi
              mAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%mAi
              VelDevFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%VelDev !Velocity deviation
            ELSE !Frontal tubes
              tAoFdown=tAiCoil
              rhAoFdown=rhAiCoil
              !wbAoFdown=wbAiCoil
              tAiFdown=tAiCoil
              rhAiFdown=rhAiCoil
              !wbAiFdown=wbAiCoil
              mAiFdown=Ckt(II)%Tube(III)%Seg(IV)%mAi
              VelDevFdown=Ckt(II)%Tube(III)%Seg(IV)%VelDev
            END IF
          ELSE
            tAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%tAo
            rhAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%rhAo
            !wbAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(K)%wbAo
            tAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%tAi
            rhAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%rhAi
            !wbAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(K)%wbAi
            mAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%mAi
            VelDevFdown=Tube(Tube(TubeNum)%Fdown)%Seg(IV)%VelDev !Velocity deviation
          END IF
        ELSE !Even tubes
          IF (Tube(Tube(TubeNum)%Fdown)%Empty .NE. 0) THEN !Empty tubes

            IF (Tube(TubeNum)%Fup .EQ. 0) THEN !Coil top, ISI - 07/29/07

              tAoFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAo
              rhAoFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAo
              tAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%tAi
              rhAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%rhAi
              mAiFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%mAi
              VelDevFdown=Tube(Tube(Tube(TubeNum)%Fdown)%Fup)%Seg(IV)%VelDev !Velocity deviation

            ELSE IF (Tube(Tube(TubeNum)%Fup)%Fdown .NE. 0) THEN
              tAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%tAo
              rhAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%rhAo
              !wbAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-K)%wbAo
              tAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%tAi
              rhAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%rhAi
              !wbAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-K)%wbAi
              mAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%mAi
              VelDevFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(NumOfMods+1-IV)%VelDev !Velocity deviation
            ELSE !Frontal tubes
              tAoFdown=tAiCoil
              rhAoFdown=rhAiCoil
              !wbAoFdown=wbAiCoil
              tAiFdown=tAiCoil
              rhAiFdown=rhAiCoil
              !wbAiFdown=wbAiCoil
              mAiFdown=Ckt(II)%Tube(III)%Seg(IV)%mAi
              VelDevFdown=Ckt(II)%Tube(III)%Seg(IV)%VelDev
            END IF
          ELSE
            tAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%tAo
            rhAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%rhAo
            !wbAoFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-K)%wbAo
            tAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%tAi
            rhAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%rhAi
            !wbAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-K)%wbAi
            mAiFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%mAi
            VelDevFdown=Tube(Tube(TubeNum)%Fdown)%Seg(NumOfMods+1-IV)%VelDev !Velocity deviation
          END IF
        END IF
      ELSE !Front row tubes
        tAoFdown=tAiCoil
        rhAoFdown=rhAiCoil
        !wbAoFdown=wbAiCoil
        tAiFdown=tAiCoil
        rhAiFdown=rhAiCoil
        !wbAiFdown=wbAiCoil
        mAiFdown=Ckt(II)%Tube(III)%Seg(IV)%mAi
        VelDevFdown=Ckt(II)%Tube(III)%Seg(IV)%VelDev
      END IF

      !Average values
      tAiFavg=(tAiFup+tAiFdown)/2
      tAoFavg=(tAoFup+tAoFdown)/2

      rhAiFavg=(rhAiFup+rhAiFdown)/2
      rhAoFavg=(rhAoFup+rhAoFdown)/2

      Ckt(II)%Tube(III)%Seg(IV)%tAi=tAiFavg-1*(tAiFavg-tAoFavg)
      Ckt(II)%Tube(III)%Seg(IV)%rhAi=rhAiFavg-1*(rhAiFavg-rhAoFavg)

      IF (Ckt(II)%Tube(III)%Seg(IV)%VelDev .LE. 0) Ckt(II)%Tube(III)%Seg(IV)%VelDev=1

      Ckt(II)%Tube(III)%Seg(IV)%mAi=(mAiFup+mAiFdown)/2
      mAiMod=Ckt(II)%Tube(III)%Seg(IV)%mAi
      Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiCoil*LmodTube/(Ltube*Nt)*Ckt(II)%Tube(III)%Seg(IV)%VelDev
      mAiMod=Ckt(II)%Tube(III)%Seg(IV)%mAi

    ELSE !Microchannel coil
      IF (I .EQ. 1) THEN !1st slab
        Slab(I)%Pass(II)%Tube(1)%Seg(IV)%tAi=tAiCoil
        Slab(I)%Pass(II)%Tube(1)%Seg(IV)%rhAi=rhAiCoil
      ELSE
        Slab(I)%Pass(II)%Tube(1)%Seg(IV)%tAi=Slab(I-1)%tAo
        Slab(I)%Pass(II)%Tube(1)%Seg(IV)%rhAi=Slab(I-1)%rhAo
      END IF
    END IF

    RETURN

  END SUBROUTINE CalcSegmentAirInletConditions

  !************************************************************************

  SUBROUTINE CalcSegmentRefInletConditions(I,II,III,IV,CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To calculate inlet refrigerant pressure and enthalpy
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !November 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: I   !Slab number
    INTEGER,INTENT(IN) :: II  !Circuit,pass number
    INTEGER,INTENT(IN) :: III !Tube number
    INTEGER,INTENT(IN) :: IV  !Segment number

    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel Evaporator

    !FLOW:

    IF (CoilType .NE. MCCONDENSER) THEN

      IF (III .EQ. 1 .AND. IV .EQ. 1) THEN !Equal to circuit inlet
        Ckt(II)%Tube(III)%Seg(IV)%pRi=Ckt(II)%pRi
        Ckt(II)%Tube(III)%Seg(IV)%hRi=Ckt(II)%hRi

      ELSE IF (K .EQ. 1) THEN !Equal to outlet of previous tube
        Ckt(II)%Tube(III)%Seg(IV)%pRi=Ckt(II)%Tube(III-1)%Seg(NumOfMods)%pRo
        Ckt(II)%Tube(III)%Seg(IV)%hRi=Ckt(II)%Tube(III-1)%Seg(NumOfMods)%hRo

      ELSE !Equal to outlet of previous module(section)
        Ckt(II)%Tube(III)%Seg(IV)%pRi=Ckt(II)%Tube(III)%Seg(IV-1)%pRo
        Ckt(II)%Tube(III)%Seg(IV)%hRi=Ckt(II)%Tube(III)%Seg(IV-1)%hRo
      END IF

    ELSE

      IF (IsParallelSlabs .GT. 0) THEN
        !ISI - 07/13/07
        IF (II .EQ. 1) THEN !1st pass
          IF (IV .EQ. 1) THEN !1st segment
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=pRiCoil
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=hRiCoil
          ELSE
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
          END IF
        ELSE
          IF (IV .EQ. 1) THEN !1st segment
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%pRo
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%hRo
          ELSE
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
            Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
          END IF
        END IF
      ELSE !Series

        IF (I .EQ. 1) THEN !1st slab
          IF (II .EQ. 1) THEN !1st pass
            IF (IV .EQ. 1) THEN !1st segment
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=pRiCoil
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=hRiCoil
            ELSE
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
            END IF
          ELSE
            IF (IV .EQ. 1) THEN !1st segment
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%hRo
            ELSE
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
            END IF
          END IF
        ELSE
          IF (II .EQ. 1) THEN !1st pass
            IF (IV .EQ. 1) THEN !1st segment
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I-1)%Pass(Slab(I-1)%Npass)%Tube(1)%Seg(NumOfMods)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I-1)%Pass(Slab(I-1)%Npass)%Tube(1)%Seg(NumOfMods)%hRo
            ELSE
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
            END IF
          ELSE
            IF (IV .EQ. 1) THEN !1st segment
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II-1)%Tube(1)%Seg(NumOfMods)%hRo
            ELSE
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%pRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%pRo
              Slab(I)%Pass(II)%Tube(1)%Seg(IV)%hRi=Slab(I)%Pass(II)%Tube(1)%Seg(IV-1)%hRo
            END IF
          END IF
        END IF

      END IF

    END IF

    RETURN

  END SUBROUTINE CalcSegmentRefInletConditions

  !************************************************************************

  SUBROUTINE CalcSegmentOutletConditions(II,III,IV,CoilType) !(I,II,III,IV,CoilType)  !RS: Debugging: Extraneous I

    !------------------------------------------------------------------------
    !Purpose:
    !To calculate segment outlet conditions
    !
    !Author
    !Ipseng Iu
    !Oklahoma State Univerity, Stillwater
    !
    !Date
    !March 2005
    !
    !Reference:
    !none
    !
    !------------------------------------------------------------------------

    USE FluidProperties_HPSim
    USE CoilCalcMod
    USE AirPropMod
    USE OilMixtureMod
    !USE DataGlobals_HPSim   !RS: Debugging: Cavallini (2/14/14)

    IMPLICIT NONE

    !INTEGER,INTENT(IN) :: I   !Slab number !RS: Debugging: Extraneous
    INTEGER,INTENT(IN) :: II  !Circuit,pass number
    INTEGER,INTENT(IN) :: III !Tube number
    INTEGER,INTENT(IN) :: IV  !Segment number

    INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel evaporator

    REAL Rtot         !Total resistance, K-m^2/W
    REAL Qsolar       !Solar radiation, kW
    REAL DPreturnbend !Pressure drop at return bend, kPa
    REAL DiffpRoMod   !Difference in pRoMod
    REAL DiffhRoMod   !Difference in hRoMod
    REAL PrevpRoMod   !Previous value of pRoMod
    REAL PrevhRoMod   !Previous value of hRoMod
    INTEGER RefBCiter             !Iteration loop counter
    LOGICAL IsTransitionSegment !Flag to indicate if it is transtion segment

    !FLOW:

    !Initialize for property iteration, to find the mean property
    hfgRoMod=0;  xRoMod=0;  vgRoMod=0;  vfRoMod=0
    muRoMod=0;  mugRoMod=0;  mufRoMod=0
    kRoMod=0;	  kfRoMod=0;  kgRoMod=0
    cpRoMod=0;  cpfRoMod=0;  cpgRoMod=0
    DTmod=0;

    PrevpRoMod=BIG
    PrevhRoMod=BIG

    IsTransitionSegment=.FALSE.

    DO RefBCiter=1, RefBCmaxIter

      !Correct quality
      IF (xRoMod .GT. 1) THEN
        xRoMod=1
      ELSEIF (xRoMod .LT. 0) THEN
        xRoMod=0
      ENDIF
      IF (xRiMod .GT. 1) THEN
        xRiMod=1
      ELSEIF (xRiMod .LT. 0) THEN
        xRiMod=0
      ENDIF

      !Calculate mean properties
      CALL CalcMeanProp(hfgRiMod,hfgRoMod,hfgRmod)
      CALL CalcMeanProp(xRiMod,xRoMod,xRmod)
      CALL CalcMeanProp(vgRiMod,vgRoMod,vgRmod)
      CALL CalcMeanProp(vfRiMod,vfRoMod,vfRmod)
      CALL CalcMeanProp(muRiMod,muRoMod,muRmod)
      CALL CalcMeanProp(mugRiMod,mugRoMod,mugRmod)
      CALL CalcMeanProp(mufRiMod,mufRoMod,mufRmod)
      CALL CalcMeanProp(kRiMod,kRoMod,kRmod)
      CALL CalcMeanProp(kfRiMod,kfRoMod,kfRmod)
      CALL CalcMeanProp(kgRiMod,kgRoMod,kgRmod)
      CALL CalcMeanProp(cpRiMod,cpRoMod,cpRmod)
      CALL CalcMeanProp(cpfRiMod,cpfRoMod,cpfRmod)
      CALL CalcMeanProp(cpgRiMod,cpgRoMod,cpgRmod)

      !Correct specific heat
      IF (cpRmod .LE. 0) THEN !ISI - 08/03/06
        IF (xRmod .LE. 0) THEN
          cpRmod = cpfRmod
        END IF
        IF (xRmod .GE. 1) THEN
          cpRmod = cpgRmod
        END IF
      END IF

      !Correct thermal conductivity
      IF (kRmod .LE. 0) THEN !ISI - 08/03/06
        IF (xRmod .LE. 0) THEN
          kRmod = kfRmod
        END IF
        IF (xRmod .GE. 1) THEN
          kRmod = kgRmod
        END IF
      END IF

      IF (muRmod .LE. 0) THEN !ISI - 08/03/06
        IF (xRmod .LE. 0) THEN
          muRmod = mufRmod
        END IF
        IF (xRmod .GE. 1) THEN
          muRmod = mugRmod
        END IF
      END IF

      LmodTPratio=0
      QmodTP=0
      LmodSHratio=0
      QmodSH=0

      !For segment covers both two phase and single phase region
      IF (RefBCiter .GT. 1 .AND. &
      ((xRiMod .GT. 0 .AND. xRiMod .LT. 1 .AND. xRoMod .LE. 0) .OR. & !Condenser outlet
      (xRiMod .GE. 1 .AND. xRoMod .LT. 1 .AND. xRoMod .GT. 0))) THEN  !Condenser inlet

      CALL CalcTransitionSegment(CoilType)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF
      IF (IsSimpleCoil .EQ. 1) THEN
        IsTransitionSegment=.TRUE.
      END IF

      !Update properties ISI - 08/03/06
      IF (cpRmod .LE. 0) THEN
        IF (xRmod .LE. 0) THEN
          cpRmod = cpfRmod
        END IF
        IF (xRmod .GE. 1) THEN
          cpRmod = cpgRmod
        END IF
      END IF

      IF (kRmod .LE. 0) THEN !ISI - 08/03/06
        IF (xRmod .LE. 0) THEN
          kRmod = kfRmod
        END IF
        IF (xRmod .GE. 1) THEN
          kRmod = kgRmod
        END IF
      END IF

      IF (muRmod .LE. 0) THEN !ISI - 08/03/06
        IF (xRmod .LE. 0) THEN
          muRmod = mufRmod
        END IF
        IF (xRmod .GE. 1) THEN
          muRmod = mugRmod
        END IF
      END IF

    END IF

    !Condenser inlet
    IF (DTmod .EQ. 0) THEN
      DTmod=(tAiMod+tRiMod)/2 !First estimate
    END IF
    !CALL hcRefside(CoilType,TubeType,IDtube,ktube,mRefMod,Qmod,AoMod,AiMod,hfgRmod, &               !Calculating the refrigerant side heat transfer coefficient
    !xRmod,xRmod,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
    !MolWeight,Psat,Pcr,Tsat,SigmaMod,DTmod,Wabsolute,EFref,hciMod)
    CALL hcRefside(CoilType,TubeType,IDtube,mRefMod,Qmod, &     !Calculating the refrigerant side heat transfer coefficient
    xRmod,xRmod,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
    Psat,Pcr,Wabsolute,EFref,hciMod)

    hciMod=hciMod*hciMultiplier

    CALL Reynolds(IDtube,mRefMod,xRmod,muRmod,mugRmod,mufRmod,ReVap,ReLiq)

    WetFlag=0
    !TsurfMod=0 !RS: Debugging: Set once never used

    !Calc. UA
    CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
    hcoMod,hciMod,AfMod,AoMod,AiMod,AmMod,UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)

    IF (xRiMod .GT. 0 .AND. xRoMod .LE. 0 .AND. LmodTPratio .LT. 1) THEN !Condenser outlet
      UA=UA*(1-LmodTPratio)
    ELSEIF (xRiMod .GE. 1 .AND. xRoMod .LT. 1 .AND. LmodSHratio .LT. 1) THEN !Condenser inlet
      UA=UA*(1-LmodSHratio)
    END IF

    !Calc. Cref
    IF (CoilType .NE. MCCONDENSER) THEN
      cRef=mRefMod*cpRmod
    ELSE
      cRef=mRefMod*cpRmod*NumOfChannels
    END IF
    IF (xRmod .LT. 1. .AND. xRmod .GT. 0.) THEN
      cRef=BIG !Phase change
    END IF

    !Calc. Cair
    !CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    cAir=mAiMod*cpAir

    !Calc. Cmin
    Cmin=MIN(cAir,cRef)

    !Calc. Epsilon
    CALL EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)

    !Calc. DT
    IF (LmodTPratio .GT. 0 .OR. LmodSHratio .GT. 0) THEN !ISI - 07/21/06
      DT=(tRmod-tAiMod)
    ELSE
      DT=(tRiMod-tAiMod)
    END IF

    !Calc. Q module
    Qmod=EPS*Cmin*DT

    !Condenser outlet
    IF (xRiMod .GT. 0 .AND. xRoMod .LE. 0) THEN
      IF (IsSimpleCoil .EQ. 1) THEN
        IF (QmodTP .NE. 0) THEN
          Qmod = QmodTP
        END IF
      ELSE
        IF (LmodTP .EQ. LmodTube) THEN
          IF (Qmod .GT. QmodTP) THEN
            Qmod = QmodTP
          END IF
        ELSE
          Qmod=Qmod+QmodTP
        END IF
      END IF
    END IF

    !Condenser inlet
    IF (xRiMod .GE. 1 .AND. xRoMod .LT. 1) THEN
      IF (IsSimpleCoil .EQ. 1) THEN
        IF (QmodSH .NE. 0) THEN
          Qmod = QmodSH
        END IF
      ELSE
        IF (LmodSH .EQ. LmodTube) THEN
          IF (Qmod .GT. QmodSH) THEN
            Qmod = QmodSH
          END IF
        ELSE
          Qmod=Qmod+QmodSH
        END IF
      END IF
    END IF

    IF (CoilType .NE. MCCONDENSER) THEN
      !Include solar radiation
      IF (IsCoolingMode .GT. 1 .AND. &
      Ckt(II)%Tube(III)%Fup .EQ. 0 .AND. Ckt(II)%Tube(III)%Fdown .EQ. 0) THEN
      Rtot=Rair*AoMod+Rrefrig*AiMod+Rtube*AmMod
      Qsolar=Rair*AoMod/Rtot*Ckt(II)%Tube(III)%Seg(IV)%Aface*SurfAbsorptivity*SolarFlux
      Qmod=Qmod-Qsolar
    END IF
  END IF

  !Calc. Outside air enthalpy
  IF (CoilType .NE. MCCONDENSER) THEN
    hRoMod=-Qmod/mRefMod+hRiMod
  ELSE
    hRoMod=-(Qmod/NumOfChannels)/mRefMod+hRiMod
  END IF

  CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
  vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
  muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)

  CALL CalcSegmentRefOutletPressure(CoilType,TubeType,pRiMod,hgRiMod,hfRiMod, & !CoilType,TubeType,tRiMod,pRiMod,hgRiMod,hfRiMod, &
  hRiMod,hRoMod,xRiMod,vRiMod,vgRiMod,vfRiMod,mRefMod, &
  muRiMod,mugRiMod,mufRiMod,LModTube,LmodTPratio, & !muRiMod,mugRiMod,mufRiMod,SigmaMod,LmodTube,LmodTPratio, &    !RS: Debugging: Extraneous SigmaMod
  Dchannel,HtCoil,Lcoil,DPrefMultiplier,pRoMod)

  IF (ErrorFlag .GT. CONVERGEERROR) THEN
    RETURN
  END IF

  CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
  vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
  muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
  IF (ErrorFlag .GT. CONVERGEERROR) THEN
    RETURN
  END IF

  IF (CoilType .NE. MCCONDENSER .AND. IsSimpleCoil .NE. 1) THEN
    !Return bend pressure drop
    IF (K .EQ. NumOfMods) THEN
      IF (J .EQ. LastTube) THEN
        IF (Ckt(II)%OutSplit .GT. 1 .OR. Ckt(II)%OutJoin .GT. 1) THEN
          CALL returnbend(CoilType,TubeType,IDtube,Pt,mRefmod,xRoMod,vRoMod,vgRoMod,vfRoMod,mugRoMod,mufRoMod,DPreturnbend)
          !CoilType,TubeType,ID,Pt,mRef,xRef,vRef,vgi,vfi,mug,muf,dPret)
          pRoMod=pRoMod-DPreturnbend
        END IF
      ELSE
        CALL returnbend(CoilType,TubeType,IDtube,Pt,mRefmod,xRoMod,vRoMod,vgRoMod,vfRoMod,mugRoMod,mufRoMod,DPreturnbend)
        pRoMod=pRoMod-DPreturnbend
      END IF

      CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
      vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
      muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
      IF (ErrorFlag .GT. CONVERGEERROR) THEN
        RETURN
      END IF

    END IF
  END IF

  IF (IsSimpleCoil .EQ. 1) THEN
    IF (IsTransitionSegment) THEN
      EXIT
    END IF
  END IF

  !Correct the equation, Sankar 2/19/2009 - 9:30pm
  DTmod=Qmod*(1/(hciMod*AiMod)+LOG(ODtube/IDtube)/(2*PI*Ktube*LmodTube))
  DiffpRoMod=ABS((pRoMod-PrevpRoMod)/PrevpRoMod)
  DiffhRoMod=ABS((hRoMod-PrevhRoMod)/PrevhRoMod)
  IF (DiffpRoMod .GT. SMALL .OR. DiffhRoMod .GT. SMALL) THEN
    PrevpRoMod=pRoMod
    PrevhRoMod=hRoMod
  ELSE
    EXIT
  END IF

END DO !end of RefBCiter

IF (RefBCiter .GT. RefBCmaxIter) THEN
  ErrorFlag=CONVERGEERROR
END IF

!Outside air temp
tAoMod=Qmod/cAir+tAiMod

!Calc. Outside air enthalpy
hAoMod=Qmod/mAiMod+hAiMod

AirPropOpt=1
AirProp%APTDB=tAoMod   !RS: Debugging: Formerly AirProp(1)
AirProp%APEnth=hAoMod   !RS: Debugging: Formerly AirProp(4)
CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
rhAoMod=AirProp%APRelHum  !RS: Debugging: Formerly AirProp(3)
wbAoMod=AirProp%APTWB  !RS: Debugging: Formerly AirProp(5)

!SigmaTest=SigmaMod  !RS: Debugging: Cavallini (2/14/14)
!hfgTest=hfgRoMod    !RS: Debugging: Cavallini (2/14/14)
!DTTest=DT           !RS: Debugging: Cavallini (2/14/14)

RETURN

END SUBROUTINE CalcSegmentOutletConditions

!************************************************************************

SUBROUTINE CalcTransitionSegment(CoilType) !,NumOfChannels)

  !------------------------------------------------------------------------
  !Purpose:
  !To calculate transition segment (both single and two phase refrigerant
  !in segment) heat transfer
  !
  !Author
  !Ipseng Iu
  !Oklahoma State University, Stillwater
  !
  !Date
  !March 2005
  !
  !Reference:
  !none
  !
  !------------------------------------------------------------------------

  USE FluidProperties_HPSim
  USE CoilCalcMod
  USE AirPropMod
  USE OilMixtureMod

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
  !3=High side interconnecting pipes;
  !4=Low side interconnecting pipes
  !5=Microchannel condenser
  !6=Microchannel evaporator
  !FLOW:

  !Condenser inlet
  IF (xRiMod .GE. 1 .AND. xRoMod .LT. 1)  THEN
    xRmod=xRiMod
  END IF

  hciMod=hciMod*hciMultiplier

  !Condenser outlet
  IF (xRiMod .GT. 0 .AND. xRoMod .LE. 0)  THEN
    IF (CoilType .NE. MCCONDENSER) THEN
      QmodTP=mRefMod*(hRiMod-hfRiMod)
    ELSE
      QmodTP=mRefMod*NumOfChannels*(hRiMod-hfRiMod)
    END IF
    cRef=BIG !Phase change
  END IF

  !Condenser inlet
  IF (xRiMod .GE. 1 .AND. xRoMod .LT. 1)  THEN
    IF (CoilType .NE. MCCONDENSER) THEN
      QmodSH=mRefMod*(hRiMod-hgRiMod)
      cRef=mRefMod*cpRmod
    ELSE
      QmodSH=mRefMod*NumOfChannels*(hRiMod-hgRiMod)
      cRef=mRefMod*NumOfChannels*cpRmod
    END IF

  END IF

  !Calc. DT
  DT=(tRiMod-tAiMod)

  !Find transition boundary
  CALL FindTransitionBoundary(CoilType)

  !Condenser outlet
  IF (xRiMod .GT. 0 .AND. xRoMod .LE. 0)  THEN
    IF (LmodTP .NE. LmodTube) THEN
      xRmod=0
    ELSE
      IF (CoilType .NE. MCCONDENSER) THEN
        hRoMod=-QmodTP/mRefMod+hRiMod !ISI - 06/18/05
      ELSE
        hRoMod=-(QmodTP/NumOfChannels)/mRefMod+hRiMod
      END IF
      Pressure=pRiMod*1000    !RS Comment: Unit Conversion
      Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion
      xRmod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)    !RS Comment: Refrigerant Quality for a Module
      IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 1328'
        ErrorFlag=REFPROPERROR
        RETURN
      END IF
      !Keeping the module qualities in decimal form
      IF (xRmod .GT. 1) THEN
        xRmod=1
      END IF
      IF (xRmod .LT. 0) THEN
        xRmod=0
      END IF
    END IF

    !ISI - 07/21/06 to update the refrigerant temperature at the transition boundary
    Pressure=pRiMod*1000    !RS Comment: Unit Conversion
    Quality=xRmod
    tRmod=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !RS Comment: Refrigerant Temperature for a Module

  END IF

  !Condenser inlet
  IF (xRiMod .GE. 1 .AND. xRoMod .LT. 1)  THEN
    IF (LmodSH .NE. LmodTube) THEN
      xRmod=(xRoMod+1)/2

      !ISI - 07/21/06 to update the refrigerant temperature at the transition boundary
      Pressure=pRiMod*1000    !RS Comment: Unit Conversion
      Quality=xRmod
      tRmod=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !RS Comment: Refrigerant Temperature for a Module
    END IF
  END IF

  RETURN

END SUBROUTINE CalcTransitionSegment

!************************************************************************

SUBROUTINE FindTransitionBoundary(CoilType)

  !------------------------------------------------------------------------
  !Purpose:
  !To find the transition boundary in transition segment
  !
  !Author
  !Ipseng Iu
  !Oklahoma State University, Stillwater
  !
  !Date
  !March 2005
  !
  !Reference:
  !none
  !
  !------------------------------------------------------------------------

  USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  USE CoilCalcMod
  USE AirPropMod

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
  !3=High side interconnecting pipes;
  !4=Low side interconnecting pipes
  !5=Microchannel condenser
  !6=Microchannel evaporator

  !FLOW:

  !Initialize
  LmodTPmin=0
  LmodTPmax=LmodTube
  LmodTP=LmodTube/2
  LmodSHmin=0
  LmodSHmax=LmodTube
  LmodSH=LmodTube/2
  DTmod=0
  DO

    !Condenser outlet
    IF (xRiMod .GT. 0 .AND. xRoMod .LE. 0)  THEN
      LmodTPratio=LmodTP/LmodTube
      LmodSHratio=1-LmodTPratio
      WetFlag=0
      !TsurfMod=0 !RS: Debugging: Set once but never used

      !ISI - 09/11/06
      IF (DTmod .EQ. 0) THEN
        DTmod=(tAiMod+tRiMod)/2 !First estimate
      END IF
      !CALL hcRefside(CoilType,TubeType,IDtube,ktube,mRefMod,Qmod,AoMod,AiMod,hfgRmod,xRiMod,0.0, &    !Calculate the refrigerant side heat transfer coefficient
      !vgRmod,vfRmod,muRmod,mugRmod,mufRmod, &
      !kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
      !MolWeight,Psat,Pcr,Tsat,SigmaMod,DTmod,Wabsolute,EFref,hciMod)
      CALL hcRefside(CoilType,TubeType,IDtube,mRefMod,Qmod, &               !Calculating the refrigerant side heat transfer coefficient
      xRimod,0.0,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
      Psat,Pcr,Wabsolute,EFref,hciMod)

      !Calc. UA
      CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
      hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
      UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)
      !CALL CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth,RowNum,tAiMod,hAiMod, &
      !hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
      !UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)

      IF (IsSimpleCoil .EQ. 1) THEN
        mAiMod=mAiCoil*LmodTP/Lcoil !ISI - 12/05/06
      END IF

      !Calc. Cair !ISI - 12/05/06
      !CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
      Cair=mAiMod*CPair

      !Calc. Cmin
      Cmin=MIN(cAir,cRef)

      !Calc. Epsilon
      CALL EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)

      !Calc. dry module heat transfer
      Qmod=EPS*Cmin*DT

      DTmod=Qmod*(1/(hciMod*AiMod*LmodTP/LmodTube)+LOG(ODtube/IDtube)/(2*PI*Ktube*LmodTP))
      IF (ABS(LmodTPratio-1) .LT. SMALL) THEN
        LmodTP=LmodTube !ISI - 06/18/05
        LmodTPratio=1 !ISI - 06/18/05
        QmodTP=Qmod
        EXIT
      END IF
      IF (ABS((Qmod-QmodTP)/QmodTP) .GT. SMALL .AND. ABS(LmodTPmax-LmodTPmin)/LmodTP .GT. SMALL) THEN
        IF (ABS(Qmod) .GT. ABS(QmodTP)) THEN
          LmodTPmax=LmodTP
        ELSE
          LmodTPmin=LmodTP
        END IF
        LmodTP=(LmodTPmax+LmodTPmin)/2
      ELSE
        IF (IsSimpleCoil .EQ. 1) THEN
          LmodTube=LmodTP
        END IF
        EXIT
      END IF
    END IF

    !Condenser inlet
    IF (xRiMod .GE. 1 .AND. xRoMod .LT. 1)  THEN
      LmodSHratio=LmodSH/LmodTube
      LmodTPratio=1-LmodSHratio
      WetFlag=0
      !TsurfMod=0 !RS: Debugging: Set once but never used

      !ISI - 09/11/06
      IF (DTmod .EQ. 0) THEN
        DTmod=(tAiMod+tRiMod)/2 !First estimate
      END IF
      !CALL hcRefside(CoilType,TubeType,IDtube,ktube,mRefMod,Qmod,AoMod,AiMod,hfgRmod,xRiMod,1.0, &    !Calculate the refrigerant side heat transfer coefficient
      !vgRmod,vfRmod,muRmod,mugRmod,mufRmod, &
      !kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
      !MolWeight,Psat,Pcr,Tsat,SigmaMod,DTmod,Wabsolute,EFref,hciMod)
      CALL hcRefside(CoilType,TubeType,IDtube,mRefMod,Qmod, &               !Calculating the refrigerant side heat transfer coefficient
      xRimod,1.0,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
      Psat,Pcr,Wabsolute,EFref,hciMod)

      !Calc. UA
      CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
      hcoMod,hciMod,AfMod*LmodSHratio,AoMod*LmodSHratio,AiMod*LmodSHratio,AmMod*LmodSHratio, &
      UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)
      !CALL CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth,RowNum,tAiMod,hAiMod, &
      !hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
      !UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)

      IF (IsSimpleCoil .EQ. 1) THEN
        mAiMod=mAiCoil*LmodSH/Lcoil !ISI - 12/05/06
      END IF

      !Calc. Cair !ISI - 12/05/06
      !CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
      Cair=mAiMod*CPair

      !Calc. Cmin
      Cmin=MIN(cAir,cRef)

      !Calc. Epsilon
      CALL EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)

      !Calc. dry module heat transfer
      Qmod=EPS*Cmin*DT

      DTmod=Qmod*(1/(hciMod*AiMod*LmodSH/LmodTube)+LOG(ODtube/IDtube)/(2*PI*Ktube*LmodSH))
      IF (ABS(LmodSHratio-1) .LT. SMALL) THEN
        LmodSH=LmodTube !ISI - 05/18/05
        LmodSHratio=1 !ISI - 05/18/05
        QmodSH=Qmod
        EXIT
      END IF
      IF (ABS((Qmod-QmodSH)/QmodSH) .GT. SMALL .AND. ABS(LmodSHmax-LmodSHmin)/LmodSH .GT. SMALL) THEN
        IF (ABS(Qmod) .GT. ABS(QmodSH)) THEN
          LmodSHmax=LmodSH
        ELSE
          LmodSHmin=LmodSH
        END IF
        LmodSH=(LmodSHmax+LmodSHmin)/2
      ELSE
        IF (IsSimpleCoil .EQ. 1) THEN
          LmodTube=LmodSH
        END IF
        EXIT
      END IF
    END IF
  END DO !End of LmodTPratio, LmodSHratio

  RETURN

END SUBROUTINE FindTransitionBoundary

!************************************************************************

SUBROUTINE CalcRefProperty(pRef,hRef,hfRef,hgRef,hfgRef,Psat,Tsat,tRef,xRef, &
  vRef,vfRef,vgRef,cpRef,cpfRef,cpgRef, &
  muRef,mufRef,mugRef,kRef,kfRef,kgRef,SigmaRef)

  !------------------------------------------------------------------------
  !Purpose:
  !To calculate refrigerant properties given pressure and enthalpy
  !
  !Author
  !Ipseng Iu
  !Oklahoma State University, Stillwater
  !
  !Date
  !March 2005
  !
  !Reference:
  !none
  !
  !------------------------------------------------------------------------

  USE FluidProperties_HPSim
  USE OilMixtureMod

  IMPLICIT NONE

  REAL, INTENT(IN)  :: pRef   !Pressure, kPa
  REAL, INTENT(IN)  :: hRef   !Enthalpy, kJ/kg
  REAL, INTENT(OUT) :: hfRef  !Liquid enthalpy, kJ/kg
  REAL, INTENT(OUT) :: hgRef  !Vapor enthalpy, kJ/kg
  REAL, INTENT(OUT) :: hfgRef !hg - hf
  REAL, INTENT(OUT) :: Psat   !Saturation pressure, kPa
  REAL, INTENT(OUT) :: Tsat   !Saturation temperature, C
  REAL, INTENT(OUT) :: tRef   !Temperature, C
  REAL, INTENT(OUT) :: xRef   !Quality
  REAL, INTENT(OUT) :: vRef   !Specific volume, m^3/kg
  REAL, INTENT(OUT) :: vfRef  !Liquid specific volume, m^3/kg
  REAL, INTENT(OUT) :: vgRef  !Vapor specific volume, m^3/kg
  REAL, INTENT(OUT) :: cpRef  !Specific heat, kJ/kg-K
  REAL, INTENT(OUT) :: cpfRef !Liquid specific heat, kJ/kg-K
  REAL, INTENT(OUT) :: cpgRef !Vapor specific heat, kJ/kg-K
  REAL, INTENT(OUT) :: muRef  !Dynamic viscoity, Pa-s
  REAL, INTENT(OUT) :: mufRef !Liquid dynamic viscoity, Pa-s
  REAL, INTENT(OUT) :: mugRef !Vapor dynamic viscoity, Pa-s
  REAL, INTENT(OUT) :: kRef   !Thermal conductivity, kW/m-K
  REAL, INTENT(OUT) :: kfRef  !Liquid thermal conductivity, kW/m-K
  REAL, INTENT(OUT) :: kgRef  !Vapor thermal conductivity, kW/m-K
  REAL, INTENT(OUT) :: SigmaRef  !Surface tension, N/m

  !LOCAL VARIABLES:
  REAL Wlocal !Local oil mass fraction

  !FLOW:

  Pressure=pRef*1000  !RS Comment: Unit Conversion
  Enthalpy=hRef*1000  !RS Comment: Unit Conversion

  tRef=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !RS Comment: Refrigerant Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3144'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  xRef=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !RS Comment: Refrigerant Quality
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3151'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  vRef=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3158'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  vRef=1/vRef !RS Comment: Refrigerant Specific Volume

  cpRef=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !RS Comment: Refrigerant Specific Heat
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3166'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  cpRef=cpRef/1000    !RS Comment: Unit Conversion

  muRef=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Dynamic Viscosity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3174'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  kRef=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !RS Comment: Refrigerant Thermal Conductivity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3181'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  kRef=kRef/1000  !RS Comment: Unit Conversion

  Temperature=tRef
  Quality=1
  IF (tRef+273.15 .GT. Tcr .OR. tRef+273.15 .LT. 0) THEN
    Psat=pRef
  ELSE
    Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !RS Comment: Saturation Pressure
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3194'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF
    Psat=Psat/1000  !RS Comment: Unit Conversion
  END IF

  SigmaRef=PQ(RefName, Pressure, Quality, 'surfacetension', RefrigIndex,RefPropErr)   !RS Comment: Surface Tension
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 4585'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  Pressure=pRef*1000  !RS Comment: Unit Conversion
  Quality=0
  hfRef=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !RS Comment: Liquid Enthalpy
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3205'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  hfRef=hfRef/1000    !RS Comment: Unit Conversion

  cpfRef=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !RS Comment: Liquid Specific Heat
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3213'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  cpfRef=cpfRef/1000  !RS Comment: Unit Conversion

  mufRef=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !RS Comment: Liquid Dynamic Viscosity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3221'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  kfRef=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !RS Comment: Liquid Thermal Conductivity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3228'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  kfRef=kfRef/1000    !RS Comment: Unit Conversion

  vfRef=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3236'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  vfRef=1/vfRef   !RS Comment: Liquid Specific Volume

  Pressure=pRef*1000  !RS Comment: Unit Conversion
  Quality=1
  hgRef=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !RS Comment: Vapor Enthalpy
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3246'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  hgRef=hgRef/1000    !RS Comment: Unit Conversion

  cpgRef=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !RS Comment: Vapor Specific Heat
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3254'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  cpgRef=cpgRef/1000  !RS Comment: Unit Conversion

  mugRef=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !RS Comment: Vapor Dynamic Viscosity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3262'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  kgRef=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !RS Comment: Vapor Thermal Conductivity
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3269'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  kgRef=kgRef/1000    !RS Comment: Unit Conversion

  vgRef=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3277'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  vgRef=1/vgRef   !RS Comment: Vapor Specific Volume

  Tsat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Saturation Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 4566'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  IF ((Tsat*1.8+32) .LT. 0) THEN !ISI - 10/11/06
    WRITE(*,*)'-- WARNING -- Condenser: Pressure drop error. Line 4933'
    ErrorFlag=DPERROR
    RETURN
  END IF

  hfgRef=hgRef-hfRef

  !Account for oil effect !ISI - 09/27/06
  IF (xRef .LT. 1 .AND. xRef .GT. 0 .AND. Wabsolute .GT. 0 .AND. LmodTPratio .EQ. 0) THEN
    Wlocal=LocalOilMassFraction(Wabsolute,xRef)
    Tsat=OilMixtureTsat(Wlocal,Pref*1e-3)
    cpfRef=OilMixtureSpecificHeat(CompManufacturer,Wlocal,cpfRef*1000,Tsat)*1e-3
    vfRef=1/OilMixtureDensity(CompManufacturer,Wlocal,1/vfRef,tRef)
    mufRef=OilMixtureViscosity(CompManufacturer,Wlocal,mufRef,MolWeight,tRef)
    SigmaRef=OilMixtureSurfaceTension(CompManufacturer,Wlocal,SigmaRef)
    kfRef=OilMixtureThermalConductivity(CompManufacturer,Wlocal,kfRef*1000)*1e-3
  END IF

  RETURN

END SUBROUTINE CalcRefProperty

!************************************************************************

SUBROUTINE CalcSegmentRefOutletPressure(CoilType,TubeType,pRi,hgRi,hfRi, & !CoilType,TubeType,tRi,pRi,hgRi,hfRi, &
  hRi,hRo,xRi,vRi,vgRi,vfRi,mRef, &
  !muRi,mugRi,mufRi,Sigma,Lsegment,LmodTPratio, & !RS: Debugging: Extraneous Sigma
  muRi,mugRi,mufRi,Lsegment,LmodTPratio, &
  IDtube,Elevation,Ltotal,DPrefMultiplier,pRo)

  !------------------------------------------------------------------------
  !Purpose:
  !To calculate segment refrigerant outlet pressure
  !
  !Author
  !Ipseng Iu
  !Oklahoma State University, Stillwater
  !
  !Date
  !March 2005
  !
  !Reference:
  !none
  !
  !------------------------------------------------------------------------

  USE FluidProperties_HPSim
  USE CoilCalcMod

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: CoilType  !1=Condenser; 2=Evaporator;
  !3=High side interconnecting pipes;
  !4=Low side interconnecting pipes
  !5=Microchannel condenser
  !6=Microchannel evaporator

  INTEGER, INTENT(IN) :: TubeType  !1=Plain; 2=General Micro Fin; 3=Herringbone;
  !4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
  !REAL, INTENT(IN) ::  tRi       !Inlet temperature, C   !RS: Debugging: Extraneous tRi
  REAL, INTENT(IN) ::  pRi       !Inlet pressure, kPa
  REAL, INTENT(IN) ::  hgRi      !Inlet vapor enthalpy, kJ/kg
  REAL, INTENT(IN) ::  hfRi      !Inlet liquid enthalpy, kJ/kg
  REAL, INTENT(IN) ::  hRi       !Inlet enthalpy, kJ/kg
  REAL, INTENT(IN) ::  hRo       !Outlet enthalpy, kJ/kg
  REAL, INTENT(IN) ::  xRi       !Inlet quality
  REAL, INTENT(IN) ::  vRi       !Inlet specific volume, m^3/kg
  REAL, INTENT(IN) ::  vgRi      !Inlet vapor specific volume, m^3/kg
  REAL, INTENT(IN) ::  vfRi      !Inlet liquid specific volume, m^3/kg
  REAL, INTENT(IN) ::  mRef      !Ref. mass flow rate, kg/s
  REAL, INTENT(IN) ::  muRi      !Inlet dynamic viscosity, Pa-s
  REAL, INTENT(IN) ::  mugRi     !Inlet vapor dynamic viscosity, Pa-s
  REAL, INTENT(IN) ::  mufRi     !Inlet liquid dynamic viscosity, Pa-s
  !REAL, INTENT(IN) ::  Sigma     !Surface tension, N/m   !RS: Debugging: Extraneous Sigma
  REAL, INTENT(IN) ::  Lsegment  !Segment length, m
  REAL, INTENT(IN) ::  LmodTPratio !Two-phase ratio
  REAL, INTENT(IN) ::  IDtube    !Tube inside diameter, m
  REAL, INTENT(IN) ::  Elevation !Elevation, m
  REAL, INTENT(IN) ::  Ltotal    !Total tube length, m
  REAL, INTENT(IN) ::  DPrefMultiplier !Pressure drop multiplier
  REAL, INTENT(OUT) :: pRo       !Outlet pressure, kPa

  REAL tRo       !Outlet temperature, C
  REAL vgRo      !Inlet vapor specific volume, m^3/kg
  REAL vfRo      !Inlet liquid specific volume, m^3/kg
  REAL vRo !Outlet specific volume, m^3/kg
  REAL xRo !Outlet quality
  REAL pRoPrev !Previous value of pRo, for iteration

  !FLOW:

  !Find outlet ref. pressure
  pRo=pRi !Initialize
  pRoPrev=BIG !Initialize
  Counter=0
  DO

    Pressure=pRo*1000   !RS Comment: Unit Conversion
    Enthalpy=hRo*1000   !RS Comment: Unit Conversion
    xRo=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Outlet Quality
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 1377'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF

    tRo=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Outlet Temperature
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3377'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF

    vRo=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 1386'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF
    vRo=1/vRo   !RS Comment: Refrigerant Outlet Specific Volume

    Pressure=pRo*1000   !RS Comment: Unit Conversion
    Quality=1
    vgRo=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
    vgRo=1/vgRo !RS Comment: Inlet Vapor Specific Volume
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3394'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF

    Pressure=pRo*1000   !RS Comment: Unit Conversion
    Quality=0
    vfRo=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
    vfRo=1/vfRo !RS Comment: Inlet Liquid Specific Volume
    IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3403'
      ErrorFlag=REFPROPERROR
      RETURN
    END IF

    CALL MODdP(CoilType,TubeType,hgRi,hfRi, &
    hRi,hRo,xRi,xRo,vRi,vRo,vgRi,vfRi,vgRo,vfRo,mRef,muRi,mugRi,mufRi, &
    Lsegment,LmodTPratio,IDtube,Elevation,Ltotal,dPfric,dPmom,dPgrav)
    !(CoilType,TubeType,tRi,tRo,pRi,hg,hf,hRi,hRo,xRi,xRo, &
    !                 vRi,vRo,vgi,vfi,vgo,vfo,mRef,muRef,mug,muf,Sigma, &
    !				 Lmod,LmodTPratio,ID,OD,HtCoil,Lcoil,dPfric,dPmom,dPgrav)

    IF (ABS(dPfric) .LT. ABS(dPmom)) THEN
      dPmom=0
    END IF

    dPmod=ABS(dPfric)+ABS(dPmom)+dPgrav
    dPmod=dPmod*DPrefMultiplier

    IF (pRi-dPmod .LT. 0) THEN
      dPmod=0
    END IF

    pRo=pRi-dPmod
    IF (ABS(pRo-pRoPrev)/pRoPrev .LT. SMALL .OR. Counter .GT. PressureMaxIter) THEN
      EXIT
    ELSE
      pRoPrev=pRo
      Counter=Counter+1
    END IF

  END DO !end of pRo

  RETURN

END SUBROUTINE CalcSegmentRefOutletPressure

!************************************************************************

SUBROUTINE MicrochannelCondenser(XIN,PAR,OUT) !(Ref$,XIN,PAR,OUT)   !RS: Debugging: Extraneous Ref$

  !-----------------------------------------------------------------------------------
  !
  !  Description:
  !  A segment-by-segment microchannel condenser model
  !  To predict coil air side and refrigerant side properties, heat transfer,
  !  and pressure drop
  !
  !  Inputs:
  !  Ref$=Refrigerant name
  !  XIN(1)=Refrigerant side mass flow rate, kg/s
  !  XIN(2)=Refrigerant side inlet (compressor outlet) pressure, kPa
  !  XIN(3)=Refrigerant side inlet (compressor outlet) enthalpy, kJ/kg
  !  XIN(4)=Air side mass flow rate, kg/s
  !  XIN(5)=Air side inlet temp. C
  !  XIN(6)=Air side inlet relative humidity
  !  XIN(7)=Evaporator outlet temperature, C
  !
  !  Parameters:
  !  PAR(1)=Barometric pressure, kPa
  !  PAR(2)=Cooling mode? 1=yes; 0=no
  !  PAR(3)=Discharge line length, m
  !  PAR(4)=Discharge line outside diameter, m
  !  PAR(5)=Discharge line tube wall thickness, m
  !  PAR(6)=Discharge line elevation, m
  !  PAR(7)=Discharge line heat loss, kW
  !  PAR(8)=Discharge line temperature change, C
  !  PAR(9)=Discharge line additional pressure drop, kPa
  !  PAR(10)=Liquid line length, m
  !  PAR(11)=Liquid line outside diameter, m
  !  PAR(12)=Liquid line tube wall thickness, m
  !  PAR(13)=Liquid line elevation, m
  !  PAR(14)=Liquid line heat loss, kW
  !  PAR(15)=Liquid line temperature change, C
  !  PAR(16)=Liquid line additional pressure drop, kPa
  !  PAR(17)=Multiplier for ref. side heat transfer correlation
  !  PAR(18)=Multiplier for ref. side pressure drop correlation
  !  PAR(19)=Multiplier for air side heat transfer correlation
  !  PAR(20)=Multiplier for air side pressure drop correlation
  !  PAR(21)=Fan power, kW
  !  PAR(22)=Fan location, 1=draw through; 2=blow through
  !  PAR(23)=Compressor heat loss, kW
  !  PAR(24)=Is compressor in air stream, 1=yes, 0=no
  !  PAR(25)=Custom air side data unit, 1=SI; 2=IP
  !  PAR(26)=Custom air heat transfer curve type, 1=Power; 2=Polynomial
  !  PAR(27)=Power coefficient for air heat transfer curve
  !  PAR(28)=Power coefficient for air heat transfer curve
  !  PAR(29)=Polynomial coefficient for air heat transfer curve
  !  PAR(30)=Polynomial coefficient for air heat transfer curve
  !  PAR(31)=Polynomial coefficient for air heat transfer curve
  !  PAR(32)=Polynomial coefficient for air heat transfer curve
  !  PAR(33)=Custom air heat transfer curve type, 1=Power; 2=Polynomial
  !  PAR(34)=Power coefficient for air heat transfer curve
  !  PAR(35)=Power coefficient for air heat transfer curve
  !  PAR(36)=Polynomial coefficient for air heat transfer curve
  !  PAR(37)=Polynomial coefficient for air heat transfer curve
  !  PAR(38)=Polynomial coefficient for air heat transfer curve
  !  PAR(39)=Polynomial coefficient for air heat transfer curve
  !
  !  Outputs:
  !  OUT(1)=Coil capacity, kW
  !  OUT(2)=Coil inlet pressure, kPa
  !  OUT(3)=Coil inlet enthalpy, kJ/kg
  !  OUT(4)=Coil inlet temperature, C
  !  OUT(5)=Coil inlet quality
  !  OUT(6)=Coil outlet pressure, kPa
  !  OUT(7)=Coil outlet enthalpy, kJ/kg
  !  OUT(8)=Coil outlet temperature, C
  !  OUT(9)=Coil outlet quality
  !  OUT(10)=Coil outlet subcooling, C
  !  OUT(11)=Liquid line outlet pressure, kPa
  !  OUT(12)=Liquid line outlet enthalpy, kJ/kg
  !  OUT(13)=Liquid line outlet temperature, C
  !  OUT(14)=Liquid line outlet quality
  !  OUT(15)=Liquid line outlet subcooling, C
  !  OUT(16)=Air side outlet temperature, C
  !  OUT(17)=Air side outlet relative humidity
  !  OUT(18)=Air side pressure drop, kPa
  !  OUT(19)=Aluminum weight, kg
  !  OUT(20)=Error flag: 0-No error
  !                      1-Condenser solution not converge
  !                      2-Refprop error
  !					   3-Circuit file error
  !  OUT(21)=Mass in discharge line, kg
  !  OUT(22)=Mass in liquid line, kg
  !  OUT(23)=Mass in coil, kg
  !  OUT(24)=Liquid mass in coil, kg
  !  OUT(25)=Vapor mass in coil, kg
  !
  !  Reference:
  !  none
  !
  !  Author:
  !  Ipseng Iu
  !  Mechanical and Aerospace Engineering
  !  Oklahoma State University, Stillwater
  !
  !  Date: November 2005
  !
  !-----------------------------------------------------------------------------------

  USE FluidProperties_HPSim
  USE CoilCalcMod
  USE AirPropMod
  USE OilMixtureMod
  USE InputProcessor    !RS: Debugging: GetObjectItem

  IMPLICIT NONE

  !Subroutine argument declarations
  !CHARACTER*80, INTENT(IN)  :: Ref$  !RS: Debugging: Extraneous Ref$
  REAL,         INTENT(IN)  :: XIN(7)
  REAL,         INTENT(IN)  :: PAR(23) !RS: Debugging: Formerly PAR(39)
  REAL,         INTENT(OUT) :: OUT(22)

  !Subroutine local variables
  INTEGER,PARAMETER :: CoilType = MCCONDENSER

  INTEGER I,II,III,IV,V !Loop counters
  LOGICAL Converged     !Convergence flag
  INTEGER RefBCiter     !Refrigerant bounadary condition iteration counter
  REAL Qpass            !Coil pass capacity, kW
  REAL QinletPass       !Inlet pass capacity, kW
  REAL pRoSlab       !Outlet Refrigerant pressure for a coil slab, kPa
  REAL hRoSlab       !Outlet Refrigerant enthalpy for a coil slab, kJ/kg
  REAL tRoSlab       !Outlet Refrigerant temperature for a coil slab, C
  REAL xRoSlab       !Outlet Refrigerant quality for a coil slab
  REAL Aface         !Coil face area, m^2
  REAL SumPro        !Sum of outlet pressures, kPa
  REAL SumMrefHro    !Sum of mdot*H (mass flow rate * enthalpy)
  REAL DPcoil, DPcoilPrev !Coil pressure drop, kPa
  REAL mdothRo !mdot x outlet enthalpy

  INTEGER, PARAMETER :: MaxNameLength = 200   !RS: Debugging: Bringing through

  CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
  INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
  REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
  INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
  INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
  INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

  !FLOW:

  mRefTot          =XIN(1)
  pRoCmp           =XIN(2)
  hRoCmp           =XIN(3)
  mAiCoil          =XIN(4)
  tAiCoil          =XIN(5)
  rhAiCoil         =XIN(6)
  !tRoEvp           =XIN(7)   !RS: Debugging: Set once but never used

  BaroPressure     =PAR(1)
  !IsCoolingMode    =PAR(2)   !RS: Debugging: Global variable now
  LdisLn           =PAR(3)
  ODdisLn          =PAR(4)
  DisLnThk         =PAR(5)
  !ElevDisLn        =PAR(6)   !RS: Debugging: Never Used?
  QdisLn           =PAR(7)
  !DTdisLn          =PAR(8)   !RS: Debugging: Never Used?
  !AddDPdisLn       =PAR(9)   !RS: Debugging: Never Used?
  LliqLn           =PAR(10)
  ODliqLn          =PAR(11)
  LiqLnThk         =PAR(12)
  ElevLiqLn        =PAR(13)
  QliqLn           =PAR(14)
  DTliqLn          =PAR(15)
  !AddDPLiqLn       =PAR(16)  !RS: Debugging: Never Used?
  hciMultiplier    =PAR(17)
  DPrefMultiplier  =PAR(18)
  hcoMultiplier    =PAR(19)
  DPairMultiplier  =PAR(20)
  PwrFan           =PAR(21)
  DrawBlow         =PAR(22)
  QlossCmp         =PAR(23)
  !IsCmpInAirStream =PAR(24)
  !CurveUnit        =PAR(25)  !RS: Debugging: Never Really Used
  !CurveTypeHTC     =PAR(26)  !RS: Debugging: Never Really Used
  !PowerAHTC        =PAR(27)  !RS: Debugging: Never Really Used
  !PowerBHTC        =PAR(28)  !RS: Debugging: Never Really Used
  !Poly1HTC         =PAR(29)  !RS: Debugging: Never Really Used
  !Poly2HTC         =PAR(30)  !RS: Debugging: Never Really Used
  !Poly3HTC         =PAR(31)  !RS: Debugging: Never Really Used
  !Poly4HTC         =PAR(32)  !RS: Debugging: Never Really Used
  !CurveTypeDP      =PAR(33)  !RS: Debugging: Never Really Used
  !PowerADP         =PAR(34)  !RS: Debugging: Never Really Used
  !PowerBDP         =PAR(35)  !RS: Debugging: Never Really Used
  !Poly1DP          =PAR(36)  !RS: Debugging: Never Really Used
  !Poly2DP          =PAR(37)  !RS: Debugging: Never Really Used
  !Poly3DP          =PAR(38)  !RS: Debugging: Never Really Used
  !Poly4DP          =PAR(39)  !RS: Debugging: Never Really Used

  !********************Refrigerant Cycle Data (Heating)***********************  !RS: Debugging: Moving: Condenser

  CALL GetObjectItem('RefrigerantCycleData(Heating)',1,Alphas,NumAlphas, &
  TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

  IsCmpInAirStream = Numbers(2) !Is Compressor in Air Stream

  IsParallelSlabs = 1

  ErrorFlag=NOERROR !Initialize

  !Coil height What is this? Is this right for a microchannel? Yes - ISI - 12/24/2009
  HtCoil=Nt*Pt

  !Coil length
  Lcoil=Nl*Nt*Ltube

  !Face area
  Aface=Ltube*Nt*Pt

  !Tube information
  LmodTube=Ltube/NumOfMods

  CPair=CPA(REAL(tAiCoil))    !RS: Replace: CPA (2/19/14)
  CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
  Cair=mAiCoil*CPAir

  AirPropOpt=2
  AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
  AirProp%APRelHum=rhAiCoil !RS: Debugging: Formerly AirProp(3)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
  hAiCoil=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)
  wbAiCoil=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)

  IF (DrawBlow .EQ. BLOWTHROUGH) THEN !Blow through
    tAiCoil=tAiCoil+PwrFan/Cair
    hAiCoil=hAiCoil+PwrFan/mAiCoil
  END IF
  IF (IsCmpInAirStream .NE. 0) THEN !Compressor in air stream
    tAiCoil=tAiCoil+QlossCmp/Cair
    hAiCoil=hAiCoil+QlossCmp/mAiCoil
  END IF

  AirPropOpt=1
  AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
  AirProp%APEnth=hAiCoil  !RS: Debugging: Formerly AirProp(4)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
  DensityIn=AirProp%APDryDens    !RS: Debugging: Formerly AirProp(7)

  !Area calculations
  CALL CalcCoilArea(CoilType,Nl,Nt,Pt,Pl,TubeDepth, &
  Ltube,IDtube,TubeHeight,Dchannel,NumOfChannels, &
  FinThk,FinSpg,Lcoil,AfCoil, &
  AoCoil,AiCoil,AmCoil)

  !****** Discharge line calculation ******
  IF (LdisLn .GT. 0) THEN
    CALL DischargeLine
    IF (ErrorFlag .GT. CONVERGEERROR) THEN
      WRITE(*,*)'-- WARNING -- DischargeLine: Refprop error. Line 4017'
      RETURN
    END IF
  ELSE
    pRiCoil=pRoCmp
    hRiCoil=hRoCmp
  END IF

  Pressure=pRiCoil*1000   !RS Comment: Unit Conversion
  Enthalpy=hRiCoil*1000   !RS Comment: Unit Conversion
  tRiCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Inlet Coil Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error. Line 3904'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF
  xRiCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Inlet Coil Quality
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error. Line 3910'
    ErrorFlag=REFPROPERROR
    RETURN
  END IF

  !Initialize boundary conditions
  DO I=1, Nl
    IF (IsParallelSlabs .GT. 0) THEN
      Slab(I)%mdot=mRefTot/Nl !ISI - 07/14/07
    ELSE
      Slab(I)%mdot=mRefTot
    END IF
    DO II=1, Slab(I)%Npass
      DO III=1, Slab(I)%Pass(II)%Ntube
        !Populating the slab segment property arrays
        Slab(I)%Pass(II)%Tube(III)%Seg%mAi=mAiCoil*LmodTube/(Ltube*Nt)
        Slab(I)%Pass(II)%Tube(III)%Seg%tAi=tAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%wbAi=wbAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%rhAi=rhAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%tAo=tAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%wbAo=wbAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%rhAo=rhAiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%pRi=pRiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%hRi=hRiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%pRo=pRiCoil
        Slab(I)%Pass(II)%Tube(III)%Seg%hRo=hRiCoil

      END DO
    END DO
  END DO

  !Calc. air side mass flow rate at the front row
  VelAvg=(mAiCoil/DensityIn)/Aface
  DO I=1, Nl
    DO II=1, Slab(I)%Npass
      DO III=1, Slab(I)%Pass(II)%Ntube
        DO IV=1, NumOfMods
          Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Aface=LmodTube/(Ltube*Nt)*Aface
          IF (Slab(I)%Pass(II)%Tube(III)%Seg(IV)%VelDev .LE. 0) THEN
            Slab(I)%Pass(II)%Tube(III)%Seg(IV)%VelDev=1
          END IF
          Slab(I)%Pass(II)%Tube(III)%Seg(IV)%mAi=mAiCoil*LmodTube/(Ltube*Nt)* &
          Slab(I)%Pass(II)%Tube(III)%Seg(IV)%VelDev
        END DO
      END DO
    END DO
  END DO

  !****** Start coil calculation ******
  Converged=.TRUE.

  !Initialize
  QmodPrev=QdisLn
  DPcoilPrev=1e20

  Qcoil=0.0

  DO RefBCiter=1, MdotMaxIter

    DO I=1, Nl !Number of slabs

      DO II=1,Slab(I)%Npass !Number of passes

        Qpass=0.0

        IF (IsSimpleCoil .EQ. 1) THEN
          IsSimpleCoil=2  !If IsSimpleCoil equals 1 then the following subroutine can't handle the microchannel
          CALL CalcCircuitRefInletConditions(I,II,CoilType)
          IsSimpleCoil=1
        ELSE
          CALL CalcCircuitRefInletConditions(I,II,CoilType)
        END IF

        IF (II .EQ. 1 .AND. Slab(I)%Ninlet .GT. 1) THEN !Multi-inlet

          SumPro=0; SumMrefHro=0;

          DO V=1, Slab(I)%Ninlet

            QinletPass=0

            Slab(I)%InletPass(V)%pRi=pRiCoil
            Slab(I)%InletPass(V)%hRi=hRiCoil

            DO III=1,1 !Slab(I)%Pass(II)%NumOfTubes !Number of tubes
              DO IV=1,NumOfMods !Number of segments

                !ref. mass flow rate for each channel
                mRefMod=Slab(I)%mdot/Slab(I)%Pass(II)%Ntube/NumOfChannels !ISI - 07/13/07

                !Surface areas
                AoMod=AoCoil*LmodTube/Lcoil
                AfMod=AfCoil*LmodTube/Lcoil
                AiMod=AiCoil*LmodTube/Lcoil
                AmMod=AmCoil*LmodTube/Lcoil

                CALL CalcCoilSegment(I,II,III,IV,CoilType)
                IF (ErrorFlag .GT. CONVERGEERROR) THEN
                  OUT(20)=ErrorFlag
                  RETURN
                END IF

                QinletPass=QinletPass+Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod

              END DO !End mod
            END DO !End tube

            Slab(I)%InletPass(V)%Qpass=QinletPass*Slab(I)%InletPass(V)%Ntube

            Slab(I)%InletPass(V)%pRo=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%pRo
            Slab(I)%InletPass(V)%hRo=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%hRo
            Slab(I)%InletPass(V)%mRef=mRefTot*Slab(I)%InletPass(V)%Ntube/Slab(I)%Pass(II)%Ntube

            Qpass=Qpass+Slab(I)%InletPass(V)%Qpass

            SumPro=SumPro+Slab(I)%InletPass(V)%pRo !Sum of outlet pressures
            SumMrefHro=SumMrefHro+Slab(I)%InletPass(V)%mRef*Slab(I)%InletPass(V)%hRo !Sum of outlet mdot*hro

          END DO !End inlet pass

          pRoCkt=SumPro/Slab(I)%Ninlet
          hRoCkt=SumMrefHro/mRefTot

        ELSE !Single inlet

          DO III=1,1 !Slab(I)%Pass(II)%NumOfTubes !Number of tubes

            DO IV=1,NumOfMods !Number of segments

              !ref. mass flow rate for each channel
              mRefMod=Slab(I)%mdot/Slab(I)%Pass(II)%Ntube/NumOfChannels !ISI - 07/13/07

              !Surface areas
              AoMod=AoCoil*LmodTube/Lcoil
              AfMod=AfCoil*LmodTube/Lcoil
              AiMod=AiCoil*LmodTube/Lcoil
              AmMod=AmCoil*LmodTube/Lcoil

              CALL CalcCoilSegment(I,II,III,IV,CoilType)
              IF (ErrorFlag .GT. CONVERGEERROR) THEN
                OUT(20)=ErrorFlag
                RETURN
              END IF

              !Calc. circuit heat transfer

            END DO !End mod

          END DO !End tube

          pRoCkt=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%pRo !Circuit outlet pressure
          hRoCkt=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%hRo !Circuit outlet enthalpy

        END IF

        Slab(I)%Pass(II)%pRo=pRoCkt
        Slab(I)%Pass(II)%hRo=hRoCkt

      END DO !End pass

      Slab(I)%pRi=Slab(I)%Pass(1)%Tube(1)%Seg(1)%pRi
      Slab(I)%hRi=Slab(I)%Pass(1)%Tube(1)%Seg(1)%hRi

      Slab(I)%tAi=Slab(I)%Pass(1)%Tube(1)%Seg(1)%tAi
      Slab(I)%rhAi=Slab(I)%Pass(1)%Tube(1)%Seg(1)%rhAi

      AirPropOpt=2
      AirProp%APTDB=Slab(I)%tAi  !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=Slab(I)%rhAi !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      Slab(I)%hAi=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)

      IF (Slab(1)%Npass .EQ. 1 .AND. Slab(1)%Ninlet .GT. 1) THEN !1-pass, multi-inlet
        pRoSlab=pRoCkt
        hRoSlab=hRoCkt
      ELSE !Single inlet
        pRoSlab=Slab(I)%Pass(Slab(I)%Npass)%pRo
        hRoSlab=Slab(I)%Pass(Slab(I)%Npass)%hRo
      END IF

      Pressure=pRoSlab*1000   !RS Comment: Unit Conversion
      Enthalpy=hRoSlab*1000   !RS Comment: Unit Conversion
      tRoSlab=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Outlet Slab Temperature
      IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(20)=ErrorFlag
        RETURN
      END IF
      xRoSlab=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Refrigerant Outlet Slab Quality
      IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(20)=ErrorFlag
        RETURN
      END IF

      Slab(I)%pRo=pRoSlab
      Slab(I)%hRo=hRoSlab
      Slab(I)%tRo=tRoSlab
      Slab(I)%xRo=xRoSlab

      Slab(I)%Qslab=Slab(I)%mdot*(Slab(I)%hRi-Slab(I)%hRo) !ISI - 01/06/08

      CALL CalcMeanProp(Slab(I)%tAi,Slab(I)%tAo,tAmod)

      !Coil air side outlet conditions
      CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
      Cair=mAicoil*CPAir

      Slab(I)%tAo=Slab(I)%tAi+Slab(I)%Qslab/Cair
      Slab(I)%hAo=Slab(I)%hAi+Slab(I)%Qslab/mAiCoil

      AirPropOpt=1
      AirProp%APTDB=Slab(I)%tAo  !RS: Debugging: Formerly AirProp(1)
      AirProp%APEnth=Slab(I)%hAo  !RS: Debugging: Formerly AirProp(4)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
      Slab(I)%rhAo=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)

    END DO !End Slabs

    !ISI - 07/14/07
    IF (Nl .GT. 1) THEN

      IF (IsParallelSlabs .GT. 0) THEN
        CALL UpdateMCrefMassFlowRate(Slab,Nl,mRefTot,DPcoil)

        IF (ABS(DPcoil-DPcoilPrev)/DPcoilPrev .GT. SMALL) THEN
          DPcoilPrev=DPcoil
        ELSE

          !Calculate outlet pressure and enthalpy
          pRoCoil=pRiCoil-DPcoil
          mdothRo=0
          DO V=1, Nl
            mdothRo=mdothRo+Slab(V)%mdot*Slab(V)%hRo
          END DO
          hRoCoil=mdothRo/mRefTot

          EXIT
        END IF

      ELSE !Series

        !Calculate outlet pressure and enthalpy
        DPcoil=Slab(1)%pRi-Slab(Nl)%pRo
        pRoCoil=Slab(Nl)%pRo
        hRoCoil=Slab(Nl)%hRo
        EXIT

      END IF
    ELSE
      !Calculate outlet pressure and enthalpy
      DPcoil=Slab(1)%pRi-Slab(1)%pRo
      pRoCoil=Slab(1)%pRo
      hRoCoil=Slab(1)%hRo
      EXIT
    END IF

  END DO !End iter

  Qcoil=mRefTot*(hRiCoil-hRoCoil)

  !Coil air side outlet conditions
  CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
  CPair=CPAirFunction(tAMod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
  Cair=mAicoil*CPAir

  tAoCoil=tAiCoil+Qcoil/Cair
  hAoCoil=hAiCoil+Qcoil/mAiCoil

  !Fan air side inlet conditions
  CPair=CPA(REAL(tAoCoil))    !RS: Replace: CPA (2/19/14)
  CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
  Cair=mAiCoil*CPAir

  IF (DrawBlow .EQ. DRAWTHROUGH) THEN !Draw through
    tAoCoil=tAoCoil+PwrFan/Cair
    hAoCoil=hAoCoil+PwrFan/mAiCoil
  END IF

  AirPropOpt=1
  AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
  AirProp%APEnth=hAiCoil  !RS: Debugging: Formerly AirProp(4)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
  rhAiCoil=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
  DensityIn=AirProp%APDryDens    !RS: Debugging: Formerly AirProp(7)

  AirPropOpt=1
  AirProp%APTDB=tAoCoil  !RS: Debugging: Formerly AirProp(1)
  AirProp%APEnth=hAoCoil  !RS: Debugging: Formerly AirProp(4)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
  rhAoCoil=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
  DensityOut=AirProp%APDryDens   !RS: Debugging: Formerly AirProp(7)

  WetFlag=0
  RowNum=0
  CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
  IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair, &
  hAoCoil)

  !CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
  !IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC, &
  !Poly1HTC,Poly2HTC,Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP, &
  !Poly1DP,Poly2DP,Poly3DP,Poly4DP,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair)

  DPair=DPair*DPairMultiplier

  Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
  Enthalpy=hRoCoil*1000   !RS Comment: Unit Conversion
  tRoCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Coil Outlet Refrigerant Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF
  xRoCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Coil Outlet Refrigerant Quality
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF

  Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
  Quality=0
  tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF

  IF (xRoCoil .LE. 0.0) THEN
    tSCoCoil=tSat-tRoCoil !Subcooling
  ELSE
    tSCoCoil=0.0
  END IF

  !****** Liquid line calculation ******
  IF (LliqLn .GT. 0) THEN
    CALL LiquidLine
    IF (ErrorFlag .GT. CONVERGEERROR) THEN
      WRITE(*,*)'LiquidLine: Refprop error.'
      OUT(20)=ErrorFlag
      RETURN
    END IF
  ELSE
    pRiExp=pRoCoil
    hRiExp=hRoCoil
  END IF

  Pressure=pRiExp*1000    !RS Comment: Unit Conversion
  Enthalpy=hRiExp*1000    !RS Comment: Unit Conversion
  tRiExp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Refrigerant Inlet Temperature
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF
  xRiExp=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Expansion Device Refrigerant Inlet Quality
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF

  Pressure=pRiExp*1000    !RS Comment: Unit Conversion
  Quality=0
  tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Pressure
  IF (RefPropErr .GT. 0) THEN
    WRITE(*,*)'-- WARNING -- MCCondenser: Refprop error.'
    ErrorFlag=REFPROPERROR
    OUT(20)=ErrorFlag
    RETURN
  END IF

  IF (xRiExp .LE. 0) THEN
    tSCiExp=tSat-tRiExp
  ELSE
    tSCiExp=0
  END IF

  OUT(1)=Qcoil
  OUT(2)=pRiCoil
  OUT(3)=hRiCoil
  OUT(4)=tRiCoil
  OUT(5)=xRiCoil
  OUT(6)=pRoCoil
  OUT(7)=hRoCoil
  OUT(8)=tRoCoil
  OUT(9)=xRoCoil
  OUT(10)=tSCoCoil
  OUT(11)=pRiExp
  OUT(12)=hRiExp
  OUT(13)=tRiExp
  OUT(14)=xRiExp
  OUT(15)=tSCiExp
  OUT(16)=tAoCoil
  OUT(17)=rhAoCoil
  OUT(18)=DPair
  OUT(19)=WeightAluminum

  OUT(21)=MassDisLn
  OUT(22)=MassLiqLn
  OUT(20)=ErrorFlag

  RETURN

END SUBROUTINE MicrochannelCondenser

!************************************************************************

SUBROUTINE LoadMicrochannelInputs(MCXIN,MCPAR) !(FTXIN,FTPAR,MCXIN,MCPAR)

  !-----------------------------------------------------------------------------------
  !
  !  Description:
  !  Transfer input data from subroutine "Condenser" to "MicrochannelCondenser"
  !
  !  Author:
  !  Ipseng Iu
  !  Mechanical and Aerospace Engineering
  !  Oklahoma State University, Stillwater
  !
  !  Date: November 2005
  !
  !-----------------------------------------------------------------------------------

  IMPLICIT NONE

  !REAL, INTENT(IN)  :: FTXIN(7)  !Fin-tube coil input data    !RS: Debugging: Formerly FTXIN(8)
  !REAL, INTENT(IN)  :: FTPAR(45) !Fin-tube coil input parameters  !RS: Debugging: Formerly FTPAR(55)
  REAL, INTENT(OUT) :: MCXIN(7)  !Microchannel coil input data
  REAL, INTENT(OUT) :: MCPAR(24) !Microchannel coil input parameters  !RS: Debugging: Formerly MCPAR(39)

  !FLOW:

  MCXIN(1)=CondIN%CInmRef !Refrigerant side mass flow rate, kg/s    !RS: Debugging: Formerly FTXIN(1)
  MCXIN(2)=CondIN%CInpRo !Refrigerant side inlet (compressor outlet) pressure, kPa !RS: Debugging: Formerly FTXIN(2)
  MCXIN(3)=CondIN%CInhRo !Refrigerant side inlet (compressor outlet) enthalpy, kJ/kg   !RS: Debugging: Formerly FTXIN(3)
  MCXIN(4)=CondIN%CInmAi !Air side mass flow rate, kg/s    !RS: Debugging: Formerly FTXIN(4)
  MCXIN(5)=CondIN%CIntAi !Air side inlet temp. C   !RS: Debugging: Formerly FTXIN(5)
  MCXIN(6)=CondIN%CInrhAi !Air side inlet relative humidity !RS: Debugging: Formerly FTXIN(6)

  MCPAR(1)=CondPAR%CondBarPress !Barometric pressure, kPa    !RS: Debugging: Formerly FTPAR(38)
  MCPAR(2)=CondPAR%CondCoolMode !Cooling mode? 1=yes; 0=no   !RS: Debugging: Formerly FTPAR(27)
  MCPAR(3)=CondPAR%CondDisLnLen  !Discharge line length, m    !RS: Debugging: Formerly FTPAR(1)
  MCPAR(4)=CondPAR%CondDisLnOD  !Discharge line outside diameter, m  !RS: Debugging: Formerly FTPAR(2)
  MCPAR(5)=CondPAR%CondDisLnTWThick  !Discharge line tube wall thickness, m   !RS: Debugging: Formerly FTPAR(3)
  MCPAR(6)=CondPAR%CondDisLnElev  !Discharge line elevation, m !RS: Debugging: Formerly FTPAR(4)
  MCPAR(7)=CondPAR%CondDisLnQLoss  !Discharge line heat loss, kW    !RS: Debugging: Formerly FTPAR(5)
  MCPAR(8)=CondPAR%CondDisLnTempChg  !Discharge line temperature change, C    !RS: Debugging: Formerly FTPAR(6)
  MCPAR(9)=CondPAR%CondDisLnAddPD  !Discharge line additional pressure drop, kPa    !RS: Debugging: Formerly FTPAR(7)
  MCPAR(10)=CondPAR%CondLiqLnLen !Liquid line length, m   !RS: Debugging: Formerly FTPAR(8)
  MCPAR(11)=CondPAR%CondLiqLnOD !Liquid line outside diameter, m !RS: Debugging: Formerly FTPAR(9)
  MCPAR(12)=CondPAR%CondLiqLnTWThick !Liquid line tube wall thickness, m     !RS: Debugging: Formerly FTPAR(10)
  MCPAR(13)=CondPAR%CondLiqLnElev !Liquid line elevation, m   !RS: Debugging: Formerly FTPAR(11)
  MCPAR(14)=CondPAR%CondLiqLnQLoss !Liquid line heat loss, kW  !RS: Debugging: Formerly FTPAR(12)
  MCPAR(15)=CondPAR%CondLiqLnTempChg !Liquid line temperature change, C  !RS: Debugging: Formerly FTPAR(13)
  MCPAR(16)=CondPAR%CondLiqLnAddPD !Liquid line additional pressure drop, kPa  !RS: Debugging: Formerly FTPAR(14)
  MCPAR(17)=CondPAR%CondMultRefQT !Multiplier for ref. side heat transfer correlation !RS: Debugging: Formerly FTPAR(30)
  MCPAR(18)=CondPAR%CondMultRefPD !Multiplier for ref. side pressure drop correlation !RS: Debugging: Formerly FTPAR(31)
  MCPAR(19)=CondPAR%CondMultAirQT !Multiplier for air side heat transfer correlation  !RS: Debugging: Formerly FTPAR(32)
  MCPAR(20)=CondPAR%CondMultAirPD !Multiplier for air side pressure drop correlation  !RS: Debugging: Formerly FTPAR(33)
  MCPAR(21)=CondPAR%CondFanPwr !Fan power, kW  !RS: Debugging: Formerly FTPAR(34)
  MCPAR(22)=CondPAR%CondFanLoc !Fan location, 1=draw through; 2=blow through   !RS: Debugging: Formerly FTPAR(35)
  MCPAR(23)=CondPAR%CondCompQLoss !Compressor heat loss, kW   !RS: Debugging: Formerly FTPAR(39)
  MCPAR(24)=CondPAR%CondPressTolConv !Is compressor in air stream, 1=yes, 0=no   !RS: Debugging: Formerly FTPAR(40)

  RETURN

END SUBROUTINE LoadMicrochannelInputs

!************************************************************************

SUBROUTINE LoadMicrochannelOutputs(MCOUT) !,FTOUT)

  !-----------------------------------------------------------------------------------
  !
  !  Description:
  !  Transfer output data from subroutine "MicrochannelCondenser" to "Condenser"
  !
  !  Author:
  !  Ipseng Iu
  !  Mechanical and Aerospace Engineering
  !  Oklahoma State University, Stillwater
  !
  !  Date: November 2005
  !
  !-----------------------------------------------------------------------------------

  IMPLICIT NONE

  REAL, INTENT(IN)  :: MCOUT(22)  !Microchannel coil output data
  !REAL, INTENT(OUT) :: FTOUT(20)  !Fin-tube coil output data  !RS: Debugging: Formerly FTOUT(29), FTOUT(24)

  !FLOW:

  CondOUT%COutpRiC=MCOUT(2)   !Coil inlet pressure, kPa
  CondOUT%COuthRiC=MCOUT(3)   !Coil inlet enthalpy, kJ/kg
  !FTOUT(3)=MCOUT(4)   !Coil inlet temperature, C
  !FTOUT(4)=MCOUT(5)   !Coil inlet quality
  CondOUT%COutpRoC=MCOUT(6)   !Coil outlet pressure, kPa
  CondOUT%COuthRoC=MCOUT(7)   !Coil outlet enthalpy, kJ/kg
  CondOUT%COuttRoC=MCOUT(8)   !Coil outlet temperature, C
  !CondOUT%8=MCOUT(9)   !Coil outlet quality  !RS: Debugging: Not really used?
  !FTOUT(9)=MCOUT(10)  !Coil outlet subcooling, C !RS: Debugging: Never used
  CondOUT%COutpRiE=MCOUT(11) !Liquid line outlet pressure, kPa
  CondOUT%COuthRiE=MCOUT(12) !Liquid line outlet enthalpy, kJ/kg
  CondOUT%COuttRiE=MCOUT(13) !Liquid line outlet temperature, C
  CondOUT%COutxRiE=MCOUT(14) !Liquid line outlet quality
  CondOUT%COuttSCiE=MCOUT(15) !Liquid line outlet subcooling, C
  CondOUT%COutQC=MCOUT(1)  !Coil capacity, kW
  CondOUT%COutMDisLn=MCOUT(21) !Mass in discharge line, kg
  CondOUT%COutMLiqLn=MCOUT(22) !Mass in liquid line, kg
  CondOUT%COutMC=0         !Mass in coil, kg
  !FTOUT(19)=0         !Liquid mass in coil, kg   !RS: Debugging: Never used
  !FTOUT(20)=0         !Vapor mass in coil, kg    !RS: Debugging: Never used
  CondOUT%COuttAoC=MCOUT(16) !Air side outlet temperature, C !RS: Debugging: Formerly FTOUT(21)
  CondOUT%COutrhAoC=MCOUT(17) !Air side outlet relative humidity   !RS: Debugging: Formerly FTOUT(22)
  CondOUT%COutDPAir=MCOUT(18) !Air side pressure drop, kPa    !RS: Debugging: Formerly FTOUT(23)
  CondOUT%COutErrFlag=MCOUT(20) !Error flag !RS: Debugging: Formerly FTOUT(24)
  !FTOUT(25)=0         !Air side heat transfer coefficients, kW/m^2-K !RS: Debugging: Never used
  !FTOUT(26)=0         !Inlet coil surface temperature, C !RS: Debugging: Never used
  !FTOUT(27)=0         !Outlet coil surface temperature, C    !RS: Debugging: Never used
  CondOUT%COutWtAl=MCOUT(19) !Aluminum weight, kg    !RS: Debugging: Formerly FTOUT(28), FTOUT(19)
  CondOUT%COutWtCu=0         !Copper weight, kg  !RS: Debugging: Formerly FTOUT(29), FTOUT(20)

  RETURN

END SUBROUTINE LoadMicrochannelOutputs

!************************************************************************

SUBROUTINE UpdateTubeDataFromCircuitData(I,J)

  !-----------------------------------------------------------------------------------
  !
  !  Description:
  !  Update the Tube data from Circuit "Ckt" data.
  !
  !  Author:
  !  Ipseng Iu
  !  Mechanical and Aerospace Engineering
  !  Oklahoma State University, Stillwater
  !
  !  Date: May 2006
  !
  !-----------------------------------------------------------------------------------

  IMPLICIT NONE

  !Subroutine passing variables
  INTEGER, INTENT(IN) :: I !Circuit number
  INTEGER, INTENT(IN) :: J !Tube number in calculation sequence

  !Subroutine local variables
  INTEGER TubeNum !Tube number in circuit diagram

  !FLOW:

  TubeNum=Ckt(I)%TubeSequence(J)
  Tube(TubeNum)=Ckt(I)%Tube(J)

  RETURN

END SUBROUTINE UpdateTubeDataFromCircuitData

SUBROUTINE GetCondProp(Out1, Out2, Out3, Out4, Out5, Out6, Out7, Out8)
  !RS: Integration: Trying to carry over the properties to output
  REAL Out1, Out2, Out3, Out4, Out5, Out6, Out7, Out8

  Out1=pRiCoil
  Out2=tRiCoil
  Out3=hRiCoil
  Out4=pRoCoil
  Out5=tRoCoil
  Out6=hRoCoil
  Out7=tAoCoil
  Out8=rhAoCoil

END SUBROUTINE

SUBROUTINE SimSimpleFan(tAiCoil,hAiCoil,MassFlow,tAoCoil,hAoCoil) !RS: Debugging: Adding in condenser fan model (9/4/14)

  ! SUBROUTINE INFORMATION:
  !       AUTHOR         Unknown
  !       DATE WRITTEN   Unknown
  !       MODIFIED       Brent Griffith, May 2009, added EMS override
  !                      Chandan Sharma, March 2011, FSEC: Added LocalTurnFansOn and LocalTurnFansOff
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS SUBROUTINE:
  ! This subroutine simulates the simple constant volume fan.

  ! METHODOLOGY EMPLOYED:
  ! Converts design pressure rise and efficiency into fan power and temperature rise
  ! Constant fan pressure rise is assumed.

  ! REFERENCES:
  ! ASHRAE HVAC 2 Toolkit, page 2-3 (FANSIM)

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! SUBROUTINE ARGUMENT DEFINITIONS:
  !Integer, Intent(IN) :: FanNum
  REAL, INTENT(IN) :: tAiCoil
  REAL, INTENT(IN) :: hAiCoil
  REAL, INTENT(IN) :: MassFlow !mAiCoil
  REAL, INTENT(OUT) :: tAoCoil
  REAL, INTENT(OUT) :: hAoCoil

  ! SUBROUTINE PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! SUBROUTINE LOCAL VARIABLE DECLARATIONS:
  !REAL(r64) RhoAir
  !REAL(r64) DeltaPress  ! [N/m2]
  !REAL(r64) FanEff
  !REAL(r64) MotInAirFrac
  !REAL(r64) MotEff
  !REAL(r64) MassFlow    ! [kg/sec]
  !unused0909      REAL(r64) Tin         ! [C]
  !unused0909      REAL(r64) Win
  REAL FanShaftPower ! power delivered to fan shaft
  REAL PowerLossToAir ! fan and motor loss to air stream (watts)
  !Integer NVPerfNum

  !NVPerfNum  = Fan(FanNum)%NVPerfNum

  !IF (NightVentOn .AND. NVPerfNum > 0) THEN
  !  DeltaPress = NightVentPerf(NVPerfNum)%DeltaPress
  !  FanEff = NightVentPerf(NVPerfNum)%FanEff
  !  MotEff = NightVentPerf(NVPerfNum)%MotEff
  !  MotInAirFrac = NightVentPerf(NVPerfNum)%MotInAirFrac
  !ELSE
  !  DeltaPress = Fan(FanNum)%DeltaPress
  !  FanEff     = Fan(FanNum)%FanEff
  !  MotEff     = Fan(FanNum)%MotEff
  !  MotInAirFrac = Fan(FanNum)%MotInAirFrac
  !END IF

  !IF (Fan(FanNum)%EMSFanPressureOverrideOn) DeltaPress = Fan(FanNum)%EMSFanPressureValue
  !IF (Fan(FanNum)%EMSFanEffOverrideOn) FanEff = Fan(FanNum)%EMSFanEffValue

  ! For a Constant Volume Simple Fan the Max Flow Rate is the Flow Rate for the fan
  !unused0909   Tin        = Fan(FanNum)%InletAirTemp
  !unused0909   Win        = Fan(FanNum)%InletAirHumRat
  !RhoAir     = Fan(FanNum)%RhoAirStdInit
  !MassFlow   = Fan(FanNum)%InletAirMassFlowRate
  !IF (Fan(FanNum)%EMSMaxMassFlowOverrideOn) MassFlow = Fan(FanNum)%EMSAirMassFlowValue
  ! MassFlow   = MIN(MassFlow,Fan(FanNum)%MaxAirMassFlowRate)
  ! MassFlow   = MAX(MassFlow,Fan(FanNum)%MinAirMassFlowRate)
  !
  !Determine the Fan Schedule for the Time step
  !If( ( GetCurrentScheduleValue(Fan(FanNum)%SchedPtr)>0.0 .or. LocalTurnFansOn) &
  !      .and. .NOT.LocalTurnFansOff  .and. Massflow>0.0) Then
  !Fan is operating
  !Fan(FanNum)%FanPower = MassFlow*DeltaPress/(FanEff*RhoAir) ! total fan power
  FanOut%HumRat=PsyWFnTdbHLocal(tAiCoil,(hAiCoil*1000))
  FanOut%RhoAir=PsyRhoAirFnPbTdbWLocal(CondPAR%CondBarPress,tAiCoil,FanOut%HumRat)
  FanOut%Power=MassFlow*FanOut%DeltaPress/(FanOut%FanEff*(FanOut%RhoAir*1000)) !RhoAir) !RS: Debugging: I think it needs to be *1000
  FanShaftPower = FanOut%MotorEff * FanOut%Power !Fan(FanNum)%FanPower  ! power delivered to shaft
  PowerLossToAir = FanShaftPower + (FanOut%Power - FanShaftPower) * FanOut%MotInAirFrac
  !Fan(FanNum)%OutletAirEnthalpy = Fan(FanNum)%InletAirEnthalpy + PowerLossToAir/MassFlow
  hAoCoil = hAiCoil +(PowerLossToAir/MassFlow)/1000
  ! This fan does not change the moisture or Mass Flow across the component
  !Fan(FanNum)%OutletAirHumRat       = Fan(FanNum)%InletAirHumRat
  !Fan(FanNum)%OutletAirMassFlowRate = MassFlow
  !Fan(FanNum)%OutletAirTemp = PsyTdbFnHW(Fan(FanNum)%OutletAirEnthalpy,Fan(FanNum)%OutletAirHumRat)
  tAoCoil = PsyTdbFnHWLocal((hAoCoil*1000),FanOut%HumRat) !Fan(FanNum)%OutletAirHumRat)

  !Else
  !  !Fan is off and not operating no power consumed and mass flow rate.
  !  Fan(FanNum)%FanPower = 0.0
  !  FanShaftPower = 0.0
  !  PowerLossToAir = 0.0
  !  Fan(FanNum)%OutletAirMassFlowRate = 0.0
  !  Fan(FanNum)%OutletAirHumRat       = Fan(FanNum)%InletAirHumRat
  !  Fan(FanNum)%OutletAirEnthalpy     = Fan(FanNum)%InletAirEnthalpy
  !  Fan(FanNum)%OutletAirTemp = Fan(FanNum)%InletAirTemp
  !  ! Set the Control Flow variables to 0.0 flow when OFF.
  !  Fan(FanNum)%MassFlowRateMaxAvail = 0.0
  !  Fan(FanNum)%MassFlowRateMinAvail = 0.0
  !
  !End If

  RETURN
END SUBROUTINE SimSimpleFan

FUNCTION PsyTdbFnHWLocal(H,dW) RESULT(TDB)  !RS: Debugging: Adding in condenser fan model (9/4/14)

  ! FUNCTION INFORMATION:
  !       AUTHOR         J. C. VanderZee
  !       DATE WRITTEN   Feb. 1994
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function provides air temperature from enthalpy and humidity ratio.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! ASHRAE HANDBOOK OF FUNDAMENTALS, 1972, P100, EQN 32
  !   by inverting function PsyHFnTdbW

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, intent(in) :: H    ! enthalpy {J/kg}
  REAL, intent(in) :: dW    ! humidity ratio
  !character(len=*), intent(in), optional :: calledfrom  ! routine this function was called from (error messages) !unused1208
  REAL        :: TDB  ! result=> dry-bulb temperature {C}

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  REAL W   ! humidity ratio

  W=MAX(dW,1.0d-5)
  TDB = (H - 2.50094d6 * W)/(1.00484d3 + 1.85895d3*W)

  RETURN
END FUNCTION PsyTdbFnHWLocal

FUNCTION PsyRhoAirFnPbTdbWLocal(pb,tdb,dw)  RESULT(rhoair)   !RS: Debugging: Adding in condenser fan model (9/4/14)

  ! FUNCTION INFORMATION:
  !       AUTHOR         G. S. Wright
  !       DATE WRITTEN   June 2, 1994
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function provides density of air as a function of barometric
  ! pressure, dry bulb temperature, and humidity ratio.

  ! METHODOLOGY EMPLOYED:
  ! ideal gas law
  !    universal gas const for air 287 J/(kg K)
  !    air/water molecular mass ratio 28.9645/18.01534

  ! REFERENCES:
  ! Wylan & Sontag, Fundamentals of Classical Thermodynamics.
  ! ASHRAE handbook 1985 Fundamentals, Ch. 6, eqn. (6),(26)

  ! USE STATEMENTS:
  ! na

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, INTENT(IN)  :: pb     ! barometric pressure (Pascals)
  REAL, INTENT(IN)  :: tdb    ! dry bulb temperature (Celsius)
  REAL, INTENT(IN)  :: dw      ! humidity ratio (kgWater/kgDryAir)
  !character(len=*), intent(in), optional :: calledfrom  ! routine this function was called from (error messages) !unused1208
  REAL        :: rhoair ! result=> density of air

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  REAL w  ! humidity ratio

  w=MAX(dw,1.0d-5)
  rhoair = pb/(287.d0*(tdb+273.15)*(1.d0+1.6077687d0*w))

  return
end function PsyRhoAirFnPbTdbWLocal

FUNCTION PsyWFnTdbHLocal(TDB,H) RESULT(W)    !RS: Debugging: Adding in condenser fan model (9/4/14)

  ! FUNCTION INFORMATION:
  !       AUTHOR         George Shih
  !       DATE WRITTEN   May 1976
  !       MODIFIED       na
  !       RE-ENGINEERED  na

  ! PURPOSE OF THIS FUNCTION:
  ! This function provides the humidity ratio from dry-bulb temperature
  ! and enthalpy.

  ! METHODOLOGY EMPLOYED:
  ! na

  ! REFERENCES:
  ! ASHRAE HANDBOOK OF FUNDAMENTALS, 1972, P100, EQN 32

  ! USE STATEMENTS:

  IMPLICIT NONE    ! Enforce explicit typing of all variables in this routine

  ! FUNCTION ARGUMENT DEFINITIONS:
  REAL, INTENT(IN) :: TDB    ! dry-bulb temperature {C}
  REAL, INTENT(IN) :: H      ! enthalpy {J/kg}
  !character(len=*), intent(in), optional :: calledfrom  ! routine this function was called from (error messages)
  REAL        :: W      ! result=> humidity ratio
  INTEGER, PARAMETER :: iPSyWFnTdbH =6

  ! FUNCTION PARAMETER DEFINITIONS:
  ! na

  ! INTERFACE BLOCK SPECIFICATIONS
  ! na

  ! DERIVED TYPE DEFINITIONS
  ! na

  ! FUNCTION LOCAL VARIABLE DECLARATIONS:
  ! na

  !CP-------- here is 1.2, 1200., 1.004, or 1004.  --------
  W=(H-1.00484d3*TDB)/(2.50094d6+1.85895d3*TDB)
  !
#ifdef EP_psych_stats
  NumTimesCalled(iPsyWFnTdbH)=NumTimesCalled(iPsyWFnTdbH)+1
#endif

  !                                      VALIDITY TEST.
  IF (W < 0.0d0) THEN
#ifdef EP_psych_errors
    IF (W < -.0001d0) THEN
      IF (.not. WarmupFlag) THEN
        IF (iPsyErrIndex(iPsyWFnTdbH) == 0) THEN
          !String=' Dry-Bulb= '//TRIM(TrimSigDigits(TDB,2))//' Enthalpy= '//TRIM(TrimSigDigits(H,3))
          !CALL ShowWarningMessage('Calculated Humidity Ratio invalid (PsyWFnTdbH)')
          if (present(calledfrom)) then
            !CALL ShowContinueErrorTimeStamp(' Routine='//trim(calledfrom)//',')
          else
            !CALL ShowContinueErrorTimeStamp(' Routine=Unknown,')
          endif
          !CALL ShowContinueError(TRIM(String))
          !String='Calculated Humidity Ratio= '//TRIM(TrimSigDigits(W,4))
          !CALL ShowContinueError(TRIM(String)//' ... Humidity Ratio set to .00001')
        ENDIF
        !CALL ShowRecurringWarningErrorAtEnd('Calculated Humidity Ratio invalid (PsyWFnTdbH)',   &
        !  iPsyErrIndex(iPsyWFnTdbH),ReportMinOf=W,ReportMaxOf=W,ReportMinUnits='[]',ReportMaxUnits='[]')
      ENDIF
    ENDIF
#endif
    W=1.d-5
  ENDIF

  ! W is the result

  RETURN
END FUNCTION PsyWFnTdbHLocal

!************************************************************************

END MODULE CondenserMod
