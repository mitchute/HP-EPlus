! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module simulates the behavior of the Evaporator component in the Heat pump.  
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
!This Component is the evaporator in the Heat Pump Vapor-Compression cycle
! A description of the component is found at:
! http://en.wikipedia.org/wiki/Evaporator
! From that website: 
!  - This component takes refrigerant from its liquid form to its natural gaseous form

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module takes the air properties entering the evaporator and using the properties of the evaporator 
! and ouputs the leaving air properties and heat transfer rates to the fluid.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! No apparent input or output files.
! The IDF file is used to read in some evaporator properties.
! 
! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There are a lot of variables defined at the module level.
! These include convergence parameters, error and location flags,
! circuiting variables, refrigerant property variables,
! coil and segment variables, and geometry properties.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 23 methods:
!    PUBLIC  Evaporator - To predict coil air side and refrigerant side properties, heat transfer, and pressure drop
!    PUBLIC  MicrochannelEvaporator -  A segment-by-segment microchannel evaporator model
!    PUBLIC  CalcEvaporatorInventory - calculate refrigerant inventory in evaporator
!    PUBLIC  PrintEvaporatorResult - To print simulation result to output file "Evaporator.csv"
!    PUBLIC  InitEvaporatorCoil - initialize evaporator geometry and circuiting
!    PUBLIC  EndEvaporatorCoil - free allocated arrays
!    PRIVATE RefrigerantParameters - calculates refrigerant parameters
!    PRIVATE SuctionLine - Calculates the condition at the suction line to the compressor
!    PRIVATE LoadMicrochannelInputs - Transfer input data from subroutine "Evaporator" to "MircochannelEvaporator"
!    PRIVATE LoadMicrochannelOutputs - Transfer output data from subroutine "MircochannelEvaporator" to "Evaporator"
!    PRIVATE InitBoundaryConditions - To initialize segment boundary conditions
!    PRIVATE CalcCircuitRefInletConditions - To calculate circuit refrigerant inlet condition according to circuitry
!    PRIVATE CalcSegmentRefInletConditions - To calculate inlet refrigerant pressure and enthalpy
!    PRIVATE CalcSegmentAirInletConditions - To calculate inlet air temp and relative humidity
!    PRIVATE CalcCoilSegment - To perform heat exchanger calculation for a segment
!    PRIVATE CalcSegmentOutletConditions - calculate segment outlet conditions
!    PRIVATE CalcTransitionSegment - To calculate transition segment (both single and two phase refrigerant
!                                    in segment) heat transfer
!    PRIVATE FindTransitionBoundary - To find the transition boundary in transition segment
!    PRIVATE CalcRefProperty - To calculate refrigerant properties given pressure and enthalpy
!    PRIVATE CalcSegmentRefOutletPressure - To calculate segment refrigerant outlet pressure
!    PRIVATE CalcWetSurfaceDing - To calculate wet surface air side heat transfer
!    PRIVATE CalcWetSurfaceBraun - To calculate wet surface air side heat transfer
!    PRIVATE CalcWetSurfaceMcQuiston - To calculate wet surface air side heat transfer
!    PRIVATE UpdateTubeDataFromCircuitData - Update the Tube data from Circuit "Ckt" data. 

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! Unit Conversions need to be changed

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-19 | JEH | Header Modification
! 2013-12-18 | RAS | Finished filling out the header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Some more documentation would be helpful. The ability of the code to switch from
! indoor circuit to outdoor and vice-versa needs to be explored; it currently seems
! to require manual changing of certain inputs and outputs for the switch.

MODULE EvaporatorMod

USE CoilCalcMod
USE DataSimulation
USE DataGlobals_HPSimIntegrated, ONLY: RefName, RefrigIndex    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)

IMPLICIT NONE

PRIVATE 

!Parameters
INTEGER,PARAMETER  :: MdotMaxIter=10      !Max. number of iterations
INTEGER,PARAMETER  :: RefBCmaxIter=20     !Max. number of iterations
INTEGER,PARAMETER  :: AirBCmaxIter=20     !Max. number of iterations
INTEGER,PARAMETER  :: PressureMaxIter=20          !Max. number of iterations
INTEGER,PARAMETER  :: WetSurfaceMaxIter=20          !Max. number of iterations
REAL,PARAMETER :: SMALL=1.0E-4  !Small number 
REAL,PARAMETER :: BIG=1.0E6     !Big number 

!Error Flags 
INTEGER,PARAMETER :: NOERROR       = 0
INTEGER,PARAMETER :: CONVERGEERROR = 1
INTEGER,PARAMETER :: REFPROPERROR  = 2
INTEGER,PARAMETER :: CKTFILEERROR  = 3
INTEGER,PARAMETER :: COILTUBEERROR = 4
INTEGER,PARAMETER :: COILFINERROR  = 5
INTEGER,PARAMETER :: AIRSIDEERROR  = 6
INTEGER,PARAMETER :: ZEROLENCOILERROR  = 7
INTEGER,PARAMETER :: CONVERGEERRORMINOR = 8

!System Types
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

!Variable for ref. mass flow rate distribution
REAL,ALLOCATABLE,DIMENSION(:),SAVE :: mRefIter !Circuit flow rate for iteration check, kg/s
INTEGER,ALLOCATABLE,DIMENSION(:),SAVE :: JoinTubes !Joined tube numbers

!Subcooling cirucits variables, ISI - 06/05/07
INTEGER :: SubcoolingTube !Subcooling tube number
INTEGER :: NumOfSubcoolingCkts !Number of subcooling circuits

!Circuitry variables
INTEGER I,J,K !Loop control
INTEGER NumOfTubes !Total number of tubes
INTEGER TubeNum    !Tube number in circuit diagram
INTEGER ErrorFlag          !0-No error
                           !1-Condenser solution not converge
						   !2-Refprop error
                           !3-Circuit file error
						   !4,5-Coil geometry misdefined
						   !6-Air side boundary condition not appropriate

!INTEGER CoilType           !1=Condenser; 2=Evaporator; 
                            !3=High side interconnecting pipes; 4=Low side interconnecting pipes

INTEGER FirstTube          !First simulation tube
INTEGER LastTube           !Last simulation tube
INTEGER EqCircuits         !1=Equivalent circuits; otherwise=no
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
REAL VelDevMod
REAL tAiMod,tAoMod,tAmod
REAL rhAiMod,rhAoMod
REAL wbAoMod
REAL hAiMod,hAoMod,hAmod    !RS: Debugging: Adding the last for !RS: Replace (2/20/14)

REAL DPmod !Pressure drop in module, kPa
REAL SigmaMod !Surface tension, N/m
REAL DTmod   !Temperature difference between saturated vapor and wall, C
REAL tRoCkt  !Circuit outlet temp. C
REAL pRiCkt  !Circuit Inlet pressure, kPa
REAL pRoCkt  !Circuit outlet pressure, kPa
REAL SumpRoCkt !Sum of outlet circuit outlet pressure, kPa
REAL hRoCkt  !Circuit outlet enthalpy, kJ/kg
REAL xRoCkt  !Circuit outlet quality

!Expansion device variables
REAL pRoExp
REAL hRoExp

!Compressor variables
REAL tRiCmp
REAL pRiCmp
REAL hRiCmp
REAL xRiCmp

!Heat transfer calc. variables
REAL mRefTot  !Refrigerant mass flow rate, kg/s

REAL mRefMod  !Module refrigerant mass flow rate, kg/s
REAL mRefCkt  !Ckt refrigerant mass flow rate, kg/s
REAL mRefJoin !Join tube total mass flow rate, kg/s
REAL Cmin     !Min. capacity rate, kW/C
REAL DT       !Temperature difference, C
REAL NTU      !Number of transfer unit
REAL Cratio   !Ratio of min to max capacity rate
REAL EPS      !Heat exchanger effectiveness
REAL Qcoil    !Total coil heat transfer, kW
REAL QcoilSens !Senible coil heat transfer, kW
REAL Qsection    !Total coil section heat transfer, kW !ISI - 09/10/07
REAL QsectionSens !Senible coil section heat transfer, kW !ISI - 09/10/07
REAL PrevQcoil !Previous value of total coil heat transfer, kW
REAL DiffQcoil !Difference of Qcoil in iteration, kW
REAL Qckt     !Circuit heat transfer, kW
REAL Qmod     !Module heat transfer, kW
REAL Qsurf    !Module heat transfer, kW
REAL QmodPrev !Previous module heat transfer for convergence calc., kW
REAL Tsurf    !Surface temperature, C
REAL TsurfNew !Surface temperature, C
REAL TsurfMax !Maximum surface temperature, C
REAL TsurfMin !Minimum surface temperature, C
REAL cAir     !Capacity rate of air, kW/C
REAL cRef     !Capacity rate of refrigerant, kW/C
REAL UA       !Overall heat transfer coefficient, kW/C
REAL hcRef    !Refrigerant film coefficent, W/m^2-C
REAL Rtube    !Thermal resistance of tube, K/W
REAL Rair     !Module air film resistance, K/W
REAL Rfrost   !Sankar
REAL Rrefrig  !Module refrigerant film resistance, K/W
REAL hco      !Air side heat tranfer coefficient, kW/m2-K
REAL hci      !Refrigerant side heat tranfer coefficient, kW/m2-K
REAL hdo      !Mass transfer coefficient, kg/m^2-s
REAL hcoMod   !Air side heat tranfer coefficient, kW/m2-K
REAL hcoDry   !Dry surface air side heat tranfer coefficient, kW/m2-K
REAL hciMod   !Refrigerant side heat tranfer coefficient, kW/m2-K
REAL EFref    !Refrigerant side heat tranfer enhancement factor
REAL Velavg   !Average face velocity, m/s
REAL ReVap    !Module Reynolds number vapor
REAL ReLiq    !Module Reynolds number liquid
REAL Const    !A constant
REAL MolWeight !Molecular weight, kg/kmol
REAL tSat     !Saturation temp., C
REAL QmodTP   !Heat transfer in two-phase region, kW 
REAL QdisTube !Distributor tube heat load, kW
REAL QmodSens !Sensible Module heat transfer, kW
REAL QmodLat  !Latent Module heat transfer, kW
REAL hciMultiplier  !Multiplier for hci
REAL hcoMultiplier  !Multiplier for hco
REAL DPrefMultiplier !Multiplier for DPref
REAL DPairMultiplier !Multiplier for DPair 
REAL AddDPSucLn !Suction line additional pressure drop, kPa
REAL DPfric !Frictional pressure drop, kPa
REAL DPgrav !Gravitational pressure drop, kPa
REAL DPmom !Momentum pressure drop, kPa
REAL FaceVel    !Face velocity, m/s
REAL DPair      !Air side pressure drop, kPa
REAL SurfAbsorptivity !Surface absorptivity
REAL SolarFlux  !Solar heat flux, kW/m2
REAL QlossCmp   !Compressor heat loss, kW
REAL IsCmpInAirStream !Is compressor in air stream, 1=yes, 0=no
!INTEGER(2) SystemType !1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat
INTEGER,SAVE :: CompManufacturer !Compressor manufacturer: 1=Copeland
                                                          !2=Bristol
														  !3=Danfoss
														  !4=Panasonic
!Properties
REAL mu       !Bulk viscosity, Pa-s
REAL muf      !Liquid viscosity, Pa-s
REAL mug      !Vapor viscosity, Pa-s
REAL kRef     !Refrigerant bulk conductivity, kW/m-C
REAL cpRef    !Ref. specific heat, kJ/(kg-C)
REAL rhoRef   !Ref. density, kg/m3
REAL CPAir                !Specific heat of air, kJ/kg-K
INTEGER(2) AirPropOpt     !Air prop calc. option
INTEGER(2) AirPropErr     !Error flag:1-error; 0-no error
!REAL AirProp(8)
REAL TwbAiMod !Inlet air wet bulb temperature, C
REAL TdpAiMod !Inlet air dewpoint temeperature, C
REAL TdbAoDry !Outlet air dry bulb temperature, dry surface, C
REAL TdbAoWet !Outlet air dry bulb temperature, wet surface, C
REAL TwbAoMod !Outlet air wet bulb temperature, C
REAL TwbAoWet !Outlet air wet bulb temperature, wet surface, C
REAL TdpAoMod !Outlet air dewpoint temperature, C
REAL wAiMod,wAoMod,wAmod
REAL hAoWet
REAL hAoDry
REAL hAiCoil,hAoCoil
REAL Cminf    !Fictitious min. capacity rate, kW/C
REAL EPSf     !Fictitious heat exchanger effectiveness
REAL EPSsDry  !Dry surface coil effectiveness
REAL EPSsWet  !Wet surface coil effectiveness
REAL QmodDry  !Dry module heat transfer, kW
REAL QmodWet  !Wet module heat transfer, kW
REAL NTUsWet  !Number of transfer unit for wet surface
REAL NTUsDry  !Number of transfer unit for dry surface
REAL TsDry    !Dry coil surface temperature, C
REAL TsWet    !Wet coil surface temperature, C
REAL HsWet    !Wet coil surface enthalpy, kJ/kg
REAL cf      !Fictitious air specific heat, kJ/kg-C
REAL cfprev  !Previous value of cf
REAL cAirf   !Fictitious capcity rate of air, kW/C
REAL UAf     !Fictitous overall heat transfer coefficient, kW/C
REAL Rairf   !Fictitous module air film resistance, K/W
REAL hcof    !Fictitious air side heat tranfer coefficient, kW/m2-K
REAL SHR     !Sensible heat ratio
REAL Rtot    !Total resistance, K-m^2/W
REAL Qsolar  !Solar radiation, kW

REAL QmodSensTot    !Total sensible heat transfer, kW
REAL QmodLatTot !Total latent heat transfer, kW
REAL QModTot    !Total heat transfer, kW

!Variables for coil
REAL mRiCoil
REAL mAiCoil
REAL tRiCoil,tRoCoil
REAL tSiCoil,tSoCoil
REAL pRiCoil,pRoCoil,pRoCoilPrev
REAL hRiCoil,hRoCoil
REAL xRiCoil,xRoCoil
REAL tAiCoil,tAoCoil
REAL rhAiCoil,rhAoCoil
REAL wAiCoil,wAoCoil
REAL wbAiCoil
REAL tSHoCoil     !Coil outlet superheat, C 
REAL tSHiCmp      !Compressor inlet superheat, C 
REAL Wabsolute    !Asolute oil mass fraction  
REAL DensityIn    !Inlet air density, kg/m3
REAL DensityOut    !Outlet air density, kg/m3

!Geometry variables
REAL Aface       !Coil face area
REAL,SAVE :: AiCoil !Inside coil surface area, m^2
REAL,SAVE :: AoCoil !Outside coil surface area, m^2
REAL,SAVE :: AfCoil !Coil fin surface area, m^2
REAL,SAVE :: AmCoil !Coil tube mean surface area, m^2
REAL AoMod       !Module outside surface area, m^2
REAL AfMod       !Module fin surface area
REAL AiMod       !Module inside surface area
REAL AiModSuc    !Module inside surface area for suction line
REAL AmMod       !Module tube mean surface area    
REAL Lcoil       !Total tube length, m
REAL LmodTube    !Module length of tube, m
REAL LmodSuc     !Module length of suction line, m
REAL,SAVE :: LmodTP     !Two-phase module length, m
REAL LmodTPmin   !Minimum two-phase module length, m
REAL LmodTPmax   !Maximum two-phase module length, m
REAL LmodTPratio !Ratio of two-phase modeul length to total module length
REAL LsucLn      !Suction line length, m
REAL ElevSucLn   !Suction line elevation, m
REAL IDsucLn     !Inside diameter of suction line, m
REAL ODsucLn     !Outside diameter of suction line, m 
REAL SucLnThk    !Suction line tube wall thickness, m
REAL DreturnBend !Return bend diameter, m
REAL LreturnBend !Return bend length, m
REAL DisLnThk    !Discharge line tube wall thickness, m
REAL LiqLnThk    !Liquid line tube wall thickness, m
REAL HtCoil      !Coil height, m
REAL FinSpg      !Fin spacing, m
REAL phi         !Parameter for fin efficiency calc.
REAL SurfEff     !Surface effectiveness
REAL FinEff      !Fin effectiveness

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
INTEGER,SAVE :: TubeOrientation !Tube orientation, 1=Horizontal; 2=Vertical
INTEGER,SAVE :: NumOfMods	  !Number of modules per tube 
INTEGER,SAVE :: NumOfChannels !Number cf channels
REAL,SAVE    :: Dchannel      !Channel diameter, m
INTEGER,SAVE :: NumOfCkts     !Number of circuits
INTEGER,SAVE :: ShiftTube     !1= last row lower than 2nd last row
                              !0= last row higher than 2nd last row

INTEGER NmodLast         !Total number of modules in the last row
INTEGER IsParallelSlabs !Parallel microchannel slabs (1=yes; 0=no)
INTEGER RowNum           !Coil row number
INTEGER Ntube            !Tube number !Loop counter
INTEGER Nckt             !Circuit number !Loop counter
INTEGER Nmod             !Module number !Loop counter
INTEGER NcktLast         !Total number of outlet circuits 
INTEGER NcktFirst        !Total number of inlet circuits 
INTEGER Nnode            !Number of split and joint nodes
INTEGER DryWet           !Dry wet flag: 1=Wet; 2=Partially wet; 3=Dry
LOGICAL,SAVE :: IsSameNumOfTubes !Flag to check if same number of tubes
                                 !in all circuit branches
REAL DrawBlow  !Fan location, 1=draw through; 2=blow through
REAL PwrFan	   !Fan power, kW
INTEGER WetFlag            !1=Wet; 0=dry
INTEGER Iter               !Iteration loop counter
INTEGER cfIter             !Iteration loop counter
INTEGER RefBCiter        !Iteration loop counter
INTEGER AirBCiter        !Iteration loop counter

INTEGER(2)       :: RefPropErr  !Error flag:1-error; 0-no error
REAL Psat,Pcr,Tcr

!Mass inventory
REAL MassSucLn    !Total refrigerant inventory in suction line, kg
REAL MassMod      !Refrigerant inventory in a module, kg
REAL MassLiqMod   !Mass in liquid phase, kg
REAL MassVapMod   !Mass in vapor phase, kg

REAL, SAVE :: WeightAluminum !Weight of aluminum, kg
REAL, SAVE :: WeightCopper   !Weight of copper, kg

REAL TestH    !RS: Debugging: Finding the entering air enthalpy hopefully

INTEGER FirstTime !Flag to indicate the first time of execution
                  !1=yes, otherwise=no
INTEGER Counter                   !Iteration loop counter
INTEGER IsSimpleCoil !Flag to indicate if it is simple coil, i.e. ignoring circuiting
                     !1=Simple coil
				     !otherwise=detailed
INTEGER NumOfSections !Number of sections, ISI - 09/10/07
INTEGER NumInletSections !Number of inlet sections, ISI - 09/10/07

TYPE (SlabInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Slab      !Coil slab pointer
TYPE (CktInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Ckt       !Circuit pointer
TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Tube     !Tube pointer
TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:,:),SAVE :: Tube2D !2-dimensional Tube pointer
TYPE (ModInfo),ALLOCATABLE,DIMENSION(:),SAVE :: SucLnSeg  !Suction line pointer
TYPE (NodeInfo),ALLOCATABLE,DIMENSION(:),SAVE :: Node     !Split or joint node
TYPE (SectionInfo),ALLOCATABLE,DIMENSION(:) :: CoilSection !Coil section, ISI - 09/10/07

PUBLIC  Evaporator
PUBLIC  MicrochannelEvaporator
PUBLIC  CalcEvaporatorInventory
PUBLIC  PrintEvaporatorResult
PUBLIC  InitEvaporatorCoil
PUBLIC  EndEvaporatorCoil
PRIVATE RefrigerantParameters
PRIVATE SuctionLine
PRIVATE LoadMicrochannelInputs
PRIVATE LoadMicrochannelOutputs
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
PRIVATE CalcWetSurfaceDing
PRIVATE CalcWetSurfaceBraun
PRIVATE CalcWetSurfaceMcQuiston
PRIVATE UpdateTubeDataFromCircuitData
PUBLIC GetQOut  !RS: Testing: Trying to integrate HPSim and EPlus
PUBLIC GetEvapProp   !RS: Integration: Trying to carry over the properties to output
PRIVATE InitEvaporatorStructures    !RS: Debugging
PUBLIC GetNodeProp  !RS: Debugging: Trying to update nodal properties

CONTAINS

!***********************************************************************************

    SUBROUTINE Evaporator(Ref$) !,XIN,PAR,OUT) !(Ref$,PureRef,XIN,PAR,OUT) !RS: Debugging: Extraneous PureRef

    !-----------------------------------------------------------------------------------
    !
    !  Description:	
    !  Ragazzi's modular coil model (Fixed length version)
    !  To predict coil air side and refrigerant side properties, heat transfer, 
    !  and pressure drop
    !
    !  HVACSIM+ Airprop subroutine required
    !  EnergyPlus REFPROP subroutines required 	
    !
    !  Inputs:
    !  Ref$=Refrigerant name
    !  PureRef=Refrigerant flag: 1=pure refrigerant
    !                            0=refrigerant mixture
    !  XIN(1)=Refrigerant side mass flow rate, kg/s
    !  XIN(2)=Refrigerant side inlet (distributor tube outlet) pressure, kPa
    !  XIN(3)=Refrigerant side inlet (distributor tube outlet) enthalpy, kJ/kg
    !  XIN(4)=Air side mass flow rate, kg/s
    !  XIN(5)=Air side inlet temp. C
    !  XIN(6)=Air side inlet relative humidity
    !  XIN(7)=Distributor tube heat load, kW
    !  XIN(8)=Solar heat flux, kW/m^2
    !  XIN(9)=Compressor Discharge Temp. C
    !
    !  Parameters:
    !  PAR(1)=Suction line length, m
    !  PAR(2)=Suction line outside diameter, m
    !  PAR(3)=Suction line tube wall thickness, m 
    !  PAR(4)=Suction line elevation, m
    !  PAR(5)=Suction line heat loss, kW
    !  PAR(6)=Suction line temperature change, C
    !  PAR(7)=Suction line additional pressure drop, kPa
    !  PAR(8)=Coil tube outside diameter, m
    !  PAR(9)=Coil tube wall thickness, m
    !  PAR(10)=Coil single tube length, m
    !  PAR(11)=Coil tube thermal conductivity, kW/m-C
    !  PAR(12)=Tube spacing in transverse direction, m (normal to air flow)
    !  PAR(13)=Row spacing in longitudinal direction, m (parallel to air flow)
    !  PAR(14)=Fin thickness, m
    !  PAR(15)=Fin pitch, fin/m
    !  PAR(16)=Fin thermal conductivity, kW/m-C
    !  PAR(17)=Number of tubes in transverse direction (normal to air flow)
    !  PAR(18)=Number of rows in longitudinal direction (parallel to air flow)
    !  PAR(19)=Number of circuits
    !  PAR(20)=Cooling mode? 1=yes; 0=no
    !  PAR(21)=Number of modules per tube
    !  PAR(22)=Fin type: 1=Plain; 2=Wavy; 3=Louver
    !  PAR(23)=Multiplier for ref. side heat transfer correlation
    !  PAR(24)=Multiplier for ref. side pressure drop correlation
    !  PAR(25)=Multiplier for air side heat transfer correlation
    !  PAR(26)=Multiplier for air side pressure drop correlation
    !  PAR(27)=Fan power, kW
    !  PAR(28)=Fan location, 1=draw through; 2=blow through
    !  PAR(29)=Surface absorptivity
    !  PAR(30)=Tube tube: 1=Smooth; 2-Microfin; 3=Herringbone; 4=Crosshatch; 
    !                     5=Herringbone w/crosshatch; 6-Turbo-A
    !  PAR(31)=Barometric pressure, kPa
    !  PAR(32)=Compressor heat loss, kW
    !  PAR(33)=System Type, 1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat
    !  PAR(34)=Pressure tolerance convergence Criteria, kPa
    !  PAR(35)=Compressor manufacturer: 1=Copeland; 2=Bristol; 
    !                                   3=Danfoss;  4=Panasonic
    !  PAR(36)=Oil mass fraction
    !  PAR(37)=Simple coil flag: 1=Simple coil; otherwise=Detailed coil
    !  PAR(38)=First time to run this model flag: 1=yes, otherwise=no
    !          for component validation, set it to 1
    !          for system validation, set it to 1 first, then zero
    !  PAR(39)= !RS: Debugging: Nothing
    !
    !  Outputs:
    !  OUT(1)=Coil outlet pressure, kPa
    !  OUT(2)=Coil outlet enthalpy, kJ/kg
    !  OUT(3)=Coil outlet temperature, C
    !  OUT(4)=Coil outlet quality
    !  OUT(5)=Coil outlet superheat, C
    !  OUT(6)=Suction line outlet pressure, kPa
    !  OUT(7)=Suction line outlet enthalpy, kJ/kg
    !  OUT(8)=Suction line outlet temperature, C
    !  OUT(9)=Suction line outlet quality
    !  OUT(10)=Suction line outlet superheat, C
    !  OUT(11)=Coil capacity, kW
    !  OUT(12)=Sensible coil capacity, kW
    !  OUT(13)=Mass in suction line, kg
    !  OUT(14)=Mass in coil, kg
    !  OUT(15)=Liquid mass in coil, kg
    !  OUT(16)=Vapor mass in coil, kg
    !  OUT(17)=Air side outlet temperature, C
    !  OUT(18)=Air side outlet relative humidity
    !  OUT(19)=Air side pressure drop, kPa
    !  OUT(20)=Error flag: 0-No error
    !                      1-Evaporator solution not converge
    !                      2-Refprop error
    !					   3-Circuit file error
    !  OUT(21)=Air side heat transfer coefficients, kW/m^2-K
    !  OUT(22)=Inlet coil surface temperature, C
    !  OUT(23)=Outlet coil surface temperature, C
    !  OUT(24)=Aluminum weight, kg
    !  OUT(25)=Copper weight, kg
    !
    !  References: 
    !  Ding, X., Eppe, J.P., Lebrun, J. and Wasacz, M. (1990). Cooling coil
    !    models to be used in transient and/or wet regimes. Theoretical analysis
    !    and experimental validation. In Proceedings of System Simulation in
    !    Buildings '90, Liege, Belgium, December 1990.
    !
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

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE CoilCalcMod
    USE AirPropMod
    USE OilMixtureMod
    USE ReversingValveMod
    USE InputProcessor   !RS: Debugging: GetObjectItem

    IMPLICIT NONE

    !Subroutine argument declarations
    CHARACTER*80,     INTENT(IN)  :: Ref$
    !INTEGER(2),       INTENT(IN)  :: PureRef   !RS: Debugging: Extraneous PureRef
    !REAL, INTENT(IN)  :: XIN(9) 
    !REAL, INTENT(IN)  :: PAR(39) !ISI - 12/21/06    !RS: Debugging: Formerly EvapPAR(54)
    !REAL, INTENT(OUT) :: OUT(17)    !RS: Debugging: Formerly OUT(20)

    !Subroutine local variables
    REAL :: MCXIN(7)  !Microchannel coil input data
    REAL :: MCPAR(33) !Microchannel coil input parameters
    REAL :: MCOUT(21) !Microchannel coil output data

    INTEGER,SAVE :: CoilType    !1=Condenser; 2=Evaporator; 
                                !3=High side interconnecting pipes; 
                                !4=Low side interconnecting pipes
                                !5=Microchannel condenser
                                !6=Microchannel evaporator

    REAL tSiSUM    !Sum of inlet surface temperature, C
    REAL tSoSUM    !Sum of outlet surface temperature, C
    REAL QcktSens !Sensible Circuit heat transfer, kW
    REAL DPair        !Air side pressure drop, kPa
    LOGICAL Converged   !Solution convergence flag
    REAL MaxResidual    !Maximum residual in iteration
    Real PTol           !Evaporator Outlet Pressure convergence criteria,kPa
    REAL tRdis          !Compressor discharge temperature, C
    REAL DPvalve        !Reversing valve pressure drop, kPa
    REAL hRsuc          !Suction enthalpy, kPa
    INTEGER NumSection  !Section number, ISI - 09/10/07
    REAL humrat !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio

    REAL :: tSAvgCoil = 0.0
    
    INTEGER, PARAMETER :: MaxNameLength = 200

    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL :: LmodTubeKeep    !RS: Debugging: Placeholder for LmodTube that keeps the previous one                                                                    
    
    CHARACTER(LEN=39),PARAMETER :: FMT_107 = "(A67,F5.6)" !VL Comment: previously !10

    !TestH=AirProp(4)    !RS: Debugging: Finding the entering air enthalpy hopefully
    
    !Flow:

    mRefTot =EvapIN%EInmRef !RS: Debugging: Formerly XIN(1)
    pRiCoil =EvapIN%EInpRi !RS: Debugging: Formerly XIN(2)
    hRiCoil =EvapIN%EInhRi !RS: Debugging: Formerly XIN(3)
    mAiCoil =EvapIN%EInmAi !RS: Debugging: Formerly XIN(4)
    tAiCoil =EvapIN%EIntAi !RS: Debugging: Formerly XIN(5)
    rhAiCoil=EvapIN%EInrhAi !RS: Debugging: Formerly XIN(6)
    SolarFlux=EvapIN%EInSolFlux    !RS: Debugging: EvapIN(8) set once to 0 !RS: Debugging: Formerly XIN(8)
    tRdis=EvapIN%EIntRdis    !RS: Debugging: Formerly XIN(9)

    !IF (QdisTube .EQ. 0) THEN  RS: Debugging: QDisTube is never used
    !    QdisTube=SMALL
    !END IF

    LsucLn    = EvapPAR%EvapSucLnLen !RS: Debugging: Formerly PAR(1)
    ODsucLn   = EvapPAR%EvapSucLnOD  !RS: Debugging: Formerly PAR(2)
    SucLnThk  = EvapPAR%EvapSucLnTWThick   !RS: Debugging: Formerly PAR(3)
    ElevSucLn = EvapPAR%EvapSucLnElev  !RS: Debugging: Formerly PAR(4)
    QsucLn    = EvapPAR%EvapSucLnQLoss  !RS: Debugging: Formerly PAR(5)
    DTsucLn   = EvapPAR%EvapSucLnTempChg  !RS: Debugging: Formerly PAR(6)
    AddDPSucLn = EvapPAR%EvapSucLnAddPD !RS: Debugging: Formerly PAR(7)

    IsSimpleCoil=EvapPAR%EvapSimpCoil !(38) !PAR(53) !ISI - 12/21/06    !RS: Debugging: Formerly PAR(37) !RS: Debugging: IsSimple
    FirstTime=EvapPAR%EvapFirstTime !(39) !PAR(54)    !ISI - 12/21/06   !RS: Debugging: Formerly PAR(38)

    IF (FirstTime .EQ. 1) THEN

        ODtube      = EvapPAR%EvapCoilTOD    !RS: Debugging: Formerly PAR(8)
        TubeThk     = EvapPAR%EvapCoilTWThick    !RS: Debugging: Formerly PAR(9)
        Ltube       = EvapPAR%EvapCoilSTLen   !RS: Debugging: Formerly PAR(10)
        Ktube       = EvapPAR%EvapCoilTThermCon   !RS: Debugging: Formerly PAR(11)
        Pt          = EvapPAR%EvapTspc   !RS: Debugging: Formerly PAR(12)
        Pl          = EvapPAR%EvapRspc   !RS: Debugging: Formerly PAR(13)
        FinThk      = EvapPAR%EvapFinThick   !RS: Debugging: Formerly PAR(14)
        FinPitch    = EvapPAR%EvapFinPitch   !RS: Debugging: Formerly PAR(15)
        Kfin        = EvapPAR%EvapFinThermCon   !RS: Debugging: Formerly PAR(16)
        Nt          = EvapPAR%EvapNt   !RS: Debugging: Formerly PAR(17)
        Nl          = EvapPAR%EvapNl   !RS: Debugging: Formerly PAR(18)
        NumOfCkts   = EvapPAR%EvapNumCkt   !RS: Debugging: Formerly PAR(19)
        NumOfMods   = EvapPAR%EvapNumMod   !RS: Debugging: Formerly PAR(21)
        FinType     = EvapPAR%EvapFinType   !RS: Debugging: Formerly PAR(22)
        TubeType    = EvapPAR%EvapTube   !RS: Debugging: Formerly PAR(30)
        CALL InitEvaporatorCoil(CoilType)
        CALL CalcMaterialWeight(CoilType,Ltube,IDtube,ODtube,TubeHeight,TubeDepth, &
        Dchannel,NumOfChannels,Pt,Pl,Nt,Nl,NumOfCkts, &
        FinThk,FinPitch,WeightAluminum,WeightCopper)
        IF (ErrorFlag .GT. CONVERGEERROR) THEN
            EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
            CALL Evaporator_Helper_1
            RETURN
        END IF
        CALL RefrigerantParameters(Ref$)
        CALL GetRefID(Ref$,RefID)
        
          !********************Refrigerant Cycle Data (Heating)***********************  !RS: Debugging: Moving: Stay here? Compressor? ORNLSolver?

        CALL GetObjectItem('RefrigerantCycleData(Heating)',1,Alphas,NumAlphas, &
                            TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

        IsCmpInAirStream = Numbers(2) !Is Compressor in Air Stream
 
    END IF

    hciMultiplier   = EvapPAR%EvapMultRefQT   !RS: Debugging: Formerly PAR(23)
    DPrefMultiplier = EvapPAR%EvapMultRefPD   !RS: Debugging: Formerly PAR(24)
    hcoMultiplier   = EvapPAR%EvapMultAirQT   !RS: Debugging: Formerly PAR(25)
    DPairMultiplier = EvapPAR%EvapMultAirPD   !RS: Debugging: Formerly PAR(26)

    PwrFan           = EvapPAR%EvapFanPwr/1000  !RS: Debugging: Formerly PAR(27) Conversion
    DrawBlow         = EvapPAR%EvapFanLoc  !RS: Debugging: Formerly PAR(28)
    SurfAbsorptivity = EvapPAR%EvapSurfAbs  !RS: Debugging: Formerly PAR(29)

    BaroPressure     = EvapPAR%EvapBarPress  !RS: Debugging: Formerly PAR(31)
    QlossCmp         = EvapPAR%EvapCompQLoss  !RS: Debugging: Formerly PAR(32)
    SystemType       = EvapPAR%EvapSysType  !RS: Debugging: Formerly PAR(33)
    PTol             = EvapPAR%EvapPressTolConv !(35) !(50)    !RS: Debugging: Formerly PAR(34)
    Wabsolute        = EvapPAR%EvapOilMassFrac !(51)    !RS: Debugging: Formerly PAR(35)
    CompManufacturer = EvapPAR%EvapCompMan !(52)    !RS: Debugging: Formerly PAR(36)

    IF (CoilType .EQ. MCEVAPORATOR) THEN
        CALL LoadMicrochannelInputs(MCXIN,MCPAR) !(XIN,PAR,MCXIN,MCPAR)
        CALL MicrochannelEvaporator(MCXIN,MCPAR,MCOUT) !(Ref$,MCXIN,MCPAR,MCOUT)    !RS: Debugging: Extraneous Ref$
        CALL LoadMicrochannelOutputs(MCOUT) !,OUT)
        RETURN
    END IF

    EqCircuits=0 !Equivalent circuit flag

    ErrorFlag=NOERROR !Initialize

    !Tube inside diameter
    IDsucLn=ODsucLn-SucLnThk*2

    !Fin spacing
    FinSpg = 1/FinPitch-FinThk

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

    CALL InitBoundaryConditions(CoilType)
    IF (ErrorFlag .GT. CONVERGEERROR) THEN
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF

    !****** Start coil calculation ******
    Qcoil=0.0; QcoilSens=0 !ISI - 09/10/07
    DO NumSection=1, NumOfSections !ISI - 09/10/07
    
        !ISI - 09/10/07
        Converged=.TRUE. 
        pRoCoilPrev=pRiCoil
        MaxResidual=0
    
        DO Iter=1, MdotMaxIter
    
            !Initialize
            mRefJoin=0; PrevQcoil=BIG; 
    
            DO AirBCiter=1, AirBCmaxIter
    
                !ISI - 09/10/07
                CoilSection(NumSection)%Qsection=0.0
                CoilSection(NumSection)%QsectionSens=0.0
                tSiSUM=0.0; tSoSUM=0.0
    
                DO I=1,CoilSection(NumSection)%NumOfCkts
    
                    Qckt=0.0; QcktSens=0
    
                    CALL CalcCircuitRefInletConditions(NumSection,I,I,CoilType) !ISI - 09/10/07
    
                    !Find first and last simulation tubes
                    IF (IsSimpleCoil .EQ. 1) THEN
                        FirstTube=1
                        LastTube=1
                    ELSE !ISI - 09/10/07
                        FirstTube=1
                        LastTube=CoilSection(NumSection)%Ckt(I)%Ntube 
                        IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
                            FirstTube=2 !Skip first tube
                        END IF 
                        IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
                            LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 !Ignore last tube
                        END IF
                    END IF
    
                    !****************** Hard code refrigerant distribution ****************** 
                    !              CoilSection(NumSection)%Ckt(1)%mRef=mRefTot*1.5/8
                    !              CoilSection(NumSection)%Ckt(2)%mRef=mRefTot*1/8
                    !              CoilSection(NumSection)%Ckt(3)%mRef=mRefTot*1.5/8
                    !              CoilSection(NumSection)%Ckt(4)%mRef=mRefTot*1.5/8
                    !              CoilSection(NumSection)%Ckt(5)%mRef=mRefTot*1/8
                    !              CoilSection(NumSection)%Ckt(6)%mRef=mRefTot*1.5/8
                    !************************************************************************
    
                    mRefMod=CoilSection(NumSection)%Ckt(I)%mRef !ISI - 09/10/07
    
                    DO J=FirstTube,LastTube
    
                        IF (IsSimpleCoil .EQ. 1) THEN
                            TubeNum=1
                        ELSE
                            TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
                        END IF
    
                        DO K=1,NumOfMods
    
                            IF (IsSimpleCoil .EQ. 1) THEN
                                !SELECT CASE(K)
                                !CASE (1)
                                !    LmodTube=Lcoil/CoilSection(NumSection)%NumOfCkts !Start with guessing the whole length
                                !CASE (2)
                                !    LmodTube=Lcoil/CoilSection(NumSection)%NumOfCkts-CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(1)%Len
                                !END SELECT
                                IF (K .EQ. 1) THEN  !RS: Debugging: Handling every module
                                    LmodTube=Lcoil/CoilSection(NumSection)%NumOfCkts
                                    LmodTubeKeep=LmodTube
                                ELSE
                                    LmodTube=LmodTubeKeep/NumOfMods 
                                !LmodTube=LModTube-CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(1)%Len/NumOfMods !RS: Debugging: Adding
        !RS: con. "NumOfMods" to account for more than one segment (2/4/14)
                                END IF
                            END IF
                            CALL CalcCoilSegment(NumSection,I,I,J,K,CoilType)  !RS: Debugging: Temporarily setting in an Epsilon-NTU method
                        
                            IF (ErrorFlag .GT. CONVERGEERROR) THEN
                                EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
                                CALL Evaporator_Helper_1
                                RETURN
                            END IF
    
                            !Calc. circuit heat transfer
                            Qckt=Qckt+CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Qmod
                            QcktSens=QcktSens+CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%QmodSens
    
                            !Calc. sum of outlet air surface temperature
                            IF (CoilSection(NumSection)%Ckt(I)%Tube(J)%Back .EQ. 0) THEN 
                                tSoSUM=tSoSUM+CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tSo
                            END IF
    
                            IF (CoilSection(NumSection)%Ckt(I)%Tube(J)%Fup .EQ. 0 .AND. &
                            CoilSection(NumSection)%Ckt(I)%Tube(J)%Fdown .EQ. 0) THEN
                                tSiSUM=tSiSUM+CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tSi
                            END IF
    
                        END DO !End mod
    
                    END DO !End tube
    
                    pRoCkt=CoilSection(NumSection)%Ckt(I)%Tube(LastTube)%Seg(NumOfMods)%pRo !Circuit outlet pressure
                    hRoCkt=CoilSection(NumSection)%Ckt(I)%Tube(LastTube)%Seg(NumOfMods)%hRo !Circuit outlet enthalpy
    
                    Pressure=pRoCkt*1000    !RS Comment: Unit Conversion
                    Enthalpy=hRoCkt*1000    !RS Comment: Unit Conversion
                    tRoCkt=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !RS Comment: Circuit Outlet Refrigerant Temperature
                    IF (RefPropErr .GT. 0) THEN
                        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2205'
                        ErrorFlag=REFPROPERROR
                        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
                        CALL Evaporator_Helper_1
                        RETURN
                    END IF
                    xRoCkt=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !RS Comment: Circuit Outlet Refrigerant Quality
                    IF (RefPropErr .GT. 0) THEN
                        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2211'
                        ErrorFlag=REFPROPERROR
                        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
                        CALL Evaporator_Helper_1
                        RETURN
                    END IF
    
                    Pressure=pRoCkt*1000    !RS Comment: Unit Conversion
                    Quality=1
                    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Saturation Temperature
                    IF (RefPropErr .GT. 0) THEN
                        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2231'
                        ErrorFlag=REFPROPERROR
                        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
                        CALL Evaporator_Helper_1
                        RETURN
                    END IF
    
                    IF (xRoCkt .GE. 1) THEN 
                        CoilSection(NumSection)%Ckt(I)%tSH=tRoCkt-tSat !Superheat
                    ELSE
                        CoilSection(NumSection)%Ckt(I)%tSH=0.0
                    END IF
    
                    CoilSection(NumSection)%Ckt(I)%tRo=tRoCkt
                    CoilSection(NumSection)%Ckt(I)%pRo=pRoCkt
                    CoilSection(NumSection)%Ckt(I)%hRo=hRoCkt
                    CoilSection(NumSection)%Ckt(I)%Qckt=Qckt          !Circuit capacity
                    CoilSection(NumSection)%Ckt(I)%QcktSens=QcktSens  !Sensible Circuit capacity
                    CoilSection(NumSection)%Qsection=CoilSection(NumSection)%Qsection+CoilSection(NumSection)%Ckt(I)%Qckt !Total coil capacity
                    CoilSection(NumSection)%QsectionSens=CoilSection(NumSection)%QsectionSens+CoilSection(NumSection)%Ckt(I)%QcktSens !Total sensible coil capacity
    
                    IF (EqCircuits .EQ. 1 .AND. IsUniformVelProfile .OR. IsSimpleCoil .EQ. 1) THEN  !Equivalent circuit and Uniform velocity profile
                        Qcoil=CoilSection(NumSection)%Qsection*CoilSection(NumSection)%NumOfCkts
                        QcoilSens=CoilSection(NumSection)%QsectionSens*CoilSection(NumSection)%NumOfCkts
                        pRoCoil=pRoCkt
                        hRoCoil=hRoCkt
                        CoilSection(NumOfSections)%pRo=pRoCkt
                        CoilSection(NumOfSections)%hRo=hRoCkt
                        EXIT
                    END IF
    
                    IF (CoilSection(NumSection)%Ckt(I)%OutSplit .LE. 1 .AND. CoilSection(NumSection)%Ckt(I)%OutJoin .LE. 1) THEN
                        mRefJoin=mRefJoin+CoilSection(NumSection)%Ckt(I)%mRef 
                    END IF
    
                END DO !End circuit
    
                IF (IsSimpleCoil .EQ. 1) THEN
                    EXIT
                END IF
    
                CALL CalcMeanProp(tAiCoil,tAoCoil,tAmod)    !RS Comment: Mean Air Coil Temperature
                CALL CalcMeanProp(hAiCoil,hAoCoil,hAmod)    !RS Comment: Mean Air Coil Temperature
                !CPair=CPA(REAL(tAmod))  !RS Comment: Finding the specific heat of air   !RS: Replace: CPA (2/19/14)
                humrat=HUMTH(tAmod,hAmod)   !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
                CPair=CPAirFunction(tAmod,humrat)  !RS: Replace: CPA (2/19/14)
                Cair=mAiCoil*CPAir  !RS Comment: Finding the capacity rate of air
    
                tAoCoil=tAiCoil+CoilSection(NumSection)%QsectionSens/Cair   !RS Comment: Air Coil Outlet Temperature
    
                DiffQcoil=ABS((CoilSection(NumSection)%Qsection-PrevQcoil)/PrevQcoil) 
                IF (DiffQcoil .GT. SMALL) THEN
                    PrevQcoil=CoilSection(NumSection)%Qsection
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
                EXIT !for equivalent circuit, no need to update mdot ref.
            END IF
    
            !Synchronize from circuits to tubes
            DO I=1, CoilSection(NumSection)%NumOfCkts
    
                FirstTube=1
                LastTube=CoilSection(NumSection)%Ckt(I)%Ntube
                IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
                    FirstTube=2 !Skip first tube
                END IF 
                IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
                    LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 !Ignore last tube
                END IF
                IF (FirstTube .GT. LastTube) THEN
                    LastTube=FirstTube !Dummy tube
                END IF
    
                DO J=FirstTube, LastTube
                    TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
                    Tube(TubeNum)=CoilSection(NumSection)%Ckt(I)%Tube(J)
                END DO
    
                IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
                    CoilSection(NumSection)%Ckt(I)%Tube(1)=Tube(CoilSection(NumSection)%Ckt(I)%TubeSequence(1))
                END IF 
                IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
                    CoilSection(NumSection)%Ckt(I)%Tube(CoilSection(NumSection)%Ckt(I)%Ntube)=Tube(CoilSection(NumSection)%Ckt(I)%TubeSequence(Ckt(I)%Ntube))
                END IF
    
            END DO !End circuit
    
            pRoCkt=0
            SumpRoCkt=0
            DO I=1, CoilSection(NumSection)%NumOfCkts
                IF (CoilSection(NumSection)%Ckt(I)%OutSplit .EQ. 0) THEN !outlet circuit
                    SumpRoCkt=SumpRoCkt+CoilSection(NumSection)%Ckt(I)%pRo  !RS Comment: Summing the outlet circuit refrigerant pressure
                END IF
            END DO !End Circuit
            CoilSection(NumSection)%pRo=SumpRoCkt/CoilSection(NumSection)%NcktLast
            IF (SumpRoCkt .EQ. 0) THEN
                SumpRoCkt=CoilSection(NumSection)%Ckt(1)%pRo !At least 1 circuit, ISI - 07/28/06
            END IF
    
            IF (ABS(CoilSection(NumSection)%pRo-pRoCoilPrev) .GT. Ptol) THEN
                MaxResidual=ABS(pRoCoilPrev-CoilSection(NumSection)%pRo)
                pRoCoilPrev=CoilSection(NumSection)%pRo
                Converged=.FALSE.   !No convergence
            END IF
    
            IF (IsSameNumOfTubes .AND. IsUniformVelProfile) THEN !RS: Debugging: Does this do anything?
                EXIT
            END IF
    
            IF (NOT(Converged) .OR. Iter .LE. 2) THEN 
                Converged=.TRUE. ! Reinitialize
    
                !Moved this subroutine to CoilCalc and share with condenser ISI - 06/05/07
                CALL UpdateRefMassFlowRate(Iter,CoilSection(NumSection)%Ckt, &
                CoilSection(NumSection)%NumOfCkts, &
                CoilSection(NumSection)%pRi,mRefTot, &
                CoilSection(NumSection)%Nnode,CoilSection(NumSection)%Node)
            ELSE
                EXIT
            END IF
    
        END DO !End iter
    
        IF (Iter .GT. MdotMaxIter) THEN
            IF (MaxResidual .GT. Ptol) THEN
                WRITE(*,FMT_107)'-- WARNING -- Evaporator: Solution not converged. Max. Residual = ',MaxResidual
                ErrorFlag=CONVERGEERROR
            END IF
        END IF
    
        IF (IsSimpleCoil .EQ. 1) THEN
            EXIT
        END IF
    
        !ISI - 09/10/07
        CoilSection(NumSection)%hRo=CoilSection(NumSection)%hRi- &
        CoilSection(NumSection)%Qsection/CoilSection(NumSection)%mRef
    
        Qcoil=Qcoil+CoilSection(NumSection)%Qsection                !RS Comment: Determining the total coil heat transfer
        Qcoilsens=QcoilSens+CoilSection(NumSection)%QsectionSens    !RS Comment: Determing the total sensible coil heat transfer
    
    END DO !Number of sections, !ISI - 09/10/07

    pRoCoil=CoilSection(NumOfSections)%pRo

    !Surface temperature
    tSiCoil=tSiSUM/NmodLast
    tSoCoil=tSoSUM/NmodLast

    tSAvgCoil = (tSiCoil+tSoCoil)/2
    CoilParams(2)%TSurfCoil=tSAvgCoil
    !Coil air side outlet condition
    !CPair=CPA(REAL(tAmod))  !RS Comment: Finding the specific heat of air   !RS: Replace: CPA (2/19/14)
    CALL CalcMeanProp(tAiCoil,tAoCoil,tAmod)    !RS Comment: Mean Air Coil Temperature
    CALL CalcMeanProp(hAiCoil,hAoCoil,hAmod)    !RS Comment: Mean Air Coil Humidity
    humrat=HUMTH(tAmod,hAmod)   !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
    CPair=CPAirFunction(tAmod,humrat) !AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAiCoil*CPAir  !RS Comment: Finding the capacity rate of air

    IF (ABS(QcoilSens) .GT. ABS(Qcoil)) THEN !Make sure sensible heat no larger than total heat, ISI - 12/07/07
        QcoilSens=Qcoil
    END IF

    tAoCoil=tAiCoil+QcoilSens/Cair  !RS Comment: Finding the Coil Outlet Air Temperature
    hAoCoil=hAiCoil+Qcoil/mAiCoil   !RS Comment: Finding the Coil Outlet Air Enthalpy

    !Fan air side inlet condition
    !CPair=CPA(REAL(tAoCoil))    !RS Comment: Finding the specific heat of air   !RS: Replace: CPA (2/19/14)
    humrat=HUMTH(tAoCoil,hAoCoil)  !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
    CPair=CPAirFunction(tAoCoil,humrat)  !RS: Replace: CPA (2/19/14)
    Cair=mAiCoil*CPAir  !RS Comment: Finding the capacity rate of air

    !IF (SystemType .NE. REHEAT) THEN !For reheat system, skip this !RS: Debugging: Already accounting for supply fan in E+ (9/3/14)
    !    IF (DrawBlow .EQ. DRAWTHROUGH) THEN !Draw through
    !        tAoCoil=tAoCoil+PwrFan/Cair
    !        hAoCoil=hAoCoil+PwrFan/mAiCoil
    !    END IF
    !END IF

    !Determining inlet and outlet air properties
    AirPropOpt=1
    AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APEnth=hAiCoil  !RS: Debugging: Formerly AirProp(4)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
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

    hRoCoil=hRiCoil-Qcoil/mRefTot   !RS Comment: Determining the Coil Outlet Refrigerant Enthalpy

    Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
    Enthalpy=hRoCoil*1000   !RS Comment: Unit Conversion
    tRoCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Coil Outlet Refrigerant Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2646'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF
    xRoCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !RS Comment: Coil Outlet Refrigerant Quality
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2652'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF

    Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
    Quality=1
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2672'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF

    IF (xRoCoil .GE. 1) THEN
        tSHoCoil=tRoCoil-tSat
    ELSE
        tSHoCoil=0
    END IF

    !****** Suction line calculation ******
    IF (LsucLn .GT. 0) THEN
        CALL SuctionLine
    ELSE
        pRiCmp=pRoCoil
        hRiCmp=hRoCoil
    END IF

    IF (SystemType .EQ. HEATPUMP) THEN !Heat pump
        !Calculate reversing valve heat transfer and pressure drop
        Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
        Enthalpy=hRiCmp*1000    !RS Comment: Unit Conversion
        tRiCmp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !RS Comment: Compressor Refrigerant Inlet Temperature
        IF (RefPropErr .GT. 0) THEN
            WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2646'
            ErrorFlag=REFPROPERROR
            EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
            CALL Evaporator_Helper_1
            RETURN
        END IF

        hRsuc=hRiCmp
        IF (IsCoolingMode .GT. 0) THEN
            !Correlation good for cooling only
            CALL SuctionHeatTransfer(mRefTot,tRdis,tRiCmp,hRsuc,hRiCmp)
            CALL SuctionPressureDrop(mRefTot,DPvalve)    	
            pRiCmp=pRiCmp-DPvalve   !Compressor Inlet Refrigerant Pressure
        END IF

    END IF

    Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
    Enthalpy=hRiCmp*1000    !RS Comment: Unit Conversion
    tRiCmp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !RS Comment: Compressor Refrigerant Inlet Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2705'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF
    xRiCmp=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !RS Comment: Compressor Refrigerant Inlet Quality
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2711'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF

    Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
    Quality=1
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !RS Comment: Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2730'
        ErrorFlag=REFPROPERROR
        EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)
        CALL Evaporator_Helper_1
        RETURN
    END IF

    IF (xRiCmp .GE. 1) THEN
        tSHiCmp=tRiCmp-tSat
    ELSE
        tSHiCmp=0
    END IF

    !Populating the OUT array
    EvapOUT%EOutpRoC=pRoCoil  !RS: Debugging: Formerly OUT(1)
    EvapOUT%EOuthRoC=hRoCoil  !RS: Debugging: Formerly OUT(2)
    EvapOUT%EOuttAoC=tAoCoil !RS: Debugging: Formerly OUT(3)
    EvapOUT%EOutrhAoC=rhAoCoil    !RS: Debugging: Only used as an output  !RS: Debugging: Formerly OUT(4)
    EvapOUT%EOutDPair=DPair   !RS: Debugging: Only used as an output  !RS: Debugging: Formerly OUT(5)
    EvapOUT%EOutpRiC=pRiCmp   !RS: Debugging: Formerly OUT(6)
    EvapOUT%EOuthRiC=hRiCmp   !RS: Debugging: Formerly OUT(7)
    EvapOUT%EOuttRiC=tRiCmp   !RS: Debugging: Formerly OUT(8)
    EvapOUT%EOutxRiC=xRiCmp   !RS: Debugging: Formerly OUT(9)
    EvapOUT%EOuttSHiC=tSHiCmp !RS: Debugging: Formerly OUT(10)
    EvapOUT%EOutQC=Qcoil   !RS: Debugging: Formerly OUT(11)
    EvapOUT%EOutQCSens=QcoilSens   !RS: Debugging: Only used as an output  !RS: Debugging: Formerly OUT(12)
    EvapOUT%EOutMSucLn=MassSucLn   !RS: Debugging: Formerly OUT(13)
    EvapOUT%EOutMC=0  !RS: Debugging: Never really used?   !RS: Debugging: Formerly OUT(14)
    EvapOUT%EOutWtAl=WeightAluminum  !RS: Debugging: Only used as an output  !RS: Debugging: Formerly OUT(15)
    EvapOUT%EOutWtCu=WeightCopper    !RS: Debugging: Only used as an output  !RS: Debugging: Formerly OUT(16)
    EvapOUT%EOutErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(17)

    CALL Evaporator_Helper_1

    RETURN

    END SUBROUTINE Evaporator

    SUBROUTINE Evaporator_Helper_1

    IF(CoolHeatModeFlag == 1) THEN
        CoilParams(1)%CoilInletRefTemp=tRiCoil
        CoilParams(1)%CoilOutletRefTemp=tRoCoil
        CoilParams(1)%CoilInletAirTemp=tAiCoil
        CoilParams(1)%CoilOutletAirTemp=tAoCoil
        CoilParams(1)%DPAir=DPair
        CoilParams(1)%DPRef=pRicoil-pRocoil
    ELSE IF(CoolHeatModeFlag == 0) THEN
        CoilParams(2)%CoilInletRefTemp=tRiCoil
        CoilParams(2)%CoilOutletRefTemp=tRoCoil
        CoilParams(2)%CoilInletAirTemp=tAiCoil
        CoilParams(2)%CoilOutletAirTemp=tAoCoil
        CoilParams(2)%DPAir=DPair
        CoilParams(2)%DPRef=pRicoil-pRocoil
    END IF

    END SUBROUTINE Evaporator_Helper_1

!************************************************************************

SUBROUTINE CalcEvaporatorInventory(MassCoil,MassLiqCoil,MassVapCoil,LiqTubeLength,VapTubeLength,TwoPhaseTubeLength,NumLiqTubes)

!------------------------------------------------------------------------
!Purpose:
!To get calculate refrigerant inventory in evaporator
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
USE CoilCalcMod

IMPLICIT NONE

REAL,INTENT(OUT) :: MassCoil           !Total refrigerant inventory in coil, kg
REAL,INTENT(OUT) ::  MassLiqCoil       !Total liquid refrigerant inventory in coil, kg
REAL,INTENT(OUT) ::  MassVapCoil       !Total vapor refrigerant inventory in coil, kg
REAL,INTENT(OUT) ::  NumLiqTubes       !Number of Liquid tubes
REAL,INTENT(OUT) :: LiqTubeLength      !Liquid tube length, m
REAL,INTENT(OUT) :: VapTubeLength      !Vapor tube length, m
REAL,INTENT(OUT) :: TwoPhaseTubeLength !Two-phase tube length, m

INTEGER :: CoilType !1=Condenser; 2=Evaporator; 
                    !3=High side interconnecting pipes; 
					!4=Low side interconnecting pipes
					!5=Microchannel condenser
					!6=Microchannel evaporator

INTEGER NumSection,I,J,K,II,III,IV !Loop Counter
REAL Lregion !Region length, m

  LiqTubeLength=0.0
  VapTubeLength=0.0
  TwoPhaseTubeLength=0.0

  MassCoil=0
  MassLiqCoil=0
  MassVapCoil=0

  IF (NumOfChannels .GT. 1) THEN
	  CoilType = MCEVAPORATOR !Microchannel coil
  ELSE
	  CoilType = EVAPORATORCOIL !Fin-tube coil
  END IF

  IF (CoilType .NE. MCEVAPORATOR) THEN !Microchannel coil

	  DO NumSection=1, NumOfSections
	  
	      DO I=1, CoilSection(NumSection)%NumOfCkts
		      CoilSection(NumSection)%Ckt(I)%Qckt=0.0
		      CoilSection(NumSection)%Ckt(I)%tSat=0.0
        
		      !Find first and last simulation tubes
		      IF (IsSimpleCoil .EQ. 1) THEN
			      FirstTube=1
			      LastTube=1
		      ELSE
			      FirstTube=1
			      LastTube=CoilSection(NumSection)%Ckt(I)%Ntube
			      IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
				      FirstTube=2 
			      END IF 

			      IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
				      LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 
			      END IF
			      IF (FirstTube .GT. LastTube) THEN
                      LastTube=FirstTube !Dummy tube
                  END IF
		      END IF

		      DO J=1,LastTube !Ckt(I)%Ntube !ISI - 10/30/06
			      DO K=1,NumOfMods
            
                      !Defining the module properties
				      pRiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%pRi  !Module Inlet Refrigerant Pressure
				      hRiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hRi  !Module Inlet Refrigerant Enthalpy
				      pRoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%pRo  !Module Outlet Refrigerant Pressure
				      hRoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hRo  !Module Outlet Refrigerant Enthalpy

				      LmodTube=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Len
				      AiMod=AiCoil*LmodTube/Lcoil

				      Pressure=pRiMod*1000  !RS Comment: Unit Conversion
				      Enthalpy=hRiMod*1000  !RS Comment: Unit Conversion
				      tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Inlet Refrigerant Temperature
				      xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Inlet Refrigerant Quality
				      vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
				      vRiMod=1/vRiMod   !Module Inlet Refrigerant Specific Volume
				      cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Inlet Refrigerant Specific Heat
				      cpRiMod=cpRiMod/1000  !RS Comment: Unit Conversion
				      muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Inlet Refrigerant Dynamic Viscosity
				      kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Inlet Refrigerant Thermal Conductivity
				      kRiMod=kRiMod/1000    !RS Comment: Unit Conversion

				      Quality=1
				      vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vgRiMod=1/vgRiMod !Module Inlet Refrigerant Vapor Specific Volume

				      Quality=0
				      vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vfRiMod=1/vfRiMod !Module Inlet Refrigerant Liquid Specific Volume

				      IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
					      cpRiMod=0
					      muRiMod=0
					      kRiMod=0
				      END IF

				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Enthalpy=hRoMod*1000  !RS Comment: Unit Conversion
				      tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Temperature
				      xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Quality
				      vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
				      vRoMod=1/vRoMod   !Module Outlet Refrigerant Specific Volume
				      cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Outlet Refrigerant Specific Heat
				      cpRoMod=cpRoMod/1000  !RS Comment: Unit Conversion
				      muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Outlet Refrigerant Dynamic Viscosity
				      kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Thermal Conductivity
				      kRoMod=kRoMod/1000    !RS Comment: Unit Conversion

				      Quality=1
				      vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vgRoMod=1/vgRoMod !Module Outlet Refrigerant Vapor Specific Volume

				      Quality=0
				      vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vfRoMod=1/vfRoMod !Module Outlet Refrigerant Liquid Specific Volume

				      Temperature=tRoMod
				      Quality=1
				      Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)    !Saturation Pressure
				      Psat=Psat/1000    !RS Comment: Unit Conversion
    		  
				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Quality=0
				      hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Enthalpy
				      hfRoMod=hfRoMod/1000  !RS Comment: Unit Conversion
				      CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Module Outlet Refrigerant Liquid Specific Heat
				      CpfRoMod=CpfRoMod/1000    !RS Comment: Unit Conversion
				      mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Module Outlet Refrigerant Liquid Dynamic Viscosity
				      kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Module Outlet Refrigerant Liquid Thermal Conductivity
				      kfRoMod=kfRoMod/1000  !RS Comment: Unit Conversion

				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Quality=1
				      hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Module Outlet Refrigerant Vapor Enthalpy
				      hgRoMod=hgRoMod/1000  !RS Comment: Unit Conversion
				      CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Module Outlet Refrigerant Vapor Specific Heat
				      CpgRoMod=CpgRoMod/1000    !RS Comment: Unit Conversion
				      mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Module Outlet Refrigerant Vapor Dynamic Viscosity
				      kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Module Outlet Refrigerant Vapor Thermal Conductivity
				      kgRoMod=kgRoMod/1000  !RS Comment: Unit Conversion

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

                      !Defining the other module values
				      mRefMod=CoilSection(NumSection)%Ckt(I)%mRef

				      tAiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tAi
				      rhAiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%rhAi

				      tAoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tAo
				      rhAoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%rhAo

				      Qmod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Qmod
				      hciMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hci
				      hcoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hco
				      ReVap=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%ReVap
				      ReLiq=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%ReLiq
				      Cair=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%cAir
				      Rair=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Rair 
				      Rtube=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Rtube

				      IF ((K .EQ. NumOfMods .AND. J .NE. LastTube) .OR. &
				          (J .EQ. LastTube .AND. (CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1 .OR. CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1))) THEN !ISI - 02/05/07
					      !Include return bend length
					      CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									     xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									     LmodTube+Lreturnbend,LmodTP,LmodTP,MassLiqMod,MassVapMod,MassMod)
                    !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
                    !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
                    !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

				      ELSE
					      CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									     xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									     LmodTube,LmodTP,LmodTP,MassLiqMod,MassVapMod,MassMod)
				      END IF
    			  
				      CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Mass=MassMod

				      IF (K .EQ. NumOfMods .OR. &
				          (J .EQ. LastTube .AND. &
				          (CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1 .OR. &
				           CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1))) THEN
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
			      EXIT
		      END IF
	      END DO !end circuit

	  END DO !end section
	  
	  NumLiqTubes=LiqTubeLength/Lcoil

  ELSE

	  DO I=1, Nl    !Number of Rows

		  DO II=1,Slab(I)%Npass

			  DO III=1,1 !NumOfTubes
				  
				  DO IV=1, NumOfMods
        
                      !Defining module properties
					  pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
					  hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi
					  pRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo
					  hRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo

					  LmodTube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len
					  AiMod=AiCoil*LmodTube/Lcoil

					  Pressure=pRiMod*1000  !RS Comment: Unit Conversion
					  Enthalpy=hRiMod*1000  !RS Comment: Unit Conversion
					  tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Temperature
					  xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Quality
					  vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
					  vRiMod=1/vRiMod   !Module Refrigerant Inlet Specific Volume
					  cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Specific Heat
					  cpRiMod=cpRiMod/1000  !RS Comment: Unit Conversion
					  muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Dynamic Viscosity
					  kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Thermal Conductivity
					  kRiMod=kRiMod/1000    !RS Comment: Unit Conversion

					  Quality=1
					  vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vgRiMod=1/vgRiMod !Module Refrigerant Inlet Vapor Specific Heat

					  Quality=0
					  vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vfRiMod=1/vfRiMod !Module Refrigerant Inlet Liquid Specific Heat

					  IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
						  cpRiMod=0
						  muRiMod=0
						  kRiMod=0
					  END IF

					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Enthalpy=hRoMod*1000  !RS Comment: Unit Conversion
					  tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Temperature
					  xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Quality
					  vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
					  vRoMod=1/vRoMod   !Module Outlet Refrigerant Specific Volume
					  cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Outlet Refrigerant Specific Heat
					  cpRoMod=cpRoMod/1000  !RS Comment: Unit Conversion
					  muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Outlet Refrigerant Dynamic Viscosity
					  kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Outlet Refrigerant Thermal Conductivity
					  kRoMod=kRoMod/1000    !RS Comment: Unit Conversion

					  Quality=1
					  vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vgRoMod=1/vgRoMod !Module Outlet Refrigerant Vapor Specific Heat

					  Quality=0
					  vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vfRoMod=1/vfRoMod !Module Outlet Refrigerant Liquid Specific Heat

					  Temperature=tRoMod
					  Quality=1
					  IF (tRoMod+273.15 .GT. Tcr) THEN
						  Psat=pRoMod
					  ELSE 
						  Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)    !Saturation Pressure
						  Psat=Psat/1000    !RS Comment: Unit Conversion
                      END IF
			  
                      !RS Comment: Module Outlet Refrigerant Liquid Properties
					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Quality=0
					  hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)
					  hfRoMod=hfRoMod/1000  !RS Comment: Unit Conversion
					  CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)
					  CpfRoMod=CpfRoMod/1000    !RS Comment: Unit Conversion
					  mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)
					  kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)
					  kfRoMod=kfRoMod/1000  !RS Comment: Unit Conversion

                      !RS Comment: Module Outlet Refrigerant Vapor Vapor Properties
					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Quality=1
					  hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)
					  hgRoMod=hgRoMod/1000  !RS Comment: Unit Conversion
					  CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)
					  CpgRoMod=CpgRoMod/1000    !RS Comment: Unit Conversion
					  mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)
					  kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)
					  kgRoMod=kgRoMod/1000  !RS Comment: Unit Conversion

					  IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
						  cpRoMod=0
						  muRoMod=0
						  kRoMod=0
                      END IF

                      !Average Values
					  mu=(muRiMod+muRoMod)/2
					  kRef=(kRiMod+kRoMod)/2
					  cpRef=(cpRiMod+cpRoMod)/2
					  rhoRef=(1/vRiMod+1/vRoMod)/2

                      !Defining the other module values
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

					  CALL Inventory(CoilType,Dchannel,mRefMod/NumOfChannels,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									 xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									 LmodTube,LmodTP,(LmodTube-LmodTP),MassLiqMod,MassVapMod,MassMod)
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

	  NumLiqTubes=LiqTubeLength/Lcoil

  END IF

  RETURN

END SUBROUTINE CalcEvaporatorInventory

!************************************************************************

SUBROUTINE PrintEvaporatorResult

!------------------------------------------------------------------------
!Purpose:
!To print simulation result to output file "evaporator.csv"
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

IMPLICIT NONE

INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

REAL :: MassCoil    !Total refrigerant inventory in coil, kg
REAL :: MassLiqCoil !Total liquid refrigerant inventory in coil, kg
REAL :: MassVapCoil !Total vapor refrigerant inventory in coil, kg

INTEGER NumSection,I,J,K,II,III,IV !Loop Counter
CHARACTER(LEN=13),PARAMETER :: FMT_100 = "(50(A12,','))"
CHARACTER(LEN=25),PARAMETER :: FMT_104 = "(3(I3,','),50(F10.3,','))"

  OPEN (17,FILE='Evaporator.csv')
  !OPEN (17,FILE='Evaporator_PlainFin.csv')    !RS: Test case evaporator output file
  QModLatTot=0
  QModSensTot=0
  QModTot=0

  MassCoil=0
  MassLiqCoil=0
  MassVapCoil=0

  IF (NumOfChannels .GT. 1) THEN
	  CoilType = MCEVAPORATOR
  ELSE
	  CoilType = EVAPORATORCOIL
  END IF

  IF (CoilType .NE. MCEVAPORATOR) THEN

	  !WRITE(17,FMT_100)'Nckt','Ntube','Nmod','tRi(C)','tRo(C)','pRi(kPa)','pRo(kPa)', &
			!	   'hRi(kJ/kg)','hRo(kJ/kg)','xRi','xRo','tAi(C)','tAo(C)', &
			!	   'rhAi','rhAo','hci(kW/m2K)','hco(kW/m2K)', &
			!	   'mu(uPa-s)','k(W/mK)','cp(kJ/kgK)','rho(kg/m3)','ReVap','ReLiq', &
			!	   'QmodTot(W)','QmodSens(W)','MassLiq(g)','MassVap(g)','Mass(g)','DryWet', &
			!	   'mdotR(kg/hr)','mdotA(kg/s)'
      WRITE(17,FMT_100)'Nckt','Ntube','Nmod','tRi(C)','tRo(C)','pRi(kPa)','pRo(kPa)', &
				   'hRi(kJ/kg)','hRo(kJ/kg)','xRi','xRo','tAi(C)','tAo(C)', &
				   'rhAi','rhAo','hci(kW/m2K)','hco(kW/m2K)', &
				   'mu(uPa-s)','k(W/mK)','cp(kJ/kgK)','rho(kg/m3)','ReVap','ReLiq', &
				   'QmodTot(W)','QmodSens(W)','QmodLat(W)','MassLiq(g)','MassVap(g)','Mass(g)','DryWet', &
				   'mdotR(kg/hr)','mdotA(kg/s)' !, 'cpAir', 'hAiMod kJ/kg', 'hAoMod' !RS: Adding in the latent heat !RS: Debugging: Adding in air specific heat and h's

	  DO NumSection=1, NumOfSections

	      DO I=1, CoilSection(NumSection)%NumOfCkts
		      CoilSection(NumSection)%Ckt(I)%Qckt=0.0
		      CoilSection(NumSection)%Ckt(I)%QcktSens=0.0
		      CoilSection(NumSection)%Ckt(I)%tSat=0.0
        
		      !Find first and last simulation tubes
		      FirstTube=1
		      LastTube=CoilSection(NumSection)%Ckt(I)%Ntube
		      IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
			      FirstTube=2 
		      END IF 
		      IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
			      LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 
		      END IF
		      IF (FirstTube .GT. LastTube) THEN
                  LastTube=FirstTube !Dummy tube
              END IF

		      DO J=1,CoilSection(NumSection)%Ckt(I)%Ntube
			      DO K=1,NumOfMods
            
                      !Defining the module values
				      pRiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%pRi
				      hRiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hRi
				      pRoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%pRo
				      hRoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hRo

				      LmodTube=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Len
				      AiMod=AiCoil*LmodTube/Lcoil

				      Pressure=pRiMod*1000  !RS Comment: Unit Conversion
				      Enthalpy=hRiMod*1000  !RS Comment: Unit Conversion
				      tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Temperature
				      xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Refrigerant Inlet Quality
				      vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
				      vRiMod=1/vRiMod   !Module Refrigerant Inlet Specific Volume
				      cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Refrigerant Inlet Specific Heat
				      cpRiMod=cpRiMod/1000  !RS Comment: Unit Conversion
				      muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Refrigerant Inlet Dynamic Viscosity
				      kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Inlet Thermal Conductivity
				      kRiMod=kRiMod/1000    !RS Comment: Unit Conversion

				      Quality=1
				      vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vgRiMod=1/vgRiMod !Module Refrigerant Inlet Vapor Specific Volume

				      Quality=0
				      vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vfRiMod=1/vfRiMod !Module Refrigerant Inlet Liquid Specific Volume

				      IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
					      cpRiMod=0
					      muRiMod=0
					      kRiMod=0
				      END IF

				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Enthalpy=hRoMod*1000  !RS Comment: Unit Conversion
				      tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Temperature
				      xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Refrigerant Outlet Quality
				      vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
				      vRoMod=1/vRoMod   !Module Refrigerant Outlet Specific Volume
				      cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Specific Heat
				      cpRoMod=cpRoMod/1000  !RS Comment: Unit Conversion
				      muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Dynamic Viscosity
				      kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Thermal Conductivity
				      kRoMod=kRoMod/1000    !RS Comment: Unit Conversion

				      Quality=1
				      vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vgRoMod=1/vgRoMod !Module Refrigerant Outlet Vapor Specific Volume

				      Quality=0
				      vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
				      vfRoMod=1/vfRoMod !Module Refrigerant Outlet Liquid Specific Volume

				      Temperature=tRoMod
				      Quality=1
				      Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)    !Saturation Pressure
				      Psat=Psat/1000    !RS Comment: Unit Conversion

				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Quality=0
				      hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Liquid Enthalpy
				      hfRoMod=hfRoMod/1000  !RS Comment: Unit Conversion
				      CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Liquid Specific Heat
				      CpfRoMod=CpfRoMod/1000    !RS Comment: Unit Conversion
				      mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Liquid Specific Heat
				      kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Thermal Conductivity
				      kfRoMod=kfRoMod/1000  !RS Comment: Unit Conversion

				      Pressure=pRoMod*1000  !RS Comment: Unit Conversion
				      Quality=1
				      hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Vapor Enthalpy
				      hgRoMod=hgRoMod/1000  !RS Comment: Unit Conversion
				      CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Module Refrigerant Outlet Vapor Specific Heat
				      CpgRoMod=CpgRoMod/1000    !RS Comment: Unit Conversion
				      mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Module Refrigerant Outlet Vapor Dynamic Viscosity
				      kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Module Refrigerant Outlet Vapor Thermal Conductivity
				      kgRoMod=kgRoMod/1000  !RS Comment: Unit Conversion

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

                      !Defining the other module values
				      mRefMod=CoilSection(NumSection)%Ckt(I)%mRef

				      Qmod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Qmod
				      hciMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hci
				      hcoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%hco
				      ReVap=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%ReVap
				      ReLiq=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%ReLiq
				      DryWet=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%DryWet
				      Cair=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%cAir
				      Rair=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Rair
				      Rtube=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Rtube

				      IF (K .EQ. NumOfMods .OR. &
				          (J .EQ. LastTube .AND. &
				          (CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1 .OR. &
				           CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1))) THEN
				           
						    !Include return bend length
					      CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									     xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									     LmodTube+Lreturnbend,LmodTP,LmodTP,MassLiqMod,MassVapMod,MassMod)
                    !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
                    !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
                    !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

				      ELSE
					      CALL Inventory(CoilType,IDtube,mRefMod,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									     xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									     LmodTube,LmodTP,LmodTP,MassLiqMod,MassVapMod,MassMod)

				      END IF
    			  
				      CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Mass=MassMod

				      !Total mass inventory, avoid duplication
				      IF (J .GE. FirstTube .AND. J .LE. LastTube) THEN
					      MassCoil=MassCoil+MassMod 
					      MassLiqCoil=MassLiqCoil+MassLiqMod
					      MassVapCoil=MassVapCoil+MassVapMod
                      END IF

                      !Defining more module values
				      mAiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%mAi
				      tAiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tAi
				      rhAiMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%rhAi

				      tAoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%tAo
				      rhAoMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%rhAo

				      MassMod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Mass

				      Qmod=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%Qmod  !RS: Commenting out because it should already have been defined
                      
                      QmodTot=Qmod+QmodTot  !RS: Debugging: Counting up the total QMod

				      QmodSens=CoilSection(NumSection)%Ckt(I)%Tube(J)%Seg(K)%QmodSens
                      QmodSensTot=QmodSens+QmodSensTot !RS: Debugging: Not sure how to actually count this up!

				      mRefCkt=CoilSection(NumSection)%Ckt(I)%mRef

                      !Keeping the module quality in the form of a decimal
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

				      !WRITE(17,FMT_104)I,J,K,tRiMod,tRoMod,pRiMod,pRoMod,hRiMod,hRoMod, &
							   !    xRiMod,xRoMod,tAiMod,tAoMod,rhAiMod,rhAoMod, &
							   !    hciMod,hcoMod,mu*1e6,kRef*1e3,cpRef,rhoRef,ReVap,ReLiq, &
							   !    Qmod*1000,QmodSens*1000,MassLiqMod*1000,MassVapMod*1000,MassMod*1000, &
							   !    FLOAT(DryWet),mRefCkt*3600,mAiMod
                     WRITE(17,FMT_104)I,J,K,tRiMod,tRoMod,pRiMod,pRoMod,hRiMod,hRoMod, &
							       xRiMod,xRoMod,tAiMod,tAoMod,rhAiMod,rhAoMod, &
							       hciMod,hcoMod,mu*1e6,kRef*1e3,cpRef,rhoRef,ReVap,ReLiq, &
							       Qmod*1000,QmodSens*1000,QmodLat*1000,MassLiqMod*1000,MassVapMod*1000,MassMod*1000, &
							       FLOAT(DryWet),mRefCkt*3600, mAiMod !, cpAIR, hAiMod, hAoMod, AirProp(4), TestH   !RS: Trying to find the latent heat !RS: Debugging: Adding in the air cp and h's

			      END DO !end Nmod

		      END DO !end Ntube
        
		      IF (EqCircuits .EQ. 1 .AND. IsUniformVelProfile) THEN
                  EXIT
              END IF

	      END DO !end circuit

      END DO !End section
  
  ELSE

	  WRITE(16,FMT_100)'Nslab','Npass','Nmod','tRi(C)','tRo(C)','pRi(kPa)','pRo(kPa)', &
				   'hRi(kJ/kg)','hRo(kJ/kg)','xRi','xRo','tAi(C)','tAo(C)', &
  				   'rhAi','rhAo','hci(W/m2K)','hco(W/m2K)', &
				   'mu(uPa-s)','k(W/mK)','cp(kJ/kgK)','rho(kg/m3)','ReVap','ReLiq', &
				   'Qmod(W)','MassLiq(g)','MassVap(g)','MassTot(g)','mdotR(kg/h)','mdotA(kg/s)' 

	  DO I=1, Nl    !Number of Rows

		  DO II=1,Slab(I)%Npass

			  DO III=1,1 !Num Of tubes
				
				  DO IV=1, NumOfMods
        
                      !Defining module properties
					  pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
					  hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi
					  pRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo
					  hRoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo

					  LmodTube=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len
					  AiMod=AiCoil*LmodTube/Lcoil

					  Pressure=pRiMod*1000  !RS Comment: Unit Conversion
					  Enthalpy=hRiMod*1000  !RS Comment: Unit Conversion
					  tRiMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Inlet Refrigerant Temperature
					  xRiMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Inlet Refrigerant Quality
					  vRiMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
					  vRiMod=1/vRiMod   !Module Inlet Refrigerant Specific Volume
					  cpRiMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Inlet Refrigerant Specific Heat
					  cpRiMod=cpRiMod/1000  !RS Comment: Unit Conversion
					  muRiMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Inlet Refrigerant Dynamic Viscosity
					  kRiMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Inlet Refrigerant Thermal Conductivity
					  kRiMod=kRiMod/1000    !RS Comment: Unit Conversion

					  Quality=1
					  vgRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vgRiMod=1/vgRiMod !Module Inlet Refrigerant Vapor Specific Volume

					  Quality=0
					  vfRiMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vfRiMod=1/vfRiMod !Module Inlet Refrigerant Liquid Specific Volume

					  IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
						  cpRiMod=0
						  muRiMod=0
						  kRiMod=0
					  END IF

					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Enthalpy=hRoMod*1000  !RS Comment: Unit Conversion
					  tRoMod=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Temperature
					  xRoMod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Module Outlet Refrigerant Quality
					  vRoMod=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
					  vRoMod=1/vRoMod   !Module Outlet Refrigerant Specific Volume
					  cpRoMod=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Module Outlet Refrigerant Specific Heat
					  cpRoMod=cpRoMod/1000  !RS Comment: Unit Conversion
					  muRoMod=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Module Outlet Refrigerant Dynamic Viscosity
					  kRoMod=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Module Outlet Refrigerant Thermal Conductivity
					  kRoMod=kRoMod/1000    !RS Comment: Unit Conversion

					  Quality=1
					  vgRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vgRoMod=1/vgRoMod !Module Outlet Refrigerant Vapor Specific Volume

					  Quality=0
					  vfRoMod=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
					  vfRoMod=1/vfRoMod !Module Outlet Refrigerant Liquid Specific Volume

					  Temperature=tRoMod
					  Quality=1
					  IF (tRoMod+273.15 .GT. Tcr) THEN  !Finding the saturation pressure
						  Psat=pRoMod
					  ELSE 
						  Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)
						  Psat=Psat/1000    !RS Comment: Unit Conversion
                      END IF
			  
                      !RS Comment: Module Outlet Refrigerant Liquid Properties
					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Quality=0
					  hfRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)
					  hfRoMod=hfRoMod/1000  !RS Comment: Unit Conversion
					  CpfRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)
					  CpfRoMod=CpfRoMod/1000    !RS Comment: Unit Conversion
					  mufRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)
					  kfRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)
					  kfRoMod=kfRoMod/1000  !RS Comment: Unit Conversion

                      !RS Comment: Module Outlet Refrigerant Vapor Properties
					  Pressure=pRoMod*1000  !RS Comment: Unit Conversion
					  Quality=1
					  hgRoMod=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)
					  hgRoMod=hgRoMod/1000  !RS Comment: Unit Conversion
					  CpgRoMod=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)
					  CpgRoMod=CpgRoMod/1000    !RS Comment: Unit Conversion
					  mugRoMod=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)
					  kgRoMod=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)
					  kgRoMod=kgRoMod/1000  !RS Comment: Unit Conversion

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

                      !Defining the other module values
					  mRefMod=mRefTot/Slab(I)%Pass(II)%Ntube !Ckt(I)%mRef

					  mAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%mAi
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

					  CALL Inventory(CoilType,Dchannel,mRefMod/NumOfChannels,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &   !RS: Debugging: Extraneous
									 xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
									 LmodTube,LmodTP,(LmodTube-LmodTP),MassLiqMod,MassVapMod,MassMod)
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
				
                      !Keeping the module quality in the decimal form
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
								   mRefMod*3600,mAiMod

				  END DO !end Nmod

			  END DO !end Ntube
  
		  END DO !end circuit

	  END DO !end Slab

  END IF

  CLOSE(17)

  RETURN

END SUBROUTINE PrintEvaporatorResult

!************************************************************************

    SUBROUTINE InitEvaporatorCoil(CoilType)

    !------------------------------------------------------------------------
    !Purpose:
    !To initialize evaporator geometry and circuiting
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
    INTEGER NumSection !Loop counter, ISI - 09/10/07
    
  INTEGER, PARAMETER :: MaxNameLength = 200

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
    
    !INTEGER,PARAMETER :: SI=1
    INTEGER,PARAMETER :: IP=2
    
    IF (ErrorFlag .EQ. 1) THEN  !RS: Carrying over from the condenser initialization
        ErrorFlag = 0 !JG set a zero error
    END IF
    
    !***** Get circuit info *****
    IF (ErrorFlag .GT. CONVERGEERROR) THEN !NE. NOERROR) THEN   !RS: Debugging: Pushing through convergence errors
        ErrorFlag=CKTFILEERROR
        CALL InitEvaporatorCoil_Helper_1
        RETURN
    END IF


    !**************************** Model *************************************

          CALL GetObjectItem('IDCcktModel',1,Alphas,NumAlphas, &
                        TmpNumbers,NumNumbers,Status)
        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
        
        ModelName = Alphas(1)
            
        IF (ModelName .EQ. 'MICROCHANNEL COIL') THEN    !Checking whether MicroChannel or Evaporator Coil
            CoilType = MCEVAPORATOR
        ELSE
            CoilType = EVAPORATORCOIL
        END IF

IF (CoilType .EQ. EVAPORATORCOIL) THEN !Fin-tube coil or MicroChannel?
        
    IF (IsCoolingMode .GT. 0) THEN  !IDC or ODC ckt?
            !IDC ckt

        !**************************** Geometry *************************************

        !Reading in the variable values
            CALL GetObjectItem('IDCcktGeometry',1,Alphas,NumAlphas, &
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
                CALL InitEvaporatorCoil_Helper_1
                RETURN

            END IF

            NumOfMods = Numbers(14)
            TempNumofMods=NumOfMods
            NumOfCkts = Numbers(15)
            NumofSections = 1   !RS Comment: Not in the input file, but needed for the code to run properly. Set to 1 as there is only one coil section here.
            !RS: Debugging: The above NumofSections should probably be made as an input in the IDF and IDD; there may be other places in
            !RS: Debugging: (con.) the code where it's hardcoded to handle only one section. (12/31/13)

            SELECT CASE (Alphas(5)(1:1))    !Tube Shift Flag
            CASE ('F','f')
                IsShift=.FALSE.
            CASE DEFAULT
                IsShift=.TRUE.
            END SELECT

            
    !*************************** Inputs ***********************
        
        !***************** Indoor coil data ***************** !RS: Debugging: Evaporator & Condenser

  CALL GetObjectItem('IndoorCoilData',1,Alphas,NumAlphas, &
                      TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
  EvapPAR%EvapFinType = Numbers(1)  !IDC_FinType    !RS: Debugging: Formerly EvapPAR(22)
  
  IDC_FinName = Alphas(1)
  
  EvapPAR%EvapFinPitch = Numbers(2)  !RS: Debugging: Formerly EvapPAR(15)
  EvapPAR%EvapFinThermCon = Numbers(3) !Fin Conductivity    !RS: Debugging: Formerly EvapPAR(16)
  EvapPAR%EvapFinThick = Numbers(4)   !Fin Thickness !RS: Debugging: Formerly EvapPAR(14)
  EvapPAR%EvapTube = Numbers(5) !Numerical Denotion of the tube type    !RS: Debugging: Formerly EvapPAR(30)
  IDC_TubeID = Numbers(6)   !Tube Inner Diameter
  EvapPAR%EvapCoilTOD = Numbers(7)   !Tube Outer Diameter    !RS: Debugging: Formerly EvapPAR(8)
  EvapPAR%EvapCoilTThermCon = Numbers(8)    !Tube Conductivity    !RS: Debugging: Formerly EvapPAR(11)
  EvapPAR%EvapRspc = Numbers(9)   !Tube Lateral Spacing  !RS: Debugging: Formerly EvapPAR(13)
  EvapPAR%EvapTspc = Numbers(10)   !Tube Vertical Spacing    !RS: Debugging: Formerly EvapPAR(12)
  EvapPAR%EvapNl = Numbers(11)  !Number of Rows    !RS: Debugging: Formerly EvapPAR(18)
  EvapPAR%EvapNt = Numbers(12)  !Number of Tubes Per Row   !RS: Debugging: Formerly EvapPAR(17)
  EvapPAR%EvapNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly EvapPAR(19)
  EvapPAR%EvapNumMod = Numbers(14)    !Number of Segments  !RS: Debugging: Formerly EvapPAR(21)
  EvapPAR%EvapCoilSTLen = Numbers(15)   !Length of Tube   !RS: Debugging: Formerly EvapPAR(10)
  EvapPAR%EvapMultRefQT = Numbers(16)   !Ref Side Heat Transfer Multiplier    !RS: Debugging: Formerly EvapPAR(23)
  EvapPAR%EvapMultRefPD = Numbers(17) !Ref Side Pressure Drop Multiplier  !RS: Debugging: Formerly EvapPAR(24)
  EvapPAR%EvapMultAirQT = Numbers(18)   !Air Side Heat Transfer Multiplier    !RS: Debugging: Formerly EvapPAR(25)
  EvapPAR%EvapMultAirPD = Numbers(19) !Air Side Pressure Drop Multiplier  !RS: Debugging: Formerly EvapPAR(26)

  !Tube wall thickness, mm or mil
  EvapPAR%EvapCoilTWThick=(EvapPAR%EvapCoilTOD-IDC_TubeID)/2  !RS: Debugging: Formerly EvapPAR(9)
  IF (Unit .EQ. IP) THEN
      EvapPAR%EvapCoilTWThick=EvapPAR%EvapCoilTWThick*1000
  END IF

	EvapPAR%EvapCoolMode=IsCoolingMode   !RS: Debugging: Formerly EvapPAR(20)

	!EvapPAR(29)=IDC_SurfAbsorptivity   !RS: Debugging: Extraneous
    
    
  !***************** Indoor fan data *****************  !RS: Debugging: Moving: Evaporator & Condenser
  
  CALL GetObjectItem('IndoorFanData',1,Alphas,NumAlphas, &
                      TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
  EvapPAR%EvapFanPwr = Numbers(1) !Fan Power   !RS: Debugging: Formerly EvapPAR(27)
  !VdotIDfan = Numbers(2)    !Fan Air Flow Rate
  EvapPAR%EvapFanLoc = Numbers(3)   !Draw Through or Blow Through  !RS: Debugging: Formerly EvapPAR(28)
  
        
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

            IF (ErrorFlag .GT. CONVERGEERROR) THEN !NE. NOERROR) THEN    !RS: Debugging: Pushing through convergence errors
                ErrorFlag=CKTFILEERROR
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF

            NumOfTubes=Nl*Nt  

            !Fin spacing
            FinSpg = 1/FinPitch-FinThk

            !For plate finned tube
            FinHeight=0 !No use for plate finned tube
            TubeDepth=ODtube
            TubeHeight=ODtube
            Dchannel=IDtube
            NumOfChannels=1

            IF (FinSpg .LT. FinThk) THEN
                ErrorFlag=COILFINERROR
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF

            IF (Pt .LT. ODtube+2*FinThk) THEN
                ErrorFlag=COILTUBEERROR
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF
            
            IF (IsSimpleCoil .EQ. 1) THEN
                IF (.NOT. ALLOCATED(Ckt)) THEN
                !NumOfMods=2
                ALLOCATE(CoilSection(NumOfSections)) 
                ALLOCATE(Ckt(NumOfCkts))
                ALLOCATE(CoilSection(1)%Ckt(NumOfCkts)) 
                ALLOCATE(Tube(NumofTubes)) 
                ALLOCATE(SucLnSeg(NumOfMods))		  
                CoilSection(1)%NumOfCkts=NumOfCkts
                DO I=1, NumOfCkts
                    Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
                    ALLOCATE(Ckt(I)%Tube(1))
                    ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
                    ALLOCATE(Ckt(I)%TubeSequence((NumofTubes/NumOfCkts)))
                !IF (I .EQ. 1) THEN
                !    CALL GetObjectItem('IDCcktCircuit1_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)     
                !ELSEIF (I .EQ. 2) THEN
                !    CALL GetObjectItem('IDCcktCircuit2_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                !ELSEIF (I .EQ. 3) THEN
                !    CALL GetObjectItem('IDCcktCircuit3_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                !ELSEIF (I .EQ. 4) THEN
                !    CALL GetObjectItem('IDCcktCircuit4_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                !ELSEIF (I .EQ. 5) THEN
                !    CALL GetObjectItem('IDCcktCircuit5_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                !ELSE
                !    CALL GetObjectItem('IDCcktCircuit6_TubeSequence',1,Alphas,NumAlphas, &
                !                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                !END IF
                !        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
                        
                    DO J=1,Ckt(I)%Ntube
                        Ckt(I)%TubeSequence(J)=I !Numbers(J)   !RS Comment: Populating the tube sequence arrays
                    END DO
                    CoilSection(1)%Ckt(I)=Ckt(I)   
                END DO
                END IF
            
            END IF
            
            IF (.NOT. ALLOCATED(Ckt)) THEN
                CALL InitEvaporatorStructures()
                ! IF (TempNumofMods .NE. 2) THEN
                !    NumOfMods=TempNumofMods
                !    DEALLOCATE(SucLnSeg)
                !    ALLOCATE(SucLnSeg(NumOfMods))
                !    DO I=1, NumOfCkts
                !        DO J=1,Ckt(I)%Ntube
                !            ALLOCATE(Ckt(I)%Tube(J)%Seg(NumOfMods))
                !        END DO
                !    END DO
                !END IF
            END IF

            IF (IsSimpleCoil .NE. 1) THEN
            
            !IF (IsSimpleCoil .EQ. 1) THEN
            !    NumOfMods=2
            !    ALLOCATE(CoilSection(NumOfSections)) 
            !    ALLOCATE(Ckt(NumOfCkts))
            !    ALLOCATE(CoilSection(1)%Ckt(NumOfCkts)) 	  
            !    CoilSection(1)%NumOfCkts=NumOfCkts
            !    DO I=1, NumOfCkts
            !        Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
            !        ALLOCATE(Ckt(I)%Tube(1))
            !        ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
            !        CoilSection(1)%Ckt(I)=Ckt(I)
            !    END DO
            !ELSE
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
            !
            !    CALL GetObjectItem('IDCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
            !                        TmpNumbers,NumNumbers,Status) 
            !    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            !
            !    DO I=1, NumOfCkts
            !        Ckt(I)%Ntube = Numbers(I)
            !        IF (ErrorFlag .NE. NOERROR) THEN 
            !            ErrorFlag=CKTFILEERROR
            !            CALL InitEvaporatorCoil_Helper_1
            !            RETURN
            !        END IF
            !        ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))          
            !        ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))  
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
                CALL InitEvaporatorCoil_Helper_1
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
                ELSEIF (I .EQ. 9) THEN
                    CALL GetObjectItem('IDCcktCircuit9_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSEIF (I .EQ. 10) THEN
                    CALL GetObjectItem('IDCcktCircuit10_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSEIF (I .EQ. 11) THEN
                    CALL GetObjectItem('IDCcktCircuit11_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSE
                    CALL GetObjectItem('IDCcktCircuit12_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                END IF
                        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
                        
                    DO J=1, Ckt(I)%Ntube
                        Ckt(I)%TubeSequence(J)=Numbers(J)   !RS Comment: Populating the tube sequence arrays
                    END DO 
                IF (ErrorFlag .NE. NOERROR) THEN 
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
                    RETURN
                END IF
            END DO
        
        !************************* Velocity Profile ********************************

            CoilSection(NumOfSections)%NumOfCkts=NumOfCkts
            IF (.NOT. ALLOCATED(CoilSection(NumOfSections)%CktNum)) THEN
		      ALLOCATE(CoilSection(NumOfSections)%CktNum(CoilSection(NumOfSections)%NumOfCkts))
		      ALLOCATE(CoilSection(NumOfSections)%mRefIter(CoilSection(NumOfSections)%NumOfCkts))
            END IF
		      DO J=1, NumOfCkts
		          CoilSection(NumOfSections)%CktNum(J)=J    !Numbering the Circuits
		      END DO

            DO I=1,2
                IF (ErrorFlag .NE. NOERROR) THEN !Tube# ,velocity Deviation from mean value
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
                    RETURN
                END IF
            END DO

            !Section data, ISI - 09/10/07
            !Initialize
            CoilSection%IsInlet = .TRUE.
            IsUniformVelProfile=.TRUE.

            CALL GetObjectItem('IDCcktVelocityProfile',1,Alphas,NumAlphas, &
                                TmpNumbers,NumNumbers,Status) 
                        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
                        
            DO I=Nt*(Nl-1)+1,Nt*Nl !last row faces air inlet (Cross flow HX)
                DO J=1, NumOfMods
                    Tube(I)%Seg(J)%VelDev = Numbers(J)  !Bringing in velocity deviation data
                END DO
                IF (ErrorFlag .NE. NOERROR) THEN 
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
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
                CALL InitEvaporatorCoil_Helper_1
                RETURN

            END IF

            !Propagate circuit info to coil section, ISI - 09/10/07
            NumInletSections=0
            DO I=1, NumOfSections
                IF (CoilSection(I)%IsInlet) THEN
                    NumInletSections=NumInletSections+1
                END IF
                !ALLOCATE(CoilSection(I)%Ckt(CoilSection(I)%NumOfCkts)) !RS: Debugging:

                DO J=1, CoilSection(I)%NumOfCkts
                    CoilSection(I)%Ckt(J)=Ckt(CoilSection(I)%CktNum(J))
                END DO
            END DO

            !Determine inlet and outlet flags, split, joint or nothing, ISI - 09/10/07
            DO NumSection=1, NumOfSections
                CoilSection(NumSection)%Nnode=1

                Nnode=1
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    !Initialize
                    CoilSection(NumSection)%Ckt(I)%InJoin=0
                    CoilSection(NumSection)%Ckt(I)%InSplit=0
                    CoilSection(NumSection)%Ckt(I)%OutJoin=0
                    CoilSection(NumSection)%Ckt(I)%OutSplit=0
                    DO J=1, CoilSection(NumSection)%NumOfCkts
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(1) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube)) THEN
                            CoilSection(NumSection)%Ckt(I)%InJoin=CoilSection(NumSection)%Ckt(I)%InJoin+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(1) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(1)) THEN
                            CoilSection(NumSection)%Ckt(I)%InSplit=CoilSection(NumSection)%Ckt(I)%InSplit+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(CoilSection(NumSection)%Ckt(I)%Ntube) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube)) THEN
                            CoilSection(NumSection)%Ckt(I)%OutJoin=CoilSection(NumSection)%Ckt(I)%OutJoin+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(CoilSection(NumSection)%Ckt(I)%Ntube) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(1)) THEN
                            CoilSection(NumSection)%Ckt(I)%OutSplit=CoilSection(NumSection)%Ckt(I)%OutSplit+1
                        END IF
                    END DO
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1 .OR. CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1) THEN
                        Nnode=Nnode+1
                    END IF
                END DO !End NumOfCkts

                CoilSection(NumSection)%Nnode=Nnode
                IF (.NOT. ALLOCATED(CoilSection(NumSection)%Node)) THEN !RS: Debugging:
                    ALLOCATE(CoilSection(NumSection)%Node(Nnode))
                END IF

                !Find split and joint tube numbers
                J=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1 ) THEN
                        J=J+1
                        CoilSection(NumSection)%Node(J)%Num=CoilSection(NumSection)%Ckt(I)%TubeSequence(1)
                    ELSEIF (CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1) THEN
                        J=J+1
                        CoilSection(NumSection)%Node(J)%Num=CoilSection(NumSection)%Ckt(I)%TubeSequence(Ckt(I)%Ntube)
                    END IF
                    IF (J .GT. Nnode) THEN
                        EXIT
                    END IF
                END DO
                CoilSection(NumSection)%Node(Nnode)%Num=0 !section outlet 

            END DO !End NumOfSections

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

            !ISI - 09/10/07
            Tube%Empty = 1 !Initialize 
            DO NumSection=1, NumOfSections

                !Find even tubes
                DO I=1, CoilSection(NumSection)%NumOfCkts

                    !Find first and last simulation tubes
                    FirstTube=1
                    LastTube=CoilSection(NumSection)%Ckt(I)%Ntube
                    IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
                        FirstTube=2 !Skip first tube
                    END IF 
                    IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
                        LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 !Ignore last tube
                    END IF
                    IF (FirstTube .GT. LastTube) THEN
                        LastTube=FirstTube !Dummy tube
                    END IF

                    DO J=FirstTube, LastTube
                        TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)  !Determining the tube number
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
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    DO J=1, CoilSection(NumSection)%Ckt(I)%Ntube
                        TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
                        Tube(TubeNum)%Empty=0
                    END DO
                END DO

                !Number of inlet circuits
                CoilSection(NumSection)%NcktFirst=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .LT. 1) THEN
                        CoilSection(NumSection)%NcktFirst=CoilSection(NumSection)%NcktFirst+1
                    END IF
                END DO
                IF (CoilSection(NumSection)%NcktFirst .EQ. 0) THEN
                    CoilSection(NumSection)%NcktFirst = 1 !At least one circuit, ISI - 07/28/06
                END IF

                !Number of outlet circuits
                CoilSection(NumSection)%NcktLast=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%OutSplit .EQ. 0) THEN
                        CoilSection(NumSection)%NcktLast=CoilSection(NumSection)%NcktLast+1
                    END IF
                END DO
                IF (CoilSection(NumSection)%NcktLast .EQ. 0) THEN
                    CoilSection(NumSection)%NcktLast = 1 !At least one circuit, ISI - 07/28/06
                END IF
            END DO !End NumOfSections
            
            END IF !RS Comment: Adding in an END IF to close the above open block (Simple Condenser or not)
            
    ELSE  !ODC or IDC circuits?
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
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF
            
            NumOfMods = Numbers(14)
            NumOfCkts = Numbers(15)
            NumofSections = 1   !RS Comment: Not in the input file, but needed for the code to run properly. Set to 1 as there is only one coil section here.

            SELECT CASE (Alphas(5)(1:1))    !Tube Shift Flag
            CASE ('F','f')
                IsShift=.FALSE.
            CASE DEFAULT
                IsShift=.TRUE.
            END SELECT
            
        !*************************** Inputs ***************************************
        
          !***************** Outdoor coil data *****************    !RS: Debugging: Moving: Evaporator & Condenser

  CALL GetObjectItem('OutdoorCoilData',1,Alphas,NumAlphas, &
                       TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
  !Fin type (1-smooth; 2-Wavy; 3-louvered)

  EvapPAR%EvapFinType = Numbers(1)  !ODC_FinType    !RS: Debugging: Formerly EvapPAR(22)
  
  ODC_FinName = Alphas(1)
  
  EvapPAR%EvapFinPitch = Numbers(2)  !ODC_FinPitch   !RS: Debugging: Formerly EvapPAR(15)
  EvapPAR%EvapFinThermCon = Numbers(3) !Conductivity of Fin !RS: Debugging: Formerly EvapPAR(16)
  EvapPAR%EvapFinThick = Numbers(4)   !Fin Thickness !RS: Debugging: Formerly EvapPAR(14)
  EvapPAR%EvapTube = Numbers(5) !Numerical Denotion of Tube Type    !RS: Debugging: Formerly EvapPAR(30)
  ODC_TubeID = Numbers(6)   !Tube Inner Diameter
  EvapPAR%EvapCoilTOD = Numbers(7)   !Tube Outer Diameter    !RS: Debugging: Formerly EvapPAR(8)
  EvapPAR%EvapCoilTThermCon = Numbers(8)    !Tube Conductivity    !RS: Debugging: Formerly EvapPAR(11)
  EvapPAR%EvapRspc = Numbers(9)   !Tube Lateral Spacing  !RS: Debugging: Formerly EvapPAR(13)
  EvapPAR%EvapTspc = Numbers(10)   !Tube Vertical Spacing    !RS: Debugging: Formerly EvapPAR(12)
  EvapPAR%EvapNl = Numbers(11)  !Number of Rows    !RS: Debugging: Formerly EvapPAR(18)
  EvapPAR%EvapNt = Numbers(12)  !Number of Tubes per Row   !RS: Debugging: Formerly EvapPAR(17)
  EvapPAR%EvapNumCkt = Numbers(13)    !Number of Circuits  !RS: Debugging: Formerly EvapPAR(19)
  EvapPAR%EvapNumMod = Numbers(14)    !Number of Segments  !RS: Debugging: Formerly EvapPAR(21)
  EvapPAR%EvapCoilSTLen = Numbers(15)   !Single Tube Length   !RS: Debugging: Formerly EvapPAR(10)
  EvapPAR%EvapMultRefQT = Numbers(16)   !Ref Side Heat Transfer Multiplier    !RS: Debugging: Formerly EvapPAR(23)
  EvapPAR%EvapMultRefPD = Numbers(17) !Ref Side Pressure Drop Multiplier  !RS: Debugging: Formerly EvapPAR(24)
  EvapPAR%EvapMultAirQT = Numbers(18)   !Air Side Heat Transfer Multiplier    !RS: Debugging: Formerly EvapPAR(25)
  EvapPAR%EvapMultAirPD = Numbers(19) !Air Side Pressure Drop Multiplier  !RS: Debugging: Formerly EvapPAR(26)

    !Tube wall thickness, mm or mil
  EvapPAR%EvapCoilTWThick=(EvapPAR%EvapCoilTOD-ODC_TubeID)/2  !RS: Debugging: Formerly EvapPAR(29)
  IF (Unit .EQ. IP) THEN
      EvapPAR%EvapCoilTWThick=EvapPAR%EvapCoilTWThick*1000
  END IF

	EvapPAR%EvapCoolMode=IsCoolingMode   !RS: Debugging: Formerly EvapPAR(20)

	!EvapPAR(29)=ODC_SurfAbsorptivity   !RS: Debugging: Extraneous
    
  !***************** Outdoor fan data ***************** !RS: Debugging: Moving: Evaporator & Condenser

  CALL GetObjectItem('OutdoorFanData',1,Alphas,NumAlphas, &
                      TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
  EvapPAR%EvapFanPwr = Numbers(1) !Fan Power   !RS: Debugging: Formerly EvapPAR(27)
  !VdotODfan = Numbers(2)    !Fan Air Flow Rate
  EvapPAR%EvapFanLoc = Numbers(3)   !Draw Through (1) or Blow Through (2)  !RS: Debugging: Formerly EvapPAR(28)

        
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

            IF (ErrorFlag .NE. NOERROR) THEN 
                ErrorFlag=CKTFILEERROR
                CALL InitEvaporatorCoil_Helper_1
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
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF

            IF (Pt .LT. ODtube+2*FinThk) THEN
                ErrorFlag=COILTUBEERROR
                CALL InitEvaporatorCoil_Helper_1
                RETURN
            END IF
            
            IF (.NOT. ALLOCATED(Ckt)) THEN
                CALL InitEvaporatorStructures
            END IF
            
            IF (IsSimpleCoil .NE. 1) THEN

            !IF (IsSimpleCoil .EQ. 1) THEN   !This is an open block currently; it will need fixing. !RS Comment: It's closed now.
            !    NumOfMods=3
            !    ALLOCATE(Ckt(NumOfCkts))	  
            !    DO I=1, NumOfCkts
            !        Ckt(I)%Ntube=1 !Initialize ISI - 12/03/06
            !        ALLOCATE(Ckt(I)%Tube(1))
            !        ALLOCATE(Ckt(I)%Tube(1)%Seg(NumOfMods))
            !    END DO
            !    ALLOCATE(CoilSection(NumOfSections))    !RS: Debugging: Allocating so it doesn't crash later...
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
            !    DO I=1, NumOfCkts
            !        CALL GetObjectItem('ODCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
            !                            TmpNumbers,NumNumbers,Status) 
            !        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
            !            
            !        Ckt(I)%Ntube=Numbers(I)
            !    IF (ErrorFlag .NE. NOERROR) THEN 
            !        ErrorFlag=CKTFILEERROR
            !        CALL InitEvaporatorCoil_Helper_1
            !        RETURN
            !    END IF
            !        ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
            !        ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
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
                CALL InitEvaporatorCoil_Helper_1
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
                ELSEIF (I .EQ. 9) THEN
                    CALL GetObjectItem('ODCcktCircuit9_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSEIF (I .EQ. 10) THEN
                    CALL GetObjectItem('ODCcktCircuit10_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSEIF (I .EQ. 11) THEN
                    CALL GetObjectItem('ODCcktCircuit11_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                ELSE
                    CALL GetObjectItem('ODCcktCircuit12_TubeSequence',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
                END IF
                    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)                
                    DO J=1, Ckt(I)%Ntube
                        Ckt(I)%TubeSequence(J)=Numbers(J)   !Populating the tube sequence array
                    END DO 
                IF (ErrorFlag .NE. NOERROR) THEN 
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
                    RETURN
                END IF
            END DO

            !************************* Velocity Profile ********************************

            CoilSection(NumOfSections)%NumOfCkts=NumOfCkts
            IF (.NOT. ALLOCATED(CoilSection(NumofSections)%CktNum)) THEN    !RS: Debugging:
		      ALLOCATE(CoilSection(NumOfSections)%CktNum(CoilSection(NumOfSections)%NumOfCkts))
		      ALLOCATE(CoilSection(NumOfSections)%mRefIter(CoilSection(NumOfSections)%NumOfCkts))
              ALLOCATE(CoilSection(NumOfSections)%Ckt(CoilSection(NumOfSections)%NumOfCkts))    !RS: Debugging:
            END IF
		      DO J=1, NumOfCkts
		          CoilSection(NumOfSections)%CktNum(J)=J
              END DO
              
            DO I=1,2
                IF (ErrorFlag .NE. NOERROR) THEN  !Tube# ,velocity Deviation from mean value
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
                    RETURN
                END IF
            END DO

            IsUniformVelProfile=.TRUE.
            CALL GetObjectItem('ODCcktVelocityProfile',1,Alphas,NumAlphas, &
                                TmpNumbers,NumNumbers,Status)
            Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)   
            
            DO I=Nt*(Nl-1)+1,Nt*Nl !last row faces air inlet (Cross flow HX)
                DO J=1,NumOfMods
                    Tube(I)%Seg(J)%VelDev = Numbers(J)  !Bringing in the velocity deviation data
                END DO
                IF (ErrorFlag .NE. NOERROR) THEN 
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
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
                CALL InitEvaporatorCoil_Helper_1
                RETURN

            END IF

            !Propagate circuit info to coil section, ISI - 09/10/07
            NumInletSections=0
            DO I=1, NumOfSections
                IF (CoilSection(I)%IsInlet) THEN
                    NumInletSections=NumInletSections+1
                END IF
                !IF (CoilSection(I)%Ckt
                !ALLOCATE(CoilSection(I)%Ckt(CoilSection(I)%NumOfCkts))

                DO J=1, CoilSection(I)%NumOfCkts
                    CoilSection(I)%Ckt(J)=Ckt(CoilSection(I)%CktNum(J))
                END DO
            END DO

            !Determine inlet and outlet flags, split, joint or nothing, ISI - 09/10/07
            DO NumSection=1, NumOfSections
                CoilSection(NumSection)%Nnode=1

                Nnode=1
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    !Initialize
                    CoilSection(NumSection)%Ckt(I)%InJoin=0
                    CoilSection(NumSection)%Ckt(I)%InSplit=0
                    CoilSection(NumSection)%Ckt(I)%OutJoin=0
                    CoilSection(NumSection)%Ckt(I)%OutSplit=0
                    DO J=1, CoilSection(NumSection)%NumOfCkts
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(1) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube)) THEN
                            CoilSection(NumSection)%Ckt(I)%InJoin=CoilSection(NumSection)%Ckt(I)%InJoin+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(1) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(1)) THEN
                            CoilSection(NumSection)%Ckt(I)%InSplit=CoilSection(NumSection)%Ckt(I)%InSplit+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(CoilSection(NumSection)%Ckt(I)%Ntube) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube)) THEN
                            CoilSection(NumSection)%Ckt(I)%OutJoin=CoilSection(NumSection)%Ckt(I)%OutJoin+1
                        END IF
                        IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(CoilSection(NumSection)%Ckt(I)%Ntube) .EQ. &
                        CoilSection(NumSection)%Ckt(J)%TubeSequence(1)) THEN
                            CoilSection(NumSection)%Ckt(I)%OutSplit=CoilSection(NumSection)%Ckt(I)%OutSplit+1
                        END IF
                    END DO
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1 .OR. CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1) THEN
                        Nnode=Nnode+1
                    END IF
                END DO !End NumOfCkts

                CoilSection(NumSection)%Nnode=Nnode
                IF (.NOT. ALLOCATED(CoilSection(NumSection)%Node)) THEN  !RS: Debugging:
                    ALLOCATE(CoilSection(NumSection)%Node(Nnode))
                END IF

                !Find split and joint tube numbers
                J=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1 ) THEN
                        J=J+1
                        CoilSection(NumSection)%Node(J)%Num=CoilSection(NumSection)%Ckt(I)%TubeSequence(1)
                    ELSEIF (CoilSection(NumSection)%Ckt(I)%OutSplit .GT. 1) THEN
                        J=J+1
                        CoilSection(NumSection)%Node(J)%Num=CoilSection(NumSection)%Ckt(I)%TubeSequence(Ckt(I)%Ntube)
                    END IF
                    IF (J .GT. Nnode) THEN
                        EXIT
                    END IF
                END DO
                CoilSection(NumSection)%Node(Nnode)%Num=0 !section outlet 

            END DO !End NumOfSections

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

            !ISI - 09/10/07
            Tube%Empty = 1 !Initialize 
            DO NumSection=1, NumOfSections

                !Find even tubes
                DO I=1, CoilSection(NumSection)%NumOfCkts

                    !Find first and last simulation tubes
                    FirstTube=1
                    LastTube=CoilSection(NumSection)%Ckt(I)%Ntube
                    IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN
                        FirstTube=2 !Skip first tube
                    END IF 
                    IF (CoilSection(NumSection)%Ckt(I)%OutJoin .GT. 1) THEN
                        LastTube=CoilSection(NumSection)%Ckt(I)%Ntube-1 !Ignore last tube
                    END IF
                    IF (FirstTube .GT. LastTube) THEN
                        LastTube=FirstTube !Dummy tube
                    END IF

                    DO J=FirstTube, LastTube
                        TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
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
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    DO J=1, CoilSection(NumSection)%Ckt(I)%Ntube
                        TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
                        Tube(TubeNum)%Empty=0
                    END DO
                END DO

                !Number of inlet circuits
                CoilSection(NumSection)%NcktFirst=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%InJoin .LT. 1) THEN
                        CoilSection(NumSection)%NcktFirst=CoilSection(NumSection)%NcktFirst+1
                    END IF
                END DO
                IF (CoilSection(NumSection)%NcktFirst .EQ. 0) THEN
                    CoilSection(NumSection)%NcktFirst = 1 !At least one circuit, ISI - 07/28/06
                END IF

                !Number of outlet circuits
                CoilSection(NumSection)%NcktLast=0
                DO I=1, CoilSection(NumSection)%NumOfCkts
                    IF (CoilSection(NumSection)%Ckt(I)%OutSplit .EQ. 0) THEN
                        CoilSection(NumSection)%NcktLast=CoilSection(NumSection)%NcktLast+1
                    END IF
                END DO
                IF (CoilSection(NumSection)%NcktLast .EQ. 0) THEN
                    CoilSection(NumSection)%NcktLast = 1 !At least one circuit, ISI - 07/28/06
                END IF

            END DO !End NumOfSections
            
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
            
            TubeThk=TubeHeight-TubeDepth    !Or, Tube OD - Tube ID

        IF (Ltube .LE. 1e-3) THEN
            ErrorFlag=ZEROLENCOILERROR
            CALL InitEvaporatorCoil_Helper_1
            RETURN
        END IF

            SELECT CASE (Alphas(5)(1:1))
            CASE ('V','v')
                TubeOrientation=VERTICAL
            CASE ('H','h')
                TubeOrientation=HORIZONTAL
            CASE DEFAULT
                TubeOrientation=HORIZONTAL
            END SELECT
            
        NumOfMods = Numbers(14) !Number of segments or modules
        NumOfChannels = Numbers(15) !Number of circuits
        Dchannel = Numbers(16)

        !*************************** Circuiting ************************************
        
        CALL GetObjectItem('ODCcktCircuiting_Slab1',1,Alphas,NumAlphas, &
                                TmpNumbers,NumNumbers,Status)
        Numbers=DBLE(TmpNumbers)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
        
        ALLOCATE(Slab(Nl))
        !Slab#,#Passes,Tubes per Pass
        DO I=1, Nl

            NumOfPasses = Numbers(2)
            ALLOCATE(Slab(I)%Pass(NumOfPasses))
            Slab(I)%Npass=NumOfPasses
            
            DO II=1, NumOfPasses

                J=2+II
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
        !IF (LineData(1:1) .EQ. 'S' .OR. LineData(1:1) .EQ. 's') THEN !Inlet pass info

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

        !ELSE
        !
        !    NumOfInlets=1
        !    DO I=1, Nl
        !
        !        ALLOCATE(Slab(I)%InletPass(NumOfInlets))
        !        Slab(I)%Ninlet=NumOfInlets
        !
        !        DO II=1, NumOfInlets
        !
        !            Slab(I)%InletPass(II)%Ntube=Slab(I)%Pass(II)%Ntube
        !
        !        END DO
        !
        !    END DO
        !END IF

        !************************* Velocity Profile ********************************

        !Tube# ,velocity Deviation from mean value

        CALL GetObjectItem('ODCcktVelocityProfile',1,Alphas,NumAlphas, &
                                TmpNumbers,NumNumbers,Status)
        Numbers=DBLE(TmpNumbers)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
        
        IsUniformVelProfile=.TRUE.
        DO II=1,Slab(Nl)%Npass
            DO III=1,Slab(Nl)%Pass(II)%Ntube
                DO IV=1, NumOfMods
                    Slab(Nl)%Pass(II)%Tube(III)%Seg(IV)%VelDev=Numbers(IV)
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

    !CLOSE(12) !Circuit file    !RS: Debugging: Removing because the file number isn't actually attached to the file (5/21/14)

    !Suction line info
    IDsucLn=ODsucLn-SucLnThk*2
    LmodSuc=Lsucln
    AiModSuc=PI*IDsucLn*LmodSuc

    !***** Allocate pointer for suction line *****
    !ALLOCATE(SucLnSeg(NumOfMods))

    CALL InitEvaporatorCoil_Helper_1

    RETURN

    END SUBROUTINE InitEvaporatorCoil
    
SUBROUTINE InitEvaporatorStructures
    
    USE InputProcessor
    
  INTEGER, PARAMETER :: MaxNameLength = 200

  CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
  INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
  REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
  INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
  INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
  INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)


                ALLOCATE(CoilSection(NumOfSections)) 
                ALLOCATE(Ckt(NumOfCkts))
                ALLOCATE(Tube(NumOfTubes))
                ALLOCATE(Tube2D(Nl,Nt))
                ALLOCATE(JoinTubes(NumOfTubes))
                
    IF (IsCoolingMode .GT. 0) THEN
      
                !NumOfMods=2
                ALLOCATE(CoilSection(NumofSections)%Ckt(NumOfCkts)) 	  
                CoilSection(NumofSections)%NumOfCkts=NumOfCkts

                DO I=1, NumOfTubes
                    ALLOCATE(Tube(I)%Seg(NumOfMods))
                END DO

                DO I=1, Nl
                    DO J=1, Nt
                        ALLOCATE(Tube2D(I,J)%Seg(NumOfMods))
                    END DO
                END DO

                CALL GetObjectItem('IDCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
                                    TmpNumbers,NumNumbers,Status) 
                Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
        
                DO I=1, NumOfCkts
                    Ckt(I)%Ntube = Numbers(I)
                    IF (ErrorFlag .NE. NOERROR) THEN 
                        ErrorFlag=CKTFILEERROR
                        CALL InitEvaporatorCoil_Helper_1
                        RETURN
                    END IF
                    ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))          
                    ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))  
                END DO
                
                ALLOCATE(SucLnSeg(NumOfMods))
      
    ELSE
                !NumOfMods=3
                
                DO I=1, NumOfTubes
                    ALLOCATE(Tube(I)%Seg(NumOfMods))
                END DO

                DO I=1, Nl
                    DO J=1, Nt
                        ALLOCATE(Tube2D(I,J)%Seg(NumOfMods))
                    END DO
                END DO

                ALLOCATE(mRefIter(NumOfCkts))
                
                DO I=1, NumOfCkts
                    CALL GetObjectItem('ODCcktCircuiting_TubeNumbers',1,Alphas,NumAlphas, &
                                        TmpNumbers,NumNumbers,Status) 
                    Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
                        
                    Ckt(I)%Ntube=Numbers(I)
                IF (ErrorFlag .NE. NOERROR) THEN 
                    ErrorFlag=CKTFILEERROR
                    CALL InitEvaporatorCoil_Helper_1
                    RETURN
                END IF
                    ALLOCATE(Ckt(I)%Tube(Ckt(I)%Ntube))
                    ALLOCATE(Ckt(I)%TubeSequence(Ckt(I)%Ntube))
                END DO
                
                ALLOCATE(SucLnSeg(NumOfMods))
    END IF
    
    END SUBROUTINE InitEvaporatorStructures
    
!************************************************************************
    
    SUBROUTINE InitEvaporatorCoil_Helper_1

    IF (ErrorFlag .NE. NOERROR) THEN
        IF (ErrorFlag .EQ. CKTFILEERROR) THEN
            WRITE(*,*)'## ERROR ## Evaporator: Circuit file error.'
        ELSEIF (ErrorFlag .EQ. COILTUBEERROR) THEN
            WRITE(*,*)'## ERROR ## Evaporator: Coil geometry misdefined.'
            WRITE(*,*)'Tube spacing is less than tube diameter.'
        ELSEIF (ErrorFlag .EQ. COILFINERROR) THEN
            WRITE(*,*)'## ERROR ## Evaporator: Coil geometry misdefined.'
            WRITE(*,*)'Fin spacing is less than fin thickness.'
        END IF
    END IF    


    END SUBROUTINE InitEvaporatorCoil_Helper_1

!************************************************************************

SUBROUTINE EndEvaporatorCoil

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

INTEGER I,II,III !Loop counters

	  !ISI - 09/10/07
  IF (ALLOCATED(CoilSection)) THEN
      DO I=1, NumOfSections
          IF (ALLOCATED(CoilSection(I)%CktNum)) THEN
              DEALLOCATE(CoilSection(I)%CktNum)
          END IF
          IF (ALLOCATED(CoilSection(I)%Node)) THEN
              DEALLOCATE(CoilSection(I)%Node)
          END IF
          IF (ALLOCATED(CoilSection(I)%mRefIter)) THEN
              DEALLOCATE(CoilSection(I)%mRefIter)
          END IF
          IF (ALLOCATED(CoilSection(I)%Ckt)) THEN
              DEALLOCATE(CoilSection(I)%Ckt)
          END IF
      END DO
      DEALLOCATE(CoilSection)
  END IF

  IF (IsSimpleCoil .EQ. 1) THEN
	  DO I=1, NumOfCkts
		!DEALLOCATE(Ckt(I)%Tube(1)%Seg)
		DEALLOCATE(Ckt(I)%Tube)
	  END DO
	  DEALLOCATE(Ckt)
	  IF (ALLOCATED(SucLnSeg)) THEN
          DEALLOCATE(SucLnSeg) !Suction line
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
			  IF (ALLOCATED(Tube(I)%Seg)) THEN
                  DEALLOCATE(Tube(I)%Seg)
              END IF
		  END DO
	  END IF
  
	  IF (ALLOCATED(Ckt)) THEN
		  DO I=1, NumOfCkts
			  IF (ALLOCATED(Ckt(I)%Tube)) THEN
                  DEALLOCATE(Ckt(I)%Tube)
              END IF
			  IF (ALLOCATED(Ckt(I)%TubeSequence)) THEN
                  DEALLOCATE(Ckt(I)%TubeSequence)
              END IF
		  END DO
		  DEALLOCATE(Ckt)
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

  IF (ALLOCATED(SucLnSeg)) THEN
      DEALLOCATE(SucLnSeg) !Suction line
  END IF

  RETURN

END SUBROUTINE EndEvaporatorCoil

!************************************************************************

SUBROUTINE SuctionLine

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
USE CoilCalcMod

IMPLICIT NONE

INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 4=Low side interconnecting pipes
INTEGER TubeType !1=Plain; 2=General Micro Fin; 3=Herringbone; 4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A

    LmodSuc=LsucLn
    AiModSuc=PI*LmodSuc

	CoilType=LOWSIDETUBE
	TubeType=SMOOTH 
	
	MassSucLn=0

    IF (DTsucLn .NE. 0) THEN !Given suction line temperature changes
		QsucLn=mRefTot*CpgRoMod*(-DTsucLn)
	END IF

    DO K=1,1

		IF (K .EQ. 1) THEN
			!Equal to coil outlet condition
			SucLnSeg(K)%pRi=pRoCoil
			SucLnSeg(K)%hRi=hRoCoil
		ELSE !Equal to outlet of previous module(section)
			SucLnSeg(K)%pRi=SucLnSeg(K-1)%pRo
			SucLnSeg(K)%hRi=SucLnSeg(K-1)%hRo
		END IF

		pRiMod=SucLnSeg(K)%pRi
		hRiMod=SucLnSeg(K)%hRi

		CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
							 vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
							 muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		IF (xRiMod .GT. 0 .AND. xRiMod .LT. 1) THEN  !Corrected by ISI 07/09/06
		  muRiMod=xRiMod*mugRiMod+(1-xRiMod)*mufRiMod
		  vRiMod=xRiMod*vgRiMod+(1-xRiMod)*vfRiMod
		END IF
		
		hRoMod=-QsucLn/mRefTot+hRiMod   !Module Refrigerant Outlet Enthalpy

		IF (xRiMod .LT. 1) THEN
			LmodTPratio=1
		ELSE
			LmodTPratio=0
		END IF

		!Find outlet ref. pressure
		CALL CalcSegmentRefOutletPressure(CoilType,TubeType,pRiMod,hgRiMod,hfRiMod, & !CoilType,TubeType,tRiMod,pRiMod,hgRiMod,hfRiMod, &
			    	                      hRiMod,hRoMod,xRiMod,vRiMod,vgRiMod,vfRiMod,mRefTot, &
										  muRiMod,mugRiMod,mufRiMod,LmodSuc,LmodTPratio, & !muRiMod,mugRiMod,mufRiMod,SigmaMod,LmodSuc,LmodTPratio, &
										  IDsucLn,ElevSucLn,LmodSuc,pRoMod)

		pRoMod=pRoMod-AddDPSucLn    !Module Refrigerant Outlet Pressure

		CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
							 vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
							 muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		IF (xRiMod .LT. 1 .AND. xRiMod .GT. 0) THEN
			cpRiMod=0
			muRiMod=0
			kRiMod=0
		END IF

		IF (xRoMod .LT. 1 .AND. xRoMod .GT. 0) THEN
			cpRoMod=0
			muRoMod=0
			kRoMod=0
        END IF

        !Average viscosities
		mu=(muRiMod+muRoMod)/2
		muf=(mufRiMod+mufRoMod)/2
		mug=(mugRiMod+mugRoMod)/2

		CALL Inventory(CoilType,IDsucLn,mRefTot,hgRoMod,hfRoMod, & !hRiMod,hRoMod, &  !RS: Debugging: Extraneous
					   xRiMod,xRoMod,vRiMod,vRoMod,vgRimod,vfRimod,vgRomod,vfRomod, &
					   LmodSuc,LmodTP,LmodTP,MassLiqMod,MassVapMod,MassMod)
                    !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
                    !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
                    !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)
		  
		SucLnSeg(K)%Mass=MassMod
		SucLnSeg(K)%pRo=pRoMod
		SucLnSeg(K)%hRo=hRoMod
		
		MassSucLn=MassSucLn+MassMod

    END DO !End Nmod

    pRiCmp=SucLnSeg(1)%pRo  !RS Comment: Compressor Refrigerant Inlet Pressure
    hRiCmp=SucLnSeg(1)%hRo  !RS Comment: Compressor Refrigerant Inlet Enthalpy

    RETURN

END SUBROUTINE SuctionLine

!************************************************************************

SUBROUTINE RefrigerantParameters(Ref$)

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

IMPLICIT NONE

  CHARACTER*80, INTENT(IN)  :: Ref$

  RefName=Ref$

  MolWeight=MW(RefName)*1000    !RS Comment: Unit Conversion

  Tcr=Tcrit(RefName)+273.15     !RS Comment: Unit Conversion, from C to K
  Pcr=Pcrit(RefName)/1000       !RS Comment: Unit Conversion

END SUBROUTINE RefrigerantParameters

!************************************************************************

SUBROUTINE InitBoundaryConditions(CoilType)

!------------------------------------------------------------------------
!Purpose:
!To initialize segment boundary conditions
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
USE AirPropMod
USE CoilCalcMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: CoilType     !1=Condenser; 2=Evaporator; 
                                   !3=High side interconnecting pipes; 
								   !4=Low side interconnecting pipes
								   !5=Microchannel condenser
INTEGER :: NumSection !Loop counter, ISI - 09/10/07

!FLOW:

  AirPropOpt=2
  AirProp%APTDB=tAiCoil    !RS: Debugging: Formerly AirProp(1)
  AirProp%APRelHum=rhAiCoil   !RS: Debugging: Formerly AirProp(3)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
  hAiCoil=AirProp%APEnth    !RS: Debugging: Formerly AirProp(4)
  wAiCoil=AirProp%APHumRat !ISI - 08/07/06    !RS: Debugging: Formerly AirProp(2)

  !air side inlet conditions
  !RS: Replace: Moving so as to update the relative humidity before calculating
  !CPair=CPA(REAL(tAiCoil))  !RS Comment: Finding the specific heat of air   !RS: Replace: CPA (2/19/14)
  CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
  Cair=mAiCoil*CPAir    !RS Comment: Finding the capacity rate of air
  
  tAoCoil=tAiCoil
  wAoCoil=wAiCoil

  !IF (DrawBlow .EQ. BLOWTHROUGH) THEN !Blow through    !RS: Debugging: Already accounting for supply fan in E+ (9/3/14)
  !    tAiCoil=tAiCoil+PwrFan/Cair
	 ! hAiCoil=hAiCoil+PwrFan/mAiCoil
  !END IF
  IF (IsCmpInAirStream .NE. 0) THEN !Compressor in air stream
	  tAiCoil=tAiCoil+QlossCmp/Cair
	  hAiCoil=hAiCoil+QlossCmp/mAiCoil
  END IF

  AirPropOpt=1
  AirProp%APTDB=tAiCoil    !RS: Debugging: Formerly AirProp(1)
  AirProp%APEnth=hAiCoil    !RS: Debugging: Formerly AirProp(4)
  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
  rhAiCoil=AirProp%APRelHum   !RS: Debugging: Formerly AirProp(3)
  DensityIn=AirProp%APDryDens  !RS: Debugging: Formerly AirProp(7)

  !****** Coil calculation ******

  !Area calculations
  CALL CalcCoilArea(CoilType,Nl,Nt,Pt,Pl,TubeDepth, &
                    Ltube,IDtube,TubeHeight,Dchannel,NumOfChannels, &
				    FinThk,FinSpg,Lcoil,AfCoil, &
				    AoCoil,AiCoil,AmCoil)

  Pressure=pRiCoil*1000 !RS Comment: Unit Conversion
  Enthalpy=hRiCoil*1000 !RS Comment: Unit Conversion
  tRiCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)    !Coil Refrigerant Inlet Temperature
  IF (RefPropErr .GT. 0) THEN
      WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 757'
      ErrorFlag=REFPROPERROR
      RETURN
  END IF

  IF (IsSimpleCoil .EQ. 1) THEN
	  !Initialize
	  DO NumSection=1, NumOfSections !ISI - 09/10/07
	      DO I=1, CoilSection(NumSection)%NumOfCkts
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%Qmod=0
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%tAi=tAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%wbAi=wbAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%rhAi=rhAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%tAo=tAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%wbAo=wbAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%rhAo=rhAiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%pRo=pRiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%hRo=hRiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%pRi=pRiCoil
		    CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%hRi=hRiCoil
            CoilSection(NumSection)%Ckt(I)%Tube(1)%Seg%wAo=wAoCoil*.9  !RS: Debugging: Trying to initialize these values

	        CoilSection(NumSection)%Ckt(I)%pRi=pRiCoil
	        CoilSection(NumSection)%Ckt(I)%pRo=pRiCoil
	        CoilSection(NumSection)%Ckt(I)%hRi=hRiCoil
	        CoilSection(NumSection)%Ckt(I)%hRo=hRiCoil
	        CoilSection(NumSection)%Ckt(I)%mRef=mRefTot/NumOfCkts
	      END DO !End NumOfCkts
	  END DO !End NumSection 
  ELSE

	  VelAvg=(mAiCoil/DensityIn)/Aface
	  DO I=1,Nt*Nl
		  DO J=1,NumOfMods
			  Tube(I)%Seg(J)%Aface=LmodTube/(Ltube*Nt)*Aface 
			  Tube(I)%Seg(J)%mAi=mAiCoil*LmodTube/(Ltube*Nt)*Tube(I)%Seg(J)%VelDev
		  END DO
	  END DO

	  !Initialize
	  DO I=1, NumOfTubes
		  Tube(I)%Seg%Qmod=0
		  Tube(I)%Seg%tAi=tAiCoil
		  Tube(I)%Seg%rhAi=rhAiCoil
		  Tube(I)%Seg%wAi=wAiCoil
		  Tube(I)%Seg%tAo=tAiCoil
		  Tube(I)%Seg%rhAo=rhAiCoil
		  Tube(I)%Seg%wAo=wAoCoil
		  Tube(I)%Seg%pRo=pRiCoil
		  Tube(I)%Seg%hRo=hRiCoil
		  Tube(I)%Seg%pRi=pRiCoil
		  Tube(I)%Seg%hRi=hRiCoil
		  Tube(I)%ID=IDtube !ISI - 06/05/07
		  Tube(I)%NumOfMods=NumOfMods !ISI - 06/05/07
	  END DO

	  !Synchronize from tube to circuits
	  DO NumSection=1, NumOfSections !ISI - 09/10/07
	      DO I=1, CoilSection(NumSection)%NumOfCkts
		      DO J=1, CoilSection(NumSection)%Ckt(I)%Ntube
			      TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
			      CoilSection(NumSection)%Ckt(I)%Tube(J)=Tube(TubeNum)
		      END DO
		      !Initialize
		      CoilSection(NumSection)%Ckt(I)%pRi=pRiCoil
		      CoilSection(NumSection)%Ckt(I)%pRo=pRiCoil
		      CoilSection(NumSection)%Ckt(I)%hRi=hRiCoil
		      CoilSection(NumSection)%Ckt(I)%hRo=hRiCoil
	      END DO
	      
	      CoilSection(NumSection)%pRi=pRiCoil
	      CoilSection(NumSection)%hRi=hRiCoil
	      CoilSection(NumSection)%pRo=pRiCoil
	      CoilSection(NumSection)%hRo=hRiCoil
	      CoilSection(NumSection)%mRef=mRefTot/NumInletSections
	      
	      CoilSection(NumSection)%mRefIter=CoilSection(NumSection)%Ckt%mRef

	      DO I=1, CoilSection(NumSection)%NumOfCkts
      			
		    IF (CoilSection(NumSection)%Ckt(I)%OutSplit .EQ. 0) THEN !Outlet circuit
			    CoilSection(NumSection)%Ckt(I)%mRef=mRefTot/CoilSection(NumSection)%NcktLast
		    ELSEIF (CoilSection(NumSection)%Ckt(I)%InJoin .EQ. 0) THEN !Inlet circuit
		       CoilSection(NumSection)%Ckt(I)%mRef=mRefTot/CoilSection(NumSection)%NcktFirst
		    ELSEIF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN !Split inlet
  			    DO J=1, CoilSection(NumSection)%NumOfCkts
  				    IF (CoilSection(NumSection)%Ckt(I)%TubeSequence(1) .EQ. CoilSection(NumSection)%Ckt(J)%TubeSequence(Ckt(J)%Ntube)) THEN
  					    CoilSection(NumSection)%Ckt(I)%mRef=CoilSection(NumSection)%Ckt(J)%mRef/CoilSection(NumSection)%Ckt(I)%InSplit
  					    EXIT !Found the split tube
  				    END IF
			    END DO
  		    ELSE IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1) THEN !Joint inlet
  			    DO J=1, CoilSection(NumSection)%NumOfCkts
  				    IF (CoilSection(NumSection)%Ckt(J)%TubeSequence(Ckt(J)%Ntube) .EQ. CoilSection(NumSection)%Ckt(I)%TubeSequence(1)) THEN
  					    CoilSection(NumSection)%Ckt(I)%mRef=CoilSection(NumSection)%Ckt(J)%mRef*CoilSection(NumSection)%Ckt(I)%InJoin
  					    EXIT !Found the joined tube
  				    END IF
  			    END DO
  		    ELSE
			    CoilSection(NumSection)%Ckt(I)%mRef=mRefTot/CoilSection(NumSection)%NumOfCkts !to take care of one tube case, ISI - 07/28/06
  		    END IF

	      END DO
    	        
      END DO !End NumSection 

  END IF
	 
  RETURN

END SUBROUTINE InitBoundaryConditions

!************************************************************************

SUBROUTINE CalcCircuitRefInletConditions(NumSection,I,II,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate circuit refrigerant inlet condition according to circuitry 
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

IMPLICIT NONE

INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: NumSection !Section number, ISI - 09/10/07
!INTEGER,INTENT(IN) :: III !Tube number
!INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel Evaporator

REAL SumMref    !Sum of mdot ref.
REAL SumMrefHri !Sum of mdot ref x hri
INTEGER J,K !Loop counters

!FLOW:

    !ISI - 09/10/07
	IF (IsSimpleCoil .EQ. 1) THEN 
		CoilSection(NumSection)%Ckt(II)%pRi=pRiCoil
		CoilSection(NumSection)%Ckt(II)%hRi=hRiCoil
		RETURN
	END IF

	IF (CoilType .NE. MCEVAPORATOR) THEN
    !ISI - 09/10/07
		IF (CoilSection(NumSection)%Ckt(I)%InSplit .GT. 1) THEN !Split inlet
			DO J=1, NumOfCkts
				IF (CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube) .EQ. &
				    CoilSection(NumSection)%Ckt(I)%TubeSequence(1)) THEN
					IF (CoilSection(NumSection)%Ckt(J)%OutJoin .LE. 1) THEN !No joint at outlet
						!Inlet conditions
						CoilSection(NumSection)%Ckt(I)%pRi=CoilSection(NumSection)%Ckt(J)%pRo
						CoilSection(NumSection)%Ckt(I)%hRi=CoilSection(NumSection)%Ckt(J)%hRo
					ELSE !Outlet has joint
						!Inlet conditions
						CoilSection(NumSection)%Ckt(I)%pRi=CoilSection(NumSection)%Ckt(J)%Tube(Ckt(J)%Ntube-1)%Seg(NumOfMods)%pRo
						CoilSection(NumSection)%Ckt(I)%hRi=CoilSection(NumSection)%Ckt(J)%Tube(Ckt(J)%Ntube-1)%Seg(NumOfMods)%hRo
					END IF
					EXIT !Found the split tube
				END IF
			END DO
		ELSE IF (CoilSection(NumSection)%Ckt(I)%InJoin .GT. 1) THEN !Joint inlet
			K=0
			CoilSection(NumSection)%Ckt(I)%pRi=0
			CoilSection(NumSection)%Ckt(I)%hRi=0
			SumMref=0
			SumMrefHri=0
			DO J=1, CoilSection(NumSection)%NumOfCkts
				IF (CoilSection(NumSection)%Ckt(J)%TubeSequence(CoilSection(NumSection)%Ckt(J)%Ntube) .EQ. &
				    CoilSection(NumSection)%Ckt(I)%TubeSequence(1)) THEN
					K=K+1
					!Inlet conditions
					CoilSection(NumSection)%Ckt(I)%pRi=CoilSection(NumSection)%Ckt(I)%pRi+CoilSection(NumSection)%Ckt(J)%pRo
					CoilSection(NumSection)%Ckt(I)%hRi=CoilSection(NumSection)%Ckt(I)%hRi+CoilSection(NumSection)%Ckt(J)%hRo
					JoinTubes(K)=J
					SumMref=SumMref+CoilSection(NumSection)%Ckt(J)%mRef
					SumMrefHri=SumMrefHri+CoilSection(NumSection)%Ckt(J)%mRef*CoilSection(NumSection)%Ckt(J)%hRo
					IF (K .EQ. CoilSection(NumSection)%Ckt(I)%InJoin) EXIT !Found all joined tubes
				END IF
			END DO
			!Calculate according to energy balance
			CoilSection(NumSection)%Ckt(I)%pRi=CoilSection(NumSection)%Ckt(I)%pRi/CoilSection(NumSection)%Ckt(I)%InJoin !Calc inlet pressure
			CoilSection(NumSection)%Ckt(I)%hRi=SumMrefHri/SumMref       !Calc inlet enthalpy

		ELSE IF (.NOT.(CoilSection(NumSection)%IsInlet) .AND. NumSection .GT. 1) THEN !ISI - 09/10/07   !RS: Debugging: Added in >1...
		    CoilSection(NumSection)%pRi=CoilSection(NumSection-1)%pRo
		    CoilSection(NumSection)%hRi=CoilSection(NumSection-1)%hRo
		    CoilSection(NumSection)%Ckt(I)%pRi=CoilSection(NumSection-1)%pRo
		    CoilSection(NumSection)%Ckt(I)%hRi=CoilSection(NumSection-1)%hRo
		ELSE !Coil inlet
			CoilSection(NumSection)%pRi=pRiCoil
			CoilSection(NumSection)%hRi=hRiCoil
			CoilSection(NumSection)%Ckt(I)%pRi=pRiCoil
			CoilSection(NumSection)%Ckt(I)%hRi=hRiCoil
		END IF

	ELSE !Microchannel evaporator

        IF (IsParallelSlabs .GT. 0) THEN
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

SUBROUTINE CalcCoilSegment(NumSection,I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To perform heat exchanger calculation for a segment
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

INTEGER,INTENT(IN) :: NumSection !Section number
INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel Evaporator

!FLOW:

	IF (CoilType .NE. MCEVAPORATOR) THEN

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Len=LmodTube

		IF (IsSimpleCoil .EQ. 1 .AND. LmodTube .LT. SMALL) THEN     !For zero length sections
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Len=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Qmod=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%QmodSens=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRo=pRoMod
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRo=hRoMod
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAo=tAoMod
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAo=rhAoMod
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%wbAo=wbAoMod

			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hci=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%EFref=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hco=0

			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReVap=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReLiq=0

			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%cAir=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rair=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rtube=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost=0

			!Surface temperature
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tSi=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tSo=0

			RETURN
		END IF

		CALL CalcSegmentRefInletConditions(NumSection,II,II,III,IV,CoilType)
							  
		CALL CalcSegmentAirInletConditions(NumSection,II,II,III,IV,CoilType)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		mAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
		tAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAi
		rhAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAi
		VelDevMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev

		WetFlag=0
		IF (IsSimpleCoil .EQ. 1) THEN
		    RowNum=0
			IF (FinType .EQ. 4) THEN
                FinType=3 !Use regular louver fin correlation, ISI - 02/12/08
            END IF
		ELSE
			RowNum=CoilSection(NumSection)%Ckt(II)%Tube(III)%RowNum
        END IF
        
        !RS: Replace: Copying this up here so that the air enthalpy will be properly updated for AirSideCalc CALL (2/21/14)
        AirPropOpt=2
		AirProp%APTDB=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAi    !RS: Debugging: Formerly AirProp(1)
		AirProp%APRelHum=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAi   !RS: Debugging: Formerly AirProp(3)
		CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
		hAiMod=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)
        
		IF (RowNum .EQ. 0) THEN
			CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiCoil,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
							 IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil, &
                             AiCoil,FaceVel,hco,DPair,hAiMod)  
            !CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,DensityIn,DensityOut,Pt,Pl,Ltube,HtCoil, &
    !IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC, &
    !Poly1HTC,Poly2HTC,Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP, &
    !Poly1DP,Poly2DP,Poly3DP,Poly4DP,Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DPair)

		ELSE
			CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiMod,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
							 IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil, &
                             AiCoil,FaceVel,hco,DPair,hAiMod)
        END IF
        !Module surface areas
	    AoMod=AoCoil*LmodTube/Lcoil
		AfMod=AfCoil*LmodTube/Lcoil
		AiMod=AiCoil*LmodTube/Lcoil
		AmMod=AmCoil*LmodTube/Lcoil

		hco=hco*hcoMultiplier
		DPair=DPair*DPairMultiplier

		hcoMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev*hco !*LmodTube/Lcoil

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hco=hcoMod

		AirPropOpt=2
		AirProp%APTDB=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAi    !RS: Debugging: Formerly AirProp(1)
		AirProp%APRelHum=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAi   !RS: Debugging: Formerly AirProp(3)
		CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
		hAiMod=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)
		TwbAiMod=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)
		TdpAiMod=AirProp%APTDP !RS: Debugging: Formerly AirProp(6)
		wAiMod=AirProp%APHumRat   !RS: Debugging: Formerly AirProp(2)

		mRefMod=CoilSection(NumSection)%Ckt(II)%mRef
		pRiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRi
		hRiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRi

		CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
							 vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
							 muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		tAoMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAo
		wAoMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%wAo 
		CALL CalcMeanProp(tAiMod,tAoMod,tAmod)  !Module Mean Air Temperature
		CALL CalcMeanProp(wAiMod,wAoMod,wAmod)  !Module Mean Air Wet Bulb Temperature

		IF (IsSimpleCoil .EQ. 1) THEN
			IF (IV .EQ. 2) THEN
				hRiMod=hgRiMod*1.001 !Perturb a little to make sure it is in single phase region
				CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
									 vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
									 muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
				IF (ErrorFlag .GT. CONVERGEERROR) RETURN
			END IF
		END IF

		CALL CalcSegmentOutletConditions(NumSection,II,II,III,IV,CoilType)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		QmodPrev=Qmod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiMod !ISI - 12/05/06
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Len=LmodTube
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Qmod=Qmod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%QmodSens=QmodSens
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRo=pRoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRo=hRoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAo=tAoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%wAo=wAoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAo=rhAoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%wbAo=wbAoMod

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hci=hciMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%EFref=EFref
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hco=hcoMod

		IF (xRmod .GE. 1) THEN
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReVap=ReVap
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReLiq=0
		ELSE IF (xRmod .LE. 0) THEN
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReVap=0
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
		ELSE
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReVap=ReVap
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%ReLiq=ReLiq
		END IF

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%cAir=cAir
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rair=Rair
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rtube=Rtube
		
		!Surface temperature
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tSi=tAiMod-ABS(Qmod)*(Rair+CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost)
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tSo=tAoMod-ABS(Qmod)*(Rair+CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost)

		IF (IsSimpleCoil .NE. 1) THEN
            CALL UpdateTubeDataFromCircuitData(NumSection,II,III)
        END IF

	ELSE !Microchannel coil

		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Len=LmodTube
					  		              
		CALL CalcSegmentRefInletConditions(NumSection,I,II,III,IV,CoilType)
					  
		CALL CalcSegmentAirInletConditions(NumSection,I,II,III,IV,CoilType)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		mAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%mAi !*Slab(I)%Pass(II)%Ntube
		tAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi
		rhAiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi
        
        !RS: Replace: Copying this up here so that the enthalpy of air will be properly updated for AirSideCalc CALL (2/21/14)
        AirPropOpt=2
		AirProp%APTDB=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAi   !RS: Debugging: Formerly AirProp(1)
		AirProp%APRelHum=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%rhAi  !RS: Debugging: Formerly AirProp(3)
		CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
		hAiMod=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

		WetFlag=0
		RowNum=0 !Ckt(I)%Tube(J)%RowNum
		CALL AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiMod,mAiCoil,DensityIn,DensityIn,Pt,Pl,Ltube,HtCoil, &
		 				 IDtube,ODtube,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg,Lcoil,AfCoil,AoCoil, &
                         AiCoil,FaceVel,hco,DPair,hAiMod)
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
		TwbAiMod=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)
		TdpAiMod=AirProp%APTDP !RS: Debugging: Formerly AirProp(6)
		wAiMod=AirProp%APHumRat   !RS: Debugging: Formerly AirProp(2)

		pRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRi
		hRiMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRi

		CALL CalcRefProperty(pRiMod,hRiMod,hfRiMod,hgRiMod,hfgRiMod,Psat,Tsat,tRiMod,xRiMod, &
						   vRiMod,vfRiMod,vgRiMod,cpRiMod,cpfRiMod,cpgRiMod, &
						   muRiMod,mufRiMod,mugRiMod,kRiMod,kfRiMod,kgRiMod,SigmaMod)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		tAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo !ISI - 12/25/06
		wAoMod=Slab(I)%Pass(II)%Tube(III)%Seg(IV)%wAo !ISI - 12/25/06
		CALL CalcMeanProp(tAiMod,tAoMod,tAmod)  !Module Mean Air Temperature
		CALL CalcMeanProp(wAiMod,wAoMod,wAmod)  !Module Mean Air Wet Bulb Temperature

		CALL CalcSegmentOutletConditions(I,I,II,III,IV,CoilType)
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF
	
		QmodPrev=Qmod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod=Qmod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%QmodSens=QmodSens
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%pRo=pRoMod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hRo=hRoMod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%tAo=tAoMod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%wAo=wAoMod
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

SUBROUTINE CalcSegmentAirInletConditions(NumSection,I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate inlet air temp and relative humidity
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

IMPLICIT NONE

INTEGER,INTENT(IN) :: NumSection !Section number
INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
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
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAi=tAiCoil
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAi=rhAiCoil
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev=1
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiCoil*Lmodtube/Lcoil !ISI - 12/05/06
		mAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi

		RETURN
	END IF

	IF (CoilType .NE. MCEVAPORATOR) THEN
		IF (Tube(TubeNum)%Fup .NE. 0) THEN !Upper front tubes
			IF (Tube(TubeNum)%Even .EQ. 0) THEN !Odd tubes
				IF (Tube(Tube(TubeNum)%Fup)%Empty .NE. 0) THEN !Empty tubes

					IF (Tube(TubeNum)%Fdown .EQ. 0) THEN !Coil bottom, ISI - 07/29/07

						IF (Tube(Tube(TubeNum)%Fup)%Fdown .EQ. 0) THEN
						    tAoFdown=tAiCoil
						    rhAoFdown=rhAiCoil
						    tAiFdown=tAiCoil
						    rhAiFdown=rhAiCoil
						    mAiFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
						    VelDevFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev						
                        ELSE
						    tAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAo
						    rhAoFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAo
						    tAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%tAi
						    rhAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%rhAi
						    mAiFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%mAi
						    VelDevFdown=Tube(Tube(Tube(TubeNum)%Fup)%Fdown)%Seg(IV)%VelDev 
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
                        mAiFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
			            VelDevFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev						
						
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
						    mAiFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
						    VelDevFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev						
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
                        mAiFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
			            VelDevFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev						
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
			mAiFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
			VelDevFup=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev
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
						mAiFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
						VelDevFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev
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
						mAiFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
						VelDevFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev
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
			mAiFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
			VelDevFdown=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev
        END IF

        !Average values
		tAiFavg=(tAiFup+tAiFdown)/2
		tAoFavg=(tAoFup+tAoFdown)/2

 		rhAiFavg=(rhAiFup+rhAiFdown)/2
		rhAoFavg=(rhAoFup+rhAoFdown)/2

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%tAi=tAiFavg-1*(tAiFavg-tAoFavg)
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%rhAi=rhAiFavg-1*(rhAiFavg-rhAoFavg)

		IF (CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev .LE. 0) THEN
            CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev=1
        END IF

		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi=(mAiFup+mAiFdown)/2
		mAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi=mAiCoil*LmodTube/(Ltube*Nt)*CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%VelDev
		mAiMod=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%mAi

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

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
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
REAL, INTENT(OUT) :: SigmaRef !Surface tension, N/m

!LOCAL VARIABLES
REAL Wlocal !Local oil mass fraction

!FLOW:

	Pressure=pRef*1000  !RS Comment: Unit Conversion
	Enthalpy=hRef*1000  !RS Comment: Unit Conversion

	tRef=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Refrigerant Temperature
    IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3144'
	    ErrorFlag=REFPROPERROR
		RETURN
    END IF
	
	xRef=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Refrigerant Quality
    IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3151'
	    ErrorFlag=REFPROPERROR
		RETURN
    END IF
	
	vRef=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr) !Refrigerant Specific Volume
    IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3158'
	    ErrorFlag=REFPROPERROR
		RETURN
    END IF	
	vRef=1/vRef

	cpRef=PH(RefName, Pressure, Enthalpy, 'specificheat', RefrigIndex,RefPropErr)   !Refrigerant Specific Volume
    IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3166'
	    ErrorFlag=REFPROPERROR
		RETURN
    END IF
	cpRef=cpRef/1000    !RS Comment: Unit Conversion
	
	muRef=PH(RefName, Pressure, Enthalpy, 'viscosity', RefrigIndex,RefPropErr)  !Refrigerant Dynamic Viscosity
    IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3174'
	    ErrorFlag=REFPROPERROR
		RETURN
    END IF
	
	kRef=PH(RefName, Pressure, Enthalpy, 'conductivity', RefrigIndex,RefPropErr)    !Refrigerant Thermal Conductivity
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3181'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	kRef=kRef/1000  !RS Comment: Unit Conversion

	Temperature=tRef
	Quality=1
	IF (tRef+273.15 .GT. Tcr .OR. tRef+273.15 .LT. 0) THEN
		Psat=pRef
	ELSE 
		Psat=TQ(RefName, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr)  !Saturation Pressure
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3194'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF
		Psat=Psat/1000  !RS Comment: Unit Conversion
	END IF

	SigmaRef=PQ(RefName, Pressure, Quality, 'surfacetension', RefrigIndex,RefPropErr)   !Refrigerant Surface Tension
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3820'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF

	Pressure=pRef*1000  !RS Comment: Unit Conversion
	Quality=0
	hfRef=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Refrigerant Liquid Enthalpy
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3205'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	hfRef=hfRef/1000    !RS Comment: Unit Conversion
	
	cpfRef=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Refrigerant Liquid Specific Heat
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3213'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	cpfRef=cpfRef/1000  !RS Comment: Unit Conversion

	mufRef=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Refrigerant Liquid Dynamic Viscosity
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3221'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	
	kfRef=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Refrigerant Liquid Thermal Conductivity
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3228'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	kfRef=kfRef/1000    !RS Comment: Unit Conversion

	vfRef=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3236'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	vfRef=1/vfRef   !Refrigerant Liquid Specific Volume

	Pressure=pRef*1000  !RS Comment: Unit Conversion
	Quality=1
	hgRef=PQ(RefName, Pressure, Quality, 'enthalpy', RefrigIndex,RefPropErr)    !Refrigerant Vapor Enthalpy
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3246'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	hgRef=hgRef/1000    !RS Comment: Unit Conversion
	
	cpgRef=PQ(RefName, Pressure, Quality, 'specificheat', RefrigIndex,RefPropErr)   !Refrigerant Vapor Specific Heat
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3254'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	cpgRef=cpgRef/1000  !RS Comment: Unit Conversion
	
	mugRef=PQ(RefName, Pressure, Quality, 'viscosity', RefrigIndex,RefPropErr)  !Refrigerant Vapor Dynamic Viscosity
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3262'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	
	kgRef=PQ(RefName, Pressure, Quality, 'conductivity', RefrigIndex,RefPropErr)    !Refrigerant Vapor Thermal Conductivity
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3269'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	kgRef=kgRef/1000    !RS Comment: Unit Conversion

	vgRef=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3277'
		ErrorFlag=REFPROPERROR
		RETURN
	END IF
	vgRef=1/vgRef   !Refrigerant Vapor Specific Volume

	Tsat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Pressure
	IF (RefPropErr .GT. 0) THEN
		WRITE(*,*)'-- WARNING -- Condenser: Refprop error. Line 3051'
		ErrorFlag=REFPROPERROR
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

SUBROUTINE CalcSegmentRefInletConditions(NumSection,I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate inlet refrigerant pressure and enthalpy
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!November 2005
!
!Reference:
!none
!
!------------------------------------------------------------------------

IMPLICIT NONE

INTEGER,INTENT(IN) :: NumSection !Section number
INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel Evaporator

!FLOW:

	IF (CoilType .NE. MCEVAPORATOR) THEN 

		IF (III .EQ. 1 .AND. IV .EQ. 1) THEN !Equal to circuit inlet
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRi=CoilSection(NumSection)%Ckt(II)%pRi
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRi=CoilSection(NumSection)%Ckt(II)%hRi

		ELSE IF (IV .EQ. 1) THEN !Equal to outlet of previous tube Changed it to IV from K
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRi=CoilSection(NumSection)%Ckt(II)%Tube(III-1)%Seg(NumOfMods)%pRo
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRi=CoilSection(NumSection)%Ckt(II)%Tube(III-1)%Seg(NumOfMods)%hRo

		ELSE !Equal to outlet of previous module(section)
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%pRi=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV-1)%pRo
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hRi=CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV-1)%hRo
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

SUBROUTINE CalcSegmentOutletConditions(NumSection,I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate segment outlet conditions
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
USE OilMixtureMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: NumSection !Section number
INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
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
INTEGER RefBCiter   !Iteration loop counter
LOGICAL IsTransitSegmentCalled !Flag to indicate if 'CalcTransitSegment' is called
LOGICAL IsTransitionSegment !Flag to indicate if it is transtion segment
REAL humrat !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio

!FLOW:

	!Initialize for property iteration, to find the mean property
	hfgRoMod=0;  xRoMod=0;  vgRoMod=0;  vfRoMod=0
	muRoMod=0;  mugRoMod=0;  mufRoMod=0
	kRoMod=0;	  kfRoMod=0;  kgRoMod=0
	cpRoMod=0;  cpfRoMod=0;  cpgRoMod=0
	DPmod=0;

	IsTransitSegmentCalled=.FALSE.
	IsTransitionSegment=.FALSE.
	
	PrevpRoMod=BIG 
	PrevhRoMod=BIG 

	tRoMod=tRiMod
	
	xRoMod=xRiMod

	hRoMod=hRiMod

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

        !Correct dynamic viscosity
		IF (muRmod .LE. 0) THEN !ISI - 08/03/06
		    IF (xRmod .LE. 0) THEN
                muRmod = mufRmod
            END IF
			IF (xRmod .GE. 1) THEN
                muRmod = mugRmod
            END IF
		END IF

		!For tube cover both two phase and single phase region
		LmodTPratio=0 !Initialize
		QmodTP=0 !Initialize
		IF (RefBCiter .GT. 1 .AND. hRiMod .LT. hgRoMod .AND. hRoMod .GE. hgRoMod) THEN 
		
			CALL CalcTransitionSegment(CoilType)
			IF (IsSimpleCoil .EQ. 1) THEN
                IsTransitionSegment=.TRUE.
            END IF
			IsTransitSegmentCalled=.TRUE.
			IF (ErrorFlag .GT. CONVERGEERROR) THEN
                RETURN
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

		IF (.NOT. IsTransitionSegment) THEN
			IF (DTmod .EQ. 0) THEN
                DTmod=(tAiMod+tRiMod)/2 !First estimate
            END IF
			!CALL hcRefside(CoilType,TubeType,IDtube,ktube,mRefMod,-Qmod,AoMod,AiMod,hfgRmod,xRmod,xRmod, &
			!			   vgRmod,vfRmod,muRmod,mugRmod,mufRmod, &
		 ! 				   kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
			!			   MolWeight,Psat,Pcr,Tsat,SigmaMod,DTmod,Wabsolute,EFref,hciMod)
            CALL hcRefside(CoilType,TubeType,IDtube,mRefMod,-Qmod, &               !Calculating the refrigerant side heat transfer coefficient
        xRmod,xRmod,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
        Psat,Pcr,Wabsolute,EFref,hciMod)
            
			hciMod=hciMod*hciMultiplier

			CALL Reynolds(IDtube,mRefMod,xRmod,muRmod,mugRmod,mufRmod,ReVap,ReLiq)
			  
			!***Dry surface calc.
			WetFlag=0
			Tsurf=0

			!Calc. UA
			CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
						hcoMod,hciMod,AfMod,AoMod,AiMod,AmMod, &
						UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)
            !CALL CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth,RowNum,tAiMod,hAiMod, &
            !hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
            !UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)

			IF (xRiMod .LT. 1 .AND. xRoMod .GE. 1 .AND. LmodTPratio .LT. 1) THEN !Evaporator outlet
				UA=UA*(1-LmodTPratio)
			END IF
			 
			!Calc. Cref
			IF (CoilType .NE. MCEVAPORATOR) THEN
				cRef=mRefMod*cpRmod
			ELSE !Microchannel coil
				cRef=mRefMod*cpRmod*NumOfChannels 
			END IF
			IF (xRmod .LT. 1. .AND. xRmod .GT. 0.) THEN
                cRef=BIG !Phase change
            END IF

			!Calc. Cair
			!CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
            CALL CalcMeanProp(TAiCoil,TAoCoil,TAmod)    !RS Comment: Mean Air Coil Temperature
            CALL CalcMeanProp(hAiCoil,hAoCoil,hAmod)    !RS Comment: Mean Air Coil Humidity
            humrat=HUMTH(tAmod,hAmod)   !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
            CPair=CPAirFunction(tAmod,humrat)  !RS: Replace: CPA (2/19/14)
			cAir=mAiMod*cpAir

			!Calc. Cmin
			Cmin=MIN(cAir,cRef)

			!Calc. Epsilon
			CALL EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)
            
			!Calc. DT
			IF (LmodTPratio .GT. 0) THEN !ISI - 07/21/06
				DT=(tRmod-tAiMod)
			ELSE
				DT=(tRiMod-tAiMod)
			END IF

			!Calc. dry module heat transfer
			QmodDry=EPS*Cmin*DT

			IF (xRiMod .LT. 1 .AND. xRoMod .GE. 1) THEN
				IF (IsSimpleCoil .EQ. 1) THEN
					IF (QmodTP .NE. 0) Qmod = QmodTP
				ELSE
					IF (LmodTP .EQ. LmodTube) THEN
						IF (QmodDry .GT. QmodTP) THEN
                            QmodDry = QmodTP
                        END IF
					ELSE
						QmodDry=QmodDry+QmodTP
					END IF
				END IF
			END IF

			!Include solar radiation
			IF (CoilType .NE. MCEVAPORATOR) THEN
				IF (IsCoolingMode .GT. 0 .AND. &
					CoilSection(NumSection)%Ckt(II)%Tube(III)%Fup .EQ. 0 .AND. &
					CoilSection(NumSection)%Ckt(II)%Tube(III)%Fdown .EQ. 0) THEN
					Rtot=Rair*AoMod+Rrefrig*AiMod+Rtube*AmMod
					Qsolar=Rair*AoMod/Rtot*CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Aface* &
					       SurfAbsorptivity*SolarFlux
					QmodDry=QmodDry-Qsolar
	  			END IF
			END IF

			!Calc. Outside air enthalpy
			hAoDry=QmodDry/mAiMod+hAiMod

			!Outside air temp
			TdbAoDry=QmodDry/cAir+tAiMod

			AirPropOpt=1
			AirProp%APTDB=TdbAoDry !RS: Debugging: Formerly AirProp(1)
			AirProp%APEnth=hAoDry   !RS: Debugging: Formerly AirProp(4)
			CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
			rhAoMod=AirProp%APRelHum  !RS: Debugging: Formerly AirProp(3)
			TwbAoMod=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)
			TdpAoMod=AirProp%APTDP !RS: Debugging: Formerly AirProp(6)
			DensityOut=AirProp%APDryDens   !RS: Debugging: Formerly AirProp(7)

			!Calc dry surface temperature
			NTUsDry=1/(cAir*Rair) 
			EPSsDry=1-EXP(-NTUsDry)

			TsDry=QmodDry/(EPSsDry*cAir)+tAiMod

			DryWet=3
			Qmod=QmodDry
			QmodSens=QmodDry
			QmodWet=0
			hAoMod=hAoDry
			tAoMod=TdbAoDry
			SHR=1
						 
		END IF

		IF (xRiMod .LT. 1 .AND. xRoMod .GE. 1) THEN
			IF (IsSimpleCoil .EQ. 1) THEN
				IF (QmodTP .NE. 0) THEN
                    Qmod = QmodTP
                END IF
			END IF
		END IF

		IF (xRmod .LT. 1 .AND. xRmod .GT. 0 .AND. NOT(IsTransitSegmentCalled)) THEN 
							  						
			!********************* Ding starts ******************************
			!Calc. temperature where moisture remove occurs

			IF (TdpAiMod .GT. tRiMod) THEN !Wet surface
				  !CALL CalcWetSurfaceDing(I,II,III,IV,CoilType)
				  CALL CalcWetSurfaceBraun(NumSection,I,II,III,IV,CoilType)
				  !CALL CalcWetSurfaceMcQuiston(I,II,III,IV,CoilType)
			END IF

		END IF

		IF (CoilType .NE. MCEVAPORATOR) THEN
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%DryWet = DryWet
		ELSE
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%DryWet = DryWet
		END IF

		AirPropOpt=1
		AirProp%APTDB=tAoMod   !RS: Debugging: Formerly AirProp(1)
		AirProp%APEnth=hAoMod   !RS: Debugging: Formerly AirProp(4)
		CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
		rhAoMod=AirProp%APRelHum  !RS: Debugging: Formerly AirProp(3)

		IF (CoilType .NE. MCEVAPORATOR) THEN
			hRoMod=-Qmod/mRefMod+hRiMod 
		ELSE !Microchannel coil
			hRoMod=-(Qmod/NumOfChannels)/mRefMod+hRiMod  
		END IF

		!To prevent enthalpy is out of range - ISI - 12/27/06
		IF (IsSimpleCoil .EQ. 1 .AND. RefBCiter .EQ. 1 .AND. hRiMod .LT. hgRoMod .AND. hRoMod .GE. hgRoMod) THEN
            CYCLE
        END IF

		CALL CalcSegmentRefOutletPressure(CoilType,TubeType,pRiMod,hgRiMod,hfRiMod, & !CoilType,TubeType,tRiMod,pRiMod,hgRiMod,hfRiMod, &
				  	                      hRiMod,hRoMod,xRiMod,vRiMod,vgRiMod,vfRiMod,mRefMod, &
										  muRiMod,mugRiMod,mufRiMod,LmodTube,LmodTPratio, & !muRiMod,mugRiMod,mufRiMod,SigmaMod,LmodSuc,LmodTPratio, &
										  IDtube,HtCoil,Lcoil,pRoMod)

		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		CALL CalcRefProperty(pRoMod,hRoMod,hfRoMod,hgRoMod,hfgRoMod,Psat,Tsat,tRoMod,xRoMod, &
							 vRoMod,vfRoMod,vgRoMod,cpRoMod,cpfRoMod,cpgRoMod, &
							 muRoMod,mufRoMod,mugRoMod,kRoMod,kfRoMod,kgRoMod,SigmaMod)
        
		IF (ErrorFlag .GT. CONVERGEERROR) THEN
            RETURN
        END IF

		IF (CoilType .NE. MCEVAPORATOR .AND. IsSimpleCoil .NE. 1) THEN
			!Return bend pressure drop
			IF (IV .EQ. NumOfMods) THEN
				IF (II .EQ. LastTube) THEN
		  			IF (CoilSection(NumSection)%Ckt(II)%OutSplit .GT. 1 .OR. &
		  			    CoilSection(NumSection)%Ckt(II)%OutJoin .GT. 1) THEN
						CALL returnbend(CoilType,TubeType,IDtube,Pt,mRefmod,xRoMod,vRoMod,vgRoMod,vfRoMod,mugRoMod,mufRoMod,DPreturnbend)
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

		IF (xRmod .GT. 0.9 .AND. RefBCiter .GE. 2) THEN
		  EXIT !doesn't converge for xRmod > 0.9.  May be due to the sharp change of the hci
		ELSE
			DTmod=-Qmod*(1/(hciMod*AiMod)+LOG(ODtube/IDtube)/(2*PI*Ktube*LmodTube))
			DiffpRoMod=ABS((pRoMod-PrevpRoMod)/PrevpRoMod)
			DiffhRoMod=ABS((hRoMod-PrevhRoMod)/PrevhRoMod)
			IF (DiffpRoMod .GT. SMALL .OR. DiffhRoMod .GT. SMALL) THEN 
				PrevpRoMod=pRoMod
				PrevhRoMod=hRoMod
			ELSE 
				EXIT
			END IF
		END IF
			
	END DO !end of RefBCiter

	IF (RefBCiter .GT. RefBCmaxIter) THEN
		!RefBCiter not converged.
		ErrorFlag=CONVERGEERROR
	END IF

	!Outside air temp
	tAoMod=QmodSens/cAir+tAiMod

	!Calc. Outside air enthalpy
	hAoMod=Qmod/mAiMod+hAiMod

	AirPropOpt=1
	AirProp%APTDB=tAoMod   !RS: Debugging: Formerly AirProp(1)
	AirProp%APEnth=hAoMod   !RS: Debugging: Formerly AirProp(4)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,
	rhAoMod=AirProp%APRelHum  !RS: Debugging: Formerly AirProp(3)
	wbAoMod=AirProp%APTWB  !RS: Debugging: Formerly AirProp(5)
	wAoMod=AirProp%APHumRat   !RS: Debugging: Formerly AirProp(2)
	 
RETURN

END SUBROUTINE CalcSegmentOutletConditions

!************************************************************************

SUBROUTINE CalcSegmentRefOutletPressure(CoilType,TubeType,pRi,hgRi,hfRi, & !CoilType,TubeType,tRi,pRi,hgRi,hfRi, &
			  	                        hRi,hRo,xRi,vRi,vgRi,vfRi,mRef, &
										muRi,mugRi,mufRi,Lsegment,LmodTPratio, & !muRi,mugRi,mufRi,Sigma,Lsegment,LmodTPratio, &
										IDtube,Elevation,Ltotal,pRo)


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

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
USE CoilCalcMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
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
REAL, INTENT(OUT) :: pRo       !Outlet pressure, kPa

REAL tRo       !Outlet temperature, C
REAL vgRo      !Outlet vapor specific volume, m^3/kg
REAL vfRo      !Outlet liquid specific volume, m^3/kg
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
		xRo=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !Refrigerant Outlet Quality
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3391'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF

	    tRo=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !Refrigerant Outlet Temperature
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3398'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF

		vRo=PH(RefName, Pressure, Enthalpy, 'density', RefrigIndex,RefPropErr)
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3405'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF
		vRo=1/vRo   !Refrigerant Outlet Specific Volume

		Pressure=pRo*1000   !RS Comment: Unit Conversion
		Quality=1
		vgRo=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3394'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF
		vgRo=1/vgRo !Refrigerant Outlet Vapor Specific Volume

		Quality=0
		vfRo=PQ(RefName, Pressure, Quality, 'density', RefrigIndex,RefPropErr)
		IF (RefPropErr .GT. 0) THEN
			WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 3403'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF
		vfRo=1/vfRo !Refrigerant Outlet Liquid Specific Volume
  
		CALL MODdP(CoilType,TubeType,hgRi,hfRi, &                                       !Calculates the module pressure drop
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

SUBROUTINE CalcTransitionSegment(CoilType)

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

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
USE CoilCalcMod
USE AirPropMod
USE OilMixtureMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel evaporator

!FLOW:

	xRmod=(xRiMod+1)/2 !Set it to first portion of the transition element

    hciMod=hciMod*hciMultiplier

	IF (CoilType .NE. MCEVAPORATOR) THEN
		QmodTP=mRefMod*(hRiMod-hgRiMod)
	ELSE
		QmodTP=mRefMod*NumOfChannels*(hRiMod-hgRiMod)
	END IF

	!Calc. Cref
	cRef=1.0E20 !Phase change 

	!Calc. DT
	DT=(tRiMod-tAiMod)

	!Find transition boundary
	CALL FindTransitionBoundary(CoilType)

	IF (LmodTP .NE. LmodTube) THEN
	    xRmod=1
	ELSE
		IF (CoilType .NE. MCEVAPORATOR) THEN
			hRoMod=-QmodTP/mRefMod+hRiMod !ISI - 06/18/05
		ELSE
			hRoMod=-(QmodTP/NumOfChannels)/mRefMod+hRiMod 
		END IF
		 
		Pressure=pRiMod*1000    !RS Comment: Unit Conversion
		Enthalpy=hRoMod*1000    !RS Comment: Unit Conversion

	    xRmod=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)    !Module Refrigerant Quality
		IF (RefPropErr .GT. 0) THEN
		    WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 1291'
			ErrorFlag=REFPROPERROR
			RETURN
		END IF
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
	tRmod=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !Module Refrigerant Temperature

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

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel evaporator

!FLOW:

	!Find module length in two-phase 
	!Initialize
	LmodTPmin=1e-6 !ISI - 13/25/06
	LmodTPmax=LmodTube
	LmodTP=LmodTube/2
	DTmod=0 
	DO
		LmodTPratio=LmodTP/LmodTube
								
		!***Dry surface calc.
		WetFlag=0
		Tsurf=0

		!ISI - 09/11/06
		IF (DTmod .EQ. 0) THEN
            DTmod=(tAiMod+tRiMod)/2 !First estimate
        END IF
		!CALL hcRefside(CoilType,TubeType,IDtube,ktube,mRefMod,-Qmod,AoMod,AiMod,hfgRmod,xRiMod,1.0, &   !Calculates the refrigerant side heat transfer coefficient
		!			   vgRmod,vfRmod,muRmod,mugRmod,mufRmod, &
		!	  		   kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
		!			   MolWeight,Psat,Pcr,Tsat,SigmaMod,DTmod,Wabsolute,EFref,hciMod)
        CALL hcRefside(CoilType,TubeType,IDtube,mRefMod,-Qmod, &               !Calculating the refrigerant side heat transfer coefficient
        xRimod,1.0,vgRmod,vfRmod,muRmod,mugRmod,mufRmod,kRmod,kfRmod,kgRmod,cpRmod,cpfRmod,cpgRmod, &
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
		CPair=CPA(REAL(tAmod))
		Cair=mAiMod*CPair

		Cmin=MIN(Cair,Cref) !ISI - 12/25/06

		!Calc. Epsilon
		CALL EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)
        
		!Calc. dry module heat transfer
		QmodDry=EPS*Cmin*DT
		Qmod=QmodDry
		IF (IsSimpleCoil .EQ. 1 .AND. DryWet .EQ. 1) THEN
			!CALL CalcWetSurfaceDing(1,1,1,1,CoilType)
			CALL CalcWetSurfaceBraun(1,1,1,1,1,CoilType)
			!CALL CalcWetSurfaceMcQuiston(1,1,1,1,CoilType)
		END IF
              
		DTmod=-Qmod*(1/(hciMod*AiMod*LmodTP/LmodTube)+LOG(ODtube/IDtube/(2*PI*Ktube*LmodTP)))
		IF (ABS(LmodTPratio-1) .LT. SMALL) THEN
			QmodTP=Qmod
			LmodTP=LmodTube
			LmodTPratio=0
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
			  
	END DO !End of LmodTP calc.

RETURN

END SUBROUTINE FindTransitionBoundary

!************************************************************************

SUBROUTINE CalcWetSurfaceDing(I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate wet surface air side heat transfer
!
!Reference:
!  Ding, X., Eppe, J.P., Lebrun, J. and Wasacz, M. (1990). Cooling coil
!    models to be used in transient and/or wet regimes. Theoretical analysis
!    and experimental validation. In Proceedings of System Simulation in
!    Buildings '90, Liege, Belgium, December 1990.
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
USE CoilCalcMod
USE AirPropMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel evaporator

!FLOW:

	!Initialize fictitious air specific heat
	hAoWet=hAoDry
	TwbAoWet=TwbAoMod
	cfprev=cpAir
	cf=0

	DO cfIter=1, WetSurfaceMaxIter

	    !Fictitious air specific heat
		cf=(hAiMod-hAoWet)/(TwbAiMod-TwbAoWet)
								   
		IF (ABS((cf-cfprev)/cfprev) .GT. SMALL .AND. cf .GT. 0 .OR. cfIter .EQ. 1) THEN
		    !Fictitious convection resistance
			hcof=(hcoMod*cf)/cPAir
			WetFlag=1
            
			IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN  !Calculates the UA
				CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
						    hcof,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio,UAf,Rairf,Rrefrig,Rtube,FinEff,SurfEff)
            !CALL CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth,RowNum,tAiMod,hAiMod, &
            !hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
            !UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)

			ELSE
				CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
						    hcof,hciMod,AfMod,AoMod,AiMod,AmMod,UAf,Rairf,Rrefrig,Rtube,FinEff,SurfEff)
			END IF
							 
			cAirf=mAiMod*cf
										 
			CALL EPScalc(cAirf,cRef,UAf,Cratio,NTU,EPSf)						 
			Cminf=MIN(cAirf,cRef)
										 
			DT=(tRiMod-TwbAiMod)
 										 
  			QmodWet=EPSf*Cminf*DT

			IF (CoilType .NE. MCEVAPORATOR) THEN
				!Include solar radiation
				IF (.NOT.(IsCoolingMode) .AND. &
					Ckt(II)%Tube(III)%Fup .EQ. 0 .AND. Ckt(II)%Tube(III)%Fdown .EQ. 0) THEN
					Rtot=Rairf*AoMod+Rrefrig*AiMod+Rtube*AmMod
					Qsolar=Rairf*AoMod/Rtot*Ckt(II)%Tube(III)%Seg(IV)%Aface*SurfAbsorptivity*SolarFlux
					QmodWet=QmodWet-Qsolar
	  			END IF
			END IF
				   
			!Calc wet surface temperature
			IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
				NTUsWet=SurfEff*AoMod*LmodTPratio*hcof/cAirf
			ELSE
				NTUsWet=SurfEff*AoMod*hcof/cAirf
			END IF

			EPSsWet=1-EXP(-NTUsWet)
			TsWet=QmodWet/(EPSsWet*cAirf)+TwbAiMod

			!Outlet air dry bulb temperature
			TdbAoWet=EPSsWet*(TsWet-tAiMod)+tAiMod

			!Outlet air wet bulb temp
			TwbAoWet=QmodWet/cAirf+TwbAiMod

			!Outlet air enthalpy
			AirPropOpt=3
			AirProp%APTDB=TdbAoWet !RS: Debugging: Formerly AirProp(1)
			AirProp%APTWB=TwbAoWet !RS: Debugging: Formerly AirProp(5)
			CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
			hAoWet=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

			!Update fictitious air specific heat
			cfprev=cf

		ELSE
		    !Calc wet surface temperature
			IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
				NTUsWet=SurfEff*AoMod*LmodTPratio*hcof/cAirf
			ELSE
				NTUsWet=SurfEff*AoMod*hcof/cAirf
			END IF
			EPSsWet=1-EXP(-NTUsWet)
			TsWet=QmodWet/(EPSsWet*cAirf)+TwbAiMod

			!Outlet air dry bulb temperature
			TdbAoWet=EPSsWet*(TsWet-tAiMod)+tAiMod

			TwbAoWet=QmodWet/cAirf+TwbAiMod

			!Outlet air enthalpy
			AirPropOpt=3
			AirProp%APTDB=TdbAoWet !RS: Debugging: Formerly AirProp(1)
			AirProp%APTWB=TwbAoWet !RS: Debugging: Formerly AirProp(5)
			CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
			hAoWet=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

			DryWet=1
			Qmod=QmodWet
			QmodSens=cAir*(TdbAoWet-tAiMod)
			SHR=QmodSens/Qmod
			hAoMod=hAoWet
			tAoMod=TdbAoWet

            !Populating the circuit or slab arrays
			IF (CoilType .NE. MCEVAPORATOR) THEN
				Ckt(II)%Tube(III)%Seg(IV)%hco = hcof
				Ckt(II)%Tube(III)%Seg(IV)%cAir=cAirf
				Ckt(II)%Tube(III)%Seg(IV)%Rair=Rairf
				Ckt(II)%Tube(III)%Seg(IV)%Rtube=Rtube
				Ckt(II)%Tube(III)%Seg(IV)%SHR=SHR
			ELSE
				Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco = hcof
				Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir=cAirf
				Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=Rairf
				Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=Rtube
				Slab(I)%Pass(II)%Tube(III)%Seg(IV)%SHR=SHR
			END IF

			EXIT
		END IF
	END DO !End of cfIter

	IF (cfIter .GT. WetSurfaceMaxIter) THEN
        WRITE(*,*)'-- WARNING -- Evaporator: Fictitious air specific heat not converged'
    END IF

		!********************* Ding ends ******************************

RETURN

END SUBROUTINE CalcWetSurfaceDing

!************************************************************************

SUBROUTINE CalcWetSurfaceBraun(NumSection,I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate wet surface air side heat transfer
!
!Reference:
!Harms, T.M.; Groll, E.A. and Braun, J.E. (2003). Accurate charge inventory  
!modeling for unitary air conditioners. HVAC&R Research, 9(1), pp. 55-78.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!July 2006
!
!Reference:
!none
!
!------------------------------------------------------------------------

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
USE CoilCalcMod
USE AirPropMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: NumSection !Section number
INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel evaporator

REAL, PARAMETER :: DeltaTemp = 0.01 !Small temperature change, C

REAL hRiSat      !Saturation air enthalpy based on refrigerant temperature, kJ/kg
REAL DeltahRiSat !Small change of hRiSat, kJ/kg
REAL CpSat       !Saturation specific heat, kJ/kg
REAL CpMoist     !Specific heat of moist air, kJ/kg
REAL Rcollar     !Tube radius with fin collar, m
REAL XM          !Heat exchanger geometric parameter, m
REAL XL          !Heat exchanger geometric parameter, m
REAL Req         !Equivalent fin radius, m
REAL mm          !Parameter for fin efficiency calculation
REAL DH          !Enthalpy difference, kJ/kg
REAL NTUwet      !Wet surface NTU
REAL hAoSatSurf  !Outlet air surface saturation enthalpy, kJ/kg
REAL tAoSat      !Outlet air saturation temperature, C
REAL FinEffwet   !Wet surface fin efficiency
REAL SurfEffwet  !Wet surface efficiency
REAL UAwet       !Wet surface UA
REAL EPSwet      !Wet surface heat exchanger efficiency
REAL NTUoWet     !Outlet wet surface NTU

REAL KFinFrost
REAL FinFrostThk
REAL KFrost
REAL FrostThk

!FLOW:

	AirPropOpt=2	
	AirProp%APTDB=tRiMod   !RS: Debugging: Formerly AirProp(1)
	AirProp%APRelHum=1    !RS: Debugging: Formerly AirProp(3)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	hRiSat=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)

	AirPropOpt=2	
	AirProp%APTDB=tRiMod+DeltaTemp !RS: Debugging: Formerly AirProp(1)
	AirProp%APRelHum=1    !RS: Debugging: Formerly AirProp(3)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	DeltahRiSat=ABS(hRisat-AirProp%APEnth)  !RS: Debugging: Formerly AirProp(4)

	CpSat=DeltahRiSat/DeltaTemp

	CpMoist=(1+0.09*wAiMod)*CpAir !www.jgsee.kmutt.ac.th/exell/JEE661/JEE661Lecture2.html

	IF(CoilType .EQ. EVAPORATORCOIL) THEN
	 Rcollar=ODtube/2+FinThk+FrostParam%Thickness
	ELSE
	 Rcollar=ODtube/2+FinThk
    END IF

	XM=Pt/2
	IF (XM .GT. Pl .AND. Pl .NE. 0) THEN
        XM = Pl
    END IF
	XL=0.5*((Pt/2)**2+Pl**2)**0.5
	Req=1.27*XM*(XL/XM-0.3)**(0.5)

	!Phi for fin efficiency calc.
	phi=(Req/Rcollar-1.0)*(1.0+0.35*LOG(Req/Rcollar))

	IF(CoilType .EQ. EVAPORATORCOIL) THEN
		Kfrost=FrostParam%Conductivity
		FrostThk=2*FrostParam%Thickness      
		KfinFrost=(Kfin*FinThk+Kfrost*FrostThk)/(FinThk+FrostThk)
		FinFrostThk=FinThk+FrostThk
		mm=SQRT(2*hcoMod/(KfinFrost*FinFrostThk)*CpSat/CpMoist)
	ELSE
		mm=SQRT(2*hcoMod/(Kfin*FinThk)*CpSat/CpMoist)
	END IF 
	
	FinEffwet=TANH(mm*Rcollar*phi)/(mm*Rcollar*phi)
	SurfEffwet=1-AfMod/AoMod*(1-FinEffwet)
	
	IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
		UAwet=1/(CpSat/(hciMod*AiMod*LmodTPratio)+CpSat*TubeThk/(Ktube*AmMod*LmodTPratio)+CpMoist/(SurfEffwet*hcoMod*AoMod*LmodTPratio))
	ELSE
		UAwet=1/(CpSat/(hciMod*AiMod)+CpSat*TubeThk/(Ktube*AmMod)+CpMoist/(SurfEffwet*hcoMod*AoMod))
	END IF

	cAir=mAiMod*CpSat
	
	NTUwet=UAwet/mAiMod
	EPSwet=1.-EXP(-NTUwet)

	DH=(hRiSat-hAiMod)
	QmodWet=EPSwet*mAiMod*DH
	
	hAoWet=QmodWet/mAiMod+hAiMod

	IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
		NTUoWet=SurfEff*hcoMod*AoMod*LmodTPratio/(mAiMod*CpMoist)
	ELSE
		NTUoWet=SurfEff*hcoMod*AoMod/(mAiMod*CpMoist)
	END IF

	hAoSatSurf=hAiMod-(hAiMod-hAoWet)/(1-EXP(-NTUoWet))

	tAoSat=TS(hAoSatSurf)   !RS: Replace: TS (2/19/14)

	TdbAoWet=tAoSat+(tAiMod-tAoSat)*EXP(-NTUoWet)

	DryWet=1
	Qmod=QmodWet
	cAir=mAiMod*CPair
	QmodSens=cAir*(TdbAoWet-tAiMod)
	SHR=QmodSens/Qmod
	hAoMod=hAoWet
	tAoMod=TdbAoWet

	IF (CoilType .NE. MCEVAPORATOR) THEN
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%hco =hcoMod
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%cAir=cAir
		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod*LmodTPratio)
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rtube=LOG(ODtube/IDtube)/(2*PI*ktube*LmodTube*LmodTPratio)
			IF(FrostParam%Conductivity==0.0) THEN
			 CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost = 0.0
			ELSE 
		     CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost=LOG((ODtube+FrostParam%Thickness)/ODtube)/(2*PI*FrostParam%Conductivity*LmodTube*LmodTPratio)
		    END IF 
		ELSE
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod)
			CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rtube=LOG(ODtube/IDtube)/(2*PI*ktube*LmodTube)
			IF(FrostParam%Conductivity==0.0) THEN
			 CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost = 0.0
			ELSE 
		     CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%Rfrost=LOG((ODtube+FrostParam%Thickness)/ODtube)/(2*PI*FrostParam%Conductivity*LmodTube)
		    END IF 
		END IF
		CoilSection(NumSection)%Ckt(II)%Tube(III)%Seg(IV)%SHR=SHR
	ELSE
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco = hcoMod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir=cAir
		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod*LmodTPratio)
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=LOG((ODtube/IDtube)/(2*PI*ktube*LmodTube*LmodTPratio))
		ELSE
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod)
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=LOG((ODtube/IDtube)/(2*PI*ktube*LmodTube))
		END IF
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%SHR=SHR
	END IF

RETURN

END SUBROUTINE CalcWetSurfaceBraun

!************************************************************************

SUBROUTINE CalcWetSurfaceMcQuiston(I,II,III,IV,CoilType)

!------------------------------------------------------------------------
!Purpose:
!To calculate wet surface air side heat transfer
!
!Reference:
!McQuiston, F. C., Parker, J. D. and Spitler, J.D. (2000). "Heating, 
!ventilating, and air conditioning analysis and design, 4th Ed." 
!New York, NY: John Wiley & Sons, Inc.
!
!Harms, T.M.; Groll, E.A. and Braun, J.E. (2003). Accurate charge inventory  
!modeling for unitary air conditioners. HVAC&R Research, 9(1), pp. 55-78.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!July 2006
!
!Reference:
!none
!
!------------------------------------------------------------------------

USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
USE CoilCalcMod
USE AirPropMod

IMPLICIT NONE

INTEGER,INTENT(IN) :: I   !Slab number
INTEGER,INTENT(IN) :: II  !Circuit,pass number
INTEGER,INTENT(IN) :: III !Tube number
INTEGER,INTENT(IN) :: IV  !Segment number

INTEGER,INTENT(IN) :: CoilType   !1=Condenser; 2=Evaporator; 
                                 !3=High side interconnecting pipes; 
								 !4=Low side interconnecting pipes
								 !5=Microchannel condenser
								 !6=Microchannel evaporator

REAL, PARAMETER :: DeltaTemp = 0.01 !Small temperature change, C

REAL hRiSat      !Saturation air enthalpy based on refrigerant temperature, kJ/kg
REAL DeltahRiSat !Small change of hRiSat, kJ/kg
REAL CpSat       !Saturation specific heat, kJ/kg
REAL CpMoist     !Specific heat of moist air, kJ/kg
REAL DH          !Enthalpy difference, kJ/kg
REAL NTUwet      !Wet surface NTU
REAL hAoSatSurf  !Outlet air surface saturation enthalpy, kJ/kg
REAL tAoSat      !Outlet air saturation temperature, C
REAL FinEffwet   !Wet surface fin efficiency
REAL SurfEffwet  !Wet surface efficiency
REAL EPSwet      !Wet surface heat exchanger efficiency
REAL NTUoWet     !Outlet wet surface NTU
REAL hcoWet      !Wet surface heat transfer coefficient, kW/m^2-K
REAL RairWet     !Wet surface resistance, kW/K
REAL wRiSat      !Saturated humidity ratio at refrigerant temperature
REAL hgRiSat     !Saturated vapor enthalpy at refrigerant temperature, kJ/kg
REAL hfRiSat     !Saturated liquid enthalpy at refrigerant temperature, kJ/kg
REAL CC          !Ratio of humidity ratio to Temperature difference, 1/C
REAL CCprev      !Previous value of CC
REAL hfgSat      !Enthalpy of vaporization at refrigerant temperature, kJ/kg
INTEGER iterCC   !Iteration counter for CC
REAL,PARAMETER :: Le=0.89 !1  !Lewis number

!FLOW:

	AirPropOpt=2	
	AirProp%APTDB=tRiMod   !RS: Debugging: Formerly AirProp(1)
	AirProp%APRelHum=1    !RS: Debugging: Formerly AirProp(3)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	hRiSat=AirProp%APEnth   !RS: Debugging: Formerly AirProp(4)
	wRiSat=AirProp%APHumRat   !RS: Debugging: Formerly AirProp(2)
	hfRiSat=hRiSat

	AirPropOpt=2	
	AirProp%APTDB=tRiMod   !RS: Debugging: Formerly AirProp(1)
	AirProp%APRelHum=0    !RS: Debugging: Formerly AirProp(3)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	hgRiSat=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)

	AirPropOpt=2	
	AirProp%APTDB=tRiMod+DeltaTemp
	AirProp%APRelHum=1    !RS: Debugging: Formerly AirProp(3)
	CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	DeltahRiSat=ABS(hRisat-AirProp%APEnth)  !RS: Debugging: Formerly AirProp(4)

	CpMoist=(1+0.09*wAiMod)*CpAir !www.jgsee.kmutt.ac.th/exell/JEE661/JEE661Lecture2.html

	!Initial guess
	wAoMod=wAiMod
	tAoMod=tAiMod
	hAoMod=hAiMod
	TwbAoMod=TwbAiMod
	CCprev=1e6

	DO iterCC=1, WetSurfaceMaxIter
	
		CC=((wAiMod-wRiSat)/(tAiMod-tRiMod)+(wAoMod-wRiSat)/(tAoMod-tRiMod))/2
		hfgSat=ABS(hgRiSat-hfRiSat)
		hcoWet=hcoMod*(1+hfgSat*CC/(Le*CpMoist)) !McQuiston
		
		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN  !Calculates the UA
			CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
						hcoWet,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio,UA,RairWet,Rrefrig,Rtube,FinEffwet,SurfEffWet)
		ELSE
			CALL CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth, &
						hcoWet,hciMod,AfMod,AoMod,AiMod,AmMod,UA,RairWet,Rrefrig,Rtube,FinEffwet,SurfEffWet)
            !CALL CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,ODtube,TubeThk,TubeDepth,RowNum,tAiMod,hAiMod, &
            !hcoMod,hciMod,AfMod*LmodTPratio,AoMod*LmodTPratio,AiMod*LmodTPratio,AmMod*LmodTPratio, &
            !UA,Rair,Rrefrig,Rtube,FinEff,SurfEff)
		END IF

		CpSat=(Le*CpMoist+hfgSat*CC)/Le !Saturation Specific Heat

		cAir=mAiMod*CpSat
		
		NTUwet=UA/cAir

		EPSwet=1.-EXP(-NTUwet)

		DH=(hRiSat-hAiMod)
		QmodWet=EPSwet*mAiMod*DH
		
		hAoWet=QmodWet/mAiMod+hAiMod

		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
			NTUoWet=SurfEff*hcoMod*AoMod*LmodTPratio/(mAiMod*CpMoist)
		ELSE
			NTUoWet=SurfEff*hcoMod*AoMod/(mAiMod*CpMoist)	
		END IF

		hAoSatSurf=hAiMod-(hAiMod-hAoWet)/(1-EXP(-NTUoWet))

		tAoSat=TS(hAoSatSurf)   !RS: Replace: TS (2/19/14)

		TdbAoWet=tAoSat+(tAiMod-tAoSat)*EXP(-NTUoWet)

		DryWet=1
		Qmod=QmodWet
		cAir=mAiMod*CPair
		QmodSens=cAir*(TdbAoWet-tAiMod)
		SHR=QmodSens/Qmod
		hAoMod=hAoWet
		tAoMod=TdbAoWet

		AirPropOpt=1	
		AirProp%APTDB=tAoMod   !RS: Debugging: Formerly AirProp(1)
		AirProp%APEnth=hAoMod   !RS: Debugging: Formerly AirProp(4)
		CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
		wAoMod=AirProp%APHumRat   !RS: Debugging: Formerly AirProp(2)
		TwbAoMod=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)

		IF (ABS((CCprev-CC)/CCprev) .LT. 1e-3) THEN
			EXIT 
		ELSE
			CCprev=CC
		END IF

	END DO 

	IF (CoilType .NE. MCEVAPORATOR) THEN
		Ckt(II)%Tube(III)%Seg(IV)%hco = hcoMod
		Ckt(II)%Tube(III)%Seg(IV)%cAir=cAir
		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
			Ckt(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod*LmodTPratio)
			Ckt(II)%Tube(III)%Seg(IV)%Rtube=LOG(ODtube/IDtube)/(2*PI*ktube*LmodTube*LmodTPratio)
		    Ckt(II)%Tube(III)%Seg(IV)%Rfrost=LOG((ODtube+FrostParam%Thickness)/ODtube)/(2*PI*FrostParam%Conductivity*LmodTube*LmodTPratio)
		ELSE
			Ckt(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod)
			Ckt(II)%Tube(III)%Seg(IV)%Rtube=LOG(ODtube/IDtube)/(2*PI*ktube*LmodTube)
		    Ckt(II)%Tube(III)%Seg(IV)%Rfrost=LOG((ODtube+FrostParam%Thickness)/ODtube)/(2*PI*FrostParam%Conductivity*LmodTube)
		END IF
		Ckt(II)%Tube(III)%Seg(IV)%SHR=SHR
	ELSE
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%hco = hcoMod
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%cAir=cAir
		IF (LmodTPratio .NE. 0 .AND. IsSimpleCoil .EQ. 1) THEN
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod*LmodTPratio)
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=LOG((ODtube/IDtube)/(2*PI*ktube*LmodTube*LmodTPratio))
		ELSE
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rair=1/(SurfEffwet*hcoMod*AoMod)
			Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Rtube=LOG((ODtube/IDtube)/(2*PI*ktube*LmodTube))
		END IF
		Slab(I)%Pass(II)%Tube(III)%Seg(IV)%SHR=SHR
	END IF

RETURN

END SUBROUTINE CalcWetSurfaceMcQuiston

!************************************************************************

SUBROUTINE MicrochannelEvaporator(XIN,PAR,OUT) !(Ref$,XIN,PAR,OUT)  !RS: Debugging: Extraneous Ref$

    !-----------------------------------------------------------------------------------
    !
    !  Description:	
    !  A segment-by-segment microchannel evaporator model
    !  To predict coil air side and refrigerant side properties, heat transfer, 
    !  and pressure drop
    !
    !  Inputs:
    !  Ref$=Refrigerant name
    !  XIN(1)=Refrigerant side mass flow rate, kg/s
    !  XIN(2)=Refrigerant side inlet pressure, kPa
    !  XIN(3)=Refrigerant side inlet enthalpy, kJ/kg
    !  XIN(4)=Air side mass flow rate, kg/s
    !  XIN(5)=Air side inlet temp. C
    !  XIN(6)=Air side inlet relative humidity
    !  XIN(7)=Compressor discharge temperature, C
    !
    !  Parameters:
    !  PAR(1)=Barometric pressure, kPa
    !  PAR(2)=Cooling mode? 1=yes; 0=no  
    !  PAR(3)=Suction line length, m
    !  PAR(4)=Suction line outside diameter, m
    !  PAR(5)=Suction line tube wall thickness, m 
    !  PAR(6)=Suction line elevation, m
    !  PAR(7)=Suction line heat loss, kW
    !  PAR(8)=Suction line temperature change, C
    !  PAR(9)=Suction line additional pressure drop, kPa
    !  PAR(10)=Multiplier for ref. side heat transfer correlation
    !  PAR(11)=Multiplier for ref. side pressure drop correlation
    !  PAR(12)=Multiplier for air side heat transfer correlation
    !  PAR(13)=Multiplier for air side pressure drop correlation
    !  PAR(14)=Fan power, kW
    !  PAR(15)=Fan location, 1=draw through; 2=blow through
    !  PAR(16)=Compressor heat loss, kW
    !  PAR(17)=Is compressor in air stream, 1=yes, 0=no
    !  PAR(18)=System Type, 1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat
    !  PAR(19)=Custom air side data unit, 1=SI; 2=IP
    !  PAR(20)=Custom air heat transfer curve type, 1=Power; 2=Polynomial
    !  PAR(21)=Power coefficient for air heat transfer curve
    !  PAR(22)=Power coefficient for air heat transfer curve
    !  PAR(23)=Polynomial coefficient for air heat transfer curve
    !  PAR(24)=Polynomial coefficient for air heat transfer curve
    !  PAR(25)=Polynomial coefficient for air heat transfer curve
    !  PAR(26)=Polynomial coefficient for air heat transfer curve
    !  PAR(27)=Custom air heat transfer curve type, 1=Power; 2=Polynomial
    !  PAR(28)=Power coefficient for air heat transfer curve
    !  PAR(29)=Power coefficient for air heat transfer curve
    !  PAR(30)=Polynomial coefficient for air heat transfer curve
    !  PAR(31)=Polynomial coefficient for air heat transfer curve
    !  PAR(32)=Polynomial coefficient for air heat transfer curve
    !  PAR(33)=Polynomial coefficient for air heat transfer curve
    !
    !  Outputs:
    !  OUT(1)=Coil capacity, kW
    !  OUT(2)=Sensible coil capacity, kW
    !  OUT(3)=Coil outlet pressure, kPa
    !  OUT(4)=Coil outlet enthalpy, kJ/kg
    !  OUT(5)=Coil outlet temperature, C
    !  OUT(6)=Coil outlet quality
    !  OUT(7)=Coil outlet superheat, C
    !  OUT(8)=Suction line outlet pressure, kPa
    !  OUT(9)=Suction line outlet enthalpy, kJ/kg
    !  OUT(10)=Suction line outlet temperature, C
    !  OUT(11)=Suction line outlet quality
    !  OUT(12)=Suction line outlet superheat, C
    !  OUT(13)=Air side outlet temperature, C
    !  OUT(14)=Air side outlet relative humidity
    !  OUT(15)=Air side pressure drop, kPa
    !  OUT(16)=Mass in suction line, kg
    !  OUT(17)=Mass in coil, kg
    !  OUT(18)=Liquid mass in coil, kg
    !  OUT(19)=Vapor mass in coil, kg
    !  OUT(20)=Aluminum weight, kg 
    !  OUT(21)=Error flag: 0-No error
    !                      1-Evaporator solution not converge
    !                      2-Refprop error
    !					   3-Circuit file error
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

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE CoilCalcMod
    USE AirPropMod
    USE OilMixtureMod
    USE ReversingValveMod
    USE InputProcessor    !RS: Debugging: GetObjectItem

    IMPLICIT NONE

    !Subroutine argument declarations
    !CHARACTER*80, INTENT(IN)  :: Ref$  !RS: Debugging: Extraneous
    REAL,         INTENT(IN)  :: XIN(7)
    REAL,         INTENT(IN)  :: PAR(18)    !RS: Debugging: Formerly PAR(33)
    REAL,         INTENT(OUT) :: OUT(21)

    !Subroutine local variables
    INTEGER,PARAMETER :: CoilType = MCEVAPORATOR 

    INTEGER I,II,III,IV,V !Loop counters
    LOGICAL Converged     !Convergence flag
    INTEGER RefBCiter     !Refrigerant bounadary condition iteration counter
    REAL Qpass			  !Coil pass capacity, kW
    REAL QpassSens		  !Coil pass sensible capacity, kW
    REAL QinletPass		  !Inlet pass capacity, kW
    REAL QinletPassSens   !Inlet pass sensible capacity, kW
    REAL pRoSlab		  !Outlet Refrigerant pressure for a coil slab, kPa
    REAL hRoSlab		  !Outlet Refrigerant enthalpy for a coil slab, kJ/kg
    REAL tRoSlab		  !Outlet Refrigerant temperature for a coil slab, C  
    REAL xRoSlab		  !Outlet Refrigerant quality for a coil slab  
    REAL Aface			  !Coil face area, m^2
    REAL SumPro			  !Sum of outlet pressures, kPa
    REAL SumMrefHro		  !Sum of mdot*H (mass flow rate * enthalpy)
    REAL tRdis            !Compressor discharge temperature, C
    REAL DPvalve          !Reversing valve pressure drop
    REAL hRsuc            !Suction enthalpy, kJ/kg
    REAL DPcoil, DPcoilPrev !Coil pressure drop, kPa
    REAL mdothRo !mdot x outlet enthalpy
    REAL humrat !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
    
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
    pRiCoil          =XIN(2)
    hRiCoil          =XIN(3)
    mAiCoil          =XIN(4)
    tAiCoil          =XIN(5)
    rhAiCoil         =XIN(6)
    tRdis            =XIN(7)

    BaroPressure     =PAR(1)
    !IsCoolingMode    =PAR(2)   !RS: Debugging: Global variable now
    LsucLn           =PAR(3)
    !ODsucLn          =PAR(4)    !RS: Debugging: Never used?
    SucLnThk         =PAR(5)
    !ElevSucLn        =PAR(6)   !RS: Debugging: Never used?
    QsucLn           =PAR(7)
    !DTsucLn          =PAR(8)   !RS: Debugging: Never used?
    AddDPsucLn       =PAR(9)
    hciMultiplier    =PAR(10)
    DPrefMultiplier  =PAR(11)
    hcoMultiplier    =PAR(12)
    DPairMultiplier  =PAR(13)
    PwrFan           =PAR(14)
    DrawBlow         =PAR(15)
    QlossCmp         =PAR(16)
    !IsCmpInAirStream =PAR(17)
    SystemType       =PAR(18)
    !CurveUnit        =PAR(19)  !RS: Debugging: Never Really Used
    !CurveTypeHTC     =PAR(20)  !RS: Debugging: Never Really Used
    !PowerAHTC        =PAR(21)  !RS: Debugging: Never Really Used
    !PowerBHTC        =PAR(22)  !RS: Debugging: Never Really Used
    !Poly1HTC         =PAR(23)  !RS: Debugging: Never Really Used
    !Poly2HTC         =PAR(24)  !RS: Debugging: Never Really Used
    !Poly3HTC         =PAR(25)  !RS: Debugging: Never Really Used
    !Poly4HTC         =PAR(26)  !RS: Debugging: Never Really Used
    !CurveTypeDP      =PAR(27)  !RS: Debugging: Never Really Used
    !PowerADP         =PAR(28)  !RS: Debugging: Never Really Used
    !PowerBDP         =PAR(29)  !RS: Debugging: Never Really Used
    !Poly1DP          =PAR(30)  !RS: Debugging: Never Really Used
    !Poly2DP          =PAR(31)  !RS: Debugging: Never Really Used
    !Poly3DP          =PAR(32)  !RS: Debugging: Never Really Used
    !Poly4DP          =PAR(33)  !RS: Debugging: Never Really Used
    
    !********************Refrigerant Cycle Data (Heating)***********************  !RS: Debugging: Moving: Condenser

  CALL GetObjectItem('RefrigerantCycleData(Heating)',1,Alphas,NumAlphas, &
                      TmpNumbers,NumNumbers,Status)
  Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

  IsCmpInAirStream = Numbers(2) !Is Compressor in Air Stream

    IsParallelSlabs = 1

    ErrorFlag=NOERROR !Initilaize

    !Coil height
    HtCoil=Nt*Pt

    !Coil length
    Lcoil=Nl*Nt*Ltube

    !Face area
    Aface=Ltube*Nt*Pt

    !Tube information
    LmodTube=Ltube/NumOfMods

    !RS: Replace: Moving this up so that Cp will be calculated with updated air properties
    AirPropOpt=2
    AirProp%APTDB=tAiCoil  !RS: Debugging: Formerly AirProp(1)
    AirProp%APRelHum=rhAiCoil !RS: Debugging: Formerly AirProp(3)
    CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
    hAiCoil=AirProp%APEnth  !RS: Debugging: Formerly AirProp(4)
    wbAiCoil=AirProp%APTWB !RS: Debugging: Formerly AirProp(5)
    
    CPair=CPA(REAL(tAiCoil))    !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAiCoil,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAiCoil*CPAir

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

    Pressure=pRiCoil*1000   !RS Comment: Unit Conversion
    Enthalpy=hRiCoil*1000   !RS Comment: Unit Conversion
    tRiCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !Coil Refrigerant Inlet Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error. Line 3904'
        ErrorFlag=REFPROPERROR
        RETURN
    END IF
    xRiCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !Coil Refrigerant Inlet Quality
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error. Line 3910'
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
    QmodPrev=1e20
    DPcoilPrev=1e20
    Qcoil=0.0; QcoilSens=0;

    DO RefBCiter=1, MdotMaxIter

        DO I=1, Nl !Number of slabs

            DO II=1,Slab(I)%Npass !Number of passes

                Qpass=0.0; QpassSens=0

                CALL CalcCircuitRefInletConditions(I,I,II,CoilType) !ISI - 09/10/07

                IF (II .EQ. 1 .AND. Slab(I)%Ninlet .GT. 1) THEN !Multi-inlet

                    SumPro=0; SumMrefHro=0;

                    DO V=1, Slab(I)%Ninlet

                        QinletPass=0; QinletPassSens=0;

                        Slab(I)%InletPass(V)%pRi=pRiCoil
                        Slab(I)%InletPass(V)%hRi=hRiCoil

                        DO III=1,1 !Slab(I)%Pass(II)%NumOfTubes !Number of tubes
                            DO IV=1,NumOfMods !Number of segments

                                !ref. mass flow rate for each channel
                                mRefMod=Slab(I)%mdot/Slab(I)%Pass(II)%Ntube/NumOfChannels !ISI - 07/13/07

                                !Module Surface Areas
                                AoMod=AoCoil*LmodTube/Lcoil 
                                AfMod=AfCoil*LmodTube/Lcoil 
                                AiMod=AiCoil*LmodTube/Lcoil 
                                AmMod=AmCoil*LmodTube/Lcoil 

                                CALL CalcCoilSegment(I,I,II,III,IV,CoilType)
                                IF (ErrorFlag .GT. CONVERGEERROR) THEN
                                    OUT(21)=ErrorFlag
                                    RETURN

                                END IF

                                QinletPass=QinletPass+Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod
                                QinletPassSens=QinletPassSens+Slab(I)%Pass(II)%Tube(III)%Seg(IV)%QmodSens

                            END DO !End mod
                        END DO !End tube

                        Slab(I)%InletPass(V)%Qpass=QinletPass*Slab(I)%InletPass(V)%Ntube
                        Slab(I)%InletPass(V)%QpassSens=QinletPassSens*Slab(I)%InletPass(V)%Ntube

                        Slab(I)%InletPass(V)%pRo=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%pRo
                        Slab(I)%InletPass(V)%hRo=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%hRo
                        Slab(I)%InletPass(V)%mRef=mRefTot*Slab(I)%InletPass(V)%Ntube/Slab(I)%Pass(II)%Ntube

                        Qpass=Qpass+Slab(I)%InletPass(V)%Qpass
                        QpassSens=QpassSens+Slab(I)%InletPass(V)%QpassSens

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

                            CALL CalcCoilSegment(I,I,II,III,IV,CoilType)
                            IF (ErrorFlag .GT. CONVERGEERROR) THEN
                                OUT(21)=ErrorFlag
                                RETURN
                            END IF

                            Qpass=Qpass+Slab(I)%Pass(II)%Tube(III)%Seg(IV)%Qmod
                            QpassSens=QpassSens+Slab(I)%Pass(II)%Tube(III)%Seg(IV)%QmodSens

                        END DO !End mod

                    END DO !End tube

                    Qpass=Qpass*Slab(I)%Pass(II)%Ntube
                    QpassSens=QpassSens*Slab(I)%Pass(II)%Ntube

                    pRoCkt=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%pRo !Circuit outlet pressure
                    hRoCkt=Slab(I)%Pass(II)%Tube(1)%Seg(NumOfMods)%hRo !Circuit outlet enthalpy

                END IF

                Slab(I)%Pass(II)%pRo=pRoCkt
                Slab(I)%Pass(II)%hRo=hRoCkt

                Qcoil=Qcoil+Qpass
                QcoilSens=QcoilSens+QpassSens

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

            SHR=QcoilSens/Qcoil

            IF (Slab(1)%Npass .EQ. 1 .AND. Slab(1)%Ninlet .GT. 1) THEN !1-pass, multi-inlet
                pRoSlab=pRoCkt
                hRoSlab=hRoCkt
            ELSE !Single inlet
                pRoSlab=Slab(I)%Pass(Slab(I)%Npass)%pRo
                hRoSlab=Slab(I)%Pass(Slab(I)%Npass)%hRo
            END IF

            Pressure=pRoSlab*1000   !RS Comment: Unit Conversion
            Enthalpy=hRoSlab*1000   !RS Comment: Unit Conversion
            tRoSlab=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !Slab Refrigerant Outlet Temperature
            IF (RefPropErr .GT. 0) THEN
                WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
                ErrorFlag=REFPROPERROR
                OUT(21)=ErrorFlag
                RETURN
            END IF
            xRoSlab=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !Slab Refrigerant Outlet Quality
            IF (RefPropErr .GT. 0) THEN
                WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
                ErrorFlag=REFPROPERROR
                OUT(21)=ErrorFlag
                RETURN
            END IF

            Slab(I)%pRo=pRoSlab
            Slab(I)%hRo=hRoSlab
            Slab(I)%tRo=tRoSlab
            Slab(I)%xRo=xRoSlab

            Slab(I)%Qslab=mRefTot*(Slab(I)%hRi-Slab(I)%hRo)

            CALL CalcMeanProp(Slab(I)%tAi,Slab(I)%tAo,tAmod)    !Calculating mean module air temperature

            !Coil air side outlet conditions
            CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
            CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
            Cair=mAicoil*CPAir

            Slab(I)%tAo=Slab(I)%tAi+Slab(I)%Qslab/Cair
            Slab(I)%hAo=Slab(I)%hAi+Slab(I)%Qslab/mAiCoil

            AirPropOpt=1
            AirProp%APTDP=Slab(I)%tAo  !RS: Debugging: Formerly AirProp(1)
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

    IF (ABS(QcoilSens) .GT. ABS(Qcoil)) THEN !Make sure sensible heat no larger than total heat, ISI - 12/07/07
        QcoilSens=Qcoil
    END IF

    !Coil air side outlet conditions
    !CPair=CPA(REAL(tAmod))  !RS: Replace: CPA (2/19/14)
    CPair=CPAirFunction(tAmod,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
    Cair=mAicoil*CPAir

    tAoCoil=tAiCoil+QcoilSens/Cair
    hAoCoil=hAiCoil+Qcoil/mAiCoil

    !Fan air side inlet conditions
    
    !CPair=CPA(REAL(tAoCoil))    !RS: Replace: CPA (2/19/14)
    humrat=HUMTH(tAoCoil,hAoCoil)  !RS: Replace: CPA (2/20/14) Finding outlet humidity ratio
    CPair=CPAirFunction(tAoCoil,humrat)  !RS: Replace: CPA (2/19/14)
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

    Pressure=pRoCoil*1000
    Enthalpy=hRoCoil*1000
    tRoCoil=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)  !Coil Refrigerant Outlet Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN

    END IF
    xRoCoil=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)  !Coil Refrigerant Outlet Quality
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN

    END IF

    Pressure=pRoCoil*1000   !RS Comment: Unit Conversion
    Quality=0
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN

    END IF

    !IF (Wabsolute .GT. 0) THEN
    !  Wlocal=LocalOilMassFraction(Wabsolute,xRoCoil)
    !  tSat=OilMixtureTsat(RefName,Wlocal,Psat/1000)
    !END IF

    IF (xRoCoil .LE. 0.0) THEN 
        tSHoCoil=tRoCoil-tSat !Superheat
    ELSE
        tSHoCoil=0.0
    END IF

    !****** Suction line calculation ******
    IF (LsucLn .GT. 0) THEN 
        CALL SuctionLine
        IF (ErrorFlag .GT. CONVERGEERROR) THEN
            WRITE(*,*)'SuctionLine: Refprop error.'
            OUT(21)=ErrorFlag
            RETURN

        END IF
    ELSE
        pRiCmp=pRoCoil
        hRiCmp=hRoCoil
    END IF

    IF (SystemType .EQ. HEATPUMP) THEN !Heat pump
        !Calculate reversing valve heat transfer and pressure drop
        Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
        Enthalpy=hRiCmp*1000    !RS Comment: Unit Conversion
        tRiCmp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Compressor Inlet Refrigerant Temperature
        IF (RefPropErr .GT. 0) THEN
            WRITE(*,*)'-- WARNING -- Evaporator: Refprop error. Line 2646'
            ErrorFlag=REFPROPERROR
            OUT(21)=ErrorFlag
            RETURN
        END IF

        hRsuc=hRiCmp
        IF (IsCoolingMode .GT. 0) THEN
            !Correlation good for cooling only
            CALL SuctionHeatTransfer(mRefTot,tRdis,tRiCmp,hRsuc,hRiCmp)
            CALL SuctionPressureDrop(mRefTot,DPvalve)
            pRiCmp=pRiCmp-DPvalve
        END IF

    END IF

    Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
    Enthalpy=hRiCmp*1000    !RS Comment: Unit Conversion
    tRiCmp=PH(RefName, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr)   !Compressor Inlet Refrigerant Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN
    END IF
    xRiCmp=PH(RefName, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Compressor Inlet Refrigerant Quality
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN
    END IF

    Pressure=pRiCmp*1000    !RS Comment: Unit Conversion
    Quality=0
    tSat=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)  !Saturation Temperature
    IF (RefPropErr .GT. 0) THEN
        WRITE(*,*)'-- WARNING -- MCEvaporator: Refprop error.'
        ErrorFlag=REFPROPERROR
        OUT(21)=ErrorFlag
        RETURN
    END IF

    !IF (Wabsolute .GT. 0) THEN
    !  Wlocal=LocalOilMassFraction(Wabsolute,xRiCmp)
    !  tSat=OilMixtureTsat(RefName,Wlocal,Psat/1000)
    !END IF

    IF (xRiCmp .GE. 1) THEN
        tSHiCmp=tRiCmp-tSat
    ELSE
        tSHiCmp=0
    END IF

    OUT(1)=Qcoil
    OUT(2)=QcoilSens
    OUT(3)=pRoCoil
    OUT(4)=hRoCoil
    OUT(5)=tRoCoil
    OUT(6)=xRoCoil
    OUT(7)=tSHoCoil
    OUT(8)=pRiCmp
    OUT(9)=hRiCmp
    OUT(10)=tRiCmp
    OUT(11)=xRiCmp
    OUT(12)=tSHiCmp
    OUT(13)=tAoCoil
    OUT(14)=rhAoCoil
    OUT(15)=DPair
    OUT(16)=MassSucLn
    OUT(17)=0
    OUT(18)=0
    OUT(19)=0
    OUT(20)=WeightAluminum

    OUT(21)=ErrorFlag

    RETURN

    END SUBROUTINE MicrochannelEvaporator

!************************************************************************

SUBROUTINE LoadMicrochannelInputs(MCXIN,MCPAR) !(FTXIN,FTPAR,MCXIN,MCPAR)

!-----------------------------------------------------------------------------------
!
!  Description:	
!  Transfer input data from subroutine "Evaporator" to "MircochannelEvaporator"
!
!  Author:
!  Ipseng Iu
!  Mechanical and Aerospace Engineering
!  Oklahoma State University, Stillwater	
!
!  Date: November 2005
!
!-----------------------------------------------------------------------------------

USE DataSimulation

IMPLICIT NONE

!REAL, INTENT(IN)  :: FTXIN(9)  !Fin-tube coil input data
!REAL, INTENT(IN)  :: FTPAR(39) !Fin-tube coil input parameters  !RS: Debugging: Formerly FTPAR(49)
REAL, INTENT(OUT) :: MCXIN(7)  !Microchannel coil input data
REAL, INTENT(OUT) :: MCPAR(19) !Microchannel coil input parameters  !RS: Debugging: Formerly MCPAR(33)

!FLOW:

  MCXIN(1)=EvapIN%EInmRef !Refrigerant side mass flow rate, kg/s  !RS: Debugging: Formerly FTXIN(1)
  MCXIN(2)=EvapIN%EInpRi !Refrigerant side inlet pressure, kPa   !RS: Debugging: Formerly FTXIN(2)
  MCXIN(3)=EvapIN%EInhRi !Refrigerant side inlet enthalpy, kJ/kg !RS: Debugging: Formerly FTXIN(3)
  MCXIN(4)=EvapIN%EInmAi !Air side mass flow rate, kg/s  !RS: Debugging: Formerly FTXIN(4)
  MCXIN(5)=EvapIN%EIntAi !Air side inlet temp. C !RS: Debugging: Formerly FTXIN(5)
  MCXIN(6)=EvapIN%EInrhAi !Air side inlet relative humidity   !RS: Debugging: Formerly FTXIN(6)
  MCXIN(7)=EvapIN%EIntRdis !Compressor discharge temperature, C    !RS: Debugging: Formerly FTXIN(7)

  MCPAR(1)=EvapPAR%EvapBarPress  !Barometric pressure, kPa !RS: Debugging: Formerly FTPAR(31)
  MCPAR(3)=EvapPAR%EvapSucLnLen   !Suction line length, m   !RS: Debugging: Formerly FTPAR(1)
  MCPAR(4)=EvapPAR%EvapSucLnOD   !Suction line outside diameter, m !RS: Debugging: Formerly FTPAR(2)
  MCPAR(5)=EvapPAR%EvapSucLnTWThick   !Suction line tube wall thickness, m  !RS: Debugging: Formerly FTPAR(3)
  MCPAR(6)=EvapPAR%EvapSucLnElev   !Suction line elevation, m    !RS: Debugging: Formerly FTPAR(4)
  MCPAR(7)=EvapPAR%EvapSucLnQLoss   !Suction line heat loss, kW   !RS: Debugging: Formerly FTPAR(5)
  MCPAR(8)=EvapPAR%EvapSucLnTempChg   !Suction line temperature change, C   !RS: Debugging: Formerly FTPAR(6)
  MCPAR(9)=EvapPAR%EvapSucLnAddPD   !Suction line additional pressure drop, kPa   !RS: Debugging: Formerly FTPAR(7)
  MCPAR(10)=EvapPAR%EvapMultRefQT !Multiplier for ref. side heat transfer correlation   !RS: Debugging: Formerly FTPAR(23)
  MCPAR(11)=EvapPAR%EvapMultRefPD !Multiplier for ref. side pressure drop correlation   !RS: Debugging: Formerly FTPAR(24)
  MCPAR(12)=EvapPAR%EvapMultAirQT !Multiplier for air side heat transfer correlation    !RS: Debugging: Formerly FTPAR(25)
  MCPAR(13)=EvapPAR%EvapMultAirPD !Multiplier for air side pressure drop correlation    !RS: Debugging: Formerly FTPAR(26)
  MCPAR(14)=EvapPAR%EvapFanPwr !Fan power, kW    !RS: Debugging: Formerly FTPAR(27)
  MCPAR(15)=EvapPAR%EvapFanLoc !Fan location, 1=draw through; 2=blow through !RS: Debugging: Formerly FTPAR(28)
  MCPAR(16)=EvapPAR%EvapCompQLoss !Compressor heat loss, kW !RS: Debugging: Formerly FTPAR(32)
  MCPAR(17)=EvapPAR%EvapSysType !Is compressor in air stream, 1=yes, 0=no !RS: Debugging: Formerly FTPAR(33)
  MCPAR(18)=EvapPAR%EvapPressTolConv !System Type, 1=A/C, 2=Heat Pump, 3=Condenser Unit, 4=Reheat  !RS: Debugging: Formerly FTPAR(34)

  RETURN

END SUBROUTINE LoadMicrochannelInputs

!************************************************************************

SUBROUTINE LoadMicrochannelOutputs(MCOUT) !,FTOUT)

!-----------------------------------------------------------------------------------
!
!  Description:	
!  Transfer output data from subroutine "MircochannelEvaporator" to "Evaporator"
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

REAL, INTENT(IN)  :: MCOUT(21)  !Microchannel coil output data
!REAL, INTENT(OUT) :: FTOUT(17)  !Fin-tube coil output data  !RS: Debugging: Formerly FTOUT(20)

!FLOW:

  EvapOUT%EOutpRoC=MCOUT(3)   !Coil outlet pressure, kPa    !RS: Debugging: Formerly FTOUT(1)
  EvapOUT%EOuthRoC=MCOUT(4)   !Coil outlet enthalpy, kJ/kg  !RS: Debugging: Formerly FTOUT(2)
  EvapOUT%EOutpRiC=MCOUT(8)   !Suction line outlet pressure, kPa    !RS: Debugging: Formerly FTOUT(6)
  EvapOUT%EOuthRiC=MCOUT(9)   !Suction line outlet enthalpy, kJ/kg  !RS: Debugging: Formerly FTOUT(7)
  EvapOUT%EOuttRiC=MCOUT(10)  !Suction line outlet temperature, C   !RS: Debugging: Formerly FTOUT(8)
  EvapOUT%EOutxRiC=MCOUT(11)  !Suction line outlet quality  !RS: Debugging: Formerly FTOUT(9)
  EvapOUT%EOuttSHiC=MCOUT(12) !Suction line outlet superheat, C !RS: Debugging: Formerly FTOUT(10)
  EvapOUT%EOutQC=MCOUT(1)  !Coil capacity, kW    !RS: Debugging: Formerly FTOUT(11)
  EvapOUT%EOutQCSens=MCOUT(2)  !Sensible coil capacity, kW   !RS: Debugging: Formerly FTOUT(12)
  EvapOUT%EOutMSucLn=MCOUT(16) !Mass in suction line, kg !RS: Debugging: Formerly FTOUT(13)
  EvapOUT%EOutMC=MCOUT(17) !Mass in coil, kg !RS: Debugging: Formerly FTOUT(14)
  EvapOUT%EOuttAoC=MCOUT(13) !Air side outlet temperature, C   !RS: Debugging: Formerly FTOUT(3)
  EvapOUT%EOutrhAoC=MCOUT(14) !Air side outlet relative humidity    !RS: Debugging: Formerly FTOUT(4)
  EvapOUT%EOutDPAir=MCOUT(15) !Air side pressure drop, kPa  !RS: Debugging: Formerly FTOUT(5)
  EvapOUT%EOutErrFlag=MCOUT(21) !Error flag   !RS: Debugging: Formerly FTOUT(17)
  EvapOUT%EOutWtAl=MCOUT(20) !Aluminum weight, kg  !RS: Debugging: Formerly FTOUT(15)
  EvapOUT%EOutWtCu=0         !Copper weight, kg    !RS: Debugging: Formerly FTOUT(16)

  RETURN

END SUBROUTINE LoadMicrochannelOutputs

!************************************************************************

SUBROUTINE UpdateTubeDataFromCircuitData(NumSection,I,J)

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
INTEGER, INTENT(IN) :: NumSection !Section number
INTEGER, INTENT(IN) :: I !Circuit number
INTEGER, INTENT(IN) :: J !Tube number in calculation sequence

!Subroutine local variables
INTEGER TubeNum !Tube number in circuit diagram

!FLOW:

  TubeNum=CoilSection(NumSection)%Ckt(I)%TubeSequence(J)
  Tube(TubeNum)=CoilSection(NumSection)%Ckt(I)%Tube(J)

  RETURN

END SUBROUTINE UpdateTubeDataFromCircuitData

SUBROUTINE GetQOut(Out1, Out2, Out3)

REAL Out1
REAL Out2
REAL Out3
!REAL WAi,HAi,TWBAi,TDPAi,RhoDAi,RhoMAi,BaroPressureAi,ErrStatAi
REAL WAi
REAL WAo,HAo,TWBAo,TDPAo,RhoDAo,RhoMAo,BaroPressureAo,ErrStatAo

    Out1=QModSensTot*1000 !Sensible Module heat transfer, W
    CALL TDB_RH (tAoMod,WAo,rhAoMod,HAo,TWBAo,TDPAo,RhoDAo,RhoMAo,BaroPressureAo,ErrStatAo)
    Out2= mAiMod*(WAo-WAi) !Latent Output, kg/s
    Out3=QModTot*1000
    !RS: Debugging: Something that should probably be looked at is if that's really the total latent output
    !RS: Debugging: Does some sort of average need to be taken?
    
END SUBROUTINE

SUBROUTINE GetNodeProp()    !RS: Debugging: Updating the air node data for E+

USE DataLoopNode
USE DXCoils !, ONLY: DXCoil, DXCoilHPSimNum
USE DataGlobals_HPSimIntegrated, ONLY: CurrentTime, BeginTimeStepFlag
USE DataHVACGlobals, ONLY: SysTimeElapsed

  INTEGER :: AirOutletNode ! air outlet node number
  INTEGER :: AirInletNode ! air inlet node number
  REAL tWBi,tDPi,RhoDi,RhoMi,ErrStati
  REAL tWBo,tDPo,RhoDo,RhoMo,ErrStato
  REAL hrAiCoil,hrAoCoil
  REAL Qair
  INTEGER :: LogFile       =153 !RS: Debugging file denotion, hopefully this works.
  CHARACTER(LEN=25),PARAMETER :: FMT_153 = "(1(A12,','),5(F10.3,','))"
  REAL, SAVE :: SysTimeTemp=0   !RS: Debugging: Keeping track of the system time (8/21/14)
  
    !OPEN(unit=LogFile,file='logfile.csv')    !RS: Debugging
    
    !AirInletNode=NodeConnectionType_OutsideAir  !RS: Debugging: Coming directly from outside
    !AirOutletNode=2 !NodeConnectionType_Internal !RS: Debugging: This is what it needs to be to connect the mass flow rate for RA-only case

    !CALL TDB_RH(tAiMod,hrAiMod,rhAiMod,hAiMod,tWBi,tDPi,RhoDi,RhoMi,BaroPressure,ErrStati)
    !CALL TDB_RH(tAoMod,hrAoMod,rhAoMod,hAoMod,tWBo,tDPo,RhoDo,RhoMo,BaroPressure,ErrStato)
    !RS: Debugging: Changed all the Mods to Coils because I think Coil might actually be correct
    CALL TDB_RH(tAiCoil,hrAiCoil,rhAiCoil,hAiCoil,tWBi,tDPi,RhoDi,RhoMi,BaroPressure,ErrStati)
    CALL TDB_RH(tAoCoil,hrAoCoil,rhAoCoil,hAoCoil,tWBo,tDPo,RhoDo,RhoMo,BaroPressure,ErrStato)
    
!Node(AirInletNode)%Enthalpy     = hAiMod
!Node(AirInletNode)%Temp         = tAiMod
!Node(AirInletNode)%HumRat       = hrAiMod
!Node(AirInletNode)%MassFlowRate = mAiMod
!Node(AirInletNode)%Press =BaroPressure
! changed outputs
!Node(AirOutletNode)%Enthalpy     = hAoMod
!Node(AirOutletNode)%Temp         = tAoMod
!Node(AirOutletNode)%HumRat       = hrAoMod
!Node(AirOutletNode)%MassFlowRate = mAiMod
! pass through outputs

  DXCoil(DXCoilHPSimNum)%InletAirMassFlowRate = mAiCoil !mAiMod !RS: Debugging: Testing to see the correct mdot air to pass back

  DXCoil(DXCoilHPSimNum)%OutletAirTemp     = tAoCoil !Mod
  DXCoil(DXCoilHPSimNum)%OutletAirHumRat   = hrAoCoil !Mod
  DXCoil(DXCoilHPSimNum)%OutletAirEnthalpy = hAoCoil !Mod
  Qair=mAiCoil*(hAoCoil-hAiCoil)
  
  DXCoil(DXCoilHPSimNum)%CndFanPwr=CondPAR%CondFanPwr   !RS: Debugging: Trying to pass through power from HPSim (9/3/14)
  DXCoil(DXCoilHPSimNum)%CompPwr=CompOUT%CmpOPwr   !RS: Debugging: Trying to pass through power from HPSim (9/3/14)
  !WRITE(LogFile,*) 'Supply Air Temp: ',tAoCoil !Mod
  !WRITE(LogFile,*) 'Return Air Temp: ',Node(3)%Temp
  !WRITE(LogFile,*) 'Qout: ',Qair
  
    !IF ((SysTimeTemp-SysTimeElapsed) .NE. 0) THEN   !RS: Debugging: Seeing if the data is passed correctly (8/20/14)
    !    WRITE(LogFile,FMT_153) 'HPSim Side',tAiCoil, hrAiCoil, tAoCoil, hrAoCoil, Qair
    !END IF
    !SysTimeTemp=SysTimeElapsed
    
END SUBROUTINE

SUBROUTINE GetEvapProp(Out1, Out2, Out3, Out4)
    !RS: Integration: Trying to carry over the properties to output
 REAL Out1, Out2, Out3, Out4
 
    Out1=tAiCoil
    Out2=hAiCoil
    Out3=tAoCoil
    Out4=hAoCoil

END SUBROUTINE

!************************************************************************

END MODULE EvaporatorMod