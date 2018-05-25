! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module contains a slew of routines used for evaporator and condenser coil calculations.
! These calculations include pressure drop, heat transfer coefficients, mass flow rates, and void factors.

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! The subroutines in this module represent evaporator and condenser coils, and the flow through them.
! Refer to the evaporator and condenser modules for information on whole components.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! Refrigerant and material properties are brought in, coil calculations run, and new refrigerant properties output.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no input or output files directly connected to this module.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! Coil, tube, slab, and circuit variables are all defined on the modular level.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 60 methods:
!    PUBLIC EPScalc -- Calculates the cross-flow heat exchanger effectiveness
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC UAcalc -- Calculates the UA for a heat exchanger segment
!      Never called
!    PUBLIC hcRefside -- Calculates the heat transfer coefficient on the refrigerant side
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC AirSideCalc -- Calculates the heat transfer resistance for the air side
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC CalcCoilArea -- Calculates the heat exchange heat transfer areas
!      Called by Condenser.f90
!      Called by Evaporator.f90
!      Called by AirSideCalc
!    PUBLIC CalcUA -- Calculates the UA value
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC MODdP -- Calculates the pressure drop in a module
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC Reynolds -- Calculates the Reynolds number
!      Called by Condenser.f90
!      Called by Evaporator.f90
!      Called by CoilCalc.f90
!    PUBLIC Prandtl -- Calculates the Prandtl number
!      Never called
!    PUBLIC manifold -- Calculates the condenser manifolds pressure drop
!      Called by Condenser.f90
!    PUBLIC returnbend -- Calculates the return bends pressure drop
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC Inventory -- Calculates refrigerant charge
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC GetRefName -- Sets refrigerant name based on refrigerant ID
!      Never called
!    PUBLIC GetRefID -- Sets refrigerant ID based on refrigerant name
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC CalcMeanProp -- Calculates the mean value of a property
!      Called by Condenser.f90
!      Called by Evaporator.f90
!      Called by CoilCalc.f90
!    PUBLIC hTPevapWattelet -- Calculates the refrigerant heat transfer coefficient in a two phase evaporating region using Wattelet method
!      Never called
!    PUBLIC MinimumFreeFlowArea -- Calculates the minimum free flow area for staggered tubes
!      Called by CoilCalc.f90
!    PUBLIC TWOPhasedPNino -- Calculates two phase pressure drop in microchannel tubes
!      Never called
!    PUBLIC TwoPhaseDPMoriyama -- Calculates two phase pressure drop using Moriyama method
!      Never called
!    PUBLIC TWOPhasedPChoi -- Calculates the pressure drop in condensing and evaporating two phase regions using Choi method
!      Never called
!    PUBLIC hSPDittus -- Calculates the refrigerant heat transfer coefficient for a single phase region using Dittus method
!      Called by hSPDittus
!    PRIVATE hSPPetukhov -- Calculates the refrigerant heat transfer coefficient for a single phase region using Petukhov method
!      Never called
!    PRIVATE hSPGnielinski -- Calculates the refrigerant heat transfer coefficient for a single phase region using Gnielinski method
!      Never called
!    PUBLIC  hTPDobson -- Calculates the refrigerant heat transfer coefficient for a two phase condensing region using Dobson method
!      Called by CoilCalc.f90
!    PRIVATE hTPCZ -- Calculates the refrigerant heat transfer coefficient for a two phase condensing region using Cavallini-Zecchin method
!      Never called
!    PUBLIC  hTPShahCond -- Calculates the refrigerant heat transfer coefficient for a two phase condensing region using Shah method
!      Called by CoilCalc
!    PRIVATE hTPShahEvap -- Calculates the refrigerant heat transfer coefficient in a two phase evaporating region using Shah method
!      Never called
!    PUBLIC  hTPCavallini -- Calculates the refrigerant heat transfer coefficient for a two phase condensing region using Cavallini method
!      Never called
!    PRIVATE hTPCavalliniAnnular -- Calculates the refrigerant heat transfer coefficient in an annular flow condensing region
!      Called by CoilCalc.f90
!    PRIVATE hTPCavalliniAnn_Strat -- Calculates the refrigerant heat transfer coefficient in an stratified annular flow condensing region
!      Called by CoilCalc.f90
!    PRIVATE hTPCavalliniStrat_Slug -- Calculates the refrigerant heat transfer coefficient in an stratified annular flow condensing region
!      Called by CoilCalc.f90
!    PRIVATE hTPconst -- Sets the refrigerant heat transfer coefficient as constant in a two phase region
!      Never called
!    PRIVATE ONEPhasedP -- Calculates the pressure drop in a single phase region
!      Called by CoilCalc.f90
!    PRIVATE FanningdP -- Calculates the friction pressure drop gradient in a single phase region
!      Called by CoilCalc.f90
!    PRIVATE FannFact -- Calculates the friction factor for a single phase region
!      Called by CoilCalc.f90
!    PUBLIC  TWOPhasedPLM -- Calculates the pressure drop in a two phase region using Lockhart-Matrinelli method
!      Never called
!    PUBLIC  TWOPhasedPSouza -- Calculates the pressure drop in condensing and evaporating two phase regions using Souza method
!      Called by CoilCalc.f90
!    PUBLIC  TWOPhasedPChang -- Calculates the pressure drop in condensing and evaporating two phase regions using Chang method
!      Never called
!    PRIVATE alphaCALC -- Calculates void fraction
!      Called by CoilCalc.f90
!    PRIVATE PHIcalc -- Calculates PHI for Martinelli friction pressure drop calculation
!      Called by CoilCalc.f90
!    PRIVATE PmomTWOphase -- Calculates the two phase momentum pressure change
!      Called by CoilCalc.f90
!    PRIVATE PmomIntegral -- Calculates intergral for evaluation of two phase momentum pressure change
!      Called by CoilCalc.f90
!    PRIVATE ElevdP -- Calculates pressure change due to elevation change
!      Called by CoilCalc.f90
!    PRIVATE Homo -- Calculates refrigerant charge using a closed homogeneous solution
!      Never called
!    PUBLIC  LockMartVoidFrac -- Calculates void fraction using Lockhart-Martinelli correlation
!      Never called
!    PUBLIC  HughmarkVoidFrac -- Calculates void fraction using Hughmark correlation
!      Never called
!    PUBLIC  GrahamVoidFrac -- Calculates void fraction using Graham et al. correlation
!      Called by CoilCalc.f90
!    PUBLIC  HarmsVoidFrac -- Calculates void fraction using Harms et al. correlation
!      Never called
!    PUBLIC  RouhanniVoidFrac -- Calculates void fraction using Rouhanni et al. correlation
!      Never called
!    PRIVATE CalcMassTP -- Calculates the charge inventory in a two phase region
!      Called by CoilCalc.f90
!    PUBLIC  CalcJfactor -- Calculates the j-factor for different fin types
!      Called by CoilCalc.f90
!    PRIVATE CalcFricfactor -- Calculates air side fraction factor for different fin types
!      Called by CoilCalc.f90
!    PRIVATE CalcFinEff -- Calculates the fin efficiency
!      Called by CalcFinEff
!    PUBLIC  CalcMaterialWeight -- Calculates the material weight in the coils
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PRIVATE CalcCustomAirHco -- Calculates the heat transfer coeffiecient on the air side using a custom curve fit
!      Never called
!    PRIVATE CalcCustomAirDP -- Calculates the pressure drop on the air side using a custom curve fit
!      Never called
!    PRIVATE CalcDPpenaltyFactor -- Calculates pressure drop penalty factor due to enhanced tube
!      Called by CoilCalc.f90
!    PRIVATE CalcHTCenhancementFactor -- Calculates heat transfer enhancement factor due to enhanced tube
!      Called by CoilCalc.f90
!    PUBLIC  UpdateRefMassFlowRate -- Updates refrigerant mass flow rate in each circuit
!      Called by Condenser.f90
!      Called by Evaporator.f90
!    PUBLIC  UpdateMCrefMassFlowRate -- Updates refrigerant mass flow rate in each microchannel slab
!      Called by Condenser.f90
!      Called by Evaporator.f90

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
! There are a number of subroutines that are never called. It would probably be a good idea to clean up anything
! that is out of date or not used. Also, some more documentation would be useful.

MODULE CoilCalcMod

!Contains all coil related calculations
USE DataSimulation

IMPLICIT NONE

PUBLIC !Share the following with Condenser and Evaporator, ISI - 06/05/07

!Data types 
TYPE ModInfo
  REAL Len             !Length, m
  REAL pRi             !Inlet ref. pressure
  REAL pRo             !Outlet ref. pressure
  REAL hRi             !Inlet ref. enthalpy
  REAL hRo             !Outlet ref. enthalpy
  REAL mAi             !air mass flow rate
  REAL tAi             !Inlet air temperature
  REAL tAo             !Outlet air temperature 
  REAL wAi             !Inlet air humidity ratio
  REAL wAo             !Outlet air humidity ratio
  REAL rhAi            !Inlet air humidity
  REAL rhAo            !Outlet air humidity
  REAL wbAi            !Inlet wet bulb air temperature
  REAL wbAo            !Outlet wet bulb air temperature
  REAL tSi             !Inlet surface temperature
  REAL tSo             !Outlet surface temperature
  REAL hci             !Module inside heat transfer coefficient
  REAL EFref           !Module inside heat transfer enhancement factor
  REAL hco             !Module outside heat transfer coefficient
  REAL ReVap           !Module Reynolds number vapor
  REAL ReLiq           !Module Reynolds number liquid
  REAL cAir            !Air specific heat
  REAL Rtube           !Tube resistance
  REAL Rair            !Air side resistance
  REAL Rfrost          !Frost resistance (Sankar)
  REAL Qmod            !Module heat transfer
  REAL QmodSens        !Module sensible heat transfer
  REAL SHR             !Sensible heat ratio
  INTEGER DryWet       !1=Wet surface; 2=Partially wet surface; 3=Dry surface
  REAL Mass            !Refrigerant mass
  REAL MassLiq         !Refrigerant mass liquid
  REAL MassVap         !Refrigerant mass vapor
  REAL VelDev          !Velocity deviation
  REAL VelAi           !air velcoity, m/s
  REAL Aface           !Face area, m^2
END TYPE ModInfo

TYPE TubeInfo
  INTEGER :: Even      !Even tubes flag
  INTEGER :: Fup       !Upper front tube#
  INTEGER :: Fdown     !Lower front tube#
  INTEGER :: Back      !Back tube flag: 1=yes; 2=no back tube
  INTEGER :: Empty     !Empty tube flag: 1=yes; 2=no
  INTEGER :: RowNum    !Row number of current tube
  INTEGER :: Nchannel  !Number of channels
  REAL    :: ID        !Tube inside diameter, m
  INTEGER :: NumOfMods !Number of tube modules/segments
  TYPE (ModInfo),ALLOCATABLE,DIMENSION(:) :: Seg  !Segment
END TYPE TubeInfo   

TYPE CktInfo
  INTEGER InJoin                   !# tubes joined at inlet
  INTEGER InSplit                  !# tubes split at inlet
  INTEGER OutJoin                  !# tubes joined at outlet
  INTEGER OutSplit                 !# tubes split at outlet
  INTEGER Ntube                    !# tubes
  INTEGER,ALLOCATABLE,DIMENSION(:) :: TubeSequence ! Tube sequence        
  REAL Conduct         !Refrigerant flow conductance
  REAL pRo             !Outlet pressure
  REAL pRi             !Inlet pressure
  REAL hRi             !Inlet enthalpy
  REAL hRo             !Outlet enthalpy
  REAL tRo             !Outlet temp.
  REAL tSat            !Saturated temp.
  REAL tSC             !Subcool temp.
  REAL tSH             !Superheat temp.
  REAL Qckt            !Circuit heat transfer
  REAL QcktSens        !Circuit sensible heat transfer
  REAL mRef            !Ref. mass flow rate 
  REAL mRefPrev        !Previous Ref. mass flow rate
  TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:) :: Tube !Tube
END TYPE CktInfo

TYPE NodeInfo
  INTEGER(2) Num       !Split or Joint node number
  REAL Pressure        !Pressure at node
END TYPE NodeInfo

TYPE PassInfo
  INTEGER Ntube        !# tubes
  REAL pRi             !Inlet pressure
  REAL pRo             !Outlet pressure
  REAL hRi             !Inlet enthalpy
  REAL hRo             !Outlet enthalpy
  REAL Qpass            !Coil pass heat transfer
  TYPE (TubeInfo),ALLOCATABLE,DIMENSION(:) :: Tube !Tube
END TYPE PassInfo

TYPE InletPassInfo
  INTEGER Ntube        !# tubes
  REAL mRef            !Ref. mass flow rate
  REAL pRi             !Inlet pressure
  REAL pRo             !Outlet pressure
  REAL hRi             !Inlet enthalpy
  REAL hRo             !Outlet enthalpy
  REAL Qpass           !Coil pass heat transfer, kW 
  REAL QpassSens       !Coil pass sensible heat transfer, kW
END TYPE InletPassInfo

TYPE SlabInfo
  INTEGER Npass        !# passes
  INTEGER Ninlet       !# inlets
  REAL tAi             !Inlet air temp
  REAL tAo             !Outlet air temp
  REAL rhAi            !Inlet air humidity
  REAL rhAo            !Outlet air humidity
  REAL hAi             !Inlet air enthalpy
  REAL hAo             !Outlet air enthalpy
  REAL mdot            !Refrigerant mass flow rate
  REAL pRi             !Inlet pressure
  REAL pRo             !Outlet pressure
  REAL hRi             !Inlet enthalpy
  REAL hRo             !Outlet enthalpy
  REAL tRi             !Inlet temp.
  REAL tRo             !Outlet temp.
  REAL xRi             !Inlet quality
  REAL xRo             !Outlet quality
  REAL tSC             !Subcool temp.
  REAL tSH             !Superheat temp.
  REAL Qslab           !Circuit heat transfer
  REAL Conduct         !Refrigerant flow conductance
  TYPE (PassInfo),ALLOCATABLE,DIMENSION(:) :: Pass !Pass 
  TYPE (InletPassInfo),ALLOCATABLE,DIMENSION(:) :: InletPass !Inlet pass 
END TYPE SlabInfo

TYPE SectionInfo !ISI - 09/10/07
  REAL Qsection        !Section capacity
  REAL QsectionSens    !Section sensible capacity
  REAL mRef            !Ref. mass flow rate
  REAL pRi             !Inlet pressure
  REAL pRo             !Outlet pressure
  REAL hRi             !Inlet enthalpy
  REAL hRo             !Outlet enthalpy
  INTEGER NumOfCkts    !Number of circuits
  INTEGER NcktLast     !Number of last circuits
  INTEGER NcktFirst    !Number of first circuits
  INTEGER Nnode        !Number of nodes
  LOGICAL IsInlet      !Inlet section flag
  INTEGER,ALLOCATABLE,DIMENSION(:) :: CktNum !Circuit number
  REAL,ALLOCATABLE,DIMENSION(:) :: mRefIter !Iteration ref. mass flow rate
  TYPE (CktInfo),ALLOCATABLE,DIMENSION(:) :: Ckt !Circuit info
  TYPE (NodeInfo),ALLOCATABLE,DIMENSION(:) :: Node !Node info
END TYPE SectionInfo

REAL,PARAMETER :: PI=3.14159265 !Pi

!Tube types
INTEGER,PARAMETER :: SMOOTH                    = 1     
INTEGER,PARAMETER :: MICROFIN                  = 2
INTEGER,PARAMETER :: HERRINGBONE               = 3
INTEGER,PARAMETER :: CROSSHATCH                = 4
INTEGER,PARAMETER :: HERRINGBONEWITHCROSSHATCH = 5

!Coil types
!INTEGER,PARAMETER :: CONDENSERCOIL  = 1
!INTEGER,PARAMETER :: EVAPORATORCOIL = 2
!INTEGER,PARAMETER :: HIGHSIDETUBE   = 3
!INTEGER,PARAMETER :: LOWSIDETUBE    = 4
!INTEGER,PARAMETER :: MCCONDENSER    = 5
!INTEGER,PARAMETER :: MCEVAPORATOR   = 6

!Fin types
INTEGER, PARAMETER :: PLAINFIN         = 1
INTEGER, PARAMETER :: WAVYFIN          = 2

!Module subroutines:
PUBLIC EPScalc
PUBLIC UAcalc
PUBLIC hcRefside 
PUBLIC AirSideCalc 
PUBLIC CalcCoilArea
PUBLIC CalcUA
PUBLIC MODdP
PUBLIC Reynolds 
PUBLIC Prandtl
PUBLIC manifold
PUBLIC returnbend 
PUBLIC Inventory
PUBLIC GetRefName
PUBLIC GetRefID
PUBLIC CalcMeanProp
PUBLIC hTPevapWattelet
PUBLIC MinimumFreeFlowArea
PUBLIC TWOPhasedPNino
PUBLIC TwoPhaseDPMoriyama
PUBLIC TWOPhasedPChoi
PUBLIC hSPDittus 
PRIVATE hSPPetukhov
PRIVATE hSPGnielinski
PUBLIC  hTPDobson
PRIVATE hTPCZ
PUBLIC  hTPShahCond
PRIVATE hTPShahEvap
PUBLIC  hTPCavallini
PRIVATE hTPCavalliniAnnular
PRIVATE hTPCavalliniAnn_Strat
PRIVATE hTPCavalliniStrat_Slug
PRIVATE hTPconst
PRIVATE ONEPhasedP
PRIVATE FanningdP
PRIVATE FannFact
PUBLIC  TWOPhasedPLM
PUBLIC  TWOPhasedPSouza
PRIVATE alphaCALC
PRIVATE PHIcalc
PRIVATE PmomTWOphase
PRIVATE PmomIntegral
PRIVATE ElevdP 
PRIVATE Homo 
PUBLIC  LockMartVoidFrac 
PUBLIC  HughmarkVoidFrac
PUBLIC  GrahamVoidFrac
PUBLIC  HarmsVoidFrac
PUBLIC  RouhanniVoidFrac
PRIVATE CalcMassTP
PUBLIC  CalcJfactor
PRIVATE CalcFricfactor
PRIVATE CalcFinEff
PUBLIC  CalcMaterialWeight
PRIVATE CalcCustomAirHco
PRIVATE CalcCustomAirDP
PRIVATE CalcDPpenaltyFactor
PRIVATE CalcHTCenhancementFactor
PUBLIC  UpdateRefMassFlowRate
PUBLIC  UpdateMCrefMassFlowRate

CONTAINS

!************************************************************************

SUBROUTINE EPScalc(cAir,cRef,UA,Cratio,NTU,EPS)

!------------------------------------------------------------------------
!Purpose:
!To calculate the cross-flow heat exchanger effectiveness.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL cAir !Air heat capacity rate, [kW/K]
REAL cRef !Refrigerant heat capacity rate, [kW/K]
REAL UA   !Overall heat transfer coefficient, [kW/m^2-K]

!Outputs:
REAL Cratio !Capacity rate ratio, = Cmin/Cmax, [-]
REAL NTU	!Number of transfer unit, [-]
REAL EPS    !Effectiveness, [-]

!Subroutine local parameters:
REAL,PARAMETER :: SMALL = 1E-6 !A small number
REAL,PARAMETER :: ZERO  = 0.0  !Zero

!Subroutine local variables:
REAL Cmin !Minimum heat capacity rate, [kW/K]
REAL Cmax !Maximum heat capacity rate, [kW/K]

!Flow:

	Cmin=MIN(cAir,cRef)
	IF (Cmin .NE. ZERO) THEN
		Cmax=MAX(cAir,cRef)
		Cratio=Cmin/Cmax
		NTU=UA/Cmin
		IF (Cratio .LT. SMALL) THEN
			EPS=1.-EXP(-NTU)
		ELSE IF (Cmin .EQ. cAir) THEN !Air side is mixed because of louver fins, ISI - 07/27/06
			!EPS=(1-EXP(-Cratio*(1-EXP(-NTU))))/Cratio !Cmin unmixed, Cmax mixed 
			EPS=1-EXP(-(1-EXP(-Cratio*NTU))/Cratio)   !Cmin mixed, Cmax unmixed
		ELSE IF (Cmin .EQ. cRef) THEN !Ref. side is unmixed because each circuit has different outlet temp., ISI - 07/27/06
			!EPS=1-EXP(-(1-EXP(-Cratio*NTU))/Cratio)   !Cmin mixed, Cmax unmixed
			EPS=(1-EXP(-Cratio*(1-EXP(-NTU))))/Cratio !Cmin unmixed, Cmax mixed 
		END IF
	ELSE
		EPS=1.0
	END IF

	RETURN

END SUBROUTINE EPScalc

!************************************************************************

SUBROUTINE UAcalc(Rtube,RairMod,RrefMod,UA)

!------------------------------------------------------------------------
!Purpose:
!To calculate the UA value for a heat exchanger segment
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE
      
!Subroutine passing variables:
!Inputs:
REAL Rtube   !Tube resistance, [(m^2-K)/kW]
REAL RairMod !Air side resistance, [(m^2-K)/kW]
REAL RrefMod !Ref side resistance, [(m^2-K)/kW]
!
!Outputs:
REAL UA !Overall heat transfer coefficient, [kW/m^2-K]

!Flow:

	UA=1./(RrefMod+Rtube+RairMod)

	RETURN

END SUBROUTINE UAcalc

!************************************************************************

SUBROUTINE hcRefside(CoilType,TubeType,ID,mRef,Qout, & !CoilType,TubeType,ID,ktube,mRef,Qout,AoMod,AiMod,hfg, &
                     xRi,xRo,vg,vf,muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV, &
					 Psat,Pcrit,Wabsolute,EF,hcRef) !MolWeight,Psat,Pcrit,Tsat,sigma,DT,Wabsolute,EF,hcRef)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant side heat transfer coefficient
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Eckels, S.J.; Doerr, T.M.; Pate, M.B. In-tube heat transfer and pressure
!drop of R-134a and ester lubricant mixtures in a smooth tube and a 
!micro-fin tube: part II - condensation. ASHRAE Transactions, v 100, n 2,  
!pp. 283-294, 1994.
!
!Eckels, S.J.; Doerr, T.M.; Pate, M.B. In-tube heat transfer and pressure 
!drop of R-134a and ester lubricant mixtures in a smooth tube and a 
!micro-fin tube: part I - evaporation. ASHRAE Transactions, v 100, n 2, 
!pp. 265-282, 1994.
!
!Outokumpu (2005). Tube side enhancement factors. Internal data.
!
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

USE OilMixtureMod
USE DataGlobals_HPSimIntegrated!RS: Debugging: Added for determination of sigma for testing (2/13/14)

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
							!7=Helical; 8=42F HXH
REAL ID			!Inside diameter, [m]
!REAL ktube		!Tube wall thermal conductivity, [kW/m-K]   !RS: Debugging: Extraneous
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL Qout		!Tube outside heat transfer, [kW]
REAL xRi		!Inlet quality of refrigerant, [-]
REAL xRo		!Outlet quality of refrigerant, [-]
REAL vg			!Specific volume of refrigerant vapor, [m^3/kg]
REAL vf			!Specific volume of refrigerant liquid, [m^3/kg]
REAL muRef		!Dynamic visocity of refrigerant, [Pa-s]
REAL mug		!Dynamic visocity of refrigerant vapor, [Pa-s]
REAL muf		!Dynamic visocity of refrigerant liquid, [Pa-s]
REAL kRef		!Thermal conductivity of refrigerant, [kW/m-K]
REAL kL			!Thermal conductivity of refrigerant liquid, [kW/m-K]
REAL kV			!Thermal conductivity of refrigerant vapor, [kW/m-K]
REAL CpRef		!Specific heat of refrigerant, [kJ/kg-K]
REAL CpL		!Specific heat of refrigerant liquid, [kJ/kg-K]
REAL CpV		!Specific heat of refrigerant vapor, [kJ/kg-K]
REAL Psat		!Saturated pressure of refrigerant, [kPa]
REAL Pcrit		!Critical pressure of refrigerant, [kPa]
REAL Wabsolute  !Absolute oil mass fraction

!Output:
REAL hcRef !Refrigerant side heat transfer coefficient, [kW/m^2-K]
REAL EF    !Heat transfer enhancement factor, [-]

!Subroutine local variables:
REAL xRef	  !Quality of refrigerant, [-]
REAL Gref     !Refrigerant mass flux, [kg/s-m^2]
!REAL Ref$   !RS: Debugging: Adding for testing of method (2/13/14)
!REAL sigma  !RS: Debugging: Cavallini (2/14/14)
!REAL hfg    !RS: Debugging: Cavallini (2/14/14)
!REAL DT     !RS: Debugging: Cavallini (2/14/14)

!RS: Debugging: Adding for testing of method (2/13/14)
!Ref$='R22'
!sigma=SigmaTest !RS: Debugging: Cavallini (2/14/14)
!hfg=hfgTest     !RS: Debugging: Cavallini (2/14/14)
!DT=DTTEST       !RS: Debugging: Cavallini (2/14/14)

!Flow:

	xRef=(xRi+xRo)/2 !ISI - 09/11/06

	Gref=mRef/(PI*ID**2/4)

	EF=1
	IF (xRef .LT. 1. .AND. xRef .GT. 0.) THEN
		IF (CoilType .EQ. CONDENSERCOIL) THEN
				CALL hTPDobson(CoilType,ID,xRef,mRef,vg,vf,mug,muf,kL,kV,CpL,CpV,hcRef)		
				!CALL hTPCavallini(CoilType,ID,xRef,mRef,vg,vf,mug,muf,kL,kV,CpL,CpV,sigma,hfg,DT,hcRef)
				!CALL hTPShahCond(ID,xRef,mRef,mug,muf,kL,kV,CpL,CpV,Psat,Pcrit,hcRef) 
				!CALL hTPCZ(ID,xRef,mRef,vg,vf,mug,muf,kL,CpL,hcRef)
		ELSEIF(CoilType .EQ. MCCONDENSER) THEN
				CALL hTPShahCond(ID,xRef,mRef,mug,muf,kL,kV,CpL,CpV,Psat,Pcrit,hcRef)
                !(ID,xRef,mRef,vgi,vfi,muRef,mug,muf,kL,kV,CpL,CpV,Psat,Pcrit,hTP)
				!CALL hTPDobson(CoilType,ID,xRef,mRef,vg,vf,mug,muf,kL,kV,CpL,CpV,hcRef)		
		ELSE
			IF (Qout .LT. 0) THEN
                Qout=0
            END IF

				hcRef=OilMixtureHTCevap(Gref,xRef,ID,muf,mug,1/vf,1/vg,kL*1e3,CpL*1e3,Wabsolute)*1e-3
	    
		END IF

	ELSE
		IF (xRef .GE. 1) THEN
		  muRef=mug
		  kRef=kV
		  CpRef=CpV
		ELSE
		  muRef=muf
		  kRef=kL
		  CpRef=CpL
		END IF
		CALL hSPDittus(CoilType,ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hcRef)
		!CALL hSPPetukhov(ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hcRef) !(ID,mRef,Quality,mu,mug,muf,k,Cp,hcRef)
		!CALL hSPGnielinski(ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hcRef)
	END IF

	IF (xRef .GE. 1 .OR. xRef .LE. 0) THEN
	    EF=1 
	ELSE !Only apply to two phase region, ISI - 09/19/20
		Gref=mRef/(ID**2*PI/4)
		CALL CalcHTCenhancementFactor(CoilType,TubeType,Gref,EF)
	END IF

	hcRef=hcRef*EF

	RETURN

END SUBROUTINE hcRefside

!************************************************************************

SUBROUTINE AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,tAiCoil,mAiCoil,rhoIn,rhoOut,Pt,Pl, &
                       Ltube,HtCoil,ID,OD,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg, &
					   Lcoil,AfCoil,AoCoil,AiCoil,FaceVel,hco,DP,hAir)

!AirSideCalc(CoilType,FinType,WetFlag,Nl,Nt,RowNum,tAiCoil,mAiCoil,rhoIn,rhoOut,Pt,Pl, &
!                       Ltube,HtCoil,ID,OD,NumOfChannels,Dchannel,TubeHeight,TubeDepth,FinThk,FinSpg, &
!					   CurveUnit,CurveTypeHTC,PowerAHTC,PowerBHTC,Poly1HTC,Poly2HTC, &
!					   Poly3HTC,Poly4HTC,CurveTypeDP,PowerADP,PowerBDP,Poly1DP,Poly2DP, &
!					   Poly3DP,Poly4DP,Lcoil,AfCoil, &
!					   AoCoil,AiCoil,FaceVel,hco,DP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the air side heat transfer resistance
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!McQuiston, F. C., Parker, J. D. and Spitler, J.D. (2000). Heating, 
!ventilating, and air conditioning - analysis and design, 4th Ed. 
!New York, NY: John Wiley & Sons, Inc.
!
!------------------------------------------------------------------------

USE AirPropMod      

IMPLICIT NONE

!Subroutine passing variables
!Inputs:
INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator
INTEGER FinType            !Fin type: 1-Plain; 2-Wavy; 3-Louver; 4-York 11-element
INTEGER WetFlag            !Wet surface flag: 1-Wet surface; Otherwise-Dry surface
INTEGER Nl                 !Number of rows
INTEGER Nt                 !Number of rows per tube
REAL tAiCoil   !Coil inlet air temperature, [C]
REAL mAiCoil   !Coil inlet air mass flow rate, [kg/s]
REAL rhoIn     !Inlet air density, [kg/m^3]
REAL rhoOut    !Outlet air density, [kg/m^3]
REAL Pt        !Tube spacing, [m]
REAL Pl        !Row spacing, [m]
REAL Ltube     !Single tube length, [m]
REAL HtCoil    !Coil height, [m]
REAL ID        !Tube inside diameter, [m]
REAL OD        !Tube outside diameter, [m]
INTEGER NumOfChannels !Number of channels
REAL Dchannel  !Channel diameter, [m]
REAL TubeHeight !Tube Height, [m]
REAL TubeDepth !Tube Depth, [m]
REAL FinThk    !Fin thickness, [m]
REAL FinSpg    !Fin spacing, [m]
REAL hAir   !RS: Replace: enthalpy of air for replacement (2/21/14)

!Outputs:
REAL Lcoil    !Coil tube length, [m]
REAL AfCoil   !Coil fin surface area, [m^2]
REAL AoCoil   !Coil outside surface area, [m^2]
REAL AiCoil   !Coil inside surface area, [m^2]
REAL FaceVel  !Face velocity, [m/s]
REAL hco      !Air side heat transfer coefficient, [kW/m^2-K]
REAL DP       !Air side pressure drop, [kPa]

!Subroutine local variables
REAL muA      !Dynamic viscosity, [Pa-s]
REAL kA       !Thermal conductivity, [kW/m-K]
REAL CPair    !Specific heat of air, [kJ/kg-K]
REAL PrAir    !Prandtl number of air, [-]
REAL rhoAvg   !Average air density, [kg/m^3]
REAL Dc       !Tube outside diameter including fin collar, [m]
REAL Rc       !Tube outside radius including fin collar, [m]
REAL HXdep    !Heat exchanger depth, [m]
REAL Jfactor  !J-factor, [-]
REAL Ffactor  !Friction factor, [-]
REAL Amin     !Minimum free flow area, [m^2]
REAL AbrCoil  !Bared tube coil surface area, [m^2]
REAL AmCoil   !Mean coil surface area, [m^2] 
REAL Acs      !Cross-sectional area, [m^2]
REAL AfrCoil  !Frontal area, [m^2]
REAL Gmax     !Maximum mass flux, [kg/s-m^2]
REAL ReDc     !Reynold number based on Dc, [-]
REAL RePl     !Reynold number based on Pl, [-]
REAL FinPitch !Fin pitch, [fins/m]
REAL sigma    !Ratio of Amin to AfrCoil, [-]
REAL Ki       !Contraction coefficient, [-]
REAL Ke       !Expansion coefficient, [-]
REAL humrat !RS: Replace: humidity ratio for replacement (2/21/14)

!Flow:

  !Obtain Air Properties
  humrat=HUMTH(tAiCoil,hAir)    !RS: Replace: humidity ratio for replacement (2/21/14)
  !muA = VISCA(REAL(tAiCoil)) !Viscosity !RS: Replace: VISCA (2/19/14)
  muA=Viscosity(tAiCoil,humrat) !RS: Replace: VISCA (2/19/14)
  !CPair=CPA(REAL(tAiCoil))  !RS: Replace: CPA (2/19/14)
  CPair=CPAirFunction(tAiCoil,humrat)  !RS: Replace: CPA (2/19/14)
  !kA = AKA(REAL(tAiCoil))   !Conductivity   !RS: Replace: AKA (2/19/14)
  kA=Conductivity(tAiCoil,humrat) !RS: Replace: VISCA (2/19/14)
  PrAir = muA*CPair/kA !Prandtl #

  !Outside radius, including collar  Sankar adjusted frost thickness
 IF(CoilType .EQ. CONDENSERCOIL) THEN      !RS: No place for the MC coil type
  Rc=OD/2+FinThk
 ELSE IF(CoilType .EQ. EVAPORATORCOIL) THEN
  Rc=OD/2+FinThk+FrostParam%Thickness
 END IF
 
  !Outside diameter, including collar
  Dc=Rc*2

  !Coil length
  Lcoil=Nl*Nt*Ltube

  !Coil frontal area
  AfrCoil=HtCoil*Ltube
  IF(CoolHeatModeFlag == 1) THEN            !Need to allow the option of a MicroChannel case
   IF(CoilType .EQ. CONDENSERCOIL) THEN
    CoilParams(2)%CoilFaceArea=AfrCoil
   ELSE IF(CoilType .EQ. EVAPORATORCOIL) THEN
    CoilParams(1)%CoilFaceArea=AfrCoil
   END IF
  ELSE IF(CoolHeatModeFlag == 0) THEN
   IF(CoilType .EQ. CONDENSERCOIL) THEN
    CoilParams(1)%CoilFaceArea=AfrCoil
   ELSE IF(CoilType .EQ. EVAPORATORCOIL) THEN
    CoilParams(2)%CoilFaceArea=AfrCoil
   END IF
  END IF    
  FinPitch=1/(FinSpg+FinThk)

  !HX depth	
  HXdep=Nl*Pl

  !Minimum flow area
  CALL MinimumFreeFlowArea(CoilType,Nl, Nt, Pl, Pt, OD,TubeHeight, FinThk, FinPitch, Ltube, Amin)

  CALL CalcCoilArea(CoilType,Nl,Nt,Pt,Pl,TubeDepth, &
                    Ltube,ID,OD,Dchannel,NumOfChannels, &
					FinThk,FinSpg,Lcoil,AfCoil, &
					AoCoil,AiCoil,AmCoil)
 
  AbrCoil=PI*Lcoil*Dc
  
  !Tube Cross sectional area
  Acs=PI*(ID**2)/4

  !Max. air mass flux 
  Gmax=mAiCoil/Amin

  !Reynolds # based on Dc
  ReDc=Gmax*Dc/muA

  !Reynolds # based on longitudinal pitch
  RePl=Gmax*Pl/muA
      
  SELECT CASE (FinType)
  CASE (1,2,3,4,6,7,8,9,10) !RS: Debugging: Allowing all fin types here
  !CASE (PLAINFIN,WAVYFIN)  !RS: No allowance for louver fins
  !CASE (PLAINFIN,WAVYFIN, 3)    !RS Comment: 3 is the Louver fin type, such as is used in the MC condenser case
  !RS Comment: Trying to implement the MC case(s), but it's unclear what's actually useful.

	FaceVel=mAiCoil/(rhoIn*AfrCoil)
	Gmax=mAiCoil/Amin
	ReDc=Gmax*Dc/muA
	
	RePl=Gmax*Pl/muA

	!J-factors
	CALL CalcJfactor(FinType,WetFlag,FinSpg,FinThk,HXdep,Nl, &
	                 Dc,Pt,Pl,Amin,AoCoil,ReDc,Jfactor)
!(CoilType,FinType,WetFlag,FinSpg,FinThk,Ltube,HXdep,Nl,Nt,RowNum, &
!                       OD,ID,Dc,Pt,Pl,TubeDepth,Aface,Amin,AoCoil,AbrCoil,ReDc,RePl,muAir, &
!					   cpAir,kAir,mAiCoil,FaceVel,Gmax,PrAir,jfactor)

	!Outside heat transfer coefficient
	hco=Jfactor*cpAir*Gmax*PrAir**(-0.667)

    !Friction factor
    CALL CalcFricfactor(FinType,WetFlag,FinSpg,FinThk,HXdep,Nl,Dc,Pt,Pl, &
					    AoCoil,ReDc,Amin,Ffactor)   !RS: Debugging: Removed extraneous variables

!(CoilType,FinType,WetFlag,FinSpg,FinThk,HXdep,Nl,Nt,Dc,OD,Pt,Pl, &
!                          Ltube,TubeDepth,AoCoil,AbrCoil,ReDc,RePl,Amin,AfrCoil,mAiCoil, &
!						  Gmax,muAir,rhoIn,rhoOut,Fricfactor)
					  
    IF (rhoOut .EQ. 0.0) THEN
	  rhoOut=rhoIn
    END IF
    rhoAvg=(rhoIn+rhoOut)/2
		  
    !Air side pressure drop 
	sigma=Amin/AfrCoil
	Ki=-0.404*sigma+0.494 !McQuiston
	Ke=-1.272*sigma+0.8726 !McQuiston
	DP=(Gmax**2)/(2*rhoIn)*((Ki+1-sigma**2)+2*(rhoIn/rhoOut-1)+Ffactor*(AoCoil/Amin)*(rhoIn/rhoAvg) &
	                        -(1-sigma**2-Ke)*rhoIn/rhoOut)*1E-3

  END SELECT
  
  !RS: Attempt at error handling (8/1/12)
  !IF (FINTYPE .NE. PLAINFIN .AND. FINTYPE .NE. WAVYFIN) THEN
  !    CALL ShowSevereError('The fin type selected cannot be handled by the program')
  !END IF
  
  RETURN

END SUBROUTINE AirSideCalc

!************************************************************************

SUBROUTINE CalcCoilArea(CoilType,Nl,Nt,Pt,Pl,TubeDepth, &
                      Ltube,ID,OD,Dchannel,NumOfChannels, &
					  FinThk,FinSpg,Lcoil,AfCoil, &
					  AoCoil,AiCoil,AmCoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate the heat exchanger heat transfer areas
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!November 2006
!
!Reference:
!none
!
!------------------------------------------------------------------------
      
IMPLICIT NONE

!Subroutine passing variables
!Inputs:
INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

INTEGER Nl                 !Number of rows
INTEGER Nt                 !Number of rows per tube
REAL Pt        !Tube spacing, [m]
REAL Pl        !Row spacing, [m]
REAL TubeDepth !Tube depth, [m]
REAL Ltube     !Single tube length, [m]
REAL ID        !Tube inside diameter, [m]
REAL OD        !Tube outside diameter, [m]
REAL Dchannel  !Microchannel diameter, [m]
INTEGER NumOfChannels !Number of channels per tube
REAL FinThk    !Fin thickness, [m]
REAL FinSpg    !Fin spacing, [m]

!Outputs:
REAL Lcoil    !Coil tube length, [m]
REAL AfCoil   !Coil fin surface area, [m^2]
REAL AoCoil   !Coil outside surface area, [m^2]
REAL AiCoil   !Coil inside surface area, [m^2]
REAL AmCoil   !Mean coil surface area, [m^2]

!Subroutine local variables
REAL Dc       !Tube outside diameter including fin collar, [m]
REAL Rc       !Tube outside radius including fin collar, [m]
REAL AbsCoil  !Coil base surface area, [m^2]
REAL AbrCoil  !Bared tube coil surface area, [m^2]
REAL AfrCoil  !Frontal area, [m^2]
REAL FinPitch !Fin pitch, [fins/m]
REAL FinHeight !Fin height, [m]
REAL TubeHeight !Tube height, [m]

REAL CoilDepth   !Coil depth, m
REAL CoilWidth   !Coil width, m
REAL CoilHeight   !Coil height, m

!Flow:

  !Coil length
  Lcoil=Nl*Nt*Ltube
  
  FinPitch=1/(FinSpg+FinThk)  

  IF (CoilType .EQ. MCCONDENSER .OR. &
      CoilType .EQ. MCEVAPORATOR) THEN !Microchannel coil

	  TubeHeight=OD
	  FinHeight=Pt-TubeHeight
	  AbrCoil=(PI*TubeHeight+TubeDepth*2)*Lcoil
	  AfCoil=2*(FinThk+TubeDepth)*FinHeight*FinPitch*Lcoil
	  AbsCoil=AbrCoil-FinThk*TubeDepth*FinPitch*Lcoil
	  AoCoil=AbsCoil+AfCoil
	  AiCoil=PI*Dchannel*NumOfChannels*Lcoil
	  AmCoil=(AiCoil+AbrCoil)/2      !Module tube mean surface area    
  ELSE
	  !Outside radius, including collar
    IF(CoilType==CONDENSERCOIL) THEN 
     Rc=OD/2+FinThk
    ELSE IF(CoilType==EVAPORATORCOIL) THEN
     Rc=OD/2+FinThk+FrostParam%Thickness
    END IF 

	  !Outside diameter, including collar
	  Dc=Rc*2

	  !Coil frontal area
	  AfrCoil=Pt*Nt*Ltube

	  !HX depth	
	  !HXdep=Nl*Pl

      !Bare tube surface area
	  !Base surface area
	  !Fin surface area
	  IF(CoilType == CONDENSERCOIL) THEN
	   AbsCoil=PI*Dc*Lcoil*(1-FinThk*FinPitch)
	   AbrCoil=PI*Lcoil*Dc   !Sankar Changed
	   AfCoil=2*(Nt*Pt*Nl*Pl-PI*Dc**2/4*Nl*Nt+FinThk*(Nt*Pt+Nl*Pl))*FinPitch*Ltube
	  ELSE IF(CoilType == EVAPORATORCOIL) THEN
	   AbsCoil=PI*Dc*Lcoil*(1-(FinThk+2*FrostParam%Thickness)*FinPitch)
	   AbrCoil=PI*Lcoil*Dc  !Sankar Changed
	   AfCoil=2*(Nt*Pt*Nl*Pl-PI*Dc**2/4*Nl*Nt+(FinThk+FrostParam%Thickness)*(Nt*Pt+Nl*Pl))*FinPitch*Ltube
      END IF
      
      IF(CoolHeatModeFlag == 1) THEN
       IF(CoilType==1) THEN
        CoilParams(2)%CoilFinArea=AfCoil
        CoilParams(2)%FinPitch=FinPitch
        CoilParams(2)%FinThickness=FinThk
       ELSE IF(CoilType == 2) THEN
        CoilParams(1)%FinPitch=FinPitch
        CoilParams(1)%FinThickness=FinThk
        CoilParams(1)%CoilFinArea=AfCoil
       END IF
      ELSE IF(CoolHeatModeFlag == 0) THEN
       IF(CoilType == 1) THEN
        CoilParams(1)%CoilFinArea=AfCoil
        CoilParams(1)%FinPitch=FinPitch
        CoilParams(1)%FinThickness=FinThk
       ELSE IF(CoilType == 2) THEN
        CoilParams(2)%CoilFinArea=AfCoil
        CoilParams(2)%FinPitch=FinPitch
        CoilParams(2)%FinThickness=FinThk
       END IF
      END IF    

	  !Total outside surface area
	  AoCoil=AbsCoil+AfCoil

	  !Weber 1991
	  CoilDepth=Nl*Pl
	  CoilWidth=Ltube
	  CoilHeight=Nt*Pt
	  !AbsCoil=PI*Dc*CoilWidth*(1-Finthk)*Nl*Nt+2*(CoilDepth*CoilHeight-PI*Dc**2/4*Nl*Nt)
	  !AfCoil=2*(CoilDepth*CoilHeight-PI*Dc**2/4*Nl*Nt+CoilHeight*Finthk)*FinPitch*CoilWidth
	  !AoCoil=AbsCoil+AfCoil

	  !Total inside surface area
	  AiCoil=PI*ID*Lcoil

	  !Tube Cross sectional area
	  !Acs=PI*(ID**2)/4

	  AmCoil=(AiCoil+AbrCoil)/2 

	  IF(CoilType==CONDENSERCOIL) THEN
	   CondTubeArea = AbsCoil
	   CondFinArea = AfCoil
	  ELSE IF(CoilType==EVAPORATORCOIL) THEN
	   EvapTubeArea = AbsCoil
	   EvapFinArea = AfCoil
	   EvapTotArea = AoCoil
	   EvapBareArea= AbrCoil
	  END IF

  END IF

RETURN

END SUBROUTINE CalcCoilArea

!************************************************************************

SUBROUTINE CalcCustomAirHco(Unit,CurveType,PowerA,PowerB, &
                            Poly1,Poly2,Poly3,Poly4,FaceVel,hco)

!------------------------------------------------------------------------
!Purpose:
!To calculate the air side heat transfer coefficient using custom 
!curve fit coefficients
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!August 2005
!
!Reference:
!none
!
!------------------------------------------------------------------------
      
IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER Unit      !Unit of of the curve fit coefficients: 1=SI, 2=IP
INTEGER CurveType !Curve fit type: 1-power; 2-polynomial
REAL PowerA  !Power fit coefficient A
REAL PowerB  !Power fit coefficient B
REAL Poly1   !Polynomial fit coefficient 1
REAL Poly2   !Polynomial fit coefficient 2
REAL Poly3   !Polynomial fit coefficient 3
REAL Poly4   !Polynomial fit coefficient 4
REAL FaceVel !Face velociy, [m/s]

!Output:
REAL hco	 !Heat transfer coefficient, [kW/m2-K]

!Flow:

	IF (CurveType .EQ. 1) THEN
		hco=PowerA*FaceVel**PowerB
	ELSE
		hco=Poly1+Poly2*FaceVel+Poly3*FaceVel**2+Poly4*FaceVel**3
	END IF

	IF (Unit .EQ. 1) THEN
		hco=hco/1000       !from W/m2-K to kW/m2-K
	ELSE
		hco=hco*5.678/1000 !from Btu/hr-ft2-K to kW/m2-K
	END IF

RETURN

END SUBROUTINE CalcCustomAirHco

!************************************************************************

SUBROUTINE CalcCustomAirDP(Unit,CurveType,PowerA,PowerB, &
                           Poly1,Poly2,Poly3,Poly4,FaceVel,DP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the air side pressure drop using custom 
!curve fit coefficients
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!August 2005
!
!Reference:
!none
!
!------------------------------------------------------------------------
      
IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER Unit      !Unit of of the curve fit coefficients: 1=SI, 2=IP
INTEGER CurveType !Curve fit type: 1-power; 2-polynomial
REAL PowerA  !Power fit coefficient A
REAL PowerB  !Power fit coefficient B
REAL Poly1   !Polynomial fit coefficient 1
REAL Poly2   !Polynomial fit coefficient 2
REAL Poly3   !Polynomial fit coefficient 3
REAL Poly4   !Polynomial fit coefficient 4
REAL FaceVel !Face velociy, [m/s]

!Output:
REAL DP		 !Air side pressure drop, kPa

!Flow:

	IF (CurveType .EQ. 1) THEN
		DP=PowerA*FaceVel**PowerB
	ELSE
		DP=Poly1+Poly2*FaceVel+Poly3*FaceVel**2+Poly4*FaceVel**3
	END IF

	IF (Unit .EQ. 1) THEN
		DP=DP/1000 !from Pa to kPa
	ELSE 
		DP=DP*249/1000 !from in.wg to kPa
	END IF

RETURN

END SUBROUTINE CalcCustomAirDP


!************************************************************************

SUBROUTINE CalcUA(CoilType,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,OD,TubeThk,TubeDepth, &
                 hco,hci,AfMod,AoMod,AiMod,AmMod,UA,Rair,Rref,Rtube,FinEff,SurfEff)

!CalcUA(CoilType,WetFlag,Kfin,FinThk,FinHeight,Ktube,Pt,Pl,OD,TubeThk,TubeDepth, &
!                  RowNum,tAiMod,hAiMod,hco,hci,AfMod,AoMod,AiMod,AmMod,UA,Rair,Rref,Rtube,FinEff,SurfEff)

!------------------------------------------------------------------------
!Purpose:
!To calculate UA value
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!
!
!------------------------------------------------------------------------
      
IMPLICIT NONE

!Subroutine passing varibles:
!Inputs:
INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator
REAL Kfin     !Fin conductivity, [kW/m-K]
REAL FinThk   !Fin thickness, [m]
REAL FinHeight !Fin height, [m]
REAL Ktube    !Tube conductivity, [kW/m-K]
REAL TubeThk  !Tube thickness, [m]
REAL OD       !Outside diameter of tube, [m]
REAL TubeDepth !Tube depth, [m]
REAL Pt       !Tube spacing, [m]
REAL Pl       !Row spacing, [m]
REAL hco      !Outside heat transfer coefficient, [kW/m^2-K]
REAL hci      !Inside heat transfer coefficient, [kW/m^2-K]
REAL AfMod    !Module fin surface area, [m^2]
REAL AoMod    !Module outside surface area, [m^2]
REAL AiMod    !Module inside surface area, [m^2]
REAL AmMod    !Module mean surface area, [m^2]

!Outputs:
REAL UA       !UA value, [kW/K]
REAL Rair     !Air side resistance, [kW/K]
REAL Rref     !Refrigerant side resistance, [kW/K]
REAL Rtube    !Tube side resistance, [kW/K]
REAL FinEff   !Fin efficiency, [-]
REAL SurfEff  !Surface efficiency, [-]

!Subroutine local variables:
REAL Rc !Tube outside radius including fin collar, [m]

REAL KFinFrost
REAL FinFrostThk
REAL KFrost
REAL FrostThk

!Flow:

  !Outside radius, including fin thickness
  IF(CoilType==EVAPORATORCOIL) THEN
   Rc=OD/2+FinThk+FrostParam%Thickness
  ELSE
   Rc=OD/2+FinThk
  END IF
  
  !Fin efficiency
  CALL CalcFinEff(CoilType,Kfin,FinThk,FinHeight,Rc,TubeDepth,Pt,Pl,hco,FinEff)

  !Surface effectiveness
  SurfEff=1-(AfMod/AoMod)*(1-FinEff)

  !Dry fin total air side resistance  
  Rair=1/(AoMod*SurfEff*hco)

  !Refrigerant side resistance
  Rref=1/(hci*AiMod)
  
  !Tube resistance   Sankar-Change the KTube
  IF(CoilType==EVAPORATORCOIL) THEN
   Kfrost=FrostParam%Conductivity
   FrostThk=2*FrostParam%Thickness      
   KfinFrost=(KTube*TubeThk+Kfrost*FrostThk)/(TubeThk+FrostThk)
   !FinFrostThk=FinThk+FrostThk !ISI - 12/24/2009
   FinFrostThk=TubeThk+FrostThk
   Rtube=FinFrostThk/(KfinFrost*AmMod) 
  ELSE
   Rtube=TubeThk/(Ktube*AmMod) 
  END IF

  !UA
  UA=1/(Rref+Rtube+Rair)

  RETURN

END SUBROUTINE CalcUA

!************************************************************************

SUBROUTINE hSPDittus(CoilType,ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hSP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in single phase
!region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Dittus, Boelter
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
REAL ID         !Tube inside diameter, [m]
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL xRef		!Refrigerant quality, [-]
REAL muRef		!Dynamic viscosity, [Pa-s]
REAL mug		!Dynamic viscosity of vapor, [Pa-s]
REAL muf		!Dynamic viscosity of liquid, [Pa-s]
REAL kRef		!Thermal conductivity, [kW/m-K] 
REAL CpRef		!Specific heat, [kJ/kg-K]

!Output:
REAL hSP        !Single phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables:
REAL ReLiq      !Liquid Reynold number
REAL ReVap      !Vapor Reynold number
REAL Re         !Reynold number
REAL Pr         !Prandtl number 
REAL Nu         !Nusselt number

!Flow:

  CALL Reynolds(ID,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
  IF (xRef .LE. 0.0) THEN
      Re=ReLiq
  ENDIF
  IF (xRef .GE. 1.0) THEN
      Re=ReVap
  ENDIF
  Pr=muRef*CpRef/kRef

  IF (CoilType .EQ. CONDENSERCOIL .OR. &
      CoilType .EQ. MCCONDENSER .OR. &
	  CoilType .EQ. HIGHSIDETUBE) THEN

	Nu=0.023*(Re**0.8)*(Pr**0.3) !Condenser
  ELSE
	Nu=0.023*(Re**0.8)*(Pr**0.4) !Evaporator
  END IF

  hSP=Nu*kRef/ID

  RETURN

END SUBROUTINE hSPDittus

!************************************************************************

SUBROUTINE hSPPetukhov(ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hSP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in single phase
!region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Pedtukhov
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID         !Tube inside diameter, [m]
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL xRef		!Refrigerant quality, [-]
REAL muRef		!Dynamic viscosity, [Pa-s]
REAL mug		!Dynamic viscosity of vapor, [Pa-s]
REAL muf		!Dynamic viscosity of liquid, [Pa-s]
REAL kRef		!Thermal conductivity, [kW/m-K] 
REAL CpRef		!Specific heat, [kJ/kg-K]

!Output:
REAL hSP        !Single phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables:
REAL ReLiq      !Liquid Reynold number
REAL ReVap      !Vapor Reynold number
REAL Re         !Reynold number
REAL Pr         !Prandtl number 
REAL Nu         !Nusselt number
REAL ff			!Intermediate variable

!Flow:

  CALL Reynolds(ID,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
  IF (xRef .LT. 0.0) THEN
      Re=ReLiq
  END IF
  IF (xRef .GT. 1.0) THEN
      Re=ReVap
  END IF
  Pr=muRef*CpRef/kRef
  ff=(1.82*LOG10(Re)-1.64)**(-2)
  Nu=ff/8*Re*Pr/(1.07+12.7*(ff/8)**0.5*(Pr**0.67-1))
  hSP=Nu*kRef/ID

  RETURN

END SUBROUTINE hSPPetukhov

!************************************************************************

SUBROUTINE hSPGnielinski(ID,mRef,xRef,muRef,mug,muf,kRef,CpRef,hSP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in single phase
!region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Gnielinski, V. New equation for heat and mass transfer in turbulent 
!pipe and channel flow. Int. Chem. Eng. 16. pp. 359-368. 1976.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID         !Tube inside diameter, [m]
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL xRef		!Refrigerant quality, [-]
REAL muRef		!Dynamic viscosity, [Pa-s]
REAL mug		!Dynamic viscosity of vapor, [Pa-s]
REAL muf		!Dynamic viscosity of liquid, [Pa-s]
REAL kRef		!Thermal conductivity, [kW/m-K] 
REAL CpRef		!Specific heat, [kJ/kg-K]

!Output:
REAL hSP        !Single phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables:
REAL ReLiq      !Liquid Reynold number
REAL ReVap      !Vapor Reynold number
REAL Re         !Reynold number
REAL Pr         !Prandtl number 
REAL Nu         !Nusselt number
REAL ff			!Intermediate variable

!Flow:

  CALL Reynolds(ID,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
  IF (xRef .LE. 0.0) THEN
      Re=ReLiq
  END IF
  IF (xRef .GE. 1.0) THEN
      Re=ReVap
  END IF
  Pr=muRef*CpRef/kRef
  IF (Re .LT. 0) THEN
      Re=1e-6
  END IF
  ff=(1.58*LOG(Re)-3.28)**(-2)
  Nu=(Re-1000)*Pr*(ff/2)/(1+12.7*SQRT(ff/2)*(Pr**(0.67)-1))

  hSP=Nu*kRef/ID

  RETURN

END SUBROUTINE hSPGnielinski

!************************************************************************

SUBROUTINE hTPDobson(CoilType,ID,xRef,mRef,vgi,vfi,mug,muf,kL,kV,CpL,CpV,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Dobson, M.K. and Chato, J.C. Condensation in smooth horizontal tubes. 
!J. of heat transfer, transactions of ASME. 120:193-312, 1998.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL vgi   !Vapor specific volume, [m^3/kg]
REAL vfi   !Liquid specific volume, [m^3/kg]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL kV    !Vapor conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL CpV   !Vapor specific heat, [kJ/kg-K]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local parameter:
REAL,PARAMETER :: xMin=0.1
REAL,PARAMETER :: xMax=0.9

!Subroutine local variables
REAL Acs    !Cross sectional area, [m^2]  
REAL Gref   !Refrigerant mass flux, [kg/s-m^2]
REAL rhog   !Vapor density, [kg/m^3]
REAL rhof   !Liquid density, [kg/m^3]
REAL Xtt    !Lockhart-Martinelli parameter
REAL PrL    !Liquid Prandtl number
REAL ReL    !Liquid Reynold number
REAL hLmin  !Liquid heat transfer coefficient for xMin, [kW/m^2-K]
REAL hLiq   !Liquid heat transfer coefficient, [kW/m^2-K]
REAL hVap   !Vapor heat transfer coefficient, [kW/m^2-K]
REAL hTPmin !Minimum two-phase heat transfer coefficient, [kW/m^2-K]
REAL hTPmax !Maximum two-phase heat transfer coefficient, [kW/m^2-K]
REAL VelGas !Wallis dimensionless gas velocity, [-]

!Flow:

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs
  rhog=1.0/vgi
  rhof=1.0/vfi

  VelGas=Gref*xRef/SQRT(9.8*ID*rhog*(rhof-rhog))

  IF (xRef .GE. xMin .AND. xRef .LE. xMax) THEN !ISI - 09/27/06
    Xtt=(rhog/rhof)**0.5*(muf/mug)**0.125*((1-xRef)/xRef)**0.875
    PrL=CpL*muf/kL
    ReL=Gref*ID*(1.0-xRef)/muf
	hLmin=0.023*ReL**0.8*PrL**0.3*kL/ID
	hTP=hLmin*2.61/Xtt**0.805   
  ELSEIF (xRef .LT. xMin) THEN
    Xtt=(rhog/rhof)**0.5*(muf/mug)**0.125*((1-xMin)/xMin)**0.875
    PrL=CpL*muf/kL
    ReL=Gref*ID*(1.0-xMin)/muf
  	hLmin=0.023*ReL**0.8*PrL**0.3*kL/ID
    hTPmin=2.61*hLmin/Xtt**0.805
   	CALL hSPDittus(CoilType,ID,mRef,0.00,muf,mug,muf,kL,CpL,hLiq)
  	hTP=xRef*(hTPmin-hLiq)/xMin+hLiq
  ELSEIF (xRef .GT. xMax) THEN
    Xtt=(rhog/rhof)**0.5*(muf/mug)**0.125*((1-xMax)/xMax)**0.875
    PrL=CpL*muf/kL
    ReL=Gref*ID*(1.0-xMax)/muf
  	hLmin=0.023*ReL**0.8*PrL**0.3*kL/ID
    hTPmax=2.61*hLmin/Xtt**0.805
  	CALL hSPDittus(CoilType,ID,mRef,1.00,mug,mug,muf,kV,CpV,hVap)
  	hTP=(xRef-1)/(xMax-1)*(hTPmax-hVap)+hVap
  END IF
	      
  RETURN

END SUBROUTINE hTPDobson

!************************************************************************

SUBROUTINE hTPCZ(ID,xRef,mRef,vgi,vfi,mug,muf,kL,CpL,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Cavallini-Zecchin
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL vgi   !Vapor specific volume, [m^3/kg]
REAL vfi   !Liquid specific volume, [m^3/kg]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables:
REAL Acs   !Cross sectional area, [m^2]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL PrL   !Liquid Prandtl number
REAL ReL   !Liquid Reynold number
REAL ReV   !Vapor Reynold number
REAL ReEq  !Equivalent Reynold number

!Flow:

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs
      
  rhog=1.0/vgi
  rhof=1.0/vfi

  PrL=CpL*muf/kL
  ReV=Gref*ID*xRef/mug  
  ReL=Gref*ID*(1.0-xRef)/muf

  ReEq=ReL+ReV*(mug/muf)*(rhof/rhog)**0.5
  hTP=0.05*ReEq**0.8*PrL**0.33*(kL/ID)

  RETURN

END SUBROUTINE hTPCZ

!************************************************************************

SUBROUTINE hTPShahCond(ID,xRef,mRef,mug,muf,kL,kV,CpL,CpV,Psat,Pcrit,hTP) !(ID,xRef,mRef,vgi,vfi,muRef,mug,muf,kL,kV,CpL,CpV,Psat,Pcrit,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Shah, M.M. (1979). A general correlation for heat transfer during film
!condensation inside pipes. International Journal of heat and mass transfer. 
!88, pp 185-196.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL kV    !Vapor conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL CpV   !Vapor specific heat, [kJ/kg-K]
REAL Psat  !Saturation pressure, [kPa]
REAL Pcrit !Critical pressure, [kPa]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local parameter:
INTEGER,PARAMETER :: CoilType=1 !Condenser
REAL,PARAMETER    :: xMax=0.95 !Maximum quality

!Subroutine local variables:
REAL hliq  !Liquid heat transfer coefficient, [kW/m^2-K] 
REAL h1    !Intermediate heat transfer coefficient, [kW/m^2-K] 
REAL ZZ    !Intermediate variable
REAL Pred  !Reduced pressure, [kPa]
REAL hTPmax
REAL hVap

!Flow:

  CALL hSPDittus(CoilType,ID,mRef,0.00,muf,mug,muf,kL,CpL,hliq)
  
  IF (xRef .LT. xMax) THEN
    h1=hliq*(1-xRef)**0.8
    Pred=Psat/Pcrit
    ZZ=(1/xRef-1)**0.8*Pred**0.4
    hTP=h1*(1+3.8/ZZ**0.95)
  ELSEIF (xRef .GE. xMax) THEN
    h1=hliq*(1-xMax)**0.8
    Pred=Psat/Pcrit
    ZZ=(1/xMax-1)**0.8*Pred**0.4
    hTPmax=h1*(1+3.8/ZZ**0.95)
	CALL hSPDittus(CoilType,ID,mRef,1.00,mug,mug,muf,kV,CpV,hVap)
	hTP=(xRef-1)/(xMax-1)*(hTPmax-hVap)+hVap
  END IF

  RETURN

END SUBROUTINE hTPShahCond

!************************************************************************

SUBROUTINE hTPCavallini(CoilType,ID,xRef,mRef,vgi,vfi,mug,muf,kL,kV,CpL,CpV,sigma,hfg,DT,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!Sept 2006
!
!Reference:
!Cavallini A, Del Col D, Doretti L, Longo G.A., Rossetto, L. "Intube
!condensation of halogenated refrigerants." ASHRAE transactions, 
!2002,108(1), pp.146-161
!
!------------------------------------------------------------------------

IMPLICIT NONE


!Subroutine passing variables:
!Inputs:
INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL vgi   !Vapor specific volume, [m^3/kg]
REAL vfi   !Liquid specific volume, [m^3/kg]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL kV    !Vapor conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL CpV   !Vapor specific heat, [kJ/kg-K]
REAL sigma !Surface tension, [N/m]  
REAL hfg   !Latent heat of condensation, [kJ/kg]
REAL DT    !Temperature difference between tube wall and saturate vapor, [C]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local parameter:
REAL,PARAMETER :: xMax=0.9

!Subroutine local variables
REAL Acs    !Cross sectional area, [m^2]  
REAL Gref   !Refrigerant mass flux, [kg/s-m^2]
REAL JG     !Dimensionless vapor velocity
REAL Xtt    !Lockhart-Martinelli parameter
REAL rhog   !Vapor density, [kg/m^3]
REAL rhof   !Liquid density, [kg/m^3]
REAL hTPmax !Two-phase heat transfer coefficient at xMax, [kW/m^2-K]
REAL hVap   !Vapor heat transfer coefficient, [kW/m^2-K] 

!Flow:

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs
  rhog=1.0/vgi
  rhof=1.0/vfi

  IF (xRef .LT. xMax) THEN
	  JG=xRef*Gref/(9.8*ID*rhog*(rhof-rhog))**0.5
	  Xtt=(muf/mug)**0.1*(rhog/rhof)**0.5*((1-xRef)/xRef)**0.9

	  IF (JG .GT. 2.5) THEN !Annular flow

		CALL hTPCavalliniAnnular(ID,xRef,Gref,rhog,rhof,mug,muf, &
								 kL,CpL,sigma,hTP)
        !hTPCavalliniAnnular(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
        !                       kL,kV,CpL,CpV,sigma,hTP)
  
	  ELSEIF (JG .LT. 2.5 .AND. Xtt .LT. 1.6) THEN !Annular-stratified flow

		CALL hTPCavalliniAnn_Strat(ID,xRef,Gref,rhog,rhof,mug,muf, &
								   kL,CpL,sigma,hfg,DT,JG,hTP)
        !(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                                 kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)
  
	  ELSEIF (JG .LT. 2.5 .AND. Xtt .GT. 1.6) THEN !Stratified and slug flow

		CALL hTPCavalliniStrat_Slug(ID,xRef,Gref,rhog,rhof,mug,muf, &
									kL,CpL,sigma,hfg,DT,JG,hTP)
        !(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
        !                          kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)
  
	  END IF

  ELSE
	  JG=xMax*Gref/(9.8*ID*rhog*(rhof-rhog))**0.5
	  Xtt=(muf/mug)**0.1*(rhog/rhof)**0.5*((1-xMax)/xMax)**0.9

	  IF (JG .GT. 2.5) THEN !Annular flow

		CALL hTPCavalliniAnnular(ID,xMax,Gref,rhog,rhof,mug,muf, &
								 kL,CpL,sigma,hTPmax)
                !hTPCavalliniAnnular(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
        !                       kL,kV,CpL,CpV,sigma,hTP)
  
	  ELSEIF (JG .LT. 2.5 .AND. Xtt .LT. 1.6) THEN !Annular-stratified flow

		CALL hTPCavalliniAnn_Strat(ID,xMax,Gref,rhog,rhof,mug,muf, &
								   kL,CpL,sigma,hfg,DT,JG,hTPmax)
        !(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                                 kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)
  
	  ELSEIF (JG .LT. 2.5 .AND. Xtt .GT. 1.6) THEN !Stratified and slug flow

		CALL hTPCavalliniStrat_Slug(ID,xMax,Gref,rhog,rhof,mug,muf, &
									kL,CpL,sigma,hfg,DT,JG,hTPmax)
                !(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
        !                          kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)
  
	  END IF

	CALL hSPDittus(CoilType,ID,mRef,1.00,mug,mug,muf,kV,CpV,hVap)
	hTP=(xRef-1)/(xMax-1)*(hTPmax-hVap)+hVap
  END IF
      
  RETURN

END SUBROUTINE hTPCavallini

!************************************************************************

SUBROUTINE hTPCavalliniAnnular(ID,xRef,Gref,rhog,rhof,mug,muf, &
                               kL,CpL,sigma,hTP)
!hTPCavalliniAnnular(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                               kL,kV,CpL,CpV,sigma,hTP)
!
!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in annular flow
!condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!Sept 2006
!
!Reference:
!Cavallini A, Del Col D, Doretti L, Longo G.A., Rossetto, L. "Intube
!condensation of halogenated refrigerants." ASHRAE transactions, 
!2002,108(1), pp.146-161
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL sigma !Surface tension, [N/m]  

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables
REAL PrL    !Liquid Prandtl number
REAL ReL    !Liquid Reynold number
REAL EE     !Intermediate variable, E in Cavallini et al.
REAL FF     !Intermediate variable, F in Cavallini et al.
REAL HH     !Intermediate variable, H in Cavallini et al.
REAL DD		!Intermediate variable, theta+ in Cavallini et al.
REAL TT	    !Intermediate variable, T+ in Cavallini et al.
REAL fliq   !Liquid friction factor 
REAL fvap   !Vapor friction factor
REAL We     !Weber number
REAL phi2Liq !Square of two-phase multiplier 
REAL DPDZliq !Liquid pressure gradient
REAL tau    !Intermediate variable
REAL mu     !Dynamic viscosity, [Pa-s]

!Flow:

  ReL=Gref*(1-xRef)*ID/muf
  PrL=CpL*muf/kL
  mu=xRef*mug+(1-xRef)*muf

  IF (ReL .LE. 1145) THEN
    DD=(ReL/2)**0.5
  ELSE
    DD=0.0504*ReL**(0.875)
  END IF

  IF (DD .LE. 5) THEN
    TT=DD*PrL
  ELSEIF (DD .GT. 5 .AND. DD .LT. 30) THEN
    TT=5*(PrL+LOG(1+PrL*(DD/5-1)))
  ELSEIF (DD .GE. 30) THEN
    TT=5*(PrL+LOG(1+5*PrL)+0.495*LOG(DD/30))
  END IF
    
  IF (Gref*ID/mu .GT. 2000) THEN
    fliq=0.046*(Gref*ID/muf)**(-0.2)
	fvap=0.046*(Gref*ID/mug)**(-0.2)
  ELSE
	fliq=16/(Gref*ID/muf)
	fvap=16/(Gref*ID/mug)
  END IF

  EE=(1-xRef)**2+xRef**2*(rhof*fvap)/(rhog*fliq)
  FF=xRef**0.6978
  HH=(rhof/rhog)**0.3278*(mug/muf)**(-1.181)*(1-mug/muf)**3.477
  We=Gref**2*ID/(rhog*sigma)

  phi2Liq=EE+(1.262*FF*HH)/We**0.1458

  DPDZliq=phi2Liq*2*fliq*Gref**2/(ID*rhof)

  tau=DPDZliq*ID/4

  hTP=rhof*CpL*(tau/rhof)**0.5/TT
	      
  RETURN

END SUBROUTINE hTPCavalliniAnnular

!************************************************************************

SUBROUTINE hTPCavalliniAnn_Strat(ID,xRef,Gref,rhog,rhof,mug,muf, &
                                 kL,CpL,sigma,hfg,DT,JG,hTP)
!(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                                 kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in annular stratified
!flow condensing region
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!Sept 2006
!
!Reference:
!Cavallini A, Del Col D, Doretti L, Longo G.A., Rossetto, L. "Intube
!condensation of halogenated refrigerants." ASHRAE transactions, 
!2002,108(1), pp.146-161
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL sigma !Surface tension, [N/m]  
REAL hfg   !Latent heat of condensation, [kJ/kg]
REAL DT    !Temperature difference between tube wall and saturate vapor, [C]
REAL JG    !Dimensionless vapor velocity
 
!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables
REAL Gref25  !Refrigerant mass flux at JG=2.5, [kg/s-m^2]
REAL hTPan25 !Annular heat transfer coefficient at JG=2.5, [kW/m^2-K]
REAL hLiq    !Liquid heat transfer coefficient, [kW/m^2-K]
REAL hTPst   !Stratifed heat transfer coefficient, [kW/m^2-K]
REAL alpha   !Zivi's void fraction
REAL theta   !Liquid phase angle

!Flow:

  Gref25=2.5*(9.8*ID*rhog*(rhof-rhog))**0.5/xRef
  
  CALL hTPCavalliniAnnular(ID,xRef,Gref25,rhog,rhof,mug,muf,&
                           kL,CpL,sigma,hTPan25)
          !hTPCavalliniAnnular(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
        !                       kL,kV,CpL,CpV,sigma,hTP)

  hliq=0.023*(kL/ID)*(Gref*ID/muf)**0.8*(CpL*muf/kL)**0.4*(1-xRef)**0.8

  alpha=xRef/(xRef+(1-xRef)*(rhog/rhof)**0.66)

  theta=PI-ACOS(2*alpha-1)

  hTPst=0.725*(1+0.82*((1-xRef)/xRef)**0.268)**(-1)* &
        (kL**3*rhof*(rhof-rhog)*9.8*hfg/(muf*ID*DT))**0.25+ &
  	    hliq*(1-theta/PI) 

  hTP=(hTPan25-hTPst)*(JG/2.5)+hTPst
	      
  RETURN

END SUBROUTINE hTPCavalliniAnn_Strat

!************************************************************************

SUBROUTINE hTPCavalliniStrat_Slug(ID,xRef,Gref,rhog,rhof,mug,muf, &
                                  kL,CpL,sigma,hfg,DT,JG,hTP)
!(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                                  kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in stratified slug
!flow condensing region
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!Sept 2006
!
!Reference:
!Cavallini A, Del Col D, Doretti L, Longo G.A., Rossetto, L. "Intube
!condensation of halogenated refrigerants." ASHRAE transactions, 
!2002,108(1), pp.146-161
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
!INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL sigma !Surface tension, [N/m]  
REAL hfg   !Latent heat of condensation, [kJ/kg]
REAL DT    !Temperature difference between tube wall and saturate vapor, [C]
REAL JG    !Dimensionless vapor velocity

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local variables
REAL x16    !Reference quality at Xtt=1.6
REAL hTP16  !Heat transfer coefficient at x16, [kW/m^2-K]
REAL hLiq   !Liquid heat transfer coefficient, [kW/m^2-K]

!Flow:

  x16=(muf/mug)**0.111*(rhog/rhof)**0.556/(1.686+(muf/mug)**0.111*(rhog/rhof)**0.556)
  
  CALL hTPCavalliniAnn_Strat(ID,x16,Gref,rhog,rhof,mug,muf, &
                             kL,CpL,sigma,hfg,DT,JG,hTP16)
  !(CoilType,ID,xRef,Gref,rhog,rhof,mug,muf, &
!                                 kL,kV,CpL,CpV,sigma,hfg,DT,JG,hTP)
  
  hLiq=0.023*(kL/ID)*(Gref*ID/muf)**0.8*(CpL*muf/kL)**0.4
  
  hTP=hLiq+xRef*(hTP16-hLiq)/x16
	      
  RETURN

END SUBROUTINE hTPCavalliniStrat_Slug

!************************************************************************

SUBROUTINE hTPShahEvap(ID,xRef,mRef,Qout,hfg,vgi,vfi,mug,muf,kL,CpL,Psat,Pcrit,hTP) !(ID,xRef,mRef,Qout,hfg,vgi,vfi,muRef,mug,muf,kL,CpL,Psat,Pcrit,hTP)
		   
!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!evaporating region
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!Shah, M.M (1977). General correlation for heat transfer during 
!subcooled boiling in pipes and annuli. ASHRAE transactions. 83(1), 
!pp. 202-217. 
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID    !Tube inside diameter, [m]
REAL xRef  !Quality of refrigerant, [-]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL Qout  !Tube outside heat transfer, [kW]
REAL hfg   !Enthalpy difference between vapor and fluid states, [kJ/kg]
REAL vgi   !Vapor specific volume, [m^3/kg]
REAL vfi   !Liquid specific volume, [m^3/kg]
REAL mug   !Vapor dynamic viscosity, [Pa-s]
REAL muf   !Liquid dynamic viscosity, [Pa-s]
REAL kL    !Liquid conductivity, [kW/m-K]
REAL CpL   !Liquid specific heat, [kJ/kg-K]
REAL Psat  !Saturation pressure, [kPa]
REAL Pcrit !Critical pressure, [kPa]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local parameter:
INTEGER,PARAMETER :: CoilType=2 !Evaporator

!Subroutine local variables:
REAL hliq    !Liquid heat transfer coefficient, [kW/m^2-K] 
REAL h1      !Intermediate heat transfer coefficient, [kW/m^2-K] 
REAL ZZ      !Intermediate variable
REAL Pred    !Reduced pressure, [kPa]
REAL rhog    !Vapor density, [kg/m^3] 
REAL rhof    !Liquid density, [kg/m^3]
REAL Acs     !Cross sectional area, [m^2]
REAL Gref    !Mass flux, [kg/s-m^2]
REAL Bo      !Boiling number
REAL Co      !Intermediate variable
REAL FrL     !Liquid Froude number
REAL NN      !Intermediate variable
REAL FF      !Intermediate variable
REAL phinb   !Intermediate variable
REAL phicb   !Intermediate variable
REAL phibs   !Intermediate variable
REAL phi     !Intermediate variable

!Flow:

  CALL hSPDittus(CoilType,ID,mRef,0.00,muf,mug,muf,kL,CpL,hliq)
  
  IF (Qout .EQ. 0) THEN
    h1=hliq*(1-xRef)**0.8
    Pred=Psat/Pcrit
    ZZ=(1/xRef-1)**0.8*Pred**0.4
    hTP=h1*(1+3.8/ZZ**0.95)
  ELSE

    rhog=1/vgi
    rhof=1/vfi

    Acs=(ID/2)**2*PI
    Gref=mRef/Acs
    
    Bo=(Qout/Acs)/(Gref*Hfg)
    Co=((1-xRef)/xRef)**0.8*(rhog/rhof)**0.5
    FrL=Gref**2/(rhof**2*9.8*ID)

    IF (FrL .LE. 0.04) THEN
      NN=0.38*FrL**(-0.3)*Co
    ELSE
      NN=Co
    END IF

    IF (Bo .GE. 0.0011) THEN
      FF=14.7
    ELSE
      FF=15.43
    END IF

    IF (NN .GT. 1) THEN
      IF (Bo .GT. 0.00003) THEN
        phinb=230*Bo**0.5
      ELSE
        phinb=1+46*Bo**0.5
      END IF
      phicb=1.8/(NN**0.8)
	ELSEIF (NN .LE. 0.1) THEN
      phibs=FF*Bo**0.5*EXP(2.74*NN**(-0.15))
    ELSE
      phibs=FF*Bo**0.5*EXP(2.74*NN**(-0.1))
    END IF

    IF (phibs .GT. phicb) THEN
      phi=phibs
    ELSE
      phi=phicb
    END IF
    
	IF (phi .LT. phinb) THEN
      phi=phinb
    END IF

    hTP = phi * hLiq

  END IF

  RETURN

END SUBROUTINE hTPShahEvap

!************************************************************************

SUBROUTINE hTPevapWattelet(ID,mRef,muf,mug,vf,vg,kL,kV,CpL,CpV,xRef,Psat,Pcrit,MolWeight,Qout,AoMod,hTP)
!(ID,mRef,muf,mug,vf,vg,kL,kV,CpL,CpV,xRef,Psat,Pcrit,MolWeight,Qout,AoMod,AiMod,hTP)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant heat transfer coefficient in two phase
!evaporating region
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!Wattelet, J.P.; Chato, J.C.; Souza, A.L. and Christoffersen, B.R. (1994). 
!Evaporative charateristics of R-12, R-134a, and a mixture at low
!mass fluxes. ASHRAE transactions, 100(1), pp. 603-615.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID         !Tube inside diameter, [m]
REAL mRef       !Refrigerant mass flow rate, [kg/s]
REAL muf        !Liquid dynamic viscosity, [Pa-s]
REAL mug        !Vapor dynamic viscosity, [Pa-s]
REAL vf         !Liquid specific volume, [m^3/kg]
REAL vg         !Vapor specific volume, [m^3/kg]
REAL kL         !Liquid conductivity, [kW/m-K]
REAL kV         !Vapor conductivity, [kW/m-K]
REAL CpL        !Liquid specific heat, [kJ/kg-K]
REAL CpV        !Vapor specific heat, [kJ/kg-K]
REAL xRef       !Quality of refrigerant, [-]
REAL Psat       !Saturation pressure, [kPa]
REAL Pcrit      !Critical pressure, [kPa]
REAL MolWeight	!Molecular weight of refrigerant, [kg/mol]      
REAL AoMod      !Outside heat transfer area, [m^2]
REAL Qout       !Tube outside heat transfer, [W]

!Output:
REAL hTP   !Two-phase heat transfer coefficient, [kW/m^2-K]

!Subroutine local parameters:
INTEGER,PARAMETER :: CoilType=2 !Evaporator
REAL,PARAMETER :: DryAngle=0 !Dry angle for annular flow
REAL,PARAMETER :: nn=2.5 !3. !Model coefficient
REAL,PARAMETER :: xMax=0.8      !Maximum quality, [-] 
REAL,PARAMETER :: xVap=1.0      !Vapor quality, [-] 

!Subroutine local variables:

REAL Acs      !Cross sectional area, [m^2]
REAL Gref     !Mass flux, [kg/s-m^2]
REAL FrL      !Liquid Froude number 
REAL Pred     !Reduced pressure
REAL RR		  !Intermediate variable
REAL Xtt	  !Lockhart-Martinelli parameter
REAL FF		  !Intermediate variable
REAL hcb	  !Heat transfer coefficient for convective boiling, [kW/m^2-K]
REAL hnb      !Heat transfer coefficient for nucleate boiling, [kW/m^2-K] 
REAL hL		  !Liquid heat transfer coefficient, [kW/m^2-K] 
REAL rhof	  !Liquid density, [kg/m^3]
REAL rhog	  !Vapor density, [kg/m^3]
REAL ReV	  !Vapor Reynold number
REAL PrV	  !Vapor Prandtl number
REAL hMax	  !Heat transfer coefficient at xMax, [kW/m^2-K] 
REAL hvap	  !Vapor heat transfer coefficient, [kW/m^2-K] 
REAL hwet	  !Wet surface heat transfer coefficient, [kW/m^2-K] 	

!Flow:

	rhof=1/vf
	rhog=1/vg

	CALL hSPDittus(CoilType,ID,mRef,0.00,muf,mug,muf,kL,CpL,hL)

	Acs=(ID/2)**2*PI
	Gref=mRef/Acs

	FrL=Gref**2/(rhof**2*9.8*ID)

	IF (FrL .LT. 0.25) THEN
		RR=1.32*FrL**0.2
	ELSE
		RR=1.0
	END IF

	Pred=Psat/Pcrit

	hnb=55.0*MolWeight**(-0.5)*(Qout/AoMod)**0.67*Pred**0.12*(-LOG10(Pred))**(-0.55)

    IF (xRef .LT. xMax) THEN
		Xtt=(rhog/rhof)**0.5*(muf/mug)**0.125*((1-xRef)/xRef)**0.875

		FF=1+1.925*Xtt**(-0.83)

		hcb=FF*hL*RR

		ReV=Gref*xRef*ID/mug

		PrV=CpV*mug/kV

		hvap=0.023*ReV**0.8*PrV**0.4*kV/ID

		hwet=(hnb**nn+hcb**nn)**(1/nn)

		htp=(hwet*(2*PI-DryAngle)+DryAngle*hvap)/(2*PI)
		htp=htp*1e-3
	ELSE
		Xtt=(rhog/rhof)**0.5*(muf/mug)**0.125*((1-xMax)/xMax)**0.875

		FF=1+1.925*Xtt**(-0.83)

		hcb=FF*hL*RR

		ReV=Gref*xMax*ID/mug

		PrV=CpV*mug/kV

		hvap=0.023*ReV**0.8*PrV**0.4*kV/ID

		hwet=(hnb**nn+hcb**nn)**(1/nn)

		hMax=(hwet*(2*PI-DryAngle)+DryAngle*hvap)/(2*PI)
		hMax=hMax*1e-3

		CALL hSPDittus(CoilType,ID,mRef,xVap,mug,mug,muf,kV*1e-3,CpV*1e-3,hVap)

		htp=(xRef-xMax)/(xVap-xMax)*(hVap-hMax)+hMax

	END IF
	
	RETURN

END SUBROUTINE hTPevapWattelet

!************************************************************************

SUBROUTINE hTPconst(Const,hTP)

!------------------------------------------------------------------------
!Purpose:
!To set the refrigerant heat transfer coefficient to constant in two phase
!region
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

!Subroutine passing variables:
!Input:
REAL Const !A constant, [kW/m^2-K]

!Output:
REAL hTP   !Heat transfer coefficient, [kW/m^2-K]

!Flow:

  hTP=Const

  RETURN

END SUBROUTINE hTPconst

!*************************************************************************************

SUBROUTINE MODdP(CoilType,TubeType,hg,hf,hRi,hRo,xRi,xRo, &
                 vRi,vRo,vgi,vfi,vgo,vfo,mRef,muRef,mug,muf, &
				 Lmod,LmodTPratio,ID,HtCoil,Lcoil,dPfric,dPmom,dPgrav)
!(CoilType,TubeType,tRi,tRo,pRi,hg,hf,hRi,hRo,xRi,xRo, &
!                 vRi,vRo,vgi,vfi,vgo,vfo,mRef,muRef,mug,muf,Sigma, &
!				 Lmod,LmodTPratio,ID,OD,HtCoil,Lcoil,dPfric,dPmom,dPgrav)

!------------------------------------------------------------------------
!Purpose:
!To calculate pressure drop in a module
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Eckels, S.J.; Doerr, T.M.; Pate, M.B. In-tube heat transfer and pressure
!drop of R-134a and ester lubricant mixtures in a smooth tube and a 
!micro-fin tube: part II - condensation. ASHRAE Transactions, v 100, n 2,  
!pp. 283-294, 1994.
!
!Eckels, S.J.; Doerr, T.M.; Pate, M.B. In-tube heat transfer and pressure 
!drop of R-134a and ester lubricant mixtures in a smooth tube and a 
!micro-fin tube: part I - evaporation. ASHRAE Transactions, v 100, n 2, 
!pp. 265-282, 1994.
!
!Outokumpu(2005). Tube side enhancement factors. Internal data.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
REAL hg         !Vapor refrigerant enthalpy, [kJ/kg]
REAL hf         !Liquid refrigerant enthalpy, [kJ/kg]
REAL hRi        !Inlet refrigerant enthalpy, [kJ/kg]
REAL hRo        !Outlet refrigerant enthalpy, [kJ/kg]
REAL xRi        !Inlet refrigerant quality, [-]
REAL xRo        !Outlet refrigerant quality, [-]
REAL vRi        !Inlet refrigerant specific volume, [m^3/kg]
REAL vRo        !Outlet refrigerant specific volume, [m^3/kg]
REAL vgi        !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi        !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo        !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo        !Liquid refrigerant specific volume, [m^3/kg]
REAL mRef       !Refrigerant mass flow rate, [kg/s]
REAL muRef      !Refrigerant dynamic viscosity, [Pa-s]
REAL mug        !Vapor refrigerant dynamic viscosity, [Pa-s] 
REAL muf        !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL Lmod       !Module length, [m]
REAL LmodTPratio !Ratio of two-phase length to the total length
REAL ID         !Tube inside diameter, [m]
REAL HtCoil     !Coil height, [m]
REAL Lcoil      !Coil length, [m]

!Outputs:
REAL dPfric    !Frictional pressure drop, [kPa]
REAL dPmom     !Momentum pressure drop, [kPa] 
REAL dPgrav    !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL dPfricSP  !Single phase frictional pressure drop, [kPa]
REAL dPmomSP   !Single phase momentum pressure drop, [kPa]
REAL dPgravSP  !Single phase gravitational pressure drop, [kPa]
REAL dPfricTP  !Two phase frictional pressure drop, [kPa]
REAL dPmomTP   !Two phase momentum pressure drop, [kPa]
REAL dPgravTP  !Two phase gravitational pressure drop, [kPa]
REAL FracTP    !Two phase fraction, [-]
REAL Gref      !Refrigerant mass flux, [kg/s-m^2]
REAL PF        !Pressure enhancement factor

!Flow:

  PF=1
  dPfricSP=0
  dPfricTP=0
  dPmomSP=0
  dPmomTP=0
  dPgravSP=0
  dPgravTP=0

  IF (xRi .GT. 1) THEN
      xRi=1
  END IF
  IF (xRo .GT. 1) THEN
      xRo=1
  ENDIF
  IF (xRi .LT. 0) THEN
      xRi=0
  ENDIF
  IF (xRo .LT. 0) THEN
      xRo=0
  END IF

  IF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. (xRo .LT. 1 .AND. xRo .GT. 0)) THEN
	IF (CoilType .EQ. MCCONDENSER .OR. &
	    CoilType .EQ. MCEVAPORATOR) THEN
		
		CALL TWOPhasedPSouza(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
							 dPfric,dPmom,dPgrav,mRef,ID,mug,muf,HtCoil,Lcoil) !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
		!CALL TWOPhasedPLM(tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
        !                  dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
		!CALL TwoPhaseDPMoriyama(CoilType,tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
        !                        dPmom,dPgrav,mRef,ID,muRef,mug,muf,sigma,HtCoil,Lcoil)
	ELSE
       !CALL TWOPhasedPChoi(hf,hg,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
							!dPfric,dPmom,dPgrav,mRef,ID,mug,muf, &
							!HtCoil,Lcoil)
                            
		CALL TWOPhasedPSouza(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
							 dPfric,dPmom,dPgrav,mRef,ID,mug,muf,HtCoil,Lcoil) !Modified by ISI 07/09/06 !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

  		!CALL TWOPhasedPLM(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
    !                          dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
	END IF
  ELSEIF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. xRo .GE. 1) THEN !Evaporator outlet

	IF (LmodTPratio .NE. 0) THEN !ISI - 10/29/06
		FracTP=LmodTPratio
	ELSE
		IF (hRi .EQ. hRo) THEN
		  FracTP=xRi/(xRo+xRi)
		ELSE
		  FracTP=ABS((hg-hRi)/(hRi-hRo))
		END IF
	END IF	

	IF (FracTP .EQ. 0) THEN
		dPfricTP=0
		dPmomTP=0
		dPgravTP=0
	ELSE

		IF (CoilType .EQ. MCCONDENSER .OR. &
		    CoilType .EQ. MCEVAPORATOR) THEN

			CALL TWOPhasedPSouza(xRi,0.9999,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !
			!CALL TwoPhaseDPMoriyama(CoilType,tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
			!						dPmom,dPgrav,mRef,ID,muRef,mug,muf,sigma,HtCoil,Lcoil)

			!CALL TWOPhasedPLM(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
            !                      dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		ELSE

			!CALL TWOPhasedPChoi(hf,hg,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
			!				dPfric,dPmom,dPgrav,mRef,ID,mug,muf, &
			!				HtCoil,Lcoil)
            
			CALL TWOPhasedPSouza(xRi,0.9999,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !Modified by ISI 07/09/06 !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
		
  		    !CALL TWOPhasedPLM(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
        !                      dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		END IF

	END IF

    IF (FracTP .NE. 1) THEN
		CALL ONEPhasedP(CoilType,1.00,vRi,vRo,vgi,vfi,Lmod*(1-FracTP), & !CoilType,tRi,pRi,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, &
						dPfricSP,dPmomSP,dPgravSP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
    END IF
	
	dPfric=dPfricSP+dPfricTP
	dPmom=dPmomSP+dPmomTP
	dPgrav=dPgravSP+dPgravTP

  ELSEIF ((xRo .LT. 1 .AND. xRo .GT. 0) .AND. xRi .GE. 1) THEN !Condenser inlet

	IF (LmodTPratio .NE. 0) THEN !ISI - 10/29/06
		FracTP=LmodTPratio
	ELSE
		IF (hRi .EQ. hRo) THEN
			FracTP=xRo/(xRo+xRi)
		ELSE
			FracTP=ABS((hg-hRo)/(hRi-hRo))
		END IF
	END IF	


	IF (FracTP .EQ. 0) THEN
		dPfricTP=0
		dPmomTP=0
		dPgravTP=0
	ELSE

		IF (CoilType .EQ. MCCONDENSER .OR. &
		    CoilType .EQ. MCEVAPORATOR) THEN

			CALL TWOPhasedPSouza(0.9999,xRo,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
			!CALL TwoPhaseDPMoriyama(CoilType,tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
			!						dPmom,dPgrav,mRef,ID,muRef,mug,muf,sigma,HtCoil,Lcoil)

			!CALL TWOPhasedPLM(tRi,tRo,pRi,0.9999,xRo,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP,dPfricTP, &
		    !                  dPmomTP,dPgravTP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		ELSE

			!CALL TWOPhasedPChoi(TubeType,tRi,tRo,hf,hg,0.9999,xRo,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, &
			!					dPfricTP,dPmomTP,dPgravTP,mRef,ID,OD,muRef,mug,muf, &
			!					HtCoil,Lcoil)

			CALL TWOPhasedPSouza(0.9999,xRo,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !Modified by ISI 07/09/06 !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
			!CALL TWOPhasedPLM(tRi,tRo,pRi,0.9999,xRo,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP,dPfricTP, &
		    !                  dPmomTP,dPgravTP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		END IF
	END IF

    IF (FracTP .NE. 1) THEN
		CALL ONEPhasedP(CoilType,1.00,vRi,vRo,vgi,vfi,Lmod*(1-FracTP), & !CoilType,tRi,pRi,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, &
						dPfricSP,dPmomSP,dPgravSP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
	END IF
	
	dPfric=dPfricSP+dPfricTP
	dPmom=dPmomSP+dPmomTP
	dPgrav=dPgravSP+dPgravTP

  ELSEIF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. xRo .LE. 0) THEN !Condenser outlet

	IF (LmodTPratio .NE. 0) THEN !ISI - 10/29/06
		FracTP=LmodTPratio
	ELSE
		IF (hRi .EQ. hRo) THEN
			 FracTP=xRi/(ABS(xRo)+xRi)
		ELSE
			 FracTP=ABS((hf-hRi)/(hRi-hRo))
		END IF
	END IF

	IF (FracTP .EQ. 0) THEN
		dPfricTP=0
		dPmomTP=0
		dPgravTP=0
	ELSE

		IF (CoilType .EQ. MCCONDENSER .OR. &
		    CoilType .EQ. MCEVAPORATOR) THEN

			CALL TWOPhasedPSouza(xRi,1.0E-6,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
			!CALL TwoPhaseDPMoriyama(CoilType,tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
            !                        dPmom,dPgrav,mRef,ID,muRef,mug,muf,sigma,HtCoil,Lcoil)
			!CALL TWOPhasedPLM(tRi,tRo,pRi,xRi,1.0E-6,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP,dPfricTP, &
		    !                  dPmomTP,dPgravTP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		ELSE

			!CALL TWOPhasedPChoi(TubeType,tRi,tRo,hf,hg,xRi,1.0E-6,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, &
			!					dPfricTP,dPmomTP,dPgravTP,mRef,ID,OD,muRef,mug,muf, &
			!					HtCoil,Lcoil)

			CALL TWOPhasedPSouza(xRi,1.0E-6,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
								 dPfricTP,dPmomTP,dPgravTP,mRef,ID,mug,muf,HtCoil,Lcoil) !Modified by ISI 07/09/06 !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
			!CALL TWOPhasedPLM(tRi,tRo,pRi,xRi,1.0E-6,vRi,vgi,vfi,vgo,vfo,Lmod*FracTP,dPfricTP, &
		    !                  dPmomTP,dPgravTP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

		END IF
	
	END IF

	IF (FracTP .NE. 1) THEN
		CALL ONEPhasedP(CoilType,0.0,vRi,vRo,vgi,vfi,Lmod*(1-FracTP), & !CoilType,tRi,pRi,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, &
						dPfricSP,dPmomSP,dPgravSP,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
	END IF
	
	dPfric=dPfricSP+dPfricTP
	dPmom=dPmomSP+dPmomTP
	dPgrav=dPgravSP+dPgravTP

  ELSE
	CALL ONEPhasedP(CoilType,xRi,vRi,vRo,vgi,vfi,Lmod, & !CoilType,tRi,pRi,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, &
                    dPfric,dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
  END IF

  Gref=mRef/(PI*ID**2/4)
  CALL CalcDPpenaltyFactor(CoilType,TubeType,Gref,PF)
  
  dPfric=dPfric*PF
  dPmom=dPmom*PF
  dPgrav=dPgrav*PF
  
  RETURN

END SUBROUTINE MODdP

!************************************************************************

SUBROUTINE ONEPhasedP(CoilType,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, & !CoilType,tRi,pRi,xRi,vRi,vRo,vgi,vfi,Lmod,dPfric, &
                      dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)
     
!------------------------------------------------------------------------
!Purpose:
!To calculate single phase pressure drop
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------
      
IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
REAL xRi      !Inlet refrigerant quality, [-]
REAL vRi      !Inlet refrigerant specific volume, [m^3/kg]
REAL vRo      !Outlet refrigerant specific volume, [m^3/kg]
REAL vgi      !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi      !Liquid refrigerant specific volume, [m^3/kg] 
REAL Lmod     !Module length, [m]
REAL mRef     !Refrigerant mass flow rate, [kg/s]
REAL ID       !Tube inside diameter, [m]
REAL muRef    !Refrigerant dynamic viscosity, [Pa-s]
REAL mug      !Vapor dynamic viscosity, [Pa-s]
REAL muf      !Liquid dynamic viscosity, [Pa-s]
REAL HtCoil   !Coil height, [m]
REAL Lcoil    !Coil length, [m]

!Outputs:
REAL dPfric   !Frictional pressure drop, [kPa]
REAL dPmom    !Momentum pressure drop, [kPa]
REAL dPgrav   !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL Gref      !Refrigerant mass flux, [kg/s-m^2]
REAL Acs       !Cross sectional area, [m^2]
REAL dPdZfLiq  !Liquid pressure drop gradient, [kPa/m]
REAL dPdZfVap  !Vapor pressure drop gradient, [kPa/m]
REAL dPdZmom   !Momentum pressure drop gradient, [kPa/m]
REAL dPdZgrav  !Gravitational pressure drop gradient, [kPa/m]
REAL dPdZfrict !Fritional pressure drop gradient, [kPa/m]

!Flow:

  CALL FanningdP(CoilType,xRi,vRi,vgi,vfi,dPdZfLiq,dPdZfVap,mRef,ID,muRef,mug,muf)

  IF (xRi .LE. 0.0) THEN
      dPdZfrict=dPdZfLiq
  END IF
  IF (xRi .GE. 1.0) THEN
      dPdZfrict=dPdZfVap
  END IF

  !Momentum change
  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs

  dPdZmom=Gref*Gref*ABS(vRi-vRo)/Lmod

  !Elevation Pchange
  CALL ElevdP(xRi,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)

  !Total Pchange
  dPfric=dPdZfrict*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3

  RETURN

END SUBROUTINE ONEPhasedP

!************************************************************************

SUBROUTINE FanningdP(CoilType,xRef,vRef,vg,vf,dPdZfLiq,dPdZfVap,mRef,ID,muRef,mug,muf) !(CoilType,tR,pR,xRef,vRef,vg,vf,dPdZfLiq,dPdZfVap,mRef,ID,muRef,mug,muf)

!------------------------------------------------------------------------
!Purpose:
!To calculate single-phase friction pressure drop gradient
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
REAL xRef  !Refrigerant quality, [-]
REAL vRef  !Refrigerant specific volume, [m^3/kg]
REAL vg    !Vapor refrigerant specific volume, [m^3/kg]
REAL vf    !Liquid refrigerant specific volume, [m^3/kg]
REAL mRef  !Refrigerant mass flow rate, [kg/s]
REAL ID    !Tube inside diameter, [m]
REAL muRef !Refrigerant dynamic viscosity, [Pa-s]
REAL mug   !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf   !Liquid refrigerant dynamic viscosity, [Pa-s]

!Outputs:
REAL dPdZfLiq !Liquid pressure drop gradient, [kPa/m]
REAL dPdZfVap !Vapor pressure drop gradient, [kPa/m]

!Subroutine local variables:
REAL ReVap    !Vapor Reynold number
REAL ReLiq    !Liquid Reynold number
REAL Acs      !Cross sectional area, [m^2]
REAL Gref     !Refrigerant mass flux, [kg/s-m^2]
REAL fLiq     !Liquid friction factor
REAL fVap     !Vapor friction factor
REAL viLiq    !Liquid specific volume, [m^3/kg]
REAL viVap    !Vapor specific volume, [m^3/kg]

!Flow:

  CALL Reynolds(ID,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)

  CALL FannFact(0.00,ReLiq,CoilType,fLiq)
  CALL FannFact(1.00,ReVap,CoilType,fVap)

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs

  IF (xRef .LE. 0.0) THEN
      viLiq=vRef
  END IF
  IF (xRef .GE. 1.0) THEN
      viVap=vRef
  END IF
  IF (xRef .GT. 0.0 .AND. xRef .LT. 1.0) THEN
    viLiq=vf
    viVap=vg
  END IF

  IF (xRef .LE. 0.0) THEN
    dPdZfVap=0.0
	dPdZfLiq=(2.*fLiq*Gref*Gref*viLiq)/ID
  ELSE IF (xRef .GE. 1.0) THEN
    dPdZfLiq=0.0
	dPdZfVap=(2.*fVap*Gref*Gref*viVap)/ID
  ELSE
    dPdZfVap=((2.*fVap*Gref*Gref*viVap)/ID)*(xRef**2.)
    dPdZfLiq=((2.*fLiq*Gref*Gref*viLiq)/ID)*((1-xRef)**2.)
  END IF

  RETURN

END SUBROUTINE FanningdP

!************************************************************************

SUBROUTINE FannFact(xRef,Re,CoilType,ff)

!------------------------------------------------------------------------
!Purpose:
!To calculate single-phase friction factor
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!Harms, T.M.; Kazmierczak, M.J. and Gerner, F.M. (1999). Developing
!convective heat transfer in deep rectangular microchannels. Int. J. Heat
!and Fluid Flow, 20, 149-157.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing parameters:
!Inputs:
REAL xRef  !Refrigerant quality
REAL Re    !Reynold number
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
!Output:
REAL ff    !Friction factor
REAL ReTransit !Reynolds number transition from laminar to turbulent

!Flow:

  ReTransit=1500 !2300.
  IF (CoilType .EQ. MCCONDENSER .OR. &
      CoilType .EQ. MCEVAPORATOR) THEN
      ReTransit=1500 !Harms et al. 1999
  END IF

  IF (Re .EQ. 0.0) THEN
      ff=0.0
  END IF
  IF (Re .GT. ReTransit) THEN
    IF (xRef .GE. 1.0) THEN
        ff=0.046/(Re**(0.20))
    END IF
    IF (xRef .LE. 0.0) THEN
        ff=0.079/(Re**(0.25))
    END IF
  ELSE
    ff=16./Re
  END IF

  RETURN

END SUBROUTINE FannFact

!************************************************************************

SUBROUTINE TWOPhasedPLM(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, & !tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
                        dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Lockhart-Martinelli
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL ID      !Tube inside diameter, [m]
REAL muRef   !Refrigerant dynamic viscosity, [Pa-s]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local parameters:

!Subroutine local variables:
REAL dPdZfLiq   !Liquid frictional pressure gradient, [kPa/m]
REAL dPdZfVap   !Vapor frictional pressure gradient, [kPa/m]
REAL dPdZfrict  !Frictional pressure gradient, [kPa/m]
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]
REAL PHIg       !Intermediate variable
REAL alphai     !Inlet void fraction
REAL alphao     !Outlet void fraction
REAL Xtt        !Martinelli parameter

!Flow:

  IF (xRo .LE. 0.0) THEN
      xRo=0.00000001
  END IF

  !Two-phase friction pressure change
  CALL FanningdP(CoilType,xRi,vRi,vgi,vfi,dPdZfLiq,dPdZfVap,mRef,ID,muRef,mug,muf)
  Xtt=(dPdZfLiq/dPdZfVap)**0.5
  CALL PHIcalc(Xtt,PHIg)
  dPdZfrict=(PHIg**2.0)*dPdZfVap

  !Chisholm
  !CALL Reynolds(ID,mRef,xRi,muRef,mug,muf,ReVap,ReLiq)
  
  !IF (ReLiq .LT. 2000 .AND. ReVap .LT. 2000) THEN
  !	  CC=5
  !ELSEIF (ReLiq .LT. 2000 .AND. ReVap .GT. 2000) THEN
  !	  CC=10
  !ELSEIF (ReLiq .GT. 2000 .AND. ReVap .LT. 2000) THEN
  !	  CC=12
  !ELSEIF (ReLiq .GT. 2000 .AND. ReVap .GT. 2000) THEN
  !	  CC=20
  !END IF

  !phif2=1+CC/Xtt+1/Xtt**2
  !phig2=1+CC*Xtt+Xtt**2
  !dPdZfrict=MAX(dPdZfLiq*phif2,dPdZfVap*phig2)

  !Two-phase momentum Pchange
  CALL alphaCALC(xRi,vgi,vfi,alphai)
  CALL alphaCALC(xRo,vgo,vfo,alphao)
  CALL PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID)

  !Two-phase elevation Pchange
  CALL ElevdP(xRi,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPdZfrict*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3
  
  IF (xRo .LE. 0.00000001) THEN
      xRo=0.0
  END IF
                  
  RETURN

END SUBROUTINE TWOPhasedPLM

!************************************************************************

SUBROUTINE TWOPhasedPSouza(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
                           dPmom,dPgrav,mRef,ID,mug,muf,HtCoil,Lcoil) !dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop in condensing, and evaporating regions
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Souza, A.L., Chato J.C., Jabardo, J.M.S., Wattelet, J.P., Paneck, L., 
!Christoffersen, B.C. and Rhines, N. Pressure drop during two-phase flow 
!of refrigerants in horizontal smooth tubes. National heat transfer 
!conf.. 1993.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL ID      !Tube inside diameter, [m]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL dPdZfric   !Frictional pressure gradient, [kPa/m]
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]
REAL alphai     !Intermediate variable
REAL alphao     !Intermediate variable
REAL Gref       !Refrigerant mass flux, [kg/s-m^2]
REAL Acs        !Cross sectional area, [m^2]
REAL DPDZliq    !Liquid only pressure graident, Pa/m
REAL rhog       !Vapor density, [kg/m^3]
REAL rhof       !Liquid density, [kg/m^3]
REAL ReLiq      !Liquid Reynold number
REAL fLiq       !Liquid friction factor
REAL FrLiq      !Liquid Froude number
REAL C1         !Intermediate variable
REAL C2         !Intermediate variable
REAL Xtt        !Lockhart Martinelli parameter
REAL phiLiq     !Liquid phase mulitplier
REAL xRef       !Refrigerant quality, [-]

!Flow:

  !Two-phase friction pressure change
  xRef=(xRi+xRo)/2
  Acs=(PI*ID**2)/4
  Gref=mRef/Acs
  rhog=1/vgi
  rhof=1/vfi

  !Souza et al. 1993
  !CALL Reynolds(ID,mRef,0.00,muf,mug,muf,ReVap,ReLiq)
  !DPDZliq=(2.0*0.079*Gref**2)/(rhof*ID*ReLiq**0.25)
  !IF (xRef .GE. 0.05) THEN
  !  Fr=Gref**2/(ID*9.8*rhof**2)
  !  XttIn=((rhog/rhof)**0.5)*((muf/mug)**0.125)*(((1-xRef)/xRef)**0.875)
  !	IF (Fr .LE. 0.7) THEN
  !	  C1=4.172+(5.48*Fr)-(1.564*Fr**2)
  !	  C2=1.773-(0.169*Fr)
  !	ELSE
  !	  C1=7.242
  !	  C2=1.655
  !	END IF
  !	phiIn=(1.376+C1*(XttIn**(-C2)))*(1-xRef)**1.75
  !	dPdZfric=DPDZliq*phiIn
  !ELSE
  !  dPdZfric=DPDZliq
  !END IF

  !Souza et al. 1993
  ReLiq=Gref*ID/muf*(1-xRef)
  fLiq=0.079/ReLiq**0.25
  DPDZliq=(2.0*fLiq*Gref**2*(1-xRef)**2)/(rhof*ID)
  DPDZfric=DPDZliq

  IF (xRef .GE. 0.05) THEN
  	  FrLiq=Gref**2/(rhof**2*9.8*ID)
  	  IF (FrLiq .GT. 0 .AND. FrLiq .LT. 0.7) THEN
  		C1=4.172+5.48*FrLiq-1.564*FrLiq**2
  		C2=1.773-0.169*FrLiq
  	  ELSE
  		C1=7.242
  		C2=1.655
  	  END IF
  	  Xtt=((rhog/rhof)**0.5)*((muf/mug)**0.125)*(((1-xRef)/xRef)**0.875)
  	  phiLiq=(1.376+C1/Xtt**C2)**0.5
  	  DPDZfric=DPDZliq*phiLiq**2
  END IF

  !Souza and Pimenta 1995 - ISI 07/07/06
  !ReLiq=Gref*ID/muf
  !fLiq=0.079/ReLiq**0.25
  !DPDZliq=(2.0*fLiq*Gref**2)/(rhof*ID)
  !gamma=(rhof/rhog)**0.5*(mug/muf)**0.125
  !Xtt=((rhog/rhof)**0.5)*((muf/mug)**0.125)*(((1-xRef)/xRef)**0.875)
  !phiLiq=(1+(gamma**2-1)*xRef**1.75*(1+0.9524*gamma*Xtt**0.4126))**0.5
  !DPDZfric=DPDZliq*phiLiq**2

  !Two-phase momentum Pchange
  CALL alphaCALC(xRi,vgi,vfi,alphai)
  CALL alphaCALC(xRo,vgo,vfo,alphao)
  CALL PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID)
	        
  !Two-phase elevation Pchange
  CALL ElevdP(xRef,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPdZfric*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3
                  
  RETURN

END SUBROUTINE TWOPhasedPSouza

!************************************************************************

SUBROUTINE TWOPhasedPChang(xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, & !tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
                           dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop in condensing, and evaporating regions
!
!Author
!Ipseng Iu
!Johnson Controls, Inc. 
!Wichita, KS
!
!Date
!August 2009
!
!Reference:
!Chang, Y.; Chiang, S.; Chung, T.; Wang, C. (2000). "Two-phase frictional
!characteristics of R-410A and air-water in a 5 mm smooth tube." 
!ASHRAE transactions, 106(1), pp. 792-797.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL ID      !Tube inside diameter, [m]
REAL muRef   !Refrigerant dynamic viscosity, [Pa-s]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL,PARAMETER :: IDReference = 0.01 !reference diameter, m
REAL,PARAMETER :: muReference = 121.33e-6 !reference viscosity, Pa-s at 25 °C, R410A
REAL dPdZfric   !Frictional pressure gradient, [kPa/m]
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]
REAL alphai     !Intermediate variable
REAL alphao     !Intermediate variable
REAL Gref       !Refrigerant mass flux, [kg/s-m^2]
REAL Acs        !Cross sectional area, [m^2]
REAL DPDZliq    !Liquid only pressure graident, Pa/m
REAL rhog       !Vapor density, [kg/m^3]
REAL rhof       !Liquid density, [kg/m^3]
REAL rhoRef     !Refrigerant density, [kg/m^3]
REAL ReVap      !Vapor Reynold number
REAL ReLiq      !Liquid Reynold number
REAL fLiq       !Liquid friction factor
REAL fvap       !Liquid friction factor
REAL Fr         !Froude number
REAL We         !Weber number
REAL phiLiq     !Liquid phase multiplier
REAL xRef       !Refrigerant quality, [-]
REAL A1,A2,A3   !Intermediate variables
REAL omega      !Small diameter correction factor
REAL ReLiqReference   !Liquid Reynolds number at 25 °C, 410A 

!Flow:

  !Two-phase friction pressure change
  xRef=(xRi+xRo)/2
  Acs=(PI*ID**2)/4
  Gref=mRef/Acs
  rhog=1/vgi
  rhof=1/vfi
  
  rhoRef = 1/(xRef/rhog+(1-xRef)/rhof)

  ReLiq=Gref*ID/muf*(1-xRef)
  fLiq=0.079/ReLiq**0.25
  DPDZliq=(2.0*fLiq*Gref**2*(1-xRef)**2)/(rhof*ID)
  DPDZfric=DPDZliq

  IF (xRef .GE. 0.05) THEN
  	  ReVap=Gref*ID/mug*xRef
  	  fvap=0.079/ReVap**0.25
  	    
  	  A1=(1-xRef)**2+xRef**2*(rhof*fvap/(rhog*fliq))
  	  A2=xRef**0.78*(1-xRef)**0.224
  	  A3=(rhof/rhog)**0.91*(mug/muf)**0.19*(1-mug/muf)**0.7
  	  Fr=Gref**2/(9.8*ID*rhoRef)
  	  We=Gref**2*ID/(rhoRef*muRef)
  	  ReLiqReference=Gref*ID/muReference*(1-xRef)
  	  omega=EXP(-(IDreference/ID)**1.7)*LOG(350/(We*(ReLiq/ReLiqReference)**1.3))

  	  phiLiq=(A1+3.24*A2*A3/(Fr**0.045*We**0.035+omega))**0.5
  	  
  	  DPDZfric=DPDZliq*phiLiq**2
  END IF

  !Two-phase momentum Pchange
  CALL alphaCALC(xRi,vgi,vfi,alphai)
  CALL alphaCALC(xRo,vgo,vfo,alphao)
  CALL PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID)
	        
  !Two-phase elevation Pchange
  CALL ElevdP(xRef,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPdZfric*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3
                  
  RETURN

END SUBROUTINE TWOPhasedPChang

!************************************************************************

SUBROUTINE TWOPhasedPChoi(hf,hg,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
                          dPfric,dPmom,dPgrav,mRef,ID,mug,muf, &
!TubeType,tRi,tRo,hf,hg,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
!                          dPfric,dPmom,dPgrav,mRef,ID,OD,muRef,mug,muf, &
						  HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop in condensing, and evaporating regions
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Choi, J.Y., Kedzierski, M. and Domanski, P.A. (1999). "A Generalized pressure
!drop correlation for evaporation and condensation of alternative
!refrigerants in smooth and microfin tubes.", NISTIR 6333. National
!Institute of Standard and Technology, Gaithersburg, MD, USA.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
!INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch
REAL hg      !Vapor refrigerant enthalpy, [kJ/kg]
REAL hf      !Liquid refrigerant enthalpy, [kJ/kg]
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL ID      !Tube inside diameter, [m]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]
REAL Gref       !Refrigerant mass flux, [kg/s-m^2]
REAL Acs        !Cross sectional area, [m^2]
REAL ReVap      !Vapor Reynold number
REAL ReLiq      !Liquid Reynold number
REAL xRef       !Refrigerant quality, [-]
REAL fN !Two-phase friction factor
REAL Kf !Pierre boiling number
REAL vo !Outlet specific volume, [m^3/kg]
REAL vi !Inlet specific volume, [m^3/kg]
REAL hfg !hg-hf, [kJ/kg]
REAL Dh !Hydraulic diameter, [m]

!Flow:

	  Dh=ID
	  Acs=PI*ID**2/4

  !Two-phase friction and momentum pressure change
  hfg=hg-hf
  Gref=mRef/Acs
  vo=xRo*vgo+(1-xRo)*vfo
  vi=xRi*vgi+(1-xRi)*vfi
  Kf=ABS(xRo-xRi)*hfg*1000/(Lmod*9.8)
  CALL Reynolds(Dh,mRef,0.00,muf,mug,muf,ReVap,ReLiq)
  
  !Choi
  fN=0.00506*ReLiq**(-0.0951)*Kf**0.1554
  dPfric=(fN*Lmod*(vo+vi)/Dh+ABS(vo-vi))*Gref**2
  dPdZmom=0

  !Pierre
  !xRef=(xRo+xRi)/2
  !vg=(vgi+vgo)/2
  !fN=0.0185*(Kf/ReLiq)**0.25
  !dPfric=(fN+(xRo-xRi)*ID/(xRef*Lmod))*(Gref**2*xRef*vg*Lmod)/ID
  !dPdZmom=0
     
  !Two-phase elevation Pchange
  CALL ElevdP(xRef,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPfric*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3
                  
  RETURN

END SUBROUTINE TWOPhasedPChoi

!************************************************************************

SUBROUTINE alphaCALC(xR,vg,vf,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction (alpha)
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xR     !Refrigerant quality
REAL vg     !Vapor refrigerant specific volume, [m^3/kg]
REAL vf     !Liquid refrigerant specific volume, [m^3/kg]

!Outputs:
REAL alpha  !Void fraction

!Flow:

  alpha=1.0/(1.0+((1.0-xR)/xR)*((vf/vg)**(2.0/3.0)))

RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE PHIcalc(Xtt,PHIg)

!------------------------------------------------------------------------
!Purpose:
!To calculate the 'PHI' factor used in the Martinelli friction pressure 
!drop calculation.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL Xtt   !Lockhart-Martinelli parameter

!Outputs:
REAL PHIg   !Phi factor

!Flow:

  PHIg=1.0+2.85*(Xtt**0.523)

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE  PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID) !(tRi,tRo,vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase momentum pressure change.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL vgi      !Vapor specific volume, [m^3/kg]
REAL vfi      !Liquid specific volume, [m^3/kg]
REAL vgo      !Vapor specific volume, [m^3/kg]
REAL vfo      !Liquid specific volume, [m^3/kg]
REAL xRi      !Inlet refrigerant quality
REAL xRo      !Outlet refrigerant quality
REAL alphai   !Inlet void fraction
REAL alphao   !Outlet void fraction
REAL Lmod     !Module length, [m]
REAL mRef     !Refrigerant mass flow rate, [kg/s]
REAL ID       !Inside tube diameter, [m]

!Outputs:
REAL dPdZmom  !Momentum pressure gradient, [kPa/m]

!Subroutine local variables
REAL dPdZmomI !Inlet pressure gradient, [kPa/m] 
REAL dPdZmomO !Outlet pressure gradient, [kPa/m]
REAL Gref     !Refrigerant mass flux, [kg/s-m^2]

!Flow:

  Gref=(4.*mRef)/(PI*ID*ID)

  CALL PmomIntegral(xRi,vgi,vfi,alphai,dPdZmomI) !(tR,xR,vg,vf,alpha,dPdZmomInt)
  CALL PmomIntegral(xRo,vgo,vfo,alphao,dPdZmomO)
  
  dPdZmom=ABS(dPdZmomI-dPdZmomO)*(Gref**2)*(1./Lmod)

  !alphai=1.0/(1.0+((1.0-xRi)/xRi)*((vfi/vgi)**(2.0/3.0)))
  !alphao=1.0/(1.0+((1.0-xRo)/xRo)*((vfo/vgo)**(2.0/3.0)))
  
  !Pin=(xRi**2*vgi/alphai+(1-xRi)**2*vfi/(1-alphai))
  !Pout=(xRo**2*vgo/alphao+(1-xRo)**2*vfo/(1-alphao))
  !Pin=Gref**2*(xRi**2*vgi/alphai+(1-xRi)**2*vfi/(1-alphai))
  !Pout=Gref**2*(xRo**2*vgo/alphao+(1-xRo)**2*vfo/(1-alphao))

  !dPdZmom=ABS(Pout-Pin)/Lmod
  IF(alphai .EQ. 1 .AND. alphao .EQ. 1) THEN    !RS: Debugging: Trying to keep dPdZmom from being NaN
    dPdZmom=0.0
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE PmomIntegral(xR,vg,vf,alpha,dPdZmomInt) !(tR,xR,vg,vf,alpha,dPdZmomInt)

!------------------------------------------------------------------------
!Purpose:
!To calculate the integral used to evaluate two-phase momentum pressure
!change at a given quality value.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xR         !Refrigerant quality
REAL vg         !Vapor specific volume, [m^3/kg]
REAL vf         !Liquid specific volume, [m^3/kg]
REAL alpha      !Void fraction

!Outputs:
REAL dPdZmomInt !Pressure gradient

!Flow:

  dPdZmomInt=(xR**2)*vg/alpha+((1.-xR)**2.)*vf/(1.-alpha)

  RETURN

END SUBROUTINE 

!************************************************************************

SUBROUTINE ElevdP(xR,vRef,vg,vf,dPdZgrav,HtCoil,Lcoil) !(xR,tR,vRef,vg,vf,dPdZgrav,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate the the pressure change due to changes in elevation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xR       !Refrigerant quality
REAL vRef     !Refrigerant specific volume, [m^3/kg]
REAL vg       !Vapor specific volume, [m^3/kg]
REAL vf       !Liquid specific volume, [m^3/kg]
REAL HtCoil   !Coil height, [m]
REAL Lcoil    !Coil length, [m]

!Outputs:
REAL dPdZgrav !Gravitational pressure gradient, [kPa/m]

!Subroutine local variables:
REAL vRef1   !Refrigerant specific volume, [m^3/kg]
REAL dZdL    !Pressure gradient, [kPa/m]

!Flow:

  IF (xR .LE. 0.0 .OR. xR .GE. 1.0) THEN
      vRef1=vRef
  END IF
  
  IF (xR .GT. 0.0 .AND. xR .LT. 1.0) THEN
      vRef1=vf+xR*(vg-vf)
  END IF

  dZdL=HtCoil/Lcoil

  dPdZgrav=(1./vRef1)*9.8*dZdL
  
  IF (dZdL .eq. 0) THEN !RS: Debugging: Trying to keep dPdZgrav from being NaN
      dPdZgrav=0
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE Reynolds(ID,mRef,xRef,muRef,muRvap,muRliq,ReVap,ReLiq)

!------------------------------------------------------------------------
!Purpose:
!To calculate the Reynolds number for both the liquid and the vapor, in
!the cases where each one flows alone in the tube
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL ID       !Tube inside diameter, [m]
REAL mRef     !Refrigerant mass flow rate, [kg/s]
REAL xRef     !Refrigerant quality
REAL muRef    !Refrigerant dynamic viscosity, [Pa-s]
REAL muRvap   !Vapor dynamic viscosity, [Pa-s]
REAL muRliq   !Liquid dynamic viscosity, [Pa-s]

!Outputs:
REAL ReVap    !Reynolds number of the vapor refrigerant
REAL ReLiq    !Reynolds number of the liquid refrigerant

!Subroutine local variables:
REAL Acs      !Cross sectional area, [m^2]
REAL Gref     !Refrigerant mass flux, [kg/s-m^2]

!Flow:

  Acs=PI*(ID**2.0)/4.0
  Gref=mRef/Acs
  IF (xRef .GE. 1.) THEN
    ReLiq=0.0
    IF (muRef .NE. 0) THEN
		ReVap=(Gref*ID/muRef)
	ELSE
		ReVap=(Gref*ID/muRvap)
	END IF
  ELSE IF (xRef .LE. 0.) THEN
    IF (muRef .NE. 0) THEN
		ReLiq=(Gref*ID/muRef)
	ELSE
		ReLiq=(Gref*ID/muRliq)
	END IF
    ReVap=0.0
  ELSE
    ReLiq=(Gref*ID/muRliq)*(1.0-xRef)
    ReVap=(Gref*ID/muRvap)*xRef
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE Prandtl(CpRef,muRef,kRef,Pr)

!------------------------------------------------------------------------
!Purpose:
!To calculate the Prandtl number for the refrigerant as a function of 
!temperature and quality
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL muRef   !Refrigerant dynamic viscosity, [Pa-s]
REAL kRef    !Refrigerant thermal conductivity, [kW/m-K]
REAL CpRef   !Refrigerant specific heat, [kJ/kg-K]

!Output:
REAL Pr      !Prandtl number

!Flow:

  Pr=CpRef*muRef/kRef

RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE manifold(CoilType,ID,mRef,xRef,vgi,vfi,mug,muf,dPman) !(CoilType,ID,mRef,xRef,vRef,vgi,vfi,mug,muf,dPman)

!------------------------------------------------------------------------
!Purpose:
!To calculate pressure drop in the condenser manifolds
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Paliwoda, A. Generalized method of pressure drop calculation across 
!pipe components containing two-phase flow of refrigerants. Int. j. of 
!refrigeration. 15(2), pp. 119-125. 1992.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER CoilType	!1-condenser; 
					!2-evaporator;
					!3-High side interconnecting pipes;
					!4-Low side interconnecting pipes
					!5-Microchannel condenser
					!6-Microchannel evaporator
REAL ID     !Inside tube diameter, [m]
REAL mRef   !Refrigerant mass flow rate, [kg/s]
REAL xRef   !Refrigerant quality
REAL vgi    !Vapor specific volume, [m^3/kg]
REAL vfi    !Liquid specific volume, [m^3/kg]
REAL mug    !Vapor dynamic viscosity, [Pa-s] 
REAL muf    !Liquid dynamic viscosity, [Pa-s]

!Outputs:
REAL dPman  !Manifold pressure drop, [kPa]

!Subroutine local parameter:
REAL, PARAMETER :: Kevap=3.15 !2.8-3.5
REAL, PARAMETER :: Kcond=2.85 !2.5-3.2
REAL, PARAMETER :: CC=0.58

!Subroutine local variables:
REAL Acs      !Cross sectional area, [m^2]
REAL Gref     !Refrigerant mass flux, [kg/s-m^2]
REAL rhog     !Vapor density, [kg/m^3]
REAL rhof     !Liquid density, [kg/m^3]
REAL phi      !Intermediate variable
REAL betac    !Two-phase multiplier for pipe component
REAL Kfactor  !K factor

!Flow:

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs
      
  rhog=1./vgi   !Vapor Density
  rhof=1./vfi   !Liquid Density
  
  IF (CoilType .EQ. HIGHSIDETUBE) THEN
	Kfactor=Kcond
  ELSEIF (CoilType .EQ. LOWSIDETUBE) THEN
	Kfactor=Kevap
  END IF

  DPman=Kfactor*Gref**2/(2*rhog)

  IF (xRef .GT. 0.0 .AND. xRef .LT. 1.0) THEN !Two-phase
	phi=rhog/rhof*(muf/mug)**0.25
	betac=(phi+CC*(1-phi)*xRef)*(1-xRef)**0.333+xRef**2.276
	DPman=DPman*betac
  END IF

  DPman=DPman/1000. !RS Comment: Unit Conversion

  RETURN
      
END SUBROUTINE
      
!************************************************************************

SUBROUTINE returnbend(CoilType,TubeType,ID,Pt,mRef,xRef,vRef,vgi,vfi,mug,muf,dPret) !CoilType,TubeType,ID,OD,Pt,mRef,xRef,vRef,vgi,vfi,muRef,mug,muf,dPret)
                          
!------------------------------------------------------------------------
!Purpose:
!To calculate pressure drop in the return bends.
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Paliwoda, A. Generalized method of pressure drop calculation across 
!pipe components containing two-phase flow of refrigerants. Int. j. of 
!refrigeration. 15(2), pp. 119-125. 1992.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER CoilType	!1-condenser; 
					!2-evaporator;
					!3-High side interconnecting pipes;
					!4-Low side interconnecting pipes
					!5-Microchannel condenser
					!6-Microchannel evaporator
INTEGER TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
					!4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
					!7=Helical; 8=42F HXH
REAL ID     !Inside tube diameter, [m]
REAL Pt     !Tube vertical pitch, [m]
REAL mRef   !Refrigerant mass flow rate, [kg/s]
REAL xRef   !Refrigerant quality
REAL vRef   !Refrigerant specific volume, [m^3/kg]
REAL vgi    !Vapor specific volume, [m^3/kg]
REAL vfi    !Liquid specific volume, [m^3/kg]
REAL mug    !Vapor dynamic viscosity, [Pa-s] 
REAL muf    !Liquid dynamic viscosity, [Pa-s]

!Output:
REAL dPret  !Return bend pressure drop, [kPa]

!Subroutine local parameter:
REAL, PARAMETER :: CC=3      !Empirical coefficient

!Subroutine local variables:
REAL Acs     !Cross sectional area, [m^2]
REAL Gref    !Refrigerant mass flux, [kg/s-m^2]
REAL rhog    !Vapor density, [kg/m^3]
REAL rhof    !Liquid density, [kg/m^3]
REAL phi     !Intermediate variable
REAL betac   !Intermediate variable
REAL rho     !Density, [kg/m^3]
REAL PF      !Pressure drop penalty factor
REAL Kfactor !K factor

!Flow:

  Acs=(PI*ID*ID)/4.0
  Gref=mRef/Acs
  rhog=1./vgi   !Vapor Density
  rhof=1./vfi   !Liquid Density
  rho=1./vRef   !Refrigerant Density
  
  !ISI - 09/12/06
  Kfactor=1/(3.426*LOG(Pt/2/ID)+3.8289) !Curve fit from table data

  IF (xRef .GT. 0.0 .AND. xRef .LT. 1.0) THEN !Two-phase
	DPret=Kfactor*Gref**2/(2*rhog)
	phi=rhog/rhof*(muf/mug)**0.25
	betac=(phi+CC*(1-phi)*xRef)*(1-xRef)**0.333+xRef**2.276
	DPret=DPret*betac
  ELSE
	DPret=Kfactor*Gref**2/(2*rho)
  END IF

  DPret=DPret/1000.0    !RS Comment: Unit Conversion

  CALL CalcDPpenaltyFactor(CoilType,TubeType,Gref,PF)
  DPret=DPret*PF

  RETURN
      
END SUBROUTINE

!************************************************************************

SUBROUTINE Inventory(CoilType,ID,mRef,hg,hf,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &  !RS: Debugging: Extraneous hRi & hRo
                    !(CoilType,TubeType,ID,ktube,mRef,Qout,hg,hf,hRi,hRo,xRi,xRo,vRi,vRo,vgi,vfi,vgo,vfo, &
                    !muRef,mug,muf,kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit,Tsat, &
                    !Cair,Const,Rair,Rtube,AiMod,Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)
                     !mug,muf, &
					 Lmod,LmodTP,LmodSP,MassLiq,MassVap,MassMod)

!------------------------------------------------------------------------
!Purpose:
!To calculate the charge in each module and the overall segment 
!refrigerant charge.
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

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
REAL ID			!Inside diameter, [m]
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL hg			!Vapor enthalpy, [kJ/kg]
REAL hf         !Liquid enthalpy, [kJ/kg]
REAL xRi        !Refrigerant inlet quality
REAL xRo        !Refrigerant outlet quality
REAL vRi        !Refrigerant inlet specific volume, [m^3/kg]
REAL vRo        !Refrigerant outlet specific volume, [m^3/kg]
REAL vgi        !Inlet vapor specific volume, [m^3/kg]
REAL vfi        !Inlet liquid specific volume, [m^3/kg]
REAL vgo        !Outlet vapor specific volume, [m^3/kg]
REAL vfo        !Outlet liquid specific volume, [m^3/kg]
REAL Lmod       !Module length, [m]
REAL LmodTP     !Two phase Module length in transition element, [m]
REAL LmodSP     !Single phase Module length in transition element, [m]

!Outputs:
REAL MassLiq    !Liquid mass, [kg]
REAL MassVap    !Vapor mass, [kg]
REAL MassMod    !Mass in module, [kg]

!Subroutine local variables:
REAL hfg        !hg-hf, [kJ/kg]
REAL Acs        !Cross sectional area, [m^2]
REAL Volmod     !Volume of module, [m^3]
REAL FracTP     !Two-phase fraction
REAL Gref       !Refrigerant mass flux, [kg/s-m^2]
REAL rhoRef     !Refrigerant density, [kg/m^3]
REAL MassLiqTP  !Liquid mass in two-phase portion, [kg]
REAL MassVapTP  !Vapor mass in two-phase portion, [kg] 
REAL MassModTP  !Mass in two-phase portion, [kg]
REAL MassModSP  !Mass in single-phase portion, [kg]

!Flow:

  Acs = PI*(ID/2)**2
  Gref=mRef/Acs

  hfg=hg-hf

  IF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. (xRo .LT. 1 .AND. xRo .GT. 0)) THEN

    CALL CalcMassTP(CoilType,ID,mRef,xRi,xRo,vgi,vfi,vgo,vfo, &
					 Lmod, &
                     MassLiq,MassVap,MassMod)
!(CoilType,TubeType,ID,ktube,mRef,Qout,hfg, &
!                      xRi,xRo,vgi,vfi,vgo,vfo,muRef,mug,muf, &
!                      kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit, &
!					  Tsat,Cair,Const,Rair,Rtube,AiMod,Lmod, &
!                      MassLiq,MassVap,MassTot)

  ELSEIF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. xRo .GE. 1) THEN !Evporator outlet

	FracTP=LmodTP/Lmod

    CALL CalcMassTP(CoilType,ID,mRef,xRi,0.9999,vgi,vfi,vgo,vfo, &
                     LmodTP, &
                     MassLiqTP,MassVapTP,MassModTP)
!(CoilType,TubeType,ID,ktube,mRef,Qout,hfg, &
!                      xRi,xRo,vgi,vfi,vgo,vfo,muRef,mug,muf, &
!                      kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit, &
!					  Tsat,Cair,Const,Rair,Rtube,AiMod,Lmod, &
!                      MassLiq,MassVap,MassTot)

    rhoRef=(1/vRo+1/vgo)/2
    Volmod=Acs*Lmod*(1-FracTP)
	MassModSP=Volmod*rhoRef

	MassLiq=MassLiqTP
	MassVap=MassModSP+MassVapTP
	MassMod=MassModSP+MassModTP

  ELSEIF ((xRi .LT. 1 .AND. xRi .GT. 0) .AND. xRo .LE. 0) THEN !Condenser outlet
	FracTP=LmodTP/Lmod

    CALL CalcMassTP(CoilType,ID,mRef,xRi,1.0E-6,vgi,vfi,vgo,vfo, &
                     LmodTP, &
                     MassLiqTP,MassVapTP,MassModTP)
!(CoilType,TubeType,ID,ktube,mRef,Qout,hfg, &
!                      xRi,xRo,vgi,vfi,vgo,vfo,muRef,mug,muf, &
!                      kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit, &
!					  Tsat,Cair,Const,Rair,Rtube,AiMod,Lmod, &
!                      MassLiq,MassVap,MassTot)

    rhoRef=(1/vRo+1/vfo)/2
    Volmod=Acs*Lmod*(1-FracTP)
	MassModSP=Volmod*rhoRef

	MassLiq=MassModSP+MassLiqTP
	MassVap=MassVapTP
	MassMod=MassModSP+MassModTP

  ELSEIF ((xRo .LT. 1 .AND. xRo .GT. 0) .AND. xRi .GE. 1) THEN !Condenser inlet

	FracTP=(Lmod-LmodSP)/Lmod

    CALL CalcMassTP(CoilType,ID,mRef,0.9999,xRo,vgi,vfi,vgo,vfo, &
                     Lmod-LmodSP, &
                     MassLiqTP,MassVapTP,MassModTP)

!(CoilType,TubeType,ID,ktube,mRef,Qout,hfg, &
!                      xRi,xRo,vgi,vfi,vgo,vfo,muRef,mug,muf, &
!                      kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit, &
!					  Tsat,Cair,Const,Rair,Rtube,AiMod,Lmod, &
!                      MassLiq,MassVap,MassTot)

    rhoRef=(1/vRi+1/vgi)/2
    Volmod=Acs*Lmod*(1-FracTP)
	MassModSP=Volmod*rhoRef

	MassLiq=MassLiqTP
	MassVap=MassModSP+MassVapTP
	MassMod=MassModSP+MassModTP

  ELSE
	!Single-phase region
    IF (xri .GE. 1) THEN
        rhoRef=1/vgi
    ELSE
        rhoRef=1/vfi
    END IF
    Volmod=Acs*Lmod
	MassMod=Volmod*rhoRef
    IF (xRi .GE. 1 ) THEN
	  MassLiq=0
	  MassVap=MassMod
	ELSE
	  MassLiq=MassMod
	  MassVap=0
	END IF
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE Homo(xRi,xRo,vgo,vfo,Lmod,ID,MassLiq,MassVap,MassTot)

!------------------------------------------------------------------------
!Purpose:
!To calculate the refrigerant inventory using the
!closed form of the homogeneous solution.
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRi   !Inlet refrigerant quality
REAL xRo   !Outlet refrigerant quality
REAL vgo   !Outlet vapor specific volume, [m^3/kg]
REAL vfo   !Outlet liquid specific volume, [m^3/kg]
REAL Lmod  !Module length, [m]
REAL ID    !Inside tube dimeter, [m]

!Outputs:
REAL MassLiq !Liquid mass, [kg]
REAL MassVap !Vapor mass, [kg]
REAL MassTot !Total mass, [kg]

!Subroutine local variables:
REAL rhog       !Vapor density, [kg/m^3]
REAL rhof       !Liquid density, [kg/m^3]
REAL Acs        !Cross sectional area, [m^2]
REAL Vol        !Volume, [m^3]
REAL c1         !Intermediate variable
REAL c2         !Intermediate variable
REAL xo         !Outlet quality
REAL xi         !Inlet quality
REAL VoidFrac   !Void fraction

!Flow:

  rhog=1.0/vgo
  rhof=1.0/vfo
  Acs=PI*(ID/2.0)**2.0
  Vol=Acs*Lmod
  c1=(rhog/rhof)**0.667
  c2=1.0-c1
  xo=xRo
  xi=xRi
  IF (xo .GT. 1.) THEN
      xo=1.0
  END IF
  IF (xi .GT. 1.) THEN
      xi=1.0
  END IF
  IF (xo .LT. 0.) THEN
      xo=0.
  END IF
  IF (xi .LT. 0.) THEN
      xi=0.
  END IF

  VoidFrac=(((xo/c2)-((c1/c2**2.)*LOG(c1+c2*xo)))- &
           ((xi/c2)-((c1/c2**2.)*LOG(c1+c2*xi))))/(xo-xi)
  
  VoidFrac=ABS(VoidFrac)

  MassLiq=Vol*rhof*(1-VoidFrac) !Liquid Mass
  MassVap=Vol*rhog*VoidFrac     !Vapor Mass
  MassTot=MassLiq+MassVap       !Total Mass
  
  RETURN

END SUBROUTINE


!************************************************************************

SUBROUTINE LockMartVoidFrac(xRef,rhog,rhof,Mug,Muf,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction using Lockhart-Martinelli correlation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Lockhart Martinelli
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRef  !Refrigerant quality
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL Mug   !Dynamic visocity of refrigerant vapor, [Pa-s]
REAL Muf   !Dynamic visocity of refrigerant liquid, [Pa-s]

!Outputs:
REAL alpha !Void fraction

!Subroutine local variables:
REAL Xtt   !Lockhart-Martinelli parameter

!Flow:

  Xtt=(((1-xRef)/xRef)**0.9)*((rhog/rhof)**0.5)*((Muf/Mug)**0.1)
  IF (Xtt .LE. 10) THEN
    alpha=(1+Xtt**0.8)**(-0.378)
  ELSE
	alpha=0.823-0.157*LOG(Xtt)
  END IF
  
  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE CalcMassTP(CoilType,ID,mRef, &
                      xRi,xRo,vgi,vfi,vgo,vfo, &
					  Lmod, &
                      MassLiq,MassVap,MassTot)

!(CoilType,TubeType,ID,ktube,mRef,Qout,hfg, &
!                      xRi,xRo,vgi,vfi,vgo,vfo,muRef,mug,muf, &
!                      kRef,kL,kV,CpRef,CpL,CpV,MolWeight,Pref,Psat,Pcrit, &
!					  Tsat,Cair,Const,Rair,Rtube,AiMod,Lmod, &
!                      MassLiq,MassVap,MassTot)

!------------------------------------------------------------------------
!Purpose:
!To calculate the charge inventory in two phase region
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!March 2005
!
!Reference:
!Ragazzi, F. and Pedersen, C.O. (1991). Modular-based computer simulation
!of an air-cooled condenser. ACRC technical report 07.
!University of Illinois, Urbana,IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
!INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch; 6=Turbo-A
REAL ID			!Inside diameter, [m]
REAL mRef		!Refrigerant mass flow rate, [kg/s]
REAL xRi        !Refrigerant inlet quality
REAL xRo        !Refrigerant outlet quality
REAL vgi        !Inlet vapor specific volume, [m^3/kg]
REAL vfi        !Inlet liquid specific volume, [m^3/kg]
REAL vgo        !Outlet vapor specific volume, [m^3/kg]
REAL vfo        !Outlet liquid specific volume, [m^3/kg]
REAL Lmod       !Module length, [m]

!Outputs:
REAL MassLiq    !Liquid mass, [kg]
REAL MassVap    !Vapor mass, [kg]
REAL MassTot    !Total mass, [kg]

!Subroutine local variables:
REAL Acs        !Cross sectional area, [m^2]
REAL Vol        !Volume, [m^3]
REAL Gref       !Mass flux, [kg/s-m^2]
REAL vg         !Vapor specific volume, [m^3/kg]
REAL vf         !Liquid specific volume, [m^3/kg]
REAL rhog       !Vapor density, [kg/m^3]
REAL rhof       !Liquid density, [kg/m^3]
REAL xo         !Outlet quality
REAL xi         !Inlet quality
REAL xRef       !Refrigerant quality
REAL alpha      !Void fraction

!Flow:

  Acs=PI*(ID/2)**2
  Vol=Acs*Lmod
  Gref=mRef/Acs

  CALL CalcMeanProp(vgo,vgi,vg) !Determining the mean vapor specific volume
  CALL CalcMeanProp(vfo,vfi,vf) !Determining the mean liquid specific volume
  rhog=1/vg     !Vapor Density
  rhof=1/vf     !Liquid Density

  xo=xRo
  xi=xRi

  IF (CoilType .EQ. HIGHSIDETUBE .OR. &
      CoilType .EQ. LOWSIDETUBE) THEN !for interconnecting pipes

    xRef=(xi+xo)/2

    !CALL HughmarkVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)
	!CALL LockMartVoidFrac(xRef,rhog,rhof,Mug,Muf,alpha)
    CALL GrahamVoidFrac(xRef,ID,Gref,rhog,alpha) !Modified by ISI 07/09/06 !(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)
    !CALL HarmsVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)

  ELSE !for condenser and evaporator
    xRef=(xi+xo)/2
    
    !CALL HughmarkVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)
    !CALL LockMartVoidFrac(x,rhog,rhof,Mug,Muf,alpha)
	CALL GrahamVoidFrac(xRef,ID,Gref,rhog,alpha) !Modified by ISI 07/09/06 !(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)
	!CALL HarmsVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)
    
  END IF

  alpha=ABS(alpha)

  MassLiq=Vol*rhof*(1-alpha)    !Liquid Mass
  MassVap=Vol*rhog*alpha        !Vapor Mass
  MassTot=MassLiq+MassVap       !Total Mass

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE HughmarkVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction using Hughmark correlation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Hughmark, G.A. (1962). Holdup in gas-liquid flow. Chemical engineering 
!progress. 58 (4), pp. 62-65.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRef  !Refrigerant quality
REAL ID    !Tube inside diameter, [m]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL Mug   !Dynamic visocity of refrigerant vapor, [Pa-s]
REAL Muf   !Dynamic visocity of refrigerant liquid, [Pa-s]

!Outputs:
REAL alpha !Void fraction

!Subroutine local parameters:
INTEGER, PARAMETER :: MaxIter = 50
REAL, PARAMETER :: alphaTol = 1E-6

!Subroutine local variables:
REAL Acs         !Cross sectional area, [m^2]
REAL alphaPrev   !Previous value of alpha
REAL beta        !Intermediate variable
REAL yL          !Intermediate variable
REAL Fr          !Intermediate variable
REAL ReAlpha     !Reynold number
REAL ZZ          !Intermediate variable
REAL KH          !Intermediate variable
REAL Error       !Error in iteration
INTEGER Iter                 !Iteration loop counter

!Flow:

  Acs = PI*(ID/2)**2
  
  alphaPrev=xRef
  DO Iter=1,MaxIter

    beta=1/(1+(1-xRef)/xRef*rhog/rhof)  
    yL=1-beta
    Fr=1/ID*(Gref*xRef/(beta*rhog))
    ReAlpha=ID*Gref/(muf+alphaPrev*(mug-muf))
	ZZ=(ID*Gref/(muf+alphaPrev*(mug-muf)))**0.167* &
	   (1/(9.8*ID)*(Gref*xRef/(rhog*beta*(1-beta)))**2)**0.125
	IF (ZZ .LT. 8) THEN
		KH=0.0017*ZZ**3-0.0393*ZZ**2+0.3258*ZZ-0.1792    
	ELSE
        KH=0.6363*ZZ**0.0887
	END IF
    alpha=KH*beta
	IF (alpha .GT. 1) THEN
        alpha=1
    END IF
  
	Error=(alphaPrev-alpha)/alphaPrev
	IF (ABS(Error) .GT. alphaTol) THEN
	  alphaPrev=alpha
	ELSE
	  EXIT
	END IF

  END DO

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE GrahamVoidFrac(xRef,ID,Gref,rhog,alpha) !(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction using Graham et al. correlation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Graham, D.M.; Kopke, H.R.; Wilson, M.J.; Yashar, D.A.; Chato, J.C.; and
!Newell, T.A. (1998). An investigation of void fraction in the stratified/
!annular/intermittent flow regions in smooth, horizontal tubes. ACRC TR-144,
!Air conditioning and refrigeration center, University of Illinois, 
!Urbana-Champaign, IL.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRef  !Refrigerant quality
REAL ID    !Tube inside diameter, [m]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]

!Outputs:
REAL alpha !Void fraction

!Subroutine local variables:
REAL Ft    !Intermediate variable

!Flow:

  Ft=(xRef**3*Gref**2/(rhog**2*9.8*ID*(1-xRef)))**0.5
  
  IF (Ft .GT. 0.01031) THEN
    alpha=1-EXP(-1-0.3*LOG(Ft)-0.0328*(LOG(Ft))**2)
  ELSE
    alpha=0
  END IF

  RETURN

END SUBROUTINE GrahamVoidFrac

!************************************************************************

SUBROUTINE HarmsVoidFrac(xRef,ID,Gref,rhog,rhof,Mug,Muf,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction using Harms et al. correlation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Harms, T.M.; Li, D.; Groll, E.A. and Braun, J.E. (2003). "A void fraction
!model for annular flow in horizontal tubes." International J. of 
!Heat and Mass Transfer. 46. pp. 4051-4057.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRef  !Refrigerant quality
REAL ID    !Tube inside diameter, [m]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL Mug   !Dynamic visocity of refrigerant vapor, [Pa-s]
REAL Muf   !Dynamic visocity of refrigerant liquid, [Pa-s]

!Outputs:
REAL alpha !Void fraction

!Subroutine local variables:
REAL Xtt   !Lockhart Martinelli parameter
REAL Ref   !Intermediate variable

!Flow:

  Xtt=((1-xRef)/xRef)**0.9*(rhog/rhof)**0.5*(mug/muf)**0.1

  Ref = (1 - xRef) * Gref * ID / muf
  
  alpha = (1 - 10.06 * Ref**(-0.875) * (1.74 + 0.104 * Ref**0.5)**2 * &
          (1.376 + 7.242 / (Xtt**1.655))**(-0.5))**2

  RETURN

END SUBROUTINE HarmsVoidFrac

!************************************************************************

SUBROUTINE RouhanniVoidFrac(xRef,ID,Gref,rhog,rhof,Sigma,alpha)

!------------------------------------------------------------------------
!Purpose:
!To calculate void fraction using Rouhanni's correlation
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!July 2006
!
!Reference:
!Wolverine Engineering Book III, Chapter 17
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
REAL xRef  !Refrigerant quality
REAL ID    !Tube inside diameter, [m]
REAL Gref  !Refrigerant mass flux, [kg/s-m^2]
REAL rhog  !Vapor density, [kg/m^3]
REAL rhof  !Liquid density, [kg/m^3]
REAL Sigma !Surface tension, [N/m]

!Outputs:
REAL alpha !Void fraction

!Subroutine local variables:
REAL,PARAMETER :: c0=0.12

REAL mRef  !Refrigerant mass flow rate, kg/s
REAL Ugu   !Drift flux, m3/s
REAL Acs   !Cross-sectional area, m2

!Flow:

  Acs=PI*ID**2/4
  mRef=Gref*Acs

  Ugu=1.18*(1-xRef)*(9.8*Sigma*(rhof-rhog)/rhof**2)**0.25
  alpha=xRef/rhog*((1+c0*(1-xRef))*(xRef/rhog+(1-xRef)/rhof)+Ugu/mRef)

  RETURN

END SUBROUTINE RouhanniVoidFrac

!************************************************************************

SUBROUTINE CalcJfactor(FinType,WetFlag,FinSpg,FinThk,HXdep,Nl, &
                       Dc,Pt,Pl,Amin,AoCoil,ReDc,jfactor)
!(CoilType,FinType,WetFlag,FinSpg,FinThk,Ltube,HXdep,Nl,Nt,RowNum, &
!                       OD,ID,Dc,Pt,Pl,TubeDepth,Aface,Amin,AoCoil,AbrCoil,ReDc,RePl,muAir, &
!					   cpAir,kAir,mAiCoil,FaceVel,Gmax,PrAir,jfactor)

!------------------------------------------------------------------------
!Purpose:
!To calculate j-factor for different fin types under wet or dry condition
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Wang, C.C., Lee, C.J., Chang, C.T. and Lin, S.P. (1999). Heat transfer
!and friction correlation for compact louvered fin-and-tube heat
!exchangers. International journal of heat and mass transfer. 42, pp.
!1945-1956.
!
!Wang, C.C., Lin, Y.T. and Lee, C.J. (2000). An airside correlation for 
!plain fin-and-tube heat exchanger in wet conditions. International
!journal of heat and mass transfer. 43, pp. 1869-1872.
!
!Wang, C.C., Lin, Y.T. and Lee, C.J. (2000). Heat and momentum transfer
!for compact louvered fin-and-tube heat exchangers in wet conditions.
!International journal of heat and mass transfer. 43, pp. 3443-3452.
!
!Wang, C.C., Hwang, Y.M. and Lin, Y.T. (2002). Empirical correlations
!for heat transfer and flow friction characteristics of herringbone
!wavy fin-and-tube heat exchangers. International journal of 
!refrigeration. 25, pp. 673-680.
!
!Wang, C.C., Lee, W.S., Sheu, W.J. (2001). A comparative study of 
!compact enhanced fin-and-tube heat exchangers. International journal
!of heat and mass transfer, 44, pp. 3565-3573.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
!INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

INTEGER FinType !1-Plain; 2-Wavy; i
INTEGER WetFlag !1=Wet; 0=Dry
REAL FinSpg     !Fin spacing, [m]
REAL FinThk     !Fin thickness, [m]
REAL HXdep      !Heat exchanger depth, [m]
INTEGER Nl      !Number of rows
REAL Dc         !Tube outside diameter including fin collar, [m]
REAL Pt         !Tube spacing, [m]
REAL Pl         !Row spacing, [m]
REAL Amin       !Mimimum free flow area, [m^2]
REAL AoCoil     !Coil outside surface area, [m^2]
REAL ReDc       !Reynold number based on Dc

!Outputs:
REAL jfactor !J-factor

!Subroutine local parameters:
REAL, PARAMETER :: theta=8 !Wavy fin angle, [deg]

!Subroutine local variables:
REAL Dh        !Hydraulic diameter, [m]
REAL FinPitch  !Fin pitch, [fins/m]
REAL P1        !Intermediate variable 
REAL P2        !Intermediate variable 
REAL P3        !Intermediate variable
REAL P4        !Intermediate variable 
REAL P5        !Intermediate variable
REAL P6        !Intermediate variable
REAL beta      !Intermediate variable
REAL j1        !Intermediate variable
REAL j2        !Intermediate variable
REAL j3        !Intermediate variable

!Flow:

  IF (WetFlag .EQ. 0) THEN

    SELECT CASE (FinType)
    !CASE (1) !Dry plain fin
    CASE (PLAINFIN)
      !McQuiston
	  !JP=ReDc**(-0.4)*(4*Pl*Pt*Amin/(PI*Dh*Dc*Aface))**(-0.15)
	  !j4=0.0014+0.2618*JP
	  !j=j4*(1-Nl*1280*RePl**(-1.2))/(1-4*1280*RePl**(-1.2))

      !Gray and Webb (1986)
	  !j4=0.14*ReDc**(-0.328)*(Pt/Pl)**(-0.502)*(FinSpg/Dc)**0.0312
	  !j=j4*0.991*(2.24*ReDc**(-0.092)*(Nl/4.)**(-0.031))**(0.607*(4.-Nl))

	  !Briggs and Young (1963)
	  !j=0.134*ReDc**(-0.319)*(FinSpg/HXdep)**0.2*(FinSpg/FinThk)**0.11

	  !Wang et al. (1996)
      !FinPitch=1/(FinSpg+FinThk)
	  !j=0.394*ReDc**(-0.392)*(FinThk/Dc)**(-0.0449)*Nl**(-0.0897)*((1/FinPitch)/Dc)**(-0.212)

	  !Wang et al. (2000)
	  !Nl=1 - 6
	  !Dc = 6.35 - 12.7 mm
	  !FinPitch = 1.19 - 8.7 mm
	  !Pt = 17.7 - 31.75 mm
      !Pl = 12.4 - 27.5 mm

  	  Dh=4*Amin*HXdep/AoCoil
	  FinPitch=1/(FinSpg+FinThk)
	  IF (Nl .EQ. 1) THEN
	    P1=1.9-0.23*LOG(ReDc)
	    P2=-0.236+0.126*LOG(ReDc)
        jfactor=0.108*ReDc**(-0.29)*(Pt/Pl)**P1*((1/FinPitch)/Dc)**(-1.084)*((1/FinPitch)/Dh)**(-0.786)*((1/FinPitch)/Pt)**P2
	  ELSE !IF (Nl .GE. 2) THEN
	    P3=-0.361-0.042*Nl/LOG(ReDc)+0.1581*LOG(Nl*((1/FinPitch)/Dc)**0.41)
	    P4=-1.224-0.076*(Pl/Dh)**1.42/LOG(ReDc)
	    P5=-0.083+0.058*Nl/LOG(ReDc)
	    P6=-5.735+1.211*LOG(ReDc/Nl)
	    jfactor=0.086*ReDc**P3*Nl**P4*((1/FinPitch)/Dc)**P5*((1/FinPitch)/Dh)**P6*((1/FinPitch)/Pt)**(-0.93)
	  END IF

	  !Rich (1975)
	  !IF (Nl .GE. 4) THEN
	  !  jfactor = 0.281 * RePl ** (-0.3941) !4-row
	  !ELSEIF (Nl .EQ. 3) THEN
	  !  jfactor = 0.5813 * RePl ** (-0.4732) !3-row
	  !ELSEIF (Nl .EQ. 2) THEN
	  !  jfactor = 1.2071 * RePl ** (-0.5465) !2-row
	  !ELSEIF (Nl .EQ. 1) THEN
	  !  jfactor = 1.9015 * RePl ** (-0.5918) !1-row
	  !END IF

    !CASE (2) !Dry wavy fin
    CASE (WAVYFIN)
	  !Nl=1-6
	  !Pl=12-33mm
	  !Pt=21-38mm
	  !Dc=7.66-16.85mm
	  !Pd=0.3-1.8mm

	  beta=PI*Dc**2/(4*Pt*Pl)
	  Dh=4*Amin*HXdep/AoCoil
	  FinPitch=1/(FinSpg+FinThk)
	  IF (ReDc .LT. 1000) THEN
        j1=0.0045-0.491*ReDc**(-0.0316-0.0171*LOG(Nl*TAND(theta))) * &
           (Pl/Pt)**(-0.109*LOG(Nl*TAND(theta)))*(Dc/Dh)**(0.542+0.0471*Nl) * &
	       (FinSpg/Dc)**0.984*(FinSpg/Pt)**(-0.349)
	    j2=-2.72+6.84*TAND(theta)
	    j3=2.66*TAND(theta)
	    jfactor=0.882*ReDc**j1*(Dc/Dh)**j2*(FinSpg/Pt)**j3*(FinSpg/Dc)**(-1.58) * &
	      (TAND(theta))**(-0.2)
      ELSE !ReDc .GE. 1000
        j1=-0.0545-0.0538*TAND(theta)-0.302*Nl**(-0.24)*(FinSpg/Pl)**(-1.3)* &
	       (Pl/Pt)**0.379*(Pl/Dh)**(-1.35)*(TAND(theta))**(-0.256)
        j2=-1.29*(Pl/Pt)**(1.77-9.43*TAND(theta))*(Dc/Dh)**(0.229-1.43*TAND(theta)) * &
	       Nl**(-0.166-1.08*TAND(theta))*(FinSpg/Pt)**(-0.174*LOG(0.5*Nl))
	    jfactor=0.0646*ReDc**j1*(Dc/Dh)**j2*(FinSpg/Pt)**(-1.03)*(Pl/Dc)**0.432 * &
	      (TAND(theta))**(-0.692)*Nl**(-0.737)
      END IF

      !RS: Code below all written by RS for testing
   ! CASE(3)   !3 means a louver
   !   Dh=4*Amin*HXdep/AoCoil    !The other cases take these into account, so maybe this one does too
	  !FinPitch=1/(FinSpg+FinThk)
   !
   !   !RS: Assuming that this is indeed a correlation for a louvered fin        
	  !!Wang et. al (1999), IJ Heat and Mass Transfer
	  !j1=-0.229+0.115*((1/FinPitch)/Dc)**0.6*(Pl/Dh)**0.54*Nl**(-0.284)*LOG(0.5*TAND(theta))
	  !j2=-0.251+0.232*Nl**1.37/(LOG(ReDc)-2.303)
	  !j3=-0.439*(FinSpg/Dh)**0.09*(Pl/Pt)**(-1.75)*Nl**(-0.93)
	  !j4=0.502*(LOG(ReDc)-2.54)
	  !jfactor=0.324*ReDc**j1*((1/FinPitch)/Pl)**j2*TAND(theta)**j3*(Pl/Pt)**j4*Nl**0.428
   !   !RS: End of test code by RS
      
	END SELECT

  ELSE 

      !RS Comment: Need a case for louvers!
    SELECT CASE (FinType)
    CASE (PLAINFIN)
      j1=0.3745-1.554*(FinSpg/Dc)**0.24*(Pl/Pt)**0.12*Nl**(-0.19)
      jfactor=19.36*ReDc**j1*(FinSpg/Dc)**1.352*(Pl/Pt)**0.6795*Nl**(-1.291)
  	END SELECT
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE CalcFricfactor(FinType,WetFlag,FinSpg,FinThk,HXdep,Nl,Dc,Pt,Pl, &
                          AoCoil,ReDc,Amin,Fricfactor)

!(CoilType,FinType,WetFlag,FinSpg,FinThk,HXdep,Nl,Nt,Dc,OD,Pt,Pl, &
!                          Ltube,TubeDepth,AoCoil,AbrCoil,ReDc,RePl,Amin,AfrCoil,mAiCoil, &
!						  Gmax,muAir,rhoIn,rhoOut,Fricfactor)

!------------------------------------------------------------------------
!Purpose:
!To calculate air side friction factor for different fin types
!under wet or dry condition
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Wang, C.C., Lee, C.J., Chang, C.T. and Lin, S.P. (1999). Heat transfer
!and friction correlation for compact louvered fin-and-tube heat
!exchangers. International journal of heat and mass transfer. 42, pp.
!1945-1956.
!
!Wang, C.C., Lin, Y.T. and Lee, C.J. (2000). An airside correlation for 
!plain fin-and-tube heat exchanger in wet conditions. International
!journal of heat and mass transfer. 43, pp. 1869-1872.
!
!Wang, C.C., Lin, Y.T. and Lee, C.J. (2000). Heat and momentum transfer
!for compact louvered fin-and-tube heat exchangers in wet conditions.
!International journal of heat and mass transfer. 43, pp. 3443-3452.
!
!Wang, C.C., Hwang, Y.M. and Lin, Y.T. (2002). Empirical correlations
!for heat transfer and flow friction characteristics of herringbone
!wavy fin-and-tube heat exchangers. International journal of 
!refrigeration. 25, pp. 673-680.
!
!Wang, C.C., Lee, W.S., Sheu, W.J. (2001). A comparative study of 
!compact enhanced fin-and-tube heat exchangers. International journal
!of heat and mass transfer, 44, pp. 3565-3573.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
!INTEGER CoilType !1=Condenser; 2=Evaporator;   !RS: Debugging: Extraneous variables
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

INTEGER FinType !1-Plain; 2-Wavy; 
INTEGER WetFlag  !1=Wet; 0=Dry
REAL FinSpg    !Fin spacing, [m]
REAL FinThk    !Fin thickness, [m]
REAL HXdep     !Heat exchanger depth, [m]
INTEGER Nl     !Number of rows
REAL Dc        !Tube outside diameter including fin collar, [m]
REAL Pt        !Tube spacing, [m]
REAL Pl        !Row spacing, [m]
REAL AoCoil    !Coil outside surface area, [m^2]
REAL ReDc      !Reynold number based on Dc
REAL Amin      !Mimimum free flow area, [m^2]

!Outputs:
REAL Fricfactor !Friction factor

!Subroutine local parameters:
REAL, PARAMETER :: theta=8 !Wavy fin angle, [deg]

!Subroutine local variables:
REAL Dh        !Hydraulic diameter, [m]
REAL FinPitch  !Fin pitch, [fins/m]
REAL F1        !Intermediate variable
REAL F2        !Intermediate variable 
REAL F3        !Intermediate variable
REAL F4        !Intermediate variable
REAL beta      !Intermediate variable

!Flow:

  IF (WetFlag .EQ. 0) THEN

    SELECT CASE (FinType)

	!insert correlations  here!

    CASE (PLAINFIN)
	         !Nl=1 - 6
			 !Dc = 6.35 - 12.7 mm
			 !FinPitch = 1.19 - 8.7 mm
			 !Pt = 17.7 - 31.75 mm
             !Pl = 12.4 - 27.5 mm

	  !McQuiston

	  !Dh=Dc*(AoCoil/AbrCoil)/(1+(Pt-Dc)/FinSpg)
	  !FP=ReDc**(-0.25)*(Dc/Dh)**0.25*((Pt-Dc)/(4*(FinSpg-FinThk)))**(-0.4)*ABS(Pt/Dh-1)**(-0.5)
      !f=0.004904+1.382*FP**2

	  !Want et al. (1996)
      !FinPitch=1/(FinSpg+FinThk)
	  !f=1.039*ReDc**(-0.418)*(FinThk/Dc)**(-0.104)*Nl**(-0.0935)*(FinPitch/Dc)**(-0.197)

	  Dh=4*Amin*HXdep/AoCoil
	  FinPitch=1/(FinSpg+FinThk)
      F1=-0.764+0.739*Pt/Pl+0.177*(1/FinPitch)/Dc-0.00758/Nl
	  F2=-15.689+64.021/LOG(ReDc)
	  F3=1.696-15.695/LOG(ReDc)
	  Fricfactor=0.0267*ReDc**F1*(Pt/Pl)**F2*((1/FinPitch)/Dc)**F3

	CASE (WAVYFIN)
	  beta=PI*Dc**2/(4*Pt*Pl)
	  Dh=2*FinSpg*(1-beta)/((1-beta)*(1/COSD(theta))+2*FinSpg*beta/Dc)
	  Dh=4*Amin*HXdep/AoCoil
	  FinPitch=1/(FinSpg+FinThk)
	  IF (ReDc .LT. 1000) THEN
	    F1=-0.574-0.137*(LOG(ReDc)-5.26)**0.245*(Pt/Dc)**(-0.765)*(Dc/Dh)**(-0.243) &
	  	   *(FinSpg/Dh)**(-0.474)*(TAND(theta))**(-0.217)*Nl**0.035
	  	F2=-3.05*TAND(theta)
	  	F3=-0.192*Nl
	  	F4=-0.646*TAND(theta)
	  	Fricfactor=4.37*ReDc**F1*(FinSpg/Dh)**F2*(Pl/Pt)**F3*(Dc/Dh)**0.2054*Nl**F4
	  ELSE
	    F1=-0.141*(FinSpg/Pl)**0.0512*(TAND(theta))**(-0.472)*(Pl/Pt)**0.35 &
	  	   *(Pt/Dh)**(0.449*TAND(theta))*Nl**(-0.049+0.237*TAND(theta))
	  	F2=-0.562*(LOG(ReDc))**(-0.0923)*Nl**0.013
	  	F3=0.302*ReDc**0.03*(Pt/Dc)**0.026
	  	F4=-0.306+3.63*TAND(theta)
	  	Fricfactor=0.228*ReDc**F1*TAND(theta)**F2*(FinSpg/Pl)**F3*(Pl/Dc)**F4*(Dc/Dh)**0.383*(Pl/Pt)**(-0.247)
      END IF
      
      !RS: Code below written by RS for testing
      !CASE (3)   !RS: Dry Louver Case
      !RS: Need Dry Louvered Fin case correlations and equations
      !Fricfactor=0.1  !Just seeing if the case works; Yep, it does
      !RS: End of test code by RS
      
    END SELECT
  ELSE 
    SELECT CASE (FinType)
    CASE (1) !Wet plain fin

    CASE DEFAULT !Wet louvered fin
        !RS Comment: Louvered fin may be case default, but what does that case actually involve? Looks like nothing currently.
        !RS: Also, this is under the ELSE, so it's totally skipped if the above IF is true.
  	END SELECT
  END IF
  
  !RS: Attempt at error handling (8/1/12)  
  !IF (FINTYPE .NE. PLAINFIN .AND. FINTYPE .NE. WAVYFIN) THEN
  !    CALL ShowFatalError('The fin type selected cannot be handled by the program.')
  !END IF

END SUBROUTINE

!************************************************************************

SUBROUTINE CalcFinEff(CoilType,Kfin,FinThk,FinHeight,Rc,TubeDepth,Pt,Pl,hco,FinEff)

!------------------------------------------------------------------------
!Purpose:
!To calculate fin efficiency
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Schmidt, ThE. Heat transfer calculations for extended surfaces. 
!Refrig. Eng., pp. 351-357. 1949.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER CoilType !1=Condenser; 2=Evaporator; 
                 !3=High side interconnecting pipes; 
				 !4=Low side interconnecting pipes
				 !5=Microchannel condenser
				 !6=Microchannel evaporator

REAL Kfin      !Fin thermal conductivity, [kW/m-K]
REAL FinThk    !Fin thickness, [m]
REAL FinHeight !Fin height, [m]
REAL Rc        !Collar radius, [m]
REAL TubeDepth !Tube depth, [m]
REAL Pt        !Tube spacing, [m]
REAL Pl        !Row spacing, [m]
REAL hco       !Heat transfer coefficent, [kW/m^2-K]

!Outputs:
REAL FinEff

!Subroutine local variables:
REAL mm		  !Parameter for fin efficiency calculation, [-]
REAL LL		  !Parameter for fin efficiency calculation, [-]
REAL phi      !Parameter for fin efficiency calculation, [-]
REAL XM       !Heat exchanger geometry parameter
REAL XL       !Heat exchanger geometry parameter
REAL Req      !Equivalent radius, [m]

REAL KFinFrost
REAL FinFrostThk
REAL KFrost
REAL FrostThk

!Flow:

	IF (CoilType .EQ. MCCONDENSER .OR. &
	    CoilType .EQ. MCEVAPORATOR) THEN !Microchannel coil

		mm=(2*hco/(Kfin*FinThk)*(1+FinThk/TubeDepth))**0.5
		LL=FinHeight/2-FinThk
		FinEff=TANH(mm*LL)/(mm*LL)

	ELSE !plate finned tube coil
		!Dry fin Schimdt
		XM=Pt/2
		IF (XM .GT. Pl .AND. Pl .NE. 0) THEN
            XM = Pl
        END IF
		XL=0.5*((Pt/2)**2+Pl**2)**0.5
		Req=1.27*XM*(XL/XM-0.3)**(0.5)

		!Phi for fin efficiency calc.
		phi=(Req/Rc-1.0)*(1.0+0.35*LOG(Req/Rc))
  !Changed by Sankar
        IF(CoilType .EQ. EVAPORATORCOIL) THEN
         Kfrost=FrostParam%Conductivity
         FrostThk=2*FrostParam%Thickness      
         KfinFrost=(Kfin*FinThk+Kfrost*FrostThk)/(FinThk+FrostThk)
         FinFrostThk=FinThk+FrostThk
		 mm=SQRT(2*hco/(KfinFrost*FinFrostThk))
		 FinEff=TANH(mm*Rc*phi)/(mm*Rc*phi)
		ELSE
		 mm=SQRT(2*hco/(Kfin*FinThk))
		 FinEff=TANH(mm*Rc*phi)/(mm*Rc*phi)
		END IF 

		!Hong and Webb
		!FinEff=TANH(m*Rc*phi)*COS(0.1*m*Rc*phi)/(m*Rc*phi) 
	END IF

	RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE CalcMeanProp(In,Out,Mean)

!------------------------------------------------------------------------
!Purpose:
!To calculate mean property
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

!Subroutine passing variables:
!Inputs:
REAL In   !Inlet property
REAL Out  !Outlet property

!Outputs:
REAL Mean !Mean property

!Flow:

  IF (In .NE. 0 .AND. Out .NE. 0) THEN
    Mean=(In+Out)/2
  ELSEIF (In .GT. 0) THEN
    Mean=In
  ELSEIF (Out .GT. 0) THEN
    Mean=Out
  ELSE
    Mean=0
  END IF

  RETURN

END SUBROUTINE

!************************************************************************

SUBROUTINE GetRefName(RefID,RefName)

!------------------------------------------------------------------------
!Purpose:
!To get refrigerant name from refrigerant ID
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

!Subroutine passing variables:
!Inputs:
INTEGER      RefID    !Refrigerant ID:
					  !1=R22; 2=R410A; 3=R407C; 4=R134A; 5=PROPANE; 
 
!Outputs:
CHARACTER*80 RefName  !Refrigerant name

!Flow:

  SELECT CASE (RefID)
  CASE (1)
    RefName='R22'
  CASE (2)
    RefName='R410A'
  CASE (3)
    RefName='R407C'
  CASE (4)
    RefName='R134A'
  CASE (5)
    RefName='PROPANE'
  CASE (6)
    RefName='R417A'
  CASE (7)
    RefName='R507A'
  CASE DEFAULT
    RefName='R22'
  END SELECT

END SUBROUTINE

!************************************************************************

SUBROUTINE GetRefID(RefName,RefID)

!------------------------------------------------------------------------
!Purpose:
!To get refrigerant ID from refrigerant name
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

!Subroutine passing variables:
!Inputs:
CHARACTER*80 RefName  !Refrigerant name
 
!Outputs:
INTEGER      RefID    !Refrigerant ID:
					  !1=R22; 2=R410A; 3=R407C; 4=R134A; 5=PROPANE; 6=R417A; 7=R507A

!Flow:

  IF (TRIM(RefName) .EQ. "R22") THEN
    RefID=1
  ELSEIF (TRIM(RefName) .EQ. "R410A") THEN
    RefID=2
  ELSEIF (TRIM(RefName) .EQ. "R407C") THEN
    RefID=3  
  ELSEIF (TRIM(RefName) .EQ. "R134A") THEN
    RefID=4  
  ELSEIF (TRIM(RefName) .EQ. "PROPANE") THEN
    RefID=5  
  ELSEIF (TRIM(RefName) .EQ. "R417A") THEN
    RefID=6
  ELSEIF (TRIM(RefName) .EQ. "R507A") THEN
    RefID=7    
  ELSE
    RefID=1
  END IF

END SUBROUTINE

!*********************************************************************************************

SUBROUTINE CalcMaterialWeight(CoilType,Ltube,IDtube,ODtube,TubeHeight,TubeDepth, &
                              Dchannel,Nchannels,Pt,Pl,Nt,Nl,NumOfCkts, &
                              FinThk,FinPitch,WeightAluminum,WeightCopper)

!------------------------------------------------------------------------
!Purpose:
!To calculate the weights of material in the coil
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!July 2005
!
!Reference:
!none
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables
INTEGER,INTENT(IN)  :: CoilType	!1-condenser 
							    !2-evaporator
							    !3-High side interconnecting pipes;
							    !4-Low side interconnecting pipes
							    !5-Microchannel condenser
								!6-Microchannel evaporator

REAL, INTENT(IN)    :: Ltube           !Single tube length, m
REAL, INTENT(IN)    :: IDtube          !Tube inside diameter, m
REAL, INTENT(IN)    :: ODtube          !Tube outside diameter, m
REAL, INTENT(IN)    :: TubeHeight      !Tube height, m
REAL, INTENT(IN)    :: TubeDepth       !Tube depth, m
REAL, INTENT(IN)    :: Dchannel		   !Channel diameter, m
INTEGER, INTENT(IN) :: Nchannels	   !Number of channels
REAL, INTENT(IN)    :: Pt              !Tube spacing, vertical, m
REAL, INTENT(IN)    :: Pl              !Row spacing, horizontal, m
INTEGER, INTENT(IN) :: Nt              !Number of tubes, vertical
INTEGER, INTENT(IN) :: Nl              !Number of rows, horizontal
INTEGER, INTENT(IN) :: NumOfCkts       !Number of circuits
REAL, INTENT(IN)    :: FinThk          !Fin thickness, m
REAL, INTENT(IN)    :: FinPitch        !Fin Pitch, fin/m
REAL, INTENT(OUT)   :: WeightAluminum  !Aluminum weight, kg
REAL, INTENT(OUT)   :: WeightCopper    !Copper weight, kg

!Subroutine local parameters
REAL,PARAMETER :: CopperDensity=8910   !Density of copper, kg/m3
REAL,PARAMETER :: AluminumDensity=2740 !Density of aluminum, kg/m3
REAL,PARAMETER :: PI=3.14159265        !Pi

!Subroutine local variables
REAL FinSpg           !Fin spacing, m
REAL FinHeight		  !Fin height, m
REAL CopperVol        !Copper volume, m3
REAL AluminumVol      !Aluminum volume, m3
REAL SingleTubeVol    !Single tube volume, m3
REAL ReturnBendVol    !Return Bend Volume, m3
REAL NumOfReturnBends !Number of return bends
REAL SingleSheetArea  !Single aluminum sheet area, m2
REAL SingleSheetVol   !Single aluminum sheet vol, m3
REAL TubeArea         !Copper tube area, m2
REAL TubeVol          !Tube volume, m3
REAL FinVol           !Fin volume, m3
REAL CollarVol        !Fin collar volume, m3
REAL Dcollar          !Collar diamter, m
REAL Lcoil			  !Coil length, m

!Flow:

	IF (CoilType .EQ. MCCONDENSER .OR. &
	    CoilType .EQ. MCEVAPORATOR) THEN

		Lcoil=Nt*Nl*Ltube
		FinHeight=Pt-TubeHeight
		TubeVol=((TubeHeight*TubeDepth+PI*TubeHeight**2/4)-PI*Dchannel**2/4*Nchannels)*Lcoil
		FinVol=FinThk*FinHeight*TubeDepth*FinPitch*Lcoil
		WeightAluminum=(TubeVol+FinVol)*AluminumDensity

	ELSE

		FinSpg=1/FinPitch-FinThk

		!Copper weight
		SingleTubeVol=PI*(ODtube**2-IDtube**2)/4*Ltube
		ReturnBendVol=PI*(ODtube**2-IDtube**2)/4*(PI*Pt/2)
		NumOfReturnBends=Nl*Nt-NumOfCkts
		CopperVol=SingleTubeVol*Nl*Nt+ReturnBendVol*NumOfReturnBends
		WeightCopper=CopperVol*CopperDensity

		!Aluminum weight
		SingleSheetArea=Pl*Nt*Pt*Nl
		TubeArea=PI*ODtube**2/4*Nl*Nt
		
		Dcollar=ODtube+FinThk*2
		CollarVol=PI*Dcollar*FinSpg*FinThk*Nl*Nt

		SingleSheetVol=(SingleSheetArea-TubeArea)*FinThk+CollarVol
		AluminumVol=FinPitch*Ltube*SingleSheetVol
		WeightAluminum=AluminumVol*AluminumDensity

	END IF

RETURN

END SUBROUTINE CalcMaterialWeight

!************************************************************************

SUBROUTINE MinimumFreeFlowArea(CoilType,Nl, Nt, Pl, Pt, OD, TubeHeight,FinThk, FinPitch, Ltube, Amin)

!------------------------------------------------------------------------
!Purpose:
!To calculate minimum free flow area, for staggered tubes
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!October 2005
!
!Reference:
!Rohsenow, W.M; Hartnett, J.P.; Ganic, E.N. (1985). Hand book of heat
!transfer fundamentals. 2nd ed. McGraw-Hill, N.Y.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables
INTEGER,INTENT(IN) :: CoilType	!1-condenser; 
					 		    !2-evaporator;
							    !3-High side interconnecting pipes;
							    !4-Low side interconnecting pipes
							    !5-Microchannel condenser

INTEGER,          INTENT(IN) :: Nl              !Number of rows, horizontal
INTEGER,          INTENT(IN) :: Nt              !Number of tubes, vertical
REAL, INTENT(IN) :: Pl              !Row spacing, horizontal, m
REAL, INTENT(IN) :: Pt              !Tube spacing, vertical, m
REAL, INTENT(IN) :: OD              !Tube outside diameter, m
REAL, INTENT(IN) :: TubeHeight      !Tube height, m
REAL, INTENT(IN) :: FinThk          !Fin thickness, m
REAL, INTENT(IN) :: FinPitch        !Fin Pitch, fin/m
REAL, INTENT(IN) :: Ltube           !Single tube length, m
REAL, INTENT(OUT) :: Amin           !Minimum free flow area, m^2

!Subroutine local variables
REAL OuterGap         !Outer free space between tubes, m
REAL InnerGap         !Inner free space between tubes, m
REAL MinGap           !Minimum of inner and outer caps, m
REAL Dcollar          !Collar diamter, m
REAL Rc

REAL FinHeight        !Fin height, m

REAL CoilDepth   !Coil depth, m
REAL CoilWidth   !Coil width, m
REAL CoilHeight   !Coil height, m
REAL s1  !half of outer free space between tubes, m
REAL s2  !Inner free space between tubes, m
REAL Smin !Minimum of s1 and s2, m

!Flow:

	IF (CoilType .EQ. MCCONDENSER .OR. &
	    CoilType .EQ. MCEVAPORATOR) THEN !Microchannel coil

		FinHeight=Pt-TubeHeight
		Amin=Nt*Pt*Ltube-TubeHeight*Ltube*Nt-FinThk*FinHeight*FinPitch*Ltube*Nt

	ELSE !Staggered plate finned tube coil

      IF(CoilType==CONDENSERCOIL) THEN 
       Rc=OD/2+FinThk
      ELSE IF(CoilType==EVAPORATORCOIL) THEN
       Rc=OD/2+FinThk+FrostParam%Thickness
      END IF 
		Dcollar=2*Rc

		OuterGap = Pt - Dcollar
		InnerGap = (((Pt / 2.0) ** 2.0 + Pl **  2.0) **  0.5 - Dcollar) * 2.0

		MinGap = OuterGap
		IF (MinGap .GT. InnerGap) THEN
		  MinGap = InnerGap
		END IF

		IF(CoilType==CONDENSERCOIL) THEN
		 IF (Nl .GT. 1) THEN
			Amin = Nt*MinGap*Ltube-FinPitch*Ltube*FinThk*MinGap*Nt
		 ELSE
			Amin = (Pt-Dcollar)*Ltube*Nt-FinPitch*Ltube*FinThk*(Pt-Dcollar)*Nt
		 END IF
		ELSE IF(CoilType==EVAPORATORCOIL) THEN
		 IF (Nl .GT. 1) THEN
			Amin = Nt*MinGap*Ltube-FinPitch*Ltube*(FinThk+2*FrostParam%Thickness)*MinGap*Nt
		 ELSE
			Amin = (Pt-(Dcollar+FrostParam%Thickness))*Ltube*Nt-FinPitch*Ltube*(FinThk+2*FrostParam%Thickness)*(Pt-(Dcollar+FrostParam%Thickness))*Nt
		 END IF
        END IF
        
		!Weber 1991
		CoilDepth=Nl*Pl
		CoilWidth=Ltube
		CoilHeight=Nt*Pt
		IF(CoilType==CONDENSERCOIL) THEN
		 s1=0.5*(1-FinThk*FinPitch)*(Pt-Dcollar)
		 s2=((Pt/2)**2+Pl**2)**0.5-Dcollar-(Pt-Dcollar)*FinThk*FinPitch
		ELSE IF(CoilType==EVAPORATORCOIL) THEN
		 s1=0.5*(1-(FinThk+2*FrostParam%Thickness)*FinPitch)*(Pt-Dcollar)
		 s2=((Pt/2)**2+Pl**2)**0.5-Dcollar-(Pt-Dcollar)*(FinThk+2*FrostParam%Thickness)*FinPitch
		END IF 
		Smin=MIN(2*s1,2*s2)
		!Amin=(CoilHeight/Pt-1)*Smin*CoilWidth-((1-FinThk*FinPitch)*(Pt-Dcollar))*CoilWidth

	END IF
    IF(CoilType==EVAPORATORCOIL) THEN
        EvapMinArea = Amin
    END IF

RETURN

END SUBROUTINE MinimumFreeFlowArea

!************************************************************************

SUBROUTINE TWOPhasedPNino(CoilType,mRef,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
                          dPfric,dPmom,dPgrav,Dh,muRef,mug,muf,Sigma,HtCoil,Lcoil)

!(CoilType,mRef,tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod, &
!                          dPfric,dPmom,dPgrav,Dh,muRef,mug,muf,Sigma,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop in microchannel tubes
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!May 2006
!
!Reference:
!Nino, V.G.; Jassim, E.W.; Hrnjak, P.S.; and Newell, T.A. (2006).  Flow-
!regime-based model for pressure drop predictions in microchannels.
!HVAC&R Research, 12(1), pp. 17-34.
!
!Garimella, S.; Agarwal, A.; and Killion, J.D. (2005). Condensation
!pressure drop in circular microchannels. Heat transfer engineering, 26(3),
!pp. 28-35.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL Dh		 !Hydraulic diameter, [m]
REAL muRef   !Refrigerant dynamic viscosity, [Pa-s]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL Sigma	 !Surface tension, [N/m]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL dPdZfric   !Friction pressure gradient, [kPa/m]
REAL dPdZfrici   !Friction pressure gradient at inlet, [kPa/m]
REAL dPdZfrico   !Friction pressure gradient at outlet, [kPa/m]
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]
REAL Gref       !Refrigerant mass flux, [kg/s-m^2]
REAL Acs        !Cross sectional area, [m^2]
REAL ReVap      !Vapor Reynold number
REAL ReLiq      !Liquid Reynold number
REAL xRef       !Refrigerant quality, [-]
REAL aa			!Intermediate parameter for intermittent flow boundary
REAL bb			!Intermediate parameter for intermittent flow boundary
REAL xInter		!Intermittent flow threshold quality 
REAL rhog		!Vapor density, [kg/m^3]
REAL rhof		!Liquid density, [kg/m^3]
REAL rho		!Density, [kg/m^3]
REAL KE			!Kinetic energy, [kPa]
REAL DPDZg		!Vapor pressure gradient, [kPa/m]
REAL Xtt		!Lockhart-Martinelli parameter
REAL WeVap		!Vapor Weber number
REAL Xann		!Nino parameter
REAL phig		!Vapor two-phase correction factor
REAL fVap		!Vapor fanning friction factor
REAL alphai     !Inlet void fraction
REAL alphao     !Outlet void fraction
!The function shows discontinuities when plotting X vs Pressure Drop
!To Resolve this, the value of Xhi-Xlow is set to a minimum value
REAL Xhi		!High Quality Value for smoothing
REAL Xlow		!High Quality Value for smoothing

!Flow:

  Acs=PI*Dh**2/4
  Gref=mRef/Acs
  
  !Threshold quality from Garimella et al. (2005)
  aa=69.57+22.6*EXP(0.259*Dh*1000)
  bb=-59.99+176.8*EXP(0.383*Dh*1000)
  xInter=aa/(Gref+bb)

  rhog=1/((vgi+vgo)/2)
  rhof=1/((vfi+vfo)/2)
 ! xRef=(xRi+xRo)/2
  IF (Xri .LT. Xro) then ! Expand range to avoid discontinuity
	Xlow=Xri
	Xhi=Xro
  else 
	Xhi=Xri
	Xlow=Xro
  endif

  If (abs(Xhi-Xlow) .LT. 0.35) then
	Xhi=(Xhi+Xlow)/2.0 + 0.35/2.0
	Xlow=(Xhi+Xlow)/2.0 - 0.35/2.0

	If (Xhi .GT. 1.0) THEN
		Xhi=1.0
		Xlow=0.65
	elseIF (Xlow .LT. 0.0) THEN
		Xlow=0.0
		Xhi=0.35
	endif
  endif

  xref=Xlow
  IF (xRef .LT. xInter) THEN !Intermittent flow
	rho=1/(xRef/rhog+(1-xRef)/rhof)
	KE=Gref**2/(2*rho)
	DPDZfrici = 0.045/Dh*KE    
  ELSE !other flow patterns
	CALL Reynolds(Dh,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
    CALL FannFact(1.00,ReVap,CoilType,fVap)
	DPDZg=2*fVap*Gref**2/(Dh*rhog)
	Xtt=((1-xRef)/xRef)**0.875*(rhog/rhof)**0.5*(muf/mug)**0.125
	WeVap=(xRef*Gref)**2*Dh/(rhog*Sigma)
	Xann=((Xtt+1/(WeVap**1.3))*(rhof/rhog)**0.9)
	phig=SQRT(EXP(-0.046*Xann)+0.22*(EXP(-0.002*Xann)-EXP(-7*Xann)))
	DPDZfrici=phig**2*DPDZg
  END IF

  xref=Xhi
  IF (xRef .LT. xInter) THEN !Intermittent flow
	rho=1/(xRef/rhog+(1-xRef)/rhof)
	KE=Gref**2/(2*rho)
	DPDZfrico = 0.045/Dh*KE    
  ELSE !other flow patterns
	CALL Reynolds(Dh,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
    CALL FannFact(1.00,ReVap,CoilType,fVap)
	DPDZg=2*fVap*Gref**2/(Dh*rhog)
	Xtt=((1-xRef)/xRef)**0.875*(rhog/rhof)**0.5*(muf/mug)**0.125
	WeVap=(xRef*Gref)**2*Dh/(rhog*Sigma)
	Xann=((Xtt+1/(WeVap**1.3))*(rhof/rhog)**0.9)
	phig=SQRT(EXP(-0.046*Xann)+0.22*(EXP(-0.002*Xann)-EXP(-7*Xann)))
	DPDZfrico=phig**2*DPDZg
  END IF

  DPDZfric=(DPDZfrici+DPDZfrico)/2.0

  xRef=(xRi+xRo)/2
  IF (xRi .LT. xInter) THEN !Intermittent flow
	rho=1/(xRi/rhog+(1-xRi)/rhof)
	KE=Gref**2/(2*rho)
	DPDZfrico = 0.045/Dh*KE    
  ELSE !other flow patterns
	CALL Reynolds(Dh,mRef,xRi,muRef,mug,muf,ReVap,ReLiq)
    CALL FannFact(1.00,ReVap,CoilType,fVap)
	DPDZg=2*fVap*Gref**2/(Dh*rhog)
	Xtt=((1-xRi)/xRi)**0.875*(rhog/rhof)**0.5*(muf/mug)**0.125
	WeVap=(xRi*Gref)**2*Dh/(rhog*Sigma)
	Xann=((Xtt+1/(WeVap**1.3))*(rhof/rhog)**0.9)
	phig=SQRT(EXP(-0.046*Xann)+0.22*(EXP(-0.002*Xann)-EXP(-7*Xann)))
	DPDZfrici=phig**2*DPDZg
  END IF

  IF (xRo .LT. xInter) THEN !Intermittent flow
	rho=1/(xRo/rhog+(1-xRo)/rhof)
	KE=Gref**2/(2*rho)
	DPDZfrico = 0.045/Dh*KE    
  ELSE !other flow patterns
	CALL Reynolds(Dh,mRef,xRo,muRef,mug,muf,ReVap,ReLiq)
    CALL FannFact(1.00,ReVap,CoilType,fVap)
	DPDZg=2*fVap*Gref**2/(Dh*rhog)
	Xtt=((1-xRo)/xRo)**0.875*(rhog/rhof)**0.5*(muf/mug)**0.125
	WeVap=(xRef*Gref)**2*Dh/(rhog*Sigma)
	Xann=((Xtt+1/(WeVap**1.3))*(rhof/rhog)**0.9)
	phig=SQRT(EXP(-0.046*Xann)+0.22*(EXP(-0.002*Xann)-EXP(-7*Xann)))
	DPDZfrico=phig**2*DPDZg
  END IF

  DPDZfric=(xRef-xRi)/(xRo-xRi)*(DPDZfrico-DPDZfrici)+DPDZfrici

  !Two-phase momentum Pchange
  CALL alphaCALC(xRi,vgi,vfi,alphai)
  CALL alphaCALC(xRo,vgo,vfo,alphao)
  CALL PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,Dh)

  !Two-phase elevation Pchange
  CALL ElevdP(xRef,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPdZfric*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3

RETURN

END SUBROUTINE TWOPhasedPNino

!************************************************************************

SUBROUTINE TwoPhaseDPMoriyama(CoilType,tRi,tRo,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
                              dPmom,dPgrav,mRef,ID,muRef,mug,muf,HtCoil,Lcoil)

!(CoilType,tRi,tRo,pRi,xRi,xRo,vRi,vgi,vfi,vgo,vfo,Lmod,dPfric, &
!                              dPmom,dPgrav,mRef,ID,muRef,mug,muf,sigma,HtCoil,Lcoil)

!------------------------------------------------------------------------
!Purpose:
!To calculate two-phase pressure drop
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Lee, H.J. and Lee, S.Y. (2001). Pressure drop correlations for two-phase
!flow within horizontal rectangular channels with small heights. Int. J.
!Multiphase Flow, 27, pp. 783-796.
!
!Moriyama, K.; Inoue, A. and Ohira, H. (1992). The thermodynamic charateristics
!of two-phase flow in extremely narrow channels (the frictional pressure drop
!void fraction of adiabatic two-component two-phase flow). Heat transfer
!Japanese Research, 21 (8): 823-837.
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
REAL tRi     !Refrigerant inlet temperature, [C]
REAL tRo     !Refrigerant Outlet temperature, [C]
REAL xRi     !Refrigerant inlet quality, [-]
REAL xRo     !Refrigerant outlet quality, [-]
REAL vRi     !Refrigerant specific volume, [m^3/kg]
REAL vgi     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfi     !Liquid refrigerant specific volume, [m^3/kg]
REAL vgo     !Vapor refrigerant specific volume, [m^3/kg]
REAL vfo     !Liquid refrigerant specific volume, [m^3/kg]
REAL Lmod    !Module length, [m]
REAL mRef    !Refrigerant mass flow rate, [kg/s]
REAL ID      !Tube inside diameter, [m]
REAL muRef   !Refrigerant dynamic viscosity, [Pa-s]
REAL mug     !Vapor refrigerant dynamic viscosity, [Pa-s]
REAL muf     !Liquid refrigerant dynamic viscosity, [Pa-s]
REAL HtCoil  !Coil height, [m]
REAL Lcoil   !Coil length, [m]

!Outputs:
REAL dPfric  !Frictional pressure drop, [kPa]
REAL dPmom   !Momentum pressure drop, [kPa]
REAL dPgrav  !Gravitational pressure drop, [kPa]

!Subroutine local variables:
REAL dPdZfLiq   !Liquid frictional pressure gradient, [kPa/m]
REAL dPdZfVap   !Vapor frictional pressure gradient, [kPa/m]
REAL dPdZfrict  !Frictional pressure gradient, [kPa/m]
REAL dPdZmom    !Momentum pressure gradient, [kPa/m]
REAL dPdZgrav   !Gravitational pressure gradient, [kPa/m]s
REAL phiLiq     !Liquid two-phase multiplier
REAL alphai     !Inlet void fraction
REAL alphao     !Outlet void fraction
REAL XX         !Martinelli parameter
REAL ReVap      !Vapor Reynold number
REAL ReLiq      !Liquid Reynold number
REAL xRef       !Refrigerant quality
REAL rhof       !Liquid density, [kg/m^3]
REAL vg         !Vapor specific volumn, [m^3/kg]
REAL vf         !Liquid specific volumn, [m^3/kg]
REAL tRef       !Refrigerant temperature, [C]
REAL KK         !Intermediate parameter

!Flow:

  xRef=(xRi+xRo)/2
  vg=(vgi+vgo)/2
  vf=(vfi+vfo)/2
  tRef=(tRi+tRo)/2
  rhof=1/vfi

  !Two-phase friction pressure change
  CALL FanningdP(CoilType,xRef,vRi,vg,vf,dPdZfLiq,dPdZfVap,mRef,ID,muRef,mug,muf)
  XX=(dPdZfLiq/dPdZfVap)**0.5
  
  CALL Reynolds(ID,mRef,xRef,muRef,mug,muf,ReVap,ReLiq)
  
  !Chisholm 1967
  !IF (ReLiq .LT. 1000 .AND. ReVap .LT. 1000) THEN !Laminar-Laminar
  !  CC=5.0
  !ELSE IF (ReLiq .GT. ReTransit .AND. ReVap .LT. 1000) THEN !Turbulent-Laminar
  !  CC=10.0
  !ELSE IF (ReLiq .LT. 1000 .AND. ReVap .GT. ReTransit) THEN !Laminar-Turbulent
  !  CC=12.0
  !ELSE IF (ReLiq .GT. ReTransit .AND. ReVap .GT. ReTransit) THEN !Turbulent-Turbulent
  !  CC=20.0
  !END IF

  !Lee 2001
  !IF (ReLiq .LT. ReTransit .AND. ReVap .LT. ReTransit) THEN !Laminar-Laminar
  !  aa=6.833e-8; qq=-1.317; rr=0.719; ss=0.557
  !ELSE IF (ReLiq .GT. ReTransit .AND. ReVap .LT. ReTransit) THEN !Turbulent-Laminar
  !  aa=6.185e-2; qq=0;      rr=0;     ss=0.726
  !ELSE IF (ReLiq .LT. ReTransit .AND. ReVap .GT. ReTransit) THEN !Laminar-Turbulent
  ! aa=3.627;    qq=0;      rr=0;     ss=0.174
  !ELSE IF (ReLiq .GT. ReTransit .AND. ReVap .GT. ReTransit) THEN !Turbulent-Turbulent
  !  aa=0.408;    qq=0;      rr=0;     ss=0.451
  !END IF
  !VelLiq=ReLiq*muf/(rhof*ID)
  !psi=muf*VelLiq/sigma
  !lambda=muf**2/(rhof*sigma*ID)
  !CC=aa*lambda**qq*psi**rr*ReLiq**ss

  !Mishima and Hibiki 1996
  !CC=21*(1-EXP(-0.319*ID)) 
  !phiLiq=(1+CC/XX+1/XX**2)**0.5

  !Moriyama and Inoue 1992
  IF (ReLiq .GT. 1.3) THEN
    KK=0.9*(ReLiq)**0.3
  ELSE
    KK=1.0
  END IF
  phiLiq=(1+KK/XX**2)**0.5

  dPdZfrict=(phiLiq**2.0)*dPdZfLiq

  !Two-phase momentum Pchange
  CALL alphaCALC(xRi,vgi,vfi,alphai)
  CALL alphaCALC(xRo,vgo,vfo,alphao)
  CALL PmomTWOphase(vgi,vfi,vgo,vfo,xRi,xRo,alphai,alphao,Lmod,dPdZmom,mRef,ID)

  !Two-phase elevation Pchange
  CALL ElevdP(xRi,vRi,vgi,vfi,dPdZgrav,HtCoil,Lcoil)
      
  !Total two-phase Pchange
  dPfric=dPdZfrict*Lmod*1E-3
  dPmom=dPdZmom*Lmod*1E-3
  dPgrav=dPdZgrav*Lmod*1E-3
                  
  RETURN

END SUBROUTINE TwoPhaseDPMoriyama

!************************************************************************

SUBROUTINE CalcDPpenaltyFactor(CoilType,TubeType,Gref,PF)

!------------------------------------------------------------------------
!Purpose:
!To calculate pressure drop penalty factor due to enhanced tube
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!June 2006
!
!Reference:
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch; 
							
REAL Gref    !Refrigerant mass flux, [kg/s-m2]

!Outputs:
REAL PF      !Penalty factor

!Flow:

  PF=1
  IF (CoilType .EQ. CONDENSERCOIL) THEN !Condenser
    SELECT CASE (TubeType) !Outokumpu data
	CASE (SMOOTH)
      PF=1
	CASE (MICROFIN) !General micro fin, Eckels et al. (1998) ASHRAE RP-630
	  PF=1
	CASE (HERRINGBONE)
	  IF (Gref .GT. 800) THEN
	    PF=3.3817
	  ELSE
	    PF = -6E-07*Gref**2 + 0.0022*Gref + 2.0057
	  END IF
	  PF=PF*1 !Modified by ISI 07/09/06
	CASE (CROSSHATCH)
	  IF (Gref .GT. 1100) THEN
	    PF=1.3081
	  ELSE
	    PF = 1E-06*Gref**2 - 0.0023*Gref + 2.6281
	  END IF
	  PF=PF*1 !Modified by ISI 07/09/06
	CASE (HERRINGBONEWITHCROSSHATCH)
	  IF (Gref .GT. 800) THEN
	    PF=2.3537
	  ELSE
	    PF = 2E-06*Gref**2 - 0.0022*Gref + 2.8337
	  END IF
	  PF=PF*1 !Modified by ISI 07/09/06
 	CASE DEFAULT
      PF=1
	END SELECT

  ELSEIF (CoilType .EQ. EVAPORATORCOIL) THEN !Evaporator
	SELECT CASE (TubeType) !Outokumpu data
	CASE (SMOOTH)
	  PF=1
	CASE (MICROFIN)
	  PF=1
	CASE (HERRINGBONE)
	  IF (Gref .GT. 350) THEN
	    PF=2.0675
	  ELSE
	    PF=6E-08*Gref**3 - 3E-05*Gref**2 + 0.0034*Gref + 1.98
	  END IF
	  PF=PF*1 !Modified by ISI 07/13/06
	CASE (CROSSHATCH)
	  IF (Gref .GT. 280) THEN ! ISI - 07/13/06
		PF=1.479996
	  ELSE
	    PF = -2E-09*Gref**3 + 7E-06*Gref**2 - 0.0035*Gref + 1.9551
	  END IF
	  PF=PF*1
	CASE (HERRINGBONEWITHCROSSHATCH)
	  IF (Gref .GT. 350) THEN
		PF=2.54705
	  ELSE
	    PF=9E-08*Gref**3 - 5E-05*Gref**2 + 0.0106*Gref + 1.1033
	  END IF
	  PF=PF*1 !Modified by ISI 07/13/06
	CASE DEFAULT
	  PF=1
	END SELECT
  END IF

RETURN

END SUBROUTINE CalcDPpenaltyFactor

!************************************************************************

SUBROUTINE CalcHTCenhancementFactor(CoilType,TubeType,Gref,EF)

!------------------------------------------------------------------------
!Purpose:
!To calculate heat transfer enhancement factor due to enhanced tube
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!June 2006
!
!Reference:
!
!------------------------------------------------------------------------

IMPLICIT NONE

!Subroutine passing variables:
!Inputs:
INTEGER          CoilType	!1-condenser; 
							!2-evaporator;
							!3-High side interconnecting pipes;
							!4-Low side interconnecting pipes
							!5-Microchannel condenser
							!6-Microchannel evaporator
INTEGER          TubeType	!1=Plain; 2=General Micro Fin; 3=Herringbone; 
							!4=Crosshatch; 5=Herringbone w/crosshatch
							
REAL Gref    !Refrigerant mass flux, [kg/s-m2]

!Outputs:
REAL EF      !Enhancement factor

!Flow:

	EF=1
	IF (CoilType .EQ. CONDENSERCOIL) THEN !Condenser

		SELECT CASE (TubeType) 
		CASE (SMOOTH) !Plain
			EF=1
		CASE (MICROFIN) !General micro fin, Eckels et al. (1998) ASHRAE RP-630
			EF=-1E-06*Gref + 2.323
		CASE (HERRINGBONE) !Herringbone, !Outokumpu data
			IF (Gref .GT. 1200) THEN
			  EF = 3.1255
			ELSE
			  EF = 3E-06*Gref**2 - 0.0069*Gref + 7.0855
			END IF
			EF=EF*1 !ISI - 07/13/06
		CASE (CROSSHATCH) !Crosshatch, !Outokumpu data
			IF (Gref .GT. 2100) THEN
			  EF = 1
			ELSE
			  EF = 4E-08*Gref**2 - 0.0008*Gref + 2.5461
			END IF
			EF=EF*1 !2 !Modified by ISI 07/09/06
		CASE (HERRINGBONEWITHCROSSHATCH) !Herringbone w/crosshatch, !Outokumpu data
			IF (Gref .GT. 1000) THEN
			  EF = 3.561
			ELSE
			  EF = 4E-06*Gref**2 - 0.0082*Gref + 7.761
			END IF
			EF=EF*1 !ISI - 07/13/06
		CASE DEFAULT
			EF=1
		END SELECT

	ELSEIF (CoilType .EQ. EVAPORATORCOIL) THEN !Evaporator
		SELECT CASE (TubeType)
		CASE (SMOOTH) !Plain
			EF=1
		CASE (MICROFIN) !General micro fin, Eckels et al. (1994)
			EF=-0.0018*Gref + 2.0332
			EF=EF*1
		CASE (HERRINGBONE)
			IF (Gref .GT. 350) THEN
			  EF=2.936117829
			ELSE
			  EF = 0.3194*Gref**0.3787
			END IF
			EF=EF*1 !ISI - 07/13/06
		CASE (CROSSHATCH)
			IF (Gref .GT. 250) THEN
			  EF=2.027555284
			ELSE
			  EF=6.383*Gref**(-0.2077)
			END IF
			EF=EF*1 !2 !ISI - 07/13/06
		CASE (HERRINGBONEWITHCROSSHATCH)
			IF (Gref .GT. 350) THEN
			  EF=3.09902775
			ELSE
			  EF=0.3185*Gref**0.3884
			END IF
			EF=EF*1 !ISI - 07/13/06
		CASE DEFAULT
			EF=1
		END SELECT
	END IF
	IF (EF .LE. 0) THEN
        EF=1
    END IF

RETURN

END SUBROUTINE CalcHTCenhancementFactor

!************************************************************************

SUBROUTINE UpdateRefMassFlowRate(Iter,Ckt,NumOfCkts,pRiCoil,mRefTot,Nnode,Node)

!------------------------------------------------------------------------
!Purpose:
!To update refrigerant mass flow rate in each circuit
!
!Author
!Ipseng Iu
!Oklahoma State Univerity, Stillwater
!
!Date
!March 2005
!
!Reference:
!Vlach, J. and Singhal, K. (1993). Computer methods for circuit 
!analysis and design. 2nd ed. Van Nostrand Reinhold. New York.
!
!------------------------------------------------------------------------

IMPLICIT NONE

INTEGER,INTENT(IN)                         :: Iter !Iteration number
TYPE (CktInfo),DIMENSION(:),INTENT(INOUT)  :: Ckt !Circuit pointer
INTEGER,INTENT(IN)                         :: NumOfCkts !Number of circuits
REAL,INTENT(IN)                            :: pRiCoil !Inlet pressure, kPa
REAL,INTENT(IN)                            :: mRefTot !Tot refrigerant mass flow rate, kg/s
INTEGER,INTENT(IN)                         :: Nnode !Number of nodes for refrigerant flow rate calculation
TYPE (NodeInfo),DIMENSION(:),INTENT(INOUT) :: Node !Split or joint node

REAL,PARAMETER :: Relax=0.6     !Relaxation

INTEGER I,J,K !Loop counter
REAL DPckt  !Refrigerant pressure drop in circuit, kPa
REAL Gckt   !Mass flux in circuit, kg/m^2
REAL pRiCkt !Circuit inlet pressure, kPa      
REAL pRoCkt !Circuit outlet pressure, kPa
INTEGER :: NumOfMods !Number of tube modules/Segments
REAL :: IDtube !Tube inside diameter, m
REAL,ALLOCATABLE,DIMENSION(:,:) :: Gmat    !matrix stores conductance
REAL,ALLOCATABLE,DIMENSION(:,:) :: GmatInv !Inverse of Gmat
REAL,ALLOCATABLE,DIMENSION(:) :: Mmat      !matrix stores mass flow rate
REAL,ALLOCATABLE,DIMENSION(:) :: Pmat      !matrix stores pressure

!FLOW:

	ALLOCATE(Gmat(Nnode,Nnode))
	ALLOCATE(GmatInv(Nnode,Nnode))
	ALLOCATE(Mmat(Nnode))
	ALLOCATE(Pmat(Nnode))

	IDtube=Ckt(1)%Tube(1)%ID
	NumOfMods=Ckt(1)%Tube(1)%NumOfMods

	!Calc. circuit conductance
	DO I=1 ,NumOfCkts

		!Circuit pressure drop
		IF (Ckt(I)%InSplit .GT. 1 .AND. Ckt(I)%OutJoin .GT. 1) THEN !Split inlet and Joint outlet
			pRiCkt=Ckt(I)%Tube(2)%Seg(1)%pRi
	  		pRoCkt=Ckt(I)%Tube(Ckt(I)%Ntube-1)%Seg(NumOfMods)%pRo 

		ELSE IF (Ckt(I)%InSplit .GT. 1) THEN !Split inlet
			pRiCkt=Ckt(I)%Tube(2)%Seg(1)%pRi
	  		pRoCkt=Ckt(I)%Tube(Ckt(I)%Ntube)%Seg(NumOfMods)%pRo

		ELSE IF (Ckt(I)%OutJoin .GT. 1) THEN !Joint outlet
			pRiCkt=Ckt(I)%Tube(1)%Seg(1)%pRi
      		pRoCkt=Ckt(I)%Tube(Ckt(I)%Ntube-1)%Seg(NumOfMods)%pRo 

		ELSE
			pRiCkt=Ckt(I)%Tube(1)%Seg(1)%pRi
			pRoCkt=Ckt(I)%Tube(Ckt(I)%Ntube)%Seg(NumOfMods)%pRo
		END IF

		DPckt=pRiCkt-pRoCkt

		Gckt=Ckt(I)%mRef/((PI*IDtube**2)/4) !Circuit mass flux
		  
		Ckt(I)%Conduct=(Ckt(I)%mRef)/DPCkt !Circuit conductance
	END DO

		  
	!***Recalc. refrigerant mass flow rates***
	Ckt%mRefPrev=Ckt%mRef !Store previous mdot

	DO I=1, Nnode

		!***Fill Gmat***
		DO J=1, Nnode
			Gmat(I,J)=0
			IF (I .EQ. J) THEN !Diagonal
				IF (Node(I)%Num .GT. 0) THEN !Node inside circuit
					DO K=1, NumOfCkts
						IF (Ckt(K)%TubeSequence(1) .EQ. Node(I)%Num .OR. &
							Ckt(K)%TubeSequence(Ckt(K)%Ntube) .EQ. Node(I)%Num) THEN
							Gmat(I,J)=Gmat(I,J)+Ckt(k)%Conduct
						END IF
					END DO !End NumOfCkts
				ELSE !Node at circuit outlet
					DO K=1, NumOfCkts 
						IF (Ckt(K)%OutSplit .LE. 1 .AND. Ckt(K)%OutJoin .LE. 1) THEN
                            Gmat(I,J)=Gmat(I,J)+Ckt(k)%Conduct
                        END IF
					END DO !End NumOfCkts
				END IF
	  
			ELSE !Off diagonal
				IF (Node(I)%Num .GT. 0) THEN !Node inside circuit
					IF (Node(J)%Num .GT. 0) THEN !Last node is circuit outlet
						DO K=1, NumOfCkts
							IF (Ckt(K)%TubeSequence(1) .EQ. Node(I)%Num .AND. &
								Ckt(K)%TubeSequence(Ckt(K)%Ntube) .EQ. Node(J)%Num) THEN
								Gmat(I,J)=Gmat(I,J)-Ckt(k)%Conduct
							ELSEIF (Ckt(K)%TubeSequence(1) .EQ. Node(J)%Num .AND. &
									Ckt(K)%TubeSequence(Ckt(K)%Ntube) .EQ. Node(I)%Num) THEN
								Gmat(I,J)=Gmat(I,J)-Ckt(k)%Conduct
							END IF
						END DO !End NumOfCkts
					ELSE
						DO K=1, NumOfCkts
							IF (Ckt(K)%TubeSequence(1) .EQ. Node(I)%Num .AND. &
								(Ckt(K)%OutSplit .LE. 1 .AND. Ckt(K)%OutJoin .LE. 1)) THEN
								Gmat(I,J)=Gmat(I,J)-Ckt(k)%Conduct
							ELSEIF (Ckt(K)%TubeSequence(1) .EQ. Node(J)%Num .AND. &
								Ckt(K)%TubeSequence(Ckt(K)%Ntube) .EQ. Node(I)%Num) THEN
								Gmat(I,J)=Gmat(I,J)-Ckt(k)%Conduct
							END IF
						END DO !End NumOfCkts
					END IF
				ELSE !Node at circuit outlet
					DO K=1, NumOfCkts 
						IF (Ckt(K)%OutSplit .LE. 1 .AND. Ckt(K)%OutJoin .LE. 1 .AND. &
						    Ckt(K)%TubeSequence(1) .EQ. Node(J)%Num) THEN
                            Gmat(I,J)=Gmat(I,J)-Ckt(k)%Conduct
                        END IF
					END DO !End NumOfCkts
				END IF
			END IF
		END DO !End Nnode J
    
		!***Fill Mmat***
		Mmat(I)=0
		DO J=1, NumOfCkts
			IF (Node(I)%Num .GT. 0) THEN !Node at circuit outlet
				IF (((Ckt(J)%InSplit .LE. 1) .AND. (Ckt(J)%InJoin .LE. 1)) .AND. &
					 ((Ckt(J)%TubeSequence(Ckt(J)%Ntube) .EQ. Node(I)%Num))) THEN
                    Mmat(I)=Mmat(I)+pRiCoil*Ckt(J)%Conduct
                END IF
				ELSE
					IF (((Ckt(J)%InSplit .LE. 1) .AND. (Ckt(J)%InJoin .LE. 1)) .AND. &
						 ((Ckt(J)%OutSplit .LE. 1) .AND. (Ckt(J)%OutJoin .LE. 1))) THEN
                        Mmat(I)=Mmat(I)+pRiCoil*Ckt(J)%Conduct
                    END IF
			END IF
		END DO
	END DO !End Node I

	Mmat(Nnode)=Mmat(Nnode)-mRefTot !Last node is the outlet of coil

	!Calc inverse of Gmat
	CALL matrix_inverse(Gmat, GmatInv, Nnode)

	!Calc node pressure
	CALL CalcMatrixMultVector(GmatInv, Nnode, Mmat, Pmat)

	!Store node pressure
	DO I=1, Nnode
		Node(I)%Pressure=Pmat(I)
	END DO

	DO I=1, NumOfCkts
		DO J=1, Nnode
			IF ((Ckt(I)%InSplit .LE. 1 .AND. Ckt(I)%InJoin .LE. 1) .AND. &
				(Ckt(I)%OutSplit .LE. 1 .AND. Ckt(I)%OutJoin .LE. 1)) THEN !Connected to coil inlet and outlet 
				
				Ckt(I)%mRef=(pRiCoil-Node(Nnode)%Pressure)*Ckt(I)%Conduct                  
			ELSEIF (Ckt(I)%InSplit .LE. 1 .AND. Ckt(I)%InJoin .LE. 1) THEN !Connected to coil inlet 
				
				IF (Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Node(J)%Num) &
				    Ckt(I)%mRef=(pRiCoil-Node(J)%Pressure)*Ckt(I)%Conduct
				ELSEIF (Ckt(I)%OutSplit .LE. 1 .AND. Ckt(I)%OutJoin .LE. 1) THEN !Connected to coil outlet
					
					IF (Ckt(I)%TubeSequence(1) .EQ. Node(J)%Num) THEN
                        Ckt(I)%mRef=(Node(J)%Pressure-Node(Nnode)%Pressure)*Ckt(I)%Conduct
                    END IF
					ELSE !Connected between circuits
						
						DO K=1, Nnode
							IF (Ckt(I)%TubeSequence(1) .EQ. Node(J)%Num .AND. &
								Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Node(K)%Num) THEN !Pressure J > Pressure K
								
								Ckt(I)%mRef=(Node(J)%Pressure-Node(K)%Pressure)*Ckt(I)%Conduct
							ELSE IF (Ckt(I)%TubeSequence(1) .EQ. Node(K)%Num .AND. &
								Ckt(I)%TubeSequence(Ckt(I)%Ntube) .EQ. Node(J)%Num) THEN !Pressure K > Pressure J
								
								Ckt(I)%mRef=(Node(K)%Pressure-Node(J)%Pressure)*Ckt(I)%Conduct
  							END IF
						END DO !End K nodes
			END IF
		END DO !End J nodes
	END DO !End circuits

	!Apply relaxation
	IF (Iter .GE. 2) THEN
		Ckt%mRef=(Ckt%mRef*(Relax)+Ckt%mRefPrev*(1-Relax))
	END IF

	IF (ALLOCATED(Gmat)) THEN
        DEALLOCATE(Gmat)
    END IF
	IF (ALLOCATED(GmatInv)) THEN
        DEALLOCATE(GmatInv)
    END IF
	IF (ALLOCATED(Mmat)) THEN
        DEALLOCATE(Mmat)
    END IF
	IF (ALLOCATED(Pmat)) THEN
        DEALLOCATE(Pmat)
    END IF

RETURN

END SUBROUTINE UpdateRefMassFlowRate

!************************************************************************

SUBROUTINE UpdateMCrefMassFlowRate(Slab,Nl,mRefTot,DPcoil)

!------------------------------------------------------------------------
!Purpose:
!To update refrigerant mass flow rate in each microchannel slab
!
!Author
!Ipseng Iu
!Oklahoma State University, Stillwater
!
!Date
!July 2007
!
!------------------------------------------------------------------------

IMPLICIT NONE

TYPE (SlabInfo),ALLOCATABLE,DIMENSION(:),INTENT(INOUT) :: Slab !Coil slab pointer
INTEGER,INTENT(IN) :: Nl !Number of slabs
REAL,INTENT(IN) :: mRefTot !Total refrigerant mass flow rate, kg/s
REAL,INTENT(OUT) :: DPcoil !Coil pressure drop, kPa

REAL,PARAMETER :: Relax=0.6     !Relaxation
INTEGER :: I !Loop counter
REAL :: TotCond !Total refrigerant flow conductance, kg/s-Pa
REAL,ALLOCATABLE,DIMENSION(:) :: mdot !kg/s

  IF (.NOT. ALLOCATED(mdot)) THEN
      ALLOCATE(mdot(Nl))
  ENDIF

  TotCond=0
  DO I=1, Nl    
    Slab(I)%Conduct=Slab(I)%mdot/(Slab(I)%pRi-Slab(I)%pRo)
    TotCond=TotCond+Slab(I)%Conduct
  END DO
  DPcoil=mRefTot/TotCond
  
  DO I=1, Nl
    mdot(I)=Slab(I)%Conduct*DPcoil
  END DO

  !Apply relaxation
  Slab%mdot=(mdot*(Relax)+Slab%mdot*(1-Relax))
  
  DEALLOCATE(mdot)
  
RETURN

END SUBROUTINE UpdateMCrefMassFlowRate

!************************************************************************

END MODULE CoilCalcMod
