! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module simulates the capillary tube in the heat pump cycle
! In previous versions, there were different models, but only one was utilized
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! A capillary tube is a passive device used to meter refrigerant in the system
! Systems using a cap tube do not have a refrigerant receiver, thus charge must be precise
! The device maintains the pressure difference at the entry to the clg coil to ensure the high pressure
!   refrigerant evaporates to provide cooling at the cooling coil
! This was from: http://inspectapedia.com/aircond/Capillary_Tubes.htm#H1

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! The simulation model takes in state points and flow rate, and outputs a (different?) flow rate and state point data

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! na

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! Two module parameters are defined: PI and a unit conversion from feet to meters

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method:
!    PUBLIC CapillaryTubeMod -- This is the only routine here, no calling tree
!      This is called by FlowRateLoop.f90 and HPdesignMod.f90 in conjunction with the rest of the simulation

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! The Choi subroutine was commented out because it wasn't being called anywhere...
!   Should it be re-employed?  Need to figure out if it is better than the original in some way

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2013-01-18 | ESL | Updated header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Why is IDDISTUBE hardwired to a particular size?
! The curve fit is specified to only be valid for R12 and R22...but the curve fit has a reference
! Perhaps the reference can provide more general cases nowadays
! General code cleanup

MODULE CapillaryTubeMod

  USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
  implicit none

  LOGICAL, EXTERNAL :: IssueRefPropError

  PUBLIC  CapillaryTubeORNL
  !PUBLIC  CapillaryTubeChoi

  PRIVATE

  REAL, PARAMETER :: PI=3.14159265
  REAL, PARAMETER :: UnitL=0.3048  !(ft * UnitL = m)

CONTAINS

  !***********************************************************************

  SUBROUTINE CapillaryTubeORNL !(Ref$) !,XIN,PAR,OUT) !(Ref$,PureRef,XIN,PAR,OUT)

    ! ----------------------------------------------------------------------
    !
    !   Description: Capillary tube model
    !
    !   Method: Semi-empirical model from ORNL Heat Pump Model
    !
    !   Inputs:
    !   Ref$=Refrigerant name
    !   PureRef=Refrigerant flag: 1=pure refrigerant
    !                             0=refrigerant mixture
    !  XIN(1)=Compressor mass flow rate, kg/s
    !  XIN(2)=Exp. device inlet pressure, kPa
    !  XIN(3)=Exp. device inlet enthalpy, kJ/kg
    !  XIN(4)=Evaporator inlet pressure, kPa
    !  XIN(5)=Evaporator outlet pressure, kPa
    !
    !  Parameters:
    !  PAR(1) = Capillary tube ID, m
    !  PAR(2) = Capillary tube length, m
    !  PAR(3) = Capillary tube coil diameter, m
    !  PAR(4)=Number of circuits in evaporator
    !  PAR(5)=Distributor tube length, m
    !
    !  Outputs:
    !  OUT(1)=Exp. device mass flow rate, kg/s
    !  OUT(2)=Exp. device outlet pressure, kPa
    !  OUT(3)=Exp. device outlet temperature, C
    !  OUT(4)=Exp. device outlet quality
    !  OUT(5)=Mass in distributor tubes, kg
    !  OUT(6)=Distributor tube capacity, kW
    !  OUT(7) = Error flag: 0-No error
    !                       1-Short tube solution error
    !                       2-Refprop error
    !
    !  Reference:
    !     ORNL Heat Pump Model
    !
    !     Input IREFC, NCAP, and XMR  are passed in COMMON FLOWBA.
    !     Other orifice geometry values are passed in COMMON CAPVAL
    !           and set in BLOCK DATA or optionally in DATAIN.
    !
    !     NOTE: Curve-Fit Equations Are Based On ASHRAE Equipment Handbook,1988,
    !           Chapt. 19, Refrigerant-Control Devices, Figures 38-41,
    !           And Are Strictly Applicable Only For R-12 And R-22.
    !     (If used for other refrigerants, they are only a rough approximation.)
    !
    !      FlowFactor, PCR, and Pcorrection Equations
    !      Adapted by C.K. Rice (April 96)
    !      from S. B. Penson,
    !           Development of a Room Air Conditioner Design Model
    !           (May 1988), and
    !      from D. L. O'Neal and S. B. Penson,
    !           An Analysis of Efficiency Improvements in Room Air Conditioners
    !           ESL/88-04, Texas A&M University, October 1988
    !
    !      Equations modified in content from O'Neal and Penson RAC code
    !      (which differs from the report documentation)
    !      are given by the modification date noted in columns 72-80
    !
    !      Solution for Capillary Tube Diameter DcapTube, Given Mass Flow Rate XMR,
    !      by C.K. Rice, April 1996.
    !
    !   Author:
    !   Ipseng Iu
    !   Building Efficiency
    !   Johnson Control, Inc.
    !   Wichita, KS
    !
    !   Date: April 2009
    !
    ! ----------------------------------------------------------------------

    USE DataSimulation !, ONLY: CTTubeID,CTTubeLen,CTTubeCoilD,CTEvapCktNum,CTDisTubeLen, &   !RS: Debugging: Replacing PAR() numbers with variables
    !                        CTIMdot, CTIPiEx, CTIHiEx, CTIPiEv, CTIPoEv, CTOMdot, CTOErrFlag, CTOToE, CTOXoE, CTOMDT, CTOPoE

    IMPLICIT NONE

    !Subroutine argument declarations
    !CHARACTER*80,     INTENT(IN) :: Ref$    !Refrigerant name
    !INTEGER(2),       INTENT(IN) :: PureRef !Refrigerant flag: 1-pure refrigerant  !RS: Debugging: Extraneous
    !0-refrigerant mixture
    !REAL, INTENT(IN) :: XIN(5)
    !REAL, INTENT(IN) :: PAR(5)
    !REAL, INTENT(OUT) :: OUT(6)   !RS: Debugging: Formerly OUT(7)

    REAL Quality,Pressure,Enthalpy

    REAL :: LcapTube    !Cap tube length, m
    REAL :: LdisTube    !Distributor tube length, m
    REAL :: LCAP        !Distributor tube length, in
    REAL :: DcapTube    !Cap tube diameter, m
    REAL :: DCAP        !Cap tube diameter, in
    REAL :: Dcoil       !Cap tube coil diameter, m
    REAL :: IDDISTUBE   !Distributor inside diameter, in
    REAL :: VolDisTube  !Distributor tube volume, m^3
    REAL :: MassDisTube !Mass in distributor tube, kg
    !REAL :: PiExp       !Inlet pressure of exp. device (Up stream pressure), kPa or psi
    !REAL :: HiExp       !Exp. device inlet enthalpy, kJ/kg
    !REAL :: TiExp       !Exp. device inlet temperature, C
    !REAL :: XiExp       !Exp. device inlet quality
    !REAL :: PoExp       !Exp. device outlet pressure, kPa or psi
    !REAL :: HoExp       !Exp. device outlet enthalpy, kJ/kg
    !REAL :: ToExp       !Exp. device outlet temperature, C
    !REAL :: XoExp       !Exp. device outlet quality
    REAL :: DPdisTot    !Total pressure drop in distributor, kPa
    !REAL :: PiEvp       !Evaporator inlet pressure, kPa
    !REAL :: HiEvp       !Evaporator inlet enthalpy, kJ/kg
    !REAL :: PoEvp       !Evaporator outlet pressure, kPa
    REAL :: mdotCmp     !Mass flow rate from compressor, kg/s
    REAL :: mdotExp     !Mass flow rate from exp. device, kg/s or lbm/hr
    REAL :: TsiExp      !Liquid saturation temperature of the upstream fluid, C
    REAL :: QdisTube    !Distributor tube capacity, kW
    REAL :: rhoiEvp     !Evaporator inlet density, kg/m^3
    REAL :: Acs         !Cross-sectional area, m^2
    REAL :: HoEvpRtd    !Evaporator outlet rated enthalpy, kJ/kg
    REAL :: MDOTSTRAIGHT !Mass flow rate of straight capillary tube, lbm/hr
    REAL :: MDOTBASE    !Baseline mass flow rate, kg/s
    REAL :: Subcooling  !Subcooling, K or R
    REAL :: aa, bb      !Empirical coefficients
    REAL :: PREF         !Reference pressure, psi
    REAL :: PCR,PCRIT1  !Critical pressure, psi
    REAL :: PINLET       !Inlet pressure, psi
    REAL :: LenRatio     !Length ratio
    REAL :: Pratio       !Pressure ratio
    REAL :: Pcorrection  !Pressure correction
    REAL :: FlowFactor   !Flow factor

    INTEGER(2)       :: Nckts       !Number of circuits in evaporator
    INTEGER          :: ErrorFlag   !0-No error
    !1-Solution error
    !2-Refprop error

    !NIST Refrigerant property variables and functions
    INTEGER(2) RefPropErr           !Error flag:1-error; 0-no error

    !Flow:

    mdotCmp = CapTubeIN%CTIMdot    !RS: Debugging: Formerly XIN(1)
    PiExp   = CapTubeIN%CTIPiEx    !RS: Debugging: Formerly XIN(2)
    HiExp   = CapTubeIN%CTIHiEx    !RS: Debugging: Formerly XIN(3)
    PiEvp   = CapTubeIN%CTIPiEv    !RS: Debugging: Formerly XIN(4)
    PoEvp   = CapTubeIN%CTIPoEv    !RS: Debugging: Formerly XIN(5)

    DcapTube = CapTubePAR%CTTubeID   !RS: Debugging: Formerly PAR(1)
    LcapTube = CapTubePAR%CTTubeLen   !RS: Debugging: Formerly PAR(2)
    Dcoil    = CapTubePAR%CTTubeCoilD   !RS: Debugging: Formerly PAR(3)
    Nckts    = CapTubePAR%CTEvapCktNum   !RS: Debugging: Formerly PAR(4)
    LdisTube = CapTubePAR%CTDisTubeLen   !RS: Debugging: Formerly PAR(5)

    ErrorFlag = 0 !Initialize

    Pressure=PiExp*1000 !RS Comment: Unit Conversion
    Enthalpy=HiExp*1000 !RS Comment: Unit Conversion
    TiExp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Inlet Temperature
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, CapTubeOUT%CTOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    XiExp=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)   !Expansion Device Inlet Quality
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag,CapTubeOUT%CTOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    Quality=0
    TsiExp=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Inlet Liquid Saturation Temperature
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, CapTubeOUT%CTOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    HiEvp=HiExp

    Pressure=PiEvp*1000 !RS Comment: Unit Conversion
    Enthalpy=HiEvp*1000 !RS Comment: Unit Conversion

    rhoiEvp=PH(Ref$,Pressure,Enthalpy,'density',RefrigIndex,RefPropErr) !Evaporator Inlet Density
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, CapTubeOUT%CTOErrFlag)) THEN       !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    IF (LdisTube .GT. 0) THEN
      LDISTUBE=LdisTube*12/UnitL !Distributor tube length, in
      CALL Distributor(Ref$,LDISTUBE,Nckts,mdotCmp,TiExp,HiExp,PoEvp, &
      HoEvpRtd,QdisTube,DPdisTot,ErrorFlag)

      IF (DPdisTot .LT. 0) THEN
        DPdisTot = 0
      END IF
      IF (DPdisTot .GT. PiEvp) THEN
        DPdisTot = 0
      END IF

      IDDISTUBE=1./4.-12.0E-3 !1/4 in OD, 12 mil-think
      VolDisTube=((IDDISTUBE/12**2)*PI/4*LDISTUBE/12*Nckts)*UnitL**3  !Distributor tube volume, m^3
      MassDisTube=rhoiEvp*VolDisTube

      PoExp=PiEvp+DPdisTot

    ELSE
      PoExp=PiEvp
      MassDisTube=0
    END IF

    HoExp=HiExp

    Pressure=PoExp*1000 !RS Comment: Unit Conversion
    Enthalpy=HoExp*1000 !RS Comment: Unit Conversion
    ToExp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Outlet Temperature
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, CapTubeOUT%CTOErrFlag)) THEN    !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    XoExp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Expansion Device Outlet Quality
    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, CapTubeOUT%CTOErrFlag)) THEN   !RS: Debugging: Formerly OUT(7)
      RETURN
    END IF

    Subcooling=TsiExp-TiExp

    IF (XiExp .LT. 0.0) THEN
      XiExp = 0.0
    END IF

    Acs=PI*(DcapTube**2)/4

    !ORNL Algorithm

    PIEXP=PiExp*0.14503798
    POEXP=PoExp*0.14503798
    SUBCOOLING=Subcooling*1.8
    LCAP=LcapTube*39.3701
    DCAP=DcapTube*39.3701

    IF (XiExp .GT. 0) THEN
      aa = 0.8184 + 0.0264*(100.*XiExp)**0.624
      bb = 534.4
      PREF = 2200.
      PCRIT1 = EXP( -1.+LOG(PIEXP)-2.27*XiExp+0.623*XiExp**2+0.235*XiExp*LOG(PIEXP) )
    ELSE
      aa = 0.4035 + 0.4175*EXP(-0.04*SUBCOOLING)
      bb  = 356.0  + 0.641*((ABS(SUBCOOLING-31.0)/10.)**3.56)
      PREF = 1500.
      PCRIT1 = EXP(-0.903+0.979*LOG(PIEXP)+0.624E-4*SUBCOOLING**2)
    END IF

    MDOTBASE= bb*(PiExp/Pref)**aa
    PINLET = EXP(6.8925-0.6109*LOG(PCRIT1))

    !**CALCULATE CAPILLARY TUBE FLOW RATE
    !  GIVEN TUBE NUMBER, DIAMETER, LENGTH, CONDITIONS, AND MASS FLOW RATE

    !CALCULATE CAPILLARY CRITICAL PRESSURE TO DETERMINE IF FLOW IS CHOKED
    LenRatio=LcapTube/(1250.*DcapTube)
    If (LenRatio .GT. 1.17) Then
      PCR = EXP(11.5292-0.4922*(LenRatio)-1.6309*LOG(PINLET))+0.0317* (LenRatio**2+0.0265*LenRatio*LOG(PINLET))
    Else
      PCR = EXP(12.2892-0.9567*LenRatio-1.7491*LOG(PINLET)+ 0.0937*LenRatio*LOG(PINLET))
    EndIf

    Pratio = (POEXP-PCR)/(PIEXP-PCR)

    !CALCULATE PRESSURE CORRECTION FACTOR (Pcorrection) IF FLOW IS NOT CHOKED
    !(OTHERWISE Pcorrection=1)

    If((Pratio.GE.0.0).and.(Pratio.LT.1.0)) Then
      Pcorrection=1.00-6.249886e-2*Pratio-7.291768e-1*Pratio**2+9.635722e-1* &
      Pratio**3-5.208687e-1*Pratio**4-6.510278e-1*Pratio**5
    Else
      Pcorrection=1.0
    EndIf

    !CALCULATE FLOW FACTOR FOR GIVEN TUBE LENGTH AND DIAMETER
    FlowFactor = EXP(7.1539+2.2088*LOG(DCAP)+0.1619*LOG(LCAP)+0.0807*   &
    LOG(LCAP)*LOG(DCAP)-0.0421*(LOG(LCAP))**2)/1.011

    !CALCULATE TOTAL FLOW THROUGH ALL CAPILLARY TUBES
    MDOTSTRAIGHT = FlowFactor*MDOTBASE*Pcorrection

    mdotExp=MDOTSTRAIGHT*(2.011*(LcapTube/DcapTube)**(-0.094)*(LcapTube/Dcoil)**(-0.0527))
    mdotExp=mdotExp*1.26E-4 !Convert from lbm/hr to kg/s
    PoExp=POEXP/0.14503798 !Convert from psi to kPa

    CapTubeOUT%CTOMdot=mdotExp  !RS: Debugging: Formerly OUT(1)
    CapTubeOUT%CTOToE=ToExp    !RS: Debugging: Formerly OUT(3)
    CapTubeOUT%CTOXoE=XoExp    !RS: Debugging: Formerly OUT(4)
    CapTubeOUT%CTOMDT=MassDisTube !RS: Debugging: Only used for output !RS: Debugging: Formerly OUT(5)
    CapTubeOUT%CTOPoE=PoExp   !RS: Debugging: Formerly the pressure
    CapTubeOUT%CTOErrFlag=ErrorFlag   !RS: Debugging: Formerly OUT(2)

    RETURN

  END SUBROUTINE CapillaryTubeORNL

  !***********************************************************************

  !    SUBROUTINE CapillaryTubeChoi(Ref$,PureRef,XIN,PAR,OUT)

  !    ! ----------------------------------------------------------------------
  !    !
  !    !   Description: Capillary tube model
  !    !
  !    !   Method: Buckingham Pi model
  !    !
  !    !   Inputs:
  !    !   Ref$=Refrigerant name
  !    !   PureRef=Refrigerant flag: 1=pure refrigerant
  !    !                             0=refrigerant mixture
  !    !  XIN(1)=Compressor mass flow rate, kg/s
  !    !  XIN(2)=Exp. device inlet pressure, kPa
  !    !  XIN(3)=Exp. device inlet enthalpy, kJ/kg
  !    !  XIN(4)=Evaporator inlet pressure, kPa
  !    !  XIN(5)=Evaporator outlet pressure, kPa
  !    !
  !    !  Parameters:
  !    !  PAR(1) = Capillary tube ID, m
  !    !  PAR(2) = Capillary tube length, m
  !    !  PAR(3) = Capillary tube coil diameter, m
  !    !  PAR(4)=Number of circuits in evaporator
  !    !  PAR(5)=Distributor tube length, m
  !    !
  !    !  Outputs:
  !    !  OUT(1)=Exp. device mass flow rate, kg/s
  !    !  OUT(2)=Exp. device outlet pressure, kPa
  !    !  OUT(3)=Exp. device outlet temperature, C
  !    !  OUT(4)=Exp. device outlet quality
  !    !  OUT(5)=Mass in distributor tubes, kg
  !    !  OUT(6)=Distributor tube capacity, kW
  !    !  OUT(7) = Error flag: 0-No error
  !    !                       1-Short tube solution error
  !    !                       2-Refprop error
  !    !
  !    !  Reference:
  !    !  Choi, J.; Kim, Y.; Chung, J.T. (2004). An empirical correlation and
  !    !  rating charts for the performance of adiabatic capillary tubes
  !    !  with alternative refrigerants. Applied thermal engineering, 24,
  !    !  pp. 29-41.
  !    !
  !    !  Author:
  !    !  Ipseng Iu
  !    !  Building Efficiency
  !    !  Johnson Control, Inc.
  !    !  Wichita, KS
  !    !
  !    !  Date: April 2009
  !    !
  !    ! ----------------------------------------------------------------------

  !    IMPLICIT NONE

  !    !Subroutine argument declarations
  !    CHARACTER*80,     INTENT(IN) :: Ref$    !Refrigerant name
  !    INTEGER(2),       INTENT(IN) :: PureRef !Refrigerant flag: 1-pure refrigerant
  !    !0-refrigerant mixture
  !    REAL, INTENT(IN) :: XIN(5)
  !    REAL, INTENT(IN) :: PAR(5)
  !    REAL, INTENT(OUT) :: OUT(7)

  !    !INTEGER         :: RefrigIndex =0
  !    REAL Temperature,Quality,Pressure,Enthalpy

  !    REAL :: LcapTube    !Cap tube length, m
  !    REAL :: LdisTube    !Distributor tube length, m
  !    REAL :: LCAP        !Distributor tube length, in
  !    REAL :: DcapTube    !Cap tube diameter, m
  !    REAL :: DCAP        !Cap tube diameter, in
  !    REAL :: Dcoil       !Cap tube coil diameter, m
  !    REAL :: IDDISTUBE   !Distributor inside diameter, in
  !    REAL :: VolDisTube  !Distributor tube volumn, m^3
  !    REAL :: MassDisTube !Mass in distributor tube, kg
  !    REAL :: PiExp       !Inlet pressure of exp. device (Up stream pressure), kPa or psi
  !    REAL :: HiExp       !Exp. device inlet enthalpy, kJ/kg
  !    REAL :: TiExp       !Exp. device inlet temperature, C
  !    REAL :: XiExp       !Exp. device inlet quality
  !    REAL :: PoExp       !Exp. device outlet pressure, kPa or psi
  !    REAL :: HoExp       !Exp. device outlet enthalpy, kJ/kg
  !    REAL :: ToExp       !Exp. device outlet temperature, C
  !    REAL :: XoExp       !Exp. device outlet quality
  !    REAL :: DPdisTot    !Total pressure drop in distributor, kPa
  !    REAL :: PiEvp       !Evaporator inlet pressure, kPa
  !    REAL :: HiEvp       !Evaporator inlet enthalpy, kJ/kg
  !    REAL :: TiEvp       !Evaporator inlet temperature, C
  !    REAL :: XiEvp       !Evaporator inlet quality
  !    REAL :: PoEvp       !Evaporator outlet pressure, kPa
  !    REAL :: mdotCmp     !Mass flow rate from compressor, kg/s
  !    REAL :: mdotExp     !Mass flow rate from exp. device, kg/s or lbm/hr
  !    REAL :: TsiExp      !Liquid saturation temperature of the upstream fluid, C
  !    REAL :: QdisTube    !Distributor tube capacity, kW
  !    REAL :: rhoiEvp     !Evaporator inlet density, kg/m^3
  !    REAL :: Acs         !Cross-sectional area, m^2
  !    REAL :: HoEvpRtd    !Evaporator outlet rated enthalpy, kJ/kg
  !    REAL :: mdotStraight !Mass flow rate of straight capillary tube, kg/hr
  !    REAL :: Subcooling  !Subcooling, K

  !    REAL :: sigmaExp !Surface tension, N/m
  !    REAL :: rhofiExp !Liquid density, kg/m3
  !    REAL :: rhogiExp !Vapor density, kg/m3
  !    REAL :: mufiExp !Liquid dynamic viscosity, Pa-s
  !    REAL :: mugiExp !Vapor dynamic viscosity, Pa-s
  !    REAL :: hfiExp !Liquid enthalpy, kJ/kg
  !    REAL :: hgiExp !Vapor enthalpy, kJ/kg
  !    REAL :: hfgExp !Latent heat, kJ/kg
  !    REAL :: Pcr !Critical pressure, kPa
  !    REAL :: Tcr !Critical temperature, C
  !    REAL :: Psat !Saturated pressure, kPa

  !    REAL :: Pi1 !Mass flow rate dimensionless group
  !    REAL :: Pi2 !Inlet pressure dimensionless group
  !    REAL :: Pi3 !Subcooling dimensionless group
  !    REAL :: Pi4 !Geometry dimensionless group
  !    REAL :: Pi5 !Density dimensionless group
  !    REAL :: Pi6 !Friction, bubble growth dimensionless group
  !    REAL :: Pi7 !Friction, bubble growth dimensionless group
  !    REAL :: Pi8 !Vaporization dimensionless group

  !    INTEGER(2)      :: Nckts       !Number of circuits in evaporator
  !    INTEGER             :: ErrorFlag   !0-No error
  !    !1-Solution error
  !    !2-Refprop error
  !    INTEGER :: I  !Iteration counter
  !    INTEGER :: II !Iteration counter

  !    !NIST Refrigerant property variables and functions
  !    INTEGER(2) RefPropOpt          !Ref prop calc. option
  !    INTEGER(2) RefPropErr          !Error flag:1-error; 0-no error
  !    REAL RefProp(28)

  !    !Flow:

  !    mdotCmp = XIN(1)
  !    PiExp   = XIN(2)
  !    HiExp   = XIN(3)
  !    PiEvp   = XIN(4)
  !    PoEvp   = XIN(5)

  !    DcapTube = PAR(1)
  !    LcapTube = PAR(2)
  !    Dcoil    = PAR(3)
  !    Nckts    = PAR(4)
  !    LdisTube = PAR(5)

  !    ErrorFlag = 0 !Initialize

  !    Pressure=PiExp*1000 !RS Comment: Unit Conversion
  !    Enthalpy=HiExp*1000 !RS Comment: Unit Conversion
  !    TiExp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Inlet Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    XiExp=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)   !Expansion Device Inlet Quality
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    sigmaExp=PH(Ref$,Pressure,Enthalpy,'surfacetension',RefrigIndex,RefPropErr) !Expansion Device Surface Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF
  !    sigmaExp=sigmaExp/1000  !RS Comment: Unit Conversion

  !    Quality=0
  !    TsiExp=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Liquid Saturation Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    rhofiExp=PQ(Ref$,Pressure,Quality,'density',RefrigIndex,RefPropErr) !Expansion Device Liquid Density
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    mufiExp=PQ(Ref$,Pressure,Quality,'viscosity',RefrigIndex,RefPropErr)    !Expansion Device Liquid Dynamic Viscosity
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    hfiExp=PQ(Ref$,Pressure,Quality,'enthalpy',RefrigIndex,RefPropErr)  !Expansion Device Liquid Enthalpy
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF
  !    hfiExp=hfiExp/1000  !RS Comment: Unit Conversion

  !    Quality=1
  !    rhogiExp=PQ(Ref$,Pressure,Quality,'density',RefrigIndex,RefPropErr) !Expansion Device Vapor Density
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    mugiExp=PQ(Ref$,Pressure,Quality,'viscosity',RefrigIndex,RefPropErr)    !Expansion Device Vapor Dynamic Viscosity
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    hgiExp=PQ(Ref$,Pressure,Quality,'enthalpy',RefrigIndex,RefPropErr)  !Expansion Device Vapor Enthalpy
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF
  !    hgiExp=hgiExp/1000  !RS Comment: Unit Conversion

  !    HiEvp=HiExp

  !    Pressure=PiEvp*1000 !RS Comment: Unit Conversion
  !    Enthalpy=HiEvp*1000 !RS Comment: Unit Conversion
  !    TiEvp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Evaporator Inlet Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    XiEvp=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)   !Evaporator Inlet Quality
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    rhoiEvp=PH(Ref$,Pressure,Enthalpy,'density',RefrigIndex,RefPropErr) !Evaporator Inlet Density
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    IF (LdisTube .GT. 0) THEN
  !        LDISTUBE=LdisTube*12/UnitL !Distributor tube length, in
  !        CALL Distributor(Ref$,LDISTUBE,Nckts,mdotCmp,TiExp,HiExp,PoEvp, &
  !        HoEvpRtd,QdisTube,DPdisTot,ErrorFlag)

  !        IF (DPdisTot .LT. 0) THEN
  !            DPdisTot = 0
  !        END IF
  !        IF (DPdisTot .GT. PiEvp) THEN
  !            DPdisTot = 0
  !        END IF

  !        IDDISTUBE=1./4.-12.0E-3 !1/4 in OD, 12 mil-think
  !        VolDisTube=((IDDISTUBE/12**2)*PI/4*LDISTUBE/12*Nckts)*UnitL**3  !Distributor tube volume, m^3
  !        MassDisTube=rhoiEvp*VolDisTube

  !        PoExp=PiEvp+DPdisTot

  !    ELSE
  !        PoExp=PiEvp
  !        MassDisTube=0
  !    END IF

  !    HoExp=HiExp

  !    Pressure=PoExp*1000 !RS Comment: Unit Conversion
  !    Enthalpy=HoExp*1000 !RS Comment: Unit Conversion
  !    ToExp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Expansion Device Outlet Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    XoExp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Expansion Device Outlet Quality
  !    IF (IssueRefPropError(RefPropErr, 'Capillary Tube', ErrorFlag, OUT(7))) THEN
  !        RETURN
  !    END IF

  !    Temperature=TiExp
  !    Quality=0
  !    Psat=TQ(Ref$, Temperature, Quality, 'pressure', RefrigIndex,RefPropErr) !Saturation Temperature
  !    IF (IssueRefPropError(RefPropErr, 'Short Tube', ErrorFlag)) THEN
  !        RETURN
  !    END IF
  !    Psat=Psat/1000  !RS Comment: Unit Conversion

  !    hfgExp=hgiExp-hfiExp
  !    Subcooling=TsiExp-TiExp

  !    Tcr=Tcrit(Ref$)
  !    Pcr=Pcrit(Ref$)/1000    !RS Comment: Unit Conversion

  !    Pi2=(PiExp-Psat)/Pcr
  !    Pi3=Subcooling/Tcr
  !    Pi4=LcapTube/DcapTube
  !    Pi5=rhofiExp/rhogiExp
  !    Pi6=(mufiExp-mugiExp)/mugiExp
  !    Pi7=sigmaExp/(DcapTube*1000*PiExp)
  !    Pi8=rhofiExp*hfgExp/(Psat)

  !    Pi1 = 0.5782E-4 * Pi2**(-0.315) * Pi3**(0.369) * Pi4**(-0.344) * Pi5**(0.034) * &
  !    Pi6**(0.04) * Pi7**(-0.458) * Pi8**(0.376)

  !    mdotStraight=Pi1*((DcapTube*1000)**2 * SQRT(rhofiExp * PiExp))
  !    mdotStraight=mdotStraight/3600 !Convert kg/hr to kg/s

  !    mdotExp=mdotStraight*(2.011*(LcapTube/DcapTube)**(-0.094)*(LcapTube/Dcoil)**(-0.0527))

  !    OUT(1)=mdotExp
  !    OUT(2)=PoExp
  !    OUT(3)=ToExp
  !    OUT(4)=XoExp
  !    OUT(5)=MassDisTube
  !    OUT(6)=QdisTube

  !    OUT(7)=ErrorFlag

  !    RETURN

  !    END SUBROUTINE CapillaryTubeChoi

  !***********************************************************************

END MODULE CapillaryTubeMod
