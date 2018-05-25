! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module contains methods and variables related to the simulation of an accumulator component in a heat pump cycle
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component represents an accumulator in a heat pump system.
! A description of the component is found at:
! http://www.e-refrigeration.com/learn-refrigeration/system-components/refrigeration-accumulator
! From that website:
!  - Accumulators store excess liquid refrig and oil that didn't boil in evap
!  - It is situated between the evap and compressor in suction line
!  - It operates as the liquid strikes a plate and reflects back into a holding tank

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! I don't want to write this right now.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! no OPEN statements found in this source file
! a single WRITE statement writes to UNIT=6 (stdout)
!...assuming no file I/O

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! A small set of variables (~15) are defined at module level, with units included for most

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains three methods:
!    PUBLIC InitAccumulator -- initialize module level data structures
!      Called by ORNLSolver.f90
!    PUBLIC CalcAccumulatorMass -- calculate the mass of refrigerant and liquid level in the accumulator, by perform pressure balance in accumulator
!      Called by HPdesignMod.f90
!    PUBLIC CalcAccumulatorDP -- calculate the pressure drop across the accumulator
!      Called by AirTempLoop.f90

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! na

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 02/03/2014 | Karthik | Variable Naming and cleanup of code

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Check if other module variables can be made method locals
! ErrorFlag values will be set from the Method IssueRefPropError. This is
! not a standard way of doing this and needs to check more on this (Karthik)
! -> Elimanate un necessary Temp Variables.
! -> Variable Names has been changed to the best knowledge, if someone would review the naming convention that would be useful

MODULE AccumulatorModule

  USE DataGlobals_HPSimIntegrated, ONLY: RefName, RefrigIndex !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
  USE DataSimulation                                !Karthik - Use the Global Variables to get data for this module.

  IMPLICIT NONE
  LOGICAL, EXTERNAL :: IssueRefPropError             !Karthik - Print out Any Errors

  PRIVATE

  REAL :: RatedPressureDrop         !Unit - kPa
  REAL :: RatedTemperatureDrop      !Unit -  K
  REAL :: CoefficientOfCurveFit_M
  REAL :: CoefficientOfCurveFit_B
  REAL :: PressureDropInAccumulator !Unit -  kPa
  REAL DiameterHole(2)              !Unit -  ft
  REAL TubeDiameter                 !Unit - ft
  REAL distanceBetweenHoles         !Unit - ft
  REAL accumulatorHeight            !Unit - ft
  REAL accumulatorDiameter          !Unit - ft
  INTEGER(2) RefPropErr             !Error flag:1-error; 0-no error
  REAL Temperature,Quality,Pressure,Enthalpy,Entropy

  !Subroutines
  PUBLIC InitAccumulator
  PUBLIC CalcAccumulatorMass
  PUBLIC CalcAccumulatorDP


CONTAINS
  ! ----------------------------------------------------------------------
  !
  !   Description: HPSIM accumulator model
  !
  !   Purpose: To calculate the mass of refrigerant and liquid level
  !            in the accumulator.
  !
  !   Method: Perform pressure balance in accumulator
  !
  !   Inputs:
  !       RefName=Refrigerant name
  !       PureRef=Refrigerant flag: 1=pure refrigerant
  !                                 0=refrigerant mixture
  !       XIN(1) = Mass flow rate, kg/s
  !       XIN(2) = Outlet refrigerant pressure, kPa
  !       XIN(3) = Outlet refrigerant enthalpy, kJ/kg
  !
  !
  !
  !   Outputs:
  !       OUT(1) = Total mass inventory, kg
  !       OUT(2) = Pressure drop, kPa
  !       OUT(3) = Error Flag
  !
  !   Author:
  !   Ipseng Iu
  !   Mechanical and Aerospace Engineering
  !   Oklahoma State University, Stillwater
  !
  !   Date: June 2005
  !
  !   NOTE - Old Comments has not been deleted.
  !
  !**** PURPOSE:
  !       TO CALCULATE REFRIGERANT LIQUID LEVEL AND MASS IN AN ACCUMULATOR
  !
  !     ADAPTED 10/85 BY C. K. RICE
  !       FROM 'ACCUM' ROUTINE WRITTEN BY PIOTR A. DOMANSKI,
  !       NATIONAL BUREAU OF STANDARDS, 6/28/82
  !
  !**** INPUT DATA:
  !       accumulatorHeight(ACCHGT)    - ACCUMULATOR HEIGHT  (FT)
  !       accumulatorDiameter(ACCDIA)    - INNER DIAMETER OF ACCUMULATOR  (FT)
  !       DiameterHole(1)        - INNER DIA. OF OIL RETURN HOLE  (FT)
  !       DiameterHole(2)        - INNER DIA. OF UPPER HOLE  (FT)
  !       TubeDiameter(ATBDIA)   - INNER DIA. OF ACCUMULATOR TUBE  (FT)
  !       distanceBetweenHoles(HOLDIS)    - VERTICAL DISTANCE BETWEEN HOLES  (FT)
  !       RMASS           - REFRIG. MASS FLOW RATE  (LBM/H)
  !       refrigerantSuperHeat            - REFRIG. SUPERHEAT  (F)
  !       volumeOfRefLeavingAccumulator             - SPECIFIC VOLUME OF REFRIGERANT LEAVING ACCUMULATOR
  !                         (FT**3/LBM)
  !       volumeOfLiquidAtAccuPressure          - SPECIFIC VOLUME OF LIQUID AT ACCUMULATOR PRESSURE
  !                         (FT**3/LBM)
  !**** OUTPUT DATA:
  !       refrigerantMassInAccumulator          - MASS OF REFRIG. IN THE ACCUMULATOR  (LBM)
  !       volumeOfAccumulatorInternal          - INTERNAL VOLUME OF ACCUMULATOR  (FT**3)
  !       liquidLevelInAccumulator          - LEVEL OF LIQUID IN ACCUMULATOR  (IN)
  !
  !        COMMON / ACCDIM / ACCHGT,ACCDIA,DiameterHole(2),ATBDIA,HOLDIS
  !
  ! ----------------------------------------------------------------------
  SUBROUTINE CalcAccumulatorMass

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    IMPLICIT NONE

    INTEGER ErrorFlag
    !0-No error
    !1-Accumulator solution not converge
    !2-Refprop error

    REAL refrigerantMassFlowRate            !Refrigerant mass flow rate, kg/s
    REAL refrigerantOutletPressure          !Outlet refrigerant pressure, kPa
    REAL refrigerantOutletEnthalpy          !Outlet refrigerant enthalpy, kJ/kg
    REAL refrigerantOutletTemperature       !Outlet refrigerant temperature, C
    REAL saturatedTemperature               !Saturated temperature, C
    REAL refrigerantOutletQuality           !Outlet refrigerant quality, -
    REAL :: vapourDensity, volumeOfRefLeavingAccumulator, liquidDensity, volumeOfLiquidAtAccuPressure, refrigerantSuperHeat
    REAL :: convergenceToleranceForAccum, volumeOfAccumulatorInternal
    REAL :: refrigerantMassInAccumulator, liquidLevelInAccumulator, areaOfAccumulator, accumulatorVerticalHeight
    REAL :: refrigerantMassFlowRateInLBM
    REAL :: refrigerantMassFlowOfLiquid, liquidDensity_Temp, areaOfTube, refrigerantMassFlowPerHoleArea
    REAL :: refrigerantLiqPressureAtHole, refrigerantVapourPressureAtTube, totalRefrigerantMassinAcc_2, maximumPressureWithOutHoleDist
    REAL :: maximumPressureWithHoleDist, refrigerantMassFlowAtTopHoleMax, pressureAtHole2, refrigerantMassFlowAtTopHole
    REAL :: differenceInMassFlow_1, totalRefrigerantMassinAcc_temp, Z1, Z2
    REAL :: indexOfHolePresent, J, refrigerantMass_Temp, I, differenceInMassFlow_2, SLOPE, ERROR
    INTEGER   :: boolTerminationFlag            !Termination Flag, change the variable name if you need to.
    CHARACTER(LEN=57),PARAMETER :: FMT_602 = "(' ACCUML DOES NOT CONVERGE, MAX.ERROR =',1PE10.3,' LBM')"
    REAL liquidColumnHeight, areaOfHole
    DIMENSION liquidColumnHeight(2), areaOfHole(2)
    REAL,PARAMETER :: PI=3.14159265

    boolTerminationFlag = .FALSE.
    refrigerantMassFlowRate=AccumIN%AccImdot   !RS: Debugging: Formerly XIN(1)
    refrigerantOutletPressure=AccumIN%AccIpRo  !RS: Debugging: Formerly XIN(2)
    refrigerantOutletEnthalpy=AccumIN%AccIhRo  !RS: Debugging: Formerly XIN(3)
    ErrorFlag=0
    Pressure=refrigerantOutletPressure*1000   !RS Comment: Unit Conversion
    Enthalpy=refrigerantOutletEnthalpy*1000   !RS Comment: Unit Conversion
    refrigerantOutletTemperature=PH(RefName,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)  !Outlet Refrigerant Temperature

    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN    !RS: Debugging: Formerly ", OUT(3)"
      RETURN
    END IF
    refrigerantOutletQuality=PH(RefName,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)  !Outlet Refrigerant Quality
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN    !RS: Debugging: Formerly ", OUT(3)"
      RETURN
    END IF

    Pressure=refrigerantOutletPressure*1000   !RS Comment: Unit Conversion
    Quality=1
    saturatedTemperature=PQ(RefName,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)  !Saturated Temperature
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN    !RS: Debugging: Formerly ", OUT(3)"
      RETURN
    END IF

    vapourDensity=PQ(RefName,Pressure,Quality,'density',RefrigIndex,RefPropErr)    !Vapor Density
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN    !RS: Debugging: Formerly ", OUT(3)"
      RETURN
    END IF

    volumeOfRefLeavingAccumulator=1/(vapourDensity*0.0625) !Convert from kg/m3 to lbm/ft3
    Quality=0
    liquidDensity=PQ(RefName,Pressure,Quality,'density',RefrigIndex,RefPropErr)    !Liquid Density
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN    !RS: Debugging: Formerly ", OUT(3)"
      RETURN
    END IF

    volumeOfLiquidAtAccuPressure=1/(liquidDensity*0.0625) !Convert from kg/m3 to lbm/ft3
    refrigerantSuperHeat=(refrigerantOutletTemperature-saturatedTemperature)*1.8

    convergenceToleranceForAccum = 0.001
    volumeOfAccumulatorInternal = 0.0
    refrigerantMassInAccumulator = 0.0
    liquidLevelInAccumulator = 0.0

    IF(accumulatorHeight.GT.0.001) THEN
      areaOfAccumulator = PI*accumulatorDiameter*accumulatorDiameter/4.
      volumeOfAccumulatorInternal = areaOfAccumulator*accumulatorHeight
      IF(.NOT. refrigerantOutletQuality.LT.1.) THEN       !Check if Only Vapor in Accumulator
        liquidColumnHeight(1) = 0.0
        refrigerantMassInAccumulator = volumeOfAccumulatorInternal/volumeOfRefLeavingAccumulator
        accumulatorVerticalHeight=accumulatorHeight
      ELSE
        !
        !       LIQUID IN ACCUMULATOR
        refrigerantMassFlowRateInLBM=refrigerantMassFlowRate/0.000126     !Convert from kg/s to lbm/hr
        refrigerantMassFlowRateInLBM = refrigerantMassFlowRateInLBM/3600. !Convert from kg/s to lbm/hr
        refrigerantMassFlowOfLiquid = (1.-refrigerantOutletQuality)*refrigerantMassFlowRateInLBM
        liquidDensity_Temp = 1./volumeOfLiquidAtAccuPressure
        areaOfHole(1) = PI*DiameterHole(1)*DiameterHole(1)/4. !LIQUID LEVEL FOR CASE OF BOTTOM HOLE ALONE
        areaOfTube = PI*TubeDiameter*TubeDiameter/4.
        refrigerantMassFlowPerHoleArea = refrigerantMassFlowOfLiquid/areaOfHole(1)
        refrigerantLiqPressureAtHole = refrigerantMassFlowPerHoleArea*refrigerantMassFlowPerHoleArea/(.585*.585*2.*liquidDensity_Temp)
        refrigerantVapourPressureAtTube = 0.5*(refrigerantOutletQuality*refrigerantMassFlowRateInLBM/areaOfTube)**2.*volumeOfRefLeavingAccumulator
        liquidColumnHeight(1) = (refrigerantLiqPressureAtHole-refrigerantVapourPressureAtTube)/(liquidDensity_Temp*32.2)

        IF(liquidColumnHeight(1).LT.0.) THEN
          liquidColumnHeight(1) = 0.
        END IF

        accumulatorVerticalHeight = accumulatorHeight-liquidColumnHeight(1)
        accumulatorVerticalHeight = AMAX1(0.,accumulatorVerticalHeight)
        totalRefrigerantMassinAcc_2 = areaOfAccumulator*(liquidColumnHeight(1)*liquidDensity_Temp+accumulatorVerticalHeight/volumeOfRefLeavingAccumulator)

        IF(liquidColumnHeight(1).LT.distanceBetweenHoles.OR.distanceBetweenHoles.EQ.0.) THEN !CHECK IF BELOW SECOND HOLE
          boolTerminationFlag = .TRUE.
        END IF
        IF (boolTerminationFlag .EQ. .FALSE.) THEN
          areaOfHole(2) = PI*DiameterHole(2)*DiameterHole(2)/4.
          !
          !       CHECK FOR FULL ACCUMULATOR
          !
          maximumPressureWithOutHoleDist = accumulatorHeight*liquidDensity_Temp*32.2 + refrigerantVapourPressureAtTube
          maximumPressureWithHoleDist = (accumulatorHeight-distanceBetweenHoles)*liquidDensity_Temp*32.2 + refrigerantVapourPressureAtTube
          refrigerantMassFlowAtTopHoleMax = 0.585*areaOfHole(1)*SQRT(2.*liquidDensity_Temp*maximumPressureWithOutHoleDist) +0.585*areaOfHole(2)*SQRT(2.*liquidDensity_Temp*maximumPressureWithHoleDist)

          IF(refrigerantMassFlowOfLiquid.GT.refrigerantMassFlowAtTopHoleMax) THEN
            accumulatorVerticalHeight = 0.0
            liquidColumnHeight(1) = accumulatorHeight
            totalRefrigerantMassinAcc_2 = areaOfAccumulator*accumulatorHeight*liquidDensity_Temp
            boolTerminationFlag = .TRUE. ! FIND LEVEL ABOVE SECOND HOLE
          END IF
        END IF

        IF (boolTerminationFlag .EQ. .FALSE.) THEN
          liquidColumnHeight(2) = liquidColumnHeight(1)-distanceBetweenHoles
          pressureAtHole2 = liquidColumnHeight(2)*liquidDensity_Temp*32.2 + refrigerantVapourPressureAtTube
          refrigerantMassFlowAtTopHole = 0.585*areaOfHole(2)*SQRT(2.*liquidDensity_Temp*pressureAtHole2)
          differenceInMassFlow_1 = refrigerantMassFlowAtTopHole
          totalRefrigerantMassinAcc_temp = totalRefrigerantMassinAcc_2
          Z1 = liquidColumnHeight(1)
          Z2 = distanceBetweenHoles
          indexOfHolePresent = 0
          DO J=1,12
            liquidColumnHeight(1) = Z2
            liquidColumnHeight(2) = Z2-distanceBetweenHoles
            refrigerantMass_Temp = 0.
            DO I=1,2
              pressureAtHole2 = liquidColumnHeight(I)*liquidDensity_Temp*32.2 + refrigerantVapourPressureAtTube
              refrigerantMassFlowAtTopHole = 0.585*areaOfHole(I)*SQRT(2.*liquidDensity_Temp*pressureAtHole2)
              refrigerantMass_Temp = refrigerantMass_Temp + refrigerantMassFlowAtTopHole
            END DO
            differenceInMassFlow_2 = refrigerantMass_Temp-refrigerantMassFlowOfLiquid
            IF(J.EQ.1.AND.differenceInMassFlow_2.GT.0.0) THEN
              indexOfHolePresent = 1
            END IF
            accumulatorVerticalHeight = accumulatorHeight - liquidColumnHeight(1)
            accumulatorVerticalHeight = AMAX1(0.,accumulatorVerticalHeight)
            totalRefrigerantMassinAcc_2 = areaOfAccumulator* (liquidColumnHeight(1)*liquidDensity_Temp + accumulatorVerticalHeight/volumeOfRefLeavingAccumulator)
            IF(indexOfHolePresent.EQ.1) THEN !SKIP OUT IF LEVEL WILL NOT RISE ABOVE SECOND HOLE
              boolTerminationFlag = .TRUE.
              EXIT
            END IF
            IF(ABS(totalRefrigerantMassinAcc_temp-totalRefrigerantMassinAcc_2).LT.convergenceToleranceForAccum) THEN  !CHECK FOR CONVERGENCE ON LEVEL ABOVE SECOND HOLE
              boolTerminationFlag = .TRUE.
              EXIT
            END IF
            SLOPE = (Z1-Z2)/(differenceInMassFlow_1-differenceInMassFlow_2)
            IF(.NOT. ABS(differenceInMassFlow_1).LT.ABS(differenceInMassFlow_2)) THEN
              totalRefrigerantMassinAcc_temp = totalRefrigerantMassinAcc_2
              differenceInMassFlow_1 = differenceInMassFlow_2
              Z1 = Z2
            END IF
            Z2 = Z1 - differenceInMassFlow_1*SLOPE
            Z2 = AMAX1(Z2,distanceBetweenHoles)
          END DO

          IF (boolTerminationFlag .EQ. .FALSE.) THEN !PRINT RESULTS
            ERROR = areaOfAccumulator*ABS(Z2-Z1)*liquidDensity_Temp
            WRITE(6,FMT_602) ERROR
          END IF

        END IF
        refrigerantMassInAccumulator = totalRefrigerantMassinAcc_2
      END IF
      liquidLevelInAccumulator = liquidColumnHeight(1)*12.
    END IF

    AccumOUT%AccOMass=refrigerantMassInAccumulator*0.4563    !Total mass, convert from lbm to kg !RS: Debugging: Formerly OUT(1)
    AccumOUT%AccODP=PressureDropInAccumulator  !RS: Debugging: Formerly OUT(2)

    RETURN

  END SUBROUTINE CalcAccumulatorMass

  !******************************************************************************

  SUBROUTINE CalcAccumulatorDP(refrigerantMassFlowRate,evaporatingTemperature,condensingTemperature,Subcooling,Superheat,liquidLineQuality,suctionLineQuality,estimatedPressureDrop)

    ! ----------------------------------------------------------------------
    !
    !   Purpose: To calculate the pressure drop.
    !
    !   Method: Catalog data curve fit
    !
    !   Author:
    !   Ipseng Iu
    !   Mechanical and Aerospace Engineering
    !   Oklahoma State University, Stillwater
    !
    !   Date: July 2005
    !
    ! ----------------------------------------------------------------------

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)

    IMPLICIT NONE

    REAL, INTENT(IN) :: refrigerantMassFlowRate !Refrigerant mass flow rate, kg/s
    REAL, INTENT(IN) :: evaporatingTemperature  !Evaporating temperature, C
    REAL, INTENT(IN) :: condensingTemperature   !Condensing temperature, C
    REAL, INTENT(IN) :: Subcooling              !Subcooling, K
    REAL, INTENT(IN) :: Superheat               !Superheat, K
    REAL, INTENT(INOUT) :: liquidLineQuality    !Liquid line quality
    REAL, INTENT(INOUT) :: suctionLineQuality   !Suction line quality
    REAL, INTENT(OUT) :: estimatedPressureDrop  !Estimated pressure drop, kPa

    INTEGER ErrorFlag
    !0-No error
    !1-Accumulator solution not converge
    !2-Refprop error

    REAL MaxCapacity !Max. system capacity, ton
    REAL systemCapacityInTon     !System capacity, ton
    REAL estimatedPressureDrop_Temp       !Estimated temperature drop, kPa
    REAL suctionPressure        !Suction presure, kPa
    REAL dischargePressure         !Discharge pressure, kPa
    REAL liquidEnthalpy         !Liquid enthalpy, kJ/kg
    REAL vaporEnthalpy         !Vapor enthalpy, kJ/kg
    REAL saturationPressure_1, saturationPressure_2!Saturation pressure, kPa
    IF (liquidLineQuality .GT. 1) THEN
      liquidLineQuality=1
    END IF
    IF (liquidLineQuality .LT. 0) THEN
      liquidLineQuality=0
    END IF
    IF (suctionLineQuality .GT. 1) THEN
      suctionLineQuality=1
    END IF
    IF (suctionLineQuality .LT. 0) THEN
      suctionLineQuality=0
    END IF

    ErrorFlag=0

    Temperature=evaporatingTemperature
    Quality=1
    suctionPressure=TQ(RefName,Temperature,Quality,'pressure',RefrigIndex,RefPropErr)  !Suction Pressure
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
      RETURN
    END IF
    suctionPressure=suctionPressure/1000  !RS Comment: Unit Conversion

    Temperature=condensingTemperature
    Quality=1
    dischargePressure=TQ(RefName,Temperature,Quality,'pressure',RefrigIndex,RefPropErr)  !Discharge Pressure
    IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
      RETURN
    END IF
    dischargePressure=dischargePressure/1000  !RS Comment: Unit Conversion

    IF (Superheat .GT. 0) THEN
      Pressure=suctionPressure*1000  !RS Comment: Unit Conversion
      Temperature=evaporatingTemperature+Superheat
      vaporEnthalpy=TP(RefName,Temperature,Pressure,'enthalpy',RefrigIndex,RefPropErr) !Vapor Enthalpy
      IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
        RETURN
      END IF
      vaporEnthalpy=vaporEnthalpy/1000  !RS Comment: Unit Conversion
    ELSE
      Pressure=suctionPressure*1000  !RS Comment: Unit Conversion
      Quality=suctionLineQuality
      vaporEnthalpy=PQ(RefName,Pressure,Quality,'enthalpy',RefrigIndex,RefPropErr) !Vapor Enthalpy
      IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
        RETURN
      END IF
      vaporEnthalpy=vaporEnthalpy/1000  !RS Comment: Unit Conversion
    END IF

    IF (Subcooling .GT. 0) THEN
      Pressure=dischargePressure*1000  !RS Comment: Unit Conversion
      Temperature=condensingTemperature-Subcooling
      liquidEnthalpy=TP(RefName,Temperature,Pressure,'enthalpy',RefrigIndex,RefPropErr) !Liquid Enthalpy
      IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
        RETURN
      END IF
      liquidEnthalpy=liquidEnthalpy/1000  !RS Comment: Unit Conversion
    ELSE
      Pressure=dischargePressure*1000  !RS Comment: Unit Conversion
      Quality=liquidLineQuality
      liquidEnthalpy=PQ(RefName,Pressure,Quality,'enthalpy',RefrigIndex,RefPropErr) !Liquid Enthalpy
      IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
        RETURN
      END IF
      liquidEnthalpy=liquidEnthalpy/1000  !RS Comment: Unit Conversion
    END IF

    systemCapacityInTon=refrigerantMassFlowRate*(vaporEnthalpy-liquidEnthalpy)*0.28435 !Convert from kW to ton

    IF (CoefficientOfCurveFit_M .NE. 0 .OR. CoefficientOfCurveFit_B .NE. 0) THEN
      MaxCapacity = CoefficientOfCurveFit_M * evaporatingTemperature + CoefficientOfCurveFit_B

      IF (RatedPressureDrop .NE. 0) THEN
        estimatedPressureDrop = systemCapacityInTon / MaxCapacity * RatedPressureDrop
      ELSE
        estimatedPressureDrop_Temp = systemCapacityInTon / MaxCapacity * RatedTemperatureDrop
        saturationPressure_1 = suctionPressure
        Temperature=evaporatingTemperature-estimatedPressureDrop_Temp
        Quality=1
        saturationPressure_2=TQ(RefName,Temperature,Quality,'pressure',RefrigIndex,RefPropErr) !Saturation Pressure 2
        IF (IssueRefPropError(RefPropErr, 'Accumulator', 2, ErrorFlag)) THEN
          RETURN
        END IF
        saturationPressure_2=saturationPressure_2/1000    !RS Comment: Unit Conversion
        RatedPressureDrop = saturationPressure_1 - saturationPressure_2 !Rated Pressure Drop
        estimatedPressureDrop = systemCapacityInTon / MaxCapacity * RatedPressureDrop
      END IF

    END IF

    PressureDropInAccumulator=estimatedPressureDrop

    RETURN

  END SUBROUTINE CalcAccumulatorDP

  !******************************************************************************

  SUBROUTINE InitAccumulator !(PAR)

    ! ----------------------------------------------------------------------
    !
    !   Purpose: To initialize accumulator parameters
    !
    !   Method: Load the parameters to public variables
    !
    !   Inputs:
    !       Ref$ = Refrigerant name
    !       PAR(1) = Accumulator inside diameter, m
    !       PAR(2) = Accumulator internal height, m
    !       PAR(3) = J-tube lower hole diameter, m
    !       PAR(4) = J-tube upper hole diameter, m
    !       PAR(5) = Distance between holes on J-tube, m
    !       PAR(6) = J-tube inside diameter, m
    !       PAR(7) = Rating pressure drop, kPa
    !       PAR(8) = Rating temperature drop, K
    !       PAR(9) = Curve fit coefficient M
    !       PAR(10) = Curve fit coefficient B

    !
    !   Author:
    !   Ipseng Iu
    !   Mechanical and Aerospace Engineering
    !   Oklahoma State University, Stillwater
    !
    !   Date: July 2005
    !
    ! ----------------------------------------------------------------------

    IMPLICIT NONE

    !REAL, INTENT(IN) :: PAR(10)

    !Flow:

    !Convert from m to ft
    accumulatorDiameter = AccumPAR%AccD/0.3048 !ACCDIA    !RS: Debugging: Formerly PAR(1)
    accumulatorHeight = AccumPAR%AccH/0.3048 !ACCHGT    !RS: Debugging: Formerly PAR(2)
    DiameterHole(1)=AccumPAR%AccD1/0.3048  !RS: Debugging: Formerly PAR(3)
    DiameterHole(2)=AccumPAR%AccD2/0.3048  !RS: Debugging: Formerly PAR(4)
    distanceBetweenHoles = AccumPAR%AccHDis/0.3048 !HOLDIS    !RS: Debugging: Formerly PAR(5)
    TubeDiameter = AccumPAR%AccDTube/0.3048 !ATBDIA   !RS: Debugging: Formerly PAR(6)

    RatedPressureDrop=AccumPAR%AccDP  !RS: Debugging: Formerly PAR(7)
    RatedTemperatureDrop=AccumPAR%AccDT  !RS: Debugging: Formerly PAR(8)
    CoefficientOfCurveFit_M=AccumPAR%AccCM   !RS: Debugging: Formerly PAR(9)
    CoefficientOfCurveFit_B=AccumPAR%AccCB  !RS: Debugging: Formerly PAR(10)

    RETURN
  END SUBROUTINE InitAccumulator
END MODULE
