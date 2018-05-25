! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This is the overarching heat pump design subroutine. It calls and manages the component modules. 

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! Since this subroutine runs the heat pump simulation, there's really no physical component.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! The input, DTVALU, appears to be a guess of temperature. There don't appear to be any direct outputs.
! This module calls the Evaporator subroutine for the first heat pump run case.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no direct input or output files connected to this routine.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! This is only a subroutine, so there are no module-level variables or structures.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method
!    PUBLIC HPDM -- Runs the entire heat pump simulation
!      Called by ChargeLoop.f90
!      Called by ORNLsolver.f90

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
! Some more documentation would probably be useful.

    SUBROUTINE HPDM(DTVALU)

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE CondenserMod
    USE EvaporatorMod
    USE ShortTubeMod
    USE CapillaryTubeMod
    USE TXVMOD
    USE AccumulatorModule
    USE DataSimulation
    USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
    USE InputProcessor    !RS: Debugging: Bringing over from GetHPSimInputs

    IMPLICIT NONE

    REAL Temperature,Quality,Pressure,Enthalpy

    REAL DTVALU

    LOGICAL PRINT

    INTEGER(2) RefPropErr			!Error flag:1-error; 0-no error

    INTEGER(2) AirPropOpt			!Air prop calc. option
    INTEGER(2) AirPropErr			!Error flag:1-error; 0-no error
    !REAL AirProp(8)		!Air properties

    INTEGER ICHRGE,IMASS,IREFC,LPRINT
    REAL TAIIEI
    INTEGER NTAMB,NCROSS
    REAL DELT2,DTVLMN
    REAL TAIIE
    INTEGER I
    REAL ERRMSG(2)
    REAL TSAT1,CONV,STEP,DIFFER,XMR
    INTEGER IERROR,IER
    REAL TAIE1,DIFF,DIFSGN,PROD,TSATSV,TSATDM,TAISV,TAIDM
    REAL TsoEvp,LsucLn
    REAL MassCoil,MassLiqCoil,MassVapCoil
    INTEGER NumIter,MaxIteration
    REAL XMRFLD,ErrXMR,TSICMPprev
    REAL Dshtb,MaxDshTb,MinDshTb
    REAL CapTubeDimension,MaxLen,MinLen
    REAL Qtxv
    REAL Subcooling, Superheat, DPtxv
    REAL ChargeCorrection !Correction charge for the charge tuning method, lbm
    REAL, EXTERNAL :: CNDNSR, EVPTR
    REAL ZeroConvergence
    REAL, PARAMETER :: Dstep=1
    REAL, PARAMETER :: CapTubeDimStep=1E-3

    LOGICAL IsSizeDiameter
    REAL DetailedQevp,DetailedDPevp
    REAL SimpleQevp,SimpleDPevp
    LOGICAL,SAVE :: IsFirstTimeEvaporator = .TRUE. !First time to call evaporator flag
    INTEGER ChargeOption !Charge option, 1=no tuning; 2=w/charge tuning

    LOGICAL :: FLAG_GOTO_950

    CHARACTER(LEN=200) :: tmpString

    LOGICAL, EXTERNAL :: IssueRefPropError
    
    INTEGER :: TimeStep1 !RS: Testing

    
    INTEGER, PARAMETER :: MaxNameLength = 200

    CHARACTER(len=MaxNameLength),DIMENSION(200) :: Alphas ! Reads string value from input file
    INTEGER :: NumAlphas               ! States which alpha value to read from a "Number" line
    REAL, DIMENSION(500) :: Numbers    ! brings in data from IP
    INTEGER :: NumNumbers              ! States which number value to read from a "Numbers" line
    INTEGER :: Status                  ! Either 1 "object found" or -1 "not found"
    INTEGER, PARAMETER :: r64=KIND(1.0D0)  !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12) 
    REAL(r64), DIMENSION(500) :: TmpNumbers !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    REAL RefSimulatedCharge     !Simulated charge at reference point, kg or lbm
    REAL SimulatedCharge2       !Simulated charge at 2nd reference point, kg or lbm
    REAL LiquidLength2          !Liquid length at 2nd reference point, m or ft
  
        
    TimeStep1=0 !Karthik - Initialize to 0
    TAISV=1.0 !Karthik - Initialize to 1
    TSATSV=1.0 !Karthik - Initialize to 1
    
    IF (EvapPAR%EvapFirstTime .EQ. 1) THEN    !RS: Debugging: Formerly EvapPAR(38)
          !*************** Charge Tuning Curve ***************  !RS: Debugging: Moving: HPDM

    CALL GetObjectItem('ChargeTuningCurve',1,Alphas,NumAlphas, &
                        TmpNumbers,NumNumbers,Status) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)     

        Numbers = DBLE(TmpNumbers) !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
  
        SELECT CASE (Alphas(1)(1:1))  !Is Charge Tuning?
            CASE ('F','f')
                IsChargeTuning=0  !RS: Debugging: If this is the case, I don't think these inputs are ever used
            CASE ('T','t')
                IsChargeTuning=1
        END SELECT
  
        RefSimulatedCharge = Numbers(1)   !Tuning Point #1 Simulated Charge
        RefLiquidLength = Numbers(2)  !Tuning Point #1 Liquid Length
        SimulatedCharge2 = Numbers(3) !Tuning Point #2 Simulated Charge
        LiquidLength2 = Numbers(4)    !Tuning Points #2 Liquid Length
  
        !store the refrigerant name in data globals
        RefName = Ref$
  
        !Calculate charge tuning curve
        IF (MODE .NE. 2 .AND. (RefLiquidLength-LiquidLength2) .NE. 0) THEN
	        IF (RefChg .GT. 0) THEN
		        ChargeCurveSlope=(SimulatedCharge2-RefSimulatedCharge)/ &
						        (LiquidLength2-RefLiquidLength)
		        ChargeCurveIntercept=RefChg-RefSimulatedCharge
	        ELSE
		        ChargeCurveSlope=0
		        ChargeCurveIntercept=0
	        END IF
        END IF
        ChargeCurveSlope=ChargeCurveSlope*UnitM/UnitL
	    ChargeCurveIntercept=ChargeCurveIntercept*UnitM
	    RefLiquidLength=RefLiquidLength*UnitL
  
    END IF

    MaxIteration=30
    ICHRGE=1
    IMASS=1
    LPRINT=1

    !for specified subcooling, set IREFC to zero
    !for specifed flow control, set IREFC to 3 
    IF (MODE .EQ. 1 .OR. MODE .EQ. 3) THEN
        IREFC=3
    ELSE
        IREFC=0
    END IF
    
    TAIIEI=TaiE

    NTAMB = 0
    NCROSS = 0
    DELT2 = 1.25

    DTVLMN = -150.0

    IF(DTVALU.LT.DTVLMN) THEN
        DTVALU = DTVLMN
    END IF

    IF(ICHRGE.NE.0) THEN
        IF(ICHRGE.NE.2) THEN
            IF (MODE .EQ. 4) THEN
                DTROC = DTVALU
                IF(DTROC.LT.0.0) THEN
                    DTROC = DTROC/200.
                END IF
            ELSE
                SUPER = DTVALU
                IF(SUPER.LE.0.0) THEN
                    SUPER = -(1.0+SUPER/500.)
                END IF
            END IF

            Temperature=(TSICMP-32)/1.8 !RS Comment: Unit Conversion, from F to C
            Quality=1
            PiCmp=TQ(Ref$,Temperature,Quality,'pressure',RefrigIndex,RefPropErr)    !Compressor Inlet Pressure
            IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                STOP
            END IF
            PiCmp=PiCmp/1000    !RS Comment: Unit Conversion

            IF (SUPER .GT. 0) THEN
                Temperature=(TSICMP+SUPER-32)/1.8   !RS Comment: Unit Conversion, from F to C
                Pressure=PiCmp*1000 !RS Comment: Unit Conversion
                HiCmp=TP(Ref$,Temperature,Pressure,'enthalpy',RefrigIndex,RefPropErr)   !Compressor Inlet Enthalpy
                IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                    STOP
                END IF
                HiCmp=HiCmp/1000    !RS Comment: Unit Conversion
                TiCmp=((TSICMP+SUPER)-32)/1.8   !RS Comment: Unit Conversion, from F to C
            ELSE
                Pressure=PiCmp*1000 !RS Comment: Unit Conversion
                Quality=-SUPER
                HiCmp=PQ(Ref$,Pressure,Quality,'enthalpy',RefrigIndex,RefPropErr)   !Compressor Inlet Enthalpy
                IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                    STOP
                END IF
                HiCmp=HiCmp/1000    !RS Comment: Unit Conversion
                TiCmp=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)    !Compressor Inlet Temperature
                IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                    STOP
                END IF
            END IF

            Xmr=CompOUT%CmpOMdot  !RS: Debugging: Formerly CompOUT(2)

            PoEvp=EvapOUT%EOutpRoC    !RS: Debugging: Formerly EvapOUT(1)

            QsucLn=EvapPAR%EvapSucLnQLoss   !RS: Debugging: Formerly EvapPAR(5)
            DTsucLn=EvapPAR%EvapSucLnTempChg  !RS: Debugging: Formerly EvapPAR(6)
            LsucLn=EvapPAR%EvapSucLnLen   !RS: Debugging: Formerly EvapPAR(1)

            IF (LsucLn .GT. 0) THEN
                IF (QsucLn .NE. 0) THEN
                    HoEvp=HiCmp-QsucLn/Xmr

                    Pressure=PoEvp*1000 !RS Comment: Unit Conversion
                    Enthalpy=HoEvp*1000 !RS Comment: Unit Conversion
                    ToEvp=PH(Ref$,Pressure,Enthalpy,'temperature',RefrigIndex,RefPropErr)   !Evaporator Outlet Temperature
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF

                    XoEvp=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)   !Evaporator Outlet Quality
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF

                    Quality=1
                    TsoEvp=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Evaporator Outlet Saturation Temperature
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF

                    SUPERE=(ToEvp-TsoEvp)*1.8

                    IF (XoEvp .LT. 1.) THEN
                        SUPERE = -XoEvp
                    END IF

                ELSEIF (DTsucLn .NE. 0) THEN

                    ToEvp=TiCmp-DTsucLn

                    Temperature=ToEvp
                    Pressure=PoEvp*1000 !RS Comment: Unit Conversion
                    HoEvp=TP(Ref$, Temperature, Pressure, 'enthalpy', RefrigIndex,RefPropErr)   !Evaporator Outlet Enthalpy
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF
                    HoEvp=HoEvp/1000    !RS Comment: Unit Conversion

                    Pressure=PoEvp*1000 !RS Comment: Unit Conversion
                    Enthalpy=HoEvp*1000 !RS Comment: Unit Conversion
                    XoEvp=PH(Ref$,Pressure,Enthalpy,'quality',RefrigIndex,RefPropErr)   !Evaporator Outlet Quality
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF

                    Quality=1
                    TsoEvp=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)   !Evaporator Outlet Saturation Temperature
                    IF (IssueRefPropError(RefPropErr, 'HPdesign')) THEN
                        STOP
                    END IF

                    SUPERE=(ToEvp-TsoEvp)*1.8

                    IF (XoEvp .LT. 1.) THEN
                        SUPERE = -XoEvp
                    END IF

                ELSE
                    SUPERE=SUPER
                END IF

            ELSE
                SUPERE=SUPER
            END IF
        ELSE
            DTROC = DTVALU
            IF(DTROC.LT.0.0) THEN
                DTROC = DTROC/200.
            END IF
        END IF

    END IF

    FLAG_GOTO_950 = .FALSE.
    DO WHILE (.TRUE.)
        
        TimeStep1 = TimeStep1+1   !RS: Testing

        CurSimTime=(TimeStep1-1)*TimeInterval  !PrevSimTime+ !RS: Testing
        
        DO I=1,2
            ERRMSG(I) = 0.0
        END DO

        !   FIND DESIRED CONDENSER SUBCOOLING
        !        OR REFRIGERANT MASS FLOW RATE BALANCE
        !   BY ADJUSTING COMPRESSOR EXIT SATURATION TEMPERATURE

        !       USE 'ZERO3' AND 'CNDNSR' TO FIND TSOCMP SUCH THAT
        !       ABS(CDTROC - DTROC)<CNDCON --  IF IREFC = 0 OR
        !       ABS(XMRFLD - XMR)<FLOCON/20 -- IF IREFC IS NOT EQUAL TO 0

        TSAT1 = TSOCMP

        CONV = CNDCON
        IF(IREFC .NE. 0) THEN
            CONV = FLOCON !/20.
        END IF
        STEP = 3

        !CALL IssueOutputMessage('|-------------------- Highside Iteration --------------------|')

        AirPropOpt=2
        AirProp%APTDB=(TaiC-32)*5/9    !RS Comment: Unit Conversion, from F to C   !RS: Debugging: Formerly AirProp(1)
        AirProp%APRelHum=RHiC !RS: Debugging: Formerly AirProp(3)
        CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
        RhoAiC=AirProp%APDryDens   !RS: Debugging: Formerly AirProp(7)

        AirPropOpt=2
        AirProp%APTDB=(TaiE-32)*5/9    !RS Comment: Unit Conversion, from F to C   !RS: Debugging: Formerly AirProp(1)
        AirProp%APRelHum=RHiE !RS: Debugging: Formerly AirProp(3)
        CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
        RhoAiE=AirProp%APDryDens   !RS: Debugging: Formerly AirProp(7)

        !Actual air mass flow rate
        XMaC=CFMcnd*RhoAiC
        XMaE=CFMevp*RhoAiE  !RS: Debugging: Uncommented since we're trying to use the optimized flow rate (10/18/14)

        CoilMode=0  !RS: Debugging: This is for a test in ZeroConvergence (11/18/13)
        TSOCMP = ZeroConvergence(TSAT1,CNDNSR,1E-3,CNDCON,STEP,DIFFER,IERROR)

        IF (IERROR .GE. 3) THEN
            CALL IssueOutputMessage('')
            CALL IssueOutputMessage('## ERROR ## Highside: Failed to find a solution.')
            CALL IssueOutputMessage('Try another condenser, compressor, or change boundary conditions.')
            !STOP  !RS: Debugging: Temporarily setting in an Epsilon-NTU method
        END IF
        FirstTimeFlowRateLoop=.FALSE.
        !ISI 05-25-05
        IF (MODE .EQ. 5) THEN !Peform the rest of the calculation for Condenser Unit
            FLAG_GOTO_950 = .TRUE.
            EXIT
        END IF

        IF (LPRINT.EQ.2) THEN 
            PRINT = .TRUE.
            DIFFER = CNDNSR(TSOCMP,IER)
        END IF
        IF (ABS(DIFFER) .GT. CONV) THEN
            IF (LPRINT .GT. 1) THEN
                IF (Unit .EQ. 1) THEN
                    CALL IssueOutputMessage('## ERROR ## Highside: Solution not converged on subcooling.')
                    WRITE(tmpString,'(F10.3)') DIFFER/1.8
                    CALL IssueOutputMessage('Difference: '//TRIM(tmpString)//DTunit)
                ELSE
                    CALL IssueOutputMessage('## ERROR ## Highside: Solution not converged on subcooling.')
                    WRITE(tmpString,'(F10.3)') DIFFER
                    CALL IssueOutputMessage('Difference: '//TRIM(tmpString)//DTunit)
                END IF  
            END IF

            ERRMSG(1) = DIFFER
        END IF
        
        EvapIN%EInmRef=MdotR           !Refrigerant side mass flow rate, kg/s    !RS: Debugging: Formerly EvapIN(1)
        EvapIN%EInhRi=CondOUT%COuthRiE     !Exp. device inlet enthalpy, kJ/kg    !RS: Debugging: Formerly EvapIN(3), CondOUT(11)
        EvapIN%EInmAi=XMaE            !Air side mass flow rate, kg/s    !RS: Debugging: Formerly EvapIN(4)
        EvapIN%EIntAi=(TAIIEI-32)/1.8 !Air side inlet temp. C   !RS: Debugging: Formerly EvapIN(5)
        EvapIN%EInrhAi=RHiE            !Air side inlet relative humidity !RS: Debugging: Formerly EvapIN(6)
        EvapIN%EIntRdis=CompOUT%CmpOTdis      !Discharge temperature, C !RS: Debugging: Formerly EvapIN(9), CompOUT(5)

        !Take compressor shell loss into account
        IF (CompPAR%CompQLossFrac .NE. 0) THEN !Shell loss in fraction    !RS: Debugging: Formerly CompPAR(21)
            EvapPAR%EvapCompQLoss=CompPAR%CompQLossFrac*CompOUT%CmpOPwr  !RS: Debugging: Formerly EvapPAR(32), CompPAR(21), CompOUT(1)
        ELSE !Shell loss in W
            EvapPAR%EvapCompQLoss=CompPAR%CompQLoss/1000    !RS: Debugging: Formerly EvapPAR(32) & CompPAR(22)
        END IF

        IF (FirstTimeHPdesignMode) THEN

            IF ((IsCoolingMode .GT. 0 .AND. IDCcoilType .EQ. MCEVAPORATOR) .OR. &
            (IsCoolingMode .LT. 1 .AND. ODCcoilType .EQ. MCEVAPORATOR)) THEN
                !Microchannel coil
                EvapPAR%EvapFirstTime=1 !First time   !RS: Debugging: Formerly EvapPAR(38)
                EvapPAR%EvapSimpCoil=0 !Detailed version !RS: Debugging: Formerly EvapPAR(37)
                CALL Evaporator(Ref$) !,EvapIN,EvapPAR,EvapOUT) !(Ref$,PureRef,EvapIN,EvapPAR,EvapOUT) !RS: Debugging: Extraneous PureRef
                EvapPAR%EvapFirstTime=0 !No longer first time !RS: Debugging: Formerly EvapPAR(38)
            ELSE
                !Plate-fin coil
                !Run both simple and detailed version to determine which one to use
                !Change the logic to reset IsFirstTimeEvaporator
                IF (IsFirstTimeEvaporator) THEN
                    EvapPAR%EvapFirstTime=1 !First time   !RS: Debugging: Formerly EvapPAR(38)
                    EvapPAR%EvapSimpCoil=0 !Detailed version !RS: Debugging: Formerly EvapPAR(37)
                    CALL Evaporator(Ref$) !EvapIN,EvapPAR,DetailedEvapOUT) !(Ref$,PureRef,EvapIN,EvapPAR,DetailedEvapOUT) !RS: Debugging: Extraneous PureRef
                    DetailedQevp=EvapOUT%EOutQC    !RS: Debugging: Formerly DetailedEvapOUT(11)
                    DetailedDPevp=EvapIN%EInpRi-EvapOUT%EOutpRiC  !RS: Debugging: Formerly EvapIN(2), DetailedEvapOUT(6)

                    EvapPAR%EvapSimpCoil=1 !Simple version   !RS: Debugging: Formerly EvapPAR(37)
                    CALL Evaporator(Ref$) !,EvapIN,EvapPAR,SimpleEvapOUT) !(Ref$,PureRef,EvapIN,EvapPAR,SimpleEvapOUT)   !RS: Debugging: Extraneous PureRef
                    SimpleQevp=EvapOUT%EOutQC    !RS: Debugging: Formerly SimpleEvapOUT(11)
                    SimpleDPevp=EvapIN%EInpRi-EvapOUT%EOutpRiC !RS: Debugging: Formerly EvapIN(2), SimpleEvapOUT(6)

                    IF (ABS((SimpleQevp-DetailedQevp)/DetailedQevp) .LT. 0.1 .AND. &
                    ABS((SimpleDPevp-DetailedDPevp)/DetailedDPevp) .LT. 0.1) THEN
                        EvapPAR%EvapSimpCoil=1 !Simple version   !RS: Debugging: Formerly EvapPAR(37)
                    ELSE
                        EvapPAR%EvapSimpCoil=0 !Detailed version !RS: Debugging: Formerly EvapPAR(37)
                    END IF
                    IsFirstTimeEvaporator=.FALSE. 

                    EvapPAR%EvapSimpCoil=1  !RS: Debugging: Always simple (6/23/14)
                    !Always detailed    !RS: Debugging: There's no need for this to be set
                    !EvapPAR%EvapSimpCoil=0 !Detailed version !RS: Debugging: Formerly EvapPAR(53) !RS: Debugging: Simple case 

                ELSE
                    CALL Evaporator(Ref$) !,EvapIN,EvapPAR,EvapOUT) !(Ref$,PureRef,EvapIN,EvapPAR,EvapOUT) !RS: Debugging: Extraneous PureRef
                    EvapPAR%EvapFirstTime=1 !0 !No longer first time !RS: Debugging: Formerly EvapPAR(38)
                END IF
            END IF

            IF (EvapOUT%EOutErrFlag .NE. 0) THEN    !RS: Debugging: Formerly EvapOUT(17)
                SELECT CASE (INT(EvapOUT%EOutErrFlag))  !RS: Debugging: Formerly EvapOUT(17)
                CASE (3,4,5)
                    STOP
                END SELECT
            END IF
            FirstTimeHPdesignMode=.FALSE.

        END IF

        IF (LPRINT .EQ. 2) THEN
            PRINT = .FALSE.
        END IF

        !   FIND DESIRED EVAPORATOR EXIT SUPERHEAT OR QUALITY
        !   BY ADJUSTING EVAPORATOR INLET AIR TEMPERATURE

        TAIE1 = TAIIEI

        STEP = 2

        !CALL IssueOutputMessage(' ')
        !CALL IssueOutputMessage('|-------------------- Lowside Iteration ---------------------|')
        IF (Unit .EQ. 1) THEN
            WRITE(tmpString, '(F10.4)') (TSICMP-32)*5/9
        ELSE
            WRITE(tmpString, '(F10.4)') TSICMP
        END IF
        !CALL IssueOutputMessage( '>> Compressor suction saturation temperature: '//TRIM(tmpString)//Tunit)

        CoilMode=1  !RS: Debugging: This is for a test in ZeroConvergence (11/18/13)
        TAIIE = ZeroConvergence(TAIE1,EVPTR,AMBCON,EVPCON,STEP,DIFFER,IERROR)

        IF (IERROR .GE. 3) THEN
            CALL IssueOutputMessage('')
            CALL IssueOutputMessage('## ERROR ## Lowside: Failed to find a solution.')
        END IF

        IF (LPRINT .GT. 2) THEN
            PRINT = .TRUE.
            DIFFER = EVPTR(TAIIE,IER)
            PRINT = .FALSE.
        END IF
        
        IF (ABS(DIFFER) .GT. EVPCON) THEN
            IF (LPRINT .GT. 1) THEN
                IF (Unit .EQ. 1) THEN
                    CALL IssueOutputMessage('## ERROR ## Lowside: Solution not converged on superheat.')
                    WRITE(tmpString,'(F10.4)') DIFFER/1.8
                    CALL IssueOutputMessage('Difference: '//tmpString//DTunit)
                ELSE
                    CALL IssueOutputMessage('## ERROR ## Lowside: Solution not converged on superheat.')
                    WRITE(tmpString,'(F10.4)') DIFFER
                    CALL IssueOutputMessage('Difference: '//tmpString//DTunit)
                END IF  
            END IF
            ERRMSG(2) = DIFFER
        END IF

        IF(LPRINT.GT.1.AND.IMASS.NE.0) THEN
            IF (AccumPAR%AccH .GT. 0) THEN !Height    !RS: Debugging: Formerly AccumPAR(2)
                AccumIN%AccImdot= MdotR    !RS: Debugging: Formerly AccumIN(1)
                AccumIN%AccIpRo=CompIN%CompInPsuc !Pressure  !RS: Debugging: Formerly CompIN(1), AccumIN(2)
                AccumIN%AccIhRo=CompIN%CompInHsuc !Enthalpy  !RS: Debugging: Formerly CompIN(3), AccumIN(3)
                CALL CalcAccumulatorMass !(AccumIN,AccumOUT)
            ELSE
                AccumOUT%AccOMass=0
            END IF

            CALL CalcCondenserInventory(MassCoil,MassLiqCoil,MassVapCoil,CondLiqTubeLength,CondVapTubeLength,CondTwoPhaseTubeLength,CondNumLiqTubes)
            CondOUT%COutMC=MassCoil    !RS: Debugging: Formerly CondOUT(18)
            CALL CalcEvaporatorInventory(MassCoil,MassLiqCoil,MassVapCoil,EvapLiqTubeLength,EvapVapTubeLength,EvapTwoPhaseTubeLength,EvapNumLiqTubes)
            EvapOUT%EOutMC=MassCoil    !RS: Debugging: Formerly EvapOUT(14)

            IF (ExpDevice .EQ. 1) THEN
                CALCHG=(CompOUT%CmpOMCmp+CondOUT%COutMDisLn+CondOUT%COutMLiqLn+CondOUT%COutMC+ &   !RS: Debugging: Formerly CompOUT(6), CondOUT(16), CondOUT(17), CondOUT(18)
                EvapOUT%EOutMSucLn+EvapOUT%EOutMC+ShTbOUT%ShTbOMDT+AccumOUT%AccOMass)/UnitM   !RS: Debugging: Formerly EvapOUT(13), EvapOUT(14), ShTbOUT(5), AccumOUT(1)
            ELSE
                CALCHG=(CompOUT%CmpOMCmp+CondOUT%COutMDisLn+CondOUT%COutMLiqLn+CondOUT%COutMC+ &   !RS: Debugging: Formerly CompOUT(6), CondOUT(16), CondOUT(17), CondOUT(18)
                EvapOUT%EOutMSucLn+EvapOUT%EOutMC+AccumOUT%AccOMass)/UnitM    !RS: Debugging: Formerly EvapOUT(13)+EvapOUT(14)+TxvOUT(5)+AccumOUT(1))/UnitM
            END IF
        END IF

        !   FIND DESIRED EVAPORATOR INLET AIR TEMPERATURE
        !   BY ADJUSTING COMPRESSOR INLET SATURATION TEMPERATURE

        DIFF = TAIIE-TAIIEI
        !VL: Previously: IF(ABS(DIFF).LE.AMBCON) GO TO 900
        IF(ABS(DIFF).LE.AMBCON) THEN
            EXIT
        END IF

        IF(NTAMB.EQ.0) THEN
            DIFSGN = DIFF
        END IF

        PROD = DIFF*DIFSGN

        IF(PROD.GT.0.0.AND.NCROSS.EQ.0) THEN

            IF(NTAMB.GT.0) THEN
                DELT2 = (TAIIE-TAIDM)/(TSICMP-TSATDM)
            END IF
            IF (ABS(DELT2) .LE. 0.05) THEN !ISI - 06/13/07
                !VL: Previously: GO TO 900 !0.05 F !ISI - 08/02/06
                EXIT
            END IF

            TSATDM = TSICMP
            TAIDM = TAIIE
            TSICMPprev=TSICMP
            TSICMP = TSICMP-DIFF/DELT2
            IF (TSICMP .GT. TAIIEI) THEN
                TSICMP=(TSICMPprev+TAIIEI)/2 !Make sure TSICMP < TAIIEI
            END IF

        ELSE    !RS: The following runs after DIFF goes negative (because it sets PROD negative) (12/19/13)

            NCROSS = 1  !RS: Ensures that this section will run until the end of the loop (12/19/13)
            
            IF(PROD.LE.0.0) THEN
                TSATSV = TSATDM
                TAISV = TAIDM
            END IF

            TSATDM = TSICMP
            TAIDM = TAIIE
            TSICMPprev=TSICMP
            TSICMP = TSICMP-(TSATSV-TSICMP)/(TAISV-TAIIE)*DIFF
            IF (ABS(TSICMPprev-TSICMP) .LE. 0.01) THEN
                !VL: Previously: GO TO 900 !0.05 F !ISI - 08/02/06
                EXIT !RS: Debugging: Seeing if it'll converge on the air inlet temp (11/13/13)
            END IF
            DIFSGN = DIFF
            IF (TSICMP .GT. TAIIEI) THEN
                TSICMP=(TSICMPprev+TAIIEI)/2 !Make sure TSICMP < TAIIEI
            END IF

        END IF

        NTAMB = NTAMB + 1
        IF(NTAMB.GT.15) THEN
            WRITE(tmpString,"('0DRIVER: ***** FAILED TO CONVERGE ON EVAPORATOR ',  'INLET AIR TEMPERATURE *****',/, '               DIFFERENCE  =',F8.3,' F')") DIFF
            CALL IssueOutputMessage(TRIM(tmpString))
            EXIT
        END IF
        IF (LPRINT .GT. 1) THEN
            WRITE(tmpString,"('0        DID NOT CONVERGE ON  EVAPORATOR INLET ',  'AIR TEMPERATURE FOR THIS SATURATION TEMPERATURE.' ,/,'         SET COMPRESSOR INLET SATURATION TEMPERATURE TO',  F8.3,' F AND GO BACK TO CONDENSER ITERATION.')")TSICMP
            CALL IssueOutputMessage(TRIM(tmpString))
        END IF

        FirstTimeAirTempLoop=.TRUE.

        IF (TSICMP .GE. TSOCMP) THEN
            CALL IssueOutputMessage('')
            CALL IssueOutputMessage('## ERROR ## HPdesign: Failed to find a solution.')
            STOP
        END IF

    END DO

    IF (FLAG_GOTO_950 .EQ. .FALSE.) THEN 

        IF (IREFC .EQ. 0) THEN
            !**************Size short tube orifice**************

            XMR=CompOUT%CmpOMdot*3600/UnitM   !RS: Debugging: Formerly CompOUT(2)

            ShTbIN%ShTbINMdotC=CompOUT%CmpOMdot  !Compressor mass flow rate, kg/s  !RS: Debugging: Formerly CompOUT(2), ShTbIN(1)
            ShTbIN%ShTbINPiE=CondOUT%COutpRiE !Exp. device inlet pressure, kPa  !RS: Debugging: Formerly CondOUT(10), ShTbIN(2)
            ShTbIN%ShTbINHiE=CondOUT%COuthRiE !Exp. device inlet enthalpy, kJ/kg    !RS: Debugging: Formerly CondOUT(11), ShTbIN(3)
            ShTbIN%ShTbINPiEv=EvapIN%EInpRi   !Evaporator inlet pressure, kPa   !RS: Debugging: Formerly EvapIN(2), ShTbIN(4)
            ShTbIN%ShTbINPoEv=EvapOUT%EOutpRoC  !Evaporator outlet pressure, kPa  !RS: Debugging: Formerly EvapOUT(1), ShTbIN(5)

            IF (ShTbPAR%ShTbTLen .LE. 0) THEN !RS: Debugging: Formerly ShTbPAR(1)
                ShTbPAR%ShTbTLen=0.0127   !RS: Debugging: Formerly ShTbPAR(1) ...and a Magic Number
                !Short Tube: Parameters not defined.
            ELSE

                !Initial guess
                NumIter=0
                MaxDshTb=0
                MinDshTb=0

                Dshtb=2.0 !1.0 !Initial guess !Short tube diameter, mm
                ShTbPAR%ShTbTID=Dshtb/1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly ShTbPAR(2)

                DO NumIter=1, MaxIteration

                    !CALL ShortTube(Ref$,PureRef,ShTbIN,ShTbPAR,ShTbOUT)
                    !CALL ShortTubePayne(Ref$,PureRef,ShTbIN,ShTbPAR,ShTbOUT)
                    CALL ShortTubePayne(Ref$) !,ShTbIN,ShTbPAR,ShTbOUT)   !RS: Debugging: Extraneous PureRef
                    IF (ShTbOUT%ShTbOErrFlag .NE. 0) THEN !RS: Debugging: Formerly ShTbOUT(7)
                        SELECT CASE (INT(ShTbOUT%ShTbOErrFlag))   !RS: Debugging: Formerly ShTbOUT(7)
                        CASE (1)
                            ShTbPAR%ShTbTID=ShTbPAR%ShTbTID*1.2   !RS: Debugging: Formerly ShTbPAR(2)
                            CYCLE
                        END SELECT
                    END IF

                    XMRFLD=ShTbOUT%ShTbOMdotE*3600/UnitM    !RS Comment: Unit Conversion, lbm/s??   !RS: Debugging: Formerly ShTbOUT(1)
                    ToExp=ShTbOUT%ShTbOToE    !RS: Debugging: Formerly ShTbOUT(3)
                    XoExp=ShTbOUT%ShTbOXoE    !RS: Debugging: Formerly ShTbOUT(4)

                    ErrXMR=ABS((XMRFLD-XMR))
                    IF (MaxDshTb .NE. 0 .AND. MinDshTb .NE. 0 .AND. ErrXMR .GT. 1E-4) THEN
                        IF (XMRFLD .GT. XMR) THEN
                            MaxDshTb=Dshtb
                        ELSE
                            MinDshTb=Dshtb
                        END IF
                        Dshtb=(MaxDshTb+MinDshTb)/2
                        ShTbPAR%ShTbTID=Dshtb/1000 !Short tube diameter, m   !RS: Debugging: Formerly ShTbPAR(2)
                    ELSEIF (ErrXMR .GT. 1E-4) THEN !Find short tube diameter by secant method
                        IF (XMRFLD .GT. XMR) THEN
                            MaxDshTb=Dshtb
                            DshTb=DshTb/(2**NumIter*Dstep)
                        ELSE
                            MinDshTb=Dshtb
                            DshTb=DshTb+2**NumIter*Dstep
                        END IF
                        IF (MaxDshTb .NE. 0 .AND. MinDshTb .NE. 0) THEN
                            Dshtb=(MaxDshTb+MinDshTb)/2
                        END IF
                        ShTbPAR%ShTbTID=Dshtb/1000 !Short tube diameter, m   !RS: Debugging: Formerly ShTbPAR(2)
                    ELSE
                        EXIT
                    END IF
                END DO
            END IF

            IF (INT(ShTbOUT%ShTbOErrFlag) .EQ. 1) THEN    !RS: Debugging: Formerly ShTbOUT(7)
                CALL IssueOutputMessage( '')
                CALL IssueOutputMessage('## ERROR ## HPdesign: Short tube solution error.')
                STOP
            END IF

            !**************Size TXV**************
            mdotr=CompOUT%CmpOMdot    !RS: Debugging: Formerly CompOUT(2)
            PiCmp=CompIN%CompInPsuc !RS: Debugging: Formerly CompIN(1)
            PoCmp=CompIN%CompInPdis !RS: Debugging: Formerly CompIN(2)
            Subcooling=CondOUT%COuttSCiE  !RS: Debugging: Formerly CondOUT(14)
            Superheat=EvapOUT%EOuttSHiC   !RS: Debugging: Formerly EvapOUT(10)
            IF (ShTbOUT%ShTbOPoE .NE. 0) THEN !RS: Debugging: Formerly ShTbOUT(2)
                DPtxv=CondOUT%COutpRiE-ShTbOUT%ShTbOPoE  !RS: Debugging: Formerly CondOUT(10), ShTbOUT(2)
            ElSE
                DPtxv=CondOUT%COutpRiE-EvapIN%EInpRi !RS: Debugging: Formerly EvapIN(2), CondOUT(10)
            END IF

            CALL TXV(mdotr,PiCmp,PoCmp,Subcooling,Superheat,DPtxv,Qtxv)    !RS: Debugging: Testing: Just commenting this out for now
            TxvOUT%TXVQ=Qtxv  !RS: Debugging: Formerly TxvPAR(1)

            !**************Size Capillary Tube**************
            CapTubeIN%CTIMdot=CompOUT%CmpOMdot  !Compressor mass flow rate, kg/s   !RS: Debugging: Formerly CapTubeIN(1), CompOUT(2)
            CapTubeIN%CTIPiEx=CondOUT%COutpRiE !Exp. device inlet pressure, kPa   !RS: Debugging: Formerly CapTubeIN(2), CondOUT(10)
            CapTubeIN%CTIHiEx=CondOUT%COuthRiE !Exp. device inlet enthalpy, kJ/kg !RS: Debugging: Formerly CapTubeIN(3), CondOUT(11)
            CapTubeIN%CTIPiEv=EvapIN%EInpRi   !Evaporator inlet pressure, kPa   !RS: Debugging: Formerly CapTubeIN(4), EvapIN(2)
            CapTubeIN%CTIPoEv=EvapOUT%EOutpRoC  !Evaporator outlet pressure, kPa   !RS: Debugging: Formerly CapTubeIN(5), EvapOUT(1)

            !Initial guess
            NumIter=0
            MaxLen=0
            MinLen=0

            IsSizeDiameter=.TRUE. !Always size diameter

            CapTubeDimension=1e-4 !1E-3 !Initial guess of capillary tube diameter
            IF (IsSizeDiameter .EQ. .TRUE.) THEN
                CapTubePAR%CTTubeID=CapTubeDimension  !RS: Debugging: Formerly CapTubePAR(1)
            ELSE
                CapTubePAR%CTTubeLen=CapTubeDimension  !RS: Debugging: Formerly CapTubePAR(2)
            END IF

            DO NumIter=1, MaxIter

                !CALL CapillaryTubeChoi(Ref$,PureRef,CapTubeIN,CapTubePAR,CapTubeOUT)  
                !CALL CapillaryTubeORNL(Ref$,PureRef,CapTubeIN,CapTubePAR,CapTubeOUT)
                CALL CapillaryTubeORNL !(Ref$) !,CapTubeIN,CapTubePAR,CapTubeOUT)    !RS: Debugging: Extraneous PureRef

                IF (CapTubeOUT%CTOErrFlag .NE. 0) THEN   !RS: Debugging: Formerly CapTubeOUT(2) 
                    SELECT CASE (INT(CapTubeOUT%CTOErrFlag))   !RS: Debugging: Formerly CapTubeOUT(2) 
                    CASE (1)
                        CapTubePAR%CTTubeID=CapTubePAR%CTTubeID*1.2 !RS: Debugging: Formerly CapTubePAR(1)
                        CYCLE
                    END SELECT
                END IF

                XMRFLD=CapTubeOUT%CTOMdot*3600/UnitM !RS Comment: Unit Conversion, lbm/s??   !RS: Debugging: Formerly CapTubeOUT(1)
                ToExp=CapTubeOUT%CTOToE !RS: Debugging: Formerly CapTubeOUT(3)
                XoExp=CapTubeOUT%CTOXoE !RS: Debugging: Formerly CapTubeOUT(4)

                ErrXMR=ABS((XMRFLD-XMR))

                IF (MaxLen .NE. 0 .AND. MinLen .NE. 0 .AND. ErrXMR .GT. 1E-2) THEN
                    IF (XMRFLD .GT. XMR) THEN
                        IF (IsSizeDiameter .EQ. .TRUE.) THEN
                            MaxLen=CapTubeDimension
                        ELSE
                            MinLen=CapTubeDimension
                        END IF
                    ELSE
                        IF (IsSizeDiameter .EQ. .TRUE.) THEN
                            MinLen=CapTubeDimension
                        ELSE
                            MaxLen=CapTubeDimension
                        END IF
                    END IF
                    CapTubeDimension=(MaxLen+MinLen)/2

                ELSEIF (ErrXMR .GT. 1E-2) THEN !Find capillary tube dimension by secant method

                    IF (XMRFLD .GT. XMR) THEN
                        IF (IsSizeDiameter .EQ. .TRUE.) THEN
                            MaxLen=CapTubeDimension
                        ELSE
                            MinLen=CapTubeDimension
                        END IF
                        CapTubeDimension=CapTubeDimension/(2**NumIter*CapTubeDimStep)
                    ELSE
                        IF (IsSizeDiameter .EQ. .TRUE.) THEN
                            MinLen=CapTubeDimension
                        ELSE
                            MaxLen=CapTubeDimension
                        END IF
                        CapTubeDimension=CapTubeDimension+2**NumIter*CapTubeDimStep
                    END IF
                    IF (MaxLen .NE. 0 .AND. MinLen .NE. 0) THEN
                        CapTubeDimension=(MaxLen+MinLen)/2
                    END IF

                ELSE
                    EXIT
                END IF

                IF (IsSizeDiameter .EQ. .TRUE.) THEN
                    CapTubePAR%CTTubeID=CapTubeDimension  !RS: Debugging: Formerly CapTubePAR(1)
                ELSE
                    CapTubePAR%CTTubeLen=CapTubeDimension  !RS: Debugging: Formerly CapTubePAR(2)
                END IF

            END DO

            IF (NumIter .GT. MaxIter) THEN
                CALL IssueOutputMessage( '')
                CALL IssueOutputMessage('## ERROR ## HPdesign: Capillary tube solution not converged.')
                STOP
            END IF

        END IF

        IF (LPRINT.LE.1.AND.IMASS.NE.0) THEN
            IF (AccumPAR%AccH .GT. 0) THEN !Height    !RS: Debugging: Formerly AccumPAR(2)
                AccumIN%AccImdot=MdotR    !RS: Debugging: Formerly AccumIN(1)
                AccumIN%AccIpRo=CompIN%CompInPsuc !Pressure  !RS: Debugging: Formerly CompIN(1), AccumIN(2)
                AccumIN%AccIhRo=CompIN%CompInHsuc !Enthalpy  !RS: Debugging: Formerly CompIN(3), AccumIN(3)
                CALL CalcAccumulatorMass !(AccumIN,AccumOUT)
            ELSE
                AccumOUT%AccOMass=0   !RS: Debugging: Formerly AccumOUT(1)
            END IF

            CALL CalcCondenserInventory(MassCoil,MassLiqCoil,MassVapCoil,CondLiqTubeLength,CondVapTubeLength,CondTwoPhaseTubeLength,CondNumLiqTubes)
            CondOUT%COutMC=MassCoil    !RS: Debugging: Formerly CondOUT(18)
            CALL CalcEvaporatorInventory(MassCoil,MassLiqCoil,MassVapCoil,EvapLiqTubeLength,EvapVapTubeLength,EvapTwoPhaseTubeLength,EvapNumLiqTubes)
            EvapOUT%EOutMC=MassCoil    !RS: Debugging: Formerly EvapOUT(14)

            CALCHG=(CompOUT%CmpOMCmp+CondOUT%COutMDisLn+CondOUT%COutMLiqLn+CondOUT%COutMC+ &   !RS: Debugging: Formerly CompOUT(6), CondOUT(16), CondOUT(17), CondOUT(18)
            EvapOUT%EOutMSucLn+EvapOUT%EOutMC+ShTbOUT%ShTbOMDT+AccumOUT%AccOMass)/UnitM   !RS: Debugging: Formerly EvapOUT(13), EvapOUT(14), ShTbOUT(5), AccumOUT(1)
        END IF

        IF(ICHRGE.EQ.0.AND.ERRMSG(1).NE.0.) THEN 
            WRITE(tmpString,"('0HPDM: **** FAILED TO CONVERGE ON SUBCOOLING *****',/,  '         DIFFERENCE  =',F8.3,' F')") ERRMSG(1)
            CALL IssueOutputMessage(TRIM(tmpString))
        END IF 
        IF(ICHRGE.EQ.0.AND.ERRMSG(2).NE.0.) THEN 
            WRITE(tmpString,"('0HPDM: **** FAILED TO CONVERGE ON SUPERHEAT *****',/,   '         DIFFERENCE  =',F8.3,' F')") ERRMSG(2)
            CALL IssueOutputMessage(TRIM(tmpString))
        END IF

        IF (IsChargeTuning .GT. 0 .AND. MODE .NE. 2) THEN !Apply charge tuning
            ChargeCorrection=(ChargeCurveIntercept+ChargeCurveSlope*(CondLiqTubeLength-RefLiquidLength))/UnitM

            IF (FirstTimeChargeLoop) THEN
                IF (CALCHG + ChargeCorrection .LT. 0) THEN
                    ChargeOption=1
                ELSE
                    ChargeOption=2
                END IF
            END IF

            IF (ChargeOption .EQ. 2) THEN
                CALCHG = CALCHG + ChargeCorrection
            ELSE
                CALCHG = CALCHG
            END IF
        END IF

        RETURN

    END IF

    CALL CalcCondenserInventory(MassCoil,MassLiqCoil,MassVapCoil,CondLiqTubeLength,CondVapTubeLength,CondTwoPhaseTubeLength,CondNumLiqTubes)
    CondOUT%COutMC=MassCoil    !RS: Debugging: Formerly CondOUT(18)

    CALCHG=(CompOUT%CmpOMCmp+CondOUT%COutMDisLn+CondOUT%COutMLiqLn+CondOUT%COutMC)/UnitM   !RS: Debugging: Formerly CompOUT(6), CondOUT(16), CondOUT(17), CondOUT(18) 

    IF (IsChargeTuning .GT. 0 .AND. MODE .NE. 2) THEN !Apply charge tuning
        ChargeCorrection=(ChargeCurveIntercept+ChargeCurveSlope*(CondLiqTubeLength-RefLiquidLength))/UnitM

        IF (FirstTimeChargeLoop) THEN
            IF (CALCHG + ChargeCorrection .LT. 0) THEN
                ChargeOption=1
            ELSE
                ChargeOption=2
            END IF
        END IF

        IF (ChargeOption .EQ. 2) THEN
            CALCHG = CALCHG + ChargeCorrection
        ELSE
            CALCHG = CALCHG
        END IF
    END IF

    RETURN

    END SUBROUTINE
