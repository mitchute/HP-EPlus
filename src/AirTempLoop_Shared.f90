! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This function calculates the air-side superheating or subcooling.  This is done by utilizing the evaporator model. 
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This function does not model any specific component of the Heat Pump.
! 
! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module takes the air inputs into the evaporator and uses them to determine 
! the superheating or subcooling temperatures in the evaporator.
!
! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! No specified input files
! No specified output files

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There are no variables called at the module level; all of the variables are called at the function level

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method:
!    REAL FUNCTION EVPTR - calculates the evaporator superheat and subcool    

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! No bugs/tickets outstanding

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-12 | JEH | Header Revisions

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! No notes/ToDo/Recommendations
!

REAL FUNCTION EVPTR(TINPUT,IERR)
    !
    !       EVPTR(TEMPERATURE) = SUPCAL(TEMPERATURE) - SUPERE
    !       EVPTR IS USED WITH THE ROOT SOLVER 'ZERO' TO FIND TAIIE SO
    !       THAT ABS(SUPCAL(TAIIE) - SUPERE) < EVPCON.
    !       'ZERO' CONTAINS ALL OF THE LOGIC NECESSARY TO ITERATE TO
    !       A ROOT, TAIIE.
    !
    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE EvaporatorMod
    USE AccumulatorModule
    USE DataSimulation
    USE DataGlobals_HPSimIntegrated, ONLY: MaxNameLength, RefrigIndex  !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
    USE UnitConvertMod, ONLY: Temperature_F2C

    IMPLICIT NONE

    REAL Quality,Pressure,Enthalpy

    REAL TINPUT
    INTEGER IERR

    INTEGER(2) RefPropErr			!Error flag:1-error; 0-no error
    
    REAL,PARAMETER :: StandardDensity=1.2 !kg/m3

    REAL TAIIE,TSATEI,SUPCAL,SUPCL,TSATCI
    REAL SUPR,TSAVG,SXIC
    REAL TsatEvp,TsatCnd,Subcooling,Superheat,AccumDP,Xliq,Xvap
    
    CHARACTER(LEN=13),PARAMETER :: FMT_800 = "(A41,F7.2,A5)"
    CHARACTER(LEN=13),PARAMETER :: FMT_804 = "(A32,F7.2,A5)"

    CHARACTER(LEN=MaxNameLength) :: PrintString ! placeholder for formatted output strings
    
    INTEGER :: PVReport=27  !RS: Debugging: Keeping track of cycle points (7/21/14)
    INTEGER, SAVE :: LowIteration=1   !RS: Debugging: Keeping track of cycle points (7/21/14)
    CHARACTER(LEN=25),PARAMETER :: FMT_127 = "(1(I3,','),2(F10.3,','))" !RS: Debugging: Keeping track of cycle points (7/21/14)
    
    OPEN(Unit=PVReport,file='PVReport.csv') !RS: Debugging: Keeping track of cycle points (7/21/14)

    ! initialize some parameters
    EVPTR = 1.0E+10
    IERR = 0
    TAIIE = TINPUT

    IF (TAIIE .LT. TSICMP) THEN
        IERR=2
        RETURN
    END IF

    IF (Unit .EQ. 1) THEN
        WRITE(PrintString, FMT_800) '>> Evaporator entering air temperature: ',Temperature_F2C(TAIIE),Tunit
    ELSE
        WRITE(PrintString, FMT_800) '>> Evaporator entering air temperature: ',TAIIE,Tunit
    END IF
    !CALL IssueOutputMessage( '')
    !CALL IssueOutputMessage(TRIM(PrintString))    

    HiEvp=EvapIN%EInhRi !RS: Debugging: Formerly EvapIN(3)

    IF (FirstTimeAirTempLoop) THEN
        EvapIN%EInpRi=CompIN%CompInPsuc+(EvapIN%EInpRi-EvapOUT%EOutpRiC)  !RS: Debugging: Formerly EvapIN(2), CompIN(1), EvapOUT(6)
        FirstTimeAirTempLoop=.FALSE.
    END IF

    PiEvp=EvapIN%EInpRi !RS: Debugging: Formerly EvapIN(2)

    Pressure=PiEvp*1000 !RS Comment: Unit Conversion
    Quality=0
    TSATEI=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)
    IF (RefPropErr .GT. 0) THEN
        CALL IssueOutputMessage( '-- WARNING -- LowSide: Refprop error.')
        IERR=1
        RETURN
    END IF
    TSATEI=TSATEI*1.8+32    !RS Comment: Unit Conversion, from C to F

    EvapIN%EInmRef=MdotR			!Refrigerant side mass flow rate, kg/s  !RS: Debugging: Formerly EvapIN(1)
    EvapIN%EInpRi=PiEvp			!Evap. inlet pressure, kPa  !RS: Debugging: Formerly EvapIN(2)
    EvapIN%EInhRi=HiEvp			!Refrigerant side inlet enthalpy, kJ/kg !RS: Debugging: Formerly EvapIN(3)
    EvapIN%EInmAi=XMaE            !Air side mass flow rate, kg/s    !RS: Debugging: Formerly EvapIN(4)
    EvapIN%EIntAi=(TAIIE-32)/1.8  !Air side inlet temp. C   !RS: Debugging: Formerly EvapIn(5)
    EvapIN%EInrhAi=RHiE            !Air side inlet relative humidity !RS: Debugging: Formerly EvapIN(6)
    EvapIN%EIntRdis=CompOUT%CmpOTdis      !Discharge temperature, C !RS: Debugging: Formerly EvapIn(9), CompOUT(5)

    !Take compressor shell loss into account
    IF (CompPAR%CompQLossFrac .NE. 0) THEN !Shell loss in fraction    !RS: Debugging: Formerly CompPAR(21)
        EvapPAR%EvapCompQLoss=CompPAR%CompQLossFrac*CompOUT%CmpOPwr  !RS: Debugging: Formerly EvapPAR(32), CompPAR(21), CompOUT(1)
    ELSE !Shell loss in W
        EvapPAR%EvapCompQLoss=CompPAR%CompQLoss/1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapPAR(32) & CompPAR(22)
    END IF

    CALL Evaporator(Ref$) !,EvapIN,EvapPAR,EvapOUT) !(Ref$,PureRef,EvapIN,EvapPAR,EvapOUT) !RS: Debugging: Extraneous PureRef
    CALL PrintEvaporatorResult 
    EvapPAR%EvapFirstTime=0 !No longer first time !RS: Debugging: Formerly EvapPAR(38)
    IF (EvapOUT%EOutErrFlag .NE. 0) THEN    !RS: Debugging: Formerly EvapOUT(17)
        SELECT CASE (INT(EvapOUT%EOutErrFlag))  !RS: Debugging: Formerly EvapOUT(17)
        CASE (2)
            IERR=1
            RETURN
        CASE (3,4,5)
            CALL IssueOutputMessage('Press return to terminate program.')
            OPEN(UNIT=19, FILE='NC.txt')    !RS: Debugging: Trying to set up a buffer program (10/9/14)
            WRITE(19,*) 'Initializing "Not Converged" file'
            CLOSE(19)
            OPEN(20, FILE='Crash.txt', STATUS='old')   !RS: Debugging: Trying to set up a buffer program (10/9/14)
            CLOSE(20, STATUS='DELETE') !RS: Debugging: Trying to set up a buffer program (10/9/14)
            WRITE(6,*) 'HPSim did not converge' !RS: Debugging: Using the log file to let wrapper program know if HPSim has crashed (12/19/14)
            WRITE(6,*) 'Fatal error recognised in EVPTR'
            CLOSE(6)
            STOP   !RS: Debugging: Pushing through for now.
        END SELECT
    END IF

    PoEvp=EvapOUT%EOutpRoC    !RS: Debugging: Formerly EvapOUT(1)
    HoEvp=EvapOUT%EOuthRoC    !RS: Debugging: Formerly EvapOUT(2)
    PiCmp=EvapOUT%EOutpRiC    !RS: Debugging: Formerly EvapOUT(6)
    HiCmp=EvapOUT%EOuthRiC    !RS: Debugging: Formerly EvapOUT(7)
    XiCmp=EvapOUT%EOutxRiC    !RS: Debugging: Formerly EvapOUT(9)

    IF (AccumPAR%AccH .GT. 0) THEN !Accumulator exists    !RS: Debugging: Formerly AccumPAR(2)
        TsatEvp=(TSICMP-32)*5/9     !RS Comment: Unit Conversion, from F to C
        TsatCnd=(TSOCMP-32)*5/9     !RS Comment: Unit Conversion, from F to C
        Subcooling=CondOUT%COuttSCiE  !RS: Debugging: Formerly CondOUT(14)
        Superheat=EvapOUT%EOuttSHiC   !RS: Debugging: Formerly EvapOUT(10)
        Xliq=CondOUT%COutxRiE    !RS: Debugging: Formerly CondOUT(13)
        Xvap=XiCmp
        CALL CalcAccumulatorDP(MdotR,TsatEvp,TsatCnd,Subcooling,Superheat, &
        Xliq,Xvap,AccumDP)

        PiCmp=PiCmp-AccumDP
        EvapOUT%EOutpRiC=PiCmp    !RS: Debugging: Formerly EvapOUT(6)

        Pressure=PiCmp*1000 !RS Comment: Unit Conversion
        Enthalpy=HiCmp*1000 !RS Comment: Unit Conversion
        XiCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr)   !Compressor Inlet Quality
        IF (RefPropErr .GT. 0) THEN
            CALL IssueOutputMessage( '-- WARNING -- LowSide: Refprop error.')
            IERR=1
            RETURN
        END IF
    END IF	
    
    WRITE(PVReport,FMT_127) LowIteration,EvapOUT%EOutpRoC,EvapOUT%EOuthRoC !RS: Debugging: Keeping track of cycle points (7/21/14)
    LowIteration=LowIteration+1 !RS: Debugging: Keeping track of cycle points (7/21/14)

    Pressure=PiCmp*1000 !RS Comment: Unit Conversion
    Quality=1
    TSATCI=PQ(Ref$,Pressure,Quality,'temperature',RefrigIndex,RefPropErr)
    IF (RefPropErr .GT. 0) THEN
        CALL IssueOutputMessage( '-- WARNING -- LowSide: Refprop error.')
        IERR=1
        RETURN
    END IF
    TSATCI=TSATCI*1.8+32    !RS Comment: Unit Conversion, from C to F

    TSAVG=(TSATEI+TSATCI)/2
    IF(TSAVG.GT.TAIIE) THEN
        IERR=2
    END IF

    IF (IERR .GE. 1) THEN
        RETURN
    END IF
    SUPCAL=EvapOUT%EOuttSHiC*1.8 !ISI - 10/07/06  !RS: Debugging: Formerly EvapOUT(10)
    SUPCL = SUPCAL

    IF (XICMP .LT. 0.) THEN !Edited for refprop table 01-15-04 - ISI
        SUPCL = SUPCAL*10 !-SUPCAL !-500.0*SUPCAL
    ELSEIF (XICMP .LT. 1. .AND. XICMP .GT. 0) THEN
        SUPCL = -500.0*(1. - XICMP)
    END IF

    SUPR = SUPER
    IF (SUPER .LT. 0.0) THEN
        SUPR = -500.*(1.0 + SUPER)
    END IF
    EVPTR = SUPCL - SUPR

    IF(SUPER.LT.0.0) THEN
        SXIC = -SUPER
        WRITE(PrintString,FMT_804) '           Desired quality = ',SXIC*100,Xunit
    ELSE
        IF (Unit .EQ. 1) THEN
            WRITE(PrintString,FMT_804) '           Desired superheat = ',SUPER/1.8,DTunit
        ELSE
            WRITE(PrintString,FMT_804) '           Desired superheat = ',SUPER,DTunit
        END IF
    END IF
    !CALL IssueOutputMessage( PrintString)

    !This IF block will always report one and only one message based on calculated "quality"
    IF (XICMP .LT. 0.0) THEN
        IF (Unit .EQ. 1) THEN
            WRITE(PrintString,FMT_804) '       Calculated subcooling = ',-SUPCAL/1.8,DTunit
        ELSE
            WRITE(PrintString,FMT_804) '       Calculated subcooling = ',-SUPCAL,DTunit
        END IF
    ELSEIF (XICMP.LT.1.0) THEN
        WRITE(PrintString,FMT_804)'        Calculated quality = ',XICMP*100,Xunit
    ELSE
        IF (Unit .EQ. 1) THEN
            WRITE(PrintString,FMT_804)'        Calculated superheat = ',SUPCAL/1.8,DTunit
        ELSE  
            WRITE(PrintString,FMT_804)'        Calculated superheat = ',SUPCAL,DTunit
        END IF
    END IF
    !CALL IssueOutputMessage(PrintString)

    RETURN

END FUNCTION
