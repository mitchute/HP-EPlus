! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This routine writes outputs to the YorkHP.out file.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component does not represent any physical component of the system

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This routine writes the file YorkHP.out.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! Output File:
!   YorkHP.out

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There are no variables or structures defined at the module level.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 methods:
!    PUBLIC DumpOutputs -- Writes the .out file outputs
!      Called by ORNLsolver.f90

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! There are no known issues.

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2013-12-17 | RAS | Filled out the header 

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! This could probably use a fair amount of clean-up and probably some revision.
! I don't think there's a need to rename all the variables brought in when they
! could be written out straight from data structures.

SUBROUTINE DumpOutputs

USE FluidProperties_HPSim
USE AirPropMod
USE DataSimulation
USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code
USE DataSimulation, ONLY: inputratio    !RS: The compressor capacity ratio for use in the new compressor model (10/3/14)

IMPLICIT NONE

INTEGER(2) RefPropErr			!Error flag:1-error; 0-no error

INTEGER(2) AirPropOpt			!Air prop calc. option
INTEGER(2) AirPropErr			!Error flag:1-error; 0-no error
!REAL AirProp(8)		!Air properties

REAL Qcnd,Qevp,QevpSens,QevpLat
REAL PwrCmp,mdot,Qtxv
REAL Dshtb,DrawBlow,CPair,Qtot,QtotSens
REAL DcapTube,LcapTube
REAL TAOCND,TAOEVP,COP,SHR,EER

REAL TSATICMP,TSUBICMP,TSUPICMP
REAL TSATOCMP,TSUBOCMP,TSUPOCMP
REAL TSATICND,TSUBICND,TSUPICND
REAL TSATOCND,TSUBOCND,TSUPOCND
REAL TSATIEXP,TSUBIEXP,TSUPIEXP
REAL TSATOEXP,TSUBOEXP,TSUPOEXP
REAL TSATIEVP,TSUBIEVP,TSUPIEVP
REAL TSATOEVP,TSUBOEVP,TSUPOEVP

REAL TDBICND,TWBICND,RHICND,HAICND
REAL TDBOCND,TWBOCND,RHOCND,DPACND,WAOCND,HAOCND
REAL TDBIEVP,TWBIEVP,RHIEVP,HAIEVP
REAL TDBOEVP,TWBOEVP,RHOEVP,DPAEVP,HAOEVP,WAOEVP

REAL MassCmp,MassCnd,MassEvp,MassSucLn,MassDisLn,MassLiqLn,MassDistube
REAL AccumDP,FilterDP,MassAccum
REAL WeightEvpAluminum,WeightEvpCopper,WeightCndAluminum,WeightCndCopper
REAL,PARAMETER :: StandardDensity=1.2 !Standard density, kg/m3
REAL,PARAMETER :: StandardSpecHeat=1.02 !Standard specific heat, kJ/kg-K
REAL SpecHeat     !Specific heat, kJ/kg-K
REAL Quality,Pressure,Enthalpy
CHARACTER (len=50) :: Title !Output file title

REAL :: CoilSurfTemp = 0.0

REAL Qtest,Qtestlat,Qtestsens !RS: Debugging: Trying to debug capacity mismatch (6/11/14)
REAL TotalPower   !RS: Debugging: Setting a total power variable for GenOpt to easily minimize (9/30/14)

CHARACTER(LEN=13),PARAMETER :: FMT_2200 = "(A32,',',A50)"
CHARACTER(LEN=3),PARAMETER :: FMT_2204 = "(A)"
CHARACTER(LEN=23),PARAMETER :: FMT_2208 = "(A33,',',F40.3,',',A15)"
CHARACTER(LEN=69),PARAMETER :: FMT_2212 = "(A24,',',A18,',',A18,',',A21,',',A29,',',A14,',',A17,',',A16,',',A26)"
CHARACTER(LEN=83),PARAMETER :: FMT_2216 = "(A24,',',F18.3,',',F18.3,',',F21.3,',',F29.3,',',F14.1,',',F17.3,',',F16.3,',',A26)"
CHARACTER(LEN=53),PARAMETER :: FMT_2220 = "(A18,',',A27,',',A27,',',A24,',',A33,',',A25,',',A20)"
CHARACTER(LEN=63),PARAMETER :: FMT_2224 = "(A22,',',F27.3,',',F27.3,',',F24.1,',',F33.3,',',F25.3,',',A22)"

  SELECT CASE(MODE)
  CASE(FIXEDORIFICESIM)
 
	  IF (SystemType .EQ. 4) THEN
		  Title='Short Tube Simulation, Reheat mode'
	  ELSEIF (IsCoolingMode .GT. 0) THEN	  
		  Title='Short Tube Simulation, Cooling Mode'
	  ELSE
		  Title='Short Tube Simulation, Heating Mode'
	  END IF
  CASE(ORIFICEANDTXVDESIGN)

	  IF (SystemType .EQ. 4) THEN
		  Title='Short Tube Simulation, Reheat mode'
	  ELSEIF (IsCoolingMode .GT. 0) THEN	  
		  Title='Fixed Subcooling Design Calculation, Cooling Mode'
	  ELSE
		  Title='Fixed Subcooling Design Calculation, Heating Mode'
	  END IF
  CASE(FIXEDSUPERHEATSIM)
   
	  IF (SystemType .EQ. 4) THEN
		  Title='Short Tube Simulation, Reheat mode'
	  ELSEIF (IsCoolingMode .GT. 0) THEN	  
		  Title='Fixed Orifice Design Calculation, Cooling Mode'
	  ELSE
		  Title='Fixed Orifice Design Calculation, Heating Mode'
	  END IF
  CASE(TXVSIMULATION)
  
	  IF (SystemType .EQ. 4) THEN
		  Title='Short Tube Simulation, Reheat mode'
	  ELSEIF (IsCoolingMode .GT. 0) THEN	  
		  Title='TXV Simulation, Cooling Mode'
	  ELSE
		  Title='TXV Simulation, Heating Mode'
	  END IF
  CASE(CONDENSERUNITSIM)  
	  Title='Condenser Unit, Cooling Mode'
  CASE(COILONLYSIM)  
	  IF (IsCoolingMode .GT. 0) THEN
	    Title='Coil Only, Cooling Mode'
	  ELSE
	    Title='Coil Only, Heating Mode'
	  END IF
  END SELECT

  !Calculate corresponding output depending on the unit type
  IF (MODE .EQ. CONDENSERUNITSIM) THEN

      !*******Compressor data*******
      PiCmp=CompIN%CompInPsuc   !RS: Debugging: Formerly CompIN(1)
      HiCmp=CompIN%CompInHsuc   !RS: Debugging: Formerly CompIN(3)

      Pressure=PiCmp*1000   !RS Comment: Unit Conversion
      Enthalpy=HiCmp*1000   !RS Comment: Unit Conversion
      TiCmp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Compressor Inlet Temperature
      XiCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Compressor Inlet Quality

      Quality=1
      TsatiCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature

      TsupiCmp=TiCmp-TsatiCmp
      IF (TsupiCmp .LT. 0 ) THEN
          TsupiCmp=0
      END IF

      Quality=0
      TsatiCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature
      TsubiCmp=TsatiCmp-TiCmp
      IF (TsubiCmp .LT. 0 ) THEN
          TsubiCmp=0
      END IF

      IF (XiCmp .GE. 1) THEN
          XiCmp=1
      END IF
      IF (XiCmp .LE. 0) THEN
          XiCmp=0
      END IF
  
      PoCmp=CompIN%CompInPdis   !RS: Debugging: Formerly CompIN(2)
      HoCmp=CompOUT%CmpOHdis  !RS: Debugging: Formerly CompOUT(3)

      Pressure=PoCmp*1000   !RS Comment: Unit Conversion
      Enthalpy=HoCmp*1000   !RS Comment: Unit Conversion
      ToCmp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Compressor Outlet Temperature
      XoCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Compressor Outlet Quality

      Quality=1
      TsatoCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature

      TsupoCmp=ToCmp-TsatoCmp
      IF (TsupoCmp .LT. 0 ) THEN
          TsupoCmp=0
      END IF

      Quality=0
      TsatoCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature
      TsuboCmp=TsatoCmp-ToCmp
      IF (TsuboCmp .LT. 0 ) THEN
          TsuboCmp=0
      END IF

      IF (XoCmp .GE. 1) THEN
          XoCmp=1
      END IF
      IF (XoCmp .LE. 0) THEN
          XoCmp=0
      END IF

      PwrCmp=CompOUT%CmpOPwr*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly CompOUT(1)
      MassCmp=CompOUT%CmpOMCmp    !RS: Debugging: Formerly CompOUT(6)

      !*******Condenser*******
      PiCnd=CondOUT%COutpRiC  !RS: Debugging: Formerly CondOUT(1)
      HiCnd=CondOUT%COuthRiC  !RS: Debugging: Formerly CondOUT(2)

      Pressure=PiCnd*1000   !RS Comment: Unit Conversion
      Enthalpy=HiCnd*1000   !RS Comment: Unit Conversion
      TiCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Inlet Temperature
      XiCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Inlet Quality

      Quality=1
      TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

      TsupiCnd=TiCnd-TsatiCnd
      IF (TsupiCnd .LT. 0 ) THEN
          TsupiCnd=0
      END IF

      Quality=0
      TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
      TsubiCnd=TsatiCnd-TiCnd
      IF (TsubiCnd .LT. 0 ) THEN
          TsubiCnd=0
      END IF

      IF (XiCnd .GE. 1) THEN
          XiCnd=1
      END IF
      IF (XiCnd .LE. 0) THEN
          XiCnd=0
      END IF

      TdbiCnd=CondIN%CIntAi !RS: Debugging: Formerly CondIN(5)
      RHiCnd=CondIN%CInrhAi  !RS: Debugging: Formerly CondIN(6)

      AirPropOpt=2
      AirProp%APTDB=TdbiCnd    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=RHiCnd !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
      hAiCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
      TwbiCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
      RhoAiC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

      CPair=CPA(TdbiCnd)    !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(TdbiCnd,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)

      PoCnd=CondOUT%COutpRoC  !RS: Debugging: Formerly CondOUT(5)
      HoCnd=CondOUT%COuthRoC  !RS: Debugging: Formerly CondOUT(6)

      Pressure=PoCnd*1000   !RS Comment: Unit Conversion
      Enthalpy=HoCnd*1000   !RS Comment: Unit Conversion
      ToCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Outlet Temperature
      XoCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Outlet Quality

      Quality=1
      TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

      TsupoCnd=ToCnd-TsatoCnd
      IF (TsupoCnd .LT. 0 ) THEN
          TsupoCnd=0
      END IF

      Quality=0
      TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
      TsuboCnd=TsatoCnd-ToCnd
      IF (TsuboCnd .LT. 0 ) THEN
          TsuboCnd=0
      END IF

      IF (XoCnd .GE. 1) THEN
          XoCnd=1
      END IF
      IF (XoCnd .LE. 0) THEN
          XoCnd=0
      END IF

      Qcnd =CondOUT%COutQC*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(15)

      TdboCnd=CondOUT%COuttAoC   !RS: Debugging: Formerly CondOUT(3)
      RHoCnd=CondOUT%COutrhAoC    !RS: Debugging: Formerly CondOUT(4)

      AirPropOpt=2
      AirProp%APTDB=TdboCnd    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=RHoCnd !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
      WaoCnd=AirProp%APHumRat !RS: Debugging: Formerly AirProp(2)
      hAoCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
      TwboCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
      RhoAoC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

      DPaCND=CondOUT%COutDPAir*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(19)

      MassCnd=CondOUT%COutMC   !RS: Debugging: Formerly CondOUT(18)
      MassDisLn=CondOUT%COutMDisLn !RS: Debugging: Formerly CondOUT(16)
      MassLiqLn=CondOUT%COutMLiqLn !RS: Debugging: Formerly CondOUT(17)

      WeightCndAluminum=CondOUT%COutWtAl !RS: Debugging: Formerly CondOUT(8)
      WeightCndCopper=CondOUT%COutWtCu   !RS: Debugging: Formerly CondOUT(9)

      !*******Evaporator*******
      PiEvp=0
      HiEvp=0
      TiEvp=0
      XiEvp=0
      TsatiEvp=0
      TsupiEvp=0
      TsubiEvp=0
      TdbiEvp=0
      RHiEvp=0
      hAiEvp=0
      TwbiEvp=0
      RhoAiE=0
      PoEvp=0
      HoEvp=0

      ToEvp=0
      XoEvp=0
      TsatoEvp=0
      TsupoEvp=0
      TsatoEvp=0
      TsuboEvp=0

      CFMevp=0
      Qevp=mdotR*(HiCmp-HoCnd)*1000
      QevpSens=0
      QevpLat=0
      
      TdboEvp=0
      RHoEvp=0
      WaoEvp=0
      hAoEvp=0
      TwboEvp=0
      RhoAoE=0

      DPaEVP=0

      MassEvp=0
      MassSucLn=0

      WeightEvpAluminum=0
      WeightEvpCopper=0      

	  !*******Exp. device*******
	  PiExp=0
	  HiExp=0
	  MassDisTube=0

	  Dshtb=0
	  DcapTube=0
	  LcapTube=0
	  Qtxv=0
    	  
	  PiExp=0
	  HiExp=0
	  MassDisTube=0
    	  	  
	  TiExp=0
	  XiExp=0

	  TsatiExp=0
      TsupiExp=0
	  TsatiExp=0
      TsubiExp=SUBCOOL/1.8
	  PoExp=0
	  HoExp=0

	  ToExp=0
	  XoExp=0

	  TsatoExp=0
	  TsupoExp=0
	  TsuboExp=0

      EER=0
      COP=0

      CalChg=0
      MassAccum=0
      AccumDP=0
      FilterDP=0
  
  ELSEIF (MODE .EQ. COILONLYSIM) THEN

      IF (IsCoolingMode .GT. 0) THEN

          !*******Condenser*******
          PiCnd=0
          HiCnd=0
          TiCnd=0
          XiCnd=0
          TsatiCnd=0
          TsupiCnd=0
          TsubiCnd=0
          TdbiCnd=0
          RHiCnd=0
          hAiCnd=0
          TwbiCnd=0
          RhoAiC=0
          PoCnd=0
          HoCnd=0

          ToCnd=CondOUT%COuttRoC  !RS: Debugging: Formerly CondOUT(7)
          XoCnd=0
          TsatoCnd=0
          TsupoCnd=0
          TsatoCnd=0
          TsuboCnd=0

          CFMcnd=0
          Qcnd =0

          TdboCnd=0
          RHoCnd=0
          WaoCnd=0
          hAoCnd=0
          TwboCnd=0
          RhoAoC=0

          DPaCND=0

          MassCnd=0
          MassDisLn=0
          MassLiqLn=0

          WeightCndAluminum=0
          WeightCndCopper=0

	      !*******Evaporator*******
	      PiEvp=EvapIN%EInpRi   !RS: Debugging: Formerly EvapIN(2)
	      HiEvp=EvapIN%EInhRi   !RS: Debugging: Formerly EvapIN(3)

	      Pressure=PiEvp*1000   !RS Comment: Unit Conversion
	      Enthalpy=HiEvp*1000   !RS Comment: Unit Conversion
	      TiEvp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Evaporator Inlet Temperature
	      XiEvp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Evaporator Inlet Quality

	      Quality=1
	      TsatiEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature

	      TsupiEvp=TiEvp-TsatiEvp
	      IF (TsupiEvp .LT. 0 ) THEN
              TsupiEvp=0
          END IF

	      Quality=0
	      TsatiEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature
	      TsubiEvp=TsatiEvp-TiEvp
	      IF (TsubiEvp .LT. 0 ) THEN
              TsubiEvp=0
          END IF

	      IF (XiEvp .GE. 1) THEN
              XiEvp=1
          END IF
	      IF (XiEvp .LE. 0) THEN
              XiEvp=0
          END IF

	      TdbiEvp=EvapIN%EIntAi !RS: Debugging: Formerly EvapIN(5)
	      RHiEvp=EvapIN%EInrhAi  !RS: Debugging: Formerly EvapIN(6)

	      AirPropOpt=2
	      AirProp%APTDB=TdbiEvp    !RS: Debugging: Formerly AirProp(1)
	      AirProp%APRelHum=RHiEvp !RS: Debugging: Formerly AirProp(3)
	      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	      hAiEvp=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
	      TwbiEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	      RhoAiE=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

	      CPair=CPA(TdbiEvp)    !RS: Replace: CPA (2/19/14)
          CPair=CPAirFunction(TdbiEvp,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)

	      PoEvp=EvapOUT%EOutpRoC  !RS: Debugging: Formerly EvapOUT(1)
	      HoEvp=EvapOUT%EOuthRoC  !RS: Debugging: Formerly EvapOUT(2)

	      Pressure=PoEvp*1000   !RS Comment: Unit Conversion
	      Enthalpy=HoEvp*1000   !RS Comment: Unit Conversion
	      ToEvp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Evaporator Outlet Temperature
	      XoEvp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Evaporator Outlet Temperature

	      Quality=1
	      TsatoEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature

	      TsupoEvp=ToEvp-TsatoEvp
	      IF (TsupoEvp .LT. 0 ) THEN
              TsupoEvp=0
          END IF

	      Quality=0
	      TsatoEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature
	      TsuboEvp=TsatoEvp-ToEvp
	      IF (TsuboEvp .LT. 0 ) THEN
              TsuboEvp=0
          END IF

	      IF (XoEvp .GE. 1) THEN
              XoEvp=1
          END IF
	      IF (XoEvp .LE. 0) THEN
              XoEvp=0
          END IF

	      DPaEvp=EvapOUT%EOutDPAir*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(5)

	      Qevp =-EvapOUT%EOutQC*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(11)
	      QevpSens=-EvapOUT%EOutQCSens*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(12)
	      IF (ABS(QevpSens) .GT. ABS(Qevp)) THEN !Make sure sensible heat is never higher than total heat. ISI - 08/02/07
	          QevpSens = Qevp
	          hAoEvp=-Qevp/1000/(CFMevp*RhoAiE)+hAiEvp
	          SpecHeat=CPA(TdbiEvp) !RS: Replace: CPA (2/19/14)
              SpecHeat=CPAirFunction(TdbiEvp,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
	          TdboEvp=-QevpSens/1000/(CFMevp*RhoAiE*SpecHeat)+TdbiEvp
	          AirPropOpt=1
	          AirProp%APTDB=TdboEvp    !RS: Debugging: Formerly AirProp(1)
	          AirProp%APEnth=hAoEvp !RS: Debugging: Formerly AirProp(4)
	          CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	          TwboEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	          RHoEvp=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
	      ELSE
    	      TdboEvp=EvapOUT%EOuttAoC   !RS: Debugging: Formerly EvapOUT(3)
	          RHoEvp=EvapOUT%EOutrhAoC   !RS: Debugging: Formerly EvapOUT(4)
	          AirPropOpt=2
	          AirProp%APTDB=TdboEvp    !RS: Debugging: Formerly AirProp(1)
	          AirProp%APRelHum=RHoEvp !RS: Debugging: Formerly AirProp(3)
	          CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	          TwboEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	          RhoAoE=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)
	      END IF

	      QevpLat=Qevp-QevpSens
	      SHR=QevpSens/Qevp

	      MassEvp=EvapOUT%EOutMC   !RS: Debugging: Formerly EvapOUT(14)
	      MassSucLn=EvapOUT%EOutMSucLn+AccumOUT%AccOMass !RS: Debugging: Formerly EvapOUT(13), AccumOUT(1)

	      WeightEvpAluminum=EvapOUT%EOutWtAl !RS: Debugging: Formerly EvapOUT(15)
	      WeightEvpCopper=EvapOUT%EOutWtCu   !RS: Debugging: Formerly EvapOUT(16)

	      !*******Exp. device*******
	      PiExp=0
	      HiExp=0
	      MassDisTube=0

	      Dshtb=0
	      DcapTube=0
	      LcapTube=0
	      Qtxv=0
        	      	  	  
	      TiExp=ToCnd
	      XiExp=1

	      TsatiExp=0

	      TsupiExp=0

          TsubiExp=0

	      IF (XiExp .GE. 1) THEN
              XiExp=1
          END IF
	      IF (XiExp .LE. 0) THEN
              XiExp=0
          END IF

	      PoExp=0
	      HoExp=0

	      ToExp=0
	      XoExp=0

	      TsatoExp=0
	      TsupoExp=0
	      TsuboExp=0

          !*******Compressor data*******
          PoCmp=0
          HoCmp=0
          ToCmp=0
          XoCmp=0
          TsatoCmp=(TSOCMP-32)*5/9  !RS Comment: Unit Conversion, from F to C
          TsupoCmp=0
          TsuboCmp=0

          PiCmp=CompIN%CompInPsuc !PoEvp    !RS: Debugging: Formerly CompIN(1)
          HiCmp=CompIN%CompInHsuc !HoEvp    !RS: Debugging: Formerly CompIN(3)

	      Pressure=PiCmp*1000   !RS Comment: Unit Conversion
	      Enthalpy=HiCmp*1000   !RS Comment: Unit Conversion
	      TiCmp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Compressor Inlet Temperature
	      XiCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Compressor Inlet Quality

          TsatiCmp=(TSICMP-32)*5/9 !TsatoEvp    !RS Comment: Unit Conversion, from F to C
          TsupiCmp=TsupoEvp
          TsubiCmp=TsuboEvp

          PwrCmp=0
          mdot=0
          MassCmp=0

      ELSE

          !*******Condenser*******
          PiCnd=CondOUT%COutpRiC  !RS: Debugging: Formerly CondOUT(1)
          HiCnd=CondOUT%COuthRiC  !RS: Debugging: Formerly CondOUT(2)

          Pressure=PiCnd*1000   !RS Comment: Unit Conversion
          Enthalpy=HiCnd*1000   !RS Comment: Unit Conversion
          TiCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Inlet Temperature
          XiCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Inlet Quality

          Quality=1
          TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

          TsupiCnd=TiCnd-TsatiCnd
          IF (TsupiCnd .LT. 0 ) THEN
              TsupiCnd=0
          END IF

          Quality=0
          TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
          TsubiCnd=TsatiCnd-TiCnd
          IF (TsubiCnd .LT. 0 ) THEN
              TsubiCnd=0
          END IF

          IF (XiCnd .GE. 1) THEN
              XiCnd=1
          END IF
          IF (XiCnd .LE. 0) THEN
              XiCnd=0
          END IF

          TdbiCnd=CondIN%CIntAi !RS: Debugging: Formerly CondIN(5)
          RHiCnd=CondIN%CInrhAi  !RS: Debugging: Formerly CondIN(6)

          AirPropOpt=2
          AirProp%APTDB=TdbiCnd    !RS: Debugging: Formerly AirProp(1)
          AirProp%APRelHum=RHiCnd !RS: Debugging: Formerly AirProp(3)
          CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
          hAiCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
          TwbiCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
          RhoAiC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

          CPair=CPA(TdbiCnd)    !RS: Replace: CPA (2/19/14)
          CPair=CPAirFunction(TdbiCnd,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)

          PoCnd=CondOUT%COutpRoC  !RS: Debugging: Formerly CondOUT(5)
          HoCnd=CondOUT%COuthRoC  !RS: Debugging: Formerly CondOUT(6)

          Pressure=PoCnd*1000   !RS Comment: Unit Conversion
          Enthalpy=HoCnd*1000   !RS Comment: Unit Conversion
          ToCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Outlet Temperature
          XoCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Outlet Quality

          Quality=1
          TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

          TsupoCnd=ToCnd-TsatoCnd
          IF (TsupoCnd .LT. 0 ) THEN
              TsupoCnd=0
          END IF

          Quality=0
          TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
          TsuboCnd=TsatoCnd-ToCnd
          IF (TsuboCnd .LT. 0 ) THEN
              TsuboCnd=0
          END IF

          IF (XoCnd .GE. 1) THEN
              XoCnd=1
          END IF
          IF (XoCnd .LE. 0) THEN
              XoCnd=0
          END IF

          Qcnd =CondOUT%COutQC*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(15)

          TdboCnd=CondOUT%COuttAoC   !RS: Debugging: Formerly CondOUT(3)
          RHoCnd=CondOUT%COutrhAoC    !RS: Debugging: Formerly CondOUT(4)

          AirPropOpt=2
          AirProp%APTDB=TdboCnd    !RS: Debugging: Formerly AirProp(1)
          AirProp%APRelHum=RHoCnd !RS: Debugging: Formerly AirProp(3)
          CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
          WaoCnd=AirProp%APHumRat !RS: Debugging: Formerly AirProp(2)
          hAoCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
          TwboCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
          RhoAoC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

          DPaCND=CondOUT%COutDPAir*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(19)

          MassCnd=CondOUT%COutMC   !RS: Debugging: Formerly CondOUT(18)
          MassDisLn=CondOUT%COutMDisLn !RS: Debugging: Formerly CondOUT(16)
          MassLiqLn=CondOUT%COutMLiqLn !RS: Debugging: Formerly CondOUT(17)

          WeightCndAluminum=CondOUT%COutWtAl !RS: Debugging: Formerly CondOUT(8)
          WeightCndCopper=CondOUT%COutWtCu   !RS: Debugging: Formerly CondOUT(9)

          !*******Evaporator*******
          PiEvp=0
          HiEvp=0
          TiEvp=0
          XiEvp=0
          TsatiEvp=0
          TsupiEvp=0
          TsubiEvp=0
          TdbiEvp=0
          RHiEvp=0
          hAiEvp=0
          TwbiEvp=0
          RhoAiE=0
          PoEvp=0
          HoEvp=0

          ToEvp=0
          XoEvp=0
          TsatoEvp=0
          TsupoEvp=0
          TsatoEvp=0
          TsuboEvp=0

          CFMevp=0
          Qevp=0
          QevpSens=0
          QevpLat=0

          TdboEvp=0
          RHoEvp=0
          WaoEvp=0
          hAoEvp=0
          TwboEvp=0
          RhoAoE=0

          DPaEVP=0
          SHR=1

          MassEvp=0
          MassSucLn=0

          WeightEvpAluminum=0
          WeightEvpCopper=0      

	      !*******Exp. device*******
	      PiExp=PoCnd
	      HiExp=HoCnd
	      MassDisTube=0

	      Dshtb=0
	      DcapTube=0
	      LcapTube=0
	      Qtxv=0
        	          	  	  
	      TiExp=ToCnd
	      XiExp=XoCnd

	      TsatiExp=TsatoCnd
          TsupiExp=TsupoCnd
	      TsatiExp=TsatoCnd
          TsubiExp=TsuboCnd
          
	      PoExp=0
	      HoExp=0

	      ToExp=0
	      XoExp=0

	      TsatoExp=0
	      TsupoExp=0
	      TsuboExp=0

          !*******Compressor data*******
          PoCmp=PiCnd
          HoCmp=HiCnd

          ToCmp=TiCnd
          XoCmp=XiCnd

          TsatoCmp=(TSOCMP-32)*5/9 !TsatiCnd    !RS Comment: Unit Conversion, from F to C
          TsupoCmp=TsupiCnd
          TsuboCmp=TsubiCnd

          IF (XoCmp .GE. 1) THEN
              XoCmp=1
          END IF
          IF (XoCmp .LE. 0) THEN
              XoCmp=0
          END IF

          PiCmp=0
          HiCmp=0
          TiCmp=0
          XiCmp=0
          TsatiCmp=(TSICMP-32)*5/9  !RS Comment: Unit Conversion, from F to C
          TsupiCmp=0
          TsubiCmp=0

          PwrCmp=0
          mdot=0
          MassCmp=0

      END IF

      EER=0
      COP=0

      CalChg=0
      MassAccum=0
      AccumDP=0
      FilterDP=0
  
  ELSE

      !*******Compressor data*******
      PiCmp=CompIN%CompInPSuc   !RS: Debugging: Formerly CompIN(1)
      HiCmp=CompIN%CompInHSuc   !RS: Debugging: Formerly CompIN(3)

      Pressure=PiCmp*1000   !RS Comment: Unit Conversion
      Enthalpy=HiCmp*1000   !RS Comment: Unit Conversion
      TiCmp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Compressor Inlet Temperature
      XiCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Compressor Inlet Quality

      Quality=1
      TsatiCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature

      TsupiCmp=TiCmp-TsatiCmp
      IF (TsupiCmp .LT. 0 ) THEN
          TsupiCmp=0
      END IF

      Quality=0
      TsatiCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature
      TsubiCmp=TsatiCmp-TiCmp
      IF (TsubiCmp .LT. 0 ) THEN
          TsubiCmp=0
      END IF

      IF (XiCmp .GE. 1) THEN
          XiCmp=1
      END IF
      IF (XiCmp .LE. 0) THEN
          XiCmp=0
      END IF
  
      PoCmp=CompIN%CompInPDis   !RS: Debugging: Formerly CompIN(2)
      HoCmp=CompOUT%CmpOHDis  !RS: Debugging: Formerly CompOUT(3)

      Pressure=PoCmp*1000   !RS Comment: Unit Conversion
      Enthalpy=HoCmp*1000   !RS Comment: Unit Conversion
      ToCmp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Compressor Outlet Temperature
      XoCmp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Compressor Outlet Quality

      Quality=1
      TsatoCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature

      TsupoCmp=ToCmp-TsatoCmp
      IF (TsupoCmp .LT. 0 ) THEN
          TsupoCmp=0
      END IF

      Quality=0
      TsatoCmp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Compressor Saturation Temperature
      TsuboCmp=TsatoCmp-ToCmp
      IF (TsuboCmp .LT. 0 ) THEN
          TsuboCmp=0
      END IF

      IF (XoCmp .GE. 1) THEN
          XoCmp=1
      END IF
      IF (XoCmp .LE. 0) THEN
          XoCmp=0
      END IF

      PwrCmp=CompOUT%CmpOPwr*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly CompOUT(1)
      MassCmp=CompOUT%CmpOMCmp    !RS: Debugging: Formerly CompOUT(6)

      !*******Condenser*******
      PiCnd=CondOUT%COutpRiC  !RS: Debugging: Formerly CondOUT(1)
      HiCnd=CondOUT%COuthRiC  !RS: Debugging: Formerly CondOUT(2)

      Pressure=PiCnd*1000   !RS Comment: Unit Conversion
      Enthalpy=HiCnd*1000   !RS Comment: Unit Conversion
      TiCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Inlet Temperature
      XiCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Inlet Quality

      Quality=1
      TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

      TsupiCnd=TiCnd-TsatiCnd
      IF (TsupiCnd .LT. 0 ) THEN
          TsupiCnd=0
      END IF

      Quality=0
      TsatiCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
      TsubiCnd=TsatiCnd-TiCnd
      IF (TsubiCnd .LT. 0 ) THEN
          TsubiCnd=0
      END IF

      IF (XiCnd .GE. 1) THEN
          XiCnd=1
      END IF
      IF (XiCnd .LE. 0) THEN
          XiCnd=0
      END IF

      TdbiCnd=CondIN%CIntAi !RS: Debugging: Formerly CondIN(5)
      RHiCnd=CondIN%CInrhAi  !RS: Debugging: Formerly CondIN(6)

      AirPropOpt=2
      AirProp%APTDB=TdbiCnd    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=RHiCnd !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
      hAiCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
      TwbiCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
      RhoAiC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

      !CPair=CPA(TdbiCnd)    !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(TdbiCnd,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)

      PoCnd=CondOUT%COutpRoC  !RS: Debugging: Formerly CondOUT(5)
      HoCnd=CondOUT%COuthRoC  !RS: Debugging: Formerly CondOUT(6)

      Pressure=PoCnd*1000   !RS Comment: Unit Conversion
      Enthalpy=HoCnd*1000   !RS Comment: Unit Conversion
      ToCnd=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Condenser Outlet Temperature
      XoCnd=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Condenser Outlet Quality

      Quality=1
      TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature

      TsupoCnd=ToCnd-TsatoCnd
      IF (TsupoCnd .LT. 0 ) THEN
          TsupoCnd=0
      END IF

      Quality=0
      TsatoCnd=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Condenser Saturation Temperature
      TsuboCnd=TsatoCnd-ToCnd
      IF (TsuboCnd .LT. 0 ) THEN
          TsuboCnd=0
      END IF

      IF (XoCnd .GE. 1) THEN
          XoCnd=1
      END IF
      IF (XoCnd .LE. 0) THEN
          XoCnd=0
      END IF

      Qtest=mdotR*(HoCnd-HiCnd)!RS: Debugging: Trying to debug capacity mismatch (6/11/14)
      Qcnd =CondOUT%COutQC*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(15)

      TdboCnd=CondOUT%COuttAoC   !RS: Debugging: Formerly CondOUT(3)
      RHoCnd=CondOUT%COutrhAoC    !RS: Debugging: Formerly CondOUT(4)
      
      Qtestsens=CPair*(TdboCnd-TdbiCnd) !RS: Debugging: Trying to debug capacity mismatch (6/11/14)

      AirPropOpt=2
      AirProp%APTDB=TdboCnd    !RS: Debugging: Formerly AirProp(1)
      AirProp%APRelHum=RHoCnd !RS: Debugging: Formerly AirProp(3)
      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
      WaoCnd=AirProp%APHumRat !RS: Debugging: Formerly AirProp(2)
      hAoCnd=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
      TwboCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
      RhoAoC=AirProp%APDryDens  !RS: Debugging: Formerly AirProp(7)

      DPaCND=CondOUT%COutDPAir*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly CondOUT(19)

      MassCnd=CondOUT%COutMC   !RS: Debugging: Formerly CondOUT(18)
      MassDisLn=CondOUT%COutMDisLn !RS: Debugging: Formerly CondOUT(16)
      MassLiqLn=CondOUT%COutMLiqLn !RS: Debugging: Formerly CondoUT(17)

      WeightCndAluminum=CondOUT%COutWtAl !RS: Debugging: Formerly CondOUT(8)
      WeightCndCopper=CondOUT%COutWtCu   !RS: Debugging: Formerly CondOUT(9)

	  !*******Exp. device*******
      IF (ExpDevice .EQ. 3) THEN
	      PiExp=CapTubeIN%CTIPiEx    !RS: Debugging: Formerly CapTubeIN(2)
	      HiExp=CapTubeIN%CTIHiEx    !RS: Debugging: Formerly CapTubeIN(3)
	      MassDisTube=CapTubeOUT%CTOMDT !RS: Debugging: Formerly CapTubeOUT(5)
      ELSE
	      PiExp=ShTbIN%ShTbINPiE   !RS: Debugging: Formerly ShTbIN(2)
	      HiExp=ShTbIN%ShTbINHiE   !RS: Debugging: Formerly ShTbIN(3)
	      MassDisTube=ShTbOUT%ShTbOMDT    !RS: Debugging: Formerly ShTbOUT(5)
      END IF
      
	  Dshtb=ShTbPAR%ShTbTID  !RS: Debugging: Formerly ShTbPAR(2)
	  Qtxv=TxvOUT%TXVQ !RS: Debugging: Formerly TxvPAR(1)
	  DcapTube=CapTubePAR%CTTubeID    !RS: Debugging: Formerly CapTubePAR(1)
	  LcapTube=CapTubePAR%CTTubeLen    !RS: Debugging: Formerly CapTubePAR(2)

	  Pressure=PiExp*1000   !RS Comment: Unit Conversion
	  Enthalpy=HiExp*1000   !RS Comment: Unit Conversion
	  TiExp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Expansion Device Inlet Temperature
	  XiExp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Expansion Device Inlet Quality

	  Quality=1
	  TsatiExp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Saturation Temperature

	  TsupiExp=TiExp-TsatiExp
	  IF (TsupiExp .LT. 0 ) THEN
          TsupiExp=0
      END IF

	  Quality=0
	  TsatiExp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Saturation Temperature
	  TsubiExp=TsatiExp-TiExp
	  IF (TsubiExp .LT. 0 ) THEN
          TsubiExp=0
      END IF

	  IF (XiExp .GE. 1) THEN
          XiExp=1
      END IF
	  IF (XiExp .LE. 0) THEN
          XiExp=0
      END IF

      IF (ExpDevice .EQ. 3) THEN
	      PoExp=CapTubeOUT%CTOPoE   !RS: Debugging: Formerly CapTubeOUT(6)???
	      HoExp=CapTubeIN%CTIHiEx    !RS: Debugging: Formerly CapTubeIN(3)
      ELSE
	      PoExp=ShTbOUT%ShTbOPoE  !RS: Debugging: Formerly ShTbOUT(2)
	      HoExp=ShTbIN%ShTbINHiE   !RS: Debugging: Formerly ShTbIN(3)
      END IF

	  Pressure=PoExp*1000   !RS Comment: Unit Conversion
	  Enthalpy=HoExp*1000   !RS Comment: Unit Conversion
	  ToExp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Expansion Device Outlet Temperature
	  XoExp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Expansion Device Outlet Quality

	  Quality=1
	  TsatoExp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Saturation Temperature

	  TsupoExp=ToExp-TsatoExp
	  IF (TsupoExp .LT. 0 ) THEN
          TsupoExp=0
      END IF

	  Quality=0
	  TsatoExp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Expansion Device Saturation Temperature

	  TsuboExp=TsatoExp-ToExp
	  IF (TsuboExp .LT. 0 ) THEN
          TsuboExp=0
      END IF

	  IF (XoExp .GE. 1) THEN
          XoExp=1
      END IF
	  IF (XoExp .LE. 0) THEN
          XoExp=0
      END IF

	  !*******Evaporator*******
	  PiEvp=EvapIN%EInpRi   !RS: Debugging: Formerly EvapIN(2)
	  HiEvp=EvapIN%EInhRi   !RS: Debugging: Formerly EvapIN(3)
      
      PwrODfan=CondPAR%CondFanPwr !*1000 !RS Comment: 1000 accounts for CondPAR conversion !RS: Debugging: Formerly CondPAR(34)
      PwrIDfan=EvapPAR%EvapFanPwr !*1000 !RS Comment: 1000 accounts for EvapPAR conversion !RS: Debugging: Formerly EvapPAR(27)

	  Pressure=PiEvp*1000   !RS Comment: Unit Conversion
	  Enthalpy=HiEvp*1000   !RS Comment: Unit Conversion
	  TiEvp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Evaporator Inlet Temperature
	  XiEvp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Evaporator Inlet Quality

	  Quality=1
	  TsatiEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature

	  TsupiEvp=TiEvp-TsatiEvp
	  IF (TsupiEvp .LT. 0 ) THEN
          TsupiEvp=0
      END IF

	  Quality=0
	  TsatiEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature
	  TsubiEvp=TsatiEvp-TiEvp
	  IF (TsubiEvp .LT. 0 ) THEN
          TsubiEvp=0
      END IF

	  IF (XiEvp .GE. 1) THEN
          XiEvp=1
      END IF
	  IF (XiEvp .LE. 0) THEN
          XiEvp=0
      END IF

	  TdbiEvp=EvapIN%EIntAi !RS: Debugging: Formerly EvapIN(5)
	  RHiEvp=EvapIN%EInrhAi  !RS: Debugging: Formerly EvapIN(6)

	  AirPropOpt=2
	  AirProp%APTDB=TdbiEvp    !RS: Debugging: Formerly AirProp(1)
	  AirProp%APRelHum=RHiEvp !RS: Debugging: Formerly AirProp(3)
	  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	  hAiEvp=AirProp%APEnth !RS: Debugging: Formerly AirProp(4)
	  TwbiEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	  RhoAiE=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

	  !CPair=CPA(TdbiEvp)    !RS: Replace: CPA (2/19/14)
      CPair=CPAirFunction(TdbiEvp,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)

	  PoEvp=EvapOUT%EOutpRoC  !RS: Debugging: Formerly EvapOUT(1)
	  HoEvp=EvapOUT%EOuthRoC  !RS: Debugging: Formerly EvapOUT(2)

	  Pressure=PoEvp*1000   !RS Comment: Unit Conversion
	  Enthalpy=HoEvp*1000   !RS Comment: Unit Conversion
	  ToEvp=PH(Ref$, Pressure, Enthalpy, 'temperature', RefrigIndex,RefPropErr) !Evaporator Outlet Temperature
	  XoEvp=PH(Ref$, Pressure, Enthalpy, 'quality', RefrigIndex,RefPropErr) !Evaporator Outlet Quality

	  Quality=1
	  TsatoEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature

	  TsupoEvp=ToEvp-TsatoEvp
	  IF (TsupoEvp .LT. 0 ) THEN
          TsupoEvp=0
      END IF

	  Quality=0
	  TsatoEvp=PQ(Ref$, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr)   !Evaporator Saturation Temperature
	  TsuboEvp=TsatoEvp-ToEvp
	  IF (TsuboEvp .LT. 0 ) THEN
          TsuboEvp=0
      END IF

	  IF (XoEvp .GE. 1) THEN
          XoEvp=1
      END IF
	  IF (XoEvp .LE. 0) THEN
          XoEvp=0
      END IF

	  DPaEvp=EvapOUT%EOutDPAir*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(5)
      
      Qtest=mdotR*(HoEvp-HiEvp)!RS: Debugging: Trying to debug capacity mismatch (6/11/14)

	  Qevp =-EvapOUT%EOutQC*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(11)
	  QevpSens=-EvapOUT%EOutQCSens*1000    !RS Comment: Unit Conversion    !RS: Debugging: Formerly EvapOUT(12)
	  IF (ABS(QevpSens) .GT. ABS(Qevp)) THEN !Make sure sensible heat is never higher than total heat. ISI - 08/02/07
	      QevpSens = Qevp
	      hAoEvp=-Qevp/1000/(CFMevp*RhoAiE)+hAiEvp
	      !SpecHeat=CPA(TdbiEvp) !RS: Replace: CPA (2/19/14)
          SpecHeat=CPAirFunction(TdbiEvp,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
	      TdboEvp=-QevpSens/1000/(CFMevp*RhoAiE*SpecHeat)+TdbiEvp
	      AirPropOpt=1
	      AirProp%APTDB=TdboEvp    !RS: Debugging: Formerly AirProp(1)
	      AirProp%APEnth=hAoEvp !RS: Debugging: Formerly AirProp(4)
	      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	      TwboEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	      RHoEvp=AirProp%APRelHum !RS: Debugging: Formerly AirProp(3)
	  ELSE
    	  TdboEvp=EvapOUT%EOuttAoC   !RS: Debugging: Formerly EvapOUT(3)
	      RHoEvp=EvapOUT%EOutrhAoC    !RS: Debugging: Formerly EvapOUT(4)
	      AirPropOpt=2
	      AirProp%APTDB=TdboEvp    !RS: Debugging: Formerly AirProp(1)
	      AirProp%APRelHum=RHoEvp !RS: Debugging: Formerly AirProp(3)
	      CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	      TwboEvp=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	      RhoAoE=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)
      END IF

      SpecHeat=CPAirFunction(TdbiEvp,AirProp%APHumRat)  !RS: Debugging: Trying to debug capacity mismatch (6/11/14)
      Qtestsens=SpecHeat*(Tdboevp-Tdbievp) !RS: Debugging: Trying to debug capacity mismatch (6/11/14)
      
      Qtestsens=SpecHeat*(Tdboevp-TsatoEvp) !RS: Debugging: Trying to debug capacity mismatch (6/11/14)
      Qtestlat=Qtest-Qtestsens !RS: Debugging: Trying to debug capacity mismatch (6/11/14)
      
	  QevpLat=Qevp-QevpSens

	  MassEvp=EvapOUT%EOutMC   !RS: Debugging: Formerly EvapOUT(14)
	  MassSucLn=EvapOUT%EOutMSucLn+AccumOUT%AccOMass !RS: Debugging: Formerly EvapOUT(13), AccumOUT(1)

	  WeightEvpAluminum=EvapOUT%EOutWtAl !RS: Debugging: Formerly EvapOUT(15)
	  WeightEvpCopper=EvapOUT%EOutWtCu   !RS: Debugging: Formerly EvapOUT(16)

	  IF (IsCoolingMode .GT. 0) THEN !ISI - 11/03/2008
	    DrawBlow=EvapPAR%EvapFanLoc    !RS: Debugging: Formerly EvapPAR(28)
	  ELSE
	    DrawBlow=CondPAR%CondFanLoc    !RS: Debugging: Formerly CondPAR(35)
	  END IF

	  IF (DrawBlow .LT. 2) THEN !draw through
		  IF (IsCoolingMode .GT. 0) THEN
			  EER=((Qevp-PwrIDfan)/UnitPwr)/(PwrCmp+PwrIDfan+PwrODfan)
			  COP=(Qevp-PwrIDfan)/(PwrCmp+PwrIDfan+PwrODfan)
			  SHR=(QevpSens-PwrIDfan)/(Qevp-PwrIDfan)
		  ELSE
			  EER=((Qcnd+PwrIDfan)/UnitPwr)/(PwrCmp+PwrIDfan+PwrODfan)
			  COP=(Qcnd+PwrIDfan)/(PwrCmp+PwrIDfan+PwrODfan)
			  SHR=(Qcnd+PwrIDfan)/(Qcnd+PwrIDfan)
		  END IF
	  ELSE
		  IF (IsCoolingMode .GT. 0) THEN
			  EER=(Qevp/UnitPwr)/(PwrCmp+PwrIDfan+PwrODfan)
			  COP=Qevp/(PwrCmp+PwrIDfan+PwrODfan)
			  SHR=QevpSens/Qevp
		  ELSE
			  EER=(Qcnd/UnitPwr)/(PwrCmp+PwrIDfan+PwrODfan)
			  COP=Qcnd/(PwrCmp+PwrIDfan+PwrODfan)
			  SHR=Qcnd/Qcnd
		  END IF
	  END IF

      CalChg=CALCHG*UnitM
      Dshtb=ShTbPAR%ShTbTID*1000 !RS Comment: Unit Conversion    !RS: Debugging: Formerly ShTbPAR(2)
	  DcapTube=CapTubePAR%CTTubeID*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly CapTubePAR(1)
	  LcapTube=CapTubePAR%CTTubeLen*1000   !RS Comment: Unit Conversion    !RS: Debugging: Formerly CapTubePAR(2)

      MassAccum=AccumOUT%AccOMass !RS: Debugging: Formerly AccumOUT(1)
      AccumDP=AccumOUT%AccODP   !RS: Debugging: Formerly AccumOUT(2)
      FilterDP=FilterOUT%FODP !RS: Debugging: Formerly FilterOUT(1)
  
  END IF
  
  IF (SystemType .EQ. 4) THEN !Reheat system
	  IF (DrawBlow .LT. 2) THEN !draw through
		  Qtot=-Qevp+Qcnd+PwrIDfan
		  QtotSens=-QevpSens+Qcnd+PwrIDfan
	  ELSE
		  Qtot=-Qevp+Qcnd
		  QtotSens=-QevpSens+Qcnd	  
	  END IF
	  	
	  !CPair=CPA(REAL(TdboCnd))  !RS: Replace: CPA (2/19/14)
      !RS: Debugging: Should update the humidity ratio (2/20/14)
      CPair=CPAirFunction(TdboCnd,AirProp%APHumRat)  !RS: Replace: CPA (2/19/14)
	  hAoCnd=Qtot/1000/(CFMevp*StandardDensity)+hAiEvp
	  TdboCnd=QtotSens/1000/(CFMevp*StandardDensity*CPair)+TdbiEvp

	  AirPropOpt=1
	  AirProp%APTDB=TdboCnd    !RS: Debugging: Formerly AirProp(1)
	  AirProp%APEnth=hAoCnd !RS: Debugging: Formerly AirProp(4)
	  CALL PsyChart(AirPropOpt,AirPropErr)  !(AirProp, ,BaroPressure,  
	  WaoCnd=AirProp%APHumRat !RS: Debugging: Formerly AirProp(2)
	  AirProp%APRelHum=RHoCnd !RS: Debugging: Formerly AirProp(3)
	  TwboCnd=AirProp%APTWB    !RS: Debugging: Formerly AirProp(5)
	  RhoAoC=AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)

	  Qevp=Qtot
	  QevpSens=QtotSens
	  SHR=QevpSens/Qevp
  END IF

  mdot=MdotR*3600   !RS Comment: Unit Conversion
      
  !Convert output data to IP unit
  IF (Unit .EQ. 2) THEN
      PICMP=PiCmp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
      HICMP=HiCmp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
      TICMP=TiCmp*1.8+32    !RS Comment: Unit Conversion, from C to F
      POCMP=PoCmp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
      HOCMP=HoCmp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
      TOCMP=ToCmp*1.8+32    !RS Comment: Unit Conversion, from C to F
      TSUPICMP=TsupiCmp*1.8
	  TSUPOCMP=TsupoCmp*1.8
      TSUBICMP=TsubiCmp*1.8
	  TSUBOCMP=TsuboCmp*1.8
	  TSATICMP=TsatiCmp*1.8+32  !RS Comment: Unit Conversion, from C to F
	  TSATOCMP=TsatoCmp*1.8+32  !RS Comment: Unit Conversion, from C to F
	  MDOT=mdot/UrefFlow

      PICND=PiCnd/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  POCND=PoCnd/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  TICND=TiCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TOCND=ToCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  HICND=HiCnd/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
	  HOCND=HoCnd/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
      TSUPICND=TsupiCnd*1.8
	  TSUPOCND=TsupoCnd*1.8
      TSUBICND=TsubiCnd*1.8
	  TSUBOCND=TsuboCnd*1.8
	  TSATICND=TsatiCnd*1.8+32  !RS Comment: Unit Conversion, from C to F
	  TSATOCND=TsatoCnd*1.8+32  !RS Comment: Unit Conversion, from C to F
	  QCND=Qcnd/UnitPwr
	  TDBICND=TdbiCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TDBOCND=TdboCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TWBICND=TwbiCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TWBOCND=TwboCnd*1.8+32    !RS Comment: Unit Conversion, from C to F
	  CFMCND=CFMcnd/UnitArFlw
	  DPACND=DPaCnd/UairPres

      PIEXP=PiExp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
      HIEXP=HiExp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
      TIEXP=TiExp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  POEXP=PoExp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  HOEXP=HoExp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
	  TOEXP=ToExp*1.8+32    !RS Comment: Unit Conversion, from C to F
      TSUPIEXP=TsupiExp*1.8
	  TSUPOEXP=TsupoExp*1.8
      TSUBIEXP=TsubiExp*1.8
	  TSUBOEXP=TsuboExp*1.8
	  TSATIEXP=TsatiExp*1.8+32  !RS Comment: Unit Conversion, from C to F
	  TSATOEXP=TsatoExp*1.8+32  !RS Comment: Unit Conversion, from C to F

      PIEVP=PiEvp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  POEVP=PoEvp/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  TIEVP=TiEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TOEVP=ToEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  HIEVP=HiEvp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
	  HOEVP=HoEvp/UnitH     !RS Comment: Unit Conversion, from kJ/kg to Btu/lbm
      TSUPIEVP=TsupiEvp*1.8
	  TSUPOEVP=TsupoEvp*1.8
      TSUBIEVP=TsubiEvp*1.8
	  TSUBOEVP=TsuboEvp*1.8
	  TSATIEVP=TsatiEvp*1.8+32  !RS Comment: Unit Conversion, from C to F
	  TSATOEVP=TsatoEvp*1.8+32  !RS Comment: Unit Conversion, from C to F

	  QEVP=Qevp/UnitPwr
	  QEVPSENS=QevpSens/UnitPwr
	  QEVPLAT=QevpLat/UnitPwr
	  TDBIEVP=TdbiEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TDBOEVP=TdboEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TWBIEVP=TwbiEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  TWBOEVP=TwboEvp*1.8+32    !RS Comment: Unit Conversion, from C to F
	  CFMEVP=CFMevp/UnitArFlw
	  DPAEVP=DPaEvp/UairPres

	  CALCHG=CalChg/UnitM       !RS Comment: Unit Conversion, from kg to lbm
	  MassCmp=MassCmp/UnitM     !RS Comment: Unit Conversion, from kg to lbm
	  MassCnd=MassCnd/UnitM     !RS Comment: Unit Conversion, from kg to lbm
	  MassEvp=MassEvp/UnitM     !RS Comment: Unit Conversion, from kg to lbm
	  MassSucLn=MassSucLn/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  MassDisLn=MassDisLn/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  MassLiqLn=MassLiqLn/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  MassDisTube=MassDisTube/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  MassAccum=MassAccum/UnitM     !RS Comment: Unit Conversion, from kg to lbm
      Dshtb=ShTbPAR%ShTbTID/UnitL*12         !RS Comment: Unit Conversion, from m to in  !RS: Debugging: Formerly ShTbPAR(2)
	  DcapTube=CapTubePAR%CTTubeID/UnitL*12   !RS Comment: Unit Conversion, from m to in  !RS: Debugging: Formerly CapTubePAR(1)
	  LcapTube=CapTubePAR%CTTubeLen/UnitL*12   !RS Comment: Unit Conversion, from m to in  !RS: Debugging: Formerly CapTubePAR(2)

      TaoEVP=TaoEVP*1.8+32  !RS Comment: Unit Conversion, from C to F
	  TaoCND=TaoCND*1.8+32  !RS Comment: Unit Conversion, from C to F
    
	  AccumDP=AccumDP/UnitP     !RS Comment: Unit Conversion, from kPa to psi
	  FilterDP=FilterDP/UnitP   !RS Comment: Unit Conversion, from kPa to psi

	  WeightEvpAluminum=WeightEvpAluminum/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  WeightEvpCopper=WeightEvpCopper/UnitM     !RS Comment: Unit Conversion, from kg to lbm

	  WeightCndAluminum=WeightCndAluminum/UnitM !RS Comment: Unit Conversion, from kg to lbm
	  WeightCndCopper=WeightCndCopper/UnitM     !RS Comment: Unit Conversion, from kg to lbm

	  WeightSucLn=WeightSucLn/UnitM             !RS Comment: Unit Conversion, from kg to lbm
	  WeightDisLn=WeightDisLn/UnitM             !RS Comment: Unit Conversion, from kg to lbm
	  WeightLiqLn=WeightLiqLn/UnitM             !RS Comment: Unit Conversion, from kg to lbm
	  WeightValveIDCLn=WeightValveIDCLn/UnitM   !RS Comment: Unit Conversion, from kg to lbm
	  WeightValveODCLn=WeightValveODCLn/UnitM   !RS Comment: Unit Conversion, from kg to lbm

	  CondLiqTubeLength=CondLiqTubeLength/UnitL             !RS Comment: Unit Conversion, from m to ft
	  CondVapTubeLength=CondVapTubeLength/UnitL             !RS Comment: Unit Conversion, from m to ft
	  CondTwoPhaseTubeLength=CondTwoPhaseTubeLength/UnitL   !RS Comment: Unit Conversion, from m to ft

	  EvapLiqTubeLength=EvapLiqTubeLength/UnitL             !RS Comment: Unit Conversion, from m to ft
	  EvapVapTubeLength=EvapVapTubeLength/UnitL             !RS Comment: Unit Conversion, from m to ft
	  EvapTwoPhaseTubeLength=EvapTwoPhaseTubeLength/UnitL   !RS Comment: Unit Conversion, from m to ft

	  !Conver Pressure to gauge basis
      !RS Comment: Unit Conversion, from kPa to psi
	  PiCmp=PiCmp-BaroPressure/UnitP
	  PoCmp=PoCmp-BaroPressure/UnitP
	  PiCnd=PiCnd-BaroPressure/UnitP
	  PoCnd=PoCnd-BaroPressure/UnitP
	  PiEvp=PiEvp-BaroPressure/UnitP
	  PoEvp=PoEvp-BaroPressure/UnitP
	  PiExp=PiExp-BaroPressure/UnitP
	  PoExp=PoExp-BaroPressure/UnitP

  END IF

  TotalPower=PWRCMP+FanOut%Power !+FanOutE%Power  !RS: Debugging: Setting a total power variable for GenOpt to easily minimize (9/30/14)
  
  !Conver Pressure to gauge basis
  IF (Unit .EQ. 1) THEN
	  PiCmp=PiCmp-BaroPressure
	  PoCmp=PoCmp-BaroPressure
	  PiCnd=PiCnd-BaroPressure
	  PoCnd=PoCnd-BaroPressure
	  PiEvp=PiEvp-BaroPressure
	  PoEvp=PoEvp-BaroPressure
	  PiExp=PiExp-BaroPressure
	  PoExp=PoExp-BaroPressure
  END IF
CoilSurfTemp = 0.0
CoilSurfTemp=CoilParams(2)%TSurfCoil*9/5+32 !RS Comment: Unit Conversion, from C to F
  !VB program report format
  WRITE(5,*)
  WRITE(5,FMT_2200)'Title (ver. 2.0 12/17/09): ',Title
  WRITE(5,*)
  WRITE(5,FMT_2204)'*** System Performance Summary ***'
  WRITE(5,*)
  WRITE(5,FMT_2208)'Evaporator gross capacity       ',QEVP,CapUnit
  WRITE(5,FMT_2208)'Condenser gross capacity        ',QCND,CapUnit
  WRITE(5,FMT_2208)'Gross sensible capacity         ',QEVPSENS,CapUnit
  WRITE(5,FMT_2208)'Gross latent capacity           ',QEVPLAT,CapUnit
  WRITE(5,FMT_2208)'Compressor power                ',PWRCMP,PwrUnit
  WRITE(5,FMT_2208)'Condenser Fan power             ',FanOut%Power,PwrUnit
  !WRITE(5,FMT_2208)'Evaporator Fan power            ',FanOutE%Power,PwrUnit    !RS: Debugging: Removing Evap Fan calls since handled on E+ side: Strictly for coupling/optimization! (11/11/14)
  WRITE(5,FMT_2208)'Total Power                     ',TotalPower,PwrUnit    !RS: Debugging: Setting a total power variable for GenOpt to easily minimize (9/30/14)
  WRITE(5,FMT_2208)'Refrigerant mass flow rate      ',MDOT,mdotUnit
  WRITE(5,FMT_2208)'COP (coefficient of performance)',COP,NoUnit
  WRITE(5,FMT_2208)'EER (energy efficiency ratio)   ',EER,EERunit
  WRITE(5,FMT_2208)'SHR (sensible heat ratio)       ',SHR,NoUnit
  WRITE(5,FMT_2208)'Condenser subcooling            ',TSUBOCND,DTunit
  WRITE(5,FMT_2208)'Expansion device subcooling     ',TSUBIEXP,DTunit
  WRITE(5,FMT_2208)'Evaporator superheat            ',TSUPOEVP,DTunit
  WRITE(5,FMT_2208)'Compressor superheat            ',TSUPICMP,DTunit
  WRITE(5,FMT_2208)'System charge                   ',CALCHG,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in compressor       ',MassCmp,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in condenser        ',MassCnd,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in evaporator       ',MassEvp,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in suction line     ',MassSucLn,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in discharge line   ',MassDisLn,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in liquid line      ',MassLiqLn,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in distributor tubes',MassDistube,MassUnit
  WRITE(5,FMT_2208)'Refrigerant in accumulator      ',MassAccum,MassUnit
  WRITE(5,FMT_2208)'Short tube diameter             ',Dshtb,MiniLunit
  WRITE(5,FMT_2208)'TXV capacity                    ',Qtxv,SysUnit
  WRITE(5,FMT_2208)'Capillary tube diameter         ',DcapTube,MiniLunit !Added for Cap. tube 04/13/2009 - ISI
  WRITE(5,FMT_2208)'Capillary tube length           ',LcapTube,MiniLunit !Added for Cap. tube 04/13/2009 - ISI 
  WRITE(5,*)
  WRITE(5,FMT_2204)'***** Refrigerant Side Data *****'
  WRITE(5,*)
  IF (Unit .EQ. 2) THEN
      WRITE(5,FMT_2212)'Location               ','Temperature (F)','Pressure (psig)','Enthalpy (Btu/lbm)','Saturation Temperature (F)','Quality (%)','Subcooling (R)','Superheat (R)','Location               '
  ELSE
      WRITE(5,FMT_2212)'Location               ','Temperature (C)','Pressure (kPa)','Enthalpy (kJ/kg)','Saturation Temperature (C)','Quality (%)','Subcooling (K)','Superheat (K)','Location               '
  END IF
  WRITE(5,FMT_2216)'Compressor suction     ',TICMP,PICMP,HICMP,TSATICMP,XICMP*100,TSUBICMP,TSUPICMP,'Compressor suction     '
  WRITE(5,FMT_2216)'Compressor discharge   ',TOCMP,POCMP,HOCMP,TSATOCMP,XOCMP*100,TSUBOCMP,TSUPOCMP,'Compressor discharge   '
  WRITE(5,FMT_2216)'Condenser inlet        ',TICND,PICND,HICND,TSATICND,XICND*100,TSUBICND,TSUPICND,'Condenser inlet        '
  WRITE(5,FMT_2216)'Condenser outlet       ',TOCND,POCND,HOCND,TSATOCND,XOCND*100,TSUBOCND,TSUPOCND,'Condenser outlet       '
  WRITE(5,FMT_2216)'Expansion device inlet ',TIEXP,PIEXP,HIEXP,TSATIEXP,XIEXP*100,TSUBIEXP,TSUPIEXP,'Expansion device inlet '
  WRITE(5,FMT_2216)'Expansion device outlet',TOEXP,POEXP,HOEXP,TSATOEXP,XOEXP*100,TSUBOEXP,TSUPOEXP,'Expansion device outlet'
  WRITE(5,FMT_2216)'Evaporator inlet       ',TIEVP,PIEVP,HIEVP,TSATIEVP,XIEVP*100,TSUBIEVP,TSUPIEVP,'Evaporator inlet       '
  WRITE(5,FMT_2216)'Evaporator outlet      ',TOEVP,POEVP,HOEVP,TSATOEVP,XOEVP*100,TSUBOEVP,TSUPOEVP,'Evaporator outlet      '
  WRITE(5,*)
  WRITE(5,FMT_2204)'********* Air Side Data *********'
  WRITE(5,*)
  IF (Unit .EQ. 2) THEN
      WRITE(5,FMT_2220)'Location         ','Dry bulb temperature (F)','Wet bulb temperature (F)','Relative Humidity (%)','Volumetric flow rate (CFM)','Pressure Drop (in-H2O)','Location         '
  ELSE
      WRITE(5,FMT_2220)'Location         ','Dry bulb temperature (C)','Wet bulb temperature (C)','Relative Humidity (%)','Volumetric flow rate (m^3/min)','Pressure Drop (Pa)','Location         '
  END IF
  WRITE(5,FMT_2224)'Condenser air inlet  ',TDBICND,TWBICND,RHICND*100,CFMCND,0.0d0,'Condenser air inlet  '
  WRITE(5,FMT_2224)'Condenser air outlet ',TDBOCND,TWBOCND,RHOCND*100,CFMCND,DPACND,'Condenser air outlet '
  WRITE(5,FMT_2224)'Evaporator air inlet ',TDBIEVP,TWBIEVP,RHIEVP*100,CFMEVP,0.0d0,'Evaporator air inlet '
  WRITE(5,FMT_2224)'Evaporator air outlet',TDBOEVP,TWBOEVP,RHOEVP*100,CFMEVP,DPAEVP,'Evaporator air outlet'
  WRITE(5,*)
  WRITE(5,FMT_2204)'***** Optimization Outputs ******'  !RS: Debugging: Setting to pass out the exiting air properties in the output file (10/29/14)
  WRITE(5,*)
  WRITE(5,FMT_2208)'Evaporator Air Enthalpy          ',EvapOUT%AirEnth,'kJ/kg' !RS: Debugging: Setting to pass out the exiting air enthalpy in the output file (10/29/14)
  WRITE(5,FMT_2208)'Evaporator Air Humidity Ratio    ',EvapOUT%AirHumRat,'' !RS: Debugging: Setting to pass out the exiting air humidity ratio in the output file (10/29/14)
  WRITE(5,FMT_2208)'Compressor Input Ratio           ',inputratio,'' !RS: Debugging: Setting to pass out the compressor ratio (11/15/14)
  WRITE(5,*)
  WRITE(5,FMT_2204)'********* Pressure Drop *********'
  WRITE(5,*)
  WRITE(5,FMT_2208)'Accumulator                     ',AccumDP,Punit
  WRITE(5,FMT_2208)'Filter Drier                    ',FilterDP,Punit   
  WRITE(5,*)
  WRITE(5,FMT_2204)'*********** Material ************'
  WRITE(5,*)
  WRITE(5,FMT_2208)'Aluminum Evaporator             ',WeightEvpAluminum,MassUnit
  WRITE(5,FMT_2208)'Aluminum Condenser              ',WeightCndAluminum,MassUnit
  WRITE(5,FMT_2208)'Copper Evaporator               ',WeightEvpCopper,MassUnit
  WRITE(5,FMT_2208)'Copper Condenser                ',WeightCndCopper,MassUnit
  WRITE(5,FMT_2208)'Copper Suction Line             ',WeightSucLn,MassUnit
  WRITE(5,FMT_2208)'Copper Discharge Line           ',WeightDisLn,MassUnit
  WRITE(5,FMT_2208)'Copper Valve-Indoor Coil Line   ',WeightValveIDCLn,MassUnit
  WRITE(5,FMT_2208)'Copper Valve-Outdoor Coil Line  ',WeightValveODCLn,MassUnit
  WRITE(5,FMT_2208)'Copper Liquid Line              ',WeightLiqLn,MassUnit
  WRITE(5,*)
  WRITE(5,FMT_2204)'*** Refrigerant Distributions ***'
  WRITE(5,*)
  WRITE(5,FMT_2208)'Condenser liquid length         ',CondLiqTubeLength,Lunit
  WRITE(5,FMT_2208)'Condenser two-phase length      ',CondTwoPhaseTubeLength,Lunit
  WRITE(5,FMT_2208)'Condenser vapor length          ',CondVapTubeLength,Lunit
  WRITE(5,FMT_2208)'Evaporator two-phase length     ',EvapTwoPhaseTubeLength,Lunit
  WRITE(5,FMT_2208)'Evaporator vapor length         ',EvapVapTubeLength,Lunit
  WRITE(5,*)
  WRITE(5,FMT_2204)'*** Frost Parameters ***'
  WRITE(5,*)
  WRITE(5,FMT_2208)'Evaporator Tube Area            ',EvapTubeArea,NoUnit
  WRITE(5,FMT_2208)'Evaporator Fin Area             ',EvapFinArea,NoUnit
  WRITE(5,FMT_2208)'Evaporator Total Area           ',EvapTotArea,NoUnit
  WRITE(5,FMT_2208)'Evaporator Bare Area            ',EvapBareArea,NoUnit
  WRITE(5,FMT_2208)'Evaporator Min Area             ',EvapMinArea,NoUnit
  WRITE(5,FMT_2208)'Condenser Tube Area             ',CondTubeArea,NoUnit
  WRITE(5,FMT_2208)'Condenser Fin Area              ',CondFinArea,NoUnit
  WRITE(5,FMT_2208)'Evaporator Surface Temperature  ',CoilSurfTemp,NoUnit
  WRITE(5,*) 'YorkHP.out successfully written'

  WRITE(6,*) 'YorkHP.out successfully written'  !RS: Debugging: Using the log file to let wrapper program know if HPSim has crashed (12/19/14)
  
RETURN

END SUBROUTINE
