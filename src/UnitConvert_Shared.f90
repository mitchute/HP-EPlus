! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module converts the input units to standard SI units 
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This does not represent a physical component of the heat pump

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module takes the parameters for each of the components and converts the inputs to SI units
! 
! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! no associated files
!
! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! Constants defined at the module level
!   UnitPwr, UnitAirFlw, UnitRfFlw, UnitP, UnitM, UnitL, UnitK
! Flags defined at module level
!   SI, IP

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains X methods:
!    PUBLIC UnitConvert
!       ORNLSolver
!    PUBLIC MicroChannelUnitConvert
!    PUBLIC FinTubeCoilUnitConvert
!       Condenser
!       Evaporator
!    PUBLIC Temperature_F2C
!       ORNLSolver
!       AirTempLoop
!    PUBLIC Temperature_C2F

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! This module is to be removed, according to Ticket #23

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-29 | JEH | Header completion

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! this module should be removed 

MODULE UnitConvertMod

IMPLICIT NONE

PRIVATE
REAL, PARAMETER :: UnitPwr   = 0.2927E-3 !(Btu/hr X UnitPwr = KW)
REAL, PARAMETER :: UnitArFlw = 0.472E-3  !(CFM X UnitArFlw = m^3/s)
REAL, PARAMETER :: UnitRfFlw = 0.4536    !(Lbm/hr X UnitRfFlw = kg/hr)
REAL, PARAMETER :: UnitP     = 6.8947453 !(psi X UnitP = kPa)
REAL, PARAMETER :: UnitM     = 0.4536    !(lbm X UnitM = kg)
REAL, PARAMETER :: UnitL     = 0.3048    !(ft X UnitL = m)
REAL, PARAMETER :: UnitK     = 0.1442E-3 !(Btu-in/hr-ft2-F X UnitK = kW/m-C)

!Unit flags !ISI - 07/14/06
INTEGER,PARAMETER :: SI=1
INTEGER,PARAMETER :: IP=2

PUBLIC UnitConvert
PUBLIC MicroChannelCoilUnitConvert
PUBLIC FinTubeCoilUnitConvert

PUBLIC Temperature_F2C  ! VL
PUBLIC Temperature_C2F  ! VL

CONTAINS

!***********************************************************************************

SUBROUTINE UnitConvert !(CFMcnd,CFMevp) !( !TaiC,TaiE,RHiC,RHiE, & !XMaC,XMaE, !Unit,& !CompPAR,CondPAR,EvapPAR,ShTbPAR,CapTubePAR, & !TxvPAR,  &
                       !AccumPAR,FilterPAR,XMaC,XMaE,TaiC,TaiE,RHiC,RHiE, &
				       !Refchg,TSOCMP,TSICMP,SUPER,SUBCOOL,BaroPressure) !, &
					   !ChargeCurveSlope,ChargeCurveIntercept,RefLiquidLength, &    !RS: Debugging: Removing these
					   !Tdis,Tliq)
! ----------------------------------------------------------------------
!
!   Description: To convert all input data into standard SI
!                Some data needs to be in IP unit for the
!                Corresponding model
!
!   Author:
!   Ipseng Iu
!   Mechanical and Aerospace Engineering
!   Oklahoma State University, Stillwater
!
!   Date:
!   December, 2002
!
! ----------------------------------------------------------------------

USE DataSimulation !, ONLY: CondDisLnOD,CondDisLnTWThick,CondDisLnQLoss, & !RS: Debugging: For replacement of array numbers with variable names
!                        CondLiqLnOD,CondLiqLnTWThick,CondLiqLnQLoss,CondCoilTOD,CondCoilTWThick,CondTspc,CondRspc, &
!                        CondFinThick,CondFanPwr,CondBarPress,CondDisLnLen,CondDisLnElev,CondDisLnTempChg,CondDisLnAddPD,CondLiqLnLen, &
!                        CondLiqLnElev,CondLiqLnTempChg,CondLiqLnAddPD,CondCoilSTLen,CondCoilTThermCon,CondFinThick,CondFinPitch, &
!                        CondFinThermCon,EvapSucLnOD,EvapSucLnTWThick,EvapSucLnQLoss,EvapCoilTOD,EvapCoilTWThick,EvapTspc,EvapRspc, &
!                        EvapFinThick,EvapFanPwr,EvapSucLnLen,EvapSucLnElev,EvapSucLnTempChg,EvapSucLnAddPD,EvapCoilTThermCon, &
!                        EvapCoilSTLen,EvapFinPitch,EvapFinThermCon,EvapBarPress,ShTbTLen,ShTbTID,ShTbChamDep,ShTbDTubeLen,CTTubeID, &
!                        CTTubeLen,CTTubeCoilD,CTDisTubeLen,AccD,AccH,AccD1,AccD2,AccHDis,AccDTube,AccDP,AccDT,CompIntVol,CompQLoss

IMPLICIT NONE

!Subroutine parameters

!Subroutine arguments
!INTEGER(2), INTENT(IN) :: Unit !Unit flag: 1=SI; 2=IP
!REAL, INTENT(INOUT) :: CompPAR(26) !Compressor model input data
!REAL, INTENT(INOUT) :: CondPAR(45) !Condenser model real number input data  !RS: Debugging: Formerly CondPAR(61)
!REAL, INTENT(INOUT) :: EvapPAR(39) !Evaporator model real number input data !RS: Debugging: Formerly EvapPAR(54)
!REAL, INTENT(INOUT) :: ShTbPAR(5)  !Short tube model input data
!REAL, INTENT(INOUT) :: CapTubePAR(5) !Capillary tube model input data
!REAL, INTENT(INOUT) :: TxvPAR(7)   !TXV model input data   !RS: Debugging: Not ever used
!REAL, INTENT(INOUT) :: AccumPAR(10) !Accumulator input data
!REAL, INTENT(INOUT) :: FilterPAR(2) !Filter drier input data
!REAL, INTENT(INOUT) :: CFMcnd      !Condenser inlet air flow rate, kg/s
!REAL, INTENT(INOUT) :: CFMevp      !Evaporator inlet air flow rate, kg/s
!REAL, INTENT(INOUT) :: TaiC      !Condenser inlet air DB temp, F
!REAL, INTENT(INOUT) :: TaiE      !Evaporator inlet air DB temp, F
!REAL, INTENT(INOUT) :: RHiC      !Condenser inlet air RH
!REAL, INTENT(INOUT) :: RHiE      !Evaporator inlet air RH
!REAL, INTENT(INOUT) :: RefChg    !Refrigerant charge, Lbm
!REAL, INTENT(INOUT) :: TSOCMP    !High side saturation temp., F
!REAL, INTENT(INOUT) :: TSICMP    !Low side saturation temp., F
!REAL, INTENT(INOUT) :: SUPER     !Superheat, F
!REAL, INTENT(INOUT) :: SUBCOOL   !Subcooling, F
!REAL, INTENT(INOUT) :: BaroPressure !Barometric pressure, kPa
!REAL, INTENT(INOUT) ::  Tdis !Discharge temperature, C !RS: Debugging: Only used for CoilOnlySim, which isn't used by us
!REAL, INTENT(INOUT) ::  Tliq !Liquid temperature, C    !RS: Debugging: Only used for CoilOnlySim, which isn't used by us

  IF (Unit .EQ. SI)THEN !SI unit inputs
    
	!CompPAR%CompIntVol=CompPAR%CompIntVol/(100**3) !Compressor internal volume, m^3   !RS: Formerly CompPAR(23)
 !
	!****Condenser input data****
    !CondPAR(1)                   !Discharge line length, m
	CondPAR%CondDisLnOD=CondPAR%CondDisLnOD/1000   !Discharge line outside diameter, m    !RS: Debugging: Formerly CondPAR(2)
	CondPAR%CondDisLnTWThick=CondPAR%CondDisLnTWThick/1000   !Discharge line tube wall thickness, m !RS: Debugging: Formerly CondPAR(3)
	!CondPAR(4)                   !Discharge line elevation, m
	CondPAR%CondDisLnQLoss=CondPAR%CondDisLnQLoss/1000   !Discharge line heat loss, kW  !RS: Debugging: Formerly CondPAR(5)
	!CondPAR(6)                   !Discharge line temperature drop, C
	!CondPAR(7)                   !Discharge line additional pressure drop, kPa
	!CondPAR(8)                   !Liquid line length, m
	CondPAR%CondLiqLnOD=CondPAR%CondLiqLnOD/1000   !Liquid line outside diameter, m   !RS: Debugging: Formerly CondPAR(9)
	CondPAR%CondLiqLnTWThick=CondPAR%CondLiqLnTWThick/1000   !Liquid line tube wall thickness, m  !RS: Debugging: Formerly CondPAR(10)
	!CondPAR(11)                  !Liquid line elevation, m
	CondPAR%CondLiqLnQLoss=CondPAR%CondLiqLnQLoss/1000 !Liquid line heat loss, kW !RS: Debugging: Formerly CondPAR(12)
	!CondPAR(13)                  !Liquid line temperature drop, C
	!CondPAR(14)                  !Liquid line additional pressure drop, kPa
	CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/1000   !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
    CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick/1000   !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
    !CondPAR(17)                   !Tube length, m 
    !CondPAR(18)                   !Tube thermal conductivity, kW/m-C 
    CondPAR%CondTspc=CondPAR%CondTspc/1000   !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
    CondPAR%CondRspc=CondPAR%CondRspc/1000   !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
    CondPAR%CondFinThick=CondPAR%CondFinThick/1000   !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
    !CondPAR(22)                   !Fin pitch, fin/m 
    !CondPAR(23)                   !Fin thermal conductivity, kW/m-C 
    !CondPAR(24)                   !Number of tubes in transverse direction (normal to air flow) 
    !CondPAR(25)                   !Number of tubes in longitudinal direction (parallel to air flow) 
    !CondPAR(26)                   !Number of circuits 
    !CondPAR(27)                   !Equivalent circuits (1=yes; 2=no) 
    !CondPAR(28)                   !Number of modules per tube 
    !CondPAR(29)                   !Fin type (1-smooth; 2-Wavy; 3-louvered) 
	CondPAR%CondFanPwr=CondPAR%CondFanPwr/1000   !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)

	!CondPAR(38)                   !Barometric pressure, kPa

	!****Evaporator input data****
	!EvapPAR(1)                   !Suction line length, m
	EvapPAR%EvapSucLnOD=EvapPAR%EvapSucLnOD/1000   !Suction line outside diameter, m  !RS: Debugging: Formerly EvapPAR(2)
	EvapPAR%EvapSucLnTWThick=EvapPAR%EvapSucLnTWThick/1000   !Suction line tube wall thickness, m   !RS: Debugging: Formerly EvapPAR(3)
	!EvapPAR(4)                   !Suction line elevation, m
	EvapPAR%EvapSucLnQLoss=EvapPAR%EvapSucLnQLoss/1000   !Suction line heat gain, kW    !RS: Debugging: Formerly EvapPAR(5)
	!EvapPAR(6)                   !Suction line temperature rise, C
	!EvapPAR(7)                   !Suction line additional pressure drop, kPa
    EvapPAR%EvapCoilTOD=EvapPAR%EvapCoilTOD/1000   !Tube outside diameter, m  !RS: Debugging: Formerly EvapPAR(8)
    EvapPAR%EvapCoilTWThick=EvapPAR%EvapCoilTWThick/1000   !Tube wall thickness, m    !RS: Debugging: Formerly EvapPAR(9)
    !EvapPAR(10)                  !Tube length, m 
    !EvapPAR(11)                 !Tube thermal conductivity, kW/m-C 
    EvapPAR%EvapTspc=EvapPAR%EvapTspc/1000 !Tube spacing in transverse direction (normal to air flow), m  !RS: Debugging: Formerly EvapPAR(12)
    EvapPAR%EvapRspc=EvapPAR%EvapRspc/1000 !Tube spacing in longitudinal direction (parallel to air flow), m  !RS: Debugging: Formerly EvapPAR(13)
    EvapPAR%EvapFinThick=EvapPAR%EvapFinThick/1000 !Fin thickness, m  !RS: Debugging: Formerly EvapPAR(14)
    !EvapPAR(15)                 !Fin pitch, fin/m 
    !EvapPAR(16)                 !Fin thermal conductivity, kW/m-C 
    !EvapPAR(17)                 !Number of tubes in transverse direction (normal to air flow) 
    !EvapPAR(18)                 !Number of tubes in longitudinal direction (parallel to air flow) 
    !EvapPAR(19)                 !Number of circuits 
    !EvapPAR(20)                 !Equivalent circuits (1=yes; 2=no) 
    !EvapPAR(21)                 !Number of modules per tube 
    !EvapPAR(22)                 !Fin type (1-smooth; 2-Wavy; 3-louvered) 
	EvapPAR%EvapFanPwr=EvapPAR%EvapFanPwr/1000 !Fan power, kW !RS: Debugging: Formerly EvapPAR(27)

	!EvapPAR(31)                   !Barometric pressure, kPa

    !****Short tube input data****
    ShTbPAR%ShTbTLen=ShTbPAR%ShTbTLen/1000 !Length, m   !RS: Debugging: Formerly ShTbPAR(1)
    ShTbPAR%ShTbTID=ShTbPAR%ShTbTID/1000 !Diameter, m !RS: Debugging: Formerly ShTbPAR(2)
    ShTbPAR%ShTbChamDep=ShTbPAR%ShTbChamDep/1000 !45 deg chamfer depth, m !RS: Debugging: Formerly ShTbPAR(3)
    !ShTbPAR(4) !Number of circuits
	!ShTbPAR(5) !Distributor tube length, m

    !****Cap. tube input data****
    CapTubePAR%CTTubeID=CapTubePAR%CTTubeID/1000 !Length, m !RS: Debugging: Formerly CapTubePAR(1)
    CapTubePAR%CTTubeLen=CapTubePAR%CTTubeLen/1000 !Diameter, m   !RS: Debugging: Formerly CapTubePAR(2)
    CapTubePAR%CTTubeCoilD=CapTubePAR%CTTubeCoilD/1000 !Coil diameter, m  !RS: Debugging: Formerly CapTubePAR(3)
    !CapTubePAR(4) !Number of circuits
	!CapTubePAR(5) !Distributor tube length, m

    !****TXV input data****
    !TxvPAR(1) !Rated TXV capacity, ton
    !TxvPAR(2) !Rated superheat, C
    !TxvPAR(3) !Static superheat, C
    !TxvPAR(4) !Bleed factor
    !TxvPAR(5) !Number of circuits in evaporator
	!TxvPAR(6) !Distributor tube length, m
	!TxvPAR(7) !Maximum effective superheat, C

	!***** Accumulator input data *****
	AccumPAR%AccD = AccumPAR%AccD/1000 !Accumulator inside diameter, m  !RS: Debugging: Formerly AccumPAR(1)
	AccumPAR%AccH = AccumPAR%AccH/1000 !Accumulator internal height, m  !RS: Debugging: Formerly AccumPAR(2)
	AccumPAR%AccD1 = AccumPAR%AccD1/1000 !J-tube lower hole diameter, m   !RS: Debugging: Formerly AccumPAR(3)
	AccumPAR%AccD2 = AccumPAR%AccD2/1000 !J-tube upper hole diameter, m   !RS: Debugging: Formerly AccumPAR(4)
	AccumPAR%AccHDis = AccumPAR%AccHDis/1000 !Distance between holes on J-tube, m !RS: Debugging: Formerly AccumPAR(5)
	AccumPAR%AccDTube = AccumPAR%AccDTube/1000 !J-tube inside diameter, m   !RS: Debugging: Formerly AccumPAR(6)
	!AccumPAR(7) = Rating pressure drop, kPa
	!AccumPAR(8) = Rating temperature drop, K
	!AccumPAR(9) = Curve fit coefficient M
	!AccumPAR(10) = Curve fit coefficient B

	!Filter drier input data
	!FilterPAR(1)=Flow capacity, ton
	!FilterPAR(2)=Rating pressure drop, kPa

    !****Refrigerant input data****
    RefChg=RefChg/UnitM !Refrigerant charge lbm, ORNL solver uses IP unit

    !Air side boundary conditions
    !XMaC             !Condenser inlet air flow rate, kg/s 
    TaiC= TaiC*1.8+32 !Condenser inlet DB temperature F, ORNL solver uses IP unit
    !RHiC             !Condenser inlet RH 
    !XMaE             !Evaporator inlet air flow rate, kg/s 
    TaiE=TaiE*1.8+32  !Evaporator inlet DB temperature F, ORNL solver uses IP unit
    !RHiE             !Evaporator inlet RH 

    !Initial guesses
    TsoCmp=TsoCmp*1.8+32 !High side saturation temp. F, ORNL solver uses IP unit
    TsiCmp=TsiCmp*1.8+32 !Low side saturation temp. F, ORNL solver uses IP unit
    Super=Super*1.8      !Superheat, F, ORNL solver uses IP unit
	Subcool=Subcool*1.8  !Subcooling, F, ORNL solver uses IP unit
	!BaroPressure !barometric pressure kPa

  ELSE !IP unit inputs
  
    !****Compressor input data****
    !CompPAR(1-10) !10 coefficients for power - DO nothing
  
    !10 coefficients for mass flow rate
    !CompPAR(11-20) !DO - nothing
      
    !Compressor shell heat loss
	!CompPAR%CompQLoss=CompPAR%CompQLoss*UnitPwr*1000 !Compressor shell heat loss W  !RS: Debugging: Formerly CompPAR(22)
   
	!CompPAR%CompIntVol=CompPAR%CompIntVol/(12**3)*(UnitL**3) !Compressor internal volume, m^3 !RS: Debugging: Formerly CompPAR(23)

	!****Condenser input data****
	CondPAR%CondDisLnLen=CondPAR%CondDisLnLen*UnitL            !Discharge line length, m    !RS: Debugging: Formerly CondPAR(1)
	CondPAR%CondDisLnOD=CondPAR%CondDisLnOD/12*UnitL         !Discharge line outside diameter, m  !RS: Debugging: Formerly CondPAR(2)
	CondPAR%CondDisLnTWThick=CondPAR%CondDisLnTWThick*0.001/12*UnitL   !Discharge line tube wall thickness, m   !RS: Debugging: Formerly CondPAR(3)
	CondPAR%CondDisLnElev=CondPAR%CondDisLnElev*UnitL            !Discharge line elevation, m !RS: Debugging: Formerly CondPAR(4)
	CondPAR%CondDisLnQLoss=CondPAR%CondDisLnQLoss*UnitPwr          !Discharge line heat loss, kW    !RS: Debugging: Formerly CondPAR(5)
	CondPAR%CondDisLnTempChg=CondPAR%CondDisLnTempChg/1.8              !Discharge line temperature drop, C  !RS: Debugging: Formerly CondPAR(6)
	CondPAR%CondDisLnAddPD=CondPAR%CondDisLnAddPD*UnitP            !Discharge line additional pressure drop, kPa    !RS: Debugging: Formerly CondPAR(7)
	CondPAR%CondLiqLnLen=CondPAR%CondLiqLnLen*UnitL          !Liquid line length, m !RS: Debugging: Formerly CondPAR(8)
	CondPAR%CondLiqLnOD=CondPAR%CondLiqLnOD/12*UnitL       !Liquid line outside diameter, m   !RS: Debugging: Formerly CondPAR(9)
	CondPAR%CondLiqLnTWThick=CondPAR%CondLiqLnTWThick*0.001/12*UnitL !Liquid line tube wall thickness, m  !RS: Debugging: Formerly CondPAR(10)
	CondPAR%CondLiqLnElev=CondPAR%CondLiqLnElev*UnitL          !Liquid line elevation, m    !RS: Debugging: Formerly CondPAR(11)
	CondPAR%CondLiqLnQLoss=CondPAR%CondLiqLnQLoss*UnitPwr        !Liquid line heat loss, kW   !RS: Debugging: Formerly CondPAR(12)
	CondPAR%CondLiqLnTempChg=CondPAR%CondLiqLnTempChg/1.8            !Liquid line temperature drop, C !RS: Debugging: Formerly CondPAR(13)
    CondPAR%CondLiqLnAddPD=CondPAR%CondLiqLnAddPD*UnitP          !Liquid line additional pressure drop, kPa   !RS: Debugging: Formerly CondPAR(14)
	CondPAR%CondCoilTOD=CondPAR%CondCoilTOD/12*UnitL       !Tube outside diameter, m    !RS: Debugging: Formerly CondPAR(15)
    CondPAR%CondCoilTWThick=CondPAR%CondCoilTWThick*0.001/12*UnitL !Tube wall thickness, m  !RS: Debugging: Formerly CondPAR(16)
    CondPAR%CondCoilSTLen=CondPAR%CondCoilSTLen/12*UnitL       !Tube length, m  !RS: Debugging: Formerly CondPAR(17)
    CondPAR%CondCoilTThermCon=CondPAR%CondCoilTThermCon*UnitK          !Tube thermal conductivity, kW/m-C   !RS: Debugging: Formerly CondPAR(18)
    CondPAR%CondTspc=CondPAR%CondTspc/12*UnitL       !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly CondPAR(19)
    CondPAR%CondRspc=CondPAR%CondRspc/12*UnitL       !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly CondPAR(20)
    CondPAR%CondFinThick=CondPAR%CondFinThick*0.001/12*UnitL !Fin thickness, m    !RS: Debugging: Formerly CondPAR(21)
    CondPAR%CondFinPitch=CondPAR%CondFinPitch*12/UnitL       !Fin pitch, fin/m    !RS: Debugging: Formerly CondPAR(22)
    CondPAR%CondFinThermCon=CondPAR%CondFinThermCon*UnitK          !Fin thermal conductivity, kW/m-C    !RS: Debugging: Formerly CondPAR(23)
    !CondPAR(24)                           !Number of tubes in transverse direction (normal to air flow) 
    !CondPAR(25)                           !Number of tubes in longitudinal direction (parallel to air flow) 
    !CondPAR(26)                           !Number of circuits 
    !CondPAR(27)                           !Equivalent circuits (1=yes; 2=no) 
    !CondPAR(28)                           !Number of modules per tube 
    !CondPAR(29)                           !Fin type (1-smooth; 2-Wavy; 3-louvered) 
	CondPAR%CondFanPwr=CondPAR%CondFanPwr*1E-3           !Fan power, kW   !RS: Debugging: Formerly CondPAR(34)
 
	CondPAR%CondBarPress=CondPAR%CondBarPress*UnitP          !Barometric pressure, kPa    !RS: Debugging: Formerly CondPAR(38)

	!****Evaporator input data****
	EvapPAR%EvapSucLnLen=EvapPAR%EvapSucLnLen*UnitL            !Suction line length, m  !RS: Debugging: Formerly EvaPAR(1)
	EvapPAR%EvapSucLnOD=EvapPAR%EvapSucLnOD/12*UnitL         !Suction line outside diameter, m    !RS: Debugging: Formerly EvapPAR(2)
	EvapPAR%EvapSucLnTWThick=EvapPAR%EvapSucLnTWThick*0.001/12*UnitL   !Suction line tube wall thickness, m !RS: Debugging: Formerly EvapPAR(3)
    EvapPAR%EvapSucLnElev=EvapPAR%EvapSucLnElev*UnitL            !Suction line elevation, m   !RS: Debugging: Formerly EvapPAR(4)
	EvapPAR%EvapSucLnQLoss=EvapPAR%EvapSucLnQLoss*UnitPwr          !Suction line heat gain, kW  !RS: Debugging: Formerly EvapPAR(5)
	EvapPAR%EvapSucLnTempChg=EvapPAR%EvapSucLnTempChg/1.8              !Suction line temperature rise, C    !RS: Debugging: Formerly EvapPAR(6)
	EvapPAR%EvapSucLnAddPD=EvapPAR%EvapSucLnAddPD*UnitP            !Suction line additional pressure drop, kPa  !RS: Debugging: Formerly EvapPAR(7)
    EvapPAR%EvapCoilTOD=EvapPAR%EvapCoilTOD/12*UnitL         !Tube outside diameter, m    !RS: Debugging: Formerly EvapPAR(8)
    EvapPAR%EvapCoilTWThick=EvapPAR%EvapCoilTWThick*0.001/12*UnitL   !Tube wall thickness, m  !RS: Debugging: Formerly EvapPAR(9)
    EvapPAR%EvapCoilSTLen=EvapPAR%EvapCoilSTLen/12*UnitL       !Tube length, m  !RS: Debugging: Formerly EvapPAR(10)
    EvapPAR%EvapCoilTThermCon=EvapPAR%EvapCoilTThermCon*UnitK          !Tube thermal conductivity, kW/m-C   !RS: Debugging: Formerly EvapPAR(11)
    EvapPAR%EvapTspc=EvapPAR%EvapTspc/12*UnitL       !Tube spacing in transverse direction (normal to air flow), m    !RS: Debugging: Formerly EvapPAR(12)
    EvapPAR%EvapRspc=EvapPAR%EvapRspc/12*UnitL       !Tube spacing in longitudinal direction (parallel to air flow), m    !RS: Debugging: Formerly EvapPAR(13)
    EvapPAR%EvapFinThick=EvapPAR%EvapFinThick*0.001/12*UnitL !Fin thickness, m    !RS: Debugging: Formerly EvapPAR(14)
    EvapPAR%EvapFinPitch=EvapPAR%EvapFinPitch*12/UnitL       !Fin pitch, fin/m    !RS: Debugging: Formerly EvapPAR(15)
    EvapPAR%EvapFinThermCon=EvapPAR%EvapFinThermCon*UnitK          !Fin thermal conductivity, kW/m-C    !RS: Debugging: Formerly EvapPAR(16)
    !EvapPAR(17)                           !Number of tubes in transverse direction (normal to air flow) 
    !EvapPAR(18)                           !Number of tubes in longitudinal direction (parallel to air flow) 
    !EvapPAR(19)                           !Number of circuits 
    !EvapPAR(20)                           !Equivalent circuits (1=yes; 2=no) 
    !EvapPAR(21)                           !Number of modules per tube 
    !EvapPAR(22)                           !Fin type (1-smooth; 2-Wavy; 3-louvered) 
    EvapPAR%EvapFanPwr=EvapPAR%EvapFanPwr*1E-3           !Fan power, kW   !RS: Debugging: EvapPAR(27)

	EvapPAR%EvapBarPress=EvapPAR%EvapBarPress*UnitP          !Barometric pressure, kPa    !RS: Debugging: EvapPAR(31)

    !****Short tube input data****
    ShTbPAR%ShTbTLen=ShTbPAR%ShTbTLen/12*UnitL !Length, m   !RS: Debugging: Formerly ShTbPAR(1)
    ShTbPAR%ShTbTID=ShTbPAR%ShTbTID/12*UnitL !Diameter, m !RS: Debugging: Formerly ShTbPAR(2)
    ShTbPAR%ShTbChamDep=ShTbPAR%ShTbChamDep/12*UnitL !45 deg chamfer depth, m !RS: Debugging: Formerly ShTbPAR(3)
    !ShTbPAR(4) !Number of circuits
	ShTbPAR%ShTbDTubeLen=ShTbPAR%ShTbDTubeLen*UnitL !Distributor tube length, m !RS: Debugging: Formerly ShTbPAR(5)

    !****Cap. tube input data****
    CapTubePAR%CTTubeLen=CapTubePAR%CTTubeLen/12*UnitL !Length, m !RS: Debugging: Formerly CapTubePAR(1)
    CapTubePAR%CTTubeID=CapTubePAR%CTTubeID/12*UnitL !Diameter, m   !RS: Debugging: Formerly CapTubePAR(2)
    CapTubePAR%CTTubeCoilD=CapTubePAR%CTTubeCoilD/12*UnitL !Coil diameter, m  !RS: Debugging: Formerly CapTubePAR(3)
    !CapTubePAR(4) !Number of circuits
	CapTubePAR%CTDisTubeLen=CapTubePAR%CTDisTubeLen*UnitL !Distributor tube length, m   !RS: Debugging: Formerly CapTubePAR(5)

    !****TXV input data****
    !TxvPAR(1)                !Rated TXV capacity, ton
    !TxvPAR(2)=TxvPAR(2)/1.8   !Rated superheat, C  !RS: Debugging: Not ever used
    !TxvPAR(3)=TxvPAR(3)/1.8   !Static superheat, C !RS: Debugging: Not ever used
    !TxvPAR(4)                !Bleed factor
    !TxvPAR(5)                !Number of circuits in evaporator
	!TxvPAR(6)=TxvPAR(6)*UnitL !Distributor tube length, m  !RS: Debugging: Not ever used
	!TxvPAR(7)=TxvPAR(7)/1.8   !Maximum effective superheat, C  !RS: Debugging: Not ever used

	!***** Accumulator input data *****
	AccumPAR%AccD = AccumPAR%AccD/12*0.3048 !Accumulator inside diameter, m !RS: Debugging: Formerly AccumPAR(1)
	AccumPAR%AccH = AccumPAR%AccH/12*0.3048 !Accumulator internal height, m !RS: Debugging: Formerly AccumPAR(2)
	AccumPAR%AccD1 = AccumPAR%AccD1/12*0.3048 !J-tube lower hole diameter, m  !RS: Debugging: Formerly AccumPAR(3)
	AccumPAR%AccD2 = AccumPAR%AccD2/12*0.3048 !J-tube upper hole diameter, m  !RS: Debugging: Formerly AccumPAR(4)
	AccumPAR%AccHDis = AccumPAR%AccHDis/12*0.3048 !Distance between holes on J-tube, m    !RS: Debugging: Formerly AccumPAR(5)
	AccumPAR%AccDTube = AccumPAR%AccDTube/12*0.3048 !J-tube inside diameter, m  !RS: Debugging: Formerly AccumPAR(6)
	AccumPAR%AccDP = AccumPAR%AccDP*UnitP !Rating pressure drop, kPa  !RS: Debugging: Formerly AccumPAR(7)
	AccumPAR%AccDT = AccumPAR%AccDT/1.8 !Rating temperature drop, K   !RS: Debugging: Formerly AccumPAR(8)
	!AccumPAR(9) = Curve fit coefficient M
	!AccumPAR(10) = Curve fit coefficient B

	!***** Filter drier input data *****
	!FilterPAR(1)=Flow capacity, ton
	!FilterPAR(2)=FilterPAR(2)*UnitP !Rating pressure drop, kPa

    !****Refrigerant input data****
    !RefChg !Refrigerant charge lbm, , ORNL solver uses IP unit

    !Air side boundary conditions
	CFMcnd=CFMcnd*UnitArFlw !XMaC=XMaC*UnitArFlw                 !Condenser inlet air flow rate, m^3/s
    !TaiC                               !Condenser inlet DB temperature F,  ORNL solver uses IP unit
    IF (RHiC .GT. 1) THEN
        RHiC=(RHiC-32)/1.8 !Condenser inlet RH or WB temp
    END IF
    CFMevp=CFMevp*UnitArFlw !XMaE=XMaE*UnitArFlw                 !Evaporator inlet air flow rate, m^3/s
    !TaiE                               !Evaporator inlet DB temperature F,  ORNL solver uses IP unit
    IF (RHiE .GT. 1) THEN
        RHiE=(RHiE-32)/1.8 !Evaporator inlet RH or WB temp
    END IF

    !Initial guesses
	BaroPressure=BaroPressure*UnitP !barometric pressure kPa

  END IF

  RETURN

END SUBROUTINE UnitConvert

!***********************************************************************************

SUBROUTINE MicroChannelCoilUnitConvert(IsSIUnit,FinPitch,Kfin,FinThk, &
								       TubeHeight,TubeDepth,TubeThk,Ktube, &
									   Pt,Ltube,Dchannel)

! ----------------------------------------------------------------------
!
!   Description: To convert microchannel coil input data into standard SI
!                Some data needs to be in IP unit for the
!                Corresponding model
!
!   Author:
!   Ipseng Iu
!   Mechanical and Aerospace Engineering
!   Oklahoma State University, Stillwater
!
!   Date:
!   November, 2005
!
! ----------------------------------------------------------------------

IMPLICIT NONE

LOGICAL,INTENT(IN)    :: IsSIunit    !SI unit flag
REAL,   INTENT(INOUT) :: FinPitch    !Fin pitch, fins/m
REAL,   INTENT(INOUT) :: Kfin        !Fin conductivity, kW/m-C 
REAL,   INTENT(INOUT) :: FinThk      !Fin thickness, m
REAL,   INTENT(INOUT) :: TubeHeight  !Tube height, m
REAL,   INTENT(INOUT) :: TubeDepth   !Tube depth, m
REAL,   INTENT(INOUT) :: TubeThk     !Tube thickness, m
REAL,   INTENT(INOUT) :: Ktube       !Tube conductivity, kW/m-C
REAL,   INTENT(INOUT) :: Pt          !Tube spacing, m
REAL,   INTENT(INOUT) :: Ltube       !Single tube length, m
REAL,   INTENT(INOUT) :: Dchannel    !Channel diameter, m

!FLOW:

  IF (IsSIunit) THEN !SI unit inputs
	  !FinPitch
	  Kfin       =Kfin*1e-3       !W/m-C to kW/m-C
	  FinThk     =FinThk*1e-3     !mm to m
	  TubeHeight =TubeHeight*1e-3 !mm to m
	  TubeDepth  =TubeDepth*1e-3  !mm to m
	  TubeThk    =TubeThk*1e-3    !mm to m
	  Ktube      =Ktube*1e-3      !W/m-C to kW/m-C
	  Pt         =Pt*1e-3         !mm to m
	  Ltube      =Ltube*1e-3      !mm to m
	  Dchannel   =Dchannel*1e-2   !mm to m
  ELSE
	  FinPitch   =FinPitch*12/UnitL     !fins/in to fins/m
	  Kfin       =Kfin*UnitK            !Btu-in/hr-ft2-F to kW/m-C
	  FinThk     =FinThk*0.001/12*UnitL !mil to m
	  TubeHeight =TubeHeight/12*UnitL   !in to m
	  TubeDepth  =TubeDepth/12*UnitL    !in to m
	  TubeThk    =TubeThk/12*UnitL      !in to m
	  Ktube      =Ktube*UnitK           !Btu-in/hr-ft2-F to kW/m-C 
	  Pt         =Pt/12*UnitL           !in to m
	  Ltube      =Ltube/12*UnitL        !in to m 
	  Dchannel   =Dchannel/12*UnitL     !in to m
  END IF

RETURN

END SUBROUTINE MicroChannelCoilUnitConvert

!***********************************************************************************

SUBROUTINE FinTubeCoilUnitConvert(IsSIUnit,FinPitch,Kfin,FinThk, &
								  ODtube,IDtube,Ktube,Pt,Pl,Ltube)
									  
! ----------------------------------------------------------------------
!
!   Description: To convert fin-tube coil input data into standard SI
!                Some data needs to be in IP unit for the
!                Corresponding model
!
!   Author:
!   Ipseng Iu
!   Mechanical and Aerospace Engineering
!   Oklahoma State University, Stillwater
!
!   Date:
!   November, 2005
!
! ----------------------------------------------------------------------

IMPLICIT NONE

LOGICAL,INTENT(IN)    :: IsSIunit    !SI unit flag
REAL,   INTENT(INOUT) :: FinPitch    !Fin pitch, fins/m
REAL,   INTENT(INOUT) :: Kfin        !Fin conductivity, kW/m-C 
REAL,   INTENT(INOUT) :: FinThk      !Fin thickness, m
REAL,   INTENT(INOUT) :: ODtube      !Tube height, m
REAL,   INTENT(INOUT) :: IDtube      !Tube depth, m
REAL,   INTENT(INOUT) :: Ktube       !Tube conductivity, kW/m-C
REAL,   INTENT(INOUT) :: Pt          !Tube spacing, m
REAL,   INTENT(INOUT) :: Pl          !Row spacing, m
REAL,   INTENT(INOUT) :: Ltube       !Single tube length, m

!FLOW:

  IF (IsSIunit) THEN !SI unit inputs
	  !FinPitch
	  Kfin       =Kfin*1e-3       !W/m-C to kW/m-C
	  FinThk     =FinThk*1e-3     !mm to m
	  ODtube     =ODtube*1e-3     !mm to m
	  IDtube     =IDtube*1e-3     !mm to m
	  Ktube      =Ktube*1e-3      !W/m-C to kW/m-C
	  Pt         =Pt*1e-3         !mm to m
	  Pl         =Pl*1e-3         !mm to m
	  Ltube      =Ltube*1e-3      !mm to m
  ELSE
	  FinPitch   =FinPitch*12/UnitL     !fins/in to fins/m
	  Kfin       =Kfin*UnitK            !Btu-in/hr-ft2-F to kW/m-C
	  FinThk     =FinThk*0.001/12*UnitL !mil to m
	  ODtube     =ODtube/12*UnitL       !in to m
	  IDtube     =IDtube/12*UnitL       !in to m
	  Ktube      =Ktube*UnitK           !Btu-in/hr-ft2-F to kW/m-C 
	  Pt         =Pt/12*UnitL           !in to m
	  Pl         =Pl/12*UnitL           !in to m
	  Ltube      =Ltube/12*UnitL        !in to m 
  END IF

RETURN

END SUBROUTINE FinTubeCoilUnitConvert

!***********************************************************************************

REAL FUNCTION Temperature_F2C(tF) 
	
!   Description: Convert temperature from Fahrenheit to Centigrade
!   Author: Venu Lolla

    IMPLICIT NONE
    REAL, INTENT(IN) :: tF

    Temperature_F2C = (tF - 32.0) / 1.8

END FUNCTION Temperature_F2C

!***********************************************************************************

REAL FUNCTION Temperature_C2F(tC) 
	
!   Description: Convert temperature from Centigrade to Fahrenheit
!   Author: Venu Lolla

    IMPLICIT NONE
    REAL, INTENT(IN) :: tC

    Temperature_C2F = tC * 1.8 + 32.0

END FUNCTION Temperature_C2F

!***********************************************************************************
END MODULE UnitConvertMod
