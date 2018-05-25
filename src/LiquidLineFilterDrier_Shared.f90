
! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This routine calculates the pressure drop over the filter drier.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component represents the filter drier in a heat pump system.
! A description of the component is found at:
! http://www.achrnews.com/articles/a-guide-to-understanding-filter-drier-functions-and-types
! From that website: 
!  - The filter absorbs moisture
!  - The filter also filters particulates out

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! It lowers the expansion device's inlet pressure.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no associated input or output files.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There are no variables or structures defined at the module level.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method:
!    PUBLIC CalcFilterDrierDP -- Calculates the pressure drop across the filter drier
!      Called by FlowRateLoop.f90

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
! The actual pressure drop calculation could use some explanation about where the equation comes from.

SUBROUTINE CalcFilterDrierDP !(XIN,PAR,OUT)

! ----------------------------------------------------------------------
!
!   Purpose: To calculate the pressure drop across the filter drier.
!
!   Method: Catalog curve fit 
!
!	Inputs:
!		RefName=Refrigerant name
!		XIN(1) = Mass flow rate, kg/s
!
!	Parameters:
!		PAR(1)=Flow capacity, ton
!		PAR(2)=Rating pressure drop, kPa
!
!	Outputs:
!		OUT(1) = Pressure drop, kPa
!
!   Author:
!   Ipseng Iu
!   Mechanical and Aerospace Engineering
!   Oklahoma State University, Stillwater
!
!   Date: July 2005
!
! ----------------------------------------------------------------------
USE DataGlobals_HPSimIntegrated, ONLY: RefName    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
USE DataSimulation !, ONLY: FilFlowCap, FilRatDP, FIDP, FODP  !RS: Debugging: Brings variable names instead of numbers through
implicit none

!Flow:

!REAL, INTENT(IN) :: XIN(1)
!REAL, INTENT(IN) :: PAR(2)
!REAL, INTENT(OUT) :: OUT(1)

!Constants from ARI std 70 
REAL,PARAMETER :: FlowRatePerTonR22 = 0.0224 !kg/s/ton '0.0064 'kg/s/kW
REAL,PARAMETER :: FlowRatePerTonR134A = 0.02345 !kg/s/ton '0.0067 'kg/s/kW
REAL,PARAMETER :: FlowRatePerTonR407C = 0.0224 !kg/s/ton '0.0064 'kg/s/kW
REAL,PARAMETER :: FlowRatePerTonR410A = 0.021 !kg/s/ton '0.006  'kg/s/kW

REAL mdot !Refrigerant mass flow rate, kg/s
REAL FlowCapacity !Flow capacity, ton
REAL RatedDP !Rating pressure drop, kPa
REAL FlowRatePerTon !Flow rate per ton, kg/s/ton
REAL DP !Pressure drop, kPa

!Flow:

	mdot=FilterIN%FIDP !RS: Debugging: Formerly XIN(1)

	FlowCapacity=FilterPAR%FilFlowCap !RS: Debugging: Formerly PAR(1)
	RatedDP=FilterPAR%FilRatDP  !RS: Debugging: Formerly PAR(2)

	IF (FlowCapacity .GT. 0) THEN
		SELECT CASE (TRIM(RefName))
		CASE('R22')
			FlowRatePerTon = FlowRatePerTonR22
		CASE('R134A')
			FlowRatePerTon = FlowRatePerTonR134A
		CASE('R410A')
			FlowRatePerTon = FlowRatePerTonR410A
		Case('R407C')
			FlowRatePerTon = FlowRatePerTonR407C
		Case DEFAULT
			FlowRatePerTon = FlowRatePerTonR22
		END SELECT

		DP = RatedDP / (FlowCapacity * FlowRatePerTon) * mdot   !Determining Pressure Drop
	ELSE
		DP=0
	END IF

	FilterOUT%FODP=DP   !RS: Debugging: Formerly OUT(1)

RETURN

END SUBROUTINE CalcFilterDrierDP




