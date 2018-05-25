! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This is a module that models a reversing valve for the Heat Pump system
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This represents a reversing valve in a heat pump system
! A description of the component is found at:
! http://en.wikipedia.org/wiki/Reversing_valve
! From that website: 
!  - It changes the direction of refrigerant flow

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module takes the suction side properties and calculates heat transfer and pressure drop through the valve 

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! no input/output files

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! nothing defined at the module level

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 2 methods:
!    PUBLIC SuctionHeatTransfer - calculates heat transfer on the suction side
!       Called by Evaporator
!    PUBLIC SuctionPressureDrop - calculates the pressure drop on the suction side
!       Called by Evaporator

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! No current issues/bugs/tickets

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-29 | JEH | Header Completed

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Nothing currently needs to be done on this module
!
MODULE ReversingValveMod

IMPLICIT NONE

!------------------------------------------------------------------------------
!
!The reversing valve model from York empirical data
!Reference:
!Damasceno G D S, Lee W N T, Rooke S P., Goldschmidt V W. 
!"Performance of a heat pump reversing valves and comparison through 
!characterising parameters." ASHRAE Transactions, 1988, Vol 94, Part 1, 
!pp. 304-317
!
!Written by IPSENG IU
!July 2007
!
!------------------------------------------------------------------------------

PUBLIC SuctionHeatTransfer
PUBLIC SuctionPressureDrop

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE SuctionHeatTransfer(mdot,Tdis,Tsuc,hsi,hso)

!Calculate heat gain at suction side

IMPLICIT NONE

REAL, INTENT(IN) :: mdot !refrigerant mass flow rate, kg/s
REAL, INTENT(IN) :: Tdis !Discharge temperature, C
REAL, INTENT(IN) :: Tsuc !Suction temperature, C
REAL, INTENT(IN) :: hsi !Suction inlet enthalpy, kJ/kg
REAL, INTENT(OUT) :: hso !Suction outlet enthalpy, kJ/kg

REAL :: UA !Overall heat transfer coefficient, kW/K

UA = 0.1541*mdot - 0.002    !RS Comment: What are the values from?
IF (UA .LT. 0) THEN
    UA=0
END IF

hso=UA*(Tdis-Tsuc)/mdot+hsi

hso = hsi   !RS: Comment: Added by John Gall (1/9/14) for the cooling-only system
!RS: This may not be necessary; it may be possible to just not call these routines (1/10/14)

END SUBROUTINE SuctionHeatTransfer

!------------------------------------------------------------------------------

SUBROUTINE SuctionPressureDrop(mdot,DP)

!Calculate suction side pressure drop

IMPLICIT NONE

REAL, INTENT(IN) :: mdot !Refrigerant mass flow rate, kg/s
REAL, INTENT(OUT) :: DP !Pressure drop kPa

DP = 327.69*mdot + 3.9633   !RS Comment: What are the values from?

IF (DP .LT. 0) THEN
    DP = 0
END IF

DP = 0  !RS: Comment: Added by John Gall (1/9/14) for the cooling-only system
!RS: This may not be necessary; it may be possible to just not call these routines (1/10/14)

END SUBROUTINE SuctionPressureDrop

!------------------------------------------------------------------------------

END MODULE ReversingValveMod
