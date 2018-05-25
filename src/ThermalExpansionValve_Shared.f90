! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module deals with the thermal expansion valve, calculating
! the Q over it.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component represents the thermal expansion valve in a heat pump system.
! A description of the component is found at:
! http://en.wikipedia.org/wiki/Thermal_expansion_valve
! http://www.swtc.edu/ag_power/air_conditioning/lecture/expansion_valve.htm
! From those websites:
!  - TXV is a flow metering device
!  - It uses fluid temperature to control how much fluid passes through
!  - This lowers the pressure and expands the fluid into a vapor

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This calculates the Q for the TXV; this has no real purpose in the
! overall simulation as the Q is never effectively used.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no input or output files directly connected to this module.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! No module level variables or structures are defined

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 4 methods:
!   PUBLIC TXV -- Calculates and returns Qtxv
!       Called once by HPdesignMod.f90
!   PRIVATE CalcQtxv -- Does the actual Qtxv calculation
!       Called once by internal TXV routine
!   PRIVATE InterpolateQtxv -- Determines the TXV capacity
!       Called once by CalcQtxv
!   PRIVATE InterpolateCF_DP -- Determines the pressure drop correction factor
!       Called once by CalcQtxv

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! This module is essentially useless in the current version of the code;
! it only returns Qtxv, which is only ever used to be printed out.

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2013-12-17 | RAS | Filled out the header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Make this module useful again by either improving what it returns or
! by allowing the rest of the code to actually use Qtxv.

MODULE TXVMOD

	USE DataGlobals_HPSimIntegrated, ONLY: RefName    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)
	implicit none

	PUBLIC TXV
	PRIVATE CalcQtxv
	PRIVATE InterpolateQtxv
	PRIVATE InterpolateCF_DP

CONTAINS
	!***********************************************************************************

	SUBROUTINE TXV(mdot,PiCmp,PoCmp,DTsub,DTsup,DP,Qtxv)

		!-----------------------------------------------------------------------------------
		!
		!  Description:
		!  To calculate the TXV capacity based on Catalog data
		!
		!  EnergyPlus REFPROP subroutines required
		!
		!  Inputs:
		!  RefName =Refrigerant name
		!  mdot = Refrigerant mass flow rate, kg/s
		!  PiCmp = Compressor suction pressure, kPa
		!  PioCmp = Compressor discharge pressure, kPa
		!  DTsub = Subcooling, K
		!  DTsup = Superheat, K
		!  DP = Pressure drop across TXV, kPa
		!
		!  Outputs:
		!  Qtxv = TXV capacity, ton
		!  ErrorFlag = Error flag: 0-No error
		!                          1-TXV solution error
		!                          2-Refprop error
		!
		!  Reference:
		!
		!  Author:
		!  Ipseng Iu
		!  Mechanical and Aerospace Engineering
		!  Oklahoma State University, Stillwater
		!
		!  Date: July 2005
		!
		!-----------------------------------------------------------------------------------

		USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
		USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code

		IMPLICIT NONE

		!Subroutine argument declarations
		REAL, INTENT(IN)  :: mdot    !Refrigerant mass flow rate, kg/s
		REAL, INTENT(IN)  :: PiCmp   !Compressor suction pressure, kPa
		REAL, INTENT(IN)  :: PoCmp   !Compressor discharge pressure, kPa
		REAL, INTENT(IN)  :: DTsub   !Subcooling, K
		REAL, INTENT(IN)  :: DTsup   !Superheat, K
		REAL, INTENT(IN)  :: DP      !Pressure drop, kPa
		REAL, INTENT(OUT) :: Qtxv    !TXV capacity, kW

		!Subroutine local variables
		REAL Temperature,Quality,Pressure
		INTEGER(2) RefPropErr			!Error flag:1-error; 0-no error

		INTEGER ErrorFlag     !0-No error
		!1-TXV solution error
		!2-Refprop error
		REAL Tliq !Liquid temperature, C
		REAL Tevp !Evaporating temperature, C
		REAL Qcmp !Compressor capacity, ton
		REAL TsatSuc !Compressor suction saturation temperature, C
		REAL TsatDis !Compressor discharge saturation temperature, C
		REAL HiCmp !Compressor suction enthalpy, kJ/kg
		REAL HiLiq !Liquid line enthalpy, kJ/kg

		!Flow:

		ErrorFlag=0

		!Calc. compressor capacity
		Pressure=PiCmp*1000   !RS Comment: Unit Conversion
		Quality=1
		TsatSuc=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !Compressor Suction Saturation Temperature
		CALL IssueRefPropError(RefPropErr, 'TXV', 2, ErrorFlag)

		Temperature=TsatSuc+DTsup
		Pressure=PiCmp*1000   !RS Comment: Unit Conversion
		HiCmp=TP(RefName, Temperature, Pressure, 'enthalpy', RefrigIndex,RefPropErr)  !Compressor Suction Enthalpy
		CALL IssueRefPropError(RefPropErr, 'TXV', 2, ErrorFlag)
		HiCmp=HiCmp/1000  !RS Comment: Unit Conversion

		Pressure=PoCmp*1000   !RS Comment: Unit Conversion
		Quality=0
		TsatDis=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !Compressor Discharge Saturation Temperature
		CALL IssueRefPropError(RefPropErr, 'TXV', 2, ErrorFlag)

		Temperature=TsatDis-DTsub
		Pressure=PoCmp*1000   !RS Comment: Unit Conversion
		HiLiq=TP(RefName, Temperature, Pressure, 'enthalpy', RefrigIndex,RefPropErr)  !Liquid Line Enthalpy
		CALL IssueRefPropError(RefPropErr, 'TXV', 2, ErrorFlag)
		HiLiq=HiLiq/1000  !RS Comment: Unit Conversion

		Qcmp=mdot*(HiCmp-HiLiq)
		Qcmp=Qcmp*0.2843

		Tevp=TsatSuc
		Tliq=TsatDis-DTsub

		CALL CalcQtxv(Qcmp,Tevp,Tliq,DP,Qtxv)

		IF (Qtxv .LT. 0) THEN
			ErrorFlag=1
			RETURN
		END IF

		RETURN

	END SUBROUTINE TXV

	!***********************************************************************************

	SUBROUTINE CalcQtxv(Qcmp,Tevp,Tliq,DP,Qtxv)

		!To calculate TXV capacity

		IMPLICIT NONE

		!Subroutine passing variables:
		REAL, INTENT(IN)  :: Qcmp    !Compressor rated capacity, ton
		REAL, INTENT(IN)  :: Tevp    !Evaporating temperature, C
		REAL, INTENT(IN)  :: Tliq    !Liquid temperature, C
		REAL, INTENT(IN)  :: DP      !Pressure drop across TXV, kPa
		REAL, INTENT(OUT) :: Qtxv    !TXV capacity, ton

		!Subroutine local variables:
		REAL TevpMax  !Maximum reference evaporating temperature, C
		REAL TevpMin  !Minimum reference evaporating temperature, C
		REAL M_Qmax   !M-Coefficient for maximum reference capacity
		REAL M_Qmin   !M-Coefficient for minimum reference capacity
		REAL B_Qmax   !B-Coefficient for maximum reference capacity
		REAL B_Qmin   !B-Coefficient for minimum reference capacity
		REAL M_Tliq   !M-Coefficient for reference liquid temperature
		REAL B_Tliq   !B-Coefficient for reference liquid temperature
		REAL a_DPmax  !a-Coefficient for maximum reference pressure drop
		REAL a_DPmin  !a-Coefficient for minimum reference pressure drop
		REAL b_DPmax  !b-Coefficient for maximum reference pressure drop
		REAL b_DPmin  !b-Coefficient for minimum reference pressure drop
		REAL CF_Tliq  !Correction factor for liquid temperature
		REAL CF_DP    !Correction factor for pressure drop

		!Flow:

		!Get capacity and liquid temperature coefficients
		SELECT CASE (TRIM(RefName))

		CASE ('R407C')
			IF (Tevp .GT. 4.44) THEN
				TevpMax=4.44;   M_Qmax=0.9277;	 B_Qmax=0.0394
				TevpMin=-6.67;  M_Qmin=0.9921;	 B_Qmin=0.0123

			ELSEIF (Tevp .LT. 4.44 .AND. Tevp .GT. -6.67) THEN
				TevpMax=4.44;   M_Qmax=0.9277;	 B_Qmax=0.0394
				TevpMin=-6.67;  M_Qmin=0.9921;	 B_Qmin=0.0123

			ELSEIF (Tevp .LT. -6.67 .AND. Tevp .GT. -17.78) THEN
				TevpMax=-6.67;    M_Qmax=0.9921;	 B_Qmax=0.0123
				TevpMin=-17.78;   M_Qmin=0.936;	 B_Qmin=-0.3312

			ELSE
				TevpMax=-6.67;    M_Qmax=0.9921;	 B_Qmax=0.0123
				TevpMin=-17.78;   M_Qmin=0.936;	 B_Qmin=-0.3312

			END IF

			M_Tliq=-0.0127; B_Tliq=1.4727

		CASE ('R134A')
			IF (Tevp .GT. 4.44) THEN
				TevpMax=4.44;   M_Qmax=1.2053;	 B_Qmax=-0.0038
				TevpMin=-6.67;  M_Qmin=1.1961;	 B_Qmin=-0.112

			ELSEIF (Tevp .LT. 4.44 .AND. Tevp .GT. -6.67) THEN
				TevpMax=4.44;   M_Qmax=1.2053;	 B_Qmax=-0.0038
				TevpMin=-6.67;  M_Qmin=1.1961;	 B_Qmin=-0.112

			ELSEIF (Tevp .LT. -6.67 .AND. Tevp .GT. -17.78) THEN
				TevpMax=-6.67;    M_Qmax=1.1961;	 B_Qmax=-0.112
				TevpMin=-17.78;   M_Qmin=0.9697;	 B_Qmin=0.0863

			ELSE
				TevpMax=-6.67;    M_Qmax=1.1961;	 B_Qmax=-0.112
				TevpMin=-17.78;   M_Qmin=0.9697;	 B_Qmin=0.0863

			END IF

			M_Tliq=-0.0127; B_Tliq=1.4782

		CASE DEFAULT

			IF (Tevp .GT. 4.44) THEN
				TevpMax=4.44;		M_Qmax=1.0085;	B_Qmax=0.0443
				TevpMin=-6.67;	M_Qmin=1.0975;	B_Qmin=0.0062

			ELSEIF (Tevp .LT. 4.44 .AND. Tevp .GT. -6.67) THEN
				TevpMax=4.44;		M_Qmax=1.0085;	B_Qmax=0.0443
				TevpMin=-6.67;	M_Qmin=1.0975;	B_Qmin=0.0062

			ELSEIF (Tevp .LT. -6.67 .AND. Tevp .GT. -17.78) THEN
				TevpMax=-6.67;    M_Qmax=1.0975;    B_Qmax=0.0062
				TevpMin=-17.78;   M_Qmin=1.0518;    B_Qmin=-0.3681

			ELSEIF (Tevp .LT. -17.78 .AND. Tevp .GT. -23.33) THEN
				TevpMax=-17.78;   M_Qmax=1.0518;	 B_Qmax=-0.3681
				TevpMin=-23.33;	M_Qmin=0.9187;	 B_Qmin=-0.6363

			ELSEIF (Tevp .LT. -23.33 .AND. Tevp .GT. -28.89) THEN
				TevpMax=-23.33;	M_Qmax=0.9187;	 B_Qmax=-0.6363
				TevpMin=-28.89;   M_Qmin=0.7775;	 B_Qmin=-0.625

			ELSEIF (Tevp .LT. -28.89 .AND. Tevp .GT. -40.0) THEN
				TevpMax=-28.89;   M_Qmax=0.7775;	 B_Qmax=-0.625
				TevpMin=-40.0;	M_Qmin=0.6189;	 B_Qmin=-0.5324

			ELSE
				TevpMax=-28.89;   M_Qmax=0.7775;	 B_Qmax=-0.625
				TevpMin=-40.0;	M_Qmin=0.6189;	 B_Qmin=-0.5324

			END IF

			M_Tliq=-0.0103; B_Tliq=1.3861

		END SELECT

		CALL InterpolateQtxv(M_Qmin,B_Qmin,TevpMin,M_Qmax,B_Qmax,TevpMax,Qcmp,Tevp,Qtxv)
		CF_Tliq=M_Tliq*Tliq+B_Tliq

		!Get pressure drop coefficients
		SELECT CASE (TRIM(RefName))

		CASE ('R134A')
			IF (Tevp .GT. 4.44) THEN
				TevpMax=4.44;   a_DPmax=0.0502;	 b_DPmax=0.4968
				TevpMin=-6.67;  a_DPmin=0.0432;	 b_DPmin=0.4975

			ELSEIF (Tevp .LT. 4.44 .AND. Tevp .GT. -6.67) THEN
				TevpMax=4.44;   a_DPmax=0.0502;	 b_DPmax=0.4968
				TevpMin=-6.67;  a_DPmin=0.0432;	 b_DPmin=0.4975

			ELSEIF (Tevp .LT. -6.67 .AND. Tevp .GT. -17.78) THEN
				TevpMax=-6.67;    a_DPmax=0.0432;	 b_DPmax=0.4975
				TevpMin=-17.78;   a_DPmin=0.0432;	 b_DPmin=0.4975

			ELSE
				TevpMax=-6.67;    a_DPmax=0.0432;	 b_DPmax=0.4975
				TevpMin=-17.78;   a_DPmin=0.0432;	 b_DPmin=0.4975

			END IF

		CASE DEFAULT
			IF (Tevp .GT. 4.44) THEN
				TevpMax=4.44;		a_DPmax=0.0389; b_DPmax=0.4971
				TevpMin=-6.67;	a_DPmin=0.034;  b_DPmin=0.5002

			ELSEIF (Tevp .LT. 4.44 .AND. Tevp .GT. -6.67) THEN
				TevpMax=4.44;		a_DPmax=0.0389; b_DPmax=0.4971
				TevpMin=-6.67;	a_DPmin=0.034;  b_DPmin=0.5002

			ELSEIF (Tevp .LT. -6.67 .AND. Tevp .GT. -17.78) THEN
				TevpMax=-6.67;    a_DPmax=0.034;  b_DPmax=0.5002
				TevpMin=-17.78;   a_DPmin=0.034;  b_DPmin=0.5002

			ELSEIF (Tevp .LT. -17.78 .AND. Tevp .GT. -23.33) THEN
				TevpMax=-17.78;   a_DPmax=0.034;   b_DPmax=0.5002
				TevpMin=-23.33;	a_DPmin=0.0321;	 b_DPmin=0.4954

			ELSEIF (Tevp .LT. -23.33 .AND. Tevp .GT. -28.89) THEN
				TevpMax=-23.33;	a_DPmax=0.0321;	 b_DPmax=0.4954
				TevpMin=-28.89;   a_DPmin=0.0321;	 b_DPmin=0.4954

			ELSEIF (Tevp .LT. -28.89 .AND. Tevp .GT. -40.0) THEN
				TevpMax=-28.89;   a_DPmax=0.0321;	 b_DPmax=0.4954
				TevpMin=-40.0;	a_DPmin=0.0278;	 b_DPmin=0.5048

			ELSE
				TevpMax=-28.89;   a_DPmax=0.0321;	 b_DPmax=0.4954
				TevpMin=-40.0;	a_DPmin=0.0278;	 b_DPmin=0.5048

			END IF

		END SELECT

		CALL InterpolateCF_DP(a_DPmin,b_DPmin,TevpMin,a_DPmax,b_DPmax,TevpMax,DP,Tevp,CF_DP)

		!Application correction factor
		Qtxv=CF_Tliq*CF_DP*Qtxv

		RETURN

	END SUBROUTINE CalcQtxv

	!***********************************************************************************

	SUBROUTINE InterpolateQtxv(Mmin,Bmin,Tmin,Mmax,Bmax,Tmax,Qcmp,Tevp,Qtxv)

		!Interpolate or extrapolate TXV capacity

		IMPLICIT NONE

		!Subroutine passing variables:
		REAL, INTENT(IN)  :: Mmin    !M-coefficient for minimum reference point
		REAL, INTENT(IN)  :: Bmin    !B-coefficient for minimum reference point
		REAL, INTENT(IN)  :: Tmin    !Minimum reference temperature, C
		REAL, INTENT(IN)  :: Mmax    !M-coefficient for maximum reference point
		REAL, INTENT(IN)  :: Bmax    !B-coefficient for maximum reference point
		REAL, INTENT(IN)  :: Tmax    !Maximum reference temperature, C
		REAL, INTENT(IN)  :: Qcmp    !Compressor rated capacity, ton
		REAL, INTENT(IN)  :: Tevp    !Evaporating temperature, C
		REAL, INTENT(OUT) :: Qtxv    !TXV capacity, ton

		!Subroutine local variables:
		REAL Qmin !Minimum capacity, ton
		REAL Qmax !Maximum capacity, ton

		!Flow:

		Qmin=Mmin*Qcmp+Bmin !Linear correlation
		Qmax=Mmax*Qcmp+Bmax !Linear correlation

		Qtxv=(Tevp-Tmin)/(Tmax-Tmin)*(Qmax-Qmin)+Qmin

		RETURN

	END SUBROUTINE InterpolateQtxv

	!***********************************************************************************

	SUBROUTINE InterpolateCF_DP(aMin,bMin,Tmin,aMax,bMax,Tmax,DP,Tevp,CF)

		!Interpolate or extrapolate pressure drop correction factor

		IMPLICIT NONE

		!Subroutine passing variables:
		REAL, INTENT(IN)  :: aMin    !a-coefficient for minimum reference point
		REAL, INTENT(IN)  :: bMin    !b-coefficient for minimum reference point
		REAL, INTENT(IN)  :: Tmin    !Minimum reference temperature, C
		REAL, INTENT(IN)  :: aMax    !a-coefficient for maximum reference point
		REAL, INTENT(IN)  :: bMax    !b-coefficient for maximum reference point
		REAL, INTENT(IN)  :: Tmax    !Maximum reference temperature, C
		REAL, INTENT(IN)  :: DP      !Pressure drop across TXV, kPa
		REAL, INTENT(IN)  :: Tevp    !Evaporating temperature, C
		REAL, INTENT(OUT) :: CF      !Pressure drop correction factor

		!Subroutine local variables:
		REAL CFmin !Minimum correction factor
		REAL CFmax !Maximum correction factor

		!Flow:

		CFmin=aMin*DP**bMin !Power fit
		CFmax=aMax*DP**bMax !Power fit

		CF=(Tevp-Tmin)/(Tmax-Tmin)*(CFmax-CFmin)+CFmin

		RETURN

	END SUBROUTINE InterpolateCF_DP

	!***********************************************************************************

END MODULE TXVMOD
