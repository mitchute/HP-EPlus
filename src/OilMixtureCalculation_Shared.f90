! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module contains property calculation routines for refrigerant-oil mixtures
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This module does not represent a physical component of the heat pump system.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This module calculates and returns fluid properties for refrigerant-oil mixtures.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no associated input or output files.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! A few Copeland and Bristol oil properties are defined, as are the integers
! denoting manufacturers.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 13 methods:
!   PUBLIC LocalOilMassFraction -- Calculates the local oil mass fraction
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureSpecificHeat -- Calculates the mixture's specific heat
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureDensity -- Calculates the mixture's density
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureViscosity -- Calculates the mixture's viscosity
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureSurfaceTension -- Calculates the mixture's surface tension
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureThermalConductivity -- Calculates the mixture's thermal conductivity
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureReynoldsNumber -- Calculates the mixture's Reynolds Number
!       Called internally only
!   PUBLIC OilMixtureTsat -- Calculates the mixture's saturation temperature
!       Called by Condenser.f90
!       Called by Evaporator.f90
!   PUBLIC OilMixtureHTCevap -- Calculates the mixture's evaporation heat transfer coefficient
!       Called by CoilCalculation.f90
!   PUBLIC OilMixtureXtt -- Calculates the mixture's thermal conductivity
!       Called internally only
!   PUBLIC OilMixtureOutletEnthalpy -- Calculates the mixture's enthalpy
!       Never called
!   PRIVATE OilViscosity -- Calculates the mixture's viscosity
!       Called internally only
!   PRIVATE OilDensity -- Calculates the mixture's density
!       Called internally only

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
! There is probably some more documentation that could be done.
! Additionally, OilMixtureOutletEnthalpy doesn't appear to ever be called.

MODULE OilMixtureMod

  !Contains refrigerant-oil mixture property calculation

  USE DataGlobals_HPSimIntegrated, ONLY: RefName    !RS Comment: Needs to be used for implementation with Energy+ currently (7/23/12)

  IMPLICIT NONE

  PUBLIC LocalOilMassFraction
  PUBLIC OilMixtureSpecificHeat
  PUBLIC OilMixtureDensity
  PUBLIC OilMixtureViscosity
  PUBLIC OilMixtureSurfaceTension
  PUBLIC OilMixtureThermalConductivity
  PUBLIC OilMixtureReynoldsNumber
  PUBLIC OilMixtureTsat
  PUBLIC OilMixtureHTCevap
  PUBLIC OilMixtureXtt
  PUBLIC OilMixtureOutletEnthalpy

  PRIVATE OilViscosity
  PRIVATE OilDensity

  PRIVATE

  !Sample oil properties from ASHRAE handbook - refrigeration or
  !Jensen, M.K.and Jackman, D.L. (1984). Prediction of nucleate
  !pool boiling heat transfer coefficient of refrigerant-oil
  !mixtures, transactions of the ASME, Journal of heat transfer,
  !106, pp. 184-190.

  REAL,PARAMETER :: CopelandPOEmoleMass=570     !POE oil molecular mass
  REAL,PARAMETER :: CopelandPOEsurfaceTension=0.0251 !POE oil surface tension, N/m
  REAL,PARAMETER :: CopelandPOEconductivity=0.1237   !POE oil thermal conductivity, W/m-K

  REAL,PARAMETER :: CopelandMOmoleMass=378     !Mineral oil oil molecular mass
  REAL,PARAMETER :: CopelandMOsurfaceTension=0.0251 !Mineral oil surface tension, N/m
  REAL,PARAMETER :: CopelandMOconductivity=0.1249   !Mineral oil thermal conductivity, W/m-K

  REAL,PARAMETER :: BristolPOEmoleMass=570     !POE oil molecular mass
  REAL,PARAMETER :: BristolPOEsurfaceTension=0.0251 !POE oil surface tension, N/m
  REAL,PARAMETER :: BristolPOEconductivity=0.1237   !POE oil thermal conductivity, W/m-K

  REAL,PARAMETER :: BristolMOmoleMass=360     !Mineral oil oil molecular mass
  REAL,PARAMETER :: BristolMOsurfaceTension=0.047 !Mineral oil surface tension, N/m
  REAL,PARAMETER :: BristolMOconductivity=0.1312   !Mineral oil thermal conductivity, W/m-K

  !Compressor Manufacturer
  INTEGER,PARAMETER :: COPELAND  = 1 !ISI - 10/05/06
  INTEGER,PARAMETER :: BRISTOL   = 2
  INTEGER,PARAMETER :: DANFOSS   = 3
  INTEGER,PARAMETER :: PANASONIC = 4

CONTAINS

  !******************************************************************************

  REAL FUNCTION LocalOilMassFraction(AsoluteOilMassFraction,MixtureQuality)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate local oil mass fraction
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    REAL, INTENT(IN) :: AsoluteOilMassFraction !Asolute oil mass fraction
    REAL, INTENT(INOUT) :: MixtureQuality         !Ref.-Oil mixture quality

    !LOCAL VARIABLES:
    REAL xMax !Maximum oil-mixture quality

    !FLOW:

    IF (MixtureQuality .LE. 0) THEN
      MixtureQuality=0
    END IF
    IF (MixtureQuality .GE. 1)THEN
      MixtureQuality=1
    END IF

    xMax=1-AsoluteOilMassFraction

    IF (MixtureQuality .GT. xMax) THEN
      LocalOilMassFraction=0
    ELSE
      LocalOilMassFraction=AsoluteOilMassFraction/(1-MixtureQuality)
    END IF

    IF (LocalOilMassFraction .GT. 1) THEN
      LocalOilMassFraction=1
    END IF

    RETURN

  END FUNCTION LocalOilMassFraction

  !******************************************************************************

  REAL FUNCTION OilMixtureTsat(Wlocal,Psat)

    USE FluidProperties_HPSim !RS Comment: Currently needs to be used for integration with Energy+ Code (6/28/12)
    USE DataGlobals_HPSimIntegrated, ONLY: RefrigIndex   !RS: Debugging: Removal of plethora of RefrigIndex definitions in the code

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate saturation temperature for refrigerant-oil mixture, C
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    REAL,         INTENT(IN) :: Psat    !Saturation pressure for refrigerant, MPa
    REAL,         INTENT(IN) :: Wlocal  !Local oil mass fraction

    !LOCAL PARAMETERS:
    REAL,PARAMETER :: a1=182.52   !Empirical coefficient
    REAL,PARAMETER :: a2=-724.21  !Empirical coefficient
    REAL,PARAMETER :: a3=3868.0   !Empirical coefficient
    REAL,PARAMETER :: a4=-5268.9  !Empirical coefficient
    REAL,PARAMETER :: b1=-0.72212 !Empirical coefficient
    REAL,PARAMETER :: b2=2.3914   !Empirical coefficient
    REAL,PARAMETER :: b3=-13.779  !Empirical coefficient
    REAL,PARAMETER :: b4=17.066   !Empirical coefficient
    REAL,PARAMETER :: DP=0.1      !Pressure perturbation, MPa

    !LOCAL VARIABLES:
    REAL AA    !Intermediate variable
    REAL BB    !Intermediate variable
    REAL a0    !Estimate coefficient
    REAL b0    !Estimate coefficient
    REAL Psat1 !1st Perturbed saturation pressure, MPA
    REAL Psat2 !2nd Perturbed saturation pressure, MPA
    REAL Tsat1 !Saturation temperature at Psat1, K
    REAL Tsat2 !Saturation temperature at Psat2, K

    !EnergyPlus Refprop variables
    INTEGER(2) RefPropErr  !Error flag:1-error; 0-no error
    REAL Quality  !0-1
    REAL Pressure !Pa

    !FLOW:

    Psat1=Psat+DP   !MPa
    Psat2=Psat-DP   !MPa

    Pressure=Psat1*1e6 !from MPa to Pa
    Quality=1
    Tsat1=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !Saturation Temperature at Psat1
    Tsat1=Tsat1+273.15 !from C to K

    Pressure=Psat2*1e6 !from MPa to Pa
    Quality=1
    Tsat2=PQ(RefName, Pressure, Quality, 'temperature', RefrigIndex,RefPropErr) !Saturation Temperature at Psat2
    Tsat2=Tsat2+273.15 !from C to K

    b0=(LOG(Psat1)*(Tsat1/Tsat2)-LOG(Psat2))/((Tsat1/Tsat2)-1)
    a0=Tsat1*(LOG(Psat1)-b0)

    AA=a0+a1*wlocal+a2*wlocal**3+a3*wlocal**5+a4*wlocal**7
    BB=b0+b1*wlocal+b2*wlocal**3+b3*wlocal**5+b4*wlocal**7

    OilMixtureTsat=AA/(LOG(Psat)-BB)-273.15

    RETURN

  END FUNCTION OilMixtureTsat

  !******************************************************************************

  SUBROUTINE OilMixtureOutletEnthalpy(CompManufacturer,CoilType,Psat,mdot,Qmod, &
    Xin,hfg,Cpf,Cpg,Wabsolute,hin,hout,Xout)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture outlet enthalpy (J/kg) and quality
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    INTEGER,INTENT(IN) :: CoilType !1=Condenser; 2=Evaporator;
    !3=High side interconnecting pipes;
    !4=Low side interconnecting pipes
    !5=Microchannel condenser
    !6=Microchannel evaporator
    REAL, INTENT(IN) :: Psat !Saturation pressure for refrigerant, MPa
    REAL, INTENT(IN) :: mdot !Mass flow rate, kg/s
    REAL, INTENT(IN) :: Qmod !Segment capacity, W
    REAL, INTENT(INOUT) :: Xin  !Inlet oil mixture quality
    REAL, INTENT(IN) :: hfg  !Latent heat, J/kg
    REAL, INTENT(IN) :: Cpf  !Liquid specific heat, J/kg-K
    REAL, INTENT(IN) :: Cpg  !Vapor specific heat, J/kg-K
    REAL, INTENT(IN) :: Wabsolute !Absolute oil mass fraction
    REAL, INTENT(IN) :: Hin  !Inlet enthalpy change, J/kg

    !OUTPUTS
    REAL, INTENT(OUT) :: Hout !Outlet enthalpy change, J/kg
    REAL, INTENT(OUT) :: Xout !Outlet quality

    !LOCAL PARAMETERS:
    REAL, PARAMETER :: SMALL = 1E-3 !Iteration tolerance
    INTEGER, PARAMETER :: MaxIter=20 !Maximum number of iterations

    !LOCAL VARIABLES:
    REAL Xmin !Minimum quality
    REAL Xmax !Maximum quality
    REAL Xmix !Mixture quality
    REAL DXmix !Mixture quality difference
    REAL TsatIn !Inlet saturation temperature, K
    REAL TsatOut !Outlet saturation temperature, K
    REAL Tsat !Saturation temperature, K
    REAL DTsat !Saturation temperature difference, K
    REAL Wlocal !Local oil mass fraction
    REAL DH !Enthalpy change, J/kg
    INTEGER I !Iteration counter
    REAL CpMix !Mixture specific heat, J/kg-K
    REAL Qguess !Guessed heat capacity, W

    !FLOW:

    Wlocal=LocalOilMassFraction(Wabsolute,Xin)

    IF (Wlocal .LE. 0 .OR. Xin .GE. 1 .OR. Xin .LE. 0) THEN
      hout=-Qmod/mdot+hin     !Outlet Enthalpy, J/kg
      DXmix=ABS(hin-hout)/hfg
      SELECT CASE (CoilType)
      CASE(1,3,5)  !High pressure side
        Xout=Xin-DXmix
      CASE DEFAULT !Low pressure side
        Xout=Xin+DXmix
      END SELECT
      RETURN
    END IF

    TsatIn=OilMixtureTsat(Wlocal,Psat)  !Inlet Saturation Temperature, K

    SELECT CASE (CoilType)
    CASE(1,3,5)  !High pressure side
      Xmin=0
      Xmax=Xin
    CASE DEFAULT !Low pressure side
      Xmin=Xin
      Xmax=1-Wabsolute
    END SELECT

    Xout=(Xmin+Xmax)/2 !Initial guess

    DO I=1, MaxIter

      Wlocal=LocalOilMassFraction(Wabsolute,Xout) !Local Oil Mass Fraction
      TsatOut=OilMixtureTsat(Wlocal,Psat)         !Outlet Saturation Temperature, K

      DTsat=ABS(TsatIn-TsatOut)
      DXmix=ABS(Xin-Xout)

      Xmix=(Xin+Xout)/2
      Tsat=(TsatIn+TsatOut)/2

      CpMix=OilMixtureSpecificHeat(CompManufacturer,Wlocal,Cpf,Tsat)  !Specific heat of the mixture

      DH=hfg*DXmix+(1-Xmix)*CpMix*DTsat+Xmix*Cpg*DTsat    !Change in enthalpy

      Qguess=mdot*DH

      IF (ABS((Qguess-Qmod)/Qmod) .GT. SMALL) THEN
        IF (ABS(Qguess) .GT. ABS(Qmod)) THEN

          SELECT CASE (CoilType)
          CASE(1,3,5)  !High pressure side
            Xmin=Xout
          CASE DEFAULT !Low pressure side
            Xmax=Xout
          END SELECT

        ELSE

          SELECT CASE (CoilType)
          CASE(1,3,5)  !High pressure side
            Xmax=Xout
          CASE DEFAULT !Low pressure side
            Xmin=Xout
          END SELECT

        END IF
        Xout=(Xmin+Xmax)/2
      ELSE
        EXIT
      END IF

    END DO

    IF (I .GT. MaxIter) THEN
      hout=-Qmod/mdot+hin         !Outlet Enthalpy, J/kg
      DXmix=ABS(hin-hout)/hfg
      SELECT CASE (CoilType)
      CASE(1,3,5)  !High pressure side
        Xout=Xin-DXmix
      CASE DEFAULT !Low pressure side
        Xout=Xin+DXmix
      END SELECT
      RETURN
    END IF

    SELECT CASE (CoilType)
    CASE(1,3,5)  !High pressure side
      hout=hin-DH
    CASE DEFAULT !Low pressure side
      hout=hin+DH
    END SELECT

    RETURN

  END SUBROUTINE OilMixtureOutletEnthalpy

  !******************************************************************************

  REAL FUNCTION OilMixtureSpecificHeat(CompManufacturer,Wlocal,Cpf,Tref)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture specific heat, J/kg-K
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Wlocal !Local oil mass fraction
    REAL, INTENT(IN) :: Cpf    !Specific heat of liquid refrigerant, J/kg-K
    REAL, INTENT(IN) :: Tref   !Refrigerant temperature, C

    !LOCAL PARAMETERS:
    REAL, PARAMETER :: Treference = 15.6 !Reference tempeature, C

    !LOCAL VARIABLES:
    REAL CpOil    !Specific heat of oil, J/kg-K
    REAL rhoOil   !Oil denstiy, kg/m3
    REAL rhoWater !Water density, kg/m3

    REAL, EXTERNAL :: WRHO !HVACSIM+ Water property routine

    !FLOW:

    rhoOil=OilDensity(CompManufacturer,Tref)

    rhoWater=WRHO(Treference)

    CpOil=4186*(0.338+0.00045*(1.8*Tref+32))/SQRT(rhoOil/rhoWater)

    OilMixtureSpecificHeat=Wlocal*CpOil+(1-Wlocal)*Cpf

    RETURN

  END FUNCTION OilMixtureSpecificHeat

  !******************************************************************************

  REAL FUNCTION OilMixtureDensity(CompManufacturer,Wlocal,rhof,Tref)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture density, kg/m3
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Wlocal !Local oil mass fraction
    REAL, INTENT(IN) :: rhof   !Density of liquid refrigerant, kg/m3
    REAL, INTENT(IN) :: Tref   !Refrigerant temperature, C

    !LOCAL VARIABLES:
    REAL rhoOil   !Oil density, kg/m3

    !FLOW:

    rhoOil=OilDensity(CompManufacturer,Tref)

    OilMixtureDensity=1/(wlocal/rhoOil+(1-wlocal)/rhof)

    RETURN

  END FUNCTION OilMixtureDensity

  !******************************************************************************

  REAL FUNCTION OilMixtureViscosity(CompManufacturer,Wlocal,muf,MoleMassRef,Tref)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture viscosity, kg/s-m
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Wlocal      !Local oil mass fraction
    REAL, INTENT(IN) :: muf         !Viscosity of liquid refrigerant, kg/s-m
    REAL, INTENT(IN) :: MoleMassRef !Molecular mass of refrigerant, kg/mol
    REAL, INTENT(IN) :: Tref   !Refrigerant temperature, C

    !LOCAL PARAMETERS:
    REAL, PARAMETER :: KK = 0.58 !Empirical constant

    !LOCAL VARIABLES:
    REAL muOil   !Oil viscosity, N/m
    REAL MoleMassOil !Molecular mass of oil, kg/mol
    REAL MoleFracOil !Oil mole fraction
    REAL MoleFracRef !Refrigerant mole fraction
    REAL XiLiq       !Liquid refrigerant Yokozeki factor
    REAL XiOil       !Oil Yokozeki factor

    !FLOW:

    muOil=OilDensity(CompManufacturer,Tref)* &
    OilViscosity(CompManufacturer,Tref)*1e-6

    IF (muOil .LE. 0) THEN
      muOil=muf
    END IF

    SELECT CASE(CompManufacturer)
    CASE (COPELAND)
      IF (TRIM(RefName) .EQ. "R22") THEN
        MoleMassOil=CopelandMOmoleMass
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        MoleMassOil=CopelandPOEmoleMass
      END IF
    CASE (BRISTOL)
      IF (TRIM(RefName) .EQ. "R22") THEN
        MoleMassOil=BristolMOmoleMass
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        MoleMassOil=BristolPOEmoleMass
      END IF
    CASE DEFAULT
      IF (TRIM(RefName) .EQ. "R22") THEN
        MoleMassOil=BristolMOmoleMass
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        MoleMassOil=BristolPOEmoleMass
      END IF
    END SELECT

    MoleFracOil=Wlocal*(MoleMassRef/MoleMassOil)/ &
    (1-Wlocal+Wlocal*(MoleMassRef/MoleMassOil))

    MoleFracRef=1-MoleFracOil

    XiLiq=MoleMassRef**KK*MoleFracRef/ &
    (MoleMassRef**KK*MoleFracRef+MoleMassOil**KK*MoleFracOil)     !Liquid Refrigerant Yokozeki Factor

    XiOil=MoleMassOil**KK*MoleFracOil/ &
    (MoleMassRef**KK*MoleFracRef+MoleMassOil**KK*MoleFracOil)     !Oil Yokozeki Factor

    OilMixtureViscosity=EXP(XiLiq*LOG(muf)+XiOil*LOG(muOil))

    RETURN

  END FUNCTION OilMixtureViscosity

  !******************************************************************************

  REAL FUNCTION OilMixtureSurfaceTension(CompManufacturer,Wlocal,sigmaLiq)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture surface tension, N/m
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Wlocal   !Local oil mass fraction
    REAL, INTENT(IN) :: sigmaLiq !Surface tension of liquid refrigerant, N/m

    !LOCAL VARIABLES:
    REAL sigmaOil !Surface tension of oil, N/m

    !FLOW:

    SELECT CASE(CompManufacturer)
    CASE (COPELAND)
      IF (TRIM(RefName) .EQ. "R22") THEN
        sigmaOil=CopelandMOsurfaceTension
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        sigmaOil=CopelandPOEsurfaceTension
      END IF
    CASE (BRISTOL)
      IF (TRIM(RefName) .EQ. "R22") THEN
        sigmaOil=BristolMOsurfaceTension
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        sigmaOil=BristolPOEsurfaceTension
      END IF
    CASE DEFAULT
      IF (TRIM(RefName) .EQ. "R22") THEN
        sigmaOil=BristolMOsurfaceTension
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        sigmaOil=BristolPOEsurfaceTension
      END IF
    END SELECT

    OilMixtureSurfaceTension=sigmaLiq+(sigmaOil-sigmaLiq)*SQRT(Wlocal)

    RETURN

  END FUNCTION OilMixtureSurfaceTension

  !******************************************************************************

  REAL FUNCTION OilMixtureThermalConductivity(CompManufacturer,Wlocal,kLiq)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture thermal conductivity, W/m-K
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Wlocal !Local oil mass fraction
    REAL, INTENT(IN) :: kLiq   !Thermal conductivity of liquid refrigerant, W/m-K

    !LOCAL VARIABLES:
    REAL kOil !Thermal conductivity of oil, W/m-K

    !FLOW:

    SELECT CASE(CompManufacturer)
    CASE (COPELAND)
      IF (TRIM(RefName) .EQ. "R22") THEN
        kOil=CopelandMOconductivity
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        kOil=CopelandPOEconductivity
      END IF
    CASE (BRISTOL)
      IF (TRIM(RefName) .EQ. "R22") THEN
        kOil=BristolMOconductivity
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        kOil=BristolPOEconductivity
      END IF
    CASE DEFAULT
      IF (TRIM(RefName) .EQ. "R22") THEN
        kOil=BristolMOconductivity
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        kOil=BristolPOEconductivity
      END IF
    END SELECT

    OilMixtureThermalConductivity=kLiq*(1-Wlocal)+kOil*Wlocal- &
    0.72*(kOil-kLiq)*(1-Wlocal)*Wlocal

    RETURN

  END FUNCTION OilMixtureThermalConductivity

  !******************************************************************************

  SUBROUTINE OilMixtureReynoldsNumber(Gtot,Xmix,ID,muMix,muVap,ReLiq,ReVap)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture thermal conductivity, W/m-K
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    REAL, INTENT(IN) :: Gtot  !Total mass flux, kg/m2-s
    REAL, INTENT(IN) :: Xmix  !Mixture quality
    REAL, INTENT(IN) :: ID    !Tube inside diameter, m
    REAL, INTENT(IN) :: muMix !Mixture viscosity, kg/s-m
    REAL, INTENT(IN) :: muVap !Vapor viscosity, kg/s-m

    !OUTPUTS
    REAL, INTENT(OUT) :: ReLiq !Liquid Reynolds number
    REAL, INTENT(OUT) :: ReVap !Vapor Reynolds number

    !LOCAL VARIABLES:
    REAL Quality

    !FLOW:

    Quality=Xmix
    IF (Xmix .GT. 1) THEN
      Quality=1
    END IF
    IF (Xmix .LT. 0) THEN
      Quality=0
    END IF

    ReLiq=Gtot*(1-Quality)*ID/muMix

    ReVap=Gtot*Quality*ID/muVap

    RETURN

  END SUBROUTINE OilMixtureReynoldsNumber

  !******************************************************************************

  REAL FUNCTION OilMixtureHTCevap(Gtot,Xmix,ID,muMix,muVap,rhoMix,rhog, &
    kLiq,Cpf,OMF)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture evaporation heat transfer
    !coefficient, W/m2-K
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    REAL, INTENT(IN) :: Gtot   !Total mass flux, kg/m2-s
    REAL, INTENT(IN) :: Xmix   !Mixutre quality
    REAL, INTENT(IN) :: ID     !Tube inside diameter, m
    REAL, INTENT(IN) :: muMix  !Mixture viscosity, kg/s-m
    REAL, INTENT(IN) :: muVap  !Vapor viscosity, kg/s-m
    REAL, INTENT(IN) :: rhoMix !Mixture density, kg/m3
    REAL, INTENT(IN) :: rhog   !Vapor density, kg/m3
    REAL, INTENT(IN) :: kLiq   !Thermal conductivity of liquid refrigerant, W/m-K
    REAL, INTENT(IN) :: Cpf    !Specific heat of liquid refrigerant, J/kg-K
    REAL, INTENT(IN) :: OMF    !Oil mass fraction

    !LOCAL VARIABLES:
    REAL CC   !Empirical coefficient
    REAL nn   !Empirical coefficient
    REAL Xtt  !Lockhart-Martinelli parameter
    REAL hLiq !Mixture heat transfer coefficient, W/m2-K

    !FLOW:

    CC = 0.1385*(OMF*100) + 4.155
    nn= -0.0036*(OMF*100) + 0.6085

    Xtt=OilMixtureXtt(Gtot,Xmix,ID,muMix,muVap,rhoMix,rhog)

    hLiq=0.023*(kLiq/ID)*(Gtot*(1-Xmix)*ID/muMix)**0.8*(Cpf*muMix/kLiq)**0.4

    OilMixtureHTCevap=CC*(1/Xtt)**nn*hLiq   !Evaporation Heat Transfer Coefficient of the Mixture

    RETURN

  END FUNCTION OilMixtureHTCevap

  !******************************************************************************

  REAL FUNCTION OilMixtureXtt(Gtot,Xmix,ID,muMix,muVap,rhoMix,rhog)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate refrigerant-oil mixture thermal conductivity, W/m-K
    !
    !REFERENCE:
    !Schwentker, R.A. (2005). Advances to computer model used in the
    !simulation and optimization of heat exchangers. MS thesis.
    !The Univeristy of Maryland, College Park.
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    REAL, INTENT(IN) :: Gtot   !Total mass flux, kg/m2-s
    REAL, INTENT(IN) :: Xmix   !Mixutre quality
    REAL, INTENT(IN) :: ID     !Tube inside diameter, m
    REAL, INTENT(IN) :: muMix  !Mixture viscosity, kg/s-m
    REAL, INTENT(IN) :: muVap  !Vapor viscosity, kg/s-m
    REAL, INTENT(IN) :: rhoMix !Mixture density, kg/m3
    REAL, INTENT(IN) :: rhog   !Vapor density, kg/m3

    !LOCAL VARIABLES:
    REAL ReLiq !Liquid Reynolds number
    REAL ReVap !Vapor Reynolds number
    REAL fmix  !Refrigerant-oil mixture friction factor
    REAL fvap  !Vapor refrigerant friction factor
    REAL DPmixDZ !Refrigerant-oil mixture pressure drop gradient, Pa/m
    REAL DPvapDZ !Vapor mixture pressure drop gradient, Pa/m

    !FLOW:

    CALL OilMixtureReynoldsNumber(Gtot,Xmix,ID,muMix,muVap,ReLiq,ReVap)

    IF (ReLiq .LT. 1500) THEN !The two equations below intercept at Re=1500
      fmix=16/ReLiq
    ELSE
      fmix=0.046/ReLiq**0.2
    END IF

    IF (ReVap .LT. 2000) THEN
      fvap=16/ReVap
    ELSE
      fvap=0.046/ReVap**0.2
    END IF

    DPmixDZ=fmix*Gtot**2*(1-Xmix)**2/(2*ID*rhoMix)
    DPvapDZ=fvap*Gtot**2*Xmix**2/(2*ID*rhog)

    OilMixtureXtt=(DPmixDZ/DPvapDZ)**0.5    !Thermal conductivity of the oil mixture

    RETURN

  END FUNCTION OilMixtureXtt

  !******************************************************************************

  REAL FUNCTION OilViscosity(CompManufacturer,Temperature)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate oil kinematic viscosity, mm2/s
    !
    !REFERENCE:
    !Copleland data
    !Assume a linear relation
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Temperature !Temperature, C

    !FLOW:

    SELECT CASE(CompManufacturer)
    CASE (COPELAND)
      !Copeland
      IF (TRIM(RefName) .EQ. "R22") THEN
        OilViscosity=-0.5433*Temperature + 59.733
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilViscosity=-0.4283*Temperature + 48.683
      END IF

    CASE (BRISTOL)
      !Bristol
      IF (TRIM(RefName) .EQ. "R22") THEN
        !OilViscosity=-1.3218*Temperature + 77.895 !Check this
        OilViscosity=-0.5433*Temperature + 59.733 !use copeland property
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilViscosity= -0.4117*Temperature + 46.467
      END IF

    CASE DEFAULT
      !Bristol
      IF (TRIM(RefName) .EQ. "R22") THEN
        !OilViscosity=-1.3218*Temperature + 77.895 !Check this
        OilViscosity=-0.5433*Temperature + 59.733 !use copeland property
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilViscosity= -0.4117*Temperature + 46.467
      END IF

    END SELECT
    RETURN

  END FUNCTION OilViscosity

  !******************************************************************************

  REAL FUNCTION OilDensity(CompManufacturer,Temperature)

    IMPLICIT NONE

    !------------------------------------------------------------------------------
    !OBJECTIVE:
    !To calculate oil density, kg/m3
    !
    !REFERENCE:
    !Copleland data
    !Assume a linear relation
    !
    !Ipseng Iu
    !September 2006
    !------------------------------------------------------------------------------

    !INPUTS:
    INTEGER, INTENT(IN) :: CompManufacturer !Compressor manufacturer
    !1=Copland; 2=Bristol
    !3=Danfoss; 4=Panasonic
    REAL, INTENT(IN) :: Temperature !Temperature, C

    !FLOW:

    SELECT CASE(CompManufacturer)
    CASE (COPELAND)
      !Copeland
      IF (TRIM(RefName) .EQ. "R22") THEN
        OilDensity=-0.62*Temperature + 876.2
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilDensity=-Temperature + 1009
      END IF

    CASE (BRISTOL)
      !Bristol
      IF (TRIM(RefName) .EQ. "R22") THEN
        OilDensity=-1.0774*Temperature + 1155.7
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilDensity=-Temperature + 1009
      END IF

    CASE DEFAULT
      !Bristol
      IF (TRIM(RefName) .EQ. "R22") THEN
        OilDensity=-1.0774*Temperature + 1155.7
      ELSE !(TRIM(RefName) .EQ. "R410A") THEN
        OilDensity=-Temperature + 1009
      END IF

    END SELECT

    RETURN

  END FUNCTION OilDensity

  !******************************************************************************

END MODULE OilMixtureMod
