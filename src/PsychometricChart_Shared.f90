! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This module calculates the properties of air dependent on various intial properties (functions as a psychrometric chart)
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This does not represent a physical component of a Heat Pump system.
! 

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! THis calculates air properties as if using a psychrometric chart
!

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! no OPEN statements or output files.
!
! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
!  Nothing defined at the module level

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains X methods:
!   TDB_H - calculates air props given dry bulb and enthalpy
!       PsyChart
!   TDB_RH = calulates air props given dry bulb and relative humididty
!       PsyChart
!   TDB_TWB - calculates air props given dry bulb and wet bulb
!       PsyChart
!   TDB_W - calculates air props given dry bulb and humidity ratio
!       PsyChart
!   DEWPOINT - calculates the dew point
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   DRYBULB - calculates the dry bulb temp
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   ENTHALPY - calculates enthalpy
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   ENTHSAT - calculates saturation enthalpy
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   HUMRATIO - calculates the humidity ratio from water vapor pressure and atmospheric pressure
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   HUMTH - calculates humidity ratio from dry bulb and enthalpy
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   RELHUM - calculates the relative humidity
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   RHODRY - calculates density of the dry air
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   RHOMOIST - calculates the density of moist air
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   SATPRESS - calculate saturation pressure
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   SATTEMP - calculate saturation temperature
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   TAIRSAT - calculate dry bulb temperature at saturation
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   WETBULB - calculate wet bulb temperature
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   XITERATE - Iterately solves for the value of X which satisfies F(X)=0
!       TAIRSAT
!       Wetbulb
!   FTWB -  AIR WET BULB TEMPERATURES AS A FUNCTION OF THE AIR DRY BULB TEMPERATURE,THE HUMIDITY RATIO AND THE BAROMETRIC PRESSURE 
!       TDB_H
!       TDB_RH
!       TDB_TWB
!       TDB_W
!   FHAIR - Enthalpy of moist air
!       FTWB
!       FHSAT
!       FTSAT
!   FHSAT - Air Saturation Enthalpy as a function of saturation temperature
!       FTWB
!       FHAIR
!       FTSAT
!   FTSAT - Saturation temperature as a function of air saturation enthalpy
!       FTWB
!       FHARI
!       FHSAT
!   FPWS - saturation pressure over liquid water
!       FHSAT
!       FTWB
!   FWPW - Humidity ratio as a function of vapor pressure and barometric pressure
!       FHSAT
!       FTWB

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! no issues at this time

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-29 | JEH | Header completion 

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Nothing to do at this point

    !******************************************************************************

    SUBROUTINE TDB_H (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    SUBROUTINE:             TDB_H
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate psychrometric properties of
    !*                            moist air given dry bulb temperature
    !*                            and enthalpy.
    !***********************************************************************
    !*    INPUT VARIABLES
    !*    TDB           Dry bulb temperature                             (C)
    !*    H             Enthalpy of moist air                         (J/kg)
    !*    BaroPressure  Barometric pressure                            (kPa)
    !*
    !*    OUTPUT VARIABLES
    !*    W             Humidity ratio                                   (-)
    !*    RH            Relative humidity                                (-)
    !*    TWB           Wet bulb temperature                             (C)
    !*    TDP           Dewpoint temperature                             (C)
    !*    RhoD          Dry air density                              (kg/m3)
    !*    RhoM          Moist air density                            (kg/m3)
    !*    ErrStat       Error flag (0=ok, 1=error)                       (-)
    !*
    !*    PROPERTIES
    !*    Patm          Atmospheric pressure                            (Pa)
    !*    CpAir         Specific heat of air                        (J/kg C)
    !*    CpVap         Specific heat of water vapor                (J/kg C)
    !*    Hfg           Reference heat of vaporization of water       (J/kg)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     Perfect gas relationships
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       SATPRESS
    !                             RELHUM
    !                             WETBULB
    !                             DEWPOINT
    !                             RHODRY
    !                             RHOMOIST
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************

    IMPLICIT NONE
    INTEGER ErrStat
    REAL Patms,Hfg,CpVap,CpAir,psat,BaroPressure
    REAL TDB,W,H,RH,TWB,TDP,RhoD,RhoM,FTWB
    REAL RELHUM,WETBULB,DEWPOINT,RHODRY,RHOMOIST,SATPRESS
    LOGICAL IsCvg

    Patms=BaroPressure*1000 !ISI 05-23-05
    Hfg=2501000.
    CpVap=1805.
    CpAir=1004.

    ErrStat = 0

    !*** Calculate the saturation pressure at a given temperature.

    psat = SATPRESS (TDB)

    !*** Calculate the humidity ratio as a function of dry bulb
    !*** temperature and enthalpy.

    W = (H-CpAir*TDB)/(Hfg+CpVap*TDB)

    !*** Calculate the relative humidity as a function of partial water
    !*** vapor pressure and saturation pressure.

    RH = RELHUM (psat,W,BaroPressure)

    !*** Calculate wet bulb temperature as a function of enthalpy,
    !*** dry bulb temperature and humidity ratio.

    TWB = WETBULB (TDB,W,BaroPressure,IsCvg)
    IF (IsCvg .EQ. .FALSE.) THEN
        TWB = FTWB(TDB,W,BaroPressure) !Use HVACSIM+ psychrometric - ISI 01/20/04
    END IF

    !*** Calculate dewpoint temperature as a function of humidity ratio.

    TDP = DEWPOINT (W,BaroPressure)

    !*** Calculate dry and moist air density.

    RhoD = RHODRY (TDB,W,BaroPressure)
    RhoM = RHOMOIST (RhoD,W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END SUBROUTINE TDB_H

    !******************************************************************************

    SUBROUTINE TDB_RH (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    SUBROUTINE:             TDB_RH
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate psychrometric properties of
    !*                            moist air with input of dry bulb
    !*                            temperature and relative humidity.
    !***********************************************************************
    !*    INPUT VARIABLES
    !*    TDB           Dry bulb temperature                             (C)
    !*    RH            Relative humidity                                (-)
    !*    BaroPressure  Barometric pressure                            (kPa)
    !*
    !*    OUTPUT VARIABLES
    !*    W             Humidity ratio                                   (-)
    !*    H             Enthalpy                                      (J/kg)
    !*    TWB           Wet bulb temperature                             (C)
    !*    TDP           Dewpoint temperature                             (C)
    !*    RhoD          Dry air density                              (kg/m3)
    !*    RhoM          Moist air density                            (kg/m3)
    !*    ErrStat       Error flag (0=ok, 1=error)                       (-)
    !*
    !*    PROPERTIES
    !*    Patm          Atmospheric pressure                            (Pa)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     None
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       SATPRESS
    !                             HUMRATIO
    !                             WETBULB
    !                             ENTHALPY
    !                             DEWPOINT
    !                             RHODRY
    !                             RHOMOIST
    !
    !     REVISION HISTORY:       None
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************
    !     INTERNAL VARIABLES:
    !     pw            Partial water vapor pressure                    (Pa)
    !***********************************************************************

    IMPLICIT NONE
    INTEGER ErrStat
    REAL TDB,Patms,psat,pw,RH,W,H,TWB,TDP,RhoD,RhoM,BaroPressure
    REAL SATPRESS,HUMRATIO,ENTHALPY,WETBULB,FTWB
    REAL DEWPOINT,RHODRY,RHOMOIST
    LOGICAL IsCvg

    Patms=BaroPressure*1000 !ISI 05-23-05

    ErrStat = 0

    !*** Calculate the saturation pressure at a given temperature.

    psat = SATPRESS (TDB)

    !*** Calculate the partial water vapor pressure as a function of the
    !*** saturation pressure and relative humidity.

    pw = RH*psat

    !*** Calculate the humidity ratio as a function of the partial water
    !*** vapor pressure.

    W = HUMRATIO (Pw,BaroPressure)

    !*** Calculate the Halpy as a function of dry bulb temperature
    !*** and humidity ratio.

    H = ENTHALPY (TDB,W)

    !*** Calculate the wet bulb temperature as a function of dry bulb
    !*** temperature and humidity ratio.

    TWB = WETBULB (TDB,W,BaroPressure,IsCvg)
    IF (IsCvg .EQ. .FALSE.) THEN
        TWB = FTWB(TDB,W,BaroPressure) !Use HVACSIM+ psychrometric - ISI 01/20/04
    END IF

    !*** Calculate the dewpoint temperature as a function of the humidity
    !*** ratio.

    TDP = DEWPOINT (W,BaroPressure)

    !*** Calculate dry and moist air densities.

    RhoD = RHODRY (TDB,W,BaroPressure)
    RhoM = RHOMOIST (RhoD,W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END SUBROUTINE TDB_RH

    !******************************************************************************

    SUBROUTINE TDB_TWB (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    SUBROUTINE:             TDB_TWB
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate psychrometric properties of
    !*                            moist air given dry bulb temperature
    !*                            and wet bulb temperature.
    !***********************************************************************
    !*    INPUT VARIABLES
    !*    TDB           Dry bulb temperature                             (C)
    !*    TWB           Wet bulb temperature                             (C)
    !*    BaroPressure  Barometric pressure                            (kPa)
    !*
    !*    OUTPUT VARIABLES
    !*    W             Humidity ratio                                   (-)
    !*    RH            Relative humidity                                (-)
    !*    H             Enthalpy of air                               (J/kg)
    !*    TDP           Dewpoint temperature                             (C)
    !*    RhoD          Dry air density                              (kg/m3)
    !*    RhoM          Moist air density                            (kg/m3)
    !*    ErrStat       Error flag (0=ok, 1=error)                       (-)
    !*
    !*    PROPERTIES
    !*    Patm          Atmospheric pressure                            (Pa)
    !*    CpAir         Specific heat of air                        (J/kg C)
    !*    CpVap         Specific heat of water vapor                (J/kg C)
    !*    CpWat         Specific heat of liquid water               (J/kg C)
    !*    Hfg           Reference heat of vaporization of water       (J/kg)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     None
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       SATPRESS
    !                             HUMRATIO
    !                             RELHUM
    !                             ENTHALPY
    !                             DEWPOINT
    !                             RHODRY
    !                             RHOMOIST
    !
    !     REVISION HISTORY:       None
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************
    !     INTERNAL VARIABLES:
    !     psat          Saturated water vapor pressure                  (Pa)
    !     wstar         Saturated humidity ratio at wet bulb temp.       (-)
    !***********************************************************************

    IMPLICIT NONE
    INTEGER ErrStat
    REAL Patms,CpWat,CpVap,Hfg,CpAir,psat,TWB,wstar,BaroPressure
    REAL W,TDB,RH,H,TDP,RhoD,RhoM
    REAL SATPRESS,HUMRATIO,RELHUM,ENTHALPY,DEWPOINT
    REAL RHODRY,RHOMOIST

    !Patms=101325.
    Patms=BaroPressure*1000 !ISI 05-23-05

    CpWat=4186.
    CpVap=1805.
    Hfg=2501000.
    CpAir=1004.

    ErrStat = 0

    !*** Calculate the saturation pressure for a given wet bulb
    !*** temperature.

    psat = SATPRESS (TWB)

    !*** Calculate the humidity ratio corresponding to saturation
    !*** pressure for a given wet bulb temperature.

    wstar =  HUMRATIO (psat,BaroPressure)

    !*** Calculate the humidity ratio as a function of dry bulb
    !*** temperature, wet bulb temperature and the humidity ratio
    !*** at saturation pressure for the given wet bulb temperature.

    W = ((Hfg+(CpVap-CpWat)*TWB)*wstar-            &
    CpAir*(TDB-TWB))/                                     &
    (Hfg+CpVap*TDB-CpWat*TWB)

    !*** Calculate the saturation pressure for a given temperature.

    psat = SATPRESS (TDB)

    !*** Calculate the relative humidity as a function of partial
    !*** water vapor pressure and water vapor pressure at saturation.

    RH = RELHUM (psat,W,BaroPressure)

    !*** Calculate enthalpy as a function of dry bulb temperature
    !*** and humidity ratio.

    H = ENTHALPY (TDB,W)

    !*** Calculate the dewpoint temperature as a function of humidity
    !*** ratio.

    TDP = DEWPOINT (W,BaroPressure)

    !*** Calculate dry and moist air densities.

    RhoD = RHODRY (TDB,W,BaroPressure)
    RhoM = RHOMOIST (RhoD,W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END SUBROUTINE TDB_TWB

    !******************************************************************************

    SUBROUTINE TDB_W (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    SUBROUTINE:             TDB_W
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate psychrometric properties of
    !*                            moist air given dry bulb temperature
    !*                            and humidity ratio.
    !***********************************************************************
    !*    INPUT VARIABLES
    !*    TDB           Dry bulb temperature                             (C)
    !*    W             Humidity ratio                                   (-)
    !*    BaroPressure  Barometric pressure                            (kPa)
    !*
    !*    OUTPUT VARIABLES
    !*    RH            Relative humidity                                (-)
    !*    H             Enthalpy of moist air                         (J/kg)
    !*    TWB           Wet bulb temperature                             (C)
    !*    TDP           Dewpoint temperature                             (C)
    !*    RhoD          Dry air density                              (kg/m3)
    !*    RhoM          Moist air density                            (kg/m3)
    !*    ErrStat       Error flag (0=ok, 1=error)                       (-)
    !*
    !*    PROPERTIES
    !*    Patm          Atmospheric pressure                            (Pa)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     None
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       SATPRESS
    !                             RELHUM
    !                             WETBULB
    !                             ENTHALPY
    !                             DEWPOINT
    !                             RHODRY
    !                             RHOMOIST
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************

    IMPLICIT NONE
    INTEGER ErrStat
    REAL Patms,psat,TDB,RH,W,H,TWB,TDP,RhoD,RhoM,BaroPressure
    REAL SATPRESS,RELHUM,ENTHALPY,WETBULB,DEWPOINT,FTWB
    REAL RHODRY,RHOMOIST
    LOGICAL IsCvg

    Patms=BaroPressure*1000 !ISI 05-23-05

    ErrStat = 0

    !*** Calculate the saturation pressure at a given temperature.

    psat = SATPRESS (TDB)

    !*** Calculate the relative humidity as a function of the partial
    !*** pressure of water vapor and the saturation pressure of water vapor.

    RH = RELHUM (psat,W,BaroPressure)

    !*** Calculate enthalpy as a function of dry bulb temperature
    !*** and humidity ratio.

    H =  ENTHALPY (TDB,W)

    !*** Calculate the wet bulb temperature as a function of enthalpy or
    !*** dry bulb temperature and humidity ratio.

    TWB = WETBULB (TDB,W,BaroPressure,IsCvg)
    IF (IsCvg .EQ. .FALSE.) THEN
        TWB = FTWB(TDB,W,BaroPressure) !Use HVACSIM+ psychrometric - ISI 01/20/04
    END IF

    !*** Calculate the dewpoint temperature as a function of the humidity
    !*** ratio.

    TDP = DEWPOINT (W,BaroPressure)

    !*** Calculate dry and moist air densities.

    RhoD = RHODRY (TDB,W,BaroPressure)
    RhoM = RHOMOIST (RhoD,W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END SUBROUTINE TDB_W

    !******************************************************************************

    REAL FUNCTION DEWPOINT (W,BaroPressure)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    FUNCTION:               DEWPOINT
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate the dewpoint temperature given
    !*                            humidity ratio
    !***********************************************************************
    !*    INPUT VARIABLES
    !*    W             Humidity ratio                                   (-)
    !*    BaroPressure  Barometric pressure                            (kPa)
    !*
    !*    OUTPUT VARIABLES
    !*    DewPoint      Dew point temperature of air                     (C)
    !*
    !*    PROPERTIES
    !*    Patm          Atmospheric pressure                            (Pa)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     None
    !
    !     DEVELOPER:              Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       None
    !
    !     REVISION HISTORY:       None
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************
    !     INTERNAL VARIABLES:
    !     pw            Partial water vapor pressure                    (Pa)
    !     small         Small number
    !***********************************************************************

    IMPLICIT NONE
    REAL Patms,small,W,pw,BaroPressure
    REAL SATTEMP
    DATA small/1.E-9/

    Patms=BaroPressure*1000 !ISI 05-23-05

    !*** Test for "dry" air

    IF (W .LT. small) THEN
        DewPoint = -999
    ELSE

        !*** Calculate the partial water vapor pressure as a function of
        !*** humidity ratio.

        pw= Patms*W/(.62198+W)

        !*** Calculate dewpoint as saturation temperature at water vapor
        !*** partial pressure

        DewPoint = SATTEMP(pw)

    ENDIF

999 RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION DEWPOINT

    !******************************************************************************

    REAL FUNCTION DRYBULB (H,W)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    FUNCTION:               DRYBULB
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate the dry bulb temperature of
    !*                            moist air from enthalpy and humidity.
    !***********************************************************************
    !*    INPUT VARIABLES:
    !*    H             Enthalpy                                      (J/kg)
    !*    W             Humidity ratio                                   (-)
    !*
    !*    OUTPUT VARIABLES:
    !*    Drybulb       Dry bulb temperature                             (C)
    !*
    !*    PROPERTIES:
    !*    CpAir         Specific heat of air                        (J/kg C)
    !*    CpVap         Specific heat of water vapor                (J/kg C)
    !*    Hfg           Reference heat of vaporization of water       (J/kg)
    !***********************************************************************
    !     MAJOR RESTRICTIONS:     Uses perfect gas relationships
    !                             Fit for enthalpy of saturated water vapor
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       None
    !
    !     REVISION HISTORY:       None
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************

    !*** Calculate the dry bulb temperature as a function of enthalpy and
    !*** humidity ratio.
    !*** hDryAir = CpAir*TDB
    !*** hSatVap = Hfg + CpVap*TDB
    !*** Enthalpy = hDryAir + W*hSatVap

    IMPLICIT NONE
    REAL CpVap,CpAir,Hfg,H,W

    CpVap = 1805.
    CpAir=1004.
    Hfg=2501000.

    Drybulb = (H-Hfg*W)/(CpAir+CpVap*W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION DRYBULB

    !******************************************************************************

    REAL FUNCTION ENTHALPY (TDB,W)

    !***********************************************************************
    !*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !***********************************************************************
    !*    FUNCTION:               ENTHALPY
    !*
    !*    LANGUAGE:               FORTRAN 90
    !*
    !*    PURPOSE:                Calculate the enthalpy of moist air.
    !***********************************************************************
    !*    INPUT VARIABLES:
    !*    TDB           Dry bulb temperature                             (C)
    !*    W             Humidity ratio                                   (-)
    !*
    !*    OUTPUT VARIABLES:
    !*    Enthalpy      Enthalpy of moist air                         (J/kg)
    !*
    !*    PROPERTIES:
    !*    CpAir         Specific heat of air                        (J/kg C)
    !*    CpVap         Specific heat of water vapor                (J/kg C)
    !*    Hfg           Reference heat of vaporization of water       (J/kg)
    !***********************************************************************
    !     MAJOR RESTRICTIONS      Uses perfect gas relationships
    !                             Fit for enthalpy of saturated water vapor
    !
    !     DEVELOPER:              Shauna Gabel
    !                             Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       None
    !
    !     REVISION HISTORY:       None
    !
    !     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !***********************************************************************
    !
    !
    !!*** Calculate the enthalpy as a function of dry bulb temperature and
    !!*** humidity ratio.

    IMPLICIT NONE
    REAL CpAir,CpVap,Hfg,hDryAir,TDB,hSatVap,W

    CpAir=1004.
    CpVap=1805.
    Hfg=2501000.

    hDryAir = CpAir*TDB
    hSatVap = Hfg + CpVap*TDB
    Enthalpy = hDryAir + W*hSatVap

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION ENTHALPY

    !******************************************************************************

    REAL FUNCTION ENTHSAT (TDB,BaroPressure)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               ENTHSAT
    !C*
    !C*    LANGUAGE:               FORTRAN 77
    !C*
    !C*    PURPOSE:                Calculate the enthalpy at saturation
    !C*                            for given dry bulb temperature
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    TDB           Dry bulb temperature                             (C)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    EnthSat       Enthalpy at saturation                        (J/kg)
    !C*
    !C*    PROPERTIES
    !C*    Patm          Atmospheric pressure                            (Pa)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       SATPRESS
    !C                             HUMRATIO
    !C                             ENTHALPY
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     psat          Saturated water vapor pressure                  (Pa)
    !C     w             Humidity ratio                                   (-)
    !C***********************************************************************
    !
    !
    IMPLICIT NONE
    REAL Patms,psat,TDB,w,BaroPressure
    REAL SATPRESS, HUMRATIO, ENTHALPY

    Patms=BaroPressure*1000 !ISI 05-23-05

    !!*** Calculate the saturation pressure at the given temperature.

    psat = SATPRESS (TDB)

    !!*** Calculate the humidity ratio from the saturation pressure

    w = HUMRATIO (psat,BaroPressure)

    !!*** Calculate the enthalpy as a function of dry bulb temperature
    !!*** and humidity ratio.

    ENTHSAT = ENTHALPY (TDB,w)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION ENTHSAT

    !******************************************************************************

    REAL FUNCTION HUMRATIO (Pw,BaroPressure)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               HUMRATIO
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate the humidity ratio from water
    !C*                            vapor pressure and atmospheric pressure
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    Patm          Atmospheric pressure                            (Pa)
    !C*    Pw            Partial water vapor pressure                    (Pa)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    HumRatio      Humidity ratio                                   (-)
    !C***********************************************************************
    !C     MAJOR RESRICTIONS:      None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C************************************************************************

    !!*** Calculate the humidity ratio.
    IMPLICIT NONE
    REAL Patms,Pw,BaroPressure

    Patms=BaroPressure*1000 !ISI 05-23-05

    HumRatio = 0.62198*Pw/(Patms-Pw)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION HUMRATIO

    !******************************************************************************

    REAL FUNCTION HUMTH (TDB,H)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               HUMTH
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate the humidity ratio of moist air
    !C*                            from dry bulb temperature and enthalpy.
    !C***********************************************************************
    !C*    INPUT VARIABLES:
    !C*    H             Enthalpy                                      (J/kg)
    !C*    TDB           Dry bulb temperature                             (C)
    !C*
    !C*    OUTPUT VARIABLES:
    !C*    HumTH         Humidity ratio                                   (-)
    !C*
    !C*    PROPERTIES:
    !C*    CpAir         Specific heat of air                        (J/kg C)
    !C*    CpVap         Specific heat of water vapor                (J/kg C)
    !C*    Hfg           Reference heat of vaporization of water       (J/kg)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     Uses perfect gas relationships
    !C                             Fit for enthalpy of saturated water vapor
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !
    !
    !!*** Calculate humidity ratio from dry bulb temperature and enthalpy
    !!*** hDryAir = CpAir*TDB
    !!*** hSatVap = Hfg + CpVap*TDB
    !!*** Enthalpy = hDryAir + W*hSatVap
    IMPLICIT NONE
    REAL CpVap,Hfg,CpAir,H,TDB

    CpVap=1805.
    Hfg=2501000.
    CpAir=1004.

    HumTH = (H-CpAir*TDB)/(Hfg+CpVap*TDB)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION HUMTH

    !******************************************************************************

    REAL FUNCTION RELHUM (Psat,HumRatio,BaroPressure)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               RELHUM
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate the relative humidity from
    !C*                            saturation and atmospheric pressures
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    Patm          Atmospheric pressure                            (Pa)
    !C*    Psat          Saturation pressure                             (Pa)
    !C*    HumRatio      Humidity ratio                                   (-)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    RelHum        Relative humidity                                (-)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     pw            Partial water vapor pressure                    (Pa)
    !C***********************************************************************

    !!*** Calculate the partial water vapor pressure as a function of
    !!*** humidity ratio.
    IMPLICIT NONE
    REAL Patms,pw,HumRatio,Psat,BaroPressure

    Patms=BaroPressure*1000 !ISI 05-23-05

    pw = Patms*HumRatio/(.62198+HumRatio)

    !!*** Calculate the relative humidity as a function of partial water
    !!*** vapor pressure and water vapor pressure at saturation.

    RelHum = pw/Psat

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION RELHUM

    !******************************************************************************

    REAL FUNCTION RHODRY (TDB,W,BaroPressure)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               RHODRY
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate dry air density.
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    TDB           Dry bulb temperature                             (C)
    !C*    W             Humidity ratio                                   (-)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    RhoDry        Density of dry air                           (kg/m3)
    !C*
    !C*    PROPERTIES
    !C*    Patm          Atmospheric pressure                            (Pa)
    !C*    RAir          Gas constant for air                        (J/kg C)
    !C*    TAbsAdd       Additive constant to convert user T to absolute T
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     Perfect gas relationships
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     pAir          Partial pressure of dry air                     (Pa)
    !C***********************************************************************
    !
    !
    IMPLICIT NONE
    REAL Patms,RAir,pAir,W,TDB,BaroPressure

    Patms=BaroPressure*1000 !ISI 05-23-05

    RAir=287.055

    !!*** Calculate the dry air density from perfect gas laws.

    pAir = 0.62198*Patms/(0.62198+W)
    RhoDry = pAir/RAir/(TDB+273.15)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION RHODRY

    !******************************************************************************

    REAL FUNCTION RHOMOIST (RhoDry,W)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               RHOMOIST
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate moist air density from dry air
    !C*                            density and humidity ratio
    !C***********************************************************************
    !C*    INPUT VARIABLES:
    !C*    RhoDry                  Dry air density                    (kg/m3)
    !C*    W                       Humidity ratio                         (-)
    !C*
    !C*    OUTPUT VARIABLES:
    !C*    RhoMoist                Density of dry air                 (kg/m3)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C************************************************************************

    !!*** Calculate the moist air density
    IMPLICIT NONE
    REAL RhoDry,W

    RhoMoist = RhoDry*(1.+W)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION RHOMOIST

    !******************************************************************************

    REAL FUNCTION SATPRESS (T)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    SUBROUTINE:             SATPRESS
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate saturation pressure of water
    !C*                            vapor as a function of temperature
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    T             Temperature                                      (C)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    SatPress      Saturation pressure                             (Pa)
    !C*
    !C*    PROPERTIES
    !C*    TKelMult      Multiplying factor to convert user T to Kelvin
    !C*    TAbsAdd       Additive factor to convert user T to absolute T
    !C*                  tKel = 1. * (T + 273.15)
    !C*    PaMult        Multiplying factor to convert user P to Pascals
    !C*    PAbsAdd       Additive factor to convert user P to absolute P
    !C*                  Pa = 1. * (P + 0.)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     173.16 K <= Temp <= 473.15 K
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       None
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C
    !C                             Hyland, R.W., and A. Wexler.  1983.
    !C                             Formulations for the thermodynamic
    !C                             properties of the saturated phases of H2O
    !C                             from 173.15 K to 473.15 K.  ASHRAE
    !C                             Transactions, Vol. 89, No. 2A, pp. 500-519
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     tKel          Temperature in Kelvin                            (K)
    !C     pascals       Saturation pressure                             (Pa)
    !C***********************************************************************

    IMPLICIT NONE
    REAL C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13
    REAL tKel,pascals,T
    DATA C1/-5674.5359/,C2/6.3925247/,C3/-0.9677843E-2/
    DATA C4/0.62215701E-6/,C5/0.20747825E-8/,C6/-0.9484024E-12/
    DATA C7/4.1635019/,C8/-5800.2206/,C9/1.3914993/,C10/-0.048640239/
    DATA C11/0.41764768E-4/,C12/-0.14452093E-7/,C13/6.5459673/

    !!*** Convert temperature from user units to Kelvin.
    tKel = 1.*(T+273.15)

    !!*** If below freezing, calculate saturation pressure over ice.

    IF (tKel .LT. 273.15) THEN

        pascals = EXP(C1/tKel+C2+tKel*(C3+tKel*(C4+tKel*(C5+C6*tKel)))+C7*ALOG(REAL(tKel)))

        !!*** If above freezing, calculate saturation pressure over liquid water.

    ELSE IF (tKel .GE. 273.15) THEN

        pascals = EXP(C8/tKel+C9+tKel*(C10+tKel*(C11+tKel*C12))+C13*ALOG(REAL(tKel)))

    ENDIF

    !!*** Convert pressure from Pascals to user units

    SatPress = pascals/1. - 0.

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION SATPRESS

    !******************************************************************************

    REAL FUNCTION SATTEMP (P)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               SATTEMP
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate the saturation (boiling)
    !C*                            temperature of water given pressure
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    P             Pressure                                        (Pa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    SatTemp       Saturation temperature of water vapor            (C)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       SATPRESS
    !C                             XITERATE
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     tSat          Water temperature guess                          (C)
    !C     pSat          Pressure corresponding to temp. guess           (Pa)
    !C     error         Deviation of dependent variable in iteration
    !C     iter          Iteration counter
    !C     icvg          Iteration convergence flag
    !C     F1,F2         Previous values of dependent variable in XITERATE
    !C     X1,X2         Previous values of independent variable in XITERATE
    !C***********************************************************************
    !
    IMPLICIT NONE
    INTEGER itmax,iter,icvg
    REAL tSat,error,P,pSat,X1,F1,X2,F2
    REAL SATPRESS,XITERATE
    DATA itmax/50/
    
    CHARACTER(LEN=117),PARAMETER :: FMT_1001 = "(/1X,'*** ERROR IN FUNCTION SatTemp ***'/ 1X,'    Saturation temperature has not converged after ',I2,' iterations'/)"

    !!*** Use an iterative process to determine the saturation temperature
    !!*** at a given pressure using a correlation of saturated water vapor
    !!*** pressure as a function of temperature

    !!*** Initial guess of boiling temperature

    tSat = 100.

    !!*** Iterate to find the saturation temperature
    !!*** of water given the total pressure

    !!*** Set iteration loop parameters

    !VL: Previously: DO 100 iter = 1,itmax
    DO iter = 1,itmax

        !!*** Calculate saturation pressure for estimated boiling temperature

        pSat = SATPRESS(tSat)

        !!*** Compare with specified pressure and update estimate of temperature

        error = P - pSat
        tSat = XITERATE (tSat,error,X1,F1,X2,F2,iter,icvg)

        !!*** If converged leave loop iteration

        !VL: Previously: IF (icvg .EQ. 1) GO TO 110
        IF (icvg .EQ. 1) THEN
            SatTemp = tSat
            RETURN
        END IF

        !!*** Water temperature not converged, repeat calculations with new
        !!*** estimate of water temperature

        !VL: Previously:100 CONTINUE
    END DO

    !!*** Saturation temperature has not converged after maximum specified
    !!*** iterations. Print error message, set return error flag, and RETURN

    WRITE(77,FMT_1001) itmax
    
    !!VL: Previously: 
!!1001 FORMAT(/1X,'*** ERROR IN FUNCTION SatTemp ***'/ 1X,'    Saturation temperature has not converged after ',I2,' iterations'/)

    !VL: Previously: 110 SatTemp = tSat
    SatTemp = tSat

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION SATTEMP

    !******************************************************************************

    REAL FUNCTION TAIRSAT (HSat,BaroPressure)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               ENTHSAT
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate the dry bulb temperature given
    !C*                            enthalpy at saturation.
    !C***********************************************************************
    !C*    INPUT VARIABLES:
    !C*    HSat          Enthalpy at saturation                        (J/kg)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES:
    !C*    TAirSat       Dry bulb temperature                             (C)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       ENTHSAT
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     error         Deviation of dependent variable in iteration
    !C     iter          Iteration counter
    !C     icvg          Iteration convergence flag
    !C     F1,F2         Previous values of dependent variable in XITERATE
    !C     X1,X2         Previous values of independent variable in XITERATE
    !C***********************************************************************

    IMPLICIT NONE
    INTEGER itmax,iter,icvg
    REAL tSat,error,HSat,X1,F1,X2,F2,BaroPressure
    REAL ENTHSAT,XITERATE
    DATA itmax/20/,tSat/50./
    
    CHARACTER(LEN=105),PARAMETER :: FMT_1001 = "(/1X,'*** ERROR IN FUNCTION TAIRSAT ***'/1X,'    Temperature has not converged after ',I2,' iterations'/)"

    !!*** Estimate saturation temperature if reasonable value not available

    IF(tSat .LT. -200. .OR. tSat .GT. 1000.) THEN
        tSat = 50.
    END IF

    !!*** Calculate saturation temperature by iteration using function to
    !!*** calculate saturation enthalpy from temperature

    !VL: Previously:DO 100 iter=1,itmax
    DO iter=1,itmax

        error = HSat - ENTHSAT(tSat,BaroPressure)
        tSat = XITERATE(tSat,error,X1,F1,X2,F2,iter,icvg)

        !!*** If converged, leave iteration loop.

        !VL: Previously: IF (icvg .EQ. 1) GO TO 110
        IF (icvg .EQ. 1)THEN
            TAirSat = tSat
            RETURN
        END IF

        !!*** Temperature not converged, repeat calculation with new
        !!*** estimate of temperature.

        !VL: Previously:100 CONTINUE
    END DO

    !!*** Temperature has not converged after maximum specified
    !!*** iterations. Print error message and RETURN

    WRITE(77,FMT_1001) itmax
    !!VL: Previously: 
!!1001 FORMAT(/1X,'*** ERROR IN FUNCTION TAIRSAT ***'/1X,'    Temperature has not converged after ',I2,' iterations'/)

    !VL: Previously: 110 CONTINUE

    TAirSat = tSat

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION TAIRSAT

    !******************************************************************************

    REAL FUNCTION WETBULB (TDB,W,BaroPressure,IsCvg)

    !C***********************************************************************
    !C*    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !C***********************************************************************
    !C*    FUNCTION:               WETBULB
    !C*
    !C*    LANGUAGE:               FORTRAN 90
    !C*
    !C*    PURPOSE:                Calculate wet bulb temperature from dry
    !C*                            bulb temperature and humidity ratio
    !C***********************************************************************
    !C*    INPUT VARIABLES
    !C*    TDB           Dry bulb temperature                             (C)
    !C*    W             Humidity ratio of air                            (-)
    !C*    BaroPressure  Barometric pressure                            (kPa)
    !C*
    !C*    OUTPUT VARIABLES
    !C*    WetBulb       Wet bulb temperature                             (C)
    !C*
    !C*    PROPERTIES:
    !C*    Patm          Atmospheric pressure                            (Pa)
    !C*    Hfg           Latent heat of vaporization of water          (J/kg)
    !C*    CpAir         Specific heat of air                        (J/kg C)
    !C*    CpVap         Specific heat of water vapor                (J/kg C)
    !C*    CpWat         Specific heat of water                      (J/kg C)
    !C***********************************************************************
    !C     MAJOR RESTRICTIONS:     None
    !C
    !C     DEVELOPER:              Shauna Gabel
    !C                             Michael J. Brandemuehl, PhD, PE
    !C                             University of Colorado at Boulder
    !C
    !C     DATE:                   January 1, 1992
    !C
    !C     SUBROUTINES CALLED:     None
    !C     FUNCTIONS CALLED:       SATPRESS
    !C                             HUMRATIO
    !C                             SATTEMP
    !C                             XITERATE
    !C
    !C     REVISION HISTORY:       None
    !C
    !C     REFERENCE:              1989 ASHRAE Handbook - Fundamentals
    !C***********************************************************************
    !C     INTERNAL VARIABLES:
    !C     tBoil         Boiling temperature of water at given pressure   (C)
    !C     psatStar      Saturation pressure at wet bulb temperature      (C)
    !C     wStar         Humidity  ratio as a function of PsatStar        (-)
    !C     newW          Humidity ratio calculated with wet bulb guess    (-)
    !C     error         Deviation of dependent variable in iteration
    !C     iter          Iteration counter
    !C     icvg          Iteration convergence flag
    !C     F1,F2         Previous values of dependent variable in XITERATE
    !C     X1,X2         Previous values of independent variable in XITERATE
    !C***********************************************************************
    !
    IMPLICIT NONE
    INTEGER itmax,iter,icvg
    REAL Patms,newW,CpWat,CpVap,Hfg,CpAir,tBoil,TDB,BaroPressure
    REAL psatStar,wStar,error,X1,F1,X2,F2,W
    REAL SATTEMP,SATPRESS,HUMRATIO,XITERATE
    LOGICAL IsCvg !Check if solution converged. Edited by ISI 01/20/04
    DATA itmax/20/
    
    CHARACTER(LEN=114),PARAMETER :: FMT_1009 = "(/1X,'*** ERROR IN FUNCTION WetBulb ***'/1X,'    Wet bulb temperature has not converged after ',I2,' iterations'/)"

    Patms=BaroPressure*1000 !ISI 05-23-05

    CpWat=4186.
    CpVap=1805.
    Hfg=2501000.
    CpAir=1004.

    IsCvg=.TRUE.

    !*** Initial temperature guess

    tBoil = SATTEMP (Patms)
    WetBulb = MAX( MIN(WetBulb,TDB,(tBoil-0.1)), 0.)

    !*** Begin iteration loop

    !VL: Previously: DO 100 iter = 1,itmax
    DO  iter = 1,itmax

        IF (WetBulb .GE. (tBoil-0.09) ) THEN
            WETBULB = tBoil-0.1
        END IF

        !*** Determine the saturation pressure for wet bulb temperature

        psatStar = SATPRESS (WetBulb)

        !*** Determine humidity ratio for given saturation pressure

        wStar = HUMRATIO (psatStar,BaroPressure)

        !*** Calculate new humidity ratio and determine difference from known
        !*** humidity ratio

        newW = ((Hfg-(CpWat-CpVap)*WetBulb)*wStar-     &
        CpAir*(TDB-WetBulb))/(Hfg+CpVap*TDB    &
        -CpWat*WetBulb)

        !*** Check error, if not satisfied, calculate new guess and iterate

        error = W-newW
        WetBulb = XITERATE(WetBulb,error,X1,F1,X2,F2,iter,icvg)

        !*** If converged, leave iteration loop.

        !VL: Previously: IF (icvg .EQ. 1) GO TO 900
        IF (icvg .EQ. 1) THEN
            IF (WetBulb .GT. TDB) WetBulb = TDB
            RETURN
        END IF

        !*** Wet bulb temperature not converged, repeat calculation with new
        !*** estimate of wet bulb temperature.

        !VL: Previously: 100 CONTINUE
    END DO

    !*** Wet bulb temperature has not converged after maximum specified
    !*** iterations. Print error message, set return error flag, and RETURN

    WRITE(77,FMT_1009) itmax
    IsCvg=.FALSE.
    
    !!VL: Previously:
!!1009 FORMAT(/1X,'*** ERROR IN FUNCTION WetBulb ***'/1X,'    Wet bulb temperature has not converged after ',I2,' iterations'/)

    !VL: Previously: 900 IF (WetBulb .GT. TDB) WetBulb = TDB
    IF (WetBulb .GT. TDB) WetBulb = TDB

    !VL: Previously: 999 RETURN
    RETURN

    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION WETBULB

    !******************************************************************************

    REAL FUNCTION XITERATE (X0,F0,X1,F1,X2,F2,ICount,ICvg)

    !**********************************************************************
    !    Copyright ASHRAE.  Toolkit for HVAC System Energy Calculations
    !**********************************************************************
    !
    !    SUBROUTINE:             XITERATE
    !
    !    LANGUAGE:               FORTRAN 90
    !
    !    PURPOSE:                Iterately solves for the value of X which
    !                            satisfies F(X)=0. Given Xi,F(Xi) pairs,
    !                            the subroutine tests for convergence and
    !                            provides a new guess for the value of the
    !                            independent variable X.
    !**********************************************************************
    !    INPUT VARIABLES
    !    F0            Current value of the function F(X)
    !    X0            Current value of X
    !    F1,F2         Two previous values of F(Xi)
    !    X1,X2         Two previous values of X
    !
    !        NOTE:     F1,X1,F2,X2 MUST BE STORED AND SAVED IN CALLING
    !                  ROUTINE.  THEY NEED NO INITIALIZATION
    !
    !    ICount        Number of iterations
    !
    !    OUTPUT VARIABLES
    !    XIterate      New estimate of X for F(X)=0
    !    ICvg          Convergence flag  ICvg = 0:  Not converged
    !                                    ICvg = 1:  Converged
    !**********************************************************************
    !     DEVELOPER:              Michael J. Brandemuehl, PhD, PE
    !                             University of Colorado at Boulder
    !
    !     DATE:                   January 1, 1992
    !
    !     INCLUDE FILES:          None
    !     SUBROUTINES CALLED:     None
    !     FUNCTIONS CALLED:       None
    !
    !     REFERENCE:              None
    !**********************************************************************
    !     INTERNAL VARIABLES
    !     small         Small number used in place of zero
    !     mode          Number of points used in fit
    !                   mode = 1:  Use XPerburb to get new X
    !                   mode = 2:  Linear equation to get new X
    !                   mode > 2:  Quadratic equation to get new X
    !     coef(i)       Coefficients for quadratic fit
    !                   F(X) = coef(1) + coef(2)*X + coef(3)*X*X
    !     check         Term under radical in quadratic solution
    !     FiQ,XiQ       REAL values of Fi,Xi
    !     slope         Slope for linear fit
    !     tolRel        Relative error tolerance
    !     xPerturb      Perturbation applied to X to initialize iteration
    !**********************************************************************

    IMPLICIT NONE
    INTEGER ICvg,ICount,mode
    REAL coef(3),check,F0Q,F1Q,F2Q,X0Q,X1Q,X2Q,xOther
    REAL tolRel,xPerturb,small,X0,X1,F0,SLOPE,F1,X2,F2
    DATA tolRel/1.E-5/,xPerturb/0.1/,small/1.E-9/

    !*** Check for convergence by comparing change in X

    IF ((ABS(X0-X1) .LT. tolRel*MAX(ABS(X0),small) .AND.               &
    ICount .NE. 1) .OR. F0 .EQ. 0.) THEN
        XIterate = X0
        ICvg=1
        RETURN
    ENDIF

    !*** Not converged.
    !*** If after the second iteration there are enough previous points to
    !    fit a quadratic for the new X.  If the quadratic fit is not
    !    applicable, mode will be set to 1 or 2 and a new X will be
    !    determined by incrementing X from xPerturb or from a linear fit.

    ICvg=0
    mode=ICount

    DO WHILE (.TRUE.)
        !VL: Previously: 10      IF (mode .EQ. 1) THEN
        IF (mode .EQ. 1) THEN

            !*** New guess is specified by xPerturb

            IF (ABS(X0) .GT. small) THEN
                XIterate = X0*(1.+xPerturb)
            ELSE
                XIterate = xPerturb
            ENDIF

        ELSEIF (mode .EQ. 2) THEN

            !*** New guess calculated from LINEAR FIT of most recent two points

            SLOPE=(F1-F0)/(X1-X0)
            IF(slope.EQ.0) THEN
                mode=1
                !VL: Previously: GO TO 10
                CYCLE
            ENDIF
            XIterate=X0-F0/SLOPE
        ELSE

            !*** New guess calculated from QUADRATIC FIT

            !*** If two Xi are equal, set mode for linear fit and return to top

            IF (X0 .EQ. X1) THEN
                X1=X2
                F1=F2
                mode=2
                !VL: Previously: GO TO 10
                CYCLE
            ELSEIF (X0 .EQ. X2) THEN
                mode=2
                !VL: Previously: GO TO 10
                CYCLE
            ENDIF

            !*** Determine quadratic coefficients from the three data points
            !*** using REAL.

            F2Q=F2
            F1Q=F1
            F0Q=F0
            X2Q=X2
            X1Q=X1
            X0Q=X0
            coef(3)=((F2Q-F0Q)/(X2Q-X0Q)-(F1Q-F0Q)/(X1Q-X0Q))/(X2Q-X1Q)
            coef(2)=(F1Q-F0Q)/(X1Q-X0Q)-(X1Q+X0Q)*coef(3)
            coef(1)=F0-(coef(2)+coef(3)*X0Q)*X0Q

            !*** If points are colinear, set mode for linear fit and return to top

            IF (ABS(coef(3)) .LT. 1.D-10) THEN
                mode=2
                !VL: Previously: GO TO 10
                CYCLE
            ENDIF

            !*** Check for precision.  If the coefficients do not accurately
            !*** predict the given data points due to round-off errors, set
            !*** mode for a linear fit and return to top.

            IF (ABS((coef(1)+(coef(2)+coef(3)*X1Q)*X1Q-F1Q)/F1Q) .GT. 1.D-4) THEN
                mode=2
                !VL: Previously: GO TO 10
                CYCLE
            ENDIF

            !*** Check for imaginary roots.  If no real roots, set mode to
            !*** estimate new X by simply incrementing by xPerturb

            check=coef(2)**2-4*coef(1)*coef(3)
            IF (check .LT. 0) THEN

                !*** Imaginary roots -- go back to linear fit

                mode=2
                !VL: Previously: GO TO 10
                CYCLE

            ELSEIF (check .GT. 0) THEN

                !*** Real unequal roots -- determine root nearest to most recent guess

                XIterate=(-coef(2)+SQRT(check))/coef(3)/2
                xOther=-XIterate-coef(2)/coef(3)
                IF (ABS(XIterate-X0) .GT. ABS(xOther-X0)) THEN
                    XIterate=xOther
                END IF
            ELSE

                !*** Real Equal Roots -- one solution

                XIterate=-coef(2)/coef(3)/2
            ENDIF

        ENDIF

        EXIT

    END DO

    !*** Set previous variable values for the next iteration

    IF (mode .LT. 3) THEN

        !*** No valid previous points to eliminate.

        X2=X1
        F2=F1
        X1=X0
        F1=F0

    ELSE

        !*** Eliminate one previous point based on sign and magnitude of F(X)
        !*** Keep the current point and eliminate one of the previous ones.

        IF (F1*F0 .GT. 0 .AND. F2*F0 .GT. 0) THEN

            !*** All previous points of same sign.  Eliminate one with biggest F(X)

            IF (ABS(F2) .GT. ABS(F1)) THEN
                X2=X1
                F2=F1
            ENDIF
        ELSE

            !*** Points of different sign.
            !*** Eliminate the previous one with the same sign as current F(X).

            IF (F2*F0 .GT. 0)  THEN
                X2=X1
                F2=F1
            ENDIF
        ENDIF
        X1=X0
        F1=F0
    ENDIF
    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************

    END FUNCTION Xiterate

    !   Notes:
    !   Subroutines in this module are adapted from ASHRAE Toolkit for Secondary 
    !   HVAC System Energy Calculations which was originally written in FORTRAN 77,
    !   and put into a Loads Toolkit FORTRAN 90 module by Leon Dong, 1998.

    !******************************************************************************
    ! Copyright as described in the ASHRAE Toolkit License Agreement.
    !******************************************************************************

    !-------------------------------------------------------------------------------

    !From HVACSIM+ Edited by ISI - 01/20/04
    ! **********************************************************************
    ! * 13. AIR WET BULB TEMPERATURES AS A FUNCTION OF THE AIR DRY BULB    *
    ! * TEMPERATURE,THE HUMIDITY RATIO AND THE BAROMETRIC PRESSURE         *
    ! *                                                                    *
    ! * DETERMINATION OF TWB BY THE SECANT METHOD                          *
    ! **********************************************************************
    !
    REAL FUNCTION FTWB(TDB,W,BaroPressure)
    !
    DATA CPW/4.187/

    REAL TDB,W,BaroPressure
    REAL PATM,TWB1,PWS1,PWS2,WS1,WS2,CPW,H1,H,TWB2
    !
    ! *********************************************************************
    ! * INITIAL VALUES - TWB2 = SATURATION TEMPERATURE CORRESPONDING TO H *
    ! *                  TWB1 = TWB2 - 1                                  *
    ! ******************************************************************** *
    !
    PATM=BaroPressure*1000 !Calvin 05-23-05

    H=FHAIR(TDB,W)
    !
    TWB1=FTSAT(H,PATM)-5.
    PWS1=FPWS(TWB1)
    WS1=FWPW(PWS1,PATM)
    H1=FHAIR(TWB1,WS1)-CPW*TWB1*(WS1-W)
    DELTH1=H1-H
    !
    TWB2=TWB1+5.
    !
    ! ******************************
    ! * STARTING OF THE ITERATIONS *
    ! ******************************
    !

    !VL: Previously: DO 10 I=1,50
    DO I=1,50
        PWS2=FPWS(TWB2)
        WS2=FWPW(PWS2,PATM)
        H2=FHAIR(TWB2,WS2)-CPW*TWB2*(WS2-W)
        DELTH2=H2-H
        !
        !VL: Previously: IF(ABS(DELTH2).LT.1.E-3) GOTO 20
        IF(ABS(DELTH2).LT.1.E-3) THEN
            FTWB=TWB2
            RETURN
        END IF

        !
        TWB=TWB1-DELTH1*(TWB2-TWB1)/(DELTH2-DELTH1)
        TWB1=TWB2
        H1=H2
        DELTH1=DELTH2
        TWB2=TWB
    END DO
    FTWB=TWB2
    RETURN
    END FUNCTION

    !-------------------------------------------------------------------------------

    ! ****************************
    ! * 10. ENTHAPY OF MOIST AIR *
    ! ****************************
    !
    REAL FUNCTION FHAIR(TDB,W)
    DATA CPA/1.006/,CPG/1.805/,HFG/2501./

    REAL TDB,W
    !
    FHAIR=CPA*TDB+W*(CPG*TDB+HFG)
    !
    RETURN
    END FUNCTION

    !-------------------------------------------------------------------------------

    ! *********************************************************************
    ! * 11. AIR SATURATION ENTHALPY AS FUNCTION OF SATURATION TEMPERATURE *
    ! *********************************************************************
    !
    REAL FUNCTION FHSAT(TSAT,PATM)

    REAL TSAT,PATM

    REAL PS,WS

    PS=FPWS(TSAT)
    WS=FWPW(PS,PATM)
    FHSAT=FHAIR(TSAT,WS)
    RETURN
    END FUNCTION

    !-------------------------------------------------------------------------------

    ! **********************************************************************
    ! * 12. SATURATION TEMPERATURE AS A FUNCTION OF THE AIR SATURATION     *
    ! * ENTHALPY                                                           *
    ! *                                                                    *
    ! * Use piecewise linear approximation if 0 < Ts < 47 and atmospheric  *
    ! * pressure within +/- 0.5% of standard atmospheric pressure          *
    ! **********************************************************************
    !
    REAL FUNCTION FTSAT(HS,PATM)
    !
    REAL HS,PATM

    REAL TS1,TS2

    DATA C0/-6.0055/,C1/0.68510/,C2/-0.0056978/,C3/3.5344E-5/, &
    C4/-1.2891E-7/,C5/2.0165E-10/
    !
    DP=ABS(PATM-101325.0)
    IF (HS.GT.9.473 .AND. HS.LT.236.759 .AND. DP.LE.500.0) THEN
        !
        !  Temperature and pressure range OK for piecewise linear approximation 
        !  produced from ASHRAE Chap 6, Table 1 by GEK
        !
        IF (HS .GE. 9.4730 .AND. HS .LT. 18.639) THEN
            !-------- ( 0 < Ts < 5 )
            FTSAT = 0.545494*HS - 5.16747
        ELSEIF (HS .GE. 18.639 .AND. HS .LT. 31.724) THEN
            !-------- ( 5 < Ts < 11 )
            FTSAT = 0.45854*HS - 3.546733
        ELSEIF (HS .GE. 31.724 .AND. HS .LT. 48.) THEN
            !-------- ( 11 < Ts < 17 )
            FTSAT = 0.36844*HS - 0.694765
        ELSEIF (HS .GE. 48. .AND. HS .LT. 68.44) THEN
            !-------- ( 17 < Ts < 23 )
            FTSAT = 0.293542*HS + 2.90998
        ELSEIF (HS .GE. 68.44 .AND. HS .LT. 94.878) THEN
            !-------- ( 23 < Ts < 29 )
            FTSAT = 0.226946*HS + 7.467811
        ELSEIF (HS .GE. 94.878 .AND. HS .LT. 129.455) THEN
            !-------- ( 29 < Ts < 35 )
            FTSAT = 0.173526*HS + 12.536224
        ELSEIF (HS .GE. 129.455 .AND. HS .LT. 175.265) THEN
            !-------- ( 35 < Ts < 41 )
            FTSAT = 0.130976*HS + 18.04453
        ELSEIF (HS .GE. 175.265 .AND. HS .LT. 236.759) THEN
            !-------- ( 41 < Ts < 47 )
            FTSAT = 0.0975705*HS + 23.89931
        ENDIF
    ELSEIF (HS .LT. 9.473) THEN
        !
        !  Outside range of correlations
        !
        FTSAT = 0.
    ELSE
        !
        !  Temperature > 47 or pressure not within 0.5% of atmospheric, use
        !  secant method to find inverse of FHSAT
        !
        ! ******************
        ! * INITIAL VALUES *
        ! *****************************************************************
        ! * FIRST VALUE : FIT ON THE ASHRAE VALUES AT BAROMETRIC PRESSURE*
        ! *****************************************************************
        !
        TS1=C0+HS*(C1+HS*(C2+HS*(C3+HS*(C4+HS*C5))))
        HS1=FHSAT(TS1,PATM)
        DELTH1=HS1-HS
        IF (ABS(DELTH1).LT.1.E-3) THEN
            FTSAT=TS1
        ELSE
            !   Iterate (why is guess temperature reduced by 5 C ?)
            TS2=TS1-5.
            !
            ! ******************************
            ! * STARTING OF THE ITERATIONS *
            ! ******************************
            !
            DO I=1,50
                HS2=FHSAT(TS2,PATM)
                DELTH2=HS2-HS
                !
                IF (ABS(DELTH2).LT.1.E-3) THEN
                    FTSAT=TS2
                    RETURN
                END IF
                !
                TS=TS1-DELTH1*(TS2-TS1)/(DELTH2-DELTH1)
                TS1=TS2
                HS1=HS2
                DELTH1=DELTH2
                TS2=TS
            END DO
            
            FTSAT=TS2

        ENDIF
    ENDIF
    RETURN
    END FUNCTION

    !-------------------------------------------------------------------------------

    ! ****************************************************************
    ! * 1. SATURATION PRESSURE OVER LIQUID WATER FOR THE TEMPERATURE *
    ! * RANGE OF 0 DEG C TO 200 DEG C                                *
    ! * REFERENCE: ASHRAE FUNDAMENTALS HANDBOOK - CHAPTER 6          *
    ! ****************************************************************
    !
    REAL FUNCTION FPWS(TDB)
    !
    REAL TDB

    DIMENSION A(10),B(10)
    !
    DATA A /52.26,71.1,95.5,126.68,166.08, 215.38,276.36,351.16,441.94,551.36/
    DATA B /611.2,517.,273.,-194.7,-982.7, -2215.2,-4044.6,-6662.6,-10293.8,-15217.7/
    !
    DATA C1/-5800.2206/,C2/1.3914993/,C3/-0.048640239/, C4/0.41764768E-4/,C5/-0.14452093E-7/,C6/6.5459673/
    !
    IT=INT(TDB/5.0) + 1
    IF (IT.GT.0 .AND. IT.LT.10) THEN
        !  Temperature in range 0 - 50 C
        !  Use linear approximations produced by GEK from ASHRAE, Chap 6, Table 1
        FPWS=A(IT)*TDB + B(IT)
    ELSEIF (IT.GE.10 .AND. IT.LT.40) THEN
        !  50 < T < 200, outside range of piece-wise linear fit, use Chap 6, eqn 4
        TD=TDB+273.15
        FPWS=EXP(C1/TD+C2+C3*TD+C4*TD**2.+C5*TD**3.+C6*ALOG(TD))
    ELSEIF (IT.LE.0) THEN
        !  T < 0 , out of range, use T=0 value and write warning
        FPWS=611.2
        !WRITE(8,*) "FUNCTION FPWS: TEMPERATURE OUT OF RANGE, T=",TDB
    ELSE
        !  T > 200 , out of range, use T=200 value and write warning
        FPWS=1555074.
    ENDIF
    RETURN
    END FUNCTION

    !-------------------------------------------------------------------------------

    ! **********************************************************************
    ! * 4. HUMIDITY RATIO AS A FUNCTION OF THE VAPOR PRESSURE AND THE      *
    ! * BAROMETRIC PRESSURE                                                *
    ! **********************************************************************
    !
    REAL FUNCTION FWPW(PW,PATM)

    REAL PW,PATM
    !
    FWPW=0.62198*PW/(PATM-PW)
    RETURN
    END FUNCTION
    
    SUBROUTINE PsyChart(AirPropOpt,AirPropErr) !(AirProp,AirPropOpt,BaroPressure,AirPropErr)

!-----------------------------------------------------------------------------------
!
!  Description:	
!  Given two moist air properties, calculates others
!
!  Input:
!  BaroPressure = Barometric pressure, kPa
!  Inputs/Outputs:
!  AirProp(1)=Dry bulb temperature, C
!  AirProp(2)=Humidity ratio
!  AirProp(3)=Relative humidity
!  AirProp(4)=Enthalpy, kJ/kg
!  AirProp(5)=Wet bulb temperature, C
!  AirProp(6)=Dew point temperature, C
!  AirProp(7)=Dry air density, kg/m^3
!  AirProp(8)=Wet air density, kg/m^3
!  AirPropOpt=Air property calc. options
!         1:Given: Dry bulb temperature, C
!                  and Enthalpy, kJ/kg                
!         2:Given: Dry bulb temperature, C
!                  and Relative humidity, from 0 to 1
!         3:Given: Dry bulb temperature, C
!                  and Wet bulb temperature, C           
!         4:Given: Dry bulb temperature, C
!                  and Humidity ratio, from 0 to 1
!        *5:Given: Wet bulb temperataure, C *-not functioning
!                  and enthalpy, kJ/kg
!  AirPropErr=Error status: 1=error; 0=No error
!
!  Reference: 
!  ASHRAE Building Loads Toolkit (2000)
!
!  Author:
!  Ipseng Iu
!  Mechanical and Aerospace Engineering
!  Oklahoma State University, Stillwater	
!
!  Date: August 2002
!
!-----------------------------------------------------------------------------------
USE DataSimulation !, ONLY: APTDB, APHumRat, APRelHum, APEnth, APTWB, APTDP, APDryDens, APWetDens

IMPLICIT NONE

!Subroutine argument declarations
!REAL, INTENT(IN) :: BaroPressure
INTEGER(2), INTENT(IN) :: AirPropOpt
INTEGER(2), INTENT(OUT) :: AirPropErr
!REAL, INTENT(INOUT) :: AirProp(8)

!Subroutine local vairables
!INTEGER, PARAMETER :: MaxIter=20 !Maximum number of iterations
REAL, PARAMETER :: Small=1E-4 !Small number 
REAL, PARAMETER :: DT=1.0 !Delta temperature, C
REAL TDB  !Dry bulb temperature, C
REAL W    !Humidity ratio
REAL RH   !Relative humidity
REAL H    !Enthalpy, kJ/kg
REAL TWB  !Wet bulb temperature, C
REAL TDP  !Dew point temperature, C
REAL RhoD !Dry air density, kg/m^3
REAL RhoM !Moist air density, kg/m^3
INTEGER ErrStat       !Error status

!Flow**
  
  OPEN (77,FILE='AirProp.err')

  TDB  = AirProp%APTDB !RS: Debugging: Formerly AirProp(1)
  W    = AirProp%APHumRat  !RS: Debugging: Formerly AirProp(2)
  RH   = AirProp%APRelHum  !RS: Debugging: Formerly AirProp(3)
  H    = AirProp%APEnth*1000   !RS: Debugging: Formerly AirProp(4)
  TWB  = AirProp%APTWB !RS: Debugging: Formerly AirProp(5)
  TDP  = AirProp%APTDP !RS: Debugging: Formerly AirProp(6)
  RhoD = AirProp%APDryDens !RS: Debugging: Formerly AirProp(7)
  RhoM = AirProp%APWetDens !RS: Debugging: Formerly AirProp(8)
  
  SELECT CASE (AirPropOpt)
  CASE (1)
	CALL TDB_H (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)
  CASE (2)
    CALL TDB_RH (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)
  CASE (3)
    CALL TDB_TWB (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)
  CASE (4)
    CALL TDB_W (TDB,W,RH,H,TWB,TDP,RhoD,RhoM,BaroPressure,ErrStat)
  CASE DEFAULT
    ErrStat=1
    WRITE(77,*)'AirProp error: Invalid input data'
  END SELECT

  IF (TWB .GE. TDB) THEN
      RH = 1
  END IF

  AirProp%APTDB=TDB    !RS: Debugging: Formerly AirProp(1)
  AirProp%APHumRat=W      !RS: Debugging: Formerly AirProp(2)
  AirProp%APRelHum=RH     !RS: Debugging: Formerly AirProp(3)
  AirProp%APEnth=H/1000 !RS: Debugging: Formerly AirProp(4)
  AirProp%APTWB=TWB    !RS: Debugging: Formerly AirProp(5)
  AirProp%APTDP=TDP    !RS: Debugging: Formerly AirProp(6)
  AirProp%APDryDens=RhoD   !RS: Debugging: Formerly AirProp(7)
  AirProp%APWetDens=RhoM   !RS: Debugging: Formerly AirProp(8)

  AirPropErr=ErrStat
  IF (AirPropErr .EQ. 0) THEN
    IF (W .LT. 0 .OR. H .LT. 0 .OR. RhoD .LT. 0 .OR. RhoM .LT. 0) THEN
      WRITE(77,*)'AirProp error: Invalid input data'
	  AirPropErr=1
	END IF
  END IF
  
  CLOSE (77)

END SUBROUTINE

