	MODULE AirPropMod
	implicit none
!=======================================================================
!
!                  AIR PROPERTIES
!
!=======================================================================
!
!  Ref:  Irvine, T.F.Jr., and Liley, P.E., "Steam and Gas Tables
!        with Computer Equations," Academic Press, Inc., 1984.
!
!***********************************************************************

	PUBLIC CPCVA
	PUBLIC CPA
	PUBLIC HA
	PUBLIC PHIA
	PUBLIC TPHIA
	PUBLIC VISCA
	PUBLIC AKA
	PUBLIC TS
    PUBLIC Conductivity
    PUBLIC Viscosity
    PUBLIC CPAirFunction
    PUBLIC HUMTH

	CONTAINS

        SUBROUTINE CPCVA(TCA,CPA,CVA,GAMMA,SONIC)
        implicit none

! ----------------------------------------------------------------------
!
!   This subroutine takes Celsius air temperature, TCA, and computes:
!     CPA: Specific heat of air at constant pressure [KJ/(kg K)]
!     CVA: Specific heat of air at constant volume [KJ/(kg K)]
!     GAMMA: The ratio CPA/CVA [dimensionless]
!     SONIC: Speed of sound in air (M/S)
!
!***********************************************************************
        REAL, INTENT(IN) :: TCA
        REAL, INTENT(OUT):: CPA, CVA, GAMMA, SONIC
        REAL A0,A1,A2,A3,A4
        REAL TCONV
        REAL R
        REAL T
        DATA A0,A1,A2/1.03409,-0.284887E-3,0.7816818E-6/
        DATA A3,A4,TCONV,R/-0.4970786E-9,0.1077024E-12,273.15,0.287040/
        T=TCA+TCONV
        CPA=A0+T*(A1+T*(A2+T*(A3+T*A4)))
        CVA=CPA-R
        GAMMA=CPA/CVA
        SONIC=SQRT(GAMMA*R*T)
        RETURN
        END SUBROUTINE

!***********************************************************************

        REAL FUNCTION CPA(TC)
        implicit none

! ----------------------------------------------------------------------
!
!   This subroutine takes Celsius air temperature, TC, and computes:
!     CP: Specific heat of air at constant pressure [KJ/(kg K)]
!
!***********************************************************************
        REAL, INTENT(IN) :: TC
        REAL A0,A1,A2,A3,A4
        REAL TCONV
        REAL R
        REAL T
        DATA A0,A1,A2/1.03409,-0.284887E-3,0.7816818E-6/
        DATA A3,A4,TCONV,R/-0.4970786E-9,0.1077024E-12,273.15,0.287040/
        T=TC+TCONV
        CPA=A0+T*(A1+T*(A2+T*(A3+T*A4)))
        RETURN
        END FUNCTION

!***********************************************************************

        REAL FUNCTION HA(TC)
        implicit none

! ----------------------------------------------------------------------
!
!   Enthalpy of air [KJ/kg] as a function of temperature [C]
!
!***********************************************************************
        REAL, INTENT(IN) :: TC
        REAL A0,A1,A2,A3
        REAL TCONV
        REAL T
        DATA A0,A1,A2,A3/12.0740,0.924502,0.115984E-3,-0.563568E-8/
        DATA TCONV/273.15/
        T=TC+TCONV
        HA=A0+T*(A1+T*(A2+T*A3))

!   Note that internal energy, U, is equal to (HA - R*T)
!   where R = 0.287040

        RETURN
        END FUNCTION
!***********************************************************************

        REAL FUNCTION PHIA(TCA)
        implicit none

! ----------------------------------------------------------------------
!
!  Entropy function, PHI(T) = [ CPA(T)/T dT ] integrated from T0 to T
!      where T0 is a reference temperature at which entropy = 0.
!      PHIA has same units as entropy [KJ/(kg K)]
!  Note that S2 - S1 = PHIA(TC2) - PHIA(TC1) - R*LOG(P2/P1)
!      where R = 0.287040 [KJ/(kg K)]
!  Other useful relationships:
!      Isentropic pressure function, LOG(P/P0) = PHIA(TCA) / R
!      Isentropic volume function, LOG(V/V0) = LOG(R*T) - LOG(P/P0)
!
!***********************************************************************
        REAL, INTENT(IN) :: TCA
        REAL A0,A1,A2
        REAL TCONV
        REAL T
        DATA A0,A1,A2,TCONV/1.386989,0.184930E-3,0.95,273.15/
        T=TCA+TCONV
        PHIA=A0+A1*T+A2*ALOG(T)
        RETURN
        END FUNCTION
!***********************************************************************

        REAL FUNCTION TPHIA(PHIA)
        implicit none

! ----------------------------------------------------------------------
!
!  Temperature [C] of air as a function of the entropy function, PHIA
!
!***********************************************************************
        REAL, INTENT(IN) :: PHIA
        REAL A0,A1,A2,A3
        REAL TCONV
        REAL R
        REAL PR
        DATA A0,A1,A2,A3/-8800.92,1269.74,-61.9391,1.03530/
        DATA R,TCONV/0.287040,273.15/
        PR=PHIA/R
        TPHIA=A0+PR*(A1+PR*(A2+PR*A3))-TCONV
        RETURN
        END FUNCTION
!***********************************************************************

        REAL FUNCTION VISCA(TC)
        implicit none

! ----------------------------------------------------------------------
!
!   Dynamic viscosity [(N S)/(M*M)] of air, from celsius temperature
!
!***********************************************************************
        REAL, INTENT(IN) :: TC
        REAL A0,A1,A2,A3,A4
        REAL B0,B1,B2,B3,B4
        REAL TCONV
        REAL T
        DATA A0,A1,A2/-0.98601,9.080125E-2,-1.17635575E-4/
        DATA A3,A4/1.2349703E-7,-5.7971299E-11/
        DATA B0,B1,B2/4.8856745,5.43232E-2,-2.4261775E-5/
        DATA B3,B4,TCONV/7.9306E-9,-1.10398E-12,273.15/
        T=TC+TCONV
        VISCA=A0+T*(A1+T*(A2+T*(A3+T*A4)))
        IF(T .GE. 600.) VISCA=B0+T*(B1+T*(B2+T*(B3+T*B4)))
        VISCA=VISCA*1.E-6
        RETURN
        END FUNCTION
!***********************************************************************

        REAL FUNCTION AKA(TC)
        implicit none

! ----------------------------------------------------------------------
!
! Thermal conductivity of air [KW/(M K)], given Celsius temperature
!
!***********************************************************************
        REAL, INTENT(IN) :: TC
        REAL C0,C1,C2,C3,C4,C5
        REAL TCONV
        REAL T
        DATA C0,C1,C2/-2.276501E-3,1.2598485E-4,-1.4815235E-7/
        DATA C3,C4,C5/1.73550646E-10,-1.066657E-13,2.47663035E-17/
        DATA TCONV/273.15/
        T=TC+TCONV
        AKA=0.001*(C0+T*(C1+T*(C2+T*(C3+T*(C4+T*C5)))))
        RETURN
        END FUNCTION
!***********************************************************************	

        REAL FUNCTION TS(HS)    !RS: Comment: Accurate between about 0 & 35 C according to John Gall (2/20/14)
        implicit none

! ----------------------------------------------------------------------
!
!  Saturation temperature [C], given saturation enthalpy [kJ/kg]
!
!***********************************************************************
        REAL, INTENT(IN) :: HS
        REAL C0,C1,C2,C3
        DATA C0,C1/1.050415E-05,-3.801049E-03/
	    DATA C2,C3/6.297500E-01,-5.549214/

	    TS=C0*HS**3+C1*HS**2+C2*HS+C3

        RETURN
        END FUNCTION
!***********************************************************************	
	

    REAL FUNCTION Viscosity(T,omega)
    !A curve-fit function for viscosity at standard atmospheric pressure using temperature and humidity ratio
    !Replacing VISCA dry air function
    !Function written by John Gall, 2/19/14, developed using EES tables
    !Function implemented by Rachel Spitler, 2/19/14
    !Viscosity returns in kg/m-s or Pa-s
    REAL T !Drybulb Temperature; enters in C
    REAL omega !Humidity Ratio; dimensionless (lbs water/lbs air)
    
    T=T+273.15 !Temperature must be in Kelvin for calculation

    Viscosity=-2.63899258E-06+1.22041696E-07*T-3.39028951E-10*T**2+9.38085074E-13*T**3-1.74727339E-15*T**4+ &
    1.88832644E-18*T**5-8.93100404E-22*T**6+3.64937140E-06*omega+5.40589192E-05*omega**2-1.75221009E-03*omega**3- &
    9.56224115E-02*omega**4+1.01721262E+00*omega**5-2.90411023E-05*omega**6-8.38510499E-08*T*omega- &
    6.66693283E-07*T*omega**2+2.87464796E-05*T*omega**3+1.08080907E-03*T*omega**4-8.56573935E-03*T*omega**5+ &
    4.17691666E-10*T**2*omega+2.91669508E-09*T**2*omega**2-1.79714662E-07*T**2*omega**3-4.45348761E-06*T**2*omega**4+ &
    1.09387919E-05*T**2*omega**5-1.05990032E-12*T**3*omega-5.26376501E-12*T**3*omega**2+5.39307579E-10*T**3*omega**3+ &
    7.54358241E-09*T**3*omega**4+8.94390631E-08*T**3*omega**5+1.56449520E-15*T**4*omega+1.46752094E-15*T**4*omega**2- &
    7.74004805E-13*T**4*omega**3-3.02639597E-12*T**4*omega**4-3.18838955E-10*T**4*omega**5-9.77624818E-19*T**5*omega+ &
    3.62228022E-18*T**5*omega**2+4.21562455E-16*T**5*omega**3-2.88448086E-15*T**5*omega**4+3.02868674E-13*T**5*omega**5

    T=T-273.15 !Converting back to C
    
    RETURN
    END FUNCTION
    
    REAL FUNCTION Conductivity(T,omega)
    !A curve-fit function for conductivity at standard atmospheric pressure using temperature and humidity ratio
    !Replacing AKA dry air function
    !Function written by John Gall, 2/19/14, developed using EES tables
    !Function implemented by Rachel Spitler, 2/19/14
    !Conductivity returns in kW/m-k
    REAL T !Drybulb Temperature; enters in C
    REAL omega !Humidity Ratio; dimensionless (lbs water/lbs air)
    
    T=T+273.15 !Temperature must be in Kelvin for calculation
    
    Conductivity=9.80682445E-04+9.08985840E-05*T-3.30678040E-08*T**2+2.36199493E-11*T**3-4.99595330E-14*T**4+ &
    6.99708802E-17*T**5-4.07364329E-20*T**6+4.84035162E-03*omega-4.61926181E-02*omega**2-2.39809182E-01*omega**3- &
    3.50193353E+01*omega**4+5.77865105E+02*omega**5-6.95234708E-02*omega**6-3.21006378E-05*T*omega+6.11486828E-04*T*omega**2+ &
    3.13721248E-03*T*omega**3+3.49142245E-01*T*omega**4-4.71971578E+00*T*omega**5+3.58018349E-08*T**2*omega- &
    3.79079053E-06*T**2*omega**2-1.19688437E-05*T**2*omega**3-1.10105993E-03*T**2*omega**4+4.55860193E-03*T**2*omega**5+ &
    4.35541522E-10*T**3*omega+1.08299931E-08*T**3*omega**2+9.21252918E-09*T**3*omega**3+5.00126319E-07*T**3*omega**4+ &
    5.78375627E-05*T**3*omega**5-8.51658226E-13*T**4*omega-1.59630736E-11*T**4*omega**2+3.79447691E-11*T**4*omega**3+ &
    3.18448537E-09*T**4*omega**4-1.94386839E-07*T**4*omega**5+5.90998540E-16*T**5*omega+9.26815281E-15*T**5*omega**2- &
    5.97258541E-14*T**5*omega**3-4.11201230E-12*T**5*omega**4+1.81436172E-10*T**5*omega**5
    
    Conductivity=Conductivity*0.001

    T=T-273.15 !Converting back to C
    
    RETURN
    END FUNCTION
    
    REAL FUNCTION CPAirFunction(T,dw)
    !Modified from the E+ function PsyCpAirFnWTdb; info on that is below
    !Replaces CPA dry air function
    !Function implemented by Rachel Spitler, 2/19/14

          ! FUNCTION INFORMATION:
          !       AUTHOR         J. C. VanderZee
          !       DATE WRITTEN   Feb. 1994
          !       MODIFIED       na
          !       RE-ENGINEERED  na

          ! PURPOSE OF THIS FUNCTION:
          ! This function provides the heat capacity of air {kJ/kg-C} as function of humidity ratio.

          ! METHODOLOGY EMPLOYED:
          ! take numerical derivative of PsyHFnTdbW function

          ! REFERENCES:
          ! see PsyHFnTdbW ref. to ASHRAE Fundamentals
          ! USAGE:  cpa = PsyCpAirFnWTdb(w,T)

          ! USE STATEMENTS:
          ! na

          ! FUNCTION ARGUMENT DEFINITIONS:
      REAL dw    ! humidity ratio {kgWater/kgDryAir}
      REAL T    ! input temperature {Celsius}
      
          ! FUNCTION LOCAL VARIABLE DECLARATIONS:
      REAL h1  ! PsyHFnTdbW result of input parameters
      REAL tt  ! input temperature (T) + .1
      REAL h2  ! PsyHFnTdbW result of input humidity ratio and tt
      REAL w  ! humidity ratio

      w=MAX(dw,1.0d-5)
      h1 = 1.00484d3*T+w*(2.50094d6+1.85895d3*T)
      tt = T + 0.1d0
      h2 = 1.00484d3*tt+w*(2.50094d6+1.85895d3*tt)
      CPAirFunction = ((h2-h1)/0.1d0)/1000 !Converting to kJ

      RETURN
    END FUNCTION CPAirFunction
    
    REAL FUNCTION HUMTH (TDB,H) !RS: Replace: Moving a copy here from PsyChart

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

    H=H*1000 !RS: Debugging: Assuming that h passed in is in kJ/kg (2/20/14)
    CpVap=1805.
    Hfg=2501000.
    CpAir=1004.

    HumTH = (H-CpAir*TDB)/(Hfg+CpVap*TDB)
    
    H=H/1000 !RS: Debugging: Converting back to kJ/kg (2/20/14)

    RETURN
    !***************************
    !   EVOLUTIONARY HISTORY:
    !   This subroutine is adapted from ASHRAE Toolkit for Secondary HVAC System 
    !   Energy Calculations which was originally written in FORTRAN 77.
    !***************************
    END FUNCTION HUMTH
    
	END MODULE AirPropMod
