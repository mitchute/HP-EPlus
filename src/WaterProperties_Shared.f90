!=======================================================================
!
!                  STEAM AND LIQUID WATER PROPERTIES
!
!=======================================================================
!
!  References:
!
!    1)  Irvine, T.F.Jr., and Liley, P.E., "Steam and Gas Tables
!        with Computer Equations," Academic Press, Inc., 1984.
!    2)  Van Wylen, G.J., and Sonntag, R.E., Fundamentals of Classical
!        Thermodynamics (SI Version).  New York: John Wiley and Sons.
!    3)  Chapman, A.J.  Heat Transfer (3rd Edition).  New York:
!        Macmillan  Publishing Co. (1974).
!    4)  Chi, J.
!    5)  Karlekar, B.V., and Desmond, R.M., Engineering Heat Transfer.
!        St. Paul: West Publishing Co., 1977.
!    6)  CRC Handbook of Chemistry and Physics, 61st Ed. (1980-1981).
!
!  This file contains the following functions:
!
!  Superheated Steam:
!       VS(PKPA,TC)     eq. from [1]
!       HS(PKPA,TC)     eq. from [1]
!       SS(PKPA,TC)     eq. from [1]
!       CPS(TC)         eq. from [2]
!       TPSS(PKPA,S)    derived using functions SS and CPS
!       VISSPH(TC)      curve fit to data in [3]
!       STEAMK(T)       curve fit to data in [3]
!       CVS(V,T)        adapted from subroutine in [4]
!
!  Saturated Steam:
!       TSATS(PKPA)     eq. from [1]
!       PSATS(TC)       eq. from [1]
!       VSATS(PKPA,TC)  eq. from [1]
!       HSATS(TC)       eq. from [1]
!       SSATS(TC)       eq. from [1]
!       VISSV(P)        curve fit to data in [3]
!
!  Two Phase:
!       HFG(TC)         eq. from [1]
!
!  Saturated Liquid Water:
!       VSATW(TC)       eq. from [1]
!       HSATW(TC)       eq. from [1]
!       SSATW(TC)       eq. from [1]
!       WMU(T)          curve fit to data in [5] for T above 100 C.
!       WCP(T)          curve fit to data in [5] for T above 100 C.
!       WK(T)           curve fit to data in [5] for T above 100 C.
!
!  Liquid Water at 1 Atmosphere:
!       WMU(T)          eq. from [6]
!       WCP(T)          eq. from [6]
!       WK(T)           curve fit to data in [6]
!       WRHO(T)         eq. from [6]
!
!***********************************************************************

      FUNCTION TSATS(PKPA)

! ----------------------------------------------------------------------
!
!  Saturation temp. of steam (C) as a function of pressure (KPA)
!
!***********************************************************************

      DATA A1,B1,C1,TCONV/42.6776,-3892.70,-9.48654,-273.15/
      DATA A2,B2,C2,PCONV/-387.592,-12587.5,-15.2578,0.001/

      P=PKPA*PCONV
      IF(P .LT. 12.33) THEN
          TSATS=TCONV+A1+B1/(ALOG(P)+C1)
      ELSE
          TSATS=TCONV+A2+B2/(ALOG(P)+C2)
      ENDIF

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION PSATS(TC)

! ----------------------------------------------------------------------
!
!  Saturation pressure of steam (KPA) as a function of temperature (C)
!
!***********************************************************************

      DATA A0,A1,A2,A3/10.4592,-0.404897E-2,-0.417520E-4,0.368510E-6/
      DATA A4,A5,A6/-0.101520E-8,0.865310E-12,0.903668E-15/
      DATA A7,A8,A9/-0.199690E-17,0.779287E-21,0.191482E-24/
      DATA A10,A11,TCONV,PCONV/-3968.06,39.5735,273.15,1000./

      T=TC+TCONV
        PLOG=A0+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*(A6+T*(A7+T*(A8+T*A9 )))))))) + A10/(T-A11)
      PSATS=PCONV*EXP(PLOG)

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION VSATS(PKPA,TC)

! ----------------------------------------------------------------------
!
!  Sat. specific volume of steam (M3/KG) given sat. T (C) and P (KPA)
!
!***********************************************************************

      DATA A,B,C,D,E1/1.,1.6351057,52.584599,-44.694653,-8.9751114/
      DATA E2,E3,E4,E5/-0.43845530,-19.179576,36.765319,-19.462437/
      DATA TCR,PCR,VCR,TCNV,PCNV/647.3,22.089,3.155E-3,273.15,0.001/

      TR=(TCR-TC-TCNV)/TCR
      Y=A+B*TR**(1./3.)+C*TR**(5./6.)+D*TR**0.875
      Y=Y+TR*(E1+TR*(E2+TR*(E3+TR*(E4+TR*E5))))
      VSATS=Y*PCR*VCR/(PKPA*PCNV)

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION VSATW(TC)

! ----------------------------------------------------------------------
!
!  Sat. specific volume of water (M3/KG) given sat. T (C)
!
!***********************************************************************

      DATA A,B,C,D,E1/1.,-1.9153882,12.015186,-7.8464025,-3.8886414/
      DATA E2,E3,E4,E5/2.0582238,-2.0829991,0.82180004,0.47549742/
      DATA TCR,VCR,TCNV/647.3,3.155E-3,273.15/

      TR=(TCR-TC-TCNV)/TCR
      Y=A+B*TR**(1./3.)+C*TR**(5./6.)+D*TR**0.875
      Y=Y+TR*(E1+TR*(E2+TR*(E3+TR*(E4+TR*E5))))
      VSATW=Y*VCR

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION HSATW(TC)

! ----------------------------------------------------------------------
!
!  Sat. enthalpy of liquid water (KJ/KG) given sat. T (C)
!
!***********************************************************************

      DATA E11,E21,E31,HFCR/624.698837,-2343.85369,-9508.12101,2099.3/
      DATA E41,E51,E61,TCNV/71628.7928,-163535.221,166531.093,273.15/
      DATA E71,A2,E12/-64785.4585,0.8839230108,-2.67172935/
      DATA E22,E32,E42/6.22640035,-13.1789573,-1.91322436/
      DATA E52,E62,E72,A3/68.7937653,-124.819906,72.1435404,1.0/
      DATA B3,C3,D3/-0.441057805,-5.52255517,6.43994847/
      DATA E13,E23,TCR/-1.64578795,-1.30574143,647.3/

      TK=TC+TCNV
      TR=(TCR-TK)/TCR
      IF(TK.LT.300)THEN
          Y=TR*(E11+TR*(E21+TR*(E31+TR*(E41+TR*(E51+TR*(E61+TR*E71))))))
      ELSE IF(TK.LT.600)THEN
          Y=TR*(E12+TR*(E22+TR*(E32+TR*(E42+TR*(E52+TR*(E62+TR*E72))))))
          Y=Y+A2
      ELSE
        Y=A3+B3*TR**(1./3.)+C3*TR**(5./6.)+D3*TR**0.875+TR*(E13+TR*E23)
      END IF
      HSATW=Y*HFCR

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION HFG(TC)

! ----------------------------------------------------------------------
!
!  Latent heat of vaporization of water (KJ/KG) given sat. T (C)
!
!***********************************************************************

      DATA E1,E2,E3,E4,E5/-3.87446,2.94553,-8.06395,11.5633,-6.02884/
      DATA B,C,D,HFGTP,TCR/0.779221,4.62668,-1.07931,2500.9,647.3/
      DATA TCNV/273.15/

      IF(TC.LT.0.) THEN
          TC=0.
      END IF
      TR=(TCR-TC-TCNV)/TCR
      IF(TR.LT.0.) THEN
          HFG=0.
          RETURN
      ENDIF
      Y=B*TR**(1./3.)+C*TR**(5./6.)+D*TR**0.875
      Y=Y+TR*(E1+TR*(E2+TR*(E3+TR*(E4+TR*E5))))
      HFG=Y*HFGTP

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION HSATS(TC)

! ----------------------------------------------------------------------
!
!  Enthalpy of saturated steam (KJ/KG) given sat. T (C)
!
!***********************************************************************

      DATA E1,E2,E3,E4/-4.81351884,2.69411792,-7.39064542,10.4961689/
      DATA E5,B,C,D/-5.46840036,0.457874342,5.08441288,-1.48513244/
      DATA A,TCR,HCR,TCNV/1.0,647.3,2099.3,273.15/

      TR=(TCR-TC-TCNV)/TCR
      Y=A+B*TR**(1./3.)+C*TR**(5./6.)+D*TR**0.875
      Y=Y+TR*(E1+TR*(E2+TR*(E3+TR*(E4+TR*E5))))
      HSATS=Y*HCR

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION SSATW(TC)

! ----------------------------------------------------------------------
!
!  Sat. entropy of liquid water [KJ/(KG K)] given sat. T (C)
!
!***********************************************************************

      DATA E11,E21,E31,SCR/-1836.92956,14706.6352,-43146.6046,4.4289/
      DATA E41,E51,E61,TCNV/48606.6733,7997.5096,-58333.9887,273.15/
      DATA E71,A2,E12,TCR/33140.0718,0.912762917,-1.75702956,647.3/
      DATA E22,E32,E42/1.68754095,5.82215341,-63.3354786/
      DATA E52,E62,E72,A3/188.076546,-252.344531,128.058531,1.0/
      DATA B3,C3,D3/-0.324817650,-2.990556709,3.2341900/
      DATA E13,E23/-0.678067859,-1.91910364/

      TK=TC+TCNV
      TR=(TCR-TK)/TCR
      IF(TK.LT.300)THEN
          Y=TR*(E11+TR*(E21+TR*(E31+TR*(E41+TR*(E51+TR*(E61+TR*E71))))))
      ELSE IF(TK.LT.600)THEN
          Y=TR*(E12+TR*(E22+TR*(E32+TR*(E42+TR*(E52+TR*(E62+TR*E72))))))
          Y=Y+A2
      ELSE
        Y=A3+B3*TR**(1./3.)+C3*TR**(5./6.)+D3*TR**0.875+TR*(E13+TR*E23)
      ENDIF
      SSATW=Y*SCR

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION SSATS(TC)

! ----------------------------------------------------------------------
!
!  Entropy of saturated steam [KJ/(KG K)] given sat. T (C)
!
!***********************************************************************

      DATA E1,E2,E3,E4,E5/-4.34839,1.34672,1.75261,-6.22295,9.99004/
      DATA A,B,C,D,TCR,SCR/1.0,0.377391,-2.78368,6.93135,647.3,4.4289/
      DATA TCNV/273.15/

      TR=(TCR-TC-TCNV)/TCR
      Y=A+B*TR**(1./3.)+C*TR**(5./6.)+D*TR**0.875
      Y=Y+TR*(E1+TR*(E2+TR*(E3+TR*(E4+TR*E5))))
      SSATS=Y*SCR

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION VS(PKPA,TC)

! ----------------------------------------------------------------------
!
!  Specific volume of superheated steam (M3/KG) given P (KPA) and T (C)
!
!***********************************************************************

      DATA R,B1,B2,B3,EM/4.61631E-4,5.27993E-2,3.75928E-3,0.022,40./
      DATA A0,A1,A2,TCNV/-3.741378,-4.7838281E-3,1.5923434E-5,273.15/
      DATA A3,C1,C2,C3,PCNV/10.,42.6776,-3892.70,-9.48654,0.001/
      DATA C4,C5,C6/-387.592,-12587.5,-15.2578/

      P=PKPA*PCNV
      T=TC+TCNV
      TS=C1+C2/(ALOG(P)+C3)
      IF(P .GE. 12.33) THEN
          TS=C4+C5/(ALOG(P)+C6)
      END IF

      VS=R*T/P-B1*EXP(-B2*T)+(B3-EXP(A0+TS*(A1+TS*A2)))/(A3*P) * EXP((TS-T)/EM)

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION HS(PKPA,TC)

! ----------------------------------------------------------------------
!
!  Enthalpy of superheated steam (KJ/KG) given P (KPA) and T (C)
!
!***********************************************************************

      DATA B11,B12,B13,B21/2041.21,-40.4002,-0.48095,1.610693/
      DATA B22,B23,B31/5.472051E-2,7.517537E-4,3.383117E-4/
      DATA B32,B33,B41,B42/-1.975736E-5,-2.87409E-7,1707.82,-16.99419/
      DATA B43,B44,B45,EM/6.2746295E-2,-1.0284259E-4,6.4561298E-8,45./
      DATA C1,C2,C3,PCNV/42.6776,-3892.70,-9.48654,0.001/
      DATA C4,C5,C6,TCNV/-387.592,-12587.5,-15.2578,273.15/

      P=PKPA*PCNV
      T=TC+TCNV
      TS=C1+C2/(ALOG(P)+C3)
      IF(P .GE. 12.33) THEN
          TS=C4+C5/(ALOG(P)+C6)
      END IF

      A0=B11+P*(B12+P*B13)
      A1=B21+P*(B22+P*B23)
      A2=B31+P*(B32+P*B33)
      A3=B41+TS*(B42+TS*(B43+TS*(B44+TS*B45)))

      HS=A0+T*(A1+T*A2)-A3*EXP((TS-T)/EM)

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION SS(PKPA,TC)

! ----------------------------------------------------------------------
!
!  Entropy of superheated steam [KJ/(KG K)] given P (KPA) and T (C)
!
!***********************************************************************

      DATA A0,A1,A2,A3/4.6162961,1.039008E-2,-9.873085E-6,5.43411E-9/
      DATA A4,B1,B2,B3,C0/-1.170465E-12,-0.4650306,0.001,10.,1.777804/
      DATA C1,C2,C3,EM/-1.802468E-2,6.854459E-5,-1.184434E-7,85./
      DATA C4,E1,E2,E3/8.142201E-11,42.6776,-3892.70,-9.48654/
      DATA E4,E5,E6,TCNV/-387.592,-12587.5,-15.2578,273.15/

      P=PKPA*B2
      T=TC+TCNV
      TS=E1+E2/(ALOG(P)+E3)
      IF(P .GE. 12.33) THEN
          TS=E4+E5/(ALOG(P)+E6)
      END IF

      SS=A0+T*(A1+T*(A2+T*(A3+T*A4)))+B1*ALOG(B2+P*B3)-EXP((TS-T)/EM)* (C0+TS*(C1+TS*(C2+TS*(C3+TS*C4))))

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION TPSS(P,S)

! ----------------------------------------------------------------------
!
!  Temperature (C) of steam,  given P (KPA) and S [KJ/(KG K)]
!
!***********************************************************************

      DATA E1,E2,E3,PCNV/42.6776,-3892.70,-9.48654,0.001/
      DATA E4,E5,E6,TABS/-387.592,-12587.5,-15.2578,273.15/
      CHARACTER(LEN=45),PARAMETER :: FMT_1 = "(' WARNING: FUNCTION TPSS FAILS TO CONVERGE')"

!  compare input entropy with saturation value

      TO=E1-TABS+E2/(ALOG(P*PCNV)+E3)
      IF(P .GE. 12330.) THEN
          TO=E4-TABS+E5/(ALOG(P*PCNV)+E6)
      END IF
      SO=SSATS(TO)
      IF(SO.GE.S)THEN
          TPSS=TO
          RETURN
      ENDIF

!  Initial guess TA is based on assumption of constant specific heat.
!  Subsequent approximations made by interpolation.

      TA=(TO+TABS)*(1.+(S-SO)/CPS(TO))-TABS
      SA=SS(P,TA)
      DO I=1,10
          T=TA+(TO-TA)*(S-SA)/(SO-SA)
          IF(ABS(T-TA).LT.0.05) THEN
              TPSS=T
              RETURN
          END IF

          TO=TA
          SO=SA
          TA=T
          SA=SS(P,TA)
      END DO

      WRITE(1,FMT_1)
      
      TPSS=T

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION CPS(T)

! ----------------------------------------------------------------------
!
!       Determine specific heat of steam, Cp, (KJ/Kg/K) given Temp. (C)
!
!       Specific heat equation from "Fundamentals of Classical
!       Thermodynamics-SI Version" by Van Wylen and Sonntag
!       Table A.9, pg. 683.
!
!       Valid for T between 300-3500 K   max error = .43%
!
!***********************************************************************

      DATA C1,C2,C3,C4,E1,E2/143.05,-183.54,82.751,-3.6989,.25,.5/
      
      CHARACTER(LEN=46),PARAMETER :: FMT_13 = "(' ',' WARNING: FUNCTION CPS: T OUT OF RANGE')"

      TK=T+273.15
      T1=TK/100.
      CPS=(C1+C2*T1**E1+C3*T1**E2+C4*T1)/18.015
      IF(TK.LT.300..OR.TK.GT.3500.) THEN
          WRITE(1,FMT_13)
      END IF
      
      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION CVS(V,T)

! ----------------------------------------------------------------------
!
!       Calculate Cv (KJ/Kg/K) given V (m3/kg) and T (C)
!
!***********************************************************************

      DIMENSION A(7)

      DATA TC,TFR,B1/1165.11,459.67,.0063101/
      DATA A/.99204818,-33.137211,416.29663,.185053,5.475,-2590.5815,113.95968/

      TR=9./5.*T+32.+TFR
      VE=(V-B1)/.062428
      CVS=A(1)+A(2)/SQRT(TR)+A(3)/TR-A(4)*A(5)**2*TR/TC**2*EXP(-A(5)*TR/TC)*(A(6)/VE+A(7)/VE**2)
      CVS=CVS*4.1868

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION VISSV(P)

! ----------------------------------------------------------------------
!
!  Calculates the dynamic viscosity (kg/m-s) of saturated
!  vapor given the pressure (kPa).  Correlation obtained from
!  a curve fit of data from 'Heat Transfer' by Alan J. Chapman, 1974.
!
!***********************************************************************

      DATA C1,C2,C3,C4/.0314,2.9675E-05,-1.60583E-08,3.768986E-12/

!       Covert pressure from kPa to psi

      PSI=P/6.894757
      VISSV=C1+C2*PSI+C3*PSI**2+C4*PSI**3

!       Convert viscosity from lbm/ft-hr to kg/m-s

      VISSV=VISSV*4.1338E-04

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION VISSPH(T)

! ----------------------------------------------------------------------
!
!   Calculates the dynamic viscosity (kg/m-s) of superheated
!   steam given the temperature (C).  The correlation is obtained
!   from a curve fit at atmospheric pressure from 'Heat Transfer'
!   by Alan J. Chapman, 1974. (Note: there is little  variation in
!   viscosity at higher pressures.)
!
!***********************************************************************

      DATA C1,C2,C3,C4/.0183161,5.7067E-05,-1.42253E-08,7.241555E-12/

!     Convert temperature from C to F

      TF=T*1.8+32.
      VISSPH=C1+C2*TF+C3*TF**2+C4*TF**3

!     Convert viscosity from lbm/ft-hr to kg/m-s

      VISSPH=VISSPH*4.1338E-04

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION STEAMK(T)

! ----------------------------------------------------------------------
!
!     Calculates thermal conductivity of superheated steam (KW/m-C)
!     given the temperature (C).  Curve fit from data in 'Heat Transfer'
!     by Alan J. Chapman, 1974.
!
!***********************************************************************

      DATA C1,C2,C3/.824272,.00254627,9.848539E-08/

!       Convert temperature from C to F

      TF=T*1.8+32.
      STEAMK=(C1+C2*TF+C3*TF**2)*.01

!       Convert K from Btu/hr-ft-F to kW/m-C

      STEAMK=STEAMK*0.0017308

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION WRHO(TW)

! ----------------------------------------------------------------------
!
!  Density eq. for water at 1 atm., from CRC Handbook of Chem. & Phys.,
!   61st Edition (1980-1981), p. F-6.  Density (kg/m3) given temp. (C).
!
!***********************************************************************

      DATA AR0,AR1,AR2,AR6/999.83952,16.945176,-.0079870401,.01687985/
      DATA AR3,AR4,AR5/-46.170461E-06,105.56302E-09,-280.54253E-12/

      WRHO=(AR0+TW*(AR1+TW*(AR2+TW*(AR3+TW*(AR4+TW*AR5)))))/(1.+AR6*TW)

      RETURN
      END FUNCTION

!***********************************************************************

      FUNCTION WMU(TW)

! ----------------------------------------------------------------------
!
!  Viscosity equations for water at 1 atm., from CRC Handbook (op.cit.),
!    page F-51.  WMU in kg/meter-sec; for centipoise, multiply by 1000.
!    For temps > 100 C, fit to data from Karlekar & Desmond (saturated).
!
!***********************************************************************

      DATA AM0,AM1,AM2,AM3,AM4/-3.30233,1301.,998.333,8.1855,.00585/
      DATA AM5,AM6,AM7,AM8/1.002,-1.3272,-0.001053,105./
      DATA A10,A11,A12,A13/.68714,-.0059231,2.1249E-05,-2.69575E-08/

      WMU=AM5*10.**((TW-20.)*(AM6+(TW-20.)*AM7)/(TW+AM8))
      IF(TW.LT.20.) THEN
          WMU=10.**(AM0+AM1/(AM2+(TW-20.)*(AM3+AM4*(TW-20.))))*100.
      END IF
      IF(TW.GT.100.) THEN
          WMU=A10+TW*(A11+TW*(A12+TW*A13))
      END IF
      WMU=0.001*WMU

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION WK(TW)

! ----------------------------------------------------------------------
!
!  Thermal conductivity equation from linear least-squares fit to data
!   in CRC Handbook (op.cit.), page E-11; temps. from 270 K to 620 K.
!    Temperature in Celsius, WK in [kW/(m K)].  Values at one atmosphere
!    for T from 0 to 100 C, at saturation for T above 100.
!
!***********************************************************************

      DATA AK0,AK1,AK2/0.560101,0.00211703,-1.05172E-05/
      DATA AK3,AK4/1.497323E-08,-1.48553E-11/

      WK=0.001*(AK0+TW*(AK1+TW*(AK2+TW*(AK3+TW*AK4))))

      RETURN
      END FUNCTION
!***********************************************************************

      FUNCTION WCP(TW)

! ----------------------------------------------------------------------
!
!  Specific heat of water at 1 atmosphere, 0 to 100 C.  Equation from
!    linear least-squares regression of data from CRC Handbook (op.cit.)
!    page D-174; in J/g-C (or kJ/kg-C).
!    For temps > 100, fit to data from Karlekar & Desmond (saturated).
!
!***********************************************************************

      DATA ACP0,ACP1,ACP2/4.21534,-0.00287819,7.4729E-05/
      DATA ACP3,ACP4/-7.79624E-07,3.220424E-09/
      DATA ACP5,ACP6,ACP7,ACP8/2.9735,.023049,-.00013953,3.092474E-07/

      WCP=ACP0+TW*(ACP1+TW*(ACP2+TW*(ACP3+TW*ACP4)))
      IF(TW.GT.100.) THEN
          WCP=ACP5+TW*(ACP6+TW*(ACP7+TW*ACP8))
      END IF

      RETURN
      END FUNCTION
