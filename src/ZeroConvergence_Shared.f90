
! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This function looks like it checks each iteration for convergence; it's a root solver.
! It's looking for the zero value. We believe this module is adapted from the ZBRENT function
! in Numerical Recipes (1996); thus it uses the Van Wijngaarden-Dekker-Brent method (Section 9.3).

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This component does not represent any physical item.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This function looks like it checks each iteration for convergence, and loops until the calculation converges.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! There are no input or output files directly connected to this function.

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! There is no module level; this is a single function.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 1 method:
!    PUBLIC ZeroConvergence -- Calls GUESS3 and iterates until convergence.
!      Called by HPdesignMod.f90

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! A better look needs to be taken at why the root bracketing is exponential.
! "upperBoundVal = upperBoundVal + 2.0**(ICOUNT-1)*DX*SIGN"

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-29 | JEH | Filled out header
! 2013-12-18 | RAS | Updated Issues & To-Do
! 2013-12-28 | Karthik | Updated the algorithm and some comments.

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! The root bracketing needs to be looked at. In addition, some documentation would be good.
!
! Brent's method is a root-finding algorithm which combines root bracketing, bisection, and inverse quadratic interpolation.
! It is sometimes known as the van Wijngaarden-Deker-Brent method.
! This link below takes you to a soft copy of the book which explains further on how method works and
! how it is implemented.
! Link is - http://www.fing.edu.uy/if/cursos/fiscomp/extras/numrec/book/f9.pdf
! I strongly recommend study the approach with Wikipedia first and then
! Once you have got a good understanding then you can study in above link. (Karthik - 12/28/2013)
! ************************************** !
! -- REFERENCES ------------------------ !
! -------------------------------------- !
! The routine is modified from the ZBRENT function given in Section 9.3
! of Numerical Recipes (1996). Some additional comments have been taken
! from or influenced by Forsythe et al. (1977).
! (List Numerical Recipes, Forsythe, etc.)

!      /*---------------------------------------------------------------------
!        |  METHOD: ZeroConvergence
!        |
!        |  PURPOSE:  Calculate the Root of a Equation.
!        |
!        |  Pre-condition: For any Given Function F(x) with the initial Interval [a0,b0], such that the conditions for the Brent Method are Satisfied.
!        |
!        |  Post-condition: [What we can expect to exist or be true after
!        |      this method has executed under the pre-condition(s).] None
!        |
!        |  PARAMETERS:
!        |            lowerBoundValue - Lower Bound Value for the Function
!        |            FunctionPointer - Address of the Function to be called
!        |            TOL1 -  Tolerence Level
!        |            TOL2 -  Tolerence Level
!        |            DX -    Step Increment for the Guess Method
!        |            FB -
!        |            IERROR - Error Code
!        |
!        |  Returns: Root Value for the Equation.
!        |
!        | LAST UPDATED - Karthik | 12/28/2013 | Updated the Algorithm, Increased the Tolerence Level (Testing)
!        *-------------------------------------------------------------------*/

REAL FUNCTION ZeroConvergence(lowerBoundValue,FunctionPointer,TOL1,TOL2,DX,FB,IERROR)
  USE DataSimulation !, ONLY: CoilMode  !RS: Debugging: This is for a test below (11/8/13)
  implicit none
  !
  LOGICAL FIRST
  INTEGER IterationsMAX !Maximum Number of Iterations before we Return from the Method
  REAL lowerBoundValue,FunctionValueForLowerBound,upperBoundValue,FunctionValueForUpperBound,FunctionPointer,TOL1,TOL2
  REAL A,B,C,D,E,EPS,FA,FB,FC,TOLX,TOLF,XM,P,Q,R,S,EPSDEFAULT
  REAL DX, ONE, EPS0
  INTEGER IERROR, IERR, iter
  EXTERNAL FunctionPointer
  PARAMETER (IterationsMAX=100,EPSDEFAULT=3.e-8)
  DATA FIRST / .TRUE. /

  !Karthik - Call Guess Method to Get the UpperBound for the selected Function Pointer Method.
  CALL GUESS3(lowerBoundValue,FunctionValueForLowerBound,upperBoundValue,FunctionValueForUpperBound,FunctionPointer,DX,TOL2,IERROR)
  !RS: Initial guesses and bounding of region
  IF (IERROR .NE. 0) THEN
    RETURN
  END IF

  !Karthik - Calculate the Machine Precision of the Current Machine.
  !This Value is used to calculate the Tolerence Later.
  !IF (FIRST) THEN
  !    EPS = 1.0
  !    EPS = EPS/2.0
  !    ONE = 1.0 + EPS
  !    DO WHILE (ONE .GT. 1.0)
  !        EPS = EPS/2.0
  !        ONE = 1.0 + EPS0
  !    END DO
  !
  !    FIRST = .FALSE.
  !END IF

  !6/3/2014 - Karthik - As per Rachel request, commenting the Machine precision calculation
  !The tolereance level for EPS has been set to 1/1000
  EPS=1.0/1000

  !Karthik - We use the Same Variable Names as with the original recipe Book, for clarity and easy understanding.

  A = lowerBoundValue
  B = upperBoundValue
  FA = FunctionValueForLowerBound
  FB = FunctionValueForUpperBound

  !Karthik - Check if the Root is Bracketed or NOT and display error message to the user.

  IF ((FA .GT.0. .AND. FB .GT. 0.) .OR. (FA .LT.0. .AND. FB .LT. 0.)) THEN
    !Karthik-Root is NOT Bracketed. One of the possibility of this is a wrong calculation of the upperBoundValue from the Guess3 Method
  END IF

  C = B
  FC = FB

  !DO WHILE (.TRUE.)
  DO iter=1,IterationsMAX !Karthik - Testing Limit the Number of Iterations to 100 to save computing time.

    IF ((FB .GT.0. .AND. FC .GT. 0.) .OR. (FB .LT.0. .AND. FC .LT. 0.)) THEN ! Rename a, b, c and adjust bounding interval d.
      C = A
      FC = FA
      D = B-A
      E = D
    END IF

    IF (ABS(FC) .LT. ABS(FB)) THEN
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
    END IF

    !TOLX = TOL1      !ISI - 05/31/05
    TOLX = 2.*EPS*ABS(B)+0.5*TOL1  !Convergence check.      !Karthik: Testing , Increasing the Tolerence Value 12/28/2013
    TOLF = TOL2      !ISI - 05/31/05
    XM = 0.5*(C - B) !ISI - 05/31/05

    !RS: Debugging: According to the original ORNL flowchart, for the evaporator case, it should
    ! only exit out if both terms are converged. Therefore, the following commented out lines are
    ! a test of this. It didn't seem to make a huge difference, so it's been commented back out
    ! for now.
    !RS: Debugging: Testing (11/8/13)
    !IF (CoilMode .EQ. 1) THEN    !RS: If this is the low-side (evaporator) loop, then
    !    IF (ABS(XM) .LE. TOLX .AND. ABS(FB) .LE. TOLF) THEN
    !        EXIT   !RS: Both convergence criteria have to be met
    !    END IF
    !ELSE   !RS: Otherwise, if it's the high-side loop, then only one of the two
    ! convergence criteria has to be met.
    !    IF (ABS(XM) .LE. TOLX) THEN
    !        EXIT
    !    END IF
    !    IF (ABS(FB) .LE. TOLF) THEN
    !        EXIT
    !    END IF
    !END IF

    !Karthik - Need to check the Criteria for TOLF.
    IF (ABS(XM) .LE. TOLX .OR. FB.EQ.0. ) THEN    !RS: Debugging: Comment out when testing the above (11/8/13)
      EXIT
    END IF
    IF (ABS(FB) .LE. TOLF) THEN
      EXIT
    END IF
    !Karthik - EXIT gives a value for ZeroConvergence and return to calling method.

    IF ((ABS(E) .GE. TOLX) .OR. (ABS(FA) .GT. ABS(FB))) THEN

      S = FB/FA
      IF (A .EQ. C) THEN  !Linear Interpolation
        P = 2.0*XM*S
        Q = 1.0 - S
      ELSE    !Attempt Inverse Quadratic Interpolation.
        Q = FA/FC
        R = FB/FC
        S = FB/FA
        P = S*(2.0*XM*Q*(Q - R) - (B - A)*(R - 1.0))
        Q = (Q - 1.0)*(R - 1.0)*(S - 1.0)
      END IF

      IF (P .GT. 0.0) THEN
        Q = -Q      !Check whether in bounds
      END IF
      P = ABS(P)

      IF (((2.0*P) .LT. (3.0*XM*Q - ABS(TOLX*Q))) .OR. (P .LT. ABS(E*Q))) THEN
        E = D       !Accept interpolation.
        D = P/Q
      ELSE
        D = XM      !Interpolation failed, use bisection.
        E = D
      END IF
    ELSE
      D = XM          !Bounds decreasing too slowly, use bisection.
      E = D
    END IF
    A = B   !Move last best guess to a.
    FA = FB
    IF (ABS(D) .GT. TOLX) THEN  !Evaluate new trial root.
      B = B + D
    ELSE
      B = B + SIGN(TOLX,XM)
    END IF

    FB = FunctionPointer(B,IERR)

    !Karthik - We have shifted this check to the begining of Loop
    !IF ((FB*(FC/ABS(FC))) .GT. 0.) THEN
    !    C = A
    !    FC = FA
    !    D = B-A
    !    E = D
    !END IF

  END DO

  !   SET ERROR CODES
  !     IERROR = 0, NORMAL RETURN
  !     IERROR = 1, TOLERANCE ON INDEPENDENT VARIABLE EXCEEDED
  !     IERROR = 2, TOLERANCE ON FUNCTION VALUE EXCEEDED
  !     IERROR = 3, TOLERANCES ON INDEPENDENT VARIABLE AND FUNCTION VALUE EXCEEDED

  ZeroConvergence = B
  IERROR = 0
  IF (ABS(XM) .GT. TOLX) THEN
    IERROR = 1
  END IF
  IF (ABS(FB).GT.TOLF) THEN
    IERROR = IERROR + 2
  END IF
  RETURN
  !
  !   COMPUTE MACHINE PRECISION "EPS"
  !
END FUNCTION
