! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This is a multi-purpose iteration routine; there is no known method used.
!
! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! This is not a physical component

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This provides an iteration scheme for other modules and components within the HP system.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! NONE

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! NONE

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
!    GUESS3
!       ZERO3
!       ZeroConvergence

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! A better look needs to be taken at why the root bracketing is exponential.
! "upperBoundValueGuess = upperBoundValueGuess + 2.0**(ICOUNT-1)*DX*SIGN"

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-29 | JEH | Filled out header
! 2013-12-18 | RAS | Updated Issues & To-Do
! 2014-01-14 | Karthik | Updated Variable Names and comments

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! The root bracketing needs to be looked at. In addition, some documentation would be good.


!      /*---------------------------------------------------------------------
!        |  METHOD: GUESS3
!        |
!        |  PURPOSE:  Guess the Other Bound Value for the Solver (Usually UpperBound Value)
!        |     
!        |  Pre-condition: 
!        |
!        |  Post-condition: 
!        |
!        |  PARAMETERS:
!        |            lowerBoundValueGuess - Lower Bound Value for the Function 
!        |            FunctionValForLowerBoundGuess - F(X) Value for the LowerBound Value
!        |            upperBoundValueGuess -  Upper Bound Value to be Guessed  
!        |            FunctionValForUpperBoundGuess -  F(X) Value for the Upper Bound Value
!        |            FunctionPointerAddress -   Address to the Method to be Called
!        |            DX - Step Value 
!                     TOL - Tolerence Value  
!        |            IERROR - Error Code
!        |
!        |  Returns: Upper and Lower Bound Value to Calling Method (ZeroConvergence). By Modifying the Input parameters values
!        | 
!        | LAST UPDATED - Karthik | 1/14/2014 | Updated the Naming of Variables and comments.
!        *-------------------------------------------------------------------*/


 
 SUBROUTINE GUESS3(lowerBoundValueGuess,FunctionValForLowerBoundGuess,upperBoundValueGuess,FunctionValForUpperBoundGuess,FunctionPointerAddress,DX,TOL,IERROR)
      implicit none

      INTEGER IERR, ICOUNT, IERROR
      REAL SIGN, SLOPE
      REAL lowerBoundValueGuess, upperBoundValueGuess, FunctionValForLowerBoundGuess, FunctionValForUpperBoundGuess, DX, TOL
      REAL YMIN, DDX
      REAL FunctionPointerAddress
        EXTERNAL FunctionPointerAddress
      CHARACTER(LEN=111),PARAMETER :: FMT_1001 = "(' GUESS3: ** FAILED TO BRACKET A SOLUTION **',/,9X,'F(',1PE12.4,') = ',1PE13.5,5X,'F(',1PE12.4,') = ',1PE13.5)"

      ICOUNT = 0
      IERROR = 0
      SIGN = 1.0
      SLOPE = 0.0

      DO WHILE (.TRUE.)
          FunctionValForLowerBoundGuess = FunctionPointerAddress(lowerBoundValueGuess,IERR)
          upperBoundValueGuess = lowerBoundValueGuess                                !Initialize ISI - 03/26/04
          ICOUNT = ICOUNT + 1
          IF (IERR .EQ. 0) THEN
              EXIT
          END IF
          IF (IERR .EQ. 1) THEN                  !ISI - 03/26/04
              lowerBoundValueGuess = lowerBoundValueGuess - DX                         !ISI - 03/26/04
          ELSE                                   !ISI - 03/26/04
              lowerBoundValueGuess = lowerBoundValueGuess + DX                         !ISI - 03/26/04
          END IF                                 !ISI - 03/26/04
          
          IF (ICOUNT .GT. 30) THEN
              IERROR = 4
              RETURN
          END IF

      END DO

      YMIN = ABS(FunctionValForLowerBoundGuess)
      IF (YMIN .LE. TOL) THEN
          RETURN              !ISI - 02/12/06
      END IF
      IF(FunctionValForLowerBoundGuess.GT.0.0) THEN
          SIGN = -1.0
      END IF
      DDX = ABS(DX)
      upperBoundValueGuess = upperBoundValueGuess + 2.0**(ICOUNT-1)*DX*SIGN      !To bracket root ISI - 03/26/04
      
      DO WHILE (.TRUE.)

          FunctionValForUpperBoundGuess = FunctionPointerAddress(upperBoundValueGuess,IERR)    !RS Comment: The compressor saturation temperature seems to be iterating the wrong way with the MC case
          ICOUNT = ICOUNT + 1
          SLOPE = (FunctionValForUpperBoundGuess - FunctionValForLowerBoundGuess)/(upperBoundValueGuess - lowerBoundValueGuess)
          IF (IERR .NE. 0) THEN 
              upperBoundValueGuess = (lowerBoundValueGuess + upperBoundValueGuess)/2.
              IF (ICOUNT .GT. 25) THEN  !15) THEN   !RS: Debugging: Trying to let it iterate more
                  IERROR = 4
                  RETURN
              END IF
              
              CYCLE
          END IF

          YMIN = AMIN1(YMIN,ABS(FunctionValForUpperBoundGuess))
          IF (YMIN .LE. TOL) THEN
              RETURN              !ISI - 02/12/06
          END IF
          DDX = ABS(DX)
          DDX=2.0**(ICOUNT-1)*ABS(DX)				!ISI - 05/10/04
!
          IF (FunctionValForLowerBoundGuess*FunctionValForUpperBoundGuess .LE. 0.) THEN
              RETURN
          END IF
          IF (ICOUNT .GT. 15) THEN
              IERROR = 4
              RETURN
          END IF
!
!	IF SLOPE IS POSITIVE AND "Y" IS GREATER THAN 0. OR
!	IF SLOPE IS NEGATIVE AND "Y" IS LESS THAN 0. MOVE THE LOWER
!	POINT TO THE LEFT
!
          SLOPE = (FunctionValForUpperBoundGuess - FunctionValForLowerBoundGuess)/(upperBoundValueGuess - lowerBoundValueGuess)
          IF (SLOPE*FunctionValForLowerBoundGuess .GE. 0.) THEN
!
!	MAKE "lowerBoundValueGuess" THE LOWER OR LEFT-HAND POINT
!
              IF (lowerBoundValueGuess .GE. upperBoundValueGuess) THEN
                  lowerBoundValueGuess = upperBoundValueGuess
                  FunctionValForLowerBoundGuess = FunctionValForUpperBoundGuess
              END IF
              
              upperBoundValueGuess = lowerBoundValueGuess - DDX
              CYCLE

          END IF

!	IF THERE IS A NEGATIVE SLOPE AND "Y" IS GREATER THAN 0. OR
!	THERE IS A POSITIVE SLOPE AND "Y" IS LESS THAN 0. THEN
!	MOVE THE RIGHT-HAND POINT FURTHER TO THE RIGHT
!
          IF (lowerBoundValueGuess .LE. upperBoundValueGuess) THEN 
!
!	MAKE lowerBoundValueGuess THE UPPER POINT
!
              lowerBoundValueGuess = upperBoundValueGuess
              FunctionValForLowerBoundGuess = FunctionValForUpperBoundGuess
          END IF

          upperBoundValueGuess = lowerBoundValueGuess + DDX
          CYCLE   !VL: Not needed ...

      END DO

 END SUBROUTINE
