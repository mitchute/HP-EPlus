! ************************************** !
! ** HEAT PUMP SIMULATION CODE HEADER ** !
! ************************************** !

! ************************************** !
! -- HIGH LEVEL OVERVIEW/DESCRIPTION --- !
! -------------------------------------- !
! This set of subroutines deals with matrix operations: matrix multiplication, inversion, and solving. 

! ************************************** !
! -- PHYSICAL DESCRIPTION -------------- !
! -------------------------------------- !
! There is no physical component representation for this module.

! ************************************** !
! -- SIMULATION DATA RESPONSIBILITIES -- !
! -------------------------------------- !
! This is purely a set of mathematical calculation subroutines.

! ************************************** !
! -- INPUT FILES/OUTPUT FILES (none) --- !
! -------------------------------------- !
! NA

! ************************************** !
! -- MODULE LEVEL VARIABLES/STRUCTURES - !
! -------------------------------------- !
! Nothing is defined at the module level; there is no module.

! ************************************** !
! -- SUMMARY OF METHODS, CALL TREE ----- !
! -------------------------------------- !
! This module contains 5 methods:
!   PUBLIC LUD -- Solves a set of linear equations and returns the Lu decomposition
!      Called by matrix_inverse
!   PUBLIC matrix_invert_step_2 -- Determines the matrix inverse by back substitution
!      Called by matrix_inverse
!   PUBLIC matrix_inverse -- Inverts a matrix
!      Called by CoilCalc.fd0
!   PUBLIC CalcMatrixMultVector -- Multiplies a single-row matrix and single-column matrix
!      Called by CoilCalc.fd0
!   PUBLIC CalcMatrixMultMatrix -- Multiplies non-singular matrices
!      Not called by any module or subroutine

! ************************************** !
! -- ISSUES/BUGS/TICKETS --------------- !
! -------------------------------------- !
! NA

! ************************************** !
! -- CHANGELOG ------------------------- !
! -------------------------------------- !
! 2012-12-11 | ESL | Initial header
! 2012-12-12 | RAS | Updated header

! ************************************** !
! -- TODO/NOTES/RECOMMENDATIONS -------- !
! -------------------------------------- !
! Some more documentation might be useful.

!*******************************************************************************

! This routine converted from a Fortran 77 routine from Wheatley
!  --------------------------------------------------------------
!
! ---------------------------------------------------------------
!
SUBROUTINE Solve(a, n, IPVT, B)

!
!  ---------------------------------------------------------------
!
!      MAKING USE OF THE LU DECOMPOSITION OF THE MATRIX A, THIS SUB-
!      ROUTINE SOLVES THE SYSTEM BY FORWARD AND BACK SUBSTITUTION.
!
!      INPUT:   A - LU MATRIX FROM SUBROUTINE LUD
!               N - NUMBER OF EQUATIONS
!               IPVT - A RECORD OF THE REARRANGEMENT OF THE ROWS
!                      OR A FROM SUBROUTINE LUD
!               B - RIGHT HAND SIDE OF THE SYSTEM OF EQUATIONS
!
!      OUTPUT:  B - THE SOLUTION VECTOR
!
!  -----------------------------------------------------------------
!
      IMPLICIT NONE

      INTEGER n
	  INTEGER, DIMENSION(n) :: IPVT
	  REAL,DIMENSION(n,n) :: a
	  REAL,DIMENSION(n) :: B

	  REAL X(n),SUM
      INTEGER IROW, JCOL, i
!
!  --------------------------------------------------------------
!
!     REARRANGE THE ELEMENTS OF THE B VECTOR. STORE THEM IN THE
!     X VECTOR.
!
      DO i = 1, n
          X(i) = B(IPVT(i))
      END DO
          
!
!  ---------------------------------------------------------------
!
!     SOLVE USING FORWARD SUBSTITUTION--LY = B
!
!  ----------------------------------------------------------------
!
      DO IROW = 2, n
          SUM = X(IROW)
          DO JCOL = 1 , (IROW - 1)
              SUM = SUM - a(IROW, JCOL) * X(JCOL)
          END DO
          X(IROW) = SUM
      END DO
!
!  ----------------------------------------------------------------
!
!     SOLVE BY BACK SUBSTITUTION--UX = Y
!
!  ----------------------------------------------------------------
!
      B(n) = X(n) / a(n, n)
      DO IROW = (n - 1), 1, -1
          SUM = X(IROW)
          DO JCOL = (IROW + 1), n
              SUM = SUM - a(IROW, JCOL) * B(JCOL)
          END DO
          B(IROW) = SUM / a(IROW, IROW)
      END DO
 
END SUBROUTINE

!*******************************************************************************

! This routine converted from a Fortran 77 routine from Wheatley

!
!             SUBROUTINE LUD
!
SUBROUTINE LUD(a, n, IPVT, det)

!
!  ---------------------------------------------------------------
!
!         SUBROUTINE LUD: THIS SUBROUTINE SOLVES A SET OF LINEAR EQUATIONS
!     RETURNS THE LU DECOMPOSITION OF THE COEFFICIENT MATRIX. THE METHOD
!     IS BASED ON THE ALGORITHM PRESENTED IN SECTION 4.
!
!     INPUT:  A  -  THE COEFFICIENT MATRIX
!             N  -  THE NUMBER OF EQUATIONS
!             NDIM  -  THE MAXIMUM ROW DIMENSION OF A
!
!     OUTPUT: A  -  THE LU DECOMPOSITION OF THE MATRIX A
!                   THE ORIGINAL MATRIX A IS LOST
!             IPVT  -  A VECTOR CONTAINING THE ORDER OF THE ROWS OF THE
!                      REARRANGED MATRIX DUE TO PIVOTING
!         DET  -  THE DETERMINANT OF THE MATRIX. IT IS SET TO
!                 0 IF ANY PIVOT ELEMENT IS LESS THAN 0.00001.
!
!  --------------------------------------------------------------------
!
      IMPLICIT NONE

	  INTEGER n
      REAL, DIMENSION(n,n) :: a
	  REAL det
	  INTEGER, DIMENSION(n) :: IPVT

      REAL Temp
      INTEGER i, IPVTMT, NLESS1, IPLUS1, j
      INTEGER KCOL, JCOL, JROW, TMPVT
!
      det = 1.
      NLESS1 = n - 1
      DO i = 1, n
          IPVT(i) = i
      END DO
!
      DO i = 1, NLESS1
          IPLUS1 = i + 1
          IPVTMT = i
!
!         FIND PIVOT ROW
!
          DO j = IPLUS1, n
              If (Abs(a(IPVTMT, i)) < Abs(a(j, i))) THEN
                  IPVTMT = j
              END IF
          END DO
!
!  ------------------------------------------------------------------
!
!     CHECK FOR SMALL PIVOT ELEMENT
!
!  ------------------------------------------------------------------
!
          If (Abs(a(IPVTMT, i)) < 0.00001) Then
              det = 0.0
          End If
!
!         INTERCHANGE ROWS IF NECESSARY
!
          If (IPVTMT /= i) Then
              TMPVT = IPVT(i)
              IPVT(i) = IPVT(IPVTMT)
              DO JCOL = 1, n
                 Temp = a(i, JCOL)
                 a(i, JCOL) = a(IPVTMT, JCOL)
                 a(IPVTMT, JCOL) = Temp
              END DO
              IPVT(IPVTMT) = TMPVT
              det = -det
          End If
!
!                  REDUCE ALL ELEMENTS BELOW THE I'TH ROW
!
            DO JROW = IPLUS1, n
                If (a(JROW, i) /= 0.0) Then
                   a(JROW, i) = a(JROW, i) / a(i, i)
                   DO KCOL = IPLUS1, n
                      a(JROW, KCOL) = a(JROW, KCOL) - a(JROW, i) * a(i, KCOL)
                   END DO
                End If
            END DO
!
       END DO  
!
       If (Abs(a(n, n)) < 0.00001) Then
            ! error...
       End If
!
!  -----------------------------------------------------------------
!
!     COMPUTE THE DETERMINANT OF THE MATRIX
!
!  ------------------------------------------------------------------
!
      DO i = 1, n
          det = det * a(i, i)
      END DO
!
End SUBROUTINE

!*******************************************************************************
!

!This subroutine written by RDD and JDS 9/1/97
!After LU decomposition is performed, this
!routine finds the matrix inverse by back substitution

SUBROUTINE matrix_invert_step_2(a, ainv, pivot, n)
implicit none

INTEGER n, i
REAL, DIMENSION(n,n) :: a, ainv
INTEGER, DIMENSION(n) :: pivot

REAL ivector(n)
INTEGER col

DO col = 1, n

    DO i = 1, n
       ivector(i) = 0.0
    END DO
       
    ivector(col) = 1
    
    Call Solve(a, n, pivot, ivector)
    
    DO i = 1, n
       ainv(i, col) = ivector(i)
    END DO
END DO

End SUBROUTINE

!*******************************************************************************

SUBROUTINE matrix_inverse(a, ainv, n)
implicit none

!
!
! Name:           matrix_inverse
!
! Developer:      Dr. J.D. Spitler
!
!Input variables:
!
!                  a - array holding the nxn matrix to be
!                      inverted
!                   n- the order of the matrix
!Output variables:
!
!                ainv- array holding the matrix inverse

INTEGER n
REAL,DIMENSION(n,n) :: a,ainv

INTEGER pivot(n)
REAL det 

Call LUD(a, n, pivot, det)
Call matrix_invert_step_2(a, ainv, pivot, n)

End SUBROUTINE

!*******************************************************************************

SUBROUTINE CalcMatrixMultVector(M1, Nrow, M2, M3)
implicit none

!PURPOSE OF THIS SUBROUTINE:
    !To calculate the product of two matrices, M3=M1*M2

!SUBROUTINE ARGUMENT DECLARATIONS:
    !Input variables:
        !M1     - Input matrix 1
        !Nrow   - Number of rows
        !M2     - Input matrix 2
    !Output variable:
        !M3 - Product of matrices M1 and M2

INTEGER Nrow
REAL M1(Nrow,Nrow),M2(Nrow),M3(Nrow)
    
!SUBROUTINE INTERNAL VARIABLE DECLARATIONS:
INTEGER i !Loop counter for matrix M1 row element
INTEGER j !Loop counter for matrix M2 column element

!FLOW:

DO i = 1, Nrow
   M3(i) = 0.0
   DO j = 1, Nrow
     M3(i) = M3(i) + M1(i, j) * M2(j)
   END DO
END DO
    
END SUBROUTINE

!*******************************************************************************

SUBROUTINE CalcMatrixMultMatrix(M1, NrowM1, NcolM1, M2, NcolM2, M3)
implicit none

!PURPOSE OF THIS SUBROUTINE:
    !To calculate the product of two matrices, M3=M1*M2

!SUBROUTINE ARGUMENT DECLARATIONS:
    !Input variables:
        !M1     - Input matrix 1
        !NrowM1 - Number of rows for matrix M1
        !NcolM1 - Number of columns for matrix M1
        !M2     - Input matrix 2
        !NcolM2 - Number of columns for matrix M2
    !Output variable:
        !M3 - Product of matrices M1 and M2

INTEGER NrowM1,NcolM1,NcolM2
REAL M1(NrowM1,NcolM1), M2(NcolM1,NcolM2), M3(NrowM1,NcolM2)
    
!SUBROUTINE INTERNAL VARIABLE DECLARATIONS:
    INTEGER i !Loop counter for matrix M1 row element
    INTEGER j !Loop counter for matrix M2 column element
    INTEGER k !Loop counter for matrix M1 column element

!FLOW:
    DO i = 1, NrowM1
        DO j = 1, NcolM2
            M3(i, j) = 0.0
            DO k = 1, NcolM1
                M3(i, j) = M3(i, j) + M1(i, k) * M2(k, j)
            END DO
        END DO
    END DO
    
END SUBROUTINE

!*******************************************************************************
