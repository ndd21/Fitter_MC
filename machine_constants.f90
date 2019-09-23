MODULE machine_constants
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(precision(1.d0),range(1.d0))

CONTAINS

FUNCTION imdcon(k) RESULT(ival)

IMPLICIT NONE
INTEGER, INTENT(IN) :: k
INTEGER             :: ival

!  ***  return integer machine-dependent constants  ***

!     ***  k = 1 means return standard output unit number.   ***
!     ***  k = 2 means return alternate output unit number.  ***
!     ***  k = 3 means return  input unit number.            ***
!          (note -- k = 2, 3 are used only by test programs.)

INTEGER :: mdcon(3) = (/ 6, 0, 5 /)

ival = mdcon(k)

RETURN
!  ***  last card of imdcon follows  ***
END FUNCTION imdcon


FUNCTION rmdcon(k) RESULT(fn_val)

!  ***  return machine dependent constants used by nl2sol  ***

IMPLICIT NONE
INTEGER, INTENT(IN) :: k
REAL (dp)           :: fn_val

!  ***  the constant returned depends upon k...

!  ***        k = 1... smallest pos. eta such that -eta exists.
!  ***        k = 2... square root of eta.
!  ***        k = 3... unit roundoff = smallest pos. no. machep such
!  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1.
!  ***        k = 4... square root of machep.
!  ***        k = 5... square root of big (see k = 6).
!  ***        k = 6... largest machine no. big such that -big exists.

SELECT CASE (k)
  CASE(1)
    fn_val = TINY(1._dp)
  CASE(2)
    fn_val = SQRT( TINY(1._dp) )
  CASE(3)
    fn_val = EPSILON(1._dp)
  CASE(4)
    fn_val = SQRT( EPSILON(1._dp) )
  CASE(5)
    fn_val = SQRT( HUGE(1._dp) )
  CASE(6)
    fn_val = HUGE(1._dp)
  CASE DEFAULT
    WRITE(*, *) '** Illegal argument, k = ', k, ' to FUNCTION rmdcon **'
END SELECT

RETURN
!  ***  last card of rmdcon follows  ***
END FUNCTION rmdcon

END MODULE machine_constants
