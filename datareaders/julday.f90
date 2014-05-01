! NAME:
!	JULDAY
!
! PURPOSE:
!	Calculate the Julian Day Number for a given month, day, and year.
!	This is the inverse of the library function CALDAT.
!	See also caldat, the inverse of this function.
!
! INPUTS:
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year. Year parameters must be valid
!              values from the civil calendar. Years B.C.E. are represented
!              as negative integers. Years in the common era are represented
!              as positive integers. In particular, note that there is no
!              year 0 in the civil calendar. 1 B.C.E. (-1) is followed by
!              1 C.E. (1).
!
!	HOUR:	Number of the hour of the day.
!
!	MINUTE:	Number of the minute of the hour.
!
!	SECOND:	Number of the second of the minute (double precision number).
!
! OUTPUTS:
!	JULDAY returns the Julian Day Number (which begins at noon) of the
!	specified calendar date.
!
! RESTRICTIONS:
!	Accuracy using IEEE double precision numbers is approximately
!      1/10000th of a second, with higher accuracy for smaller (earlier)
!      Julian dates.
!
! MODIFICATION HISTORY:
!	Translated from "Numerical Recipies in C", by William H. Press,
!	Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
!	Cambridge University Press, 1988 (second printing).
!
!	AB, September, 1988
!	DMS, April, 1995, Added time of day.
!
!      Translated into Fortran90 by Vijay Natraj, JPL, February 24 2012

function JULDAY(MONTH, DAY, YEAR, Hour, Minute, Second)

implicit none

!  Inputs

integer(kind=4), intent(in)  :: MONTH
integer(kind=4), intent(in)  :: DAY
integer(kind=4), intent(in)  :: YEAR
integer(kind=4), intent(in)  :: HOUR
integer(kind=4), intent(in)  :: MINUTE
real(kind=8),    intent(in)  :: SECOND

!  Outputs

real(kind=8)                 :: JULDAY

!  Local variables

integer(kind=4)              :: GREG
integer(kind=4)              :: min_calendar
integer(kind=4)              :: max_calendar
integer(kind=4)              :: bc
integer(kind=4)              :: L_YEAR
integer(kind=4)              :: inJanFeb
integer(kind=4)              :: JY
integer(kind=4)              :: JM
integer(kind=4)              :: JA
integer(kind=4)              :: JUL
real(kind=8)                 :: eps

! Gregorian Calendar was adopted on Oct. 15, 1582
! skipping from Oct. 4, 1582 to Oct. 15, 1582

GREG = 2299171  ! incorrect Julian day for Oct. 25, 1582

! check if date is within allowed range

min_calendar = -4716
max_calendar = 5000000
IF ((YEAR .LT. min_calendar) .OR. (YEAR .GT. max_calendar)) &
   write(0,*) 'Value of Julian date is out of allowed range'

IF (YEAR .LT. 0) THEN
   bc = 1
ELSE
   bc = 0
ENDIF
L_YEAR = YEAR + bc

IF (MONTH .LE. 2) THEN
   inJanFeb = 1
ELSE
   inJanFeb = 0
ENDIF

JY = YEAR - inJanFeb
JM = MONTH + 1 + 12*inJanFeb

JUL = FLOOR(365.25d0 * JY) + FLOOR(30.6001d0 * JM) + DAY + 1720995

! Test whether to change to Gregorian Calendar.

IF (JUL .GE. GREG) THEN
   JA = FLOOR(0.01d0 * JY)
   JUL = JUL + 2 - JA + FLOOR(0.25d0 * JA)
ENDIF

! Add a small offset so we get the hours, minutes, & seconds back correctly
! if we convert the Julian dates back. This offset is proportional to the
! Julian date, so small dates (a long, long time ago) will be "more" accurate.

! eps = (MACHAR(/DOUBLE)).eps

eps = 2.2204460d-16 ! For Ganesha, calculated from IDL output of above statement

IF (ABS(JUL) .GE. 1) THEN
   eps = eps*ABS(JUL)
ENDIF

! For Hours, divide by 24, then subtract 0.5, in case we have unsigned integers.

JULDAY = real(JUL,kind=8) + real(Hour,kind=8)/24.d0 - &
         0.5d0 + real(Minute,kind=8)/1440.d0 + &
         Second/86400.d0 + eps

RETURN
END FUNCTION JULDAY
