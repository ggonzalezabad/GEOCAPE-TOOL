! NAME:
!	CALDAT
!
! PURPOSE:
!	Return the calendar date and time given julian date.
!	This is the inverse of the function JULDAY.
!
! INPUTS:
!	JULIAN contains the Julian Day Number (which begins at noon) of the
!	specified calendar date.  It should be a long integer.
!
! OUTPUTS:
!	(Trailing parameters may be omitted if not required.)
!	MONTH:	Number of the desired month (1 = January, ..., 12 = December).
!
!	DAY:	Number of day of the month.
!
!	YEAR:	Number of the desired year.
!
!	HOUR:	Hour of the day
!
!	Minute: Minute of the day
!
!	Second: Second (and fractions) of the day.
!
! RESTRICTIONS:
!	Accuracy using IEEE double precision numbers is approximately
!	1/10000th of a second.
!
! MODIFICATION HISTORY:
!	Translated from "Numerical Recipies in C", by William H. Press,
!	Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling.
!	Cambridge University Press, 1988 (second printing).
!
!	DMS, July 1992.
!	DMS, April 1996, Added HOUR, MINUTE and SECOND keyword
!	AB, 7 December 1997, Generalized to handle array input.
!	AB, 3 January 2000, Make seconds output as DOUBLE in array output.
!      CT, Nov 2006: For Hour/Min/Sec, tweak the input to make sure hours
!      and minutes are correct. Restrict hours to 0-23 & min to 0-59.
!
!      Translated into Fortran90 by Vijay Natraj, JPL, February 24 2012

SUBROUTINE CALDAT(julian, month, day, year, hour, minute, second)

implicit none

!  Inputs
  
real(kind=8), intent(in)      :: julian

!  Outputs

integer(kind=4), intent(out)  :: month
integer(kind=4), intent(out)  :: day
integer(kind=4), intent(out)  :: year
integer(kind=4), intent(out)  :: hour
integer(kind=4), intent(out)  :: minute
integer(kind=4), intent(out)  :: second
   
!  Local variables

real(kind=8)                  :: min_julian
real(kind=8)                  :: max_julian
integer(kind=4)               :: igreg
integer(kind=4)               :: julLong
integer(kind=4)               :: jalpha
integer(kind=4)               :: ja
integer(kind=4)               :: jb
integer(kind=4)               :: jc
integer(kind=4)               :: jd
integer(kind=4)               :: je
real(kind=8)                  :: fraction
real(kind=8)                  :: eps

min_julian = -1095.d0
max_julian = 1827933925.d0
IF ((julian .LT. min_julian) .OR. (julian .GT. max_julian)) THEN
   write(0,*) 'Value of Julian date is out of allowed range.'
ENDIF

igreg = 2299161                 ! Beginning of Gregorian calendar
julLong = FLOOR(julian + 0.5d0) ! Better be long

IF (julLong .GE. igreg) THEN    ! Gregorian
   jalpha = FLOOR(((julLong - 1867216) - 0.25d0) / 36524.25d0)
   ja = julLong + 1 + jalpha - FLOOR(0.25d0 * jalpha)
ELSE IF (julLong .LT. 0) THEN
   ja = julLong + 36525 * (1-julLong/36525)
ELSE
   ja = julLong
ENDIF

jb = ja + 1524
jc = FLOOR(6680.d0 + ((jb-2439870)-122.1d0)/365.25d0)
jd = FLOOR(365.d0 * jc + (0.25d0 * jc))
je = FLOOR((jb - jd) / 30.6001d0)

day = jb - jd - FLOOR(30.6001d0 * je)
month = je - 1
month = MOD(month-1,12) + 1
year = jc - 4715
IF (month .GT. 2) year = year - 1
IF (year .LE. 0) year = year - 1
IF (julLong .LT. 0) year = year - 100 * (1-julLong/36525)

!  hours, minutes, seconds

fraction = julian + 0.5d0 - julLong
IF (ABS(julLong) .GE. 1) THEN
  eps = 1.d-12*ABS(julLong)
ELSE
  eps = 1.d-12
ENDIF
hour = FLOOR(fraction * 24.d0 + eps)
fraction = fraction - hour/24.d0
minute = FLOOR(fraction*1440.d0 + eps)
second = (fraction - minute/1440.d0)*86400.d0

RETURN
END
