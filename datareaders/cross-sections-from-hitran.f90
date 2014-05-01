PROGRAM cross_sections_from_hitran

IMPLICIT none

INTEGER, PARAMETER :: maxlines  = 60000        ! spectral lines
INTEGER, PARAMETER :: maxmols   = 39           ! number of molecules in hitran
INTEGER, PARAMETER :: maxiso    = 8            ! maximum # of isotopes
INTEGER, PARAMETER :: maxpoints = 600000       ! number of points in spectrum

! constants
REAL(KIND=8), PARAMETER :: pi = 3.14159265358979d0
REAL(KIND=8), PARAMETER :: c = 2.99792458d10
REAL(KIND=8), PARAMETER :: p0 = 1013.25d0
REAL(KIND=8), PARAMETER :: t0 = 296.d0          ! hitran standard
! codata 2002 constants
REAL(KIND=8), PARAMETER :: h = 6.6260693d-27
REAL(KIND=8), PARAMETER :: an = 6.0221415d23
REAL(KIND=8), PARAMETER :: r = 82.057463d0      ! derived
REAL(KIND=8), PARAMETER :: rk = 1.3806505d-16
REAL(KIND=8), PARAMETER :: du = 2.6867773d16    ! derived
REAL(KIND=8), PARAMETER :: c2 = 1.4387752d0

INTEGER, DIMENSION(maxlines)             :: mol, iso
REAL(KIND=8), DIMENSION(maxlines)        :: sigma0, strnth, einstein, alpha, &
     elow, coeff, selfbrdn
REAL(KIND=8), DIMENSION(maxmols)         :: qpower
REAL(KIND=8), DIMENSION(maxmols, maxiso) :: amu
REAL(KIND=8), DIMENSION(maxpoints)       :: pos, spec, specmod, voigtx, v
CHARACTER(LEN=132)                       :: input, output

! Input parameters
INTEGER      :: molnum, npoints, nvoigt, nmod
REAL(KIND=8) :: start, step, press, temp, hw1e

INTEGER      :: i, j, mol_temp, iso_temp, nlines, nvlo, nvhi, idx, MEND
REAL(KIND=8) :: sigma0_temp, strnth_temp, einstein_temp, alpha_temp, &
     selfbrdn_temp, elow_temp, coeff_temp
REAL(KIND=8) :: vg, voigta, ratio1, ratio2, ratio, vnorm, rt0t, rc2t, rc2t0
CHARACTER(LEN=2) :: molc

INTEGER, EXTERNAL :: ibin

write (*, '(5x, a)') 'enter hitran input file'
read (*, '(a)') input
write (*, '(5x, a)') 'enter cross section output file'
read (*, '(a)') output

! read logicals and parameters controlling calculations
!
! molnum is the hitran molecule number.
!
! start, step, and npoints define the calculation grid. in order to avoid
! undersampling the grid should be at least as fine as 1/3 of the smallest
! gaussian fwhm of a spectral line, or 0.5550 times the smallest gaussian hw1e.
! this insures that the maximum sampling error is <1.24e-4 of the full-scale
! line-shape. later: add automatic check for undersampling wrt voigt widths
!
! press is the pressure in atmospheres (millibars / 1013.25). temp is the
! temperature in degrees kelvin.
!
! nvoigt is the number of grid points to each side of a spectral line for
! performing the voigt calculation. assuming calculation to <= 1e-4 is desired,
! nvoigt should be the greater of (1) 100 * hwhm_max / step, where hwhm_max is
! the largest lorentzian hwhm of a line and step is the grid spacing;
! (2) 3.035 * hw1e_max / step, where hw1e_max is the largest gaussian hw1e of a
! line.
!
! hw1e is the gaussian slit with at 1/e intensity. hw1e = the full width at
! half-maximum, fwhm, * 0.60056 = fwhm / (2 * sqrt (ln 2)).
!
! nmod gives the option not to write out humongous spectral files by printing
! every nmod^th spectral value
write (*, '(a)') &
  'read molnum, start, step, npoints, press, temp, nvoigt, hw1e, nmod'
read (*, *) molnum, start, step, npoints, press, temp, nvoigt, hw1e, nmod

! initialize cross sections and calculate the spectrum grid.
spec = 0.d0
specmod = 0.d0
do i = 1, npoints
  pos (i) = start + (i - 1) * step
end do
if (npoints > maxpoints) pause 'npoints > maxpoints!!!'

! setup hitran
CALL hitran_setup (maxmols, maxiso, qpower, amu)
!WRITE(molc, '(I2.2)') molnum
!input='../geocape_data/HITRAN/' // molc // '_hit08.par'
open (unit = 1, file = input, status = 'old')
open (unit = 2, file = output, status = 'unknown')

! read lines

i = 1
DO 
   READ (1, '(i2, i1, f12.6, 2e10.3, 2f5.4, f10.4, f4.2)', IOSTAT=MEND) mol_temp, &
        iso_temp, sigma0_temp, strnth_temp, einstein_temp, alpha_temp, selfbrdn_temp, &
        elow_temp, coeff_temp
   IF (MEND < 0 .OR. sigma0_temp > pos(npoints) + 2.0 * nvoigt * step) EXIT
   !IF (MEND < 0 ) EXIT
   !write (*, *) mol_temp, molnum
   ! only count lines for the specified molecule
   IF (mol_temp == molnum .AND. sigma0_temp > pos(1) - 2.0 * nvoigt * step) then
   !IF (mol_temp == molnum ) then
      mol(i) = mol_temp
      iso(i) = iso_temp
      sigma0(i) = sigma0_temp
      strnth(i) = strnth_temp
      einstein(i) = einstein_temp
      alpha(i) = alpha_temp
      selfbrdn(i) = selfbrdn_temp
      elow(i) = elow_temp
      coeff(i) = coeff_temp
      !WRITE(*, '(I5,3F12.4)') i, sigma0(i), pos(1), pos(npoints)
      i = i + 1
   ENDIF
ENDDO
nlines = i - 1
WRITE (*, *) 'nlines = ', nlines
IF (nlines > maxlines) PAUSE 'nlines > maxlines!!!'

! loop over lines to fill out cross section array
rt0t = t0 / temp; rc2t = c2 / temp; rc2t0 = c2 / t0
DO i = 1, nlines
   vg = 4.30140d-7 * sigma0 (i) * dsqrt (temp / amu (mol (i), iso (i)))
   voigta = press * alpha(i) * (rt0t)**coeff(i) / vg
   ratio1 = dexp(-elow(i) * rc2t) - dexp(-(sigma0(i) + elow(i)) * rc2t)
   ratio2 = dexp(-elow(i) * rc2t0) - dexp(-(sigma0(i) + elow(i)) * rc2t0)
   ratio = ratio1 / ratio2 * ((rt0t)**qpower (mol(i)))
   vnorm = ratio * strnth(i) / vg
   idx = ibin (sigma0(i), pos, npoints)
   nvlo = MAX(1, idx - nvoigt)
   nvhi = MIN(npoints, idx + nvoigt)
   IF (idx == 0 .AND. sigma0(i) < pos(1) - nvoigt * step .OR. sigma0(i) > pos(npoints) + nvoigt * step) CYCLE
   voigtx(nvlo:nvhi) = (pos(nvlo:nvhi) - sigma0(i)) / vg
   CALL voigt (voigtx, voigta, v, npoints, nvlo, nvhi)
   WRITE(*, '(4I6,6D14.5)') i, idx, nvlo, nvhi, voigta*vg, vg!, MINVAL(v(nvlo:nvhi))*vnorm, MAXVAL(v(nvlo:nvhi))*vnorm
   spec(nvlo:nvhi) = spec(nvlo:nvhi) + vnorm * v(nvlo:nvhi)
ENDDO

! convolve with instrument function
CALL gauss (pos, spec, specmod, npoints, hw1e)

! output results
DO i = 1, npoints
   WRITE (2, '(f10.3, 1p2e13.5)') pos(i), spec(i), specmod(i)
ENDDO

CLOSE(unit = 1)
CLOSE(unit = 2)

END PROGRAM cross_sections_from_hitran
!

SUBROUTINE hitran_setup (maxmols, maxiso, qpower, amu)

IMPLICIT NONE
INTEGER, INTENT(IN) :: maxmols, maxiso
REAL(KIND=8), DIMENSION(maxmols), INTENT(OUT)         :: qpower
REAL(KIND=8), DIMENSION(maxmols, maxiso), INTENT(OUT) :: amu

INTEGER :: i

qpower = 1.5d0
amu = 0.d0

! hitran numbers:
!   h2o (1)
!   co2 (2)
!    o3 (3)
!   n2o (4)
!    co (5)
!   ch4 (6)
!    o2 (7)
!    no (8)
!   so2 (9)
!   no2 (10)
!   nh3 (11)
!  hno3 (12)
!    oh (13)
!    hf (14)
!   hcl (15)
!   hbr (16)
!    hi (17)
!   clo (18)
!   ocs (19)
!  h2co (20)
!  hocl (21)
!    n2 (22)
!   hcn (23)
! ch3cl (24)
!  h2o2 (25)
!  c2h2 (26)
!  c2h6 (27)
!   ph3 (28)
!  cof2 (29)
!   sf6 (30)
!   h2s (31)
! hcooh (32)
!   ho2 (33)
!     o (34)
!clono2 (35)
!   no+ (36)
!  hobr (37)
!  c2h4 (38)
! ch3oh (39)
! ch3br (40)
! ch3cn (41)
!   cf4 (42)

do i = 1, maxmols
  if (i .eq. 34) qpower (i) = 0.d0
  if (i .eq. 2 .or. i .eq. 4 .or. i .eq. 5 .or. i .eq. 7 .or. i .eq. 8 &
  .or. i .eq. 13 .or. i .eq. 14 .or. i .eq. 15 .or. i .eq. 16 .or. i &
  .eq. 17 .or. i .eq. 18 .or. i .eq. 19 .or. i .eq. 22 .or. i .eq. 23 &
  .or. i .eq. 26 .or. i .eq. 36) qpower (i) = 1.d0
end do
amu (1, 1) = 18.010565  
amu (1, 2) = 20.014811  
amu (1, 3) = 19.014780
amu (1, 4) = 19.016740  
amu (1, 5) = 21.020985  
amu (1, 6) = 20.020956  
amu (2, 1) = 43.989830  
amu (2, 2) = 44.993185  
amu (2, 3) = 45.994076  
amu (2, 4) = 44.994045  
amu (2, 5) = 46.997431  
amu (2, 6) = 45.997400  
amu (2, 7) = 47.998322  
amu (2, 8) = 46.998291  
amu (3, 1) = 47.984745  
amu (3, 2) = 49.988991  
amu (3, 3) = 49.988991  
amu (3, 4) = 48.988960  
amu (3, 5) = 48.988960  
amu (4, 1) = 44.001062  
amu (4, 2) = 44.998096  
amu (4, 3) = 44.998096  
amu (4, 4) = 46.005308  
amu (4, 5) = 45.005278  
amu (5, 1) = 27.994915  
amu (5, 2) = 28.998270  
amu (5, 3) = 29.999161  
amu (5, 4) = 28.999130  
amu (5, 5) = 31.002516  
amu (5, 6) = 30.002485  
amu (6, 1) = 16.031300  
amu (6, 2) = 17.034655  
amu (6, 3) = 17.037475  
amu (7, 1) = 31.989830  
amu (7, 2) = 33.994076  
amu (7, 3) = 32.994045  
amu (8, 1) = 29.997989  
amu (8, 2) = 30.995023  
amu (8, 3) = 32.002234  
amu (9, 1) = 63.961901  
amu (9, 2) = 65.957695  
amu (10, 1) = 45.992904  
amu (11, 1) = 17.026549  
amu (11, 2) = 18.023583  
amu (12, 1) = 62.995644  
amu (13, 1) = 17.002740  
amu (13, 2) = 19.006986  
amu (13, 3) = 18.008915  
amu (14, 1) = 20.006229  
amu (15, 1) = 35.976678  
amu (15, 2) = 37.973729  
amu (16, 1) = 79.926160  
amu (16, 2) = 81.924115  
amu (17, 1) = 127.912297  
amu (18, 1) = 50.963768  
amu (18, 2) = 52.960819  
amu (19, 1) = 59.966986  
amu (19, 2) = 61.962780  
amu (19, 3) = 60.970341  
amu (19, 4) = 60.966371  
amu (19, 5) = 61.971231  
amu (20, 1) = 30.010565  
amu (20, 2) = 31.013920  
amu (20, 3) = 32.014811  
amu (21, 1) = 51.971593  
amu (21, 2) = 53.968644  
amu (22, 1) = 28.006147  
amu (23, 1) = 27.010899  
amu (23, 2) = 28.014254  
amu (23, 3) = 28.007933  
amu (24, 1) = 49.992328  
amu (24, 2) = 51.989379  
amu (25, 1) = 34.005480  
amu (26, 1) = 26.015650  
amu (26, 2) = 27.019005  
amu (27, 1) = 30.046950  
amu (28, 1) = 33.997238  
amu (29, 1) = 65.991722  
amu (30, 1) = 145.962492
amu (31, 1) = 33.987721
amu (31, 2) = 35.983515
amu (31, 3) = 34.987105
amu (32, 1) = 46.005480
amu (33, 1) = 32.997655
amu (34, 1) = 15.994915
amu (35, 1) = 96.956672
amu (35, 2) = 98.953723
amu (26, 1) = 29.997989
amu (37, 1) = 95.921076
amu (37, 2) = 97.919027
amu (38, 1) = 28.031300
amu (38, 2) = 29.034655
amu (39, 1) = 32.026215

RETURN
END SUBROUTINE hitran_setup

!
SUBROUTINE voigt (x, a, v, ndim, nvlo, nvhi)

! the following calculated voigt values at all grid values for each point
! subroutine voigt (x, a, v, nx)
! voigt first and second derivatives commented out
! subroutine voigt (x, a, v, dv, d2v, nx)

IMPLICIT NONE
INTEGER, INTENT(IN)      :: ndim, nvlo, nvhi
REAL(KIND=8), INTENT(IN) :: a
REAL(KIND=8), DIMENSION(ndim), INTENT(IN)  :: x
REAL(KIND=8), DIMENSION(ndim), INTENT(OUT) :: v ! , dv, d2v 

REAL (KIND=8), PARAMETER  :: sqrtpi =1.77245385090551d0, &
     twooverpi = 0.63661977236758d0,  fouroverpi = 1.27323954473516d0

INTEGER      :: capn, i, nu, n, in, np1
REAL(KIND=8) :: lamda, sfac, absx, s, h, h2, r1, r2, s1, s2, t1, t2, c !, c2v
LOGICAL      :: b

! a = 0.
IF (a < 1.0d-8) THEN
   v(nvlo:nvhi) = dexp(-x(nvlo:nvhi)**2) / sqrtpi
   !dv(nvlo:nvhi) = -2.0d0 * x(nvlo:nvhi) * v(nvlo:nvhi)
   !d2v(nvlo:nvhi) = (4.0d0 * x(nvlo:nvhi) ** 2 - 2.d0) * v(nvlo:nvhi)

   ! add lorentzian check here, for speed
ELSE
  ! coefficient for second derivative
  ! c2v = 4.0d0 * a * a + 2.d0

  sfac = 1.0d0 - a / 4.29d0
  DO i = nvlo, nvhi
  ! do i = 1, nx
     absx = dabs (x (i))
     IF ((a < 4.29d0) .AND. (absx < 5.33d0)) THEN
        s = sfac * dsqrt (1.d0 - (x (i) / 5.33d0)**2)
        h = 1.6d0 * s
        h2 = 2.0d0 * h
        capn = 6.d0 + 23.0d0 * s
        lamda = h2**capn
        nu = 9.0d0 + 21.0d0 * s
     ELSE
        h = 0.0d0
        capn = 0
        nu = 8
     ENDIF
     b = (h == 0.0d0) .or. (lamda == 0.0d0)
     r1 = 0.0d0
     r2 = 0.0d0
     s1 = 0.0d0
     s2 = 0.0d0
     n = nu
     DO in = 1, nu + 1
        np1 = n + 1
        t1 = a + h + dfloat (np1) * r1
        t2 = absx - dfloat (np1) * r2
        c = .5d0 / (t1 * t1 + t2 * t2)
        r1 = c * t1
        r2 = c * t2
        IF ((h > 0.0d0) .AND. (n <= capn)) THEN
           t1 = lamda + s1
           s1 = r1 * t1 - r2 * s2
           s2 = r2 * t1 + r1 * s2
           lamda = lamda / h2
        ENDIF
        n = n - 1
     ENDDO
     IF (b) THEN
        v (i) = twooverpi * r1
        !     dv (i) = fouroverpi * (a * r2 - absx * r1)
     ELSE
        v (i) = twooverpi * s1
        !     dv (i) = fouroverpi * (a * s2 - absx * s1)
     ENDIF
     !   dv (i) = -dsign (dv (i), x (i))
     !   d2v (i) = fouroverpi * a - (c2v + 4.d0 * x (i) * x (i)) * &
     !   v (i) - 4.d0 * x (i) * dv (i)
  ENDDO
ENDIF

RETURN
END SUBROUTINE voigt

!
FUNCTION ibin (vtarget, array, nentries) RESULT(idx)

! binary search in an array of real numbers in increasing order.
! returned is the number of the last entry which is less than target, or
! 0 if not within array. (this was written to find values enclosing
! target for a linear interpolation scheme.) 4/9/84 john lavagnino;
! adapted from jon bentley, cacm february 1984, vol. 27, no. 2, p. 94.

IMPLICIT NONE
INTEGER, INTENT(IN)      :: nentries
REAL(KIND=8), INTENT(IN) :: vtarget
REAL(KIND=8), DIMENSION(nentries), INTENT(IN) :: array

INTEGER :: upper, lower, middle
INTEGER :: idx

lower = 0
upper = nentries + 1

DO WHILE (lower + 1 /= upper)
  middle = (lower + upper) / 2
  IF (array(middle) < vtarget) THEN
     lower = middle
  ELSE
     upper = middle
  ENDIF
ENDDO

! at this point, either array (lower) <= target <= array (upper), or
! lower = 0, or upper = nentries + 1 (initial values).
IF (lower > 0 .and. upper /= nentries + 1) THEN
   idx = lower
ELSE
   idx = 0
ENDIF

END FUNCTION ibin

!
SUBROUTINE gauss (pos, spec, specmod, npoints, hw1e)

! convolves input spectrum with gaussian slit function of specified hw1e
! (half-width at 1/e intensity). Assumes input spectrum has constant
! spacing (will need to modify or replace with nm version when and if
! needed).


INTEGER,                            INTENT (IN)  :: npoints
REAL (KIND=8),                      INTENT (IN)  :: hw1e
REAL (KIND=8), DIMENSION (npoints), INTENT (IN)  :: pos, spec
REAL (KIND=8), DIMENSION (npoints), INTENT (OUT) :: specmod

INTEGER                            :: nhi, nlo, i, j, nslit
REAL (KIND=8)                      :: emult, delwvl, slitsum, slit0
REAL (KIND=8), DIMENSION (npoints) :: slit

specmod = spec
IF ( hw1e == 0.0 ) RETURN

emult = - 1.d0 / (hw1e**2)
delpos = pos (2) - pos (1)

! calculate slit function values out to 0.001 times x0 value, normalize
! so that sum = 1.
slitsum = 1.0 ;  slit0 = 1.0
i = 1  ;  nslit = 0
DO WHILE ( nslit <= npoints )
   slit (i) = EXP (emult * (delpos * i)**2)
   slitsum = slitsum + 2.0 * slit(i)
   IF (slit (i) <= 0.001 ) EXIT
   i = i + 1
ENDDO
nslit = i

slit0 = slit0 / slitsum
slit(1:nslit) = slit(1:nslit) / slitsum

! convolve spectrum. don't reflect at endpoints (for now).
specmod = slit0 * spec
DO i = 1, npoints
   DO j = 1, nslit
      nlo = i - j
      nhi = i + j

      ! Refect at endpoints
      !IF (nlo < 1) nlo = -nlo + 2
      !IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)

      IF (nlo >= 1) specmod (i) = specmod (i) + slit (j) * spec (nlo)
      IF (nhi <= npoints) specmod (i) = specmod (i) + slit (j) * spec (nhi)
   ENDDO
ENDDO

RETURN
END SUBROUTINE gauss

!------------------------------------------------------------------------------
!S+
! NAME:
!       StrUpCase
!
! PURPOSE:
!       Function to convert an input string to upper case.
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       Result = StrUpCase( String )
!
! INPUT ARGUMENTS:
!       String:  Character string to be converted to upper case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( * )
!                DIMENSION:  Scalar
!                ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result:  The input character string converted to upper case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( LEN(String) )
!                DIMENSION:  Scalar
!
! CALLS:
!       None.
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! EXAMPLE:
!       string = 'this is a string'
!       WRITE( *, '( a )' ) StrUpCase( string )
!   THIS IS A STRING
!
! PROCEDURE:
!       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!       1995 Springer-Verlag, New York.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String )
  
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  
  ! -- Local variables
  INTEGER :: i, n
  
  ! -- Copy input string
  Output_String = Input_String
  
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
     
     ! -- Find location of letter in lower case constant string
     n = INDEX( LOWER_CASE, Output_String( i:i ) )
     
     ! -- If current substring is a lower case letter, make it upper case
     IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n )
     
  END DO
  
END FUNCTION StrUpCase


!------------------------------------------------------------------------------
!S+
! NAME:
!       StrLowCase
!
! PURPOSE:
!       Function to convert an input string to lower case.
!
! CATEGORY:
!       Utility
!
! LANGUAGE:
!       Fortran-95
!
! CALLING SEQUENCE:
!       Result = StrLowCase( String )
!
! INPUT ARGUMENTS:
!       String: Character string to be converted to lower case.
!               UNITS:      N/A
!               TYPE:       CHARACTER( * )
!               DIMENSION:  Scalar
!               ATTRIBUTES: INTENT( IN )
!
! OPTIONAL INPUT ARGUMENTS:
!       None.
!
! OUTPUT ARGUMENTS:
!       None.
!
! OPTIONAL OUTPUT ARGUMENTS:
!       None.
!
! FUNCTION RESULT:
!       Result:  The input character string converted to lower case.
!                UNITS:      N/A
!                TYPE:       CHARACTER( LEN(String) )
!                DIMENSION:  Scalar
!
! CALLS:
!       None.
!
! SIDE EFFECTS:
!       None.
!
! RESTRICTIONS:
!       None.
!
! EXAMPLE:
!       string = 'THIS IS A STRING'
!       WRITE( *, '( a )' ) StrLowCase( string )
!   this is a string
!
! PROCEDURE:
!       Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
!       1995 Springer-Verlag, New York.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 18-Oct-1999
!                       paul.vandelst@ssec.wisc.edu
!S-
!------------------------------------------------------------------------------

FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String )
  
  ! -- Argument and result
  CHARACTER( * ), INTENT( IN )     :: Input_String
  CHARACTER( LEN( Input_String ) ) :: Output_String
  
  CHARACTER( * ), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
  CHARACTER( * ), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
  
  ! -- Local variables
  INTEGER :: i, n
  
  
  ! -- Copy input string
  Output_String = Input_String
  
  ! -- Loop over string elements
  DO i = 1, LEN( Output_String )
     
     ! -- Find location of letter in upper case constant string
     n = INDEX( UPPER_CASE, Output_String( i:i ) )
     
     ! -- If current substring is an upper case letter, make it lower case
     IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n )
     
  END DO
  
END FUNCTION StrLowCase

