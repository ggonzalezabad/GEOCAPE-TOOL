!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    solar_angles
!
! PURPOSE
!   Calculate Elevation and Azimuth Angles of Sun at given location and time
!
! INPUT
!    yr  = 4-digit year
!    mon = month of year, 1-12
!    day = day of month
!    hr  = hour, 0-23
!    min = minute, 0-59
!    sec = second, [0.d0,60.d0)
!    tz  = local time zone (0.d0 = UTC, -8.d0 = PST)
!    lat = Latitude of location in degrees, positive North, negative South
!    lon = Longitude of location in degrees, positive East, negative West
!    
! OUTPUT
!    solar_angles is returned array of length 2.
!    solar_angles(1) is Solar Azimuth Angle in Degrees, 0-360.
!    solar_angles(2) is Solar Elevation Angle in Degrees, 0-90. Negative is Sun below Horizon.
!
! REFERENCES
!    NOAA Solar Calculator at http://www.esrl.noaa.gov/gmd/grad/solcalc/
!      (Decl, Eqtime, Az, and EL agree with values returned by above.)
!    General info about calculation details is at 
!      http://www.srrb.noaa.gov/highlights/sunrise/calcdetails.html
!         "The calculations in the NOAA Sunrise/Sunset and Solar Position Calculators
!          are based on equations from Astronomical Algorithms, by Jean Meeus."
!    Code used by NOAA Solar Calculator to get Declination and Equation of Time is at
!       http://www.srrb.noaa.gov/highlights/sunrise/program.txt
!
! HISTORY
!    written February 2010 by Joyce Wolf for Annmarie Eldering and Susan Kulawik (JPL)
!    translated into Fortran90 24 Feb 2012 by Vijay Natraj (JPL)
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

FUNCTION solar_angles(yr,mo,day,hr,min,sec,tz,lat,lon)

implicit none

!  Inputs
  
integer(kind=4), intent(in)  :: yr
integer(kind=4), intent(in)  :: mo
integer(kind=4), intent(in)  :: day
integer(kind=4), intent(in)  :: hr
integer(kind=4), intent(in)  :: min
real(kind=8),    intent(in)  :: sec
real(kind=8),    intent(in)  :: tz
real(kind=8),    intent(in)  :: lat
real(kind=8),    intent(in)  :: lon

!  Outputs
   
real(kind=8)                 :: solar_angles(2)
   
!  Local variables

real(kind=8)                 :: d2r, r2d
real(kind=8)                 :: julianday, julday
real(kind=8)                 :: t, meanlon, mean_anom
real(kind=8)                 :: mrad, sinm, sin2m, sin3m, c
real(kind=8)                 :: truelon, true_anom
real(kind=8)                 :: omega, app_lon, e, seconds, eps0
real(kind=8)                 :: epsilon, sint, decl
real(kind=8)                 :: m, lz, y, sin2lz, cos2lz, sin4lz
real(kind=8)                 :: eqtime, utc, ha
real(kind=8)                 :: declr, latr, har
real(kind=8)                 :: elr, yaz, xaz, azr, az, el

d2r = acos(-1.d0)/180.d0
r2d = 180.d0/acos(-1.d0)

! t = Julian centuries past J2000.0
! Adjust if time is local

julianday = JULDAY(mo,day,yr,hr,min,sec)
t = (julianday-2451545.d0)/36525.d0 - tz/24.d0

! Geometric Mean Longitude of Sun in degrees

meanlon = 280.46646d0 + t * (36000.76983d0 + 0.0003032d0 * t)
do while (meanlon .GT. 360.d0)
   meanlon = meanlon-360.d0
end do
do while (meanlon .LT. 0.d0)
   meanlon = meanlon + 360.d0
end do

! Geometric Mean Anomaly of Sun in degrees

mean_anom = 357.52911d0 + t * (35999.05029d0 - 0.0001537d0 * t)

! Sun Equation of Center

mrad = d2r*mean_anom
sinm = sin(mrad)
sin2m = sin(2.d0*mrad)
sin3m = sin(3.d0*mrad)
c = sinm * (1.914602d0 - t * (0.004817d0 + 0.000014d0 * t)) + &
    sin2m * (0.019993d0 - 0.000101d0 * t) + sin3m * 0.000289d0

! Sun's True Longitude in degrees

truelon = meanlon + c

! Sun's True Anomaly in degrees

true_anom = mean_anom + c

! Sun's Apparent Longitude in degrees

omega = 125.04d0 - 1934.136d0 * t
app_lon = truelon - 0.00569d0 - 0.00478d0 * sin(d2r*omega)

! Eccentricity of Earth Orbit

e = 0.016708634d0 - t * (0.000042037d0 + 0.0000001267d0 * t)

! Mean Obliquity of Ecliptic in degrees

seconds = 21.448d0 - t*(46.8150d0 + t*(0.00059d0 - t*(0.001813d0)))
eps0 = 23.d0 + (26.d0 + (seconds/60.d0)) / 60.d0

! Corrected Obliquity

epsilon = eps0 + 0.00256d0 * cos(d2r*omega)

! Sun's Declination in degrees

sint = sin(d2r*epsilon) * sin(d2r*app_lon)
decl = r2d*asin(sint)

! Equation of Time in minutes

m = mean_anom
lz = meanlon
y = (tan(d2r*epsilon/2.d0))**2
sin2lz = sin(2.d0 * d2r*lz)
sinm   = sin(d2r*m)
cos2lz = cos(2.d0 * d2r*lz)
sin4lz = sin(4.d0 * d2r*lz)
sin2m  = sin(2.d0 * d2r*m)
eqtime = y * sin2lz - 2.d0 * e * sinm &
          + 4.d0 * e * y * sinm * cos2lz &
          - 0.5d0 * y * y * sin4lz - 1.25d0 * e * e * sin2m
eqtime = r2d * eqtime * 4.d0

! Now use DECL and EQTIME to calculate Solar Elevation and Azimuth at (Lat,Lon) 
! for time specified by yr,mo,day,hr,min,sec

! UTC time in Decimal Hours

utc = real(hr,kind=8) + (real(min,kind=8)+sec/60.d0) / 60.d0 - tz

! Hour Angle in degrees, negative is before noon.

ha = 180.d0 - 15.d0*utc - lon - eqtime/4.d0

! Convert to Radians

declr = d2r*decl
latr = d2r*lat
har = d2r*ha

! Solar Elevation Angle

elr = asin(sin(latr)*sin(declr)+cos(latr)*cos(declr)*cos(har))

! Solar Azimuthal Angle
! atan2(y,x) = arctan(y,x) is in range (-pi,pi].

yaz = -sin(har)*cos(declr)
xaz = cos(har)*cos(declr)*sin(latr) - sin(declr)*cos(latr)

! Above formulas are for az measured from due South
! So rotate by 180 to to get az measured from due North

azr = atan2(-yaz,-xaz)
az = r2d* azr
if (az .LT. 0.d0) az = az + 360.d0
el = r2d*elr

solar_angles(1) = az
solar_angles(2) = el

RETURN
end FUNCTION solar_angles
