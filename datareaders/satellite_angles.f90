!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    satellite_angles
!
! PURPOSE
!    For geostationary satellite GEOCAPE, at a given target location,
!    calculate spacecraft Azimuth Viewing Angle VAA
!    and spacecraft Zenith Viewing Angle VZA;
!    also given UTC time calculate
!    Solar Zenith Angle (SZA) and Solar Azimuth Angle(SAA).
!    (User may then calculate Relative Azimuth Angle RAA = VAA-SAA)
!
! INPUTS
!    Yr,Mo,Day,Hr,Min,Sec = UTC Time
!    tz    = local time zone (0.d0 = UTC, -8.d0 = PST)
!    Latp  = latitude of target in degrees (-90 to 90)
!    Lonp  = longitude of target in degrees (-180 to 180)
!    Latss = latitude of satellite in degrees (-90 to 90)
!    Lonss = longitude of satellite in degrees (-180 to 180)
!    hgtss = altitude of satellite in km
!
! OUTPUTS
!   satellite_angles    = returned array of length 4
!   satellite_angles(1) = Spacecraft Viewing Azimuth Angle in deg (0-360)
!   satellite_angles(2) = Spacecraft Viewing Zenith Angle in deg (0-90)
!   satellite_angles(3) = Solar Azimuth Angle in degrees (0-360)
!   satellite_angles(4) = Solar Zenith Angles in degrees (0-180, >90 is below horizon)
!
! REQUIRED
!   solar_angles.f90
!   view_angles.f90
!   julday.f90
!
! HISTORY
!   written 29 April 2010 by Joyce Wolf for Annmarie Eldering & S. Kulawik
!   translated into Fortran90 24 Feb 2012 by Vijay Natraj
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION satellite_angles(yr,mo,day,hr,min,sec,tz,latp,lonp,latss,lonss,hgtss)

implicit none

!  Inputs

integer(kind=4), intent(in) :: yr
integer(kind=4), intent(in) :: mo
integer(kind=4), intent(in) :: day
integer(kind=4), intent(in) :: hr
integer(kind=4), intent(in) :: min
real(kind=8),    intent(in) :: sec
real(kind=8),    intent(in) :: tz
real(kind=8),    intent(in) :: latp
real(kind=8),    intent(in) :: lonp
real(kind=8),    intent(in) :: latss
real(kind=8),    intent(in) :: lonss
real(kind=8),    intent(in) :: hgtss

!  Outputs

real(kind=8)                :: satellite_angles(4)

!  Local variables

real(kind=8)                :: vangles(2), phiv, thetav
real(kind=8)                :: azel(2), saa, sza

interface
function view_angles(latp,lonp,latss,lonss,hgtss)
real(kind=8),    intent(in)  :: latp
real(kind=8),    intent(in)  :: lonp
real(kind=8),    intent(in)  :: latss
real(kind=8),    intent(in)  :: lonss
real(kind=8),    intent(in)  :: hgtss
real(kind=8)                 :: view_angles(2)
end function view_angles
end interface

interface
function solar_angles(yr,mo,day,hr,min,sec,tz,lat,lon) 
integer(kind=4), intent(in)  :: yr
integer(kind=4), intent(in)  :: mo
integer(kind=4), intent(in)  :: day
integer(kind=4), intent(in)  :: hr
integer(kind=4), intent(in)  :: min
real(kind=8),    intent(in)  :: sec
real(kind=8),    intent(in)  :: tz
real(kind=8),    intent(in)  :: lat
real(kind=8),    intent(in)  :: lon
real(kind=8)                 :: solar_angles(2)
end function solar_angles
end interface

vangles = view_angles(latp,lonp,latss,lonss,hgtss)
phiv    = vangles(1)
thetav  = 90.d0 - vangles(2)

azel = solar_angles(yr,mo,day,hr,min,sec,tz,latp,lonp)
saa = azel(1)
sza = 90.d0 - azel(2)

satellite_angles(1) = phiv
satellite_angles(2) = thetav
satellite_angles(3) = saa
satellite_angles(4) = sza

return
end function satellite_angles
