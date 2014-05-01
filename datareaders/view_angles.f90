!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    view_angles
!
! PURPOSE
!    Calculate viewing angles of satellite at a given target location P, 
!    subsatellite point SSP, and satellite altitude h.
!
! INPUT
!    latp  = latitude of target in degrees (-90 to 90)
!    lonp  = longitude of target in degrees (-180 to 180)
!    latss = latitude of subsatellite point in degrees (-90 to 90)
!    lonss = longitude of subsatellite point in degrees (-180 to 180)
!    h     = altitude above Earth in km
!
! OUTPUT
!    view_angles is returned array of length 2.
!    view_angles(1) is Satellite Azimuth Angle in Degrees, 0-360, clockwise from North at Target.
!    view_angles(2) is Satellite Elevation Angle in Degrees, 0-90, measured at target between
!    satellite and local horizontal (Angle EPSILON in Reference).
!    Example in Reference: Given LATP,LONP = 22, -160; LATSSP,LONSSP = 10, -175; h=1000 km;
!    view_angles = [232.5, 14.42] NOTE: We need view angles for satellite as viewed from target, 
!    hence results are for this scenario and not that in the reference.
!    
! REFERENCE
!    Space Mission Analysis and Design, 
!        by James R. Wertz, Wiley J. Larson. 
!        Chapter 5, Space Mission Geometry, pp 112-114. 
!        (These pages are available online in
!         http://astrobooks.com/files/SMAD3Err3rd.pdf)
!    See also http://www.aoe.vt.edu/~cdhall/courses/aoe4140/missa.pdf
!    and http://en.wikipedia.org/wiki/Spherical_law_of_cosines
!
! HISTORY
!    written 1 April 2010 by Joyce Wolf for Annemarie Eldering and Susan Kulawik (JPL)
!    updated 30 April 2010 to add warning if target is outside 
!    region visible to spacecraft (spacecraft is below horizon at target).
!    translated into Fortran90 24 Feb 2012 by Vijay Natraj
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

FUNCTION view_angles(latp, lonp, latss, lonss, h)

implicit none

!  Inputs

real(kind=8),    intent(in)  :: latp
real(kind=8),    intent(in)  :: lonp
real(kind=8),    intent(in)  :: latss
real(kind=8),    intent(in)  :: lonss
real(kind=8),    intent(in)  :: h

!  Outputs
   
real(kind=8)                 :: view_angles(2)

!  Local variables

real(kind=8)                 :: d2r, re, srho
real(kind=8)                 :: deltaL, cdeltaL, clatss, slatss, clatp, slatp
real(kind=8)                 :: clambda, slambda, cphiv, phiv
real(kind=8)                 :: taneta, eta, ceps, eps, r2d
real(kind=8)                 :: lambda, lambda0

d2r = acos(-1.d0)/180.d0

! Earth radius in km

re = 6378.d0

! sin(angular radius of Earth at height h)

srho = re/(re+h)

deltaL = abs(lonss-lonp)
cdeltaL = cos(d2r*deltaL)
clatss = cos(d2r*latss)
slatss = sin(d2r*latss)
clatp = cos(d2r*latp)
slatp = sin(d2r*latp)

! Find lambda, central angle of great circle arc connecting P and SSP.
! use Law of Cosines for spherical triangle with vertices P, SSP, and North Pole.
! sides are central angles (great circle arcs) 90-LatP, 90-LatSS, and Lambda.
! Law of Cosines: cos(c) = cos(a) cos(b) +sin(a) sin(b) cos(C),
! where a, b, c are the sides and C is the corner angle opposite side c.

! cos(lambda)

clambda = slatss*slatp + clatss*clatp*cdeltaL

! sin(lambda) 

slambda = sin(acos(clambda))

! cos phiv (Phi_V is azimuth of satellite measured from North at target P).
! Use Law of Cosines on Spherical Triangle formed by P, North Pole, SSP.

cphiv = (slatss - slatp*clambda) / ( clatp*slambda)
if (cphiv .GT. 1.d0) cphiv = 1.d0
if (cphiv .LT. -1.d0) cphiv = -1.d0
phiv = acos(cphiv)

! tan eta

taneta = srho*slambda / (1.d0-srho*clambda)
eta = atan(taneta)

! cos epsilon

ceps = sin(eta)/srho
eps = acos(ceps)

r2d = 180.d0/acos(-1.d0)
view_angles(1) = r2d*phiv
if (lonp-lonss .GT. 0.d0) view_angles(1) = 360.d0 - view_angles(1)
view_angles(2) = r2d*eps

! Check for spacecraft below horizon at target

lambda = r2d*acos(clambda)
lambda0 = r2d*acos(srho)
if (lambda .GT. lambda0) then
   write(0,*) 'WARNING: SPACECRAFT BELOW HORIZON AT TARGET'
   write(0,*) 'Lambda  (Central Angle between Target and SSP) = ', lambda
   write(0,*) 'Lambda0 (Central Angle Visible to Spacecraft)  = ', lambda0
   view_angles(2) = -view_angles(2)
endif

return
end function view_angles
