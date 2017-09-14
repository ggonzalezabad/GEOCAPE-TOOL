! This tool is taken from Rob's Combo-VLIDORT tool, used to generate aerosol
! or trace gas profiles/plumes using generalized distribution function (GDF) 
subroutine generate_plume &
     ( max_layers, nlayers, level_heights, & ! input
     profile_type, z_upperlimit, z_lowerlimit,  & ! input
     z_peakheight, total_column, half_width,    & ! input
     relaxation,                                & ! input
     profile, d_profile_dcol,                   & ! output
     d_profile_dpkh, d_profile_dhfw,            & ! output
     d_profile_drel, fail, message)              ! output
  
! ------
! Inputs
! ------

! Dimensioning
  INTEGER, INTENT(IN) :: max_layers, nlayers

! H level input
  REAL(KIND=8), DIMENSION(0:MAX_LAYERS), INTENT(IN) :: level_heights   

! Profile shape type
  CHARACTER(LEN=3) :: profile_type

! Height levels controlling profile shape. All in [km].

! z_upperlimit = uppermost height value at start  of GDF profile
! z_lowerlimit = lowest    height value at bottom of GDF profile
! z_peakheight = height at which peak value of GDF profile occurs

! Notes:
!   (1) Must have  z_u > z_p >/= z_l. This is checked internally
!   (2) Must have  z_l > heightgrid(nlayers). Also checked internally
  REAL(KIND=8), INTENT(IN) :: z_upperlimit, z_lowerlimit,  z_peakheight

! Half width
  REAL(KIND=8), INTENT(IN) :: half_width

! Total column between upper and lower limits
  REAL(KIND=8), INTENT(IN) :: total_column

! Exponential relaxation
  REAL(KIND=8), INTENT(IN) :: relaxation

! -------
! Outputs
! -------
  REAL(KIND=8), DIMENSION(MAX_LAYERS), INTENT(OUT)  :: profile
! Umkehr profile derivatives (3 of them)
!   Ashplume:  (w.r.t. total tau and peak height and half width)
  REAL(KIND=8), DIMENSION(MAX_LAYERS), INTENT(OUT) :: d_profile_dcol
  REAL(KIND=8), DIMENSION(MAX_LAYERS), INTENT(OUT) :: d_profile_dpkh
  REAL(KIND=8), DIMENSION(MAX_LAYERS), INTENT(OUT) :: d_profile_dhfw
  REAL(KIND=8), DIMENSION(MAX_LAYERS), INTENT(OUT) :: d_profile_drel

! Error handling
  LOGICAL, INTENT(OUT)         :: FAIL
  CHARACTER*(*), INTENT(INOUT) :: MESSAGE
  CHARACTER(LEN=256)           :: ACTION

! ---------------
! Local variables
! ---------------

! Control for calculating derivatives
  logical :: do_derivatives

! Other variables
  integer :: i

! Exception handling. Check physical limits of plume
  fail = .false.
  message = ' '
  action  = ' '
  if ( level_heights(nlayers) .gt. z_lowerlimit ) then
     fail = .true.
     message = 'adjust Aerosol Plume to be all above ground level'
     action  = 'Set z_lowerlimit > level_heights(nlayers)'
     return
  endif
  
! Initialize
  do i = 1, nlayers
     profile(i)        = 0.0d0
     d_profile_dcol(i) = 0.0d0
     d_profile_dpkh(i) = 0.0d0
     d_profile_dhfw(i) = 0.0d0
     d_profile_drel(i) = 0.0d0
  enddo
  
! Call
  do_derivatives = .true.
  call profiles_gdfone  & 
       ( max_layers, nlayers, level_heights, do_derivatives, & ! Input
       z_upperlimit, z_peakheight, z_lowerlimit,             & ! Input
       half_width, total_column,                             & ! Input
       profile, d_profile_dcol,                              & ! Output
       d_profile_dpkh, d_profile_dhfw,                       & ! Output
       fail, message, action )                                 ! Output
              
  RETURN
END SUBROUTINE generate_plume

! 

SUBROUTINE profiles_gdfone                               & 
     ( maxlayers, nlayers, heightgrid, do_derivatives,   & ! Input
     z_upperlimit, z_peakheight, z_lowerlimit,           & ! Input
     half_width, total_column,                           & ! Input
     profile, d_profile_dcol,                            & ! Output
     d_profile_dpkh, d_profile_dhfw,                     & ! Output
     fail, message, action ) 

!  Inputs
!  ======

!  Height grid in [km]

  INTEGER, INTENT(in)                               :: maxlayers, nlayers
  REAL (kind=8), DIMENSION(0:maxlayers), INTENT(IN) :: heightgrid
  
!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of GDF profile
!     z_lowerlimit = lowest    height value at bottom of GDF profile
!     z_peakheight = height at which peak value of GDF profile occurs

!  Notes:
!    (1) Must have  z_u > z_p >/= z_l. This is checked internally
!    (2) Must have  z_l > heightgrid(nlayers). Also checked internally

  REAL(KIND=8), INTENT(IN) :: z_upperlimit, z_lowerlimit,  z_peakheight
  
  !  Half width
  
  REAL(KIND=8), INTENT(IN) :: half_width
  
!  Total column between upper and lower limits
  
  REAL(KIND=8), INTENT(IN) :: total_column

!  Control for calculating derivatives w.r.t. quantities

  LOGICAL, INTENT(IN)      :: do_derivatives

!  Output
!  ======

!  Umkehr profile 

  REAL(KIND=8), DIMENSION(MAXLAYERS), INTENT(OUT) :: profile

!  Umkehr profile derivatives (3 of them)
!    Ashplume:  (w.r.t. total tau and peak height and half width)

  REAL(KIND=8), DIMENSION(MAXLAYERS), INTENT(OUT):: d_profile_dcol
  REAL(KIND=8), DIMENSION(MAXLAYERS), INTENT(OUT):: d_profile_dpkh
  REAL(KIND=8), DIMENSION(MAXLAYERS), INTENT(OUT):: d_profile_dhfw

!  Exception handling

  LOGICAL, INTENT(OUT)         :: fail
  CHARACTER*(*), INTENT(INOUT) :: message, action
  
!  Local variables
!  ===============

!  Local commons

  REAL(KIND=8) :: omega,z1,z0,z2,h
  common /constraints /omega,z1,z0,z2,h
  
!  help variables
     
  INTEGER      :: n, nupper, nlower
  REAL(KIND=8) :: a,zt1,zt2,z,q,r,za,zb
  REAL(KIND=8) :: a_d,a_p,a_w,term,q1,q2,r1,r2,rsq1,rsq2
  
!  Debug

  LOGICAL :: do_debug_profile
  INTEGER :: m,msav
  INTEGER, PARAMETER :: mdebug = 1001
  REAL(KIND=8) :: sum, sumd
  REAL(KIND=8), DIMENSION(mdebug) :: pp,zp
  
!  Initial checks
!  ==============

!  set debug output

!      do_debug_profile = .true.
  do_debug_profile = .false.

!  Initialise exception handling

  fail = .false.
  message = ' '
  action  = ' '

!  Check input

  if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
     fail = .true.
     message = 'Bad input: GDF Lower limit is below ground'
     action  = ' increase z_lowerlimit > heightgrid(nlayers)'
     return
  endif
  
  if ( z_lowerlimit .gt. z_peakheight ) then
     fail = .true.
     message = 'Bad input: GDF Lower limit is above peak height'
     action  = ' increase z_peakheight > z_lowerlimit'
     return
  endif
  
  if ( z_peakheight .gt. z_upperlimit ) then
     fail = .true.
     message = 'Bad input: GDF transition height above upper limit'
     action  = ' increase z_upperlimit > z_peakheight'
     return
  endif
  
!  Find intercept layers
!  ---------------------

  n = 0
  do while ( heightgrid(n).ge.z_upperlimit  )
     n = n + 1
  enddo
  nupper = n
  
  do while (heightgrid(n).gt.z_lowerlimit)
     n = n + 1
  enddo
  nlower = n
  
  !  Zero output between the intercept layers

  do n = nupper, nlower
     profile(n) = 0.0d0
  enddo
  if ( do_derivatives ) then
     do n = nupper, nlower
        d_profile_dcol(n) = 0.0d0
        d_profile_dpkh(n) = 0.0d0
        d_profile_dhfw(n) = 0.0d0
     enddo
  endif
  
!  Basic settings for local commons
  !  --------------------------------
  
  h  = half_width
  z1 = z_upperlimit
  z0 = z_peakheight
  z2 = z_lowerlimit
  omega = total_column
  
!  set profile parameters

  zt1 = z1 - z0
  zt2 = z0 - z2
  q1 = dexp ( - h * zt1)
  q2 = dexp ( - h * zt2)
  r1 = 1.0d0 / ( 1.0d0 + q1 )
  r2 = 1.0d0 / ( 1.0d0 + q2 )
  term = ( r1 + r2 - 1.0d0 )
  a  = omega * h / term
  
!  Set profile derivative parameters

  rsq1 = r1 * r1 * q1
  rsq2 = r2 * r2 * q2
  a_d = a / omega
  a_p  = - a * h * ( rsq2 - rsq1 ) / term
  a_w = ( omega - a * ( rsq2*zt2 + rsq1*zt1 ) ) / term
  
!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit. Partly EXP

  n = nupper
  za = z1
  zb = heightgrid(nupper)
  if ( do_derivatives ) then
     call set_profile_plus (za,zb,a,a_d,a_p,a_w, &
             profile(n),     d_profile_dcol(n),  &
             d_profile_dpkh(n),d_profile_dhfw(n))
  else
     call set_profile_only(za,zb,a,profile(n))
  endif
  
  !  Intermediate Layers

  do n = nupper+1,nlower-1
     za = heightgrid(n-1)
     zb = heightgrid(n)
     if ( do_derivatives ) then
        call set_profile_plus (za,zb,a,a_d,a_p,a_w,  & 
                profile(n),     d_profile_dcol(n),   &
                d_profile_dpkh(n),d_profile_dhfw(n))
     else
        call set_profile_only(za,zb,a,profile(n))
     endif
  enddo
  
  !  Layer containing z_lowerlimit
  !   ( Will not be done if nupper = nlower)
  
  if ( nupper .lt. nlower ) then
     n = nlower
     za = heightgrid(n-1)
     zb = z2
     if ( do_derivatives ) then
        call set_profile_plus (za,zb,a,a_d,a_p,a_w, &
                profile(n),     d_profile_dcol(n), &
              d_profile_dpkh(n),d_profile_dhfw(n))
     else
        call set_profile_only(za,zb,a,profile(n))
     endif
  endif
  
  !  Check integrated profile
  !  ------------------------
  
  sum = 0.0d0
  sumd = 0.0d0
  do n = nupper, nlower
     sum = sum + profile(n)
     if(do_derivatives)sumd = sumd + d_profile_dcol(n)
  enddo
  
  !  Debug profile itself
  !  --------------------
  
  if ( do_debug_profile ) then
     m = 0
     do z = z1, z2, -0.01
        m = m + 1
        if(m.gt.mdebug)stop'debug dimension needs to be increased'
        zp(m) = z
        q = dexp ( - h * dabs(( z - z0 )) )
        r = 1.0d0 + q
        pp(m) = a * q / r / r
     enddo
     msav = m
     do m = 1, msav
        write(81,'(13f10.4)')zp(m),pp(m)
     enddo
  endif
  
  !  Finish
  
  return
end SUBROUTINE profiles_gdfone

subroutine set_profile_only(za,zb,a,prof)
  
  !  Profile assignation
  
  !  Arguments
  
  REAL(KIND=8), INTENT(IN)  :: za,zb,a
  REAL(KIND=8), INTENT(OUT) :: prof
  
  REAL(KIND=8) :: omega,z1,z0,z2,h
  common /constraints /omega,z1,z0,z2,h
  
  !  Local variables
  REAL(KIND=8) :: fac,ra,rb
  
  !  initialise
  
  prof = 0.0d0
  
  fac = a / h
  ra = 1.0d0 / ( 1.0d0 + dexp ( - h * dabs(za-z0)) ) 
  rb = 1.0d0 / ( 1.0d0 + dexp ( - h * dabs(zb-z0)) )
  if ( za.gt.z0 ) then
     if ( zb .ge. z0 ) then
        prof = fac * ( ra - rb )
     else 
        prof = fac * ( ra + rb - 1.0d0 )
     endif
  else
     prof = fac * ( rb - ra )
  endif
  
  !  Finish
  
  return
end subroutine set_profile_only


subroutine set_profile_plus ( za,zb,a,a_d,a_p,a_w, &
     prof,deriv1,deriv2,deriv3)
  
!  Profile and derivative assignation

!  Arguments

  REAL(KIND=8), INTENT(IN)  :: za,zb,a,a_d,a_p,a_w
  REAL(KIND=8), INTENT(OUT) :: prof,deriv1,deriv2,deriv3

  REAL(KIND=8) :: omega,z1,z0,z2,h
  common /constraints /omega,z1,z0,z2,h

!  Local variables

  REAL(KIND=8) :: fac,qa,qb,ra,rb,mza,mzb,rsqa,rsqb,rab,fac_d
  REAL(KIND=8) :: fac_p,ra_p,rb_p,rab_p
  REAL(KIND=8) :: fac_w,ra_w,rb_w,rab_w

!  initialise

  prof   = 0.0d0
  deriv1 = 0.0d0
  deriv2 = 0.0d0
  deriv3 = 0.0d0
  
  !  Type 1

  fac = a / h
  mza =  dabs(za-z0)
  mzb =  dabs(zb-z0)
  qa = dexp ( - h * mza)
  qb = dexp ( - h * mzb)
  ra = 1.0d0 / ( 1.0d0 + qa ) 
  rb = 1.0d0 / ( 1.0d0 + qb )
  fac_d = a_d/h
  fac_p = a_p/h
  fac_w = (a_w/h) - (fac/h)
  rsqa =  qa * ra * ra
  rsqb =  qb * rb * rb
  ra_p = rsqa * h
  rb_p = rsqb * h
  ra_w = rsqa * mza
  rb_w = rsqb * mzb
  if ( za.gt.z0 ) then
     if ( zb .ge. z0 ) then
        rab   =   ra   - rb
        rab_p = - ra_p + rb_p
        rab_w =   ra_w - rb_w
     else 
        rab   =   ra   + rb   - 1.0d0
        rab_p = - ra_p + rb_p
        rab_w =   ra_w + rb_w
     endif
  else
     rab    = rb   - ra
     rab_p  = rb_p - ra_p
     rab_w  = rb_w - ra_w
  endif
  prof   = fac   * rab
  deriv1 = fac_d * rab
  deriv2 = fac_p * rab + fac * rab_p
  deriv3 = fac_w * rab + fac * rab_w
  
!  Finish

  return
end subroutine set_profile_plus
