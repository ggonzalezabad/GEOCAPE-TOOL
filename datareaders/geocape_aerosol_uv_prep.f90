subroutine geocape_aerosol_setter_uv                   &
      ( GC_maxlayers, maxlambdas, profile_data,        &  ! Input
        GC_nlayers, nlambdas, lambdas,                 &  ! Input
        aer_flags, aer_opdeps, aer_ssalbs, aer_gpars,  &  ! Output
        fail, message )                                   ! Exception handling

!  inputs
!  ------

!  Dimensioning

   integer, intent(in)  :: GC_maxlayers, maxlambdas

!  Data from file-read

   real(kind=8), dimension(69,27),  intent(in) :: profile_data

!  number of layers has 'GC_' prefix (distinguish with VLIDORT variable)

   integer,  intent(in) :: GC_nlayers

!  Window of wavelengths

   integer,                                intent(in) :: nlambdas
   real(kind=8), dimension(maxlambdas),    intent(in) :: lambdas

!  Output
!  ------

!  Aerosol Stuff (Output from "geocape_aerosol_setter_1")

   logical,      dimension(GC_maxlayers), intent(INOUT)             :: aer_flags
   real(kind=8), dimension(GC_maxlayers,maxlambdas), intent(INOUT)  :: aer_opdeps
   real(kind=8), dimension(GC_maxlayers,maxlambdas), intent(INOUT)  :: aer_ssalbs
   real(kind=8), dimension(GC_maxlayers,maxlambdas), intent(INOUT)  :: aer_gpars

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  GEOCAPE profile, received from J. Joiner, 3 June 2009.
!   .... has the following 27-column attributes.....

!    1.   Pressure at bottom of layer (Pa),
!    2.   Pressure at midpoint of layer (Pa),
!    3.   Geopotential height at bottom of layer (m),
!    4.   Geopotential height at midpoint of layer (m),
!    5.   Latitude,
!    6.   Longitude,
!    7.   Temperature (K),
!    8.   Water Vapor mixing ratio (kg/kg),
!    9.   O3 (ppm),
!    10.  NO2 (ppm),
!    11.  HCHO (ppm),
!    12.  SO2 (ppm),
!    13.  CO (ppm),
!    14.  cloud optical thickness of ice,
!    15.  cloud optical thickness for water,
!    16.  aerosol optical depth at 300nm,
!    17.  aerosol optical depth at 400nm,
!    18.  aerosol optical depth at 600nm,
!    19.  aerosol optical depth at 999nm,
!    20.  asymmetry parameter at 300nm,
!    21.  asymmetry parameter at 400nm,
!    22.  asymmetry parameter at 600nm,
!    23.  asymmetry parameter at 999nm,
!    24.  single scattering albedo at 300nm,
!    25.  single scattering albedo at 400nm,
!    26.  single scattering albedo at 600nm,
!    27.  single scattering albedo at 999nm

!  Local variables
!  ---------------

!  Array of data wavelengths

   real(kind=8), dimension(4) :: aerwavs
   data aerwavs / 300.0, 400.0, 600.0, 999.0 /

!  help variables

   integer       :: n, n1, w, w1, w2, ww
   real(kind=8)  :: diff, f1, f2
   logical       :: entries

!  initialize

   fail    = .false.
   message = ' '

   aer_flags(:)    = .false.
   aer_opdeps(:,:) = 0.0d0
   aer_ssalbs(:,:) = 0.0d0
   aer_gpars(:,:)  = 0.0d0
  
!  Set the flags

   n = GC_nlayers
   entries = .true.
   do while (entries.and.n.gt.1)
     n1 = GC_nlayers - n + 1
     entries = (dabs(profile_data(n1,16)+999.00d0).gt.1.0d-06)
     if (entries) aer_flags(n) = .true.
     n = n - 1
   enddo

!  Check wavelengths are within range

   do w = 1, nlambdas
      if ( lambdas(w).lt.aerwavs(1) ) then
         message = 'Some wavelengths less than 300 nm (lower limit): Reset!'
         fail = .true.
      endif
      if ( lambdas(w).gt.aerwavs(4) ) then
         message = 'Some wavelengths greater than 999 nm (upper limit): Reset!'
         fail = .true.
      endif
   enddo
   if ( fail) return

!  Set properties (linearly interpolate with wavelength)

   do w = 1, nlambdas
      do ww = 1, 3
         if ( lambdas(w).ge.aerwavs(ww).and.lambdas(w).lt.aerwavs(ww+1) ) w1 = ww
      enddo
      w2 = w1 + 1
      diff = aerwavs(w2) - aerwavs(w1)
      f2 = (   lambdas(w) - aerwavs(w1) ) / diff
      f1 = ( - lambdas(w) + aerwavs(w2) ) / diff
      do n = 1, GC_nlayers
         n1 = GC_nlayers-n+1
         if ( aer_flags(n1) ) then
            aer_opdeps(n1,w) = f1*profile_data(n,w1+16) + f2*profile_data(n,w2+16)
            aer_ssalbs(n1,w) = f1*profile_data(n,w1+24) + f2*profile_data(n,w2+24)
            aer_gpars(n1,w)  = f1*profile_data(n,w1+20) + f2*profile_data(n,w2+20)
         endif
      enddo
   enddo

!  Finish

   return
end subroutine geocape_aerosol_setter_uv

