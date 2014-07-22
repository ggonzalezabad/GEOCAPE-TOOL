! Read generic atmospheric profiles 
! Here is the typical header of the data, lines marked with *** will be read
! The line below the End_of_Headers stores the standard ID for each column
! Here are the list of 30 recognizable standard IDs (no space inside an ID)
!    Lev: level or layer
!    PRESSURE: hPa
!    ALTITUDE: km   (not used)
!    TATM:     K
!    O3, NO2, HCHO, SO2, H2O, GLYX, BRO, OCLO, IO, CO, CO2, N2O, CH4, O2, NO, HNO3, OCS (VMR)
!    OSU, OBC, OOC, OSF, OSC, ODU, OWC, OIC, RH
! Note the input profiles must meet the following conditions:
! 1. From surface to TOA
! 2. The first 4 fields got to be Level, pressure, altitude, and TATM in order(always at level)
!    The number of fields and the order after TATM do not matter
! 3. The units for gases is VMR
! 4. The units for aerosol, cloud optical depth is dimensionless
! 5. If the data is specified for layer, Gas is for layer mean VMR between l and l + 1
! 6. Aerosol and cloud data is always specified at each layer
    
!Geos-Chem, GEOS-5, 2.5x2, Hourly data
!Data_Size = 48 x 21
!ProfileIsLevel = 0
!YYYY_MM_DD_HH = 2007 07 15 19
!Lon_Lat =   -97.50  38.00
!SZA_VZA_AZA =    17.2902   44.0700  147.4352
!SurfaceTemperature(K) =   303.3730
!SurfacePressure(hPa) =   961.8620
!IsOcean =  0
!WindSpeed(m/s) =     2.9121
!CloudFraction =     0.1065
!CloudTopPressure =   607.1039
!****** End_Of_Header ******
!Lev Pressure Altitude   TATM       O3         NO2       HCHO        SO2       GLYX        H2O        BRO       OCLO       OSU        OBC        OOC       OSF       OSC       ODU       OWC       OIC        RH    
!      hPa      km        K        VMR        VMR        VMR        VMR        VMR        VMR        VMR        VMR                                                                                            % 
subroutine geocape_profile_reader_1 &
           ( filename, maxl, maxc,                                                          &          ! input
             profile_data, footprint_data, profids, profile_is_level, GC_nlayers, pnc, message, fail ) ! output

!  inptu filename

   character(len=*),                intent(in)  :: filename

!  Main output
   integer, intent(in)                             :: maxl, maxc
   integer,                        intent(out)     :: GC_nlayers, pnc
   real(kind=8), dimension(maxl,maxc), intent(out) :: profile_data
   real(kind=8), dimension(15),     intent(out)    :: footprint_data
   logical, intent(out)                            :: profile_is_level  !T: level  F: layer
   character(len=256), intent(out)                 :: profids

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables

   integer      :: i, j, nr, nc, np, itmp
   real(kind=8), dimension(maxl,maxc+1) :: tmp_data
   character*80 :: dummy

!  initialize

   fail    = .false.
   message = ' '
   do i = 1, maxl
      do j = 1,maxc
         profile_data(i,j) = -999.
      enddo
   enddo

   np = 1
   do while (filename(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Open and read file

   open(1,file=filename(1:np),err=90, status='old')
   read(1,*)
   read(1,*) dummy,dummy,nr,dummy,nc
   if (nr .ge. maxl .or. nc .ge. maxc+1) then
      message = 'Need to increase maxl or maxc in geocape_profile_reader_1'
      fail = .true.
      return
   endif
   read(1,*) dummy, dummy, itmp
   if (itmp .eq. 0) then
      profile_is_level = .FALSE.
   else
      profile_is_level = .TRUE.
   endif
   
   read(1,*) dummy, dummy, footprint_data(1:4)  ! Year, mon, day, utc
   read(1,*) dummy, dummy, footprint_data(5:6)  ! longitude, latitude
   read(1,*) dummy, dummy, footprint_data(7:9)  ! SZA, VZA, AZA
   read(1,*) dummy, dummy, footprint_data(10)   ! Surface temperature
   read(1,*) dummy, dummy, footprint_data(11)   ! Surface pressure
   read(1,*) dummy, dummy, footprint_data(12)   ! Is ocean? (1: Ocean, for ocean BRDF)
   read(1,*) dummy, dummy, footprint_data(13)   ! Wind speed
   read(1,*) dummy, dummy, footprint_data(14)   ! Cloud fraction     
   read(1,*) dummy, dummy, footprint_data(15)   ! Cloud top pressure (for Lambertian clouds)
   read(1,*)
   read(1,'(A256)') profids
   read(1,*)

   do i = 1, nr
     read(1,*) (tmp_data(i,j),j=1,nc)
     do j = 1,3
        profile_data(i,j) = tmp_data(i,j+1)
     enddo
     do j = 4,nc-1
        profile_data(i,j) = tmp_data(i,j+1)
     enddo
   enddo

   GC_nlayers = nr-1
   pnc = nc - 1
   
   close(1)

   return

! error return

90 continue
   fail = .true.
   message = 'Open failure for profile file = '//filename(1:LEN(filename))
   return

end subroutine geocape_profile_reader_1

subroutine geocape_profile_setter_1                        &
    ( maxl, maxc, GC_nlayers, pnc, ngases, profile_data,   &  ! Input
      profids, profile_is_level, which_gases,              &  ! Input
      heights, temperatures, pressures,                    &  ! Output
      aircolumns, daircolumns_dT,                          &  ! Output
      gas_partialcolumns, gas_totalcolumns,                &  ! Output
      fail, message )                                         ! Output

!  inputs
!  ------

!  Dimensioning

   integer, intent(in) :: maxl, maxc, GC_nlayers, pnc, ngases

!  Data from file-read

   real(kind=8), dimension(maxl,maxc),  intent(in) :: profile_data
   character(len=256), intent(in)                  :: profids
   logical, intent(in)                             :: profile_is_level

!  Trace gas control

   character(Len=4), dimension ( ngases ), intent(in) :: which_gases

!  Output
!  ------

!  Atmospheric quantities (PTH)

   real(kind=8), dimension ( 0:GC_nlayers ), intent(out)  :: heights
   real(kind=8), dimension ( 0:GC_nlayers ), intent(out)  :: temperatures
   real(kind=8), dimension ( 0:GC_nlayers ), intent(out)  :: pressures

!  Air density Partial columns and T-derivative

   real(kind=8), dimension (GC_nlayers), intent(out)          :: aircolumns
   real(kind=8), dimension (GC_nlayers), intent(out)          :: daircolumns_dT

!  Trace gas partial columns (profile) and total columns

   real(kind=8), dimension (GC_nlayers,ngases), intent(out) :: gas_partialcolumns
   real(kind=8), dimension (ngases), intent(out)            :: gas_totalcolumns

!  Exception handling

   logical, intent(INOUT)       :: fail
   character*(*), intent(INOUT) :: message

!  Local variables
!  ---------------

!  Array of level temperatures

   real(kind=8), dimension(1:GC_nlayers) :: layertemp

!  Array of derived gas constants

   real(kind=8), dimension(GC_nlayers)   :: gasconstants

! Array of data fields
   character(len=4), dimension(0:maxc)   :: datanames  ! datanames(0) should not be used
   integer, dimension(ngases)            :: gasidxs
   character(len=4)                      :: gasname
   real(kind=8), dimension(GC_nlayers)   :: gasconcen

!  help variables

   integer       :: i, n, n1, g, ngas_check
   real(kind=8)  :: rho1, rho2, col, pp, delp, temp

!  Parameters: Loschmidt's number (particles/cm3), STP parameters

   real(kind=8), parameter ::  RHO_STAND = 2.68675D+19
   real(kind=8), parameter ::  PZERO     = 1013.25D0
   real(kind=8), parameter ::  TZERO     = 273.15D0
   real(kind=8), parameter ::  RHO_ZERO  = RHO_STAND * TZERO / PZERO
   real(kind=8), parameter ::  CONST     = 1.0D+05 * RHO_ZERO
   real(kind=8), parameter ::  DU_TO_CM2 = 2.68668D16
   real(kind=8), parameter ::  O2RATIO   = 0.2095D0

!  initialize

   fail    = .false.
   message = ' '
   
!  Set level data for P, T and H
!    --- Convert to [km] units (height), hPa (pressures)

   do n = 0, GC_nlayers
      n1 = GC_nlayers - n + 1
      temperatures(n) = profile_data(n1,3)
      pressures(n) = profile_data(n1,1)
      heights(n) = profile_data(n1, 2)
   enddo

!  re-calculate heights by Hydrostatic Eqn. (includes TOA)
!  Assumes Surface is Zero. TOA is calculated automatically
     
!   heights(GC_nlayers) = 0.0d0
!   ccon = - 9.81d0 * 28.9d0 / 8314.0d0 * 500.0d0
!  do n = GC_nlayers, 1, -1
!     avit = (1.0d0/temperatures(n-1))+(1.0d0/temperatures(n))
!      heights(n-1) = heights(n) - dlog(pressures(n)/pressures(n-1))/avit/ccon
!   enddo

!  develop air density
!  -------------------

!    Fiddle "pressure-difference method)" 
!         for derivative of Air density w.r.t Temperature

   do n = 1, GC_nlayers
      n1 = n - 1
      rho1 = pressures(n1)/ temperatures(n1)
      rho2 = pressures(n)/ temperatures(n)
      temp = 0.5d0 * (temperatures(n1)+temperatures(n))
      aircolumns(n)     = 0.5d0 * const * ( rho1 + rho2 ) * (heights(n1)-heights(n))
      delp              = pressures(n) - pressures(n1)
      gasconstants(n)   = aircolumns(n) * temp / delp
      layertemp(n)      = temp
      daircolumns_dT(n) = - aircolumns(n) / layertemp(n)
   enddo

! parse profile fields: datanames(1:pnc) corresponds to profile_data(:, 1:pnc)
   read (profids, *) datanames(0:pnc)
  
! find indices for each gases (note, o2 and o4 are special, which are derived for air columns)
   do g = 1, ngases
      gasidxs(g) = -1   ! initialize to -1
      gasname = which_gases(g)
      if (gasname == 'O4  ') gasname = 'O2  '
      do i = 1, pnc
         !if ( TRIM(StrUpCase(datanames(i))) == TRIM(StrUpCase(gasname)) ) then
         if ( TRIM(datanames(i)) == TRIM(gasname) ) then
            gasidxs(g) = i
            CYCLE
         endif
      enddo

      !write(*, '(I4,2A6,I4,A6)') g, which_gases(g), gasname, gasidxs(g), datanames(gasidxs(g))
   enddo
      
!  Develop gas partial columns
!  ---------------------------

!    First default, 3 June 2009. 4 UV gases (O3, NO2, HCHO, SO2)
!    T-derivatives not required explicitly, handled by air column T-deriv above.

   pp = 1.0d0  !1.0d-6
   ngas_check = 0
   do g = 1, ngases
      
      ! Get gas concentration in terms of ppm
      if (gasidxs(g) .ge. 1) then
         do n = 1, GC_nlayers
            n1 = GC_nlayers + 1 - n
            if (profile_is_level) then
               gasconcen(n) = ( profile_data(n1, gasidxs(g)) + profile_data(n1+1, gasidxs(g)) ) /2.
            else
               gasconcen(n) = profile_data(n1, gasidxs(g))
            endif
         enddo
      endif
      
      if ( which_gases(g) .eq. 'O2  ' ) then
         ngas_check = ngas_check + 1
         do n = 1, GC_nlayers
            if (gasidxs(g) == -1) then
               gas_partialcolumns(n,g) = aircolumns(n) * O2RATIO  
            else 
               gas_partialcolumns(n,g) = pp * aircolumns(n) * gasconcen(n) 
            endif
         enddo
      else if ( which_gases(g) .eq. 'O4  ' ) then
         ngas_check = ngas_check + 1
         do n = 1, GC_nlayers
            if (gasidxs(g) == -1) then
               gas_partialcolumns(n,g) = (aircolumns(n) * O2RATIO) ** 2.0 / (heights(n-1)-heights(n)) / 1.0D5  
            else 
               gas_partialcolumns(n,g) = (pp * aircolumns(n) * gasconcen(n) ) ** 2.0 &
                    / (heights(n-1)-heights(n)) / 1.0D5 
            endif
         enddo
      else if (gasidxs(g) .ge. 1) then    ! Other gases
         ngas_check = ngas_check + 1
         do n = 1, GC_nlayers
            gas_partialcolumns(n,g) = pp * aircolumns(n) * gasconcen(n)  
         enddo
      endif
   enddo
   
!  Check that All input gases have been found

   if ( ngas_check .ne. ngases ) then
      message = 'Not all desired trace gases are present in data set: Reduce choice!'
      fail = .true.
      return
   endif

!  Set non-physical entries to zero.

   do g = 1, ngases
      do n = 1, GC_nlayers
         if (gas_partialcolumns(n,g).lt.0.0d0)gas_partialcolumns(n,g)=0.0d0
      enddo
   enddo

!  Develop total columns in [DU]

   do g = 1, ngases
      col = 0.0d0
      do n = 1, GC_nlayers
         col = col + gas_partialcolumns(n,g)
      enddo
      gas_totalcolumns(g) = col / du_to_cm2
   enddo

!  Finish

   return
end subroutine geocape_profile_setter_1


subroutine get_cldprof_from_atmosprof(maxl, maxc, nlayers, pnc, profids, profile_data, &
     maxcld, ncld, cld_types, cld_profile, cld_lowers, cld_uppers)

  implicit none
  
  ! ========================
  ! Input/output parameters
  ! ========================
  integer, intent(in)    :: maxl, maxc, nlayers, pnc, maxcld
  integer, intent(out)   :: ncld
  real(kind=8), dimension(maxl, maxc),  intent(in)         :: profile_data
  character(len=256), intent(in)                           :: profids
  character(len=2), dimension(maxcld), intent(out)         :: cld_types
  integer, dimension (maxcld), intent(out)                 :: cld_lowers, cld_uppers
  real(kind=8), dimension(maxcld, nlayers),  intent(out)   :: cld_profile
  
  ! Local variables
  integer :: i, n, n1
  character(len=4), dimension(0:maxc)   :: datanames  ! datanames(0) should not be used

  ! parse profile fields: datanames(1:pnc) corresponds to profile_data(:, 1:pnc)
  read (profids, *) datanames(0:pnc)

  !OSU, OBC, OOC, OSF, OSC, ODU, OWC, OIC, RH
  ncld = 0
  cld_profile = 0.0d0
  do i = 1, pnc
     if (datanames(i) == 'OWC ' .or. datanames(i) == 'OIC ') then
        if (SUM(profile_data(1:nlayers, i)) .gt. 0.0d0) then
           ncld = ncld + 1
           if (datanames(i) == 'OWC ') then 
              cld_types(ncld) = 'WC'
           else
              cld_types(ncld) = 'IC'
           endif
           
           do n = 1, nlayers
              n1 = nlayers + 1 - n
              cld_profile(ncld, n) = profile_data(n1, i)
           enddo
           
           ! Find cld_uppers, cld_lowers
           do n = 1, nlayers
              if (cld_profile(ncld, n) .gt. 0.0d0) then 
                 cld_uppers(ncld) = n
                 exit
              endif
           enddo
           
           do n = nlayers, 1, -1
              if (cld_profile(ncld, n) .gt. 0.0d0) then 
                 cld_lowers(ncld) = n + 1
                 exit
              endif
           enddo
        endif           
     endif
  enddo

  return
end subroutine get_cldprof_from_atmosprof

subroutine get_aerprof_from_atmosprof(maxl, maxc, nlayers, pnc, profids, profile_data, &
     maxaer, naer, aer_types, aer_profile)

  implicit none
  
  ! ========================
  ! Input/output parameters
  ! ========================
  integer, intent(in)    :: maxl, maxc, nlayers, pnc, maxaer
  integer, intent(out)   :: naer
  real(kind=8), dimension(maxl, maxc),  intent(in)         :: profile_data
  character(len=256), intent(in)                           :: profids
  character(len=2), dimension(maxaer), intent(out)         :: aer_types
  real(kind=8), dimension(maxaer, nlayers),  intent(out)   :: aer_profile
  
  ! Local variables
  integer :: i, n, n1
  character(len=4), dimension(0:maxc)   :: datanames  ! datanames(0) should not be used

  ! parse profile fields: datanames(1:pnc) corresponds to profile_data(:, 1:pnc)
  read (profids, *) datanames(0:pnc)

  !OSU, OBC, OOC, OSF, OSC, ODU, OWC, OIC, RH
  naer = 0
  aer_profile = 0.0d0
  do i = 1, pnc
     if ( datanames(i) == 'OSU ' .or. datanames(i) == 'OBC ' .or. &
          datanames(i) == 'OOC ' .or. datanames(i) == 'OSF ' .or. &
          datanames(i) == 'OSC ' .or. datanames(i) == 'ODU ') then
        if (SUM(profile_data(1:nlayers, i)) .gt. 0.0d0) then
           naer = naer + 1
           aer_types(naer) = datanames(i)(2:3)
           
           do n = 1, nlayers
              n1 = nlayers + 1 - n
              aer_profile(naer, n) = profile_data(n1, i)
           enddo
           
        endif           
     endif
  enddo

  return
end subroutine get_aerprof_from_atmosprof

subroutine insert_clouds(maxlayers, nlayers, heights, pressures, temperatures,    &
     aircolumns, daircolumns_dT, gas_partialcolumns, maxgases, ngases,            &
     do_lambertian_cld, maxcld, ncld, cld_bots, cld_tops, cld_taus, &
     cld_lowers, cld_uppers, cld_opdeps, maxaer, naer, aer_opdeps, fail, message)

  implicit none

  ! ========================
  ! Input/output parameters
  ! ========================
  integer, intent(IN)    :: maxlayers, maxgases, ngases, maxcld, ncld, maxaer, naer
  logical, intent(IN)    :: do_lambertian_cld
  integer, intent(INOUT) :: nlayers
  real(kind=8), dimension (maxcld), intent(IN)               :: cld_taus
  real(kind=8), dimension (maxcld), intent(INOUT)            :: cld_bots, cld_tops
  real(kind=8), dimension (0:maxlayers), intent(INOUT)         :: heights, pressures, temperatures
  real(kind=8), dimension (maxlayers), intent(INOUT)           :: aircolumns, daircolumns_dT
  real(kind=8), dimension (maxlayers, maxgases), intent(INOUT) :: gas_partialcolumns
  character*(*), intent(INOUT)                                 :: message

  logical, INTENT(OUT)                             :: fail
  integer, dimension (maxcld), intent(out)       :: cld_lowers, cld_uppers
  real(kind=8), dimension (maxcld, maxlayers), intent(out)   :: cld_opdeps
  real(kind=8), dimension (maxaer, maxlayers), intent(inout) :: aer_opdeps

  ! ================
  ! Local variables
  ! ================
  integer      :: icld, i, istart
  real(kind=8) :: ext, frac, presfrac

  fail = .false.
  message = ' '
  cld_opdeps = 0.0d0
  cld_lowers = -1
  cld_uppers = -1
  
  ! Use linear inteprolation
  if (do_lambertian_cld) then
     do i = 1, nlayers
        if (abs(cld_tops(1) - heights(i)) < 1.0E-1) then ! merge levels within 100 m
           cld_tops(1) = heights(i)
           cld_uppers(1) = i + 1
           exit
        else if (cld_tops(1) > heights(i) ) then
           
           heights(i+1:nlayers + 1) = heights(i:nlayers)
           pressures(i+1:nlayers + 1) = pressures(i:nlayers)
           temperatures(i + 1:nlayers+1) = temperatures(i:nlayers)
           heights(i) = cld_tops(1)
           cld_uppers(1) = i + 1
           
           frac = (heights(i) - heights(i-1)) / (heights(i + 1) - heights(i - 1))
           temperatures(i) = temperatures(i - 1) * (1.0 - frac) + temperatures(i + 1) * frac
           pressures(i) = EXP(frac * (LOG(pressures(i + 1)) - LOG(pressures(i - 1))) + LOG(pressures(i - 1)))
           presfrac = (pressures(i) - pressures(i - 1)) / (pressures(i+1) - pressures(i-1))
           aircolumns(i + 1: nlayers + 1) = aircolumns(i : nlayers)
           aircolumns(i) =  aircolumns(i + 1)  * presfrac
           aircolumns(i + 1) =  aircolumns(i + 1)  * (1.0 - presfrac)
           gas_partialcolumns(i + 1: nlayers + 1, 1:ngases) = gas_partialcolumns(i : nlayers, 1:ngases)
           gas_partialcolumns(i, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * presfrac
           gas_partialcolumns(i + 1, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * (1.0 - presfrac)
           aer_opdeps(1:naer, i + 1: nlayers + 1) = aer_opdeps(1:naer, i : nlayers)
           aer_opdeps(1:naer, i) =  aer_opdeps(1:naer, i + 1)  * presfrac
           aer_opdeps(1:naer, i + 1) =  aer_opdeps(1:naer, i + 1)  * (1.0 - presfrac)
           daircolumns_dT(i) = -aircolumns(i) / (temperatures(i-1) + temperatures(i)) * 2.0
           daircolumns_dT(i + 1) = -aircolumns(i + 1) / (temperatures(i) + temperatures(i+1)) * 2.0  
           nlayers = nlayers + 1

           if (nlayers > maxlayers) then
              message = 'Need to increase maxlayers!!!'
              fail = .true.; return
           endif
           exit
        endif
     enddo
  else

     istart = 1
     do icld = ncld, 1, -1

        ! insert cloud top
        do i = istart, nlayers
           if (abs(cld_tops(icld) - heights(i)) < 1.0E-2) then
              cld_tops(icld) = heights(i)
              cld_uppers(icld) = i + 1
              exit
           else if (cld_tops(icld) > heights(i) ) then             
              heights(i + 1: nlayers + 1) = heights(i: nlayers)
              pressures(i + 1: nlayers + 1) = pressures(i: nlayers)
              temperatures(i + 1: nlayers + 1) = temperatures(i: nlayers)
              heights(i) = cld_tops(icld)
              
              frac = (heights(i) - heights(i-1)) / (heights(i + 1) - heights(i - 1))
              temperatures(i) = temperatures(i - 1) * (1.0 - frac) + temperatures(i + 1) * frac
              pressures(i) = EXP(frac * (LOG(pressures(i + 1)) - LOG(pressures(i - 1))) + LOG(pressures(i - 1)))
              presfrac = (pressures(i) - pressures(i - 1)) / (pressures(i+1) - pressures(i-1))
              aircolumns(i + 1: nlayers + 1) = aircolumns(i : nlayers)
              aircolumns(i) =  aircolumns(i + 1)  * presfrac
              aircolumns(i + 1) =  aircolumns(i + 1)  * (1.0 - presfrac)
              aer_opdeps(1:naer, i + 1: nlayers + 1) = aer_opdeps(1:naer, i : nlayers)
              aer_opdeps(1:naer, i) =  aer_opdeps(1:naer, i + 1)  * presfrac
              aer_opdeps(1:naer, i + 1) =  aer_opdeps(1:naer, i + 1)  * (1.0 - presfrac)
              gas_partialcolumns(i + 1: nlayers + 1, 1:ngases) = gas_partialcolumns(i : nlayers, 1:ngases)
              gas_partialcolumns(i, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * presfrac
              gas_partialcolumns(i + 1, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * (1.0 - presfrac)
              daircolumns_dT(i) = -aircolumns(i) / (temperatures(i-1) + temperatures(i)) * 2.0
              daircolumns_dT(i + 1) = -aircolumns(i + 1) / (temperatures(i) + temperatures(i+1)) * 2.0   

              nlayers = nlayers + 1
              cld_uppers(icld) = i + 1
              
              if (nlayers > maxlayers) then
                 message = 'Need to increase maxlayers!!!'
                 fail = .true.; return
              endif
              exit
           endif
           
        enddo

        ! insert cloud bottom
        istart = i
        do i = istart, nlayers
           if (abs(cld_bots(icld) - heights(i)) < 1.0E-2) then
              cld_bots(icld) = heights(i)
              cld_lowers(icld) = i 
              exit
           else if (cld_bots(icld) > heights(i) ) then
              
              heights(i + 1: nlayers + 1) = heights(i: nlayers)
              pressures(i + 1: nlayers + 1) = pressures(i: nlayers)
              temperatures(i + 1: nlayers + 1) = temperatures(i: nlayers)
              heights(i) = cld_bots(icld)            
              
              frac = (heights(i) - heights(i-1)) / (heights(i + 1) - heights(i - 1))
              temperatures(i) = temperatures(i - 1) * (1.0 - frac) + temperatures(i + 1) * frac
              pressures(i) = EXP(frac * (LOG(pressures(i + 1)) - LOG(pressures(i - 1))) + LOG(pressures(i - 1)))
              presfrac = (pressures(i) - pressures(i - 1)) / (pressures(i+1) - pressures(i-1))
              aircolumns(i + 1: nlayers + 1) = aircolumns(i : nlayers)
              aircolumns(i) =  aircolumns(i + 1)  * presfrac
              aircolumns(i + 1) =  aircolumns(i + 1)  * (1.0 - presfrac)
              aer_opdeps(1:naer, i + 1: nlayers + 1) = aer_opdeps(1:naer, i : nlayers)
              aer_opdeps(1:naer, i) =  aer_opdeps(1:naer, i + 1)  * presfrac
              aer_opdeps(1:naer, i + 1) =  aer_opdeps(1:naer, i + 1)  * (1.0 - presfrac)
              gas_partialcolumns(i + 1: nlayers + 1, 1:ngases) = gas_partialcolumns(i : nlayers, 1:ngases)
              gas_partialcolumns(i, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * presfrac
              gas_partialcolumns(i + 1, 1:ngases) =  gas_partialcolumns(i + 1, 1:ngases)  * (1.0 - presfrac)
              daircolumns_dT(i) = -aircolumns(i) / (temperatures(i-1) + temperatures(i)) * 2.0
              daircolumns_dT(i + 1) = -aircolumns(i + 1) / (temperatures(i) + temperatures(i+1)) * 2.0   

              nlayers = nlayers + 1
              cld_lowers(icld) = i 
              
              if (nlayers > maxlayers) then
                 message = 'Need to increase maxlayers!!!'
                 fail = .true.; return
              endif
              exit
           endif
           
        enddo
        istart = i
        
        ext = cld_taus(icld) / (heights(cld_uppers(icld)-1) - heights(cld_lowers(icld)))
        cld_opdeps(icld, cld_uppers(icld):cld_lowers(icld)) = &
             (heights(cld_uppers(icld)-1:cld_lowers(icld)-1) - heights(cld_uppers(icld):cld_lowers(icld))) * ext
        
     enddo
     
  endif
  
  return
end subroutine insert_clouds

subroutine  convert_cldspec_p2z(nl, ps, zs, do_lambertian_cld, use_cldprof, mc, nc, ctops, cbots)
  
  implicit none

  ! ========================
  ! Input/output parameters
  ! ========================
  integer, intent(IN)    :: nl, nc, mc
  logical, intent(IN)    :: do_lambertian_cld, use_cldprof
  real(kind=8), dimension (0:nl), intent(IN)    :: zs, ps
  real(kind=8), dimension (mc), intent(INOUT)   :: cbots, ctops

  ! ================
  ! Local variables
  ! ================
  integer                       :: i
  real(kind=8)                  :: ptmp, ztmp
  real(kind=8), dimension(0:nl) :: y2, plog
  
  if (.not. do_lambertian_cld .and. use_cldprof) return

  plog = log10(ps)
  call spline(plog(0:nl), zs(0:nl), nl+1,0.0d0,0.0d0,y2)
  
  if (do_lambertian_cld) then
     ptmp = log10(ctops(1))
     call splint(plog(0:nl), zs(0:nl), y2, nl+1, ptmp, ztmp) 
     ctops(1) = ztmp
  else 
     do i = 1, nc
        ptmp = log10(ctops(i))
        call splint(plog(0:nl), zs(0:nl), y2, nl+1, ptmp, ztmp) 
        ctops(i) = ztmp

        ptmp = log10(cbots(i))
        call splint(plog(0:nl), zs(0:nl), y2, nl+1, ptmp, ztmp) 
        cbots(i) = ztmp        
     enddo
  endif

return
end subroutine convert_cldspec_p2z
