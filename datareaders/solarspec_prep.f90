!xliu, 5/1/2013, remove nr=49951 to nwmax=80093 
! as one of the solar irradiance has 80093 lines
! nr will be determined during the reading process
subroutine solarspec_prep_newkur &
           ( filename, maxwavnums, nwavnums, wavnums,  & ! input
             solar_spec_data, message, fail )            ! output

!  The reader to read in solar spectrum data in a .dat file

!  input filename

   character(len=*),               intent(in)  :: filename

!  other inputs

!  Dimensioning

   integer, intent(in)       :: maxwavnums

!  wave numbers

   integer, intent(in)       :: nwavnums
   real(kind=8), intent(in)  :: wavnums ( maxwavnums )

!  Main output

   real(kind=8), dimension(maxwavnums), intent(out) :: solar_spec_data

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables
   integer, parameter             :: nwmax = 80095
   integer                        :: i, nr, np, io
   character(LEN=10)              :: charnw
   real(kind=8), dimension(nwmax) :: solar_spec_raw_data   
   real(kind=8), dimension(nwmax) :: wn_raw_data           
   real(kind=8)                   :: grad, wav
   integer                        :: m1, m2, n    

!  initialize

   fail    = .false.
   message = ' '
   !nr = 49951  !xliu
   solar_spec_data(:) = 0.d0

   np = 1
   do while (filename(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Open and read file
   open(1,file=filename(1:np),err=90, status='old')

   do i = 1,2
      read(1,*)
   enddo

! xliu, disable the use of fixed nr here as file may have different # of wavelengths
!   do i = 1, nr
!     read(1,*) wn_raw_data(i),solar_spec_raw_data(i)
!   enddo
   io = 0
   i = 1
   do 
      read(1, *, iostat=io) wn_raw_data(i),solar_spec_raw_data(i)
      if (io .lt. 0) then
         exit
      else if (io .gt. 0) then
         go to 90
      else
         i = i + 1
      endif
   enddo
   nr = i - 1
   close(1)

   if (nr .gt. nwmax) then 
      fail = .true.
      write(charnw, '(I10)') nr
      message = 'Increase nwmax in solarspec_prep_newkur to ' // charnw
      return
   endif  
!xliu */

! Check that the high resolution solar spectra covers the spectral
! region of interest. If not print warning message to screen and
! continue
   IF (wn_raw_data(1) .GT. wavnums(nwavnums) .OR. &
        wn_raw_data(nr) .LT. wavnums(1) ) THEN

      print*, '!!!Warning: Spectral window not covered by high resolution solar spectrum:'
      print 1000, '    High resolution solar spectrum:', wn_raw_data(1),' nm -',  &
           wn_raw_data(nr),' nm'
      print 1000, '    Spectral window calculation:', wavnums(nwavnums),' nm -', &
           wavnums(1),' nm'
      
      1000 FORMAT (A,F9.2,A,F9.2,A)
   ENDIF

!  Wavelength assignments (linear interpolation)
   do n = 1, nwavnums
      wav = wavnums(n)  
      if ( wav.le.wn_raw_data(1) ) then
         solar_spec_data(n) = solar_spec_raw_data(1)
      else if ( wav.ge.wn_raw_data(nr) ) then
         solar_spec_data(n) = solar_spec_raw_data(nr)
      else
         m2 = 1
         do while (wav .gt. wn_raw_data(m2) )
            m2 = m2 + 1
         enddo 
         m1 = m2 - 1
         grad = ( solar_spec_raw_data(m2) - solar_spec_raw_data(m1) ) / &
                ( wn_raw_data(m2) - wn_raw_data(m1) )
         solar_spec_data(n) = solar_spec_raw_data(m1) + grad * ( wav - wn_raw_data(m1) )
      endif
   enddo

   ! Convert it to W/ m^2 / cm^-1, added by xliu
   solar_spec_data(1:nwavnums) = solar_spec_data(1:nwavnums) * 1.0d4

   return

! error return

90 continue
   fail = .true.
   message = 'Open/read failure for solarspec file = '//filename(1:LEN(filename))
   return

 end subroutine solarspec_prep_newkur

subroutine solarspec_prep_chance &
           ( filename, maxlambdas, nlambdas, lambdas,  & ! input
             solar_spec_data, message, fail )            ! output

!  The reader to read in solar spectrum data in a .dat file

!  input filename

   character(len=*),               intent(in)  :: filename

!  other inputs

!  Dimensioning

   integer, intent(in)       :: maxlambdas

!  wave numbers

   integer, intent(in)       :: nlambdas
   real(kind=8), intent(in)  :: lambdas ( maxlambdas )

!  Main output

   real(kind=8), dimension(maxlambdas), intent(out) :: solar_spec_data

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables

   integer      :: i, nr, np
   real(kind=8), dimension(80093) :: solar_spec_raw_data
   real(kind=8), dimension(80093) :: wn_raw_data
   real(kind=8)       :: grad, wav
   integer            :: m1, m2, n    

!  initialize

   fail    = .false.
   message = ' '
   nr = 80093

   solar_spec_data(:) = 0.d0

   np = 1
   do while (filename(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Open and read file

   open(1,file=filename(1:np),err=90, status='old')

   do i = 1, nr
     read(1,*) wn_raw_data(i),solar_spec_raw_data(i)
   enddo

   close(1)

!  Wavelength assignments (linear interpolation)
   do n = 1, nlambdas
      wav = lambdas(n)  
      if ( wav.le.wn_raw_data(1) ) then
         solar_spec_data(n) = solar_spec_raw_data(1)
      else if ( wav.ge.wn_raw_data(nr) ) then
         solar_spec_data(n) = solar_spec_raw_data(nr)
      else
         m2 = 1
         do while (wav .gt. wn_raw_data(m2) )
            m2 = m2 + 1
         enddo 
         m1 = m2 - 1
         grad = ( solar_spec_raw_data(m2) - solar_spec_raw_data(m1) ) / &
                ( wn_raw_data(m2) - wn_raw_data(m1) )
         solar_spec_data(n) = solar_spec_raw_data(m1) + grad * ( wav - wn_raw_data(m1) )
      endif
   enddo

   return

! error return

90 continue
   fail = .true.
   message = 'Open failure for solarspec file = '//filename(1:LEN(filename))
   return

 end subroutine solarspec_prep_chance


! Calculate sun-earth distance correction factor ratio between daily value 
! of the sun-earth distance and its mean yearly value (1 AU)
subroutine sun_earth_distance (jday, rdis)

  implicit none
  !=================
  ! input variables
  !=================
  real(kind=8), intent (in)  :: jday
  
  ! ================
  ! output variables
  ! ================
  real(kind=8), intent (out) :: rdis
  
  !local variables
  real(KIND=8), PARAMETER :: pi = 3.14159265358979d0
  real(kind=8)            :: phi
  
  phi = 2.0 * pi * (jday - 1.d0) / 365.0
  rdis = sqrt(1.0d0/(1.00011d0+0.034221d0 * cos(phi) + 0.00128d0 * sin(phi) + &
       0.000719d0 * cos(2.d0*phi) + 0.000077d0 * sin(2.0 * phi)))
  
  return
  
end subroutine sun_earth_distance

