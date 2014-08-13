
subroutine geocape_xsec_setter_1                                            &
    ( maxlambdas, maxgases, maxlayers, xsec_data_path, xsec_data_filenames, & ! Input
      nlambdas, lambdas, ngases, which_gases,                               & ! Input
      gas_xsecs_type, nlayers, mid_pressures, temperatures,                 & ! Input
      gas_xsecs, o3c1_xsecs, o3c2_xsecs,                                    & ! Output
      Rayleigh_xsecs, Rayleigh_depols,                                      & ! Output
!      message, fail, gas_xsecsdt)                                            ! Errors
      message, fail) 

  implicit none

!  Inputs
!  ------

!  Dimensioning

   integer,          intent(in) :: maxlambdas, maxgases, maxlayers

!  File paths and names

   character(LEN=*), intent(in) :: xsec_data_path
   character(LEN=*), intent(in) :: xsec_data_filenames(maxgases)

!  gases

   integer,                                intent(in) :: ngases, nlayers
   character(Len=4), dimension (maxgases), intent(in) :: which_gases
   integer, dimension (maxgases), intent(in)          :: gas_xsecs_type

!  Window of wavelengths

   integer,                                intent(in) :: nlambdas
   real(kind=8), dimension(maxlambdas),    intent(in) :: lambdas

! Mean pressures and temperatures
   real(kind=8), dimension(maxlayers), intent(in)     :: mid_pressures, temperatures   

!  output
!  ------

!  Trace gases

   real(kind=8), dimension(maxlambdas, maxlayers, maxgases), intent(out) :: gas_xsecs
!   real(kind=8), dimension(maxlambdas, maxlayers, maxgases), intent(out), optional :: gas_xsecsdt
   real(kind=8), dimension(maxlambdas),          intent(out) :: o3c1_xsecs
   real(kind=8), dimension(maxlambdas),          intent(out) :: o3c2_xsecs

!  Rayleigh

   real(kind=8), dimension(maxlambdas),          intent(out) :: Rayleigh_xsecs
   real(kind=8), dimension(maxlambdas),          intent(out) :: Rayleigh_depols

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables

   integer, parameter :: maxspline = 100000
   character(LEN=256) :: filename, xsec_file_name

   logical       :: reading, is_wavenum
   integer       :: nbuff, ndata, n, g, np, nf, ngas_check, errstat, io
   real(kind=8)  :: lamb1, lamb2, wav, val, c1, c2,conversion, xsec, fwhm
   real(kind=8)  :: CO2_PPMV_MIXRATIO
   real(kind=8)  :: x(maxspline), y(maxspline), y2(maxspline)
   real(kind=8)  :: yc1(maxspline), yc2(maxspline)
   real(kind=8)  :: pressures(maxlayers)
   character(LEN=6) :: the_molecule

!  initialize exception handling

   fail    = .false.
   message = ' '

!  Initialize output

   gas_xsecs = 0.0d0;  o3c1_xsecs = 0.0d0;  o3c2_xsecs = 0.0d0
   Rayleigh_xsecs = 0.0d0; Rayleigh_depols = 0.0d0

!  Length of directory name

   np = 1
   do while (xsec_data_path(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  gas check starter

   ngas_check = 0

!  Adopt a simple buffering system, windows 0.5 nm extended
!    This can be varied to suit the particular data set
!   lamb1 = lambdas(1)        - 0.5d0
!   lamb2 = lambdas(nlambdas) + 0.5d0

!  ##################################################################

!  Start the Gas loop

   do g = 1, ngases

      !print *, g, which_gases(g), gas_xsecs_type(g), TRIM(ADJUSTL(xsec_data_filenames(g)))

      if (gas_xsecs_type(g) .eq. 1 .or. gas_xsecs_type(g) .eq. 2) then

!  Check off this gas
         ngas_check = ngas_check + 1

!  Buffer control, initialization
         lamb1 = lambdas(1)        - 0.05d0
         lamb2 = lambdas(nlambdas) + 0.05d0
         conversion = 1.0d+20

!  Prepare filename
         xsec_file_name = xsec_data_filenames(g)
         nf = 1
         do while (xsec_file_name(nf:nf).ne. ' ')
           nf = nf + 1 
         enddo
         nf = nf - 1
         filename = xsec_data_path(1:np)//xsec_file_name(1:nf)

         nbuff = 0
         ndata = 0
      endif

!  temperature-independent cross section from file
!  ===

      if ( gas_xsecs_type(g) .eq. 1) then

!  Read the data, check that window is covered

         open(1, file = filename, err = 90, status='old')
!  .. Read first line, and check window is covered
         read(1,*, iostat=io) wav, val
         if (lamb1 .lt. wav ) then
            message = which_gases(g) // ' Xsec data does not cover input window at lower end'
!            print *, which_gases(g) // ' Xsec data does not cover input window at lower end'
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. io .eq. 0 )
            ndata = ndata + 1
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff) = wav
                  y(nbuff) = val * conversion
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
            read (1,*,iostat=io) wav, val
         enddo

!  .. Check if last data line is present, then window not covered
         if (reading) then
            message = which_gases(g) // ' Xsec data does not cover input window at upper end'
!            print *, which_gases(g) // ' Xsec data does not cover input window at upper end'
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in geocape_xsec_setter_1: '
            go to 91
         endif
                          
!  Spline the data. Replace with BSPLINES ???????????????

         call spline(x,y,nbuff,0.0d0,0.0d0,y2)
         do n = 1, nlambdas
            if (lambdas(n) >= x(1) .and. lambdas(n) <= x(nbuff)) then 
               call splint(x,y,y2,nbuff,lambdas(n),xsec)
               if (xsec .gt. 0.0) gas_xsecs(n,1,g) = xsec / conversion
            else
               gas_xsecs(n, 1, g) = 0.d0
            endif
         enddo

         !print *, nbuff, nlambdas, g, which_gases(g)
         !print *, x(1), y(1)
         !print *, x(2), y(2)
         !print *, x(nbuff-1), y(nbuff-1)
         !print *, x(nbuff), y(nbuff)

!  For O3
!  ==

      else if ( gas_xsecs_type(g) .eq. 2) then

!  Read the data, check that window is covered

         open(1, file = filename, err = 90, status='old')
!  .. First line is a dummy
         read(1,*)
!  .. Read second line, and check window is covered
         read(1,*, iostat = io) wav, val, c1, c2
         if (lamb1 .lt. wav ) then
!            print *, 'O3 Xsec data does not cover input window at lower end'
         endif
!  .. Read more lines, saving to spline buffer
         reading = .true.
         do while (reading .and. io .eq. 0 )
            ndata = ndata + 1
            if ( wav .gt. lamb1 ) then
               if ( wav .lt. lamb2 ) then
                  nbuff = nbuff + 1
                  x(nbuff)   = wav
                  y(nbuff)   = val
                  yc1(nbuff) = c1
                  yc2(nbuff) = c2
               endif
               if ( wav .ge. lamb2 ) reading = .false.
            endif
            read(1,*, iostat = io ) wav, val, c1, c2
         enddo

!  .. Check if last data line is present, then window not covered
         if (reading) then
!            print *, 'O3 Xsec data does not cover input window at upper end'
         endif
         close(1)

         if (nbuff .gt. maxspline) then
            fail = .true.
            message = 'need to increase maxspline in geocape_xsec_setter_1: '
            go to 91
         endif
         
!  Spline the data
         call spline(x,y,nbuff,0.0d0,0.0d0,y2)

         do n = 1, nlambdas
            if (lambdas(n) >= x(1) .and. lambdas(n) <= x(nbuff)) then 
               call splint(x,y,y2,nbuff,lambdas(n),xsec)
               gas_xsecs(n,1,g) = xsec / conversion
            else
               gas_xsecs(n,1,g) = 0.0d0            
            endif

         enddo

         call spline(x,yc1,nbuff,0.0d0,0.0d0,y2)
         do n = 1, nlambdas
            if (lambdas(n) >= x(1) .and. lambdas(n) <= x(nbuff)) then 
               call splint(x,yc1,y2,nbuff,lambdas(n),xsec)
               o3c1_xsecs(n) = xsec / conversion
            else
               o3c1_xsecs(n) = 0.d0
            endif
         enddo

         call spline(x,yc2,nbuff,0.0d0,0.0d0,y2)
         do n = 1, nlambdas
            if (lambdas(n) >= x(1) .and. lambdas(n) <= x(nbuff)) then 
               call splint(x,yc2,y2,nbuff,lambdas(n),xsec)
               o3c2_xsecs(n) = xsec / conversion
            else
               o3c2_xsecs(n) = 0.0d0
            endif
   
         enddo
         gas_xsecs(1:nlambdas, 2, g) = o3c1_xsecs(1:nlambdas)
         gas_xsecs(1:nlambdas, 3, g) = o3c2_xsecs(1:nlambdas)
         
! Use HITRAN database
      else if (gas_xsecs_type(g) .eq. 3 .and. trim(adjustl(xsec_data_filenames(g))) .eq. 'HITRAN') then

         the_molecule = which_gases(g) // '  '
         is_wavenum = .FALSE.
         fwhm=0.d0
         pressures(1:nlayers) = mid_pressures(1:nlayers) / 1013.25
                  
         !if (PRESENT(gas_xsecsdt)) then  ! Currently, optional arguments does not work
         !
         !   call get_hitran_crs(the_molecule, nlambdas, lambdas(1:nlambdas), is_wavenum, nlayers, &
         !        pressures(1:nlayers), temperatures(1:nlayers), fwhm, gas_xsecs(1:nlambdas, 1:nlayers, g), &
         !        errstat, gas_xsecsdt(1:nlambdas, 1:nlayers, g))
         !else
         call get_hitran_crs(the_molecule, nlambdas, lambdas(1:nlambdas), is_wavenum, nlayers, &
              pressures(1:nlayers), temperatures(1:nlayers), fwhm, gas_xsecs(1:nlambdas, 1:nlayers, g), errstat)
         !endif

         if (errstat == 1) then
            message='Error in getting HITRAN cross sections for molecule ' // TRIM(the_molecule)
            fail = .true.
            go to 91
         else
            ! Check off this gas           
            ngas_check = ngas_check + 1
         endif
!  Finish cross - sections

      endif

   enddo

!  CHeck that all gases have been accounted for

!   print *, 'ngas_check = ', ngas_check, ' ngases = ', ngases

   if ( ngas_check .ne. ngases ) then
      fail    = .true.
      message = ' Not all gases accounted for; Reduce choices!'
      go to 91
   endif

!  RAYLEIGH
!  ========

!  Set CO2 mixing ratio

   CO2_PPMV_MIXRATIO = 385.0d0

!  Call

   call RAYLEIGH_FUNCTION                       &
          ( MAXLAMBDAS, CO2_PPMV_MIXRATIO,      &
            NLAMBDAS,   LAMBDAS,                &
            RAYLEIGH_XSECS, RAYLEIGH_DEPOLS )

!  normal return

   return

!  error returns

90 continue
   fail    = .true.
   message = 'Open failure for Xsec file = '//filename(1:LEN(filename))
   return

91 continue
   return

end subroutine geocape_xsec_setter_1


SUBROUTINE RAYLEIGH_FUNCTION                       &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS,   &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.
!     Wavelengths in NM

!  Input arguments
!  ---------------

!  wavelength
 
      INTEGER     :: FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      real(kind=8), dimension ( FORWARD_MAXLAMBDAS ) :: FORWARD_LAMBDAS 

!  CO2 mixing ratio

      real(kind=8) :: CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      real(kind=8), dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_XSEC 
      real(kind=8), dimension ( FORWARD_MAXLAMBDAS ) :: RAYLEIGH_DEPOL

!  Local variables
!  ---------------

      INTEGER      :: W
      real(kind=8) :: MASS_DRYAIR
      real(kind=8) :: NMOL, PI, CONS
      real(kind=8) :: MO2,MN2,MARG,MCO2,MAIR
      real(kind=8) :: FO2,FN2,FARG,FCO2,FAIR
      real(kind=8) :: LAMBDA_C,LAMBDA_M,LPM2,LP2
      real(kind=8) :: N300M1,NCO2M1,NCO2
      real(kind=8) :: NCO2SQ, NSQM1,NSQP2,TERM,WAV

!  data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      real(kind=8), PARAMETER ::        S0_A = 15.0556D0
      real(kind=8), PARAMETER ::        S0_B = 28.9595D0

      real(kind=8), PARAMETER ::        S1_A = 8060.51D0
      real(kind=8), PARAMETER ::        S1_B = 2.48099D+06
      real(kind=8), PARAMETER ::        S1_C = 132.274D0
      real(kind=8), PARAMETER ::        S1_D = 1.74557D+04
      real(kind=8), PARAMETER ::        S1_E = 39.32957D0

      real(kind=8), PARAMETER ::        S2_A = 0.54D0

      real(kind=8), PARAMETER ::        S3_A = 1.034D0
      real(kind=8), PARAMETER ::        S3_B = 3.17D-04
      real(kind=8), PARAMETER ::        S3_C = 1.096D0
      real(kind=8), PARAMETER ::        S3_D = 1.385D-03
      real(kind=8), PARAMETER ::        S3_E = 1.448D-04

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop

      DO W = 1, FORWARD_NLAMBDAS

!  Convert to Angstroms

      WAV = FORWARD_LAMBDAS(W) * 10.0d0

!  wavelength in micrometers

      LAMBDA_M = 1.0D-04 * WAV
      LAMBDA_C = 1.0D-08 * WAV
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

!  step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  end loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION 


  SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER            ::  n
      REAL(KIND=8)       :: yp1,ypn,x(n),y(n),y2(n)
!xl      INTEGER, PARAMETER :: NMAX=6094
      INTEGER            :: i,k
!xl     REAL(Kind=8)       :: p,qn,sig,un,u(NMAX)
      REAL(Kind=8)       :: p,qn,sig,un,u(n)

!  Check dimension

!xl      if (n.gt.NMAX) then
!xl        write(*,*)'Error in spline routine: too small NMAX =',NMAX
!xl        stop
!xl      endif

      if (yp1.gt..99e30) then
        y2(1)=0.0d0
        u(1)=0.0d0
      else
        y2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.0d0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - &
                     (y(i)-y(i-1)) / (x(i)-x(i-1))   &
                   ) / (x(i+1)-x(i-1)) - sig*u(i-1)  &
            )/p
      enddo

      if (ypn.gt..99d30) then
        qn=0.0d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo

      return
  END SUBROUTINE SPLINE


  SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER            ::  n
      REAL(KIND=8)       :: x, y, xa(n),ya(n),y2a(n)
      INTEGER            ::  k,khi,klo
      REAL(KIND=8)       :: a,b,h

      klo=1
      khi=n

 1    if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.0d0) stop
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+ &
        ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0

      return
  END SUBROUTINE SPLINT


