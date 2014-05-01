
subroutine geocape_surface_setter_1                  &
    ( surface_data_path, latitude, longitude, month, & ! inputs
      maxlambdas, nlambdas, lambdas,                 & ! inputs
      ground_ler,                                    & ! Output
      message, fail )                                  ! Exception handling

!  Based on Latitude (range [-90,90], Longitude [-180,180], and month
!  This routine Extracts the TEMIS  surface albedo for a range of
!  specified wavelengths

!  inputs
!  ------

!  Directory path

   character*(*), intent(in) :: surface_data_path

!  Dimensioning

   integer, intent(in)       :: maxlambdas

!  geo-information

   real(kind=8), intent(in)  :: latitude, longitude
   integer, intent(in)       :: month

!  wavelengths

   integer, intent(in)       :: nlambdas
   real(kind=8), intent(in)  :: lambdas ( maxlambdas )

!  output
!  ------

   real(kind=8), intent(out) :: ground_ler ( maxlambdas )

   logical,       intent(inout)  :: fail
   character*(*), intent(inout)  :: message

!  local
!  -----

   real(kind=8)       :: lams(11), start, finis
   real(kind=8)       :: x1, x2, a1, a2, ler1(11), grad, wav
   integer            :: ms, mf, m, m1, m2, p1, p2, q1, q2, n, np
   character(len=100) :: filename
   character(len=75)  :: header(3)

   character(len=2)   :: cmonth(12)
   character(len=3)   :: cwavs(11)
   integer            :: iler(180,360)

   data   cmonth / '01', '02', '03', '04', '05', '06', &
                   '07', '08', '09', '10', '11', '12' /
   data   cwavs  / '335', '380', '416', '440', '463', '494', &
                   '555', '610', '670', '758', '772' /
   data lams     / 335.0, 380.0, 416.0, 440.0, 463.0, 494.5, &
                   555.0, 610.0, 670.0, 758.0, 772.0 /

!  Initialise

   fail    = .false.
   message = ' '

!  Determine which data you need

   start = lambdas(1)
   finis = lambdas (nlambdas)
   if ( start .le. lams(1) )then
      ms = 1
   else
      ms = 1
      do while (start .gt. lams(ms))
         ms = ms + 1
      enddo
      ms = ms - 1
   endif
   if ( finis .ge. lams(11) ) then
      mf = 11
   else
      mf = ms
      do while (finis .gt. lams(mf))
         mf = mf + 1
      enddo
   endif

!  interpolation points

   if ( latitude.ge.89.5d0 ) then
      p1 = 180
      p2 = 180
   else if ( latitude.le.-89.5d0 ) then
      p1 = 1
      p2 = 1
   else
      p2 = 1
      do while (latitude .gt. dble(p2-90)+0.5d0 )
         p2 = p2 + 1
      enddo
      p1 = p2 - 1
   endif

   if ( longitude.le.-179.5d0 .or. longitude.ge.179.5d0) then
      q1 = 360
      q2 = 1
   else
      q2 = 1
      do while ( longitude .gt.  dble(q2-180)+0.5d0 )
         q2 = q2 + 1
      enddo
      q1 = q2 - 1
   endif

!  Bilinear coefficients

   x2 = dble(p2-90)+0.5d0 - latitude
   x1 = latitude - dble(p1-90)-0.5d0
   a2 = dble(q2-180)+0.5d0 - longitude
   a1 = longitude - dble(q1-180)-0.5d0

! debug
!      write(*,*)p1,p2,q1,q2
!      write(*,*)x1,x2
!      write(*,*)a1,a2

!  Length of directory name

   np = 1
   do while (surface_data_path(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Geographical interpolation (bilinear)

   do m = ms, mf
      filename = surface_data_path(1:np)//&
         'sacspecTOTL'//cmonth(month)//'_'//cwavs(m)//'.dat'
      call rd_cscdb(filename,header,iler,fail,message)
      if ( fail ) return
      ler1(m) = x2*a2*dble(iler(p1,q1))+x1*a2*dble(iler(p2,q1)) + &
                x2*a1*dble(iler(p1,q2))+x1*a1*dble(iler(p2,q2))
      ler1(m) = ler1(m) * 1.0d-03
   enddo

!  Wavelength assignments (linear interpolation)

   do n = 1, nlambdas
      wav = lambdas(n)
      if ( wav.le.lams(1) ) then
         m1 = 1
         ground_ler(n) = ler1(m1)
      else if ( wav.ge.lams(11) ) then
         m2 = mf
         ground_ler(n) = ler1(m2)
      else
         m2 = 1
         do while (wav .gt. lams(m2) )
            m2 = m2 + 1
         enddo
         m1 = m2 - 1
         grad = ( ler1(m2) - ler1(m1) ) / (lams(m2) - lams(m1))
         ground_ler(n) = ler1(m1) + grad * ( wav - lams(m1) )
      endif
   enddo

!  Finish

   return
end subroutine geocape_surface_setter_1

subroutine rd_cscdb(filename,header,iler,fail,message)
   integer              :: i,j
   character(len=*)     :: filename
   character(len=75)    :: header(3)
   character(len=14)    :: string
   integer              :: iler(180,360)
   logical              :: fail
   character*(*)        :: message

!  Initialise

   fail    = .false.
   message = ' '
 
!  Open surface TEMIS file (adapted code from Koelemeijer)

   open(10,file=filename(1:LEN(filename)),err = 90, status = 'old')
   do i=1,3
      read(10,'(a)') header(i)
   enddo
   do i=1,180
      read(10,91) (iler(i,j),j=1,360),string
   enddo
   close(10)
91 format(14(25i3/),10i3,a14)
   return

!  Error return

90 continue
   fail    = .true.
   message = 'Open failure for Surface file = '//filename(1:LEN(filename))
   return

end subroutine rd_cscdb


subroutine geocape_surface_setter_2                            &
    ( albspectra_fname, maxlambdas, nlambdas, lambdas,         & ! inputs
      ground_ler,  message, fail )                               ! outputs & Exception handling

  implicit none

!  inputs
!  ------

!  filename

   character*(*), intent(in) :: albspectra_fname

!  Dimensioning

   integer, intent(in)       :: maxlambdas

!  wavelengths

   integer, intent(in)       :: nlambdas
   real(kind=8), intent(in)  :: lambdas ( maxlambdas )

!  output
!  ------

   real(kind=8), intent(out) :: ground_ler ( maxlambdas )

   logical,       intent(inout)  :: fail
   character*(*), intent(inout)  :: message

!  local
!  -----
   integer, parameter :: maxnalb = 3000
   integer :: nalb, fidx, lidx, i, ntmp, np, errstat
   real(kind=8), dimension(maxnalb) :: albspectra_lams, albspectra

   fail = .false.; errstat = 0
   message = ' '
   
   np = 1
   do while (albspectra_fname(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

   open(1, file=albspectra_fname(1:np), err =90, status='old')
   read(1, *, err=90) nalb
   do i = 1, nalb
      read(1, *, err=90) albspectra_lams(i), albspectra(i)
   enddo
   close(1)

   albspectra_lams(1:nalb) = albspectra_lams(1:nalb) * 1000.0
   albspectra(1:nalb) = albspectra(1:nalb) / 100.0
   
   fidx = minval(minloc(lambdas(1:nlambdas), mask=(lambdas(1:nlambdas) >= albspectra_lams(1))))
   if (fidx > 1) ground_ler(1:fidx-1) = albspectra(1)
   
   lidx = minval(maxloc(lambdas(1:nlambdas), mask=(lambdas(1:nlambdas) <= albspectra_lams(nalb))))
   if (lidx > 0 .and. lidx < nlambdas) ground_ler(lidx+1:nlambdas) = albspectra(nalb)
   
   if (fidx > 0 .and. lidx > 0 .and. fidx < lidx) then
      ntmp = lidx - fidx + 1
      call bspline(albspectra_lams(1:nalb), albspectra(1:nalb), nalb, &
           lambdas(fidx:lidx), ground_ler(fidx:lidx), ntmp, errstat)
   endif

   if (errstat < 0) go to 90
   return

! Finish

90 continue
   
   message = 'Error in surface albedo interpolation!!!'
   fail = .true.

   return
end subroutine geocape_surface_setter_2
