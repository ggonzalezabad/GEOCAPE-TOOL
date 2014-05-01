! -------------------------------------------------------------------

subroutine geocape_target_reader &
           ( filename,  &     ! input
             geometry_data, time_data, message, fail ) ! output

!  The reader to read in TES geometry information in a .asc file

!  input: filename
!  output: geometry_data(3), time_data(3)

!  TES geometry_data(3) has following data
!    1.   SZA (degrees),
!    2.   VZA (degrees),
!    3.   AZM (degrees)

!  TES time_data(3) has following data
!    1.   Year,
!    2.   Month,
!    3.   Date

!  inptu filename

   character(len=*),               intent(in)  :: filename

!  Main output

   real(kind=8), dimension(3),     intent(out) :: geometry_data
   integer, dimension(3),          intent(out) :: time_data

!  Exception handling

   logical,       intent(INOUT) :: fail
   character*(*), intent(INOUT) :: message

!  Local variables

   integer      :: i, nr, nc, np
   character*80 :: dummy, time_string

!  initialize

   fail    = .false.
   message = ' '

   np = 1
   do while (filename(np:np).ne. ' ')
     np = np + 1 
   enddo
   np = np - 1

!  Open and read file

   open(1,file=filename(1:np),err=90, status='old')
   read(1,*)
   read(1,*) dummy,dummy,nr,dummy,nc
   read(1,*)
   read(1,*) dummy,dummy,time_string
   do i = 1,3 
      read(1,*)
   enddo
   read(1,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy, &
             geometry_data

!  Extract year/month/date from the "time_string" variable

   read(time_string(1:4),'(I4)') time_data(1)
   read(time_string(6:7),'(I2)') time_data(2)
   read(time_string(9:10),'(I2)') time_data(3)

   close(1)

   return

! error return

90 continue
   fail = .true.
   message = 'Open failure for profile file = '//filename(1:LEN(filename))
   return

end subroutine geocape_target_reader
