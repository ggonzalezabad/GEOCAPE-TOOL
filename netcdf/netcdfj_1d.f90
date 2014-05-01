      SUBROUTINE NETCDFJ_1D(NWAV,WAV,TEMP,DUMMY,DLEN,FNAME,FLEN,UNITNAME,UNITLEN, &
			    LAT,LON,DATE,MON,YEAR,NSZA,SZA,NVZA,VZA,NPHI,PHI,NLAY, &
			    HGT)
   
!=============================================================================
!   NWAV    Number of wavelengths
!   TEMP    Variable to be output to NETCDF file
!   DUMMY   Description of variable
!=============================================================================

      IMPLICIT NONE

      INCLUDE '../netcdf/netcdf.inc'

      INTEGER :: NWAV, NSZA, NVZA, NPHI, NLAY
      INTEGER :: ICNTN, NCID, RCODE, TIMESTEP
      INTEGER :: LATID, LONID, DATEID, MONID, YEARID, SZAID, VZAID, PHIID
      INTEGER :: HGTID, WAVID, TEMPID
      INTEGER :: LATDIM, LONDIM, DATEDIM, MONDIM, YEARDIM, SZADIM, VZADIM, PHIDIM
      INTEGER :: HGTDIM, WAVDIM
      INTEGER :: DLEN
      INTEGER :: FLEN
      INTEGER :: UNITLEN

      LOGICAL :: NCEXISTS

      CHARACTER(LEN=DLEN) :: DUMMY
      CHARACTER(LEN=FLEN) :: FNAME
      CHARACTER(LEN=UNITLEN) :: UNITNAME

! The number of dimensions we will have: Wav

      INTEGER, PARAMETER :: NDIMS = 1

! A vector containing all of the spacetime coordinates:

      INTEGER, DIMENSION(NDIMS) :: SPACETIME

! Dummy variable that holds the valid_min and valid_max values for variables

      INTEGER, DIMENSION(2) :: MINMAX

!=============================================================================
! Vector inputs which are allocated as AIPTs, and passed in by reference
! to the subroutine as arguments.  I took these definitions from the file
! allocation.f

      DOUBLE PRECISION, DIMENSION(1), INTENT(IN)      :: LAT
      DOUBLE PRECISION, DIMENSION(1), INTENT(IN)      :: LON
      INTEGER, DIMENSION(1), INTENT(IN)               :: DATE
      INTEGER, DIMENSION(1), INTENT(IN)               :: MON
      INTEGER, DIMENSION(1), INTENT(IN)               :: YEAR
      DOUBLE PRECISION, DIMENSION(NSZA), INTENT(IN)   :: SZA
      DOUBLE PRECISION, DIMENSION(NVZA), INTENT(IN)   :: VZA
      DOUBLE PRECISION, DIMENSION(NPHI), INTENT(IN)   :: PHI
      DOUBLE PRECISION, DIMENSION(0:NLAY), INTENT(IN) :: HGT
      DOUBLE PRECISION, DIMENSION(NWAV), INTENT(IN)   :: WAV
      DOUBLE PRECISION, DIMENSION(NWAV), INTENT(IN)   :: TEMP

! For the variables that vary in time and 3D space, we need a generic "START"
! and "COUNT" variable - beginning at the first entry in a data slice, and
! going for a length of one data slice

      INTEGER, DIMENSION(NDIMS) :: NDIMSTART
      INTEGER, DIMENSION(NDIMS) :: NDIMCOUNT

      DATA ICNTN /0/

!===================***********************************=======================
!===================*  Okay, let's get going already! *=======================
!===================***********************************=======================

! Increment counter.

      ICNTN = ICNTN+1

!=============================================================================
! First we will try and open the NetCDF file - IF the NCOPN fails, then we
! will go ahead and create a new one, and define all the dimensions and
! variables, and initialize the non-time-varying variables
! NOTE: the filename we're trying to open should not be hardcoded.  It needs
! to be dynamically read from *somewhere*, I just don't know where.
!=============================================================================

      INQUIRE (FILE=TRIM(FNAME),EXIST=NCEXISTS)

      IF((NCEXISTS).AND.(ICNTN.GT.1)) THEN
        NCID = NCOPN(TRIM(FNAME), NCWRITE, RCODE)

!=============================================================================
! If we're not defining the NetCDF structure this time through - because it
! was pre-existing, we need to read some information from the file first,
! so we know where to put the data.
! NCDID gets a dimension ID
! NCVID gets a variable ID
!=============================================================================

        TIMESTEP = TIMESTEP+1

!=============================================================================
! These variables are output for every run
!=============================================================================

        LATID       = NCVID (NCID, 'Latitude',  RCODE)
        LONID       = NCVID (NCID, 'Longitude',  RCODE)
        DATEID      = NCVID (NCID, 'Date',  RCODE)
        MONID       = NCVID (NCID, 'Month',  RCODE)
        YEARID      = NCVID (NCID, 'Year',  RCODE)
        SZAID      = NCVID (NCID, 'SolarZenithAngle',  RCODE)
        VZAID      = NCVID (NCID, 'ViewingZenithAngle',  RCODE)
        PHIID      = NCVID (NCID, 'RelativeAzimuthAngle',  RCODE)
        HGTID      = NCVID (NCID, 'Altitude',  RCODE)
        WAVID       = NCVID (NCID, 'Wavelength',  RCODE)
        TEMPID      = NCVID (NCID, TRIM(DUMMY),   RCODE)

!*****

      ELSE

!=============================================================================
! This is the first time we've accessed the file - the timestep is 1
!=============================================================================

        TIMESTEP = 1

!=============================================================================
! Create a NetCDF structure, and point it at the output file.
! Store the filehandle (unit number) in the NCID variable.
! NOTE: In reality, we should have a variable somewhere that describes the
! filename we are outputting to... but I don't know it off the top of my head
! NOTE: When we're going to be appending to the file, we don't want this to
! be NCCLOB.  But for now we're just trying to get one timestep working -
! the *first* one :)  We'd also put this inside a conditional and check to see
! if the file was already around - and if so just open it up and start writing
! instead of bothering to set everything up again.
!  If apppending, use NCNOCLOB
!  If not appending, use NCCLOB
!  cgm - 7/12/2001 - if we are in this section, we are either writing to
!  a new file .NOT.(NCEXISTS) or we are overwriting a file (NETCDF.NE.-1)
!  and not appending.  Therefore, set the flag for NCCRE to NCCLOB.
!=============================================================================

        NCID = NCCRE(TRIM(FNAME), NCCLOB, RCODE)
        IF (RCODE.EQ.-1) THEN
          WRITE(0,*) ' ERROR in NETCDF_OUT: NCCRE failed'
          WRITE(6,*) ' ERROR in NETCDF_OUT: NCCRE failed'
          RETURN
        END IF

! Create the dimensions of the dataset:

        LATDIM  = NCDDEF (NCID, 'Latitude', 1, RCODE)
        LONDIM  = NCDDEF (NCID, 'Longitude', 1, RCODE)
        DATEDIM = NCDDEF (NCID, 'Date', 1, RCODE)
        MONDIM  = NCDDEF (NCID, 'Month', 1, RCODE)
        YEARDIM = NCDDEF (NCID, 'Year', 1, RCODE)
        SZADIM  = NCDDEF (NCID, 'SolarZenithAngle', NSZA, RCODE)
        VZADIM  = NCDDEF (NCID, 'ViewingZenithAngle', NVZA, RCODE)
        PHIDIM  = NCDDEF (NCID, 'RelativeAzimuthAngle', NPHI, RCODE)
        HGTDIM  = NCDDEF (NCID, 'Altitude', NLAY+1, RCODE)
        WAVDIM  = NCDDEF (NCID,  'Wavelength',  NWAV,  RCODE)

!=============================================================================
! Create the coordinate (aka independent) variables:
!=============================================================================

        LATID  = NCVDEF (NCID, 'Latitude', NCDOUBLE, 1, LATDIM, &
                         RCODE)
        LONID  = NCVDEF (NCID, 'Longitude', NCDOUBLE, 1, LONDIM, &
                         RCODE)  
        DATEID = NCVDEF (NCID, 'Date', NCLONG, 1, DATEDIM, &
                         RCODE)
        MONID  = NCVDEF (NCID, 'Month', NCLONG, 1, MONDIM, &
                         RCODE)
        YEARID = NCVDEF (NCID, 'Year', NCLONG, 1, YEARDIM, &
                         RCODE)
        SZAID  = NCVDEF (NCID, 'SolarZenithAngle', NCDOUBLE, 1, SZADIM,  RCODE)
        VZAID  = NCVDEF (NCID, 'ViewingZenithAngle', NCDOUBLE, 1, VZADIM, RCODE)
        PHIID  = NCVDEF (NCID, 'RelativeAzimuthAngle', NCDOUBLE, 1, PHIDIM, RCODE)
        HGTID  = NCVDEF (NCID, 'Altitude', NCDOUBLE, 1, HGTDIM,  RCODE)
        WAVID  = NCVDEF (NCID, 'Wavelength', NCDOUBLE, 1, WAVDIM,  &
     		         RCODE)

!=============================================================================
! Create a vector containing the often referred to spacetime coordinates:
!=============================================================================

        SPACETIME(1) = WAVDIM

!=============================================================================
! Create the dependent variables:
!=============================================================================

! Radiance/Surface Albedo Jacobian:

        TEMPID=NCVDEF(NCID,TRIM(DUMMY),NCDOUBLE,NDIMS,SPACETIME,RCODE)

!=============================================================================
! Assign attributes (meta-data) to the various variables:
!=============================================================================

        CALL NCAPTC (NCID, NCGLOBAL, 'title', NCCHAR, DLEN, &
                     TRIM(DUMMY), RCODE)

        CALL NCAPTC (NCID, LATID, 'units', NCCHAR, 7, 'Degrees', &
                     RCODE)
        CALL NCAPTC (NCID, LONID, 'units', NCCHAR, 7, 'Degrees', &
                     RCODE)
        CALL NCAPTC (NCID, DATEID, 'units', NCCHAR, 8, 'Unitless', &
                     RCODE)
        CALL NCAPTC (NCID, MONID, 'units', NCCHAR, 8, 'Unitless', &
                     RCODE)
        CALL NCAPTC (NCID, YEARID, 'units', NCCHAR, 8, 'Unitless', &
                     RCODE)
        CALL NCAPTC (NCID, SZAID, 'units', NCCHAR, 7, 'Degrees', &
                     RCODE)
        CALL NCAPTC (NCID, VZAID, 'units', NCCHAR, 7, 'Degrees', &
                     RCODE)
        CALL NCAPTC (NCID, PHIID, 'units', NCCHAR, 7, 'Degrees', &
                     RCODE)
        CALL NCAPTC (NCID, HGTID, 'units', NCCHAR, 2, 'km', &
                     RCODE)

        CALL NCAPTC (NCID, TEMPID, 'units', NCCHAR, UNITLEN, TRIM(UNITNAME), &
                     RCODE)
        MINMAX = (/MINVAL(TEMP), MAXVAL(TEMP)/)
        CALL NCAPT (NCID, WAVID,'valid_range', NCLONG, 2, MINMAX, RCODE)

        CALL NCAPTC(NCID,WAVID,'units',NCCHAR,2,'nm',RCODE)
        MINMAX = (/MINVAL(WAV), MAXVAL(WAV)/)
        CALL NCAPT (NCID, WAVID,'valid_range', NCLONG, 2, MINMAX, RCODE)

!=============================================================================

!=============================================================================
! Get out of 'define' mode, and into 'data' mode
!=============================================================================

        CALL NCENDF (NCID, RCODE)

!=============================================================================
! Fill in the values of the non time varying variables.
! Remember, we're still in "initialization" here - this is only done the
! first time through.
!=============================================================================

        CALL NCVPT (NCID, LATID, 1, 1, LAT, RCODE)
        CALL NCVPT (NCID, LONID, 1, 1, LON, RCODE)
        CALL NCVPT (NCID, DATEID, 1, 1, DATE, RCODE)
        CALL NCVPT (NCID, MONID, 1, 1, MON, RCODE)
        CALL NCVPT (NCID, YEARID, 1, 1, YEAR, RCODE)
        CALL NCVPT (NCID, SZAID, 1, NSZA, SZA, RCODE)
        CALL NCVPT (NCID, VZAID, 1, NVZA, VZA, RCODE)
        CALL NCVPT (NCID, PHIID, 1, NPHI, PHI, RCODE)
        CALL NCVPT (NCID, HGTID, 1, NLAY+1, HGT, RCODE)
        CALL NCVPT (NCID, WAVID, 1, NWAV, WAV, RCODE)

      END IF

!=============================================================================
!====*********************************************************************====
!====*        From here on out, the code gets executed every time        *====
!====*********************************************************************====
!=============================================================================
!=============================================================================
! Define the START and COUNT arrays for each of the array variables.
!=============================================================================

      NDIMSTART = (/ 1 /)
      NDIMCOUNT = (/ NWAV /)

!=============================================================================
! Fill in the values for other variables
!=============================================================================

      CALL NCVPT (NCID, TEMPID, NDIMSTART, NDIMCOUNT, TEMP, RCODE)

!==============================================================================
! CLOSE the NetCDF file
!==============================================================================

      CALL NCCLOS (NCID, RCODE)

      return
      end
