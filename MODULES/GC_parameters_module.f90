!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_parameters_module.f90
!
! !DESCRIPTION: This module contains dimensioning parameters,
! angle conversion and PI definition.
!\\
!\\
! !INTERFACE: 
!
MODULE GC_parameters_module
!
  IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
!
!  Dimensioning
!  ------------
! ----------------------------------------------------------------
! Dimensions with 'GC_' prefix, distinguish from VLIDORT variables
! ----------------------------------------------------------------
  INTEGER, PARAMETER :: GC_maxlayers      = 51
  INTEGER, PARAMETER :: GC_maxgeometries  = 2
  INTEGER, PARAMETER :: GC_maxuserlevels  = 2

! ----------------------------  
! Wavelengths, gases, messages
! ----------------------------
  INTEGER, PARAMETER :: maxgases      = 10
  INTEGER, PARAMETER :: maxaer        = 6
  INTEGER, PARAMETER :: maxcld        = 3
  INTEGER, PARAMETER :: maxflambdas   = 62001
  INTEGER, PARAMETER :: maxlambdas    = 62001
  INTEGER, PARAMETER :: maxmessages   = 100
  INTEGER, PARAMETER :: maxmoms       = 32
  INTEGER, PARAMETER :: maxgksec      = 6, maxgkmatc = 8
  INTEGER, PARAMETER :: maxscatter    = 3  ! Molecules, aerosols, clouds
  INTEGER, DIMENSION(maxgkmatc), PARAMETER :: &
           greekmat_idxs = (/1, 2, 5, 6, 11, 12, 15, 16/), &
           phasmoms_idxs = (/1, 5, 5, 2, 3, 6, 6, 4/)

! --------------
! Some constants
! --------------
  REAL(KIND=8), PARAMETER :: pi      = 3.14159265358979d0
  REAL(KIND=8), PARAMETER :: deg2rad = pi / 180.d0,    &
                             rad2deg = 180.d0 / pi
  INTEGER,      PARAMETER :: max_ch_len = 256

! ------------------------------------
! File units for different input files
! ------------------------------------
  INTEGER, PARAMETER :: ctrunit = 11
  INTEGER, PARAMETER :: aerunit = 12
  INTEGER, PARAMETER :: cldunit = 13
  INTEGER, PARAMETER :: errunit = 15

! ------------------
! Profile PARAMETERS
! ------------------
  INTEGER, PARAMETER         :: pmaxl = GC_maxlayers+1, pmaxc = 30
!
! !REVISION HISTORY:
!  April 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!------------------------------------------------------------------------------
!EOC

END MODULE GC_parameters_module
