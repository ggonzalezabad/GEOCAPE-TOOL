
! ###############################################################
! #                                                             #
! #                    THE VECTOR LIDORT MODEL                  #
! #                                                             #
! #  (Vector LInearized Discrete Ordinate Radiative Transfer)   #
! #   -      --         -        -        -         -           #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7                             #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   May 2014       (2.7)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS          (2.4)            #
! #       NEW: BPDF Land-surface KERNELS       (2.4R)           #
! #       NEW: Thermal Emission Treatment      (2.4RT)          #
! #       Consolidated BRDF treatment          (2.4RTC)         #
! #       f77/f90 Release                      (2.5)            #
! #       External SS / New I/O Structures     (2.6)            #
! #                                                             #
! #       Surface-leaving, BRDF Albedo-scaling (2.7)            # 
! #       Taylor series, Black-body Jacobians  (2.7)            #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

module VBRDF_Sup_Inputs_def

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!     New OG inputs are :
!       Observation-Geometry New dimensioning.    MAX_USER_OBSGEOMS
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control
!     Added Overall Exact flag for better control

!  This module contains the following structures:

!  VBRDF_Sup_Inputs - Intent(In) for VBRDFSup

use VLIDORT_PARS

implicit none

! #####################################################################
! #####################################################################

type VBRDF_Sup_Inputs

!  Top level flags (same as VLIDORT)
!  --------------------------------

!  Stream angle flag

      LOGICAL :: BS_DO_USER_STREAMS

!  BRDF surface flag

      LOGICAL :: BS_DO_BRDF_SURFACE

!  Surface emission

      LOGICAL :: BS_DO_SURFACE_EMISSION

!   !@@ Solar sources + Observational Geometry flag !@@

      LOGICAL :: BS_DO_SOLAR_SOURCES
      LOGICAL :: BS_DO_USER_OBSGEOMS

!  Numbers and Geometry (same as VLIDORT)
!  --------------------------------------

!  Number of Stokes components

      INTEGER :: BS_NSTOKES

!  Number of discrete ordinate streams

      INTEGER :: BS_NSTREAMS

!  number of solar beams to be processed

      INTEGER :: BS_NBEAMS

!  Bottom-of-atmosphere solar zenith angles, DEGREES

      REAL(fpk), dimension (MAXBEAMS) :: BS_BEAM_SZAS

!  user-defined relative azimuths

      INTEGER                                 :: BS_N_USER_RELAZMS
      REAL(fpk), dimension (MAX_USER_RELAZMS) :: BS_USER_RELAZMS

!  User-defined zenith angle input 

      INTEGER                                 :: BS_N_USER_STREAMS
      REAL(fpk), dimension (MAX_USER_STREAMS) :: BS_USER_ANGLES_INPUT

!  !@@ Observational geometry inputs

      INTEGER                                    :: BS_N_USER_OBSGEOMS
      REAL(fpk), dimension (MAX_USER_OBSGEOMS,3) :: BS_USER_OBSGEOMS

!  BRDF-specific inputs
!  --------------------

!   Number and index-list of bidirectional functions

      INTEGER                                            :: BS_N_BRDF_KERNELS
      CHARACTER (LEN=10), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_NAMES
      INTEGER, dimension ( MAX_BRDF_KERNELS )            :: BS_WHICH_BRDF

!  Parameters required for Kernel families

      INTEGER  , dimension ( MAX_BRDF_KERNELS ) :: &
        BS_N_BRDF_PARAMETERS
      REAL(fpk), dimension ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS ) :: &
        BS_BRDF_PARAMETERS

!  Lambertian Surface control

      LOGICAL, dimension ( MAX_BRDF_KERNELS )   :: BS_LAMBERTIAN_KERNEL_FLAG

!  Input kernel amplitude factors

      REAL(fpk), dimension ( MAX_BRDF_KERNELS ) :: BS_BRDF_FACTORS

!  Number of azimuth quadrature streams for BRDF

      INTEGER :: BS_NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL :: BS_DO_SHADOW_EFFECT

!  Exact flag (!@@) and Exact only flag --> no Fourier term calculations

      LOGICAL :: BS_DO_EXACT
      LOGICAL :: BS_DO_EXACTONLY

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

      LOGICAL   :: BS_DO_WSA_SCALING
      LOGICAL   :: BS_DO_BSA_SCALING
      REAL(fpk) :: BS_WSA_VALUE, BS_BSA_VALUE

!  Multiple-scattering Glitter options
!  -----------------------------------

!  Multiple reflectance corrections for GLITTER kernels (All of them!)

      LOGICAL :: BS_DO_GLITTER_MSRCORR

!  Multiple reflectance correction for exact-term Glitter kernels only

      LOGICAL :: BS_DO_GLITTER_MSRCORR_EXACTONLY

!  Correction order for the Multiple reflectance computations
!    ( = 0, no correction, 1, 2, 3   ec.)
!  Warning; using S > 0 can increase CPU dramatically

      INTEGER :: BS_GLITTER_MSRCORR_ORDER

!  Quadrature orders for MSRCORRemacs

      INTEGER :: BS_GLITTER_MSRCORR_NMUQUAD
      INTEGER :: BS_GLITTER_MSRCORR_NPHIQUAD

end type VBRDF_Sup_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VBRDF_Sup_Inputs

end module VBRDF_Sup_Inputs_def

