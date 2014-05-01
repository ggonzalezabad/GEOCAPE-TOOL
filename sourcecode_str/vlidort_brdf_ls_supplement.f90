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
! #  Address :      RT Solutions, inc.                          #
! #            9 Channing Street                                #
! #             Cambridge, MA 02138, USA                        #
! #            Tel: (617) 492 1183                              #
! #                                                             #
! #  Email :      rtsolutions@verizon.net                       #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4RT, 2.5            #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.5)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       THREADED/OPTIMIZED F90 code         (2.5)             #
! #                                                             #
! ###############################################################

!    #####################################################
!    #                                                   #
!    #   This Version of VLIDORT comes with a GNU-style  #
!    #   license. Please read the license carefully.     #
!    #                                                   #
!    #####################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VLIDORT_BRDF_INPUTS_PLUS                         #
! #                                                             #
! #            VLIDORT_BRDF_LS_MASTER (master), calling         #
! #                                                             #
! #              VLIDORT_BRDF_MAKER_PLUS                        #
! #              VLIDORT_GCMCRI_MAKER_PLUS                      #
! #              BRDF_QUADRATURE_Gaussian                       #
! #              BRDF_QUADRATURE_Trapezoid (not used)           #
! #              VLIDORT_BRDF_LS_FOURIER                        #
! #                                                             #
! ###############################################################

SUBROUTINE VLIDORT_BRDF_LS_MASTER                                          &
        ( DO_USER_STREAMS, DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,              & ! Inputs
          DO_SURFACE_EMISSION,  DO_DEBUG_RESTORATION, NSTOKES,             & ! Inputs
          NMOMENTS_INPUT, N_BRDF_KERNELS, WHICH_BRDF,                      & ! Inputs
          LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,             & ! Inputs
          N_BRDF_PARAMETERS, BRDF_PARAMETERS, BRDF_NAMES,                  & ! Inputs
          DO_KERNEL_PARAMS_WFS, DO_KERNEL_FACTOR_WFS, DO_KPARAMS_DERIVS,   & ! Inputs
          N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS, N_SURFACE_WFS,         & ! Inputs
          NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,                & ! Inputs
          BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,                      & ! Inputs
          BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,                    & ! Outputs
          LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,        & ! Outputs
          EXACTDB_BRDFUNC,  LS_EXACTDB_BRDFUNC,                            & ! Outputs
          EMISSIVITY, USER_EMISSIVITY, LS_EMISSIVITY, LS_USER_EMISSIVITY )   ! Outputs

!  Prepares (linearizations of) the bidirectional reflectance functions
!  necessary for VLIDORT.

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Input arguments
!  ===============

!  User stream Control

      LOGICAL, intent(in) :: DO_USER_STREAMS

!  Surface emission

      LOGICAL, intent(in) :: DO_SURFACE_EMISSION

!  Debug flag for restoration

      LOGICAL, intent(in) :: DO_DEBUG_RESTORATION

!  number of Stokes components

      INTEGER, intent(in) :: NSTOKES

!  Input number of moments (only used for restoration debug)

      INTEGER, intent(in) :: NMOMENTS_INPUT

!   Number and index-list of bidirectional functions

      INTEGER, intent(in) :: N_BRDF_KERNELS
      INTEGER, intent(in) :: WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER     , intent(inout) :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(kind=8), intent(inout) :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF names

      CHARACTER*10, intent(in)    :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Lambertian Surface control

      LOGICAL, intent(in) :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(kind=8), intent(in) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER, intent(in) :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, intent(in) :: DO_SHADOW_EFFECT

!  Multiple reflectance correction for direct beam flag
!              (only for GLITTER type kernels)

      LOGICAL, intent(in) :: DO_COXMUNK_DBMS

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL, intent(in) :: DO_KERNEL_FACTOR_WFS ( MAX_BRDF_KERNELS )
      LOGICAL, intent(in) :: DO_KERNEL_PARAMS_WFS ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, intent(in) :: DO_KPARAMS_DERIVS  ( MAX_BRDF_KERNELS )

!  number of surfaceweighting functions

      INTEGER, intent(inout) :: N_SURFACE_WFS
      INTEGER, intent(in)    :: N_KERNEL_FACTOR_WFS
      INTEGER, intent(in)    :: N_KERNEL_PARAMS_WFS

!  Local angle control

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: N_USER_RELAZMS

!  Angles

      REAL(kind=8), intent(in) ::  BEAM_SZAS         (MAXBEAMS)
      REAL(kind=8), intent(in) ::  USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Output arguments
!  ================

!  Exact (direct bounce) BRDF (same all threads)
    
      REAL(kind=8), intent(out) :: EXACTDB_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Fourier components of BRDF, in the following order (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(kind=8), intent(out) :: BRDF_F_0 ( MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAXSTREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: BRDF_F   ( MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: USER_BRDF_F_0 ( MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: USER_BRDF_F   ( MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAX_USER_STREAMS, MAXSTREAMS )

!  Linearized Exact (direct bounce) BRDF (same all threads)
    
      REAL(kind=8), intent(out) :: LS_EXACTDB_BRDFUNC &
          ( MAX_SURFACEWFS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Linearized Fourier components of BRDF, (same all threads)

!    incident solar directions,   reflected quadrature streams
!    incident quadrature streams, reflected quadrature streams
!    incident solar directions,   reflected user streams
!    incident quadrature streams, reflected user streams

      REAL(kind=8), intent(out) :: LS_BRDF_F_0 ( MAX_SURFACEWFS, MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAXSTREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: LS_BRDF_F   ( MAX_SURFACEWFS, MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: LS_USER_BRDF_F_0 &
                                   ( MAX_SURFACEWFS, MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAX_USER_STREAMS, MAXBEAMS )
      REAL(kind=8), intent(out) :: LS_USER_BRDF_F   &
                                   ( MAX_SURFACEWFS, MAXSTOKES_SQ, 0:MAXMOMENTS_INPUT, MAX_USER_STREAMS, MAXSTREAMS )

!  Fourier components of emissivity

      REAL(kind=8), intent(out) :: USER_EMISSIVITY(MAXSTOKES,MAX_USER_STREAMS)
      REAL(kind=8), intent(out) :: EMISSIVITY     (MAXSTOKES,MAXSTREAMS)


!  Fourier components of emissivity

      REAL(kind=8), intent(out) :: LS_USER_EMISSIVITY ( MAX_SURFACEWFS,MAXSTOKES,MAX_USER_STREAMS)
      REAL(kind=8), intent(out) :: LS_EMISSIVITY      ( MAX_SURFACEWFS,MAXSTOKES,MAXSTREAMS      )

!  BRDF functions
!  --------------

!  BRDF External functions
!  =======================

!  lambertian

      EXTERNAL       LAMBERTIAN_VFUNCTION

!  Modis-type kernels

      EXTERNAL       ROSSTHIN_VFUNCTION
      EXTERNAL       ROSSTHICK_VFUNCTION
      EXTERNAL       LISPARSE_VFUNCTION
      EXTERNAL       LIDENSE_VFUNCTION
      EXTERNAL       HAPKE_VFUNCTION
      EXTERNAL       ROUJEAN_VFUNCTION
      EXTERNAL       RAHMAN_VFUNCTION

!  Cox-munk types

      EXTERNAL       COXMUNK_VFUNCTION
      EXTERNAL       COXMUNK_VFUNCTION_DB
      EXTERNAL       GISSCOXMUNK_VFUNCTION
      EXTERNAL       GISSCOXMUNK_VFUNCTION_DB

!  GCM CRI is not an external call
!      EXTERNAL       GCMCRI_VFUNCTION
!      EXTERNAL       GCMCRI_VFUNCTION_DB

!  Two new Land BRDFs, introduced 1 July 2008

      EXTERNAL       RHERMAN_VFUNCTION
      EXTERNAL       BREON_VFUNCTION

!  new for Version 2.4R, introduced 30 April 2009, 6 May 2009
!    2008 Veg/Soil functions based on Breon work 2008 as supplied to OCO
!    2009 function is final Kernel supplied by Breon, May 5 2009.

      EXTERNAL       BPDF2008VEG_VFUNCTION
      EXTERNAL       BPDF2008SOIL_VFUNCTION
      EXTERNAL       BPDF2009_VFUNCTION

!  Hapke old uses exact DISORT code
!      EXTERNAL       HAPKE_VFUNCTION_OLD

!  BRDFs with derivatives

      EXTERNAL       LISPARSE_VFUNCTION_PLUS
      EXTERNAL       LIDENSE_VFUNCTION_PLUS
      EXTERNAL       HAPKE_VFUNCTION_PLUS
      EXTERNAL       RAHMAN_VFUNCTION_PLUS

      EXTERNAL       COXMUNK_VFUNCTION_PLUS
      EXTERNAL       COXMUNK_VFUNCTION_DB_PLUS
      EXTERNAL       GISSCOXMUNK_VFUNCTION_PLUS
      EXTERNAL       GISSCOXMUNK_VFUNCTION_DB_PLUS

!  OTHERS


!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: BRDFUNC   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: BRDFUNC_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8)  :: USER_BRDFUNC   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: USER_BRDFUNC_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  DB Kernel values

      REAL(kind=8)  :: DBKERNEL_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8)  :: EBRDFUNC      ( MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: USER_EBRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8)  :: D_USER_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_USER_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Linearized Exact DB values

      REAL(kind=8)  :: D_DBKERNEL_BRDFUNC ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8)  :: D_EBRDFUNC      ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8)  :: D_USER_EBRDFUNC ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      REAL(kind=8)  :: PHIANG(MAX_USER_RELAZMS)
      REAL(kind=8)  :: COSPHI(MAX_USER_RELAZMS)
      REAL(kind=8)  :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      REAL(kind=8)  :: SZASURCOS(MAXBEAMS)
      REAL(kind=8)  :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates

      REAL(kind=8)  :: QUAD_STREAMS(MAXSTREAMS)
      REAL(kind=8)  :: QUAD_WEIGHTS(MAXSTREAMS)
      REAL(kind=8)  :: QUAD_SINES  (MAXSTREAMS)
      REAL(kind=8)  :: QUAD_STRMWTS(MAXSTREAMS)

!  Viewing zenith streams

      REAL(kind=8)  :: USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=8)  :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      INTEGER       :: NBRDF_HALF
      REAL(kind=8)  :: X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8)  :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations

      REAL(kind=8)  :: BAX_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8)  :: CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8)  :: SXE_BRDF ( MAXSTHALF_BRDF )

!  Azimuth factors

      REAL(kind=8)  :: BRDF_COSAZMFAC(MAXSTREAMS_BRDF)
      REAL(kind=8)  :: BRDF_SINAZMFAC(MAXSTREAMS_BRDF)

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: LOCAL_BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8)  :: LOCAL_BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8)  :: LOCAL_USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8)  :: LOCAL_USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8)  :: LOCAL_EMISSIVITY      ( MAXSTOKES, MAXSTREAMS       )
      REAL(kind=8)  :: LOCAL_USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  Local Derivative-kernel Fourier components
!  ==========================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8)  :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8)  :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8)  :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8)  :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8)  :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAXSTREAMS       )
      REAL(kind=8)  :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAX_USER_STREAMS )

!  Other local variables
!  =====================

!  Spherical albedo

      REAL(kind=8)       :: SPHERICAL_ALBEDO(MAX_BRDF_KERNELS)

!  help

      INTEGER            :: K, W, P, I, J, IB, UI, UM, IA, M, O1, Q, NSTOKESSQ
      INTEGER            :: WOFFSET ( MAX_BRDF_KERNELS)
      INTEGER            :: BRDF_NPARS, NMOMENTS
      REAL(kind=8)       :: PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL            :: DERIVS ( MAX_BRDF_PARAMETERS )
      REAL(kind=8)       :: MUX, DELFAC, HELP_A, SUM, FF, ARGUMENT
      LOGICAL            :: ADD_FOURIER

!  default, use Gaussian quadrature

      LOGICAL, parameter :: DO_BRDFQUAD_GAUSSIAN = .true.

!  Set up Quadrature streams

      CALL GAULEG(0.0d0,1.0d0, QUAD_STREAMS, QUAD_WEIGHTS, NSTREAMS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = DSQRT(1.0d0-QUAD_STREAMS(I)*QUAD_STREAMS(I))
        QUAD_STRMWTS(I) = QUAD_STREAMS(I) * QUAD_WEIGHTS(I)
      enddo

!  Number of Stokes components squared
!    ** Bookkeeping for surface kernel Cox-Munk types
!    ** Only the Giss CoxMunk kernel is vectorized (as of 19 January 2006)

!   Additional code for complex RI Giss Cox-Munk, 15 march 2010.
!     3 parameters are PARS(1) = sigma_sq
!                      PARS(2) = Real (RI)
!                      PARS(3) = Imag (RI)

      NSTOKESSQ  = 1
      DO K = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(K) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(K) .EQ. 'GissCoxMnk' ) THEN
            N_BRDF_PARAMETERS(K) = 3
            IF ( DO_SHADOW_EFFECT ) THEN
              BRDF_PARAMETERS(K,3) = ONE
            ELSE
              BRDF_PARAMETERS(K,3) = ZERO
            ENDIF
         ELSE IF ( BRDF_NAMES(K) .EQ. 'GCMcomplex' ) THEN
            N_BRDF_PARAMETERS(K) = 3
         ENDIF
          IF ( BRDF_NAMES(K) .EQ. 'GissCoxMnk'  .OR. &
               BRDF_NAMES(K) .EQ. 'GCMcomplex' ) THEN
            NSTOKESSQ = NSTOKES * NSTOKES
         ENDIF
      ENDDO

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams
!    Warning, this shoudl be the BOA angle. OK for the non-refractive case.
!        
      DO IB = 1, NBEAMS
        MUX =  DCOS(BEAM_SZAS(IB)*DEG_TO_RAD)
        SZASURCOS(IB) = MUX
        SZASURSIN(IB) = DSQRT(ONE-MUX*MUX)
      ENDDO

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = DCOS(USER_ANGLES_INPUT(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = DSQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

      DO IA = 1, N_USER_RELAZMS
        PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
        COSPHI(IA) = DCOS(PHIANG(IA))
        SINPHI(IA) = DSIN(PHIANG(IA))
      ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian                                          &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid                                         &
        ( DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,                      & ! inputs
          X_BRDF, CX_BRDF, SX_BRDF, A_BRDF, BAX_BRDF, CXE_BRDF, SXE_BRDF )       ! Outputs
      ENDIF

!  Number of weighting functions, and offset

      W = 0
      WOFFSET(1) = 0
      DO K = 1, N_BRDF_KERNELS
        IF ( DO_KERNEL_FACTOR_WFS(K) ) W = W + 1
        DO P = 1, N_BRDF_PARAMETERS(K)
          IF ( DO_KERNEL_PARAMS_WFS(K,P) ) W = W + 1
        ENDDO
        IF ( K.LT.N_BRDF_KERNELS ) WOFFSET(K+1) = W
      ENDDO

      N_SURFACE_WFS = N_KERNEL_FACTOR_WFS + N_KERNEL_PARAMS_WFS
      IF ( W .ne. N_SURFACE_WFS ) stop'bookkeeping wrong'

!  Initialise outputs
!  ------------------

!  Not necessary to initialize derivative outputs

!  zero the BRDF Fourier components

      DO O1 = 1, NSTOKESSQ
        DO M = 0, NMOMENTS
          DO I = 1, NSTREAMS
            DO IB = 1, NBEAMS
              BRDF_F_0(O1,M,I,IB) = ZERO
            ENDDO
            DO J = 1, NSTREAMS
              BRDF_F(O1,M,I,J) = ZERO
            ENDDO
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UM = 1, N_USER_STREAMS
              DO IB = 1, NBEAMS
                USER_BRDF_F_0(O1,M,UM,IB) = ZERO
              ENDDO
              DO J = 1, NSTREAMS
                  USER_BRDF_F(O1,M,UM,J) = ZERO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

!  Zero Exact Direct Beam BRDF

      DO O1 = 1, NSTOKESSQ
        DO IA = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UM = 1, N_USER_STREAMS
              EXACTDB_BRDFUNC(O1,UM,IA,IB) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Initialize surface emissivity

      IF ( DO_SURFACE_EMISSION ) THEN
        DO O1 = 1, NSTOKES
          DO I = 1, NSTREAMS
            EMISSIVITY(O1,I) = ONE
          ENDDO
          IF ( DO_USER_STREAMS ) THEN
            DO UI = 1, N_USER_STREAMS
              USER_EMISSIVITY(O1,UI) = ONE
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  Fill BRDF arrays
!  ----------------

      DO K = 1, N_BRDF_KERNELS

!  Copy parameter variables into local quantities

        BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO P = 1, MAX_BRDF_PARAMETERS
          PARS(P) = BRDF_PARAMETERS(K,P)
        ENDDO
        IF ( DO_KPARAMS_DERIVS(K) ) THEN
          DO P = 1, MAX_BRDF_PARAMETERS
            DERIVS(P) = DO_KERNEL_PARAMS_WFS(K,P)
          ENDDO
        ENDIF

!  Lambertian kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. LAMBERTIAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( LAMBERTIAN_VFUNCTION, LAMBERTIAN_VFUNCTION,           & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( ROSSTHIN_VFUNCTION, ROSSTHIN_VFUNCTION,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( ROSSTHICK_VFUNCTION, ROSSTHICK_VFUNCTION,             & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( LISPARSE_VFUNCTION_PLUS, LISPARSE_VFUNCTION_PLUS,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( LISPARSE_VFUNCTION, LISPARSE_VFUNCTION,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( LIDENSE_VFUNCTION_PLUS, LIDENSE_VFUNCTION_PLUS,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( LIDENSE_VFUNCTION, LIDENSE_VFUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( HAPKE_VFUNCTION_PLUS, HAPKE_VFUNCTION_PLUS,                     & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( HAPKE_VFUNCTION, HAPKE_VFUNCTION,                      & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( ROUJEAN_VFUNCTION, ROUJEAN_VFUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Rahman kernel: (3 free parameters) 

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( RAHMAN_VFUNCTION_PLUS, RAHMAN_VFUNCTION_PLUS,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( RAHMAN_VFUNCTION, RAHMAN_VFUNCTION,                   & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
        ENDIF

!  Scalar-only original Cox-Munk kernel: (2 free parameters, Shadow = Third).
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
         IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
         IF (  DO_COXMUNK_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( COXMUNK_VFUNCTION_PLUS, COXMUNK_VFUNCTION_DB_PLUS,              & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
         ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( COXMUNK_VFUNCTION, COXMUNK_VFUNCTION_DB,              & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( COXMUNK_VFUNCTION_PLUS, COXMUNK_VFUNCTION_PLUS,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( COXMUNK_VFUNCTION, COXMUNK_VFUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
         ENDIF
        ENDIF

!  GISS Vector Cox-Munk kernel: (2 free parameters, Shadow = Third). Real Refractive Index.
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX ) THEN
         IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
         IF (  DO_COXMUNK_DBMS ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( GISSCOXMUNK_VFUNCTION_PLUS, GISSCOXMUNK_VFUNCTION_DB_PLUS,      & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION_DB,      & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
         ENDIF
         ELSE
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_BRDF_MAKER_PLUS &
           ( GISSCOXMUNK_VFUNCTION_PLUS, GISSCOXMUNK_VFUNCTION_PLUS,         & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_BRDF_MAKER &
           ( GISSCOXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION,         & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
          ENDIF
         ENDIF
        ENDIF

!  New code . Giss Cox_munk with Complex RI : (3 Free parameters). Shadow is optional.

        IF ( WHICH_BRDF(K) .EQ. GCMCRI_IDX ) THEN
          IF ( DO_KPARAMS_DERIVS(K) ) THEN
            CALL VLIDORT_GCMCRI_MAKER_PLUS &
           ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,                              & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,     & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs
          ELSE
            CALL VLIDORT_GCMCRI_MAKER &
           ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,                    & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
         ENDIF
        ENDIF

!  Rondeaux-Herman (vegetation) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RHERMAN_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( RHERMAN_VFUNCTION, RHERMAN_VFUNCTION,                 & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Breon et al (desert) model (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BREON_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( BREON_VFUNCTION, BREON_VFUNCTION,                     & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  BPDF 2008 Vegetation kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2008VEG_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( BPDF2008VEG_VFUNCTION, BPDF2008VEG_VFUNCTION,         & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  BPDF 2008 Soil kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2008SOIL_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( BPDF2008SOIL_VFUNCTION, BPDF2008SOIL_VFUNCTION,       & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  BPDF 2009 kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. BPDF2009_IDX ) THEN
          CALL VLIDORT_BRDF_MAKER &
           ( BPDF2009_VFUNCTION, BPDF2009_VFUNCTION,               & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                 & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,     & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,     & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,   & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,   & ! Inputs
             X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,         & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,              & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC )    ! Outputs
        ENDIF

!  Exact BRDFUNC
!  ------------

!  factor

        FF = BRDF_FACTORS(K)

!  Compute Exact Direct Beam BRDF

        DO O1 = 1, NSTOKESSQ
          DO IA = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                EXACTDB_BRDFUNC(O1,UM,IA,IB) = EXACTDB_BRDFUNC(O1,UM,IA,IB) + FF * DBKERNEL_BRDFUNC(O1,UM,IA,IB)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Linearization w.r.t Kernel Factor

        W  = WOFFSET(K)
        IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
          W = W + 1
          DO O1 = 1, NSTOKESSQ
            DO IA = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  LS_EXACTDB_BRDFUNC(W,O1,UM,IA,IB) = DBKERNEL_BRDFUNC(O1,UM,IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Linearization w.r.t Kernel parameters

        DO P = 1, BRDF_NPARS
          IF ( DERIVS(P) ) THEN
            W = W + 1
            DO O1 = 1, NSTOKESSQ
              DO IA = 1, N_USER_RELAZMS
                DO IB = 1, NBEAMS
                  DO UM = 1, N_USER_STREAMS
                    LS_EXACTDB_BRDFUNC(W,O1,UM,IA,IB) = FF * D_DBKERNEL_BRDFUNC(P,O1,UM,IA,IB)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO

!  Fourier Work now
!  ================

        DO M = 0, NMOMENTS

!  Fourier addition flag

          ADD_FOURIER = ( .not.LAMBERTIAN_KERNEL_FLAG(K) .or. &
                          (LAMBERTIAN_KERNEL_FLAG(K).AND.M.EQ.0) )

!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            DELFAC   = ONE
            DO I = 1, NSTREAMS_BRDF
              BRDF_COSAZMFAC(I) = A_BRDF(I)
              BRDF_SINAZMFAC(I) = ZERO
            ENDDO
          ELSE
            DELFAC   = TWO
            DO I = 1, NSTREAMS_BRDF
              ARGUMENT = DBLE(M) * X_BRDF(I)
              BRDF_COSAZMFAC(I) = A_BRDF(I) * DCOS ( ARGUMENT )
              BRDF_SINAZMFAC(I) = A_BRDF(I) * DSIN ( ARGUMENT )
            ENDDO
          ENDIF

!  Call

          CALL VLIDORT_BRDF_FOURIER                                                      &
          ( DO_USER_STREAMS, DO_SURFACE_EMISSION, LAMBERTIAN_KERNEL_FLAG(K), M, NSTOKES, & ! Inputs
            NSTOKESSQ, NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,      & ! Inputs
            DELFAC, BRDF_FACTORS(K), BRDF_COSAZMFAC, BRDF_SINAZMFAC, A_BRDF, BAX_BRDF,   & ! Inputs
            BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,   & ! Inputs
            LOCAL_BRDF_F, LOCAL_BRDF_F_0, LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,        & ! Outputs
            LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                                      ! Outputs

!  Linear call

          IF ( BRDF_NPARS .gt. 0 ) then
            CALL VLIDORT_BRDF_LS_FOURIER                                                           &                
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION, LAMBERTIAN_KERNEL_FLAG(K), M, NSTOKES, NSTOKESSQ, & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,  BRDF_NPARS, DERIVS,       & ! Inputs
           DELFAC, BRDF_FACTORS(K), BRDF_COSAZMFAC, BRDF_SINAZMFAC, A_BRDF, BAX_BRDF,              & ! Inputs
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC,  & ! Inputs
           D_LOCAL_BRDF_F, D_LOCAL_BRDF_F_0, D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,           & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                                            ! Outputs
          ENDIF

!  Spherical albedo (debug only)

          IF ( M .EQ. 0 ) THEN
            Q = 1
            IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, NSTREAMS
               SUM = ZERO
               DO J = 1, NSTREAMS
                SUM = SUM + LOCAL_BRDF_F(Q,I,J) * QUAD_STRMWTS(J)
               ENDDO
               HELP_A = HELP_A + SUM * QUAD_STRMWTS(I)
              ENDDO
              SPHERICAL_ALBEDO(K) = HELP_A*FOUR
             ENDIF
          ENDIF

!  Start Fourier addition

          IF ( ADD_FOURIER ) THEN

!  Kernel combinations (for quadrature reflectance)
!  ------------------------------------------------

!  factor

            FF = BRDF_FACTORS(K)

!  Basic Kernel sum

            DO Q = 1, NSTOKESSQ
              DO I = 1, NSTREAMS
                DO IB = 1, NBEAMS
                  BRDF_F_0(Q,M,I,IB) = BRDF_F_0(Q,M,I,IB) + FF * LOCAL_BRDF_F_0(Q,I,IB)
                ENDDO
                DO J = 1, NSTREAMS
                  BRDF_F(Q,M,I,J) = BRDF_F(Q,M,I,J) + FF * LOCAL_BRDF_F(Q,I,J)
                ENDDO
              ENDDO
            ENDDO

!  Linearization w.r.t Kernel Factor

            W  = WOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              W = W + 1
              DO Q = 1, NSTOKESSQ
                DO I = 1, NSTREAMS
                  DO IB = 1, NBEAMS
                    LS_BRDF_F_0(W,Q,M,I,IB) = LOCAL_BRDF_F_0(Q,I,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    LS_BRDF_F(W,Q,M,I,J) =  LOCAL_BRDF_F(Q,I,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, BRDF_NPARS
              IF ( DERIVS(P) ) THEN
                W = W + 1
                DO Q = 1, NSTOKESSQ
                  DO I = 1, NSTREAMS
                    DO IB = 1, NBEAMS
                      LS_BRDF_F_0(W,Q,M,I,IB) = FF*D_LOCAL_BRDF_F_0(P,Q,I,IB)
                    ENDDO
                    DO J = 1, NSTREAMS
                      LS_BRDF_F(W,Q,M,I,J) = FF*D_LOCAL_BRDF_F(P,Q,I,J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Kernel combinations (for user-stream reflectance)
!  -------------------------------------------------

!  Basic kernel summation

            IF ( DO_USER_STREAMS ) THEN
              DO Q = 1, NSTOKESSQ
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    USER_BRDF_F_0(Q,M,UM,IB) = USER_BRDF_F_0(Q,M,UM,IB) + FF * LOCAL_USER_BRDF_F_0(Q,UM,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    USER_BRDF_F(Q,M,UM,J) = USER_BRDF_F(Q,M,UM,J) + FF * LOCAL_USER_BRDF_F(Q,UM,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel Factor

            W  = WOFFSET(K)
            IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
              W = W + 1
              DO Q = 1, NSTOKESSQ
                DO UM = 1, N_USER_STREAMS
                  DO IB = 1, NBEAMS
                    LS_USER_BRDF_F_0(W,Q,M,UM,IB) = LOCAL_USER_BRDF_F_0(Q,UM,IB)
                  ENDDO
                  DO J = 1, NSTREAMS
                    LS_USER_BRDF_F(W,Q,M,UM,J) =  LOCAL_USER_BRDF_F(Q,UM,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Linearization w.r.t Kernel parameters

            DO P = 1, BRDF_NPARS
              IF ( DERIVS(P) ) THEN
                W = W + 1
                DO Q = 1, NSTOKESSQ
                  DO UM = 1, N_USER_STREAMS
                    DO IB = 1, NBEAMS
                      LS_USER_BRDF_F_0(W,Q,M,UM,IB) = FF * D_LOCAL_USER_BRDF_F_0(P,Q,UM,IB)
                    ENDDO
                    DO J = 1, NSTREAMS
                      LS_USER_BRDF_F(W,Q,M,UM,J) = FF * D_LOCAL_USER_BRDF_F(P,Q,UM,J)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO

!  Total emissivities
!  ------------------

!  only if flagged

            IF ( DO_SURFACE_EMISSION.and. M.eq.0 ) THEN

!  Basci kernel contributions

              DO Q = 1, NSTOKES
                DO I = 1, NSTREAMS
                  EMISSIVITY(Q,I) = EMISSIVITY(Q,I) - LOCAL_EMISSIVITY(Q,I)
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  DO UI = 1, N_USER_STREAMS
                    USER_EMISSIVITY(Q,UI) = USER_EMISSIVITY(Q,UI) - LOCAL_USER_EMISSIVITY(Q,UI)
                  ENDDO
                ENDIF
              ENDDO

!  Linearization w.r.t Kernel Factor

              W  = WOFFSET(K)
              IF ( DO_KERNEL_FACTOR_WFS(K) ) THEN
                W = W + 1
                DO Q = 1, NSTOKES
                  DO I = 1, NSTREAMS
                    LS_EMISSIVITY(W,Q,I) = - LOCAL_EMISSIVITY(Q,I) / FF
                  ENDDO
                  IF ( DO_USER_STREAMS ) THEN
                    DO UI = 1, N_USER_STREAMS
                      LS_USER_EMISSIVITY(W,Q,UI) = - LOCAL_USER_EMISSIVITY(Q,UI) / FF
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF

!  Linearization w.r.t Kernel parameters

              DO P = 1, BRDF_NPARS
                IF ( DERIVS(P) ) THEN
                  W = W + 1
                  DO Q = 1, NSTOKES
                    DO I = 1, NSTREAMS
                      LS_EMISSIVITY(W,Q,I) = - D_LOCAL_EMISSIVITY(P,Q,I)
                    ENDDO
                    IF ( DO_USER_STREAMS ) THEN
                      DO UI = 1, N_USER_STREAMS
                        LS_USER_EMISSIVITY(W,Q,UI) = - D_LOCAL_USER_EMISSIVITY(P,Q,UI)
                      ENDDO
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

!  End emissivity clause

            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  End kernel loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE VLIDORT_BRDF_LS_MASTER

!

SUBROUTINE VLIDORT_BRDF_MAKER_PLUS                                           & 
           ( BRDF_VFUNCTION_PLUS, BRDF_VFUNCTION_DB_PLUS,                    & ! Inputs
             DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,               & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, CX_BRDF,  & ! Inputs
             SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,                      & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs

      implicit none

!  Prepares the bidirectional reflectance scatter matrices

!  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_VFUNCTION_PLUS
      EXTERNAL         BRDF_VFUNCTION_DB_PLUS

!  Local flags

      LOGICAL, intent(in) :: DO_USER_STREAMS
      LOGICAL, intent(in) :: DO_SURFACE_EMISSION

!  Number of Azimuth waudrature streams

      INTEGER, intent(in) :: NSTREAMS_BRDF
      INTEGER, intent(in) :: NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER, intent(in) :: NSTOKESSQ 

!  Local number of Kernel parameters

      INTEGER, intent(in) :: BRDF_NPARS
      
!  Local angle control

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: N_USER_RELAZMS

!  Local angles

      REAL(kind=8), intent(in) ::  PHIANG(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  COSPHI(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  SINPHI(MAX_USER_RELAZMS)

      REAL(kind=8), intent(in) ::  SZASURCOS(MAXBEAMS)
      REAL(kind=8), intent(in) ::  SZASURSIN(MAXBEAMS)

      REAL(kind=8), intent(in) ::  QUAD_STREAMS(MAXSTREAMS)
      REAL(kind=8), intent(in) ::  QUAD_SINES  (MAXSTREAMS)

      REAL(kind=8), intent(in) ::  USER_STREAMS(MAX_USER_STREAMS)
      REAL(kind=8), intent(in) ::  USER_SINES  (MAX_USER_STREAMS)

!  Local parameter array

      REAL(kind=8), intent(in) ::  PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL     , intent(in) ::  DERIVS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      REAL(kind=8), intent(in) ::  X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  CXE_BRDF ( MAXSTHALF_BRDF )
      REAL(kind=8), intent(in) ::  SXE_BRDF ( MAXSTHALF_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: BRDFUNC   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: BRDFUNC_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: USER_BRDFUNC   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_BRDFUNC_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: DBKERNEL_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: EBRDFUNC      ( MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_EBRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: D_USER_BRDFUNC     &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_BRDFUNC_0   &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: D_DBKERNEL_BRDFUNC  &
                                  ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: D_EBRDFUNC         &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_EBRDFUNC    &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER          :: W, I, UI, J, K, KE, IB, Q
      REAL(kind=8)     :: DFUNC ( 16, MAX_BRDF_PARAMETERS ), FUNC(16)

!  Exact DB calculation
!  --------------------

      DO K = 1, N_USER_RELAZMS
         DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              CALL BRDF_VFUNCTION_DB_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),            & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),           & ! Inputs
                 FUNC, DFUNC )                                                ! Output
              DO Q = 1, NSTOKESSQ
                DBKERNEL_BRDFUNC(Q,UI,K,IB) = FUNC(Q)
                DO W = 1, BRDF_NPARS
                  D_DBKERNEL_BRDFUNC(W,Q,UI,K,IB)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
         ENDDO
      ENDDO

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                 SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),             &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),          &
                 FUNC, DFUNC ) 
            DO Q = 1, NSTOKESSQ
              BRDFUNC_0(Q,I,IB,K) = FUNC(Q)
              DO W = 1, BRDF_NPARS
                D_BRDFUNC_0(W,Q,I,IB,K)  = DFUNC(Q,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),           &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),          &
                 FUNC, DFUNC ) 
            DO Q = 1, NSTOKESSQ
              BRDFUNC(Q,I,J,K) = FUNC(Q)
              DO W = 1, BRDF_NPARS
                D_BRDFUNC(W,Q,I,J,K)  = DFUNC(Q,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),               &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),          &
                 FUNC, DFUNC ) 
              DO Q = 1, NSTOKESSQ
                EBRDFUNC(Q,I,KE,K) = FUNC(Q)
                DO W = 1, BRDF_NPARS
                  D_EBRDFUNC(W,Q,I,KE,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
          DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),            &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         &
                 FUNC, DFUNC )
              DO Q = 1, NSTOKESSQ
                USER_BRDFUNC_0(Q,UI,IB,K) = FUNC(Q)
                DO W = 1, BRDF_NPARS
                  D_USER_BRDFUNC_0(W,Q,UI,IB,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),          &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         &
                 FUNC, DFUNC )
              DO Q = 1, NSTOKESSQ
                USER_BRDFUNC(Q,UI,J,K) = FUNC(Q)
                DO W = 1, BRDF_NPARS
                  D_USER_BRDFUNC(W,Q,UI,J,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_VFUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, PARS, DERIVS, NSTOKESSQ,  & 
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),              &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),         &
                   FUNC, DFUNC )
                DO Q = 1, NSTOKESSQ
                  USER_EBRDFUNC(Q,UI,KE,K) = FUNC(Q)
                  DO W = 1, BRDF_NPARS
                    D_USER_EBRDFUNC(W,Q,UI,KE,K)  = DFUNC(Q,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE VLIDORT_BRDF_MAKER_PLUS

!

SUBROUTINE VLIDORT_GCMCRI_MAKER_PLUS                                         & 
           ( DO_USER_STREAMS, DO_SURFACE_EMISSION,                           & ! Inputs
             DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,                              & ! Inputs
             NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, NPARS,                    & ! Inputs
             NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,               & ! Inputs
             QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,             & ! Inputs
             SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, X_BRDF, CX_BRDF,  & ! Inputs
             SX_BRDF, CXE_BRDF, SXE_BRDF, PARS, DERIVS,                      & ! Inputs
             DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                        & ! Outputs
             BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,             & ! Outputs
             D_DBKERNEL_BRDFUNC, D_BRDFUNC, D_USER_BRDFUNC,                  & ! Outputs
             D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC )      ! Outputs

      implicit none

!  Prepares the bidirectional reflectance scatter matrices

!  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Input arguments
!  ===============

!  Local flags

      LOGICAL, intent(in) :: DO_USER_STREAMS
      LOGICAL, intent(in) :: DO_SURFACE_EMISSION
      LOGICAL, intent(in) :: DO_SHADOW_EFFECT
      LOGICAL, intent(in) :: DO_COXMUNK_DBMS

!  Number of Azimuth waudrature streams

      INTEGER, intent(in) :: NSTREAMS_BRDF
      INTEGER, intent(in) :: NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER, intent(in) :: NSTOKESSQ 

!  Local number of Kernel parameters

      INTEGER, intent(in) :: NPARS
      
!  Local angle control

      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: N_USER_RELAZMS

!  Local angles

      REAL(kind=8), intent(in) ::  PHIANG(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  COSPHI(MAX_USER_RELAZMS)
      REAL(kind=8), intent(in) ::  SINPHI(MAX_USER_RELAZMS)

      REAL(kind=8), intent(inout) ::  SZASURCOS(MAXBEAMS) !xliu, change from in to inout
      REAL(kind=8), intent(in) ::  SZASURSIN(MAXBEAMS)

      REAL(kind=8), intent(inout) ::  QUAD_STREAMS(MAXSTREAMS) !xliu, change from in to inout
      REAL(kind=8), intent(in) ::  QUAD_SINES  (MAXSTREAMS)

      REAL(kind=8), intent(inout) ::  USER_STREAMS(MAX_USER_STREAMS) !xliu, change from in to inout
      REAL(kind=8), intent(in) ::  USER_SINES  (MAX_USER_STREAMS)

!  Local parameter array

      REAL(kind=8), intent(in) ::  PARS   ( MAX_BRDF_PARAMETERS )
      LOGICAL     , intent(in) ::  DERIVS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      REAL(kind=8), intent(in) ::  X_BRDF  ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  CX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  SX_BRDF ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(inout) ::  CXE_BRDF ( MAXSTHALF_BRDF ) !xliu, change from in to inout
      REAL(kind=8), intent(in) ::  SXE_BRDF ( MAXSTHALF_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: BRDFUNC   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: BRDFUNC_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: USER_BRDFUNC   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_BRDFUNC_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: DBKERNEL_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: EBRDFUNC      ( MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: USER_EBRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output Linearizations of BRDF functions (parameter derivatives)
!  ===============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: D_USER_BRDFUNC     &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_BRDFUNC_0   &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      REAL(kind=8), intent(out) :: D_DBKERNEL_BRDFUNC  &
                                  ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity

      REAL(kind=8), intent(out) :: D_EBRDFUNC         &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(out) :: D_USER_EBRDFUNC    &
                                   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER          :: W, I, UI, J, K, KE, IB, Q
      REAL(kind=8)     :: DFUNC ( 16, MAX_BRDF_PARAMETERS ), FUNC(16)
      LOGICAL          :: DOSHADOW

!  Local
!  -----

      DOSHADOW = DO_SHADOW_EFFECT

!  Exact DB calculation
!  --------------------

      IF ( DO_COXMUNK_DBMS ) THEN  
        DO K = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              CALL GCMCRI_VFUNCTION_DB_PLUS                                       &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW,   & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),                  & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),                 & ! Inputs
                 FUNC, DFUNC )                                                      ! Output
              DO Q = 1, NSTOKESSQ
                DBKERNEL_BRDFUNC(Q,UI,K,IB) = FUNC(Q)
                DO W = 1, NPARS
                  D_DBKERNEL_BRDFUNC(W,Q,UI,K,IB)  = DFUNC(Q,W)
                ENDDO
              ENDDO  
            ENDDO
          ENDDO
        ENDDO
      ELSE   
        DO K = 1, N_USER_RELAZMS
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              CALL GCMCRI_VFUNCTION_PLUS                                          &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW,   & ! Inputs
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),                  & ! Inputs
                 USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K),                 & ! Inputs
                 FUNC, DFUNC )                                                      ! Output
              DO Q = 1, NSTOKESSQ
                DBKERNEL_BRDFUNC(Q,UI,K,IB) = FUNC(Q)
                DO W = 1, NPARS
                  D_DBKERNEL_BRDFUNC(W,Q,UI,K,IB)  = DFUNC(Q,W)
                ENDDO
              ENDDO  
            ENDDO
          ENDDO
        ENDDO
      ENDIF  

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam

      DO IB = 1, NBEAMS 
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GCMCRI_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW,  & 
                 SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),                  &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               &
                 FUNC, DFUNC ) 
            DO Q = 1, NSTOKESSQ
              BRDFUNC_0(Q,I,IB,K) = FUNC(Q)
              DO W = 1, NPARS
                D_BRDFUNC_0(W,Q,I,IB,K)  = DFUNC(Q,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GCMCRI_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW, & 
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),               &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              &
                 FUNC, DFUNC ) 
            DO Q = 1, NSTOKESSQ
              BRDFUNC(Q,I,J,K) = FUNC(Q)
              DO W = 1, NPARS
                D_BRDFUNC(W,Q,I,J,K)  = DFUNC(Q,W)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

      IF ( DO_SURFACE_EMISSION ) THEN
        DO I = 1, NSTREAMS
          DO KE = 1, NBRDF_HALF
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW, & 
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I),                   &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              &
                 FUNC, DFUNC ) 
              DO Q = 1, NSTOKESSQ
                EBRDFUNC(Q,I,KE,K) = FUNC(Q)
                DO W = 1, NPARS
                  D_EBRDFUNC(W,Q,I,KE,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam

        DO IB = 1, NBEAMS
          DO UI = 1, N_USER_STREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW, & 
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),                &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),             &
                 FUNC, DFUNC )
              DO Q = 1, NSTOKESSQ
                USER_BRDFUNC_0(Q,UI,IB,K) = FUNC(Q)
                DO W = 1, NPARS
                  D_USER_BRDFUNC_0(W,Q,UI,IB,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  incident quadrature directions

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GCMCRI_VFUNCTION_PLUS &
               ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW, & 
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),              &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),             &
                 FUNC, DFUNC )
              DO Q = 1, NSTOKESSQ
                USER_BRDFUNC(Q,UI,J,K) = FUNC(Q)
                DO W = 1, NPARS
                  D_USER_BRDFUNC(W,Q,UI,J,K)  = DFUNC(Q,W)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions

        IF ( DO_SURFACE_EMISSION ) THEN
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NBRDF_HALF
              DO K = 1, NSTREAMS_BRDF
                CALL GCMCRI_VFUNCTION_PLUS &
                 ( MAX_BRDF_PARAMETERS, NPARS, PARS, DERIVS, NSTOKESSQ, DOSHADOW, & 
                   CXE_BRDF(KE), SXE_BRDF(KE), USER_STREAMS(UI),                  &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),             &
                   FUNC, DFUNC )
                DO Q = 1, NSTOKESSQ
                  USER_EBRDFUNC(Q,UI,KE,K) = FUNC(Q)
                  DO W = 1, NPARS
                    D_USER_EBRDFUNC(W,Q,UI,KE,K)  = DFUNC(Q,W)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
END SUBROUTINE VLIDORT_GCMCRI_MAKER_PLUS

!

SUBROUTINE VLIDORT_BRDF_LS_FOURIER                                                                 &
         ( DO_USER_STREAMS, DO_SURFACE_EMISSION, LAMBERTIAN_FLAG, M, NSTOKES, NSTOKESSQ,           & ! Inputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,  BRDF_NPARS, BRDF_DERIVS,  & ! Inputs
           DELFAC, FACTOR, BRDF_COSAZMFAC, BRDF_SINAZMFAC, A_BRDF, BAX_BRDF,                       & ! Inputs
           D_BRDFUNC, D_USER_BRDFUNC, D_BRDFUNC_0, D_USER_BRDFUNC_0, D_EBRDFUNC, D_USER_EBRDFUNC,  & ! Inputs
           D_LOCAL_BRDF_F, D_LOCAL_BRDF_F_0, D_LOCAL_USER_BRDF_F, D_LOCAL_USER_BRDF_F_0,           & ! Outputs
           D_LOCAL_EMISSIVITY,  D_LOCAL_USER_EMISSIVITY )                                            ! Outputs

!  Prepares Fourier component of the bidirectional reflectance functions

      IMPLICIT NONE

!  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Input arguments
!  ===============

!  Control

      LOGICAL, intent(in) :: LAMBERTIAN_FLAG
      LOGICAL, intent(in) :: DO_USER_STREAMS
      LOGICAL, intent(in) :: DO_SURFACE_EMISSION

!  Local numbers

      INTEGER, intent(in) :: M, NSTOKES, NSTOKESSQ
      INTEGER, intent(in) :: NSTREAMS
      INTEGER, intent(in) :: NBEAMS
      INTEGER, intent(in) :: N_USER_STREAMS
      INTEGER, intent(in) :: NSTREAMS_BRDF, NBRDF_HALF

!  linearization Control

      INTEGER, intent(in)      :: BRDF_NPARS
      LOGICAL, intent(in)      :: BRDF_DERIVS ( MAX_BRDF_PARAMETERS )

!  Surface factors

      REAL(kind=8), intent(in) :: DELFAC, FACTOR

!  Azimuth cosines/sines and weights

      REAL(kind=8), intent(in) ::  BRDF_COSAZMFAC ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  BRDF_SINAZMFAC ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  A_BRDF         ( MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) ::  BAX_BRDF       ( MAXSTHALF_BRDF  )

!  Local Linearizations of BRDF functions (parameter derivatives)
!  ==============================================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(in) :: D_BRDFUNC   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF ) 
      REAL(kind=8), intent(in) :: D_BRDFUNC_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      REAL(kind=8), intent(in) :: D_USER_BRDFUNC   &
                                  ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) :: D_USER_BRDFUNC_0 &
                                  ( MAX_BRDF_PARAMETERS ,MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Values for Emissivity

      REAL(kind=8), intent(in) :: D_EBRDFUNC       &
                                  ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS,       MAXSTHALF_BRDF, MAXSTREAMS_BRDF )
      REAL(kind=8), intent(in) :: D_USER_EBRDFUNC  &
                                  ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTHALF_BRDF, MAXSTREAMS_BRDF )

!  Output: Derivative-kernel Fourier components
!  ============================================

!  at quadrature (discrete ordinate) angles

      REAL(kind=8), intent(out) :: D_LOCAL_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: D_LOCAL_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      REAL(kind=8), intent(out) :: D_LOCAL_USER_BRDF_F   ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      REAL(kind=8), intent(out) :: D_LOCAL_USER_BRDF_F_0 ( MAX_BRDF_PARAMETERS, MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      REAL(kind=8), intent(out) :: D_LOCAL_EMISSIVITY      ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAXSTREAMS       )
      REAL(kind=8), intent(out) :: D_LOCAL_USER_EMISSIVITY ( MAX_BRDF_PARAMETERS, MAXSTOKES, MAX_USER_STREAMS )

!  local variables
!  ===============

      INTEGER         :: I, UI, J, K, KPHI, IB, Q, O1, O2, W
      REAL(kind=8)    :: SUM, REFL, HELP, EMISS(16)
      INTEGER         :: COSSIN_MASK(16)
      DATA COSSIN_MASK / 1,1,2,0,1,1,2,0,2,2,1,0,0,0,0,1 /

!  surface factor

      HELP = HALF * DELFAC

!  Quadrature outgoing directions
!  ------------------------------

!  Zeroing

      D_LOCAL_BRDF_F        = ZERO
      D_LOCAL_BRDF_F_0      = ZERO
      D_LOCAL_USER_BRDF_F   = ZERO
      D_LOCAL_USER_BRDF_F_0 = ZERO

!  Incident Solar beam (direct beam reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO W = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(W) ) THEN
            DO IB = 1, NBEAMS
              DO I = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  SUM = ZERO
                  IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC_0(W,Q,I,IB,K)*BRDF_COSAZMFAC(K)
                    ENDDO
                  ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC_0(W,Q,I,IB,K)*BRDF_SINAZMFAC(K)
                    ENDDO
                  ENDIF
                  D_LOCAL_BRDF_F_0(W,Q,I,IB) = SUM * HELP
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  incident quadrature directions (surface multiple reflections)

      IF ( .NOT. LAMBERTIAN_FLAG ) THEN
        DO W = 1, BRDF_NPARS
          IF ( BRDF_DERIVS(W) ) THEN
            DO I = 1, NSTREAMS
              DO J = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  SUM = ZERO
                  IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC(W,Q,I,J,K)*BRDF_COSAZMFAC(K)
                    ENDDO
                  ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                    DO K = 1, NSTREAMS_BRDF
                      SUM  = SUM + D_BRDFUNC(W,Q,I,J,K)*BRDF_SINAZMFAC(K)
                    ENDDO
                  ENDIF
                  D_LOCAL_BRDF_F(W,Q,I,J) = SUM * HELP
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam (direct beam reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                  DO Q = 1, NSTOKESSQ
                    SUM = ZERO
                    IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC_0(W,Q,UI,IB,K)*BRDF_COSAZMFAC(K)
                      ENDDO
                    ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC_0(W,Q,UI,IB,K)*BRDF_SINAZMFAC(K)
                      ENDDO
                    ENDIF
                    D_LOCAL_USER_BRDF_F_0(W,Q,UI,IB) = SUM * HELP
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

!  incident quadrature directions (surface multiple reflections)

        IF ( .NOT. LAMBERTIAN_FLAG ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO UI = 1, N_USER_STREAMS
                DO J = 1, NSTREAMS
                  DO Q = 1, NSTOKESSQ
                    SUM = ZERO
                    IF ( COSSIN_MASK(Q).EQ.1 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC(W,Q,UI,J,K)*BRDF_COSAZMFAC(K)
                      ENDDO
                    ELSE IF ( COSSIN_MASK(Q).EQ.2 ) THEN
                      DO K = 1, NSTREAMS_BRDF
                        SUM  = SUM + D_USER_BRDFUNC(W,Q,UI,J,K)*BRDF_SINAZMFAC(K)
                      ENDDO
                    ENDIF
                    D_LOCAL_USER_BRDF_F(W,Q,UI,J) = SUM * HELP
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF

      ENDIF

!  Emissivity
!  ----------

!  Assumed to exist only for the total intensity
!        (first element of Stokes Vector) - is this right ??????

      IF ( DO_SURFACE_EMISSION ) THEN

!  Lambertian case

        IF ( LAMBERTIAN_FLAG.and.M.EQ.0 ) THEN
          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO I = 1, NSTREAMS
                D_LOCAL_EMISSIVITY(W,1,I) = ZERO
              ENDDO
              IF ( DO_USER_STREAMS ) THEN
                DO UI = 1, N_USER_STREAMS
                  D_LOCAL_USER_EMISSIVITY(W,1,UI) = ZERO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF

!  bidirectional reflectance

        IF ( .not. LAMBERTIAN_FLAG ) THEN

!  Inserted Polarization sum here.   Still to be checked.....!!!!!!!

!  Quadrature polar directions

          DO W = 1, BRDF_NPARS
            IF ( BRDF_DERIVS(W) ) THEN
              DO I = 1, NSTREAMS
                DO Q = 1, NSTOKESSQ
                  REFL = ZERO
                  DO KPHI= 1, NSTREAMS_BRDF
                    SUM = ZERO
                    DO K = 1, NBRDF_HALF
                      SUM = SUM + D_EBRDFUNC(W,Q,I,K,KPHI) * BAX_BRDF(K)
                    ENDDO
                    REFL = REFL + A_BRDF(KPHI) * SUM
                  ENDDO
                  EMISS(Q) = REFL
                ENDDO
                DO O1 = 1, NSTOKES
                  REFL = ZERO
                  DO O2 = 1, NSTOKES
                    Q = 4 * ( O1 - 1 ) + O2
                    REFL = REFL + EMISS(Q) 
                  ENDDO
                  D_LOCAL_EMISSIVITY(W,O1,I) = REFL * FACTOR
                ENDDO
              ENDDO
            ENDIF
          ENDDO

!   user-defined polar directions

          IF ( DO_USER_STREAMS ) THEN
            DO W = 1, BRDF_NPARS
              IF ( BRDF_DERIVS(W) ) THEN
                DO UI = 1, N_USER_STREAMS
                  DO Q = 1, NSTOKESSQ
                    REFL = ZERO
                    DO KPHI= 1, NSTREAMS_BRDF
                      SUM = ZERO
                      DO K = 1, NBRDF_HALF
                        SUM = SUM + D_USER_EBRDFUNC(W,Q,UI,K,KPHI) * BAX_BRDF(K)
                      ENDDO
                      REFL = REFL + A_BRDF(KPHI) * SUM
                    ENDDO
                    EMISS(Q) = REFL
                  ENDDO
                  DO O1 = 1, NSTOKES
                    REFL = ZERO
                    DO O2 = 1, NSTOKES
                      Q = 4 * ( O1 - 1 ) + O2
                      REFL = REFL + EMISS(Q) 
                    ENDDO
                    D_LOCAL_USER_EMISSIVITY(W,O1,UI) = REFL * FACTOR
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF

!  Not lambertian

        ENDIF

!  end emissivity clause

      ENDIF

!  Finish

      RETURN
END SUBROUTINE VLIDORT_BRDF_LS_FOURIER

!

SUBROUTINE VLIDORT_BRDF_INPUTS_PLUS ( FILNAM,                             & ! Inputs
           DO_USER_STREAMS, DO_BRDF_SURFACE, DO_SURFACE_EMISSION,         & ! Outputs
           NSTOKES, N_BRDF_KERNELS, WHICH_BRDF, BRDF_NAMES,               & ! Outputs
           LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,           & ! Outputs
           N_BRDF_PARAMETERS, BRDF_PARAMETERS,                            & ! Outputs
           DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,                             & ! Outputs
           DO_KERNEL_PARAMS_WFS, DO_KERNEL_FACTOR_WFS, DO_KPARAMS_DERIVS, & ! Outputs
           N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS, N_SURFACE_WFS,       & ! Outputs
           NBEAMS, NSTREAMS, N_USER_STREAMS, N_USER_RELAZMS,              & ! Outputs
           BEAM_SZAS, USER_ANGLES_INPUT, USER_RELAZMS,                    & ! Outputs
           STATUS, NMESSAGES, MESSAGES, ACTIONS )                           ! Outputs

!  Input routine for BRDF program

      implicit none

!  Include file of Dimensions

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Module arguments (input filename)

      CHARACTER*(*), intent(in)  :: FILNAM

!  stream angle flag

      LOGICAL, intent(out) :: DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL, intent(out) :: DO_BRDF_SURFACE

!  Surface emission

      LOGICAL, intent(out) :: DO_SURFACE_EMISSION

!  number of Stokes components

      INTEGER, intent(out) :: NSTOKES

!   Number and index-list and names of bidirectional functions

      INTEGER, intent(out)      :: N_BRDF_KERNELS
      INTEGER, intent(out)      :: WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER*10, intent(out) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER     , intent(out) :: N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      REAL(kind=8), intent(out) ::  BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL, intent(out) :: LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors

      REAL(kind=8), intent(out) :: BRDF_FACTORS ( MAX_BRDF_KERNELS )

!  Number of azimuth quadrature streams for BRDF

      INTEGER, intent(out) :: NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL, intent(out) :: DO_SHADOW_EFFECT

!  Multiple reflectance correction for direct beam flag
!              (only for GLITTER type kernels)

      LOGICAL, intent(out) :: DO_COXMUNK_DBMS

!   Flags for WF of bidirectional function parameters and factors

      LOGICAL, intent(out) :: DO_KERNEL_FACTOR_WFS  ( MAX_BRDF_KERNELS )
      LOGICAL, intent(out) :: DO_KERNEL_PARAMS_WFS  ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  derived quantity (tells you when to do BRDF derivatives)

      LOGICAL, intent(out) :: DO_KPARAMS_DERIVS  ( MAX_BRDF_KERNELS )

!  number of surface weighting functions

      INTEGER, intent(out) :: N_SURFACE_WFS
      INTEGER, intent(out) :: N_KERNEL_FACTOR_WFS
      INTEGER, intent(out) :: N_KERNEL_PARAMS_WFS

!  Number of discrete ordinate streams

      INTEGER, intent(out) :: NSTREAMS

!  Local angle control

      INTEGER, intent(out) :: NBEAMS
      INTEGER, intent(out) :: N_USER_STREAMS
      INTEGER, intent(out) :: N_USER_RELAZMS

!  Angles

      REAL(kind=8), intent(out) :: BEAM_SZAS         (MAXBEAMS)
      REAL(kind=8), intent(out) :: USER_RELAZMS      (MAX_USER_RELAZMS)
      REAL(kind=8), intent(out) :: USER_ANGLES_INPUT (MAX_USER_STREAMS)

!  Exception handling. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER, intent(out)        :: STATUS
      INTEGER, intent(out)        :: NMESSAGES
      CHARACTER*(*), intent(out)  :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(out)  :: ACTIONS (0:MAX_MESSAGES)

!  local variables
!  ===============

      CHARACTER(Len=9), parameter ::  PREFIX = 'VLIDORT -' 

      INTEGER           :: DUM_INDEX, DUM_NPARS
      CHARACTER(Len=10) :: DUM_NAME
      LOGICAL           :: ERROR
      CHARACTER(Len=80) :: PAR_STR
      LOGICAL           :: GFINDPAR
      INTEGER           :: I, J, K, L, FILUNIT, LEN_STRING, NM

      EXTERNAL             GFINDPAR
      EXTERNAL             LEN_STRING

!  Check list of Kernel names

      CHARACTER*10 BRDF_CHECK_NAMES ( MAXBRDF_IDX )
      DATA BRDF_CHECK_NAMES / &
          'Lambertian', &
          'Ross-thin ', &
          'Ross-thick', &
          'Li-sparse ', &
          'Li-dense  ', &
          'Hapke     ', &
          'Roujean   ', &
          'Rahman    ', &
          'Cox-Munk  ', &
          'GissCoxMnk', &
          'GCMcomplex', &
          'RondHerman', &
          'Breon     ', &
          'BPDF08Veg ', &
          'BPDF08Soil', &
          'BPDF2009  '/  

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize Angle control
!  ========================

      DO_USER_STREAMS = .FALSE.
      NSTREAMS = 0

      NBEAMS   = 0
      DO I = 1, MAXBEAMS
        BEAM_SZAS(I) = ZERO
      ENDDO
      N_USER_STREAMS = 0
      DO I = 1, MAX_USER_STREAMS
        USER_ANGLES_INPUT(I) = ZERO
      ENDDO
      N_USER_RELAZMS = 0
      DO I = 1, MAX_USER_RELAZMS
        USER_RELAZMS(I) = ZERO
      ENDDO
      NSTOKES = 0

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      DO_SHADOW_EFFECT    = .FALSE.
      DO_COXMUNK_DBMS     = .FALSE.
      DO_SURFACE_EMISSION = .FALSE.

      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        BRDF_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

      N_SURFACE_WFS  = 0
      N_KERNEL_FACTOR_WFS = 0
      N_KERNEL_PARAMS_WFS = 0
      DO K = 1, MAX_BRDF_KERNELS
        DO_KPARAMS_DERIVS(K) = .false.
        DO_KERNEL_FACTOR_WFS(K) = .FALSE.
        DO L = 1, MAX_BRDF_PARAMETERS
          DO_KERNEL_PARAMS_WFS(K,L) = .FALSE.
        ENDDO
      ENDDO

!  number of Stokes components
!  ===========================

      PAR_STR = 'Number of Stokes components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Read Angle stuff
!  ================

!  user-defined Stream angle

      PAR_STR = 'User-defined stream angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
         READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  Solar beams
!  ===========

!  number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        RETURN
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'TOA solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Azimuth angles
!  ==============

!  Number of angles

      PAR_STR = 'Number of user-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  check not exceeding dimensioned number

      IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
        NM = NM + 1
        MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
        ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
        STATUS       = VLIDORT_SERIOUS
        NMESSAGES    = NM
        RETURN
      ENDIF

!  Angles

      PAR_STR = 'User-defined relative azimuth angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, N_USER_RELAZMS
          READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User defined stream angles (should be positive)
!  ==========================

      IF ( DO_USER_STREAMS ) THEN

!  Number of angles

        PAR_STR = 'Number of user-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check dimension

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  Angles

        PAR_STR = 'User-defined stream angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES_INPUT(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      ENDIF

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Surface emission flag

      PAR_STR = 'Do surface emission?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  BRDF inputs

      IF ( DO_BRDF_SURFACE ) THEN

!  number of kernels, check this value

        PAR_STR = 'Number of bidirectional reflectance kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension (=3)'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of bidirectional reflectance streams'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  =  'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  Main kernel input

        PAR_STR = 'Kernel names, indices, amplitudes, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
               BRDF_NAMES(I), WHICH_BRDF(I), BRDF_FACTORS(I), &
              N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,3)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 56     FORMAT( A10, I2, F6.2, I2, 3F12.6 )

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do shadow effect for Cox-Munk kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
          ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Multiple reflectance DB correction (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do multiple reflectance for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_COXMUNK_DBMS
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Linearized input
!  ----------------

        PAR_STR = 'Kernels, indices, # pars, Jacobian flags'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,57,ERR=998) &
            DUM_NAME, DUM_INDEX,DUM_NPARS,DO_KERNEL_FACTOR_WFS(I), &
                (DO_KERNEL_PARAMS_WFS(I,J),J=1,3)

            IF ( DUM_NAME .NE. BRDF_NAMES(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Kernel name not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel name'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

            IF ( DUM_INDEX .NE. WHICH_BRDF(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input BRDF Index name not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of BRDF kernel Index'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

            IF ( DUM_NPARS .NE. N_BRDF_PARAMETERS(I) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Input Number of BRDF parameters not same as earlier list'
              ACTIONS(NM)  = 'Check second occurence of N_BRDF_PARAMETERS'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF

!  Compute total number of pars

            IF ( DO_KERNEL_FACTOR_WFS(I) ) THEN
              N_KERNEL_FACTOR_WFS = N_KERNEL_FACTOR_WFS  + 1
            ENDIF
            DO J = 1, N_BRDF_PARAMETERS(I)
              IF ( DO_KERNEL_PARAMS_WFS(I,J) ) THEN
                N_KERNEL_PARAMS_WFS = N_KERNEL_PARAMS_WFS + 1
              ENDIF
            ENDDO
            DO_KPARAMS_DERIVS(I) = (N_KERNEL_PARAMS_WFS.GT.0)

          ENDDO
          N_SURFACE_WFS = N_KERNEL_FACTOR_WFS+N_KERNEL_PARAMS_WFS
        ENDIF

        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
 57     FORMAT( A10, I3, I2, 1X, L2, 2X, 3L2 ) 

!  Check total number of BRDF weighting functions is not out of bounds

        IF ( N_SURFACE_WFS .GT. MAX_SURFACEWFS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of Surface WFs > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_SURFACEWFS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          RETURN
        ENDIF

!  Check Kernel indices are within bounds, Check BRDF name is on Accepted list

        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
            NM = NM + 1
            MESSAGES(NM) =  'Bad input: BRDF Index not on list of indices'
            ACTIONS(NM)  =  'Re-set input value: Look in VLIDORT.PARS for correct index'
            STATUS = VLIDORT_SERIOUS
            NMESSAGES = NM
            RETURN
          ELSE
            IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
              NM = NM + 1
              MESSAGES(NM) =  'Bad input: BRDF kernel name not one of Accepted list'
              ACTIONS(NM)  =  'Re-set input value: Look in VLIDORT.PARS for correct name'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              RETURN
            ENDIF
          ENDIF
        ENDDO

!  End BRDF clause
   
      ENDIF

!  Successful finish

      CLOSE(FILUNIT)
      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//FILNAM(1:LEN_STRING(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      RETURN

!  line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//PAR_STR(1:LEN_STRING(PAR_STR))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)

!  Finish

      RETURN
END SUBROUTINE VLIDORT_BRDF_INPUTS_PLUS

