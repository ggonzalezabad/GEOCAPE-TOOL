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
! #            LISPARSE_VFUNCTION_PLUS                          #
! #            LIDENSE_VFUNCTION_PLUS                           #
! #            HAPKE_VFUNCTION_PLUS                             #
! #            RAHMAN_VFUNCTION_PLUS                            #
! #            COXMUNK_VFUNCTION_PLUS                           #
! #            GISSCOXMUNK_VFUNCTION_PLUS                       #
! #            COXMUNK_VFUNCTION_DB_PLUS                        #
! #            GISSCOXMUNK_VFUNCTION_DB_PLUS                    #
! #            RHERMAN_VFUNCTION_PLUS                           #
! #            BREON_VFUNCTION_PLUS                             #
! # Additional code for complex RI Giss Cox-Munk, 15 march 2010.#
! #                                                             #
! #            GCMCRI_VFUNCTION_PLUS                            #
! #            GCMCRI_VFUNCTION_DB_PLUS                         #
! #                                                             #
! #                                                             #
! ###############################################################

!

SUBROUTINE HAPKE_VFUNCTION_PLUS                                     &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              HAPKE_VKERNEL, HAPKE_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      REAL(kind=8), intent(out) :: HAPKE_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: HAPKE_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  Hapke Kernel function.
!    - New version, Fresh Coding
!    - Old version uses DISORT code; for validation.

!  input variables:

!    XI, SXI  : Cosine/Sine of angle of reflection (positive)
!    XJ, SXJ  : Cosine/Sine of angle of incidence (positive)
!    XPHI     : Difference of azimuth angles of incidence and reflection
!    PARS(1)  : single scattering albedo in Hapke's BDR model
!    PARS(2)  : angular width parameter of opposition effect in Hapke's model
!    PARS(3)  : Empirical hot spot multiplier

!  local variables
!    B0_EMPIR : empirical factor to account for the finite size of
!               particles in Hapke's BDR model
!    B_HOT    : term that accounts for the opposition effect
!               (retroreflectance, hot spot) in Hapke's BDR model
!    CTHETA   : cosine of phase angle in Hapke's BDR model
!    GAMMA    : albedo factor in Hapke's BDR model
!    PHASE    : scattering phase function in Hapke's BDR model
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model

!  Local variables

      INTEGER       :: J, O1
      REAL(kind=8)  :: CTHETA, THETA, PHASE
      REAL(kind=8)  :: HOTSPOT, B0_EMPIR, HELP_HOT, B_HOT
      REAL(kind=8)  :: SSALBEDO, GAMMA, REFLEC, FUNCTION
      REAL(kind=8)  :: HELP_J, GHELP_J, TERM_J
      REAL(kind=8)  :: HELP_I, GHELP_I, TERM_I
      REAL(kind=8)  :: TI_TJ, DT1, DT2
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ
        HAPKE_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          HAPKE_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI

!  (1,1) function (scalar form)

!  geometrical part

!  This is the code that is in DISORT - not right, I think.
!        CTHETA = XI * XJ + DABS(SXI) *  DABS(SXJ) * CKPHI

      CTHETA = XI * XJ + SXI * SXJ * CKPHI
      THETA  = DACOS( CTHETA )
      PHASE  = ONE + HALF * CTHETA

!  hot spot parameterization

      HOTSPOT  = PARS(2)
      B0_EMPIR = PARS(3)
      HELP_HOT = HOTSPOT + DTAN ( HALF * THETA )
      B_HOT    = B0_EMPIR * HOTSPOT / HELP_HOT

!  Albedo parameterization

      SSALBEDO = PARS(1)
      GAMMA    = DSQRT ( ONE - SSALBEDO )
      HELP_J   = TWO * XJ
      GHELP_J  = ( ONE + HELP_J * GAMMA )
      TERM_J   = ( ONE + HELP_J ) / GHELP_J
      HELP_I   = TWO * XI
      GHELP_I  = ( ONE + HELP_I * GAMMA )
      TERM_I   = ( ONE + HELP_I ) / GHELP_I
      TI_TJ    = TERM_J * TERM_I

!  Function

      REFLEC           = SSALBEDO * QUARTER / ( XI + XJ )
      FUNCTION         = ( ONE + B_HOT ) * PHASE + TI_TJ - ONE
      HAPKE_VKERNEL(1) = REFLEC * FUNCTION

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  ssalbedo derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        DT1 = HAPKE_VKERNEL(1) / SSALBEDO
        DT2 = ( HELP_J / GHELP_J ) + ( HELP_I / GHELP_I )
        DT2 = DT2 * TI_TJ * HALF / GAMMA
        HAPKE_VDERIVATIVES(1,1) = DT1 + DT2 * REFLEC
      ENDIF

!  Hotspot  derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        DT1 = B_HOT * ( B0_EMPIR - B_HOT ) / B0_EMPIR / HOTSPOT
        HAPKE_VDERIVATIVES(1,2) = DT1 * REFLEC * PHASE
      ENDIF

!  empirical factor derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        DT1 = B_HOT / B0_EMPIR 
        HAPKE_VDERIVATIVES(1,3) = DT1 * REFLEC * PHASE
      ENDIF

!  Finish

      RETURN
END SUBROUTINE HAPKE_VFUNCTION_PLUS

!

SUBROUTINE LISPARSE_VFUNCTION_PLUS                                  &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              LISPARSE_VKERNEL, LISPARSE_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      REAL(kind=8), intent(out) :: LISPARSE_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: LISPARSE_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables

      INTEGER       :: J, O1
      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, QR
      REAL(kind=8)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_QR, DY_Q
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      DO O1 = 1, N_BRDF_STOKESSQ
        LISPARSE_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          LISPARSE_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      Q = ONE
      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .LE. ONE ) THEN
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      QR = Q * R 
      LISPARSE_VKERNEL(1) = HALF * P - QR

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LISPARSE_VDERIVATIVES(1,1) = - R * DY_Q
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_QR = DX_R * Q + DX_Q * R
        LISPARSE_VDERIVATIVES(1,2) = HALF * DX_P - DX_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LISPARSE_VFUNCTION_PLUS

!

SUBROUTINE LIDENSE_VFUNCTION_PLUS                                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              LIDENSE_VKERNEL, LIDENSE_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      REAL(kind=8), intent(out) :: LIDENSE_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: LIDENSE_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables

      INTEGER       :: J, O1
      REAL(kind=8)  :: X_INC, X_REF, SX_INC, SX_REF, ANG_P, TX
      REAL(kind=8)  :: T_INC, T_REF, T_INC_SQ, T_REF_SQ
      REAL(kind=8)  :: CKSI, DELTA, T, COST, SINT, DSQ, SINTCOST
      REAL(kind=8)  :: A, B, H, R, P, Q, DT1, DT2, DT2SQ, P_QR
      REAL(kind=8)  :: A2, R2, DX_H, DX_R, DX_Q, DX_P, DX_P_QR, DY_Q
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise
!    -- Return for special case

      DO O1 = 1, N_BRDF_STOKESSQ
        LIDENSE_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          LIDENSE_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( ( XI .EQ. XJ ) .AND. ( CKPHI.EQ.ONE ) ) RETURN

!  Function
!  ========

!  .. incidence

      TX       = SXJ / XJ
      T_INC    = PARS(2) * TX
      T_INC_SQ = T_INC * T_INC
      ANG_P    = DATAN ( T_INC )
      X_INC    = DCOS(ANG_P)
      SX_INC   = DSIN(ANG_P)

!  .. reflection

      TX       = SXI / XI
      T_REF    = PARS(2) * TX
      T_REF_SQ = T_REF * T_REF
      ANG_P    = DATAN ( T_REF )
      X_REF    = DCOS(ANG_P)
      SX_REF   = DSIN(ANG_P)

!  ksi cosine

      CKSI = X_INC  * X_REF + SX_INC * SX_REF * CKPHI

!  contributions P and R and derivatives (if flagged)

      P = ( ONE + CKSI ) / X_REF
      A = ( ONE / X_INC )
      B = ( ONE / X_REF )
      R = A + B

      DX_R = ZERO
      DX_P = ZERO
      IF ( DO_DERIV_PARS(2) ) THEN
        DX_R = R * ( ONE - ( ONE / A / B ) ) / PARS(2)
        R2   = R * R
        A2   = A * A
        DX_P = ( P * ( ONE + A2 ) - ( R2 / B ) ) / PARS(2) / A2
      ENDIF

!  evaluate cos(t)

      DT1   = T_REF_SQ + T_INC_SQ
      DT2   = T_INC * T_REF
      DT2SQ = DT2 * DT2
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      DSQ   = DELTA * DELTA
      H     = DSQRT ( DSQ + SKPHI * SKPHI * DT2SQ )
      COST  = PARS(1) * H / R

!  set Q function and its derivatives if flagged

      Q    = ONE
      DX_Q = ZERO
      DY_Q = ZERO
      IF ( COST .LE. ONE ) THEN
        T        = DACOS(COST)
        SINT     = DSQRT ( ONE - COST * COST )
        SINTCOST = SINT * COST
        Q = ONE -  ( ( T - SINTCOST ) / PIE )
        IF ( DO_DERIV_PARS(2) ) THEN
          DX_H = ( TWO * H * H - DSQ ) / H / PARS(2)
          DX_Q = TWO * SINTCOST * ( (DX_H/H) - (DX_R/R) ) / PIE
        ENDIF
        IF ( DO_DERIV_PARS(1) ) THEN
          DY_Q = TWO * SINTCOST / PIE / PARS(1)
        ENDIF
      ENDIF

!  set the kernel
!  --------------

      P_QR = P / Q / R 
      LIDENSE_VKERNEL(1) = P_QR - TWO

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

      IF ( DO_DERIV_PARS(1) ) THEN
        LIDENSE_VDERIVATIVES(1,1) = - P_QR * DY_Q / Q 
      ENDIF

      IF ( DO_DERIV_PARS(2) ) THEN
        DX_P_QR = ( DX_P / P ) - ( DX_R / R ) - ( DX_Q / Q )
        LIDENSE_VDERIVATIVES(1,2) = P_QR * DX_P_QR
      ENDIF

!  Finish

      RETURN
END SUBROUTINE LIDENSE_VFUNCTION_PLUS

!

SUBROUTINE RAHMAN_VFUNCTION_PLUS                                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              RAHMAN_VKERNEL, RAHMAN_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      REAL(kind=8), intent(out) :: RAHMAN_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: RAHMAN_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables

      INTEGER       :: J, O1
      REAL(kind=8)  :: T_INC, T_REF, DT1, DT2 
      REAL(kind=8)  :: CXI, DELTA, K1_SQ, FACT, K0, K1, K2
      REAL(kind=8)  :: HELPM, HELPR, HELPG, D_HELPM, D_FACT
      REAL(kind=8)  :: GEOM, PHASE, RFAC, RFAC1, D_K0, D_K1, D_K2
      REAL(kind=8)  :: XPHI, CKPHI

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ
        RAHMAN_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          RAHMAN_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO
      XPHI  = PIE - PHI
      CKPHI = - CPHI
      IF ( XI.EQ.ZERO .OR. XJ.EQ.ZERO ) RETURN

!  parameters

      K0 = PARS(1)
      K1 = PARS(2)
      K2 = PARS(3)

!  geometrical angle xi

      CXI = XI * XJ + SXI * SXJ * CKPHI
      IF ( CXI .GT. ONE ) CXI = ONE

!  Phase function

      K1_SQ = K1 * K1
      HELPM = ONE - K1_SQ 
      FACT  = ONE + K1_SQ + TWO * K1 * CXI
      PHASE = HELPM / ( FACT ** ONEP5 )

!  Delta and R-factor

      T_INC = SXI / XI
      T_REF = SXJ / XJ
      DT1   = T_INC*T_INC + T_REF*T_REF
      DT2   = T_INC * T_REF
      DELTA = DSQRT ( DT1 - TWO * DT2 * CKPHI )
      HELPR = ONE / ( ONE + DELTA )
      RFAC  = ( ONE - K0 ) * HELPR
      RFAC1 = ONE + RFAC

!  Geom factor and kernel

      HELPG = XI * XJ * ( XI + XJ )
      GEOM  = HELPG ** ( K2 - ONE)
      RAHMAN_VKERNEL(1) = K0 * PHASE * RFAC1 * GEOM

!  Other vector functions

!    P L A C E H O L D E R

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  K0 derivative

      IF ( DO_DERIV_PARS(1) ) THEN
        D_K0   = ( ONE / K0 ) - ( HELPR / RFAC1 )
        RAHMAN_VDERIVATIVES(1,1) = RAHMAN_VKERNEL(1) * D_K0
      ENDIF

!  Phase function derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        D_FACT  =   TWO * K1 + TWO * CXI
        D_HELPM = - TWO * K1
        D_K1    = ( D_HELPM / HELPM ) - ONEP5 * ( D_FACT / FACT )
        RAHMAN_VDERIVATIVES(1,2) = RAHMAN_VKERNEL(1) * D_K1
      ENDIF

!  K2 derivative

      IF ( DO_DERIV_PARS(3) ) THEN
        D_K2 = DLOG ( HELPG )
        RAHMAN_VDERIVATIVES(1,3) = RAHMAN_VKERNEL(1) * D_K2
      ENDIF

!  Finish

      RETURN
END SUBROUTINE RAHMAN_VFUNCTION_PLUS

!

SUBROUTINE COXMUNK_VFUNCTION_PLUS                                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CKPHI, SKPHI,                  &
              COXMUNK_VKERNEL, COXMUNK_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CKPHI, SKPHI

      REAL(kind=8), intent(out) :: COXMUNK_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: COXMUNK_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  Critical exponent taken out

      REAL(kind=8),PARAMETER    ::  CRITEXP = 88.0D0 

!  Local variables

      INTEGER       :: J, O1
      REAL(kind=8)  :: Z, Z1, Z2, Z2_SQ_M1, H1, H2, RP, RL, XMP
      REAL(kind=8)  :: A, B, TA, ARGUMENT, PROB, FAC1, FAC2, CKPHI_NEG
      REAL(kind=8)  :: S1, S2, S3, XXI, XXJ
      REAL(kind=8)  :: T1_I, T2_I, DCOT_I
      REAL(kind=8)  :: T1_R, T2_R, DCOT_R
      REAL(kind=8)  :: SHADOWI, SHADOWR, SHADOW

      REAL(kind=8)  :: H1H2, H2Z2, TA_SQ, DFAC2, DH1, DH2, DRP, DRL
      REAL(kind=8)  :: D_S1, D_S2, D_T1, D_T2
      REAL(kind=8)  :: D_SHADOWI, D_SHADOWR, D_SHADOW
      REAL(kind=8), external :: DERFC

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ
        COXMUNK_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          COXMUNK_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO

!  Comment. 18 January 2006.
!  We have found in comparisons with the Giss Cox-Munk code that
!  the input COSPHI (CKPHI) here is the negative of what we actually need,
!  so introduce local variable which takes care of this

!  Also removed factor of PIE in the kernel denominator
!   This makes the output exactly same as GISS model for R(1,1)

      CKPHI_NEG = - CKPHI

!  (1,1) function (scalar form)

!  ..Scatter angles

! old   Z = - XI * XJ + SXI * SXJ * CKPHI   
! old   IF ( Z .LT. MINUS_ONE) Z = MINUS_ONE
! old   Z1 = DACOS(-Z)
! old   Z2 = DCOS(Z1*HALF)

      Z = XI * XJ + SXI * SXJ * CKPHI_NEG
      IF ( Z .GT. ONE) Z = ONE
      Z1 = DACOS(Z)
      Z2 = DCOS(Z1*HALF)

!  .. Fresnel coefficients

      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
      H1 = PARS(2) * Z2
      H2 = DSQRT ( PARS(2) + Z2_SQ_M1 )
      H1H2 = H1 + H2
      RP = ( H1 - H2 ) / H1H2
      H2Z2 = Z2 + H2
      RL = ( Z2 - H2 ) / H2Z2
      XMP = HALF * ( RP*RP + RL*RL )

!  Coxmunk Function

      A = TWO * Z2                 
      B = ( XI + XJ ) / A                  
      IF ( B .GT. ONE ) B = ONE           
      A  = PIO2 - DASIN(B)
      TA = DTAN(A)
      TA_SQ = TA * TA            
      ARGUMENT = TA_SQ  / PARS(1)
      IF ( ARGUMENT .LT. CRITEXP ) THEN
        PROB = DEXP ( - ARGUMENT )
        FAC1 = PROB / PARS(1)
        FAC2 = QUARTER / XI / ( B ** FOUR )
        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
      ENDIF

!  inverse slope-squared derivative (regular term)

      IF ( DO_DERIV_PARS(1) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DFAC2 = ( PARS(1) - TA_SQ ) / PARS(1) / PARS(1)
          COXMUNK_VDERIVATIVES(1,1) = - COXMUNK_VKERNEL(1) * DFAC2
        ENDIF
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code

      S1 = DSQRT(PARS(1)/PIE)
      S3 = ONE/(DSQRT(PARS(1)))
      S2 = S3*S3

      XXI  = XI*XI
      DCOT_I = XI/DSQRT(ONE-XXI)
      T1_I   = DEXP(-DCOT_I*DCOT_I*S2)
      T2_I   = DERFC(DCOT_I*S3)
      SHADOWI = HALF * ( S1*T1_I/DCOT_I - T2_I )

      XXJ  = XJ*XJ
      DCOT_R = XJ/DSQRT(ONE-XXJ)
      T1_R   = DEXP(-DCOT_R*DCOT_R*S2)
      T2_R   = DERFC(DCOT_R*S3)
      SHADOWR = HALF * ( S1*T1_R/DCOT_R - T2_R )

      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      COXMUNK_VKERNEL(1) = COXMUNK_VKERNEL(1) * SHADOW

!  Scalar derivatives (1,1) matrix entry
!  -------------------------------------

!  inverse slope-squared derivative

      IF ( DO_DERIV_PARS(1) ) THEN

!  add the shadow

        D_S1 = HALF / PIE / S1
        D_S2 = - S2 * S2
        D_T1 = - T1_I * DCOT_I * DCOT_I * D_S2
        D_T2 = S2 * S2 * DCOT_I * S1 * T1_I 
        D_SHADOWI = HALF * ( D_S1*T1_I/DCOT_I + S1*D_T1/DCOT_I - D_T2 )

        D_T1 = - T1_R * DCOT_R * DCOT_R * D_S2
        D_T2 = S2 * S2 * DCOT_R * S1 * T1_R
        D_SHADOWR = HALF * ( D_S1*T1_R/DCOT_R + S1*D_T1/DCOT_R - D_T2 )
        D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )

        COXMUNK_VDERIVATIVES(1,1) = COXMUNK_VDERIVATIVES(1,1) * SHADOW + &
                                    COXMUNK_VKERNEL(1) * D_SHADOW / SHADOW

      ENDIF

!  square refractive index derivative

      IF ( DO_DERIV_PARS(2) ) THEN
        IF ( ARGUMENT .LT. CRITEXP ) THEN
          DH1 = Z2
          DH2 = HALF / H2
          DRP = ( DH1 * ( ONE - RP ) - DH2 * ( ONE + RP ) ) / H1H2
          DRL =  - DH2 * ( ONE + RL ) / H2Z2
          DFAC2 = ( RP*DRP + RL*DRL ) / XMP
          COXMUNK_VDERIVATIVES(1,2) = COXMUNK_VKERNEL(1) * DFAC2
        ENDIF
      ENDIF

!  Finish

      RETURN
END SUBROUTINE COXMUNK_VFUNCTION_PLUS

!

SUBROUTINE GISSCOXMUNK_VFUNCTION_PLUS                               &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, XPHI_REF, CKPHI_REF, SKPHI_REF,     &
              GISSCOXMUNK_VKERNEL, GISSCOXMUNK_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)     :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)     :: PARS ( MAXPARS )
      LOGICAL     , intent(in)     :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(inout)  :: XI, XJ
      REAL(kind=8), intent(in)     :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF

      REAL(kind=8), intent(out)    :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out)    :: GISSCOXMUNK_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  Critical exponent taken out

      REAL(kind=8),PARAMETER    ::  CRITEXP = 88.0D0 

!  Local variables

      INTEGER       :: J, O1, KERNELMASK(10), I, IM
      REAL(kind=8)  :: XPHI_INC, CKPHI_INC, SKPHI_INC
      REAL(kind=8)  :: VI1, VI2, VI3, VR1, VR2, VR3
      REAL(kind=8)  :: unit1, unit2, unit3, fact1, factor
      REAL(kind=8)  :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      REAL(kind=8)  :: TI1, TI2, TI3, TR1, TR2, TR3
      REAL(kind=8)  :: PI1, PII2, PI3, PR1, PR2, PR3
      REAL(kind=8)  :: PIKR, PRKI, TIKR, TRKI
      REAL(kind=8)  :: E1, E2, E3, E4, SIGMA2
      REAL(kind=8)  :: CF11, CF12, CF21, CF22, RDZ2, RDZ4
      REAL(kind=8)  :: VP1, VP2, VP3, DMOD, DEX, DCOEFF
      REAL(kind=8)  :: AF, AF11, AF12, AF21, AF22
      REAL(kind=8)  :: C21, C22, CTTTP, CTTPT, CTTPP
      REAL(kind=8)  :: CTPPT, CTPPP, CPTPP
      REAL(kind=8)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(kind=8)  :: SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT
      DATA KERNELMASK / 1,2,3,5,6,7,9,10,11,16 /

      REAL(kind=8)  :: D_S1, D_S2, D_T1, D_T2, DERFAC
      REAL(kind=8)  :: D_SHADOWI, D_SHADOWR, D_SHADOW
      REAL(kind=8)  :: D_DCOEFF_0, D_DCOEFF, D_DEX, D_ARGUMENT
      REAL(kind=8), external :: DERFC

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        GISSCOXMUNK_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          GISSCOXMUNK_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER) 

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = ONE
      CN2 = PARS(2)

!  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF
    
      RDZ2 = UNIT3*UNIT3
      RDZ4 = RDZ2*RDZ2

      DCOEFF_0   = ONE/(8.0D0*XI*XJ*DMOD*RDZ4*SIGMA2)
      ARGUMENT  = (UNIT1*UNIT1 + UNIT2*UNIT2) / (TWO*SIGMA2*RDZ2)
      IF ( ARGUMENT .GT. CRITEXP ) THEN
        DEX = ZERO
      ELSE
        DEX = DEXP(-ARGUMENT)
      ENDIF

      DCOEFF = DCOEFF_0 * FACT1 * FACT1 * DEX

!  Derivative of DCOEFF w.r.t. SIGMA2

      D_DCOEFF = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( ARGUMENT .LT. CRITEXP ) THEN
        D_ARGUMENT = - ARGUMENT / SIGMA2
        D_DCOEFF_0 = - DCOEFF_0 / SIGMA2
        D_DEX      = - D_ARGUMENT * DEX
        D_DCOEFF = FACT1*FACT1 * ( DCOEFF_0 * D_DEX + D_DCOEFF_0 * DEX )
       ENDIF
      ENDIF

      AF  = HALF * DCOEFF
      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

!  original code
!      R(1,1)=(AF11+AF12+AF21+AF22)*AF
!      R(1,2)=(AF11-AF12+AF21-AF22)*AF
!      R(2,1)=(AF11-AF22+AF12-AF21)*AF
!      R(2,2)=(AF11-AF12-AF21+AF22)*AF

!  Transcribed code
      GISSCOXMUNK_VKERNEL(1) = (AF11+AF12+AF21+AF22)*AF
      GISSCOXMUNK_VKERNEL(2) = (AF11-AF12+AF21-AF22)*AF
      GISSCOXMUNK_VKERNEL(5) = (AF11-AF22+AF12-AF21)*AF
      GISSCOXMUNK_VKERNEL(6) = (AF11-AF12-AF21+AF22)*AF

!  Original code
!      CI=(0D0, -1D0)
!      C21=DCONJG(CF21)
!      C22=DCONJG(CF22)
!      CTTTP=CF11*DCONJG(CF12)
!      CTTPT=CF11*C21
!      CTTPP=CF11*C22
!      CTPPT=CF12*C21
!      CTPPP=CF12*C22
!      CPTPP=CF21*C22

!  replica for real variables only
      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

!  original code
!      R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
!      R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
!      R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
!      R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
!      R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
!      R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
!      R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
!      R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
!      R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
!      R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
!      R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
!      R(4,4)=    ( CTTPP-CTPPT)*DCOEFF

!  New code (several entries are zero)

      GISSCOXMUNK_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
      GISSCOXMUNK_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF
      GISSCOXMUNK_VKERNEL(16) =    ( CTTPP-CTPPT)*DCOEFF

!  Derivative before shadow effect

      DERFAC = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( DCOEFF .NE. ZERO ) THEN
        DERFAC = D_DCOEFF / DCOEFF
       ENDIF
       DO I = 1, 10
        IM = KERNELMASK(I)
        GISSCOXMUNK_VDERIVATIVES(IM,1) = GISSCOXMUNK_VKERNEL(IM)*DERFAC
       ENDDO
      ENDIF

!  No Shadow code if not flagged

      IF ( PARS(3) .EQ. ZERO ) RETURN

!  Shadow code (includes derivative if flagged)

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3
      D_S1 = ONE / PIE / S1
      D_S2 = - S2 / SIGMA2

      SHADOWI   = ZERO
      D_SHADOWI = ZERO
      IF ( XI .NE. ONE ) THEN
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWI = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF

      SHADOWR   = ZERO
      D_SHADOWR = ZERO
      IF ( XJ .NE. ONE ) THEN
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWR = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF
      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, 10
       IM = KERNELMASK(I)
       GISSCOXMUNK_VKERNEL(IM) = GISSCOXMUNK_VKERNEL(IM) * SHADOW
      ENDDO

      IF ( DO_DERIV_PARS(1) ) THEN
       D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
       DO I = 1, 10
        IM = KERNELMASK(I)
        GISSCOXMUNK_VDERIVATIVES(IM,1) =  GISSCOXMUNK_VDERIVATIVES(IM,1) * SHADOW + &
                                          GISSCOXMUNK_VKERNEL(IM) * D_SHADOW / SHADOW
       ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GISSCOXMUNK_VFUNCTION_PLUS

!

SUBROUTINE GCMCRI_VFUNCTION_PLUS                                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, NSSQ, DO_SHADOW, &
              XJ, SXJ, XI, SXI, XPHI_REF, CKPHI_REF, SKPHI_REF,     &
              GCMCRI_VKERNEL, GCMCRI_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)     :: MAXPARS, NPARS, NSSQ
      LOGICAL     , intent(in)     :: DO_SHADOW
      REAL(kind=8), intent(in)     :: PARS ( MAXPARS )
      LOGICAL     , intent(in)     :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(inout)  :: XI, XJ
      REAL(kind=8), intent(in)     :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF

      REAL(kind=8), intent(out)    :: GCMCRI_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out)    :: GCMCRI_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  Critical exponent taken out

      REAL(kind=8),PARAMETER    ::  CRITEXP = 88.0D0 

!  Local variables

      INTEGER       :: J, O1, KERNELMASK(16), I, IM
      REAL(kind=8)  :: XPHI_INC, CKPHI_INC, SKPHI_INC
      REAL(kind=8)  :: VI1, VI2, VI3, VR1, VR2, VR3
      REAL(kind=8)  :: unit1, unit2, unit3, fact1, factor
      REAL(kind=8)  :: XI1, TI1, TI2, TI3, TR1, TR2, TR3
      REAL(kind=8)  :: PI1, PII2, PI3, PR1, PR2, PR3
      REAL(kind=8)  :: PIKR, PRKI, TIKR, TRKI
      REAL(kind=8)  :: E1, E2, E3, E4, SIGMA2
      REAL(kind=8)  :: RDZ2, RDZ4
      REAL(kind=8)  :: VP1, VP2, VP3, DMOD, DEX, DCOEFF
      REAL(kind=8)  :: AF, AF11, AF12, AF21, AF22
      REAL(kind=8)  :: S1, S2, S3, XXI, XXJ, T1, T2, DCOT
      REAL(kind=8)  :: SHADOWI, SHADOWR, SHADOW, DCOEFF_0, ARGUMENT

      COMPLEX (kind=8)  :: CI, CF11, CF12, CF21, CF22, C21, C22
      COMPLEX (kind=8)  :: CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      COMPLEX (kind=8)  :: CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

      DATA KERNELMASK / 1,2,3,5,6,7,9,10,11,4,8,12,13,14,15,16 /

      REAL(kind=8)  :: D_S1, D_S2, D_T1, D_T2, DERFAC
      REAL(kind=8)  :: D_SHADOWI, D_SHADOWR, D_SHADOW
      REAL(kind=8)  :: D_DCOEFF_0, D_DCOEFF, D_DEX, D_ARGUMENT
      REAL(kind=8), external :: DERFC

!  Initialise

      DO O1 = 1, NSSQ 
        GCMCRI_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          GCMCRI_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER) 

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!               Remark on Use of shadow effect
!               ------------------------------
!  Shadow effect is controlled by the third parameter. That is, if
!  PARS(3) not equal to then shadow effect will be included.
!    --- NPARS should always be 3 for this Kernel.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Slope square is PARS(1)

      SIGMA2 = PARS(1)

!  refractive Index, Real + Imaginary, are PARS 2 and 3

      CN1 = ( ONE, ZERO )
      CN2 = CMPLX ( PARS(2), PARS(3) )

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS
!   -------------------------------

!  this is the original code, but now C-variables are complex

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = CDSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  CALCULATION OF THE STOKES REFLECTION MATRIX
!  -------------------------------------------

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF
    
      RDZ2 = UNIT3*UNIT3
      RDZ4 = RDZ2*RDZ2

      DCOEFF_0   = ONE/(8.0D0*XI*XJ*DMOD*RDZ4*SIGMA2)
      ARGUMENT  = (UNIT1*UNIT1 + UNIT2*UNIT2) / (TWO*SIGMA2*RDZ2)
      IF ( ARGUMENT .GT. CRITEXP ) THEN
        DEX = ZERO
      ELSE
        DEX = DEXP(-ARGUMENT)
      ENDIF

      DCOEFF = DCOEFF_0 * FACT1 * FACT1 * DEX

!  Derivative of DCOEFF w.r.t. SIGMA2

      D_DCOEFF = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( ARGUMENT .LT. CRITEXP ) THEN
        D_ARGUMENT = - ARGUMENT / SIGMA2
        D_DCOEFF_0 = - DCOEFF_0 / SIGMA2
        D_DEX      = - D_ARGUMENT * DEX
        D_DCOEFF = FACT1*FACT1 * ( DCOEFF_0 * D_DEX + D_DCOEFF_0 * DEX )
       ENDIF
      ENDIF

      AF  = HALF * DCOEFF
      AF11 = CDABS(CF11)
      AF12 = CDABS(CF12)
      AF21 = CDABS(CF21)
      AF22 = CDABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

!  original code
!      R(1,1)=(AF11+AF12+AF21+AF22)*AF
!      R(1,2)=(AF11-AF12+AF21-AF22)*AF
!      R(2,1)=(AF11-AF22+AF12-AF21)*AF
!      R(2,2)=(AF11-AF12-AF21+AF22)*AF

!  Transcribed code

      GCMCRI_VKERNEL(1) = (AF11+AF12+AF21+AF22)*AF
      GCMCRI_VKERNEL(2) = (AF11-AF12+AF21-AF22)*AF
      GCMCRI_VKERNEL(5) = (AF11-AF22+AF12-AF21)*AF
      GCMCRI_VKERNEL(6) = (AF11-AF12-AF21+AF22)*AF

!  Key Debug statement to track down DMOD
!      write(*,*)'(1,1) Fresnel = ',(AF11+AF12+AF21+AF22)*HALF/DMOD

!  Original code

      CI  = ( ZERO, MINUS_ONE )
      C21=DCONJG(CF21)
      C22=DCONJG(CF22)
      CTTTP=CF11*DCONJG(CF12)
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

!  replica for real variables only
!      C21 = CF21
!      C22 = CF22
!      CTTTP=CF11*CF12
!      CTTPT=CF11*C21
!      CTTPP=CF11*C22
!      CTPPT=CF12*C21
!      CTPPP=CF12*C22
!      CPTPP=CF21*C22

!  original code
!      R(1,3)=    (-CTTTP-CPTPP)*DCOEFF
!      R(1,4)=-CI*( CTTTP+CPTPP)*DCOEFF
!      R(2,3)=    (-CTTTP+CPTPP)*DCOEFF
!      R(2,4)=-CI*( CTTTP-CPTPP)*DCOEFF
!      R(3,1)=    (-CTTPT-CTPPP)*DCOEFF
!      R(3,2)=    (-CTTPT+CTPPP)*DCOEFF
!      R(3,3)=    ( CTTPP+CTPPT)*DCOEFF
!      R(3,4)= CI*( CTTPP-CTPPT)*DCOEFF
!      R(4,1)= CI*( CTTPT+CTPPP)*DCOEFF
!      R(4,2)= CI*( CTTPT-CTPPP)*DCOEFF
!      R(4,3)=-CI*( CTTPP+CTPPT)*DCOEFF
!      R(4,4)=    ( CTTPP-CTPPT)*DCOEFF

!  New code

      GCMCRI_VKERNEL(3)  =    (-CTTTP-CPTPP)*DCOEFF
      GCMCRI_VKERNEL(7)  =    (-CTTTP+CPTPP)*DCOEFF
      GCMCRI_VKERNEL(9)  =    (-CTTPT-CTPPP)*DCOEFF
      GCMCRI_VKERNEL(10) =    (-CTTPT+CTPPP)*DCOEFF
      GCMCRI_VKERNEL(11) =    ( CTTPP+CTPPT)*DCOEFF
      GCMCRI_VKERNEL(16) =    ( CTTPP-CTPPT)*DCOEFF

      IF ( NSSQ.EQ.16) THEN
        GCMCRI_VKERNEL(4)  = - CI * ( CTTTP+CPTPP)*DCOEFF
        GCMCRI_VKERNEL(8)  = - CI * ( CTTTP-CPTPP)*DCOEFF
        GCMCRI_VKERNEL(12) =   CI * ( CTTPP-CTPPT)*DCOEFF
        GCMCRI_VKERNEL(13) =   CI * ( CTTPT+CTPPP)*DCOEFF
        GCMCRI_VKERNEL(14) =   CI * ( CTTPT-CTPPP)*DCOEFF
        GCMCRI_VKERNEL(15) = - CI * ( CTTPP+CTPPT)*DCOEFF
        GCMCRI_VKERNEL(16) =        ( CTTPP-CTPPT)*DCOEFF
      ENDIF

!  Derivative before shadow effect

      DERFAC = ZERO
      IF ( DO_DERIV_PARS(1) ) THEN
       IF ( DCOEFF .NE. ZERO ) THEN
        DERFAC = D_DCOEFF / DCOEFF
       ENDIF
       DO I = 1, NSSQ
        IM = KERNELMASK(I)
        GCMCRI_VDERIVATIVES(IM,1) = GCMCRI_VKERNEL(IM)*DERFAC
       ENDDO
      ENDIF

!  No Shadow code if not flagged

      IF ( .not. DO_SHADOW ) RETURN

!  Shadow code (includes derivative if flagged)

      S1 = DSQRT(TWO*SIGMA2/PIE)
      S3 = ONE/(DSQRT(TWO*SIGMA2))
      S2 = S3*S3
      D_S1 = ONE / PIE / S1
      D_S2 = - S2 / SIGMA2

      SHADOWI   = ZERO
      D_SHADOWI = ZERO
      IF ( XI .NE. ONE ) THEN
       XXI  = XI*XI
       DCOT = XI/DSQRT(ONE-XXI)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWI = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWI = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF

      SHADOWR   = ZERO
      D_SHADOWR = ZERO
      IF ( XJ .NE. ONE ) THEN
       XXJ  = XJ*XJ
       DCOT = XJ/DSQRT(ONE-XXJ)
       T1   = DEXP(-DCOT*DCOT*S2)
       T2   = DERFC(DCOT*S3)
       SHADOWR = HALF*(S1*T1/DCOT-T2)
       IF ( DO_DERIV_PARS(1) ) THEN
        D_T1 = - T1 * DCOT * DCOT * D_S2
        D_T2 = TWO * S2 * DCOT * T1 / PIE / S1
        D_SHADOWR = HALF * ( D_S1*T1/DCOT + S1*D_T1/DCOT - D_T2 )
       ENDIF
      ENDIF
      SHADOW = ONE/(ONE+SHADOWI+SHADOWR)

      DO I = 1, NSSQ
       IM = KERNELMASK(I)
       GCMCRI_VKERNEL(IM) = GCMCRI_VKERNEL(IM) * SHADOW
      ENDDO

      IF ( DO_DERIV_PARS(1) ) THEN
       D_SHADOW = - SHADOW * SHADOW * ( D_SHADOWI + D_SHADOWR )
       DO I = 1, NSSQ
        IM = KERNELMASK(I)
        GCMCRI_VDERIVATIVES(IM,1) =  GCMCRI_VDERIVATIVES(IM,1) * SHADOW + &
                                     GCMCRI_VKERNEL(IM) * D_SHADOW / SHADOW
       ENDDO
      ENDIF

!  Finish

      RETURN
END SUBROUTINE GCMCRI_VFUNCTION_PLUS

!^L

SUBROUTINE COXMUNK_VFUNCTION_DB_PLUS                                &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              COXMUNK_VKERNEL, COXMUNK_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(in)  :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI

      REAL(kind=8), intent(out) :: COXMUNK_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: COXMUNK_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer       :: n, k, i, i1, q, N_phiquad_HALF
      REAL(kind=8)  :: XM, SXM, sum_pr, pr, dfunc(16,3)
      REAL(kind=8)  :: sumr, w_p, reflec_0(16), reflec_1(16)
      double precision d_reflec_0(16,3), d_reflec_1(16,3)
      double precision phi_sub1, cphi_sub1, sphi_sub1
      double precision phi_sub2, cphi_sub2, sphi_sub2

!  Local quadrature stuff

      integer, parameter :: max_msrs_muquad  = 40
      integer, parameter :: max_msrs_phiquad = 100
      integer            :: n_muquad, n_phiquad

!  arrays

      REAL(kind=8)  :: X_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: W_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: SX_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: WXX_MUQUAD(max_msrs_muquad)

      REAL(kind=8)  :: X_PHIQUAD (max_msrs_phiquad)
      REAL(kind=8)  :: W_PHIQUAD (max_msrs_phiquad)

      REAL(kind=8)  :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: D_R0_QUAD_IN  (max_msrs_muquad,max_msrs_phiquad,2)
      REAL(kind=8)  :: D_R0_OUT_QUAD (max_msrs_muquad,max_msrs_phiquad,2)

!  Safety first zeroing, only of (1,1) component

      REFLEC_0(1) = ZERO
      REFLEC_1(1) = ZERO
      COXMUNK_VKERNEL(1) = ZERO

      DO Q = 1, NPARS
        IF ( DO_DERIV_PARS(Q) ) THEN
          COXMUNK_VDERIVATIVES(1,Q) = ZERO
          D_REFLEC_0(1,Q) = ZERO
          D_REFLEC_1(1,Q) = ZERO
        ENDIF
      ENDDO

!  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 100
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

!  Single scattering (zero order), Phi is in degrees here!

      CALL COXMUNK_VFUNCTION_PLUS                                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              REFLEC_0, D_REFLEC_0 )

!  Quadrature output for first order R/T calculations 

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          SPHI_SUB1 = DSIN(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          SPHI_SUB2 = DSIN(PHI_SUB2)
          CALL COXMUNK_VFUNCTION_PLUS                               &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XM, SXM, XI, SXI, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,     &
              R0_OUT_QUAD(1,K,N), DFUNC )
          D_R0_OUT_QUAD(K,N,1) = DFUNC(1,1)
          D_R0_OUT_QUAD(K,N,2) = DFUNC(1,2)
          CALL COXMUNK_VFUNCTION_PLUS                               &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XM, SXM, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,     &
              R0_QUAD_IN(1,K,N),  DFUNC )
          D_R0_QUAD_IN(K,N,1) = DFUNC(1,1)
          D_R0_QUAD_IN(K,N,2) = DFUNC(1,2)
        ENDDO
      ENDDO

!  compute the next order

      SUMR = ZERO
      DO K = 1, n_muquad
        SUM_PR = ZERO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          SUM_PR = SUM_PR + W_P * R0_QUAD_IN(1,K,N) * R0_OUT_QUAD(1,K,N)
        ENDDO
        SUMR = SUMR + SUM_PR * WXX_MUQUAD(K)
      ENDDO
      REFLEC_1(1) = SUMR

!  Compute total

      COXMUNK_VKERNEL(1) = REFLEC_0(1) + REFLEC_1(1)

!  Derivatives

      DO Q = 1, 2
       IF ( DO_DERIV_PARS(Q) ) THEN
        SUMR = ZERO
        DO K = 1, n_muquad
         PR = ZERO
         DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          PR = PR + W_P *   R0_QUAD_IN(1,K,N)   * D_R0_OUT_QUAD(K,N,Q) &
                  + W_P * D_R0_QUAD_IN(K,N,Q) *   R0_OUT_QUAD(1,K,N)
         ENDDO
         SUMR = SUMR + PR * WXX_MUQUAD(K)
        ENDDO
        D_REFLEC_1(1,Q) = SUMR
        COXMUNK_VDERIVATIVES(1,Q) = D_REFLEC_0(1,Q) + D_REFLEC_1(1,Q)
       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE COXMUNK_VFUNCTION_DB_PLUS

!

SUBROUTINE GISSCOXMUNK_VFUNCTION_DB_PLUS                            &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
              GISSCOXMUNK_VKERNEL, GISSCOXMUNK_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(inout)  :: XI, XJ !xliu, add XI,XJ
      REAL(kind=8), intent(in)  :: SXI, SXJ, PHI, CPHI, SKPHI !xliu: remove XI,XJ

      REAL(kind=8), intent(out) :: GISSCOXMUNK_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: GISSCOXMUNK_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer       :: n, k, i, i1, q, N_phiquad_HALF,o1,o2,o3
      REAL(kind=8)  :: XM, SXM, sum_pr(16), dfunc(16,3), sum
      REAL(kind=8)  :: sumr(16), w_p, reflec_0(16), reflec_1(16)
      double precision d_reflec_0(16,3), d_reflec_1(16,3)
      double precision phi_sub1, cphi_sub1, sphi_sub1, s1
      double precision phi_sub2, cphi_sub2, sphi_sub2, s2

!  Local quadrature stuff

      integer, parameter :: max_msrs_muquad  = 40
      integer, parameter :: max_msrs_phiquad = 100
      integer            :: n_muquad, n_phiquad

!  arrays

      REAL(kind=8)  :: X_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: W_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: SX_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: WXX_MUQUAD(max_msrs_muquad)

      REAL(kind=8)  :: X_PHIQUAD (max_msrs_phiquad)
      REAL(kind=8)  :: W_PHIQUAD (max_msrs_phiquad)

      REAL(kind=8)  :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: D_R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad,2)
      REAL(kind=8)  :: D_R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad,2)

!  Indices

      INTEGER       :: MASKIT(3,3),LNS,NS,NSS,M,M12,M13,M32

!  Safety first zeroing

      DO O1 = 1, N_BRDF_STOKESSQ
        REFLEC_0(O1) = ZERO
        REFLEC_1(O1) = ZERO
        GISSCOXMUNK_VKERNEL(O1) = ZERO
        DO Q = 1, NPARS
         IF ( DO_DERIV_PARS(Q) ) THEN
          GISSCOXMUNK_VDERIVATIVES(O1,Q) = ZERO
          D_REFLEC_0(O1,Q) = ZERO
          D_REFLEC_1(O1,Q) = ZERO
         ENDIF
       ENDDO
      ENDDO

!  Masking limits

      NSS = N_BRDF_STOKESSQ
      IF ( N_BRDF_STOKESSQ.EQ.1  ) LNS = 1
      IF ( N_BRDF_STOKESSQ.EQ.4  ) LNS = 2
      IF ( N_BRDF_STOKESSQ.EQ.9  ) LNS = 3
      NS = LNS
      IF ( N_BRDF_STOKESSQ.EQ.16 ) THEN
        LNS = 3
        NS = LNS + 1 
      ENDIF

!  masking array

      MASKIT(1,1) = 1
      MASKIT(1,2) = 2
      MASKIT(1,3) = 3
      MASKIT(2,1) = 5
      MASKIT(2,2) = 6
      MASKIT(2,3) = 7
      MASKIT(3,1) = 9
      MASKIT(3,2) = 10
      MASKIT(3,3) = 11

!  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 100
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

!  Single scattering (zero order), Phi is in degrees here!

      CALL GISSCOXMUNK_VFUNCTION_PLUS                                &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ,  &
              XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI,                    &
              REFLEC_0, D_REFLEC_0 )

!  Quadrature output for first order R/T calculations 

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          SPHI_SUB1 = DSIN(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          SPHI_SUB2 = DSIN(PHI_SUB2)
          CALL GISSCOXMUNK_VFUNCTION_PLUS                            &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ,  &
              XI, SXI, XM, SXM, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,      &
              R0_OUT_QUAD(1,K,N), DFUNC )
          DO O1 = 1, N_BRDF_STOKESSQ
            D_R0_OUT_QUAD(O1,K,N,1) = DFUNC(O1,1)
            D_R0_OUT_QUAD(O1,K,N,2) = DFUNC(O1,2)
          ENDDO
          CALL GISSCOXMUNK_VFUNCTION_PLUS                            &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ,  &
              XM, SXM, XJ, SXJ, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,      &
              R0_QUAD_IN(1,K,N),  DFUNC )
          DO O1 = 1, N_BRDF_STOKESSQ
            D_R0_QUAD_IN(O1,K,N,1) = DFUNC(O1,1)
            D_R0_QUAD_IN(O1,K,N,2) = DFUNC(O1,2)
          ENDDO
        ENDDO
      ENDDO

!  compute the next order

      DO M = 1, N_BRDF_STOKESSQ
        SUMR(M) = ZERO
      ENDDO
      DO K = 1, n_muquad
        DO M = 1, N_BRDF_STOKESSQ
          SUM_PR(M) = ZERO
        ENDDO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          DO O1 = 1, LNS
            DO O2 = 1, LNS
              M12 = MASKIT(O1,O2)
              SUM = ZERO
              DO O3 = 1, LNS
                M13 = MASKIT(O1,O3)
                M32 = MASKIT(O3,O2)
                SUM = SUM + R0_QUAD_IN(M13,K,N) * R0_OUT_QUAD(M32,K,N)
              ENDDO  
              SUM_PR(M12) = SUM_PR(M12) + W_P * SUM
            ENDDO
          ENDDO
          IF (NS.EQ.4) THEN
            SUM = R0_QUAD_IN(NSS,K,N) * R0_OUT_QUAD(NSS,K,N)
            SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
          ENDIF
        ENDDO
        DO M = 1, N_BRDF_STOKESSQ
          SUMR(M) =  SUMR(M) + SUM_PR(M) * WXX_MUQUAD(K)
        ENDDO
      ENDDO
      DO M = 1, N_BRDF_STOKESSQ
        REFLEC_1(M) = SUMR(M)
      ENDDO

!  Compute total

      DO M = 1, N_BRDF_STOKESSQ
        GISSCOXMUNK_VKERNEL(M) = REFLEC_0(M) + REFLEC_1(M)
      ENDDO

!  Derivatives

      DO Q = 1, NPARS
       IF ( DO_DERIV_PARS(Q) ) THEN
        DO M = 1, N_BRDF_STOKESSQ
         SUMR(M) = ZERO
        ENDDO
        DO K = 1, n_muquad
         DO M = 1, N_BRDF_STOKESSQ
          SUM_PR(M) = ZERO
         ENDDO
         DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          DO O1 = 1, LNS
            DO O2 = 1, LNS
              M12 = MASKIT(O1,O2)
              S1 = ZERO
              S2 = ZERO
              DO O3 = 1, LNS
                M13 = MASKIT(O1,O3)
                M32 = MASKIT(O3,O2)
                S1 = S1 + R0_QUAD_IN(M13,K,N) * D_R0_OUT_QUAD(M32,K,N,Q)
                S2 = S2 + D_R0_QUAD_IN(M13,K,N,Q) * R0_OUT_QUAD(M32,K,N)
              ENDDO  
              SUM_PR(M12) = SUM_PR(M12) + W_P * ( S1 + S2 )
            ENDDO
          ENDDO
          IF (NS.EQ.4) THEN
            SUM =   R0_QUAD_IN(NSS,K,N)   * D_R0_OUT_QUAD(NSS,K,N,Q) + &
                  D_R0_QUAD_IN(NSS,K,N,Q) *   R0_OUT_QUAD(NSS,K,N) 
            SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
          ENDIF
         ENDDO
         DO M = 1, N_BRDF_STOKESSQ
          SUMR(M) =  SUMR(M) + SUM_PR(M) * WXX_MUQUAD(K)
         ENDDO
        ENDDO
        DO M = 1, N_BRDF_STOKESSQ
          D_REFLEC_1(M,Q) = SUMR(M)
          GISSCOXMUNK_VDERIVATIVES(M,Q) =  D_REFLEC_0(M,Q) + D_REFLEC_1(M,Q)
        ENDDO
       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE GISSCOXMUNK_VFUNCTION_DB_PLUS

!

SUBROUTINE GCMCRI_VFUNCTION_DB_PLUS                        &
   ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, NSSQ, DO_SHADOW, &
     XJ, SXJ, XI, SXI, PHI, CPHI, SKPHI,                   &
     GCMCRI_VKERNEL, GCMCRI_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)  :: MAXPARS, NPARS, NSSQ
      REAL(kind=8), intent(in)  :: PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_DERIV_PARS ( MAXPARS )
      LOGICAL     , intent(in)  :: DO_SHADOW
      REAL(kind=8), intent(inout)  :: XI, XJ ! xliu, change xi, xj from in to inout
      REAL(kind=8), intent(in)  :: SXI, SXJ, PHI, CPHI, SKPHI ! xliu

      REAL(kind=8), intent(out) :: GCMCRI_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out) :: GCMCRI_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables
!  ---------------

!  help variables

      integer       :: n, k, i, i1, q, N_phiquad_HALF,o1,o2,o3
      REAL(kind=8)  :: XM, SXM, sum_pr(16), dfunc(16,3), sum
      REAL(kind=8)  :: sumr(16), w_p, reflec_0(16), reflec_1(16)
      double precision d_reflec_0(16,3), d_reflec_1(16,3)
      double precision phi_sub1, cphi_sub1, sphi_sub1, s1
      double precision phi_sub2, cphi_sub2, sphi_sub2, s2

!  Local quadrature stuff

      integer, parameter :: max_msrs_muquad  = 40
      integer, parameter :: max_msrs_phiquad = 100
      integer            :: n_muquad, n_phiquad

!  arrays

      REAL(kind=8)  :: X_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: W_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: SX_MUQUAD (max_msrs_muquad)
      REAL(kind=8)  :: WXX_MUQUAD(max_msrs_muquad)

      REAL(kind=8)  :: X_PHIQUAD (max_msrs_phiquad)
      REAL(kind=8)  :: W_PHIQUAD (max_msrs_phiquad)

      REAL(kind=8)  :: R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad)
      REAL(kind=8)  :: D_R0_QUAD_IN  (16,max_msrs_muquad,max_msrs_phiquad,2)
      REAL(kind=8)  :: D_R0_OUT_QUAD (16,max_msrs_muquad,max_msrs_phiquad,2)

!  Indices

      INTEGER       :: MASKIT(4,4),LNS,NS,NSS,M,M12,M13,M32

!  Safety first zeroing

      DO O1 = 1, NSSQ
        REFLEC_0(O1) = ZERO
        REFLEC_1(O1) = ZERO
        GCMCRI_VKERNEL(O1) = ZERO
        DO Q = 1, NPARS
         IF ( DO_DERIV_PARS(Q) ) THEN
          GCMCRI_VDERIVATIVES(O1,Q) = ZERO
          D_REFLEC_0(O1,Q) = ZERO
          D_REFLEC_1(O1,Q) = ZERO
         ENDIF
       ENDDO
      ENDDO

!  Masking limits

      NSS = NSSQ
      IF ( NSSQ.EQ.1  ) LNS = 1
      IF ( NSSQ.EQ.4  ) LNS = 2
      IF ( NSSQ.EQ.9  ) LNS = 3
      NS = LNS
      IF ( NSSQ.EQ.16 ) THEN
        LNS = 3
        NS = LNS + 1 
      ENDIF

!  masking array

      MASKIT(1,1) = 1
      MASKIT(1,2) = 2
      MASKIT(1,3) = 3
      MASKIT(2,1) = 5
      MASKIT(2,2) = 6
      MASKIT(2,3) = 7
      MASKIT(3,1) = 9
      MASKIT(3,2) = 10
      MASKIT(3,3) = 11

      MASKIT(1,4) = 4
      MASKIT(2,4) = 8
      MASKIT(3,4) = 12
      MASKIT(4,1) = 13
      MASKIT(4,2) = 14
      MASKIT(4,3) = 15
      MASKIT(4,4) = 16

!  Air to water, Polar quadrature

      n_muquad = 40
      CALL GAULEG ( ZERO, ONE, X_muquad, W_muquad, n_muquad )
      DO I = 1, N_MUQUAD
        XM = X_MUQUAD(I)
        SX_MUQUAD(I) = DSQRT(ONE-XM*XM)
        WXX_MUQUAD(I) = XM * XM * W_MUQUAD(I)
      ENDDO

!  Azimuth quadrature

      n_phiquad = 100
      N_phiquad_HALF = N_PHIQUAD / 2
      CALL GAULEG ( ZERO, ONE, X_PHIQUAD, W_PHIQUAD, N_PHIQUAD_HALF )
      DO I = 1, N_PHIQUAD_HALF
        I1 = I + N_PHIQUAD_HALF
        X_PHIQUAD(I1) = - X_PHIQUAD(I)
        W_PHIQUAD(I1) =   W_PHIQUAD(I)
      ENDDO
      DO I = 1, N_PHIQUAD
        X_PHIQUAD(I)  = PIE * X_PHIQUAD(I)
      ENDDO

!  Single scattering (zero order), Phi is in degrees here!

      CALL GCMCRI_VFUNCTION_PLUS                                     &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, NSSQ, DO_SHADOW,  &
              XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI,                    &
              REFLEC_0, D_REFLEC_0 )

!  Quadrature output for first order R/T calculations 

      DO K = 1, n_muquad
        XM  = X_MUQUAD(K)
        SXM = SX_MUQUAD(K)
        DO N = 1, N_PHIQUAD
          PHI_SUB1  = X_PHIQUAD(N)
          CPHI_SUB1 = DCOS(PHI_SUB1)
          SPHI_SUB1 = DSIN(PHI_SUB1)
          PHI_SUB2  = PHI*DEG_TO_RAD - X_PHIQUAD(N)
          CPHI_SUB2 = DCOS(PHI_SUB2)
          SPHI_SUB2 = DSIN(PHI_SUB2)
          CALL GCMCRI_VFUNCTION_PLUS                                 &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, NSSQ, DO_SHADOW,  &
              XI, SXI, XM, SXM, PHI_SUB2, CPHI_SUB2, SPHI_SUB2,      &
              R0_OUT_QUAD(1,K,N), DFUNC )
          DO O1 = 1, NSSQ
            D_R0_OUT_QUAD(O1,K,N,1) = DFUNC(O1,1)
            D_R0_OUT_QUAD(O1,K,N,2) = DFUNC(O1,2)
          ENDDO
          CALL GCMCRI_VFUNCTION_PLUS                                 &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, NSSQ, DO_SHADOW,  &
              XM, SXM, XJ, SXJ, PHI_SUB1, CPHI_SUB1, SPHI_SUB1,      &
              R0_QUAD_IN(1,K,N),  DFUNC )
          DO O1 = 1, NSSQ
            D_R0_QUAD_IN(O1,K,N,1) = DFUNC(O1,1)
            D_R0_QUAD_IN(O1,K,N,2) = DFUNC(O1,2)
          ENDDO
        ENDDO
      ENDDO

!  compute the next order

      DO M = 1, NSSQ
        SUMR(M) = ZERO
      ENDDO
      DO K = 1, n_muquad
        DO M = 1, NSSQ
          SUM_PR(M) = ZERO
        ENDDO
        DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          DO O1 = 1, LNS
            DO O2 = 1, LNS
              M12 = MASKIT(O1,O2)
              SUM = ZERO
              DO O3 = 1, LNS
                M13 = MASKIT(O1,O3)
                M32 = MASKIT(O3,O2)
                SUM = SUM + R0_QUAD_IN(M13,K,N) * R0_OUT_QUAD(M32,K,N)
              ENDDO  
              SUM_PR(M12) = SUM_PR(M12) + W_P * SUM
            ENDDO
          ENDDO
          IF (NS.EQ.4) THEN
            SUM = R0_QUAD_IN(NSS,K,N) * R0_OUT_QUAD(NSS,K,N)
            SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
          ENDIF
        ENDDO
        DO M = 1, NSSQ
          SUMR(M) =  SUMR(M) + SUM_PR(M) * WXX_MUQUAD(K)
        ENDDO
      ENDDO
      DO M = 1, NSSQ
        REFLEC_1(M) = SUMR(M)
      ENDDO

!  Compute total

      DO M = 1, NSSQ
        GCMCRI_VKERNEL(M) = REFLEC_0(M) + REFLEC_1(M)
      ENDDO

!  Derivatives

      DO Q = 1, NPARS
       IF ( DO_DERIV_PARS(Q) ) THEN
        DO M = 1, NSSQ
         SUMR(M) = ZERO
        ENDDO
        DO K = 1, n_muquad
         DO M = 1, NSSQ
          SUM_PR(M) = ZERO
         ENDDO
         DO N = 1, N_PHIQUAD
          W_P  = W_PHIQUAD(N)
          DO O1 = 1, LNS
            DO O2 = 1, LNS
              M12 = MASKIT(O1,O2)
              S1 = ZERO
              S2 = ZERO
              DO O3 = 1, LNS
                M13 = MASKIT(O1,O3)
                M32 = MASKIT(O3,O2)
                S1 = S1 + R0_QUAD_IN(M13,K,N) * D_R0_OUT_QUAD(M32,K,N,Q)
                S2 = S2 + D_R0_QUAD_IN(M13,K,N,Q) * R0_OUT_QUAD(M32,K,N)
              ENDDO  
              SUM_PR(M12) = SUM_PR(M12) + W_P * ( S1 + S2 )
            ENDDO
          ENDDO
          IF (NS.EQ.4) THEN
            SUM =   R0_QUAD_IN(NSS,K,N)   * D_R0_OUT_QUAD(NSS,K,N,Q) + &
                  D_R0_QUAD_IN(NSS,K,N,Q) *   R0_OUT_QUAD(NSS,K,N) 
            SUM_PR(NSS) = SUM_PR(NSS) + W_P * SUM
          ENDIF
         ENDDO
         DO M = 1, NSSQ
          SUMR(M) =  SUMR(M) + SUM_PR(M) * WXX_MUQUAD(K)
         ENDDO
        ENDDO
        DO M = 1, NSSQ
          D_REFLEC_1(M,Q) = SUMR(M)
          GCMCRI_VDERIVATIVES(M,Q) =  D_REFLEC_0(M,Q) + D_REFLEC_1(M,Q)
        ENDDO
       ENDIF
      ENDDO

!  Finish

      RETURN
END SUBROUTINE GCMCRI_VFUNCTION_DB_PLUS

!

SUBROUTINE RHERMAN_VFUNCTION_PLUS                                   &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, XPHI_REF, CKPHI_REF, SKPHI_REF,     &
              RHERMAN_VKERNEL, RHERMAN_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)     :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)     :: PARS ( MAXPARS )
      LOGICAL     , intent(in)     :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(inout)  :: XI, XJ
      REAL(kind=8), intent(in)     :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF

      REAL(kind=8), intent(out)    :: RHERMAN_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out)    :: RHERMAN_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables

      INTEGER       :: J, O1

      REAL(kind=8)  :: RAHMAN_VKERNEL      (MAXSTOKES_SQ)
      REAL(kind=8)  :: HFUNCTION

      REAL(kind=8)  :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
      REAL(kind=8)  :: VI1, VI2, VI3, VR1, VR2, VR3
      REAL(kind=8)  :: VP1, VP2, VP3, DMOD
      REAL(kind=8)  :: unit1, unit2, unit3, fact1, factor
      REAL(kind=8)  :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      REAL(kind=8)  :: TI1, TI2, TI3, TR1, TR2, TR3
      REAL(kind=8)  :: PI1, PII2, PI3, PR1, PR2, PR3
      REAL(kind=8)  :: PIKR, PRKI, TIKR, TRKI
      REAL(kind=8)  :: E1, E2, E3, E4
      REAL(kind=8)  :: CF11, CF12, CF21, CF22
      REAL(kind=8)  :: AF11, AF12, AF21, AF22
      REAL(kind=8)  :: C21, C22, CTTTP, CTTPT, CTTPP
      REAL(kind=8)  :: CTPPT, CTPPP, CPTPP

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ
        RHERMAN_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          RHERMAN_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO

!#####################################################################
!#####################################################################
!  COXMUNK scalar stuff...................
!  Also removed factor of PIE in the kernel denominator
!   This makes the output exactly same as GISS model for R(1,1)
!      CKPHI_NEG = - CKPHI
!      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
!      IF ( Z .GT. ONE) Z = ONE
!      Z1 = DACOS(Z)
!      Z2 = DCOS(Z1*HALF)
!  .. Fresnel coefficients
!      REFIDX = 1.5d0
!      REFIDX_SQ = REFIDX * REFIDX
!      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
!      H1 = REFIDX_SQ * Z2
!      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
!      RP = ( H1 - H2 ) / ( H1 + H2 )
!      RL = ( Z2 - H2 ) / ( Z2 + H2 )
!  XMP = r(i) = (1,1) reflection coefficient for natural lights
!      XMP = HALF * ( RP*RP + RL*RL )
!  Comment. This is the Lenoble SOS code
!      xind=1.50
!      xx=cos(pi*thd/360.)
!      yy=sqrt(xind*xind+xx*xx-1.0)
!      zz=xx*xind*xind
!      rl=(zz-yy)/(zz+yy)
!      rr=(xx-yy)/(xx+yy)
!  where
!     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
!     xx   = Z2
!     yy   = H2
!     zz   = H1
!     xind = REFIDX
!     rl   = RP
!     rr   = RL
!  Coxmunk Function
!      A = TWO * Z2
!      B = ( XI + XJ ) / A
!      IF ( B .GT. ONE ) B = ONE
!      A = PIO2 - DASIN(B)
!      TA = DTAN(A)
!      ARGUMENT = TA * TA  / PARS(1)
!      IF ( ARGUMENT .LT. CRITEXP ) THEN
!        PROB = DEXP ( - ARGUMENT )
!        FAC1 = PROB / PARS(1)
!        FAC2 = QUARTER / XI / ( B ** FOUR )
!        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
!      ENDIF
!#####################################################################
!#####################################################################

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR A SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY.

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   RHERMAN_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

!  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

!  Set the H-function

      HFUNCTION = 0.25d0 / ( XI + XJ )

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      RHERMAN_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      RHERMAN_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      RHERMAN_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      RHERMAN_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

!  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      RHERMAN_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      RHERMAN_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      RHERMAN_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      RHERMAN_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      RHERMAN_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      RHERMAN_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

!  Add the Diffuse term to R(1,1)
!  ------------------------------

!   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
!   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

!  This is just the Rahman Kernel.........different name !!

      XPHI_RAH = 180.0d0 - XPHI_REF
      CALL RAHMAN_VFUNCTION_PLUS                                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, XPHI_RAH, CKPHI_REF, SKPHI_REF,     &
              RAHMAN_VKERNEL, RHERMAN_VDERIVATIVES )

!  Add to the specular term

      RHERMAN_VKERNEL(1) = RHERMAN_VKERNEL(1) + RAHMAN_VKERNEL(1)

!  Here is the original code from France........
!    cs: cosinus of the solar zenith angle
!    cv: cosinus of the viewing zenith angle
!    phi: azimuth between the solar and observation vertical planes
!      xx=cs**(vk-1.)
!      yy=cs**(vk-1.)
!      zz=(cs+cv)**(1.-vk)
!      FF1=rho*xx*yy/zz
!      xx=sqrt(1-cs*cs)
!      yy=sqrt(1-cv*cv)
!      ww=cs*cv-xx*yy*cos(pi*phi/180.)
!      aa=1+gt*gt+2*gt*ww
!      FF2=(1-gt*gt)/(aa**1.5)
!      vv=xx/cs
!      ww=yy/cv
!      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
!      FF3=1+(1-rho)/(1+G)
!      RBD=FF1*FF2*FF3
! RBD...........................IS THE REFLECTION COEFFICIENT

!  Finish

      RETURN
END SUBROUTINE RHERMAN_VFUNCTION_PLUS

!

SUBROUTINE BREON_VFUNCTION_PLUS                                     &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, XPHI_REF, CKPHI_REF, SKPHI_REF,     &
              BREON_VKERNEL, BREON_VDERIVATIVES )

      implicit none

!  include file of constants

      INCLUDE '../includes/VLIDORT.PARS_F90'

!  Subroutine arguments

      INTEGER     , intent(in)     :: MAXPARS, NPARS, N_BRDF_STOKESSQ
      REAL(kind=8), intent(in)     :: PARS ( MAXPARS )
      LOGICAL     , intent(in)     :: DO_DERIV_PARS ( MAXPARS )
      REAL(kind=8), intent(inout)  :: XI, XJ
      REAL(kind=8), intent(in)     :: SXI, SXJ, XPHI_REF, CKPHI_REF, SKPHI_REF


      REAL(kind=8), intent(out)    :: BREON_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8), intent(out)    :: BREON_VDERIVATIVES ( MAXSTOKES_SQ, MAXPARS )

!  local variables

      INTEGER       :: J, O1

      REAL(kind=8)  :: RAHMAN_VKERNEL(MAXSTOKES_SQ)
      REAL(kind=8)  :: HFUNCTION

      REAL(kind=8)  :: XPHI_INC, CKPHI_INC, SKPHI_INC, XPHI_RAH
      REAL(kind=8)  :: VI1, VI2, VI3, VR1, VR2, VR3
      REAL(kind=8)  :: VP1, VP2, VP3, DMOD
      REAL(kind=8)  :: unit1, unit2, unit3, fact1, factor
      REAL(kind=8)  :: XI1, CN1, CN2, CXI2, C2, C1, CRPER, CRPAR
      REAL(kind=8)  :: TI1, TI2, TI3, TR1, TR2, TR3
      REAL(kind=8)  :: PI1, PII2, PI3, PR1, PR2, PR3
      REAL(kind=8)  :: PIKR, PRKI, TIKR, TRKI
      REAL(kind=8)  :: E1, E2, E3, E4
      REAL(kind=8)  :: CF11, CF12, CF21, CF22
      REAL(kind=8)  :: AF11, AF12, AF21, AF22
      REAL(kind=8)  :: C21, C22, CTTTP, CTTPT, CTTPP
      REAL(kind=8)  :: CTPPT, CTPPP, CPTPP

!  Initialise

      DO O1 = 1, N_BRDF_STOKESSQ 
        BREON_VKERNEL(O1) = ZERO
        DO J = 1, NPARS
          BREON_VDERIVATIVES(O1,J) = ZERO
        ENDDO
      ENDDO

!#####################################################################
!#####################################################################
!  COXMUNK scalar stuff...................
!  Also removed factor of PIE in the kernel denominator
!   This makes the output exactly same as GISS model for R(1,1)
!      CKPHI_NEG = - CKPHI
!      Z = XI * XJ + SXI * SXJ * CKPHI_NEG 
!      IF ( Z .GT. ONE) Z = ONE
!      Z1 = DACOS(Z)
!      Z2 = DCOS(Z1*HALF)
!  .. Fresnel coefficients
!      REFIDX = 1.5d0
!      REFIDX_SQ = REFIDX * REFIDX
!      Z2_SQ_M1 = Z2 * Z2 + MINUS_ONE
!      H1 = REFIDX_SQ * Z2
!      H2 = DSQRT ( REFIDX_SQ + Z2_SQ_M1 )
!      RP = ( H1 - H2 ) / ( H1 + H2 )
!      RL = ( Z2 - H2 ) / ( Z2 + H2 )
!  XMP = r(i) = (1,1) reflection coefficient for natural lights
!      XMP = HALF * ( RP*RP + RL*RL )
!  Comment. This is the Lenoble SOS code
!      xind=1.50
!      xx=cos(pi*thd/360.)
!      yy=sqrt(xind*xind+xx*xx-1.0)
!      zz=xx*xind*xind
!      rl=(zz-yy)/(zz+yy)
!      rr=(xx-yy)/(xx+yy)
!  where
!     thd  = z1 * 180 / pi  (in degrees, Z1 is in radians)
!     xx   = Z2
!     yy   = H2
!     zz   = H1
!     xind = REFIDX
!     rl   = RP
!     rr   = RL
!  Coxmunk Function
!      A = TWO * Z2
!      B = ( XI + XJ ) / A
!      IF ( B .GT. ONE ) B = ONE
!      A = PIO2 - DASIN(B)
!      TA = DTAN(A)
!      ARGUMENT = TA * TA  / PARS(1)
!      IF ( ARGUMENT .LT. CRITEXP ) THEN
!        PROB = DEXP ( - ARGUMENT )
!        FAC1 = PROB / PARS(1)
!        FAC2 = QUARTER / XI / ( B ** FOUR )
!        COXMUNK_VKERNEL(1) = XMP * FAC1 * FAC2 / XJ
!      ENDIF
!#####################################################################
!#####################################################################

!  Transcription of the RMATR subroutine from Mishchenko/Travis code.

!   CALCULATION OF THE STOKES REFLECTION MATRIX FOR
!   ILLUMINATION FROM ABOVE FOR
!   A STATISTICALLY ROUGH SURFACE SEPARATING TWO HALF-SPACES
!   WITH REFRACTIVE INDICES OF THE UPPER AND LOWER HALF-SPACES EQUAL TO
!   CN1 AND CN2, RESPECTIVELY. THE EFFECT OF SHADOWING IS NOT
!   INCLUDED IN THIS SUBROUTINE BUT IS ADDED IN THE MAIN PROGRAM.

!   SIGMA2 = s**2 = MEAN SQUARE SURFACE SLOPE (EQ. (18) IN THE JGR PAPER) 

!   XI = ABS(COSINE OF THE INCIDENT ZENITH ANGLE)
!   XJ = ABS(COSINE OF THE REFLECTION ZENITH ANGLE).
!   SXI and SXJ are the respective SINES (input)
!   XPHI_REF = REFLECTION AZIMUTH ANGLE
!   GISSCOXMUNK_VKERNEL(16-elements) = (4X4) REFLECTION MATRIX

!  For real case, incident azimuth taken to be zero

      XPHI_INC  = ZERO
      CKPHI_INC = ONE
      SKPHI_INC = ZERO

!  Check for limiting cases

      IF(DABS(XI-1D0).LT.1d-9) XI = 0.999999999999d0
      IF(DABS(XJ-1D0).LT.1d-9) XJ = 0.999999999999d0

!  help variables (coordinate transformations)

      VI1 = SXI * CKPHI_INC
      VI2 = SXI * SKPHI_INC
      VI3 = -XI
      VR1 = SXJ * CKPHI_REF
      VR2 = SXJ * SKPHI_REF
      VR3 = XJ

!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION (normalized to 1)

      UNIT1  = VI1-VR1
      UNIT2  = VI2-VR2
      UNIT3  = VI3-VR3
      FACT1  = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      FACTOR = DSQRT(ONE/FACT1)

!   FRESNEL REFLECTION COEFFICIENTS, assume only real for now
!   ---------------------------------------------------------

      CN1 = ONE
      CN2 = 1.5d0

!  this is the original code, but now C-variables are real

      XI1 =  FACTOR*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      CXI2 = ONE - (ONE-XI1*XI1)*CN1*CN1/(CN2*CN2)
      CXI2 = DSQRT(CXI2)
      C1 = CN1*XI1
      C2 = CN2*CXI2
      CRPER = (C1-C2)/(C1+C2)
      C1 = CN2*XI1
      C2 = CN1*CXI2
      CRPAR = (C1-C2)/(C1+C2)

!  CALCULATION OF THE AMPLITUDE SCATTERING MATRIX
!  ----------------------------------------------

      TI1 = - XI * CKPHI_INC
      TI2 = - XI * SKPHI_INC
      TI3 = - SXI

      TR1 = + XJ * CKPHI_REF
      TR2 = + XJ * SKPHI_REF
      TR3 = -SXJ

      PI1  = - SKPHI_INC
      PII2 = + CKPHI_INC
      PI3  = ZERO

      PR1 = - SKPHI_REF
      PR2 = + CKPHI_REF
      PR3 = ZERO

      PIKR = PI1*VR1 + PII2*VR2 + PI3*VR3
      PRKI = PR1*VI1 + PR2*VI2  + PR3*VI3
      TIKR = TI1*VR1 + TI2*VR2  + TI3*VR3
      TRKI = TR1*VI1 + TR2*VI2  + TR3*VI3

      E1 = PIKR*PRKI
      E2 = TIKR*TRKI
      E3 = TIKR*PRKI
      E4 = PIKR*TRKI

!  Set the H-function

      HFUNCTION = 0.25d0 / ( XI * XJ )

!  Settting (1,1), (1,2), (2,1) and (2,2) components
!  -----------------------------------------------

      CF11 =  E1*CRPER+E2*CRPAR
      CF12 = -E3*CRPER+E4*CRPAR
      CF21 = -E4*CRPER+E3*CRPAR
      CF22 =  E2*CRPER+E1*CRPAR

!  Not to forget the normalization

      VP1 = VI2*VR3-VI3*VR2
      VP2 = VI3*VR1-VI1*VR3
      VP3 = VI1*VR2-VI2*VR1
      DMOD = VP1*VP1+VP2*VP2+VP3*VP3
      DMOD = DMOD*DMOD

!  if DMOD = 0, that is | n x n_0 | ^ 4 = 0 in M-T formula)
!    Then we need to set the ratio CF11 / DMOD

      IF ( DMOD .EQ. ZERO ) THEN
        CF11 = CRPAR
        CF22 = CRPER
        DMOD = 1.0d0
      ENDIF

!  Final H-function needs to be normalized

      HFUNCTION = HFUNCTION / DMOD

!  Compute the Electromagnetic contributions

      AF11 = DABS(CF11)
      AF12 = DABS(CF12)
      AF21 = DABS(CF21)
      AF22 = DABS(CF22)
      AF11 = AF11*AF11
      AF12 = AF12*AF12
      AF21 = AF21*AF21
      AF22 = AF22*AF22

      FACTOR = HALF * HFUNCTION
      BREON_VKERNEL(1) = (AF11+AF12+AF21+AF22) * FACTOR
      BREON_VKERNEL(2) = (AF11-AF12+AF21-AF22) * FACTOR
      BREON_VKERNEL(5) = (AF11-AF22+AF12-AF21) * FACTOR
      BREON_VKERNEL(6) = (AF11-AF12-AF21+AF22) * FACTOR

!  Setting (1,3), (2,3), (3,1), (3,2), (3,3) and (4,4) components
!  --------------------------------------------------------------

      C21 = CF21
      C22 = CF22
      CTTTP=CF11*CF12
      CTTPT=CF11*C21
      CTTPP=CF11*C22
      CTPPT=CF12*C21
      CTPPP=CF12*C22
      CPTPP=CF21*C22

      FACTOR = HFUNCTION
      BREON_VKERNEL(3)  =    (-CTTTP-CPTPP) * FACTOR
      BREON_VKERNEL(7)  =    (-CTTTP+CPTPP) * FACTOR
      BREON_VKERNEL(9)  =    (-CTTPT-CTPPP) * FACTOR
      BREON_VKERNEL(10) =    (-CTTPT+CTPPP) * FACTOR
      BREON_VKERNEL(11) =    ( CTTPP+CTPPT) * FACTOR
      BREON_VKERNEL(16) =    ( CTTPP-CTPPT) * FACTOR

!  Add the Diffuse term to R(1,1)
!  ------------------------------

!   IN THE PREVIOUS EQUATION, THE BRDF MODEL WE USE IS THE SINYUK-ET-AL MODEL
!   WITH FREE PARAMETERS rho,vk AND gt (see Sinyuk et al. paper).

!  This is just the Rahman Kernel.........different name !!

      XPHI_RAH = 180.0d0 - XPHI_REF
      CALL RAHMAN_VFUNCTION_PLUS                                    &
            ( MAXPARS, NPARS, PARS, DO_DERIV_PARS, N_BRDF_STOKESSQ, &
              XJ, SXJ, XI, SXI, XPHI_RAH, CKPHI_REF, SKPHI_REF,     &
              RAHMAN_VKERNEL, BREON_VDERIVATIVES )

!  Add to the specular term

      BREON_VKERNEL(1) = BREON_VKERNEL(1) + RAHMAN_VKERNEL(1)

!  Here is the original code from France........
!    cs: cosinus of the solar zenith angle
!    cv: cosinus of the viewing zenith angle
!    phi: azimuth between the solar and observation vertical planes
!      xx=cs**(vk-1.)
!      yy=cs**(vk-1.)
!      zz=(cs+cv)**(1.-vk)
!      FF1=rho*xx*yy/zz
!      xx=sqrt(1-cs*cs)
!      yy=sqrt(1-cv*cv)
!      ww=cs*cv-xx*yy*cos(pi*phi/180.)
!      aa=1+gt*gt+2*gt*ww
!      FF2=(1-gt*gt)/(aa**1.5)
!      vv=xx/cs
!      ww=yy/cv
!      G=sqrt(vv*vv+ww*ww+2*vv*ww*cos(pi*phi/180))
!      FF3=1+(1-rho)/(1+G)
!      RBD=FF1*FF2*FF3
! RBD...........................IS THE REFLECTION COEFFICIENT

!  Finish

      RETURN
END SUBROUTINE BREON_VFUNCTION_PLUS
