module L_IQGreenfunc_m

   USE VLIDORT_PARS

   public

   contains

SUBROUTINE VLIDORT_L_IQGF_SOLUTION (                  &
        LAYER, FOURIER, IBEAM, NSTOKES, NSTREAMS,     &
        NMOMENTS, FLUX_FACTOR, DO_VARY, N_PARAMETERS, &
        DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF, DFLUX, &
        QUAD_STRMWTS, QUAD_WEIGHTS, L_OMEGA_GREEK,    &
        PI_XQP, PI_X0P, PI_XQM_POST, K_REAL, SOLA_XPOS, L_SOLA_XPOS,  &
        NORM_SAVED, DMI, DPI, ATERM_SAVE, BTERM_SAVE, & ! Input
        L_ATERM_SAVE, L_BTERM_SAVE )       ! Output

!  Linearization of the Green's function solution (non-multipliers)

!  module, dimensions and numbers

      USE VLIDORT_PARS

      IMPLICIT NONE

!  subroutine input arguments
!  ==========================

      INTEGER, INTENT (IN) ::           LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR

      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS

      INTEGER, INTENT (IN) ::  LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::  DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  L_SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

!  Saved quantities for the Green function solution

      DOUBLE PRECISION, intent(in)  :: NORM_SAVED(MAXLAYERS,MAXSTREAMS_2)
      DOUBLE PRECISION, intent(in)  :: ATERM_SAVE(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: BTERM_SAVE(MAXSTREAMS_2,MAXLAYERS)
      DOUBLE PRECISION, intent(in)  :: DMI(MAXSTREAMS,2), DPI(MAXSTREAMS,2)

!  subroutine output arguments
!  ===========================

!  Linearized Saved quantities for the Green function solution

      DOUBLE PRECISION, intent(inout) :: L_ATERM_SAVE(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)
      DOUBLE PRECISION, intent(inout) :: L_BTERM_SAVE(MAXSTREAMS_2,MAXLAYERS,MAX_ATMOSWFS)

!  Local variables
!  ---------------

!  linearizations of component variables

      DOUBLE PRECISION :: L_DMI(MAXSTREAMS,2,MAX_ATMOSWFS)
      DOUBLE PRECISION :: L_DPI(MAXSTREAMS,2,MAX_ATMOSWFS)

!  help variables

      LOGICAL   :: DO_FIRST
      INTEGER   :: I, J, I1, J1, L, N, M, K, AA, O1, O2, O3, Q, IB
      INTEGER   :: NSTKS_NSTRMS_IQ,NSTOKES_IQ
      LOGICAL   :: AAMASK(MAXEVALUES)

      DOUBLE PRECISION :: L_SUM_LA, L_SUM_LB, L_NORM, L_ATERM, L_BTERM, TMA, TPA, SUM, GSUM
      DOUBLE PRECISION :: TA1, TA2, TB1, TB2, F1, NORM, S_TPA, S_TMA, T1, T2
      DOUBLE PRECISION :: L_HELP_Q1(0:MAXMOMENTS,MAXSTOKES)
!  Linearized norms for the Green function solution

      DOUBLE PRECISION  :: L_NORM_SAVED(MAXSTREAMS_2,MAX_ATMOSWFS)

!  initialise indices

      N  = LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  Check existence
!  ===============

!  This section added by R. Spurr, RTSOLUTIONS Inc., 3/4/06.

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. &
                      DO_LAYER_SCATTERING(M,N)

!  Filter

      NSTOKES_IQ = 1   ;   IF ( NSTOKES.GT.1) NSTOKES_IQ = 2
      NSTKS_NSTRMS_IQ = NSTREAMS * NSTOKES_IQ

!  second shot. 3 July 2012, Works.
      aamask = .false.
      aamask(1:NSTOKES_IQ*NSTREAMS) = .true.

!  If there is nothing varying or if there is no solution then
!        Zero the output values and exit

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) THEN
        DO AA = 1, NSTREAMS * NSTOKES_IQ
          DO Q = 1, N_PARAMETERS
            L_ATERM_SAVE(AA,N,Q) = ZERO
            L_BTERM_SAVE(AA,N,Q) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Form quantities independent of optical depth
!  ============================================

!  For all layers, save the norms

      IF ( DO_VARY ) THEN
        DO Q = 1, N_PARAMETERS
          AA = 0
          DO K = 1, K_REAL(N)
            IF ( AAMASK(K) ) then
              AA = AA + 1 ; NORM = ZERO
              DO J = 1, NSTREAMS
                J1 = J + NSTREAMS
                T1 = ZERO ; T2 = ZERO
                DO O1 = 1, NSTOKES_IQ
                  T1 = T1 + L_SOLA_XPOS(J,O1,K,N,Q)  * SOLA_XPOS(J,O1,K,N) 
                  T2 = T2 + L_SOLA_XPOS(J1,O1,K,N,Q) * SOLA_XPOS(J1,O1,K,N) 
                ENDDO
                NORM = NORM + QUAD_STRMWTS(J)*(T1-T2)
              ENDDO
              L_NORM_SAVED(AA,Q) = TWO * NORM
            ENDIF
          ENDDO
        ENDDO
      ENDIF

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  first function
!    (linearization of primary scattering of beam)

        IF ( DO_VARY ) THEN
          DO L = M, NMOMENTS
            DO O1 = 1, NSTOKES_IQ
              SUM = ZERO
              DO O2 = 1, NSTOKES_IQ
                GSUM = ZERO
                DO O3 = 1, NSTOKES_IQ
                  GSUM = GSUM + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
                ENDDO
                SUM = SUM + L_OMEGA_GREEK(L,N,O1,O2,Q) * GSUM
              ENDDO
              L_HELP_Q1(L,O1) = SUM
            ENDDO
          ENDDO
        ENDIF

!  set up linearizations of help arrays (independent of eigenvector)

        DO I = 1, NSTREAMS
          DO O1 = 1, NSTOKES_IQ
            S_TPA = ZERO
            S_TMA = ZERO
            DO L = M, NMOMENTS
              TPA = ZERO
              TMA = ZERO
              DO O2 = 1, NSTOKES_IQ
                TPA = TPA + PI_XQP(L,I,O1,O2)      * L_HELP_Q1(L,O2)
                TMA = TMA + PI_XQM_POST(L,I,O1,O2) * L_HELP_Q1(L,O2)
              ENDDO
              S_TPA = S_TPA + TPA
              S_TMA = S_TMA + TMA
            ENDDO
            L_DPI(I,O1,Q) = F1 * S_TPA
            L_DMI(I,O1,Q) = F1 * S_TMA
          ENDDO
        ENDDO

!  linearize quantities independent of TAU (L_ATERM_SAVE, L_BTERM_SAVE)

        AA = 0
        DO K = 1, K_REAL(N)
          IF ( AAMASK(K) ) then
            AA = AA + 1 
            L_SUM_LA = ZERO
            L_SUM_LB = ZERO
            DO I = 1, NSTREAMS
              I1 = I + NSTREAMS
              TA1 = ZERO ; TA2 = ZERO
              TB1 = ZERO ; TB2 = ZERO
              DO O1 = 1, NSTOKES
                TA1 = TA1 + L_DPI(I,O1,Q) *   SOLA_XPOS(I,O1,K,N) &
                          +   DPI(I,O1)   * L_SOLA_XPOS(I,O1,K,N,Q) 
                TA2 = TA2 + L_DMI(I,O1,Q) *   SOLA_XPOS(I1,O1,K,N) &
                          +   DMI(I,O1)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
                TB1 = TB1 + L_DMI(I,O1,Q) *   SOLA_XPOS(I,O1,K,N) &
                          +   DMI(I,O1)   * L_SOLA_XPOS(I,O1,K,N,Q) 
                TB2 = TB2 + L_DPI(I,O1,Q) *   SOLA_XPOS(I1,O1,K,N) &
                          +   DPI(I,O1)   * L_SOLA_XPOS(I1,O1,K,N,Q) 
              ENDDO
              L_SUM_LA  = L_SUM_LA + QUAD_WEIGHTS(I) * ( TA1 + TA2 )
              L_SUM_LB  = L_SUM_LB + QUAD_WEIGHTS(I) * ( TB1 + TB2 )
            ENDDO
            L_NORM = L_NORM_SAVED(AA,Q)
            L_ATERM = ( L_SUM_LA / ATERM_SAVE(AA,N) ) - L_NORM
            L_BTERM = ( L_SUM_LB / BTERM_SAVE(AA,N) ) - L_NORM
            L_ATERM_SAVE(AA,N,Q) = L_ATERM / NORM_SAVED(N,AA)
            L_BTERM_SAVE(AA,N,Q) = L_BTERM / NORM_SAVED(N,AA)
          ENDIF
        ENDDO

!   End parameter loop

      ENDDO

!  Finish

      RETURN
END SUBROUTINE VLIDORT_L_IQGF_SOLUTION

!

SUBROUTINE VLIDORT_L_IQGF_USERSOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, DO_VARY, N_PARAMETERS, &
        DO_UPWELLING, DO_DNWELLING, NSTOKES,                &
        DO_OBSERVATION_GEOMETRY, FLUX_FACTOR,               &
        LAYER_PIS_CUTOFF, NMOMENTS, N_USER_STREAMS,         &
        DO_LAYER_SCATTERING, LOCAL_UM_START,                &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN,             &
        DFLUX, L_OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P,       &
        L_UPAR_DN_1, L_UPAR_UP_1 )

      USE VLIDORT_PARS

      IMPLICIT NONE

!  Regular input (strictly in)

      INTEGER, INTENT (IN) ::          GIVEN_LAYER
      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          IBEAM
      LOGICAL, INTENT (IN) ::          DO_VARY
      INTEGER, INTENT (IN) ::          N_PARAMETERS
      LOGICAL, INTENT (IN) ::          DO_UPWELLING
      LOGICAL, INTENT (IN) ::          DO_DNWELLING
      INTEGER, INTENT (IN) ::          NSTOKES
      LOGICAL, INTENT (IN) ::          DO_OBSERVATION_GEOMETRY
      DOUBLE PRECISION, INTENT (IN) :: FLUX_FACTOR
      INTEGER, INTENT (IN) ::          LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::          NMOMENTS
      INTEGER, INTENT (IN) ::          N_USER_STREAMS
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::          LOCAL_UM_START
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::          STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) :: DFLUX ( MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) :: L_OMEGA_GREEK &
         ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  Linearized (Ostensibly OUTPUT)

      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION, INTENT (INOUT) :: L_UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER, DO_FIRST
      INTEGER ::          IB, UM, L, N, M, O1, O2, O3
      INTEGER ::          Q

      DOUBLE PRECISION :: SUM1, F1, T1, SP, SPS

      DOUBLE PRECISION :: L_HELP_Q1(0:MAXMOMENTS,MAXSTOKES)

!  Layer and Fourier

      N  = GIVEN_LAYER
      M  = FOURIER
      IB = IBEAM
      F1 = FLUX_FACTOR / PI4

!  Linearizations for quantities varying in Layer N

!  Check existence
!  ---------------

!  Only a solution if the layer is active and not below Cutoff.

      DO_FIRST = ( N .LE. LAYER_PIS_CUTOFF(IB) ) .AND. &
                      DO_LAYER_SCATTERING(M,N)

!  If no solution or no variation, zero output and go to Part 2

      IF ( .NOT. DO_VARY .OR. .NOT. DO_FIRST ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
          DO Q = 1, N_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                L_UPAR_UP_1(UM,O1,N,Q)       = ZERO
              ENDDO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO O1 = 1, NSTOKES
                L_UPAR_DN_1(UM,O1,N,Q)       = ZERO
              ENDDO
            ENDDO
           ENDIF
          ENDDO
        ELSE
          DO Q = 1, N_PARAMETERS
           IF ( DO_UPWELLING ) THEN
            DO O1 = 1, NSTOKES
              L_UPAR_UP_1(IB,O1,N,Q)        = ZERO
            ENDDO
           ENDIF
           IF ( DO_DNWELLING ) THEN
            DO O1 = 1, NSTOKES
              L_UPAR_DN_1(IB,O1,N,Q)        = ZERO
            ENDDO
           ENDIF
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. &
                   STERM_LAYERMASK_DN(N) )

!  start parameter loop

      DO Q = 1, N_PARAMETERS

!  first function
!    (linearization of primary scattering of beam)

       IF ( DO_LAYER ) THEN
        DO L = M, NMOMENTS
         DO O1 = 1, NSTOKES
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SP = ZERO
            DO O3 = 1, NSTOKES
              SP = SP + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SPS = SPS + L_OMEGA_GREEK(L,N,O1,O2,Q) * SP
          ENDDO
          L_HELP_Q1(L,O1) = SPS
         ENDDO
        ENDDO
       ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integrated contributions

       IF ( DO_UPWELLING .AND. STERM_LAYERMASK_UP(N) ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           T1 = ZERO
           DO L = M, NMOMENTS
            SUM1 = ZERO
            DO O2 = 1, NSTOKES
             SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUM(L,UM,O1,O2)
            ENDDO
            T1 = T1 + SUM1
           ENDDO
           L_UPAR_UP_1(UM,O1,N,Q)       = F1 * T1
          ENDDO
         ENDDO
        ELSE
          DO O1 = 1, NSTOKES
           T1 = ZERO
           DO L = M, NMOMENTS
            SUM1 = ZERO
            DO O2 = 1, NSTOKES
             SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUM(L,IB,O1,O2)
            ENDDO
            T1 = T1 + SUM1
           ENDDO
           L_UPAR_UP_1(IB,O1,N,Q)       = F1 * T1
          ENDDO
        ENDIF
       ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integrated contributions

       IF ( DO_DNWELLING .AND. STERM_LAYERMASK_DN(N) ) THEN
        IF (.not. DO_OBSERVATION_GEOMETRY ) THEN
         DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
           T1 = ZERO
           DO L = M, NMOMENTS
            SUM1 = ZERO
            DO O2 = 1, NSTOKES
             SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUP(L,UM,O1,O2)
            ENDDO
            T1 = T1 + SUM1
           ENDDO
           L_UPAR_DN_1(UM,O1,N,Q)       = F1 * T1
          ENDDO
         ENDDO
        ELSE
          DO O1 = 1, NSTOKES
           T1 = ZERO
           DO L = M, NMOMENTS
            SUM1 = ZERO
            DO O2 = 1, NSTOKES
             SUM1 = SUM1 + L_HELP_Q1(L,O2)*PI_XUP(L,IB,O1,O2)
            ENDDO
            T1 = T1 + SUM1
           ENDDO
           L_UPAR_DN_1(IB,O1,N,Q)       = F1 * T1
          ENDDO
        ENDIF
       ENDIF

!  end parameter loop

      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_L_IQGF_USERSOLUTION


!  End

end module L_IQGreenfunc_m
