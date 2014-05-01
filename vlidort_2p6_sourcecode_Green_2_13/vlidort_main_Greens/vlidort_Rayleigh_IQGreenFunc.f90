module IQGreenfunc_m

   USE VLIDORT_PARS

   public

   contains

      SUBROUTINE VLIDORT_IQGF_SOLUTION (              &
        GIVEN_LAYER, FOURIER, IBEAM,                  &
        NSTOKES, NSTREAMS, NMOMENTS, FLUX_FACTOR,     &
        DO_LAYER_SCATTERING, LAYER_PIS_CUTOFF, DFLUX, &
        QUAD_STRMWTS, QUAD_WEIGHTS, OMEGA_GREEK,      &
        INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,  &
        PI_XQP, PI_X0P, PI_XQM_POST,                  &
        K_REAL, KEIGEN, SOLA_XPOS, T_DELT_EIGEN,      &
        NORM_SAVED, DPI, DMI, ATERM_SAVE, BTERM_SAVE, &
        GAMMA_M, GAMMA_P, GFUNC_UP, GFUNC_DN,         &
        WUPPER, WLOWER, STATUS, MESSAGE, TRACE )

!  This is the Green's function solution for IQ only

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM

      INTEGER, INTENT (IN) ::           NSTOKES
      INTEGER, INTENT (IN) ::           NSTREAMS
      INTEGER, INTENT (IN) ::           NMOMENTS
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR

      INTEGER, INTENT (IN) ::  LAYER_PIS_CUTOFF ( MAXBEAMS )
      LOGICAL, INTENT (IN) ::  DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) ::  QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  QUAD_WEIGHTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_MUBAR ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  INITIAL_TRANS ( MAXLAYERS, MAXBEAMS )
      DOUBLE PRECISION, INTENT (IN) ::  AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )

      DOUBLE PRECISION, INTENT (IN) ::  PI_XQP &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XQM_POST &
          ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

      INTEGER, INTENT (IN) ::           K_REAL ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  KEIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  SOLA_XPOS &
          ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  T_DELT_EIGEN &
          ( MAXEVALUES, MAXLAYERS )

!  Output
!  ======

!  Saved quantities for the Green function solution

      double precision, INTENT (INOUT) :: ATERM_SAVE(MAXSTREAMS_2,MAXLAYERS)
      double precision, INTENT (INOUT) :: BTERM_SAVE(MAXSTREAMS_2,MAXLAYERS)

!  Holding arrays for Multiplier coefficients

      double precision, INTENT (INOUT) :: GAMMA_M(MAXSTREAMS_2,MAXLAYERS)
      double precision, INTENT (INOUT) :: GAMMA_P(MAXSTREAMS_2,MAXLAYERS)

!  Green's function GFUNC multipliers

      double precision, INTENT (INOUT) :: GFUNC_UP(MAXSTREAMS_2,MAXLAYERS)
      double precision, INTENT (INOUT) :: GFUNC_DN(MAXSTREAMS_2,MAXLAYERS)

!  Norm stuff (need for linearization)

      double precision, INTENT (INOUT) :: DMI(MAXSTREAMS,2), DPI(MAXSTREAMS,2)
      double precision, INTENT (INOUT) :: NORM_SAVED(MAXLAYERS,MAXSTREAMS_2)

!  Particular integrals at layer boundaries

      DOUBLE PRECISION, INTENT (INOUT) ::  WUPPER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) ::  WLOWER &
          ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

!  Excpetion handling

      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (INOUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (INOUT) :: TRACE

!  Local variables
!  ---------------

!  help variables

      INTEGER ::          I, J, I1, J1, L, N, M, K, AA, O1, O2, O3
      INTEGER ::          NSTKS_NSTRMS_IQ,NSTOKES_IQ
      LOGICAL ::          AAMASK(MAXEVALUES)

      DOUBLE PRECISION :: RHO_M, RHO_P, SUM, GSUM, F1, S_P_U, S_P_L, S_M_U, S_M_L
      DOUBLE PRECISION :: TPA, TMA, S_TPA, S_TMA, SUM_LA, SUM_LB, del, RHOM
      DOUBLE PRECISION :: CONST, SECBAR, WDEL, ZDEL, ZWDEL, NORM, T1, T2
      DOUBLE PRECISION :: HELP_QFUNC ( MAXLAYERS, 0:MAXMOMENTS, MAXSTOKES )
      double precision :: CFUNC(MAXSTREAMS_2)
      double precision :: DFUNC(MAXSTREAMS_2)

!  Initial section
!  ---------------

!  Initialize status

      STATUS  = VLIDORT_SUCCESS
      MESSAGE = ' '
      TRACE   = ' '

!  initialise indices

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No result if conditions not right

      IF ( M.ne.0.and.NSTOKES.eq.4) then
        status = vlidort_serious
        message = 'Wrong conditions for IQ_Greensfunc calculation'
        trace = 'Beginning of routine VLIDORT_IQGREENSFUNC_SOLUTION'
        return
      endif

!  No particular solution beyond the cutoff layer
!  OR No scattering in this layer then:
!    ---> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
             .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        DO I = 1, 2*NSTREAMS
          DO O1 = 1, NSTOKES
            WUPPER(I,O1,N) = ZERO
            WLOWER(I,O1,N) = ZERO
          ENDDO
        ENDDO
        RETURN
      ENDIF

!  Filter

      NSTOKES_IQ = 1   ;   IF ( NSTOKES.GT.1) NSTOKES_IQ = 2
      NSTKS_NSTRMS_IQ = NSTREAMS * NSTOKES_IQ

!  second shot. 3 July 2012, Works.
      aamask = .false.
      aamask(1:NSTOKES_IQ*NSTREAMS) = .true.

!  First shot
!      aamask(1:NSTOKES*NSTREAMS) = .true.
!      IF ( NSTOKES.eq.3) then
!         do i = 1, nstreams
!            aamask(3*i) = .false.
!         enddo
!      ENDIF

!  Norms

      AA = 0
      DO K = 1, K_REAL(N)
        IF ( AAMASK(K) ) then
          AA = AA + 1 ; NORM = ZERO
          DO J = 1, NSTREAMS
            J1 = J + NSTREAMS
            T1 = ZERO ; T2 = ZERO
            DO O1 = 1, NSTOKES_IQ
              T1 = T1 + SOLA_XPOS(J,O1,K,N)  * SOLA_XPOS(J,O1,K,N) 
              T2 = T2 + SOLA_XPOS(J1,O1,K,N) * SOLA_XPOS(J1,O1,K,N) 
            ENDDO
            NORM = NORM + QUAD_STRMWTS(J)*(T1-T2)
          ENDDO
          NORM_SAVED(N,AA) = NORM
        ENDIF
      ENDDO

!  constants for the layer

      SECBAR = AVERAGE_SECANT(N,IBEAM)
      CONST  = INITIAL_TRANS(N,IBEAM)

!  Gamma constants

      AA = 0
      DO K = 1, K_REAL(N)
        IF ( AAMASK(K) ) then
          AA = AA + 1 
          RHO_M = SECBAR - KEIGEN(K,N)
          RHO_P = SECBAR + KEIGEN(K,N)
          GAMMA_P(AA,N) = ONE / RHO_P
          GAMMA_M(AA,N) = ONE / RHO_M
!          if (n.eq.62)write(*,*)N,K,AA,SECBAR,KEIGEN(K,N),GAMMA_M(AA,N)
        ENDIF
      ENDDO

!  Optical depth integrations for the discrete ordinate solution

      AA = 0
      WDEL    = T_DELT_MUBAR(N,IBEAM)
      DO K = 1, K_REAL(N)
        IF ( AAMASK(K) ) then
          AA = AA + 1 ; rhom = one / GAMMA_M(AA,N)
          ZDEL  = T_DELT_EIGEN(K,N)
          ZWDEL = ZDEL * WDEL
          if ( ABS(RHOM).lt.1.0d-03 ) then
            del = -LOG(WDEL)/SECBAR
            CALL LIMIT_GCFUNC  ( RHOM, DEL, ZDEL, CFUNC(AA) )

!            write(776,*)n,aa,RHOM,ZDEL,del,CFUNC(AA),( ZDEL - WDEL ) * GAMMA_M(AA,N)

          else
             CFUNC(AA)  = ( ZDEL - WDEL ) * GAMMA_M(AA,N)
          endif
          DFUNC(AA)  = ( ONE - ZWDEL ) * GAMMA_P(AA,N)
        ENDIF
      ENDDO

!  Set up sum and difference vectors for Beam source terms
!  ( sum vector may be required again in linearization )
!  Auxiliary matrix for Q functions

      DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES_IQ
          SUM = ZERO
          DO O2 = 1, NSTOKES_IQ
            GSUM = ZERO
            DO O3 = 1, NSTOKES_IQ
              GSUM = GSUM + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SUM = SUM + OMEGA_GREEK(L,N,O1,O2) * GSUM
          ENDDO
          HELP_QFUNC(N,L,O1) = SUM
        ENDDO
      ENDDO

      DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES_IQ
          S_TPA = ZERO
          S_TMA = ZERO
          DO L = M, NMOMENTS
            TPA = ZERO
            TMA = ZERO
            DO O2 = 1, NSTOKES_IQ
              TPA = TPA + PI_XQP(L,I,O1,O2)      * HELP_QFUNC(N,L,O2)
              TMA = TMA + PI_XQM_POST(L,I,O1,O2) * HELP_QFUNC(N,L,O2)
            ENDDO
            S_TPA = S_TPA + TPA
            S_TMA = S_TMA + TMA
          ENDDO
          DPI(I,O1) = F1 * S_TPA
          DMI(I,O1) = F1 * S_TMA
        ENDDO
      ENDDO

!  For each eigenstream, get the terms ATERM_SAVE and BTERM_SAVE

      AA = 0
      DO K = 1, K_REAL(N)
        IF ( AAMASK(K) ) then
          AA = AA + 1 
          SUM_LA = ZERO
          SUM_LB = ZERO
          DO I = 1, NSTREAMS
            I1 = I + NSTREAMS
            TPA = ZERO ; TMA = ZERO
            DO O1 = 1, NSTOKES
              TPA = TPA + DPI(I,O1) * SOLA_XPOS(I,O1,K,N) &
                        + DMI(I,O1) * SOLA_XPOS(I1,O1,K,N)
              TMA = TMA + DMI(I,O1) * SOLA_XPOS(I,O1,K,N) &
                        + DPI(I,O1) * SOLA_XPOS(I1,O1,K,N)
            ENDDO
            SUM_LA  = SUM_LA + QUAD_WEIGHTS(I)*TPA
            SUM_LB  = SUM_LB + QUAD_WEIGHTS(I)*TMA
          ENDDO
          ATERM_SAVE(AA,N) = SUM_LA / NORM_SAVED(N,AA)
          BTERM_SAVE(AA,N) = SUM_LB / NORM_SAVED(N,AA)
        ENDIF
      ENDDO

!  Green function multipliers For each eigenstream

      DO AA = 1, NSTKS_NSTRMS_IQ
        GFUNC_DN(AA,N) = CFUNC(AA) * ATERM_SAVE(AA,N) * CONST
        GFUNC_UP(AA,N) = DFUNC(AA) * BTERM_SAVE(AA,N) * CONST
      ENDDO

!  Set particular integral from Green function expansion

      DO I = 1, NSTREAMS
        I1 = I + NSTREAMS
        DO O1 = 1, NSTOKES
          S_P_U = ZERO
          S_P_L = ZERO
          S_M_U = ZERO
          S_M_L = ZERO
          AA = 0
          DO K = 1, K_REAL(N)
            IF ( AAMASK(K) ) then
              AA = AA + 1 
              S_P_U = S_P_U + GFUNC_UP(AA,N)*SOLA_XPOS(I1,O1,K,N)
              S_M_U = S_M_U + GFUNC_UP(AA,N)*SOLA_XPOS(I,O1,K,N)
              S_P_L = S_P_L + GFUNC_DN(AA,N)*SOLA_XPOS(I,O1,K,N)
              S_M_L = S_M_L + GFUNC_DN(AA,N)*SOLA_XPOS(I1,O1,K,N)
            ENDIF
          ENDDO
          WUPPER(I,O1,N)  = S_P_U
          WUPPER(I1,O1,N) = S_M_U
          WLOWER(I1,O1,N) = S_M_L
          WLOWER(I,O1,N)  = S_P_L
        ENDDO
      ENDDO

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_IQGF_SOLUTION

!

      SUBROUTINE VLIDORT_IQGF_USERSOLUTION ( &
        GIVEN_LAYER, FOURIER, IBEAM, &
        DO_UPWELLING, DO_DNWELLING, &
        DO_OBSERVATION_GEOMETRY, &
        NSTOKES, FLUX_FACTOR, &
        LAYER_PIS_CUTOFF, NMOMENTS, N_USER_STREAMS, &
        DO_LAYER_SCATTERING, LOCAL_UM_START, &
        STERM_LAYERMASK_UP, STERM_LAYERMASK_DN, DFLUX, &
        OMEGA_GREEK, PI_XUP, PI_XUM, PI_X0P, &
        UPAR_DN_1, UPAR_UP_1 )

      USE VLIDORT_PARS

      IMPLICIT NONE

      INTEGER, INTENT (IN) ::           GIVEN_LAYER
      INTEGER, INTENT (IN) ::           FOURIER
      INTEGER, INTENT (IN) ::           IBEAM
      LOGICAL, INTENT (IN) ::           DO_OBSERVATION_GEOMETRY
      LOGICAL, INTENT (IN) ::           DO_UPWELLING
      LOGICAL, INTENT (IN) ::           DO_DNWELLING
      INTEGER, INTENT (IN) ::           NSTOKES
      DOUBLE PRECISION, INTENT (IN) ::  FLUX_FACTOR
      INTEGER, INTENT (IN) ::           LAYER_PIS_CUTOFF ( MAXBEAMS )
      INTEGER, INTENT (IN) ::           NMOMENTS
      INTEGER, INTENT (IN) ::           N_USER_STREAMS
      LOGICAL, INTENT (IN) ::           DO_LAYER_SCATTERING &
          ( 0:MAXMOMENTS, MAXLAYERS )
      INTEGER, INTENT (IN) ::           LOCAL_UM_START
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_UP ( MAXLAYERS )
      LOGICAL, INTENT (IN) ::           STERM_LAYERMASK_DN ( MAXLAYERS )
      DOUBLE PRECISION, INTENT (IN) ::  DFLUX ( MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  OMEGA_GREEK &
          ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (IN) ::  PI_XUP &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_XUM &
          ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) ::  PI_X0P &
          ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_DN_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION, INTENT (INOUT) :: UPAR_UP_1 &
          ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )

!  Local variables
!  ---------------

      LOGICAL ::          DO_LAYER
      INTEGER ::          UM, L, N, M, O1, O2, O3
      DOUBLE PRECISION :: SUM1, T1, F1, SP, SPS
      DOUBLE PRECISION :: HELP_Q1 ( 0:MAXMOMENTS, MAXSTOKES )

!  Layer and Fourier

      N = GIVEN_LAYER
      M = FOURIER
      F1 = FLUX_FACTOR / PI4

!  No particular solution beyond the cutoff layer.
!  OR No scattering in this layer
!    --> Zero the boundary layer values and exit

      IF ( N .GT. LAYER_PIS_CUTOFF(IBEAM) .OR. &
             .NOT. DO_LAYER_SCATTERING(M,N) ) THEN
        IF ( DO_UPWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              UPAR_UP_1(UM,O1,N) = ZERO
            ENDDO
          ENDDO
        ENDIF
        IF ( DO_DNWELLING ) THEN
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              UPAR_DN_1(UM,O1,N) = ZERO
            ENDDO
          ENDDO
        ENDIF
        RETURN
      ENDIF

!  existence flag

      DO_LAYER = ( STERM_LAYERMASK_UP(N) .OR. &
                   STERM_LAYERMASK_DN(N) )

!  first function

      IF ( DO_LAYER ) THEN
       DO L = M, NMOMENTS
        DO O1 = 1, NSTOKES
          SPS = ZERO
          DO O2 = 1, NSTOKES
            SP = ZERO
            DO O3 = 1, NSTOKES
              SP = SP + PI_X0P(L,IBEAM,N,O2,O3) * DFLUX(O3)
            ENDDO
            SPS = SPS + OMEGA_GREEK(L,N,O1,O2) * SP
          ENDDO
          HELP_Q1(L,O1) = F1 * SPS
        ENDDO
       ENDDO
      ENDIF

!  Now sum over all harmonic contributions (Upwelling)
!    Direct and integarated contributions

      IF ( DO_UPWELLING ) THEN
       IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUM(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
          ENDDO
          UPAR_UP_1(UM,O1,N) = T1
         ENDDO
        ENDDO
       ELSE
         DO O1 = 1, NSTOKES
          T1 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUM(L,IBEAM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
          ENDDO
          UPAR_UP_1(IBEAM,O1,N) = T1
         ENDDO
       ENDIF
      ENDIF

!  Now sum over all harmonic contributions (Downwelling)
!    Direct and integarated contributions

      IF ( DO_DNWELLING ) THEN
       IF ( .not. DO_OBSERVATION_GEOMETRY ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO O1 = 1, NSTOKES
          T1 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUP(L,UM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
          ENDDO
          UPAR_DN_1(UM,O1,N) = T1
         ENDDO
        ENDDO
       ELSE
         DO O1 = 1, NSTOKES
          T1 = ZERO
          DO L = M, NMOMENTS
           SUM1 = ZERO
           DO O2 = 1, NSTOKES
            SUM1 = SUM1 + HELP_Q1(L,O2)*PI_XUP(L,IBEAM,O1,O2)
           ENDDO
           T1 = T1 + SUM1
          ENDDO
          UPAR_DN_1(IBEAM,O1,N) = T1
         ENDDO
       ENDIF
      ENDIF

!  debug stokes = 1

!      DO UM = LOCAL_UM_START, N_USER_STREAMS
!        write(*,*)n,um,UPAR_UP_1(UM,1,N),UPAR_UP_2(UM,1,N)
!      ENDDO
!      DO UM = LOCAL_UM_START, N_USER_STREAMS
!        write(*,*)n,um,UPAR_DN_1(UM,1,N),UPAR_DN_2(UM,1,N)
!      ENDDO
!      PAUSE

!  Finish

      RETURN
      END SUBROUTINE VLIDORT_IQGF_USERSOLUTION

      SUBROUTINE LIMIT_GCFUNC &
          ( EPS, DELTA, ZDEL, CFUNC )

!  Green's function solution multiplier CFUNC (whole layer)
!    Small number expansion to second order

      IMPLICIT NONE

!  input arguments

      DOUBLE PRECISION, INTENT(IN) :: eps, delta, zdel

!  output arguments

      DOUBLE PRECISION, INTENT(OUT) :: cfunc

!  local declarations

      DOUBLE PRECISION :: power, power2

!  initialise

      cfunc = 0.0d0

!  evaluate

      power = delta * eps
      power2 = power * power
      cfunc = zdel*delta*(1.0d0-0.5d0*power+power2/6.0d0)

!  Finish

      return
      END SUBROUTINE LIMIT_GCFUNC

end module IQGreenfunc_m

