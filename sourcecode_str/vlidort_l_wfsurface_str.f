C ###############################################################
C #                                                             #
C #                    THE VLIDORT  MODEL                       #
C #                                                             #
C #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
C #  -          --         -        -        -         -        #
C #                                                             #
C ###############################################################

C ###############################################################
C #                                                             #
C #  Author :      Robert. J. D. Spurr                          #
C #                                                             #
C #  Address :      RT Solutions, Inc.                          #
C #            9 Channing Street                                #
C #             Cambridge, MA 02138, USA                        #
C #            Tel: (617) 492 1183                              #
C #                                                             #
C #  Email :      rtsolutions@verizon.net                       #
C #                                                             #
C #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT           #
C #  Release Date :   December 2005  (2.0)                      #
C #  Release Date :   March 2007     (2.2)                      #
C #  Release Date :   October 2007   (2.3)                      #
C #  Release Date :   December 2008  (2.4)                      #
C #  Release Date :   April/May 2009 (2.4R)                     #
C #  Release Date :   July 2009      (2.4RT)                    #
C #                                                             #
C #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
C #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
C #       NEW: Thermal Emission Treatment     (2.4RT)           #
C #                                                             #
C ###############################################################

C    #####################################################
C    #                                                   #
C    #   This Version of VLIDORT comes with a GNU-style  #
C    #   license. Please read the license carefully.     #
C    #                                                   #
C    #####################################################

C ##########################################################
C #                                                        #
C # Subroutines in this Module                             #
C #                                                        #
C #     Top level routines--------------                   #
C #            VLIDORT_SURFACE_WFS  (master)               #
C #                                                        #
C #     Linearized BVP Column, surface WFs ---------       #
C #            LS_BVP_COLUMN_SETUP                         #
C #                                                        #
C #     BOA surface source terms ---------                 #
C #            GET_LS_BOA_SOURCE                           #
C #                                                        #
C #     Recursion relations ---------                      #
C #            UPUSER_SURFACEWF                            #
C #            DNUSER_SURFACEWF                            #
C #                                                        #
C #     Post-processing at user angles --------            #
C #            LS_WHOLELAYER_STERM_UP                      #
C #            LS_WHOLELAYER_STERM_DN                      #
C #            LS_PARTLAYER_STERM_UP                       #
C #            LS_PARTLAYER_STERM_DN                       #
C #                                                        #
C ##########################################################

      SUBROUTINE VLIDORT_SURFACE_WFS
     I       (  DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          N_SURFACE_WFS, FOURIER_COMPONENT, IBEAM,
     I          SURFACE_FACTOR, FLUX_MULTIPLIER,
     O          STATUS )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include file of solution variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  Module arguments
C  ----------------

C  local control flags

      LOGICAL          DO_INCLUDE_DIRECTBEAM
      LOGICAL          DO_INCLUDE_SURFEMISS
      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Number of weighting functions

      INTEGER          N_SURFACE_WFS

C  Fourier surface factor, Fourier number, Beam index, Flux multiplier

      DOUBLE PRECISION SURFACE_FACTOR
      INTEGER          FOURIER_COMPONENT, IBEAM
      DOUBLE PRECISION FLUX_MULTIPLIER

C  Output status

      INTEGER          STATUS

C  Local variables
C  ---------------

C  error tracing variables

      INTEGER          INFO
      CHARACTER*3      CI
      CHARACTER*70     MAIL, TRACE

C  Linearized BOA terms

      DOUBLE PRECISION LS_BOA_SOURCE
     &       ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE
     &       ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES)

C  Other local variables

      INTEGER          N, Q, K, K0, K1, K2, KO1, C0
      INTEGER          IROW, IROW1, IROW_S, IROW1_S 

C  Initialise status

       STATUS = VLIDORT_SUCCESS

C  Regular BVP Solution
C  ====================

C   NO TELESCOPING HERE

C  BV solution for perturbed integration constants
C  -----------------------------------------------

C  Compute the main column B' where AX = B'

C       write(*,*)DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, 
C     I      FOURIER_COMPONENT, N_SURFACE_WFS 

      CALL LS_BVP_COLUMN_SETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, 
     I      SURFACE_FACTOR, FOURIER_COMPONENT, N_SURFACE_WFS, IBEAM )

C  BVP back-substitution: With compression (multilayers)
C  -----------------------------------------------------

      IF ( NLAYERS .GT. 1 ) THEN

C  LAPACK substitution (DGBTRS) using RHS column vector COL2_WF
C  BV solution for perturbed integration constants
C    ( call to LAPACK solver routine for back substitution )

        CALL DGBTRS
     &     ( 'n', NTOTAL, N_SUBDIAG, N_SUPDIAG, N_SURFACE_WFS,
     &        BANDMAT2, MAXBANDTOTAL, IPIVOT,
     &        COL2_WFALB, MAXTOTAL, INFO )

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (multilayer) in VLIDORT_SURFACE WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, all layers

        DO Q = 1, N_SURFACE_WFS
          DO N = 1, NLAYERS
            C0 = (N-1)*NSTKS_NSTRMS_2
            DO K = 1, K_REAL(N)
              IROW = K
              IROW1 = IROW + NSTKS_NSTRMS
              NCON_ALB(Q,K,N) = COL2_WFALB(C0+IROW,Q)
              PCON_ALB(Q,K,N) = COL2_WFALB(C0+IROW1,Q)
            ENDDO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2 * K - 2
              K1 = KO1 + K0
              K2 = K1  + 1
              IROW    = K + K_REAL(N)
              IROW1   = IROW + NSTKS_NSTRMS
              IROW_S  = IROW + K_COMPLEX(N)
              IROW1_S = IROW_S + NSTKS_NSTRMS
              NCON_ALB(Q,K1,N) = COL2_WFALB(C0+IROW,   Q)
              NCON_ALB(Q,K2,N) = COL2_WFALB(C0+IROW_S, Q)
              PCON_ALB(Q,K1,N) = COL2_WFALB(C0+IROW1,  Q)
              PCON_ALB(Q,K2,N) = COL2_WFALB(C0+IROW1_S,Q)
            ENDDO
          ENDDO
        ENDDO

C  Solve the boundary problem: No compression, Single Layer only
C  -------------------------------------------------------------

      ELSE IF ( NLAYERS .EQ. 1 ) THEN

C  LAPACK substitution (DGETRS) using RHS column vector SCOL2_WFALB

        CALL DGETRS
     &     ( 'N', NTOTAL, N_SURFACE_WFS, SMAT2, MAXSTRMSTKS_2, 
     &        SIPIVOT, SCOL2_WFALB, MAXSTRMSTKS_2, INFO )

C  (error tracing)

        IF ( INFO .LT. 0 ) THEN
          WRITE(CI, '(I3)' ) INFO
          MAIL  = 'argument i illegal value, for i = '//CI
          TRACE = 'DGBTRS call (Reg. 1 layer) in VLIDORT_SURFACE_WFS'
          STATUS = VLIDORT_SERIOUS
          CALL VLIDORT_ERROR_TRACE ( MAIL, TRACE, STATUS )
          RETURN
        ENDIF

C  Set Linearized integration constants NCON_ALB and PCON_ALB, 1 layer

        DO Q = 1, N_SURFACE_WFS
          N = 1
          DO K = 1, K_REAL(N)
            IROW = K
            IROW1 = IROW + NSTKS_NSTRMS
            NCON_ALB(Q,K,N) = SCOL2_WFALB(IROW,Q)
            PCON_ALB(Q,K,N) = SCOL2_WFALB(IROW1,Q)
          ENDDO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            IROW    = K + K_REAL(N)
            IROW1   = IROW + NSTKS_NSTRMS
            IROW_S  = IROW + K_COMPLEX(N)
            IROW1_S = IROW_S + NSTKS_NSTRMS
            NCON_ALB(Q,K1,N) = SCOL2_WFALB(IROW,   Q)
            NCON_ALB(Q,K2,N) = SCOL2_WFALB(IROW_S, Q)
            PCON_ALB(Q,K1,N) = SCOL2_WFALB(IROW1,  Q)
            PCON_ALB(Q,K2,N) = SCOL2_WFALB(IROW1_S,Q)
          ENDDO
        ENDDO

C  end clause
 
      ENDIF

C  debug------------------------------------------
c        if ( do_debug_write.and.fourier_component.eq.0 ) then
c         DO N = 1, NLAYERS
c          DO K = 1, K_REAL(N)
c           write(86,'(3i2,1p6e13.5)')FOURIER_COMPONENT,N,K,
c     &                LCON(K,N), MCON(K,N),
c     &                NCON_ALB(1,K,N),PCON_ALB(1,K,N),
c     &                NCON_ALB(1,K,N),PCON_ALB(1,K,N)
c          ENDDO
c         ENDDO
c        ENDIF


C  Get the Post-processed weighting functions
C  ==========================================

C  Upwelling weighting functions
C  -----------------------------

      IF ( DO_UPWELLING ) THEN

C  Get the surface term (L_BOA_SOURCE). External Function

        CALL GET_LS_BOA_SOURCE
     I        ( DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, N_SURFACE_WFS,
     I          FOURIER_COMPONENT, IBEAM, 
     O          LS_BOA_SOURCE,
     O          LS_BOA_THTONLY_SOURCE )

        CALL UPUSER_SURFACEWF
     I      ( FLUX_MULTIPLIER, IBEAM, N_SURFACE_WFS,
     I        LS_BOA_SOURCE )

      ENDIF

C  Downwelling Albedo weighting functions
C  --------------------------------------

      IF ( DO_DNWELLING ) THEN
        CALL DNUSER_SURFACEWF
     I       ( FLUX_MULTIPLIER, IBEAM, N_SURFACE_WFS )
      ENDIF

C  mean value output
C  -----------------

      IF ( DO_INCLUDE_MVOUTPUT.OR.DO_QUAD_OUTPUT ) THEN
        CALL VLIDORT_LS_INTEGRATED_OUTPUT
     I  ( DO_INCLUDE_MVOUTPUT, 
     I    FLUX_MULTIPLIER, IBEAM, N_SURFACE_WFS,
     I    LS_BOA_THTONLY_SOURCE )
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_BVP_COLUMN_SETUP
     I   (  DO_INCLUDE_DIRECTBEAM, DO_INCLUDE_SURFEMISS, 
     I      SURFACE_FACTOR, FOURIER_COMPONENT,
     I      N_SURFACE_WFS, IBEAM_INDEX )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of linearized solution variables (output stored here)

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'

C  inputs
C  ------

C  .. direct beam inputs

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  .. surface emission input

      LOGICAL          DO_INCLUDE_SURFEMISS

C  Surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Fourier component

      INTEGER          FOURIER_COMPONENT

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  Beam index

      INTEGER          IBEAM_INDEX

C  local variables
C  ---------------

      INTEGER          N, IB, Q, M
      INTEGER          K, KO1, K0, K1, K2
      INTEGER          I, J, O1, O2, OM, IR, IROW, C0, CM
      DOUBLE PRECISION L_BEAM, L_HOM_R, L_HOM_CR, H1R, H1I
      DOUBLE PRECISION FACTOR_BRDF, REFL_ATTN, AWF_DIRECT, AWF_EMISS
      DOUBLE PRECISION H_1,   H_2,   H_1_S,   H_2_S, H1, H2, DBR
      DOUBLE PRECISION H_1_CR,   H_2_CR,   H_1_CI,   H_2_CI
      DOUBLE PRECISION H_1_S_CR, H_2_S_CR, H_1_S_CI, H_2_S_CI

C  Help arrays

      DOUBLE PRECISION PV_W  ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION HV_P  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION HV_M  ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

C  Ground level boundary condition
C  -------------------------------

C  Initialise

      IB = IBEAM_INDEX
      M = FOURIER_COMPONENT
      FACTOR_BRDF = SURFACE_FACTOR

C  initialise. Vitally necessary
C    We found that when commented out, the whole thing didn't work.
C    R. Spurr and V. Natraj, 20 january 2006

      DO Q = 1, N_SURFACE_WFS
        DO I = 1, NTOTAL
          COL2_WFALB(I,Q) = ZERO
        ENDDO
      ENDDO

C  If this is the surface

      N  = NLAYERS
      C0 = N*NSTKS_NSTRMS_2 - NSTKS_NSTRMS

C  Lambertian case
C  ===============

C  Skip if not flagged

      IF ( .not. DO_LAMBERTIAN_SURFACE ) goto 998

C  Only 1 weighting function. Q = 1.

      DO Q = 1, N_SURFACE_WFS
        DO I = 1, NSTREAMS
          IR = NSTOKES*(I-1)
          DO O1 = 1, NSTOKES
            IROW = IR + O1
            CM   = C0 + IROW

C  Beam contribution for this Kernel

            L_BEAM = R2_BEAM(I,O1)

C  Real homogeneous solutions for this Kernel

            L_HOM_R = ZERO
            DO K = 1, K_REAL(N)
              L_HOM_R = L_HOM_R
     &          + LCON(K,N) * R2_HOMP(I,O1,K) * T_DELT_EIGEN(K,N) 
     &          + MCON(K,N) * R2_HOMM(I,O1,K)
            ENDDO

C  Complex homogeneous solutions for this Kernel

            L_HOM_CR = ZERO
            KO1 = K_REAL(N) + 1
            DO K = 1, K_COMPLEX(N)
              K0 = 2*K-2
              K1 = KO1 + K0
              K2 = K1 + 1
              H1R =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K1,N)
     &              - R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K2,N)
              H1I =   R2_HOMP(I,O1,K1) * T_DELT_EIGEN(K2,N)
     &              + R2_HOMP(I,O1,K2) * T_DELT_EIGEN(K1,N)
              L_HOM_CR = L_HOM_CR
     &                  + LCON(K1,N) * H1R  - LCON(K2,N) * H1I 
     &                  + MCON(K1,N) * R2_HOMM(I,O1,K1)
     &                  - MCON(K2,N) * R2_HOMM(I,O1,K2)
            ENDDO

C  Final contribution

            COL2_WFALB(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

C  End Streams and Stokes loops

          ENDDO
        ENDDO

C  Add direct beam variation of albedo

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              COL2_WFALB(CM,Q) = 
     &         COL2_WFALB(CM,Q) + DIRECT_BEAM(I,IB,O1)
            ENDDO
          ENDDO
        ENDIF

C  If surface emission, include emissivity variation
C    This code added for Version 2.4RT

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          O1   = 1
          DO I = 1, NSTREAMS
            IR   = NSTOKES*(I-1)
            IROW = IR + O1
            CM   = C0 + IROW
            COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) - SURFBB
          ENDDO
        ENDIF

C  Copy for the single layer case

        IF ( NLAYERS .EQ. 1 ) THEN
          DO N = 1, NTOTAL
            SCOL2_WFALB(N,Q) = COL2_WFALB(N,Q)
          ENDDO
        ENDIF

C  End parameter loop

      ENDDO

C  Return after finishing Lambertian case

      RETURN

C  BRDF boundary conditions
C  ========================

C  Continuation point

998   continue

C  Save some quantities
C  --------------------

C  ( This is a repetition of earlier code and could be stored)

C  start loops

      DO J = 1, NSTREAMS
        DO O1 = 1, NSTOKES

C  Beam

          PV_W(J,O1) = WLOWER(J,O1,N) * QUAD_STRMWTS(J)

C  real homogeneous solution contributions

          DO K = 1, K_REAL(N)
            H1 = SOLA_XPOS(J,O1,K,N)
            H2 = SOLB_XNEG(J,O1,K,N)
            HV_P(J,O1,K) = QUAD_STRMWTS(J)*H1
            HV_M(J,O1,K) = QUAD_STRMWTS(J)*H2
          ENDDO

C  Complex homogeneous solution contributions

          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            HV_P(J,O1,K1) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K1,N)
            HV_P(J,O1,K2) = QUAD_STRMWTS(J)* SOLA_XPOS(J,O1,K2,N)
            HV_M(J,O1,K1) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K1,N)
            HV_M(J,O1,K2) = QUAD_STRMWTS(J)* SOLB_XNEG(J,O1,K2,N)
          ENDDO

C  End loops

        ENDDO
      ENDDO

C  Diffuse scatter contributions
C  -----------------------------

C  Start weighting function loop

      DO Q = 1, N_SURFACE_WFS

C  start loops

        DO I = 1, NSTREAMS
         IR = NSTOKES*(I-1)
         DO O1 = 1, NSTOKES
          IROW = IR + O1
          CM   = C0 + IROW

C  Beam contribution for this Kernel

          REFL_B = ZERO
          DO J = 1, NSTREAMS
            DO O2 = 1, NSTOKES
              OM = MUELLER_INDEX(O1,O2)
              REFL_B = REFL_B + PV_W(J,O2) * LS_BRDF_F(Q,M,OM,J,I)
            ENDDO
          ENDDO
          L_BEAM = REFL_B * FACTOR_BRDF

C  Real homogeneous solutions contribution to this Kernel

          L_HOM_R = ZERO
          DO K = 1, K_REAL(N)
            H_1 = ZERO
            H_2 = ZERO
            DO J = 1, NSTREAMS
              H_1_S = ZERO
              H_2_S = ZERO
              DO O2 = 1, NSTOKES
                OM  = MUELLER_INDEX(O1,O2)
                H_1_S = H_1_S + HV_P(J,O2,K) * LS_BRDF_F(Q,M,OM,J,I)
                H_2_S = H_2_S + HV_M(J,O2,K) * LS_BRDF_F(Q,M,OM,J,I)
              ENDDO
              H_1 = H_1 + H_1_S
              H_2 = H_2 + H_2_S
            ENDDO
            H_1 = FACTOR_BRDF * H_1
            H_2 = FACTOR_BRDF * H_2
            L_HOM_R = L_HOM_R + LCON(K,N) * H_1 * T_DELT_EIGEN(K,N) 
     &                        + MCON(K,N) * H_2
          ENDDO

C  homogeneous complex solutions

          L_HOM_CR = ZERO
          KO1 = K_REAL(N) + 1
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1 + 1
            H_1_CR = ZERO
            H_2_CR = ZERO
            H_1_CI = ZERO
            H_2_CI = ZERO
            DO J = 1, NSTREAMS
              H_1_S_CR = ZERO
              H_2_S_CR = ZERO
              H_1_S_CI = ZERO
              H_2_S_CI = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                DBR = LS_BRDF_F(Q,M,OM,J,I)
                H_1_S_CR = H_1_S_CR + HV_P(J,O2,K1) * DBR
                H_2_S_CR = H_2_S_CR + HV_M(J,O2,K1) * DBR
                H_1_S_CI = H_1_S_CI + HV_P(J,O2,K2) * DBR
                H_2_S_CI = H_2_S_CI + HV_M(J,O2,K2) * DBR
              ENDDO
              H_1_CR = H_1_CR + H_1_S_CR
              H_2_CR = H_2_CR + H_2_S_CR
              H_1_CI = H_1_CI + H_1_S_CI
              H_2_CI = H_2_CI + H_2_S_CI
            ENDDO
            H_1_CR = FACTOR_BRDF  * H_1_CR
            H_1_CI = FACTOR_BRDF  * H_1_CI
            H_2_CR = FACTOR_BRDF  * H_2_CR
            H_2_CI = FACTOR_BRDF  * H_2_CI
            H1R =   H_1_CR * T_DELT_EIGEN(K1,N)
     &            - H_1_CI * T_DELT_EIGEN(K2,N)
            H1I =   H_1_CR * T_DELT_EIGEN(K2,N)
     &            + H_1_CI * T_DELT_EIGEN(K1,N)
            L_HOM_CR = L_HOM_CR
     &        + LCON(K1,N) *   H1R  - LCON(K2,N) *   H1I 
     &        + MCON(K1,N) * H_2_CR - MCON(K2,N) * H_2_CI
          ENDDO

C  Final contribution

          COL2_WFALB(CM,Q) = L_BEAM + L_HOM_R + L_HOM_CR

C  End loops

         ENDDO
        ENDDO

C  Direct beam reflection

        IF ( DO_INCLUDE_DIRECTBEAM ) THEN
          REFL_ATTN = ATMOS_ATTN(IB)
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              OM = MUELLER_INDEX(O1,1)
              AWF_DIRECT = REFL_ATTN * LS_BRDF_F_0(Q,M,OM,I,IB)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_DIRECT
            ENDDO
          ENDDO
        ENDIF

C  If surface emission, include emissivity variation, BRDF surface

        IF ( DO_INCLUDE_SURFEMISS ) THEN
          DO I = 1, NSTREAMS
            IR = NSTOKES*(I-1)
            DO O1 = 1, NSTOKES
              IROW = IR + O1
              CM   = C0 + IROW
              AWF_EMISS = SURFBB * LS_EMISSIVITY(Q,O1,I)
              COL2_WFALB(CM,Q) = COL2_WFALB(CM,Q) + AWF_EMISS
            ENDDO
          ENDDO
        ENDIF

C  Copy for the single layer case

        IF ( NLAYERS .EQ. 1 ) THEN
           DO N = 1, NTOTAL
            SCOL2_WFALB(N,Q) = COL2_WFALB(N,Q)
          ENDDO
        ENDIF

C  End parameter loop

      ENDDO

C  debug

c      if ( do_debug_write ) then
c        DO N = 1, NTOTAL
c          write(85,'(2i4,1p4e17.9)')IBEAM_INDEX, N, COL2_WFALB(N,1)
c        ENDDO
c        pause
c      ENDIF
           
C  Finish

      RETURN
      END

c

      SUBROUTINE GET_LS_BOA_SOURCE
     I        ( DO_INCLUDE_DIRECTBEAM,
     I          DO_INCLUDE_SURFEMISS,
     I          DO_INCLUDE_MVOUTPUT,
     I          SURFACE_FACTOR, N_SURFACE_WFS,
     I          FOURIER_COMPONENT, IBEAM, 
     O          LS_BOA_SOURCE,
     O          LS_BOA_THTONLY_SOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables

      INCLUDE '../includes/VLIDORT_BRDF.VARS'
      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include file of Linearized solution variables

      INCLUDE '../includes/VLIDORT_L_BRDF.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_THERMALSUP.VARS'

C  Module arguments
C  ----------------

C  directbeam inclusion flag

      LOGICAL          DO_INCLUDE_DIRECTBEAM

C  surface emissivity inclusion flag

      LOGICAL          DO_INCLUDE_SURFEMISS

C  MV output inclusion flag

      LOGICAL          DO_INCLUDE_MVOUTPUT

C  Fourier surface factor

      DOUBLE PRECISION SURFACE_FACTOR

C  Number of weighting functions

      INTEGER          N_SURFACE_WFS

C  Fourier number

      INTEGER          FOURIER_COMPONENT

C  Beam index

      INTEGER          IBEAM

C  output

      DOUBLE PRECISION LS_BOA_SOURCE
     &       ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_BOA_THTONLY_SOURCE
     &       ( MAX_SURFACEWFS, MAXSTREAMS, MAXSTOKES)

C  Local variables
C  ---------------

      LOGICAL          DO_QTHTONLY
      INTEGER          UM, I, J, N, IB, O1, O2, M, OM, Q
      INTEGER          K, KO1, K0, K1, K2    
      DOUBLE PRECISION INTEGRAND(MAX_SURFACEWFS,MAXSTREAMS,MAXSTOKES)
      DOUBLE PRECISION SUM_R, SUM_CR, REFLEC, S_REFLEC, REFL_ATTN
      DOUBLE PRECISION H1, H2, NXR, PXR, NXR1, NXR2, PXR1

C  Initialise
C  ----------

C  Special flag

      DO_QTHTONLY = ( DO_THERMAL_TRANSONLY ) .AND.
     &      ( DO_QUAD_OUTPUT .OR. DO_INCLUDE_MVOUTPUT )

C  Shorthand

      N   = NLAYERS
      KO1 = K_REAL(N) + 1
      IB  = IBEAM
      M   = FOURIER_COMPONENT

C  initialise Derivative of BOA source function

      DO Q = 1, N_SURFACE_WFS
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO O1 = 1, NSTOKES
            LS_BOA_SOURCE(Q,UM,O1) = ZERO
          ENDDO
        ENDDO
      ENDDO

C  Thermal tranmsittance only, special term

      IF ( DO_QTHTONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO I = 1, NSTREAMS
            DO O1 = 1, NSTOKES
              LS_BOA_THTONLY_SOURCE(Q,I,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Skip diffuse-field variation for thermal transmittance-only

      IF ( DO_THERMAL_TRANSONLY .OR. .NOT. DO_USER_STREAMS ) GO TO 599

C  Contribution due to derivatives of BV constants
C  -----------------------------------------------

C  First compute derivative of downward intensity Integrand at stream angles 
C        .. reflectance integrand  = a(j).x(j).I_DOWN(-j)

C  start loops

      DO Q = 1, N_SURFACE_WFS
       DO I = 1, NSTREAMS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

         SUM_R = ZERO
         DO K = 1, K_REAL(N)
          NXR = NCON_ALB(Q,K,N) * SOLA_XPOS(I,O1,K,N)
          PXR = PCON_ALB(Q,K,N) * SOLB_XNEG(I,O1,K,N)
          SUM_R = SUM_R + NXR*T_DELT_EIGEN(K,N) + PXR
         ENDDO

C  Complex solutions

         SUM_CR = ZERO
         DO K = 1, K_COMPLEX(N)
          K0 = 2 * K - 2
          K1 = KO1 + K0
          K2 = K1  + 1
          NXR1 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K1,N)
     &           - NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K2,N)
          NXR2 =   NCON_ALB(Q,K1,N) * SOLA_XPOS(I,O1,K2,N)
     &           + NCON_ALB(Q,K2,N) * SOLA_XPOS(I,O1,K1,N)
          PXR1 =   PCON_ALB(Q,K1,N) * SOLB_XNEG(I,O1,K1,N)
     &           - PCON_ALB(Q,K2,N) * SOLB_XNEG(I,O1,K2,N)
          H1 =  NXR1 * T_DELT_EIGEN(K1,N)
     &         -NXR2 * T_DELT_EIGEN(K2,N)
          H2 =  PXR1
          SUM_CR = SUM_CR + H1 + H2
         ENDDO

C  Final result

         INTEGRAND(Q,I,O1) = QUAD_STRMWTS(I) * ( SUM_R + SUM_CR )

C  end loops

        ENDDO
       ENDDO
      ENDDO

C  integrated reflectance term
C  ---------------------------

C  integrate Lambertian case, same for all user-streams

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
       IF ( FOURIER_COMPONENT.EQ.0 ) THEN
        O1 = 1
        Q = 1
        REFLEC = ZERO
        DO J = 1, NSTREAMS
           REFLEC = REFLEC + INTEGRAND(Q,J,O1)
        ENDDO
        REFLEC = SURFACE_FACTOR * REFLEC * LAMBERTIAN_ALBEDO
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          LS_BOA_SOURCE(Q,UM,O1) = REFLEC
        ENDDO
       ENDIF
      ENDIF

C  BRDF case

      IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + INTEGRAND(Q,J,O2) * 
     &                           USER_BRDF_F(M,OM,UM,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              LS_BOA_SOURCE(Q,UM,O1) = REFLEC * SURFACE_FACTOR
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Contributions due to direct variation of kernel parameter
C  ---------------------------------------------------------

C  Lambertian (this is the albedo)

      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER_COMPONENT .EQ. 0 ) THEN
          O1 = 1
          Q = 1
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            REFLEC = ZERO
            DO J = 1, NSTREAMS
              REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
            ENDDO
            REFLEC = REFLEC * LAMBERTIAN_ALBEDO * SURFACE_FACTOR
            LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
          ENDDO
        ENDIF
      ENDIF


C   Non-Lambertian (need derivative of BRDF Fourier term)

      IF ( .not.DO_LAMBERTIAN_SURFACE ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO J = 1, NSTREAMS
                S_REFLEC = ZERO
                DO O2 = 1, NSTOKES
                  OM = MUELLER_INDEX(O1,O2)
                  S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) *
     &                              LS_USER_BRDF_F(Q,M,OM,UM,J)
                ENDDO
                REFLEC = REFLEC + S_REFLEC
              ENDDO
              REFLEC = REFLEC * SURFACE_FACTOR
              LS_BOA_SOURCE(Q,UM,O1)=LS_BOA_SOURCE(Q,UM,O1) + REFLEC
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Continuation point for avoiding diffuse field computation

 599  continue

C  Thermal tranmsittance and Integrated output (quadrature terms)
C     This is the reflectance term

      IF ( DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q = 1
          REFLEC = ZERO
          DO J = 1, NSTREAMS
            REFLEC = REFLEC + STOKES_DOWNSURF(J,O1)
          ENDDO
          REFLEC = SURFACE_FACTOR * REFLEC
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = 
     &          LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO I = 1, NSTREAMS
              DO O1 = 1, NSTOKES
                REFLEC = ZERO
                DO J = 1, NSTREAMS
                  S_REFLEC = ZERO
                  DO O2 = 1, NSTOKES
                    OM = MUELLER_INDEX(O1,O2)
                    S_REFLEC = S_REFLEC + STOKES_DOWNSURF(J,O2) *
     &                              LS_BRDF_F(Q,M,OM,I,J)
                  ENDDO
                  REFLEC = REFLEC + S_REFLEC
                ENDDO
                REFLEC = REFLEC * SURFACE_FACTOR
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = 
     &            LS_BOA_THTONLY_SOURCE(Q,I,O1)  + REFLEC
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Add linearization term for variation of direct beam reflectance 

      IF ( DO_INCLUDE_DIRECTBEAM.and..not.DO_DBCORRECTION ) THEN
        REFL_ATTN = ATMOS_ATTN(IB)
        IF ( DO_LAMBERTIAN_SURFACE.AND.M.EQ.0) THEN
          O1 = 1
          Q = 1
          REFLEC = SURFACE_FACTOR * REFL_ATTN * LAMBERTIAN_ALBEDO
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
          ENDDO
        ELSE IF ( .not. DO_LAMBERTIAN_SURFACE ) THEN
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              REFLEC = ZERO
              DO O2 = 1, NSTOKES
                OM = MUELLER_INDEX(O1,O2)
                REFLEC = REFLEC + 
     *             LS_USER_BRDF_F_0(Q,M,OM,UM,IB) * FLUXVEC(O2)
              ENDDO
c              REFLEC = SURFACE_FACTOR * REFLEC
              REFLEC = SURFACE_FACTOR * REFL_ATTN * REFLEC ! V. Natraj, 11/15/2010
              LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) + REFLEC
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Add emissivity variation at user defined angles.
C    Apparenly only present for Fourier zero
C  (expression for emissivity variation follows from Kirchhoff's law)

      IF ( DO_INCLUDE_SURFEMISS .and. M.eq.0 ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q = 1
          DO UM = LOCAL_UM_START, N_USER_STREAMS
            LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) - SURFBB
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O 1= 1, NSTOKES
              DO UM = LOCAL_UM_START, N_USER_STREAMS
                LS_BOA_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1) 
     &                      - SURFBB * LS_USER_EMISSIVITY(Q,O1,UM)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Thermal Transmittance only (For Integrated product)
C    Surface emissivity term

      IF ( DO_INCLUDE_SURFEMISS .and. DO_QTHTONLY ) THEN
        IF ( DO_LAMBERTIAN_SURFACE ) THEN
          O1 = 1
          Q = 1
          DO I = 1, NSTREAMS
            LS_BOA_THTONLY_SOURCE(Q,I,O1) = 
     &       LS_BOA_THTONLY_SOURCE(Q,I,O1) -  SURFBB
          ENDDO
        ELSE
          DO Q = 1, N_SURFACE_WFS
            DO O 1= 1, NSTOKES
              DO I = 1, NSTREAMS
                LS_BOA_THTONLY_SOURCE(Q,I,O1) = 
     &                LS_BOA_THTONLY_SOURCE(Q,I,O1) 
     &                   - SURFBB * LS_EMISSIVITY(Q,O1,I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF

C  Finish

      RETURN
      END

C

      SUBROUTINE UPUSER_SURFACEWF
     I   ( FLUX_MULTIPLIER, IBEAM, N_SURFACE_WFS,
     I     LS_BOA_SOURCE )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include files of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Flux multiplier

      DOUBLE PRECISION FLUX_MULTIPLIER

C  Beam index

      INTEGER          IBEAM

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  derivatives of reflected surface upwelling intensity

      DOUBLE PRECISION LS_BOA_SOURCE
     &    ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, Q, UT, IB

      DOUBLE PRECISION LS_CUMUL_SOURCE
     &        ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_LAYER_SOURCE
     &        ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_FINAL_SOURCE

C  index

      IB = IBEAM

C  Zero all Fourier components - New rule, better for safety
C    Only did this for components close to zenith (formerly)

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, N_USER_STREAMS
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

      IF ( DO_USER_STREAMS ) THEN

C  Set the cumulative source term equal to the BOA sum

        DO UM = LOCAL_UM_START, N_USER_STREAMS
         DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
           LS_CUMUL_SOURCE(Q,UM,O1) = LS_BOA_SOURCE(Q,UM,O1)
          ENDDO
         ENDDO
        ENDDO

      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NUT = 0
      NSTART = NLAYERS
      NUT_PREV = NSTART + 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = N_USER_LEVELS, 1, -1

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_UP(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term (MSST) output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL + 1
          DO N = NSTART, NUT, -1
            CALL LS_WHOLELAYER_STERM_UP
     I       ( N_SURFACE_WFS, IB, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1)
     &               + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_UP
     I        ( N_SURFACE_WFS, IB, UT, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1)
     &          + T_UTUP_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                  SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) =
     &                   FLUX_MULTIPLIER * LS_FINAL_SOURCE
                ENDDO
              ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,UM,IB,O1,UPIDX) =
     &                     FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT .NE. NUT_PREV ) NSTART = NUT - 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE DNUSER_SURFACEWF
     I    ( FLUX_MULTIPLIER,
     I      IBEAM, N_SURFACE_WFS )

C  include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of setup, solution and reflectance variables (input)

      INCLUDE '../includes/VLIDORT_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_REFLECTANCE.VARS'

C  include files of linearized setup and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SETUPS.VARS'
      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  include file of result variables (module output stored here)

      INCLUDE '../includes/VLIDORT_L_RESULTS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Beam index

      INTEGER          IBEAM

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  Flux multiplier = F/4.pi

      DOUBLE PRECISION FLUX_MULTIPLIER

C  local variables
C  ---------------

      INTEGER          N, NUT, NSTART, NUT_PREV, NLEVEL, O1
      INTEGER          UTA, UM, Q, UT, IB

      DOUBLE PRECISION LS_CUMUL_SOURCE
     &        ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_LAYER_SOURCE
     &        ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_TOA_SOURCE
     &        ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )
      DOUBLE PRECISION LS_FINAL_SOURCE

C  Initialise

      IB = IBEAM

C  Zero all Fourier component output

      IF ( DO_USER_STREAMS ) THEN
        DO UTA = 1, N_USER_LEVELS
          DO UM = 1, LOCAL_UM_START 
            DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) = ZERO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Initialize post-processing recursion
C  ====================================

C  Get the linearized TOA source terms

      IF ( DO_USER_STREAMS ) THEN
        DO UM = LOCAL_UM_START, N_USER_STREAMS
          DO Q = 1, N_SURFACE_WFS
            DO O1 = 1, NSTOKES
              LS_TOA_SOURCE(Q,UM,O1)   = ZERO
              LS_CUMUL_SOURCE(Q,UM,O1) = LS_TOA_SOURCE(Q,UM,O1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C  Recursion Loop for linearized Post-processing
C  =============================================

C  initialise cumulative source term loop

      NUT = 0
      NSTART = 1
      NUT_PREV = NSTART - 1

C  loop over all output optical depths
C  -----------------------------------

      DO UTA = 1, N_USER_LEVELS

C  Layer index for given optical depth

        NLEVEL = UTAU_LEVEL_MASK_DN(UTA)

C  Cumulative source terms to layer NUT (user-defined stream angles only)
C    1. Get layer source terms
C    2. Find cumulative source term
C    3. Set multiple scatter source term output if flagged

        IF ( DO_USER_STREAMS ) THEN
          NUT = NLEVEL
          DO N = NSTART, NUT
            CALL LS_WHOLELAYER_STERM_DN
     &       ( N_SURFACE_WFS, IB, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                LS_CUMUL_SOURCE(Q,UM,O1) = LS_LAYER_SOURCE(Q,UM,O1)
     &               + T_DELT_USERM(N,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
              ENDDO
             ENDDO
            ENDDO
          ENDDO
        ENDIF

C  Offgrid output
C  --------------

        IF ( PARTLAYERS_OUTFLAG(UTA) ) THEN

          UT = PARTLAYERS_OUTINDEX(UTA)
          N  = PARTLAYERS_LAYERIDX(UT)

C  User-defined stream output, add additional partial layer source term

          IF ( DO_USER_STREAMS ) THEN
            CALL LS_PARTLAYER_STERM_DN
     I        ( N_SURFACE_WFS, IB, UT, N, LS_LAYER_SOURCE )
            DO UM = LOCAL_UM_START, N_USER_STREAMS
             DO Q = 1, N_SURFACE_WFS
              DO O1 = 1, NSTOKES
                LS_FINAL_SOURCE = LS_LAYER_SOURCE(Q,UM,O1)
     &           + T_UTDN_USERM(UT,UM) * LS_CUMUL_SOURCE(Q,UM,O1)
                SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) =
     &                FLUX_MULTIPLIER * LS_FINAL_SOURCE
              ENDDO
             ENDDO
            ENDDO
          ENDIF

C  Ongrid output
C  -------------

        ELSE

C  User-defined stream output, just set to the cumulative source term

          IF ( DO_USER_STREAMS ) THEN
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              DO Q = 1, N_SURFACE_WFS
                DO O1 = 1, NSTOKES
                  SURFACEWF_F(Q,UTA,UM,IB,O1,DNIDX) =
     &                    FLUX_MULTIPLIER * LS_CUMUL_SOURCE(Q,UM,O1)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

C  Check for updating the recursion

        IF ( DO_USER_STREAMS ) THEN
          IF ( NUT. NE. NUT_PREV ) NSTART = NUT + 1
          NUT_PREV = NUT
        ENDIF

C  end loop over optical depth

      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_WHOLELAYER_STERM_UP
     I       ( N_SURFACE_WFS, IBEAM, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE
     &     ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  Homogeneous solutions
C  =====================

C  Loops over user angles, weighting functions and Stokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N)*UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N)*UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_2(K,UM,N)
            H2 =  PUXR * HMULT_1(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N)
     &              - NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N)
     &              + NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N)
     &              - PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N)
     &              + PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_2(K1,UM,N)
     &           - NUXR2 * HMULT_2(K2,UM,N)
            H2 =   PUXR1 * HMULT_1(K1,UM,N)
     &           - PUXR2 * HMULT_1(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_WHOLELAYER_STERM_DN
     I       ( N_SURFACE_WFS, IBEAM, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  Indices

      INTEGER          IBEAM, GIVEN_LAYER

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE
     &     ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      IB  = IBEAM

C  Homogeneous solutions
C  =====================

C  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N)*UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N)*UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * HMULT_1(K,UM,N)
            H2 =  PUXR * HMULT_2(K,UM,N)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N)
     &              - NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N)
     &              + NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N)
     &              - PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N)
     &              + PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * HMULT_1(K1,UM,N)
     &           - NUXR2 * HMULT_1(K2,UM,N)
            H2 =   PUXR1 * HMULT_2(K1,UM,N)
     &           - PUXR2 * HMULT_2(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_PARTLAYER_STERM_UP
     I       ( N_SURFACE_WFS, IBEAM, OFFGRID_INDEX, GIVEN_LAYER,
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  layer and beam index

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE
     &     ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N) * UHOM_UPDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N) * UHOM_UPUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_UD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_UU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K1,N)
     &              - NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_UPDN(UM,O1,K2,N)
     &              + NCON_ALB(Q,K2,N) * UHOM_UPDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K1,N)
     &              - PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_UPUP(UM,O1,K2,N)
     &              + PCON_ALB(Q,K2,N) * UHOM_UPUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_UD(K1,UM,N)
     &           - NUXR2 * UT_HMULT_UD(K2,UM,N)
            H2 =   PUXR1 * UT_HMULT_UU(K1,UM,N)
     &           - PUXR2 * UT_HMULT_UU(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

C  Finish

      RETURN
      END

C

      SUBROUTINE LS_PARTLAYER_STERM_DN
     I       ( N_SURFACE_WFS, IBEAM, OFFGRID_INDEX, GIVEN_LAYER, 
     O         LS_LAYERSOURCE )

C  Include files
C  -------------

C  include file of dimensions and numbers

      INCLUDE '../includes/VLIDORT.PARS'

C  include files of input and bookkeeping variables

      INCLUDE '../includes/VLIDORT_INPUTS.VARS'
      INCLUDE '../includes/VLIDORT_BOOKKEEP.VARS'

C  include files of ssolution and multiplier variables (input)

      INCLUDE '../includes/VLIDORT_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_MULTIPLIERS.VARS'

C  include files of linearized multiplier and solution variables (input)

      INCLUDE '../includes/VLIDORT_L_SOLUTION.VARS'
      INCLUDE '../includes/VLIDORT_L_MULTIPLIERS.VARS'

C  Subroutine input arguments
C  --------------------------

C  Number of surface weighting functions

      INTEGER          N_SURFACE_WFS

C  layer and beam indices

      INTEGER          GIVEN_LAYER, IBEAM

C  offgrid optical depth index

      INTEGER          OFFGRID_INDEX

C  Subroutine output arguments
C  ---------------------------

      DOUBLE PRECISION LS_LAYERSOURCE
     &     ( MAX_SURFACEWFS, MAX_USER_STREAMS, MAXSTOKES )

C  local variables
C  ---------------

      INTEGER          N, UM, O1, IB, UT, K, KO1, K0, K1, K2, Q
      DOUBLE PRECISION SHOM_R, SHOM_CR, H1, H2
      DOUBLE PRECISION NUXR, PUXR, NUXR1, NUXR2, PUXR1, PUXR2

C  Thermal transmittance only

      IF ( DO_THERMAL_TRANSONLY ) THEN
        DO Q = 1, N_SURFACE_WFS
          DO O1 = 1, NSTOKES
            DO UM = LOCAL_UM_START, N_USER_STREAMS
              LS_LAYERSOURCE(Q,UM,O1) = ZERO
            ENDDO
          ENDDO
        ENDDO 
        RETURN
      ENDIF

C  local indices

      N   = GIVEN_LAYER
      KO1 = K_REAL(N) + 1
      UT  = OFFGRID_INDEX
      IB  = IBEAM

C  Partial layer source function ( Homogeneous/constants variation )
C  =================================================================

C  Loop over user angles, weighting functions and STokes

      DO UM = LOCAL_UM_START, N_USER_STREAMS
       DO Q = 1, N_SURFACE_WFS
        DO O1 = 1, NSTOKES

C  Real homogeneous solutions

          SHOM_R = ZERO
          DO K = 1, K_REAL(N)
            NUXR = NCON_ALB(Q,K,N) * UHOM_DNDN(UM,O1,K,N)
            PUXR = PCON_ALB(Q,K,N) * UHOM_DNUP(UM,O1,K,N)
            H1 =  NUXR * UT_HMULT_DD(K,UM,UT)
            H2 =  PUXR * UT_HMULT_DU(K,UM,UT)
            SHOM_R = SHOM_R + H1 + H2
          ENDDO

C  Complex homogeneous solutions

          SHOM_CR = ZERO
          DO K = 1, K_COMPLEX(N)
            K0 = 2 * K - 2
            K1 = KO1 + K0
            K2 = K1  + 1
            NUXR1 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K1,N)
     &              - NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K2,N)
            NUXR2 =   NCON_ALB(Q,K1,N) * UHOM_DNDN(UM,O1,K2,N)
     &              + NCON_ALB(Q,K2,N) * UHOM_DNDN(UM,O1,K1,N)
            PUXR1 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K1,N)
     &              - PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K2,N)
            PUXR2 =   PCON_ALB(Q,K1,N) * UHOM_DNUP(UM,O1,K2,N)
     &              + PCON_ALB(Q,K2,N) * UHOM_DNUP(UM,O1,K1,N)
            H1 =   NUXR1 * UT_HMULT_DD(K1,UM,N)
     &           - NUXR2 * UT_HMULT_DD(K2,UM,N)
            H2 =   PUXR1 * UT_HMULT_DU(K1,UM,N)
     &           - PUXR2 * UT_HMULT_DU(K2,UM,N)
            SHOM_CR = SHOM_CR + H1 + H2
          ENDDO

C  homogeneous contribution

          LS_LAYERSOURCE(Q,UM,O1) = SHOM_R + SHOM_CR

C  End loops over Q, O1 and UM

        ENDDO
       ENDDO
      ENDDO

C  Finish

      RETURN
      END

