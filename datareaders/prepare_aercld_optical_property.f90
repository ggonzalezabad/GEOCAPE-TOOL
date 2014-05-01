SUBROUTINE prepare_aercld_optical_property(database_dir, mtype, ntype, thetypes, nlayers, &
     profile, flags, lambdas, nlambda, lambda0, nmoms, ngksec, &
     opdeps, ssas,  qexts, phfcn, message, fail) 

  IMPLICIT NONE

  ! ========================
  ! Input/output parameters
  ! ========================
  INTEGER, INTENT(IN) :: nlambda, nmoms, ngksec, mtype, ntype, nlayers
  CHARACTER (LEN=2), DIMENSION(mtype), INTENT(IN)      :: thetypes
  REAL(KIND=8), DIMENSION (nlambda), INTENT(IN)        :: lambdas
  LOGICAL, DIMENSION(nlayers), INTENT(IN)              :: flags
  REAL(KIND=8), DIMENSION (mtype, nlayers), INTENT(IN) :: profile
  REAL(KIND=8), INTENT(IN)                             :: lambda0
  CHARACTER(LEN=256),INTENT(IN)                        :: database_dir
  LOGICAL, INTENT(OUT)                                 :: fail
  CHARACTER*(*), INTENT(INOUT)                         :: message
  REAL(KIND=8), DIMENSION (nlambda, nlayers), INTENT(OUT) :: opdeps, ssas, qexts
  REAL(KIND=8), DIMENSION (nlambda, nlayers, 0:nmoms, ngksec), INTENT(OUT) :: phfcn
  
  INTEGER, PARAMETER :: naer=8, mwl=401, maxmoms=1000, maxgksec=6
  ! Black carbon, organic carbon, dust, sulfate, sea salt fine mode, 
  ! sea salt coarse mode, water clouds, ice clouds
  CHARACTER(LEN=2), dimension(naer), parameter :: aertypes = &
       (/'BC', 'OC', 'DU', 'SU', 'SF', 'SC', 'WC', 'IC'/)
  INTEGER, DIMENSION(naer), SAVE           :: nwls
  REAL(KIND=8), DIMENSION(naer, mwl), SAVE :: wls0, qext0, ssa0
  REAL(KIND=8), DIMENSION(naer, mwl, 0:maxmoms, maxgksec), SAVE :: phfcn0
  LOGICAL, SAVE :: first = .TRUE. 

  ! Local variable
  INTEGER       :: i, it, j, n, k, idx0, nw0, aeridx
  INTEGER, DIMENSION (nlambda) :: idxs
  LOGICAL                      :: bound
  REAL(kind=8)                 :: dum, aqext, tsca
  REAL(KIND=8), DIMENSION(1)   :: tmp_lambda0, fs0
  REAL(kind=8), DIMENSION (mtype)   :: scaods
  REAL(KIND=8), DIMENSION (nlambda) :: fs
  REAL(KIND=8), DIMENSION (mtype, nlambda) :: ssas1, qexts1
  REAL(KIND=8), DIMENSION (mtype, nlambda, 0:nmoms, ngksec) :: phfcn1

  fail    = .FALSE.
  message = ' '

  IF (first) THEN

     ! Read aerosol data
     IF (nmoms > maxmoms) THEN    
        message = 'Not enough phase moments, reduce nmoms!!!'
        fail = .true.; return
     ENDIF

     DO n = 1, naer
        OPEN(1, FILE = TRIM(ADJUSTL(database_dir)) // '/AerCldProp/' // aertypes(n) // '_vprop1000.dat', STATUS='old')
        READ(1, *)      ! Header 
        READ(1, *)      ! Aerosol label
        READ(1, *) nwls(n)
        IF (nwls(n) > mwl) THEN
           message = ' Aerosol input wavelengths exceeds MWL. Increase MWL!!!'
           fail = .true.; return
        ENDIF
        
        DO i = 1, nwls(n)
           READ(1, *) wls0(n, i), qext0(n, i), dum, ssa0(n, i), &
                (phfcn0(n, i, j, 1), j = 0, nmoms)
           
           DO k = 2, maxgksec
              READ(1, *) (phfcn0(n, i, j, k), j = 0, nmoms)
           ENDDO
        ENDDO
        wls0(n, 1:nwls(n)) = wls0(n, 1:nwls(n)) * 1000.d0 ! convert to nm
        
        CLOSE(1) 
     ENDDO

     first = .FALSE.
  ENDIF

  ! Find which aerosol type
  DO it = 1, ntype 
     DO i = 1, naer
        IF (thetypes(it) == aertypes(i)) THEN
           aeridx = i; EXIT
        ENDIF
     ENDDO
     
     IF (i == naer + 1) THEN
        message = 'No aerosol properties found for aerosol/cloud ' // thetypes(it)
        fail = .true.; return
     ENDIF
     
     bound = .TRUE.
     nw0 = nwls(aeridx)
     
     call interpol(lambdas, nlambda, wls0(aeridx, 1:nw0), nw0, bound, fs)
     idxs = INT(fs); fs = fs - idxs
     WHERE (idxs > nw0) 
        idxs = nw0 - 1; fs = 1
     ENDWHERE
     WHERE (idxs < 1) 
        idxs = 1; fs = 0
     ENDWHERE
     ssas1 (it, 1:nlambda)  = ssa0(aeridx, idxs) * (1.d0 - fs) + ssa0(aeridx, idxs + 1) * fs
     qexts1(it, 1:nlambda)  = qext0(aeridx, idxs) * (1.d0 - fs) + qext0(aeridx, idxs + 1) * fs
     
     DO i = 1, ngksec
        DO j = 0, nmoms
           phfcn1(it, 1:nlambda, j, i) = phfcn0(aeridx, idxs, j, i) * (1.d0 - fs) + &
                phfcn0(aeridx, idxs+1, j, i) * fs
        ENDDO
     ENDDO
     
     tmp_lambda0(1) = lambda0
     call interpol(tmp_lambda0, 1, wls0(aeridx, 1:nw0), nw0, bound, fs0)
     idx0 = INT(fs0(1));  fs0(1) = fs0(1) - idx0
     IF (idx0 > nw0) THEN
        idx0 = nw0 - 1; fs0(1) = 1.d0
     ENDIF
     IF (idx0 < 1) THEN
        idx0 = 1; fs0(1) = 0.d0
     ENDIF
     aqext = qext0(aeridx, idx0) * (1.d0 - fs0(1)) + qext0(aeridx, idx0 + 1) * fs0(1)
     qexts1(it, 1:nlambda) = qexts1(it, 1:nlambda) / aqext
     
  ENDDO

  opdeps = 0.0d0
  ssas   = 0.0d0
  qexts  = 0.0d0
  phfcn  = 0.0d0

  DO i = 1, nlayers
     IF (flags(i)) THEN
        DO j = 1, nlambda
           opdeps(j, i) = SUM(profile(1:ntype, i) * qexts1(1:ntype, j))
           qexts (j, i) = opdeps(j, i) / SUM(profile(1:ntype, i))
           scaods(1:ntype) = profile(1:ntype, i) * qexts1(1:ntype, j) * ssas1(1:ntype, j)
           tsca = SUM(scaods(1:ntype))
           ssas  (j, i)   = tsca / opdeps(j, i)

           DO n = 1, ngksec
              DO k = 0, nmoms
                 phfcn(j, i, k, n) = SUM(phfcn1(1:ntype, j, k, n) * scaods(1:ntype)) / tsca
              ENDDO
           ENDDO
           
        ENDDO
     ENDIF
  ENDDO

  RETURN
  
END SUBROUTINE prepare_aercld_optical_property
