! =========================================================================
!
! Convolves input spectrum with an asymmetric Gaussian slit function of
! specified HW1E (half-width at 1/e intensity)
!
! The :symetric Gaussian g(x) is defined as
!                   _                _
!                  |         x^2      |
!      g(x) =  EXP | - -------------- |
!                  |_     hw1e ^2    _|
!  FWHM = 2.0 * sqrt(ln(2.0)) * hw1e = 1.66551 * hw1e
! =========================================================================

! Assume wavelength grid is evenly spaced and output/input has the same spectral grid
SUBROUTINE gauss (wvlarr, specarr, specmod, npoints, fwhm)
  
  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                            INTENT (IN)    :: npoints
  REAL (KIND=8),                      INTENT (IN)    :: fwhm
  REAL (KIND=8), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr
  REAL (KIND=8), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                            :: nhi, nlo, i, j, num_slit
  REAL (KIND=8)                      :: emult, delwvl, slitsum, slit0, hw1e
  REAL (KIND=8), DIMENSION (npoints) :: slit

  hw1e = fwhm / 1.66551

  ! ----------------------------------------------------------------
  ! Initialization of output variable (default for "no convolution")
  ! ----------------------------------------------------------------
  specmod(1:npoints) = specarr(1:npoints)
  
  ! --------------------------------------
  ! No convolution if halfwidth @ 1/e is 0
  ! --------------------------------------
  IF ( hw1e == 0.0 ) RETURN

  emult  = -1.0 / ( hw1e * hw1e )
  delwvl = wvlarr(2) - wvlarr(1)
 
  !  Calculate slit function values out to 0.001 times x0 value,
  !     normalize so that sum = 1.

  slitsum = 1.0 ;  slit0 = 1.0
  i = 1  ;  num_slit = 0
  DO WHILE ( num_slit <= npoints )
     slit (i) = EXP (emult * (delwvl * i)**2)
     slitsum = slitsum + 2.0 * slit (i)
     IF (slit (i) <= 0.001 ) EXIT 
     i = i + 1
  ENDDO
  num_slit = i
  
  slit0 = slit0 / slitsum
  slit(1:num_slit) = slit(1:num_slit) / slitsum

  ! Convolve spectrum.  reflect at endpoints.
  ! Doesn't look right
  specmod(1:npoints) = slit0 * specarr(1:npoints)
  DO i = 1, npoints
     DO j = 1, num_slit
        nlo = i - j 
        IF (nlo < 1) nlo = -nlo + 2 
        nhi = i + j
        IF ( nhi > npoints ) nhi = npoints - MOD(nhi, npoints)
        specmod(i) = specmod(i) + slit(j) * ( specarr(nlo) + specarr(nhi) )
     END DO
  END DO

  RETURN
END SUBROUTINE gauss

SUBROUTINE gaussio (wvlarr, specarr, i0, specmod, npoints, fwhm, scalex)

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                            INTENT (IN)    :: npoints
  REAL (KIND=8),                      INTENT (IN)    :: fwhm, scalex
  REAL (KIND=8), DIMENSION (npoints), INTENT (IN)    :: wvlarr, specarr, i0
  REAL (KIND=8), DIMENSION (npoints), INTENT (OUT)   :: specmod

  ! ======================
  ! Local variables
  ! ======================
  REAL (KIND=8), DIMENSION (npoints) :: i0mod, abspec, abspecmod

  specmod = specarr
  if (fwhm == 0.0) RETURN

  CALL gauss(wvlarr, i0, i0mod, npoints, fwhm)
  abspec = i0 * exp(-specarr * scalex)
  CALL gauss(wvlarr, abspec, abspecmod, npoints, fwhm)

  specmod = -LOG(abspecmod / i0mod) / scalex

  RETURN

END SUBROUTINE gaussio

! convolve high-resolution spectra to low resolution spectra
! fwave: high-resolution wavelength grid
! fspec: high-resolution reference spectra
! nf:    number of wavelengths at high-resolution
! nspec: number of spectra
! fwhm:  slit function in terms of FWHM (nm)
! cwave: low-resolution wavelength grid (could be the same as high-resolution grid)
! cspec: spectra at fwhm
! nc:    number of wavelengths for the low-resolution grid
!SUBROUTINE gauss_f2c (fwave, fspec, nf, nspec, fwhm, cwave, cspec, nc)
!
!  IMPLICIT NONE
!
!  ! =======================
!  ! Input/Output variables
!  ! =======================
!  INTEGER,                               INTENT (IN) :: nc, nf, nspec
!  REAL (KIND=8),                         INTENT (IN) :: fwhm
!  REAL (KIND=8), DIMENSION (nf), INTENT (IN)         :: fwave
!  REAL (KIND=8), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
!  REAL (KIND=8), DIMENSION (nc), INTENT (IN)         :: cwave
!  REAL (KIND=8), DIMENSION (nc, nspec), INTENT (OUT) :: cspec
!
!  ! ===============
!  ! Local variables
!  ! ===============
!  INTEGER                       :: i, j, midx, sidx, eidx, nhalf
!  REAL (KIND=8)                 :: hw1esq, dfw, ssum, hw1e
!  REAL (KIND=8), DIMENSION (nf) :: slit
!
!  if (fwhm == 0.0) then
!     return
!  endif
!
!  dfw  = fwave(2) - fwave(1)
!  hw1e = fwhm / 1.66551; hw1esq = hw1e ** 2
!  nhalf  = hw1e / ABS(dfw) * 2.65
!
!  DO i = 1, nc
!     ! Find the closest pixel
!     if (dfw > 0) then
!        midx = MINVAL(MAXLOC(fwave, MASK=(fwave <= cwave(i)))) + 1
!     else
!        midx = MINVAL(MINLOC(fwave, MASK=(fwave <= cwave(i)))) + 1
!     endif
!
!     sidx = MAX(midx - nhalf, 1)
!     eidx = MIN(nf, midx + nhalf)
!     slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / hw1esq )
!
!     ssum = SUM(slit(sidx:eidx))
!     DO j = 1, nspec
!        cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
!     ENDDO
!  ENDDO
!
!  RETURN
!
!END SUBROUTINE gauss_f2c

! Updates on April 13, 2015
! Make gauss_f2c more generic to handle special cases: points outside the fine spectral grid
SUBROUTINE gauss_f2c (fwave, fspec, nf, nspec, fwhm, cwave, cspec, nc)
  
  IMPLICIT NONE
  
  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                               INTENT (IN) :: nc, nf, nspec
  REAL (KIND=8),                         INTENT (IN) :: fwhm  
  REAL (KIND=8), DIMENSION (nf), INTENT (IN)         :: fwave
  REAL (KIND=8), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=8), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=8), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                       :: i, j, midx, sidx, eidx, nhalf
  REAL (KIND=8)                 :: hw1esq, dfw, ssum, hw1e
  REAL (KIND=8), PARAMETER      :: slit_trunc = 2.65 ! Truncate slit values < 1/1000th of maximum
  REAL (KIND=8), DIMENSION (nf) :: slit

  if (fwhm == 0.0) then
     return
  endif

  dfw  = MIN(fwave(2) - fwave(1), fwave(nf) - fwave(nf-1)) 
  hw1e = fwhm / 1.66551; hw1esq = hw1e ** 2
  nhalf  = hw1e / ABS(dfw) * slit_trunc

  DO i = 1, nc
     ! Find the closest pixel
     if (dfw > 0) then
        midx = MINVAL(MAXLOC(fwave, MASK=(fwave <= cwave(i))))
     else
        midx = MINVAL(MINLOC(fwave, MASK=(fwave <= cwave(i))))
     endif
     IF (midx < 0) midx = 1
     
     IF ( ABS(cwave(i)-fwave(midx)) < slit_trunc * hw1e ) THEN
        sidx = MAX(midx - nhalf, 1)
        eidx = MIN(nf, midx + nhalf)
        slit(sidx:eidx) = EXP(-(cwave(i) - fwave(sidx:eidx))**2 / hw1esq )     
        ssum = SUM(slit(sidx:eidx))
        DO j = 1, nspec
           cspec(i, j) = SUM(fspec(sidx:eidx, j) * slit(sidx:eidx)) / ssum
        ENDDO
     ELSE  ! cwave(i) is too far away from find spectral points
        cspec(i, 1:nspec) = 0.d0
     ENDIF
  ENDDO
  
  RETURN
  
END SUBROUTINE gauss_f2c

! Same as above, except for correcting solar I0 effect using high resolution
! solar reference spectrum following the Beer's Law
! i0: high resolution solar reference
! scalex: scaling, normally number of molecules

! April 13, 2015: does not allow division by zero (i.e., when ci0 = 0.0)
SUBROUTINE gauss_f2ci0(fwave, fspec, i0, nf, nspec, scalex, fwhm, cwave, cspec, nc)

  IMPLICIT NONE

  ! =======================
  ! Input/Output variables
  ! =======================
  INTEGER,                               INTENT (IN) :: nc, nf, nspec
  REAL (KIND=8),                         INTENT (IN) :: fwhm, scalex  
  REAL (KIND=8), DIMENSION (nf), INTENT (IN)         :: fwave, i0
  REAL (KIND=8), DIMENSION (nf, nspec), INTENT (IN)  :: fspec
  REAL (KIND=8), DIMENSION (nc), INTENT (IN)         :: cwave
  REAL (KIND=8), DIMENSION (nc, nspec), INTENT (OUT) :: cspec

  ! ===============
  ! Local variables
  ! ===============
  INTEGER                             :: i, ntemp
  REAL (KIND=8), DIMENSION(nc)        :: ci0
  REAL (KIND=8), DIMENSION(nc, nspec) :: cabspec
  REAL (KIND=8), DIMENSION(nf, nspec) :: abspec

  if (fwhm == 0.0) then
     return
  endif

  ntemp = 1
  CALL gauss_f2c (fwave, i0, nf, ntemp, fwhm, cwave, ci0, nc)
  
  ! Follow Beer's Law
  DO i = 1, nspec
     abspec(:, i) = i0 * EXP(-fspec(:, i) * scalex)
  ENDDO
  
  CALL gauss_f2c (fwave, abspec, nf, nspec, fwhm, cwave, cabspec, nc)

  DO i = 1, nspec
     WHERE (ci0 /= 0.d0)
        cspec(:, i) = - LOG(cabspec(:, i) / ci0) / scalex
     ENDWHERE
  ENDDO

  RETURN

END SUBROUTINE gauss_f2ci0

