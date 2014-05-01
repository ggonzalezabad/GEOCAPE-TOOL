! ==================================================================================!
! Purpose: This routine provides an interface to Christopher E. Sioris's single     !
! scattering first-order rotational raman scattering model (Sioris and Evans, 2001) !
! Xiong Liu, 07/08/2007                                                             !
!                                                                                   !
! Inputs/Outputs:                                                                   !
! nz: number of atmospheric layers                                                  !
! nw: number of wavelengths                                                         !
! sza: solar zenith angle  (degree)                                                 !
! vza: viewing zenith angle (degree)                                                !
! sca: scattering angle (degree, backscattering: 180 )                              !
! albedo: Lambertian surface albedo (only for satellite observations)               !
! do_upwelling: true for satellite observations else: ground-based measurements     !
! ts:  temperature profile fro TOA to BOS (K)                                       !
! wave: wavelengths (nm), increasing order                                          !
! sol:  solar irradiance                                                            !
! rhos: number of air molecules at each layer (molecules cm^-2)                     !
! taus: optical thickness at each wavelength and each layer                         !
! rspec: Ring spectrum at wave (i.e., relative diff. with and without Ring effect)  !
! problems: if true, errors occur during computing the Ring spectrum                !
!                                                                                   !
! Notes:                                                                            !
! (1) Ring effect is not computed (i.e., 0) for the few nms (~3nm) at each end. So  !
!     if you want to calculate Ring effect for 310-330 nm, you need to provide      !
!     sol/taus for 307-333 nm.                                                      !
! (2) If taus and sol are at high spectral resolution, the Ring effect is also at   !
!     high spectral resolution and you need to convolve and sample it to low        !
!     resolution. On the other hand, if taus and sol are at low resolution          !
!     (e.g. OMI), so is the ring spectrum. You just need to sample/interpolate it   !
!     to the grid of you want.                                                      !
! ==================================================================================! 

SUBROUTINE GET_RAMAN (nz, nw, sza, vza, sca, albedo, do_upwelling, ts, rhos, &
     wave, sol, taus, rspec, problems)
 
  IMPLICIT NONE
  INTEGER, PARAMETER         :: dp = KIND(1.0D0) 

  ! ========================
  ! Input/output variables
  ! ========================
  INTEGER,        INTENT(IN) :: nz, nw
  LOGICAL,        INTENT(IN) :: do_upwelling
  LOGICAL,       INTENT(OUT) :: problems
  REAL (KIND=dp), INTENT(IN) :: sza, sca, vza, albedo
  REAL (KIND=dp), DIMENSION(nz), INTENT(IN)     :: ts, rhos
  REAL (KIND=dp), DIMENSION(nw), INTENT(IN)     :: wave, sol
  REAL (KIND=dp), DIMENSION(nw, nz), INTENT(IN) :: taus
  REAL (KIND=dp), DIMENSION(nw),    INTENT(OUT) :: rspec

  ! ========================
  ! Local Variables
  ! ========================
  INTEGER,        PARAMETER :: MAXNU   = 75000, NedgePos = 218
  REAL (KIND=dp), PARAMETER :: pi      = 3.14159265358979_dp
  REAL (KIND=dp), PARAMETER :: deg2rad = pi / 180.0_dp

  INTEGER                              :: nuhi, nulo, nu, i, j, errstat
  REAL (KIND=dp)                       :: scl, cosvza, cossza, tmpalb
  REAL (KIND=dp), DIMENSION(nz)        :: ctau
  REAL (KIND=dp), DIMENSION(nw, nz)    :: strans, vtrans
  REAL (KIND=dp), DIMENSION(MAXNU, nz) :: st, vt
  REAL (KIND=dp), DIMENSION(MAXNU)     :: ring, ramanwav

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=9), PARAMETER :: modulename = 'GET_RAMAN'
  
  problems = .FALSE.
  IF ( .NOT. do_upwelling) THEN
     ! Ignore the surface contribution for ground-based observations.
     tmpalb = 0.0
  ELSE
     tmpalb = albedo
  ENDIF
  
  ! Get position for raman calculation
  nuhi = INT(1.0D7 / wave(1))
  nulo = INT(1.0D7 / wave(nw)) + 1
  nu = 0
  DO i = nulo, nuhi
     nu = nu + 1
     ramanwav(nu) = 1.0D7 / REAL(i, KIND=dp)
  ENDDO
  ramanwav(nu) = wave(1); ramanwav(1) = wave(nw)
  
  IF (nuhi - nulo + 1 > MAXNU) THEN
     WRITE(*, *) modulename, ': Need to increase MAXNU!!!'
     errstat = .TRUE.; RETURN
  ELSE IF (nuhi <= nulo) THEN
     WRITE(*, *) modulename, ': nulo>=nuhi, should never happen!!!'
     errstat = .TRUE.; RETURN
  ENDIF

  cossza = COS(sza * deg2rad); cosvza = COS(vza * deg2rad)

  ! Compute optical depth
  DO i = 1, nw 
     ctau = 0.0
     DO j = 2, nz
        ctau(j) = ctau(j-1) + taus(i, j)
     ENDDO
     
     scl = sol(i) * ((wave(i) / wave(1)) ** (-4.0))

     ! Assume plane parallel. In spherical geometry, you can replace 
     ! ctau/cossza with slant optical thickness
     strans(i, 1:nz) = scl * EXP(-ctau(1:nz) / cossza)
     IF (do_upwelling) THEN
        vtrans(i, 1:nz) = EXP(-ctau(1:nz) / cosvza)
     ELSE
        vtrans(i, 1:nz) = EXP(-(ctau(nz) - ctau(1:nz)) / cosvza)
     ENDIF
  ENDDO

  ! Interpolate to raman grid in wavenumber
  DO i = 1, nz
     CALL BSPLINE(wave, strans(:, i), nw, ramanwav(1:nu), st(1:nu, i), nu, errstat)     
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
        problems = .TRUE.; RETURN
     ENDIF

     CALL BSPLINE(wave, vtrans(:, i), nw, ramanwav(1:nu), vt(1:nu, i), nu, errstat)
     IF (errstat < 0) THEN
        WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
        problems = .TRUE.; RETURN
     ENDIF
  ENDDO
  
  ! Call raman program
  CALL RAMAN(nulo, nuhi, nu, nz, sca, tmpalb, ts, rhos, st(1:nu,1:nz), vt(1:nu,1:nz), ring(1:nu))
 
  ! Interpolate calculated ring back to gome radiance grids
  CALL BSPLINE(ramanwav(1:nu), ring(1:nu), nu, wave(1:nw), rspec(1:nw), nw, errstat)
    
  ! Set edge pixels to zero
  DO i = 1, nw
     IF (wave(i) >= ramanwav(nu-NedgePos) ) THEN
        EXIT
     ELSE
        rspec(i) = 0.D0
     ENDIF
  ENDDO

  DO i = nw, 1, -1
     IF (wave(i) <= ramanwav(NedgePos+1) ) THEN
        EXIT
     ELSE
        rspec(i) = 0.D0
     ENDIF
  ENDDO

  IF (errstat < 0) THEN
     WRITE(*, *) modulename, ': BSPLINE error, errstat = ', errstat
     problems = .TRUE.; RETURN
  ENDIF
  
  RETURN
  
END SUBROUTINE GET_RAMAN



