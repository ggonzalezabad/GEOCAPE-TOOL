PROGRAM test_get_hitran

INTEGER, PARAMETER :: maxz = 60, maxlambda=50001
CHARACTER(LEN=6)   :: the_molecule
INTEGER            :: nlambda, nz, errstat
LOGICAL            :: is_wavenum
REAL(KIND=8)       :: fwhm
REAL(KIND=8), DIMENSION(maxlambda)       :: lambda
REAL(KIND=8), DIMENSION(maxz)            :: ps, ts
REAL(KIND=8), DIMENSION(maxlambda, maxz) :: crs

the_molecule ='H2O   '
nz = 2
ps(1) = 1.0; ts(1) = 300.0
ps(2) = 1.0; ts(2) = 300.0
is_wavenum = 0

nlambda = 5057
DO i = 1, nlambda
   lambda(i) = 757.5 + 0.003 * (i - 1)
ENDDO
fwhm = 0.0
CALL get_hitran_crs(the_molecule, nlambda, lambda(1:nlambda), is_wavenum, &
     nz, ps(1:2), ts(1:2), fwhm, crs(1:nlambda, 1:2), errstat)

IF (errstat == 0) THEN
   WRITE(*, *) 'Get HITRAN Cross Section sucessfully!!!'
ELSE
   WRITE(*, *) 'Error in getting HITRAN Cross Section!!!'
ENDIF

STOP


END PROGRAM test_get_hitran
