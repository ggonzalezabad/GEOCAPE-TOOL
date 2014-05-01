! This subroutine combines spline and splint function
! in Numberical Recipes by Press et al., 1997.
SUBROUTINE BSPLINE(xa, ya, n, x, y, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y

  
  REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
  REAL (KIND=dp), DIMENSION(np)             :: xcpy
  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
  
  errstat = 0
  IF (n < 3) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF

  xmin = MINVAL(x); xmax = MAXVAL(x)
  xamin = MINVAL(xa); xamax = MAXVAL(xa)
  IF (xmin < xamin .OR. xmax > xamax) THEN
     errstat =  -3; RETURN
  ENDIF
  
  IF (xa(1) < xa(n)) THEN
     CALL SPLINE1(xa, ya, n, y2a)
     CALL SPLINT1(xa, ya, y2a, n, x, y, np)
  ELSE
     xacpy = -xa; xcpy = -x
     CALL SPLINE1(xacpy, ya, n, y2a)
     CALL SPLINT1(xacpy, ya, y2a, n, xcpy, y, np)
  ENDIF

  RETURN
END SUBROUTINE BSPLINE

! modified to always use "natural" boundary conditions
SUBROUTINE SPLINE1 (x, y, n, y2)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: x, y  
  REAL (KIND=dp), DIMENSION(n), INTENT(OUT) :: y2
  
  REAL (KIND=dp), DIMENSION(n)     		 :: u
  INTEGER       :: i, k
  REAL(KIND=dp) :: sig, p, qn, un
  
  y2 (1) = 0.0
  u (1) = 0.0
  
  DO i = 2, n - 1
     sig = (x (i) - x (i - 1)) / (x (i + 1) -x (i - 1))
     p = sig * y2 (i - 1) + 2.D0
     y2 (i) = (sig - 1.) / p
     u (i) = (6._dp * ((y (i + 1) - y (i)) / (x (i + 1) - x (i)) -  & 
          (y (i) - y (i - 1)) / (x (i) - x (i - 1))) / (x (i + 1) - &
          x (i - 1)) - sig * u (i - 1)) / p
  ENDDO
  
  qn = 0.0
  un = 0.0
  y2 (n) = (un - qn * u (n - 1)) / (qn * y2 (n - 1) + 1.0)
  DO k = n - 1, 1, -1
     y2 (k) = y2 (k) * y2 (k + 1) + u (k)
  ENDDO
  
  RETURN
END SUBROUTINE SPLINE1

! This code could be optimized if x is in increasing/descreasing order
SUBROUTINE SPLINT1 (xa, ya, y2a, n, x, y, m)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n, m
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
  REAL (KIND=dp), DIMENSION(m), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(m), INTENT(OUT):: y
  
  INTEGER        :: ii, klo, khi, k 
  REAL (KIND=dp) :: h, a, b

  !klo = 1; khi = n
  DO ii = 1, m 
     klo = 1; khi = n
    
     !IF ( khi - klo == 1) THEN
     !   IF (x(ii) > xa(khi) )   THEN
     !       khi = n
     !   ENDIF
     !ENDIF

     DO WHILE (khi - klo > 1)
        k = (khi + klo) / 2
        IF (xa (k) > x(ii)) THEN
           khi = k
        ELSE
           klo = k
        ENDIF
     ENDDO
     
     h = xa (khi) - xa (klo)
     IF (h == 0.0) STOP 'Bad xa input in: splint!!!'
     a = (xa (khi) - x(ii)) / h
     b = (x(ii) - xa (klo)) / h
     
     y(ii) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + &
          (b**3 - b) * y2a (khi)) * (h**2) / 6.0
  ENDDO
  
  RETURN
END SUBROUTINE SPLINT1
