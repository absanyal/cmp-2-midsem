REAL(8) FUNCTION indicator(width, x)
  IMPLICIT NONE
  REAL(8) :: x, width
  IF ( ABS(x) .LE. width / 2 ) THEN
     indicator = 0.0
  ELSE
     indicator = 1.0
  END IF
END FUNCTION indicator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM HF
  USE constants
  IMPLICIT NONE
  REAL(8) :: indicator, harmonic, square
  INTEGER :: i, j, k, l, m
  CALL assign_values()

  OPEN(unit = 15, file = 'initial.dat', status = 'unknown')
  OPEN(unit = 25, file = 'final.dat', status = 'unknown')

  !create mesh
  DO i = 0, no_grid
     x(i) = xmin + dx * i
  END DO

  !potential
  DO i = 0, no_grid
    v(i) = 10**8 * indicator(width,x(i)) !square well
    !v(i) = 0.5 * mass * omega**2.0 * (x(i))**2.0 !harmonic
  END DO

  !form hamiltonian
  vhf = 0.0

  DO k = 1, no_loop !main loop for iteration

    h(:,:) = 0.0
    DO i = 0, no_grid
      DO j = 0, no_grid
        IF (i .EQ. j) h(i, j) = (1.0 / (dx ** 2.0))+ v(i) + vhf(i)
        IF((j == i + 1).OR.(j == i - 1)) &
        h(i, j)= -1.0 / ( 2.0 * dx ** 2.0)
        h(0, no_grid)= -1.0 / ( 2.0 * dx ** 2.0)
        h(no_grid, 0)= -1.0 / ( 2.0 * dx ** 2.0)
      END DO
    END DO

    !diagonalization
    CALL dsyevd('V','U',size,h,size,e,work,lwork,iwork,liwork,info)
    IF(info /= 0) STOP 'matrix diagonalization failed'

    !rho_h
    rho = 0.0
    DO m=0,no_p - 1
      DO j=1,no_grid
        rho(j) = rho(j) + h(j,m) * h(j,m)
      END DO
    END DO

    !rho_hf
    rhohf=0.00
    DO i = 0, no_grid
      DO j = 0, no_grid
        DO l = 0, no_p - 1
          DO m = 0, no_p - 1
            IF( rho(i) .NE. 0.0000) THEN
              rhohf(i,j) = &
              rhohf(i, j) + h(j, m) * h(j, l) * h(i, l) * h(i, m) / (rho(i))
            END IF
          END DO
        END DO
      END DO
    END DO

    !HF potential
    vhf = 0.0
    DO i = 0, no_grid
      DO j = 0, no_grid
        vhf(i) = &
        vhf(i) + (rho(j) - rhohf(i, j))* dx / (ABS(x(i) - x(j)) + (dx / 2.0) )
      END DO
    END DO

    !file output
    IF (k .EQ. 1) THEN
      DO j = 0, no_grid
        WRITE(15, *) x(j), h(j,0)**2.0, h(j,1)**2.0
      END DO
    ELSE
      IF (k .EQ. no_loop) THEN
        DO j = 0, no_grid
          WRITE(25, *) x(j), h(j,0)**2.0, h(j,1)**2.0
        END DO
      END IF
    END IF
    PRINT *, k, e(0), e(1)

  END DO !main loop ends

END PROGRAM HF
