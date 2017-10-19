REAL(8) FUNCTION indicator(width, x)
  IMPLICIT NONE
  REAL(8) :: x, width
  IF ( ABS(x) .LE. 0.5 ) THEN
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
  INTEGER :: i, j, k
  CALL assign_values()

  !create mesh
  DO i = 0, no_grid
     x(i) = xmin + dx * i
  END DO

  !potential
  DO i = 0, no_grid
    v(i) = 10**8 * indicator(width,x(i))
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
      END DO
    END DO

  END DO !main loop ends

END PROGRAM HF
