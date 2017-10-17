RECURSIVE FUNCTION factorial(x) RESULT(f)
  IMPLICIT NONE
  INTEGER :: x
  INTEGER :: f
  IF (x == 0) THEN
     f = 1
  ELSE
     f = x * factorial(x-1)
  END IF
END FUNCTION factorial

REAL FUNCTION HermitePoly(x, n)
  IMPLICIT NONE
  REAL :: x
  INTEGER :: n
  IF ( n .EQ. 0) THEN
     HermitePoly = 1.0
  END IF
  IF ( n .EQ. 1 ) THEN
     HermitePoly = 2.0 * x
  END IF
END FUNCTION HermitePoly

REAL FUNCTION HOS(x, n, omega)
  IMPLICIT NONE
  REAL :: mass, hbar, alpha, HermitePoly, x, omega
  INTEGER :: n
  REAL, PARAMETER :: pi = 3.14159265
  INTEGER :: factorial
  mass = 1.0
  hbar = 1.0
  alpha = mass * omega / hbar
  HOS = ((1.0/( (2.0**n) * factorial(n) ) ) ** 0.5) * ( (alpha/pi) ** (0.25) )  &
   * EXP( -alpha * (x ** 2.0) / 2.0) * HermitePoly((alpha ** 0.5) * x, n)
END FUNCTION HOS

REAL FUNCTION psi(x1, x2, omega)
  IMPLICIT NONE
  REAL :: x1, x2, omega, HOS
  psi =  (1.0/(2.0 ** (0.5))) * ( HOS(x1,0,omega) * HOS(x2,1,omega) &
  -HOS(x2,0,omega) * HOS(x1,1,omega) )
END FUNCTION psi

REAL FUNCTION E_c(x1, x2, omega)
  IMPLICIT NONE
  REAL :: x1, x2, omega, psi
  REAL, PARAMETER :: A = 1.0 / 2000.0
  E_c = psi(x1, x2, omega) * (1.0 / ( ABS(x1 - x2) + A )) * psi(x1, x2, omega)
END FUNCTION E_c

REAL FUNCTION integrate_E_c(lowlimx, uplimx, lowlimy, uplimy, no_of_gridx, no_of_gridy, omega)
  IMPLICIT NONE
  REAL :: s, hx, hy, lowlimx, uplimx, lowlimy, uplimy, xi, yj, E_c, omega
  INTEGER :: i, j, no_of_gridx, no_of_gridy
  hx = (uplimx - lowlimx) / no_of_gridx
  hy = (uplimy - lowlimy) / no_of_gridy
  s = 0
  DO i = 0, no_of_gridx
    DO j = 0, no_of_gridy
     xi = lowlimx + hx/2 + i * hx
     yj = lowlimy + hy/2 + j * hy
     s = s + hx * hy * E_c(xi, yj, omega)
   END DO
  END DO
  integrate_E_c = s
END FUNCTION integrate_E_c

PROGRAM main
  IMPLICIT NONE
  REAL :: omega, integrate_E_c
  REAL, PARAMETER :: lowlim = -5.0, uplim = 5.0
  INTEGER, PARAMETER :: grid = 1000
  INTEGER, PARAMETER :: out_unit = 20
  REAL :: omega_start, omega_stop, omega_step
  PRINT *, 'Enter the start, stop and step-size for omega:'
  READ *, omega_start, omega_stop, omega_step
  omega = omega_start
  OPEN (unit = out_unit, file = "Q2_pert.dat", &
   action = "write", status = "unknown")
  DO WHILE ( omega .LE. omega_stop)
    PRINT *, omega, &
     integrate_E_c(lowlim, uplim, lowlim, uplim, grid, grid, omega)
    WRITE (out_unit,*) omega, &
     integrate_E_c(lowlim, uplim, lowlim, uplim, grid, grid, omega)
    omega = omega + omega_step
  END DO
  CLOSE (out_unit)
END PROGRAM main
