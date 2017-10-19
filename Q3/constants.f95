MODULE constants
  IMPLICIT NONE

  !global constants
  INTEGER :: no_grid = 1000, no_p = 2, no_loop = 5
  INTEGER :: size !size of matrix
  REAL(8) :: xmin = -5.0, xmax = 5.0, dx !interval related
  INTEGER, ALLOCATABLE::iwork(:)
  REAL(8), ALLOCATABLE:: work(:)
  REAL(8), ALLOCATABLE :: h(:,:), e(:), v(:), x(:), &
       rho(:),rhohf(:,:),vhf(:)
  INTEGER :: lwork,info,liwork

  !physical properties
  REAL(8), PARAMETER :: mass = 1.0, hbar = 1.0
  REAL(8) :: omega = 1.0
  REAL(8) :: width = 1

CONTAINS

  SUBROUTINE assign_values()
    dx = (xmax - xmin) / no_grid
    size = no_grid + 1
    liwork = 3 + 5 * size
    lwork= 1 + 6 * size + 2 * size ** 2
    ALLOCATE (&
         work(1 + 6 * size + 2 * size**2),&
         e(0:no_grid), h(0:no_grid, 0:no_grid), iwork(3 + 5 * size), &
         v(0:no_grid), x(0:no_grid), rho(0:no_grid), &
         rhohf(0:no_grid, 0:no_grid), vhf(0:no_grid) &
         )
  END SUBROUTINE assign_values

END MODULE constants
