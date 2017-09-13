
!> Computes the inverse cumulative density function (CDF), i.e., the quantile,
! of the standard normal distribution given u uniform on the unit hypercube.
FUNCTION beasley_springer_moro(u) RESULT(x)
    IMPLICIT NONE

    INTEGER :: j
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u

    DOUBLE PRECISION :: r
    DOUBLE PRECISION, DIMENSION(SIZE(u)) :: x, y


    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: a = (/ &
            2.50662823884, &
            -18.61500062529, &
            41.39119773534, &
            -25.44106049637 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: b = (/ &
            -8.47351093090, &
            23.08336743743, &
            -21.06224101826, &
            3.13082909833 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:8) :: c = (/ &
            0.3374754822726147, &
            0.9761690190917186, &
            0.1607979714918209, &
            0.0276438810333863, &
            0.0038405729373609, &
            0.0003951896511919, &
            0.0000321767881768, &
            0.0000002888167364, &
            0.0000003960315187 /)


    y = u - 0.5D0

    DO j = 1, SIZE(u)
        IF (ABS(y(j)) < 0.42) THEN
            r = y(j)*y(j)
            x(j) = y(j)*(((a(3)*r + a(2))*r + a(1))*r + a(0))/((((b(3)*r + b(2))*r + b(1))*r + b(0))*r + 1)
        ELSE
            IF (y(j) > 0) THEN
                r = LOG(-LOG(1-u(j)))
            ELSE IF (y(j) < 0) THEN
                r = LOG(-LOG(u(j)))
            END IF
            x(j) = c(0) + r*(c(1) + r*(c(2) + r*(c(3) + r*(c(4) + r*(c(5) + r*(c(6) + r*(c(7) + r*c(8))))))))
            IF (y(j) < 0) THEN
                x(j) = -x(j)
            END IF
        END IF
    END DO

END FUNCTION beasley_springer_moro

!==============================================================================!


!> Returns a d-dimensional Sobol sequence of p points following a standard
!  normal distribution
SUBROUTINE sobol_stdnormal(d, skip, x_stdnormal)
  USE sobol
  IMPLICIT NONE
  !> dimension
  INTEGER(kind = 4), INTENT(IN) :: d

  !> number of initial points to be skipped
  INTEGER(kind = 8), INTENT(IN) :: skip   

  !> return an array of doubles, standard normal
  DOUBLE PRECISION, DIMENSION(d), INTENT(OUT) :: x_stdnormal     

  INTERFACE
     FUNCTION beasley_springer_moro(u)
       DOUBLE PRECISION :: u(:)
       DOUBLE PRECISION :: beasley_springer_moro(SIZE(u))
     END FUNCTION beasley_springer_moro
    END INTERFACE

    x_stdnormal = beasley_springer_moro(i8_sobol(INT(d, 8), skip))

END SUBROUTINE sobol_stdnormal
