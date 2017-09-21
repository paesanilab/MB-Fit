
!> Computes the inverse cumulative density function (CDF), i.e., the quantile,
! of the standard normal distribution given u uniform on the unit hypercube.
FUNCTION beasley_springer_moro_r(u) RESULT(x)
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

END FUNCTION beasley_springer_moro_r

!==============================================================================!

SUBROUTINE init_random_seed()
  ! "Initialize a pseudo-random number sequence"
  ! source:  http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED

  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: seed(:)
  INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
  INTEGER(8) :: count, tms

  CALL random_SEED(size = n)
  ALLOCATE(seed(n))
  ! First try if the OS provides a random number generator
  OPEN(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  IF (istat == 0) THEN
     READ(un) seed
     CLOSE(un)
  ELSE
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     CALL system_CLOCK(count)
     IF (count /= 0) THEN
        t = TRANSFER(count, t)
     ELSE
        CALL date_and_TIME(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = TRANSFER(tms, t)
     END IF
     s = IEOR(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = IEOR(s, pid)
     IF (n >= 3) THEN
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        IF (n > 3) THEN
           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        END IF
     ELSE
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     END IF
  END IF
  CALL random_SEED(put=seed)
END SUBROUTINE init_random_seed

!==============================================================================!

!> Returns a d-dimensional Sobol sequence of p points following a standard
!  normal distribution
SUBROUTINE random_stdnormal(d, lseed, x_stdnormal)
  !use sobol
  IMPLICIT NONE
  !> dimension
  INTEGER(kind = 4), INTENT(IN) :: d

  !> boolean set to .TRUE. to initialize seed for random_number
  LOGICAL :: lseed

  !> pseudorandom d-dimensional point to be transformed to standard normal
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r

  !> return an array of doubles, standard normal
  DOUBLE PRECISION, DIMENSION(d), INTENT(OUT) :: x_stdnormal     

  INTERFACE
     FUNCTION beasley_springer_moro_r(u)
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u
       DOUBLE PRECISION :: beasley_springer_moro_r(SIZE(u))
     END FUNCTION beasley_springer_moro_r
  END INTERFACE

  ALLOCATE(r(d))

  IF ( lseed .EQV. .TRUE.) THEN
     CALL init_random_seed()
  ENDIF
  CALL random_NUMBER(r)
  x_stdnormal = beasley_springer_moro_r(r)

  DEALLOCATE(r)


END SUBROUTINE random_stdnormal
