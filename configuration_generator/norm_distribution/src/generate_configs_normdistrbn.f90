
MODULE generate_configs_module

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: deg=180/dacos(-1.0d0)
  DOUBLE PRECISION, PARAMETER :: pi=dacos(-1.0d0)
  DOUBLE PRECISION, PARAMETER :: hbar = 1.0d0
  DOUBLE PRECISION, PARAMETER :: bohr = 0.52917721092
  DOUBLE PRECISION, PARAMETER :: autoeV = 27.211385
  DOUBLE PRECISION, PARAMETER :: eVtoau = 3.674932379d-2
  DOUBLE PRECISION, PARAMETER :: autocm = 2.194746313d5
  DOUBLE PRECISION, PARAMETER :: autokcalmol = 627.5096
  DOUBLE PRECISION, PARAMETER :: Ktoau = 3.1668114d-6
  DOUBLE PRECISION, PARAMETER :: cmtoau = 4.5563352527d-6
  DOUBLE PRECISION, PARAMETER :: melectron=1822.88839

  DOUBLE PRECISION, PARAMETER :: Hmass = 1.00782503223*melectron
  DOUBLE PRECISION, PARAMETER :: Dmass = 2.01410177812*melectron
  DOUBLE PRECISION, PARAMETER :: Tmass = 3.0160492779*melectron
  DOUBLE PRECISION, PARAMETER :: Bmass = 11.00930536*melectron
  DOUBLE PRECISION, PARAMETER :: Cmass = 12.0000000*melectron
  DOUBLE PRECISION, PARAMETER :: Nmass = 14.00307400443*melectron
  DOUBLE PRECISION, PARAMETER :: Omass = 15.99491461957*melectron
  DOUBLE PRECISION, PARAMETER :: Fmass = 18.99840316273*melectron
  DOUBLE PRECISION, PARAMETER :: Pmass = 30.97376199842*melectron
  DOUBLE PRECISION, PARAMETER :: Smass = 31.9720711744*melectron
  DOUBLE PRECISION, PARAMETER :: Clmass = 34.968852682*melectron
  DOUBLE PRECISION, PARAMETER :: Brmass = 78.9183376*melectron
  DOUBLE PRECISION, PARAMETER :: Imass = 126.9044719*melectron
  DOUBLE PRECISION, PARAMETER :: Limass = 7.0160034366*melectron
  DOUBLE PRECISION, PARAMETER :: Namass = 22.989769282*melectron
  DOUBLE PRECISION, PARAMETER :: Kmass  = 38.9637064864*melectron
  DOUBLE PRECISION, PARAMETER :: Rbmass = 84.9114977282*melectron
  DOUBLE PRECISION, PARAMETER :: Csmass = 132.905429*melectron

  DOUBLE PRECISION, ALLOCATABLE :: sqrt_mass(:)
  CHARACTER(len=2), ALLOCATABLE :: atom_type(:)

CONTAINS

  FUNCTION atom_mass(atom)

    IMPLICIT NONE

    CHARACTER(len=2), INTENT(in) :: atom
    DOUBLE PRECISION :: atom_mass


    IF (atom == 'H' .OR. atom == '1') THEN
       atom_mass = Hmass
    ELSE IF (atom == 'D') THEN
       atom_mass = Dmass
    ELSE IF (atom == 'T') THEN
       atom_mass = Tmass
    ELSE IF (atom == 'B' .OR. atom == '5') THEN
       atom_mass = Bmass
    ELSE IF (atom == 'C' .OR. atom == '6') THEN
       atom_mass = Cmass
    ELSE IF (atom == 'N' .OR. atom == '7') THEN
       atom_mass = Nmass
    ELSE IF (atom == 'O' .OR. atom == '8') THEN
       atom_mass = Omass
    ELSE IF (atom == 'F' .OR. atom == '9') THEN
       atom_mass = Fmass
    ELSE IF (atom == 'P' .OR. atom == '15') THEN
       atom_mass = Pmass
    ELSE IF (atom == 'S' .OR. atom == '16') THEN
       atom_mass = Smass
    ELSE IF (atom == 'Cl' .OR. atom == '17') THEN
       atom_mass = Clmass
    ELSE IF (atom == 'Br' .OR. atom == '35') THEN
       atom_mass = Brmass
    ELSE IF (atom == 'I' .OR. atom == '53') THEN
       atom_mass = Imass
    ELSE IF (atom == 'Li' .OR. atom == '3') THEN
       atom_mass = Limass
    ELSE IF (atom == 'Na' .OR. atom == '11') THEN
       atom_mass = Namass
    ELSE IF (atom == 'K' .OR. atom == '19') THEN
       atom_mass = Kmass
    ELSE IF (atom == 'Rb' .OR. atom == '37') THEN
       atom_mass = Rbmass
    ELSE IF (atom == 'Cs' .OR. atom == '55') THEN
       atom_mass = Csmass
    ELSE
       WRITE(*,*) 'atom ', atom, ' is not recognized'
       STOP 
    END IF

  END FUNCTION atom_mass

END MODULE generate_configs_module

!=====================================================================!

PROGRAM generate_configs_normdistrbn

  USE generate_configs_module

  IMPLICIT NONE

  INTEGER :: i, j, k, m, n
  INTEGER :: dim, dimnull                     ! full dimension of coordinate space (dim = 3*natoms); dimension of nullspace
  INTEGER :: nconfigs                         ! number of configurations to be generated
  INTEGER(kind=8) :: skip                     ! number of initial points to skip when using Sobol's sequence
  INTEGER :: maxindexT, maxindexA, nmax       ! for incremental progression of T, A
  DOUBLE PRECISION :: Tmin, Tmax, factorT     ! for incremental progression of T
  DOUBLE PRECISION :: Amin, Amax, factorA     ! for incremental progression of A
  DOUBLE PRECISION, ALLOCATABLE :: T(:)       ! temperatures, indexed over configurations
  DOUBLE PRECISION, ALLOCATABLE :: A(:)       ! A := kBT/hw, indexed over configurations; kB is Boltzmann's constant, h is hbar
  DOUBLE PRECISION, ALLOCATABLE :: q(:)       ! reference geometry coordinates (Cartesian)
  DOUBLE PRECISION, ALLOCATABLE :: r(:)       ! displacment from reference geometry (Cartesian)
  DOUBLE PRECISION, ALLOCATABLE :: rn(:)      ! distorted/displaced coordinates (Cartesian)
  DOUBLE PRECISION, ALLOCATABLE :: w(:)       ! eigenfrequencies, ordered from least to greatest
  DOUBLE PRECISION, ALLOCATABLE :: d(:)       ! d = d(w,T); defines mass-scaled covariance matrix D=d*U*U^T
  DOUBLE PRECISION, ALLOCATABLE :: U(:,:)     ! columns are the (orthonormal) eigenvectors/normal modes, ordered;
                                              ! mass-weighted, Cartesian coordinates
  DOUBLE PRECISION, ALLOCATABLE :: G(:,:)     ! square root of the mass-scaled covariance matrix D
  DOUBLE PRECISION, ALLOCATABLE :: C(:,:)     ! to confirm orthonormality of the normal modes; C := U^T*U

  DOUBLE PRECISION :: freq_cutoff=10*cmtoau   ! ensures finite breadth of sampling

  CHARACTER(len=100) :: input_xyz             ! file containing reference geometry, xyz file format
  CHARACTER(len=100) :: input_normal_modes    ! file containing eigenfrequencies, normal modes (ordered)
  CHARACTER(len=100) :: output_xyz            ! file containing generated configurations
  CHARACTER(len=1) :: sobolchar               ! 'P' for pseudorandom sequence,'Q' for quasirandom (Sobol) sequence
  LOGICAL :: lseed = .TRUE.                   ! '.TRUE.' for new seed for random number generator for each run
  LOGICAL :: geometric = .FALSE.              ! use a geometric progression of T, A
  LOGICAL :: linear = .FALSE.                 ! use a linear progression of T, A
  LOGICAL :: verbose = .FALSE.                ! '.TRUE.' for verbose printing containing d(w,T) values

  REAL :: initial_time, final_time

  CALL CPU_TIME(initial_time)



  ! read main input
  READ(*,*) input_xyz
  READ(*,*) input_normal_modes
  READ(*,*) dim, dimnull
  READ(*,*) sobolchar, nconfigs
  READ(*,*) output_xyz
  READ(*,*) geometric, linear
  READ(*,*) verbose

  IF (dimnull >= dim) THEN
     WRITE(*,*) 'Dimension of nullspace must be less than total dimension.'
     STOP
  END IF
  IF (nconfigs < 2) THEN
     WRITE(*,*) 'Number of configurations to generate must be greater than 1.'
     STOP
  ENDIF

  ALLOCATE(atom_type(dim/3), sqrt_mass(dim), q(dim), r(dim), rn(dim))
  ALLOCATE(w(dim), d(dim))
  ALLOCATE(G(dim, dim), U(dim, dim-dimnull), C(dim-dimnull, dim-dimnull))


  ! read reference geometry/configuration
  OPEN(unit = 76, file = input_xyz, action = 'read', status = 'old')
  READ(76,*)
  READ(76,*)
  DO i = 1, dim/3
     READ(76,*) atom_type(i), q(3*i-2), q(3*i-1), q(3*i) ! Angstrom
     sqrt_mass(3*i-2:3*i) = SQRT(atom_mass(atom_type(i)))
  END DO
  CLOSE(unit=76)

  ! read eigenfrequencies, normal modes
  U(:,:) = 0.0d0
  OPEN(unit = 54, file = input_normal_modes, status = 'old', action = 'read')
  DO k = 1, dim-dimnull   ! vibrational degrees of freedom
     READ(54,*)           ! normal mode: k
     READ(54,*) w(k)      ! frequency in cm^(-1)
     READ(54,*)           ! reduced mass
     DO i = 1, dim/3
        READ(54,*) U(i*3-2,k), U(i*3-1,k), U(i*3,k) ! mass-weighted Cartesian coordinates
     ENDDO
     READ(54,*)
  ENDDO
  CLOSE(unit=54)

  ! generate linear or geometric progression of T, A:
  maxindexA = nconfigs/2-1
  maxindexT = maxindexA + MOD(nconfigs,2)

  ALLOCATE(T(0:maxindexT))
  ALLOCATE(A(0:maxindexA))

  IF (geometric .EQV. .TRUE.) THEN
     Tmin = w(1)/autocm                 ! assume lowest frequency is first
     Tmax = 2.0d0*w(dim-dimnull)/autocm ! assume highest frequency is last
     Amin = 1.0d0
     Amax = 2.0d0
     factorT = (Tmin/Tmax)**(-1.0d0/maxindexT)
     factorA = (Amin/Amax)**(-1.0d0/maxindexA)
     WRITE(*,*) 'factorT: ', factorT
     WRITE(*,*) 'factorA: ', factorA
     T(0) = Tmin
     A(0) = Amin
     DO n = 1, maxindexT
        T(n) = T(n-1)*factorT
     END DO
     DO n = 1, maxindexA
        A(n) = A(n-1)*factorA
     END DO
  ELSE IF (linear .EQV. .TRUE.) THEN
     Tmin = 0.0d0                 ! (computed d(k) values will be subject to 'freq_cutoff' variable)
     Tmax = w(dim-dimnull)/autocm ! assume highest frequency is last
     Amin = 0.0d0                 ! (computed d(k) values will be subject to 'freq_cutoff' variable)
     Amax = 2.0d0
     factorT = (Tmax-Tmin)/maxindexT
     factorA = (Amax-Amin)/maxindexA
     T(0) = Tmin
     A(0) = Amin
     DO n = 1, maxindexT
        T(n) = T(n-1) + factorT
     END DO
     DO n = 1, maxindexA
        A(n) = A(n-1) + factorA
     END DO
  ELSE
     WRITE(*,*) 'Choose geometric or linear temperature progression.'
     STOP
  END IF

  ! convert to atomic units:
  !T = T*Ktoau ! for later: adding single-temp option
  q(:) = q(:)/bohr
  w(:) = w(:)/autocm

  ! unscale the normal modes:
  DO k=1, dim-dimnull
     U(:,k)=U(:,k)*sqrt_mass(:)    
     U(:,k)=U(:,k)/SQRT(SUM(U(:,k)**2))    
  ENDDO

  ! confirm that the normal modes are orthonormal:
  C(:,:) = 0.0d0
  WRITE(*,*) "U^T*U: "
  C=MATMUL(TRANSPOSE(U),U)
  DO i=1,dim-dimnull
     WRITE(*,'(*(F11.4))') C(i,:)
  END DO
  WRITE(*,*)

  ! write eigenfrequencies in wavenumbers and equivalent temperature T=hw/kB:
  DO k = 1, dim-dimnull
     WRITE(*,'(A,I0,A,F11.4,A,F11.4,A)') 'w[',k,'] = ', w(k)*autocm, ' cm^{-1} ~ ', w(k)/Ktoau, ' kelvin'
  ENDDO
  WRITE(*,*)

  ! for Sobol's sequence, skip largest power of 2 less than or equal to nconfigs:
  skip = 2_8
  DO WHILE (skip*2 <= INT(nconfigs,8))
     skip = skip*2
  ENDDO


  ! generate configurations:
  w(:) = dabs(w(:))                                    ! to accomodate imaginary frequencies
  OPEN(unit = 89, file = output_xyz, STATUS='UNKNOWN') ! will overwrite file if it exists
  DO m = 0, 1                                          ! sample quantum harmonic distributions (m = 0 case)
                                                       ! followed by "fixed-A" distributions   (m = 1 case)
     IF (m==0) nmax = maxindexT
     IF (m==1) nmax = maxindexA
     DO n = 0, nmax
        ! classical case: D = kB*T*(K^{-1})
        ! quantum case: d(k) = 0.5d0/(DTANH(w(k)/(2*T))*w(k))
        ! determine G := D^{1/2}; generate Sobol points corresponding to 
        ! normal distribution N(q,D) at temperature T
        G(:,:) = 0.0d0
        DO k = 1, dim-dimnull
           IF(w(k) < freq_cutoff) THEN
              d(k) = 0.0d0
           ELSE
              IF (T(n) < -0.0d-8) THEN  ! negative values of temperature used solely 
                                        ! to indicate use of the "fixed-A" distribution
                                        ! within the same loop used to sample the 
                                        ! harmonic distribution.
                 !d(k) = 1/w(k)                                ! classical
                 d(k) = 0.5d0/(DTANH(0.5d0/A(n))*w(k))         ! quantum
                 IF (verbose .AND. k==1) WRITE(*,'(A,F11.4)') 'A = ', A(n) 
                 IF (verbose) WRITE(*,'(A,I0,A,F11.4)') 'd(w[', k, '],T) = ', d(k)
              ELSE IF(T(n) > 1.0d-8) THEN
                 !d(k) = T/w(k)**2                             ! classical
                 d(k) = 0.5d0/(DTANH(w(k)/(2*T(n)))*w(k))      ! quantum
                 IF (verbose .AND. k==1) WRITE(*,'(A,F11.4,A)') 'T = ', T(n)/Ktoau, ' kelvin'
                 IF (verbose) WRITE(*,'(A,I0,A,F11.4)') 'd(w[', k, '],T) = ', d(k)
              ELSE
                 !d(k) = T/w(k)**2                             ! classical
                 d(k) = 0.5d0/w(k)                             ! quantum
                 IF (verbose .AND. k==1 .AND. m==0) WRITE(*,*) 'T =      0.0000 kelvin '
                 IF (verbose .AND. k==1 .AND. m==1) WRITE(*,*) 'A =      0.0000 '
                 IF (verbose) WRITE(*,'(A,I0,A,F11.4)') 'd(w[', k, '],T) = ', d(k)
              ENDIF
              DO i = 1, dim
                 DO j = 1, dim
                    G(i,j) = G(i,j) + DSQRT(d(k))*U(i,k)*U(j,k) ! D^{1/2} of the mass-scaled D, atomic units
                 ENDDO
              ENDDO
           ENDIF
        ENDDO  ! k = 1, dim-dimnull
        ! generate configurations according to the normal distribution N(q,D)
        r(:)  = 0.0d0
        rn(:) = 0.0d0
        IF (sobolchar == 'P' .OR. sobolchar == 'p' ) THEN
           CALL random_stdnormal(dim, lseed, rn)
        ELSE IF (sobolchar == 'Q' .OR. sobolchar == 'q') THEN
           CALL sobol_stdnormal(dim, skip, rn)
        ELSE
           WRITE(*,*) 'Invalid sequence option: choose "P" for pseudorandom or "Q" for quasirandom.'
           STOP
        END IF
        DO i = 1,dim
           r(i) = dot_PRODUCT(G(:,i),rn(:))
        ENDDO
        DO i = 1,dim
           r(i)  = r(i)/sqrt_mass(i) ! unscale
           rn(i) = r(i)+q(i)         ! atomic units
        ENDDO
        ! output the configuration, xyz file format:
        WRITE(89,'(I2)') dim/3 ! natoms
        WRITE(89,'(I9)') m*(maxindexT+1)+n+1 ! print index in place of energy
        DO j=1,dim/3
           WRITE(89,'(A2, 3F13.8)') atom_type(j), rn(j*3-2)*bohr, rn(j*3-1)*bohr, rn(j*3)*bohr ! angstrom
        ENDDO

     END DO  ! n = 0, nmax
     T(:) = -1.0d0*T(:) ! Move on to "einstein temperature" distributions.
                        ! Negative values of temperature are introduced solely
                        ! to indicate use of the "fixed-A" distribution.
  END DO  ! m = 0,1 

  CLOSE(unit=89)
  DEALLOCATE(atom_type, sqrt_mass)
  DEALLOCATE(q, r, rn)
  DEALLOCATE(w, d, G, U, C)
  DEALLOCATE(T, A)

  CALL CPU_TIME(final_time)
  WRITE(*,*) 'nconfigs: ', nconfigs
  WRITE(*,*) 'CPU TIME (seconds): ', final_time - initial_time


END PROGRAM generate_configs_normdistrbn
