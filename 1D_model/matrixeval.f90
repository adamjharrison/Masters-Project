!gfortran matrixeval.f90 -o matrixeval to compile
!./matrixeval 10.0 output10.csv to run for rsig=10.0
program matrixeval !Evaluation of matrix elements in Equation 24 of Granucci, J. Chem. Phys., 2012, 137, 22A501
use, intrinsic :: iso_fortran_env, only: REAL64, INT32

implicit none

!Define double and integers
integer, parameter :: dp = REAL64
integer, parameter :: DefInt = INT32

!Parameters for PESs
real(kind=dp), parameter :: alpha_1=0.35_dp
real(kind=dp), parameter :: alpha_2=0.25_dp
real(kind=dp), parameter :: a_2=0.5_dp
real(kind=dp), parameter :: dE=0.04_dp

!Parameters for SOC sigmoid function
real(kind=dp), parameter :: r_c=10.0_dp
real(kind=dp), parameter :: d_sig=2.0_dp
real(kind=dp), parameter :: c0=0.001_dp

real(kind=dp), parameter :: r_start = 1.0_dp
real(kind=dp) :: i = 1.0_dp
complex(kind=dp), parameter :: c1 = cmplx(0.0005_dp, 0.0005_dp, kind=dp) !
complex(kind=dp), parameter :: im = cmplx(0.0_dp, 1.0_dp, kind=dp) !Define imaginary complex

real(kind=dp), parameter :: n_points = 20.0_dp !Number of points
integer (kind=DefInt), parameter :: n = 4
real(kind=dp) :: E_s, E_t, b, r, sig, a_1, r_sig
complex(kind=dp) :: z

complex(kind=dp), DIMENSION(4, 4) :: V4 !4x4 Matrix
integer(kind=DefInt)     ::  row, col
integer :: file_unit

!Get command line arguments for given model variant
character(len=256) :: arg_rsig, arg_file
call get_command_argument(1, arg_rsig)
call get_command_argument(2, arg_file)
read(arg_rsig, *) r_sig

file_unit = 10
OPEN(UNIT=file_unit, FILE= trim(arg_file), STATUS='replace', ACTION='write')
write(file_unit,'(A)') "bohrs, E_s, E_t, SOC_Real, SOC_Imag" !Write column names to file

a_1 = (a_2*exp(-alpha_2*r_c)-dE)/(exp(-alpha_1*r_c)) !Determine a_1 so that S_1 crosses T_1 at r_c=10 bohr

DO WHILE(r .lt. n_points) !Evaluate the matrix elements over a one dimensional grid
    r = r_start + real(i-1,dp)

    !Calculate SOC between S_1 and T_1 at given position
    sig = calc_sigma(r,r_sig,d_sig)
    b = c0*sig
    z = c1*sig

    !Calculate energy of states at given position
    E_s = a_1*exp(-alpha_1*r)+dE
    E_t = a_2*exp(-alpha_2*r)


    V4 = cmplx(0.0_dp,0.0_dp,kind=dp)

    !Add energies to the diagonal columns of the matrix
    V4(1,1) = cmplx(E_s,0.0_dp,kind=dp)
    V4(2,2) = cmplx(E_t,0.0_dp,kind=dp)
    V4(3,3) = cmplx(E_t,0.0_dp,kind=dp)
    V4(4,4) = cmplx(E_t,0.0_dp,kind=dp)
    
    !Off-diagonal SOC terms between singlet and the triplet sublevels
    V4(1,2) = z
    V4(4,1) = z
    V4(2,1) = conjg(z)
    V4(1,4) = conjg(z)
    V4(1,3) = im*b
    V4(3,1) = -im*b

    !Write matrix elements to file
    write(file_unit,'(G0, 4( ", ", G0 ))') r,real(V4(1, 1)), real(V4(2, 2)), real(V4(1, 2)), aimag(V4(1, 3))
    i = i + 0.1_dp
END DO
CLOSE(file_unit)
contains

    function calc_sigma(r, rs, drs) result(sigma) !Determine the SOC coupling over 
        real(kind=dp), intent(in) :: r, rs, drs
        real(kind=dp) :: sigma, term, rmd, rpd
        term = (r-rs)/drs
        rmd = rs-(drs/2.0_dp)
        rpd = rs+(drs/2.0_dp)

        !Sigmoid function conditions
        if (r.le.rmd) then
            sigma = 1.0_dp
        else if (r.gt.rmd .and. r.lt.rpd) then
            sigma = 4*((term)**3)-3*term
        else
            sigma= -1.0_dp
        end if
    end function calc_sigma

end program matrixeval
