! Iterative linear solver for dense matrices using Jacobi and CG modules. See
! the subroutine usage() below for a list of valid command-line options.
!
! Reads one or more input groups of a matrix A and vectors b, y formatted as
!
!   N
!   a11 a12 ... a1N
!     ...
!   aN1 aN2 ... aNN
!   b11 b21 ... bN1
!   x11 x12 ... xN1
!
! and prints the solution to the linear system Ax = b using input vector x
! as an initial estimate. To terminate the program enter N = 0.
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 14-Apr-2014  SJB  Initial version
! 15-Apr-2014  SJB  Use iterative_utils
! 19-Apr-2014  SJB  Use callback, implement CG / SD
! 21-Mar-2015  SJB  Add --version option

program dense_solver

    ! Local libraries for Jacobi/GS/SOR and CG iterative methods
    use iterative_utils
    use jacobi
    use cg

    implicit none
    include "git_revision.fi"

    integer, parameter :: use_jacobi = 1, use_gs = 2, use_sor = 3,          &
        use_cg = 4, use_sd = 5

    integer :: N, i, method, max_iterations, iterations
    real(kind=DP) :: omega, tolerance, norm_B
    logical :: verbose
    real(kind=DP), dimension(:, :), allocatable :: A
    real(kind=DP), dimension(:), allocatable :: B, X
    character(len=30) :: fmt

    ! Set options from command line
    call process_options
    do
        read *, N
        if (N < 1) exit

        ! Allocate storage and read input group
        allocate(A(N, N), B(N), X(N))
        do i = 1, N
            read *, A(i, :)
        end do
        read *, B
        read *, X

        ! Compute norm of RHS for use in stopping test
        norm_B = norm2(B)

        ! Set output format for verbose output
        write(fmt, "(A, I8, A)") "(i8,", N + 1, "f12.6)"

        ! Call appropriate solver
        select case (method)
            case (use_jacobi)
                call dense_jacobi(A, B, X, callback, iterations)
            case (use_gs)
                call dense_gs(A, B, X, callback, iterations)
            case (use_sor)
                call dense_sor(A, B, X, omega, callback, iterations)
            case (use_cg)
                call dense_cg(A, B, X, callback, iterations)
            case (use_sd)
                call dense_sd(A, B, X, callback, iterations)
            case default
                call fatal("Unknown method specified")
        end select

        ! Print iteration count and result
        print *, "iterations =", iterations
        print *, "solution ="
        do i = 1, N
            print "(f20.10)", X(i)
        end do

        ! Clean up for next input group
        deallocate(A, B, X)
    end do

contains

    ! Called for each iteration of the chosen solver, where "iteration 0"
    ! is the unmodified initial estimate. This function handles per-iteration
    ! printing if the user has requested verbose output, and implements the
    ! stopping test norm(residual) / norm(RHS) < tolerance. It returns true
    ! if the solver should continue or false if it should stop.
    logical function callback(iteration, estimate)
        integer, intent(in) :: iteration
        real(kind=DP), dimension(:), intent(in) :: estimate

        real(kind=DP) :: norm_residual

        norm_residual = residual2(A, B, estimate)
        if (verbose) print fmt, iteration, norm_residual, estimate
        if (iteration >= max_iterations) then
            callback = .false.
        else
            callback = (norm_residual / norm_B >= tolerance)
        end if
    end function callback

    ! Set default program options then process command-line arguments to
    ! modify them as directed by the user
    subroutine process_options
        integer :: i
        character(len=30) :: arg

        method = use_gs
        omega = 1.0_DP
        max_iterations = 1000
        tolerance = 1e-6_DP
        verbose = .false.
        i = 1
        do while (i <= command_argument_count())
            call get_command_argument(i, arg)
            select case (arg)
                case ("-j")
                    method = use_jacobi
                case ("-g")
                    method = use_gs
                case ("-s")
                    method = use_sor
                    i = i + 1
                    if (i > command_argument_count()) call usage
                    call get_command_argument(i, arg)
                    read (arg, *) omega
                case ("-c")
                    method = use_cg
                case ("-d")
                    method = use_sd
                case ("-m")
                    i = i + 1
                    if (i > command_argument_count()) call usage
                    call get_command_argument(i, arg)
                    read (arg, *) max_iterations
                case ("-t")
                    i = i + 1
                    if (i > command_argument_count()) call usage
                    call get_command_argument(i, arg)
                    read (arg, *) tolerance
                case ("-v")
                    verbose = .true.
                case ("--version")
                    print *, 'dense_solver (git ', git_revision, ')'
                    print *, 'built with ', trim(cg_version()), ', ',       &
                        trim(jacobi_version()), ','
                    print *, '    ', trim(iterative_utils_version())
                    stop
                case default
                    call usage
            end select
            i = i + 1
        end do
    end subroutine process_options

    ! Print a help message on how to invoke the program, then quit
    subroutine usage()
        print *, "Usage: dense_solver [-j | -g | -s omega | -c | -d] " //   &
            "[-m max_iterations]"
        print *, "              [-t tolerance] [-v]"
        print *, "  -j                 Use Jacobi iterative solver"
        print *, "  -g                 Use Gauss-Seidel iterative solver"
        print *, "  -s omega           Use SOR with given value of omega"
        print *, "  -c                 Use conjugate gradient solver"
        print *, "  -d                 Use steepest descent solver"
        print *, "  -m max_iterations  Set maximum iterations"
        print *, "  -t tolerance       Set stopping tolerance"
        print *, "  -v                 Display verbose output"
        print *, "  --version          Print version and exit"
        stop
    end subroutine usage

end program dense_solver
