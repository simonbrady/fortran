! Solve 2-D parabolic systems with dense solution matrix and zero Dirichlet
! boundary conditions.
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 20-May-2014  SJB  Initial version
! 22-May-2014  SJB  Implement ADI
! 24-May-2014  SJB  Use init_tridiag_const in do_adi
! 27-May-2014  SJB  Add option to include boundary in output
! 21-Mar-2015  SJB  Add --version option

program dense_parabolic

    ! Get real kind DP (portable double precision) from module supplied with
    ! Intel Math Kernel Library
    use f95_precision

    ! Local module for tridiagonal matrix type
    use tridiag

    implicit none
    include "git_revision.fi"

    integer, parameter :: use_explicit = 1, use_adi = 2

    ! Numbers smaller than this are considered to be zero
    real(kind=DP), parameter :: eps = 1e-15_DP

    integer :: Nx, Ny, Nt, method, old, new
    real(kind=DP) :: dx, dy, dt, mu_x, mu_y
    logical :: print_boundary, verbose
    real(kind=DP), dimension(:, :, :), allocatable :: U
    real(kind=DP), dimension(:), allocatable :: zeros
    character(len=30) :: row_fmt

    ! Set options from command line
    call process_options
    dx = 1.0_DP / Nx
    dy = 1.0_DP / Ny
    mu_x = dt / dx**2
    mu_y = dt / dy**2
    allocate(U(Nx - 1, Ny - 1, 2), zeros(Nx + 1))
    zeros = 0.0
    write(row_fmt, "(A, I8, A)") "(", Nx + 1, "f10.5)"

    ! Set initial condition
    old = 1
    new = 2
    call init_cond(Nx, Ny, dx, dy, U(:, :, old))
    call print_soln(U(:, :, old), 0)

    ! Call appropriate solver
    select case (method)
        case (use_explicit)
            call do_explicit
        case (use_adi)
            call do_adi
        case default
            call fatal("Unknown method specified")
    end select

contains

    ! Set matrix U to initial condition across mesh, excluding boundary values
    ! which are always assumed to be zero
    subroutine init_cond(Nx, Ny, dx, dy, U)
        integer, intent(in) :: Nx, Ny
        real(kind=DP), intent(in) :: dx, dy
        real(kind=DP), dimension(:, :), intent(out) :: U

        integer :: r, s
        real(kind=DP) :: x, y

        do s = 1, Ny - 1
            do r = 1, Nx - 1
                x = r * dx
                y = s * dy
                U(r, s) = (1 - 4 * (x - 0.5)**2) * (1 - 4 * (y - 0.5)**2)
            end do
        end do
    end subroutine init_cond

    ! 2-D explicit scheme for zero Dirichlet boundary condition (Morton and
    ! Mayers, pp62-64)
    subroutine do_explicit
        integer :: n, r, s
        real(kind=DP) :: rhs

        do n = 1, Nt
            do s = 1, Ny - 1
                do r = 1, Nx - 1
                    rhs = (1 - 2 * mu_x - 2 * mu_y) * U(r, s, old)
                    if (r > 1) rhs = rhs + mu_x * U(r - 1, s, old)
                    if (r < Nx - 1) rhs = rhs + mu_x * U(r + 1, s, old)
                    if (s > 1) rhs = rhs + mu_x * U(r, s - 1, old)
                    if (s < Ny - 1) rhs = rhs + mu_x * U(r, s + 1, old)
                    U(r, s, new) = rhs
                end do
            end do
            call print_soln(U(:, :, new), n)
            ! Swap old and new indices for next iteration
            old = 3 - old
            new = 3 - new
        end do
    end subroutine do_explicit

    ! ADI scheme for zero Dirichlet boundary condition (Morton and Mayers,
    ! pp64-66)
    subroutine do_adi
        type(tridiag_matrix) :: A1, A2, B1, B2
        integer :: n, r, s

        ! Use tridiagonal matrices to solve ADI equations in matrix-vector form:
        ! (3.15a) A1 * U^{n + 1/2} = B1 * U^n
        ! (3.15b) A2 * U^{n + 1} = B2 * U^{n + 1/2}
        call init_tridiag_const(A1, Nx - 1, -0.5 * mu_x, 1 + mu_x, -0.5 * mu_x)
        call init_tridiag_const(A2, Ny - 1, -0.5 * mu_y, 1 + mu_y, -0.5 * mu_y)
        call init_tridiag_const(B1, Ny - 1, 0.5 * mu_y, 1 - mu_y, 0.5 * mu_y)
        call init_tridiag_const(B2, Nx - 1, 0.5 * mu_x, 1 - mu_x, 0.5 * mu_x)
        do n = 1, Nt
            ! Explicit half-step for each column (RHS of eqn 3.15a)
            do r = 1, Nx - 1
                call tridiag_mult(B1, U(r, :, old), U(r, :, new))
            end do
            ! Implicit half-step for each row (LHS of eqn 3.15a)
            do s = 1, Ny - 1
                call tridiag_solve(A1, U(:, s, old), U(:, s, new))
            end do
            ! Explicit half-step for each row (RHS of eqn 3.15b)
            do s = 1, Ny - 1
                call tridiag_mult(B2, U(:, s, old), U(:, s, new))
            end do
            ! Implicit half-step for each column (LHS of eqn 3.15b)
            do r = 1, Nx - 1
                call tridiag_solve(A2, U(r, :, old), U(r, :, new))
            end do
            call print_soln(U(:, :, old), n)
        end do
    end subroutine do_adi

    ! Print approximate solution in matrix form, optionally including boundary
    ! values
    subroutine print_soln(U, n)
        real(kind=DP), dimension(:, :), intent(in) :: U
        integer, intent(in) :: n

        integer :: s

        if (verbose) then
            print *, ''
            print '(A, I6)', 'n =', n
            print *, ''
        end if
        if (verbose .or. n == Nt) then
            if (print_boundary) then
                print row_fmt, zeros
                do s = 1, Ny - 1
                    print row_fmt, 0.0, U(:, s), 0.0
                end do
                print row_fmt, zeros
            else
                do s = 1, Ny - 1
                    print row_fmt, U(:, s)
                end do
            end if
        end if
    end subroutine

    ! Set default program options then process command-line arguments to
    ! modify them as directed by the user
    subroutine process_options
        integer :: i
        character(len=30) :: arg

        method = use_explicit
        Nt = 15
        Nx = 10
        Ny = 10
        dt = 0.001_DP
        print_boundary = .false.
        verbose = .false.
        i = 1
        do while (i <= command_argument_count())
            call get_command_argument(i, arg)
            select case (arg)
                case ("-a")
                    method = use_adi
                case ("-e")
                    method = use_explicit
                case ("-n")
                    Nt = get_int_arg(i)
                    if (Nt < 0) call usage
                case ("-x")
                    Nx = get_int_arg(i)
                    if (Nx < 1) call usage
                case ("-y")
                    Ny = get_int_arg(i)
                    if (Ny < 1) call usage
                case ("-t")
                    dt = get_real_arg(i)
                    if (dt < eps) call usage
                case ("-b")
                    print_boundary = .true.
                case ("-v")
                    verbose = .true.
                case ("--version")
                    print *, 'dense_parabolic (git ', git_revision, ')'
                    print *, 'built with ', trim(tridiag_version())
                    stop
                case default
                    call usage
            end select
            i = i + 1
        end do
    end subroutine process_options

    ! Read the next command-line argument as an integer-valued parameter
    integer function get_int_arg(arg_num)
        integer, intent(in out) :: arg_num

        character(len=30) :: arg

        arg_num = arg_num + 1
        if (arg_num > command_argument_count()) call usage
        call get_command_argument(arg_num, arg)
        read (arg, *) get_int_arg
    end function get_int_arg

    ! Read the next command-line argument as a real-valued parameter
    real(kind=DP) function get_real_arg(arg_num)
        integer, intent(in out) :: arg_num

        character(len=30) :: arg

        arg_num = arg_num + 1
        if (arg_num > command_argument_count()) call usage
        call get_command_argument(arg_num, arg)
        read (arg, *) get_real_arg
    end function get_real_arg

    ! Print a help message on how to invoke the program, then quit
    subroutine usage
        print *, "Usage: dense_parabolic [-a | -e] [-n time_steps] " //     &
            "[-x x_steps] [-y y_steps]"
        print *, "              [-t time_step] [-b] [-v]"
        print *, "  -a                 Use ADI method"
        print *, "  -e                 Use explicit method"
        print *, "  -n time_steps      Number of steps in time"
        print *, "  -x x_steps         Number of steps in x-direction"
        print *, "  -y y_steps         Number of steps in y-direction"
        print *, "  -t time_step       Length of each time step"
        print *, "  -b                 Include boundary in output"
        print *, "  -v                 Verbose output"
        print *, "  --version          Print version and exit"
        stop
    end subroutine usage

end program dense_parabolic
