program omp_test

    use iterative_utils
    use jacobi
    use cg
    use mkl_wrapper

    implicit none
    include "git_revision.fi"

    ! Model parameters
    real(kind=DP), parameter :: lx = 1, ly = 1

    integer, parameter :: use_jacobi = 1, use_gs = 2, use_cg = 3

    ! Local variables
    type(sparse_matrix) :: A
    real(kind = DP), dimension(:), allocatable :: B, X, soln
    real(kind = DP) :: hx, hy, tolerance, start_time, elapsed_time, norm_A, &
        norm_B
    integer :: Nx, Ny, method, max_iterations, iterations, verbose,         &
        stop_cond, problem, thread_count
    logical :: show_input
    character(len=30) :: fmt

    call process_options
    hx = lx / Nx
    hy = ly / Ny
    call omp_set_num_threads(thread_count)

    select case (problem)
        case (1)
            call poisson_five_pt(Nx, Ny, hx, hy, prob1_f1, zero, zero, A, B)
            call exact_solution(Nx, Ny, hx, hy, prob1_u, soln)
        case (2)
            call poisson_five_pt(Nx, Ny, hx, hy, prob2_f1, zero, prob2_u,   &
                A, B)
            call exact_solution(Nx, Ny, hx, hy, prob2_u, soln)
        case default
            call fatal("Unknown problem specified")
    end select
    allocate(X(A%cols))
    X = 0
    if (show_input) then
        call dump_input(Nx, Ny, A, B)
    end if
    write(fmt, "(A, I8, A)") "(i8, f24.15,", size(X), "f20.15)"
    norm_A = sparse_norm(A)
    norm_B = norm(B)
    start_time = dsecnd()
    select case (method)
        case (use_jacobi)
            call sparse_jacobi(A, B, X, callback, iterations)
        case (use_gs)
            call sparse_gs(A, B, X, callback, iterations)
        case (use_cg)
            call sparse_cg(A, B, X, callback, iterations)
        case default
            call fatal("Unknown method specified")
    end select
    elapsed_time = dsecnd() - start_time
    if (verbose > 1) then
        print *, "Approximation ="
        call print_solution(Nx, Ny, X)
        print *, "Error ="
        call print_solution(Nx, Ny, X - soln)
    end if
    print "(3i8,f20.15,f15.3)", Nx, Ny, iterations, norm(X - soln), elapsed_time

contains

    real(kind=DP) function prob1_f1(x, y)
        real(kind=DP), intent(in) :: x, y
        prob1_f1 = 32 * (y**2 - y + x**2 - x)
    end function prob1_f1

    real(kind=DP) function prob1_u(x, y)
        real(kind=DP), intent(in) :: x, y
        prob1_u = 16 * x * (1 - x) * y * (1 - y)
    end function prob1_u

    real(kind=DP) function prob2_f1(x, y)
        real(kind=DP), intent(in) :: x, y
        prob2_f1 = -2 * pi**2 * sin(pi * x) * cos(pi * y)
    end function prob2_f1

    real(kind=DP) function prob2_u(x, y)
        real(kind=DP), intent(in) :: x, y
        prob2_u = sin(pi * x) * cos(pi * y)
    end function prob2_u

    ! Utility function of two variables that always returns 0, for use as a
    ! placeholder when calling gen_five_pt_diag. The multiplication by x and
    ! y is to suppress "unused variable" compiler warnings.
    real(kind=DP) function zero(x, y)
        real(kind=DP), intent(in) :: x, y
        zero = 0 * x * y
    end function zero

    ! Called for each iteration of the chosen solver, where "iteration 0"
    ! is the unmodified initial estimate. This function handles per-iteration
    ! printing if the user has requested verbose output, and implements the
    ! stopping test norm(residual) / norm(RHS) < tolerance. It returns true
    ! if the solver should continue or false if it should stop.
    logical function callback(iteration, estimate)
        integer, intent(in) :: iteration
        real(kind=DP), dimension(:), intent(in) :: estimate

        real(kind=DP) :: norm_residual
        real(kind=DP), save :: init_residual, prev_residual

        norm_residual = sparse_residual(A, B, estimate)
        if (iteration == 0) then
            init_residual = norm_residual
            prev_residual = 0
        end if
        select case (verbose)
            case (1, 2)
                print fmt, iteration, norm_residual
            case (3)
                print fmt, iteration, norm_residual, estimate
            case default
                ! Do nothing
        end select
        if (iteration >= max_iterations) then
            callback = .false.
        else
            select case (stop_cond)
                case (1)
                    callback = (abs(norm_residual - prev_residual) >=       &
                        tolerance)
                    prev_residual = norm_residual
                case (2)
                    callback = (norm_residual / init_residual >= tolerance)
                case (3)
                    callback = (norm_residual / norm_B >= tolerance)
                case (4)
                    callback = (norm_residual / (norm_A * norm(estimate) +  &
                        norm_B) >= tolerance)
                case (5)
                    callback = (norm(estimate - soln) >= tolerance)
                case default
                    call fatal("Invalid stopping condition")
            end select
        end if
    end function callback

    subroutine exact_solution(Nx, Ny, hx, hy, u, A)
        integer, intent(in) :: Nx, Ny
        real(kind=DP), intent(in) :: hx, hy
        interface
            real(kind=DP) function u(x, y)
                import :: DP
                real(kind=DP), intent(in) :: x, y
            end function u
        end interface
        real(kind=DP), dimension(:), allocatable, intent(out) :: A

        integer :: i, j

        allocate(A((Nx - 1) * (Ny - 1)))
        do i = 1, Nx - 1
            do j = 1, Ny - 1
                A(i + (j - 1) * (Nx - 1)) = u(i * hx, j * hy)
            end do
        end do
    end subroutine exact_solution

    subroutine poisson_five_pt(Nx, Ny, hx, hy, f1, f2, g, A, B)
        integer, intent(in) :: Nx, Ny
        real(kind=DP), intent(in) :: hx, hy
        interface
            real(kind=DP) function f1(x, y)
                import :: DP
                real(kind=DP), intent(in) :: x, y
            end function f1
            real(kind=DP) function f2(x, y)
                import :: DP
                real(kind=DP), intent(in) :: x, y
            end function f2
            real(kind=DP) function g(x, y)
                import :: DP
                real(kind=DP), intent(in) :: x, y
            end function g
        end interface
        type(sparse_matrix), intent(out) :: A
        real(kind=DP), dimension(:), allocatable, intent(out) :: B

        real(kind = DP) :: x, y, lx, ly, rhs
        integer N, NNZ, i, j, k, row

        lx = Nx * hx
        ly = Ny * hy
        N = (Nx - 1) * (Ny - 1)
        ! Number of non-zero elements in the discretisation matrix:
        ! 5 for each interior point not adjacent to the boundary
        ! 4 for each interior point adjacent to the boundary on one side
        ! 3 for each interior point adjacent to the boundary on two sides
        NNZ = 5 * (Nx - 3) * (Ny - 3) + 4 * (2 * (Nx - 3 + Ny - 3)) + 3 * 4
        call init_sparse(A, N, N, NNZ)
        allocate(B(N))
        k = 1
        do j = 1, Ny - 1
            do i = 1, Nx - 1
                x = i * hx
                y = j * hy
                row = i + (j - 1) * (Nx - 1)
                A%row_offsets(row) = k
                rhs = -f1(x, y)
                if (j > 1) then
                    A%values(k) = -1 / hy**2
                    A%col_offsets(k) = row - (Nx - 1)
                    k = k + 1
                else
                    rhs = rhs + g(x, 0.0_DP) / hy**2
                end if
                if (i > 1) then
                    A%values(k) = -1 / hx**2
                    A%col_offsets(k) = row - 1
                    k = k + 1
                else
                    rhs = rhs + g(0.0_DP, y) / hx**2
                end if
                A%values(k) = 2 / hx**2 + 2 / hy**2 + f2(x, y)
                A%col_offsets(k) = row
                k = k + 1
                if (i < Nx - 1) then
                    A%values(k) = -1 / hx**2
                    A%col_offsets(k) = row + 1
                    k = k + 1
                else
                    rhs = rhs + g(lx, y) / hx**2
                end if
                if (j < Ny - 1) then
                    A%values(k) = -1 / hy**2
                    A%col_offsets(k) = row + (Nx - 1)
                    k = k + 1
                else
                    rhs = rhs + g(x, ly) / hy**2
                end if
                B(row) = rhs
            end do
        end do
    end subroutine poisson_five_pt

    subroutine process_options
        integer :: i
        character(len=30) :: arg

        thread_count = 1
        show_input = .false.
        method = use_cg
        verbose = 0
        max_iterations = 1000000
        tolerance = 1e-6_DP
        stop_cond = 3
        problem = 1
        Nx = 4
        Ny = Nx
        i = 1
        do while (i <= command_argument_count())
            call get_command_argument(i, arg)
            select case (arg)
                case ("-o")
                    thread_count = get_int_arg(i)
                    if (thread_count < 1) call usage
                case ("-i")
                    show_input = .true.
                case ("-g")
                    method = use_gs
                case ("-j")
                    method = use_jacobi
                case ("-c")
                    method = use_cg
                case ("-v")
                    verbose = get_int_arg(i)
                    if (verbose < 0 .or. verbose > 3) call usage
                case ("-n")
                    Nx = get_int_arg(i)
                    if (Nx < 3) call usage
                    Ny = Nx
                case ("-m")
                    max_iterations = get_int_arg(i)
                    if (max_iterations < 1) call usage
                case ("-t")
                    tolerance = get_real_arg(i)
                    if (tolerance < 0) call usage
                case ("-s")
                    stop_cond = get_int_arg(i)
                    if (stop_cond < 1 .or. stop_cond > 5) call usage
                case ("-p")
                    problem = get_int_arg(i)
                    if (problem < 1 .or. problem > 2) call usage
                case ("--version")
                    print *, 'omp_test (git ', git_revision, ')'
                    stop
                case default
                    call usage
            end select
            i = i + 1
        end do
    end subroutine process_options

    integer function get_int_arg(arg_num)
        integer, intent(in out) :: arg_num

        character(len=30) :: arg

        arg_num = arg_num + 1
        if (arg_num > command_argument_count()) call usage
        call get_command_argument(arg_num, arg)
        read (arg, *) get_int_arg
    end function get_int_arg

    real(kind=DP) function get_real_arg(arg_num)
        integer, intent(in out) :: arg_num

        character(len=30) :: arg

        arg_num = arg_num + 1
        if (arg_num > command_argument_count()) call usage
        call get_command_argument(arg_num, arg)
        read (arg, *) get_real_arg
    end function get_real_arg

    subroutine usage
        print *, "Usage: omp_test [-o threads] [-i] [-j | -g | -c] [-v] [-n mesh_size]"
        print *, "                [-m max_iterations] [-s stop_cond] " //   &
            "[-t tolerance]"
        print *, "  -o threads         Use OpenMP with given number of threads"
        print *, "  -i                 Show input to solver"
        print *, "  -j                 Use Jacobi solver"
        print *, "  -g                 Use Gauss-Seidel solver"
        print *, "  -c                 Use Conjugate Gradient solver"
        print *, "  -v level           Verbosity (0-3)"
        print *, "  -n mesh_size       Set mesh size (Nx, Ny)"
        print *, "  -m max_iterations  Set maximum iterations"
        print *, "  -t tolerance       Tolerance for stopping condition"
        print *, "  -s stop_cond       Stopping condition:"
        print *, "                       1 - Difference between iterates"
        print *, "                       2 - Residual vs initial residual"
        print *, "                       3 - Residual vs norm of RHS"
        print *, "                       4 - Residual vs norms of system"
        print *, "                       5 - Error from true solution"
        print *, "  -p problem         Model problem (1-2)"
        print *, "  --version          Print version and exit"
        stop
    end subroutine usage

    subroutine dump_input(Nx, Ny, A, B)
        integer, intent(in) :: Nx, Ny
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B

        integer :: N
        character(len=30) :: realfmt

        N = (Nx - 1) * (Ny - 1)
        if (A%rows /= N .or. A%cols /= N) then
            call fatal("dump_input: mesh dimensions do not match size " //  &
                "of input A")
        end if
        if (size(B) /= N) then
            call fatal("dump_input: mesh dimensions do not match size " //  &
                "of input B")
        end if
        write(realfmt, "(A, I8, A)") "(", N, "f16.6)"
        print *, N
        call dump_sparse(A, "f16.6")
        call print_sparse(A, "f16.6")
        print realfmt, B
    end subroutine dump_input

    subroutine print_solution(Nx, Ny, B)
        integer, intent(in) :: Nx, Ny
        real(kind=DP), dimension(:), intent(in) :: B

        integer :: i
        character(len=30) :: realfmt

        write(realfmt, "(A, I8, A)") "(", Nx - 1, "f15.10)"
        do i = 1, Ny - 1
            print realfmt, B((i - 1) * (Nx - 1) + 1 : i * (Nx - 1))
        end do
    end subroutine print_solution

end program omp_test
