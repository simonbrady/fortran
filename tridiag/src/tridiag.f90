! Tridiagonal matrix type and supporting functions
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 22-May-2014  SJB  Initial version
! 24-May-2014  SJB  Implement init_tridiag_const
! 21-Mar-2015  SJB  Add tridiag_version

module tridiag

    ! Get real kind DP (portable double precision) from module supplied with
    ! Intel Math Kernel Library
    use f95_precision

    implicit none

    ! Only export entities explicity marked as public
    private

    ! Subroutines exported from this module
    public :: init_tridiag, init_tridiag_const, clear_tridiag, dump_tridiag
    public :: tridiag_set, tridiag_get, tridiag_to_dense
    public :: tridiag_mult, tridiag_solve
    public :: tridiag_version

    ! Derived type for tridiagonal N*N matrices with each diagonal stored as a
    ! dense column vector of an N*3 matrix
    type, public :: tridiag_matrix
        ! Dimension of matrix
        integer :: N
        ! values(:, 1) holds first subdiagonal, values(:, 2) holds main
        ! diagonal, and values(:, 3) holds first superdiagonal
        real(kind=DP), dimension(:, :), allocatable :: values
    end type tridiag_matrix

contains

    ! Initialise storage for a new tridiagonal matrix of dimension N*N
    subroutine init_tridiag(A, N)
        type(tridiag_matrix), intent(out) :: A
        integer, intent(in) :: N

        A%N = N
        allocate(A%values(N, 3))
    end subroutine init_tridiag

    ! Initialise storage for a new tridiagonal matrix of dimension N*N and set
    ! diagonal elements to constant values:
    !  p - first subdiagonal
    !  q - main diagonal
    !  r - first superdiagonal
    subroutine init_tridiag_const(A, N, p, q, r)
        type(tridiag_matrix), intent(out) :: A
        integer, intent(in) :: N
        real(kind=DP), intent(in) :: p, q, r

        integer :: i

        call init_tridiag(A, N)
        do i = 1, N
            A%values(i, 1) = p
            A%values(i, 2) = q
            A%values(i, 3) = r
        end do
        ! Zero out undefined entries
        A%values(1, 1) = 0
        A%values(N, 3) = 0
    end subroutine init_tridiag_const

    ! Deallocate storage for a tridiagonal matrix so it can be reinitialised
    subroutine clear_tridiag(A)
        type(tridiag_matrix), intent(in out) :: A

        A%N = 0
        deallocate(A%values)
    end subroutine clear_tridiag

    ! Diagnostic routine to dump the internal structure of a tridiagonal matrix
    subroutine dump_tridiag(A, entry_fmt)
        type(tridiag_matrix), intent(in) :: A
        character(len=*), intent(in) :: entry_fmt

        character(len=30) :: row_fmt
        integer :: i

        write(row_fmt, "(A, I8, A)") "(", A%N, entry_fmt // ")"
        print *, A%N
        do i = 1, 3
            print row_fmt, A%values(:, i)
        end do
    end subroutine dump_tridiag

    ! Set (row, col) element of tridiagonal matrix A to val
    subroutine tridiag_set(A, row, col, val)
        type(tridiag_matrix), intent(in out) :: A
        integer, intent(in) :: row, col
        real(kind=DP), intent(in) :: val

        integer :: i

        i = col - row + 2
        if (row < 1 .or. row > A%N .or. i < 1 .or. i > 3) then
            call fatal("tridiag_set: invalid index specified")
        else
            A%values(row, i) = val
        end if
    end subroutine tridiag_set

    ! Get (row, col) element of tridiagonal matrix A
    real(kind=DP) function tridiag_get(A, row, col)
        type(tridiag_matrix), intent(in) :: A
        integer, intent(in) :: row, col

        integer :: i

        if (row < 1 .or. row > A%N .or. col < 1 .or. col > A%N) then
            call fatal("tridiag_get: invalid index specified")
        end if
        i = col - row + 2
        if (i < 1 .or. i > 3) then
            tridiag_get = 0
        else
            tridiag_get = A%values(row, i)
        end if
    end function tridiag_get

    ! Copy triagonal matrix A to a new dense matrix B of the same dimension
    subroutine tridiag_to_dense(A, B)
        type(tridiag_matrix), intent(in) :: A
        real(kind=DP), dimension(:, :), allocatable, intent(out) :: B

        integer :: i, j

        if (A%N < 1) call fatal("tridiag_to_dense: invalid matrix dimension")
        allocate(B(A%N, A%N))
        do i = 1, A%N
            do j = 1, A%N
                B(i, j) = tridiag_get(A, i, j)
            end do
        end do
    end subroutine tridiag_to_dense

    ! Compute matrix-vector product y = Ax
    subroutine tridiag_mult(A, X, Y)
        type(tridiag_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: X
        real(kind=DP), dimension(:), intent(out) :: Y

        integer :: i

        if (A%N < 2) call fatal("tridiag_mult: invalid matrix dimension")
        if (size(X) /= A%N) call fatal("tridiag_mult: length mismatch for X")
        if (size(Y) /= A%N) call fatal("tridiag_mult: length mismatch for Y")

        Y(1) = A%values(1, 2) * X(1) + A%values(1, 3) * X(2)
        do i = 2, A%N - 1
            Y(i) = A%values(i, 1) * X(i - 1) + A%values(i, 2) * X(i) +      &
                A%values(i, 3) * X(i + 1)
        end do
        Y(A%N) = A%values(A%N, 1) * X(A%N - 1) + A%values(A%N, 2) * X(A%N)
    end subroutine tridiag_mult

    ! Solve tridiagonal system Ax = b using Thomas algorithm. Based on Morton
    ! and Mayers, "Numerical Solution of PDEs" (2nd ed.), pp24-25, with the
    ! following mapping between their notation and our variables:
    !
    ! a_j : -A%values(j, 1)    b_j : A%values(j, 2)     c_j : -A%values(j, 3)
    ! d_j : B(j)               e_j : alpha(j)           f_j : beta(j)
    !
    ! Note that a and c are negated in the book's formulation of the problem.
    subroutine tridiag_solve(A, X, B)
        type(tridiag_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(out) :: X
        real(kind=DP), dimension(:), intent(in) :: B

        real(kind=DP), dimension(:), allocatable :: alpha, beta
        integer :: i

        if (A%N < 2) call fatal("tridiag_solve: invalid matrix dimension")
        if (size(X) /= A%N) call fatal("tridiag_solve: length mismatch for X")
        if (size(B) /= A%N) call fatal("tridiag_solve: length mismatch for B")
        allocate(alpha(A%N), beta(A%N))
        ! Forward pass
        alpha(1) = -A%values(1, 3) / A%values(1, 2)
        beta(1) = B(1) / A%values(1, 2)
        do i = 2, A%N
            alpha(i) = -A%values(i, 3) / (A%values(i, 2) + A%values(i, 1) * &
                alpha(i - 1))
            beta(i) = (B(i) - A%values(i, 1) * beta(i - 1)) /               &
                (A%values(i, 2) + A%values(i, 1) * alpha(i - 1))
        end do
        ! Backward pass
        X(A%N) = beta(A%N)
        do i = A%N - 1, 1, -1
            X(i) = beta(i) + alpha(i) * X(i + 1)
        end do
    end subroutine tridiag_solve

    ! Library version string
    character(len=50) function tridiag_version()
        include "git_revision.fi"
        tridiag_version = 'tridiag (git ' // git_revision // ')'
    end function tridiag_version

end module tridiag
