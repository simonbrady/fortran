! Jacobi, Gauss-Seidel and SOR iterative solvers
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 14-Apr-2014  SJB  Initial version
! 15-Apr-2014  SJB  Use iterative_utils
! 19-Apr-2014  SJB  Use callback
! 21-Mar-2015  SJB  Add jacobi_version

module jacobi

    ! Local library of common utility functions
    use iterative_utils

    implicit none

    ! Only export entities explicity marked as public
    private

    ! Subroutines exported from this module
    public :: dense_jacobi, dense_gs, dense_sor
    public :: sparse_jacobi, sparse_gs
    public :: jacobi_version

contains

    ! Jacobi method for dense real matrices
    !
    ! Input parameters:
    !
    ! A - N-by-N matrix for system Ax = b
    !
    ! B - N-vector of RHS
    !
    ! X - N-vector containing initial guess on input and overwritten with
    !     computed solution on output
    !
    ! callback - caller-supplied function to implement stopping test and
    !            return true/false depending on whether the solver should
    !            keep iterating or not. Called with the current iteration
    !            number and estimate (an N-vector), where iteration 0 is
    !            the initial estimate before any steps have been taken.
    !
    ! Output parameters:
    !
    ! iterations - set to the number of iterations taken

    subroutine dense_jacobi(A, B, X, callback, iterations)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B
        real(kind=DP), dimension(:), intent(in out) :: X
        interface
            logical function callback(iteration, estimate)
                import :: DP
                integer, intent(in) :: iteration
                real(kind=DP), dimension(:), intent(in) :: estimate
            end function callback
        end interface
        integer, intent(out) :: iterations

        ! U contains the current and previous estimates, indexed by variables
        ! old and new. These indices swap back and forth as we proceed so only
        ! the two most recent estimates are needed.
        real(kind=DP), dimension(:, :), allocatable :: U
        integer :: N, i, k, old, new

        ! Validate sizes of input matrices
        N = get_dim(A, B, X)

        ! Set working storage for iterations and initial estimate
        allocate(U(N, 2))
        new = 1
        old = 2
        U(:, new) = X
        k = 0

        ! Main iteration loop
        do while (callback(k, U(:, new)))
            k = k + 1
            old = 3 - old           ! Swap between 1 and 2
            new = 3 - new
            do i = 1, N
                U(i, new) = (B(i) -                                         &
                    dot_product(A(i, 1 : i - 1), U(1 : i - 1, old)) -       &
                    dot_product(A(i, i + 1 : N), U(i + 1 : N, old))) / A(I, I)
            end do
        end do

        ! Set output parameters
        X = U(:, new)
        iterations = k
    end subroutine dense_jacobi

    ! Gauss-Seidel method for dense real matrices
    !
    ! Input/output parameters same as jacobi

    subroutine dense_gs(A, B, X, callback, iterations)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B
        real(kind=DP), dimension(:), intent(in out) :: X
        interface
            logical function callback(iteration, estimate)
                import :: DP
                integer, intent(in) :: iteration
                real(kind=DP), dimension(:), intent(in) :: estimate
            end function callback
        end interface
        integer, intent(out) :: iterations

        ! Gauss-Seidel doesn't require any working storage for estimates since
        ! we can update X in place
        integer :: N, i, k

        N = get_dim(A, B, X)
        k = 0
        do while (callback(k, X))
            k = k + 1
            do i = 1, N
                X(i) = (B(i) - dot_product(A(i, 1 : i - 1), X(1 : i - 1)) - &
                    dot_product(A(i, i + 1 : N), X(i + 1 : N))) / A(I, I)
            end do
        end do

        ! Set output parameter
        iterations = k
    end subroutine dense_gs

    ! SOR for dense real matrices
    !
    ! Input/output parameters same as dense_jacobi, with the addition of
    ! input omega which is the SOR weighting parameter

    subroutine dense_sor(A, B, X, omega, callback, iterations)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B
        real(kind=DP), dimension(:), intent(in out) :: X
        real(kind=DP), intent(in) :: omega
        interface
            logical function callback(iteration, estimate)
                import :: DP
                integer, intent(in) :: iteration
                real(kind=DP), dimension(:), intent(in) :: estimate
            end function callback
        end interface
        integer, intent(out) :: iterations

        ! SOR performs a weighted combination of the previous estimate and a
        ! Gauss-Seidel update, so we require three vectors for working storage.
        ! Like jacobi the old and new variables swap back and forth between 1
        ! and 2, while row 3 is dedicated to storing the intermediate estimate.
        integer, parameter :: temp = 3
        real(kind=DP), dimension(:, :), allocatable :: U
        integer :: N, i, k, old, new

        ! Validate sizes of input matrices
        N = get_dim(A, B, X)

        ! Set working storage for iterations and initial estimate
        allocate(U(N, 3))
        new = 1
        old = 2
        U(:, new) = X
        k = 0

        ! Main iteration loop
        do while (callback(k, U(:, new)))
            k = k + 1
            old = 3 - old
            new = 3 - new
            do i = 1, N
                U(i, temp) = (B(i) -                                        &
                    dot_product(A(i, 1 : i - 1), U(1 : i - 1, temp)) -      &
                    dot_product(A(i, i + 1 : N), U(i + 1 : N, old))) / A(I, I)
            end do
            U(:, new) = (1 - omega) * U(:, old) + omega * U(:, temp)
        end do

        ! Set output parameters
        X = U(:, new)
        iterations = k
    end subroutine dense_sor

    ! Jacobi method for sparse real matrices
    !
    ! Input/output parameters same as dense_jacobi, except that A is an N*N
    ! sparse matrix using the sparse_matrix type defined in interative_utils.

    subroutine sparse_jacobi(A, B, X, callback, iterations)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B
        real(kind=DP), dimension(:), intent(in out) :: X
        interface
            logical function callback(iteration, estimate)
                import :: DP
                integer, intent(in) :: iteration
                real(kind=DP), dimension(:), intent(in) :: estimate
            end function callback
        end interface
        integer, intent(out) :: iterations

        ! U contains the current and previous estimates, indexed by variables
        ! old and new. These indices swap back and forth as we proceed so only
        ! the two most recent estimates are needed.
        real(kind=DP), dimension(:, :), allocatable :: U
        integer :: N, i, k, old, new

        ! Validate sizes of input matrices
        N = get_sparse_dim(A, B, X)

        ! Set working storage for iterations and initial estimate
        allocate(U(N, 2))
        new = 1
        old = 2
        U(:, new) = X
        k = 0

        ! Main iteration loop
        do while (callback(k, U(:, new)))
            k = k + 1
            old = 3 - old           ! Swap between 1 and 2
            new = 3 - new
            do i = 1, N
                U(i, new) = (B(i) -                                         &
                    sparse_partial_row_mult(A, i, 1, i - 1, U(:, old), 1) - &
                    sparse_partial_row_mult(A, i, i + 1, N, U(:, old),      &
                        i + 1)) / sparse_entry(A, i, i)
            end do
        end do

        ! Set output parameters
        X = U(:, new)
        iterations = k
    end subroutine sparse_jacobi

    ! Gauss-Seidel method for sparse real matrices
    !
    ! Input/output parameters same as jacobi

    subroutine sparse_gs(A, B, X, callback, iterations)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B
        real(kind=DP), dimension(:), intent(in out) :: X
        interface
            logical function callback(iteration, estimate)
                import :: DP
                integer, intent(in) :: iteration
                real(kind=DP), dimension(:), intent(in) :: estimate
            end function callback
        end interface
        integer, intent(out) :: iterations

        ! Gauss-Seidel doesn't require any working storage for estimates since
        ! we can update X in place
        integer :: N, i, k

        N = get_sparse_dim(A, B, X)
        k = 0
        do while (callback(k, X))
            k = k + 1
            do i = 1, N
                X(i) = (B(i) -                                              &
                    sparse_partial_row_mult(A, i, 1, i - 1, X, 1) -         &
                    sparse_partial_row_mult(A, i, i + 1, N, X, i + 1)) /    &
                    sparse_entry(A, i, i)
            end do
        end do

        ! Set output parameter
        iterations = k
    end subroutine sparse_gs

    ! Library version string
    character(len=50) function jacobi_version()
        include "git_revision.fi"
        jacobi_version = 'jacobi (git ' // git_revision // ')'
    end function jacobi_version

end module jacobi
