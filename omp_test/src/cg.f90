! Conjugate gradient and steepest descent iterative solvers
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 19-Apr-2014  SJB  Initial version

module cg

    ! Local library of common utility functions
    use iterative_utils

    implicit none

    ! Only export entities explicity marked as public
    private

    ! Subroutines exported from this module
    public :: dense_cg, sparse_cg, dense_sd

contains

    ! Conjugate gradient method for dense real matrices, taken almost verbatim
    ! from Golub and Ortega, figure 9.1.2, p374
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

    subroutine dense_cg(A, B, X, callback, iterations)
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

        integer :: N, k
        real(kind=DP) :: rdotr, rdotr_new, alpha, beta, pdotAp
        real(kind=DP), dimension(:), allocatable :: R, P, AP

        N = get_dim(A, B, X)
        allocate(R(N), P(N), AP(N))

        k = 0
        R = B - matmul(A, X)
        rdotr = dot_product(R, R)
        P = R
        do while (callback(k, X))
            k = k + 1
            AP = matmul(A, P)
            pdotAp = dot_product(P, AP)
            alpha = -rdotr / pdotAP
            X = X - alpha * P
            R = R + alpha * AP
            rdotr_new = dot_product(R, R)
            beta = rdotr_new / rdotr
            rdotr = rdotr_new
            P = R + beta * P
        end do
        iterations = k
    end subroutine dense_cg

    ! Conjugate Gradient method for sparse matrices. Same parameters as
    ! dense_cg except that A is of type sparse_matrix declared in module
    ! iterative_utils.

    subroutine sparse_cg(A, B, X, callback, iterations)
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

        integer :: N, k
        real(kind=DP) :: rdotr, rdotr_new, alpha, beta, pdotAp
        real(kind=DP), dimension(:), allocatable :: R, P, AP

        N = get_sparse_dim(A, B, X)
        allocate(R(N), P(N), AP(N))

        k = 0
        call sparse_mv_mult(A, X, P)         ! Use P as temp here
        R = B - P
        rdotr = dot_product(R, R)
        P = R
        do while (callback(k, X))
            k = k + 1
            call sparse_mv_mult(A, P, AP)
            pdotAp = dot_product(P, AP)
            alpha = -rdotr / pdotAP
            X = X - alpha * P
            R = R + alpha * AP
            rdotr_new = dot_product(R, R)
            beta = rdotr_new / rdotr
            rdotr = rdotr_new
            P = R + beta * P
        end do
        iterations = k
    end subroutine sparse_cg

    ! Steepest descent method for dense real matrices, based on CG with P = R
    !
    ! Input/output parameters same as cg

    subroutine dense_sd(A, B, X, callback, iterations)
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

        integer :: N, k
        real(kind=DP) :: rdotr, rdotAr, alpha
        real(kind=DP), dimension(:), allocatable :: R, AR

        N = get_dim(A, B, X)
        allocate(R(N), AR(N))

        k = 0
        R = B - matmul(A, X)
        do while (callback(k, X))
            k = k + 1
            rdotr = dot_product(R, R)
            AR = matmul(A, R)
            rdotAr = dot_product(R, AR)
            alpha = -rdotr / rdotAr
            X = X - alpha * R
            R = R + alpha * AR
        end do
        iterations = k
    end subroutine dense_sd

end module cg
