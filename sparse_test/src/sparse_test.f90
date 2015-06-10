! Test sparse_matrix type in iterative utils
!
! Expected output:
!   5-row tridiagonal (-1, 4, -1) matrix
!   scalars 0, 4, 0
!   row vector (2, 4, 6, 8, 16)
!   column vector (-3, 8, 8, 16, -5)

program sparse_test
    use iterative_utils

    implicit none
    include "git_revision.fi"

    integer, parameter :: N = 5
    type(sparse_matrix) :: A
    real(kind=DP), dimension(N) :: X, Y
    integer :: i, offset

    print *, 'sparse_test (git ', git_revision, ')'
    print *, 'built with ', trim(iterative_utils_version())
    call init_sparse(A, N, N, 3 * N - 2)
    offset = 1
    do i = 1, N
        X(i) = i
        A%row_offsets(i) = offset
        if (i > 1) then
            A%values(offset) = -1
            A%col_offsets(offset) = i - 1
            offset = offset + 1
        end if
        A%values(offset) = 4
        A%col_offsets(offset) = i
        offset = offset + 1
        if (i < N) then
            A%values(offset) = -1
            A%col_offsets(offset) = i + 1
            offset = offset + 1
        end if
    end do
    call dump_sparse(A, "f5.1")
    print *, ""
    call print_sparse(A, "f5.1")
    print *, ""
    print *, sparse_entry(A, 1, 3)
    print *, sparse_entry(A, 3, 3)
    print *, sparse_entry(A, 5, 3)
    print *, ""
    call sparse_mv_mult(A, X, Y)
    print *, Y
    print *, ""
    do i = 1, N
        print *, sparse_partial_row_mult(A, i, 2, 4, X, 3)
    end do
    call clear_sparse(A)

end program sparse_test
