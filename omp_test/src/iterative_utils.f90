! Common utility functions for iterative solvers - OpenMP version
!
! Simon Brady <simon@hikari.org.nz>
!
! Revision history:
!
! 15-Apr-2014  SJB  Initial version
! 19-Apr-2014  SJB  Add get_dim, residual2

module iterative_utils

    ! Get real kind DP (portable double precision) from module supplied with
    ! Intel Math Kernel Library
    use f95_precision

    ! Fortran 95 interface to the BLAS, also from MKL, used for IAMAX function
    use blas95

    implicit none

    ! Only export entities explicity marked as public
    private

    ! Re-export double precision kind parameter from f95_precision
    public :: DP

    ! Treat magnitudes smaller than this as equivalent to zero
    real(kind=DP), parameter, public :: eps = epsilon(1.0_DP)

    real(kind=DP), parameter, public :: pi = 3.141592653589793238462643_DP

    ! Subroutines exported from this module
    public :: residual, residual2, norm, get_dim
    public :: init_sparse, clear_sparse, dump_sparse, print_sparse
    public :: get_sparse_dim, sparse_entry
    public :: sparse_norm, sparse_residual, sparse_residual2
    public :: sparse_mv_mult, sparse_row_mult, sparse_partial_row_mult

    ! Derived type for sparse matrices in CSR format (see Saad, "Iterative
    ! Methods for Sparse Linear Systems", 2nd ed, p90)
    type, public :: sparse_matrix
        integer :: rows, cols
        real(kind=DP), dimension(:), allocatable :: values
        integer, dimension(:), allocatable :: col_offsets, row_offsets
    end type sparse_matrix

contains

    ! Compute infinty-norm of residual Ax - b
    real(kind=DP) function residual(A, B, X)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        residual = norm(matmul(A, X) - B)
    end function residual

    ! Compute 2-norm of residual Ax - b
    real(kind=DP) function residual2(A, B, X)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        residual2 = norm2(matmul(A, X) - B)
    end function residual2

    ! Infinity norm - IAMAX is BLAS95 interface to BLAS IxAMAX (e.g. IDAMAX
    ! for double precision reals)
    real(kind=DP) function norm(X)
        real(kind=DP), dimension(:), intent(in) :: X

        norm = abs(X(IAMAX(X)))
    end function norm

    ! Get dimension of system and validate that all components are
    ! consistently sized
    integer function get_dim(A, B, X)
        real(kind=DP), dimension(:, :), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        integer, dimension(2) :: input_shape

        input_shape = shape(A)
        if (input_shape(1) /= input_shape(2)) call fatal("A is not square")
        if (size(B) /= input_shape(2)) call fatal("size mismatch " //       &
            "between A and B")
        if (size(X) /= input_shape(2)) call fatal("size mismatch " //       &
            "between A and X")
        get_dim = input_shape(2)
    end function get_dim

    ! Initialise an M*N sparse matrix with space for NNZ non-zero entries,
    ! and set sentinel entry at end of row_offsets
    subroutine init_sparse(A, M, N, NNZ)
        type(sparse_matrix), intent(out) :: A
        integer, intent(in) :: M, N, NNZ

        A%rows = M
        A%cols = N
        allocate(A%values(NNZ), A%col_offsets(NNZ), A%row_offsets(M + 1))
        A%row_offsets(M + 1) = NNZ + 1
    end subroutine init_sparse

    ! Delete storage associated with a sparse matrix so it can be reinitialised
    subroutine clear_sparse(A)
        type(sparse_matrix), intent(in out) :: A

        A%rows = 0
        A%cols = 0
        deallocate(A%values, A%row_offsets, A%col_offsets)
    end subroutine clear_sparse

    ! Print internal structure of sparse matrix A using given format for
    ! entries
    subroutine dump_sparse(A, entry_fmt)
        type(sparse_matrix), intent(in) :: A
        character(len=*), intent(in) :: entry_fmt

        character(len=30) :: intfmt, realfmt

        write(intfmt, "(A, I8, A)") "(", max(size(A%row_offsets),           &
            size(A%col_offsets)), "i8)"
        write(realfmt, "(A, I8, A)") "(", size(A%values), entry_fmt // ")"
        print *, A%rows, A%cols
        print intfmt, A%row_offsets
        print intfmt, A%col_offsets
        print realfmt, A%values
    end subroutine dump_sparse

    ! Print sparse matrix A using given format for entries
    subroutine print_sparse(A, entry_fmt)
        type(sparse_matrix), intent(in) :: A
        character(len=*), intent(in) :: entry_fmt

        real(kind=DP), dimension(:), allocatable :: row
        integer :: i, j
        character(len=30) :: row_fmt

        allocate(row(A%cols))
        write(row_fmt, "(A, I8, A)") "(", A%cols, entry_fmt // ")"
        do i = 1, A%rows
            row = 0
            do j = A%row_offsets(i), A%row_offsets(i + 1) - 1
                row(A%col_offsets(j)) = A%values(j)
            end do
            print row_fmt, row
        end do
    end subroutine print_sparse

    ! Get dimension of system with sparse matrix and validate that all
    ! components are consistently sized
    integer function get_sparse_dim(A, B, X)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        if (A%rows /= A%cols) call fatal("A is not square")
        if (size(B) /= A%cols) call fatal("size mismatch between A and B")
        if (size(X) /= A%cols) call fatal("size mismatch between A and X")
        get_sparse_dim = A%cols
    end function get_sparse_dim

    ! Get entry of sparse matrix A at (row, col)
    real(kind=DP) function sparse_entry(A, row, col)
        type(sparse_matrix), intent(in) :: A
        integer, intent(in) :: row, col

        integer :: i
        real(kind=DP) :: e

        if (row < 1 .or. row > A%rows) call fatal("sparse_entry: invalid row")
        if (col < 1 .or. col > A%cols) call fatal("sparse_entry: invalid col")
        e = 0
        do i = A%row_offsets(row), A%row_offsets(row + 1) - 1
            if (A%col_offsets(i) == col) then
                e =  A%values(i)
                exit
            end if
        end do
        sparse_entry = e
    end function sparse_entry

    ! Infinity-norm (maximum absolute row sum) of sparse matrix A
    real(kind=DP) function sparse_norm(A)
        type(sparse_matrix), intent(in) :: A

        real(kind=DP), dimension(:), allocatable :: row_sum
        integer :: i, p, q

        allocate(row_sum(A%rows))
        do i = 1, A%rows
            p = A%row_offsets(i)
            q = A%row_offsets(i + 1) - 1
            row_sum(i) = sum(A%values(p : q))
        end do
        sparse_norm = norm(row_sum)
    end function sparse_norm

    ! Infinity-norm of residual Ax - b where A is sparse
    real(kind=DP) function sparse_residual(A, B, X)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        real(kind=DP), dimension(:), allocatable :: Y

        allocate(Y(size(X)))
        call sparse_mv_mult(A, X, Y)
        sparse_residual = norm(Y - B)
    end function sparse_residual

    ! 2-norm of residual Ax - b where A is sparse
    real(kind=DP) function sparse_residual2(A, B, X)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: B, X

        real(kind=DP), dimension(:), allocatable :: Y

        allocate(Y(size(X)))
        call sparse_mv_mult(A, X, Y)
        sparse_residual2 = norm2(Y - B)
    end function sparse_residual2

    ! Sparse matrix-vector multiply: compute Y = AX where A is a sparse M*N
    ! matrix, X is a dense N-vector and Y is a dense M-vector
    subroutine sparse_mv_mult(A, X, Y)
        type(sparse_matrix), intent(in) :: A
        real(kind=DP), dimension(:), intent(in) :: X
        real(kind=DP), dimension(:), intent(out) :: Y

        integer :: i

        if (size(X) /= A%cols) call fatal("sparse_mv_mult: size " //        &
            "mismatch between A and X")
        if (size(Y) /= A%rows) call fatal("sparse_mv_mult: size " //        &
            "mismatch between A and Y")
        !$omp parallel do
        do i = 1, A%rows
            Y(i) = sparse_row_mult(A, i, X)
        end do
        !$omp end parallel do
    end subroutine sparse_mv_mult

    ! Multiply vector by single row of a sparse matrix, returning a scalar
    real(kind=DP) function sparse_row_mult(A, row, X)
        type(sparse_matrix), intent(in) :: A
        integer, intent(in) :: row
        real(kind=DP), dimension(:), intent(in) :: X

        integer :: p, q

        p = A%row_offsets(row)
        q = A%row_offsets(row + 1) - 1
        sparse_row_mult = dot_product(A%values(p : q), X(A%col_offsets(p : q)))
    end function sparse_row_mult

    ! Multiply (part of a) vector by partial row of a sparse matrix, returning
    ! a scalar. start_col and end_col specifiy the column range of the matrix
    ! row to consider: non-zero entries in this range will be multiplied by
    ! corresponding elements of the vector, where column start_col of the
    ! matrix corresponds to element start_elem of the vector. For example,
    ! sparse_partial_row_mult(A, row, 1, A%cols, X, 1) is equivalent to
    ! sparse_row_mult(A, row, X).
    real(kind=DP) function sparse_partial_row_mult(A, row, start_col,       &
            end_col, X, start_elem)
        type(sparse_matrix), intent(in) :: A
        integer, intent(in) :: row, start_col, end_col, start_elem
        real(kind=DP), dimension(:), intent(in) :: X

        integer :: p, q, i

        if (start_col > end_col .or. A%row_offsets(row) ==                  &
                A%row_offsets(row + 1)) then
            sparse_partial_row_mult = 0
        else
            p = A%row_offsets(row + 1)
            q = 0
            do i = A%row_offsets(row), A%row_offsets(row + 1) - 1
                if (A%col_offsets(i) >= start_col) then
                    p = i
                    exit
                end if
            end do
            do i = A%row_offsets(row + 1) - 1, p, -1
                if (A%col_offsets(i) <= end_col) then
                    q = i
                    exit
                end if
            end do
            if (p <= q) then
                sparse_partial_row_mult = dot_product(A%values(p : q),      &
                    X(A%col_offsets(p : q) + start_elem - start_col))
            else
                sparse_partial_row_mult = 0
            end if
        end if
    end function sparse_partial_row_mult

end module iterative_utils
