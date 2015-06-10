! Wrap MKL include file in a module to speed up compilation for programs
! that use it.

module mkl_wrapper
    include "mkl.fi"
end module mkl_wrapper
