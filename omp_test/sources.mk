TARGET=omp_test

# Override default compilation order
SOURCES=iterative_utils cg jacobi fatal omp_test

commonpath=../..
include $(commonpath)/exec.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
omp_test.o: $(srcdir)/omp_test.f90 git_revision.fi

FCFLAGS += -openmp

# Additional libraries, e.g. -lfoo
LIBS += -lmkl_blas95_lp64
