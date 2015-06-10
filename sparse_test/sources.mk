TARGET=sparse_test

commonpath=../..
include $(commonpath)/exec.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
sparse_test.o: $(srcdir)/sparse_test.f90 git_revision.fi

# Additional libraries, e.g. -lfoo
LIBS += -literative_utils -lmkl_blas95_lp64
