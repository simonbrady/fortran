TARGET=dense_solver

commonpath=../..
include $(commonpath)/exec.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
dense_solver.o: $(srcdir)/dense_solver.f90 git_revision.fi

# Additional libraries, e.g. -lfoo
LIBS += -ljacobi -lcg -literative_utils -lmkl_blas95_lp64
