TARGET=elliptic

commonpath=../..
include $(commonpath)/exec.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
elliptic.o: $(srcdir)/elliptic.f90 git_revision.fi

FCFLAGS += -fpp

# Additional libraries, e.g. -lfoo
LIBS += -ljacobi -lcg -literative_utils -lmkl_blas95_lp64
