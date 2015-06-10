TARGET=dense_parabolic

commonpath=../..
include $(commonpath)/exec.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
dense_parabolic.o: $(srcdir)/dense_parabolic.f90 git_revision.fi

# Additional libraries, e.g. -lfoo
LIBS += -ltridiag
