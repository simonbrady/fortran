TARGET=libiterative_utils.a

commonpath=../..
include $(commonpath)/lib.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
iterative_utils.o: $(srcdir)/iterative_utils.f90 git_revision.fi

# Additional libraries, e.g. -lfoo
LIBS +=
