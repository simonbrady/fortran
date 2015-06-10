TARGET=libtridiag.a

commonpath=../..
include $(commonpath)/lib.mk
include $(commonpath)/git_revision_fi.mk

# Additional dependencies
tridiag.o: $(srcdir)/tridiag.f90 git_revision.fi

# Additional libraries, e.g. -lfoo
LIBS +=
