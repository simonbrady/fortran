TARGET=libmkl_wrapper.a

commonpath=../..
include $(commonpath)/lib.mk

FCFLAGS += -nowarn

# Additional libraries, e.g. -lfoo
LIBS +=
