# Common makefile fragment

all: $(TARGET)

srcdir=../src

SOURCES ?= $(patsubst $(srcdir)/%.f90, %, $(wildcard $(srcdir)/*.f90))
OBJECTS := $(addsuffix .o, $(SOURCES))

# Tools and default flags

FC=ifort
FCFLAGS=-mkl -warn all -traceback -axAVX -heap-arrays 
FLFLAGS=-mkl
LIBS=
AR=ar
ARFLAGS=-cru
INSTALL=install

ifeq ($(DEBUG), 1)
FCFLAGS += -g -O0 -check all -DDEBUG=$(DEBUG)
endif

# Target directories for make install

prefix=$(HOME)
exec_prefix=$(prefix)
includedir=$(exec_prefix)/include
moduledir=$(includedir)
ifeq ($(DEBUG), 0)
bindir=$(exec_prefix)/bin
libdir=$(exec_prefix)/lib
else
bindir=$(exec_prefix)/bin/debug
libdir=$(exec_prefix)/lib/debug
endif

FCFLAGS += -I$(moduledir)
FLFLAGS += -L$(libdir)

# Default rules

%.o: $(srcdir)/%.f90
	$(FC) $(FCFLAGS) -c -o $@ $<
	
.PHONY: clean distclean install

clean:
	-rm *.o *.mod git_revision.fi

distclean: clean
	-rm $(TARGET)
