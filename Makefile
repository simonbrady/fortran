# For a complete clean build, do:
#
# make distclean && make libinstall && make

LIBDIRS=mkl_wrapper iterative_utils cg jacobi tridiag
EXEDIRS=dense_parabolic dense_solver elliptic omp_test sparse_test
SUBDIRS=$(LIBDIRS) $(EXEDIRS)

MARKDOWN=markdown2

.PHONY: $(SUBDIRS) install clean distclean libinstall

all: $(SUBDIRS) doc

lib: $(LIBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

install clean distclean:
	for f in $(SUBDIRS); do \
		$(MAKE) -C $$f $@; \
	done

libinstall:
	for f in $(LIBDIRS); do \
		$(MAKE) -C $$f install; \
	done

doc: README.html

%.html: %.md
	$(MARKDOWN) $< > $@
