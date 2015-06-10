all: debug release

.PHONY: debug release install clean distclean

debugdir=Debug_Intel64
releasedir=Release_Intel64

debug:
	@$(MAKE) -C $(debugdir)

release:
	@$(MAKE) -C $(releasedir)

SUBDIRS=$(debugdir) $(releasedir)

install: debug release
	@for d in $(SUBDIRS); do $(MAKE) -C $$d install; done

clean distclean:
	@for d in $(SUBDIRS); do $(MAKE) -C $$d $@; done
