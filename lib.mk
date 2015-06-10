# Common makefile fragment for projects that produce libraries

include $(commonpath)/common.mk

$(TARGET): $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $?

install: $(patsubst %.o, %.mod, $(OBJECTS))
	for m in $^; do echo $$m; $(INSTALL) -m 644 $$m $(moduledir); done
	$(INSTALL) -m 644 $(TARGET) $(libdir)
