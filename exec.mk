# Common makefile fragment for projects that produce executables

include $(commonpath)/common.mk

$(TARGET): $(OBJECTS)
	$(FC) $(FLFLAGS) -o $@ $^ $(LIBS)

install:
	$(INSTALL) -m 755 $(TARGET) $(bindir)
 