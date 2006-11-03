PLATFORM = $(shell PATH=/software/Modules/bin:$(PATH) platform)
curdir = .

# should point to build directory of interactive2 (where configure was called):
INT2BUILDDIR = ../../interactive2/build-$(PLATFORM)
include $(INT2BUILDDIR)/config/Makefile.include

#CXXFLAGS += -Wall -g
CXXFLAGS += -Wall
CPPFLAGS += -I.
LIBCPPFLAGS += -I.
TARGET = hourglassmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

OBJS = \
    hourglassmodule.lo \
	polygon.lo \
	statistics.lo \
	delaunay.lo \
	cppmap.lo \
# Ulli added the following to OBJS, uncomment if you need this:
#    $(INT2BUILDDIR)/src/vigracmodule/.libs/pythonimage.o

IMPEX_LIB = `vigra-config --impex_lib`

$(TARGET): $(OBJS)
	$(LINKCXXMODULE) -o $(TARGET) $(OBJS) $(BOOST_PYTHON_LIB) $(dynmoddir)/vigracmodule.la
	test -s $(MODULE) || ln -s .libs/$(MODULE) .

install: $(TARGET)
	$(INSTALL) -d $(dynmoddir)
	$(LIBINSTALL) $(TARGET) $(dynmoddir)

hourglass: hourglass.o
	$(CXX) -o $@ $< $(IMPEX_LIB)

delaunay: delaunay.o deltest.o
	$(CXX) -o $@ delaunay.o deltest.o

$(INT2BUILDDIR)/config/Makefile.include:
	@echo
	@echo "*** Check if your vigra/interactive2 build directory is $(INT2BUILDDIR)."
	@echo "*** Otherwise, change the definition of INT2BUILDDIR in the Makefile."

clean:
	rm -f $(OBJS) $(TARGET) hourglass delaunay *.o *.ld

ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
include delaunay.ld
include $(patsubst %.lo,%.ld,$(OBJS))
endif
