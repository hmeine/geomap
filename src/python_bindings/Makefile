PLATFORM = $(shell /software/Modules/bin/platform)
#PLATFORM = gentoo
curdir = .

INT2BUILDDIR = ../../interactive2/build-$(PLATFORM)
include $(INT2BUILDDIR)/config/Makefile.include

CPPFLAGS += -I.
TARGET = hourglassmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

OBJS = \
    hourglassmodule.lo \
	polygon.lo \
	statistics.lo \
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
	$(CXX) -o hourglass hourglass.o $(IMPEX_LIB)

$(INT2BUILDDIR)/config/Makefile.include:
	@echo
	@echo "*** Check if your vigra/interactive2 build directory is $(INT2BUILDDIR)."
	@echo "*** Otherwise, change the definition of INT2BUILDDIR in the Makefile."

ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
include hourglassmodule.ld
include polygon.ld
include statistics.ld
endif
