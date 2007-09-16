PLATFORM = $(shell PATH=/software/Modules/bin:$(PATH) platform)
curdir = .

# should point to build directory of interactive2 (where configure was called):
INT2BUILDDIR = ../../interactive2/build-$(PLATFORM)
include $(INT2BUILDDIR)/config/Makefile.include

# at KOGS, try: module add pkgconfig libsigc++
LOCAL_CPPFLAGS = -I. -I../cellimage `pkg-config --cflags sigc++-2.0`

# define this if you have installed John Maddock's math_toolkit from
# the boost vault (which is needed for proper t-test support), see
# http://johnmaddock.co.uk/ or just do the following here in hourglass/:
#  svn co http://svn.boost.org/svn/boost/sandbox/math_toolkit/boost
#LOCAL_CPPFLAGS += -DHAVE_MATH_TOOLKIT

CXXFLAGS += -Wall -Wno-strict-aliasing
#CXXFLAGS += -g -W
CXXFLAGS += -O3 -DNDEBUG
CPPFLAGS := $(LOCAL_CPPFLAGS) $(CPPFLAGS)
LIBCPPFLAGS := $(LOCAL_CPPFLAGS) $(LIBCPPFLAGS)
TARGET = hourglassmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

OBJS = \
	shapecontext.lo \
    hourglassmodule.lo \
	polygon.lo \
	statistics.lo \
	delaunay.lo \
	cppmap.lo \
	cppmapmodule.lo \
	cppmapmodule_stats.lo \
	dsl.lo \

IMPEX_LIB = `vigra-config --impex_lib`

$(TARGET): $(OBJS)
	$(LINKCXXMODULE) -o $(TARGET) $(OBJS) $(BOOST_PYTHON_LIB) `pkg-config --libs sigc++-2.0` `vigra-interactive-config --pythonimage-lib`
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
	rm -f $(OBJS) $(TARGET) hourglass delaunay *.o *.ld *.d

ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
include $(patsubst %.lo,%.ld,$(OBJS))
endif
