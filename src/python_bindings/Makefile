
PLATFORM = $(shell PATH=/software/Modules/bin:$(PATH) platform)
curdir = .

include ../Makefile.include

LOCAL_CPPFLAGS = -I. -I../cellimage \
	-I"`python -c 'import numpy; print numpy.get_include()'`"

# define this if you have installed John Maddock's math_toolkit from
# the boost vault (which is needed for proper t-test support), see
# http://johnmaddock.co.uk/ or just do the following here in hourglass/:
#  svn co http://svn.boost.org/svn/boost/sandbox/math_toolkit/boost
#LOCAL_CPPFLAGS += -DHAVE_MATH_TOOLKIT

CXXFLAGS += -Wall -Wno-strict-aliasing
#CXXFLAGS += -g -W
CXXFLAGS += -pipe -DNDEBUG
CPPFLAGS := $(LOCAL_CPPFLAGS) $(CPPFLAGS)
LIBCPPFLAGS := $(LOCAL_CPPFLAGS) $(LIBCPPFLAGS)
TARGET = geomapmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

#	channelrepresentation.lo \

OBJS = \
	vectorconv.lo \
	pythondiff2d.lo \
	shapecontext.lo \
	geomapmodule.lo \
	polygon.lo \
	statistics.lo \
	features.lo \
	ellipse.lo \
	cppmapmodule_utils.lo \
	cppmapmodule_stats.lo \
	cppmapmodule.lo \
	cppmap.lo \
	cppmap_utils.lo \
	crackedgemap.lo \
	dsl.lo \

$(TARGET): $(OBJS)
	$(LINKCXXMODULE) -o $(TARGET) $(OBJS) $(BOOST_PYTHON_LIB) -lboost_signals vigranumpycore.so
	test -s $(MODULE) || ln -s .libs/$(MODULE) .

test: $(TARGET)
	python test.py

install: $(TARGET)
	$(INSTALL) -d $(dynmoddir)
	$(LIBINSTALL) $(TARGET) $(dynmoddir)

$(INT2BUILDDIR)/config/Makefile.include:
	@echo
	@echo "*** Check if your vigra/interactive2 build directory is $(INT2BUILDDIR)."
	@echo "*** Otherwise, change the definition of INT2BUILDDIR in the Makefile."
	@exit 1

clean:
	rm -f $(OBJS) $(TARGET) *.o *.ld *.d *.so

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.lo,%.ld,$(OBJS))
endif
