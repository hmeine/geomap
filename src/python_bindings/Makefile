PLATFORM = $(shell /software/Modules/bin/platform)
curdir = .

include ../../interactive2/build.$(PLATFORM)/config/Makefile.include

CPPFLAGS += -I. 
TARGET = hourglassmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

OBJS = \
    hourglassmodule.lo \
    ../../$(PLATFORM)/lib/pythonimage.o


IMPEX_LIB = `vigra-config --impex_lib`

$(TARGET): $(OBJS)
	$(LINKCXXMODULE) -o $(TARGET) $(OBJS) $(BOOST_PYTHON_LIB)
	test -s $(MODULE) || ln -s .libs/$(MODULE) .

hourglass: hourglass.o
	$(CXX) -o hourglass hourglass.o $(IMPEX_LIB)

        
ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
include hourglassmodule.ld
endif

