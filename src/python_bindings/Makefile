include ../../interactive2/config/Makefile.include

CPPFLAGS += -I. 
TARGET = hourglassmodule.la
MODULE = $(patsubst %.la,%.so,$(TARGET))

OBJS = \
    hourglassmodule.lo \
    ../../interactive2/src/vigracmodule/pythonimage.lo


IMPEX_LIB = `vigra-config --impex_lib`

hourglass: hourglass.o
	$(CXX) -o hourglass hourglass.o $(IMPEX_LIB)

        
$(TARGET): $(OBJS)
	$(LINKCXXMODULE) -o $(TARGET) $(OBJS) $(BOOST_PYTHON_LIB)
	test -s $(MODULE) || ln -s .libs/$(MODULE) .

ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
include hourglassmodule.ld
endif

