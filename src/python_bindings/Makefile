include ../../current/config/Makefile.include

CPPFLAGS += \
 -I.

IMPEX_LIB = `vigra-config --impex_lib`

hourglass: hourglass.o
	$(CXX) -o hourglass hourglass.o $(IMPEX_LIB)

        
ifneq "$(MAKECMDGOALS)" "clean"
include hourglass.d
endif

