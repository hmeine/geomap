include ../../current/config/Makefile.include

CXXFLAGS := -g -I. -I../../test/include

all: testneighborhood watershed

watershed: watershed.o
	$(CXX) -o watershed watershed.o `vigra-config --impex-lib`

testneighborhood: testneighborhood.o
	$(CXX) -o testneighborhood testneighborhood.o

testconfigurations: testconfigurations.o
	$(CXX) -o testconfigurations testconfigurations.o

test: test.o
	$(CXX) -o test test.o `vigra-config --impex-lib`
	./test

clean::
	rm -f *.o watershed testneighborhood test core res.xv *.plg *ilk
	rm -rf Debug ii_files

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, watershed.o test.o testneighborhood.o)
endif
