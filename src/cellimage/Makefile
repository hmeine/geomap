include ../../current/config/Makefile.include

CXXFLAGS := -g -I. -I../../test/include

all: testneighborhood testfoureight watershed

testfoureight: testfoureight.o
	$(CXX) -o testfoureight testfoureight.o `vigra-config --impex-lib`

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
	rm -f *.o watershed testneighborhood test core *.plg *ilk
	rm -f seeds.xv res.xv labelImage.xv
	rm -rf Debug ii_files

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, watershed.o test.o testfoureight.o testneighborhood.o)
endif
