include ../../current/config/Makefile.include

CXXFLAGS := -g -I. -I../../test/include -I/usr/local/boost_1_29_0/src -Wall -Wno-deprecated

all: testneighborhood testfoureight watershed

testfoureight: testfoureight.o
	$(CXX) -o testfoureight testfoureight.o `vigra-config --impex-lib`

testrect: testrect.o
	$(CXX) -o testrect testrect.o

watershed: watershed.o
	$(CXX) -o watershed watershed.o `vigra-config --impex-lib`

testit: testit.o
	$(CXX) -o testit testit.o

testconfigurations: testconfigurations.o
	$(CXX) -o testconfigurations testconfigurations.o

test: test.o
	$(CXX) -o test test.o `vigra-config --impex-lib`
	./test

clean::
	rm -f *.o watershed testit test core *.plg *ilk
	rm -f seeds.xv res.xv labelImage.xv
	rm -rf Debug ii_files

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, watershed.o test.o testfoureight.o testrect.o testit.o)
endif
