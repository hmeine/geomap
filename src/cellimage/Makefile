include ../../current/config/Makefile.include

CXXFLAGS := -O -g -I. -I../../test/include

watershed:: watershed.o
	$(CXX) -o watershed watershed.o `vigra-config --impex-lib`

test:  test.o
	$(CXX) -o test test.o `vigra-config --impex-lib`
	./test

clean::
	rm -f *.o test core res.xv *.plg *ilk
	rm -rf Debug ii_files

watershed.o: watershed.cxx foureightsegmentation.hxx configurations.h
test.o: foureightsegmentation.hxx configurations.h pixelneighborhood.hxx
