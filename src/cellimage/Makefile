include $(VIGRA_PATH)/config/Makefile.include

		   CXXFLAGS = -O -g -I. -I../../test/include -I../include
	                LIB = -lvigraimpex -ltiff -ljpeg -lm
		 
watershed:: watershed.o
	$(CXX) $(CXXFLAGS) $(CXXINC) -o watershed watershed.o \
         $(CXXLIBDIR) $(LIB)

test:  test.o 
	$(CXX) $(CXXFLAGS) $(CXXINC) -o test test.o $(CXXLIBDIR) $(LIB)
	test
        
save::
	cp  *.cxx *.hxx *.h Makefile \
               /db/versant/save/vigraCD/experiments/cellimage
	
clean::
	rm -f *.o test core res.xv *.plg *ilk
	rm -rf Debug ii_files
	

watershed.o: watershed.cxx foureightsegmentation.hxx configurations.h
test.o: foureightsegmentation.hxx configurations.h pixelneighborhood.hxx
