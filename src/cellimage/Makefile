include ../../interactive/config/Makefile.include

CXXFLAGS = -g -Wall -Wno-deprecated
CPPFLAGS += -I. -I../../test/include

CELLIMAGE_OBJS = \
	cellimage_module.lo \
	cellimage_nodes.lo \
	cellimage_edges.lo \
	cellimage_faces.lo \
	cellimage_pyramid.lo \
	cellimage_dart+info.lo \
	foureightsegmentation.lo
CELLIMAGE_TARGET = cellimage.la

ALL_TARGETS = testfoureight watershed testiterator testrect $(CELLIMAGE_TARGET)

all: $(ALL_TARGETS)

install: $(CELLIMAGE_TARGET)
	$(INSTALL) -d $(dynmoddir)
	$(LIBINSTALL) $(CELLIMAGE_TARGET) $(dynmoddir)

$(CELLIMAGE_TARGET): $(CELLIMAGE_OBJS)
	$(LINKCXXMODULE) -o $(CELLIMAGE_TARGET) $(CELLIMAGE_OBJS) $(BOOST_PYTHON_LIB)

testfoureight: testfoureight.o foureightsegmentation.o
	$(CXX) -o testfoureight testfoureight.o foureightsegmentation.o `vigra-config --impex-lib`

watershed: watershed.o foureightsegmentation.o
	$(CXX) -o watershed watershed.o foureightsegmentation.o `vigra-config --impex-lib`

testiterator: testiterator.o
	$(CXX) -o testiterator testiterator.o

testrect: testrect.o
	$(CXX) -o testrect testrect.o

testconfigurations: testconfigurations.o
	$(CXX) -o testconfigurations testconfigurations.o

test: test.o
	$(CXX) -o test test.o `vigra-config --impex-lib`
	./test

clean::
	rm -f $(ALL_TARGETS)
	rm -f *.o test core *.plg *ilk
	rm -f seeds.xv res.xv labelImage.xv
	rm -rf Debug ii_files

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, testfoureight.o testrect.o watershed.o testiterator.o testconfigurations.o test.o)
include $(patsubst %.lo, %.ld, $(CELLIMAGE_OBJS))
endif
