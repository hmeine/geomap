PLATFORM = $(shell PATH=/software/Modules/bin:$(PATH) platform)
curdir = .

# should point to build directory of interactive2 (where configure was called):
INT2BUILDDIR = ../../vigra/interactive2/build-$(PLATFORM)
include $(INT2BUILDDIR)/config/Makefile.include

CXXFLAGS = -Wall -Wno-deprecated -O2 -DNDEBUG
CPPFLAGS += -I../../vigra/current/test/include

LIBCPPFLAGS += -I"`python -c 'import numpy; print numpy.get_include()'`"

CELLIMAGE_OBJS = \
	cellimage_module.lo \
	cellimage_nodes.lo \
	cellimage_edges.lo \
	cellimage_faces.lo \
	cellimage_pyramid.lo \
	cellimage_dart+info.lo \
	foureightsegmentation.lo \
	cellstatistics.lo
CELLIMAGE_TARGET = cellimage.la

ALL_TARGETS = $(CELLIMAGE_TARGET) testfoureight

all: $(ALL_TARGETS)

install: $(CELLIMAGE_TARGET)
	$(INSTALL) -d $(dynmoddir)
	$(LIBINSTALL) $(CELLIMAGE_TARGET) $(dynmoddir)

$(CELLIMAGE_TARGET): $(CELLIMAGE_OBJS)
	$(LINKCXXMODULE) -o $(CELLIMAGE_TARGET) $(CELLIMAGE_OBJS) $(BOOST_PYTHON_LIB) vigranumpycore.so

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
	rm -f *.[do] *.l[doa] test core *.plg *ilk
	rm -f seeds.xv res.xv labelImage.xv
	rm -rf Debug ii_files

$(INT2BUILDDIR)/config/Makefile.include:
	@echo
	@echo "*** Check if your vigra/interactive2 build directory is $(INT2BUILDDIR)."
	@echo "*** Otherwise, change the definition of INT2BUILDDIR in the Makefile."
	@exit 1

ifneq "$(MAKECMDGOALS)" "clean"
include $(patsubst %.o, %.d, testfoureight.o testrect.o watershed.o testiterator.o testconfigurations.o test.o)
include $(patsubst %.lo, %.ld, $(CELLIMAGE_OBJS))
endif
