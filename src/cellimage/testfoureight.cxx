/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/impex.hxx>
#include "foureightsegmentation.hxx"
#include "mydebug.hxx"
#include "debugimage.hxx"
#include "crop.hxx"
#include "cellimage.hxx"
#include <unittest.hxx>

using namespace vigra;
using namespace vigra::functor;
using namespace vigra::cellimage;

typedef GeoMap::DartTraverser DartTraverser;

template<class ARRAY>
std::ostream &outputArray(std::ostream &s, ARRAY &d, unsigned int len)
{
    s << "{";
    for(unsigned int i= 0; i < len - 1; ++i)
        s << d[i] << ", ";
    if(len)
        s << d[len - 1];
    s << "}";
    return s;
}

template<class T>
std::ostream &operator <<(std::ostream &s, std::vector<T> d)
{
    outputArray(s, d, d.size());
    return s;
}

void testNodeDegrees(GeoMap &seg,
                     unsigned int expected[], unsigned int sizeofExpected)
{
    unsigned int lenExpected = sizeofExpected / sizeof(unsigned int);
    shouldEqual(lenExpected, seg.nodeCount());
    bool goodLen = lenExpected == seg.nodeCount();
    if(!goodLen)
        std::cerr << "wrong number (" << lenExpected
                  << ") of node degrees expected, nodeCount is "
                  << seg.nodeCount() << "\n";
    std::vector<unsigned int> real(seg.nodeCount());
    unsigned int pos = 0;
    bool goodData = true;
    for(GeoMap::NodeIterator
            node= seg.nodesBegin(); node.inRange(); ++node, ++pos)
    {
        real[pos] = node->degree;
        if(goodLen && (expected[pos] != real[pos]))
        {
            std::cerr << "wrong degree for node " << node->label
                      << " (position " << pos << "): got "
                      << real[pos] << ", expected " << expected[pos] << "\n";
            goodData = false;
        }
    }
    if(!(goodLen && goodData))
    {
        std::cerr << "EXPECTED ";
        outputArray(std::cerr, expected, lenExpected) << ",\n"
            "GOT " << real << std::endl;
    }
    should(goodData);
}

void validateAnchors(GeoMap &seg)
{
    for(GeoMap::NodeIterator
            node= seg.nodesBegin(); node.inRange(); ++node)
    {
        shouldEqual(node->anchor.isSingular(), node->degree == 0);
    }
}

struct GeoMapTest
{
    cellimage::GeoMap *segmentation;

    GeoMapTest()
    {
        BImage image;
        ImageImportInfo info("testboundaries.png");
        image.resize(info.size());
        importImage(info, destImage(image));

        std::cerr << "creating GeoMap from testboundaries.png..";
        segmentation = new GeoMap(srcImageRange(image), 0,
												 CellTypeVertex);
        std::cerr << "done.\n";

        debugImage(srcImageRange(segmentation->cellImage), std::cerr, 2);
    }

    void test()
    {
        validateAnchors(*segmentation);

        unsigned int nodeDegrees[] =
            {2, 2, 2, 2, 3, 2, 4, 3, 2, 3, 3, 2, 4, 2, 2, 2, 2, 2, 2};
        testNodeDegrees(*segmentation, nodeDegrees, sizeof(nodeDegrees));
    }
};

struct ConsistencyTest
{
    cellimage::GeoMap *segmentation;

    ConsistencyTest()
    {
        IImage image;
        ImageImportInfo info("labels.xv");
        image.resize(info.size());
        importImage(info, destImage(image));

        std::cerr << "creating GeoMap from labels.xv..";
        segmentation = new GeoMap(srcImageRange(image), 0,
												 CellTypeVertex);
        std::cerr << "done.\n";
    }

	void debugRange()
	{
		debugImage(crop(srcImageRange(segmentation->cellImage),
						Rect2D(309, 136, 323, 157)), std::cerr, 4);
	}

	void test()
	{
		segmentation->mergeFaces(segmentation->edge(1357).start);
		segmentation->mergeFaces(segmentation->edge(1423).start);
		segmentation->mergeEdges(segmentation->node(873).anchor);
		segmentation->mergeEdges(segmentation->node(883).anchor);
		segmentation->mergeFaces(segmentation->edge(1478).start);
		segmentation->mergeEdges(segmentation->node(934).anchor);
		segmentation->mergeEdges(segmentation->node(911).anchor);
		segmentation->mergeFaces(segmentation->edge(1521).start);
		segmentation->mergeEdges(segmentation->node(917).anchor);
	}
};

struct GeoMapTestSuite
: public vigra::test_suite
{
    GeoMapTestSuite()
    : vigra::test_suite("GeoMapTestSuite")
    {
        add(testCase(&GeoMapTest::test));
        add(testCase(&ConsistencyTest::test));
    }
};

int main(int argc, char ** argv)
{
    GeoMapTestSuite geoMapTestSuite;
    int failed = geoMapTestSuite.run();
    std::cout << geoMapTestSuite.report() << std::endl;

    return failed;
}
