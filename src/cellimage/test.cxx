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
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "foureightsegmentation.hxx"

struct FourEightSegmentationTest
{
    BImage image;
    BImage cells;
    IImage labels;
    FourEightSegmentation segmentation;

    FourEightSegmentationTest()
    : image(15, 13), cells(19, 17), labels(19, 17)
    {
        int i;

        char oimage[] = "110111011011011"
                        "010110101101100"
                        "101010101101101"
                        "111100011010011"
                        "111100111011100"
                        "111011010111101"
                        "110110000001011"
                        "001101110110110"
                        "110011000101001"
                        "110010111101111"
                        "101101111100110"
                        "101110000011001"
                        "011110110110101";

        BImage::ScanOrderIterator img = image.begin();

        for(i=0; i<15*13; ++i)   img[i] = (oimage[i] == '1') ? 1 : 0;

        char cimage[] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 0,
            0, 1, 0, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
            0, 2, 1, 0, 2, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 2, 1, 2, 0,
            0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0,
            0, 1, 0, 0, 0, 0, 2, 2, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0,
            0, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 1, 0, 0, 0, 2, 1, 2, 0,
            0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0,
            0, 1, 0, 0, 1, 0, 0, 1, 2, 1, 2, 1, 1, 0, 1, 0, 0, 1, 0,
            0, 2, 1, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 2, 0, 0, 1, 2, 0,
            0, 1, 0, 0, 2, 2, 0, 0, 1, 1, 2, 0, 1, 0, 1, 1, 0, 1, 0,
            0, 1, 0, 0, 2, 2, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
            0, 1, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 2, 1, 0, 0, 1, 2, 0,
            0, 1, 0, 1, 0, 0, 0, 2, 1, 1, 2, 1, 0, 0, 2, 2, 0, 1, 0,
            0, 2, 2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0,
            0, 2, 2, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 2, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

        BImage::ScanOrderIterator cimg = cells.begin();

        for(i=0; i<19*17; ++i)   cimg[i] = cimage[i];

        int limage[] = {
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  1,  1,  1,  3,  3,  3,  2,  6,  6,  3,  8,  8,  4, 10, 10,  5,  0,
            0,  0,  1,  1,  2,  2,  2,  2,  2,  3,  3,  7,  4,  4,  9,  5,  5, 11,  0,
            0,  6, 13,  1,  7,  2,  2,  4,  6,  5,  3,  3,  7,  4,  4,  8, 16,  9,  0,
            0, 12,  7, 13,  7, 14,  2,  4,  6,  5,  3,  3, 10,  4,  4, 15,  8, 17,  0,
            0, 12,  7,  7,  7,  7, 11, 11,  5,  3,  3, 18,  9, 19, 12,  8,  8, 17,  0,
            0, 12,  7,  7,  7,  7, 11, 11,  3,  3,  3, 18,  9,  9,  9, 12, 23, 13,  0,
            0, 12,  7,  7,  7, 20, 10, 10, 21,  3, 18,  9,  9,  9,  9, 22, 11, 24,  0,
            0, 12,  7,  7, 20, 10, 10, 25, 14, 26, 15, 28, 28,  9, 22, 11, 11, 24,  0,
            0, 16, 30, 17, 10, 10, 25, 12, 12, 12, 27, 13, 13, 18, 11, 11, 36, 19,  0,
            0, 29, 14, 14, 17, 17, 12, 12, 34, 34, 20, 13, 35, 15, 36, 36, 15, 37,  0,
            0, 29, 14, 14, 17, 17, 12, 34, 13, 13, 13, 13, 35, 15, 15, 15, 15, 37,  0,
            0, 29, 14, 31, 16, 16, 17, 13, 13, 13, 13, 13, 21, 39, 15, 15, 40, 22,  0,
            0, 29, 14, 31, 16, 16, 16, 17, 33, 33, 23, 38, 17, 17, 24, 24, 18, 41,  0,
            0, 25, 25, 16, 16, 16, 16, 32, 19, 19, 42, 17, 17, 43, 20, 44, 18, 41,  0,
            0, 25, 25, 45, 45, 45, 45, 26, 46, 46, 27, 47, 47, 28, 48, 29, 49, 30,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
        };

        IImage::ScanOrderIterator limg = labels.begin();

        for(i=0; i<19*17; ++i)   limg[i] = limage[i];

        segmentation.init(srcImageRange(image));
    }

    void cellClassificationTest()
    {
        FourEightSegmentation::CellImageIterator c1 = segmentation.cellsUpperLeft();
        BImage::Iterator c2 = cells.upperLeft()+Dist2D(2,2);

        int x,y;
        for(y=-2; y<segmentation.height()+2; ++y)
        {
            for(x=-2; x<segmentation.width()+2; ++x) should(c1(x,y) == c2(x,y));
        }
    }

    void cellLabelingTest()
    {
        FourEightSegmentation::LabelImageIterator l1 = segmentation.labelsUpperLeft();
        IImage::Iterator l2 = labels.upperLeft()+Dist2D(2,2);

        int x,y;
        for(y=-2; y<segmentation.height()+2; ++y)
        {
            for(x=-2; x<segmentation.width()+2; ++x) should(l1(x,y) == l2(x,y));
        }
    }

    void nodeIteratorTest()
    {
        float epsilon = 0.00001;

        FourEightSegmentation::NodeAccessor node;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();
        int i = 0;

        should(node.label(n) == i);
        should(node.x(n) == -1);
        should(node.y(n) == -1);
        should(node.degree(n) == 2);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 2);
        should(node.y(n) == -1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 6);
        should(node.y(n) == -0.5);
        should(node.degree(n) == 4);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 9);
        should(node.y(n) == -1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 12);
        should(node.y(n) == -1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == -1);
        should(node.degree(n) == 2);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == -1);
        should(node.y(n) == 1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 2);
        should(node.y(n) == 1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 13);
        should(node.y(n) == 1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == 1);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 10);
        should(node.y(n) == 2);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 4.5);
        should(node.y(n) == 3.5);
        should(node.degree(n) == 5);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 12.5);
        should(node.y(n) == 3.5);
        should(node.degree(n) == 4);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == 4);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 6);
        should(node.y(n) == 6);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 8);
        should(node.y(n) == 6);
        should(node.degree(n) == 4);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == -1);
        should(node.y(n) == 7);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(abs(node.x(n) - 20.0 / 7.0) < epsilon);
        should(abs(node.y(n) - 62.0 / 7.0) < epsilon);
        should(node.degree(n) == 7);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 11);
        should(node.y(n) == 7);
        should(node.degree(n) == 4);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == 7);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 8);
        should(node.y(n) == 8);
        should(node.degree(n) == 2);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 10);
        should(node.y(n) == 10);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == 10);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 8);
        should(node.y(n) == 11);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 12.5);
        should(node.y(n) == 11);
        should(node.degree(n) == 4);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == -0.5);
        should(node.y(n) == 12.5);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 5);
        should(node.y(n) == 13);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 8);
        should(node.y(n) == 13);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 11);
        should(node.y(n) == 13);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 13);
        should(node.y(n) == 13);
        should(node.degree(n) == 3);
        ++i; ++n;
        should(node.label(n) == i);
        should(node.x(n) == 15);
        should(node.y(n) == 13);
        should(node.degree(n) == 2);
        ++i; ++n;

        should(n == end);
    }

    void rayIteratorTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::NodeAtStartAccessor startnode;
        FourEightSegmentation::EdgeAccessor edge;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::CellImageIterator c = segmentation.cellsUpperLeft();
        FourEightSegmentation::LabelImageIterator l = segmentation.labelsUpperLeft();

        std::vector<int> nodeCount(segmentation.numberOfNodes(), 0);
        std::vector<int> edgeCount(segmentation.numberOfEdges(), 0);

        for(; n != end; ++n)
        {
            FourEightSegmentation::RayCirculator ray = node.rayCirculator(n);
            FourEightSegmentation::RayCirculator end = ray;

            should(startnode.x(ray) == node.x(n));
            should(startnode.y(ray) == node.y(n));

            int degree = 0;

            do
            {
                ++degree;
                ++edgeCount[edge.label(ray)];
                ray.jumpToOpposite();
                ++nodeCount[startnode.label(ray)];
                ray.jumpToOpposite();
            }
            while(++ray != end);

            should(degree == node.degree(n));

            nodeCount[startnode.label(ray)] -= degree;
        }

        int i;
        for(i=0; i<nodeCount.size(); ++i)
        {
            should(nodeCount[i] == 0);
        }
        for(i=0; i<edgeCount.size(); ++i)
        {
            should(edgeCount[i] == 2);
        }
    }

    void rayIteratorTest1()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::NodeAtStartAccessor startnode;
        FourEightSegmentation::EdgeAccessor edge;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::CellImageIterator c = segmentation.cellsUpperLeft();
        FourEightSegmentation::LabelImageIterator l = segmentation.labelsUpperLeft();

        std::vector<int> nodeCount(segmentation.numberOfNodes(), 0);
        std::vector<int> edgeCount(segmentation.numberOfEdges(), 0);

        for(; n != end; ++n)
        {
            FourEightSegmentation::RayCirculator ray = node.rayCirculator(n);
            FourEightSegmentation::RayCirculator end = ray;

            should(startnode.x(ray) == node.x(n));
            should(startnode.y(ray) == node.y(n));

            int degree = 0;

            do
            {
                ++degree;
                ++edgeCount[edge.label(ray)];
                ray.jumpToOpposite();
                ++nodeCount[startnode.label(ray)];
                ray.jumpToOpposite();
            }
            while(--ray != end);

            should(degree == node.degree(n));

            nodeCount[startnode.label(ray)] -= degree;
        }

        int i;
        for(i=0; i<nodeCount.size(); ++i)
        {
            should(nodeCount[i] == 0);
        }
        for(i=0; i<edgeCount.size(); ++i)
        {
            should(edgeCount[i] == 2);
        }
    }

    void edgeIteratorTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::NodeAtStartAccessor startnode;
        FourEightSegmentation::NodeAtEndAccessor endnode;
        FourEightSegmentation::EdgeAccessor edge;

        FourEightSegmentation::EdgeIterator e = segmentation.edgesBegin();
        FourEightSegmentation::EdgeIterator end = segmentation.edgesEnd();
        int i = 0;

        for(; e != end; ++e, ++i)
        {
            should(edge.label(e) == i);

            // check coordinates of start and end nodes
            FourEightSegmentation::NodeIterator n = startnode.nodeIterator(e);

            should(node.x(n) == startnode.x(e));
            should(node.y(n) == startnode.y(e));

            n = endnode.nodeIterator(e);

            should(node.x(n) == endnode.x(e));
            should(node.y(n) == endnode.y(e));

            // check self-loop
            if(startnode.label(e) == endnode.label(e))
            {
                should(startnode.rayCirculator(e) != endnode.rayCirculator(e));
            }

            // check ray circulator at start node
            FourEightSegmentation::RayCirculator r = startnode.rayCirculator(e);
            FourEightSegmentation::NodeIterator  ne = startnode.nodeIterator(e);
            FourEightSegmentation::RayCirculator rt = node.rayCirculator(ne);
            FourEightSegmentation::RayCirculator rtend = rt;

            bool found = false;
            do
            {
                if(r == rt)
                {
                    found = true;
                    break;
                }
            }
            while(++rt != rtend);

            should(found == true);

            // check ray circulator at end node
            r = endnode.rayCirculator(e);
            ne = endnode.nodeIterator(e);
            rt = node.rayCirculator(ne);
            rtend = rt;

            found = false;
            do
            {
                if(r == rt)
                {
                    found = true;
                    break;
                }
            }
            while(++rt != rtend);

            should(found == true);
        }

        should(i == 50);
    }

    void faceIteratorTest()
    {
        FourEightSegmentation::FaceAccessor face;
        FourEightSegmentation::EdgeAccessor edge;
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;
        FourEightSegmentation::FaceAtRightAccessor rightface;
        FourEightSegmentation::NodeAtStartAccessor startnode;

        FourEightSegmentation::FaceIterator f = segmentation.facesBegin();
        FourEightSegmentation::FaceIterator fend = segmentation.facesEnd();

        std::vector<int> nodeCount(segmentation.numberOfNodes(), 0);
        std::vector<int> edgeCount(segmentation.numberOfEdges(), 0);
        std::vector<int> faceCount(segmentation.numberOfFaces(), 0);

        int i = 0;
        for(; f != fend; ++f, ++i)
        {
            should(face.label(f) == i);
            should(face.countBoundaryComponents(f) == 1);

            FourEightSegmentation::BoundaryComponentsIterator b =
                              face.beginBoundaryComponentsIterator(f);
            FourEightSegmentation::BoundaryComponentsIterator bend =
                              face.endBoundaryComponentsIterator(f);
            FourEightSegmentation::ContourCirculator c =
                              face.contourCirculator(b);
            FourEightSegmentation::ContourCirculator cend = c;

            int boundary_length = 0;

            do
            {
                ++boundary_length;

                should(leftface.label(c) == i);

                should(rightface.label(c) < segmentation.numberOfFaces());
                ++faceCount[rightface.label(c)];
                should(startnode.label(c) < segmentation.numberOfNodes());
                ++nodeCount[startnode.label(c)];
                should(edge.label(c) < segmentation.numberOfEdges());
                ++edgeCount[edge.label(c)];
            }
            while(++c != cend);

            faceCount[leftface.label(c)] -= boundary_length;
        }

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator nend = segmentation.nodesEnd();
        for(; n != nend; ++n)
        {
            should(nodeCount[node.label(n)] == node.degree(n));
        }
        FourEightSegmentation::EdgeIterator e = segmentation.edgesBegin();
        FourEightSegmentation::EdgeIterator eend = segmentation.edgesEnd();
        for(; e != eend; ++e)
        {
            should(edgeCount[edge.label(e)] == 2);
        }
        f = segmentation.facesBegin();
        for(; f != fend; ++f)
        {
            should(faceCount[face.label(f)] == 0);
        }
    }
};

struct FourEightSegmentationBadNodeTest
{
    BImage image;

    FourEightSegmentationBadNodeTest()
    : image(5,5)
    {
        int i;

        char oimage[] = "11011"
                        "11011"
                        "00100"
                        "11011"
                        "11011";

        BImage::ScanOrderIterator img = image.begin();

        for(i=0; i<5*5; ++i)   img[i] = (oimage[i] == '1') ? 1 : 0;

    }

    void badNodeTest()
    {
        FourEightSegmentation segmentation;
        bool exceptionCaught = false;
        try
        {
            segmentation.init(srcImageRange(image));
        }
        catch(PreconditionViolation &)
        {
            exceptionCaught = true;
        }

        should(exceptionCaught);
    }

};

struct FourEightSegmentationLabelCirclesTest
{
    BImage image;
    IImage labels;
    FourEightSegmentation segmentation;

    FourEightSegmentationLabelCirclesTest()
    : image(5,5), labels(9,9)
    {
        int i;

        char oimage[] = "11111"
                        "11011"
                        "10101"
                        "11011"
                        "11111";
        BImage::ScanOrderIterator img = image.begin();

        for(i=0; i<5*5; ++i)   img[i] = (oimage[i] == '1') ? 1 : 0;

        int lab[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 1, 1, 1, 1, 1, 1, 0,
                      0, 0, 1, 1, 1, 1, 1, 2, 0,
                      0, 0, 1, 1, 4, 1, 1, 2, 0,
                      0, 0, 1, 4, 2, 4, 1, 2, 0,
                      0, 0, 1, 1, 4, 1, 1, 2, 0,
                      0, 0, 1, 1, 1, 1, 1, 2, 0,
                      0, 2, 3, 3, 3, 3, 3, 3, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0 };

       IImage::ScanOrderIterator limg = labels.begin();

        for(i=0; i<9*9; ++i)   limg[i] = lab[i];

        segmentation.init(srcImageRange(image));
    }

    void labelCirclesTest()
    {
        should(segmentation.numberOfFaces() == 3);
        should(segmentation.numberOfNodes() == 5);
        should(segmentation.numberOfEdges() == 5);

        FourEightSegmentation::LabelImageIterator l1 = segmentation.labelsUpperLeft();
        IImage::Iterator l2 = labels.upperLeft()+Dist2D(2,2);

        int x,y;
        for(y=-2; y<segmentation.height()+2; ++y)
        {
            for(x=-2; x<segmentation.width()+2; ++x) should(l1(x,y) == l2(x,y));
        }
    }

    void leftFaceTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::RayCirculator r = node.rayCirculator(n);
        should(node.degree(n) == 2);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 0);
        ++r;
        should(leftface.label(r) == 1);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 2);
        ++r;
        should(leftface.label(r) == 1);
        ++n;
        should(n == end);
    }

    void contourCirculatorTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::RayCirculator r = node.rayCirculator(n);
        FourEightSegmentation::ContourCirculator c =
                                            leftface.contourCirculator(r);
        FourEightSegmentation::ContourCirculator cend = c;

        int i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 1);
        ++n;
        should(n == end);
    }

    void faceIteratorTest()
    {
        FourEightSegmentation::FaceAccessor face;
        FourEightSegmentation::FaceAtLeftAccessor leftface;
        FourEightSegmentation::FaceAtRightAccessor rightface;
        FourEightSegmentation::NodeAtStartAccessor startnode;
        FourEightSegmentation::NodeAtEndAccessor endnode;

        FourEightSegmentation::FaceIterator f = segmentation.facesBegin();
        FourEightSegmentation::FaceIterator end = segmentation.facesEnd();

        should(face.label(f) == 0);
        should(face.countBoundaryComponents(f) == 1);

        FourEightSegmentation::BoundaryComponentsIterator b =
                          face.beginBoundaryComponentsIterator(f);
        FourEightSegmentation::BoundaryComponentsIterator bend =
                          face.endBoundaryComponentsIterator(f);
        FourEightSegmentation::ContourCirculator c =
                          face.contourCirculator(b);
        FourEightSegmentation::ContourCirculator cend = c;

        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 0);
        should(endnode.label(c) == 1);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 1);
        should(endnode.label(c) == 3);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 3);
        should(endnode.label(c) == 2);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 2);
        should(endnode.label(c) == 0);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(face.label(f) == 1);
        should(face.countBoundaryComponents(f) == 2);

        b = face.beginBoundaryComponentsIterator(f);
        bend = face.endBoundaryComponentsIterator(f);
        c = face.contourCirculator(b);
        cend = c;

        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 0);
        should(endnode.label(c) == 2);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 2);
        should(endnode.label(c) == 3);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 3);
        should(endnode.label(c) == 1);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 1);
        should(endnode.label(c) == 0);
        ++c;
        should(c == cend);
        ++b;
        c = face.contourCirculator(b);
        cend = c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 2);
        should(startnode.label(c) == 4);
        should(endnode.label(c) == 4);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(face.label(f) == 2);
        should(face.countBoundaryComponents(f) == 1);

        b = face.beginBoundaryComponentsIterator(f);
        bend = face.endBoundaryComponentsIterator(f);
        c = face.contourCirculator(b);
        cend = c;
        should(leftface.label(c) == 2);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 4);
        should(endnode.label(c) == 4);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(f == end);
    }
};

struct FourEightSegmentationLabelSelfLoopTest
{
    BImage image;
    IImage labels;
    FourEightSegmentation segmentation;

    FourEightSegmentationLabelSelfLoopTest()
    : image(5,5), labels(9,9)
    {
        int i;

        char oimage[] = "11111"
                        "10011"
                        "10101"
                        "11011"
                        "11111";

        BImage::ScanOrderIterator img = image.begin();

        for(i=0; i<5*5; ++i)   img[i] = (oimage[i] == '1') ? 1 : 0;

        int lab[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 1, 1, 1, 1, 1, 1, 0,
                      0, 0, 1, 1, 1, 1, 1, 2, 0,
                      0, 0, 1, 2, 3, 1, 1, 2, 0,
                      0, 0, 1, 3, 2, 3, 1, 2, 0,
                      0, 0, 1, 1, 3, 1, 1, 2, 0,
                      0, 0, 1, 1, 1, 1, 1, 2, 0,
                      0, 3, 4, 4, 4, 4, 4, 4, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0 };

        IImage::ScanOrderIterator limg = labels.begin();

        for(i=0; i<9*9; ++i)   limg[i] = lab[i];

        segmentation.init(srcImageRange(image));
    }

    void labelCirclesTest()
    {
        should(segmentation.numberOfFaces() == 3);
        should(segmentation.numberOfNodes() == 5);
        should(segmentation.numberOfEdges() == 5);

        FourEightSegmentation::LabelImageIterator l1 = segmentation.labelsUpperLeft();
        IImage::Iterator l2 = labels.upperLeft()+Dist2D(2,2);

        int x,y;
        for(y=-2; y<segmentation.height()+2; ++y)
        {
            for(x=-2; x<segmentation.width()+2; ++x) should(l1(x,y) == l2(x,y));
        }
    }

    void leftFaceTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::RayCirculator r = node.rayCirculator(n);
        should(node.degree(n) == 2);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 2);
        ++r;
        should(leftface.label(r) == 1);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 1);
        ++r;
        should(leftface.label(r) == 0);
        ++n;
        should(node.degree(n) == 2);
        r = node.rayCirculator(n);
        should(leftface.label(r) == 0);
        ++r;
        should(leftface.label(r) == 1);
        ++n;
        should(n == end);
    }

    void contourCirculatorTest()
    {
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator end = segmentation.nodesEnd();

        FourEightSegmentation::RayCirculator r = node.rayCirculator(n);
        FourEightSegmentation::ContourCirculator c =
                                            leftface.contourCirculator(r);
        FourEightSegmentation::ContourCirculator cend = c;

        int i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 1);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        r = node.rayCirculator(n);
        c = leftface.contourCirculator(r);
        cend = c;
        i = 0;
        do { ++i; } while(++c != cend);
        should(i == 4);
        ++n;
        should(n == end);
    }

    void faceIteratorTest()
    {
        FourEightSegmentation::FaceAccessor face;
        FourEightSegmentation::FaceAtLeftAccessor leftface;
        FourEightSegmentation::FaceAtRightAccessor rightface;
        FourEightSegmentation::NodeAtStartAccessor startnode;
        FourEightSegmentation::NodeAtEndAccessor endnode;

        FourEightSegmentation::FaceIterator f = segmentation.facesBegin();
        FourEightSegmentation::FaceIterator end = segmentation.facesEnd();

        should(face.label(f) == 0);
        should(face.countBoundaryComponents(f) == 1);

        FourEightSegmentation::BoundaryComponentsIterator b =
                          face.beginBoundaryComponentsIterator(f);
        FourEightSegmentation::BoundaryComponentsIterator bend =
                          face.endBoundaryComponentsIterator(f);
        FourEightSegmentation::ContourCirculator c =
                          face.contourCirculator(b);
        FourEightSegmentation::ContourCirculator cend = c;

        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 0);
        should(endnode.label(c) == 1);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 1);
        should(endnode.label(c) == 4);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 4);
        should(endnode.label(c) == 3);
        ++c;
        should(leftface.label(c) == 0);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 3);
        should(endnode.label(c) == 0);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(face.label(f) == 1);
        should(face.countBoundaryComponents(f) == 2);

        b = face.beginBoundaryComponentsIterator(f);
        bend = face.endBoundaryComponentsIterator(f);
        c = face.contourCirculator(b);
        cend = c;

        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 0);
        should(endnode.label(c) == 3);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 3);
        should(endnode.label(c) == 4);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 4);
        should(endnode.label(c) == 1);
        ++c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 0);
        should(startnode.label(c) == 1);
        should(endnode.label(c) == 0);
        ++c;
        should(c == cend);
        ++b;
        c = face.contourCirculator(b);
        cend = c;
        should(leftface.label(c) == 1);
        should(rightface.label(c) == 2);
        should(startnode.label(c) == 2);
        should(endnode.label(c) == 2);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(face.label(f) == 2);
        should(face.countBoundaryComponents(f) == 1);

        b = face.beginBoundaryComponentsIterator(f);
        bend = face.endBoundaryComponentsIterator(f);
        c = face.contourCirculator(b);
        cend = c;
        should(leftface.label(c) == 2);
        should(rightface.label(c) == 1);
        should(startnode.label(c) == 2);
        should(endnode.label(c) == 2);
        ++c;
        should(c == cend);
        ++b;
        should(b == bend);

        ++f;
        should(f == end);
    }
};

struct FourEightSegmentationRealisticTest
{
    BImage image;
    FourEightSegmentation segmentation;

    FourEightSegmentationRealisticTest()
    {
        ImageImportInfo info("testedges.xv");

        int w = info.width();
        int h = info.height();
        image.resize(w, h);

        importImage(info, destImage(image));


        segmentation.init(srcImageRange(image));
    }

    void consistencyCheck()
    {
        FourEightSegmentation::FaceAccessor face;
        FourEightSegmentation::EdgeAccessor edge;
        FourEightSegmentation::NodeAccessor node;
        FourEightSegmentation::FaceAtLeftAccessor leftface;
        FourEightSegmentation::FaceAtRightAccessor rightface;
        FourEightSegmentation::NodeAtStartAccessor startnode;

        FourEightSegmentation::FaceIterator f = segmentation.facesBegin();
        FourEightSegmentation::FaceIterator fend = segmentation.facesEnd();

        std::vector<int> nodeCount(segmentation.numberOfNodes(), 0);
        std::vector<int> edgeCount(segmentation.numberOfEdges(), 0);
        std::vector<int> faceCount(segmentation.numberOfFaces(), 0);

        int i = 0;
        for(; f != fend; ++f, ++i)
        {
            should(face.label(f) == i);

            FourEightSegmentation::BoundaryComponentsIterator b =
                              face.beginBoundaryComponentsIterator(f);
            FourEightSegmentation::BoundaryComponentsIterator bend =
                              face.endBoundaryComponentsIterator(f);

            for(; b != bend; ++b)
            {
                FourEightSegmentation::ContourCirculator c =
                                  face.contourCirculator(b);
                FourEightSegmentation::ContourCirculator cend = c;

                int boundary_length = 0;

                do
                {
                    ++boundary_length;

                    should(leftface.label(c) == i);

                    should(rightface.label(c) < segmentation.numberOfFaces());
                    ++faceCount[rightface.label(c)];
                    should(startnode.label(c) < segmentation.numberOfNodes());
                    ++nodeCount[startnode.label(c)];
                    should(edge.label(c) < segmentation.numberOfEdges());
                    ++edgeCount[edge.label(c)];
               }
                while(++c != cend);

                faceCount[leftface.label(c)] -= boundary_length;
            }
        }

        FourEightSegmentation::NodeIterator n = segmentation.nodesBegin();
        FourEightSegmentation::NodeIterator nend = segmentation.nodesEnd();
        for(; n != nend; ++n)
        {
            should(nodeCount[node.label(n)] == node.degree(n));
        }
        FourEightSegmentation::EdgeIterator e = segmentation.edgesBegin();
        FourEightSegmentation::EdgeIterator eend = segmentation.edgesEnd();
        for(; e != eend; ++e)
        {
            should(edgeCount[edge.label(e)] == 2);
        }
        f = segmentation.facesBegin();
        for(; f != fend; ++f)
        {
            should(faceCount[face.label(f)] == 0);
        }
    }
};

struct FourEightSegmentationTestSuite
: public TestSuite
{
    FourEightSegmentationTestSuite()
    : TestSuite("FourEightSegmentationTestSuite")
    {
        add( testCase( &FourEightSegmentationTest::cellClassificationTest));
        add( testCase( &FourEightSegmentationTest::cellLabelingTest));
        add( testCase( &FourEightSegmentationTest::nodeIteratorTest));
        add( testCase( &FourEightSegmentationTest::rayIteratorTest));
        add( testCase( &FourEightSegmentationTest::rayIteratorTest1));
        add( testCase( &FourEightSegmentationTest::edgeIteratorTest));
        add( testCase( &FourEightSegmentationTest::faceIteratorTest));
        add( testCase( &FourEightSegmentationBadNodeTest::badNodeTest));
        add( testCase( &FourEightSegmentationLabelCirclesTest::labelCirclesTest));
        add( testCase( &FourEightSegmentationLabelCirclesTest::leftFaceTest));
        add( testCase( &FourEightSegmentationLabelCirclesTest::contourCirculatorTest));
        add( testCase( &FourEightSegmentationLabelCirclesTest::faceIteratorTest));
        add( testCase( &FourEightSegmentationLabelSelfLoopTest::labelCirclesTest));
        add( testCase( &FourEightSegmentationLabelSelfLoopTest::leftFaceTest));
        add( testCase( &FourEightSegmentationLabelSelfLoopTest::contourCirculatorTest));
        add( testCase( &FourEightSegmentationLabelSelfLoopTest::faceIteratorTest));
        add( testCase( &FourEightSegmentationRealisticTest::consistencyCheck));
    }
};

int main()
{
    FourEightSegmentationTestSuite test;

    test.run();

    std::cout << test.report() << std::endl;

    return 0;
}


