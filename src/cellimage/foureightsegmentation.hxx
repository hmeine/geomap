#ifndef VIGRA_FOUREIGHTSEGMENTATION_HXX
#define VIGRA_FOUREIGHTSEGMENTATION_HXX

#include <vigra/error.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <vigra/contourcirculator.hxx>
//#include "circulatoradaptor.hxx"
#include "celltypes.hxx"
#include "filteriterator.hxx"
#include "cellimage.hxx"
#include "rect2d.hxx"

// temporary,debug
#include <vigra/impex.hxx>
#include <iostream>
#include "debugimage.hxx"
#include "mydebug.hxx"

namespace vigra {

namespace CellImage {

class FourEightSegmentation;

// -------------------------------------------------------------------
//                            LabelScanIterator
// -------------------------------------------------------------------
template<class LabelImageIterator, class ImageIterator = Diff2D>
class LabelScanIterator
{
    LabelImageIterator cellLR_, cellIter_;
    typename LabelImageIterator::value_type cellPixelValue_;
    ImageIterator imageIter_;
    unsigned int width_;

public:
        /** the iterator's value type
        */
    typedef typename ImageIterator::value_type value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename ImageIterator::reference reference;

        /** the iterator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename ImageIterator::pointer pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    LabelScanIterator()
    {}

    LabelScanIterator(LabelImageIterator cellUL, LabelImageIterator cellLR,
                      typename LabelImageIterator::value_type cellPixelValue,
                      ImageIterator imageIter = ImageIterator())
        : cellLR_(cellLR), cellIter_(cellUL),
          cellPixelValue_(cellPixelValue),
          imageIter_(imageIter),
          width_(cellLR.x - cellUL.x)
    {
        if(cellIter_ != cellLR_ && *cellIter_ != cellPixelValue_)
            operator++();
    }

    LabelScanIterator & operator++()
    {
        ++cellIter_.x, ++imageIter_.x;
        while((cellIter_.x != cellLR_.x) && (*cellIter_ != cellPixelValue_))
            ++cellIter_.x, ++imageIter_.x;

        if(cellIter_.x == cellLR_.x)
        {
            cellIter_.x -= width_, imageIter_.x -= width_;
            ++cellIter_.y, ++imageIter_.y;

            if(cellIter_.y != cellLR_.y)
                operator++();
            else
                cellIter_ = cellLR_;
        }
        return *this;
    }

    LabelScanIterator operator++(int)
    {
        LabelScanIterator ret(*this);
        operator++();
        return ret;
    }

	bool atEnd() const
	{
		return cellIter_ == cellLR_;
	}

	bool inRange() const
	{
		return cellIter_ != cellLR_;
	}

	bool operator==(LabelScanIterator const &other) const
    {
        return cellIter_ == other.cellIter_;
    }

    bool operator!=(LabelScanIterator const &other) const
    {
        return cellIter_ != other.cellIter_;
    }

    reference operator*() const
    {
        return *imageIter_;
    }

    pointer operator->() const
    {
        return imageIter_.operator->();
    }
};

// -------------------------------------------------------------------
//                             EdgelIterator
// -------------------------------------------------------------------
/**
 * An EdgelIterator starts walking in the direction of the Circulator
 * given on construction and walks along CellTypeLine type'ed
 * pixels. atEnd() will become true if the iterator steps on a
 * CellTypeVertex pixel or tries to pass one diagonally.
 *
 * It is used by the DartTraverser for nextAlpha() and in labelEdge().
 */
class EdgelIterator
{
    CellImageEightCirculator neighborCirc_;
    bool atEnd_;

public:
    EdgelIterator(CellImageEightCirculator const & n)
        : neighborCirc_(n), atEnd_(false)
    {}

    CellImageEightCirculator::reference operator*() const
    {
        return *neighborCirc_;
    }

    CellImageEightCirculator::pointer operator->() const
    {
        return neighborCirc_.operator->();
    }

    bool atEnd() const
    {
        return atEnd_;
    }

    bool inRange() const
    {
        return !atEnd_;
    }

    EdgelIterator & operator++()
    {
        neighborCirc_.moveCenterToNeighbor();
        neighborCirc_.turnRight();

        while(1)
        {
            if(neighborCirc_->type() == CellTypeVertex)
            {
                atEnd_ = true;
                break;
            }
            if(neighborCirc_->type() == CellTypeLine)
            {
                break;
            }
            ++neighborCirc_;
        }

        if(neighborCirc_.isDiagonal() &&
           neighborCirc_[1].type() == CellTypeVertex)
        {
            ++neighborCirc_;
            atEnd_ = true;
        }

        return *this;
    }

    EdgelIterator & jumpToOpposite()
    {
        while(!atEnd())
            operator++();
        neighborCirc_.swapCenterNeighbor();
        return *this;
    }

    operator CellImageEightCirculator()
    {
        return neighborCirc_;
    }
};

// -------------------------------------------------------------------
//                         FourEightSegmentation
// -------------------------------------------------------------------
class FourEightSegmentation
{
public:
    struct NodeInfo;
    struct EdgeInfo;
    struct FaceInfo;

    // -------------------------------------------------------------------
    //                  FourEightSegmentation::DartTraverser
    // -------------------------------------------------------------------
    class DartTraverser
    {
        CellImageEightCirculator neighborCirc_;
        FourEightSegmentation * segmentation_;
        bool isSingular_;

            // prevent double stop at a line pixel from different source
            // vertex pixels
        bool badDiagonalConfig() const
        {
            return (neighborCirc_->type() == CellTypeLine &&
                    (neighborCirc_[1].type() == CellTypeVertex ||
                     neighborCirc_[-1].type() == CellTypeVertex));
        }

    public:
        DartTraverser() : segmentation_(0L) {}

        DartTraverser(FourEightSegmentation * segmentation,
                      CellImageEightCirculator const & circ)
            : neighborCirc_(circ),
              segmentation_(segmentation)
        {
            vigra_precondition(neighborCirc_.center()->type() == CellTypeVertex,
            "FourEightSegmentation::DartTraverser(): center is not a node");

            vigra_precondition(neighborCirc_->type() != CellTypeVertex,
            "FourEightSegmentation::DartTraverser(): neighbor is a node");

            while(badDiagonalConfig())
            {
                ++neighborCirc_;
                // pointing from vertex to vertex pixel now
                neighborCirc_.swapCenterNeighbor();
                // still pointing from vertex to vertex pixel
                ++neighborCirc_;
                // finally, pointing from vertex to line pixel
            }

            isSingular_ = false;
            CellImageEightCirculator nend = neighborCirc_;
            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                }
                tryNextSigma();
                if(neighborCirc_ == nend)
                {
                    // did not find any adjacent line pixel
                    isSingular_ = true;
                    break;
                }
            }
            /*CellImageEightCirculator n = neighborCirc_;
              isSingular_ = true;
              do
              {
              if(n->type() != CellTypeRegion)
              {
              isSingular_ = false;
              break;
              }
              }
              while(++n != neighborCirc_);

              if(neighborCirc_->type() != CellTypeLine)
              operator++();*/
        }

        DartTraverser & nextAlpha()
        {
            if(isSingular_)
                return *this;

            EdgelIterator line(neighborCirc_);
            line.jumpToOpposite();

            neighborCirc_ = line;

            return *this;
        }

        DartTraverser & prevAlpha()
        {
            return nextAlpha();
        }

        DartTraverser & nextPhi()
        {
            return nextAlpha().prevSigma();
        }

        DartTraverser & prevPhi()
        {
            return nextSigma().prevAlpha();
        }

        DartTraverser & nextSigma()
        {
            if(isSingular_)
                return *this;

            tryNextSigma();

            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                }
                tryNextSigma();
            }
            return *this;
        }

        DartTraverser & prevSigma()
        {
            if(isSingular_)
                return *this;

            tryPrevSigma();

            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                }
                tryPrevSigma();
            }
            return *this;
        }

        bool isSingular() const
        {
            return isSingular_;
        }

        CellPixel::LabelType edgeLabel() const
        {
            return neighborCirc_->label();
        }

        CellPixel::LabelType leftFaceLabel() const
        {
            vigra_precondition(neighborCirc_[1].type() == CellTypeRegion,
                "insufficient algorithm for DartTraverser::rightFace()");
            return neighborCirc_[1].label();
        }

        CellPixel::LabelType rightFaceLabel() const
        {
            vigra_precondition(neighborCirc_[-1].type() == CellTypeRegion,
                "insufficient algorithm for DartTraverser::rightFace()");
            return neighborCirc_[-1].label();
        }

        NodeInfo & startNode()
        {
            return segmentation_->nodeList[neighborCirc_.center()->label()];
        }

        NodeInfo & endNode()
        {
            DartTraverser reverseDart_(*this);
            reverseDart_.nextAlpha();
            return reverseDart_.startNode();
        }

        EdgeInfo & edge()
        {
            return segmentation_->edgeList[edgeLabel()];
        }

        FaceInfo & leftFace()
        {
            return segmentation_->faceList[leftFaceLabel()];
        }

        FaceInfo & rightFace()
        {
            return segmentation_->faceList[rightFaceLabel()];
        }

        bool operator==(DartTraverser const & o) const
        {
            return neighborCirc_ == o.neighborCirc_;
        }

        bool operator!=(DartTraverser const & o) const
        {
            return neighborCirc_ != o.neighborCirc_;
        }

        FourEightSegmentation * segmentation() const
        {
            return segmentation_;
        }

        const CellImageEightCirculator &neighborCirculator() const
        {
            return neighborCirc_;
        }

    private:
        void tryNextSigma()
        {
            ++neighborCirc_;

            if(badDiagonalConfig())
                ++neighborCirc_;
        }

        void tryPrevSigma()
        {
            --neighborCirc_;

            if(badDiagonalConfig())
                --neighborCirc_;
        }
    };

    // -------------------------------------------------------------------

    struct CellInfo
    {
        CellPixel::LabelType label;
        Rect2D bounds;

        CellInfo() : label(NumericTraits<CellPixel::LabelType>::max()) {}

        bool initialized() const
        { return label != NumericTraits<CellPixel::LabelType>::max(); }

        void uninitialize()
        { label = NumericTraits<CellPixel::LabelType>::max(); }
    };

    struct NodeInfo : public CellInfo
    {
        DartTraverser anchor;
        // remove all following members?
        float centerX, centerY;
        int size;
    };

    struct EdgeInfo : public CellInfo
    {
        DartTraverser start, end;
    };

    struct FaceInfo : public CellInfo
    {
        std::vector<DartTraverser> contours;
    };

    typedef std::vector<NodeInfo> NodeList;
    typedef std::vector<EdgeInfo> EdgeList;
    typedef std::vector<FaceInfo> FaceList;

    struct FilterInitialized
    {
        bool operator()(CellInfo const &c) const
        {
            return c.initialized();
        }
    };

    typedef FilterIterator<NodeList::iterator, FilterInitialized> NodeIterator;
    typedef FilterIterator<EdgeList::iterator, FilterInitialized> EdgeIterator;
    typedef FilterIterator<FaceList::iterator, FilterInitialized> FaceIterator;

    typedef std::vector<DartTraverser>::iterator BoundaryComponentsIterator;

public:
    template<class SrcIter, class SrcAcc>
    void init(SrcIter ul, SrcIter lr, SrcAcc src)
    {
        width_ = lr.x - ul.x;
        height_ = lr.y - ul.y;
        unsigned int totalwidth = width_ + 4;
        unsigned int totalheight = height_ + 4;

        nodeCount_ = edgeCount_ = faceCount_ = 0;

        cellImage.resize(totalwidth, totalheight);
        cellImage = CellPixel(CellTypeRegion, 0);

        cells = cellImage.upperLeft() + Diff2D(2,2);

        // extract contours in input image and put frame around them
        BImage contourImage(totalwidth, totalheight);
        initFourEightSegmentationContourImage(ul, lr, src, contourImage);
        initCellImage(contourImage);

        std::cerr << "FourEightSegmentation::label0Cells()\n";
        CellPixel::LabelType maxNodeLabel = label0Cells();

        std::cerr << "FourEightSegmentation::label1Cells(maxNodeLabel= "
                  << maxNodeLabel << ")\n";
        CellPixel::LabelType maxEdgeLabel = label1Cells(maxNodeLabel);

        std::cerr << "FourEightSegmentation::label2Cells()\n";
        CellPixel::LabelType maxFaceLabel = label2Cells(contourImage);

        std::cerr << "FourEightSegmentation::labelCircles(maxNodeLabel= "
                  << maxNodeLabel << ", maxEdgeLabel= " << maxEdgeLabel << ")\n";
        labelCircles(maxNodeLabel, maxEdgeLabel);

        exportImage(srcImageRange(cellImage,
                                  CellImageTypeAccessor<unsigned char>()),
                    ImageExportInfo("test.png"));
        std::cerr << "** test.png saved. **\n";

        std::cerr << "FourEightSegmentation::initNodeList(maxNodeLabel= "
                  << maxNodeLabel << ")\n";
        initNodeList(maxNodeLabel);
        std::cerr << "FourEightSegmentation::initEdgeList(maxEdgeLabel= "
                  << maxEdgeLabel << ")\n";
        initEdgeList(maxEdgeLabel);
        std::cerr << "FourEightSegmentation::initFaceList(maxFaceLabel= "
                  << maxFaceLabel << ")\n";
        initFaceList(contourImage, maxFaceLabel);
    }

    template<class SrcIter, class SrcAcc>
    void init(triple<SrcIter, SrcIter, SrcAcc> src)
    {
        init(src.first, src.second, src.third);
    }

    unsigned int width() const { return width_; }
    unsigned int height() const { return height_; }

    // the fooCount()s tell how many fooList elements are initialized()
    unsigned int nodeCount() const { return nodeCount_; }
    unsigned int edgeCount() const { return edgeCount_; }
    unsigned int faceCount() const { return faceCount_; }

    NodeIterator nodesBegin()
        { return NodeIterator(nodeList.begin(), nodeList.end()); }
    NodeIterator nodesEnd()
        { return NodeIterator(nodeList.end(), nodeList.end()); }
    NodeIterator findNode(unsigned int node)
        { return NodeIterator(nodeList.begin() + node, nodeList.end()); }
    NodeInfo & node(unsigned int node)
        { return nodeList[node]; }

    EdgeIterator edgesBegin()
        { return EdgeIterator(edgeList.begin(), edgeList.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edgeList.end(), edgeList.end()); }
    EdgeIterator findEdge(unsigned int edge)
        { return EdgeIterator(edgeList.begin() + edge, edgeList.end()); }
    EdgeInfo & edge(unsigned int edge)
        { return edgeList[edge]; }

    FaceIterator facesBegin()
        { return FaceIterator(faceList.begin(), faceList.end()); }
    FaceIterator facesEnd()
        { return FaceIterator(faceList.end(), faceList.end()); }
    FaceIterator findFace(unsigned int face)
        { return FaceIterator(faceList.begin() + face, faceList.end()); }
    FaceInfo & face(unsigned int face)
        { return faceList[face]; }

    CellImage cellImage;
    CellImage::traverser cells;

    template<class ImageIterator>
    inline LabelScanIterator<CellImage::traverser, ImageIterator>
    nodeScanIterator(int node, ImageIterator const &upperLeft)
    {
        return cellScanIterator(nodeList[node], CellTypeVertex, upperLeft);
    }

    template<class ImageIterator>
    inline LabelScanIterator<CellImage::traverser, ImageIterator>
    edgeScanIterator(int edge, ImageIterator const &upperLeft)
    {
        return cellScanIterator(edgeList[edge], CellTypeLine, upperLeft);
    }

    template<class ImageIterator>
    inline LabelScanIterator<CellImage::traverser, ImageIterator>
    faceScanIterator(int face, ImageIterator const &upperLeft)
    {
        return cellScanIterator(faceList[face], CellTypeRegion, upperLeft);
    }

    typedef LabelScanIterator<CellImage::traverser, CellImage::traverser>
    CellScanIterator;

    void mergeFaces(DartTraverser & dart)
    {
        EdgeInfo &edge= dart.edge();
        FaceInfo &face1= dart.leftFace();
        FaceInfo &face2= dart.rightFace();

        // relabel cells in cellImage:
        for(CellScanIterator it= edgeScanIterator(edge.label,
                                                  cells + edge.bounds.upperLeft());
            it.inRange(); ++it)
            *it= CellPixel(CellTypeRegion, face1.label);
        for(CellScanIterator it= faceScanIterator(face2.label,
                                                  cells + face2.bounds.upperLeft());
            it.inRange(); ++it)
            *it= CellPixel(CellTypeRegion, face1.label);
        
        // update bounds:
        face1.bounds |= edge.bounds;
        face1.bounds |= face2.bounds;

        // TODO: update contours, anchor

        edge.uninitialize();
        face2.uninitialize();
    }

private:
    unsigned int nodeCount_, edgeCount_, faceCount_;

    friend struct DartTraverser;

    NodeList nodeList;
    EdgeList edgeList;
    FaceList faceList;

    void initCellImage(BImage & contourImage);
    CellPixel::LabelType label0Cells();
    CellPixel::LabelType label1Cells(CellPixel::LabelType maxNodeLabel);
    CellPixel::LabelType label2Cells(BImage & contourImage);
    void labelCircles(CellPixel::LabelType & maxNodeLabel,
                      CellPixel::LabelType & maxEdgeLabel);

    void labelEdge(CellImageEightCirculator rayAtStart,
                   CellPixel::LabelType newLabel);

    void initNodeList(CellPixel::LabelType maxNodeLabel);
    void initEdgeList(CellPixel::LabelType maxEdgeLabel);
    void initFaceList(BImage & contourImage, CellPixel::LabelType maxFaceLabel);
    void initBoundingBoxes(CellPixel::LabelType maxNodeLabel,
                           CellPixel::LabelType maxEdgeLabel,
                           CellPixel::LabelType maxFaceLabel);

    template<class ImageIterator>
    inline LabelScanIterator<CellImage::traverser, ImageIterator>
    cellScanIterator(CellInfo cell, CellType cellType,
					 ImageIterator const &upperLeft);

private:
    unsigned int width_, height_;
};

// -------------------------------------------------------------------
//                    FourEightSegmentation functions
// -------------------------------------------------------------------
template<class ImageIterator>
LabelScanIterator<CellImage::traverser, ImageIterator>
FourEightSegmentation::cellScanIterator(
    CellInfo cell, CellType cellType, ImageIterator const &upperLeft)
{
    return LabelScanIterator<CellImage::traverser, ImageIterator>
        (cells + cell.bounds.upperLeft(), cells + cell.bounds.lowerRight(),
         CellPixel(cellType, cell.label),
         upperLeft + cell.bounds.upperLeft());
}

template<class SrcIter, class SrcAcc>
void initFourEightSegmentationContourImage(SrcIter ul, SrcIter lr, SrcAcc src,
                                           BImage & contourImage)
{
    int w = lr.x - ul.x;
    int h = lr.y - ul.y;
    int x,y;

    initImageBorder(destImageRange(contourImage),
                    1, 0);
    initImageBorder(srcIterRange(contourImage.upperLeft()+Diff2D(1,1),
                                 contourImage.lowerRight()-Diff2D(1,1),
                                 contourImage.accessor()),
                    1, 1);

    typedef typename SrcAcc::value_type SrcType;
    const SrcType zero = NumericTraits<SrcType>::zero();
    for(y=0; y<h; ++y, ++ul.y)
    {
        SrcIter sx = ul;
        for(x=0; x<w; ++x, ++sx.x)
        {
            if(src(sx) == zero)
                contourImage(x+2, y+2) = 1;
        }
    }
}

void FourEightSegmentation::initCellImage(BImage & contourImage)
{
    BImage::traverser rawLine = contourImage.upperLeft() + Diff2D(1,1);
    CellImage::traverser cellLine = cellImage.upperLeft() + Diff2D(1,1);
    for(int y=-1; y<=(int)height_; ++y, ++rawLine.y, ++cellLine.y)
    {
        BImage::traverser raw = rawLine;
        CellImage::traverser cell = cellLine;
        for(int x=-1; x<=(int)width_; ++x, ++raw.x, ++cell.x)
        {
            if(*raw == 0)
            {
                cell->setType(CellTypeRegion);
            }
            else
            {
                vigra::NeighborhoodCirculator<BImage::traverser, EightNeighborCode>
                    neighbors(raw, EightNeighborCode::SouthEast);
                vigra::NeighborhoodCirculator<BImage::traverser, EightNeighborCode>
                    end = neighbors;
                int conf = 0;
                do
                {
                    conf = (conf << 1) | *neighbors;
                }
                while(--neighbors != end);

                if(cellConfigurations[conf] == CellTypeError)
                {
                    char message[200];
                    sprintf(message, "FourEightSegmentation::init(): "
                            "Configuration at (%d, %d) must be thinned further",
                            x, y);

                    vigra_precondition(0, message);
                }

                cell->setType(cellConfigurations[conf]);
            }
        }
    }
}

// -------------------------------------------------------------------

CellPixel::LabelType FourEightSegmentation::label0Cells()
{
    BImage nodeImage(width_+4, height_+4);
    BImage::traverser nodes = nodeImage.upperLeft() + Diff2D(2,2);

    for(int y=-2; y<(int)height_+2; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-2, y);

        for(int x=-2; x<(int)width_+2; ++x, ++cell.x)
        {
            if(cell->type() == CellTypeVertex)
            {
                nodes(x,y) = 1;

                // test for forbidden configuration
                CellImageEightCirculator n(cell);
                CellImageEightCirculator nend = n;

                do
                {
                    if(n->type() == CellTypeLine && n[1].type() == CellTypeLine)
                    {
                        char msg[200];
                        sprintf(msg, "initFourEightSegmentation(): "
                                "Node at (%d, %d) has two incident edgels from the same edge (direction: %d)",
                                x, y, n - nend);
                        vigra_precondition(0, msg);
                    }
                }
                while(++n != nend);
            }
            else
            {
                nodes(x,y) = 0;
            }
        }
    }

    return labelImageWithBackground(
        srcImageRange(nodeImage),
        destImage(cellImage, CellImageLabelWriter<CellTypeVertex>()), true, 0);
}

// -------------------------------------------------------------------

CellPixel::LabelType FourEightSegmentation::label1Cells(
    CellPixel::LabelType maxNodeLabel)
{
    std::vector<bool> nodeProcessed(maxNodeLabel + 1, false);

    CellPixel::LabelType maxEdgeLabel = 0;

    for(int y=-1; y<=(int)height_; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, y);

        for(int x=-1; x<=(int)width_; ++x, ++cell.x)
        {
            if(cell->type() != CellTypeVertex)
                continue;
            if(nodeProcessed[cell->label()])
                continue;

            nodeProcessed[cell->label()] = true;

            DartTraverser rayAtStart(
                this, CellImageEightCirculator(cell, EightNeighborCode::West));
            DartTraverser rayEnd = rayAtStart;

            do
            {
                if(rayAtStart.edgeLabel() != 0)
                    continue;

                labelEdge(rayAtStart.neighborCirculator(), ++maxEdgeLabel);
            }
            while(rayAtStart.nextSigma() != rayEnd);
        }
    }

    return maxEdgeLabel;
}

// -------------------------------------------------------------------

CellPixel::LabelType FourEightSegmentation::label2Cells(BImage & contourImage)
{
    // labelImageWithBackground() starts with label 1, so don't
    // include outer border (infinite regions shall have label 0)
    return labelImageWithBackground(
        srcIterRange(contourImage.upperLeft() + Diff2D(1,1),
                     contourImage.lowerRight() - Diff2D(1,1),
                     contourImage.accessor()),
        destIter(cellImage.upperLeft() + Diff2D(1,1),
                 CellImageLabelWriter<CellTypeRegion>()),
        false, 1);
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelCircles(
    CellPixel::LabelType & maxNodeLabel, CellPixel::LabelType & maxEdgeLabel)
{
    for(int y=-1; y<=(int)height_; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, y);

        for(int x=-1; x<=(int)width_; ++x, ++cell.x)
        {
            if(cell->label() != 0)
                continue;

            // found a circle (not labeled by previous steps)

            // mark first point as node
            (*cell) = CellPixel(CellTypeVertex, ++maxNodeLabel);

            CellImageEightCirculator rayAtStart(cell);
            CellImageEightCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart->type() != CellTypeLine)
                    continue;
                if(rayAtStart->label() != 0)
                    continue;

                labelEdge(rayAtStart, ++maxEdgeLabel);
            }
            while(++rayAtStart != rayEnd);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelEdge(CellImageEightCirculator rayAtStart,
                                      CellPixel::LabelType newLabel)
{
    EdgelIterator edge(rayAtStart);

    // follow the edge and relabel it
    for(; !edge.atEnd(); ++edge)
    {
        edge->setLabel(newLabel, CellTypeLine);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initNodeList(CellPixel::LabelType maxNodeLabel)
{
    nodeList.resize(maxNodeLabel + 1);
    std::vector<int> crackCirculatedAreas(maxNodeLabel + 1, 0);

    for(Point2D pos= Point2D(-1, -1); pos.y<=(int)height_; ++pos.y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, pos.y);

        for(pos.x=-1; pos.x<=(int)width_; ++pos.x, ++cell.x)
        {
            if(cell->type() != CellTypeVertex)
                continue;

            CellPixel::LabelType index = cell->label();
            vigra_precondition(index < nodeList.size(),
                               "nodeList must be large enough!");

            if(!nodeList[index].initialized())
            {
                nodeList[index].label = index;
                ++nodeCount_;

                nodeList[index].centerX = pos.x;
                nodeList[index].centerY = pos.y;
                nodeList[index].size = 1;
                nodeList[index].anchor = DartTraverser(
                    this, CellImageEightCirculator(cell,
                                                   EightNeighborCode::West));

                // calculate area from following the outer contour of the node
                CrackContourCirculator<CellImage::traverser> crack(cell);
                CrackContourCirculator<CellImage::traverser> crackend(crack);
                do
                {
                    crackCirculatedAreas[index] += crack.diff().x * crack.pos().y -
                                                   crack.diff().y * crack.pos().x;
                }
                while(++crack != crackend);

                crackCirculatedAreas[index] /= 2;
            }
            else
            {
                nodeList[index].centerX += pos.x;
                nodeList[index].centerY += pos.y;

                // calculate area from counting the pixels of the node
                nodeList[index].size += 1;
            }
            nodeList[index].bounds |= pos;
        }
    }

    for(NodeIterator node= nodesBegin(); node.inRange(); ++node)
	{
		node->centerX /= node->size;
        node->centerY /= node->size;

        // methods to calculate the area must yield identical values
        if(crackCirculatedAreas[node->label] != node->size)
        {
            std::cerr << "FourEightSegmentation::initNodeList(): "
                      << "Node " << node->label << " at ("
                      << node->centerX << ", "
                      << node->centerY << ") has a hole:\n";
            Point2D center((int)node->centerX, (int)node->centerY);
            debugImage(crop(srcImageRange(cellImage),
                            Rect2D(center.x, center.y, center.x+5, center.y+5)),
                       std::cerr, 4);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initEdgeList(CellPixel::LabelType maxEdgeLabel)
{
    edgeList.resize(maxEdgeLabel + 1);

    NodeIterator n(nodeList.begin(), nodeList.end());

    for(; n.inRange(); ++n)
    {
        DartTraverser dart = n->anchor;
        if(dart.isSingular())
            continue;
        DartTraverser dartEnd = dart;

        do
        {
            CellPixel::LabelType label = dart.edgeLabel();
            vigra_precondition(label < edgeList.size(),
                               "edgeList must be large enough!");
            if(!edgeList[label].initialized())
            {
                edgeList[label].label = label;
                ++edgeCount_;
                edgeList[label].start = dart;
                edgeList[label].end = dart;
                edgeList[label].end.nextAlpha();
            }
        }
        while(dart.nextSigma() != dartEnd);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initFaceList(
    BImage & contourImage, CellPixel::LabelType maxFaceLabel)
{
    faceList.resize(maxFaceLabel + 1);

    IImage contourLabelImage(width_ + 4, height_ + 4);
    contourLabelImage = 0;
    int contourComponentsCount =
        labelImageWithBackground(srcImageRange(contourImage),
                                 destImage(contourLabelImage), true, 0);
    IImage::traverser contourLabel =
        contourLabelImage.upperLeft() + Diff2D(2, 2);

    std::vector<bool> contourProcessed(contourComponentsCount + 1, false);

    // process outer face
    faceList[0].label= 0;
    ++faceCount_;
    faceList[0].bounds |= Point2D(-1, -1) + cellImage.size();
    //faceList[0].anchor = Diff2D(-2, -2); FIXME: are the contour-anchors enough?
    DartTraverser anchor(this, CellImageEightCirculator(cells + Diff2D(-1, -1),
                                                        EightNeighborCode::West));
    faceList[0].contours.push_back(anchor.prevSigma());
    contourProcessed[contourLabel(-1, -1)] = true;

    for(Point2D pos= Point2D(-1, -1); pos.y <= (int)height_; ++pos.y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, pos.y);
        CellImage::traverser leftNeighbor = cell - Diff2D(1, 0);

        for(pos.x= -1; pos.x <= (int)width_; ++pos.x, ++cell.x, ++leftNeighbor.x)
        {
            if(cell->type() != CellTypeRegion)
            {
                if(cell->type() == CellTypeLine)
                    edgeList[cell->label()].bounds |= pos;
                continue;
            }

            CellPixel::LabelType index = cell->label();
            vigra_precondition(index < faceList.size(),
                               "faceList must be large enough!");

            if(!faceList[index].initialized())
            {
                std::cerr << "found face " << index << " at "
                          << pos.x << "," << pos.y << "\n";
                faceList[index].label = index;
                ++faceCount_;
                //faceList[index].anchor = pos;
                faceList[index].bounds |= pos;

                // find incident node
                if(leftNeighbor->type() == CellTypeVertex)
                {
                    vigra_precondition(leftNeighbor->type() == CellTypeVertex,
                                       "leftNeighbor expected to be a vertex");

                    DartTraverser anchor(this,
                                         CellImageEightCirculator(leftNeighbor));
                    anchor.prevSigma();

                    vigra_invariant(anchor.leftFaceLabel() == index,
                                    "FourEightSegmentation::initFaceList()");

                    faceList[index].contours.push_back(anchor);
                }
                else
                {
                    vigra_precondition(leftNeighbor->type() == CellTypeLine,
                                       "leftNeighbor expected to be an edge");

                    CellPixel::LabelType edgeIndex = leftNeighbor->label();

                    vigra_precondition(edgeList[edgeIndex].initialized(),
                                       "EdgeInfo expected to be initialized");

                    DartTraverser anchor = edgeList[edgeIndex].start;
                    if(anchor.leftFaceLabel() != index)
                        anchor.nextAlpha();

                    vigra_invariant(anchor.leftFaceLabel() == index,
                                    "FourEightSegmentation::initFaceList()");

                    faceList[index].contours.push_back(anchor);
                }
            }
            else
            {
                faceList[index].bounds |= pos;

                // look for inner contours
                CellImageEightCirculator neighbor(cell);
                CellImageEightCirculator nend = neighbor;

                do
                {
                    int boundaryIndex = contourLabel[neighbor.base() - cells];
                    if(boundaryIndex == 0 || contourProcessed[boundaryIndex])
                        continue;

                    // found an inner contour
                    contourProcessed[boundaryIndex] = true;

                    // find incident node
                    if(neighbor->type() == CellTypeVertex)
                    {
                        // this is the node
                        CellImageEightCirculator n = neighbor;
                        n.swapCenterNeighbor();

                        DartTraverser anchor(this, n);
                        anchor.prevSigma();

                        vigra_invariant(anchor.leftFaceLabel() == index,
                                        "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(anchor);
                    }
                    else
                    {
                        vigra_precondition(neighbor->type() == CellTypeLine,
                                           "neighbor expected to be an edge");

                        CellPixel::LabelType edgeIndex = neighbor->label();

                        vigra_precondition(edgeList[edgeIndex].initialized(),
                                           "EdgeInfo should be initialized");

                        DartTraverser anchor = edgeList[edgeIndex].start;
                        if(anchor.leftFaceLabel() != index)
                            anchor.nextAlpha();

                        vigra_invariant(anchor.leftFaceLabel() == index,
                                        "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(anchor);
                    }
                }
                while(++neighbor != nend);
            }
        }
    }
}

} // namespace CellImage

} // namespace vigra

#endif /* VIGRA_FOUREIGHTSEGMENTATION_HXX */
