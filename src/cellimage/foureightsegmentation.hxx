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
#include <vigra/rect2d.hxx>

#include "filteriterator.hxx"
#include "cellimage.hxx"

#include <functional>

namespace vigra {

namespace cellimage {

template<class Array, class Result = typename Array::value_type,
         class Index = unsigned int>
struct LookupFunctor : public std::unary_function<Index, Result>
{
    const Array &array_;

    LookupFunctor(const Array &array)
    : array_(array)
    {}

    Result operator()(Index i) const
    { return array_[i]; }
};

template<class Operation1, class Operation2>
class ComposeFunctor
: public std::unary_function<typename Operation2::argument_type,
                             typename Operation1::result_type>
{
  protected:
    Operation1 op1;
    Operation2 op2;
  public:
    ComposeFunctor(const Operation1& x = Operation1(),
                   const Operation2& y = Operation2())
    : op1(x), op2(y) {}

    typename Operation1::result_type
    operator()(const typename Operation2::argument_type& x) const
    {
        return op1(op2(x));
    }
};

// -------------------------------------------------------------------
//                            LabelScanIterator
// -------------------------------------------------------------------
template<class LabelTraverser, class SrcTraverser = Diff2D>
class LabelScanIterator
{
    LabelTraverser cellLR_, cellIter_;
    typename LabelTraverser::value_type cellPixelValue_;
    SrcTraverser imageIter_;
    unsigned int width_;

public:
        /** the iterator's value type
        */
    typedef typename SrcTraverser::value_type value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename SrcTraverser::reference reference;

        /** the iterator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename SrcTraverser::pointer pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    LabelScanIterator()
    {}

    LabelScanIterator(LabelTraverser cellUL, LabelTraverser cellLR,
                      typename LabelTraverser::value_type cellPixelValue,
                      SrcTraverser imageIter = SrcTraverser())
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
                {
                    if(*cellIter_ != cellPixelValue_)
                        operator++();
                }
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

            recheckSingularity();
        }

            /**
             * Makes calls to isSingular() return a sane value again after
             * the underlying cell complex has changed. Returns new, valid
             * value of isSingular.
             *
             * If (still) pointing to a CellTypeLine pixel, returns with
             * !isSingular(). Else, performs operation like nextSigma()
             * (returning eventually pointing to a new edge with
             * !isSingular()) but stops with (isSingular()==true) if it
             * finds no adjacent CellTypeLine pixel.
             */
        bool recheckSingularity()
        {
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

            return isSingular_;
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

        CellLabel startNodeLabel() const
        {
            return neighborCirc_.center()->label();
        }

        CellLabel endNodeLabel() const
        {
            DartTraverser reverseDart_(*this);
            reverseDart_.nextAlpha();
            return reverseDart_.startNodeLabel();
        }

        CellLabel edgeLabel() const
        {
            return neighborCirc_->label();
        }

        CellLabel leftFaceLabel() const
        {
            vigra_precondition(neighborCirc_[1].type() == CellTypeRegion,
                "insufficient algorithm for DartTraverser::rightFace()");
            return neighborCirc_[1].label();
        }

        CellLabel rightFaceLabel() const
        {
            vigra_precondition(neighborCirc_[-1].type() == CellTypeRegion,
                "insufficient algorithm for DartTraverser::rightFace()");
            return neighborCirc_[-1].label();
        }

        NodeInfo & startNode() const
        {
            return segmentation_->nodeList_[startNodeLabel()];
        }

        NodeInfo & endNode() const
        {
            return segmentation_->nodeList_[endNodeLabel()];
        }

        EdgeInfo & edge() const
        {
            return segmentation_->edgeList_[edgeLabel()];
        }

        FaceInfo & leftFace() const
        {
            return segmentation_->faceList_[leftFaceLabel()];
        }

        FaceInfo & rightFace() const
        {
            return segmentation_->faceList_[rightFaceLabel()];
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

        void reparent(FourEightSegmentation * segmentation)
        {
            neighborCirc_ = CellImageEightCirculator(
                segmentation->cells +
                (neighborCirc_.center() - segmentation_->cells),
                neighborCirc_.direction());
            segmentation_ = segmentation;
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

    typedef std::vector<DartTraverser>
        ContourComponents;
    typedef std::vector<DartTraverser>::iterator
        ContourComponentsIterator;
    typedef std::vector<DartTraverser>::const_iterator
        ConstContourComponentsIterator;

    struct CellInfo
    {
        CellLabel label;
        Rect2D bounds;
        unsigned int size;

        CellInfo() : label(NumericTraits<CellLabel>::max()) {}

        bool initialized() const
        { return label != NumericTraits<CellLabel>::max(); }

        void uninitialize()
        { label = NumericTraits<CellLabel>::max(); }
    };

    struct NodeInfo : CellInfo
    {
        DartTraverser anchor;
        // remove following members?
        float centerX, centerY;
    };

    struct EdgeInfo : CellInfo
    {
        DartTraverser start, end;
    };

    struct FaceInfo : CellInfo
    {
        ContourComponents contours;
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

    typedef FilterIterator<NodeList::iterator, FilterInitialized>
        NodeIterator;
    typedef FilterIterator<EdgeList::iterator, FilterInitialized>
        EdgeIterator;
    typedef FilterIterator<FaceList::iterator, FilterInitialized>
        FaceIterator;
    typedef FilterIterator<NodeList::const_iterator, FilterInitialized>
        ConstNodeIterator;
    typedef FilterIterator<EdgeList::const_iterator, FilterInitialized>
        ConstEdgeIterator;
    typedef FilterIterator<FaceList::const_iterator, FilterInitialized>
        ConstFaceIterator;

  public:
    FourEightSegmentation(const FourEightSegmentation &other)
    {
        deepCopy(other);
    }

    template<class SrcIter, class SrcAcc>
    FourEightSegmentation(SrcIter ul, SrcIter lr, SrcAcc src)
    {
        init(ul, lr, src);
    }

    template<class SrcIter, class SrcAcc>
    FourEightSegmentation(triple<SrcIter, SrcIter, SrcAcc> src)
    {
        init(src.first, src.second, src.third);
    }

  protected:
    template<class SrcIter, class SrcAcc>
    void init(SrcIter ul, SrcIter lr, SrcAcc src)
    {
        nodeCount_ = edgeCount_ = faceCount_ = 0;

        cellImage.resize(lr.x - ul.x + 4, lr.y - ul.y + 4);
        cellImage = CellPixel(CellTypeRegion, 0);

        cells = cellImage.upperLeft() + Diff2D(2, 2);

        // extract contours in input image and put frame around them
        BImage contourImage(cellImage.size());
        initFourEightSegmentationContourImage(ul, lr, src, contourImage);
        initCellImage(contourImage);

        CellLabel maxNodeLabel = label0Cells();
        CellLabel maxEdgeLabel = label1Cells(maxNodeLabel);
        CellLabel maxFaceLabel = label2Cells(contourImage);
        labelCircles(maxNodeLabel, maxEdgeLabel);
        initNodeList(maxNodeLabel);
        initEdgeList(maxEdgeLabel);
        initFaceList(contourImage, maxFaceLabel);
    }

    template<class SrcIter, class SrcAcc>
    void init(triple<SrcIter, SrcIter, SrcAcc> src)
    {
        init(src.first, src.second, src.third);
    }

  public:
    FourEightSegmentation &operator=(const FourEightSegmentation &other)
    {
        return deepCopy(other);
    }

    // the fooCount()s tell how many fooList_ elements are initialized()
    unsigned int nodeCount() const { return nodeCount_; }
    CellLabel maxNodeLabel() const { return nodeList_.size(); }
    unsigned int edgeCount() const { return edgeCount_; }
    CellLabel maxEdgeLabel() const { return edgeList_.size(); }
    unsigned int faceCount() const { return faceCount_; }
    CellLabel maxFaceLabel() const { return faceList_.size(); }

    NodeIterator nodesBegin()
        { return NodeIterator(nodeList_.begin(), nodeList_.end()); }
    NodeIterator nodesEnd()
        { return NodeIterator(nodeList_.end(), nodeList_.end()); }
    NodeIterator findNode(unsigned int node)
        { return NodeIterator(nodeList_.begin() + node, nodeList_.end()); }
    NodeInfo & node(unsigned int node)
        { return nodeList_[node]; }

    ConstNodeIterator nodesBegin() const
        { return ConstNodeIterator(nodeList_.begin(), nodeList_.end()); }
    ConstNodeIterator nodesEnd() const
        { return ConstNodeIterator(nodeList_.end(), nodeList_.end()); }
    ConstNodeIterator findNode(unsigned int node) const
        { return ConstNodeIterator(nodeList_.begin() + node, nodeList_.end()); }
    const NodeInfo & node(unsigned int node) const
        { return nodeList_[node]; }

    EdgeIterator edgesBegin()
        { return EdgeIterator(edgeList_.begin(), edgeList_.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edgeList_.end(), edgeList_.end()); }
    EdgeIterator findEdge(unsigned int edge)
        { return EdgeIterator(edgeList_.begin() + edge, edgeList_.end()); }
    EdgeInfo & edge(unsigned int edge)
        { return edgeList_[edge]; }

    ConstEdgeIterator edgesBegin() const
        { return ConstEdgeIterator(edgeList_.begin(), edgeList_.end()); }
    ConstEdgeIterator edgesEnd() const
        { return ConstEdgeIterator(edgeList_.end(), edgeList_.end()); }
    ConstEdgeIterator findEdge(unsigned int edge) const
        { return ConstEdgeIterator(edgeList_.begin() + edge, edgeList_.end()); }
    const EdgeInfo & edge(unsigned int edge) const
        { return edgeList_[edge]; }

    FaceIterator facesBegin()
        { return FaceIterator(faceList_.begin(), faceList_.end()); }
    FaceIterator facesEnd()
        { return FaceIterator(faceList_.end(), faceList_.end()); }
    FaceIterator findFace(unsigned int face)
        { return FaceIterator(faceList_.begin() + face, faceList_.end()); }
    FaceInfo & face(unsigned int face)
        { return faceList_[face]; }

    ConstFaceIterator facesBegin() const
        { return ConstFaceIterator(faceList_.begin(), faceList_.end()); }
    ConstFaceIterator facesEnd() const
        { return ConstFaceIterator(faceList_.end(), faceList_.end()); }
    ConstFaceIterator findFace(unsigned int face) const
        { return ConstFaceIterator(faceList_.begin() + face, faceList_.end()); }
    const FaceInfo & face(unsigned int face) const
        { return faceList_[face]; }

    CellImage cellImage;
    CellImage::traverser cells;

    template<class SrcTraverser>
    inline LabelScanIterator<CellImage::traverser, SrcTraverser>
    nodeScanIterator(int node, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(nodeList_[node], CellTypeVertex, upperLeft,
                                cropToBaseImage);
    }

    template<class SrcTraverser>
    inline LabelScanIterator<CellImage::traverser, SrcTraverser>
    edgeScanIterator(int edge, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(edgeList_[edge], CellTypeLine, upperLeft,
                                cropToBaseImage);
    }

    template<class SrcTraverser>
    inline LabelScanIterator<CellImage::traverser, SrcTraverser>
    faceScanIterator(int face, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(faceList_[face], CellTypeRegion, upperLeft,
                                cropToBaseImage);
    }

    typedef LabelScanIterator<CellImage::traverser, CellImage::traverser>
        CellScanIterator;

  protected:
    unsigned char findContourComponent(const ContourComponents &contours,
                                       const DartTraverser & dart);

    void removeNodeFromContours(ContourComponents &contours,
                                CellLabel nodeLabel);

  public:
    FaceInfo &removeIsolatedNode(const DartTraverser & dart);

    FaceInfo &mergeFaces(const DartTraverser & dart);

    FaceInfo &removeBridge(const DartTraverser & dart);

    EdgeInfo &mergeEdges(const DartTraverser & dart);

  private:
    FourEightSegmentation &deepCopy(const FourEightSegmentation &other)
    {
        cellImage = other.cellImage;
        cells = cellImage.upperLeft() + Diff2D(2,2);

        nodeList_ = other.nodeList_;
        nodeCount_ = other.nodeCount_;
        edgeList_ = other.edgeList_;
        edgeCount_ = other.edgeCount_;
        faceList_ = other.faceList_;
        faceCount_ = other.faceCount_;

        for(NodeIterator it= nodesBegin(); it.inRange(); ++it)
        {
            it->anchor.reparent(this);
        }
        
        for(EdgeIterator it= edgesBegin(); it.inRange(); ++it)
        {
            it->start.reparent(this);
            it->end.reparent(this);
        }

        for(FaceIterator it= facesBegin(); it.inRange(); ++it)
        {
            for(ContourComponentsIterator contour= it->contours.begin();
                contour != it->contours.end(); ++contour)
            {
                contour->reparent(this);
            }
        }

        return *this;
    }

    unsigned int nodeCount_, edgeCount_, faceCount_;

    NodeList nodeList_;
    EdgeList edgeList_;
    FaceList faceList_;

    friend struct DartTraverser;

    void initCellImage(BImage & contourImage);
    CellLabel label0Cells();
    CellLabel label1Cells(CellLabel maxNodeLabel);
    CellLabel label2Cells(BImage & contourImage);
    void labelCircles(CellLabel & maxNodeLabel,
                      CellLabel & maxEdgeLabel);

    void labelEdge(CellImageEightCirculator rayAtStart,
                   CellLabel newLabel);

    void initNodeList(CellLabel maxNodeLabel);
    void initEdgeList(CellLabel maxEdgeLabel);
    void initFaceList(BImage & contourImage, CellLabel maxFaceLabel);
    void initBoundingBoxes(CellLabel maxNodeLabel,
                           CellLabel maxEdgeLabel,
                           CellLabel maxFaceLabel);

    template<class SrcTraverser>
    inline LabelScanIterator<CellImage::traverser, SrcTraverser>
    cellScanIterator(CellInfo cell, CellType cellType,
					 SrcTraverser const &upperLeft,
                     bool cropToBaseImage) const;
};

// -------------------------------------------------------------------
//                    FourEightSegmentation functions
// -------------------------------------------------------------------
template<class SrcTraverser>
LabelScanIterator<CellImage::traverser, SrcTraverser>
FourEightSegmentation::cellScanIterator(
    CellInfo cell, CellType cellType, SrcTraverser const &upperLeft,
    bool cropToBaseImage) const
{
    //std::cerr << "cellScanIterator for " << CellPixel(cellType, cell.label)
    //          << " begins at " << cell.bounds.upperLeft() << std::endl;
    /*Rect2D cellBounds = cropToBaseImage ?
                        cell.bounds & Rect2D(Diff2D(0, 0),
											 cellImage.size() - Diff2D(4, 4)) :
                                             cell.bounds;*/
    return LabelScanIterator<CellImage::traverser, SrcTraverser>
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

} // namespace cellimage

} // namespace vigra

#endif /* VIGRA_FOUREIGHTSEGMENTATION_HXX */
