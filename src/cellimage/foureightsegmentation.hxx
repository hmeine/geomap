#ifndef VIGRA_GEOMAP_HXX
#define VIGRA_GEOMAP_HXX

#include <vigra/error.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <vigra/contourcirculator.hxx>
//#include "circulatoradaptor.hxx"

#include "filteriterator.hxx"
#include "cellimage.hxx"

#include <functional>

namespace vigra {

namespace cellimage {

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
        if(!((cellLR.x > cellUL.x) && (cellLR.y > cellUL.y)))
            cellIter_ = cellLR_;
        else
            if(cellIter_ != cellLR_ && *cellIter_ != cellPixelValue_)
                operator++();
    }

    LabelScanIterator &operator++()
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

    /*explicit*/ operator SrcTraverser()
    {
        return imageIter_;
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
    EdgelIterator(CellImageEightCirculator const &n)
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

    EdgelIterator &operator++()
    {
        neighborCirc_.moveCenterToNeighbor();
        neighborCirc_.turnRight();

        while(1)
        {
            if(neighborCirc_->type() != CellTypeRegion)
                break;
            ++neighborCirc_;
        }

        if(neighborCirc_.isDiagonal() &&
           (neighborCirc_[1].type() != CellTypeRegion))
            ++neighborCirc_;

        atEnd_ = (neighborCirc_->type() == CellTypeVertex);

        return *this;
    }

    EdgelIterator &jumpToOpposite()
    {
        atEnd_ = false;
        while(!atEnd())
            operator++();
        neighborCirc_.swapCenterNeighbor();
        return *this;
    }

    const CellImageEightCirculator &neighborCirculator() const
    {
        return neighborCirc_;
    }

    const CellImageEightCirculator::base_type &base() const
    {
        return neighborCirc_.base();
    }
};

template <class IMAGEITERATOR>
class CrackEdgeIterator : public CrackContourCirculator<IMAGEITERATOR>
{
    bool atEnd_;

  public:
    typedef CrackContourCirculator<IMAGEITERATOR> Base;

        /** the circulator tag
        */
        //    typedef forward_iterator_tag iterator_category;

    typedef typename Base::NEIGHBORHOODCIRCULATOR NeighborhoodCirculator;

    CrackEdgeIterator(NeighborhoodCirculator const & circ)
    : Base(circ), atEnd_(false)
    {}

    bool atEnd() const
    {
        return atEnd_;
    }

    bool inRange() const
    {
        return !atEnd_;
    }

    EdgelIterator &operator++()
    {
        unsigned char edgeCount = 1;
        typename NeighborhoodCirculator::value_type label(label_);
        if(*neighborCirc_ != label)
        {
            label = *neighborCirc_;
            ++edgeCount;
        }
        if(neighborCirc_[1] != label)
        {
            label = neighborCirc_[1];
            ++edgeCount;
        }
        if(neighborCirc_[2] != label)
        {
            label = neighborCirc_[2];
            ++edgeCount;
        }
        atEnd_ = (edgeCount > 2);

        if(!atEnd_)
            Base::operator++();
        else
        {
            // turn around:
            neighborCirc_.turnRight();
            neighborCirc_.moveCenterToNeighbor();
            neighborCirc_.turnRight();
        }

        return *this;
    }

    EdgelIterator &jumpToOpposite()
    {
        atEnd_ = false;
        while(!atEnd())
            operator++();
        neighborCirc_.swapCenterNeighbor();
        return *this;
    }
};

/********************************************************************/
/*                                                                  */
/*                              GeoMap                              */
/*                                                                  */
/********************************************************************/
class GeoMap
{
public:
    struct NodeInfo;
    struct EdgeInfo;
    struct FaceInfo;

    /********************************************************************/
    /*                                                                  */
    /*                      GeoMap::DartTraverser                       */
    /*                                                                  */
    /********************************************************************/
    class DartTraverser
    {
        CellImageEightCirculator neighborCirc_;

        GeoMap *segmentation_;

    public:
        DartTraverser() : segmentation_(0L) {}

        DartTraverser(GeoMap *segmentation,
                      CellImageEightCirculator const &circ)
            : neighborCirc_(circ),
              segmentation_(segmentation)
        {
            vigra_precondition(neighborCirc_.center()->type() == CellTypeVertex,
            "GeoMap::DartTraverser(): center is not a node");

            vigra_precondition(neighborCirc_->type() != CellTypeVertex,
            "GeoMap::DartTraverser(): neighbor is a node");

            while(badDiagonalConfig())
            {
                ++neighborCirc_;
                // pointing from vertex to vertex pixel now
                neighborCirc_.swapCenterNeighbor();
                // still pointing from vertex to vertex pixel
                ++neighborCirc_;
                // finally, pointing from vertex to line pixel
            }

            if(neighborCirc_->type() != CellTypeLine)
                carefulNextSigma();
        }

    private:
            // the "careful" sigma operations guarantee that at least
            // one step is made (they work even if isSingular()), and
            // do not loop infinitely in case no edge is attached
            // (i.e. if the dart is indeed singular)        
        DartTraverser &carefulNextSigma() throw ()
        {
            CellImageEightCirculator nend = neighborCirc_;
            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                    ++neighborCirc_;
                }
                tryNextSigma();
                if(neighborCirc_ == nend)
                {
                    // did not find any adjacent line pixel
                    break;
                }
            }
            return *this;
        }

        DartTraverser &carefulPrevSigma() throw ()
        {
            CellImageEightCirculator nend = neighborCirc_;
            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                    --neighborCirc_;
                }
                tryPrevSigma();
                if(neighborCirc_ == nend)
                {
                    // did not find any adjacent line pixel
                    break;
                }
            }
            return *this;
        }

        friend class GeoMap;

    public:
        DartTraverser &nextAlpha() throw ()
        {
            if(isSingular())
                return *this;

            if(segmentation_->initialized())
            {
                EdgeInfo &cEdge(edge());
                return operator=(
                    operator==(cEdge.start) ? cEdge.end : cEdge.start);
            }

            EdgelIterator line(neighborCirc_);
            line.jumpToOpposite();

            neighborCirc_ = line.neighborCirculator();

            return *this;
        }

        DartTraverser &prevAlpha() throw ()
        {
            return nextAlpha();
        }

        DartTraverser &nextPhi() throw ()
        {
            return nextAlpha().prevSigma();
        }

        DartTraverser &prevPhi() throw ()
        {
            return nextSigma().prevAlpha();
        }

        DartTraverser &nextSigma() throw ()
        {
            if(isSingular())
                return *this;

            tryNextSigma();

            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                    ++neighborCirc_;
                }
                tryNextSigma();
            }
            return *this;
        }

        DartTraverser &prevSigma() throw ()
        {
            if(isSingular())
                return *this;

            tryPrevSigma();

            while(neighborCirc_->type() != CellTypeLine)
            {
                if(neighborCirc_->type() == CellTypeVertex)
                {
                    neighborCirc_.swapCenterNeighbor();
                    --neighborCirc_;
                }
                tryPrevSigma();
            }
            return *this;
        }

        bool isSingular() const throw ()
        {
            return neighborCirc_->type() != CellTypeLine;
        }

        CellLabel startNodeLabel() const throw ()
        {
            return neighborCirc_.center()->label();
        }

        CellLabel endNodeLabel() const throw ()
        {
            DartTraverser reverseDart_(*this);
            reverseDart_.nextAlpha();
            return reverseDart_.startNodeLabel();
        }

        CellLabel edgeLabel() const throw ()
        {
            return neighborCirc_->label();
        }

        CellLabel leftFaceLabel() const throw ()
        {
//             vigra_precondition(neighborCirc_[1].type() == CellTypeRegion,
//                 "insufficient algorithm for DartTraverser::leftFaceLabel()");
            return neighborCirc_[1].label();
        }

        CellLabel rightFaceLabel() const throw ()
        {
//             vigra_precondition(neighborCirc_[-1].type() == CellTypeRegion,
//                 "insufficient algorithm for DartTraverser::rightFaceLabel()");
            return neighborCirc_[-1].label();
        }

        NodeInfo &startNode() const throw ()
        {
            return segmentation_->nodeList_[startNodeLabel()];
        }

        NodeInfo &endNode() const throw ()
        {
            return segmentation_->nodeList_[endNodeLabel()];
        }

        EdgeInfo &edge() const throw ()
        {
            return segmentation_->edgeList_[edgeLabel()];
        }

        FaceInfo &leftFace() const throw ()
        {
            return segmentation_->faceList_[leftFaceLabel()];
        }

        FaceInfo &rightFace() const throw ()
        {
            return segmentation_->faceList_[rightFaceLabel()];
        }

        bool operator==(DartTraverser const &o) const throw ()
        {
            // FIXME: think about singular darts - I would like to
            // save this comparison for efficiency, however I don't
            // know where to turn singular DTs to a canonical
            // direction nor do I know a good name for another
            // comparison function.
            // ATM we don't need correct equality checks for singular
            // darts.

            if(!isSingular())
                return neighborCirc_ == o.neighborCirc_;
            else
                return startNodeLabel() == o.startNodeLabel();
        }

        bool operator!=(DartTraverser const &o) const throw ()
        {
            // see op==
            return neighborCirc_ != o.neighborCirc_;
        }

        const CellImageEightCirculator &neighborCirculator() const throw ()
        {
            return neighborCirc_;
        }

        GeoMap *segmentation() const throw ()
        {
            return segmentation_;
        }

        void reparent(GeoMap *segmentation) throw ()
        {
            neighborCirc_ = CellImageEightCirculator(
                segmentation->cells +
                (neighborCirc_.center() - segmentation_->cells),
                neighborCirc_.direction());
            segmentation_ = segmentation;
        }

        typedef unsigned int Serialized;

        Serialized serialize() const throw ()
        {
            Diff2D cPos(neighborCirc_.center() -
                        segmentation_->cellImage.upperLeft());
            return (int)neighborCirc_.direction() << 28 | // needs 3 bits
                cPos.x << 14 | cPos.y;
        }

        DartTraverser(GeoMap *segmentation,
                      Serialized serialized)
        : neighborCirc_(segmentation->cellImage.upperLeft() +
                        Diff2D((serialized >> 14) & 0x3fff,
                               serialized & 0x3fff),
                        (CellImageEightCirculator::Direction)(serialized >> 28)),
          segmentation_(segmentation)
        {}

    private:
        void tryNextSigma() throw ()
        {
            ++neighborCirc_;

            if(badDiagonalConfig())
                ++neighborCirc_;
        }

        void tryPrevSigma() throw ()
        {
            --neighborCirc_;

            if(badDiagonalConfig())
                --neighborCirc_;
        }

            // prevent double stop at a line pixel from different source
            // vertex pixels
        bool badDiagonalConfig() const
        {
            return (neighborCirc_->type() == CellTypeLine &&
                    (neighborCirc_[1].type() == CellTypeVertex ||
                     neighborCirc_[-1].type() == CellTypeVertex));
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
        unsigned short degree;
    };

    struct EdgeInfo : CellInfo
    {
        DartTraverser start, end;

        bool isLoop() const
        { return start.startNodeLabel() == end.startNodeLabel(); }

        bool isBridge() const
        { return start.leftFaceLabel() == start.rightFaceLabel(); }
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
    GeoMap(const GeoMap &other)
    : initialized_(false)
    {
        deepCopy(other);
    }

    template<class SrcIter, class SrcAcc>
    GeoMap(SrcIter ul, SrcIter lr, SrcAcc src,
           typename SrcAcc::value_type boundaryValue,
           CellType cornerType)
    : initialized_(false)
    {
        init(ul, lr, src, boundaryValue, cornerType);
    }

    template<class SrcIter, class SrcAcc>
    GeoMap(triple<SrcIter, SrcIter, SrcAcc> src,
           typename SrcAcc::value_type boundaryValue,
           CellType cornerType)
    : initialized_(false)
    {
        init(src.first, src.second, src.third, boundaryValue, cornerType);
    }

    GeoMap(const CellImage &importImage);

  protected:
    template<class SrcIter, class SrcAcc>
    void init(SrcIter ul, SrcIter lr, SrcAcc src,
              typename SrcAcc::value_type boundaryValue,
              CellType cornerType)
    {
        // extract contours in input image and put frame around them
        BImage contourImage(lr.x - ul.x + 4, lr.y - ul.y + 4);
        initGeoMapContourImage(ul, lr, src, contourImage,
                               boundaryValue);

        cellImage.resize(contourImage.size());
        cells = cellImage.upperLeft() + Diff2D(2, 2);
        initCellImage(contourImage, cornerType);

        CellLabel maxNodeLabel = label0Cells();
        CellLabel maxEdgeLabel = label1Cells(maxNodeLabel);
        CellLabel maxFaceLabel = label2Cells(contourImage);
        labelCircles(maxNodeLabel, maxEdgeLabel);

        nodeCount_ = edgeCount_ = faceCount_ = 0;

        initNodeList(maxNodeLabel);
        initEdgeList(maxEdgeLabel);
        initFaceList(contourImage, maxFaceLabel);

        initialized_ = true;
    }

  public:
    GeoMap &operator=(const GeoMap &other)
    {
        return deepCopy(other);
    }

    bool initialized() const { return initialized_; }

        // the fooCount()s tell how many fooList_ elements are initialized()
    unsigned int nodeCount() const { return nodeCount_; }
    CellLabel maxNodeLabel() const { return nodeList_.size(); }
    unsigned int edgeCount() const { return edgeCount_; }
    CellLabel maxEdgeLabel() const { return edgeList_.size(); }
    unsigned int faceCount() const { return faceCount_; }
    CellLabel maxFaceLabel() const { return faceList_.size(); }

        // node access
    NodeIterator nodesBegin()
        { return NodeIterator(nodeList_.begin(), nodeList_.end()); }
    NodeIterator nodesEnd()
        { return NodeIterator(nodeList_.end(), nodeList_.end()); }
    NodeIterator findNode(CellLabel node)
        { return NodeIterator(nodeList_.begin() + node, nodeList_.end()); }
    NodeInfo &node(CellLabel node)
        { return nodeList_[node]; }

    ConstNodeIterator nodesBegin() const
        { return ConstNodeIterator(nodeList_.begin(), nodeList_.end()); }
    ConstNodeIterator nodesEnd() const
        { return ConstNodeIterator(nodeList_.end(), nodeList_.end()); }
    ConstNodeIterator findNode(CellLabel node) const
        { return ConstNodeIterator(nodeList_.begin() + node, nodeList_.end()); }
    const NodeInfo &node(CellLabel node) const
        { return nodeList_[node]; }

        // edge access
    EdgeIterator edgesBegin()
        { return EdgeIterator(edgeList_.begin(), edgeList_.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edgeList_.end(), edgeList_.end()); }
    EdgeIterator findEdge(CellLabel edge)
        { return EdgeIterator(edgeList_.begin() + edge, edgeList_.end()); }
    EdgeInfo &edge(CellLabel edge)
        { return edgeList_[edge]; }

    ConstEdgeIterator edgesBegin() const
        { return ConstEdgeIterator(edgeList_.begin(), edgeList_.end()); }
    ConstEdgeIterator edgesEnd() const
        { return ConstEdgeIterator(edgeList_.end(), edgeList_.end()); }
    ConstEdgeIterator findEdge(CellLabel edge) const
        { return ConstEdgeIterator(edgeList_.begin() + edge, edgeList_.end()); }
    const EdgeInfo &edge(CellLabel edge) const
        { return edgeList_[edge]; }

        // face access
    FaceIterator facesBegin()
        { return FaceIterator(faceList_.begin(), faceList_.end()); }
    FaceIterator facesEnd()
        { return FaceIterator(faceList_.end(), faceList_.end()); }
    FaceIterator findFace(CellLabel face)
        { return FaceIterator(faceList_.begin() + face, faceList_.end()); }
    FaceInfo &face(CellLabel face)
        { return faceList_[face]; }

    ConstFaceIterator facesBegin() const
        { return ConstFaceIterator(faceList_.begin(), faceList_.end()); }
    ConstFaceIterator facesEnd() const
        { return ConstFaceIterator(faceList_.end(), faceList_.end()); }
    ConstFaceIterator findFace(CellLabel face) const
        { return ConstFaceIterator(faceList_.begin() + face, faceList_.end()); }
    const FaceInfo &face(CellLabel face) const
        { return faceList_[face]; }

        // cellimage access
    CellImage cellImage;
    CellImage::traverser cells; // points to cellImage[2, 2] which is coord (0, 0)

    template<class ScanTargetTraverser>
    struct ScanIterator
    {
        typedef LabelScanIterator<
            CellImage::const_traverser, ScanTargetTraverser>
            type;
    };

    template<class SrcTraverser>
    inline typename ScanIterator<SrcTraverser>::type
    nodeScanIterator(int node, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(nodeList_[node], CellTypeVertex, upperLeft,
                                cropToBaseImage);
    }

    template<class SrcTraverser>
    inline typename ScanIterator<SrcTraverser>::type
    edgeScanIterator(int edge, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(edgeList_[edge], CellTypeLine, upperLeft,
                                cropToBaseImage);
    }

    template<class SrcTraverser>
    inline typename ScanIterator<SrcTraverser>::type
    faceScanIterator(int face, SrcTraverser const &upperLeft,
                     bool cropToBaseImage = true) const
    {
        return cellScanIterator(faceList_[face], CellTypeRegion, upperLeft,
                                cropToBaseImage);
    }

  protected:
    unsigned int findComponentAnchor(const FaceInfo &face,
                                     const DartTraverser &dart);

    void removeNodeFromContours(ContourComponents &contours,
                                CellLabel nodeLabel);

  public:
    FaceInfo &removeIsolatedNode(const DartTraverser &dart);

    FaceInfo &mergeFaces(const DartTraverser &dart);

    FaceInfo &removeBridge(const DartTraverser &dart);

    EdgeInfo &mergeEdges(const DartTraverser &dart);

  protected:
    void checkConsistency();

    GeoMap &deepCopy(const GeoMap &other);

    unsigned int nodeCount_, edgeCount_, faceCount_;

    NodeList nodeList_;
    EdgeList edgeList_;
    FaceList faceList_;

    bool initialized_;

    friend struct DartTraverser;

    void initCellImage(BImage &contourImage, CellType cornerType);
    CellLabel label0Cells();
    CellLabel label1Cells(CellLabel maxNodeLabel);
    CellLabel label2Cells(BImage &contourImage);
    void labelCircles(CellLabel &maxNodeLabel,
                      CellLabel &maxEdgeLabel);

    void labelEdge(CellImageEightCirculator rayAtStart,
                   CellLabel newLabel);

    void initNodeList(CellLabel maxNodeLabel);
    void initEdgeList(CellLabel maxEdgeLabel);
    void initFaceList(BImage &contourImage, CellLabel maxFaceLabel);
    void initBoundingBoxes(CellLabel maxNodeLabel,
                           CellLabel maxEdgeLabel,
                           CellLabel maxFaceLabel);

    typedef ScanIterator<CellImage::traverser>::type
        CellScanIterator;

    template<class SrcTraverser>
    inline typename ScanIterator<SrcTraverser>::type
    cellScanIterator(const CellInfo &cell, CellType cellType,
                     SrcTraverser const &upperLeft,
                     bool cropToBaseImage) const;
};

// -------------------------------------------------------------------
//                    GeoMap functions
// -------------------------------------------------------------------
template<class SrcTraverser>
inline typename GeoMap::ScanIterator<SrcTraverser>::type
GeoMap::cellScanIterator(
    const CellInfo &cell, CellType cellType, SrcTraverser const &upperLeft,
    bool cropToBaseImage) const
{
    Rect2D cellBounds = cropToBaseImage ?
                        cell.bounds & Rect2D(cellImage.size() - Diff2D(4, 4)) :
                        cell.bounds;
    return LabelScanIterator<CellImage::const_traverser, SrcTraverser>
        (cells + cellBounds.upperLeft(), cells + cellBounds.lowerRight(),
         CellPixel(cellType, cell.label),
         upperLeft + cellBounds.upperLeft());
}

template<class SrcIter, class SrcAcc>
void initGeoMapContourImage(SrcIter ul, SrcIter lr, SrcAcc src,
                            BImage &contourImage,
                            typename SrcAcc::value_type boundaryValue)
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

    for(y=0; y<h; ++y, ++ul.y)
    {
        SrcIter sx = ul;
        for(x=0; x<w; ++x, ++sx.x)
        {
            contourImage(x+2, y+2) = (src(sx) == boundaryValue ? 1 : 0);
        }
    }
}

void debugDart(const GeoMap::DartTraverser &dart);

} // namespace cellimage

} // namespace vigra

std::ostream &
operator<<(std::ostream & out,
           const vigra::cellimage::GeoMap::DartTraverser & d);

#endif /* VIGRA_GEOMAP_HXX */
