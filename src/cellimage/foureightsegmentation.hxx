#ifndef VIGRA_FOUREIGHTSEGMENTATION_HXX
#define VIGRA_FOUREIGHTSEGMENTATION_HXX

#include "vigra/error.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/labelimage.hxx"
//#include "circulatoradaptor.hxx"
#include "configurations.h"
#include "pixelneighborhood.hxx"

namespace vigra {

namespace FourEightSegmentation {

class FourEightSegmentation;

// -------------------------------------------------------------------
//                       OldNeighborhoodCirculator
// -------------------------------------------------------------------

    // sits on center(), points to location() (diff is diff()), can
    // turn/jumpToOpposite, gives cells/labels of center(), location()
    // and forward/backward neighbors
class OldNeighborhoodCirculator
{
protected:
    FourEightSegmentation *segmentation_;
    Diff2D center_;
    EightNeighborOffsetCirculator neighbor_;

public:
    // FIXME: where is this default constructor needed?
    // what should it (not) do?
    OldNeighborhoodCirculator()
        : segmentation_(0)
    {}

    OldNeighborhoodCirculator(FourEightSegmentation * segmentation,
                              Diff2D const & center,
                              EightNeighborOffsetCirculator::Direction dir = EightNeighborOffsetCirculator::East)
        : segmentation_(segmentation), center_(center), neighbor_(dir)
    {}

    // accessor functions for cell types and labels from
    // segmentation_
    inline unsigned char cell() const;
    inline unsigned char neighborCell() const;
    inline unsigned char forwardNeighborCell() const;
    inline unsigned char backwardNeighborCell() const;
    inline int label() const;
    inline int neighborLabel() const;
    inline int forwardNeighborLabel() const;
    inline int backwardNeighborLabel() const;

    OldNeighborhoodCirculator & operator++()
    {
        ++neighbor_;
        return *this;
    }
    OldNeighborhoodCirculator operator++(int)
    {
        OldNeighborhoodCirculator ret(*this);
        operator++();
        return ret;
    }
    OldNeighborhoodCirculator & operator+=(int d)
    {
        neighbor_ += d;
        return *this;
    }

    OldNeighborhoodCirculator & operator--()
    {
        --neighbor_;
        return *this;
    }
    OldNeighborhoodCirculator operator--(int)
    {
        OldNeighborhoodCirculator ret(*this);
        operator--();
        return ret;
    }
    OldNeighborhoodCirculator & operator-=(int d)
    {
        neighbor_ -= d;
        return *this;
    }

    OldNeighborhoodCirculator & jumpToOpposite()
    {
        center_ += diff();
        neighbor_.turnRound();

        return *this;
    }

    OldNeighborhoodCirculator & turnRight()
    {
        neighbor_.turnRight();
        return *this;
    }
    OldNeighborhoodCirculator & turnLeft()
    {
        neighbor_.turnLeft();
        return *this;
    }
    OldNeighborhoodCirculator & turnRound()
    {
        neighbor_.turnRound();
        return *this;
    }

    OldNeighborhoodCirculator & translate(Diff2D const & d)
    {
        center_ += d;
        return *this;
    }

    bool operator==(OldNeighborhoodCirculator const & o) const
    {
        return center_ == o.center_ && neighbor_ == o.neighbor_;
    }

    bool operator!=(OldNeighborhoodCirculator const & o) const
    {
        return center_ != o.center_ || neighbor_ != o.neighbor_;
    }

    Diff2D const & center() const
    {
        return center_;
    }

    Diff2D const & diff() const
    {
        return neighbor_.diff();
    }

    Diff2D location() const
    {
        return center_ + diff();
    }

    EightNeighborOffsetCirculator::Direction directionCode() const
    {
        return neighbor_.direction();
    }

    bool isDiagonal() const
    {
        return neighbor_.isDiagonal();
    }

    FourEightSegmentation * segmentation() const
    {
        return segmentation_;
    }
};

// -------------------------------------------------------------------
//                            CrackCirculator
// -------------------------------------------------------------------
class CrackCirculator
{
protected:
    OldNeighborhoodCirculator neighborCirc_;
    int label;
    Diff2D pos_;

public:
    CrackCirculator(OldNeighborhoodCirculator const & n)
        : neighborCirc_(n),
          label(n.cell()),
          pos_(0, 0)
    {
        neighborCirc_.turnLeft();
    }

    CrackCirculator & operator++()
    {
        pos_ += neighborCirc_.diff();

        neighborCirc_--;
        neighborCirc_.translate(neighborCirc_.diff());

        if(neighborCirc_.cell() == label)
        {
            --neighborCirc_;
        }
        else
        {
            neighborCirc_ += 3;
            neighborCirc_.translate(neighborCirc_.diff());
            if(neighborCirc_.cell() == label)
            {
                neighborCirc_.turnRight();
            }
            else
            {
                neighborCirc_ += 2;
                neighborCirc_.translate(neighborCirc_.diff());
                neighborCirc_.turnRight();
            }
        }

        return *this;
    }

    bool operator==(CrackCirculator const & o) const
    {
        return neighborCirc_ == o.neighborCirc_;
    }

    bool operator!=(CrackCirculator const & o) const
    {
        return neighborCirc_ != o.neighborCirc_;
    }

    Diff2D const & diff() const { return neighborCirc_.diff(); }

    Diff2D const & pos() const { return pos_; }
};

// -------------------------------------------------------------------
//                             EdgelIterator
// -------------------------------------------------------------------
class EdgelIterator
{
    OldNeighborhoodCirculator neighborCirc_;
    bool isEnd_;

public:
    EdgelIterator(OldNeighborhoodCirculator const & n)
        : neighborCirc_(n), isEnd_(false)
    {}

    Diff2D location() const
    {
        return neighborCirc_.location();
    }

    bool isEnd() const
    {
        return isEnd_;
    }

    EdgelIterator & operator++()
    {
        neighborCirc_.jumpToOpposite();
        neighborCirc_.turnLeft();

        while(1)
        {
            if(neighborCirc_.neighborCell() == CellConfigurationVertex)
            {
                isEnd_ = true;
                break;
            }
            if(neighborCirc_.neighborCell() == CellConfigurationLine)
            {
                break;
            }
            ++neighborCirc_;
        }

        if(neighborCirc_.isDiagonal() &&
           neighborCirc_.forwardNeighborCell() == CellConfigurationVertex)
        {
            ++neighborCirc_;
            isEnd_ = true;
        }

        return *this;
    }

    EdgelIterator & jumpToOpposite()
    {
        while(!isEnd())
            operator++();
        neighborCirc_.jumpToOpposite();
        return *this;
    }

    operator OldNeighborhoodCirculator()
    {
        return neighborCirc_;
    }
};

// -------------------------------------------------------------------
//                             RayCirculator
// -------------------------------------------------------------------
struct RayCirculator
{
private:
    OldNeighborhoodCirculator neighborCirc_;
    bool isSingular_;

public:
    // default constructor needed for NodeInfo/EdgeInfo, init() must not be called!
    RayCirculator() {}

    RayCirculator(FourEightSegmentation * seg, Diff2D const & loc,
                  EightNeighborOffsetCirculator::Direction dir = EightNeighborOffsetCirculator::East)
        : neighborCirc_(seg, loc, dir)
    {
        init();
    }

    RayCirculator(OldNeighborhoodCirculator const & n)
        : neighborCirc_(n)
    {
        init();
    }

    RayCirculator & operator++()
    {
        if(isSingular_) return *this;

        tryNext();

        while(neighborCirc_.neighborCell() != CellConfigurationLine)
        {
            if(neighborCirc_.neighborCell() == CellConfigurationVertex)
            {
                neighborCirc_.jumpToOpposite();
            }
            tryNext();
        }
        return *this;
    }

    RayCirculator operator++(int)
    {
        RayCirculator ret(*this);
        operator++();
        return ret;
    }

    RayCirculator & operator--()
    {
        if(isSingular_) return *this;

        tryPrev();

        while(neighborCirc_.neighborCell() != CellConfigurationLine)
        {
            if(neighborCirc_.neighborCell() == CellConfigurationVertex)
            {
                neighborCirc_.jumpToOpposite();
            }
            tryPrev();
        }
        return *this;
    }

    RayCirculator operator--(int)
    {
        RayCirculator ret(*this);
        operator--();
        return ret;
    }

    RayCirculator & jumpToOpposite()
    {
        if(isSingular_) return *this;

        EdgelIterator line(neighborCirc_);
        line.jumpToOpposite();

        neighborCirc_ = line;

        return *this;
    }

    bool operator==(RayCirculator const & o) const
    {
        return neighborCirc_ == o.neighborCirc_;
    }

    bool operator!=(RayCirculator const & o) const
    {
        return neighborCirc_ != o.neighborCirc_;
    }

    FourEightSegmentation * segmentation() const
    {
        return neighborCirc_.segmentation();
    }

    Diff2D const & center() const { return neighborCirc_.center(); }

    int nodeLabel() const { return neighborCirc_.label(); }
    int edgeLabel() const { return neighborCirc_.neighborLabel(); }
    int leftFaceLabel() const { return neighborCirc_.forwardNeighborLabel(); }
    int rightFaceLabel() const { return neighborCirc_.backwardNeighborLabel(); }

    inline int degree() const;
    inline float x() const;
    inline float y() const;

    const OldNeighborhoodCirculator &neighborhoodCirculator() const
    {
        return neighborCirc_;
    }

private:
    void init()
    {
        vigra_precondition(neighborCirc_.cell() == CellConfigurationVertex,
        "FourEightSegmentation::RayCirculator(): center is not a node");

        vigra_precondition(neighborCirc_.neighborCell() != CellConfigurationVertex,
        "FourEightSegmentation::RayCirculator(): neighbor is a node");

        OldNeighborhoodCirculator n = neighborCirc_;
        isSingular_ = true;
        do
        {
            if(n.neighborCell() != CellConfigurationRegion)
            {
                isSingular_ = false;
                break;
            }
        }
        while(++n != neighborCirc_);

        if(neighborCirc_.neighborCell() != CellConfigurationLine)
            operator++();
    }

    void tryNext()
    {
        ++neighborCirc_;

        if(badDiagonalConfig()) ++neighborCirc_;
    }

    void tryPrev()
    {
        --neighborCirc_;

        if(badDiagonalConfig()) --neighborCirc_;
    }

    // prevent double stop at a line pixel from different source
    // vertex pixels
    bool badDiagonalConfig()
    {
        return (neighborCirc_.neighborCell() == CellConfigurationLine &&
                (neighborCirc_.forwardNeighborCell() == CellConfigurationVertex ||
                 neighborCirc_.backwardNeighborCell() == CellConfigurationVertex));
    }
};

// -------------------------------------------------------------------
//                           ContourCirculator
// -------------------------------------------------------------------
struct ContourCirculator
{
    RayCirculator ray_;

    ContourCirculator(RayCirculator r)
        : ray_(r)
    {}

    ContourCirculator & operator++()
    {
        ray_.jumpToOpposite();
        --ray_;
        return *this;
    }

    ContourCirculator operator++(int)
    {
        ContourCirculator ret(*this);
        operator++();
        return ret;
    }

    ContourCirculator & operator--()
    {
        ++ray_;
        ray_.jumpToOpposite();
        return *this;
    }

    ContourCirculator operator--(int)
    {
        ContourCirculator ret(*this);
        operator--();
        return ret;
    }

    ContourCirculator & jumpToOpposite()
    {
        ray_.jumpToOpposite();
        return *this;
    }

    bool operator==(ContourCirculator const & o) const
    {
        return ray_ == o.ray_;
    }

    bool operator!=(ContourCirculator const & o) const
    {
        return ray_ != o.ray_;
    }

    FourEightSegmentation * segmentation() const
    {
        return ray_.segmentation();
    }

    int nodeLabel() const { return ray_.nodeLabel(); }
    int edgeLabel() const { return ray_.edgeLabel(); }
    int leftFaceLabel() const { return ray_.leftFaceLabel(); }
    int rightFaceLabel() const { return ray_.rightFaceLabel(); }

    int degree() const { return ray_.degree(); }
    float x() const { return ray_.x(); }
    float y() const { return ray_.y(); }

    RayCirculator const & ray() const { return ray_; }
};

// -------------------------------------------------------------------
//                         FourEightSegmentation
// -------------------------------------------------------------------
class FourEightSegmentation
{
public:
    typedef BImage::Iterator CellImageIterator;
    typedef IImage::Iterator LabelImageIterator;

private:
    struct NodeInfo
    {
        int label;
        float x, y;
        int size;
        int degree;
        RayCirculator ray;

        NodeInfo() : label(-1) {}
        bool initialized() const { return label >= 0; }
    };

    struct EdgeInfo
    {
        int label;
        RayCirculator start, end;

        EdgeInfo() : label(-1) {}
        bool initialized() const { return label >= 0; }
    };

    struct FaceInfo
    {
        int label;
        Diff2D anchor;
        std::vector<ContourCirculator> contours;

        FaceInfo() : label(-1) {}
        bool initialized() const { return label >= 0; }
    };

    typedef std::vector<NodeInfo> NodeList;
    typedef std::vector<EdgeInfo> EdgeList;
    typedef std::vector<FaceInfo> FaceList;

public:
    typedef NodeList::iterator NodeIterator;
    typedef EdgeList::iterator EdgeIterator;
    typedef FaceList::iterator FaceIterator;

    typedef std::vector<ContourCirculator>::iterator BoundaryComponentsIterator;

    // -------------------------------------------------------------------
    //                  FourEightSegmentation::NodeAccessor
    // -------------------------------------------------------------------
    struct NodeAccessor
    {
        int degree(NodeIterator & i) const
        {
            return (*i).degree;
        }

        float x(NodeIterator & i) const
        {
            return (*i).x;
        }

        float y(NodeIterator & i) const
        {
            return (*i).y;
        }

        int label(NodeIterator & i) const
        {
            return (*i).label;
        }

        RayCirculator rayCirculator(NodeIterator & i) const
        {
            return (*i).ray;
        }
    };

    // -------------------------------------------------------------------
    //               FourEightSegmentation::NodeAtStartAccessor
    // -------------------------------------------------------------------
    struct NodeAtStartAccessor
    {
        int degree(RayCirculator & i) const
        {
            return i.degree();
        }

        int degree(ContourCirculator & i) const
        {
            return i.degree();
        }

        int degree(EdgeIterator & i) const
        {
            return degree((*i).start);
        }

        float x(RayCirculator & i) const
        {
            return i.x();
        }

        float y(RayCirculator & i) const
        {
            return i.y();
        }

        float x(ContourCirculator & i) const
        {
            return i.x();
        }

        float y(ContourCirculator & i) const
        {
            return i.y();
        }

        float x(EdgeIterator & i) const
        {
            return x((*i).start);
        }

        float y(EdgeIterator & i) const
        {
            return y((*i).start);
        }

        int label(RayCirculator & i) const
        {
            return i.nodeLabel();
        }

        int label(ContourCirculator & i) const
        {
            return i.nodeLabel();
        }

        int label(EdgeIterator & i) const
        {
            return label((*i).start);
        }

        RayCirculator rayCirculator(ContourCirculator & i) const
        {
            return i.ray();
        }

        RayCirculator rayCirculator(EdgeIterator & i) const
        {
            return (*i).start;
        }

        NodeIterator nodeIterator(RayCirculator & i) const
        {
            return i.segmentation()->findNode(i.nodeLabel());
        }

        NodeIterator nodeIterator(ContourCirculator & i) const
        {
            return i.segmentation()->findNode(i.nodeLabel());
        }

        NodeIterator nodeIterator(EdgeIterator & i) const
        {
            return nodeIterator((*i).start);
        }
    };

    // -------------------------------------------------------------------
    //                FourEightSegmentation::NodeAtEndAccessor
    // -------------------------------------------------------------------
    struct NodeAtEndAccessor
    {
        int degree(RayCirculator i) const
        {
            return i.jumpToOpposite().degree();
        }

        int degree(ContourCirculator i) const
        {
            return i.jumpToOpposite().degree();
        }

        int degree(EdgeIterator & i) const
        {
            return (*i).end.degree();
        }

        float x(RayCirculator i) const
        {
            return i.jumpToOpposite().x();
        }

        float y(RayCirculator i) const
        {
            return i.jumpToOpposite().y();
        }

        float x(ContourCirculator i) const
        {
            return i.jumpToOpposite().x();
        }

        float y(ContourCirculator i) const
        {
            return i.jumpToOpposite().y();
        }

        float x(EdgeIterator & i) const
        {
            return (*i).end.x();
        }

        float y(EdgeIterator & i) const
        {
            return (*i).end.y();
        }

        int label(RayCirculator i) const
        {
            return i.jumpToOpposite().nodeLabel();
        }

        int label(ContourCirculator i) const
        {
            return i.jumpToOpposite().nodeLabel();
        }

        int label(EdgeIterator & i) const
        {
            return (*i).end.nodeLabel();
        }

        RayCirculator rayCirculator(EdgeIterator & i) const
        {
            return (*i).end;
        }

        NodeIterator nodeIterator(RayCirculator & i) const
        {
            return i.segmentation()->findNode(label(i));
        }

        NodeIterator nodeIterator(ContourCirculator & i) const
        {
            return i.segmentation()->findNode(label(i));
        }

        NodeIterator nodeIterator(EdgeIterator & i) const
        {
            return (*i).end.segmentation()->findNode(label(i));
        }
    };

    // -------------------------------------------------------------------
    //                  FourEightSegmentation::EdgeAccessor
    // -------------------------------------------------------------------
    struct EdgeAccessor
    {
        int label(RayCirculator & i) const
        {
            return i.edgeLabel();
        }

        int label(ContourCirculator & i) const
        {
            return i.edgeLabel();
        }

        int label(EdgeIterator & i) const
        {
            return (*i).label;
        }
    };

    // -------------------------------------------------------------------
    //               FourEightSegmentation::FaceAtLeftAccessor
    // -------------------------------------------------------------------
    struct FaceAtLeftAccessor
    {
        int label(RayCirculator & i) const
        {
            return i.leftFaceLabel();
        }

        int label(ContourCirculator & i) const
        {
            return i.leftFaceLabel();
        }

        ContourCirculator contourCirculator(RayCirculator & i) const
        {
            return ContourCirculator(i);
        }

        ContourCirculator contourCirculator(EdgeIterator & i) const
        {
            return ContourCirculator((*i).start);
        }
    };

    // -------------------------------------------------------------------
    //               FourEightSegmentation::FaceAtRightAccessor
    // -------------------------------------------------------------------
    struct FaceAtRightAccessor
    {
        int label(RayCirculator & i) const
        {
            return i.rightFaceLabel();
        }

        int label(ContourCirculator & i) const
        {
            return i.rightFaceLabel();
        }

        ContourCirculator contourCirculator(RayCirculator & i) const
        {
            return ContourCirculator(i).jumpToOpposite();
        }

        ContourCirculator contourCirculator(EdgeIterator & i) const
        {
            return ContourCirculator((*i).end);
        }
    };

    // -------------------------------------------------------------------
    //                  FourEightSegmentation::FaceAccessor
    // -------------------------------------------------------------------
    struct FaceAccessor
    {
        int label(FaceIterator & i) const
        {
            return (*i).label;
        }

        int countBoundaryComponents(FaceIterator & i) const
        {
            return (*i).contours.size();
        }

        BoundaryComponentsIterator beginBoundaryComponentsIterator(FaceIterator & i) const
        {
            return (*i).contours.begin();
        }

        BoundaryComponentsIterator endBoundaryComponentsIterator(FaceIterator & i) const
        {
            return (*i).contours.end();
        }

        ContourCirculator contourCirculator(BoundaryComponentsIterator & i) const
        {
            return *i;
        }
    };

    // -------------------------------------------------------------------

    template <class SrcIter, class SrcAcc>
    void init(SrcIter ul, SrcIter lr, SrcAcc src)
    {
        width_ = lr.x - ul.x;
        height_ = lr.y - ul.y;
        int totalwidth = width_ + 4;
        int totalheight = height_ + 4;

        cellImage.resize(totalwidth, totalheight);
        cellImage = CellConfigurationRegion;

        labelImage.resize(totalwidth, totalheight);
        labelImage = 0;

        cells = cellImage.upperLeft() + Diff2D(2,2);
        labels = labelImage.upperLeft() + Diff2D(2,2);

        // extract contours in input image and put frame around them
        BImage contourImage(totalwidth, totalheight);
        initFourEightSegmentationContourImage(ul, lr, src, contourImage);

        initCellImage(contourImage);

        int nodeCount = label0Cells();
        int edgeCount = label1Cells(nodeCount);
        int faceCount = label2Cells(contourImage);

        labelCircles(nodeCount, edgeCount);

        // decrement labels:
        IImage::ScanOrderIterator i = labelImage.begin();
        IImage::ScanOrderIterator iend = labelImage.end();
        for(; i != iend; ++i)
            --(*i);

        initNodeList(nodeCount);
        initEdgeList(edgeCount);
        initFaceList(contourImage, faceCount);
    }

    template <class SrcIter, class SrcAcc>
    void init(triple<SrcIter, SrcIter, SrcAcc> src)
    {
        init(src.first, src.second, src.third);
    }

    NodeIterator findNode(int const & i) const { return const_cast<NodeList &>(nodeList).begin() + i; }
    EdgeIterator findEdge(int const & i) const { return const_cast<EdgeList &>(edgeList).begin() + i; }
    FaceIterator findFace(int const & i) const { return const_cast<FaceList &>(faceList).begin() + i; }

    NodeIterator nodesBegin() const { return const_cast<NodeList &>(nodeList).begin(); }
    NodeIterator nodesEnd() const { return const_cast<NodeList &>(nodeList).end(); }
    EdgeIterator edgesBegin() const { return const_cast<EdgeList &>(edgeList).begin(); }
    EdgeIterator edgesEnd() const { return const_cast<EdgeList &>(edgeList).end(); }
    FaceIterator facesBegin() const { return const_cast<FaceList &>(faceList).begin(); }
    FaceIterator facesEnd() const { return const_cast<FaceList &>(faceList).end(); }

    CellImageIterator cellsUpperLeft() const { return cells; }
    CellImageIterator cellsLowerRight() const { return cells + Diff2D(width_, height_); }

    LabelImageIterator labelsUpperLeft() const { return labels; }
    LabelImageIterator labelsLowerRight() const { return labels + Diff2D(width_, height_); }

    int width() const { return width_; }
    int height() const { return height_; }

    int nodeCount() const { return nodeList.size(); }
    int edgeCount() const { return edgeList.size(); }
    int faceCount() const { return faceList.size(); }

    BImage cellImage;
    BImage::Iterator cells;
    IImage labelImage;
    IImage::Iterator labels;

    NodeInfo const & node(int i) const { return nodeList[i]; }
    EdgeInfo const & edge(int i) const { return edgeList[i]; }
    FaceInfo const & face(int i) const { return faceList[i]; }

private:
    void initCellImage(BImage & contourImage);
    int label0Cells();
    int label1Cells(int nodeCount);
    int label2Cells(BImage & contourImage);
    void labelCircles(int & nodeCount, int & edgeCount);

    void labelLine(OldNeighborhoodCirculator rayAtStart, int newLabel);

    void initNodeList(int nodeCount);
    void initEdgeList(int edgeCount);
    void initFaceList(BImage & contourImage, int faceCount);

private:
    int width_, height_;

    NodeList nodeList;
    EdgeList edgeList;
    FaceList faceList;
};

// -------------------------------------------------------------------
//                  OldNeighborhoodCirculator functions
// -------------------------------------------------------------------
inline unsigned char OldNeighborhoodCirculator::cell() const
{
    return segmentation_->cells[center()];
}

inline unsigned char OldNeighborhoodCirculator::neighborCell() const
{
    return segmentation_->cells[location()];
}

inline unsigned char OldNeighborhoodCirculator::forwardNeighborCell() const
{
    return segmentation_->cells[center() + neighbor_[1]];
}

inline unsigned char OldNeighborhoodCirculator::backwardNeighborCell() const
{
    return segmentation_->cells[center() + neighbor_[-1]];
}

inline int OldNeighborhoodCirculator::label() const
{
    return segmentation_->labels[center()];
}

inline int OldNeighborhoodCirculator::neighborLabel() const
{
    return segmentation_->labels[location()];
}

inline int OldNeighborhoodCirculator::forwardNeighborLabel() const
{
    return segmentation_->labels[center() + neighbor_[1]];
}

inline int OldNeighborhoodCirculator::backwardNeighborLabel() const
{
    return segmentation_->labels[center() + neighbor_[-1]];
}

// -------------------------------------------------------------------
//                        RayCirculator functions
// -------------------------------------------------------------------
inline int RayCirculator::degree() const
{
    return segmentation()->node(nodeLabel()).degree;
}

inline float RayCirculator::x() const
{
    return segmentation()->node(nodeLabel()).x;
}

inline float RayCirculator::y() const
{
    return segmentation()->node(nodeLabel()).y;
}

// -------------------------------------------------------------------
//                    FourEightSegmentation functions
// -------------------------------------------------------------------
template <class SrcIter, class SrcAcc>
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
    SrcType zero = NumericTraits<SrcType>::zero();
    for(y=0; y<h; ++y, ++ul.y)
    {
        SrcIter sx = ul;
        for(x=0; x<w; ++x, ++sx.x)
        {
            if(src(sx) == zero)  contourImage(x+2, y+2) = 1;
        }
    }
}

void FourEightSegmentation::initCellImage(BImage & contourImage)
{
    BImage::Iterator raw = contourImage.upperLeft() + Diff2D(1,1);

    int x,y;

    for(y=-1; y<=height_; ++y, ++raw.y)
    {
        BImage::Iterator rx = raw;
        for(x=-1; x<=width_; ++x, ++rx.x)
        {
            if(*rx == 0)
            {
                cells(x,y) = CellConfigurationRegion;
            }
            else
            {
                EightNeighborOffsetCirculator neighbors(EightNeighborOffsetCirculator::SouthEast);
                EightNeighborOffsetCirculator end = neighbors;

                int conf = 0;
                do
                {
                    conf = (conf << 1) | rx[neighbors.diff()];
                }
                while(--neighbors != end);

                if(cellConfigurations[conf] == CellConfigurationError)
                {
                    char message[200];
                    sprintf(message, "FourEightSegmentation::init(): "
                            "Configuration at (%d, %d) must be thinned further",
                            x, y);

                    vigra_precondition(0, message);
                }

                cells(x,y) = cellConfigurations[conf];
            }
        }
    }
}

// -------------------------------------------------------------------

int FourEightSegmentation::label0Cells()
{
    BImage nodeImage(width_+4, height_+4);
    BImage::Iterator nodes = nodeImage.upperLeft() + Diff2D(2,2);

    int x,y;

    for(y=-2; y<height_+2; ++y)
    {
        for(x=-2; x<width_+2; ++x)
        {
            if(cells(x,y) == CellConfigurationVertex)
            {
                nodes(x,y) = 1;

                // test for forbidden configuration
                OldNeighborhoodCirculator n(this, Diff2D(x,y));
                OldNeighborhoodCirculator nend = n;

                do
                {
                    if(n.neighborCell() == CellConfigurationLine &&
                       n.forwardNeighborCell() == CellConfigurationLine)
                    {
                        char msg[200];
                        sprintf(msg,"initFourEightSegmentation(): "
                                    "Node at (%d, %d) has two incident edgels form the same edge",
                                x,y);
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

    return labelImageWithBackground(srcImageRange(nodeImage), destImage(labelImage), true, 0);
}

// -------------------------------------------------------------------

int FourEightSegmentation::label1Cells(int nodeCount)
{
    int x,y;

    std::vector<int> nodeProcessed(nodeCount + 1, 0);

    int edgeCount = 0;

    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(cells(x,y) != CellConfigurationVertex)
                continue;
            if(nodeProcessed[labels(x,y)])
                continue;

            nodeProcessed[labels(x,y)] = 1;

            RayCirculator rayAtStart(this, Diff2D(x,y), EightNeighborOffsetCirculator::West);
            RayCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart.edgeLabel() != 0)
                    continue;

                labelLine(rayAtStart.neighborhoodCirculator(), ++edgeCount);
            }
            while(++rayAtStart != rayEnd);
        }
    }

    return edgeCount;
}

// -------------------------------------------------------------------

int FourEightSegmentation::label2Cells(BImage & contourImage)
{
    return labelImageWithBackground(srcImageRange(contourImage), destImage(labelImage), false, 1);
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelCircles(int & nodeCount, int & edgeCount)
{
    int x,y;

    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(labels(x,y) != 0)
                continue;

            // found a circle

            // mark first point as node
            cells(x,y) = CellConfigurationVertex;
            labels(x,y) = ++nodeCount;

            OldNeighborhoodCirculator rayAtStart(this, Diff2D(x,y));
            OldNeighborhoodCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart.neighborCell() != CellConfigurationLine)
                    continue;
                if(rayAtStart.neighborLabel() != 0)
                    continue;

                labelLine(rayAtStart, ++edgeCount);
            }
            while(++rayAtStart != rayEnd);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelLine(OldNeighborhoodCirculator rayAtStart,
                                      int newLabel)
{
    EdgelIterator line(rayAtStart);

    // follow the line and relabel it
    for(;!line.isEnd(); ++line)
    {
        labels[line.location()] = newLabel;
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initNodeList(int nodeCount)
{
    nodeList.resize(nodeCount, NodeInfo());
    std::vector<int> areas(nodeCount, 0);

    int x,y;

    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(cells(x,y) != CellConfigurationVertex)
                continue;

            int index = labels(x,y);

            if(!nodeList[index].initialized())
            {
                nodeList[index].label = labels(x,y);

                nodeList[index].x = x;
                nodeList[index].y = y;
                nodeList[index].size = 1;
                nodeList[index].ray = RayCirculator(this, Diff2D(x,y), EightNeighborOffsetCirculator::West);

                // calculate degree of the node
                RayCirculator r = nodeList[index].ray;
                RayCirculator rend = nodeList[index].ray;
                nodeList[index].degree = 0;
                do
                {
                    ++nodeList[index].degree;
                }
                while(++r != rend);

                // calculate area from following the outer contour of the node
                OldNeighborhoodCirculator
                    neighbor(this, Diff2D(x,y), EightNeighborOffsetCirculator::West);
                CrackCirculator crack(neighbor);
                CrackCirculator crackend(crack);

                do
                {
                    areas[index] += crack.diff().x * crack.pos().y -
                                    crack.diff().y * crack.pos().x;
                }
                while(++crack != crackend);

                areas[index] /= 2;
            }
            else
            {
                nodeList[index].x += x;
                nodeList[index].y += y;

                // calculate area from counting the pixels of the node
                nodeList[index].size += 1;
            }
        }
    }

    int i;
    for(i=0; i < nodeList.size(); ++i)
    {
        nodeList[i].x /= nodeList[i].size;
        nodeList[i].y /= nodeList[i].size;

        // methods to calculate the area must yield identical values
        if(areas[i] != nodeList[i].size)
        {
            char msg[200];
            sprintf(msg, "initFourEightSegmentation(): "
                         "Node at (%d, %d) has a hole",
                    nodeList[i].ray.center().x, nodeList[i].ray.center().y);
            vigra_precondition(0, msg);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initEdgeList(int edgeCount)
{
    edgeList.resize(edgeCount, EdgeInfo());

    NodeAccessor node;
    EdgeAccessor edge;

    NodeIterator n = nodeList.begin();
    NodeIterator nend = nodeList.end();

    for(; n != nend; ++n)
    {
        RayCirculator r = node.rayCirculator(n);
        RayCirculator rend = r;

        do
        {
            int index = edge.label(r);
            if(!edgeList[index].initialized())
            {
                edgeList[index].label = index;
                edgeList[index].start = r;
                edgeList[index].end = r;
                edgeList[index].end.jumpToOpposite();
            }
        }
        while(++r != rend);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initFaceList(BImage & contourImage, int number_of_faces)
{
    faceList.resize(number_of_faces, FaceInfo());

    IImage contourLabelImage(width_+4, height_+4);
    IImage::Iterator contourLabel = contourLabelImage.upperLeft() + Diff2D(2,2);
    contourLabelImage = 0;
    int countContourComponents =
        labelImageWithBackground(srcImageRange(contourImage),
                                 destImage(contourLabelImage), true, 0);

    std::vector<bool> contourProcessed(countContourComponents + 1, false);

    // process outer face
    faceList[0].label= 0;
    faceList[0].anchor = Diff2D(-2, -2);
    RayCirculator ray(this, Diff2D(-1, -1), EightNeighborOffsetCirculator::West);
    --ray;
    faceList[0].contours.push_back(ContourCirculator(ray));
    contourProcessed[contourLabel(-1, -1)] = true;

    FaceAtLeftAccessor leftface;

    int x,y;

    for(y=0; y<height_; ++y)
    {
        for(x=0; x<width_; ++x)
        {
            if(cells(x,y) != CellConfigurationRegion)
                continue;

            int index = labels(x,y);

            if(!faceList[index].initialized())
            {
                faceList[index].label = index;
                faceList[index].anchor = Diff2D(x,y);

                // find incident node
                if(cells(x-1,y) == CellConfigurationVertex)
                {
                    // this is the node
                    RayCirculator ray(this, Diff2D(x-1, y), EightNeighborOffsetCirculator::East);
                    --ray;

                    vigra_invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");

                    faceList[index].contours.push_back(ContourCirculator(ray));
                }
                else
                {
                    // its an edge
                    int lineindex = labels(x-1,y);

                    ContourCirculator c(edgeList[lineindex].start);
                    if(leftface.label(c) != index)
                    {
                        c.jumpToOpposite();
                    }

                    vigra_invariant(leftface.label(c) == index, "FourEightSegmentation::initFaceList()");

                    faceList[index].contours.push_back(c);
                }
            }
            else
            {
                // look for inner contours
                OldNeighborhoodCirculator
                    neighbor(this, Diff2D(x,y), EightNeighborOffsetCirculator::East);
                OldNeighborhoodCirculator
                    nend = neighbor;

                do
                {
                    int bindex = contourLabel[neighbor.location()];
                    if(bindex == 0 || contourProcessed[bindex])
                        continue;

                    // found an inner contour
                    contourProcessed[bindex] = true;

                    // find incident node
                    if(cells[neighbor.location()] == CellConfigurationVertex)
                    {
                        // this is the node
                        OldNeighborhoodCirculator n = neighbor;
                        n.jumpToOpposite();
                        RayCirculator ray(n);
                        --ray;

                        vigra_invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(ContourCirculator(ray));
                    }
                    else
                    {
                        // its an edge
                        int lineindex = neighbor.neighborLabel();

                        ContourCirculator c(edgeList[lineindex].start);
                        if(leftface.label(c) != index)
                        {
                            c.jumpToOpposite();
                        }

                        vigra_invariant(leftface.label(c) == index, "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(c);
                    }
                }
                while(++neighbor != nend);
            }
        }
    }
}

} // namespace FourEightSegmentation

} // namespace vigra

#endif /* VIGRA_FOUREIGHTSEGMENTATION_HXX */
