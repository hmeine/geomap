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

struct CellPixel
{
    CellConfiguration type_;
    int label_;

    CellPixel() {}
    CellPixel(CellConfiguration type, int label = 0)
        : type_(type), label_(label)
    {}

    inline int label() const { return label_; }
    inline void setLabel(int label) { label_= label; }
    inline void setLabel(int label, CellConfiguration) { label_= label; }
    inline CellConfiguration type() const { return type_; }
    inline void setType(CellConfiguration type) { type_ = type; }
};

typedef BasicImage<CellPixel> CellImage;
typedef vigra::NeighborhoodCirculator<CellImage::Iterator, EightNeighborOffsetCirculator>
    CellImageEightCirculator;

struct CellImageLabelAccessor
{
    typedef int value_type;

    template<class Iterator>
    int operator()(const Iterator &it) const
    {
        return it->label();
    }

    template<class Iterator>
    void set(int label, const Iterator &it) const
    {
        it->setLabel(label);
    }
};

template<CellConfiguration type>
struct CellImageLabelWriter
{
    typedef int value_type;

    template<class Iterator>
    void set(int label, const Iterator &it) const
    {
        it->setLabel(label, type);
    }
};

// -------------------------------------------------------------------
//                            CrackCirculator
// -------------------------------------------------------------------
class CrackCirculator
{
protected:
    CellImageEightCirculator neighborCirc_;
    int typeLabel_;
    Diff2D pos_;

public:
    CrackCirculator(CellImageEightCirculator const & n)
        : neighborCirc_(n),
          typeLabel_(n.center()->type()),
          pos_(0, 0)
    {
        neighborCirc_.turnLeft();
    }

    CrackCirculator & operator++()
    {
        pos_ += neighborCirc_.diff();

        neighborCirc_--;

        if(neighborCirc_->type() == typeLabel_)
        {
            neighborCirc_.moveCenterToNeighbor(); // TODO: simplify moveCenterToNeighbor()s
            --neighborCirc_;
        }
        else
        {
            neighborCirc_.moveCenterToNeighbor(); // jump out
            neighborCirc_ += 3;
            if(neighborCirc_->type() == typeLabel_)
            {
                neighborCirc_.moveCenterToNeighbor();
                neighborCirc_.turnRight();
            }
            else
            {
                neighborCirc_.moveCenterToNeighbor();
                neighborCirc_.turnLeft();
                neighborCirc_.moveCenterToNeighbor();
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

    Diff2D pos() const { return pos_; }
};

// -------------------------------------------------------------------
//                             EdgelIterator
// -------------------------------------------------------------------
/**
 * An EdgelIterator starts walking in the direction of the Circulator
 * given on construction and walks along CellConfigurationLine type'ed
 * pixels. isEnd() will become true if the iterator steps on a
 * CellConfigurationVertex pixel or tries to pass one diagonally.
 *
 * It is used by the RayCirculator for jumpToOpposite() and in
 * labelLine().
 */
class EdgelIterator
{
    CellImageEightCirculator neighborCirc_;
    bool isEnd_;

public:
    EdgelIterator(CellImageEightCirculator const & n)
        : neighborCirc_(n), isEnd_(false)
    {}

    CellImageEightCirculator::base_type location() const
    {
        return neighborCirc_.base();
    }

    bool isEnd() const
    {
        return isEnd_;
    }

    EdgelIterator & operator++()
    {
        neighborCirc_.moveCenterToNeighbor();
        neighborCirc_.turnRight();

        while(1)
        {
            if(neighborCirc_->type() == CellConfigurationVertex)
            {
                isEnd_ = true;
                break;
            }
            if(neighborCirc_->type() == CellConfigurationLine)
            {
                break;
            }
            ++neighborCirc_;
        }

        if(neighborCirc_.isDiagonal() &&
           neighborCirc_[1].type() == CellConfigurationVertex)
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
        neighborCirc_.swapCenterNeighbor();
        return *this;
    }

    operator CellImageEightCirculator()
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
    CellImageEightCirculator neighborCirc_;
    FourEightSegmentation * segmentation_;
    bool isSingular_;

public:
    // default constructor needed for NodeInfo/EdgeInfo, init() must not be called!
    //RayCirculator() : {}

    RayCirculator(FourEightSegmentation * segmentation,
                  CellImageEightCirculator const & circ)
        : neighborCirc_(circ),
          segmentation_(segmentation)
    {
        vigra_precondition(neighborCirc_.center()->type() == CellConfigurationVertex,
        "FourEightSegmentation::RayCirculator(): center is not a node");

        vigra_precondition(neighborCirc_->type() != CellConfigurationVertex,
        "FourEightSegmentation::RayCirculator(): neighbor is a node");

        CellImageEightCirculator n = neighborCirc_;
        isSingular_ = true;
        do
        {
            if(n->type() != CellConfigurationRegion)
            {
                isSingular_ = false;
                break;
            }
        }
        while(++n != neighborCirc_);

        if(neighborCirc_->type() != CellConfigurationLine)
            operator++();
    }

    RayCirculator & operator++()
    {
        if(isSingular_) return *this;

        tryNext();

        while(neighborCirc_->type() != CellConfigurationLine)
        {
            if(neighborCirc_->type() == CellConfigurationVertex)
            {
                neighborCirc_.swapCenterNeighbor();
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

        while(neighborCirc_->type() != CellConfigurationLine)
        {
            if(neighborCirc_->type() == CellConfigurationVertex)
            {
                neighborCirc_.swapCenterNeighbor();
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
        return segmentation_;
        //return neighborCirc_.segmentation();
    }

    CellImage::Iterator center() const { return neighborCirc_.center(); }

    int nodeLabel() const { return neighborCirc_.center()->label(); }
    int edgeLabel() const { return neighborCirc_->label(); }
    int leftFaceLabel() const { return neighborCirc_[1].label(); }
    int rightFaceLabel() const { return neighborCirc_[-1].label(); }

    inline int degree() const;
    inline float x() const;
    inline float y() const;

    const CellImageEightCirculator &neighborCirculator() const
    {
        return neighborCirc_;
    }

private:
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
        return (neighborCirc_->type() == CellConfigurationLine &&
                (neighborCirc_[1].type() == CellConfigurationVertex ||
                 neighborCirc_[-1].type() == CellConfigurationVertex));
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
    struct NodeInfo
    {
        int label;
        float x, y;
        int size;
        int degree;
        RayCirculator ray;

        NodeInfo(RayCirculator r) : label(-1), ray(r) {}
        bool initialized() const { return label >= 0; }
    };

    struct EdgeInfo
    {
        int label;
        RayCirculator start, end;

        EdgeInfo(RayCirculator s, RayCirculator e) : label(-1), start(s), end(e) {}
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
        cellImage = CellPixel(CellConfigurationRegion, 0);

        cells = cellImage.upperLeft() + Diff2D(2,2);

        // extract contours in input image and put frame around them
        BImage contourImage(totalwidth, totalheight);
        initFourEightSegmentationContourImage(ul, lr, src, contourImage);

        initCellImage(contourImage);

        int nodeCount = label0Cells();
        int edgeCount = label1Cells(nodeCount);
        int faceCount = label2Cells(contourImage);

        labelCircles(nodeCount, edgeCount);

        // decrement labels:
        CellImage::ScanOrderIterator i = cellImage.begin();
        CellImage::ScanOrderIterator iend = cellImage.end();
        for(; i != iend; ++i)
            i->setLabel(i->label() - 1);

        std::cerr << "FourEightSegmentation::initNodeList(nodeCount= " << nodeCount << ")\n";
        initNodeList(nodeCount);
        std::cerr << "FourEightSegmentation::initEdgeList(edgeCount= " << edgeCount << ")\n";
        initEdgeList(edgeCount);
        std::cerr << "FourEightSegmentation::initFaceList(faceCount= " << faceCount << ")\n";
        initFaceList(contourImage, faceCount);
        std::cerr << "FourEightSegmentation::~init()\n";
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

    CellImage::Iterator cellsUpperLeft() const { return cells; }
    CellImage::Iterator cellsLowerRight() const { return cells + Diff2D(width_, height_); }

    int width() const { return width_; }
    int height() const { return height_; }

    int nodeCount() const { return nodeList.size(); }
    int edgeCount() const { return edgeList.size(); }
    int faceCount() const { return faceList.size(); }

    CellImage cellImage;
    CellImage::Iterator cells;

    NodeInfo const & node(int i) const { return nodeList[i]; }
    EdgeInfo const & edge(int i) const { return edgeList[i]; }
    FaceInfo const & face(int i) const { return faceList[i]; }

private:
    void initCellImage(BImage & contourImage);
    int label0Cells();
    int label1Cells(int nodeCount);
    int label2Cells(BImage & contourImage);
    void labelCircles(int & nodeCount, int & edgeCount);

    void labelLine(CellImageEightCirculator rayAtStart, int newLabel);

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

    for(int y=-1; y<=height_; ++y, ++raw.y)
    {
        BImage::Iterator rx = raw;
        for(int x=-1; x<=width_; ++x, ++rx.x)
        {
            if(*rx == 0)
            {
                cells(x,y).setType(CellConfigurationRegion);
            }
            else
            {
                vigra::NeighborhoodCirculator<BImage::Iterator, EightNeighborOffsetCirculator>
                    neighbors(rx, EightNeighborOffsetCirculator::SouthEast);
                vigra::NeighborhoodCirculator<BImage::Iterator, EightNeighborOffsetCirculator>
                    end = neighbors;

                int conf = 0;
                do
                {
                    conf = (conf << 1) | *neighbors;
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

                cells(x,y).setType(cellConfigurations[conf]);
            }
            std::cerr << cells(x,y).type();
        }
        std::cerr << std::endl;
    }
}

// -------------------------------------------------------------------

int FourEightSegmentation::label0Cells()
{
    std::cerr << "FourEightSegmentation::label0Cells()\n";

    BImage nodeImage(width_+4, height_+4);
    BImage::Iterator nodes = nodeImage.upperLeft() + Diff2D(2,2);

    for(int y=-2; y<height_+2; ++y)
    {
        CellImage::Iterator cell = cells + Diff2D(-2, y);
        for(int x=-2; x<width_+2; ++x, ++cell.x)
        {
            if(cell->type() == CellConfigurationVertex)
            {
                std::cout << "node found at " << x << ", " << y << " " << std::endl;

                nodes(x,y) = 1;

                // test for forbidden configuration
                CellImageEightCirculator n(cell);
                CellImageEightCirculator nend = n;

                do
                {
                    if(n->type() == CellConfigurationLine && n[1].type() == CellConfigurationLine)
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

    return labelImageWithBackground(srcImageRange(nodeImage),
                                    destImage(cellImage, CellImageLabelWriter<CellConfigurationVertex>()), true, 0);
}

// -------------------------------------------------------------------

int FourEightSegmentation::label1Cells(int nodeCount)
{
    std::cerr << "FourEightSegmentation::label1Cells(" << nodeCount << ")\n";

    std::vector<bool> nodeProcessed(nodeCount + 1, false);

    int edgeCount = 0;

    for(int y=-1; y<=height_; ++y)
    {
        CellImage::Iterator cell = cells + Diff2D(-1, y);
        for(int x=-1; x<=width_; ++x, ++cell.x)
        {
            if(cell->type() != CellConfigurationVertex)
                continue;
            if(nodeProcessed[cell->label()])
                continue;

            std::cerr << "unprocessed node found at " << x << ", " << y
                      << "(label: " << cell->label() << ")\n";

            nodeProcessed[cell->label()] = true;

            RayCirculator rayAtStart(this, CellImageEightCirculator(cell, EightNeighborOffsetCirculator::West));
            RayCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart.edgeLabel() != 0)
                    continue;

                std::cerr << "labelling line.." << std::endl;
                labelLine(rayAtStart.neighborCirculator(), ++edgeCount);
                std::cerr << "labelling line done.." << std::endl;
            }
            while(++rayAtStart != rayEnd);
        }
    }

    return edgeCount;
}

// -------------------------------------------------------------------

int FourEightSegmentation::label2Cells(BImage & contourImage)
{
    std::cerr << "FourEightSegmentation::label2Cells()\n";
    return labelImageWithBackground(srcImageRange(contourImage),
                                    destImage(cellImage, CellImageLabelWriter<CellConfigurationRegion>()),
                                    false, 1);
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelCircles(int & nodeCount, int & edgeCount)
{
    std::cerr << "FourEightSegmentation::labelCircles(" << nodeCount
              << " nodes, " << nodeCount << " edges)\n";
    int x,y;

    for(y=-1; y<=height_; y++)
    {
        CellImage::Iterator cell = cells + Diff2D(-1, y);
        for(x=-1; x<=width_; x++, cell.x++)
        {
            if(cell->label() != 0)
                continue;

            // found a circle

            // mark first point as node
            (*cell) = CellPixel(CellConfigurationVertex, ++nodeCount);

            CellImageEightCirculator rayAtStart(cell);
            CellImageEightCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart->type() != CellConfigurationLine)
                    continue;
                if(rayAtStart->label() != 0)
                    continue;

                labelLine(rayAtStart, ++edgeCount);
            }
            while(++rayAtStart != rayEnd);
        }
    }
    std::cerr << "FourEightSegmentation::~labelCircles()\n";
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelLine(CellImageEightCirculator rayAtStart,
                                      int newLabel)
{
    EdgelIterator line(rayAtStart);

    // follow the line and relabel it
    for(;!line.isEnd(); ++line)
    {
        line.location()->setLabel(newLabel, CellConfigurationLine);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initNodeList(int nodeCount)
{
    nodeList.resize(nodeCount, NodeInfo(RayCirculator(this, CellImageEightCirculator(cells-Diff2D(1,1)))));
    std::vector<int> crackCirculatedAreas(nodeCount, 0);

    int x,y;

    for(y=-1; y<=height_; ++y)
    {
        CellImage::Iterator cell = cells + Diff2D(-1, y);
        for(x=-1; x<=width_; ++x, ++cell.x)
        {
            if(cell->type() != CellConfigurationVertex)
                continue;

            int index = cell->label();
            std::cerr << "found node at " << x << "," << y << " with label " << cell->label() << "\n";

            if(!nodeList[index].initialized())
            {
                std::cerr << "initializing node " << index << "\n";
                nodeList[index].label = index;

                nodeList[index].x = x;
                nodeList[index].y = y;
                nodeList[index].size = 1;
                std::cerr << "creating RayCirculator 1\n";
                nodeList[index].ray = RayCirculator(this, CellImageEightCirculator(cell, EightNeighborOffsetCirculator::West));

                // calculate degree of the node
                std::cerr << "creating RayCirculator 2&3\n";
                RayCirculator r = nodeList[index].ray;
                RayCirculator rend = nodeList[index].ray;
                nodeList[index].degree = 0;
                do
                {
                    ++nodeList[index].degree;
                }
                while(++r != rend);

                // calculate area from following the outer contour of the node
                CellImageEightCirculator
                    neighbor(cell, EightNeighborOffsetCirculator::West);

                CrackCirculator crack(neighbor);
                CrackCirculator crackend(crack);
                do
                {
                    std::cerr << "crack.diff(): (" << crack.diff().x << ", " << crack.diff().y << ") "
                              << "crack.pos(): (" << crack.pos().x << ", " << crack.pos().y << ")\n";
                    crackCirculatedAreas[index] += crack.diff().x * crack.pos().y -
                                                   crack.diff().y * crack.pos().x;
                }
                while(++crack != crackend);

                crackCirculatedAreas[index] /= 2;
                std::cerr << "calculated crackCirculatedAreas.\n";
            }
            else
            {
                std::cerr << "node " << index << " has label " << nodeList[index].label << "\n";
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
        if(crackCirculatedAreas[i] != nodeList[i].size)
        {
            std::cerr << "crackCirculatedAreas[i]==" << crackCirculatedAreas[i] << ", "
                      << "nodeList[i].size==" << nodeList[i].size << std::endl;
            char msg[200];
            sprintf(msg, "FourEightSegmentation::initNodeList(): "
                    "Node %d at (%d, %d) has a hole", i,
                    nodeList[i].ray.center().x, nodeList[i].ray.center().y);
            vigra_precondition(0, msg);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initEdgeList(int edgeCount)
{
    edgeList.resize(edgeCount, EdgeInfo(RayCirculator(this, CellImageEightCirculator(cells-Diff2D(1,1))),
                                        RayCirculator(this, CellImageEightCirculator(cells-Diff2D(1,1)))));

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
    RayCirculator ray(this, CellImageEightCirculator(cells + Diff2D(-1, -1), EightNeighborOffsetCirculator::West));
    --ray;
    faceList[0].contours.push_back(ContourCirculator(ray));
    contourProcessed[contourLabel(-1, -1)] = true;

    FaceAtLeftAccessor leftface;

    int x,y;

    for(y=0; y<height_; ++y)
    {
        for(x=0; x<width_; ++x)
        {
            if(cells(x,y).type() != CellConfigurationRegion)
                continue;

            int index = cells(x,y).label();

            if(!faceList[index].initialized())
            {
                faceList[index].label = index;
                faceList[index].anchor = Diff2D(x,y);

                // find incident node
                if(cells(x-1,y).type() == CellConfigurationVertex)
                {
                    // this is the node
                    RayCirculator ray(this, CellImageEightCirculator(cells + Diff2D(x-1, y)));
                    --ray;

                    vigra_invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");

                    faceList[index].contours.push_back(ContourCirculator(ray));
                }
                else
                {
                    // its an edge
                    int lineindex = cells(x-1, y).label();

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
                CellImageEightCirculator neighbor(cells + Diff2D(x,y));
                CellImageEightCirculator nend = neighbor;

                do
                {
                    int bindex = contourLabel[neighbor.base()-cells];
                    if(bindex == 0 || contourProcessed[bindex])
                        continue;

                    // found an inner contour
                    contourProcessed[bindex] = true;

                    // find incident node
                    if(neighbor->type() == CellConfigurationVertex)
                    {
                        // this is the node
                        CellImageEightCirculator n = neighbor;
                        n.swapCenterNeighbor();
                        RayCirculator ray(this, n);
                        --ray;

                        vigra_invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(ContourCirculator(ray));
                    }
                    else
                    {
                        // its an edge
                        int lineindex = neighbor->label();

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
