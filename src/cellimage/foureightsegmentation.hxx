#ifndef VIGRA_FOUREIGHTSEGMENTATION_HXX
#define VIGRA_FOUREIGHTSEGMENTATION_HXX

#include "vigra/error.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/labelimage.hxx"
#include "configurations.h"
#include "pixelneighborhood.hxx"

class FourEightSegmentation;

/***********************************************************************/
/*                                                                     */
/*            FourEightSegmentationNeighborhoodCirculator              */
/*                                                                     */
/***********************************************************************/

struct FourEightSegmentationNeighborhoodCirculator
{
    Dist2D center_;
    FourEightSegmentation * segmentation_;
    EightNeighborCoding neighbor;
    
    FourEightSegmentationNeighborhoodCirculator() {}
    
    FourEightSegmentationNeighborhoodCirculator(FourEightSegmentation * seg, 
           Dist2D const & loc, 
           EightNeighborCoding::Directions dir = EightNeighborCoding::East)
    : segmentation_(seg), center_(loc), neighbor(dir)
    {}
    
    FourEightSegmentationNeighborhoodCirculator & operator++()
    {
        ++neighbor;
        return *this;
    }
    
    FourEightSegmentationNeighborhoodCirculator operator++(int)
    {
        FourEightSegmentationNeighborhoodCirculator ret(*this);
        operator++();
        return ret;
    }
    
    FourEightSegmentationNeighborhoodCirculator & operator+=(int d)
    {
        neighbor += d;
        return *this;
    }
    
    FourEightSegmentationNeighborhoodCirculator & operator--()
    {
        --neighbor;
        return *this;
    }
    
    FourEightSegmentationNeighborhoodCirculator operator--(int)
    {
        FourEightSegmentationNeighborhoodCirculator ret(*this);
        operator--();
        return ret;
    }
    
    FourEightSegmentationNeighborhoodCirculator & operator-=(int d)
    {
        neighbor -= d;
        return *this;
    }
    
    FourEightSegmentationNeighborhoodCirculator & jumpToOpposite()
    {
        center_ += diff();
        neighbor.turnRound();
        
        return *this;
    }
    
    FourEightSegmentationNeighborhoodCirculator & turnRight() { 
        neighbor.turnRight(); return *this; 
    }
    FourEightSegmentationNeighborhoodCirculator & turnLeft()  { 
        neighbor.turnLeft(); return *this; 
    }
    FourEightSegmentationNeighborhoodCirculator & turnRound() { 
        neighbor.turnRound(); return *this; 
    }
    FourEightSegmentationNeighborhoodCirculator & translate(Dist2D const & d) { 
        center_ += d; return *this; 
    }
    
    bool operator==(FourEightSegmentationNeighborhoodCirculator const & o) const 
    {
        return center_ == o.center_ && neighbor == o.neighbor;
    }
    
    bool operator!=(FourEightSegmentationNeighborhoodCirculator const & o) const 
    {
        return center_ != o.center_ || neighbor != o.neighbor;
    }

    Dist2D const & center() const {
        return center_;
    }
    
    Dist2D const & diff() const
    {
        return neighbor.diff();
    }
    
    Dist2D location() const {
        return center_ + diff();
    }
    
    EightNeighborCoding::Directions directionCode() const { return *neighbor; }
    
    inline int cell() const;
    inline int label() const;
    inline int neighborCell() const; 
    inline int neighborLabel() const;
    inline int forwardNeighborLabel() const;
    inline int backwardNeighborLabel() const;
    inline int forwardNeighborCell() const;
    inline int backwardNeighborCell() const;
    
    bool isDiagonal() const { return neighbor.isDiagonal(); }
    
    FourEightSegmentation * segmentation() const { return segmentation_; }
};


/***********************************************************************/
/*                                                                     */
/*              FourEightSegmentationCrackCirculator                   */
/*                                                                     */
/***********************************************************************/

struct FourEightSegmentationCrackCirculator
{
    FourEightSegmentationNeighborhoodCirculator neighbor;
    int label;
    int x_, y_;
    
    FourEightSegmentationCrackCirculator(
        FourEightSegmentationNeighborhoodCirculator const & n)
    : neighbor(n), label(n.cell()),
      x_(0), y_(0)
    {
        neighbor.turnLeft();
    }
    
    FourEightSegmentationCrackCirculator & operator++()
    {
        x_ += neighbor.diff().width;
        y_ += neighbor.diff().height;
        
        --neighbor;
        neighbor.translate(neighbor.diff());
        
        if(neighbor.cell() == label) 
        {
            --neighbor;
        }
        else
        {
            neighbor += 3;
            neighbor.translate(neighbor.diff());
            if(neighbor.cell() == label) 
            {
                neighbor.turnRight();
            }
            else
            {
                neighbor += 2;
                neighbor.translate(neighbor.diff());
                neighbor.turnRight();
            }
        }
        
	return *this;
    }
    
    bool operator==(FourEightSegmentationCrackCirculator const & o) const
    {
        return neighbor == o.neighbor;
    }
    
    bool operator!=(FourEightSegmentationCrackCirculator const & o) const
    {
        return neighbor != o.neighbor;
    }
    
    Dist2D const & diff() const { return neighbor.diff(); }
    
    int x() const { return x_; }
    int y() const { return y_; }
    
};

/***********************************************************************/
/*                                                                     */
/*                FourEightSegmentationEdgelIterator                   */
/*                                                                     */
/***********************************************************************/

struct FourEightSegmentationEdgelIterator
{
    FourEightSegmentationNeighborhoodCirculator neighbor;
    bool is_end;
    
    FourEightSegmentationEdgelIterator(
        FourEightSegmentationNeighborhoodCirculator const & n)
    : neighbor(n), is_end(false)
    {}
    
    FourEightSegmentationEdgelIterator & operator++()
    {
        neighbor.jumpToOpposite();
        neighbor.turnLeft();

        while(1)
        {
            if(neighbor.neighborCell() == CellConfigurationsVertex) 
            {
                is_end = true;
                break;
            }
            if(neighbor.neighborCell() == CellConfigurationsLine) 
            {
                break;
            }
            ++neighbor;
        }

        if(neighbor.isDiagonal() &&
           neighbor.forwardNeighborCell() == CellConfigurationsVertex) 
        {
            ++neighbor;
            is_end = true;
        }

	return *this;
    }
    
    bool isEnd() const { return is_end; }
    
    Dist2D location() const { return neighbor.location(); }
};

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentationRayCirculator                  */
/*                                                                     */
/***********************************************************************/

struct FourEightSegmentationRayCirculator
{
  private:
    FourEightSegmentationNeighborhoodCirculator neighbor;
    bool is_singular;
    
  public:
    FourEightSegmentationRayCirculator() {}
    
    FourEightSegmentationRayCirculator(FourEightSegmentation * seg, Dist2D const & loc,
           EightNeighborCoding::Directions dir = EightNeighborCoding::East)
    : neighbor(seg, loc, dir)
    {
        init();
    }
    
    FourEightSegmentationRayCirculator(FourEightSegmentationNeighborhoodCirculator const & n)
    : neighbor(n)
    {
        init();
    }
    
    FourEightSegmentationRayCirculator & operator++()
    {
        if(is_singular) return *this;
        
        tryNext();

        while(neighbor.neighborCell() != CellConfigurationsLine)
        {
            if(neighbor.neighborCell() == CellConfigurationsVertex)
            {
                neighbor.jumpToOpposite();
            }
            tryNext();
        }
        return *this;
    }
    
    FourEightSegmentationRayCirculator operator++(int)
    {
        FourEightSegmentationRayCirculator ret(*this);
        operator++();
        return ret;
    }
    
    FourEightSegmentationRayCirculator & operator--()
    {
        if(is_singular) return *this;
        
        tryPrev();
        
        while(neighbor.neighborCell() != CellConfigurationsLine)
        {
            if(neighbor.neighborCell() == CellConfigurationsVertex)
            {
                neighbor.jumpToOpposite();
            }
            tryPrev();
        }
        return *this;
    }
    
    FourEightSegmentationRayCirculator operator--(int)
    {
        FourEightSegmentationRayCirculator ret(*this);
        operator--();
        return ret;
    }
    
    FourEightSegmentationRayCirculator &  jumpToOpposite()
    {
        if(is_singular) return *this;

        FourEightSegmentationEdgelIterator line(neighbor);
        
        while(!line.isEnd()) 
        {
            ++line;
        }
        
        neighbor = line.neighbor.jumpToOpposite();
        
        return *this;
    }
    
    bool operator==(FourEightSegmentationRayCirculator const & o) const 
    {
        return neighbor == o.neighbor;
    }
    
    bool operator!=(FourEightSegmentationRayCirculator const & o) const 
    {
        return neighbor != o.neighbor;
    }
    
    FourEightSegmentation * segmentation() const { return neighbor.segmentation(); }
    
    Dist2D const & center() const { return neighbor.center(); }
    
    int nodeLabel() const { return neighbor.label(); }
    int edgeLabel() const { return neighbor.neighborLabel(); }
    int leftFaceLabel() const { return neighbor.forwardNeighborLabel(); }
    int rightFaceLabel() const { return neighbor.backwardNeighborLabel(); }
    
    inline int degree() const;
    inline float x() const;
    inline float y() const;
    
  private:
  
    friend class FourEightSegmentation;
    
    void init()
    {
        precondition(neighbor.cell() == CellConfigurationsVertex, 
        "FourEightSegmentationRayCirculator(): center is not a node");
        
        precondition(neighbor.neighborCell() != CellConfigurationsVertex, 
        "FourEightSegmentationRayCirculator(): neighbor is a node");
        
        FourEightSegmentationNeighborhoodCirculator n = neighbor;
        is_singular = true;
        do
        {
            if(n.neighborCell() != CellConfigurationsRegion) 
            {
                is_singular = false;
                break;
            }
        }
        while(++n != neighbor);       

        if(neighbor.neighborCell() != CellConfigurationsLine) operator++();
    }

    void tryNext()
    {
        ++neighbor;
        
        if(badDiagonalConfig()) ++neighbor;
    }
    
    void tryPrev()
    {
        --neighbor;
        
        if(badDiagonalConfig()) --neighbor;
    }
    
    bool badDiagonalConfig() 
    {
        return (neighbor.neighborCell() == CellConfigurationsLine &&
                (neighbor.forwardNeighborCell() == CellConfigurationsVertex ||
                 neighbor.backwardNeighborCell() == CellConfigurationsVertex));
    }
};

/***********************************************************************/
/*                                                                     */
/*               FourEightSegmentationContourCirculator                */
/*                                                                     */
/***********************************************************************/

struct FourEightSegmentationContourCirculator
{
    FourEightSegmentationRayCirculator ray_;
    
    FourEightSegmentationContourCirculator(FourEightSegmentationRayCirculator r)
    : ray_(r)
    {}
    
    FourEightSegmentationContourCirculator & operator++()
    {
        ray_.jumpToOpposite();
        --ray_;
        return *this;
    }
    
    FourEightSegmentationContourCirculator operator++(int)
    {
        FourEightSegmentationContourCirculator ret(*this);
        operator++();
        return ret;
    }
    
    FourEightSegmentationContourCirculator & operator--()
    {
        ++ray_;
        ray_.jumpToOpposite();
        return *this;
    }
    
    FourEightSegmentationContourCirculator operator--(int)
    {
        FourEightSegmentationContourCirculator ret(*this);
        operator--();
        return ret;
    }
    
    FourEightSegmentationContourCirculator & jumpToOpposite()
    {
        ray_.jumpToOpposite();
        return *this;
    }
    
    bool operator==(FourEightSegmentationContourCirculator const & o) const
    {
        return ray_ == o.ray_;
    }
    
    bool operator!=(FourEightSegmentationContourCirculator const & o) const
    {
        return ray_ != o.ray_;
    }
    
    FourEightSegmentation * segmentation() const { return ray_.segmentation(); }
    
    int nodeLabel() const { return ray_.nodeLabel(); }
    int edgeLabel() const { return ray_.edgeLabel(); }
    int leftFaceLabel() const { return ray_.leftFaceLabel(); }
    int rightFaceLabel() const { return ray_.rightFaceLabel(); }
    
    int degree() const { return ray_.degree(); }
    float x() const { return ray_.x(); }
    float y() const { return ray_.y(); }

    FourEightSegmentationRayCirculator const & ray() const { return ray_; }
};

/***********************************************************************/
/*                                                                     */
/*                      FourEightSegmentation                          */
/*                                                                     */
/***********************************************************************/

class FourEightSegmentation
{
  public:
  
    typedef BImage::Iterator CellImageIterator;
    typedef IImage::Iterator LabelImageIterator;
    typedef FourEightSegmentationRayCirculator RayCirculator;
    typedef FourEightSegmentationContourCirculator ContourCirculator;
    
    friend class FourEightSegmentationNeighborhoodCirculator;
    friend class FourEightSegmentationRayCirculator;
    
  private:
    
    typedef FourEightSegmentationNeighborhoodCirculator NeighborhoodCirculator;

    struct NodeInfo
    {
        int label;
        float x, y;
        int size;
        int degree;
        RayCirculator ray; 
        
        NodeInfo() : label(-1) {}
        
        bool notInitialized() const { return label < 0; }
    };
  
    struct EdgeInfo
    {
        int label;
        RayCirculator start, end; 
        
        EdgeInfo() : label(-1) {}
        
        bool notInitialized() const { return label < 0; }
    };
  
    struct FaceInfo
    {
        int label;
        Dist2D anchor;
        vector<ContourCirculator> contours;
        
        FaceInfo() : label(-1) {}
        
        bool notInitialized() const { return label < 0; }
    };
    
    typedef std::vector<NodeInfo> NodeList;
    typedef std::vector<EdgeInfo> EdgeList;
    typedef std::vector<FaceInfo> FaceList;
    
  public:
    typedef NodeList::iterator NodeIterator;
    typedef EdgeList::iterator EdgeIterator;
    typedef FaceList::iterator FaceIterator;
    
    typedef vector<ContourCirculator>::iterator BoundaryComponentsIterator;

    /***********************************************************************/
    /*                                                                     */
    /*              FourEightSegmentation::NodeAccessor                    */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*          FourEightSegmentation::NodeAtStartAccessor                 */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*           FourEightSegmentation::NodeAtEndAccessor                  */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*              FourEightSegmentation::EdgeAccessor                    */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*            FourEightSegmentation::FaceAtLeftAccessor                */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*            FourEightSegmentation::FaceAtRightAccessor               */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*              FourEightSegmentation::FaceAccessor                    */
    /*                                                                     */
    /***********************************************************************/

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
  
    /***********************************************************************/
    /*                                                                     */
    /*                 FourEightSegmentation functions                     */
    /*                                                                     */
    /***********************************************************************/

    FourEightSegmentation()
    {}
    
    template <class SrcIter, class SrcAcc>
    void init(SrcIter ul, SrcIter lr, SrcAcc src)
    {
        resize(lr.x - ul.x, lr.y - ul.y);

        BImage rawborders(totalwidth_, totalheight_);
        BImage::Iterator tmp = rawborders.upperLeft() + Dist2D(2,2);

        initFourEightSegmentationRawBorders(ul, lr, src, rawborders);

        initCellImage(rawborders);
        
        int nodeCount, edgeCount, faceCount;

        nodeCount = label0Cells();
        edgeCount = label1Cells(nodeCount);
        faceCount = label2Cells(rawborders);

        labelCircles(nodeCount, edgeCount);

        decrementLabels();

        initNodeList(nodeCount);
        initEdgeList(edgeCount);
        initFaceList(rawborders, faceCount);
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
    CellImageIterator cellsLowerRight() const { return cells + Dist2D(width_, height_); }
    
    LabelImageIterator labelsUpperLeft() const { return labels; }
    LabelImageIterator labelsLowerRight() const { return labels + Dist2D(width_, height_); }
    
    int width() const { return width_; }
    int height() const { return height_; }
    
    int numberOfNodes() const { return nodeList.size(); }
    int numberOfEdges() const { return edgeList.size(); }
    int numberOfFaces() const { return faceList.size(); }
    
  private:
    
    NodeInfo const & node(int i) const { return nodeList[i]; }
    EdgeInfo const & edge(int i) const { return edgeList[i]; }
    FaceInfo const & face(int i) const { return faceList[i]; }
    
    void resize(int w, int h)
    {
        width_ = w;
        height_ = h,
        totalwidth_ = w + 4;
        totalheight_ = h + 4;
        
        cellimage.resize(totalwidth_, totalheight_);
        cellimage = CellConfigurationsRegion;
        
        labelimage.resize(totalwidth_, totalheight_);
        labelimage = 0;
        
        cells = cellimage.upperLeft() + Dist2D(2,2);
        labels = labelimage.upperLeft() + Dist2D(2,2);
    }
    
    void initCellImage(BImage & rawborders);
    int label0Cells();
    int label1Cells(int nodeCount);
    int label2Cells(BImage & rawborders);
    void labelCircles(int & nodeCount, int & edgeCount);
    void labelLine(FourEightSegmentationNeighborhoodCirculator rayAtStart, int new_label);
    void decrementLabels();
    void initNodeList(int nodeCount);
    void initEdgeList(int edgeCount);
    void initFaceList(BImage & rawborders, int faceCount);
    
    BImage cellimage;
    BImage::Iterator cells;
    IImage labelimage;
    IImage::Iterator labels;
    int width_, height_, totalwidth_, totalheight_;
    
    NodeList nodeList;
    EdgeList edgeList;
    FaceList faceList;
};

/***********************************************************************/
/*                                                                     */
/*           FourEightSegmentationNeighborhoodCirculator               */
/*                                                                     */
/***********************************************************************/

inline int FourEightSegmentationNeighborhoodCirculator::cell() const
{
    return segmentation_->cells[center()];
}

inline int FourEightSegmentationNeighborhoodCirculator::label() const 
{
    return segmentation_->labels[center()];
}

inline int FourEightSegmentationNeighborhoodCirculator::neighborCell() const 
{
    return segmentation_->cells[location()];
}

inline int FourEightSegmentationNeighborhoodCirculator::neighborLabel() const 
{
    return segmentation_->labels[location()];
}

inline int FourEightSegmentationNeighborhoodCirculator::forwardNeighborLabel() const 
{
    return segmentation_->labels[center() + neighbor.nextDiff()];
}

inline int FourEightSegmentationNeighborhoodCirculator::backwardNeighborLabel() const 
{
    return segmentation_->labels[center() + neighbor.prevDiff()];
}

inline int FourEightSegmentationNeighborhoodCirculator::forwardNeighborCell() const 
{
    return segmentation_->cells[center() + neighbor.nextDiff()];
}

inline int FourEightSegmentationNeighborhoodCirculator::backwardNeighborCell() const 
{
    return segmentation_->cells[center() + neighbor.prevDiff()];
}

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentationRayCirculator                  */
/*                                                                     */
/***********************************************************************/

inline int FourEightSegmentationRayCirculator::degree() const
{
    return segmentation()->node(nodeLabel()).degree;
}

inline float FourEightSegmentationRayCirculator::x() const 
{
    return segmentation()->node(nodeLabel()).x;
}

inline float FourEightSegmentationRayCirculator::y() const 
{
    return segmentation()->node(nodeLabel()).y;
}

/***********************************************************************/
/*                                                                     */
/*                 initFourEightSegmentationRawBorders                 */
/*                                                                     */
/***********************************************************************/

template <class SrcIter, class SrcAcc>
void initFourEightSegmentationRawBorders(SrcIter ul, SrcIter lr, SrcAcc src,
                               BImage & rawborders)
{
    int w = lr.x - ul.x;
    int h = lr.y - ul.y;
    int x,y;
    
    typedef typename SrcAcc::value_type SrcType;
    SrcType zero = NumericTraits<SrcType>::zero();
    
    rawborders = 0;
    initImageBorder(srcIterRange(rawborders.upperLeft()+Dist2D(1,1),
                                 rawborders.lowerRight()-Dist2D(1,1),
                                 rawborders.accessor()),
                                 1,1);

    for(y=0; y<h; ++y, ++ul.y)
    {
        SrcIter sx = ul;
        for(x=0; x<w; ++x, ++sx.x)
        {
            if(src(sx) == zero)  rawborders(x+2, y+2) = 1;
        }
    }
}

/***********************************************************************/
/*                                                                     */
/*                FourEightSegmentation::initCellImage                 */
/*                                                                     */
/***********************************************************************/

void 
FourEightSegmentation::initCellImage(BImage & rawborders)
{
    BImage::Iterator raw = rawborders.upperLeft() + Dist2D(1,1);
    
    int x,y;
    
    for(y=-1; y<=height_; ++y, ++raw.y)
    {
        BImage::Iterator rx = raw;
        for(x=-1; x<=width_; ++x, ++rx.x)
        {
            if(*rx == 0) 
            {
                cells(x,y) = CellConfigurationsRegion;
            }
            else
            {
                EightNeighborCoding neighbors(EightNeighborCoding::SouthEast);
                EightNeighborCoding end = neighbors;

                int conf = 0;
                do
                {
                    conf = (conf << 1) | rx[neighbors.diff()];
                }
                while(--neighbors != end);

                if(cellConfigurations[conf] == CellConfigurationsError)
                {
                    char message[200];
                    sprintf(message, "FourEightSegmentation::init(): "
                            "Configuration at (%d, %d) must be thinned further",
                            x, y);

                    precondition(0, message);
                }

                cells(x,y) = cellConfigurations[conf];
            }
        }
    }
}

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentation::label2Cells                  */
/*                                                                     */
/***********************************************************************/

int FourEightSegmentation::label2Cells(BImage & rawborders)
{
    return labelImageWithBackground(srcImageRange(rawborders), destImage(labelimage), false, 1);
}

/***********************************************************************/
/*                                                                     */
/*                FourEightSegmentation::extract0Cells                 */
/*                                                                     */
/***********************************************************************/

int FourEightSegmentation::label0Cells()
{
    BImage nodeimage(totalwidth_, totalheight_);
    BImage::Iterator nodes = nodeimage.upperLeft() + Dist2D(2,2);
    
    int x,y;
    
    for(y=-2; y<height_+2; ++y)
    {
        for(x=-2; x<width_+2; ++x)
        {
            if(cells(x,y) == CellConfigurationsVertex)   
            {
                nodes(x,y) = 1;
                
                // test for forbidden configuration
                FourEightSegmentationNeighborhoodCirculator n(this, Dist2D(x,y));
                FourEightSegmentationNeighborhoodCirculator nend = n;

                do
                {
                    if(n.neighborCell() == CellConfigurationsLine &&
                       n.forwardNeighborCell() == CellConfigurationsLine)
                    {
                        char msg[200];
                        sprintf(msg,"initFourEightSegmentation(): "
                                    "Node at (%d, %d) has two incident edgels form the same edge", 
                                x,y);
                        precondition(0, msg);
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
    
    return labelImageWithBackground(srcImageRange(nodeimage), destImage(labelimage), true, 0);
}

/***********************************************************************/
/*                                                                     */
/*                  FourEightSegmentation::labelLine                   */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::labelLine(
       FourEightSegmentationNeighborhoodCirculator rayAtStart,
       int new_label)
{
    FourEightSegmentationEdgelIterator line(rayAtStart);

    // follow the line and relabel it
    for(;!line.isEnd(); ++line)
    {
        labels[line.location()] = new_label;
    }
}

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentation::label1Cells                  */
/*                                                                     */
/***********************************************************************/

int FourEightSegmentation::label1Cells(int number_of_nodes)
{
    int x,y;
    
    std::vector<int> nodeProcessed(number_of_nodes + 1, 0);
    
    int number_of_edges = 0;
    
    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(cells(x,y) != CellConfigurationsVertex)  continue;
            if(nodeProcessed[labels(x,y)])  continue;
            
            nodeProcessed[labels(x,y)] = 1;
            
            RayCirculator rayAtStart(this, Dist2D(x,y), EightNeighborCoding::West);
            RayCirculator rayEnd = rayAtStart;
            
            do
            {
                if(rayAtStart.edgeLabel() != 0) continue;
                
                labelLine(rayAtStart.neighbor, ++number_of_edges);
            }
            while(++rayAtStart != rayEnd);
        }
    }
    
    return number_of_edges;
}

/***********************************************************************/
/*                                                                     */
/*                FourEightSegmentation::labelCircles                  */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::labelCircles(int & number_of_nodes, int & number_of_edges)
{
    int x,y;
        
    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(labels(x,y) != 0)  continue;
            
            // found a circle
            
            // mark first point as node
            cells(x,y) = CellConfigurationsVertex;
            labels(x,y) = ++number_of_nodes;
            
            FourEightSegmentationNeighborhoodCirculator rayAtStart(this, Dist2D(x,y));
            FourEightSegmentationNeighborhoodCirculator rayEnd = rayAtStart;
            
            do
            {
                if(rayAtStart.neighborCell() != CellConfigurationsLine) continue;
                if(rayAtStart.neighborLabel() != 0) continue;
                
                labelLine(rayAtStart, ++number_of_edges);
            }
            while(++rayAtStart != rayEnd);
        }
    }
}


/***********************************************************************/
/*                                                                     */
/*              FourEightSegmentation::decrementLabels                 */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::decrementLabels()
{
    IImage::ScanOrderIterator i = labelimage.begin();
    IImage::ScanOrderIterator iend = labelimage.end();
    
    for(; i != iend; ++i) --(*i);
}

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentation::initNodeList                 */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::initNodeList(int number_of_nodes)
{
    nodeList.resize(number_of_nodes, NodeInfo());
    std::vector<int> areas(number_of_nodes, 0);

    int x,y;
    
    for(y=-1; y<=height_; ++y)
    {
        for(x=-1; x<=width_; ++x)
        {
            if(cells(x,y) != CellConfigurationsVertex)  continue;
            
            int index = labels(x,y);
            
            if(nodeList[index].notInitialized())
            {
                nodeList[index].label = labels(x,y);

                nodeList[index].x = x;
                nodeList[index].y = y;
                nodeList[index].size = 1;
                nodeList[index].ray = RayCirculator(this, Dist2D(x,y), EightNeighborCoding::West);
                
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
                FourEightSegmentationNeighborhoodCirculator 
                         neighbor(this, Dist2D(x,y), EightNeighborCoding::West);
                FourEightSegmentationCrackCirculator crack(neighbor);
                FourEightSegmentationCrackCirculator crackend(crack);

                do
                {
                    areas[index] += crack.diff().width*crack.y() - 
                                    crack.diff().height*crack.x();
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
                    nodeList[i].ray.center().width, nodeList[i].ray.center().height);
            precondition(0, msg);
        }
    }
}

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentation::initEdgeList                 */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::initEdgeList(int number_of_edges)
{
    edgeList.resize(number_of_edges, EdgeInfo());
    
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
            if(edgeList[index].notInitialized())
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

/***********************************************************************/
/*                                                                     */
/*                 FourEightSegmentation::initFaceList                 */
/*                                                                     */
/***********************************************************************/

void FourEightSegmentation::initFaceList(BImage & rawborders, int number_of_faces)
{
    faceList.resize(number_of_faces, FaceInfo());    

    IImage borderlabels(totalwidth_, totalheight_);
    IImage::Iterator borderlab = borderlabels.upperLeft() + Dist2D(2,2);
    borderlabels = 0;
    int countBorderComponents = 
        labelImageWithBackground(srcImageRange(rawborders), destImage(borderlabels), true, 0);
        
    vector<int> borderProcessed(countBorderComponents + 1, 0);
    
    // process outer face
    faceList[0].label= 0;
    faceList[0].anchor = Dist2D(-2, -2);
    RayCirculator ray(this, Dist2D(-1, -1), EightNeighborCoding::West);
    --ray;
    faceList[0].contours.push_back(ContourCirculator(ray));
    borderProcessed[borderlab(-1, -1)] = 1;
    
    FaceAtLeftAccessor leftface;
                    
    int x,y;
    
    for(y=0; y<height_; ++y)
    {
        for(x=0; x<width_; ++x)
        {
            if(cells(x,y) != CellConfigurationsRegion)  continue;
            
            int index = labels(x,y);
            
            if(faceList[index].notInitialized())
            {
                faceList[index].label = index;
                faceList[index].anchor = Dist2D(x,y);
                
                // find incident node
                if(cells(x-1,y) == CellConfigurationsVertex)
                {
                    // this is the node
                    RayCirculator ray(this, Dist2D(x-1, y), EightNeighborCoding::East);
                    --ray;
                    
                    invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");
                    
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
                    
                    invariant(leftface.label(c) == index, "FourEightSegmentation::initFaceList()");
                    
                    faceList[index].contours.push_back(c);
                }
            }
            else 
            {
                // look for inner contours
                FourEightSegmentationNeighborhoodCirculator 
                           neighbor(this, Dist2D(x,y), EightNeighborCoding::East);
                FourEightSegmentationNeighborhoodCirculator nend = neighbor;
                
                do
                {
                    int bindex = borderlab[neighbor.location()];
                    if(bindex == 0 || borderProcessed[bindex]) continue;
                    
                    // found an inner contour
                    borderProcessed[bindex] = 1;
                    
                    // find incident node
                    if(cells[neighbor.location()] == CellConfigurationsVertex)
                    {
                        // this is the node
                        FourEightSegmentationNeighborhoodCirculator n = neighbor;
                        n.jumpToOpposite();
                        RayCirculator ray(n);
                        --ray;

                        invariant(leftface.label(ray) == index, "FourEightSegmentation::initFaceList()");

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

                        invariant(leftface.label(c) == index, "FourEightSegmentation::initFaceList()");

                        faceList[index].contours.push_back(c);
                    }
                    
                }
                while(++neighbor != nend);
            }
        }
    }
}

#endif /* VIGRA_FOUREIGHTSEGMENTATION_HXX */
