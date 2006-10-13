#include "cppmap.hxx"
#include <boost/python/detail/api_placeholder.hpp>
#include <vigra/tinyvector.hxx>
#include <vigra/pythonimage.hxx>
#include <vigra/pythonutil.hxx>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "exporthelpers.hxx"

const CellLabel UNINITIALIZED_CELL_LABEL =
    vigra::NumericTraits<CellLabel>::max();

class GeoMap::Node
{
  protected:
    typedef std::vector<int> DartLabels;

    GeoMap        *map_;
    CellLabel      label_;
    vigra::Vector2 position_;
    DartLabels     darts_;

    friend class Dart; // give access to darts_
    friend class GeoMap; // give access to darts_

  public:
    Node(GeoMap *map, const vigra::Vector2 &position)
    : map_(map),
      label_(map->nodes_.size()),
      position_(position)
    {
        map_->nodes_.push_back(GeoMap::Nodes::value_type(this));
        ++map_->nodeCount_;
        map_->nodeMap_.insert(PositionedNodeLabel(position_, label_));
    }

    bool initialized() const
    {
        return map_ != NULL;
    }

    void uninitialize()
    {
        GeoMap *map = map_; // local copy (prevent 2nd uninitialize() through destructor)
        map_ = NULL;
        --map->nodeCount_;
        map->nodeMap_.erase(
            map->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                                  vigra::NumericTraits<double>::epsilon()));
        RESET_PTR(map->nodes_[label_]); // may have effect like "delete this;"
#ifdef USE_INSECURE_CELL_PTRS
        delete this;
#endif
    }

    CellLabel label() const
    {
        return label_;
    }

    const vigra::Vector2 &position() const
    {
        return position_;
    }

    void setPosition(const vigra::Vector2 &p);

    inline Dart anchor() const;

    unsigned int degree() const
    {
        return darts_.size();
    }

    inline bool operator==(const GeoMap::Node &other)
    {
        return label() == other.label() && map_ == other.map_;
    }

    inline bool operator!=(const GeoMap::Node &other)
    {
        return !operator==(other);
    }

    GeoMap *map() const
    {
        return map_;
    }

  private:
    Node(const Node &) {} // disallow copying
    Node &operator=(const Node &) { return *this; }
};

class GeoMap::Edge
: public vigra::BBoxPolygon<vigra::Vector2>
{
  public:
    typedef vigra::BBoxPolygon<vigra::Vector2> Base;

  protected:
    GeoMap   *map_;
    CellLabel label_;
    CellLabel startNodeLabel_, endNodeLabel_;
    CellLabel leftFaceLabel_, rightFaceLabel_;
    int       protection_;

    friend class Dart; // allow setLeftFaceLabel

  public:
    template<class POINTS>
    Edge(GeoMap *map, CellLabel startNodeLabel, CellLabel endNodeLabel,
         const POINTS &p)
    : Base(p),
      map_(map),
      label_(map->edges_.size()),
      startNodeLabel_(startNodeLabel),
      endNodeLabel_(endNodeLabel),
      leftFaceLabel_(UNINITIALIZED_CELL_LABEL),
      rightFaceLabel_(UNINITIALIZED_CELL_LABEL),
      protection_(0)
    {
        map_->edges_.push_back(GeoMap::Edges::value_type(this));
        ++map_->edgeCount_;
    }

    bool initialized() const
    {
        return map_ != NULL;
    }

    void uninitialize()
    {
        GeoMap *map = map_; // local copy (prevent 2nd uninitialize() through destructor)
        map_ = NULL;
        --map->edgeCount_;
        RESET_PTR(map->edges_[label_]); // may have effect like "delete this;"
#ifdef USE_INSECURE_CELL_PTRS
        delete this;
#endif
    }

    CellLabel label() const
    {
        return label_;
    }

    inline Dart dart() const;

    CellLabel startNodeLabel() const
    {
        return startNodeLabel_;
    }

    GeoMap::Nodes::value_type startNode() const
    {
        vigra_precondition(initialized(), "startNode() of uninitialized edge!");
        return map_->node(startNodeLabel_);
    }

    CellLabel endNodeLabel() const
    {
        return endNodeLabel_;
    }

    GeoMap::Nodes::value_type endNode() const
    {
        vigra_precondition(initialized(), "endNode() of uninitialized edge!");
        return map_->node(endNodeLabel_);
    }

    CellLabel leftFaceLabel() const
    {
        return leftFaceLabel_;
    }

    GeoMap::Faces::value_type leftFace() const
    {
        vigra_precondition(initialized(), "leftFace() of uninitialized edge!");
        return map_->face(leftFaceLabel_);
    }

    CellLabel rightFaceLabel() const
    {
        return rightFaceLabel_;
    }

    GeoMap::Faces::value_type rightFace() const
    {
        vigra_precondition(initialized(), "rightFace() of uninitialized edge!");
        return map_->face(rightFaceLabel_);
    }

    bool isBridge() const
    {
        return leftFaceLabel_ == rightFaceLabel_;
    }

    bool isLoop() const
    {
        return startNodeLabel_ == endNodeLabel_;
    }

    inline bool operator==(const GeoMap::Edge &other)
    {
        return label() == other.label() && map_ == other.map_;
    }

    inline bool operator!=(const GeoMap::Edge &other)
    {
        return !operator==(other);
    }

    GeoMap *map() const
    {
        return map_;
    }

  private:
    Edge(const Edge &) {} // disallow copying
    Edge &operator=(const Edge &) { return *this; }
};

class DartPointIter
{
    CELL_PTR(GeoMap::Edge) edge_;
    int index_, inc_, end_;

  public:
        /** the iterator's value type
        */
    typedef GeoMap::Edge::value_type value_type;

        /** the iterator's reference type (return type of <tt>*iter</tt>)
        */
    typedef value_type & reference;

        /** the iterator's pointer type (return type of <tt>operator-></tt>)
        */
    typedef value_type * pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    DartPointIter(GeoMap::Dart const &dart);

    DartPointIter & operator++()
    {
        index_ += inc_;
        return *this;
    }

    DartPointIter operator++(int)
    {
        DartPointIter ret(*this);
        operator++();
        return ret;
    }

    /**
     * the opposite of inRange(); true if this iterator is behind the
     * range and should not be dereferenced any more
     */
    bool atEnd() const
    {
        return index_ == end_;
    }

    /**
     * the opposite of atEnd(); true if this iterator is dereferencable
     */
    bool inRange() const
    {
        return index_ != end_;
    }

    reference operator*() const
    {
        return (*edge_)[index_];
    }

    pointer operator->() const
    {
        return &(operator*());
    }
};

class GeoMap::Dart
{
  protected:
    GeoMap *map_;
    int     label_;

    void setLeftFaceLabel(CellLabel label)
    {
        if(label_ > 0)
            guaranteedEdge()->leftFaceLabel_ = label;
        else
            guaranteedEdge()->rightFaceLabel_ = label;
    }

    friend class Face; // allow setLeftFaceLabel in Face constructor

  public:
    Dart(GeoMap *map, int label)
    : map_(map),
      label_(label)
    {}

    Dart clone() const
    {
        return Dart(map_, label_);
    }

    int label() const
    {
        return label_;
    }

    CellLabel edgeLabel() const
    {
        return abs(label_);
    }

    CellLabel startNodeLabel() const
    {
        if(label_ > 0)
            return guaranteedEdge()->startNodeLabel();
        else
            return guaranteedEdge()->endNodeLabel();
    }

    CellLabel endNodeLabel() const
    {
        if(label_ > 0)
            return guaranteedEdge()->endNodeLabel();
        else
            return guaranteedEdge()->startNodeLabel();
    }

//         def _setStartNode(self, node):
//             """changes corresponding start/end node of this dart's
//             edge and the first/last point of its' polygon, too"""
//             if self._label > 0:
//                 self.edge()._startNodeLabel = node._label
//                 self.edge()[0] = node.position()
//                 #self.edge().invalidateProperties()
//             else:
//                 self.edge()._endNodeLabel = node._label
//                 self.edge()[-1] = node.position()
//                 #self.edge().invalidateProperties()

    CellLabel leftFaceLabel() const
    {
        if(label_ > 0)
            return guaranteedEdge()->leftFaceLabel();
        else
            return guaranteedEdge()->rightFaceLabel();
    }

    CellLabel rightFaceLabel() const
    {
        if(label_ > 0)
            return guaranteedEdge()->rightFaceLabel();
        else
            return guaranteedEdge()->leftFaceLabel();
    }


    GeoMap::Edges::value_type edge() const
    {
        return map_->edge(edgeLabel());
    }

    GeoMap::Edges::value_type guaranteedEdge() const
    {
        GeoMap::Edges::value_type result(edge());
        vigra_precondition(result,
            "Cannot operate on invalid dart belonging to removed edge!");
        return result;
    }

    GeoMap::Nodes::value_type startNode() const
    {
        return map_->node(startNodeLabel());
    }

    GeoMap::Nodes::value_type endNode() const
    {
        return map_->node(endNodeLabel());
    }

    GeoMap::Faces::value_type leftFace() const
    {
        return map_->face(leftFaceLabel());
    }

    GeoMap::Faces::value_type rightFace() const
    {
        return map_->face(rightFaceLabel());
    }

    double partialArea() const
    {
        if(label_ > 0)
            return guaranteedEdge()->partialArea();
        else
            return -guaranteedEdge()->partialArea();
    }

    DartPointIter pointIter() const
    {
        return DartPointIter(*this);
    }

    typedef GeoMap::Edge::value_type value_type;

    const value_type &operator[](int index) const
    {
        if(label_ > 0)
            return (*guaranteedEdge())[index];
        else
            return (*guaranteedEdge())[size()-1-index];
    }

    GeoMap::Edge::size_type size() const
    {
        return guaranteedEdge()->size();
    }

    Dart &nextAlpha()
    {
        label_ = -label_;
        return *this;
    }

    Dart &nextSigma(int times = 1)
    {
        Node::DartLabels &darts(startNode()->darts_);
        int i = 0;
        for(; i < (int)darts.size(); ++i)
            if(darts[i] == label_)
                break;
        vigra_precondition(i < (int)darts.size(),
                           "Dart not attached to its startnode??");
        i = (i + times) % (int)darts.size();
        if(i < 0)
            i += darts.size();
        label_ = darts[i];
        return *this;
    }

    Dart &prevSigma(int times = 1)
    {
        return nextSigma(-times);
    }

    Dart &nextPhi()
    {
        return nextAlpha().prevSigma();
    }

    Dart &prevPhi()
    {
        return nextSigma().nextAlpha();
    }

    bool operator==(const Dart &other) const
    {
        return label_ == other.label_;
    }

    bool operator!=(const Dart &other) const
    {
        return label_ != other.label_;
    }

//         def alphaOrbit(self):
//             return self._orbit("nextAlpha")

//         def sigmaOrbit(self):
//             return self._orbit("nextSigma")

//         def phiOrbit(self):
//             return self._orbit("nextPhi")

//         def _orbit(self, opName):
//             dart = self.clone()
//             op = getattr(dart, opName)
//             while True:
//                 yield dart.clone()
//                 op()
//                 if dart == self:
//                     break

};

DartPointIter::DartPointIter(GeoMap::Dart const &dart)
: edge_(dart.guaranteedEdge())
{
    if(dart.label() > 0)
    {
        index_ = 0;
        inc_ = 1;
        end_ = dart.size();
    }
    else
    {
        index_ = dart.size() - 1;
        inc_ = -1;
        end_ = -1;
    }
}

/*
 * Note: This code is based on the assumption that a dart must always
 * have a least two points!
 */
class ContourPointIter
{
    DartPointIter dpi_;
    GeoMap::Dart dart_, end_;

  public:
        /** the iterator's value type
        */
    typedef GeoMap::Edge::value_type value_type;

        /** the iterator's reference type (return type of <tt>*iter</tt>)
        */
    typedef value_type & reference;

        /** the iterator's pointer type (return type of <tt>operator-></tt>)
        */
    typedef value_type * pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    ContourPointIter(GeoMap::Dart const &dart, bool firstTwice = false)
    : dpi_(dart),
      dart_(dart),
      end_(dart)
    {
        if(!firstTwice)
            ++dpi_;
    }

    ContourPointIter & operator++()
    {
        ++dpi_;
        if(dpi_.atEnd())
        {
            if(dart_.nextPhi() != end_)
            {
                dpi_ = DartPointIter(dart_);
                ++dpi_;
            }
        }
        return *this;
    }

    ContourPointIter operator++(int)
    {
        ContourPointIter ret(*this);
        operator++();
        return ret;
    }

    /**
     * the opposite of inRange(); true if this iterator is behind the
     * range and should not be dereferenced any more
     */
    bool atEnd() const
    {
        return dpi_.atEnd();
    }

    /**
     * the opposite of atEnd(); true if this iterator is dereferencable
     */
    bool inRange() const
    {
        return dpi_.inRange();
    }

    reference operator*() const
    {
        return *dpi_;
    }

    pointer operator->() const
    {
        return &(operator*());
    }
};

double contourArea(const GeoMap::Dart &dart)
{
    double result = 0.0;
    GeoMap::Dart d(dart);
    do
    {
        if(!d.guaranteedEdge()->isBridge())
            result += d.partialArea();
    }
    while(d.nextPhi() != dart);
    return result;
}

Polygon contourPoly(const GeoMap::Dart &dart)
{
    Polygon result;
    GeoMap::Dart d(dart);
    do
    {
        if(d.label() < 0)
        {
            Polygon rev(*d.edge());
            rev.reverse();
            result.extend(rev);
        }
        else
            result.extend(*d.edge());
    }
    while(d.nextPhi() != dart);
    return result;
}

typedef vigra::MultiArray<2, int>::difference_type IVector2;

inline IVector2 intVPos(const Vector2 &p)
{
    return IVector2((int)floor(p[0]+0.5), (int)floor(p[1]+0.5));
}

class GeoMap::Face
{
  public:
    typedef Edge::BoundingBox BoundingBox;

  protected:
    GeoMap             *map_;
    CellLabel           label_;
    std::vector<Dart>   anchors_;
    mutable BoundingBox boundingBox_;
    mutable bool        boundingBoxValid_;
    mutable double      area_;
    mutable bool        areaValid_;
    int                 pixelArea_;

    friend class GeoMap; // give access to pixelArea_

  public:
    Face(GeoMap *map, Dart anchor)
    : map_(map),
      label_(map->faces_.size()),
      boundingBoxValid_(false),
      areaValid_(false),
      pixelArea_(0)
    {
        map_->faces_.push_back(GeoMap::Faces::value_type(this));
        ++map_->faceCount_;

        if(label_)
        {
            anchors_.push_back(anchor);

            for(; anchor.leftFaceLabel() == UNINITIALIZED_CELL_LABEL;
                anchor.nextPhi())
            {
                // don't calculate area on-the-fly here; we want to
                // exclude bridges from the area!
                anchor.setLeftFaceLabel(label_);
            }
        }
    }

    bool initialized() const
    {
        return map_ != NULL;
    }

    void uninitialize()
    {
        GeoMap *map = map_; // local copy (prevent 2nd uninitialize() through destructor)
        map_ = NULL;
        --map->faceCount_;
        RESET_PTR(map->faces_[label_]); // may have effect like "delete this;"
#ifdef USE_INSECURE_CELL_PTRS
        delete this;
#endif
    }

    CellLabel label() const
    {
        return label_;
    }

    const BoundingBox &boundingBox() const
    {
        vigra_precondition(label_, "infinite face has no boundingBox()!");

        if(!boundingBoxValid_)
        {
            boundingBox_ = BoundingBox();
            Dart anchor(anchors_[0]), dart(anchor);
            do
            {
                boundingBox_ |= dart.edge()->boundingBox();
            }
            while(dart.nextPhi() != anchor);
            boundingBoxValid_ = true;
        }
        return boundingBox_;
    }

    bool contains(const Vector2 &point) const
    {
        vigra_precondition(initialized(), "contains() of uninitialized face!");
        if(map_->labelImage_ && (*map_->labelImage_)[intVPos(point)] == (int)label_)
            return true;
        unsigned int i = 0;
        if(label_)
        {
            if(!boundingBox().contains(point))
                return false;
            if(!contourPoly(anchors_[0]).contains(point))
                return false;
            ++i;
        }
        for(; i < anchors_.size(); ++i)
            if(contourPoly(anchors_[i]).contains(point))
                return false;
        return true;
    }

    double area() const
    {
        if(!areaValid_)
        {
            area_ = 0.0;
            for(unsigned int i = 0; i < anchors_.size(); ++i)
            {
                area_ += contourArea(anchors_[i]);
            }
            areaValid_ = true;
        }
        return area_;
    }

    const Dart &contour(unsigned int index = 0)
    {
        return anchors_[index];
    }

    void embedContour(const Dart &anchor)
    {
        anchors_.push_back(anchor);

        Dart dart(anchor); // we need a non-const reference
        for(; dart.leftFaceLabel() != label_; dart.nextPhi())
            dart.setLeftFaceLabel(label_);

        if(areaValid_)
            area_ += contourArea(dart);

        vigra_postcondition(dart == anchor,
                            "contour labeled partially?!");
    }

    inline bool operator==(const GeoMap::Face &other)
    {
        return label() == other.label() && map_ == other.map_;
    }

    inline bool operator!=(const GeoMap::Face &other)
    {
        return !operator==(other);
    }

    GeoMap *map() const
    {
        return map_;
    }

  private:
    Face(const Face &) {} // disallow copying
    Face &operator=(const Face &) { return *this; }
};

void GeoMap::Node::setPosition(const vigra::Vector2 &p)
{
    vigra_precondition(initialized(), "setPosition() of uninitialized node!");
    map_->nodeMap_.erase(
        map_->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                               vigra::NumericTraits<double>::epsilon()));
    position_ = p;
    for(unsigned int i = 0; i < darts_.size(); ++i)
    {
        if(i > 0)
        {
            (*map_->edge( i))[ 0] = p;
        }
        else
        {
            GeoMap::Edge &edge(*map_->edge(-i));
            edge[edge.size()-1] = p;
        }
    }
    map_->nodeMap_.insert(PositionedNodeLabel(p, label_));
}

inline GeoMap::Dart GeoMap::Node::anchor() const
{
    vigra_precondition(initialized(), "anchor() of uninitialized node!");
    return Dart(map_, darts_[0]);
}

GeoMap::GeoMap(bp::list nodePositions,
               bp::list edgeTuples, vigra::Size2D imageSize)
: nodeCount_(0),
  edgeCount_(0),
  faceCount_(0),
  imageSize_(imageSize),
  labelImage_(NULL)
{
    if(len(nodePositions))
        std::cerr << "initalizing nodes...\n";
    nodes_.push_back(NULL_PTR(Node));
    for(int i = 1; i < len(nodePositions); ++i)
    {
        bp::extract<Vector2> ve(nodePositions[i]);
        if(ve.check())
            addNode(ve());
        else
            nodes_.push_back(NULL_PTR(Node));
    }

    if(len(edgeTuples))
        std::cerr << "initalizing edges...\n";
    edges_.push_back(NULL_PTR(Edge));
    for(int i = 1; i < len(edgeTuples); ++i)
    {
        bp::extract<bp::tuple> ete(edgeTuples[i]);
        if(ete.check())
        {
            bp::tuple edgeTuple(ete());
            bp::extract<Vector2Array> pe(edgeTuple[2]);
            if(!pe.check())
            {
                std::cerr << "why, oh why, do I have to die??\n";
                bp::throw_type_error(
                    "GeoMap.__init__: edge geometry not convertable to Vector2Array");
            }
            CellLabel startNodeLabel = bp::extract<CellLabel>(edgeTuple[0])();
            CellLabel endNodeLabel   = bp::extract<CellLabel>(edgeTuple[1])();
            addEdge(startNodeLabel, endNodeLabel, pe());
        }
        else
            edges_.push_back(NULL_PTR(Edge));
    }

    if(edgeCount_)
    {
        std::cerr << "initalizing sigma orbits...\n";
        sortEdgesDirectly();
        std::cerr << "initalizing contours...\n";
        initContours();
        std::cerr << "initalizing faces (" << faceCount() << " contours)...\n";
        embedFaces();
    }
}

GeoMap::~GeoMap()
{
    // make sure the cells' destructors don't access this map!
    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
}

inline GeoMap::Dart GeoMap::dart(int label)
{
    return GeoMap::Dart(this, label);
}

double angleTheta(double dy, double dx); // implemented in polygon.cxx

CELL_PTR(GeoMap::Face) GeoMap::faceAt(const vigra::Vector2 &position)
{
    if(labelImage_)
    {
        GeoMap::LabelImage::difference_type p(intVPos(position));
        if(labelImage_->isInside(p))
        {
            int faceLabel = (*labelImage_)[p];
            if(faceLabel > 0)
                return face(faceLabel);
        }
    }

    FaceIterator it = facesBegin();
    for(++it; it.inRange(); ++it)
        if((*it)->contains(position))
            return *it;

    return face(0);
}

CELL_PTR(GeoMap::Node) GeoMap::addNode(
    const vigra::Vector2 &position)
{
    GeoMap::Node *result = new GeoMap::Node(this, position);
    return node(result->label());
}

CELL_PTR(GeoMap::Edge) GeoMap::addEdge(
    CellLabel startNodeLabel, CellLabel endNodeLabel,
    const Vector2Array &points, CellLabel label)
{
    if(label > edges_.size())
        edges_.resize(label, NULL_PTR(GeoMap::Edge));
    CELL_PTR(GeoMap::Node)
        startNode = node(startNodeLabel),
        endNode = node(endNodeLabel);
    vigra_precondition(
        startNode && endNode, "addEdge(): invalid start- or endNodeLabel!");
    GeoMap::Edge *result = new GeoMap::Edge(
        this, startNodeLabel, endNodeLabel, points);
    startNode->darts_.push_back( (int)result->label());
    endNode  ->darts_.push_back(-(int)result->label());
    return edge(result->label());
}

void GeoMap::sortEdgesDirectly()
{
    typedef std::pair<double, int> DartAngle;

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        std::vector<DartAngle> dartAngles;

        GeoMap::Node::DartLabels &dartLabels((*it)->darts_);

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
        {
            GeoMap::Dart d(dart(dartLabels[i]));
            vigra_precondition(
                d.size() >= 2, "cannot measure angle of darts with < 2 points!");
            dartAngles.push_back(
                DartAngle(angleTheta(-d[1][1] + d[0][1],
                                      d[1][0] - d[0][0]),
                          dartLabels[i]));
        }

        std::sort(dartAngles.begin(), dartAngles.end());

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
        {
            dartLabels[i] = dartAngles[i].second;
        }
    }
}

class DartPosition
{
  public:
    DartPosition(const GeoMap::Dart &dart)
    : dart_(dart),
      hitEnd_(false),
      pointIter_(DartPointIter(dart)),
      segmentIndex_(0),
      arcLength_(0.0),
      partialLength_(0.0),
      position_(*pointIter_)
    {
        p1_ = *pointIter_;
        p2_ = *++pointIter_;
    }

    bool atEnd() const
    {
        return hitEnd_;
    }

    const vigra::Vector2 &operator()() const
    {
        return position_;
    }

    void leaveCircle(const vigra::Vector2 &center, double radius2)
    {
        while((p2_ - center).squaredMagnitude() < radius2)
            if(!nextSegment())
                break;

        position_ = p2_;
    }

    void intersectCircle(const vigra::Vector2 &center, double radius2)
    {
        while((p2_ - center).squaredMagnitude() < radius2)
        {
            if(!nextSegment())
            {
                partialLength_ = 0.0;
                position_ = p2_;
                return;
            }
        }

        vigra::Vector2 diff(p2_ - p1_);
        double dist2 = diff.squaredMagnitude();
        double lambda = (
            (std::sqrt(radius2 * dist2
                       - vigra::sq(p2_[0]*p1_[1] - p1_[0]*p2_[1]
                                   + center[0]*diff[1] - diff[0]*center[1]))
             - dot(diff, p1_ - center))
            / dist2);
        diff *= lambda;
        partialLength_ += diff.magnitude();
        position_ = p1_ + diff;
    }

    int dartLabel() const
    {
        return dart_.label();
    }

  protected:
    bool nextSegment()
    {
        vigra_precondition(!hitEnd_, "DartPosition: trying to proceed from end");
        arcLength_ += (p2_ - p1_).magnitude();
        p1_ = p2_;
        ++pointIter_;
        if(pointIter_.atEnd())
        {
            hitEnd_ = true;
            return false;
        }
        p2_ = *pointIter_;
        ++segmentIndex_;
        return true;
    }

    GeoMap::Dart dart_;
    bool hitEnd_;
    DartPointIter pointIter_;
    unsigned int segmentIndex_;
    double arcLength_, partialLength_;
    vigra::Vector2 p1_, p2_, position_;
};

struct DartPositionAngle
{
    DartPosition dp;
    double absAngle, angle;

    DartPositionAngle(const GeoMap::Dart &dart)
    : dp(dart)
    {}

    bool operator<(const DartPositionAngle &other) const
    {
        return angle < other.angle;
    }
};

typedef std::vector<DartPositionAngle> DartPositionAngles;
typedef DartPositionAngles::iterator DPAI;

template<class Iterator>
void rotateArray(Iterator begin, Iterator newBegin, Iterator end)
{
    typedef std::vector<typename Iterator::value_type> TempArray;
    TempArray temp(begin, end);
    typename Iterator::difference_type pos(newBegin - begin);
    std::copy(temp.begin() + pos, temp.end(), begin);
    std::copy(temp.begin(), temp.begin() + pos, begin + (end - newBegin));
}

inline double normAngle(double diff)
{
    if(diff < -M_PI)
        diff += 2*M_PI;
    if(diff >= M_PI)
        diff -= 2*M_PI;
    return diff;
}

void sortEdgesInternal(const vigra::Vector2 &currentPos,
                       double referenceAngle,
                       DPAI dpBegin, DPAI dpEnd,
                       double stepDist2, double minAngle)
{
    if(dpEnd - dpBegin < 2)
        return;

    bool unsortableState = true;
    for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
    {
        if(!dpi->dp.atEnd())
        {
            unsortableState = false;
            dpi->dp.intersectCircle(currentPos, stepDist2);
            //dpi->dp.leaveCircle(currentPos, stepDist2);
        }

        dpi->absAngle =
            std::atan2(-dpi->dp()[1] + currentPos[1],
                        dpi->dp()[0] - currentPos[0]);

        dpi->angle = normAngle(dpi->absAngle - referenceAngle);
    }

    if(unsortableState)
    {
        vigra_fail("Unsortable group of edges occured and not handled yet!");
        return;
    }

    std::sort(dpBegin, dpEnd);

    // handle cyclicity of array first (by rotation if necessary):
    DPAI firstGroupStart = dpEnd;
    bool needRotation = false;
    while((--firstGroupStart)->angle + minAngle < dpBegin->angle)
    {
        needRotation = true;
        if(firstGroupStart == dpBegin)
        {
            needRotation = false; // whole array unsortable
            break;
        }
    }

    if(needRotation)
        rotateArray(dpBegin, firstGroupStart, dpEnd);

    // look for groups of parallel edges
    DPAI groupStart = dpBegin,
          groupLast = groupStart, // for convenience, this is always groupEnd - 1
           groupEnd = groupLast + 1;
    for(; true; ++groupLast, ++groupEnd)
    {
        // group ending?
        if((groupEnd == dpEnd) || // last group
           (groupEnd->angle >= groupLast->angle + minAngle)) // decision here
        {
            // recursion needed if > one dart in group:
            if(groupLast != groupStart)
            {
                // determine mean position of dart positions in subgroup:
                vigra::Vector2 meanPos(groupLast->dp());
                for(DPAI dpi = groupStart; dpi != groupLast; ++dpi)
                    meanPos += dpi->dp();
                meanPos /= (groupEnd - groupStart);

                // sort subgroup recursively:
                sortEdgesInternal(meanPos, normAngle(
                                      groupStart->absAngle +
                                      normAngle(groupLast->absAngle -
                                                groupStart->absAngle) / 2),
                                  groupStart, groupEnd,
                                  stepDist2, minAngle);
            }

            if(groupEnd == dpEnd)
                break; // loop end

            groupStart = groupEnd;
        }
    }
}

void GeoMap::sortEdgesEventually(double stepDist, double minDist)
{
    double minAngle = std::atan2(minDist, stepDist),
          stepDist2 = vigra::sq(stepDist);

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        DartPositionAngles dartPositions;

        GeoMap::Node::DartLabels &dartLabels((*it)->darts_);

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
            dartPositions.push_back(
                DartPositionAngle(dart(dartLabels[i])));

        sortEdgesInternal((*it)->position(), 0.0,
                          dartPositions.begin(), dartPositions.end(),
                          stepDist2, minAngle);

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
        {
            dartLabels[i] = dartPositions[i].dp.dartLabel();
        }
    }
}

void GeoMap::initContours()
{
    vigra_precondition(faces_.size() == 0,
        "initContours() called when contours were already present");

    new Face(this, Dart(this, 0)); // create infinite face, dart will be ignored

    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        if((*it)->leftFaceLabel() == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart( (*it)->label()));
        if((*it)->rightFaceLabel() == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart(-(*it)->label()));
    }
}

struct AbsAreaCompare
{
        // FIXME: actually, const pointers would suffice:
    bool operator()(CELL_PTR(GeoMap::Face) f1, CELL_PTR(GeoMap::Face) f2) const
    {
        double a1 = f1->area(), a2 = f2->area();
        double absdiff = fabs(a1) - fabs(a2);
        if(fabs(absdiff) < 1e-2 && ((a1 < 0) != (a2 < 0)))
            return (a1 < 0); // for faces with equal area, prefer the exterior one
        return absdiff > 0; // else, prefer face with larger absolute area
    }
};

void GeoMap::embedFaces(bool initLabelImage)
{
    vigra_precondition(!labelImage_,
        "embedFaces() called with already-initialized labelImage");

    if(initLabelImage)
        labelImage_ = new LabelImage(
            LabelImage::size_type(imageSize().width(), imageSize().height()), 0);

    // copy and remove all preliminary contours except the infinite one:
    GeoMap::Faces contours(faces_.begin() + 1, faces_.end());
    std::sort(contours.begin(), contours.end(), AbsAreaCompare());
    std::fill(faces_.begin() + 1, faces_.end(), NULL_PTR(Face));

    for(unsigned int i = 0; i < contours.size(); ++i)
    {
        GeoMap::Face &contour(*contours[i]); // FIXME: const

        GeoMap::Dart anchor(contour.contour(0));

        bool isExterior = contour.area() <= 0;

        if(!isExterior)
        {
            faces_[contour.label()] = contours[i];

            if(initLabelImage)
            {
                std::auto_ptr<vigra::Scanlines> scanlines(
                    scanPoly(contourPoly(anchor), imageSize().height()));
                contour.pixelArea_ =
                    fillScannedPoly(*scanlines, (int)contour.label(),
                                    labelImage_->traverser_begin(),
                                    labelImage_->size(),
                                    vigra::StandardValueAccessor<int>());
                drawScannedPoly(*scanlines, -1,
                                labelImage_->traverser_begin(),
                                labelImage_->size(),
                                vigra::StandardValueAccessor<int>());
            }
        }
        else
        {
            // contour is a hole, determine parent face
            CELL_PTR(GeoMap::Face) parent = NULL_PTR(GeoMap::Face);

            if(initLabelImage)
            {
                ContourPointIter cpi(anchor);
                while(cpi.inRange())
                {
                    GeoMap::LabelImage::difference_type p(intVPos(*cpi++));
                    if(labelImage_->isInside(p))
                    {
                        int parentLabel = (*labelImage_)[p];
                        if(parentLabel >= 0)
                        {
                            parent = face(parentLabel);
                            break;
                        }
                    }
                }
            }

            if(!parent)
            {
                ContourPointIter cpi(anchor);
                while(cpi.inRange())
                {
                    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
                    {
                        if((*it)->contains(*cpi++))
                        {
                            parent = *it;
                            goto parent_found; // double break
                        }
                    }
                }
            }

            if(!parent)
            {
                parent = face(0);
//                 vigra_postcondition(
//                     parent->contains(anchor[0]),
//                     "contour could not be embedded (parent not found)");
            }

        parent_found:
            std::cerr << "  embedding contour " << contour.label() << " in face "
                      << parent->label() << "\n";
            parent->embedContour(anchor);
            contour.uninitialize();
        }
    }
}

CELL_PTR(GeoMap::Node) GeoMap::nearestNode(
    const vigra::Vector2 &position,
    double maxSquaredDist)
{
    NodeMap::iterator n(
        nodeMap_.nearest(PositionedNodeLabel(position, 0), maxSquaredDist));
    if(n != nodeMap_.end())
        return node(n->second.payload);
    return NULL_PTR(GeoMap::Node);
}

bool GeoMap::checkConsistency()
{
    bool result = true;
    std::cerr << "GeoMap[" << this << "].checkConsistency()\n";
    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if((*it)->map() != this)
        {
            std::cerr << "  Node " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }
    }
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        if((*it)->map() != this)
        {
            std::cerr << "  Edge " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }
    }
    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
    {
        if((*it)->map() != this)
        {
            std::cerr << "  Face " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }
    }
    return result;
}

/********************************************************************/

template<class Iterator>
struct RangeIterWrapper
: bp::class_<Iterator>
{
    RangeIterWrapper(const char *name)
    : bp::class_<Iterator>(name, bp::no_init)
    {
        def("__iter__", (Iterator &(*)(Iterator &))&returnSelf,
            bp::return_internal_reference<>());
        def("next", &nextIterPos, CELL_RETURN_POLICY());
    }

    static Iterator &returnSelf(Iterator &v)
    {
        return v;
    }

    static typename Iterator::value_type nextIterPos(Iterator &v)
    {
        if(!v.inRange())
        {
            PyErr_SetString(PyExc_StopIteration, "cells iterator exhausted");
            bp::throw_error_already_set();
        }
        return *v++;
    }

};

using namespace boost::python;

void defMap()
{
    {
        CELL_RETURN_POLICY crp;

        scope geoMap(
            class_<GeoMap, boost::noncopyable>(
                "GeoMap", init<bp::list, bp::list, vigra::Size2D>())
            .def("node", &GeoMap::node, crp)
            .def("nodeIter", &GeoMap::nodesBegin)
            .def("edge", &GeoMap::edge, crp)
            .def("edgeIter", &GeoMap::edgesBegin)
            .def("face", &GeoMap::face, crp)
            .def("faceIter", &GeoMap::facesBegin)
            .def("dart", &GeoMap::dart)
            .add_property("nodeCount", &GeoMap::nodeCount)
            .add_property("edgeCount", &GeoMap::edgeCount)
            .add_property("faceCount", &GeoMap::faceCount)
            .def("maxNodeLabel", &GeoMap::maxNodeLabel,
                 "Returns an upper bound on the node labels.\n"
                 "Actually, this is the max. node label + 1, so that you can use it as LUT size.")
            .def("maxEdgeLabel", &GeoMap::maxEdgeLabel,
                 "Returns an upper bound on the edge labels.\n"
                 "Actually, this is the max. edge label + 1, so that you can use it as LUT size.")
            .def("maxFaceLabel", &GeoMap::maxFaceLabel,
                 "Returns an upper bound on the face labels.\n"
                 "Actually, this is the max. face label + 1, so that you can use it as LUT size.")
            .def("imageSize", &GeoMap::imageSize,
                 return_value_policy<copy_const_reference>())
            .def("addNode", &GeoMap::addNode, crp,
                 args("position"))
            .def("addEdge", &GeoMap::addEdge, crp,
                 (arg("startNodeLabel"), arg("endNodeLabel"), arg("points"),
                  arg("label") = 0))
            .def("sortEdgesDirectly", &GeoMap::sortEdgesDirectly)
            .def("sortEdgesEventually", &GeoMap::sortEdgesEventually,
                 args("stepDist", "minDist"))
            .def("initContours", &GeoMap::initContours)
            .def("embedFaces", &GeoMap::embedFaces, (arg("initLabelImage") = true))
            .def("nearestNode", &GeoMap::nearestNode, crp,
                 (arg("position"), arg(
                     "maxSquaredDist") = vigra::NumericTraits<double>::max()))
            .def("checkConsistency", &GeoMap::checkConsistency)
            );

        RangeIterWrapper<GeoMap::NodeIterator>("NodeIterator");
        RangeIterWrapper<GeoMap::EdgeIterator>("EdgeIterator");
        RangeIterWrapper<GeoMap::FaceIterator>("FaceIterator");

        class_<GeoMap::Node, boost::noncopyable>("Node", init<GeoMap *, const vigra::Vector2 &>())
            .def("initialized", &GeoMap::Node::initialized)
            .def("label", &GeoMap::Node::label)
            .def("position", &GeoMap::Node::position,
                 return_value_policy<copy_const_reference>())
            .def("setPosition", &GeoMap::Node::setPosition)
            .def("degree", &GeoMap::Node::degree)
            .def("anchor", &GeoMap::Node::anchor)
            .def(self == self)
            .def(self != self)
        ;

        class_<GeoMap::Edge, bases<Polygon>, boost::noncopyable>("Edge", no_init)
            .def(init<GeoMap *, CellLabel, CellLabel, GeoMap::Edge::Base>())
            .def("initialized", &GeoMap::Edge::initialized)
            .def("label", &GeoMap::Edge::label)
            .def("startNodeLabel", &GeoMap::Edge::startNodeLabel)
            .def("startNode", &GeoMap::Edge::startNode, crp)
            .def("endNodeLabel", &GeoMap::Edge::endNodeLabel)
            .def("endNode", &GeoMap::Edge::endNode, crp)
            .def("leftFaceLabel", &GeoMap::Edge::leftFaceLabel)
            .def("leftFace", &GeoMap::Edge::leftFace, crp)
            .def("rightFaceLabel", &GeoMap::Edge::rightFaceLabel)
            .def("rightFace", &GeoMap::Edge::rightFace, crp)
            .def("isBridge", &GeoMap::Edge::isBridge)
            .def("isLoop", &GeoMap::Edge::isLoop)
            .def(self == self)
            .def(self != self)
        ;

        class_<GeoMap::Face, boost::noncopyable>("Face", init<GeoMap *, GeoMap::Dart>())
            .def("initialized", &GeoMap::Face::initialized)
            .def("label", &GeoMap::Face::label)
            .def("boundingBox", &GeoMap::Face::boundingBox,
                 return_value_policy<copy_const_reference>())
            .def("contains", &GeoMap::Face::contains)
            .def("area", &GeoMap::Face::area)
            .def("contour", &GeoMap::Face::contour,
                 return_internal_reference<>())
            .def(self == self)
            .def(self != self)
        ;

        return_internal_reference<> rself; // "return self" policy

        RangeIterWrapper<DartPointIter>("DartPointIter");

        class_<GeoMap::Dart>("Dart", no_init)
            .def(init<GeoMap *, int>())
            .def("clone", &GeoMap::Dart::clone)
            .def("label", &GeoMap::Dart::label)
            .def("edgeLabel", &GeoMap::Dart::edgeLabel)
            .def("edge", &GeoMap::Dart::edge, crp)
            .def("guaranteedEdge", &GeoMap::Dart::guaranteedEdge, crp)
            .def("startNodeLabel", &GeoMap::Dart::startNodeLabel)
            .def("startNode", &GeoMap::Dart::startNode, crp)
            .def("endNodeLabel", &GeoMap::Dart::endNodeLabel)
            .def("endNode", &GeoMap::Dart::endNode, crp)
            .def("leftFaceLabel", &GeoMap::Dart::leftFaceLabel)
            .def("leftFace", &GeoMap::Dart::leftFace, crp)
            .def("rightFaceLabel", &GeoMap::Dart::rightFaceLabel)
            .def("rightFace", &GeoMap::Dart::rightFace, crp)
            .def("partialArea", &GeoMap::Dart::partialArea, crp)
            .def("__getitem__", &Array__getitem__<GeoMap::Dart>)
            .def("__iter__", &GeoMap::Dart::pointIter)
            .def("__len__", &GeoMap::Dart::size)
            .def("nextAlpha", &GeoMap::Dart::nextAlpha, rself)
            .def("nextSigma", &GeoMap::Dart::nextSigma, (arg("times") = 1), rself)
            .def("prevSigma", &GeoMap::Dart::prevSigma, (arg("times") = 1), rself)
            .def("nextPhi", &GeoMap::Dart::nextPhi, rself)
            .def("prevPhi", &GeoMap::Dart::prevPhi, rself)
            .def(self == self)
            .def(self != self)
        ;

#ifndef USE_INSECURE_CELL_PTRS
        register_ptr_to_python< CELL_PTR(GeoMap::Node) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Edge) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Face) >();
#endif
    }

    RangeIterWrapper<ContourPointIter>("ContourPointIter")
        .def(init<GeoMap::Dart, bool>((arg("dart"), arg("firstTwice") = false)));

    def("contourArea", &contourArea,
        "contourArea(anchor) -> float\n\n"
        "Returns the area of contourPoly(anchor) (is however much faster than\n"
        "using that function, since it simply sums up all partialArea()s of the\n"
        "darts in the phi orbit.");
    def("contourPoly", &contourPoly,
        "contourPoly(anchor) -> Polygon\n\n"
        "Returns a Polygon composed by traversing anchor's phi orbit once.");
}
