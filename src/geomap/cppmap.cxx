#include "cppmap.hxx"
#include <vigra/tinyvector.hxx>
#include <vigra/pythonimage.hxx>
#include <vigra/pythonutil.hxx>
#include <vigra/copyimage.hxx>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "exporthelpers.hxx"

#ifdef _MSC_VER
template <class T>
inline bool isnan(T t) { return _isnan(t); }
#endif

template<class Container>
void removeAll(Container &container,
               const typename Container::value_type &element)
{
    typename Container::iterator it;
    while((it = std::find(container.begin(), container.end(), element))
          != container.end())
        container.erase(it);
}

template<class Container>
void removeOne(Container &container,
               const typename Container::value_type &element)
{
    container.erase(std::find(container.begin(), container.end(), element));
}

/********************************************************************/

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
    friend class GeoMap; // give access to darts_ (add edge, sort edges, Euler..)

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
        --map_->nodeCount_;
        map_->nodeMap_.erase(
            map_->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                                   vigra::NumericTraits<double>::epsilon()));
        RESET_PTR(map_->nodes_[label_]); // may have effect like "delete this;"
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
    GeoMap      *map_;
    CellLabel    label_;
    CellLabel    startNodeLabel_, endNodeLabel_;
    CellLabel    leftFaceLabel_, rightFaceLabel_;
    unsigned int protection_;

    mutable std::auto_ptr<vigra::Scanlines> scanLines_;

    friend class Dart; // allow setLeftFaceLabel
    friend CELL_PTR(GeoMap::Edge) GeoMap::mergeEdges(Dart &);
    friend CELL_PTR(GeoMap::Edge) GeoMap::splitEdge(
        Edge &, unsigned int, const vigra::Vector2 &, bool);
    friend void GeoMap::splitParallelEdges();

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
      protection_(0),
      scanLines_(NULL)
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
        --map_->edgeCount_;
        RESET_PTR(map_->edges_[label_]); // may have effect like "delete this;"
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

    unsigned int protection() const
    {
        return protection_;
    }

    void protect(unsigned int flag, bool onoff)
    {
        if(onoff)
            protection_ |= flag;
        else
            protection_ &= ~flag;
    }

    GeoMap *map() const
    {
        return map_;
    }

    const vigra::Scanlines &scanLines() const
    {
        if(!scanLines_.get())
            scanLines_ = scanPoly(*this);
        return *scanLines_;
    }

  private:
    Edge(const Edge &) : vigra::BBoxPolygon<vigra::Vector2>() {} // disallow copying
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
         * Change the direction of traversal, without changing the
         * current position.  (Thus, atEnd()/inRange() will also not
         * change.)
         */
    void reverse()
    {
        if(inc_ < 0)
        {
            inc_ = 1;
            end_ = edge_->size();
        }
        else
        {
            inc_ = -1;
            end_ = -1;
        }
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
    friend CELL_PTR(GeoMap::Face) GeoMap::mergeFaces(Dart &);

    Dart &turnSigma(int times)
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

    GeoMap *map() const
    {
        return map_;
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
        if(!result)
        {
            std::stringstream s;
            s << "Cannot operate on invalid dart " << label()
              << " belonging to removed edge!";
            vigra_precondition(result, s.str());
        }
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

    Dart &nextSigma()
    {
        return turnSigma(1);
    }

    Dart &prevSigma()
    {
        return turnSigma(-1);
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
};

inline GeoMap::Dart GeoMap::Edge::dart() const
{
    return map_->dart(label());
}

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
    typedef std::vector<Dart> Contours;
    typedef Contours::const_iterator ContourIterator;

  protected:
    GeoMap             *map_;
    CellLabel           label_;
    std::vector<Dart>   anchors_;
    mutable BoundingBox boundingBox_;
    mutable bool        boundingBoxValid_;
    mutable double      area_;
    mutable bool        areaValid_;
    unsigned int        pixelArea_;

    friend class GeoMap; // give access to pixelArea_ and anchors_ (Euler ops...)

    unsigned int findComponentAnchor(const GeoMap::Dart &dart);

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
        --map_->faceCount_;
        RESET_PTR(map_->faces_[label_]); // may have effect like "delete this;"
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
        if(map_->labelImage_)
        {
            int l = (*map_->labelImage_)[intVPos(point)];
            if(l > 0 && (map_->faceLabelLUT_[l] == label_))
                return true;
        }
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

    unsigned int pixelArea() const
    {
        return pixelArea_;
    }

    const Dart &contour(unsigned int index = 0)
    {
        return anchors_[index];
    }

    ContourIterator contoursBegin() const
    {
        return anchors_.begin();
    }

    ContourIterator contoursEnd() const
    {
        return anchors_.end();
    }

    ContourIterator holesBegin() const
    {
        ContourIterator result(anchors_.begin());
        if(label())
            ++result;
        return result;
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
            GeoMap::Edge &edge(*map_->edge(-(int)i));
            edge[edge.size()-1] = p;
        }
    }
    map_->nodeMap_.insert(PositionedNodeLabel(p, label_));
}

inline GeoMap::Dart GeoMap::Node::anchor() const
{
    vigra_precondition(initialized(), "anchor() of uninitialized node!");
    vigra_precondition(degree(), "anchor() of degree 0 node!");
    return Dart(map_, darts_[0]);
}

GeoMap::GeoMap(vigra::Size2D imageSize)
: nodeCount_(0),
  edgeCount_(0),
  faceCount_(0),
  imageSize_(imageSize),
  labelImage_(NULL),
  edgesSorted_(false)
{
    edges_.push_back(NULL_PTR(Edge));
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
                return face(faceLabelLUT_[faceLabel]);
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

CELL_PTR(GeoMap::Node) GeoMap::addNode(
    const vigra::Vector2 &position, CellLabel label)
{
    if(label > nodes_.size())
        nodes_.resize(label, NULL_PTR(GeoMap::Node));
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

CELL_PTR(GeoMap::Edge) GeoMap::addEdge(
    GeoMap::Dart startNeighbor, GeoMap::Dart endNeighbor,
    const Vector2Array &points)
{
    CELL_PTR(GeoMap::Node)
        startNode = startNeighbor.startNode(),
        endNode = endNeighbor.startNode();
    GeoMap::Edge *result = new GeoMap::Edge(
        this, startNode->label(), endNode->label(), points);

    GeoMap::Node::DartLabels::iterator
        pos1(std::find(startNode->darts_.begin(),
                       startNode->darts_.end(), startNeighbor.label())),
        pos2(std::find(endNode  ->darts_.begin(),
                       endNode  ->darts_.end(), endNeighbor  .label()));
    startNode->darts_.insert(pos1,  result->label());
    endNode  ->darts_.insert(pos2, -(int)result->label());
    return edge(result->label());
}

void GeoMap::removeEdge(GeoMap::Dart &dart)
{
    if(mapInitialized())
    {
        if(dart.edge()->isBridge())
            removeBridge(dart);
        else
            mergeFaces(dart);
    }
    else
    {
        // this is just a graph -> detach from nodes & uninitialize()

        GeoMap::Node &node1(*dart.startNode());
        GeoMap::Node &node2(*dart.endNode());

        removeOne(node1.darts_,  dart.label());
        removeOne(node2.darts_, -dart.label());

        dart.edge()->uninitialize();
    }
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

    edgesSorted_ = true;
}

class DartPosition
{
  public:
    DartPosition(const GeoMap::Dart &dart)
    : dart_(dart),
      pointIter_(dart),
      segmentIndex_(0),
      arcLength_(0.0),
      partialArcLength_(0.0),
      position_(*pointIter_)
    {
        p1_ = *pointIter_;
        p2_ = *++pointIter_;
    }

    bool atEnd() const
    {
        return pointIter_.atEnd();
    }

    const vigra::Vector2 &operator()() const
    {
        return position_;
    }

    GeoMap::Dart dart() const
    {
        return dart_;
    }

    int dartLabel() const
    {
        return dart_.label();
    }

    unsigned int segmentIndex() const
    {
        return segmentIndex_;
    }

    double arcLength() const
    {
        return arcLength_ + partialArcLength_;
    }

    const vigra::Vector2 &segmentStart() const
    {
        return p1_;
    }

    const vigra::Vector2 &segmentEnd() const
    {
        return p2_;
    }

    double segmentLength() const
    {
        return (p2_ - p1_).magnitude();
    }

    bool gotoArcLength(double arcLength)
    {
        while(arcLength < arcLength_)
            if(!prevSegmentInternal())
                return false;
        do
        {
            double rest = arcLength - arcLength_;
            Vector2 diff(p2_ - p1_);
            if(diff.squaredMagnitude() > rest*rest)
            {
                position_ = p1_ + diff*rest/diff.magnitude();
                partialArcLength_ = rest;
                return true;
            }
        }
        while(nextSegmentInternal());
        partialArcLength_ = 0.0;
        return false;
    }

    bool gotoNextSegment()
    {
        bool result = nextSegmentInternal();
        position_ = p1_;
        partialArcLength_ = 0.0;
        return result;
    }

    bool gotoPrevSegment()
    {
        bool result = prevSegmentInternal();
        position_ = p1_;
        partialArcLength_ = 0.0;
        return result;
    }

    bool leaveCircle(const vigra::Vector2 &center, double radius2)
    {
        while((p2_ - center).squaredMagnitude() < radius2)
            if(!nextSegmentInternal())
                break;

        position_ = p2_;
        partialArcLength_ = (p2_ - p1_).magnitude();
        return !atEnd();
    }

    bool intersectCircle(const vigra::Vector2 &center, double radius2)
    {
        // unfortunately, this prevents larger steps:
//         if((p1_ - center).squaredMagnitude() >= radius2)
//         {
//             std::cerr << "intersectCircle: we are already outside!\n";
//             position_ = p1_;
//             return;
//         }
        while((p2_ - center).squaredMagnitude() < radius2)
        {
            if(!nextSegmentInternal())
            {
                position_ = p2_;
                return false;
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
        if(!isnan(lambda))
            diff *= lambda;
        else
        {
            std::cerr << "intersectCircle: error interpolating between " << p1_ << " and " << p2_ << " to a squared distance of " << radius2 << " from " << center << "!\n";
        }
        position_ = p1_ + diff;
        partialArcLength_ = diff.magnitude();
        return true;
    }

  protected:
    bool nextSegmentInternal()
    {
        if(atEnd())
            return false;
        arcLength_ += (p2_ - p1_).magnitude();
        p1_ = p2_;
        ++pointIter_;
        if(pointIter_.atEnd())
            return false;
        p2_ = *pointIter_;
        ++segmentIndex_;
        return true;
    }

    bool prevSegmentInternal()
    {
        if(!segmentIndex_)
            return false;
        // now I assume that pointIter_ can step backwards and still
        // is inRange():
        pointIter_.reverse();
        p2_ = p1_;
        ++pointIter_;
        pointIter_.reverse();
        p1_ = *pointIter_;
        arcLength_ -= (p2_ - p1_).magnitude();
        --segmentIndex_;
        return true;
    }

    GeoMap::Dart dart_;
    DartPointIter pointIter_;
    unsigned int segmentIndex_;
    double arcLength_, partialArcLength_;
    vigra::Vector2 p1_, p2_, position_;
};

struct DartPositionAngle
{
    struct EdgePosition
    {
        unsigned int segmentIndex;
        double arcLength;
        vigra::Vector2 position;
    };

    struct CommonPos : public EdgePosition
    {
        CommonPos()
        : isSet(false)
        {}

        const vigra::Vector2 &set(const DartPosition &dp)
        {
            position = dp();
            segmentIndex = dp.segmentIndex();
            arcLength = dp.arcLength();
            isSet = true;
            return position;
        }

        bool isSet;
    };

    DartPosition dp;
    double absAngle, angle;
    CommonPos commonPos;

    DartPositionAngle(const GeoMap::Dart &dart)
    : dp(dart)
    {}

    bool operator<(const DartPositionAngle &other) const
    {
        return angle < other.angle;
    }

    struct SplitPos : public EdgePosition
    {
        int dartLabel, sigmaPos;
        unsigned int splitGroup, newEdgeLabel;

        SplitPos(const EdgePosition &ep, int dl, unsigned int sg)
        : EdgePosition(ep),
          dartLabel(dl),
          splitGroup(sg)
        {}

        bool operator<(const SplitPos &other) const
        {
            return arcLength > other.arcLength;
        }
    };

    SplitPos splitPos(unsigned int group) const
    {
        vigra_precondition(
            commonPos.isSet, "splitPos() called with uninitialized commonPos");
        SplitPos result(commonPos, dp.dartLabel(), group);
        if(dp.dartLabel() < 0)
        {
            GeoMap::Edge &edge(*dp.dart().edge());
            result.segmentIndex = edge.size()-2 - result.segmentIndex;
            result.arcLength = edge.length() - result.arcLength;
        }
        return result;
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

class PlannedSplits : public std::vector<DartPositionAngle::SplitPos>
{
  public:
    PlannedSplits()
    : splitGroupCount(0)
    {}

    unsigned int splitGroupCount;
};

void sortEdgesInternal(const vigra::Vector2 &currentPos,
                       double referenceAngle,
                       DPAI dpBegin, DPAI dpEnd,
                       double stepDist2, double minAngle,
                       GeoMap::UnsortableGroups &unsortable,
                       PlannedSplits *splitInfo,
                       bool parallel)
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
        GeoMap::UnsortableGroups::value_type unsortableGroup;
        for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
            unsortableGroup.push_back(dpi->dp.dartLabel());
        unsortable.push_back(unsortableGroup);
        return;
    }

    std::sort(dpBegin, dpEnd);

    // handle cyclicity of array first (by rotation if necessary):
    DPAI firstGroupStart = dpEnd;
    if((--firstGroupStart)->angle + minAngle - 2*M_PI > dpBegin->angle)
    {
        // first and last dart are less than minAngle apart,..
        DPAI prev(firstGroupStart);
        while(--firstGroupStart > dpBegin)
        {
            // ..determine last decision point..
            if(prev->angle - firstGroupStart->angle >= minAngle)
            {
                // ..and rotate array to have whole dart group together
                // (start with "prev" == last correct firstGroupStart as newBegin)
                rotateArray(dpBegin, prev, dpEnd);
                break;
            }
            prev = firstGroupStart;
        }
    }

    bool storedSplitPos = false;
    PlannedSplits::size_type storedSplitsOffset = 0;

    // look for groups of parallel edges
    DPAI groupStart = dpBegin,
          groupLast = groupStart, // for convenience; this is always groupEnd - 1
           groupEnd = groupLast + 1;
    for(; true; ++groupLast, ++groupEnd)
    {
        // group ending?
        if((groupEnd == dpEnd) || // last group
           (groupEnd->angle >= groupLast->angle + minAngle)) // decision here
        {
            if(splitInfo && groupEnd != dpEnd && parallel && !storedSplitPos)
            {
                // build group of edges to be split:
                // (since we have been given group of parallel darts and could
                // decide about the order at the current position, but not in
                // the last common position)
                unsigned int splitGroup = splitInfo->splitGroupCount++;
                storedSplitsOffset = splitInfo->size();
                for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
                    splitInfo->push_back(dpi->splitPos(splitGroup));
                storedSplitPos = true;
            }

            // recursion needed if > one dart in group:
            if(groupLast != groupStart)
            {
                // determine mean position of dart positions in subgroup:
                vigra::Vector2 meanPos(0, 0);
                for(DPAI dpi = groupStart; dpi != groupEnd; ++dpi)
                    meanPos += dpi->commonPos.set(dpi->dp);
                meanPos /= (groupEnd - groupStart);

                // sort parallel subgroup recursively:
                sortEdgesInternal(meanPos, normAngle(
                                      groupStart->absAngle +
                                      normAngle(groupLast->absAngle -
                                                groupStart->absAngle) / 2),
                                  groupStart, groupEnd,
                                  stepDist2, minAngle,
                                  unsortable, splitInfo,
                                  true);
            }

            if(groupEnd == dpEnd)
                break; // loop end

            groupStart = groupEnd;
        }
    }

    if(storedSplitPos)
    {
        // Later, in order to get the merging of the split nodes
        // right, we need to know the sigma order we just found out
        // here.  We do not do it then, because it becomes much more
        // complicated after the splitting.

        PlannedSplits::iterator
            storedSplitsBegin(splitInfo->begin() + storedSplitsOffset),
            storedSplitsEnd(storedSplitsBegin + (dpEnd - dpBegin));

        // Unfortunately, this is O(n^2) - I thought long about this,
        // but I cannot see any shortcut.  Maybe I am too tired, but
        // right now I have to get it working.
        int sigmaPos = 0;
        for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi, ++sigmaPos)
        {
            int dartLabel = dpi->dp.dartLabel();
            for(PlannedSplits::iterator splIt = // wordplay ;-)
                    storedSplitsBegin; splIt != storedSplitsEnd; ++splIt)
            {
                if(splIt->dartLabel == dartLabel)
                {
                    splIt->sigmaPos = sigmaPos;
                    break;
                }
            }
        }
    }
}

void GeoMap::sortEdgesEventually(double stepDist, double minDist,
                                 UnsortableGroups &unsortable,
                                 bool splitEdges)
{
    double minAngle = std::atan2(minDist, stepDist),
          stepDist2 = vigra::sq(stepDist);

    if(splitEdges)
        splitInfo_ = std::auto_ptr<PlannedSplits>(new PlannedSplits());

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        GeoMap::Node::DartLabels &dartLabels((*it)->darts_);

        DartPositionAngles dartPositions;

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
            dartPositions.push_back(
                DartPositionAngle(dart(dartLabels[i])));

        //std::cerr << "sorting darts of node " << (*it)->label() << ":\n";
        sortEdgesInternal((*it)->position(), 0.0,
                          dartPositions.begin(), dartPositions.end(),
                          stepDist2, minAngle,
                          unsortable, splitInfo_.get(), false);

        for(unsigned int i = 0; i < dartLabels.size(); ++i)
            dartLabels[i] = dartPositions[i].dp.dartLabel();
    }

    edgesSorted_ = true;
}

struct MergeDart
{
    int dartLabel, sigmaPos;
    bool turnLater;

    MergeDart(int dl, bool tl, int sp)
    : dartLabel(dl), sigmaPos(sp), turnLater(tl)
    {}

        // default constructor - mark as uninitialized through zero dartLabel
    MergeDart()
    : dartLabel(0)
    {}

    bool operator<(MergeDart const &other) const
    {
        return !dartLabel || (other.dartLabel &&
                              sigmaPos < other.sigmaPos);
    }
};

std::string Dart__repr__(GeoMap::Dart const &dart);

void GeoMap::splitParallelEdges()
{
    vigra_precondition(splitInfo_.get(), "splitParallelEdges(): no planned splits (set splitEdges parameter of sortEdgesEventually?)");

    std::vector<PlannedSplits::difference_type> groupPositions;
    for(PlannedSplits::iterator it = splitInfo_->begin();
        it != splitInfo_->end(); ++it)
    {
        if(it->splitGroup == groupPositions.size())
            groupPositions.push_back(it - splitInfo_->begin());
    }

    std::sort(splitInfo_->begin(), splitInfo_->end());

    typedef std::vector<MergeDart> MergeDarts;
    MergeDarts mergeDarts(splitInfo_->size());
    for(PlannedSplits::iterator it = splitInfo_->begin();
        it != splitInfo_->end(); ++it)
    {
        CELL_PTR(Edge) newEdge(splitEdge(
                                   *edge(abs(it->dartLabel)),
                                   it->segmentIndex, it->position));

        PlannedSplits::difference_type &pos(
            groupPositions[it->splitGroup]);

        mergeDarts[pos] =
            MergeDart((int)newEdge->label(), it->dartLabel > 0, it->sigmaPos);

        ++pos;
    }

    splitInfo_.reset(); // splitting finished, free memory

    // for each split group, merge the resulting nodes:
    const double checkSurvivorDist  = 0.5;
    const double checkSurvivorDist2 = checkSurvivorDist*checkSurvivorDist;

    MergeDarts::iterator mergeDartsGroupEnd = mergeDarts.begin();
    for(unsigned int i = 0; i < groupPositions.size(); ++i)
    {
        MergeDarts::iterator mergeDartsGroupBegin = mergeDartsGroupEnd;
        mergeDartsGroupEnd = mergeDarts.begin() + (int)groupPositions[i];

        std::sort(mergeDartsGroupBegin, mergeDartsGroupEnd);

        // skip deleted edges (e.g. unsortable)
        while(mergeDartsGroupBegin->dartLabel == 0)
            ++mergeDartsGroupBegin;

        // now search for the best continuation to choose the survivor
        double bestContinuationValue = 0.0;
        int bestContinuationIndex = 0;

        std::vector<DartPosition> dps;
        for(MergeDarts::iterator it = mergeDartsGroupBegin;
            it != mergeDartsGroupEnd; ++it)
        {
            GeoMap::Dart d(dart(it->dartLabel));
            vigra::Vector2 nodePos(d.startNode()->position());

            DartPosition dp1(d);
            dp1.leaveCircle(nodePos, checkSurvivorDist2);
            d.nextSigma();
            DartPosition dp2(d);
            dp2.leaveCircle(nodePos, checkSurvivorDist2);

            vigra::Vector2
                v1(dp1() - nodePos),
                v2(nodePos - dp2());

            double cont = dot(v1, v2)/(v1.magnitude()*v2.magnitude());
            if(cont > bestContinuationValue)
            {
                bestContinuationValue = cont;
                bestContinuationIndex = it - mergeDartsGroupBegin;
            }
        }

        Dart survivor(
            dart(mergeDartsGroupBegin[bestContinuationIndex].dartLabel));
        if(mergeDartsGroupBegin[bestContinuationIndex].turnLater)
            survivor.nextSigma();

        GeoMap::Node &survivingNode(*survivor.startNode());
        int sigmaNeighborLabel(survivor.clone().nextSigma().label());
        GeoMap::Node::DartLabels::iterator
            targetSigmaPos(
                std::find(survivingNode.darts_.begin(), survivingNode.darts_.end(),
                          sigmaNeighborLabel));

        if(bestContinuationIndex)
        {
            // loop from survivor-1 to "begin" of sigma orbit
            MergeDarts::iterator it = mergeDartsGroupBegin + bestContinuationIndex;
            do
            {
                --it;

                // determine mergeDart and relocateDart
                GeoMap::Dart mergeDart(dart(it->dartLabel));
                if(it->turnLater)
                    mergeDart.nextSigma();

                vigra_invariant(mergeDart.startNode()->degree() == 2,
                                "merge nodes are expected to have degree 2");

                GeoMap::Dart relocateDart(mergeDart);
                relocateDart.nextSigma();

                // re-attach relocateDart to surviving node
                CELL_PTR(GeoMap::Edge) relocateEdge(relocateDart.edge());
                if(relocateDart.label() < 0)
                {
                    relocateEdge->endNodeLabel_ = survivingNode.label();
                    (*relocateEdge)[relocateEdge->size() - 1] =
                        survivingNode.position();
                }
                else
                {
                    relocateEdge->startNodeLabel_ = survivingNode.label();
                    (*relocateEdge)[0] = survivingNode.position();
                }
                targetSigmaPos =
                    survivingNode.darts_.insert(targetSigmaPos, relocateDart.label());

                // remove mergeDart and its startNode:
                GeoMap::Node &mergedNode(*mergeDart.startNode());
                removeEdge(mergeDart);
                removeIsolatedNode(mergedNode);
            }
            while(it != mergeDartsGroupBegin);
        }

        targetSigmaPos = std::find(
            survivingNode.darts_.begin(), survivingNode.darts_.end(),
            sigmaNeighborLabel);
        ++targetSigmaPos;

        for(MergeDarts::iterator it =
                mergeDartsGroupBegin + bestContinuationIndex + 1;
            it != mergeDartsGroupEnd; ++it)
        {
            // determine mergeDart and relocateDart
            GeoMap::Dart mergeDart(dart(it->dartLabel));
            if(it->turnLater)
                mergeDart.nextSigma();

            vigra_invariant(mergeDart.startNode()->degree() == 2,
                            "merge nodes are expected to have degree 2");

            GeoMap::Dart relocateDart(mergeDart);
            relocateDart.nextSigma();

            // re-attach relocateDart to surviving node
            CELL_PTR(GeoMap::Edge) relocateEdge(relocateDart.edge());
            if(relocateDart.label() < 0)
            {
                relocateEdge->endNodeLabel_ = survivingNode.label();
                (*relocateEdge)[relocateEdge->size() - 1] =
                    survivingNode.position();
            }
            else
            {
                relocateEdge->startNodeLabel_ = survivingNode.label();
                (*relocateEdge)[0] = survivingNode.position();
            }
            targetSigmaPos =
                survivingNode.darts_.insert(targetSigmaPos, relocateDart.label());
            ++targetSigmaPos;

            // remove mergeDart and its startNode:
            GeoMap::Node &mergedNode(*mergeDart.startNode());
            removeEdge(mergeDart);
            removeIsolatedNode(mergedNode);
        }
    }
}

void GeoMap::setSigmaMapping(SigmaMapping const &sigmaMapping, bool sorted)
{
    vigra_precondition(sigmaMapping.size() >= (2*edges_.size() - 1),
                       "setSigmaMapping: sigmaMapping too small!");
    SigmaMapping::const_iterator sigma(
        sigmaMapping.begin() + (sigmaMapping.size() / 2));

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if(!(*it)->degree())
            continue;

        GeoMap::Node::DartLabels &dartLabels((*it)->darts_);

        GeoMap::Node::DartLabels singleOrbit;
        int label = (*it)->anchor().label();
        do
        {
            singleOrbit.push_back(label);
            label = sigma[label];
        }
        while(label != singleOrbit[0]);

        vigra_precondition(singleOrbit.size() == dartLabels.size(),
                           "setSigmaMapping: sigma orbit has wrong size");
        std::swap(singleOrbit, dartLabels);
    }

    edgesSorted_ = sorted;
}

std::auto_ptr<GeoMap::SigmaMapping> GeoMap::sigmaMapping()
{
    std::auto_ptr<GeoMap::SigmaMapping> result(
        new GeoMap::SigmaMapping(2*edges_.size() - 1, 0));
    GeoMap::SigmaMapping::iterator
        sigma(result->begin() + (result->size() / 2));

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if(!(*it)->degree())
            continue;

        GeoMap::Dart anchor((*it)->anchor()), d(anchor);
        do
        {
            int dl = d.label();
            sigma[dl] = d.nextSigma().label();
        }
        while(d != anchor);
    }

    return result;
}

void GeoMap::initializeMap(bool initLabelImage)
{
    vigra_precondition(!mapInitialized(),
                       "initializeMap() called more than once");
    if(!edgesSorted())
        sortEdgesDirectly();

    initContours();
    embedFaces(initLabelImage);
}

void GeoMap::initContours()
{
    new Face(this, Dart(this, 0)); // create infinite face, dart will be ignored

    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        if((*it)->leftFaceLabel() == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart( (int)(*it)->label()));
        if((*it)->rightFaceLabel() == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart(-(int)(*it)->label()));
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

typedef vigra::MultiArray<2, int> LabelImage;

void markEdgeInLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage);

void GeoMap::embedFaces(bool initLabelImage)
{
    vigra_precondition(!labelImage_,
        "embedFaces() called with already-initialized labelImage");

    if(initLabelImage)
    {
        vigra_precondition(imageSize_.area() > 0,
                           "initLabelImage: non-zero imageSize must be given!");
        labelImage_ = new LabelImage(
            LabelImage::size_type(imageSize().width(), imageSize().height()), 0);
        faceLabelLUT_.resize(faces_.size());
    }

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
                // no need for rawAddEdgeToLabelImage here, since we
                // work with darts anyways, and there's no easy way to
                // ensure that the negative counts will not be wrong
                // (esp. also for interior bridges etc.)
                drawScannedPoly(*scanlines, -1,
                                labelImage_->traverser_begin(),
                                labelImage_->size(),
                                vigra::StandardValueAccessor<int>());
                faceLabelLUT_[contour.label()] = contour.label();
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
//             std::cerr << "  embedding contour " << contour.label()
//                       << " in face " << parent->label() << "\n";
            parent->embedContour(anchor);
            contour.uninitialize();
        }
    }

    if(initLabelImage)
    {
        // in the case of holes, the face's pixelArea()s will be wrong:
        for(FaceIterator it = facesBegin(); it.inRange(); ++it)
            (*it)->pixelArea_ = 0;

        // interior bridges may not have been set to negative labels yet:
        for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
            if((*it)->isBridge())
                drawScannedPoly((*it)->scanLines(), -1,
                                labelImage_->traverser_begin(),
                                labelImage_->size(),
                                vigra::StandardValueAccessor<int>());

        // remove temporary edge markings and fix pixelAreas:
        for(GeoMap::LabelImage::traverser lrow = labelImage_->traverser_begin();
            lrow != labelImage_->traverser_end(); ++lrow)
        {
            for(GeoMap::LabelImage::traverser::next_type lit = lrow.begin();
                lit != lrow.end(); ++lit)
            {
                int label = *lit;
                if(label < 0)
                    *lit = 0;
                else
                    ++face(label)->pixelArea_;
            }
        }

        // redo all edge markings correctly:
        for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
            markEdgeInLabelImage((*it)->scanLines(),
                                 *labelImage_);
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
    //std::cerr << "GeoMap[" << this << "].checkConsistency()\n";
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
        GeoMap::Face &face(**it);
        if(face.map() != this)
        {
            std::cerr << "  Face " << face.label() << " has wrong map()!\n";
            result = false;
            break;
        }

        for(GeoMap::Face::ContourIterator ci = face.contoursBegin();
            ci != face.contoursEnd(); ++ci)
        {
            Dart anchor(*ci), dart(anchor);
            do
            {
                if(dart.leftFaceLabel() != face.label())
                {
                    std::cerr << "  Dart " << dart.label()
                              << " has leftFaceLabel() " << dart.leftFaceLabel()
                              << " != " << face.label() << "!\n";
                    result = false;
                    break;
                }
            }
            while(dart.nextPhi() != anchor);
        }
    }
    return result;
}

/********************************************************************/

bool GeoMap::removeIsolatedNode(GeoMap::Node &node)
{
    if(!removeNodeHook(node))
        return false;

    node.uninitialize();
    return true;
}

void rawAddEdgeToLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage, int diff)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(labelImage.size()[1], scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size()[0])
                end = labelImage.size()[0];

            for(int x = begin; x < end; ++x)
                labelImage[LabelImage::difference_type(x, y)] += diff;
        }
    }
}

void markEdgeInLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(labelImage.size()[1], scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size()[0])
                end = labelImage.size()[0];

            for(int x = begin; x < end; ++x)
            {
                LabelImage::difference_type pos(x, y);
                int label = labelImage[pos];
                labelImage[pos] = (label >= 0 ? -1 : label-1);
            }
        }
    }
}

void removeEdgeFromLabelImage(
    const vigra::Scanlines &scanlines,
    LabelImage &labelImage,
    LabelImage::value_type substituteLabel,
    PixelList &outputPixels)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(labelImage.size()[1], scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size()[0])
                end = labelImage.size()[0];

            for(int x = begin; x < end; ++x)
            {
                LabelImage::reference old(labelImage(x, y));
                if(old != -1)
                {
                    old += 1;
                }
                else
                {
                    old = substituteLabel;
                    outputPixels.push_back(vigra::Point2D(x, y));
                }
            }
        }
    }
}

CELL_PTR(GeoMap::Edge) GeoMap::mergeEdges(GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "mergeEdges called on removed dart!");
    Dart d1(dart);
    d1.nextSigma();
    vigra_precondition(d1.edgeLabel() != dart.edgeLabel(),
                       "mergeEdges called on self-loop!");

    Dart d2(d1);
    d2.nextSigma();
    vigra_precondition(d2 == dart,
                       "mergeEdges cannot remove node with degree > 2!");

    vigra_assert(d1.leftFaceLabel() == d2.rightFaceLabel(),
                 "mergeEdges: broken map");
    vigra_assert(d2.leftFaceLabel() == d1.rightFaceLabel(),
                 "mergeEdges: broken map");

    GeoMap::Face *faces[2];
    faces[0] = &(*dart.leftFace());
    faces[1] = &(*dart.rightFace());
    for(GeoMap::Face **faceIt = faces; faceIt != faces+2; ++faceIt)
    {
        GeoMap::Face::Contours::iterator cEnd = (*faceIt)->anchors_.end();
        for(GeoMap::Face::Contours::iterator it = (*faceIt)->anchors_.begin();
            it != cEnd; ++it)
        {
            if(it->edgeLabel() == d2.edgeLabel())
            {
                it->nextPhi();
                break;
            }
        }
    }

    GeoMap::Node &mergedNode(*d1.startNode());
    GeoMap::Edge &survivor(*d1.edge());
    GeoMap::Edge &mergedEdge(*d2.edge());

    if(!preMergeEdgesHook(d1))
        return NULL_PTR(GeoMap::Edge);
    if(!removeNodeHook(mergedNode))
        return NULL_PTR(GeoMap::Edge);

    // TODO: history append?

    GeoMap::Dart d(d2);
    d.nextAlpha();
    GeoMap::Node &changedEndNode(*d.startNode());
    unsigned int cenDartIndex = 0;
    for(; cenDartIndex < changedEndNode.darts_.size(); ++cenDartIndex)
        if(changedEndNode.darts_[cenDartIndex] == d.label())
            break;
    vigra_postcondition(cenDartIndex < changedEndNode.darts_.size(),
                        "changing dart not found at node");

    if(labelImage_)
    {
        rawAddEdgeToLabelImage(mergedEdge.scanLines(), *labelImage_, 1);
        rawAddEdgeToLabelImage(survivor.scanLines(), *labelImage_, 1);
    }

    if(survivor.startNodeLabel() != mergedNode.label())
    {
        if(mergedEdge.startNodeLabel() != mergedNode.label())
            mergedEdge.reverse();
        survivor.extend(mergedEdge);
        survivor.endNodeLabel_ = changedEndNode.label();
    }
    else
    {
        survivor.reverse();
        if(mergedEdge.startNodeLabel() != mergedNode.label())
            mergedEdge.reverse();
        survivor.extend(mergedEdge);
        survivor.reverse();
        survivor.startNodeLabel_ = changedEndNode.label();
    }

    changedEndNode.darts_[cenDartIndex] = d1.label();

    if(labelImage_)
    {
        (*survivor.scanLines_) += mergedEdge.scanLines();
        rawAddEdgeToLabelImage(survivor.scanLines(), *labelImage_, -1);
    }

    mergedNode.uninitialize();
    mergedEdge.uninitialize();

    postMergeEdgesHook(survivor);

    return this->edge(survivor.label());
}

CELL_PTR(GeoMap::Edge) GeoMap::splitEdge(
    GeoMap::Edge &edge, unsigned int segmentIndex)
{
    return splitEdge(edge, segmentIndex, edge[0], false);
}

CELL_PTR(GeoMap::Edge) GeoMap::splitEdge(
    GeoMap::Edge &edge, unsigned int segmentIndex,
    const vigra::Vector2 &newPoint, bool insertPoint)
{
    vigra_precondition(segmentIndex < edge.size() - 1,
                       "splitEdge: invalid segmentIndex");

    preSplitEdgeHook(edge, segmentIndex, newPoint, insertPoint);

    CELL_PTR(GeoMap::Node) newNode(
        addNode(insertPoint ? newPoint : edge[segmentIndex]));

    if(labelImage_)
        rawAddEdgeToLabelImage(edge.scanLines(), *labelImage_, 1);

    if(insertPoint)
    {
        edge.insert(edge.begin() + segmentIndex + 1, newPoint);
        ++segmentIndex;
    }

    GeoMap::Edge *result = new GeoMap::Edge(
        this, newNode->label(), edge.endNodeLabel(), edge.split(segmentIndex));
    result->leftFaceLabel_ = edge.leftFaceLabel_;
    result->rightFaceLabel_ = edge.rightFaceLabel_;
    result->protection_ = edge.protection_;

    newNode->darts_.push_back(-(int)edge.label());
    newNode->darts_.push_back( (int)result->label());

    edge.endNodeLabel_ = newNode->label();

    GeoMap::Node &changedNode(*result->endNode());
    *std::find(changedNode.darts_.begin(),
               changedNode.darts_.end(), -(int)edge.label())
        = -(int)result->label();

    if(labelImage_)
    {
        edge.scanLines_.reset();
        rawAddEdgeToLabelImage(edge.scanLines(), *labelImage_, -1);
        rawAddEdgeToLabelImage(result->scanLines(), *labelImage_, -1);
    }

    postSplitEdgeHook(edge, *result);

    return this->edge(result->label());
}

unsigned int GeoMap::Face::findComponentAnchor(const GeoMap::Dart &dart)
{
    for(unsigned int i = 0; i < anchors_.size(); ++i)
        if(anchors_[i] == dart)
            return i;

    for(unsigned int i = 0; i < anchors_.size(); ++i)
    {
        GeoMap::Dart d(anchors_[i]);
        while(d.nextPhi() != anchors_[i])
            if(d == dart)
                return i;
    }

    vigra_fail("findComponentAnchor failed: dart not found in face contours!");
    return 42; // never reached
}

void GeoMap::associatePixels(GeoMap::Face &face, const PixelList &pixels)
{
    face.pixelArea_ += pixels.size();
//     for(unsigned int i = 0; i < pixels.size(); ++i)
//         pixelBounds_ |= pixels[i];
    associatePixelsHook(face, pixels);
}

CELL_PTR(GeoMap::Face) GeoMap::removeBridge(GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "removeBridge called on removed dart!");

    GeoMap::Edge &edge(*dart.edge());
    GeoMap::Face &face(*dart.leftFace());
    vigra_precondition(face.label() == dart.rightFace()->label(),
                       "removeBridge needs a bridge dart!");
    GeoMap::Node &node1(*dart.startNode());
    GeoMap::Node &node2(*dart.endNode());
    vigra_precondition(node1.label() != node2.label(),
                       "Inconsistent map: bridge to be removed is also a self-loop!?");

    if(!preRemoveBridgeHook(dart))
        return NULL_PTR(GeoMap::Face);

    // TODO: history append?

    Dart newAnchor1(dart), newAnchor2(dart);
    newAnchor1.prevSigma();
    newAnchor2.nextAlpha().prevSigma();
    unsigned int contourIndex = face.findComponentAnchor(dart);

    removeOne(node1.darts_,  dart.label());
    removeOne(node2.darts_, -dart.label());

    if(contourIndex == 0)
    {
        // determine outer anchor, swap if necessary:
        if(newAnchor1.edgeLabel() == dart.edgeLabel() ||
           newAnchor2.edgeLabel() != dart.edgeLabel() &&
           contourArea(newAnchor1) < contourArea(newAnchor2))
            std::swap(newAnchor1, newAnchor2);
    }

    face.anchors_[contourIndex] = newAnchor1;
    face.anchors_.push_back(newAnchor2);

    PixelList associatedPixels;
    if(labelImage_)
        removeEdgeFromLabelImage(
            edge.scanLines(), *labelImage_, face.label(), associatedPixels);

    // remove singular nodes
    if(newAnchor1.edgeLabel() == dart.edgeLabel())
    {
        removeIsolatedNode(*newAnchor1.startNode());
        face.anchors_.erase(face.anchors_.begin() + contourIndex);
    }
    if(newAnchor2.edgeLabel() == dart.edgeLabel())
    {
        removeIsolatedNode(*newAnchor2.startNode());
        face.anchors_.erase(face.anchors_.end() - 1);
    }

    edge.uninitialize();

    postRemoveBridgeHook(face);

    if(associatedPixels.size())
        associatePixels(face, associatedPixels);

    return this->face(face.label());
}

CELL_PTR(GeoMap::Face) GeoMap::mergeFaces(GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "mergeFaces called on removed dart!");

    GeoMap::Dart removedDart(dart);

    if(dart.leftFace()->area() < dart.rightFace()->area())
        removedDart.nextAlpha();
    if(!removedDart.rightFaceLabel()) // face 0 shall stay face 0
        removedDart.nextAlpha();

    GeoMap::Edge &mergedEdge(*removedDart.edge());
    GeoMap::Face &survivor(*removedDart.leftFace());
    GeoMap::Face &mergedFace(*removedDart.rightFace());
    GeoMap::Node &node1(*removedDart.startNode());
    GeoMap::Node &node2(*removedDart.endNode());

    vigra_precondition(survivor.label() != mergedFace.label(),
                       "mergeFaces(): dart belongs to a bridge!");

    unsigned int contour1 = survivor.findComponentAnchor(removedDart);
    unsigned int contour2 = mergedFace.findComponentAnchor(
        GeoMap::Dart(removedDart).nextAlpha());

    if(!preMergeFacesHook(dart))
        return NULL_PTR(GeoMap::Face);

    // TODO: history append?

    // remember bounding box of merged face for later updating
    GeoMap::Face::BoundingBox mergedBBox;
    if(survivor.boundingBoxValid_)
        mergedBBox = mergedFace.boundingBox();

    // relabel contour's leftFaceLabel
    for(unsigned int i = 0; i < mergedFace.anchors_.size(); ++i)
    {
        GeoMap::Dart d(mergedFace.anchors_[i]);
        while(d.nextPhi().leftFaceLabel() != survivor.label())
            d.setLeftFaceLabel(survivor.label());
    }

    // re-use an old anchor for the merged contour
    if(survivor.anchors_[contour1].edgeLabel() == mergedEdge.label())
    {
        survivor.anchors_[contour1].nextPhi();
        if(survivor.anchors_[contour1].edgeLabel() == mergedEdge.label())
        {
            survivor.anchors_[contour1] = mergedFace.anchors_[contour2];
            if(survivor.anchors_[contour1].edgeLabel() == mergedEdge.label())
                survivor.anchors_[contour1].nextPhi();
        }
    }

    // check validity of found anchor
    if(survivor.anchors_[contour1].edgeLabel() == mergedEdge.label())
    {
        vigra_precondition(node1.label() == node2.label(),
                           "special-case: merging a self-loop");
        // results in an isolated node:
        survivor.anchors_.erase(survivor.anchors_.begin() + contour1);
    }

    // copy all remaining anchors into survivor's list:
    for(unsigned int i = 0; i < mergedFace.anchors_.size(); ++i)
    {
        if(i != contour2)
            survivor.anchors_.push_back(mergedFace.anchors_[i]);
    }

    // relabel region in image
    PixelList associatedPixels;
    if(labelImage_)
    {
//         relabelImage(map.labelImage.subImage(mergedFace.pixelBounds_),
//                      mergedFace.label(), survivor.label())
        for(unsigned int i = 0; i < faceLabelLUT_.size(); ++i)
            if(faceLabelLUT_[i] == mergedFace.label())
                faceLabelLUT_[i] = survivor.label();

        removeEdgeFromLabelImage(
            mergedEdge.scanLines(),
            *labelImage_, survivor.label(), associatedPixels);

//         survivor.pixelBounds_ |= mergedFace.pixelBounds_;
    }

    removeOne(node1.darts_,  removedDart.label());
    removeOne(node2.darts_, -removedDart.label());

    // remove singular nodes
    bool removeNode1 = !node1.degree();
    if(!node2.degree() && node2.label() != node1.label())
        removeIsolatedNode(node2);
    if(removeNode1)
        removeIsolatedNode(node1);

    if(survivor.areaValid_)
        survivor.area_ += mergedFace.area();
    survivor.pixelArea_ += mergedFace.pixelArea_;

    if(survivor.boundingBoxValid_)
        survivor.boundingBox_ |= mergedBBox;

    mergedEdge.uninitialize();
    mergedFace.uninitialize();

    postMergeFacesHook(survivor);

    if(associatedPixels.size())
        associatePixels(survivor, associatedPixels);

    return this->face(survivor.label());
}

/********************************************************************/
/*                                                                  */
/*                Python export code & API wrappers                 */
/*                                                                  */
/********************************************************************/

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
namespace bp = boost::python;

CELL_PTR(GeoMap::Edge) pySplitEdge(
    GeoMap &geomap, GeoMap::Edge &edge, unsigned int segmentIndex,
    bp::object newPoint)
{
    bp::extract<vigra::Vector2> insertPoint(newPoint);
    if(insertPoint.check())
        return geomap.splitEdge(edge, segmentIndex, insertPoint());
    return geomap.splitEdge(edge, segmentIndex);
}

/********************************************************************/

template<GeoMap::Dart &(GeoMap::Dart::*nextMethod)(void)>
class OrbitIterator
{
    GeoMap::Dart dart_, end_;
    bool atEnd_;

  public:
    typedef const GeoMap::Dart        value_type;
    typedef value_type               &reference;
    typedef value_type               *pointer;
    typedef std::forward_iterator_tag iterator_category;

    OrbitIterator(GeoMap::Dart const &dart)
    : dart_(dart),
      end_(dart),
      atEnd_(false)
    {
    }

    OrbitIterator & operator++()
    {
        if((dart_.*nextMethod)() == end_)
            atEnd_ = true;
        return *this;
    }

    OrbitIterator operator++(int)
    {
        OrbitIterator ret(*this);
        operator++();
        return ret;
    }

    bool atEnd() const
    {
        return atEnd_;
    }

    bool inRange() const
    {
        return !atEnd_;
    }

    reference operator*() const
    {
        return dart_;
    }

    pointer operator->() const
    {
        return &(operator*());
    }
};

typedef OrbitIterator<&GeoMap::Dart::nextPhi> PhiOrbitIterator;
typedef OrbitIterator<&GeoMap::Dart::nextSigma> SigmaOrbitIterator;
typedef OrbitIterator<&GeoMap::Dart::nextAlpha> AlphaOrbitIterator;

std::auto_ptr<PhiOrbitIterator>
phiOrbit(const GeoMap::Dart &dart)
{
    return std::auto_ptr<PhiOrbitIterator>(new PhiOrbitIterator(dart));
}

std::auto_ptr<SigmaOrbitIterator>
sigmaOrbit(const GeoMap::Dart &dart)
{
    return std::auto_ptr<SigmaOrbitIterator>(new SigmaOrbitIterator(dart));
}

std::auto_ptr<AlphaOrbitIterator>
alphaOrbit(const GeoMap::Dart &dart)
{
    return std::auto_ptr<AlphaOrbitIterator>(new AlphaOrbitIterator(dart));
}

template<class Iterator>
struct RangeIterAdapter
{
    typedef typename Iterator::value_type value_type;
    typedef typename Iterator::reference  reference;
    typedef typename Iterator::pointer    pointer;
    typedef std::forward_iterator_tag     iterator_category;

    RangeIterAdapter(const Iterator &begin, const Iterator &end)
    : iter_(begin),
      end_(end)
    {}

    bool inRange() const
    {
        return iter_ != end_;
    }

    bool atEnd() const
    {
        return iter_ == end_;
    }

    RangeIterAdapter &operator++()
    {
        ++iter_;
        return *this;
    }

    RangeIterAdapter operator++(int)
    {
        RangeIterAdapter ret(*this);
        operator++();
        return ret;
    }

    reference operator*() const
    {
        return *iter_;
    }

    pointer operator->() const
    {
        return &(operator*());
    }

    Iterator iter_, end_;
};

#ifdef USE_INSECURE_CELL_PTRS
// This is quite dangerous, *but*: The real lifetime of
// the referenced objects / cells are unknown, since any Euler
// operation might invalidate them.
#  define CELL_RETURN_POLICY bp::return_value_policy<reference_existing_object>
#else
#  define CELL_RETURN_POLICY bp::default_call_policies
#endif

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

/********************************************************************/

class SimpleCallback
{
  public:
    virtual ~SimpleCallback()
    {
        disconnect();
    }

    void disconnect()
    {
        for(unsigned int i = 0; i < connections_.size(); ++i)
            if(connections_[i].connected())
                connections_[i].disconnect();
        connections_.clear();
    }

  protected:
    std::vector<sigc::connection> connections_;
};

class RemoveNodeCallback : public SimpleCallback
{
  public:
    RemoveNodeCallback(GeoMap *geomap, bp::object callback)
    : callback_(callback)
    {
        connections_.push_back(
            geomap->removeNodeHook.connect(
                sigc::mem_fun(this, &RemoveNodeCallback::operator())));
    }

    bool operator()(GeoMap::Node &node)
    {
        return bp::extract<bool>(callback_(boost::ref(node)))();
    }

  protected:
    bp::object callback_;
};

class MergeEdgesCallbacks : public SimpleCallback
{
  public:
    MergeEdgesCallbacks(GeoMap *geomap,
                        bp::object preOpCallback,
                        bp::object postOpCallback)
    : preOpCallback_(preOpCallback),
      postOpCallback_(postOpCallback)
    {
        if(preOpCallback)
            connections_.push_back(
                geomap->preMergeEdgesHook.connect(
                    sigc::mem_fun(this, &MergeEdgesCallbacks::preMergeEdges)));
        if(postOpCallback)
            connections_.push_back(
                geomap->postMergeEdgesHook.connect(
                    sigc::mem_fun(this, &MergeEdgesCallbacks::postMergeEdges)));
    }

    bool preMergeEdges(const GeoMap::Dart &dart)
    {
        return bp::extract<bool>(
            preOpCallback_(boost::ref(dart)))();
    }

    void postMergeEdges(GeoMap::Edge &edge)
    {
        postOpCallback_(boost::ref(edge));
    }

  protected:
    bp::object preOpCallback_, postOpCallback_;
};

class SplitEdgeCallbacks : public SimpleCallback
{
  public:
    SplitEdgeCallbacks(GeoMap *geomap,
                       bp::object preOpCallback,
                       bp::object postOpCallback)
    : preOpCallback_(preOpCallback),
      postOpCallback_(postOpCallback)
    {
        if(preOpCallback)
            connections_.push_back(
                geomap->preSplitEdgeHook.connect(
                    sigc::mem_fun(this, &SplitEdgeCallbacks::preSplitEdge)));
        if(postOpCallback)
            connections_.push_back(
                geomap->postSplitEdgeHook.connect(
                    sigc::mem_fun(this, &SplitEdgeCallbacks::postSplitEdge)));
    }

    void preSplitEdge(GeoMap::Edge &edge, unsigned int segmentIndex,
                      vigra::Vector2 const &newPoint, bool insertPoint)
    {
        preOpCallback_(boost::ref(edge), segmentIndex,
                       insertPoint ? bp::object(newPoint) : bp::object());
    }

    void postSplitEdge(GeoMap::Edge &edge, GeoMap::Edge &newEdge)
    {
        postOpCallback_(boost::ref(edge), boost::ref(newEdge));
    }

  protected:
    bp::object preOpCallback_, postOpCallback_;
};

class RemoveBridgeCallbacks : public SimpleCallback
{
  public:
    RemoveBridgeCallbacks(GeoMap *geomap,
                          bp::object preOpCallback,
                          bp::object postOpCallback)
    : preOpCallback_(preOpCallback),
      postOpCallback_(postOpCallback)
    {
        if(preOpCallback)
            connections_.push_back(
                geomap->preRemoveBridgeHook.connect(
                    sigc::mem_fun(this, &RemoveBridgeCallbacks::preRemoveBridge)));
        if(postOpCallback)
            connections_.push_back(
                geomap->postRemoveBridgeHook.connect(
                    sigc::mem_fun(this, &RemoveBridgeCallbacks::postRemoveBridge)));
    }

    bool preRemoveBridge(const GeoMap::Dart &dart)
    {
        return bp::extract<bool>(
            preOpCallback_(boost::ref(dart)))();
    }

    void postRemoveBridge(GeoMap::Face &face)
    {
        postOpCallback_(boost::ref(face));
    }

  protected:
    bp::object preOpCallback_, postOpCallback_;
};

class MergeFacesCallbacks : public SimpleCallback
{
  public:
    MergeFacesCallbacks(GeoMap *geomap,
                        bp::object preOpCallback,
                        bp::object postOpCallback)
    : preOpCallback_(preOpCallback),
      postOpCallback_(postOpCallback)
    {
        if(preOpCallback)
            connections_.push_back(
                geomap->preMergeFacesHook.connect(
                    sigc::mem_fun(this, &MergeFacesCallbacks::preMergeFaces)));
        if(postOpCallback)
            connections_.push_back(
                geomap->postMergeFacesHook.connect(
                    sigc::mem_fun(this, &MergeFacesCallbacks::postMergeFaces)));
    }

    bool preMergeFaces(const GeoMap::Dart &dart)
    {
        return bp::extract<bool>(
            preOpCallback_(boost::ref(dart)))();
    }

    void postMergeFaces(GeoMap::Face &face)
    {
        postOpCallback_(boost::ref(face));
    }

  protected:
    bp::object preOpCallback_, postOpCallback_;
};

class AssociatePixelsCallback : public SimpleCallback
{
  public:
    AssociatePixelsCallback(GeoMap *geomap, bp::object callback)
    : callback_(callback)
    {
        connections_.push_back(
            geomap->associatePixelsHook.connect(
                sigc::mem_fun(this, &AssociatePixelsCallback::operator())));
    }

    void operator()(GeoMap::Face &face, const PixelList &pixels)
    {
        callback_(boost::ref(face), pixels);
    }

  protected:
    bp::object callback_;
};

std::auto_ptr<SimpleCallback>
addRemoveNodeCallback(GeoMap *geomap, bp::object callback)
{
    return std::auto_ptr<SimpleCallback>(
        new RemoveNodeCallback(geomap, callback));
}

std::auto_ptr<SimpleCallback>
addMergeEdgesCallbacks(GeoMap *geomap,
                       bp::object preOpCallback,
                       bp::object postOpCallback)
{
    return std::auto_ptr<SimpleCallback>(
        new MergeEdgesCallbacks(geomap, preOpCallback, postOpCallback));
}

std::auto_ptr<SimpleCallback>
addSplitEdgeCallbacks(GeoMap *geomap,
                       bp::object preOpCallback,
                       bp::object postOpCallback)
{
    return std::auto_ptr<SimpleCallback>(
        new SplitEdgeCallbacks(geomap, preOpCallback, postOpCallback));
}

std::auto_ptr<SimpleCallback>
addRemoveBridgeCallbacks(GeoMap *geomap,
                         bp::object preOpCallback,
                         bp::object postOpCallback)
{
    return std::auto_ptr<SimpleCallback>(
        new RemoveBridgeCallbacks(geomap, preOpCallback, postOpCallback));
}

std::auto_ptr<SimpleCallback>
addMergeFacesCallbacks(GeoMap *geomap,
                       bp::object preOpCallback,
                       bp::object postOpCallback)
{
    return std::auto_ptr<SimpleCallback>(
        new MergeFacesCallbacks(geomap, preOpCallback, postOpCallback));
}

std::auto_ptr<SimpleCallback>
addAssociatePixelsCallback(GeoMap *geomap,
                           bp::object callback)
{
    return std::auto_ptr<SimpleCallback>(
        new AssociatePixelsCallback(geomap, callback));
}

/********************************************************************/

// template<class OriginalImage>
// class FaceColorStatistics
// {
//     FaceColorStatistics(GeoMap &map, const OriginalImage &originalImage);

//     GeoMap &map_;
//     const OriginalImage &originalImage_;
// };

// template<class OriginalImage>
// FaceColorStatistics::FaceColorStatistics(
//     GeoMap &map, const OriginalImage &originalImage)
// : map_(map),
//   originalImage_(originalImage)
// {
// }

/********************************************************************/

using namespace boost::python;

std::auto_ptr<GeoMap>
createGeoMap(bp::list nodePositions,
             bp::list edgeTuples,
             vigra::Size2D imageSize)
{
    std::auto_ptr<GeoMap> result(new GeoMap(imageSize));

    for(int i = 0; i < len(nodePositions); ++i)
    {
        bp::extract<Vector2> ve(nodePositions[i]);
        if(ve.check())
            result->addNode(ve(), i);
    }

    if(edgeTuples)
        vigra_precondition(edgeTuples[0] == bp::object(),
            "GeoMap.__init__: edgeTuples[0] must be None (valid labels start at 1)");

    for(int i = 1; i < len(edgeTuples); ++i)
    {
        bp::extract<bp::tuple> ete(edgeTuples[i]);
        if(ete.check())
        {
            bp::tuple edgeTuple(ete());
            // FIXME: check length of tuple
            bp::extract<Vector2Array> pe(edgeTuple[2]);
            if(!pe.check())
                bp::throw_type_error(
                    "GeoMap.__init__: edge geometry not convertable to Vector2Array");
            CellLabel startNodeLabel = bp::extract<CellLabel>(edgeTuple[0])();
            CellLabel endNodeLabel   = bp::extract<CellLabel>(edgeTuple[1])();
            result->addEdge(startNodeLabel, endNodeLabel, pe(), i);
        }
    }

    return result;
}

GeoMap::FaceIterator faceIter(GeoMap &geomap, bool skipInfinite)
{
    GeoMap::FaceIterator result(geomap.facesBegin());
    if(skipInfinite)
        ++result;
    return result;
}

typedef RangeIterAdapter<GeoMap::Face::ContourIterator> ContourRangeIterator;

ContourRangeIterator
faceContours(const GeoMap::Face &face)
{
    return ContourRangeIterator(face.contoursBegin(), face.contoursEnd());
}

ContourRangeIterator
faceHoleContours(const GeoMap::Face &face)
{
    return ContourRangeIterator(face.holesBegin(), face.contoursEnd());
}

bp::object
labelImage(const GeoMap &map)
{
    if(!map.hasLabelImage())
        return bp::object();
    vigra::PythonGrayImage result(map.imageSize());
    copyImage(srcIterRange(map.labelsUpperLeft(), map.labelsLowerRight(),
                           map.labelAccessor()),
              destImage(result));
    return bp::object(result);
}

bp::list GeoMap_sortEdgesEventually(
    GeoMap &map, double stepDist, double minDist, bool splitEdges)
{
    GeoMap::UnsortableGroups unsortable;
    map.sortEdgesEventually(stepDist, minDist, unsortable, splitEdges);
    bp::list result;
    for(GeoMap::UnsortableGroups::iterator it = unsortable.begin();
        it != unsortable.end(); ++it)
    {
        std::cerr << "copying unsortable group (" << it->size() << " darts): ";
        bp::list unsortableGroup;
        for(GeoMap::UnsortableGroups::value_type::iterator dit = it->begin();
            dit != it->end(); ++dit)
        {
            std::cerr << ".";
            unsortableGroup.append(*dit);
        }
        std::cerr << "\n";
        result.append(unsortableGroup);
    }
    return result;
}

struct GeoMapPickleSuite : bp::pickle_suite
{
    static bp::tuple getinitargs(GeoMap &map)
    {
        bp::list nodePositions;
        for(GeoMap::NodeIterator it = map.nodesBegin(); it.inRange(); ++it)
        {
            for(unsigned int i = 0; i < (len(nodePositions)- (*it)->label()); ++i)
                nodePositions.append(bp::object());
            nodePositions.append((*it)->position());
        }

        bp::list edgeTuples;
        for(GeoMap::EdgeIterator it = map.edgesBegin(); it.inRange(); ++it)
        {
            for(unsigned int i = 0; i < (len(edgeTuples)- (*it)->label()); ++i)
                edgeTuples.append(bp::object());
            edgeTuples.append(
                bp::make_tuple(
                    (*it)->startNodeLabel(),
                    (*it)->endNodeLabel(),
                    Polygon(**it)));
        }

        return bp::make_tuple(nodePositions, edgeTuples, map.imageSize());
    }

    static bp::tuple getstate(GeoMap &map)
    {
        bp::list pySigmaMapping;
        std::auto_ptr<GeoMap::SigmaMapping> sigmaMapping(map.sigmaMapping());
        for(GeoMap::SigmaMapping::iterator it = sigmaMapping->begin();
            it != sigmaMapping->end(); ++it)
        {
            pySigmaMapping.append(*it);
        }
        return bp::make_tuple(
            pySigmaMapping,
            map.edgesSorted(),
            map.mapInitialized(),
            map.hasLabelImage());
    }

    static void setstate(GeoMap &map, bp::tuple state)
    {
        bp::list pySigmaMapping((
            bp::extract<bp::list>(state[0])()));
        bool edgesSorted = bp::extract<bool>(state[1])();
        bool mapInitialized = bp::extract<bool>(state[2])();
        bool initLabelImage = bp::extract<bool>(state[3])();

        GeoMap::SigmaMapping sigmaMapping(bp::len(pySigmaMapping));
        unsigned int i = 0;
        for(GeoMap::SigmaMapping::iterator it = sigmaMapping.begin();
            it != sigmaMapping.end(); ++it, ++i)
        {
            *it = bp::extract<int>(pySigmaMapping[i])();
        }

        map.setSigmaMapping(sigmaMapping, edgesSorted);
        if(mapInitialized)
            map.initializeMap(initLabelImage);
    }
};

std::string Node__repr__(GeoMap::Node const &node)
{
    std::stringstream s;
    s.precision(3);
    s << "<GeoMap.Node " << node.label() << " at " << node.position()
      << ", degree " << node.degree() << ">";
    return s.str();
}

std::string Edge__repr__(GeoMap::Edge const &edge)
{
    std::stringstream s;
    s.unsetf(std::ios::scientific);
    s.precision(1);
    s << "<GeoMap.Edge " << edge.label()
      << ", node " << edge.startNodeLabel() << " -> " << edge.endNodeLabel()
      << ", partial area " << edge.partialArea() << ", length " << edge.length()
      << ", " << edge.size() << " points>";
    return s.str();
}

std::string Face__repr__(GeoMap::Face const &face)
{
    std::stringstream s;
    s << "<GeoMap.Face " << face.label() << ", "
      << (face.contoursEnd() - face.holesBegin())
      << " holes, area " << face.area() << ">";
    return s.str();
}

std::string Dart__repr__(GeoMap::Dart const &dart)
{
    std::stringstream s;
    s << "<GeoMap.Dart " << dart.label();
    if(!dart.edge())
        s << ", invalid!>";
    else
        s << ", node " << dart.startNodeLabel() << " -> " << dart.endNodeLabel()
          << ", faces " << dart.leftFaceLabel() << "(l), " << dart.rightFaceLabel()
          << "(r)>";
    return s.str();
}

template<class T>
T returnCopy(const T &v)
{
    return v;
}

void defMap()
{
    CELL_RETURN_POLICY crp;

    {
        scope geoMap(
            class_<GeoMap, boost::noncopyable>(
                "GeoMap", init<vigra::Size2D>(arg("imageSize")))
            .def("__init__", make_constructor(
                     &createGeoMap, default_call_policies(),
                     (arg("nodePositions") = list(),
                      arg("edgeTuples") = list(),
                      arg("imageSize") = vigra::Size2D(0, 0))))
            .def("node", &GeoMap::node, crp)
            .def("nodeIter", &GeoMap::nodesBegin)
            .def("edge", &GeoMap::edge, crp)
            .def("edgeIter", &GeoMap::edgesBegin)
            .def("face", &GeoMap::face, crp)
            .def("faceIter", &faceIter, arg("skipInfinite") = false)
            .def("dart", &GeoMap::dart)
            .def("faceAt", &GeoMap::faceAt, crp)
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
            .def("addNode", (CELL_PTR(GeoMap::Node) (GeoMap::*)(const vigra::Vector2 &))&GeoMap::addNode, crp, args("position"))
            .def("addNode", (CELL_PTR(GeoMap::Node) (GeoMap::*)(const vigra::Vector2 &, CellLabel))&GeoMap::addNode, crp, args("position", "label"))
            .def("addEdge", (CELL_PTR(GeoMap::Edge)(GeoMap::*)(CellLabel, CellLabel, const Vector2Array &, CellLabel))&GeoMap::addEdge, crp,
                 (arg("startNodeLabel"), arg("endNodeLabel"), arg("points"),
                  arg("label") = 0))
            .def("addEdge", (CELL_PTR(GeoMap::Edge)(GeoMap::*)(GeoMap::Dart, GeoMap::Dart, const Vector2Array &))&GeoMap::addEdge, crp,
                 (arg("startNeighbor"), arg("endNeighbor"), arg("points")))
            .def("removeEdge", &GeoMap::removeEdge, crp)
            .def("splitEdge", &pySplitEdge,
                 (arg("edge"), arg("segmentIndex"), arg("newPoint") = object()))

            .def("removeIsolatedNode", &GeoMap::removeIsolatedNode, arg("dart"))
            .def("mergeEdges", &GeoMap::mergeEdges, arg("dart"), crp)
            .def("removeBridge", &GeoMap::removeBridge, arg("dart"), crp)
            .def("mergeFaces", &GeoMap::mergeFaces, arg("dart"), crp)
            .def("sortEdgesDirectly", &GeoMap::sortEdgesDirectly)
            .def("sortEdgesEventually", &GeoMap_sortEdgesEventually,
                 (arg("stepDist"), arg("minDist"), arg("splitEdges") = true))
            .def("edgesSorted", &GeoMap::edgesSorted)
            .def("splitParallelEdges", &GeoMap::splitParallelEdges)
            .def("initializeMap", &GeoMap::initializeMap, (arg("initLabelImage") = true))
            .def("mapInitialized", &GeoMap::mapInitialized)
            .def("hasLabelImage", &GeoMap::hasLabelImage)
            .def("labelImage", &labelImage)
            .def("nearestNode", &GeoMap::nearestNode, crp,
                 (arg("position"), arg(
                     "maxSquaredDist") = vigra::NumericTraits<double>::max()))
            .def("checkConsistency", &GeoMap::checkConsistency)
            .def_pickle(GeoMapPickleSuite())
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
            .def("__repr__", &Node__repr__)
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
            .def("dart", &GeoMap::Edge::dart)
            .def("isBridge", &GeoMap::Edge::isBridge)
            .def("isLoop", &GeoMap::Edge::isLoop)
            .def("scanLines", &GeoMap::Edge::scanLines,
                 return_value_policy<copy_const_reference>())
            .def("protection", &GeoMap::Edge::protection)
            .def("protect", &GeoMap::Edge::protect,
                 (arg("flag"), arg("onoff") = true))
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Edge__repr__)
        ;

        class_<GeoMap::Face, boost::noncopyable>("Face", init<GeoMap *, GeoMap::Dart>())
            .def("initialized", &GeoMap::Face::initialized)
            .def("label", &GeoMap::Face::label)
            .def("boundingBox", &GeoMap::Face::boundingBox,
                 return_value_policy<copy_const_reference>())
            .def("contains", &GeoMap::Face::contains)
            .def("area", &GeoMap::Face::area)
            .def("pixelArea", &GeoMap::Face::pixelArea)
            .def("contour", &GeoMap::Face::contour,
                 return_internal_reference<>(),
                 arg("index") = 0)
            .def("contours", &faceContours)
            .def("holeContours", &faceHoleContours)
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Face__repr__)
        ;

        RangeIterWrapper<ContourRangeIterator>("ContourIterator");

        return_internal_reference<> rself; // "return self" policy

        RangeIterWrapper<DartPointIter>("DartPointIter");

        class_<GeoMap::Dart>("Dart", no_init)
            .def(init<GeoMap *, int>())
            .def("clone", &GeoMap::Dart::clone)
            .def("label", &GeoMap::Dart::label)
            .def("map", &GeoMap::Dart::map,
                 return_value_policy<reference_existing_object>())
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
            .def("nextSigma", &GeoMap::Dart::nextSigma, rself)
            .def("prevSigma", &GeoMap::Dart::prevSigma, rself)
            .def("nextPhi", &GeoMap::Dart::nextPhi, rself)
            .def("prevPhi", &GeoMap::Dart::prevPhi, rself)
            .def("phiOrbit", &phiOrbit)
            .def("sigmaOrbit", &sigmaOrbit)
            .def("alphaOrbit", &alphaOrbit)
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Dart__repr__)
        ;

        RangeIterWrapper<PhiOrbitIterator>("PhiOrbitIterator");
        register_ptr_to_python< std::auto_ptr<PhiOrbitIterator> >();
        RangeIterWrapper<SigmaOrbitIterator>("SigmaOrbitIterator");
        register_ptr_to_python< std::auto_ptr<SigmaOrbitIterator> >();
        RangeIterWrapper<AlphaOrbitIterator>("AlphaOrbitIterator");
        register_ptr_to_python< std::auto_ptr<AlphaOrbitIterator> >();

#ifndef USE_INSECURE_CELL_PTRS
        register_ptr_to_python< CELL_PTR(GeoMap::Node) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Edge) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Face) >();
#endif

        class_<SimpleCallback>("SimpleCallback", no_init)
            .def("disconnect", &SimpleCallback::disconnect)
        ;

        def("addRemoveNodeCallback", &addRemoveNodeCallback);
        def("addMergeEdgesCallbacks", &addMergeEdgesCallbacks);
        def("addSplitEdgeCallbacks", &addSplitEdgeCallbacks);
        def("addRemoveBridgeCallbacks", &addRemoveBridgeCallbacks);
        def("addMergeFacesCallbacks", &addMergeFacesCallbacks);
        def("addAssociatePixelsCallback", &addAssociatePixelsCallback);
        register_ptr_to_python< std::auto_ptr<SimpleCallback> >();
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

    class_<DartPosition>("DartPosition", init<GeoMap::Dart>())
        .def("atEnd", &DartPosition::atEnd)
        .def("__call__", &DartPosition::operator(),
             return_value_policy<copy_const_reference>())
        .def("__copy__", &returnCopy<DartPosition>)
        .def("dartLabel", &DartPosition::dartLabel)
        .def("segmentIndex", &DartPosition::segmentIndex)
        .def("arcLength", &DartPosition::arcLength)
        .def("segmentStart", &DartPosition::segmentStart,
             return_value_policy<copy_const_reference>())
        .def("segmentEnd", &DartPosition::segmentEnd,
             return_value_policy<copy_const_reference>())
        .def("segmentLength", &DartPosition::segmentLength)
        .def("gotoArcLength", &DartPosition::gotoArcLength)
        .def("gotoNextSegment", &DartPosition::gotoNextSegment)
        .def("gotoPrevSegment", &DartPosition::gotoPrevSegment)
        .def("leaveCircle", &DartPosition::leaveCircle)
        .def("intersectCircle", &DartPosition::intersectCircle)
    ;
}
