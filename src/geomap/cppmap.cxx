#include "cppmap.hxx"
#include <vigra/tinyvector.hxx>
#include <vigra/copyimage.hxx>
#include <vigra/crackconnections.hxx>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
inline int isnan(double t) { return _isnan(t); }
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

class GeoMap::Node : boost::noncopyable
{
  protected:
    GeoMap        *map_;
    CellLabel      label_;
    vigra::Vector2 position_;
    int            anchor_;

    friend class GeoMap; // give access to anchor_ (add edge, sort edges, Euler..)
    friend class SigmaAnchor; // give access to anchor_

  public:
    Node(GeoMap *map, const vigra::Vector2 &position)
    : map_(map),
      label_(map->nodes_.size()),
      position_(position),
      anchor_(0)
    {
        map_->nodes_.push_back(GeoMap::Nodes::value_type(this));
        ++map_->nodeCount_;
        map_->nodeMap_.insert(PositionedNodeLabel(position_, label_));
    }

    bool initialized() const
    {
        return map_ != NULL;
    }

  protected:
    void uninitialize()
    {
        GeoMap *map = map_;
        map_ = NULL; // DON'T MESS WITH THIS!
        --map->nodeCount_;
        map->nodeMap_.erase(
            map->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                                   vigra::NumericTraits<double>::epsilon()));
        RESET_PTR(map->nodes_[label_]); // may have effect like "delete this;"!
    }

  public:
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

    bool isIsolated() const
    {
        return !anchor_;
    }

    inline unsigned int degree() const;

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
: public vigra::BBoxPolygon<vigra::Vector2>, boost::noncopyable
{
  public:
    typedef vigra::BBoxPolygon<vigra::Vector2> Base;

  protected:
    GeoMap      *map_;
    CellLabel    label_;
    CellLabel    startNodeLabel_, endNodeLabel_;
    CellLabel    leftFaceLabel_, rightFaceLabel_;
    unsigned int flags_;

    mutable std::auto_ptr<vigra::Scanlines> scanLines_;

    friend class Dart; // allow setLeftFaceLabel
    friend class GeoMap;

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
      flags_(0),
      scanLines_(NULL)
    {
        map_->edges_.push_back(GeoMap::Edges::value_type(this));
        ++map_->edgeCount_;
    }

    bool initialized() const
    {
        return map_ != NULL;
    }

  protected:
    void uninitialize()
    {
        GeoMap *map = map_;
        map_ = NULL;
        --map->edgeCount_;
        RESET_PTR(map->edges_[label_]); // may have effect like "delete this;"
    }

  public:
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

    unsigned int flags() const
    {
        return flags_;
    }

    unsigned int flag(unsigned int which) const
    {
        return flags_ & which;
    }

    void setFlag(unsigned int flag, bool onoff = true)
    {
        if(onoff)
            flags_ |= flag;
        else
            flags_ &= ~flag;
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
        label_ = map_->sigmaMapping_[label_];
        return *this;
    }

    Dart &prevSigma()
    {
        label_ = map_->sigmaInverseMapping_[label_];
        return *this;
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

inline unsigned int GeoMap::Node::degree() const
{
    if(!anchor_)
        return 0;

    int result = 0;
    GeoMap::Dart d(map_, anchor_);
    do
    {
        ++result;
    }
    while(d.nextSigma().label() != anchor_);
    return result;
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

class GeoMap::Face : boost::noncopyable
{
  public:
    typedef Edge::BoundingBox BoundingBox;
    typedef std::vector<Dart> Contours;
    typedef Contours::const_iterator ContourIterator;

  protected:
    GeoMap              *map_;
    CellLabel            label_;
    std::vector<Dart>    anchors_;
    mutable unsigned int flags_;
    mutable BoundingBox  boundingBox_;
    mutable double       area_;
    unsigned int         pixelArea_;

    static const unsigned int BOUNDING_BOX_VALID = 0x8000000;
    static const unsigned int AREA_VALID         = 0x4000000;

    friend class GeoMap; // give access to pixelArea_ and anchors_ (Euler ops...)

    unsigned int findComponentAnchor(const GeoMap::Dart &dart);

  public:
    Face(GeoMap *map, Dart anchor)
    : map_(map),
      label_(map->faces_.size()),
      flags_(0),
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

  protected:
    void uninitialize()
    {
        GeoMap *map = map_;
        map_ = NULL;
        --map->faceCount_;
        RESET_PTR(map->faces_[label_]); // may have effect like "delete this;"
    }

  public:
    CellLabel label() const
    {
        return label_;
    }

    const BoundingBox &boundingBox() const
    {
        vigra_precondition(label_, "infinite face has no boundingBox()!");

        if(!flag(BOUNDING_BOX_VALID))
        {
            boundingBox_ = BoundingBox();
            Dart anchor(anchors_[0]), dart(anchor);
            do
            {
                boundingBox_ |= dart.edge()->boundingBox();
            }
            while(dart.nextPhi() != anchor);
            flags_ |= BOUNDING_BOX_VALID;
        }
        return boundingBox_;
    }

    bool contains(const Vector2 &point) const
    {
        vigra_precondition(initialized(), "contains() of uninitialized face!");
        if(map_->labelImage_)
        {
            IVector2 iPos(intVPos(point));
            if(map_->labelImage_->isInside(iPos))
            {
                int l = (*map_->labelImage_)[iPos];
                if(l > 0 && (map_->faceLabelLUT_[l] == label_))
                    return true;
            }
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
        if(!flag(AREA_VALID))
        {
            area_ = 0.0;
            for(unsigned int i = 0; i < anchors_.size(); ++i)
            {
                area_ += contourArea(anchors_[i]);
            }
            flags_ |= AREA_VALID;
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

        if(flag(AREA_VALID))
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

    unsigned int flags() const
    {
        return flags_;
    }

    unsigned int flag(unsigned int which) const
    {
        return flags_ & which;
    }

    void setFlag(unsigned int flag, bool onoff = true)
    {
        if(onoff)
            flags_ |= flag;
        else
            flags_ &= ~flag;
    }

    GeoMap *map() const
    {
        return map_;
    }

  private:
    Face(const Face &) {} // disallow copying
    Face &operator=(const Face &) { return *this; }
};

class GeoMap::SigmaAnchor
{
  public:
    SigmaAnchor(const GeoMap::Node &node)
    : isSingular_(node.isIsolated()),
      dartLabel_(node.anchor_),
      nodeLabel_(node.label()),
      map_(node.map())
    {
        vigra_precondition(isSingular_ || !node.map()->mapInitialized(),
            "sigma position of initialized GeoMap not fully specified");
    }

    SigmaAnchor(const GeoMap::Dart &dart)
    : isSingular_(false),
      dartLabel_(dart.label()),
      nodeLabel_(dart.startNodeLabel()),
      map_(dart.map())
    {
    }

    bool isSingular() const
    {
        return isSingular_;
    }

    int dartLabel() const
    {
        return dartLabel_;
    }

    CellLabel nodeLabel() const
    {
        return nodeLabel_;
    }

    bool operator==(const GeoMap::SigmaAnchor &other) const
    {
        if(isSingular() != other.isSingular())
            return false;
        if(isSingular())
            return nodeLabel_ == other.nodeLabel_;
        else
            return dartLabel_ == other.dartLabel_;
    }

  private:
    bool isSingular_;
    int dartLabel_;
    CellLabel nodeLabel_;
    GeoMap *map_;
};

void GeoMap::Node::setPosition(const vigra::Vector2 &p)
{
    vigra_precondition(initialized(), "setPosition() of uninitialized node!");
    map_->nodeMap_.erase(
        map_->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                               vigra::NumericTraits<double>::epsilon()));
    position_ = p;

    GeoMap::Dart d(map_, anchor_);
    do
    {
        if(d.label() > 0)
        {
            (*map_->edge(d.label()))[ 0] = p;
        }
        else
        {
            GeoMap::Edge &edge(*map_->edge(-d.label()));
            edge[edge.size()-1] = p;
        }
    }
    while(d.nextSigma().label() != anchor_);

    map_->nodeMap_.insert(PositionedNodeLabel(p, label_));
}

inline GeoMap::Dart GeoMap::Node::anchor() const
{
    vigra_precondition(initialized(), "anchor() of uninitialized node!");
    vigra_precondition(anchor_ != 0, "anchor() of degree 0 node!");
    return Dart(map_, anchor_);
}

GeoMap::GeoMap(vigra::Size2D imageSize)
: sigmaMappingArray_(101, 0),
  sigmaInverseMappingArray_(101, 0),
  sigmaMapping_(sigmaMappingArray_.begin() + 50),
  sigmaInverseMapping_(sigmaInverseMappingArray_.begin() + 50),
  nodeCount_(0),
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
    vigra_precondition(mapInitialized(),
        "faceAt() called on graph (mapInitialized() == false)!");

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
    const GeoMap::SigmaAnchor &startNeighbor,
    const GeoMap::SigmaAnchor &endNeighbor,
    const Vector2Array &points, CellLabel label)
{
    if(label > edges_.size())
        edges_.resize(label, NULL_PTR(GeoMap::Edge));
    GeoMap::Edge *result = new GeoMap::Edge(
        this, startNeighbor.nodeLabel(),  endNeighbor.nodeLabel(), points);

    if(startNeighbor.isSingular())
    {
        insertSigmaPredecessor(result->startNode()->anchor_, (int)result->label());
        result->startNode()->anchor_ = (int)result->label();
    }
    else
        insertSigmaPredecessor(startNeighbor.dartLabel(), (int)result->label());

    if((endNeighbor == startNeighbor) && result->partialArea() < 0)
    {
        insertSigmaPredecessor((int)result->label(), -(int)result->label());
    }
    else if(endNeighbor.isSingular())
    {
        insertSigmaPredecessor(result->endNode()->anchor_, -(int)result->label());
        result->endNode()->anchor_ = -(int)result->label();
    }
    else
        insertSigmaPredecessor(endNeighbor.dartLabel(), -(int)result->label());

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

        detachDart( dart.label());
        detachDart(-dart.label());

        dart.edge()->uninitialize();
    }
}

void GeoMap::sortEdgesDirectly()
{
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    typedef std::pair<double, int> DartAngle;

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if((*it)->isIsolated())
            continue;

        std::vector<DartAngle> dartAngles;

        GeoMap::Dart anchor((*it)->anchor()), d(anchor);
        do
        {
            vigra_precondition(
                d.size() >= 2, "cannot measure angle of darts with < 2 points!");
            dartAngles.push_back(
                DartAngle(angleTheta(-d[1][1] + d[0][1],
                                      d[1][0] - d[0][0]),
                          d.label()));
        }
        while(d.nextSigma() != anchor);

        std::sort(dartAngles.begin(), dartAngles.end());

        int predecessor = dartAngles[0].second;
        for(unsigned int i = 1; i < dartAngles.size(); ++i)
        {
            vigra_precondition(
                dartAngles[i-1].first != dartAngles[i].first,
                "sortEdgesDirectly: edges leave node at identical angles!");

            sigmaMapping_[predecessor] = dartAngles[i].second;
            sigmaInverseMapping_[dartAngles[i].second] = predecessor;
            predecessor = dartAngles[i].second;
        }
        sigmaMapping_[predecessor] = dartAngles[0].second;
        sigmaInverseMapping_[dartAngles[0].second] = predecessor;
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
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    double minAngle = std::atan2(minDist, stepDist),
          stepDist2 = vigra::sq(stepDist);

    if(splitEdges)
        splitInfo_ = std::auto_ptr<PlannedSplits>(new PlannedSplits());

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if((*it)->isIsolated())
            continue;

        DartPositionAngles dartPositions;

        GeoMap::Dart anchor((*it)->anchor()), d(anchor);
        do
        {
            dartPositions.push_back(DartPositionAngle(d));
        }
        while(d.nextSigma() != anchor);

        sortEdgesInternal((*it)->position(), 0.0,
                          dartPositions.begin(), dartPositions.end(),
                          stepDist2, minAngle,
                          unsortable, splitInfo_.get(), false);

        int predecessor = dartPositions[0].dp.dartLabel();
        for(unsigned int i = 1; i < dartPositions.size(); ++i)
        {
            int successor = dartPositions[i].dp.dartLabel();
            sigmaMapping_[predecessor] = successor;
            sigmaInverseMapping_[successor] = predecessor;
            predecessor = successor;
        }
        sigmaMapping_[predecessor] = dartPositions[0].dp.dartLabel();
        sigmaInverseMapping_[dartPositions[0].dp.dartLabel()] = predecessor;
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
        CellLabel newEdgeLabel =
            splitEdge(*edge(abs(it->dartLabel)),
                      it->segmentIndex, it->position)->label();

        PlannedSplits::difference_type &pos(
            groupPositions[it->splitGroup]);

        mergeDarts[pos] =
            MergeDart((int)newEdgeLabel, it->dartLabel > 0, it->sigmaPos);

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

            // intersect checkSurvivorDist-circle with dart
            DartPosition dp1(d);
            dp1.leaveCircle(nodePos, checkSurvivorDist2);
            d.nextSigma();
            DartPosition dp2(d);
            dp2.leaveCircle(nodePos, checkSurvivorDist2);

            // determine vectors between node pos. & intersections..
            vigra::Vector2
                v1(dp1() - nodePos),
                v2(nodePos - dp2());

            // ..and choose dart with smallest enclosed angle:
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
        int sigmaCenterLabel(survivor.clone().nextSigma().label()),
            sigmaNeighborLabel(sigmaCenterLabel);

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

                checkConsistency(); // no modification should've happened so far

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

                detachDart(relocateDart.label());
                insertSigmaPredecessor(sigmaNeighborLabel, relocateDart.label());
                sigmaNeighborLabel = relocateDart.label();

                // remove mergeDart and its startNode:
                GeoMap::Node &mergedNode(*mergeDart.startNode());
                removeEdge(mergeDart);
                removeIsolatedNode(mergedNode);
            }
            while(it != mergeDartsGroupBegin);
        }

        sigmaNeighborLabel = sigmaMapping_[sigmaCenterLabel];

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
            detachDart(relocateDart.label());
            insertSigmaPredecessor(sigmaNeighborLabel, relocateDart.label());

            // remove mergeDart and its startNode:
            GeoMap::Node &mergedNode(*mergeDart.startNode());
            removeEdge(mergeDart);
            removeIsolatedNode(mergedNode);
        }
    }
}

void GeoMap::setSigmaMapping(SigmaMapping const &sigmaMapping, bool sorted)
{
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    vigra_precondition(sigmaMapping.size() >= (2*edges_.size() - 1),
                       "setSigmaMapping: sigmaMapping too small!");
    SigmaMapping::const_iterator sigma(
        sigmaMapping.begin() + (sigmaMapping.size() / 2));

    if(sigmaMappingArray_.size() < 2*edges_.size() - 1)
        resizeSigmaMapping(2*edges_.size() - 1);

    std::copy(sigma - (edges_.size() - 1), sigma + edges_.size(),
              sigmaMapping_ - (edges_.size() - 1));
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        int label = (int)(*it)->label();
        sigmaInverseMapping_[sigmaMapping_[ label]] =  label;
        sigmaInverseMapping_[sigmaMapping_[-label]] = -label;
    }

    edgesSorted_ = sorted;
}

void GeoMap::initializeMap(bool initLabelImage)
{
    vigra_precondition(!mapInitialized(),
                       "initializeMap() called more than once");
    if(!edgesSorted())
        sortEdgesDirectly();

    initContours();
    //std::cerr << faceCount_ << " contours found, embedding...\n";
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

void GeoMap::resizeSigmaMapping(SigmaMapping::size_type newSize)
{
    SigmaMapping
        newSigma(newSize, 0),
        newSigmaInverse(newSize, 0);
    std::copy(sigmaMappingArray_.begin(), sigmaMappingArray_.end(),
              newSigma.begin() + sigmaMappingArray_.size()/2);
    std::copy(sigmaInverseMappingArray_.begin(),
              sigmaInverseMappingArray_.end(),
              newSigmaInverse.begin() + sigmaMappingArray_.size()/2);

    std::swap(sigmaMappingArray_, newSigma);
    sigmaMapping_ = sigmaMappingArray_.begin() + newSize/2;

    std::swap(sigmaInverseMappingArray_, newSigmaInverse);
    sigmaInverseMapping_ = sigmaInverseMappingArray_.begin() + newSize/2;
}

void GeoMap::insertSigmaPredecessor(int successor, int newPredecessor)
{
    while(sigmaMappingArray_.size() < 2*(unsigned)abs(newPredecessor)+1)
        resizeSigmaMapping(2*sigmaMappingArray_.size()-1);

    if(!successor)
    {
        sigmaMapping_[newPredecessor] = newPredecessor;
        sigmaInverseMapping_[newPredecessor] = newPredecessor;
        return;
    }

    int oldPredecessor = sigmaInverseMapping_[successor];
    sigmaMapping_[oldPredecessor] = newPredecessor;
    sigmaMapping_[newPredecessor] = successor;
    sigmaInverseMapping_[successor] = newPredecessor;
    sigmaInverseMapping_[newPredecessor] = oldPredecessor;
}

void GeoMap::detachDart(int dartLabel)
{
    GeoMap::Node &node(*dart(dartLabel).startNode());

    int successor = sigmaMapping_[dartLabel];
    int predecessor = sigmaInverseMapping_[dartLabel];

    sigmaMapping_[predecessor] = successor;
    sigmaInverseMapping_[successor] = predecessor;

    if(node.anchor_ == dartLabel)
    {
        if(successor != dartLabel)
            node.anchor_ = successor;
        else
            node.anchor_ = 0;
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
    //std::cerr << "GeoMap[" << this << "].checkConsistency()\n";
    bool result = true;

    if((sigmaMapping_ !=
        sigmaMappingArray_.begin() + (sigmaMappingArray_.size()-1)/2) ||
       (sigmaInverseMapping_ !=
        sigmaInverseMappingArray_.begin() + (sigmaMappingArray_.size()-1)/2) ||
       !(sigmaMappingArray_.size() & 1) ||
       (sigmaMappingArray_.size() != sigmaInverseMappingArray_.size()))
    {
        std::cerr << "  Sigma mapping arrays not correctly setup!\n";
        std::cerr << "      sigma: size " << sigmaMappingArray_.size()
                  << ", center @" << (sigmaMapping_ - sigmaMappingArray_.begin())
                  << "\n";
        std::cerr << "    inverse: size " << sigmaInverseMappingArray_.size()
                  << ", center @" << (
                      sigmaInverseMapping_ - sigmaInverseMappingArray_.begin())
                  << "\n";
        result = false;
    }
    else
    {
        if(sigmaMappingArray_.size() < maxEdgeLabel()*2-1)
        {
            std::cerr << "  Sigma mapping arrays not large enough!\n";
            result = false;
        }

        int dist = sigmaMappingArray_.size()/2;
        for(int label = -dist; label <= dist; ++label)
        {
            if((unsigned)abs(sigmaMapping_[label]) >= maxEdgeLabel() ||
               (unsigned)abs(sigmaInverseMapping_[label]) >= maxEdgeLabel())
            {
                std::cerr << "  Sigma mapping arrays contain junk (sigma["
                          << label << "] = " << sigmaMapping_[label]
                          << ", sigma^-1[" << label << "] = "
                          << sigmaInverseMapping_[label] << ")!\n";
                result = false;
                continue;
            }
            if(sigmaMapping_[label] &&
               (CellLabel)abs(label) < maxEdgeLabel() &&
               edges_[abs(label)].get() &&
               sigmaInverseMapping_[sigmaMapping_[label]] != label)
            {
                std::cerr << "  Sigma inverse is not correct ("
                          << label << " -> " << sigmaMapping_[label]
                          << " -> " << sigmaInverseMapping_[sigmaMapping_[label]]
                          << ")!\n";
                result = false;
            }
        }
        if(!result)
            return result; // function may not terminate in this state
    }

    unsigned int actualNodeCount = 0, actualEdgeCount = 0, actualFaceCount = 0;

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        ++actualNodeCount;

        if((*it)->map() != this)
        {
            std::cerr << "  Node " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }

        if((*it)->isIsolated())
            continue;

        Dart anchor((*it)->anchor()), dart(anchor);
        do
        {
            if(!dart.edge().get())
            {
                std::cerr << "  Node " << (*it)->label()
                          << " has broken sigma orbit: dart "
                          << dart.label() << " does not exist!\n";

                std::cerr << "    orbit so far: "
                          << anchor.label() << " (anchor)";
                Dart d(anchor);
                while(d != dart)
                {
                    d.nextSigma();
                    std::cerr << " -> " << d.label();
                }
                std::cerr << "\n";

                result = false;
                break;
            }

            if(dart.startNodeLabel() != (*it)->label())
            {
                std::cerr << "  Node " << (*it)->label()
                          << " has broken sigma orbit:\n"
                          << "  contains Dart " << dart.label()
                          << " from node " << dart.startNodeLabel()
                          << " -> " << dart.endNodeLabel() << "\n";
                result = false;
                break;
            }
        }
        while(dart.nextSigma() != anchor);
    }
    if(actualNodeCount != nodeCount())
    {
        std::cerr << "  Node count wrong (" << nodeCount()
                  << ", should be " << actualNodeCount << ")!\n";
            result = false;
    }

    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        ++actualEdgeCount;

        if((*it)->map() != this)
        {
            std::cerr << "  Edge " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }

        if(mapInitialized())
        {
            if((*it)->isBridge() && (*it)->isLoop())
            {
                std::cerr << "  Edge " << (*it)->label()
                          << " is both loop and bridge!?\n";
                result = false;
            }

            if(!(*it)->leftFace() || !(*it)->rightFace())
            {
                std::cerr << "  Edge " << (*it)->label() << " has invalid faces:\n";
                if(!(*it)->leftFace())
                    std::cerr << "     left face (" << (*it)->leftFaceLabel()
                              << ") does not exist!\n";
                if(!(*it)->rightFace())
                    std::cerr << "    right face (" << (*it)->rightFaceLabel()
                              << ") does not exist!\n";
                result = false;
                continue; // no need to check geometry
            }
        }

        if((*it)->size() < 2)
        {
            std::cerr << "  Edge " << (*it)->label() << " is too short ("
                      << (*it)->size() << " points)!\n";
            result = false;
        }
        else
        {
            if((**it)[0] != (*it)->startNode()->position() ||
                (**it)[(*it)->size()-1] != (*it)->endNode()->position())
            {
                std::cerr << "  Edge " << (*it)->label()
                          << " has non-matching end positions:\n"

                          << "    start node " << (*it)->startNodeLabel() << " at "
                          << (*it)->startNode()->position() << " is "
                          << ((**it)[0] - (*it)->startNode()->position()).magnitude()
                          << " pixels from " << (**it)[0] << "\n"

                          << "    end node " << (*it)->endNodeLabel() << " at "
                          << (*it)->endNode()->position() << " is "
                          << ((**it)[(*it)->size()-1]
                              - (*it)->endNode()->position()).magnitude()
                          << " pixels from " << (**it)[(*it)->size()-1] << "\n";
                result = false;
            }
        }
    }
    if(actualEdgeCount != edgeCount())
    {
        std::cerr << "  Edge count wrong (" << edgeCount()
                  << ", should be " << actualEdgeCount << ")!\n";
            result = false;
    }

    if(!mapInitialized())
        return result; // cannot test Faces yet (do not exist)

    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
    {
        ++actualFaceCount;

        GeoMap::Face &face(**it);
        if(face.map() != this)
        {
            std::cerr << "  Face " << face.label() << " has wrong map()!\n";
            result = false;
            break;
        }

        typedef std::map<int, unsigned int> SeenAnchors;
        SeenAnchors seenAnchors;

        bool outerContour = face.label() > 0;
        for(GeoMap::Face::ContourIterator ci = face.contoursBegin();
            ci != face.contoursEnd(); ++ci, outerContour = false)
        {
            double area = 0.0;
            int canonicalAnchor(ci->label());
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

                if(!dart.edge()->isBridge())
                    area += dart.partialArea();

                if(ci->label() < canonicalAnchor)
                    canonicalAnchor = ci->label();
            }
            while(dart.nextPhi() != anchor);

            if(dart != anchor)
                continue; // loop canceled, area not valid

            if((area <= 0) == outerContour)
            {
                std::cerr << "  Face " << face.label() << " contains an "
                          << (outerContour
                              ? "outer anchor with negative area"
                              : "inner anchor with positive area")
                          << ":\n  contour from dart " << ci->label()
                          << " has area " << area << "!\n";
                result = false;
            }

            std::pair<SeenAnchors::iterator, bool>
                seenAnchor(seenAnchors.insert(
                               std::make_pair(canonicalAnchor,
                                              ci - face.contoursBegin())));
            if(!seenAnchor.second)
            {
                std::cerr << "  Face " << face.label()
                          << " contains duplicate anchors at (contour indices "
                          << (ci - face.contoursBegin()) << " and "
                          << seenAnchor.first->second << ")!\n";
                result = false;
            }
        }
    }
    if(actualFaceCount != faceCount())
    {
        std::cerr << "  Face count wrong (" << faceCount()
                  << ", should be " << actualFaceCount << ")!\n";
            result = false;
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

    vigra_assert((d1.leftFaceLabel() == d2.rightFaceLabel()) &&
                 (d2.leftFaceLabel() == d1.rightFaceLabel()),
                 "mergeEdges: broken map (left/rightFaceLabel)");

    if(d1.label() > 0 && d2.label() < 0)
        std::swap(d1, d2); // minimize number of reverse()s necessary

    if(mapInitialized())
    {
        GeoMap::Face *faces[2];
        faces[0] = &(*d2.leftFace());
        faces[1] = &(*d2.rightFace());
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
    }

    GeoMap::Edge &survivor(*d1.edge());
    GeoMap::Node &mergedNode(*d2.startNode());
    GeoMap::Edge &mergedEdge(*d2.edge());

    if(!preMergeEdgesHook(d1))
        return NULL_PTR(GeoMap::Edge);
    if(!removeNodeHook(mergedNode))
        return NULL_PTR(GeoMap::Edge);

    // TODO: history append?

    int successor = sigmaMapping_[-d2.label()];
    int predecessor = sigmaInverseMapping_[-d2.label()];

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

        survivor.endNodeLabel_ = d2.endNodeLabel();
    }
    else
    {
        survivor.reverse();
        if(mergedEdge.startNodeLabel() != mergedNode.label())
            mergedEdge.reverse();
        survivor.extend(mergedEdge);
        survivor.reverse();

        survivor.startNodeLabel_ = d2.endNodeLabel();
    }

    // replace -d2 with d1 within orbits / anchor:
    if(successor == -d2.label())
    {
        // no other dart at that end node
        sigmaMapping_[d1.label()] = d1.label();
        sigmaInverseMapping_[d1.label()] = d1.label();
        d2.endNode()->anchor_ = d1.label();
    }
    else
    {
        sigmaMapping_[predecessor] = d1.label();
        sigmaMapping_[d1.label()] = successor;
        sigmaInverseMapping_[successor] = d1.label();
        sigmaInverseMapping_[d1.label()] = predecessor;
        if(d2.endNode()->anchor_ == -d2.label())
            d2.endNode()->anchor_ = d1.label();
    }

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

    GeoMap::Node
        &newNode(*addNode(insertPoint ? newPoint : edge[segmentIndex])),
        &changedNode(*edge.endNode());

    if(labelImage_)
        rawAddEdgeToLabelImage(edge.scanLines(), *labelImage_, 1);

    if(insertPoint)
    {
        edge.insert(edge.begin() + segmentIndex + 1, newPoint);
        ++segmentIndex;
    }

    GeoMap::Edge *result = new GeoMap::Edge(
        this, newNode.label(), changedNode.label(), edge.split(segmentIndex));
    result->leftFaceLabel_ = edge.leftFaceLabel_;
    result->rightFaceLabel_ = edge.rightFaceLabel_;
    result->flags_ = edge.flags_;

    if(sigmaMappingArray_.size() < 2*result->label()+1)
        resizeSigmaMapping(2*sigmaMappingArray_.size()-1);

    int successor = sigmaMapping_[-(int)edge.label()];
    int predecessor = sigmaInverseMapping_[-(int)edge.label()];

    // setup newNode (fixed degree 2)
    sigmaMapping_[-(int)edge.label()] = (int)result->label();
    sigmaMapping_[(int)result->label()] = -(int)edge.label();
    sigmaInverseMapping_[(int)result->label()] = -(int)edge.label();
    sigmaInverseMapping_[-(int)edge.label()] = (int)result->label();
    newNode.anchor_ = (int)result->label();

    // edge now ends in newNode
    edge.endNodeLabel_ = newNode.label();

    if(successor == -(int)edge.label())
    {
        // changedNode has degree 1, replace -edgeLabel with -resultLabel
        sigmaMapping_[-(int)result->label()] = -(int)result->label();
        sigmaInverseMapping_[-(int)result->label()] = -(int)result->label();
        changedNode.anchor_ = -(int)result->label();
    }
    else
    {
        sigmaMapping_[predecessor] = -(int)result->label();
        sigmaMapping_[-(int)result->label()] = successor;
        sigmaInverseMapping_[successor] = -(int)result->label();
        sigmaInverseMapping_[-(int)result->label()] = predecessor;
        if(changedNode.anchor_ == -(int)edge.label())
            changedNode.anchor_ = -(int)result->label();
    }

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

    // remove both darts from both sigma orbits:
    detachDart( dart.label());
    detachDart(-dart.label());

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
    if(survivor.flag(GeoMap::Face::BOUNDING_BOX_VALID))
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

    // remove both darts from both sigma orbits:
    detachDart( removedDart.label());
    detachDart(-removedDart.label());

    // remove singular nodes
    if(node1.isIsolated() && node2.label() != node1.label())
        removeIsolatedNode(node1);
    if(node2.isIsolated())
        removeIsolatedNode(node2);

    if(survivor.flag(GeoMap::Face::AREA_VALID))
        survivor.area_ += mergedFace.area();
    survivor.pixelArea_ += mergedFace.pixelArea_;

    if(survivor.flag(GeoMap::Face::BOUNDING_BOX_VALID))
        survivor.boundingBox_ |= mergedBBox;

    mergedEdge.uninitialize();
    mergedFace.uninitialize();

    postMergeFacesHook(survivor);

    if(associatedPixels.size())
        associatePixels(survivor, associatedPixels);

    return this->face(survivor.label());
}

/********************************************************************/

// template<class OriginalImage>
// class FaceColorStatistics
// {
//     template<int SplineOrder>
//     FaceColorStatistics(GeoMap &map, const OriginalImage &originalImage,
//                         int minSampleCount = 1);

//     GeoMap &map_;
//     const OriginalImage &originalImage_;
// };

// template<class OriginalImage>
// template<int SplineOrder>
// FaceColorStatistics<OriginalImage>::FaceColorStatistics(
//     GeoMap &map, const OriginalImage &originalImage,
//     int minSampleCount)
// : map_(map),
//   originalImage_(originalImage)
// {

// }
