#ifndef VIGRA_CPPMAP_HXX
#define VIGRA_CPPMAP_HXX

#include "filteriterator.hxx"
#include "labellut.hxx"
#include "vigra/positionedmap.hxx"
#include "vigra/polygon.hxx"
#include <vector>
#include <list>
#include <vigra/multi_array.hxx>
#include <sigc++/sigc++.h>
#include <boost/utility.hpp> // boost::noncopyable

#include <cfloat>

#ifdef _MSC_VER
inline int isnan(double t) { return _isnan(t); }
#else
#include <math.h>
// isnan is C99, Linux manpage says:
// "Compile with -std=c99; link with -lm."
// this is a workaround according to P.J. Plauger (Dinkumware):
// # define cppmap_isnan_workaround(x) ((x) != (x))
#endif

// The define USE_INSECURE_CELL_PTRS can be used to switch between
// "safe" cell handling e.g. for Python and a possibly faster C++ way.

#ifdef USE_INSECURE_CELL_PTRS
#  define CELL_PTR(Type) Type *
#  define NULL_PTR(Type) (Type *)NULL
#  define RESET_PTR(ptr) delete ptr; ptr = NULL
#else
#  include <boost/shared_ptr.hpp>
#  define CELL_PTR(Type) boost::shared_ptr<Type>
#  define NULL_PTR(Type) boost::shared_ptr<Type>()
#  define RESET_PTR(ptr) ptr.reset()
#endif

typedef unsigned int CellLabel;

// functor for FilterIterator to skip NULL cells
template<class POINTER>
struct NotNull
{
    bool operator()(const POINTER &p) const
    {
        return p;
    }
};

typedef vigra::TinyVector<double, 2>       Vector2;
typedef vigra::PointArray<Vector2>         Vector2Array;
typedef vigra::BBoxPolygon<vigra::Vector2> Polygon;

typedef vigra::PointArray<vigra::Point2D> PixelList;

// "accumulator" for libsigc++, which calls pre-operation
// callbacks in order but cancels whenever a callback returns false
struct interruptable_accumulator
{
    typedef bool result_type;

    template<typename T_iterator>
    result_type operator()(T_iterator first, T_iterator last) const
    {
        for(; first != last; ++first)
            if(!*first)
                return false;
        return true;
    }
};

namespace detail {

typedef vigra::MultiArray<2, int>::difference_type IVector2;

inline IVector2 intVPos(const Vector2 &p)
{
    return IVector2((int)floor(p[0]+0.5), (int)floor(p[1]+0.5));
}

class PlannedSplits;

} // namespace detail

/********************************************************************/
/*                                                                  */
/*                              GeoMap                              */
/*                                                                  */
/********************************************************************/

class GeoMap
{
  public:
    class Node;
    class Edge;
    class Face;
    class Dart;
    class SigmaAnchor;

    typedef std::vector< CELL_PTR(Node) > Nodes;
    typedef std::vector< CELL_PTR(Edge) > Edges;
    typedef std::vector< CELL_PTR(Face) > Faces;

    typedef vigra::SafeFilterIterator<Nodes::iterator, NotNull<Nodes::value_type> >
        NodeIterator;
    typedef vigra::SafeFilterIterator<Edges::iterator, NotNull<Edges::value_type> >
        EdgeIterator;
    typedef vigra::SafeFilterIterator<Faces::iterator, NotNull<Faces::value_type> >
        FaceIterator;
    typedef vigra::SafeFilterIterator<Nodes::const_iterator, NotNull<Nodes::value_type> >
        ConstNodeIterator;
    typedef vigra::SafeFilterIterator<Edges::const_iterator, NotNull<Edges::value_type> >
        ConstEdgeIterator;
    typedef vigra::SafeFilterIterator<Faces::const_iterator, NotNull<Faces::value_type> >
        ConstFaceIterator;

    typedef std::vector<int> SigmaMapping;
    typedef std::vector<double> EdgePreferences;

    typedef vigra::ConstImageIterator<int> LabelImageIterator;
    struct LabelImageAccessor {
        typedef int value_type;

        template<class Iterator>
        value_type operator()(Iterator it) const
        {
            int label = *it;
            if(label >= 0)
                label = faceLabelLUT_[label];
            return label;
        }

        template<class Iterator>
        value_type operator()(Iterator it,
                              typename Iterator::difference_type diff) const
        {
            int label = it[diff];
            if(label >= 0)
                label = faceLabelLUT_[label];
            return label;
        }

      protected:
        friend class GeoMap;

        LabelImageAccessor(LabelLUT const &faceLabelLUT)
        : faceLabelLUT_(faceLabelLUT)
        {}

        LabelLUT const &faceLabelLUT_;
    };

  protected:
    SigmaMapping
        sigmaMappingArray_,
        sigmaInverseMappingArray_;
    SigmaMapping::iterator
        sigmaMapping_,
        sigmaInverseMapping_;

    Nodes nodes_;
    Edges edges_;
    Faces faces_;

    unsigned int nodeCount_;
    unsigned int edgeCount_;
    unsigned int faceCount_;

    typedef vigra::PositionedObject<vigra::Vector2, CellLabel> PositionedNodeLabel;
    typedef vigra::Map2D<PositionedNodeLabel> NodeMap;
    NodeMap nodeMap_;

    typedef vigra::MultiArray<2, int> LabelImage;

    vigra::Size2D imageSize_;
    LabelImage   *labelImage_;
    LabelLUT      faceLabelLUT_;

    bool edgesSorted_;
    std::auto_ptr<detail::PlannedSplits> splitInfo_;
    std::auto_ptr<EdgePreferences> edgePreferences_;

  public:
    GeoMap(vigra::Size2D imageSize);
    GeoMap(const GeoMap &other);
    ~GeoMap();

    NodeIterator nodesBegin()
        { return NodeIterator(nodes_.begin(), nodes_.end()); }
    NodeIterator nodesEnd()
        { return NodeIterator(nodes_.end(), nodes_.end()); }
    CELL_PTR(Node) node(CellLabel label)
    {
        vigra_precondition(label < nodes_.size(), "invalid node label!");
        return nodes_[label];
    }
    ConstNodeIterator nodesBegin() const
        { return ConstNodeIterator(nodes_.begin(), nodes_.end()); }
    ConstNodeIterator nodesEnd() const
        { return ConstNodeIterator(nodes_.end(), nodes_.end()); }
    CELL_PTR(const Node) node(CellLabel label) const
    {
        vigra_precondition(label < nodes_.size(), "invalid node label!");
        return nodes_[label];
    }

    EdgeIterator edgesBegin()
        { return EdgeIterator(edges_.begin(), edges_.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edges_.end(), edges_.end()); }
    CELL_PTR(Edge) edge(CellLabel label)
    {
        vigra_precondition(label < edges_.size(), "invalid edge label!");
        return edges_[label];
    }
    ConstEdgeIterator edgesBegin() const
        { return ConstEdgeIterator(edges_.begin(), edges_.end()); }
    ConstEdgeIterator edgesEnd() const
        { return ConstEdgeIterator(edges_.end(), edges_.end()); }
    CELL_PTR(const Edge) edge(CellLabel label) const
    {
        vigra_precondition(label < edges_.size(), "invalid edge label!");
        return edges_[label];
    }

    FaceIterator facesBegin()
        { return FaceIterator(faces_.begin(), faces_.end()); }
    FaceIterator facesEnd()
        { return FaceIterator(faces_.end(), faces_.end()); }
    CELL_PTR(Face) face(CellLabel label)
    {
        vigra_precondition(label < faces_.size(), "invalid face label!");
        return faces_[label];
    }
    ConstFaceIterator facesBegin() const
        { return ConstFaceIterator(faces_.begin(), faces_.end()); }
    ConstFaceIterator facesEnd() const
        { return ConstFaceIterator(faces_.end(), faces_.end()); }
    CELL_PTR(const Face) face(CellLabel label) const
    {
        vigra_precondition(label < faces_.size(), "invalid face label!");
        return faces_[label];
    }

    inline Dart dart(int label);
    CELL_PTR(Face) faceAt(const vigra::Vector2 &position);

    CellLabel nodeCount() const { return nodeCount_; }
    CellLabel maxNodeLabel() const { return nodes_.size(); }
    CellLabel edgeCount() const { return edgeCount_; }
    CellLabel maxEdgeLabel() const { return edges_.size(); }
    CellLabel faceCount() const { return faceCount_; }
    CellLabel maxFaceLabel() const { return faces_.size(); }

    const vigra::Size2D &imageSize() const
    {
        return imageSize_;
    }

    CELL_PTR(Node) addNode(const vigra::Vector2 &position);
    CELL_PTR(Node) addNode(const vigra::Vector2 &position, CellLabel label);
    CELL_PTR(Edge) addEdge(const SigmaAnchor &startNeighbor,
                           const SigmaAnchor &endNeighbor,
                           const Vector2Array &points, CellLabel label = 0);
    // will always return NULL if !mapInitialized():
    CELL_PTR(Face) removeEdge(Dart &dart);

    void sortEdgesDirectly();

    typedef std::list<std::list<int> > UnsortableGroups;
    void sortEdgesEventually(double stepDist, double minDist,
                             UnsortableGroups &unsortable,
                             bool splitEdges);
    void setEdgePreferences(std::auto_ptr<EdgePreferences> edgePreferences)
    {
        vigra_precondition(
            splitInfo_.get(),
            "setting edge preferences futile - splitting impossible");
        edgePreferences_ = edgePreferences;
    }
    void splitParallelEdges();
        // for debugging / paper writing only:
    detail::PlannedSplits *internalSplitInfo()
        { return splitInfo_.get(); }

    void setSigmaMapping(SigmaMapping const &sigmaMapping, bool sorted = true);
    const SigmaMapping &sigmaMapping()
        { return sigmaMappingArray_; }

    bool edgesSorted() const { return edgesSorted_; }

    void initializeMap(bool initLabelImage = true);
    bool mapInitialized() const  { return faces_.size() > 0; }
    bool hasLabelImage() const { return labelImage_; }

    LabelImageIterator labelsUpperLeft() const
    {
        return LabelImageIterator(labelImage_->data(), labelImage_->shape(0));
    }
    LabelImageIterator labelsLowerRight() const
    {
        return labelsUpperLeft() + imageSize_;
    }
    LabelImageAccessor labelAccessor() const
    {
        return LabelImageAccessor(faceLabelLUT_);
    }

    vigra::triple<LabelImageIterator, LabelImageIterator, LabelImageAccessor>
    srcLabelRange() const
    {
        vigra_precondition(
            hasLabelImage(), "trying to access labelImage of GeoMap w/o label image!");
        return srcIterRange(labelsUpperLeft(),
                            labelsLowerRight(),
                            labelAccessor());
    }

        // maxFaceLabel is to be +1 here, too:
    void changeFaceLabels(const std::vector<CellLabel> &newFaceLabels,
                          CellLabel maxFaceLabel);

  protected:
    void initContours();
    void embedFaces(bool initLabelImage);
    void resizeSigmaMapping(SigmaMapping::size_type newSize);
    void insertSigmaPredecessor(int successor, int newPredecessor);
    void detachDart(int dartLabel);

  public:
    CELL_PTR(Node) nearestNode(
        const vigra::Vector2 &position,
        double maxSquaredDist = vigra::NumericTraits<double>::max());

    bool checkConsistency();

  protected:
    void associatePixels(Face &face, const PixelList &pixels);

  public:
    bool removeIsolatedNode(Node &node);
    CELL_PTR(Edge) mergeEdges(Dart &dart);
    CELL_PTR(Edge) splitEdge(Edge &edge, unsigned int segmentIndex);
    CELL_PTR(Edge) splitEdge(Edge &edge, unsigned int segmentIndex,
                             const vigra::Vector2 &newPoint,
                             bool insertPoint = true);
    CELL_PTR(Face) removeBridge(Dart &dart);
    CELL_PTR(Face) mergeFaces(Dart &dart);

        // callbacks using libsigc++ <http://libsigc.sourceforge.net/>:
    sigc::signal<bool, Node &>::accumulated<interruptable_accumulator>
        removeNodeHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preMergeEdgesHook;
    sigc::signal<void, Edge &>
        postMergeEdgesHook;
    sigc::signal<void, Edge &, unsigned int, vigra::Vector2 const &, bool>
        preSplitEdgeHook;
    sigc::signal<void, Edge &, Edge &>
        postSplitEdgeHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preRemoveBridgeHook;
    sigc::signal<void, Face &>
        postRemoveBridgeHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preMergeFacesHook;
    sigc::signal<void, Face &>
        postMergeFacesHook;
    sigc::signal<void, const Face &, const PixelList &>
        associatePixelsHook;
};

/********************************************************************/
/*                                                                  */
/*                           GeoMap::Node                           */
/*                                                                  */
/********************************************************************/

class GeoMap::Node : boost::noncopyable
{
  protected:
    GeoMap        *map_;
    CellLabel      label_;
    vigra::Vector2 position_;
    int            anchor_;

    friend class GeoMap; // give access to anchor_ (add edge, sort edges, Euler..)
    friend class SigmaAnchor; // give access to anchor_

    inline void uninitialize();

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

  public:
    bool initialized() const
    {
        return map_ != NULL;
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

    bool isIsolated() const
    {
        return !anchor_;
    }

    inline unsigned int degree() const;

    inline bool hasMinDegree(unsigned int minDegree) const;

    inline bool hasDegree(unsigned int exactDegree) const;

    bool operator==(const GeoMap::Node &other)
    {
        return label() == other.label();
    }

    bool operator!=(const GeoMap::Node &other)
    {
        return !operator==(other);
    }

    GeoMap *map() const
    {
        return map_;
    }
};

/********************************************************************/
/*                                                                  */
/*                           GeoMap::Edge                           */
/*                                                                  */
/********************************************************************/

const CellLabel UNINITIALIZED_CELL_LABEL =
    vigra::NumericTraits<CellLabel>::max();

class GeoMap::Edge
: public Polygon, boost::noncopyable
{
  public:
    typedef vigra::BBoxPolygon<vigra::Vector2> Base;

    enum { ALL_PROTECTION = 0xff };

  protected:
    GeoMap      *map_;
    CellLabel    label_;
    CellLabel    startNodeLabel_, endNodeLabel_;
    CellLabel    leftFaceLabel_, rightFaceLabel_;
    unsigned int flags_;

    mutable std::auto_ptr<vigra::Scanlines> scanLines_;

    friend class Dart; // allow setLeftFaceLabel
    friend class GeoMap;

    inline void uninitialize();

    void concatenate(Edge &other, bool atEnd, bool reverse);

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

  public:
    bool initialized() const
    {
        return map_ != NULL;
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
        return map_->face(leftFaceLabel());
    }

    CellLabel rightFaceLabel() const
    {
        return rightFaceLabel_;
    }

    GeoMap::Faces::value_type rightFace() const
    {
        vigra_precondition(initialized(), "rightFace() of uninitialized edge!");
        return map_->face(rightFaceLabel());
    }

    bool isBridge() const
    {
        return leftFaceLabel() == rightFaceLabel();
    }

    bool isLoop() const
    {
        return startNodeLabel_ == endNodeLabel_;
    }

    bool operator==(const GeoMap::Edge &other)
    {
        return label() == other.label();
    }

    bool operator!=(const GeoMap::Edge &other)
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

    CellLabel &internalLeftFaceLabel()
    {
        if(label_ > 0)
            return guaranteedEdge()->leftFaceLabel_;
        else
            return guaranteedEdge()->rightFaceLabel_;
    }

    friend class Face; // allow internalLeftFaceLabel in Face constructor
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

double contourArea(const GeoMap::Dart &dart);
double contourLength(const GeoMap::Dart &dart);
double isoperimetricQuotient(const GeoMap::Dart &dart);
Polygon contourPoly(const GeoMap::Dart &dart);

/********************************************************************/
/*                                                                  */
/*                           GeoMap::Face                           */
/*                                                                  */
/********************************************************************/

class GeoMap::Face : boost::noncopyable
{
  public:
    typedef Edge::BoundingBox BoundingBox;
    // Interestingly, std::vector seems to be faster than
    // std::list here, although we frequently need
    // erase/concatenation (cf. splice in cppmap.cxx):
    typedef std::vector<Dart> Contours;
    typedef Contours::const_iterator ContourIterator;

  protected:
    GeoMap              *map_;
    CellLabel            label_;
    Contours             anchors_;
    mutable unsigned int flags_;
    mutable BoundingBox  boundingBox_;
    mutable double       area_;
    unsigned int         pixelArea_;

    enum {
        BOUNDING_BOX_VALID = 0x80000000U,
        AREA_VALID         = 0x40000000U,
        INTERNAL_FLAGS     = 0xf0000000U,
    };

    friend class GeoMap; // give access to pixelArea_ and anchors_ (Euler ops...)

    inline void uninitialize();
    typedef Contours::iterator AnchorIterator; // non-const ContourIterator
    AnchorIterator findComponentAnchor(const GeoMap::Dart &dart);

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

            for(; anchor.internalLeftFaceLabel() == UNINITIALIZED_CELL_LABEL;
                anchor.nextPhi())
            {
                anchor.internalLeftFaceLabel() = label_;

                // don't calculate area on-the-fly here; we want to
                // exclude bridges from the area!
            }
        }
    }

  public:
    bool initialized() const
    {
        return map_ != NULL;
    }

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
            Dart anchor(*anchors_.begin()), dart(anchor);
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
            detail::IVector2 iPos(detail::intVPos(point));
            if(map_->labelImage_->isInside(iPos))
            {
                int l = (*map_->labelImage_)[iPos];
                if(l > 0)
                    return map_->faceLabelLUT_[l] == label_;
            }
        }
        ContourIterator it = contoursBegin();
        if(label_)
        {
            if(!boundingBox().contains(point))
                return false;
            if(!contourPoly(*it).contains(point))
                return false;
            ++it;
        }
        for(; it != contoursEnd(); ++it)
            if(contourPoly(*it).contains(point))
                return false;
        return true;
    }

    double area() const
    {
        if(!flag(AREA_VALID))
        {
            area_ = 0.0;
            for(ContourIterator it = contoursBegin();
                it != contoursEnd(); ++it)
            {
                area_ += contourArea(*it);
            }
            flags_ |= AREA_VALID;
        }
        return area_;
    }

    std::auto_ptr<vigra::Scanlines> scanLines() const
    {
        int startIndex = (int)floor(boundingBox().begin()[1] + 0.5);
        int endIndex = (int)floor(boundingBox_.end()[1] + 0.5);

        std::auto_ptr<vigra::Scanlines> result(
            new vigra::Scanlines(startIndex, endIndex - startIndex + 1));

        Dart anchor(*anchors_.begin()), dart(anchor);
        do
        {
            if(dart.label() > 0)
            {
                result->merge(dart.edge()->scanLines());
            }
            else
            {
                vigra::Scanlines sl(dart.edge()->scanLines());
                sl.reverse();
                result->merge(sl);
            }
        }
        while(dart.nextPhi() != anchor);

        result->normalize();

        return result;
    }

    unsigned int pixelArea() const
    {
        return pixelArea_;
    }

    const Dart &contour()
    {
        return *contoursBegin();
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

    unsigned int holeCount() const
    {
        return anchors_.size() - (label() ? 1 : 0);
    }

    void embedContour(const Dart &anchor);

    bool operator==(const GeoMap::Face &other)
    {
        return label() == other.label();
    }

    bool operator!=(const GeoMap::Face &other)
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
};

/********************************************************************/

inline GeoMap::Dart GeoMap::dart(int label)
{
    return GeoMap::Dart(this, label);
}

inline void GeoMap::Node::uninitialize()
{
    GeoMap *map = map_;
    map_ = NULL; // DON'T MESS WITH THIS!
    --map->nodeCount_;
    map->nodeMap_.erase(
        map->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                              vigra::NumericTraits<double>::epsilon()));
    RESET_PTR(map->nodes_[label_]); // may have effect like "delete this;"!
}

inline void GeoMap::Edge::uninitialize()
{
    GeoMap *map = map_;
    map_ = NULL;
    --map->edgeCount_;
    RESET_PTR(map->edges_[label_]);
}

inline void GeoMap::Face::uninitialize()
{
    GeoMap *map = map_;
    map_ = NULL;
    --map->faceCount_;
    RESET_PTR(map->faces_[label_]);
}

inline GeoMap::Dart GeoMap::Node::anchor() const
{
    vigra_precondition(initialized(), "anchor() of uninitialized node!");
    vigra_precondition(anchor_ != 0, "anchor() of degree 0 node!");
    return Dart(map_, anchor_);
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

inline bool GeoMap::Node::hasMinDegree(unsigned int minDegree) const
{
    if(!anchor_)
        return minDegree == 0;

    if(!minDegree)
        return true;

    GeoMap::Dart d(map_, anchor_);
    do
    {
        if(--minDegree == 0)
            return true;
    }
    while(d.nextSigma().label() != anchor_);

    return false;
}

inline bool GeoMap::Node::hasDegree(unsigned int exactDegree) const
{
    if(!anchor_)
        return exactDegree == 0;

    if(!exactDegree)
        return false;

    GeoMap::Dart d(map_, anchor_);
    do
    {
        if(d.nextSigma().label() == anchor_)
            return exactDegree == 1;
    }
    while(--exactDegree);
    return false;
}

inline GeoMap::Dart GeoMap::Edge::dart() const
{
    return map_->dart(label());
}

/********************************************************************/

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

/********************************************************************/

#include <iostream> // FIXME: not here, please!

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

/********************************************************************/

namespace detail {

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
        unsigned int splitGroup;

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

class PlannedSplits : public std::vector<DartPositionAngle::SplitPos>
{
  public:
    PlannedSplits()
    : splitGroupCount(0)
    {}

    unsigned int splitGroupCount;
};

} // namespace detail

#endif // VIGRA_CPPMAP_HXX
