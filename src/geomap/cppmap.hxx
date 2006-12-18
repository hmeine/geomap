#ifndef CPPMAP_HXX
#define CPPMAP_HXX

#include "filteriterator.hxx"
#include "vigra/positionedmap.hxx"
#include "vigra/polygon.hxx"
#include <vector>
#include <boost/python.hpp> // FIXME: separate this from plain C++ interface
#include <vigra/multi_array.hxx>
#include <sigc++/sigc++.h>
namespace bp = boost::python;

//#define USE_INSECURE_CELL_PTRS

// The define USE_INSECURE_CELL_PTRS can be used to switch between
// "safe" cell handling e.g. for Python and a possibly faster C++ way.

#ifdef USE_INSECURE_CELL_PTRS
#  define CELL_PTR(Type) Type *
#  define NULL_PTR(Type) (Type *)NULL
#  define RESET_PTR(ptr) ptr = NULL
// This is quite dangerous, *but*: The real lifetime of
// the referenced objects / cells are unknown, since any Euler
// operation might invalidate them.
#  define CELL_RETURN_POLICY bp::return_value_policy<reference_existing_object>
#else
#  include <boost/shared_ptr.hpp>
#  define CELL_PTR(Type) boost::shared_ptr<Type>
#  define NULL_PTR(Type) boost::shared_ptr<Type>()
#  define RESET_PTR(ptr) ptr.reset()
#  define CELL_RETURN_POLICY bp::default_call_policies
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

/********************************************************************/
/*                                                                  */
/*                              GeoMap                              */
/*                                                                  */
/********************************************************************/

class PlannedSplits;

class GeoMap
{
  public:
    class Node;
    class Edge;
    class Face;
    class Dart;

    typedef std::vector< CELL_PTR(Node) > Nodes;
    typedef std::vector< CELL_PTR(Edge) > Edges;
    typedef std::vector< CELL_PTR(Face) > Faces;

    typedef vigra::SafeFilterIterator<Nodes::iterator, NotNull<Nodes::value_type> >
        NodeIterator;
    typedef vigra::SafeFilterIterator<Edges::iterator, NotNull<Edges::value_type> >
        EdgeIterator;
    typedef vigra::SafeFilterIterator<Faces::iterator, NotNull<Faces::value_type> >
        FaceIterator;

    typedef vigra::ConstImageIterator<int> LabelImageIterator;
    struct LabelImageAccessor {
        typedef int value_type;

        template<class Iterator>
        value_type operator()(Iterator it)
        {
            int label = *it;
            if(label >= 0)
                label = faceLabelLUT_[label];
            return label;
        }

        template<class Iterator>
        value_type operator()(Iterator it,
                              typename Iterator::difference_type diff)
        {
            int label = it[diff];
            if(label >= 0)
                label = faceLabelLUT_[label];
            return label;
        }

      protected:
        friend class GeoMap;

        LabelImageAccessor(std::vector<CellLabel> const &faceLabelLUT)
        : faceLabelLUT_(faceLabelLUT)
        {}

        std::vector<CellLabel> const &faceLabelLUT_;
    };

  protected:
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

    vigra::Size2D          imageSize_;
    LabelImage            *labelImage_;
    std::vector<CellLabel> faceLabelLUT_;

    bool edgesSorted_;
    std::auto_ptr<PlannedSplits> splitInfo_;

  public:
    GeoMap(bp::list nodePositions,
           bp::list edgeTuples, vigra::Size2D imageSize);

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

    EdgeIterator edgesBegin()
        { return EdgeIterator(edges_.begin(), edges_.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edges_.end(), edges_.end()); }
    CELL_PTR(Edge) edge(CellLabel label)
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
    CELL_PTR(Edge) addEdge(CellLabel startNodeLabel, CellLabel endNodeLabel,
                           const Vector2Array &points, CellLabel label = 0);
    CELL_PTR(Edge) addEdge(Dart startNeighbor, Dart endNeighbor,
                           const Vector2Array &points);
    void removeEdge(Dart &dart);

    void sortEdgesDirectly();

    typedef std::list<std::list<int> > UnsortableGroups;
    void sortEdgesEventually(double stepDist, double minDist,
                             UnsortableGroups &unsortable,
                             bool splitEdges);
    void splitParallelEdges();

    typedef std::vector<int> SigmaMapping;
    void setSigmaMapping(SigmaMapping const &sigmaMapping, bool sorted = true);
    std::auto_ptr<SigmaMapping > sigmaMapping();

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

  protected:
    void initContours();
    void embedFaces(bool initLabelImage);

  public:
    CELL_PTR(Node) nearestNode(
        const vigra::Vector2 &position,
        double maxSquaredDist = vigra::NumericTraits<double>::max());

    bool checkConsistency();

  protected:
    void associatePixels(Face &face, const PixelList &pixels);

  public:
    bool removeIsolatedNode(Node &node);
    Edge &mergeEdges(Dart &dart);
    Edge &splitEdge(Edge &edge, unsigned int segmentIndex);
    Edge &splitEdge(Edge &edge, unsigned int segmentIndex,
                    const vigra::Vector2 &newPoint, bool insertPoint = true);
    Face &removeBridge(Dart &dart);
    Face &mergeFaces(Dart &dart);

        // callbacks using libsigc++ <http://libsigc.sourceforge.net/>:
    sigc::signal<bool, Node &>::accumulated<interruptable_accumulator>
        removeNodeHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preMergeEdgesHook;
    sigc::signal<void, Edge &>
        postMergeEdgesHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preRemoveBridgeHook;
    sigc::signal<void, Face &>
        postRemoveBridgeHook;
    sigc::signal<bool, const Dart &>::accumulated<interruptable_accumulator>
        preMergeFacesHook;
    sigc::signal<void, Face &>
        postMergeFacesHook;
    sigc::signal<void, Face &, const PixelList &>
        associatePixelsHook;

  private:
    GeoMap(const GeoMap &) {} // disallow copying
    GeoMap &operator=(const GeoMap &) { return *this; }
};

#endif // CPPMAP_HXX
