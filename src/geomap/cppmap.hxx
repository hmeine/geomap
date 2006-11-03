#ifndef CPPMAP_HXX
#define CPPMAP_HXX

#include "filteriterator.hxx"
#include "vigra/positionedmap.hxx"
#include "vigra/polygon.hxx"
#include <vector>
#include <boost/python.hpp> // FIXME: separate this from plain C++ interface
#include <vigra/multi_array.hxx>
namespace bp = boost::python;

//#define USE_INSECURE_CELL_PTRS

// The define USE_INSECURE_CELL_PTRS can be used to switch between
// "safe" cell handling e.g. for Python and a possibly faster C++ way.
// Actually, it's not 100% safe anyways ATM, since the cell
// destructors do access the map, so they must not be deleted after
// the Map. ;-(

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

typedef std::vector<vigra::Point2D> PixelList;

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

    class ModificationCallback;

    typedef std::vector< CELL_PTR(Node) > Nodes;
    typedef std::vector< CELL_PTR(Edge) > Edges;
    typedef std::vector< CELL_PTR(Face) > Faces;

    typedef vigra::FilterIterator<Nodes::iterator, NotNull<Nodes::value_type> >
        NodeIterator;
    typedef vigra::FilterIterator<Edges::iterator, NotNull<Edges::value_type> >
        EdgeIterator;
    typedef vigra::FilterIterator<Faces::iterator, NotNull<Faces::value_type> >
        FaceIterator;

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

    typedef std::vector<ModificationCallback *> ModificationCallbacks;
    typedef ModificationCallbacks::iterator MCIterator;

    ModificationCallbacks removeNodeHooks_;
    ModificationCallbacks mergeEdgesHooks_;
    ModificationCallbacks removeBridgeHooks_;
    ModificationCallbacks mergeFacesHooks_;
    ModificationCallbacks associatedPixelsHooks_;

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
    void sortEdgesDirectly();
    void sortEdgesEventually(double stepDist, double minDist);
    void initContours();
    void embedFaces(bool initLabelImage = true);

    CELL_PTR(Node) nearestNode(
        const vigra::Vector2 &position,
        double maxSquaredDist = vigra::NumericTraits<double>::max());

    bool checkConsistency();

  protected:
    void associatePixels(Face &face, const PixelList &pixels);

  public:
    void removeIsolatedNode(Node &node);
    Edge &mergeEdges(Dart &dart);
    Face &removeBridge(Dart &dart);
    Face &mergeFaces(Dart &dart);

  private:
    GeoMap(const GeoMap &) {} // disallow copying
    GeoMap &operator=(const GeoMap &) { return *this; }
};

#endif // CPPMAP_HXX
