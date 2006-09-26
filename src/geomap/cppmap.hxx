#ifndef CPPMAP_HXX
#define CPPMAP_HXX

#include "filteriterator.hxx"
#include "vigra/positionedmap.hxx"
#include <vector>
#include <boost/python.hpp> // FIXME: separate this from plain C++ interface
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
#  define CELL_RETURN_POLICY return_value_policy<reference_existing_object>
#else
#  include <boost/shared_ptr.hpp>
#  define CELL_PTR(Type) boost::shared_ptr<Type>
#  define NULL_PTR(Type) boost::shared_ptr<Type>()
#  define RESET_PTR(ptr) ptr.reset()
#  define CELL_RETURN_POLICY default_call_policies
#endif

typedef unsigned int CellLabel;

// functor for FilterIterator to skip NULL cells
template<class POINTER>
struct NotNull
{
    bool operator()(const POINTER &p) const
    {
        return p.get();
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
    class Node;
    class Edge;
    class Face;
    class Dart;

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

    PositionedMap nodeMap_;

    vigra::Size2D imageSize_;

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
        return nodes_[label];
    }

    EdgeIterator edgesBegin()
        { return EdgeIterator(edges_.begin(), edges_.end()); }
    EdgeIterator edgesEnd()
        { return EdgeIterator(edges_.end(), edges_.end()); }
    CELL_PTR(Edge) edge(CellLabel label)
    {
        return edges_[label];
    }

    FaceIterator facesBegin()
        { return FaceIterator(faces_.begin(), faces_.end()); }
    FaceIterator facesEnd()
        { return FaceIterator(faces_.end(), faces_.end()); }
    CELL_PTR(Face) face(CellLabel label)
    {
        return faces_[label];
    }

    inline Dart dart(int label);

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

    void sortEdgesDirectly();
};

#endif // CPPMAP_HXX
