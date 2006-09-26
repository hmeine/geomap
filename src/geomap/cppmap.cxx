#include "cppmap.hxx"
#include <boost/python/detail/api_placeholder.hpp>
#include <vigra/tinyvector.hxx>
#include <vigra/pythonimage.hxx>
#include "vigra/polygon.hxx"
#include <iostream>
#include "exporthelpers.hxx"

typedef vigra::TinyVector<double, 2> Vector2;
typedef vigra::PointArray<Vector2>   Vector2Array;
typedef vigra::BBoxPolygon<vigra::Vector2> Polygon;

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
        map_->nodeMap_.insert(position_, bp::object(label_));
    }

    ~Node()
    {
        RESET_PTR(map_->nodes_[label_]);
        --map_->nodeCount_;
        map_->nodeMap_.remove(position_);
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
      leftFaceLabel_(0),
      rightFaceLabel_(0),
      protection_(0)
    {
        map_->edges_.push_back(GeoMap::Edges::value_type(this));
        ++map_->edgeCount_;
    }

    ~Edge()
    {
        RESET_PTR(map_->edges_[label_]);
        --map_->edgeCount_;
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
        return map_->node(startNodeLabel_);
    }

    CellLabel endNodeLabel() const
    {
        return endNodeLabel_;
    }

    GeoMap::Nodes::value_type endNode() const
    {
        return map_->node(endNodeLabel_);
    }

    CellLabel leftFaceLabel() const
    {
        return leftFaceLabel_;
    }

    GeoMap::Faces::value_type leftFace() const
    {
        return map_->face(leftFaceLabel_);
    }

    CellLabel rightFaceLabel() const
    {
        return rightFaceLabel_;
    }

    GeoMap::Faces::value_type rightFace() const
    {
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

//         def _setLeftFaceLabel(self, label):
//             if self.label() > 0:
//                 self.edge()._leftFaceLabel = label
//             else:
//                 self.edge()._rightFaceLabel = label

//         def edge(self):
//             """Returns corresponding edge or None if that edge has already
//             been removed."""
//             result = self._map().edge(self.edgeLabel())
//             return result

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

//         def __iter__(self):
//             if self._label > 0:
//                 return self.guaranteedEdge().__iter__()
//             else:
//                 return self.guaranteedEdge().__reviter__()

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

            for(; !anchor.leftFaceLabel(); anchor.nextPhi())
            {
                // don't calculate area on-the-fly here; we want to
                // exclude bridges from the area!
                anchor.setLeftFaceLabel(label_);
            }
        }
    }

    ~Face()
    {
        RESET_PTR(map_->faces_[label_]);
        --map_->faceCount_;
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
        if(!boundingBox().contains(point))
            return false;
        unsigned int i = 0;
        if(label_)
        {
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
        for(; !dart.leftFaceLabel(); dart.nextPhi())
            dart.setLeftFaceLabel(label_);

        if(areaValid_)
            area_ += contourArea(dart);

        vigra_postcondition(dart == anchor,
                            "contour labeled partially?!");
    }
};

void GeoMap::Node::setPosition(const vigra::Vector2 &p)
{
    map_->nodeMap_.remove(position_);
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
    map_->nodeMap_.insert(p, bp::object(label_));
}

inline GeoMap::Dart GeoMap::Node::anchor() const
{
    return Dart(map_, darts_[0]);
}

GeoMap::GeoMap(bp::list nodePositions,
               bp::list edgeTuples, vigra::Size2D imageSize)
: nodeCount_(0),
  edgeCount_(0),
  faceCount_(0),
  imageSize_(imageSize)
{
    for(int i = 0; i < len(nodePositions); ++i)
    {
        bp::extract<Vector2> ve(nodePositions[i]);
        if(ve.check())
            new Node(this, ve());
        else
            nodes_.push_back(NULL_PTR(Node));
    }

    for(int i = 0; i < len(edgeTuples); ++i)
    {
        bp::extract<bp::tuple> ete(edgeTuples[i]);
        if(ete.check())
        {
            bp::tuple edgeTuple(ete());
            bp::extract<Vector2Array> pe(edgeTuple[2]);
            if(!pe.check())
            {
                std::cerr << "why, oh why, do I have to die??\n";
                PyErr_SetString(PyExc_TypeError,
                    "GeoMap.__init__: edge geometry not convertable to Vector2Array");
                bp::throw_error_already_set();
            }
            CellLabel startNodeLabel = bp::extract<CellLabel>(edgeTuple[0])();
            CellLabel endNodeLabel   = bp::extract<CellLabel>(edgeTuple[1])();
            Edge *edge = new Edge(this, startNodeLabel, endNodeLabel, pe());
            node(startNodeLabel)->darts_.push_back(edge->label());
            node(endNodeLabel)->darts_.push_back(-edge->label());
        }
        else
            edges_.push_back(NULL_PTR(Edge));
    }

    sortEdgesDirectly();
}

GeoMap::~GeoMap()
{
    // make sure these objects are deleted first
    // (their destructors access this map!)
    nodes_.clear();
    edges_.clear();
    faces_.clear();
}

inline GeoMap::Dart GeoMap::dart(int label)
{
    return GeoMap::Dart(this, label);
}

double angleTheta(double dy, double dx); // implemented in polygon.cxx

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

template<class T>
T &returnSelf(T &v)
{
    return v;
}

template<class Iterator>
typename Iterator::value_type nextIterPos(Iterator &v)
{
    if(!v.inRange())
    {
        PyErr_SetString(PyExc_StopIteration, "cells iterator exhausted");
        bp::throw_error_already_set();
    }
    return *v++;
}

using namespace boost::python;

void defMap()
{
    {
        CELL_RETURN_POLICY crp;

        scope geoMap(
            class_<GeoMap>("GeoMap", init<
                           bp::list, bp::list, vigra::Size2D>())
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
            .def("sortEdgesDirectly", &GeoMap::sortEdgesDirectly)
            );

        // return_internal_reference is a simplification, since the
        // true owner would be the GeoMap object (however, the
        // iterator's lifetime is expected to be long enough):
        class_<GeoMap::NodeIterator>("NodeIterator")
            .def("__iter__", (GeoMap::NodeIterator &
                              (*)(GeoMap::NodeIterator &))&returnSelf,
                 return_internal_reference<>())
            .def("next", &nextIterPos<GeoMap::NodeIterator>, crp);

        class_<GeoMap::EdgeIterator>("EdgeIterator")
            .def("__iter__", (GeoMap::EdgeIterator &
                              (*)(GeoMap::EdgeIterator &))&returnSelf,
                 return_internal_reference<>())
            .def("next", &nextIterPos<GeoMap::EdgeIterator>, crp);

        class_<GeoMap::FaceIterator>("FaceIterator")
            .def("__iter__", (GeoMap::FaceIterator &
                              (*)(GeoMap::FaceIterator &))&returnSelf,
                 return_internal_reference<>())
            .def("next", &nextIterPos<GeoMap::FaceIterator>, crp);

        class_<GeoMap::Node>("Node", init<GeoMap *, const vigra::Vector2 &>())
            .def("label", &GeoMap::Node::label)
            .def("position", &GeoMap::Node::position,
                 return_value_policy<copy_const_reference>())
            .def("setPosition", &GeoMap::Node::setPosition)
            .def("degree", &GeoMap::Node::degree)
            .def("anchor", &GeoMap::Node::anchor)
        ;

        class_<GeoMap::Edge, bases<Polygon> >("Edge", no_init)
            .def(init<GeoMap *, CellLabel, CellLabel, GeoMap::Edge::Base>())
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
        ;

        class_<GeoMap::Face>("Face", init<GeoMap *, GeoMap::Dart>())
            .def("label", &GeoMap::Face::label)
            .def("boundingBox", &GeoMap::Face::boundingBox,
                 return_value_policy<copy_const_reference>())
            .def("area", &GeoMap::Face::area)
        ;

        return_internal_reference<> rself; // "return self" policy

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

    def("contourArea", &contourArea,
        "contourArea(anchor) -> float\n\n"
        "Returns the area of contourPoly(anchor) (is however much faster than\n"
        "using that function, since it simply sums up all partialArea()s of the\n"
        "darts in the phi orbit.");
    def("contourPoly", &contourPoly,
        "contourPoly(anchor) -> Polygon\n\n"
        "Returns a Polygon composed by traversing anchor's phi orbit once.");
}

