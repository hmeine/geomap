#include "cppmap.hxx"
#include "exporthelpers.hxx"
#include <vigra/pythonimage.hxx>
#include <vigra/pythonutil.hxx>
#include <vigra/copyimage.hxx>
#include <iostream>

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

/********************************************************************/

class SimpleCallback
{
  public:
    virtual ~SimpleCallback()
    {
        disconnect();
    }

    bool connected() const
    {
        return connections_.size() && connections_[0].connected();
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

    void operator()(const GeoMap::Face &face, const PixelList &pixels)
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
            CELL_PTR(GeoMap::Node) startNode(result->node(startNodeLabel));
            CELL_PTR(GeoMap::Node) endNode(result->node(endNodeLabel));
            vigra_precondition(
                startNode && endNode, "invalid start- or endNodeLabel!");
            result->addEdge(*startNode, *endNode, pe(), i);
        }
    }

    return result;
}

bp::object
GeoMap__copy__(bp::object map)
{
    GeoMap *newMap(new GeoMap(bp::extract<const GeoMap &>(map)));
    bp::object result(bp::detail::new_reference(bp::managingPyObject(newMap)));

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        map.attr("__dict__"));

    return result;
}

bp::object
GeoMap__deepcopy__(bp::object map, bp::dict memo)
{
    bp::object copyMod = bp::import("copy");
    bp::object deepcopy = copyMod.attr("deepcopy");

    GeoMap *newMap(new GeoMap(bp::extract<const GeoMap &>(map)));
    bp::object result(bp::detail::new_reference(bp::managingPyObject(newMap)));

    // HACK: mapId shall be the same as the result of id(map) in Python -
    // please tell me that there is a better way! (and which ;-p)
    int mapId = (int)(map.ptr());
    memo[mapId] = result;

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        deepcopy(bp::extract<bp::dict>(map.attr("__dict__"))(), memo));

    return result;
}

CELL_PTR(GeoMap::Edge) addEdgeBackwardCompatibility(
    GeoMap &geomap,
    CellLabel startNodeLabel, CellLabel endNodeLabel,
    const Vector2Array &points, CellLabel label)
{
    std::cerr << "API warning: addEdge() now takes Node objects, not labels!\n";
    CELL_PTR(GeoMap::Node)
        startNode(geomap.node(startNodeLabel)),
        endNode(geomap.node(endNodeLabel));
    vigra_precondition(
        startNode && endNode, "invalid start- or endNodeLabel!");
    return geomap.addEdge(*startNode, *endNode, points, label);
}

GeoMap::FaceIterator faceIter(GeoMap &geomap, bool skipInfinite)
{
    vigra_precondition(geomap.mapInitialized(),
                       "faceIter() called on graph (!mapInitialized())");
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
    copyImage(map.srcLabelRange(), destImage(result));
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

bp::list GeoMap_sigmaMapping(GeoMap &map)
{
    bp::list result;
    const GeoMap::SigmaMapping &sigmaMapping(map.sigmaMapping());
    GeoMap::SigmaMapping::const_iterator
        first = sigmaMapping.begin(),
        last = sigmaMapping.end();
    while(*first == 0 && last[-1] == 0 && last >= first + 3)
    {
        ++first;
        --last;
    }
    for(GeoMap::SigmaMapping::const_iterator it = first;
        it != last; ++it)
    {
        if(map.edge(abs(*it)))
            result.append(*it);
        else
            result.append(0);
    }
    return result;
}

void GeoMap_setSigmaMapping(GeoMap &map, bp::list pySigmaMapping, bool edgesSorted = true)
{
    GeoMap::SigmaMapping sigmaMapping(bp::len(pySigmaMapping));

    unsigned int i = 0;
    for(GeoMap::SigmaMapping::iterator it = sigmaMapping.begin();
        it != sigmaMapping.end(); ++it, ++i)
    {
        *it = bp::extract<int>(pySigmaMapping[i])();
    }

    map.setSigmaMapping(sigmaMapping, edgesSorted);
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

    static bp::tuple getstate(bp::object pyMap)
    {
        GeoMap &map((bp::extract<GeoMap &>(pyMap)()));

        bp::list pySigmaMapping = GeoMap_sigmaMapping(map);

        bp::list edgeFlags;
        for(GeoMap::EdgeIterator it = map.edgesBegin(); it.inRange(); ++it)
        {
            edgeFlags.append((*it)->flags());
        }

        bp::list faceFlags, faceAnchors, faceLabels;
        for(GeoMap::FaceIterator it = map.facesBegin(); it.inRange(); ++it)
        {
            faceFlags.append((*it)->flags() & ~0xf0000000U);
            faceAnchors.append((*it)->contour().label());
            faceLabels.append((*it)->label());
        }

        return bp::make_tuple(
            pySigmaMapping,
            map.edgesSorted(),
            map.mapInitialized(),
            map.hasLabelImage(),
            edgeFlags,
            faceFlags, faceAnchors, make_tuple(faceLabels, map.maxFaceLabel()),
            pyMap.attr("__dict__"));
    }

    static bool getstate_manages_dict() { return true; }

    static void setstate(bp::object pyMap, bp::tuple state)
    {
        GeoMap &map((bp::extract<GeoMap &>(pyMap)()));

        bp::list pySigmaMapping((
            bp::extract<bp::list>(state[0])()));
        bool edgesSorted = bp::extract<bool>(state[1])();
        bool mapInitialized = bp::extract<bool>(state[2])();
        bool initLabelImage = bp::extract<bool>(state[3])();
        bp::list edgeFlags = bp::extract<bp::list>(state[4])();
        bp::list faceFlags = bp::extract<bp::list>(state[5])();
        bp::list faceAnchors = bp::extract<bp::list>(state[6])();
        bool hasLabels = len(state) > 8; // backward compat.
        bp::list faceLabels;
        CellLabel newMaxFaceLabel = 0;
        if(hasLabels)
        {
            faceLabels = bp::extract<bp::list>(state[7][0])();
            newMaxFaceLabel = bp::extract<CellLabel>(state[7][1])();
        }
        bp::object __dict__ = state[-1];

        GeoMap_setSigmaMapping(map, pySigmaMapping, edgesSorted);
        if(mapInitialized)
            map.initializeMap(initLabelImage);

        unsigned int i = 0;
        for(GeoMap::EdgeIterator it = map.edgesBegin(); it.inRange(); ++it, ++i)
        {
            (*it)->setFlag(bp::extract<unsigned int>(edgeFlags[i])());
        }

        std::vector<CellLabel> newFaceLabels(map.maxFaceLabel());
        for(i = 0; i < map.faceCount(); ++i)
        {
            GeoMap::Dart anchor =
                map.dart(bp::extract<unsigned int>(faceAnchors[i])());
            anchor.leftFace()->setFlag(
                bp::extract<unsigned int>(faceFlags[i])() & ~0xf0000000U);
            if(hasLabels)
                newFaceLabels[anchor.leftFaceLabel()] =
                    bp::extract<CellLabel>(faceLabels[i])();
        }

        if(hasLabels)
            map.changeFaceLabels(newFaceLabels, newMaxFaceLabel);

        bp::extract<bp::dict>(pyMap.attr("__dict__"))().update(__dict__);
    }
};

void GeoMap_setEdgePreferences(GeoMap &geomap,
                               bp::list edgePreferences)
{
    std::auto_ptr<GeoMap::EdgePreferences> cppep(
        new GeoMap::EdgePreferences(len(edgePreferences)));

    for(unsigned int i = 0; i < cppep->size(); ++i)
    {
        double pref = bp::extract<double>(edgePreferences[i])();
        (*cppep)[i] = pref;
    }

    geomap.setEdgePreferences(cppep);
}

bp::object GeoMap_internalSplitInfo(GeoMap &geoMap)
{
    const detail::PlannedSplits *splitInfo = geoMap.internalSplitInfo();
    if(!splitInfo)
        return bp::object();
    bp::list result;
    for(unsigned int i = 0; i < splitInfo->size(); ++i)
    {
        const detail::PlannedSplits::value_type &si((*splitInfo)[i]);
        result.append(bp::make_tuple(
                          si.segmentIndex, si.arcLength, si.position,
                          si.dartLabel, si.sigmaPos, si.splitGroup));
    }
    return result;
}

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
      //<< ", node " << edge.startNodeLabel() << " -> " << edge.endNodeLabel()
      << ", faces " << edge.leftFaceLabel() << "(l), " << edge.rightFaceLabel()
      << ", partial area " << edge.partialArea() << ", length " << edge.length()
      << ", " << edge.size() << " points>";
    return s.str();
}

std::string Face__repr__(GeoMap::Face const &face)
{
    std::stringstream s;
    s << "<GeoMap.Face " << face.label() << ", "
      << face.holeCount() << " holes, area " << face.area();
    if(face.area() != face.pixelArea()) // prevent redundant output for crack-edge maps
        s << " (" << face.pixelArea() << " px)";
    s << ">";
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

std::string SimpleCallback__repr__(SimpleCallback const &cb)
{
    std::stringstream s;
    s << "<SimpleCallback, ";
    if(cb.connected())
        s << "active>";
    else
        s << "detached>";
    return s.str();
}

template<class T>
T returnCopy(const T &v)
{
    return v;
}

/********************************************************************/

void defMapStats();
void defMapUtils();

void defMap()
{
    using namespace boost::python;

    CELL_RETURN_POLICY crp;

    {
        scope geoMap(
            class_<GeoMap, boost::noncopyable>(
                "GeoMap",
                "The GeoMap class manages a topologically consistent set of nodes,\n"
                "edges, and faces.\n\n"
                "You can get the number of nodes via the `nodeCount` property, and\n"
                "iterate over all existing nodes with ::\n\n"
                "  for node in amap.nodeIter():\n"
                "    ...\n\n"
                "The same goes for edges and faces.\n\n"
                "You can monitor changes on the map by adding your own callbacks via:\n\n"
                "addMergeEdgesCallbacks(preOpCallback, postOpCallback)\n"
                "  called before/after merging two edges with a node of degree two\n"
                "addMergeFacesCallbacks(preOpCallback, postOpCallback)\n"
                "  called before/after merging two faces via a common edge\n"
                "addRemoveBridgeCallbacks(preOpCallback, postOpCallback)\n"
                "  called before/after removing a bridge (\"loose\" edge within a face)\n"
                "addRemoveNodeCallback(callback)\n"
                "  called before removing a node (callbacks are passed the Node object\n"
                "  being removed)\n"
                "addAssociatePixelsCallback(callback)\n"
                "  called when the removal of edges results in pixels being\n"
                "  associated with the surrounding face (called with the face as\n"
                "  first, and a list of Point2Ds as second parameter)\n\n"
                "GeoMap objects can be pickled, at will guarantee persistence of edges\n"
                "and nodes with their geometry and labels (face labels may change!).",
                init<vigra::Size2D>(arg("imageSize"),
                    "GeoMap(nodePositions, edges, imageSize)\n\n"
                    "Creates a GeoMap of the given `imageSize`.\n\n"
                    "`nodePositions`\n"
                    "  is a list of node positions (``None`` values are ignored).\n"
                    "`edges`\n"
                    "  is a list of edge triples (startNodeLabel, endNodeLabel, points)\n"
                    "  (index 0 is ignored, so are ``None`` values inside the list).\n"
                    "`imageSize`\n"
                    "  is the size of the labelImage to be created\n"
                    "  (usually, the original image size).\n\n"
                    "If `nodePositions` and `edges` are given (they can be omitted),\n"
                    "adds these to the map as if `addNode`/`addEdge` were called.\n"
                    "``None`` is allowed inside these two lists,\n"
                    "and all created Node/Edge objects will have labels that\n"
                    "point to their corresponding entries in the lists (that is, ``None``\n"
                    "will lead to gaps in the labels).\n"
                    "Since edges must have non-zero labels, ``edges[0]`` *must* be ``None``."))
            .def("__init__", make_constructor(
                     &createGeoMap, default_call_policies(),
                     (arg("nodePositions") = list(),
                      arg("edgeTuples") = list(),
                      arg("imageSize") = vigra::Size2D(0, 0))))
            .def("__copy__", &GeoMap__copy__)
            .def("__deepcopy__", &GeoMap__deepcopy__)
            .def("node", (CELL_PTR(GeoMap::Node)(GeoMap::*)(CellLabel))&GeoMap::node, crp,
                 "node(label) -> Node\n\n"
                 "Return Node object for the given label.")
            .def("nodeIter", (GeoMap::NodeIterator(GeoMap::*)())&GeoMap::nodesBegin,
                 "Iterates over all existing nodes.\n\n"
                 ">>> for node in amap.nodeIter():\n"
                 "...     print node.label(), node.degree(), node.anchor()")
            .def("edge", (CELL_PTR(GeoMap::Edge)(GeoMap::*)(CellLabel))&GeoMap::edge, crp,
                 "edge(label) -> Edge\n\n"
                 "Return Edge object for the given label.")
            .def("edgeIter", (GeoMap::EdgeIterator(GeoMap::*)())&GeoMap::edgesBegin,
                 "Iterates over all existing edges.\n\n"
                 ">>> for edge in amap.edgeIter():\n"
                 "...     print \"Edge %d has %d points\" % len(edge)")
            .def("face", (CELL_PTR(GeoMap::Face)(GeoMap::*)(CellLabel))&GeoMap::face, crp,
                 "face(label) -> Face\n\n"
                 "Return Face object for the given label.")
            .def("faceIter", &faceIter, arg("skipInfinite") = false,
                 "Iterates over all existing faces.\n\n"
                 ">>> for face in amap.faceIter():\n"
                 "...     for dart in face.contours():\n"
                 "...         print dart\n\n"
                 "For your convenience, it supports skipping the infinite face:\n\n"
                 ">>> assert amap.faceIter(skipInfinite = True).next().label() > 0\n")
            .def("dart", &GeoMap::dart,
                 "dart(label) -> Dart\n\n"
                 "Return `Dart` object for the given `label`.\n"
                 "The dart will point to the `edge` labelled ``abs(label)``, and\n"
                 "negative dart labels correspond to the opposite half-edge\n"
                 "(meaning that Darts with negative labels start at the end of\n"
                 "the corresponding edge).")
            .def("faceAt", &GeoMap::faceAt, crp)
            .add_property("nodeCount", &GeoMap::nodeCount,
                          "Return the number of nodes in this graph/map.")
            .add_property("edgeCount", &GeoMap::edgeCount,
                          "Return the number of edges in this graph/map.")
            .add_property("faceCount", &GeoMap::faceCount,
                          "Return the number of faces in this graph/map.")
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
                 return_value_policy<copy_const_reference>(),
                 "imageSize() -> Size2D\n\n"
                 "Return the image size (as passed to __init__).\n"
                 "If the map has a labelImage, this is its size.")
            .def("addNode", (CELL_PTR(GeoMap::Node) (GeoMap::*)(const vigra::Vector2 &))&GeoMap::addNode, crp, args("position"),
                 "addNode(position) -> Node\n\n"
                 "Add node at the given position and return the new Node\n"
                 "object.")
            .def("addNode", (CELL_PTR(GeoMap::Node) (GeoMap::*)(const vigra::Vector2 &, CellLabel))&GeoMap::addNode, crp, args("position", "label"))
            .def("addEdge", &addEdgeBackwardCompatibility, crp,
                 (arg("startNodeLabel"), arg("endNodeLabel"),
                  arg("points"), arg("label") = 0),
                 "addEdge(startNeighbor, endNeighbor, points, label = None) -> Edge\n\n"
                 "Add edge new edge with given geometry to the map and return\n"
                 "the new `Edge` object.  If start/endNeighbor are Node objects,\n"
                 "the edge is inserted at an arbitrary point in the sigma orbit,\n"
                 "if Dart objects are given, the new edge becomes the sigma\n"
                 "predecessor of the given start and/or edge darts.  (For\n"
                 "isolated nodes, you have to pass a Node object.)\n\n"
                 "A label >= maxEdgeLabel() may be given in order to\n"
                 "force the given label (the edges with lower labels will then\n"
                 "stay unused ATM).")
            .def("addEdge", &GeoMap::addEdge, crp,
                 (arg("startNeighbor"), arg("endNeighbor"),
                  arg("points"), arg("label") = 0))
            .def("removeEdge", &GeoMap::removeEdge, crp,
                 "removeEdge(dart)\n\n"
                 "Equivalent to either `removeBridge` or `mergeFaces`, depending on\n"
                 "the type of edge.  May even be called if not `mapInitialized()`\n"
                 "in order to remove an edge from the graph.")
            .def("splitEdge", &pySplitEdge,
                 (arg("edge"), arg("segmentIndex"), arg("newPoint") = object()),
                 "splitEdge(self, edge, segmentIndex, newPoint = None) -> Edge\n\n"
                 "Splits the given `edge` by cutting it into two at the given position.\n"
                 "`newPoint` should be a point in the polygon segment with the given index.\n"
                 "If None or not given, the edge will be split exactly at the beginning\n"
                 "of the segment given by `segmentIndex`.\n\n"
                 "splitEdge() returns the new resulting edge whose `Edge.startNode()` is a new\n"
                 "node that has become the new `Edge.endNode()` of the `edge` that was split.")

            .def("removeIsolatedNode", &GeoMap::removeIsolatedNode, arg("dart"),
                 "removeIsolatedNode(node)\n\n"
                 "Euler operation removing an isolated node.  In contrast to the\n"
                 "other operations, this one cannot be parametrized with a dart,\n"
                 "since by definition no dart is attached to an isolated node.\n\n"
                 "ATM, this operation is called internally (from the other Euler\n"
                 "operations) whenever a node becomes isolated, and does not have to\n"
                 "be called manually.  (This is the reason why it does not appear in\n"
                 "map.history.)")
            .def("mergeEdges", &GeoMap::mergeEdges, arg("dart"), crp,
                 "mergeEdges(dart)\n\n"
                 "Euler operation merging two edges around a node of degree 2.\n"
                 "Return the surviving (merged) edge (which might either be\n"
                 "dart.edge() or dart.nextSigma().edge()).  The resulting edge\n"
                 "has the same direction as before.")
            .def("removeBridge", &GeoMap::removeBridge, arg("dart"), crp)
            .def("mergeFaces", &GeoMap::mergeFaces, arg("dart"), crp)
            .def("sortEdgesDirectly", &GeoMap::sortEdgesDirectly,
                 "sortEdgesDirectly()\n\n"
                 "Sort edges around nodes by taking into account the direction\n"
                 "of the first polyline segment attached to a node.  This simple\n"
                 "method is not suitable for subpixel watersheds, see\n"
                 "sortEdgesEventually instead.")
            .def("sortEdgesEventually", &GeoMap_sortEdgesEventually,
                 (arg("stepDist"), arg("minDist"), arg("splitEdges") = true),
                 "sortEdgesEventually(stepDist, minDist, splitEdges = True)\n\n"
                 "Sort edges by following them until the eventually part.\n"
                 "minDist is the minimal distance the flowlines must be apart,\n"
                 "stepDist is the (arclength) step size the sigma sorting\n"
                 "algorithm will use to walk along groups of parallel edges.\n\n"
                 "If splitEdges is True, the positions of divergence are stored\n"
                 "for later calling of splitEdges() - note that the edges are\n"
                 "not split if you do not call the latter method.")
            .def("sigmaMapping", &GeoMap_sigmaMapping,
                 "Returns a list containing the sigma permutation.\n"
                 "The list is a mapping from dart labels to the labels of their\n"
                 "sigma successors, where the mapping must be seen relative to\n"
                 "the middle of the list (its length is always odd).")
            .def("setSigmaMapping", &GeoMap_setSigmaMapping,
                 (arg("sigmaMapping"), arg("sorted") = true),
                 "setSigmaMapping(sigmaMapping, sorted = True)\n"
                 "Initialize the sigma permutation.\n"
                 "The argument is expected to be a mapping from dart labels to\n"
                 "the labels of their sigma successors, see `sigmaMapping()`.\n"
                 "Usually, the GeoMap will thus become oriented,\n"
                 "i.e. `edgesSorted()` will be True afterwards, as if\n"
                 "e.g. `sortEdgesDirectly()` was called.\n"
                 "The optional parameter `sorted` can be used to specify whether\n"
                 "the given sigmaMapping shall make the GeoMap sorted.")
            .def("edgesSorted", &GeoMap::edgesSorted,
                 "edgesSorted() -> bool\n\n"
                 "Return whether the edges have already been sorted.\n"
                 "Normally this is the case when one of\n"
                 "`sortEdgesDirectly`/`sortEdgesEventually`/`setSigmaMapping`\n"
                 "has been used.")
            .def("splitParallelEdges", &GeoMap::splitParallelEdges)
            .def("setEdgePreferences", &GeoMap_setEdgePreferences,
                 "setEdgePreferences(list)\n\n"
                 "Set edge preferences (one float per edge).  This is used for\n"
                 "deciding upon the survivor when merging parallel darts in\n"
                 "`splitParallelEdges()`; edges with higher values are more likely\n"
                 "to be preserved.  (If no preferences are given, the dart with\n"
                 "the smallest curvature around the split node survives.)")
            .def("_internalSplitInfo", &GeoMap_internalSplitInfo,
                 "for debugging / paper writing only\n"
                 "list of (segmentIndex, arcLength, position, dartLabel, sigmaPos, splitGroup)")
            .def("initializeMap", &GeoMap::initializeMap, (arg("initLabelImage") = true),
                 "initializeMap(initLabelImage = True) -> None\n\n"
                 "This finishes the initialization of a GeoMap.  Call this after\n"
                 "setting up the geometry (adding nodes/edges via\n"
                 "`addNode`/`addEdge`) and initializing the sigma orbits\n"
                 "(e.g. `sortEdgesDirectly`/`sortEdgesEventually`/`setSigmaMapping`).\n"
                 "Given the embedded graph, this will initialize the faces and\n"
                 "their embedding.\n\n"
                 "Subsequently, `mapInitialized()` will return True and\n"
                 "`initializeMap()` must not be called again.\n\n"
                 "If initLabelImage is False, `labelImage()` will not\n"
                 "return a label image, and `hasLabelImage()` will return\n"
                 "False.  This is useful e.g. for GeoMaps with Delaunay edges,\n"
                 "where most pixels facets in the label image are crossed by an\n"
                 "edge (and thus don't contain face labels) anyways.")
            .def("mapInitialized", &GeoMap::mapInitialized,
                 "mapInitialized() -> bool\n\n"
                 "Return whether initializeMap() has already been\n"
                 "called.  Otherwise, someMap is considered to be a graph.  See\n"
                 "also edgesSorted().")
            .def("hasLabelImage", &GeoMap::hasLabelImage,
                 "hasLabelImage() -> bool\n\n"
                 "Return True when someMap has a labelImage,\n"
                 "i.e. initializeMap(...) has already been called with\n"
                 "its initLabelImage parameter set to True (default).  Then,\n"
                 "labelImage() returns a GrayImage, else None.")
            .def("labelImage", &labelImage,
                 "labelImage() -> GrayImage/None\n\n"
                 "Return a GrayImage where all pixels that are entirely inside\n"
                 "a region are assigned that regions' label.  All pixels whose\n"
                 "facet is crossed by an edge are assigned negative numbers\n"
                 "indicating how many edges intersect this pixel.\n\n"
                 "If hasLabelImage() is False, this method may return\n"
                 "None.")
            .def("nearestNode", &GeoMap::nearestNode, crp,
                 (arg("position"), arg(
                      "maxSquaredDist") = vigra::NumericTraits<double>::max()),
                 "nearestNode(position[, maxSquaredDist]) -> Node\n\n"
                 "Return the nearest node to the given position.  If\n"
                 "`maxSquaredDist` dist is given and no Node is within range, ``None``\n"
                 "is returned instead.")
            .def("checkConsistency", &GeoMap::checkConsistency,
                 "checkConsistency() -> bool\n\n"
                 "Performs a series of consistency/sanity checks and returns\n"
                 "True if someMap seems to be OK.")
            .def_pickle(GeoMapPickleSuite())
            );

        RangeIterWrapper<GeoMap::NodeIterator, CELL_RETURN_POLICY>("_NodeIterator");
        RangeIterWrapper<GeoMap::EdgeIterator, CELL_RETURN_POLICY>("_EdgeIterator");
        RangeIterWrapper<GeoMap::FaceIterator, CELL_RETURN_POLICY>("_FaceIterator");

        class_<GeoMap::Node, boost::noncopyable>(
            "Node",
            "Represents a node of the GeoMap.  You can get an `anchor`\n"
            "Dart that is attached to the Node (if not `isIsolated()`), and\n"
            "you can query its `position()` and `degree()`.",
            no_init)
            .def("initialized", &GeoMap::Node::initialized)
            .def("label", &GeoMap::Node::label,
                 "label() -> int\n\n"
                 "Return the label of this Node.  Node labels are >= 0 and\n"
                 "you can query the maximal value by `GeoMap.maxNodeLabel()`.")
            .def("position", &GeoMap::Node::position,
                 return_value_policy<copy_const_reference>())
            .def("setPosition", &GeoMap::Node::setPosition)
            .def("degree", &GeoMap::Node::degree)
            .def("hasMinDegree", &GeoMap::Node::hasMinDegree, arg("minDegree"))
            .def("hasDegree", &GeoMap::Node::hasDegree, arg("exactDegree"))
            .def("isIsolated", &GeoMap::Node::isIsolated,
                 "isIsolated() -> bool\n\n"
                 "Return True iff there is no edge attached to this node.\n"
                 "You must not call `anchor` on such Darts.")
            .def("anchor", &GeoMap::Node::anchor)
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Node__repr__)
        ;

        class_<GeoMap::Edge, bases<Polygon>, boost::noncopyable>(
            "Edge",
            "Represents an edge within the GeoMap.  This class is\n"
            "derived from `Polygon`, so a GeoMap.Edge *is* a `Polygon` and\n"
            "you can iterate it, ``len(edge)`` gives the number of points,\n"
            "``edge[5]`` returns its 6th support point etc.\n\n"
            "Furthermore, it gives access to the left/right face, the\n"
            "start/end node and you can easily create a `dart` which points\n"
            "to this edge.\n\n"
            "Edges (and Faces) allow you to use `setFlag` in order to tag\n"
            "edges with bit-flags.  `flags` returns all flags or'ed\n"
            "together and `flag` can be used to simply test for a flag.\n"
            "There is a module ``flag_constants`` that is used to note\n"
            "which flags are already used by some applications/framework\n"
            "parts.",
            no_init)
            .def("initialized", &GeoMap::Edge::initialized)
            .def("label", &GeoMap::Edge::label,
                 "label() -> int\n\n"
                 "Return the label of this Edge.  Edge labels are >= 1 and\n"
                 "you can query the maximal value by `GeoMap.maxEdgeLabel()`.")
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
            .def("flags", &GeoMap::Edge::flags)
            .def("flag", &GeoMap::Edge::flag, arg("which"))
            .def("setFlag", &GeoMap::Edge::setFlag,
                 (arg("flag"), arg("onoff") = true))
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Edge__repr__)
        ;
        scope().attr("Edge").attr("ALL_PROTECTION") =
            (unsigned int)GeoMap::Edge::ALL_PROTECTION;

        class_<GeoMap::Face, boost::noncopyable>(
            "Face",
            "Represents a face within the GeoMap.  This is used to store\n"
            "anchors for the `contours` (i.e. one for the outer `contour`,\n"
            "and one for each of the `holeContours`).  Furthermore, it\n"
            "caches information like the `boundingBox` and/or the `area` of\n"
            "the whole face.\n\n"
            "Faces (and Edges) allow you to use `setFlag` in order to tag\n"
            "edges with bit-flags.  `flags` returns all flags or'ed\n"
            "together and `flag` can be used to simply test for a flag.\n"
            "There is a module ``flag_constants`` that is used to note\n"
            "which flags are already used by some applications/framework\n"
            "parts.",
            no_init)
            .def("initialized", &GeoMap::Face::initialized)
            .def("label", &GeoMap::Face::label,
                 "label() -> int\n\n"
                 "Return the label of this Face.  Face labels are >= 0 and\n"
                 "you can query the maximal value by `GeoMap.maxFaceLabel()`.\n"
                 "Face 0 is always the surrounding infinite face.")
            .def("boundingBox", &GeoMap::Face::boundingBox,
                 return_value_policy<copy_const_reference>(),
                 "boundingBox() -> BoundingBox\n\n"
                 "Return the (axis-parallel) bounding box of this face.")
            .def("contains", &GeoMap::Face::contains,
                 "contains(point) -> bool\n\n"
                 "Return whether this face contains the given point, i.e.\n"
                 "whether it is within the exterior contour, but not within\n"
                 "any hole contour.")
            .def("area", &GeoMap::Face::area,
                 "area() -> float\n\n"
                 "Return area of this region (without holes, i.e.\n"
                 "area(exterior)-sum_i(area(interior_i)) ).  The return\n"
                 "value for the map's infinite exterior face (label 0) is\n"
                 "not defined.")
            .def("pixelArea", &GeoMap::Face::pixelArea)
            .def("contour", &GeoMap::Face::contour,
                 return_value_policy<copy_const_reference>())
            .def("contours", &faceContours)
            .def("holeContours", &faceHoleContours)
            .def("holeCount", &GeoMap::Face::holeCount)
            .def("scanLines", &GeoMap::Face::scanLines)
            .def("flags", &GeoMap::Face::flags)
            .def("flag", &GeoMap::Face::flag, arg("which"))
            .def("setFlag", &GeoMap::Face::setFlag,
                 (arg("flag"), arg("onoff") = true))
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Face__repr__)
        ;

        RangeIterWrapper<ContourRangeIterator>("_ContourIterator");

        return_internal_reference<> rself; // "return self" policy

        RangeIterWrapper<DartPointIter>("_DartPointIter");

        class_<GeoMap::Dart>("Dart",
                             "Points to a dart within the GeoMap.  Dart objects are used\n"
                             "like iterators, for traversing a GeoMap and analyzing its\n"
                             "topology and/or geometry.  A Dart behaves like a `Polygon`,\n"
                             "i.e. like a sequence of points (you can iterate it,\n"
                             "``len(dart)`` gives the number of points, ``dart[5]`` returns\n"
                             "the 6th support point of the directed polygon, etc.)",
                             no_init)
            .def(init<GeoMap *, int>())
            .def("clone", &GeoMap::Dart::clone,
                 "Return a copy of this Dart.  This is especially useful\n"
                 "if you want to call nextFoo() without changing the\n"
                 "original object.")
            .def("label", &GeoMap::Dart::label,
                 "label() -> int\n\n"
                 "Return the label of this Dart.")
            .def("map", &GeoMap::Dart::map,
                 return_value_policy<reference_existing_object>(),
                 "map() -> GeoMap\n\n"
                 "Return the GeoMap this Dart belongs to.")
            .def("edgeLabel", &GeoMap::Dart::edgeLabel,
                 "edgeLabel() -> int\n\n"
                 "Return the label of the edge.  This is equivalent to\n"
                 "abs(dart.label()).")
            .def("edge", &GeoMap::Dart::edge, crp,
                "edge() -> Edge/None\n\n"
                 "Returns corresponding edge or None if that edge has\n"
                 "already been removed.")
            .def("guaranteedEdge", &GeoMap::Dart::guaranteedEdge, crp,
                "guaranteedEdge() -> Edge\n\n"
                 "Equivalent to `edge()`, but throws an exception if the\n"
                 "edge is not valid anymore (i.e. has been removed via an\n"
                 "Euler operation.")
            .def("startNodeLabel", &GeoMap::Dart::startNodeLabel,
                 "startNodeLabel() -> int\n\n"
                 "Return the label of the node that is at the start of this dart.")
            .def("startNode", &GeoMap::Dart::startNode, crp,
                 "startNode() -> Node\n\n"
                 "Return the node that is at the start of this dart.")
            .def("endNodeLabel", &GeoMap::Dart::endNodeLabel,
                 "endNodeLabel() -> int\n\n"
                 "Return the label of the node that is at the end of this dart.")
            .def("endNode", &GeoMap::Dart::endNode, crp,
                 "endNode() -> Node\n\n"
                 "Return the node that is at the end of this dart.\n"
                 "(Actually, this can be seen as a convenience method for\n"
                 "dart.clone().nextAlpha().startNode().")
            .def("leftFaceLabel", &GeoMap::Dart::leftFaceLabel,
                 "leftFaceLabel() -> int\n\n"
                 "Return the label of the face to the left of this dart.")
            .def("leftFace", &GeoMap::Dart::leftFace, crp,
                 "leftFace() -> Face\n\n"
                 "Return the face to the left of this dart.")
            .def("rightFaceLabel", &GeoMap::Dart::rightFaceLabel,
                 "rightFaceLabel() -> int\n\n"
                 "Return the label of the face to the right of this dart.")
            .def("rightFace", &GeoMap::Dart::rightFace, crp,
                 "rightFace() -> Face\n\n"
                 "Return the face to the right of this dart.")
            .def("partialArea", &GeoMap::Dart::partialArea, crp,
                 "partialArea() -> float\n\n"
                 "Return the (signed) partial area contributed by this dart.\n"
                 "This is the enclosed area if the dart is a self-loop.")
            .def("__getitem__", &Array__getitem__<GeoMap::Dart>)
            .def("__iter__", &GeoMap::Dart::pointIter)
            .def("__len__", &GeoMap::Dart::size)
            .def("nextAlpha", &GeoMap::Dart::nextAlpha, rself,
                 "Jump to the opposite dart.  (The opposite dart is the\n"
                 "successor in the alpha permutation.)  Since alpha is\n"
                 "an involution, two calls to nextAlpha() cancel each other.\n"
                 "Return the modified dart, so that you can combine multiple\n"
                 "operations like dart.nextAlpha().prevSigma()")
            .def("nextSigma", &GeoMap::Dart::nextSigma, rself,
                 "Turn counter-clockwise around the `startNode`.\n"
                 "(The cyclic order is encoded by the sigma permutation.)\n"
                 "The modified dart is returned, so that you can combine\n"
                 "multiple operations like dart.nextSigma().nextAlpha()")
            .def("prevSigma", &GeoMap::Dart::prevSigma, rself,
                 "Turn clockwise around the `startNode`.\n"
                 "(The cyclic order is encoded by the sigma permutation.)\n"
                 "The modified dart is returned, so that you can combine\n"
                 "multiple operations like dart.prevSigma().nextAlpha()")
            .def("nextPhi", &GeoMap::Dart::nextPhi, rself,
                 "Follow this contour of the face to the left.\n"
                 "(The phi permutation is composed by alpha and sigma^-1.)\n"
                 "This is a shortcut for `nextAlpha()`.`prevSigma()`.\n"
                 "The modified dart is returned, so that you can combine\n"
                 "multiple operations like dart.nextPhi().nextSigma()")
            .def("prevPhi", &GeoMap::Dart::prevPhi, rself,
                 "Follow this contour of the face to the left (backwards).\n"
                 "(The phi permutation is composed by alpha and sigma^-1.)\n"
                 "This is a shortcut for `nextSigma()`.`nextAlpha()`.\n"
                 "The modified dart is returned, so that you can combine\n"
                 "multiple operations like dart.prevPhi().nextSigma()")
            .def("phiOrbit", &phiOrbit,
                 "phiOrbit() -> iterator\n\n"
                 "Return an iterator over the darts in the phi orbit.\n"
                 "The phi orbit contains all darts in this contour of the\n"
                 "leftFace, starting with dart itself.")
            .def("sigmaOrbit", &sigmaOrbit,
                 "sigmaOrbit() -> iterator\n\n"
                 "Return an iterator over the darts in the sigma orbit.\n"
                 "i.e. a counter-clockwise ordering of all darts attached to\n"
                 "dart's startNode, starting with dart itself.")
            .def("alphaOrbit", &alphaOrbit,
                 "alphaOrbit() -> iterator\n\n"
                 "Return an iterator over the darts in the alpha orbit.\n"
                 "(Since alpha is an involution, there are always exactly\n"
                 "two darts in this orbit).")
            .def(self == self)
            .def(self != self)
            .def("__repr__", &Dart__repr__)
        ;

        RangeIterWrapper<PhiOrbitIterator>("_PhiOrbitIterator");
        register_ptr_to_python< std::auto_ptr<PhiOrbitIterator> >();
        RangeIterWrapper<SigmaOrbitIterator>("_SigmaOrbitIterator");
        register_ptr_to_python< std::auto_ptr<SigmaOrbitIterator> >();
        RangeIterWrapper<AlphaOrbitIterator>("_AlphaOrbitIterator");
        register_ptr_to_python< std::auto_ptr<AlphaOrbitIterator> >();

#ifndef USE_INSECURE_CELL_PTRS
        register_ptr_to_python< CELL_PTR(GeoMap::Node) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Edge) >();
        register_ptr_to_python< CELL_PTR(GeoMap::Face) >();
#endif

        class_<SimpleCallback>("SimpleCallback",
            "Internal class to manage Euler operation callbacks.\n"
            "Removes the callbacks from the caller on __del__ or `disconnect()`.",
                               no_init)
            .def("disconnect", &SimpleCallback::disconnect,
                 "Disconnect this callback.  This operation is not\n"
                 "reversible, you have to reconnect in the same way as before to\n"
                 "aquire a new `SimpleCallback`.")
            .def("connected", &SimpleCallback::connected)
            .def("__repr__", &SimpleCallback__repr__)
        ;

        def("addRemoveNodeCallback", &addRemoveNodeCallback,
            "addRemoveNodeCallback(callback) -> SimpleCallback\n\n"
            "Add a callback to be called called before removing an\n"
            "isolated node.  The callback will be called with the node to\n"
            "be removed as parameter.  If any callback does not return True,\n"
            "the operation will be canceled.");
        def("addMergeEdgesCallbacks", &addMergeEdgesCallbacks,
            "addMergeEdgesCallbacks(preOpCallback, postOpCallback)\n\n"
            "Add callbacks to be called called before/after merging two\n"
            "edges with a node of degree two.  The dart passed to\n"
            "preOpCallback belongs to the surviving edge and starts at the\n"
            "merged node, the postOpCallback gets the resulting Edge as\n"
            "parameter.  If any preOpCallback does not return True,\n"
            "the operation will be canceled.");
        def("addSplitEdgeCallbacks", &addSplitEdgeCallbacks,
            "addSplitEdgeCallbacks(preOpCallback, postOpCallback) -> SimpleCallback\n\n"
            "Add callbacks to be called called before/after splitting an\n"
            "edge.  The preOpCallback is called with the three arguments\n"
            "(edge, segmentIndex, newPoint) given to splitEdge, and the\n"
            "postOpCallback is given the resulting two edges (edge,\n"
            "newEdge).  A splitEdge() operation cannot be canceled through\n"
            "preOpCallbacks (their return values are ignored).");
        def("addRemoveBridgeCallbacks", &addRemoveBridgeCallbacks,
            "addRemoveBridgeCallbacks(preOpCallback, postOpCallback) -> SimpleCallback\n\n"
            "Add callbacks to be called called before/after removing a\n"
            "bridge (\"loose\" edge within a face).  The preOpCallback is\n"
            "called with the dart to be removed, the postOpCallback gets\n"
            "the surrounding Face the Edge was merged into as parameter.\n"
            "If any preOpCallback does not return True, the operation will\n"
            "be canceled.");
        def("addMergeFacesCallbacks", &addMergeFacesCallbacks,
            "addMergeFacesCallbacks(preOpCallback, postOpCallback) -> SimpleCallback\n\n"
            "Add callbacks to be called called before/after merging two\n"
            "faces via a common edge.  The rightFace() of the dart passed\n"
            "to preOpCallback will be merged into its leftFace(), the\n"
            "postOpCallback gets the resulting Face as parameter.  If any\n"
            "preOpCallback does not return True, the operation will be\n"
            "canceled.");
        def("addAssociatePixelsCallback", &addAssociatePixelsCallback,
            "addAssociatePixelsCallback(callback) -> SimpleCallback\n\n"
            "Add a callback to be called whenever the removal of edges\n"
            "results in pixels being associated with the surrounding face\n"
            "(i.e. after a mergeFaces or removeBridge operation).  The\n"
            "callback will be called with the face as first, and a list of\n"
            "Point2Ds as second parameter.  All internal GeoMap structures\n"
            "are updated before the callbacks happen.");
        register_ptr_to_python< std::auto_ptr<SimpleCallback> >();

        geoMap.attr("BYTES_PER_NODE") = sizeof(GeoMap::Node);
        geoMap.attr("BYTES_PER_EDGE") = sizeof(GeoMap::Edge);
        geoMap.attr("BYTES_PER_FACE") = sizeof(GeoMap::Face);
        geoMap.attr("BYTES_PER_MAP") = sizeof(GeoMap);
    }

    RangeIterWrapper<ContourPointIter>("ContourPointIter")
        .def(init<GeoMap::Dart, bool>((arg("dart"), arg("firstTwice") = false)));

    def("contourArea", &contourArea,
        "contourArea(anchor) -> float\n\n"
        "Returns the area of contourPoly(anchor) (is however much faster than\n"
        "using that function, since it simply sums up all partialArea()s of the\n"
        "darts in the phi orbit.");
    def("contourLength", &contourLength,
        "contourLength(anchor) -> float\n\n"
        "Returns the length of contourPoly(anchor) (is however much faster than\n"
        "using that function, since it simply sums up all length()s of the\n"
        "darts in the phi orbit.");
    def("isoperimetricQuotient", &isoperimetricQuotient,
        "isoperimetricQuotient(anchor) -> length\n\n"
        "Returns the isoperimetric quotient for contourPoly(anchor).\n"
        "This is defined by sq(contourLength(anchor))/(4*pi*contourArea(anchor)).");
    def("contourPoly", &contourPoly,
        "contourPoly(anchor) -> Polygon\n\n"
        "Returns a Polygon composed by traversing anchor's phi orbit once.");

    class_<DartPosition>("DartPosition",
                         "Helper class for traversing a Dart's geometry.",
                         init<GeoMap::Dart>())
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

    implicitly_convertible<GeoMap::Node, GeoMap::SigmaAnchor>();
    implicitly_convertible<GeoMap::Dart, GeoMap::SigmaAnchor>();

    defMapStats();
    defMapUtils();
}
