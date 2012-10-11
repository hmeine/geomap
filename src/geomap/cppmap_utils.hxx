#ifndef CPPMAP_UTILS_HXX
#define CPPMAP_UTILS_HXX

#include "cppmap.hxx"
#include "polygon.hxx"
#include <boost/utility.hpp>
#include <vector>
#include <list>

class EdgeProtection : boost::noncopyable
{
  protected:
    boost::shared_ptr<GeoMap> map_;
    std::vector<sigc::connection> connections_;

  public:
    EdgeProtection(boost::shared_ptr<GeoMap> map = boost::shared_ptr<GeoMap>())
    {
        if(map)
            attachHooks(map);
    }

    const boost::shared_ptr<GeoMap> map() const
    {
        return map_;
    }

    void attachHooks(boost::shared_ptr<GeoMap> map)
    {
        vigra_precondition(!map_, "trying to attach to more than once?!");
        connections_.push_back(
            map->preMergeFacesHook.connect(
                sigc::mem_fun(this, &EdgeProtection::preRemoveEdge)));
        connections_.push_back(
            map->preRemoveBridgeHook.connect(
                sigc::mem_fun(this, &EdgeProtection::preRemoveEdge)));
        connections_.push_back(
            map->preMergeEdgesHook.connect(
                sigc::mem_fun(this, &EdgeProtection::preMergeEdges)));
        map_ = map;
    }

    void detachHooks()
    {
        for(unsigned int i = 0; i < connections_.size(); ++i)
            if(connections_[i].connected())
                connections_[i].disconnect();
        connections_.clear();
        map_.reset();
    }

        /// do not allow removal of protected edges
    bool preRemoveEdge(const GeoMap::Dart &dart)
    {
        return !dart.edge()->flag(GeoMap::Edge::ALL_PROTECTION);
    }

        /// only allow edge merging if the edges carry the same flags
    bool preMergeEdges(const GeoMap::Dart &dart)
    {
        return (dart.edge()->flags() ==
                dart.clone().nextSigma().edge()->flags());
    }
};

/********************************************************************/

GeoMap::FacePtr mergeFacesCompletely(GeoMap::Dart &dart, bool mergeDegree2Nodes);

/// Removes all isolated nodes with map.removeIsolatedNode(...) and
/// returns the number of successful operations (= nodes removed).
unsigned int removeIsolatedNodes(GeoMap &map);

/// Removes all degree-2-nodes with map.mergeEdges(...) and
/// returns the number of successful operations (= nodes removed).
unsigned int mergeDegree2Nodes(GeoMap &map);

unsigned int removeBridges(GeoMap &map, std::list<CellLabel> &bridges);

unsigned int removeBridges(GeoMap &map);

/// Removes all edges whose labels are in `edgeLabels`.
/// Uses an optimized sequence of basic Euler operations.
template<class ITERATOR>
unsigned int removeEdges(
    GeoMap &map, ITERATOR edgeLabelsBegin, ITERATOR edgeLabelsEnd)
{
    unsigned int result = 0;

    typedef std::list<CellLabel> Bridges;
    Bridges bridges;

    for(; edgeLabelsBegin != edgeLabelsEnd; ++edgeLabelsBegin)
    {
        GeoMap::EdgePtr
            edge = map.edge(*edgeLabelsBegin);

        vigra_precondition(edge, "removeEdges: illegal edge label");

        if(edge->isBridge())
        {
            edge->setFlag(GeoMap::Edge::REMOVE_BRIDGE);
            bridges.push_back(edge->label());
        }
        else if(map.mergeFaces(edge->dart()))
            ++result;
    }

    result += removeBridges(map, bridges);
    result += removeIsolatedNodes(map); // FIXME: depend on allowIsolatedNodes
    result += mergeDegree2Nodes(map);
    return result;
}

template<class DestIterator, class SizeType, class DestAccessor>
void drawLabelImage(const GeoMap &geomap,
                    DestIterator dul, SizeType ds, DestAccessor a,
                    bool negativeEdgeLabels = true)
{
    for(GeoMap::ConstFaceIterator
            it = geomap.finiteFacesBegin(); it.inRange(); ++it)
    {
        std::auto_ptr<vigra::Scanlines> scanlines =
            (*it)->scanLines();
        fillScannedPoly(*scanlines, (int)(*it)->label(),
                        dul, ds, a);
        if(negativeEdgeLabels)
            vigra::drawScannedPoly(*(*it)->scanLines(), -1,
                                   dul, ds, a);
    }

    if(!negativeEdgeLabels)
    {
        vigra::Rect2D imageROI(vigra::Size2D(ds[0], ds[1]));

        for(GeoMap::ConstEdgeIterator
                it = geomap.edgesBegin(); it.inRange(); ++it)
        {
            for(vigra::ScanlinesIter sit((*it)->scanLines());
                sit.inRange(); ++sit)
            {
                vigra::Point2D p(*sit);
                if(!imageROI.contains(p))
                    continue;

                // TODO: unconditionally write negative labels above
                // and only call faceAt for these?
                a.set(geomap.faceAt(Vector2(p.x, p.y))->label(),
                      dul, SizeType(p.x, p.y));
            }
        }
    }
}

template<class DestIterator, class SizeType, class DestAccessor>
void drawLabelImage(const GeoMap &geomap,
                    vigra::triple<DestIterator, SizeType, DestAccessor> d,
                    bool negativeEdgeLabels = true)
{
    return drawLabelImage(geomap, d.first, d.second, d.third, negativeEdgeLabels);
}

#endif // CPPMAP_UTILS_HXX
