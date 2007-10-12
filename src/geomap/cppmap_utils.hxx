#ifndef CPPMAP_UTILS_HXX
#define CPPMAP_UTILS_HXX

#include "cppmap.hxx"
#include <boost/utility.hpp>

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

CELL_PTR(GeoMap::Face) mergeFacesCompletely(
    GeoMap::Dart &dart, bool mergeDegree2Nodes)
{
    vigra_precondition(!dart.edge()->isBridge(),
                       "mergeFacesCompletely(): dart belongs to a bridge!");
    GeoMap *map = dart.map();
    CellLabel rightLabel = dart.rightFaceLabel();

    std::vector<int> commonDarts;

    GeoMap::Dart d(dart);
    do
    {
        if(d.rightFaceLabel() == rightLabel)
        {
            if(d.edge()->flag(GeoMap::Edge::ALL_PROTECTION))
                return NULL_PTR(GeoMap::Face);
            commonDarts.push_back(d.label());
        }
    }
    while(d.nextPhi() != dart);

    std::vector<CellLabel> affectedNodes;

    CELL_PTR(GeoMap::Face) survivor = NULL_PTR(GeoMap::Face);
    for(std::vector<int>::iterator it = commonDarts.begin();
        it != commonDarts.end(); ++it)
    {
        d = map->dart(*it);
        affectedNodes.push_back(d.startNodeLabel());
        affectedNodes.push_back(d.endNodeLabel());
        if(!survivor)
            survivor = map->mergeFaces(d); // first common edge
        else
            map->removeBridge(d);
    }

    for(std::vector<CellLabel>::iterator it = affectedNodes.begin();
        it != affectedNodes.end(); ++it)
    {
        CELL_PTR(GeoMap::Node) node = map->node(*it);
        if(!node)
            continue;
        if(node->isIsolated())
            map->removeIsolatedNode(*node);
        if(mergeDegree2Nodes && node->hasDegree(2))
        {
            d = node->anchor();
            if(d.endNodeLabel() != node->label())
                map->mergeEdges(d);
        }
    }

    return survivor;
}

#endif // CPPMAP_UTILS_HXX
