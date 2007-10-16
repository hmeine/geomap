#ifndef CPPMAP_UTILS_HXX
#define CPPMAP_UTILS_HXX

#include "cppmap.hxx"
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

/********************************************************************/

/// Removes all isolated nodes with map.removeIsolatedNode(...) and
/// returns the number of successful operations (= nodes removed).
unsigned int removeIsolatedNodes(GeoMap &map)
{
    unsigned int result = 0;
    for(GeoMap::NodeIterator it = map.nodesBegin(); it.inRange(); ++it)
        if((*it)->isIsolated() && map.removeIsolatedNode(**it))
            ++result;
    return result;
}

/// Removes all degree-2-nodes with map.mergeEdges(...) and
/// returns the number of successful operations (= nodes removed).
unsigned int mergeDegree2Nodes(GeoMap &map)
{
    unsigned int result = 0;
    for(GeoMap::NodeIterator it = map.nodesBegin(); it.inRange(); ++it)
        if((*it)->hasDegree(2))
        {
            GeoMap::Dart dart((*it)->anchor());
            if(!(*it)->anchor().edge()->isLoop() &&
               map.mergeEdges(dart))
                ++result;
        }
    return result;
}

unsigned int removeBridges(GeoMap &map, std::list<CellLabel> &bridges)
{    
    unsigned int result = 0;

    typedef std::list<CellLabel> Bridges;

    while(bridges.size())
    {
        bool degree1Found = false;

        // search for bridge with endnode of degree 1:
        for(Bridges::iterator it = bridges.begin(); it != bridges.end(); ++it)
        {
            GeoMap::Dart dart(map.dart((int)*it));
            if(!dart.edge())
                continue;
            // degree 1 endnodes?
            if(dart.clone().nextAlpha().nextSigma().label() == -dart.label())
                dart.nextAlpha();
            else if(dart.clone().nextSigma() != dart)
                continue;

            degree1Found = true;
            while(true)
            {
                GeoMap::Dart next = dart.clone().nextPhi();
                if(map.removeBridge(dart))
                    ++result;

                if(next.edgeLabel() == dart.edgeLabel())
                    break;
                // degree != 1?
                if(next.clone().nextSigma() != next)
                    break;

                dart = next;
                if(!dart.edge()->flag(GeoMap::Edge::REMOVE_BRIDGE))
                    break;
            }
        }

        Bridges::iterator next = bridges.begin();
        do
        {
            Bridges::iterator it = next;
            ++next;
            if(!map.edge(*it))
                bridges.erase(it);
        }
        while(next != bridges.end());

        if(!degree1Found)
        {
            if(map.removeBridge(map.dart(bridges.front())))
                ++result;
            bridges.erase(bridges.begin());
        }
    }

    return result;
}

unsigned int removeBridges(GeoMap &map)
{
    std::list<CellLabel> bridgeLabels;

    for(GeoMap::EdgeIterator it = map.edgesBegin(); it.inRange(); ++it)
        if((*it)->isBridge())
        {
            (*it)->setFlag(GeoMap::Edge::REMOVE_BRIDGE);
            bridgeLabels.push_back((*it)->label());
        }

    return removeBridges(map, bridgeLabels);
}

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
        CELL_PTR(GeoMap::Edge)
            edge = map.edge(*edgeLabelsBegin);
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

#endif // CPPMAP_UTILS_HXX
