#include "cppmap_utils.hxx"

GeoMap::FacePtr mergeFacesCompletely(
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

    GeoMap::FacePtr survivor = NULL_PTR(GeoMap::Face);
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
        GeoMap::NodePtr node = map->node(*it);
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

unsigned int removeIsolatedNodes(GeoMap &map)
{
    unsigned int result = 0;
    for(GeoMap::NodeIterator it = map.nodesBegin(); it.inRange(); ++it)
        if((*it)->isIsolated() && map.removeIsolatedNode(**it))
            ++result;
    return result;
}

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

    unsigned int const UNREMOVABLE_BRIDGE = 0x80000;

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
                else
                    dart.edge()->setFlag(UNREMOVABLE_BRIDGE);

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
            GeoMap::EdgePtr bridge = map.edge(*it);
            if(!bridge)
                bridges.erase(it);
            else if(bridge->flag(UNREMOVABLE_BRIDGE))
            {
                bridge->setFlag(UNREMOVABLE_BRIDGE, false); // clean up
                bridges.erase(it);
            }
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
