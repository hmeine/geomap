import hourglass

BORDER_PROTECTION = 1

# --------------------------------------------------------------------
#                             simple proxies
# --------------------------------------------------------------------

def removeIsolatedNode(node):
    return node.anchor().map().removeIsolatedNode(node)

def mergeEdges(dart):
    return dart.map().mergeEdges(dart)

def removeBridge(dart):
    return dart.map().removeBridge(dart)

def mergeFaces(dart):
    return dart.map().mergeFaces(dart)

# --------------------------------------------------------------------
#                       composed Euler operations
# --------------------------------------------------------------------

def removeEdge(dart):
    """removeEdge(dart)

    Composed Euler operation which calls either removeBridge or
    mergeFaces, depending on whether the dart belongs to an edge or a
    bridge.

    Returns the surviving face."""

    map = dart.map()
    if dart.leftFaceLabel() == dart.rightFaceLabel():
        return map.removeBridge(dart)
    else:
        return map.mergeFaces(dart)

def mergeFacesCompletely(dart, doRemoveDegree2Nodes = True):
    """mergeFacesCompletely(dart, doRemoveDegree2Nodes = True)

    In contrast to the Euler operation mergeFaces(), this function
    removes all common edges of the two faces, not only the single
    edge belonging to dart.

    Furthermore, if the optional parameter doRemoveDegree2Nodes is
    True (default), all nodes whose degree is reduced to two will be
    merged into their surrounding edges.

    Returns the surviving face."""
    
    #print "mergeFacesCompletely(%s, %s)" % (dart, doRemoveDegree2Nodes)
    if dart.edge().isBridge():
        raise TypeError("mergeFacesCompletely(): dart belongs to a bridge!")
    map = dart.map()
    rightLabel = dart.rightFaceLabel()
    commonEdgeList = []
    for contourIt in dart.phiOrbit():
        if contourIt.rightFaceLabel() == rightLabel:
            if contourIt.edge().protection():
                return None
            commonEdgeList.append(contourIt)

    assert commonEdgeList, "mergeFacesCompletely(): no common edges found!"
    affectedNodes = []

    survivor = None
    for dart in commonEdgeList:
        affectedNodes.append(dart.startNodeLabel())
        affectedNodes.append(dart.endNodeLabel())
        assert dart.edge().protection() == 0
        if survivor == None:
            survivor = map.mergeFaces(dart) # first common edge
        else:
            assert survivor == map.removeBridge(dart)

    for nodeLabel in affectedNodes:
        node = map.node(nodeLabel)
        if not node: continue
        if node.degree == 0:
            map.removeIsolatedNode(node)
        if doRemoveDegree2Nodes and node.degree() == 2:
            d = node.anchor()
            if d.endNodeLabel() != node.label():
                map.mergeEdges(d)

    return survivor

def mergeFacesByLabel(map, label1, label2, doRemoveDegree2Nodes = True):
    """mergeFacesByLabel(map, label1, label2, doRemoveDegree2Nodes = True)

    Similar to mergeFacesCompletely() (which is called to perform the
    action), but is parametrized with two face labels and finds a
    common dart of these labels.

    Returns the surviving face (or None if no common edge was found)."""
    
    face1 = map.face(label1)
    if not face1:
        sys.stderr.write("mergeFacesByLabel: face with label1 = %d does not exist!\n"
                         % (label1, ))
        return
    if not map.face(label2):
        sys.stderr.write("mergeFacesByLabel: face with label2 = %d does not exist!\n"
                         % (label2, ))
        return
    for dart in face1.contours():
        for contourIt in dart.phiOrbit():
            if contourIt.rightFaceLabel() == label2:
                return mergeFacesCompletely(contourIt, doRemoveDegree2Nodes)

# --------------------------------------------------------------------

def thresholdMergeCost(map, mergeCostFunctor, maxCost, costs = None, q = None):
    """thresholdMergeCost(map, mergeCostFunctor, maxCost, costs = None, q = None)

    Merges faces of the given map in order of increasing costs until
    no more merges are assigned a cost <= maxCost by the
    mergeCostFunctor.

    The mergeCostFunctor is called on a dart for each edge to
    determine the cost of removing that edge.  This is used to
    initialize a DynamicCostQueue, which is used internally in order
    to always remove the edge with the lowest cost assigned.  The
    function returns a pair of the number of operations performed and
    the DynamicCostQueue, and you may pass the latter as optional
    argument 'd' into a subsequent call of thresholdMergeCost (with
    the same map and an increased maxCost) in order to re-use it
    (mergeCostFunctor is not used then).

    If the optional argument costs is given, it should be a mutable
    sequence which is append()ed the costs of each performed
    operation."""
    
    result = 0

    if q == None:
        q = hourglass.DynamicCostQueue(map.maxEdgeLabel()+1)
        for edge in map.edgeIter():
            q.insert(edge.label(), mergeCostFunctor(edge.dart()))
        
    while not q.empty():
        edgeLabel, cost = q.pop()
        if cost > maxCost:
            break

        edge = map.edge(edgeLabel)
        if not edge or edge.protection():
            continue
        d = edge.dart()
        if edge.isBridge():
            survivor = removeBridge(d)
        else:
            survivor = mergeFacesCompletely(d)
            if survivor:
                for anchor in survivor.contours():
                    for dart in anchor.phiOrbit():
                        q.setCost(dart.edgeLabel(),
                                  mergeCostFunctor(dart))
        if survivor:
            result += 1
            if costs:
                costs.append(cost)
    
    return result, q

# --------------------------------------------------------------------

def holeComponent(dart, includeExterior = False):
    """Returns list of all faces of the connected combinatorial map of
    the given dart.  If includeExterior is True, the first element of
    the result is dart.leftFace() (which is expected to be the face
    the above map is embedded in).  By default, this face is not
    returned.

    This is supposed to be called with e.g. the inner (hole) anchor
    darts of a face."""

    result = [dart.leftFace()]
    seen   = {dart.leftFaceLabel() : None}
    border = dict.fromkeys([d.rightFace() for d in dart.phiOrbit()])
    while border:
        face = border.popitem()[0]
        if seen.has_key(face.label()):
            continue
        result.append(face)
        seen[face.label()] = None
        for dart in face.contour().phiOrbit():
            border[dart.rightFace()] = None
    if not includeExterior:
        del result[0]
    return result

def showHomotopyTree(face, indentation = ""):
    """showHomotopyTree(root)

    Prints homotopy tree to stdout, starting with a given face
    as root.  You can also pass a GeoMap as parameter, in which case
    the tree will start with its infinite face."""
    
    if not hasattr(face, "label"):
        face = face.face(0)
    print indentation + str(face)
    for anchor in face.holeContours():
        for hole in holeComponent(anchor):
            showHomotopyTree(hole, indentation + "  ")
