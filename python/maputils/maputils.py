import hourglass, sys, math, time
from vigra import * # FIXME?

import flag_constants

# --------------------------------------------------------------------
#                            edge protection
# --------------------------------------------------------------------

from weakref import ref

class EdgeProtection(object):
    def __init__(self, map):
        self._attach(map)
        # prevent cycles if this is an attribute of the map:
        self._map = ref(map) # only needed for pickle support

    def _attach(self, map):
        assert not hasattr(self, "_attachedHooks"), \
               "trying to attach to more than one GeoMap?!"
        self._attachedHooks = (
            map.addMergeFacesCallbacks(self.preRemoveEdge, None),
            map.addRemoveBridgeCallbacks(self.preRemoveEdge, None),
            map.addMergeEdgesCallbacks(self.preMergeEdges, None))

    def detachHooks(self):
        for cb in self._attachedHooks:
            cb.disconnect()

    def preRemoveEdge(self, dart):
        "do not allow removal of protected edges"
        return not dart.edge().flag(flag_constants.ALL_PROTECTION)

    def preMergeEdges(self, dart):
        "only allow edge merging if the edges carry the same flags"
        return (dart.edge().flags() ==
                dart.clone().nextSigma().edge().flags())

    def __getinitargs__(self):
        return (self._map(), )

def protectFace(face, protect = True):
    """protectFace(face, protect = True)

    Sets the PROTECTED_FACE of 'face' according to 'protect'.
    Subsequently, sets the CONTOUR_PROTECTION of each edge in the
    contours of 'face' iff either of the adjacent faces is
    protected."""
    
    face.setFlag(flag_constants.PROTECTED_FACE, protect)
    for dart in contourDarts(face):
        dart.edge().setFlag(flag_constants.CONTOUR_PROTECTION,
                            protect or dart.rightFace().flag(
            flag_constants.PROTECTED_FACE))

# --------------------------------------------------------------------
#                      subpixel watershed functions
# --------------------------------------------------------------------

def filterSaddlePoints(rawSaddles, biSIV, threshold, maxDist):
    """filterSaddlePoints(rawSaddles, biSIV, threshold, maxDist) -> list

    Filters saddle points first by checking the threshold in the
    boundary indicator SplineImageView 'biSIV', then by removing
    duplicates (points which are nearer to previous points than
    maxDist - IOW: Points that come first win.)

    Returns a list of the remaining points together with their indices
    in the original array (in enumerate() fashion)."""

    maxSquaredDist = math.sq(maxDist)
    result = []
    sorted = []
    for k, saddle in enumerate(rawSaddles):
        if k == 0:
            assert not saddle, "rawSaddles[0] expected to be None"
            continue
        v = biSIV[saddle]
        if threshold and v < threshold:
            continue
        sorted.append((-v, k, saddle))

    sorted.sort()

    knownSaddles = hourglass.PositionedMap()
    for _, k, saddle in sorted:
        if not knownSaddles(saddle, maxSquaredDist):
            result.append((k, saddle))
            knownSaddles.insert(saddle, saddle)

    return result

def subpixelWatershedData(spws, biSIV = None, threshold = None, mask = None,
                          maxSaddleDist = 0.12, # sensible: ssMinDist
                          perpendicularDistEpsilon = 0.1, maxStep = 0.1):
    """subpixelWatershedData(spws, biSIV, threshold) -> tuple

    Calculates and returns a pair with a list of maxima and a list of
    edges (composed flowline pairs) from the SubPixelWatersheds object
    'spws'.  Each edge is represented with a quadrupel of
    startNodeIndex, endNodeIndex, edge polygon and the index of the
    original saddle within that polygon.  The first edge is None.
    That output is suitable and intended for addFlowLinesToMap().

    Gives verbose output during operation and uses filterSaddlePoints
    to filter out saddle points below the threshold within the
    SplineImageView 'biSIV' which should contain the boundary
    indicator and to filter out duplicate saddlepoints (pass the
    optional argument 'maxSaddleDist' to change the default of 0.1
    here).

    Each found edge polygon is simplified using simplifyPolygon with
    perpendicularDistEpsilon = 0.1 and maxStep = 0.1 (use the optional
    parameters with the same names to change the default)."""
    
    sys.stdout.write("- finding critical points..")
    c = time.clock()

    if mask:
        spws.findCriticalPoints(mask)

    rawSaddles = spws.saddles()
    maxima     = spws.maxima()
    sys.stdout.write("done. (%ss, %d maxima, %d saddles, %d minima)\n" % (
        time.clock()-c,
        len(maxima)-1, len(rawSaddles)-1, len(spws.minima())-1))

    def calculateEdge(index):
        flowline = spws.edge(index)
        if not flowline:
            return flowline
        sn, en, poly, saddleIndex = flowline
        if perpendicularDistEpsilon:
            simple = hourglass.simplifyPolygon(
                poly[:saddleIndex+1], perpendicularDistEpsilon, maxStep)
            newSaddleIndex = len(simple)-1
            simple.extend(hourglass.simplifyPolygon(
                poly[saddleIndex:], perpendicularDistEpsilon, maxStep))
            poly = simple
            saddleIndex = newSaddleIndex
        return (sn, en, poly, saddleIndex, index)

    saddles = filterSaddlePoints(rawSaddles, biSIV, threshold, maxSaddleDist)
    prefix = "- following %d/%d edges.." % (len(saddles), len(rawSaddles))
    c = time.clock()
    flowlines = [None]
    percentGranularity = len(saddles) / 260 + 1
    for k, _ in saddles:
        edgeTuple = calculateEdge(k)
        if edgeTuple:
            flowlines.append(edgeTuple)
        if k % percentGranularity == 0:
            sys.stderr.write("%s %d%%\r" % (
                prefix, 100 * (len(flowlines)-1) / len(saddles), ))
    print "%s done (%ss)." % (prefix, time.clock()-c)

    if len(flowlines)-1 < len(saddles):
        print "  (%d edges parallel to border removed)" % (
            len(saddles) - (len(flowlines)-1))

    if perpendicularDistEpsilon:
        print "  (simplified edges with perpendicularDistEpsilon = %s, maxStep = %s)" % (
            perpendicularDistEpsilon, maxStep)

    return maxima, flowlines

def _handleUnsortable(map, unsortable):
    # remove unsortable self-loops (most likely flowlines which
    # did not get far and were connected to the startNode by
    # addFlowLinesToMap())
    for group in unsortable:
        i = 0
        while i < len(group):
            dart = map.dart(group[i])
            if dart.edge().isLoop() and -dart.label() in group:
                group.remove(dart.label())
                group.remove(-dart.label())
            else:
                i += 1
    try:
        while True:
            unsortable.remove([])
    except ValueError:
        pass
    assert not unsortable, "unhandled unsortable edges occured"

def subpixelWatershedMap(maxima, flowlines, imageSize,
                         perpendicularDistEpsilon = 0.1, maxStep = 0.1,
                         borderConnectionDist = 0.1,
                         ssStepDist = 0.2, ssMinDist = 0.12,
                         performBorderClosing = True,
                         performEdgeSplits = True,
                         wsStatsSpline = None,
                         minima = None,
                         Map = hourglass.GeoMap):

    spmap = Map(maxima, [], imageSize)

    deleted = addFlowLinesToMap(flowlines, spmap)
    if deleted:
        print "  skipped %d broken flowlines." % len(deleted)

    if performBorderClosing:
        print "  adding border edges and EdgeProtection..."
        connectBorderNodes(spmap, borderConnectionDist)
        spmap.edgeProtection = EdgeProtection(spmap)

    print " ",; removeCruft(spmap, 1) # FIXME: removeIsolatedNodes
    unsortable = spmap.sortEdgesEventually(
        ssStepDist, ssMinDist, performEdgeSplits)
    _handleUnsortable(spmap, unsortable)

    if wsStatsSpline:
        c = time.clock()
        print "  initializing watershed statistics...",
        from statistics import WatershedStatistics
        spmap.wsStats = WatershedStatistics(
            spmap, flowlines, wsStatsSpline)
        print " (%ss)" % (time.clock()-c, )

    if performEdgeSplits:
        print "  splitting & joining parallel edges..."
        spmap.splitParallelEdges()

    spmap.initializeMap()

    if minima:
        assert wsStatsSpline, \
            "minima given, but basin statistics need a wsStatsSpline, too"
        c = time.clock()
        print "  initializing watershed basin statistics...",
        from statistics import WatershedBasinStatistics
        spmap.wsBasinStats = WatershedBasinStatistics(
            spmap, minima[1:], wsStatsSpline)
        print " (%ss)" % (time.clock()-c, )

    return spmap

def addFlowLinesToMap(edges, map):
    """addFlowLinesToMap(edges, map)

    This function expects `edges` to be a list of
    SubPixelWatersheds.edge() return values and adds edges for each
    flowline.  It contains some special handling of flowlines:

    * Additional Nodes may be added for flowlines that did not end in
      a maximum - extra care is taken not to insert multiple Nodes at
      (nearly) the same position.

    * Self-loops with area zero are not added.

    Returns edgeTuples that could not be added to the GeoMap."""
    
    # Node 0 conflicts with our special handling of 0 values:
    assert not map.node(0), \
           "addFlowLinesToMap: Node with label zero should not exist!"

    result = []
    for edgeLabel, edgeTuple in enumerate(edges):
        if not edgeTuple:
            continue

        startNodeLabel = edgeTuple[0]
        endNodeLabel = edgeTuple[1]
        if startNodeLabel == -2 and endNodeLabel == -2:
            continue # unwanted edge parallel to border

        # be careful not to modify the original 'edges' passed:
        points = hourglass.Polygon(edgeTuple[2])

        assert len(points) >= 2, "edges need to have at least two (end-)points"

        # FIXME: now, we allow jumping back, replacing the last edge
        # point with the node position.  For this, a distance of 0.25
        # could be too much - remove more points?

        # if flowlines do not end in a maximum (their end labels
        # are negative), we assign end nodes to them:
        if startNodeLabel <= 0:
            pos = points[0]
            nearestNode = map.nearestNode(pos, 0.25)
            # if there is already a node nearby and in the right
            # direction (prevents loops), attach edge to it:
            if nearestNode:
                diff = nearestNode.position() - pos
                # we are not interested in creating short, unsortable self-loops:
                if nearestNode.label() != endNodeLabel:
                    startNode = nearestNode
                    if dot(diff, pos - points[1]) >= 0: # don't jump back
                        if diff.squaredMagnitude():
                            # include node position if not present
                            points.insert(0, nearestNode.position())
                    else:
                        points[0] = nearestNode.position()
            else: # no suitable Node found -> add one
                startNode = map.addNode(pos)
            startNodeLabel = startNode.label()
        else:
            startNode = map.node(startNodeLabel)
            assert startNode, "invalid startNodeLabel!"

        # ..handle Edge end the same as the start:
        if endNodeLabel <= 0:
            pos = points[-1]
            nearestNode = map.nearestNode(pos, 0.25)
            if nearestNode:
                diff = nearestNode.position() - pos
                if nearestNode.label() != startNodeLabel:
                    endNode = nearestNode
                    if dot(diff, pos - points[-2]) >= 0:
                        if diff.squaredMagnitude():
                            points.append(nearestNode.position())
                    else:
                        points[-1] = nearestNode.position()
            else:
                endNode = map.addNode(pos)
            endNodeLabel = endNode.label()
        else:
            endNode = map.node(endNodeLabel)
            assert endNode, "invalid endNodeLabel!"

        if len(points) == 2 and points[0] == points[1]:
            # don't add self-loops with area zero
            result.append((startNodeLabel, endNodeLabel, points))
            continue
        
        if startNodeLabel > 0 and endNodeLabel > 0:
            map.addEdge(startNode, endNode, points, edgeLabel)
        else:
            result.append((startNodeLabel, endNodeLabel, points))

    return result

# --------------------------------------------------------------------
#                        map creation helpers
# --------------------------------------------------------------------

def connectBorderNodes(map, epsilon,
                       samePosEpsilon = 1e-6, aroundPixels = False):
    """connectBorderNodes(map, epsilon,
                       samePosEpsilon = 1e-6, aroundPixels = False)

    Inserts border edges around (see below for details) the image into
    `map`; all Nodes which are less than `epsilon` pixels away from
    that border are connected to it.  If the distance is larger than
    samePosEpsilon, an extra perpendicular edge is used for the
    connection, otherwise the border will run through the existing
    node.

    The optional parameter aroundPixels determines the position of the
    border: If it is set to False (default), the border will run
    through the pixel centers, i.e. from 0/0 to w-1/h-1 (where w,h =
    map.imageSize()).  Otherwise, the border will run around the pixel
    facets, i.e. be 0.5 larger in each direction."""
    
    dist = aroundPixels and 0.5 or 0.0
    x1, y1 = -dist, -dist
    x2, y2 = map.imageSize()[0] - 1 + dist, map.imageSize()[1] - 1 + dist

    left = []
    right = []
    top = []
    bottom = []
    for node in map.nodeIter():
        p = node.position()
        if   p[0] > x2 - epsilon:
            right.append(node)
        elif p[0] < x1 + epsilon:
            left.append(node)
        elif p[1] > y2 - epsilon:
            bottom.append(node)
        elif p[1] < y1 + epsilon:
            top.append(node)

    def XPosCompare(node1, node2):
        return cmp(node1.position()[0], node2.position()[0])

    def YPosCompare(node1, node2):
        return cmp(node1.position()[1], node2.position()[1])

    left.sort(YPosCompare); left.reverse()
    right.sort(YPosCompare)
    top.sort(XPosCompare)
    bottom.sort(XPosCompare); bottom.reverse()

    borderEdges = []
    lastNode = None

    lastPoints = [Vector2(x1, y1)]
    for node in top:
        thisPoints = []
        if node.position()[1] > y1 + samePosEpsilon:
            thisPoints.append(Vector2(node.position()[0], y1))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(Vector2(x2, y1))
    for node in right:
        thisPoints = []
        if node.position()[0] < x2 - samePosEpsilon:
            thisPoints.append(Vector2(x2, node.position()[1]))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(Vector2(x2, y2))
    for node in bottom:
        thisPoints = []
        if node.position()[1] < y2 - samePosEpsilon:
            thisPoints.append(Vector2(node.position()[0], y2))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(Vector2(x1, y2))
    for node in left:
        thisPoints = []
        if node.position()[0] > x1 + samePosEpsilon:
            thisPoints.append(Vector2(x1, node.position()[1]))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    if not borderEdges:
        return

    lastPoints.extend(borderEdges[0][0])
    borderEdges[0] = (lastPoints, borderEdges[0][1])

    endNode = lastNode
    for points, node in borderEdges:
        startNode = endNode
        endNode = node
        map.addEdge(startNode, endNode, points) \
                        .setFlag(flag_constants.BORDER_PROTECTION)

# --------------------------------------------------------------------
#                         consistency checks
# --------------------------------------------------------------------

class _CLCFaceLookup(object):
    def __init__(self, map):
        self._map = map
        self.errorCount = 0
        self.errorLabels = []
        self.pixelAreas = [0] * map.maxFaceLabel()

    def __call__(self, faceLabel):
        if faceLabel < 0:
            return
        faceLabel = int(faceLabel)
        self.pixelAreas[faceLabel] += 1
        try:
            assert self._map.face(faceLabel) != None
        except:
            self.errorCount += 1
            if not faceLabel in self.errorLabels:
                self.errorLabels.append(faceLabel)

def checkLabelConsistency(map):
    """Checks that no unknown positive/face labels occur in the
    labelimage, and that each face's pixelArea() is correct."""
    fl = _CLCFaceLookup(map)
    labelImage = map.labelImage()
    inspectImage(labelImage, fl)
    if fl.errorCount:
        sys.stderr.write("labelImage contains %d pixels with unknown faces!\n" % (
            fl.errorCount, ))
        sys.stderr.write("  unknown face labels found: %s\n" % (fl.errorLabels, ))
        if fl.errorCount < 40:
            for p in labelImage.size():
                if int(labelImage[p]) in fl.errorLabels:
                    print "   label %d at %s" % (int(labelImage[p]), p)
    result = (fl.errorCount == 0)
    for face in map.faceIter(skipInfinite = True):
        if face.pixelArea() != fl.pixelAreas[face.label()]:
            sys.stderr.write(
                "pixel area of Face %d is wrong (%d, should be %d)\n" %
                (face.label(), face.pixelArea(), fl.pixelAreas[face.label()]))
            result = False
    return result

# --------------------------------------------------------------------

def showMapStats(map):
    """showMapStats(map)
    Dumps number of cells and points to stdout (returns nothing.)"""
    pointCount = 0
    totalLength = 0.0
    for edge in map.edgeIter():
        pointCount += len(edge)
        totalLength += edge.length()

    if map.edgeCount == 0:
        print "empty map!"
        return

    print ("%d nodes, %d edges, and %d faces with a total of %d points\n  " + \
          "(mean p/edge: %.2f, mean dist.: %.2f, density: %.2f p/px)") % (
        map.nodeCount, map.edgeCount, map.faceCount, pointCount,
        float(pointCount) / map.edgeCount,
        totalLength / (pointCount - map.edgeCount),
        float(pointCount)/map.imageSize().area())
    
    if hasattr(map, "deleted") and map.deleted:
        print "%d edges were deleted (e.g. at image border)." % (
            len(map.deleted), )
    if hasattr(map, "unsortable") and map.unsortable:
        print "%d unsortable groups of edges occured." % (
            len(map.unsortable), )
    if hasattr(map, "unembeddableContours") and map.unembeddableContours:
        print "%d outer contours could not be embedded into " \
              "their surrounding faces!" % (len(map.unembeddableContours), )

def degree2Nodes(map):
    return [node for node in map.nodeIter() if node.degree() == 2]

# --------------------------------------------------------------------
#                    basic map cleanup operations
# --------------------------------------------------------------------

def removeCruft(map, what = 3, doChecks = False):
    """removeCruft(map, what = 3, doChecks = False)
    
    what is a bit-combination of
    1: for the removal of degree 0-nodes (default)
    2: removal of degree 2-nodes (default)
    4: removal of bridges
    8: removal of edges (i.e. all non-protected)

    if doChecks is True, consistency checks are performed after every
    operation.  As soon as that fails, removeCruft returns False.

    After normal operation, removeCruft returns the number of
    operations performed."""

    class OperationCounter(object):
        def __init__(self):
            self.count = 0

        def perform(self, op, dart):
            if op(dart):
                self.count += 1
            return True

    class CarefulCounter(OperationCounter):
        def perform(self, op, dart):
            OperationCounter.perform(self, op, dart)
            return checkConsistency(map)

    if doChecks:
        result = CarefulCounter()
    else:
        result = OperationCounter()

    if what & 8:
        for edge in map.edgeIter():
            if edge.leftFaceLabel() != edge.rightFaceLabel():
                if not result.perform(map.mergeFaces, edge.dart()):
                    return False

    if what & 4:
        for edge in map.edgeIter():
            if edge.leftFaceLabel() == edge.rightFaceLabel():
                if not result.perform(map.removeBridge, edge.dart()):
                    return False

    if what & 2:
        for node in map.nodeIter():
            if node.degree() == 2 and \
                   (node.anchor().endNode() != node):
                if not result.perform(map.mergeEdges, node.anchor()):
                    return False

    if what & 1:
        for node in map.nodeIter():
            if node.degree() == 0:
                if not result.perform(map.removeIsolatedNode, node):
                    return False

    print "removeCruft(): %d operations performed." % result.count
    return result.count

def removeSmallRegions(map, minArea):
    result = 0
    for face in map.faceIter(skipInfinite = True):
        if face.area() < minArea:
            po = list(face.contour().phiOrbit())
            if len(po) > 1:
                sys.stderr.write("WARNING: removeSmallRegions() has to decide about neighbor\n  which %s is to be merged into!\n" % face)
            for dart in po:
                if mergeFacesCompletely(dart):
                    result += 1
                    break
    return result

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
            if contourIt.edge().flag(flag_constants.ALL_PROTECTION):
                return None
            commonEdgeList.append(contourIt)

    assert commonEdgeList, "mergeFacesCompletely(): no common edges found!"
    affectedNodes = []

    survivor = None
    for dart in commonEdgeList:
        affectedNodes.append(dart.startNodeLabel())
        affectedNodes.append(dart.endNodeLabel())
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
        if not edge or edge.flag(flag_constants.ALL_PROTECTION):
            continue
        d = edge.dart()
        if edge.isBridge():
            survivor = removeBridge(d)
        else:
            survivor = mergeFacesCompletely(d)
            if survivor:
                for dart in contourDarts(survivor):
                    q.setCost(dart.edgeLabel(),
                              mergeCostFunctor(dart))
        if survivor:
            result += 1
            if costs != None:
                costs.append(cost)
    
    return result, q

# --------------------------------------------------------------------
#                   topological utility functions
# --------------------------------------------------------------------

def contourDarts(face):
    """contourDarts(face)

    Generator function which iterates over all darts within all
    contours, i.e. substitues the following standard loop construct:

      for anchor in face.contours():
          for dart in anchor.phiOrbit():
              doSth(dart)

    Thus, you can now write:

      for dart in contourDarts(face):
          doSth(dart)"""

    for anchor in face.contours():
        for dart in anchor.phiOrbit():
            yield dart

def neighborFaces(face):
    """neighborFaces(face)

    Generator function which yields exactly one dart for each adjacent
    face (even in case of multiple common edges)."""
    
    seen = {face.label(): True}
    for dart in contourDarts(face):
        fl = dart.rightFaceLabel()
        if not fl in seen:
            seen[fl] = True
            yield dart

def holeComponent(dart, includeExterior = False):
    """holeComponent(dart, includeExterior = False) -> list

    Returns list of all faces of the connected combinatorial map of
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

def edgeAtBorder(edge):
    return edge.flag(flag_constants.BORDER_PROTECTION)

def nodeAtBorder(node):
    """nodeAtBorder(node) -> bool

    Returns True iff the Node is adjacent to at least one edge
    marked with BORDER_PROTECTION."""
    if node.isIsolated():
        return False
    for dart in node.anchor().sigmaOrbit():
        if edgeAtBorder(dart.edge()):
            return True
    return False

# --------------------------------------------------------------------

from heapq import heappush, heappop

def minimumSpanningTree(map, edgeCosts):
    """minimumSpanningTree(map, edgeCosts)

    Given a cost associated with each edge of the map, this function
    finds the minimum spanning tree of the map's boundary graph (None
    in edgeCosts is allowed and is handled as if the corresponding
    edge was missing).  The result is a modified copy of the edgeCosts
    list, with all non-MST-edges set to None.  This can be used for
    the waterfall algorithm by Meyer and Beucher, see waterfall().

    The complexity is O(edgeCount*faceCount), but the performance
    could be improved a little."""

    faceLabels = [None] * map.maxFaceLabel()
    for face in map.faceIter():
        faceLabels[face.label()] = face.label()

    print "- initializing priority queue for MST..."
    heap = []
    for edge in map.edgeIter():
        edgeLabel = edge.label()
        cost = edgeCosts[edgeLabel]
        if cost != None:
            heappush(heap, (cost, edgeLabel))

    print "- building MST..."
    result = list(edgeCosts)
    while heap:
        _, edgeLabel = heappop(heap)
        edge = map.edge(edgeLabel)
        lfl = faceLabels[edge.leftFaceLabel()]
        rfl = faceLabels[edge.rightFaceLabel()]
        if lfl != rfl:
            for i in range(len(faceLabels)): # this could be optimized
                if faceLabels[i] == rfl:
                    faceLabels[i] = lfl
        else:
            result[edgeLabel] = None

    return result

def seededRegionGrowingStatic(map, faceLabels, edgeCosts):
    """seededRegionGrowingStatic(map, faceLabels, edgeCosts)

    Given a set of seed labels != None, extend labels to neighbor
    faces via edges with costs != None.  Stop when no more such
    growing is possible and return the resulting faceLabels.

    The 'static' postfix indicates that the edgeCosts do not change
    during the process.  (Furthermore, the GeoMap is not changed.)"""

    faceLabels = list(faceLabels)

    heap = []
    for face in map.faceIter():
        if not faceLabels[face.label()]:
            continue
        for dart in contourDarts(face):
            if not faceLabels[dart.rightFaceLabel()]:
                cost = edgeCosts[dart.edgeLabel()]
                if cost != None:
                    heappush(heap, (cost, dart.label()))
    while heap:
        _, dartLabel = heappop(heap)
        dart = map.dart(dartLabel)
        if faceLabels[dart.rightFaceLabel()]:
            continue # labelled in the meantime
        faceLabels[dart.rightFaceLabel()] = faceLabels[dart.leftFaceLabel()]
        for dart in contourDarts(dart.rightFace()):
            if not faceLabels[dart.rightFaceLabel()]:
                cost = edgeCosts[dart.edgeLabel()]
                if cost != None:
                    heappush(heap, (cost, dart.label()))
    
    return faceLabels

def _addNeighborsToHeap(face, heap, mergeCostMeasure):
    for dart in neighborFaces(face):
        neighbor = dart.rightFace()
        if not neighbor.flag(flag_constants.SRG_SEED |
                             flag_constants.SRG_BORDER):
            cost = mergeCostMeasure(dart)
            heappush(heap, (cost, neighbor.label())) # FIXME: cost not updated
            neighbor.setFlag(flag_constants.SRG_BORDER) # prevents second heappush

def seededRegionGrowing(map, mergeCostMeasure):
    """seededRegionGrowing(map, mergeCostMeasure)"""

    for face in map.faceIter():
        face.setFlag(flag_constants.SRG_BORDER, False)

    heap = []
    for face in map.faceIter():
        if face.flag(flag_constants.SRG_SEED):
            _addNeighborsToHeap(face, heap, mergeCostMeasure)

    while heap:
        _, faceLabel = heappop(heap)
        face = map.face(faceLabel)

        best = None
        for dart in neighborFaces(face):
            neighbor = dart.rightFace()
            if neighbor.flag(flag_constants.SRG_SEED):
                cost = mergeCostMeasure(dart)
                if not best or cost < best[0]:
                    best = (cost, dart)

        survivor = mergeFacesCompletely(best[1])
        if survivor:
            survivor.setFlag(flag_constants.SRG_BORDER, False)
            survivor.setFlag(flag_constants.SRG_SEED)
            _addNeighborsToHeap(survivor, heap, mergeCostMeasure)

def regionalMinima(map, mst):
    """regionalMinima(map, mst)

    Returns an array with face labels, where only 'regional minima' of
    the given MSt are labelled.  See the article 'Fast Implementation
    of Waterfall Based on Graphs' by Marcotegui and Beucher."""

    faceLabels = [None] * map.maxFaceLabel()
    for edge in map.edgeIter():
        edgeCost = mst[edge.label()]
        if edgeCost == None:
            continue
        isMinimum = True
        for start in edge.dart().alphaOrbit():
            for dart in contourDarts(face):
#               if dart == start: # not needed for strict comparison below
#                   continue
                otherCost = mst[dart.edgeLabel()]
                if otherCost != None and otherCost < edgeCost:
                    isMinimum = False
                    break
            if not isMinimum:
                break
        if isMinimum:
            faceLabels[edge.leftFaceLabel()] = edge.label()
            faceLabels[edge.rightFaceLabel()] = edge.label()
    return faceLabels

def waterfallLabels(map, edgeCosts, mst = None):
    """waterfallLabels(map, edgeCosts, mst = None)

    Perform a single iteration of the waterfall algorithm by
    Marcotegui and Beucher and return an array with face labels."""

    if not mst:
        mst = minimumSpanningTree(map, edgeCosts)

    print "- finding regional minima..."
    faceLabels = regionalMinima(map, mst)

    print "- extending faceLabels via MST..."
    return seededRegionGrowingStatic(map, faceLabels, mst)

def waterfall(map, edgeCosts, mst = None):
    """waterfall(map, edgeCosts, mst = None)

    Perform a single iteration of the waterfall algorithm by
    Marcotegui and Beucher and perform the face merging in the given
    map."""

    c = time.clock()
    faceLabels = waterfallLabels(map, edgeCosts, mst)
    print "- merging regions according to labelling..."
    for edge in map.edgeIter():
        if edge.flag(flag_constants.ALL_PROTECTION):
            continue
        if edge.isBridge():
            print "ERROR: waterfall should not be called on maps with bridges!"
            removeBridge(edge.dart())
            continue
        lfl = faceLabels[edge.leftFaceLabel()]
        rfl = faceLabels[edge.rightFaceLabel()]
        if lfl == rfl:
            mergeFacesCompletely(edge.dart())
    print "  total waterfall() time: %ss." % (time.clock() - c, )

def mst2map(mst, map):
    """Visualize a minimumSpanningTree() result as a GeoMap.  Note
    that the nodes are put in the bbox centers, which is wildly
    inaccurate, especially without edge splitting."""
    nodePositions = [None] * map.maxFaceLabel()
    for face in map.faceIter(skipInfinite = True):
        bbox = face.boundingBox()
        nodePositions[face.label()] = bbox.begin() + bbox.size() / 2
    edgeTuples = [None]
    for edgeLabel in mst:
        edge = map.edge(edgeLabel)
        snl = edge.leftFaceLabel()
        enl = edge.rightFaceLabel()
        edgeTuples.append((snl, enl, [nodePositions[snl], nodePositions[enl]]))
    result = hourglass.GeoMap(nodePositions, edgeTuples, map.imageSize())
    result.initializeMap(False)
    return result

# --------------------------------------------------------------------
#                               tests
# --------------------------------------------------------------------

# if __name__ == '__main__':
    
