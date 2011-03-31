"""maputils - general (i.e.topological) GeoMap utilities and segmentation algorithms"""

import vigra, geomap, sys, math, time, weakref, copy

import flag_constants, progress, sivtools

# --------------------------------------------------------------------
#                            edge protection
# --------------------------------------------------------------------

from geomap import EdgeProtection

def protectFace(face, protect = True, flag = flag_constants.PROTECTED_FACE):
    """Sets the PROTECTED_FACE flag of `face` according to `protect`.
    Subsequently, sets the CONTOUR_PROTECTION of each edge in the
    contours of `face` iff either of the adjacent faces is
    protected."""
    
    face.setFlag(flag, protect)
    for dart in contourDarts(face):
        dart.edge().setFlag(flag_constants.CONTOUR_PROTECTION,
                            protect or dart.rightFace().flag(flag))

# --------------------------------------------------------------------
#                      subpixel watershed functions
# --------------------------------------------------------------------

class PassValueFilter(object):
    def __init__(self, biSIV, threshold):
        self._biSIV = biSIV
        self._threshold = threshold

    def __call__(self, position):
        return self._biSIV[position] >= self._threshold

    def __str__(self):
        return "pass value >= %s" % self._threshold

class SaddleTensorFeature(object):
    def __init__(self, gmSIV, stSIV):
        self._hessian = sivtools.HessianSIVProxy(gmSIV)
        self._stSIV = stSIV

    def saddleDir(self, position):
        """Return tangent direction of flowline starting at the given
        saddle position.  (Direction of larger eigenvector.)"""
        gxx, gxy, gyy = self._hessian[position]
        return 0.5 * math.atan2(-2.0*gxy, gxx-gyy)

    def _orthoGrad2(self, position, onlyDir = False):
        # tangent angle of flowline:
        ta_fl = self.saddleDir(position)
        c = math.cos(ta_fl)
        s = math.sin(ta_fl)

        # component of structure tensor in normal direction:
        s11, s12, s22 = self._stSIV[position]
        og2 = max(0, s11*vigra.sq(s) - 2*s12*c*s + s22*vigra.sq(c))

        if onlyDir:
            # large eigenvalue of structure tensor:
            ev = 0.5 * (s11 + s22 + math.sqrt(vigra.sq(s11 - s22)
                                              + 4.0*vigra.sq(s12)))
            return og2 / ev

        return og2

    def orthogonalGradient(self, position):
        return math.sqrt(self._orthoGrad2(position))

    def directionMatch(self, position):
        return self._orthoGrad2(position, True)

class SaddleOrthogonalGradientFilter(SaddleTensorFeature):
    def __init__(self, gmSIV, stSIV, threshold):
        SaddleTensorFeature.__init__(self, gmSIV, stSIV)
        self._threshold = threshold

    def __call__(self, position):
        return math.sqrt(self._orthoGrad2(position)) >= self._threshold

    def __str__(self):
        return "orthogonal gradient >= %s" % self._threshold

class SaddleDirectionMatchFilter(SaddleTensorFeature):
    def __init__(self, gmSIV, stSIV, threshold):
        SaddleTensorFeature.__init__(self, gmSIV, stSIV)
        self._threshold = threshold

    def __call__(self, position):
        return self._orthoGrad2(position, True) >= self._threshold

    def __str__(self):
        return "direction match >= %s" % self._threshold

class FilterGroup(tuple):
    def __call__(self, position):
        for filter in self:
            if not filter(position):
                return False
        return True

    def __str__(self):
        return " and ".join(map(str, self))

def filterSaddlePoints(rawSaddles, biSIV, filter, maxDist):
    """Filter saddle points first by checking the passValue in the
    boundary indicator SplineImageView `biSIV`, then by removing
    duplicates (points which are nearer to previous points than
    `maxDist` - IOW: Points that come first win.)

    Return a list of the remaining points together with their indices
    in the original array (in enumerate() fashion)."""

    maxSquaredDist = math.sq(maxDist)
    result = []
    sorted = []
    for k, saddle in enumerate(rawSaddles):
        if k == 0:
            assert not saddle, "rawSaddles[0] expected to be None"
            continue
        if filter and not filter(saddle):
            continue
        sorted.append((-biSIV[saddle], k, saddle))

    sorted.sort()

    knownSaddles = geomap.PositionedMap()
    for _, k, saddle in sorted:
        if not knownSaddles(saddle, maxSquaredDist):
            result.append((k, saddle))
            knownSaddles.insert(saddle, saddle)

    return result

def subpixelWatershedData(spws, biSIV, filter = None, mask = None,
                          minSaddleDist = 0.1, # sensible: ssMinDist
                          perpendicularDistEpsilon = 0.1, maxStep = 0.1):
    """Calculates and returns a pair with a list of maxima and a list of
    edges (composed flowline pairs) from the SubPixelWatersheds object
    `spws`.  Each edge is represented with a quadrupel of
    startNodeIndex, endNodeIndex, edge polygon and the index of the
    original saddle within that polygon.  The first edge is None.
    That output is suitable and intended for addFlowLinesToMap().

    Gives verbose output during operation and uses filterSaddlePoints
    to filter out saddle points using the given filter criterion and
    to filter out duplicate saddlepoints (pass the optional argument
    `minSaddleDist` to change the default of 0.1 here).

    Each found edge polygon is simplified using simplifyPolygon with
    perpendicularDistEpsilon = 0.1 and maxStep = 0.1 (use the optional
    parameters with the same names to change the default)."""
    
    sys.stdout.write("- finding critical points..\n")
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
            simple = geomap.simplifyPolygon(
                poly[:saddleIndex+1], perpendicularDistEpsilon, maxStep)
            newSaddleIndex = len(simple)-1
            simple.extend(geomap.simplifyPolygon(
                poly[saddleIndex:], perpendicularDistEpsilon, maxStep))

            poly = simple
            saddleIndex = newSaddleIndex

        return (sn, en, poly, saddleIndex, index)

    saddles = filterSaddlePoints(rawSaddles, biSIV, filter, minSaddleDist)
    prefix = "- following %d/%d edges.." % (len(saddles), len(rawSaddles)-1)
    if filter:
        prefix += " (%s)" % (filter, )
    else:
        prefix += " (no saddle threshold)"

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
    loopCount = 0
    for group in unsortable:
        i = 0
        while i < len(group):
            dart = map.dart(group[i])
            if dart.edge().isLoop() and -dart.label() in group:
                group.remove(dart.label())
                group.remove(-dart.label())
                loopCount += 1
            else:
                i += 1
    if loopCount:
        print "*** BIG FAT WARNING: %d unsortable self-loops found. ***" % loopCount
    try:
        while True:
            unsortable.remove([])
    except ValueError:
        pass
    if unsortable:
        raise RuntimeError("unhandled unsortable edges occured (%s)" % unsortable)

def subpixelWatershedMapFromData(
    maxima, flowlines, imageSize,
    borderConnectionDist = 0.1,
    ssStepDist = 0.2, ssMinDist = 0.12,
    performBorderClosing = True,
    performEdgeSplits = True,
    wsStatsSpline = None,
    minima = None,
    cleanup = True):
    """`performEdgeSplits` may be a callable that is called before
    splitting edges, with the complete locals() dict as parameter
    (useful keys are 'spmap' and all parameters of this function).

    `cleanup` determines whether bridges, isolated nodes, and degree 2
    nodes are removed (default: True).  If False, the result also gets
    an additional attribute 'deleted' which contains the flowlines
    skipped by addFlowLinesToMap.
    """

    # copy flowline data (may be modified by addFlowLinesToMap for
    # proper WatershedStatistics):
    flowlines = list(flowlines)
    if flowlines[0] is not None:
        flowlines.insert(0, None)

    print "- initializing GeoMap from %d flowlines..." % (len(flowlines) - 1, )
    spmap = geomap.GeoMap(maxima, [], imageSize)

    deleted = addFlowLinesToMap(flowlines, spmap, imageSize,
                                minMaxBorderDist = ssMinDist)
    if deleted:
        print "  skipped %d flowlines (at border / degenerate loops), %d left" \
              % (len(deleted), spmap.edgeCount)

    if performBorderClosing:
        print "  adding border edges and EdgeProtection...",
        connectBorderNodes(spmap, borderConnectionDist)
        spmap.edgeProtection = EdgeProtection(spmap)
        print "(now %d edges)" % spmap.edgeCount

    if cleanup:
        removeIsolatedNodes(spmap)
    else:
        spmap.deleted = deleted

    unsortable = spmap.sortEdgesEventually(
        ssStepDist, ssMinDist, bool(performEdgeSplits))
    _handleUnsortable(spmap, unsortable)

    if wsStatsSpline:
        p = progress.StatusMessage("  initializing watershed statistics")
        from statistics import WatershedStatistics
        spmap.wsStats = WatershedStatistics(
            spmap, flowlines, wsStatsSpline)
        p.finish()

    if performEdgeSplits:
        if not isinstance(performEdgeSplits, int):
            performEdgeSplits(locals())
        p = progress.StatusMessage("  splitting & joining parallel edges")
        spmap.splitParallelEdges()
        p.finish()

    spmap.initializeMap()

    if minima:
        assert wsStatsSpline, \
            "minima given, but basin statistics need a wsStatsSpline, too"
        p = progress.StatusMessage("  initializing watershed basin statistics")
        from statistics import WatershedBasinStatistics
        spmap.wsBasinStats = WatershedBasinStatistics(
            spmap, minima[1:], wsStatsSpline)
        spmap.wsStats.setBasins(spmap.wsBasinStats)
        p.finish()

    if cleanup:
        removeBridges(spmap)
        mergeDegree2Nodes(spmap)

    return spmap

def subpixelWatershedMap(
    boundaryIndicator, splineOrder = 5,
    saddleThreshold = None, mask = None,
    saddleOrthoGrad = None, saddleDirMatch = None, biTensor = None,
    perpendicularDistEpsilon = 0.1, maxStep = 0.1,
    borderConnectionDist = 0.1,
    ssStepDist = 0.2, ssMinDist = 0.12,
    performBorderClosing = True,
    performEdgeSplits = True,
    initWSStats = True,
    cleanup = True):

    """`boundaryIndicator` should be an image with e.g. a gradient
    magnitude.
    
    If `mask` is False, all pixels are searched for critical
    points.  By default (mask == None), pixels below saddleThreshold/2
    are skipped.

    If `initWSStats` is True, the resulting GeoMap will have an
    attribute 'wsStats' with a `statistics.WatershedStatistics`
    instance.

    For an explanation of the parameters, look for the documentation
    of the two main internally used functions `subpixelWatershedData`
    (which creates the subpixel watershed data, i.e. critical points
    and flowlines) and `subpixelWatershedMapFromData` (which
    initializes the GeoMap and uses the sigma sorting and statistics
    parameters)."""

    SPWS = getattr(geomap, "SubPixelWatersheds%d" % splineOrder)
    spws = SPWS(boundaryIndicator)

    if hasattr(boundaryIndicator, "siv"):
        siv = boundaryIndicator.siv
    else:
        SIV = sivtools.sivByOrder(splineOrder)
        siv = SIV(boundaryIndicator)

    if mask is None and saddleThreshold or saddleOrthoGrad:
        threshold = saddleThreshold or saddleOrthoGrad
        mask = vigra.transformImage(
            boundaryIndicator, "\l x: x > %s ? 1 : 0" % (threshold/2, ))

    filters = []
    if saddleThreshold:
        filters.append(PassValueFilter(siv, saddleThreshold))
    if saddleOrthoGrad or saddleDirMatch:
        assert biTensor, "for saddle filtering using SOG/SDM, biTensor is needed!"
        if not hasattr(biTensor, "siv"):
            biTensor.siv = sivtools.TensorSIVProxy(biTensor)
        if saddleOrthoGrad:
            filters.append(
                SaddleOrthogonalGradientFilter(siv, biTensor.siv, saddleOrthoGrad))
        if saddleDirMatch:
            filters.append(
                SaddleDirectionMatchFilter(siv, biTensor.siv, saddleDirMatch))

    filter = None
    if len(filters) == 1:
        filter = filters[0]
    elif filters:
        filter = FilterGroup(filters)

    maxima, flowlines = subpixelWatershedData(
        spws, siv, filter, mask,
        minSaddleDist = ssMinDist,
        perpendicularDistEpsilon = perpendicularDistEpsilon, maxStep = maxStep)
    
    spwsMap = subpixelWatershedMapFromData(
        maxima, flowlines, boundaryIndicator.size(),
        borderConnectionDist = borderConnectionDist,
        ssStepDist = ssStepDist, ssMinDist = ssMinDist,
        performBorderClosing = performBorderClosing,
        performEdgeSplits = performEdgeSplits,
        wsStatsSpline = initWSStats and siv or None,
        minima = initWSStats and spws.minima() or None,
        cleanup = cleanup)

    return spwsMap

def addFlowLinesToMap(edges, map, imageSize = None,
                      minSaddleBorderDist = 1e-4,
                      minMaxBorderDist = 1e-2):
    """addFlowLinesToMap(edges, map, imageSize = None)

    This function expects `edges` to be a list of
    SubPixelWatersheds.edge() return values and adds edges for each
    flowline.  The labels of the resulting edges will be valid indices
    into the `edges` sequence for the corresponding source data.

    It contains some special handling of flowlines:

    * Some artifact flowlines will be recognized and not added:

      * Self-loops with area zero are not added.
    
      * Flowlines whose saddles are within `minSaddleBorderDist` of the
        image border are presumed to be parallel to the border and are
        not added at all.

      * Flowlines that are *entirely* within `minMaxBorderDist` of the
        image border are handled the same (i.e. also discarded).
        
    * Flowlines that partially run outside the image are clipped (at a
      `minSaddleBorderDist` border within the image range), and only
      the part containing the saddle point will be added to the map.
      Additionally, the given `edges` sequence is modified in-place(!)
      in order to allow later WatershedStatistics to know about the
      modified saddle indices.

    * Additional Nodes may be added for flowlines that did not end in
      a maximum - extra care is taken not to insert multiple Nodes at
      (nearly) the same position.

    Returns list of edgeTuples that were not added to the GeoMap.
    This includes clipped parts of flowlines, with node labels
    etc. set to None."""
    
    # Node 0 conflicts with our special handling of 0 values:
    if map.maxNodeLabel():
        assert not map.node(0), \
               "addFlowLinesToMap: Node with label zero should not exist!"

    if edges[0]:
        raise ValueError("Cannot add edge with label 0 to map")

    # we don't want to add edges that run along the border to the Map
    # - a border can be added by connectBorderNodes if desired, and
    # parallel, double border edges lead to unsortable edges:
    if imageSize:
        imageBox = geomap.BoundingBox(imageSize - (1,1))
        clipBox  = copy.copy(imageBox)
        innerBox = copy.copy(imageBox)

        clipBox.addBorder(-minSaddleBorderDist)
        innerBox.addBorder(-minMaxBorderDist)

    import polytools

    result = []
    for edgeLabel, edgeTuple in enumerate(edges):
        if not edgeTuple:
            continue

        startNodeLabel = edgeTuple[0]
        endNodeLabel = edgeTuple[1]
        points = geomap.Polygon(edgeTuple[2])
        #saddleIndex = edgeTuple[3]
        #saddleLabel = edgeTuple[4]

        if imageSize and not innerBox.contains(points.boundingBox()):
            saddleIndex = edgeTuple[3]
            # don't add edges that run along the border at all,
            # i.e. exclude edges whose saddle point is not within
            # the boundingBox
            if not clipBox.contains(points[saddleIndex]):
                result.append(edgeTuple)
                continue

            if (innerBox & points.boundingBox()).isEmpty():
                result.append(edgeTuple)
                continue

            for cp in polytools.clipPoly(points, clipBox):
                try:
                    newSaddleIndex = list(cp).index(points[saddleIndex])
                except ValueError:
                    result.append((None, None, cp, None, None))
                    continue
                #sys.stderr.write("WARNING: flowline (edge label %d) left image range, clipped from %d to %d points (saddle index %d->%d)\n" % (edgeLabel, len(edgeTuple[2]), len(cp), saddleIndex, newSaddleIndex))

                # correct edgeTuple - check which end nodes are still connected:
                if points[0] != cp[0]:
                    startNodeLabel = -1
                if points[-1] != cp[-1]:
                    endNodeLabel = -1
                points = cp
                edgeTuple = (
                    startNodeLabel, endNodeLabel, points, newSaddleIndex, ) + edgeTuple[4:]
                edges[edgeLabel] = edgeTuple # write back for statistics (see docstring)
                break

        # meaning of node labels:
        #   >0: index of maximum
        # =  0: step size below threshold
        # = -1: dito, at image border
        # = -2: unwanted edge parallel to border (CODE REMOVED)
        # = -3: too many steps
        # = -4: dito, at image border

        if startNodeLabel == -2 and endNodeLabel == -2:
            # FIXME: check that partialArea of closed edge is near zero
            result.append(edgeTuple)
            continue # unwanted edge parallel to border

        # be careful not to modify the original 'edges' passed:
        points = geomap.Polygon(edgeTuple[2])

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
                startNode = nearestNode
                if vigra.dot(diff, pos - points[1]) >= 0: # don't jump back
                    if diff.squaredMagnitude():
                        # include node position if not present
                        points.insert(0, nearestNode.position())
                else:
                    points[0] = nearestNode.position()
            else: # no suitable Node found -> add one
                startNode = map.addNode(pos)
        else:
            startNode = map.node(startNodeLabel)
            assert startNode, "invalid startNodeLabel!"

        # ..handle Edge end the same as the start:
        if endNodeLabel <= 0:
            pos = points[-1]
            nearestNode = map.nearestNode(pos, 0.25)
            if nearestNode:
                diff = nearestNode.position() - pos
                endNode = nearestNode
                if vigra.dot(diff, pos - points[-2]) >= 0:
                    if diff.squaredMagnitude():
                        points.append(nearestNode.position())
                else:
                    points[-1] = nearestNode.position()
            else:
                endNode = map.addNode(pos)
        else:
            endNode = map.node(endNodeLabel)
            assert endNode, "invalid endNodeLabel!"

        # we need to avoid creating short, unsortable
        # self-loops with area zero:
        if startNode == endNode and (
            startNodeLabel <= 0 or endNodeLabel <= 0 or
            len(points) == 2):
            result.append(edgeTuple)
            continue
        
        assert startNode and endNode
        map.addEdge(startNode, endNode, points, edgeLabel)

    return result

# --------------------------------------------------------------------
#                        map creation helpers
# --------------------------------------------------------------------

def mapFromEdges(edges, imageSize):
    """Quickly create a GeoMap from a set of edges, i.e. derive nodes
    from edge endpoints.  Returns a GeoMap that is not yet
    initialized, i.e. `GeoMap.edgesSorted()` will return False."""
    
    def getNode(map, position):
        node = map.nearestNode(position, 0.001)
        if not node:
            node = map.addNode(position)
        return node
    
    result = geomap.GeoMap([], [], imageSize)
    for edge in edges:
        sn = getNode(result, edge[0])
        en = getNode(result, edge[-1])
        result.addEdge(sn, en, edge)
    
    return result

def gridMap(gridSize = (10, 10), firstPos = vigra.Vector2(0.5, 0.5),
            dist = vigra.Vector2(1, 1), imageSize = None):
    """Create a GeoMap representing a rectangular grid."""
    
    xDist = vigra.Vector2(dist[0], 0)
    yDist = vigra.Vector2(0, dist[1])
    if not imageSize:
        imageSize = (int(math.ceil(gridSize[0] * dist[0] + 2*(firstPos[0]+0.5))),
                     int(math.ceil(gridSize[1] * dist[1] + 2*(firstPos[1]+0.5))))

    map = geomap.GeoMap(imageSize)

    def addEdge(n1, n2):
        map.addEdge(n1, n2, [n1.position(), n2.position()])

    prevRow = None
    for y in range(gridSize[1]+1):
        row = []
        for x in range(gridSize[0]+1):
            pos = firstPos + x * xDist + y * yDist
            row.append(map.addNode(pos))

        for x in range(gridSize[0]):
            addEdge(row[x], row[x+1])

        if prevRow:
            for x in range(gridSize[0]+1):
                addEdge(row[x], prevRow[x])

        prevRow = row

    return map

def clipMapEdgesAtBorder(map):
    """Clip all edges of the given map at its image border.  Return
    new, uninitialized GeoMap object as created by `mapFromEdges` and
    modified by `connectBorderNodes`."""
    import polytools
    imageBox = geomap.BoundingBox(map.imageSize() - (1,1))

    edges = sum((polytools.clipPoly(edge, imageBox)
                 for edge in map.edgeIter()), [])

    result = mapFromEdges(edges, map.imageSize())
    connectBorderNodes(result, 1e-8)
    #result.initializeMap()
    return result

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

    lastPoints = [vigra.Vector2(x1, y1)]
    for node in top:
        thisPoints = []
        if node.position()[1] > y1 + samePosEpsilon:
            thisPoints.append(vigra.Vector2(node.position()[0], y1))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(vigra.Vector2(x2, y1))
    for node in right:
        thisPoints = []
        if node.position()[0] < x2 - samePosEpsilon:
            thisPoints.append(vigra.Vector2(x2, node.position()[1]))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(vigra.Vector2(x2, y2))
    for node in bottom:
        thisPoints = []
        if node.position()[1] < y2 - samePosEpsilon:
            thisPoints.append(vigra.Vector2(node.position()[0], y2))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    lastPoints.append(vigra.Vector2(x1, y2))
    for node in left:
        thisPoints = []
        if node.position()[0] > x1 + samePosEpsilon:
            thisPoints.append(vigra.Vector2(x1, node.position()[1]))
        thisPoints.append(node.position())
        lastPoints.extend(thisPoints)
        borderEdges.append((lastPoints, node))
        thisPoints.reverse()
        lastPoints = thisPoints
        lastNode = node

    if not borderEdges:
        cornerNode = map.addNode(lastPoints[0])
        lastPoints.append(lastPoints[0]) # close loop
        map.addEdge(cornerNode, cornerNode, lastPoints) \
                        .setFlag(flag_constants.BORDER_PROTECTION)
        return

    lastPoints.extend(borderEdges[0][0])
    borderEdges[0] = (lastPoints, borderEdges[0][1])

    endNode = lastNode
    for points, node in borderEdges:
        startNode = endNode
        endNode = node
        if not map.edgesSorted():
            map.addEdge(startNode, endNode, points) \
                                   .setFlag(flag_constants.BORDER_PROTECTION)
        else: # ouch! map already sorted!
            start = startNode
            if not startNode.isIsolated():
                start = startNode.anchor()
            end = endNode
            if not endNode.isIsolated():
                end = endNode.anchor()
            map.addEdge(start, end, points) \
                               .setFlag(flag_constants.BORDER_PROTECTION)
    if map.edgesSorted():
        sys.stderr.write("*** BIG FAT WARNING! ***\n  connectBorderNodes() called on already sorted map - resorting necessary!\n")

def copyMapContents(sourceMap, destMap = None, edgeTransform = None):
    """Add all nodes and edges of sourceMap to the destMap.  The sigma
    order is preserved.  Thus, it is an error to copy the contents of
    a graph - not sourceMap.edgesSorted() - into a map with
    destMap.edgesSorted().

    You can pass an optional edgeTransform to map edges to other
    Polygons: edgeTransform may be a callable that is called with Edge
    instances and must return an object suitable as the geometry
    parameter of GeoMap.addEdge.  If None is returned, the edge is
    skipped.

    Returns a tuple (destMap, nodes, edges), where nodes and edges are
    lists that map source nodes/edges onto their new counterparts in
    destMap.

    If destMap is None, a new GeoMap of the same size is created and
    returned (but not initialized!)."""

    if destMap == None:
        destMap = sourceMap.__class__(sourceMap.imageSize())
        if sourceMap.edgesSorted():
            # mark as sorted (sigma order will be copied from sourceMap):
            destMap.sortEdgesDirectly()

    assert sourceMap.edgesSorted() or not destMap.edgesSorted(), \
           "refraining from inserting unsorted edges into sorted map!"
    
    nodes = [None] * sourceMap.maxNodeLabel()
    # for better edgeTransform support, nodes are added on-the-fly now
    
    edges = [None] * sourceMap.maxEdgeLabel()
    for edge in sourceMap.edgeIter():
        geometry = edge
        if edgeTransform:
            geometry = edgeTransform(geometry)
            if geometry is None:
                continue

        startNeighbor = nodes[edge.startNodeLabel()]
        if startNeighbor == None:
            startNeighbor = destMap.addNode(geometry[0])
            nodes[edge.startNodeLabel()] = startNeighbor
        elif not startNeighbor.isIsolated():
            if edgeTransform:
                assert geometry[0] == startNeighbor.position(), "copyMapContents: edgeTransform leads to inconsistent node positions!"
            
            neighbor = edge.dart()
            while neighbor.nextSigma().edgeLabel() > edge.label() \
                      or neighbor.label() == -edge.label() \
                      or (edges[neighbor.edgeLabel()] is None
                          and neighbor.edgeLabel() != edge.label()):
                # (last cond. because edgeTransform may throw away edges)
                pass
            startNeighbor = edges[neighbor.edgeLabel()].dart()
            if neighbor.label() < 0:
                startNeighbor.nextAlpha()

        endNeighbor = nodes[edge.endNodeLabel()]
        if endNeighbor == None:
            endNeighbor = destMap.addNode(geometry[-1])
            nodes[edge.endNodeLabel()] = endNeighbor
        elif not endNeighbor.isIsolated():
            if edgeTransform:
                assert geometry[-1] == endNeighbor.position(), "copyMapContents: edgeTransform leads to inconsistent node positions!"
            
            neighbor = edge.dart().nextAlpha()
            while neighbor.nextSigma().edgeLabel() > edge.label() \
                      or neighbor.label() == edge.label() \
                      or (edges[neighbor.edgeLabel()] is None
                          and neighbor.edgeLabel() != edge.label()):
                # (last cond. because edgeTransform may throw away edges)
                pass
            endNeighbor = edges[neighbor.edgeLabel()].dart()
            if neighbor.label() < 0:
                endNeighbor.nextAlpha()

        newEdge = destMap.addEdge(startNeighbor, endNeighbor, geometry)
        newEdge.setFlag(edge.flags()) # retain Edge flags
        edges[edge.label()] = newEdge

    # now add isolated nodes:
    for node in sourceMap.nodeIter():
        if not nodes[node.label()]:
            nodes[node.label()] = destMap.addNode(node.position())

    return destMap, nodes, edges

def copyMap(sourceMap, edgeTransform = None):
    """Copy a complete GeoMap.  This builds upon `copyMapContents`,
    but also initializes the resulting map and copies all Face flags.

    There are only two use cases for this:

    1) Compressing all cell labels after many edges or nodes have been
       removed.

    2) Transforming edges, i.e. passing an `edgeTransform` that
       changes (e.g. smoothes) each Edge's geometry.
    
    This function cannot be used with an edgeTransform that returns
    None for any edge, since there must be a bijection between source
    and target Faces (it might be possible to lift this requirement
    w.r.t. bridges in the future)."""

    assert sourceMap.mapInitialized(), \
           "use copyMapContents for copying node/edge data only"

    # store edge mapping 'em':
    result, _, em = copyMapContents(sourceMap, None, edgeTransform)
    result.initializeMap(sourceMap.hasLabelImage())

    # transfer Face flags using 'em':
    for face in sourceMap.faceIter():
        anchor = face.contour()
        newAnchor = em[abs(anchor.label())].dart()
        if anchor.label() < 0:
            newAnchor.nextAlpha()
        newAnchor.leftFace().setFlag(face.flags() & 0x0fffffff)

    return result

def detachMapStats(map, verbose = False):
    """Calls 'detachHooks' on all properties of the given map.
    (Use this to prevent leaking statistics when discarding maps.)
    For convenience, this function acts as a no-op when called with
    map = None."""

    if map is None:
        return
    
    for a in map.__dict__:
        o = getattr(map, a)
        if hasattr(o, "detachHooks"):
            if verbose:
                print "detaching hooks of", o
            o.detachHooks()

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
    vigra.inspectImage(labelImage, fl)
    if fl.errorCount:
        sys.stderr.write("labelImage contains %d pixels with unknown faces!\n" % (
            fl.errorCount, ))
        sys.stderr.write("  unknown face labels found: %s\n" % (fl.errorLabels, ))
        if fl.errorCount < 40:
            for p in vigra.meshIter(labelImage.size()):
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

def drawLabelImage(aMap, scale = 1, verbose = True):
    """Return freshly drawn label image.  Should be identical to
    aMap.labelImage() if aMap.hasLabelImage() and scale == 1.  scale
    may be used to draw a label image with a different resolution;
    e.g. scale = 2 means that the returned image is (2w x 2h)
    where w,h = aMap.imageSize()."""

    total = aMap.faceCount - 1
    done = 1
    result = vigra.GrayImage(aMap.imageSize()*scale)
    holes = list(aMap.face(0).holeContours())
    # shift sampling points from middle of pixel (0.5, 0.5)
    # to middle of new, scaled pixel (scale/2, scale/2):
    offset = vigra.Vector2(scale / 2.0 - 0.5, scale / 2.0 - 0.5)
    for contour in holes:
        for hole in holeComponent(contour):
            if verbose and done % 23 == 0:
                sys.stdout.write("\r[%d%%] Face %d/%d" % (
                    done*100/total, done, total))
                sys.stdout.flush()
            poly = geomap.Polygon(
                geomap.contourPoly(hole.contour()) * scale + offset)
            sl = geomap.scanPoly(poly, result.height())
            geomap.fillScannedPoly(sl, result, hole.label())
            holes.extend(hole.holeContours())
            done += 1
    return result

def checkLabelConsistencyThoroughly(aMap):
    """Call `checkLabelConsistency` and additionally check that the
    labels are at the correct position (using `drawLabelImage`)."""
    
    class AssertRightLabel(object):
        def __init__(self):
            self.correct = True
        
        def __call__(self, label, shouldBe):
            if label >= 0:
                self.correct = self.correct and (label == shouldBe)

    return checkLabelConsistency(aMap) and \
           vigra.inspectImage(aMap.labelImage(),
                              drawLabelImage(aMap, verbose = False),
                              AssertRightLabel()).correct

def checkCachedPropertyConsistency(aMap):
    """Check whether edge and face bounds and areas are valid."""
    result = True
    realPolys = mapValidEdges(lambda edge: geomap.Polygon(list(edge)), aMap)
    for edge in aMap.edgeIter():
        poly = realPolys[edge.label()]
        if abs(poly.length() - edge.length()) > 1e-6:
            sys.stderr.write("Edge %d has cached length %s instead of %s!\n" % (
                edge.label(), edge.length(), poly.length()))
            result = False
        if poly.boundingBox() != edge.boundingBox():
            sys.stderr.write("Edge %d has cached boundingBox %s instead of %s!\n" % (
                edge.label(), edge.boundingBox(), poly.boundingBox()))
            result = False
        if abs(poly.partialArea() - edge.partialArea()) > 1e-5:
            sys.stderr.write("Edge %d has cached partialArea %s instead of %s!\n" % (
                edge.label(), edge.partialArea(), poly.partialArea()))
            result = False
    for face in aMap.faceIter():
        bbox = geomap.BoundingBox()
        area = 0.0
        for dart in face.contour().phiOrbit():
            edge = realPolys[dart.edgeLabel()]
            bbox |= edge.boundingBox()
            if dart.edge().isBridge():
                continue
            if dart.label() > 0:
                area += edge.partialArea()
            else:
                area -= edge.partialArea()
        
        if face.label() and bbox != face.boundingBox():
            sys.stderr.write("Face %d has cached bounding box %s instead of %s!\n" % (
                face.label(), face.boundingBox(), bbox))
            result = False
        if abs(area - face.area()) > 1e-5:
            sys.stderr.write("Face %d has cached area %s instead of %s!\n" % (
                face.label(), face.area(), area))
            result = False
    return result

def checkAllConsistency(aMap):
    """Return whether aMap.checkConsistency(),
    `checkCachedPropertyConsistency(aMap)`, and
    `checkLabelConsistencyThoroughly(aMap)` succeed (the latter in turn
    calls `checkLabelConsistency(aMap)` and is skipped if not
    aMap.hasLabelImage())."""
    return aMap.checkConsistency() and \
           checkCachedPropertyConsistency(aMap) and \
           (not aMap.hasLabelImage() or checkLabelConsistencyThoroughly(aMap))

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

    # (the following attributes are no longer attached to the map)
    if hasattr(map, "deleted") and map.deleted:
        print "%d edges were deleted (e.g. at image border)." % (
            len(map.deleted), )
    if hasattr(map, "unsortable") and map.unsortable:
        print "%d unsortable groups of edges occured." % (
            len(map.unsortable), )
    if hasattr(map, "unembeddableContours") and map.unembeddableContours:
        print "%d outer contours could not be embedded into " \
              "their surrounding faces!" % (len(map.unembeddableContours), )

    thisMap = "for GeoMap @0x%8x>" % id(map)
    for key, value in map.__dict__.items():
        print "  %s: %s" % (key, str(value).replace(thisMap, "for this GeoMap>"))

def degree2Nodes(map):
    return [node for node in map.nodeIter() if node.hasDegree(2)]

# --------------------------------------------------------------------
#                    basic map cleanup operations
# --------------------------------------------------------------------

def removeCruft(map, what = 3, doChecks = False):
    """removeCruft(map, what = 3, doChecks = False)
    
    `what` is a bit-combination of
    1: for the removal of degree 0-nodes (default)
    2: removal of degree 2-nodes (default)
    4: removal of bridges
    8: removal of edges (i.e. all non-protected)

    If `doChecks` is True, consistency checks are performed after
    every operation.  As soon as that fails, removeCruft returns
    False.

    After normal operation, removeCruft returns the number of
    operations performed.

    Consider using the following specialized functions with more
    meaningful names instead:
    
    - `removeIsolatedNodes`
    - `mergeDegree2Nodes`
    - `removeBridges`
    - `removeUnProtectedEdges`"""

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
            return map.checkConsistency()

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
            if node.hasDegree(2) and \
                   (node.anchor().endNode() != node):
                if not result.perform(map.mergeEdges, node.anchor()):
                    return False

    if what & 1:
        for node in map.nodeIter():
            if node.isIsolated():
                if not result.perform(map.removeIsolatedNode, node):
                    return False

    print "removeCruft(): %d operations performed." % result.count
    return result.count

from geomap import \
     removeIsolatedNodes, mergeDegree2Nodes, removeBridges, removeEdges

def removeUnProtectedEdges(map):
    return removeEdges(map, [
        edge.label() for edge in map.edgeIter()
        if not edge.flag(flag_constants.ALL_PROTECTION)])

def mergeFaceToBestNeighbor(face, costMeasure):
    best = None
    for dart in contourDarts(face):
        if dart.edge().flag(flag_constants.ALL_PROTECTION):
            continue
        cost = costMeasure(dart)
        if best is None or cost < best[0]:
            best = cost, dart
    return mergeFacesCompletely(best[1])

def removeSmallRegions(
    map, minArea = None, minPixelArea = None, costMeasure = None):
    """Merge faces whose area is < minArea with any neighor.
    Alternatively, pixelArea() is compared with minPixelArea.
    
    If there is more than one neighbor sharing an unprotected edge
    with the face, it must be decided which neighbor to merge into.
    If no `costMeasure` has been passed for this decision, an
    arbitrary neighbor is chosen and a warning is issued.  (Otherwise,
    the `costMeasure` is called for one dart per neighbored face and
    the lowest value decides about the merged faces.)"""

    assert minArea or minPixelArea, "either give minArea or minPixelArea!"

    result = 0
    for face in map.faceIter(skipInfinite = True):
        if minArea and face.area() >= minArea:
            continue
        if minPixelArea and face.pixelArea() >= minPixelArea:
            continue

        # similar to neighborFaces, except for the protection check:
        neighbors = []
        seen = {face.label(): True}
        for dart in contourDarts(face):
            if dart.edge().flag(flag_constants.ALL_PROTECTION):
                continue
            fl = dart.rightFaceLabel()
            if not fl in seen:
                seen[fl] = True
                neighbors.append(dart)
        
        if len(neighbors) > 1:
            if not costMeasure:
                sys.stderr.write("WARNING: removeSmallRegions() makes random decision about neighbor\n  which %s is to be merged into!\n" % face)
            else:
                neighbors.sort(key = lambda dart: costMeasure(dart))

        for dart in neighbors:
            if mergeFacesCompletely(dart):
                result += 1
                break

    if result:
        # the survivors from above may still be too small..
        result += removeSmallRegions(map, minArea, minPixelArea, costMeasure)

    return result

# --------------------------------------------------------------------
#                       composed Euler operations
# --------------------------------------------------------------------

from geomap import mergeFacesCompletely

def findCommonDart(face1, face2):
    """Find a dart with leftFace() == face1 and rightFace() == face2."""
    for contour in face1.contours():
        for contourIt in contour.phiOrbit():
            if contourIt.rightFaceLabel() == face2.label():
                return contourIt

def mergeFacesByLabel(map, label1, label2, mergeDegree2Nodes = True):
    """mergeFacesByLabel(map, label1, label2, mergeDegree2Nodes = True)

    Similar to mergeFacesCompletely() (which is called to perform the
    action), but is parametrized with two face labels and finds a
    common dart of these labels.

    Returns the surviving face (or None if no common edge was found)."""
    
    face1 = map.face(label1)
    face2 = map.face(label2)
    assert face1, "mergeFacesByLabel: face with label1 = %d does not exist!" % (label1, )
    assert face2, "mergeFacesByLabel: face with label2 = %d does not exist!" % (label2, )
    dart = findCommonDart(face1, face2)
    return dart and mergeFacesCompletely(dart, mergeDegree2Nodes)

# --------------------------------------------------------------------

class History(list):
    """List of (operation name, dart label) pairs.  Used to manage a
    list of operations that happend on a GeoMap.  For example, you may
    do:

    >>> map1 = ...
    >>> map2 = copy.copy(map1)
    >>> history = LiveHistory(map1)
    >>> someExpensiveAnalysisMethod(map1) # modify map1 via Euler ops
    >>> history.replay(map2)"""
    
    @staticmethod
    def load(filename):
        return History(eval(file(filename).read()))
    
    def save(self, filename):
        import pprint
        f = file(filename, "w")
        pprint.pprint(self, stream = f)
        f.close()

    def __getslice__(self, b, e):
        return History(list.__getslice__(self, b, e))

    def commands(self, mapName = "map"):
        """Return string with python commands replaying this history.
        (Very useful for creating code for unittests or the like.)"""
        result = ""
        for opName, label in self:
            result += "%s.%s(%s.dart(%d))\n" % (mapName, opName,
                                                mapName, label)
        return result

    def replay(self, map, careful = False, verbose = True):
        result = 0
        
#         counter = 20
        for opName, param in self:
            if careful and not map.checkConsistency():# or not checkLabelConsistency(map)
                sys.stderr.write(
                    "History.replay: map inconsistent, will not replay anything!\n")
                return result

            if verbose:
                print "replaying %s(dart %d)" % (opName, param)
            op = getattr(map, opName)
            assert op(map.dart(param)) != None, \
                   "History.replay failed (history does not match GeoMap anymore?)"
            result += 1

#             counter -= 1
#             if not counter:
#                 QtCore.qApp.processEvents()
#                 counter = 20

        if careful:
            map.checkConsistency()
        return result

class LiveHistory(History):
    """`History` list which attaches to a GeoMap and updates itself
    via Euler operation callbacks."""

    __slots__ = ("_map", "_op")
    
    _mergeFaces = 'mergeFaces'
    _removeBridge = 'removeBridge'
    _mergeEdges = 'mergeEdges'

    def __init__(self, map):
        self._map = weakref.ref(map) # only needed for pickle support
        self._attachHooks()
    
    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeFacesCallbacks(self.preMergeFaces, self.confirm),
            self._map().addRemoveBridgeCallbacks(self.preRemoveBridge, self.confirm),
            self._map().addMergeEdgesCallbacks(self.preMergeEdges, self.confirm),
            )

    def __getstate__(self):
        return (self._map(), )

    def __setstate__(self, (map, )): # cf. __init__
        self._map = weakref.ref(map)
        self._attachHooks()

    def detachHooks(self):
        for cb in self._attachedHooks:
            cb.disconnect()
    
    def preMergeFaces(self, dart):
        self._op = (self._mergeFaces, dart.label())
        return True
    
    def preRemoveBridge(self, dart):
        self._op = (self._removeBridge, dart.label())
        return True
    
    def preMergeEdges(self, dart):
        self._op = (self._mergeEdges, dart.label())
        return True
    
    def confirm(self, *args):
        self.append(self._op)

# --------------------------------------------------------------------
#                      Automatic Region Merger
# --------------------------------------------------------------------

class AutomaticMethodBase(object):
    __slots__ = ["_step", "_queue"]

    def nextCost(self):
        """Returns the cost of the merge operation that would be
        performed next by `mergeStep()`.

        Note that `mergeStep` does not always perform a step (i.e. the
        edge at the front of the internal priority queue might have
        been removed already); see `ensureValidNext()`."""

        return self._queue.top()[1]

    def ensureValidNext(self):
        """Ensure that the next `mergeStep()` really performs a merge
        operation.  This also means that `nextCost()` /
        `nextEdgeLabel()` return the corresponding valid values."""
        
        q = self._queue
        map = self._map
        while q and not map.edge(q.top()[0]):
            q.pop()

    def nextLabel(self):
        """Returns the label of the dart that would be the parameter
        of the merge operation that would be performed next by
        `mergeStep()` (for SRG, return the label of the face being
        associated next).
        
        Note that `mergeStep` does not always perform a step (i.e. the
        edge at the front of the internal priority queue might have
        been removed already); see `ensureValidNext()`."""
        
        self._ensureValidNext()
        return self._queue.top()[0]

    def step(self):
        """Returns the number of steps performed so far."""
        return self._step

    def merge(self):
        oldStep = self._step
        while self._queue:
            self.mergeStep()
        return self._step - oldStep

    def mergeSteps(self, count):
        return self.mergeToStep(self._step + count)

    def mergeToStep(self, targetStep):
        oldStep = self._step
        while self._queue and self._step < targetStep:
            self.mergeStep()
        return self._step - oldStep

    def mergeToCost(self, maxCost):
        oldStep = self._step
        while self._queue and self.nextCost() <= maxCost:
            self.mergeStep()
        return self._step - oldStep

class AutomaticRegionMerger(AutomaticMethodBase):
    """Merges faces in order of increasing costs.  The given
    mergeCostMeasure is called on a dart for each edge to determine
    the cost of removing that edge / merging the adjacent regions.

    Internally, a DynamicCostQueue is used in order to always remove
    the edge with the lowest assigned cost.  After each operation, the
    costs of all edges around the surviving face are recalculated."""

    __slots__ = ["_map", "_mergeCostMeasure", "_costLog",
                 "_mergeOperation", "_updateNeighborHood"]

    def __init__(self, map, mergeCostMeasure, q = None,
                 completeMerge = True, updateNeighborHood = True):
        """Set `completeMerge` to False to turn the automatic region
        merger into a simpler automatic edge remover tool.

        Set `updateNeighborHood` to True (default) if the
        mergeCostMeasure depends on the face statistics and the
        internal cost queue should be updated accordingly after a
        merge, which possibly changed the statistics."""
        
        self._map = map
        self._mergeCostMeasure = mergeCostMeasure
        self._step = 0
        self._updateNeighborHood = updateNeighborHood

        if completeMerge:
            self._mergeOperation = mergeFacesCompletely
        else:
            self._mergeOperation = self._map.mergeFaces

        if q == None:
            # FIXME: (why) is the +1 needed?
            q = geomap.DynamicCostQueue(map.maxEdgeLabel()+1)
            for edge in map.edgeIter():
                if edge.flag(flag_constants.ALL_PROTECTION):
                    continue
                cost = mergeCostMeasure(edge.dart())
                if cost is not None:
                    q.insert(edge.label(), cost)
        
        self._queue = q
        self._costLog = None

    def mergeStep(self):
        """Fetch next edge from cost queue and remove it from the map.
        
        If the edge is protected or nonexistent, do nothing and return
        None (i.e. not every call results in a merge step!).  Else,
        return the surviving Face and increment the step counter
        (cf. step())."""
        
        edgeLabel, mergeCost = self._queue.pop()
        edge = self._map.edge(edgeLabel)
        if not edge or edge.flag(flag_constants.ALL_PROTECTION):
            return

#         assert not self._updateNeighborHood or \
#                mergeCost == self._mergeCostMeasure(edge.dart())

        d = edge.dart()
        if edge.isBridge():
            survivor = self._map.removeBridge(d)
        else:
            survivor = self._mergeOperation(d)
            if survivor and self._updateNeighborHood:
                q = self._queue
                mcm = self._mergeCostMeasure
                for dart in contourDarts(survivor):
                    if dart.edge().flag(flag_constants.ALL_PROTECTION):
                        continue
                    cost = mcm(dart)
                    if cost is not None:
                        q.setCost(dart.edgeLabel(), cost)

        if survivor:
            if self._costLog is not None:
                self._costLog.append(mergeCost)
            self._step += 1

        return survivor

def thresholdMergeCost(map, mergeCostMeasure, maxCost, costs = None, q = None):
    """thresholdMergeCost(map, mergeCostMeasure, maxCost, costs = None, q = None)

    Merges faces of the given map in order of increasing costs until
    no more merges are assigned a cost <= maxCost by the
    mergeCostMeasure.

    The mergeCostMeasure is called on a dart for each edge to
    determine the cost of removing that edge.  This is used to
    initialize a DynamicCostQueue, which is used internally in order
    to always remove the edge with the lowest cost assigned.  The
    function returns a pair of the number of operations performed and
    the DynamicCostQueue, and you may pass the latter as optional
    argument 'd' into a subsequent call of thresholdMergeCost (with
    the same map and an increased maxCost) in order to re-use it
    (mergeCostMeasure is not used then).

    If the optional argument costs is given, it should be a mutable
    sequence which is append()ed the costs of each performed
    operation."""
    
    arm = AutomaticRegionMerger(map, mergeCostMeasure, q)
    arm._costLog = costs

    return arm.mergeToCost(maxCost), arm._queue

# --------------------------------------------------------------------

def extractFaceClasses(map, classes = None, skipInfinite = True):
    """Return dictionary mapping each of the specified `classes` to
    the corresponding list of faces within `map`.

    `classes` may be a sequence of face flags that specifies the exact
    classes to be searched for.  The resulting dict will have these
    classes as keys.

    `classes` may also be an integer, which is interpreted as a flags
    *mask*, i.e. the set of classes will be the set of disjoint values
    after applying this mask to the map's faces' flags.

    If `classes` is None, all occuring non-internal flags are
    collected (i.e. it is the same as the mask FACE_NONINTERNAL =
    0x0fffffff).

    If `skipInfinite` is not set to false, the infinite face will not
    be considered."""

    if classes is None:
        classes = flag_constants.FACE_NONINTERNAL
    if isinstance(classes, int):
        classes = set(face.flag(classes)
                      for face in map.faceIter(skipInfinite))

    classes = list(classes)
    classes.sort(reverse = True)

    result = dict((flag, []) for flag in classes)

    for face in map.faceIter(skipInfinite):
        for flag in classes:
            if flag and face.flag(flag) == flag:
                result[flag].append(face)
                break

    return result

# --------------------------------------------------------------------

def classifyFacesFromLabelImage(map, labelImage):
    """Given a labelImage, returns the label of each Face within that
    image.  This assumes that each region contains only pixels of the
    same label; otherwise an assertion will be triggered."""
    
    import statistics
    faceStats = statistics.FaceColorStatistics(map, labelImage)
    faceStats.detachHooks()

    result = [None] * map.maxFaceLabel()
    for face in map.faceIter(skipInfinite = True):
        l = faceStats[face.label()]
        assert float(int(l)) == l, "each Face must be entirely within one region of the labelImage"
        result[face.label()] = int(l)

    return result

def applyFaceClassification(map, faceClasses, ignoreNone = True):
    """Removes all edges between faces with the same class/label.
    Uses `removeEdges()` internally.  `faceClasses` should be a
    mapping from face labels to labels/classes that can be compared
    for equality.  E.g. to faces face1 and face2 will be merged iff
    faceClasses[face1.label()] == faceClasses[face2.label()].

    If `ignoreNone` is set (default), faces classified with `None`
    values will be ignored (i.e. not merged with neighbors that are
    also assigned None).
    
    `classifyFacesFromLabelImage` creates a suitable sequence (but
    there are better, more direct ways)."""

    return removeEdges(map, [
        edge.label() for edge in map.edgeIter()
        if faceClasses[edge.leftFaceLabel()] == faceClasses[edge.rightFaceLabel()]
        and (not ignoreNone or faceClasses[edge.leftFaceLabel()] is not None)])

def extractContractionKernel(map):
    """Returns a face classification suitable for
    `applyFaceClassification`, which represents the transformation
    from the original segmentation into the given `map`.  Internally,
    this is based on the LabelLUT that is maintained for the
    labelImage, i.e. it does not work for maps without labelImage
    ATM."""

    result = [None] * map.maxFaceLabel()
    labelLUT = map.faceLabelLUT()
    for face in map.faceIter(skipInfinite = True):
        for mergedLabel in labelLUT.merged(face.label()):
            result[mergedLabel] = face.label()
    return result

# --------------------------------------------------------------------
#                   geometrical utility functions
# --------------------------------------------------------------------

from geomap import centroid, contourPoly

def pointInFace(face, level = 2):
    """Given a face, return a point for which face.contains(point)
    returns True.  One possible use for this may be to save a
    classification of GeoMap Faces such that it can be applied to
    another GeoMap with different face labels, or even with slightly
    changed boundary geometry."""
    
    sl = face.scanLines()
    for y in range(sl.startIndex(), sl.endIndex()):
        inside = 0
        x = 0
        for seg in sl[y]:
            if inside and x < seg.begin:
                return vigra.Vector2(x, y)
            x = seg.end
            inside += seg.direction

    result = geomap.centroid(geomap.contourPoly(face.contour()))
    if face.contains(result):
        return result

    import numpy

    def subsampledPoint(face, level = 2):
        result = []
        bbox = face.boundingBox()
        midPoint = (bbox.begin() + bbox.end())/2
        xRange = numpy.arange(bbox.begin()[0], bbox.end()[0], 1.0/level)
        for y in numpy.arange(bbox.begin()[1], bbox.end()[1], 1.0/level):
            for x in xRange:
                p = vigra.Vector2(x, y)
                if face.contains(p):
                    result.append(((p-midPoint).squaredMagnitude(), p))
        if not result:
            return _pointInHole(face, level+1)
        result.sort()
        return result[0][1]

    return subsampledPoint(face)

def facePixels(face): # TODO: add ROI parameter to restrict range?
    """Generator function that iterates over the pixels within a face
    (those assigned the face label by fillScannedPoly)."""
    scanLines = face.scanLines()
    for y in range(scanLines.startIndex(), scanLines.endIndex()):
        inside = 0
        x = 0
        scanLine = scanLines[y]
        for segment in scanLine:
            begin = segment.begin
            end = segment.end
            if inside > 0:
                while x < begin:
                    yield vigra.Point2D(x, y)
                    x += 1
            x = end
            inside += segment.direction

# --------------------------------------------------------------------
#                   topological utility functions
# --------------------------------------------------------------------

def contourDarts(face):
    """contourDarts(face)

    Generator function which iterates over all darts within all
    contours, i.e. substitues the following standard loop construct:

      for contour in face.contours():
          for dart in contour.phiOrbit():
              doSth(dart)

    Thus, you can now write:

      for dart in contourDarts(face):
          doSth(dart)"""

    for contour in face.contours():
        for dart in contour.phiOrbit():
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
    for contour in face.holeContours():
        for hole in holeComponent(contour):
            showHomotopyTree(hole, indentation + "  ")

def edgeAtBorder(edge):
    """Return True iff edge has the BORDER_PROTECTION flag set."""
    return edge.flag(flag_constants.BORDER_PROTECTION)

def nonBorderEdges(map):
    """Iterate over all edges that do not have the BORDER_PROTECTION
    flag."""
    
    for edge in map.edgeIter():
        if not edge.flag(flag_constants.BORDER_PROTECTION):
            yield edge

def nonBackgroundFaces(map):
    """Iterate over all faces that do not have the BACKGROUND_FACE (=
    OUTER_FACE) flag."""
    
    for face in map.faceIter():
        if not face.flag(flag_constants.BACKGROUND_FACE):
            yield face

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

def mapValidEdges(function, geomap, default = None, skipBorder = False):
    """Similar to map(function, geomap.edgeIter()), but preserves
    labels.  I.e., equivalent to [default] * geomap.maxEdgeLabel() and
    a successive filling in of function results for valid edges.

    Note that due to the analogy to the builtin `map`, the order of the
    arguments is different to most other functions within this module,
    which usually have the geomap as first argument."""
    
    result = [default] * geomap.maxEdgeLabel()
    for edge in \
            skipBorder and nonBorderEdges(geomap) or geomap.edgeIter():
        result[edge.label()] = function(edge)
    return result

def mapValidFaces(function, geomap, default = None):
    """Similar to map(function, geomap.faceIter()), but preserves
    labels.  I.e., equivalent to [default] * geomap.maxFaceLabel() and
    a successive filling in of function results for valid faces.

    Note that due to the analogy to the builtin `map`, the order of the
    arguments is different to most other functions within this module,
    which usually have the geomap as first argument."""
    
    result = [default] * geomap.maxFaceLabel()
    for face in geomap.faceIter():
        result[face.label()] = function(face)
    return result

def mapValidDarts(function, geomap, default = None):
    """Similar to `mapValidEdges`, but calls function for each Dart
    with positive label.  Especially useful to create edgeCosts for
    passing e.g. to `minimumSpanningTree`."""
    
    result = [default] * geomap.maxEdgeLabel()
    for edge in geomap.edgeIter():
        result[edge.label()] = function(edge.dart())
    return result

# --------------------------------------------------------------------

def cpuCount(fallback = 1):
    '''
    Determine the number of CPU-cores via multiprocessing.cpu_count(),
    falling back to the given value (default = single core).
    '''
    import multiprocessing
    try:
        return multiprocessing.cpu_count()
    except NotImplementedError:
        print "WARNING: Could not determine core count. Using fallback (%d)." % fallback
    return fallback

def chop(l, n):
    '''
    Chop list l into n chunks
    '''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def _mapValidFacesWorker(function, geoMap, faceLabels, resultQueue):
    '''
    Worker function targeted by a single process during multi-core processing of a
     geomap's faces.
    @param function: (inspecting) function to be applied to the faces
    @param geoMap: GeoMap to be processed
    @param faceLabels: Subset of the geoMaps face labels to be processed
    @param resultQueue: Queue shared by all processes working on the geoMap object
    '''
    for faceLabel in faceLabels:
        result = function(geoMap.face(faceLabel))
        resultQueue.put((faceLabel, result))

def mapValidFacesParallel(function, geoMap, default = None, workerCount = None):
    '''
    API-compatible variant of mapValidFaces that applies the given
    function to batches of faces in parallel (using the
    multiprocessing module) and collects the results.  By default, one
    worker process is spawned for each cpu core (cf. cpuCount()).
    '''
    from multiprocessing import Queue, Process

    if workerCount is None:
        workerCount = cpuCount()

    faceLabels = chop([face.label() for face in geoMap.faceIter()], workerCount)

    resultQueue = Queue()

    processes = [Process(target = _mapValidFacesWorker,
                         args = (function, geoMap, workerLabels, resultQueue))
                 for workerLabels in faceLabels]

    for process in processes:
        process.start()
    for process in processes:
        process.join()

    result = [default] * geoMap.maxFaceLabel()
    while not resultQueue.empty():
        faceLabel, value = resultQueue.get()
        result[faceLabel] = value

    return result

# --------------------------------------------------------------------
#                       Seeded Region Growing
# --------------------------------------------------------------------

from heapq import heappush, heappop

def seededRegionGrowingStatic(map, faceLabels, edgeCosts):
    """seededRegionGrowingStatic(map, faceLabels, edgeCosts)

    Given a set of seed labels != None, extend labels to neighbor
    faces via edges with costs != None.  Stop when no more such
    growing is possible and return the resulting faceLabels.

    The 'static' postfix indicates that the edgeCosts do not change
    during the process.  (Furthermore, the GeoMap is not changed.)"""

    if not hasattr(edgeCosts, "__getitem__"):
        edgeCosts = mapValidDarts(edgeCosts, map)

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
            continue # labeled in the meantime
        faceLabels[dart.rightFaceLabel()] = faceLabels[dart.leftFaceLabel()]
        for dart in contourDarts(dart.rightFace()):
            if not faceLabels[dart.rightFaceLabel()]:
                cost = edgeCosts[dart.edgeLabel()]
                if cost != None:
                    heappush(heap, (cost, dart.label()))
    
    return faceLabels

class StandardCostQueue(object):
    """Wrapper around the standard heappush/pop, implementing the
    DynamicCostQueue API."""
    
    def __init__(self):
        self.heap = []

    def insert(self, index, cost):
        heappush(self.heap, (cost, index))

    setCost = insert

    def empty(self):
        return not self.heap

    def __nonzero__(self):
        return bool(self.heap)

    def top(self):
        cost, index = self.heap[0]
        return index, cost

    def pop(self):
        cost, index = heappop(self.heap)
        return index, cost

class SeededRegionGrowing(AutomaticMethodBase):
    """Implements Seeded Region Growing [Adams, Bischof].  Expects seed
    faces to be marked with the SRG_SEED flag.  Modifies the map by
    merging all non-seed faces into the growing seed regions until
    only these are left (or no more merges are possible, e.g. due to
    edge protection)."""
    
    __slots__ = ["_map", "_mergeCostMeasure", "_costLog",
                 "_neighborSkipFlags"]
    
    def __init__(self, map, mergeCostMeasure,
                 dynamic = True, stupidInit = False):
        self._map = map
        self._mergeCostMeasure = mergeCostMeasure
        self._step = 0
        
        for face in map.faceIter():
            face.setFlag(flag_constants.SRG_BORDER, False)

        self._neighborSkipFlags = flag_constants.SRG_SEED
        if dynamic:
            self._queue = geomap.DynamicCostQueue(map.maxFaceLabel())
        else:
            self._queue = StandardCostQueue()
            if stupidInit:
                self._neighborSkipFlags |= flag_constants.SRG_BORDER

        for face in map.faceIter():
            if face.flag(flag_constants.SRG_SEED):
                self._addNeighborsToQueue(face)

        assert self._queue, "SeededRegionGrowing: No seeds found (mark Faces with SRG_SEED)!"

        if not dynamic and not stupidInit:
            self._neighborSkipFlags |= flag_constants.SRG_BORDER

        self._costLog = None

    def _addNeighborsToQueue(self, face):
        for dart in neighborFaces(face):
            neighbor = dart.rightFace()
            if not neighbor.flag(self._neighborSkipFlags):
                cost = self._mergeCostMeasure(dart)
                self._queue.setCost(neighbor.label(), cost)
                neighbor.setFlag(flag_constants.SRG_BORDER)

    def mergeStep(self):
        # fetch candidate Face from queue:
        while True:
            faceLabel, mergeCost = self._queue.pop()
            face = self._map.face(faceLabel)
            if face and not face.flag(flag_constants.SRG_SEED):
                break

        assert face, "no more candidate faces found"

        # look for neighbor with lowest merge cost:
        best = None
        for dart in neighborFaces(face):
            neighbor = dart.rightFace()
            if neighbor.flag(flag_constants.SRG_SEED):
                cost = self._mergeCostMeasure(dart)
                if not best or cost < best[0]:
                    best = (cost, dart)

        # grow region:
        survivor = mergeFacesCompletely(best[1])
        if survivor:
            survivor.setFlag(flag_constants.SRG_BORDER, False)
            survivor.setFlag(flag_constants.SRG_SEED)
            self._addNeighborsToQueue(survivor)
            if self._costLog is not None:
                self._costLog.append(mergeCost)
            self._step += 1

        return survivor

    growStep   = mergeStep
    grow       = AutomaticMethodBase.merge
    growSteps  = AutomaticMethodBase.mergeSteps
    growToStep = AutomaticMethodBase.mergeToStep
    growToCost = AutomaticMethodBase.mergeToCost

def seededRegionGrowing(map, mergeCostMeasure, dynamic = False, stupidInit = False):
    """seededRegionGrowing(map, mergeCostMeasure, dynamic = False, stupidInit = False)

    Implements Seeded Region Growing [Adams, Bischof].  Expects seed
    faces to be marked with the SRG_SEED flag.  Modifies the map by
    merging all non-seed faces into the growing seed regions until
    only these are left (or no more merges are possible, e.g. due to
    edge protection).

    The optional parameter `dynamic` decides whether the costs of each
    possible merge shall be updated in each step.

    If `stupidInit` is set, the initialization is done exactly as in the
    original paper - the difference is that a face adjacent to two
    seeds will be associated the cost of merging it into the neighbor
    with the lower label (processing order) instead of the minimum of
    all costs.  Note that this streamlines the algorithm, but is
    dangerous when using superpixels, which makes the mentioned case
    happen regularly (instead of being an exception as in the original
    pixel-based SRG application).

    Equivalent to::

      srg = SeededRegionGrowing(map, mergeCostMeasure, dynamic, stupidInit)
      srg.grow()"""

    srg = SeededRegionGrowing(map, mergeCostMeasure, dynamic, stupidInit)
    srg.grow()

# --------------------------------------------------------------------
#                          MST / waterfall
# --------------------------------------------------------------------

def minimumSpanningTree(map, edgeCosts):
    """minimumSpanningTree(map, edgeCosts)

    Given a cost associated with each edge of the map, this function
    finds the minimum spanning tree of the map's region adjacency
    graph (None in edgeCosts is allowed and is handled as if the
    corresponding edge was missing).  The result is a modified copy of
    the edgeCosts list, with all non-MST-edges set to None.  This can
    be used for the waterfall algorithm by Meyer and Beucher, see
    waterfall().

    The complexity is O(edgeCount*faceCount), but the performance
    could be improved a little."""

    if not hasattr(edgeCosts, "__getitem__"):
        edgeCosts = mapValidDarts(edgeCosts, map)

    print "- initializing priority queue for MST..."
    heap = []
    for edge in map.edgeIter():
        edgeLabel = edge.label()
        cost = edgeCosts[edgeLabel]
        if cost != None:
            heappush(heap, (cost, edgeLabel))

    print "- building MST..."
    faceLabels = geomap.LabelLUT(map.maxFaceLabel())
    result = list(edgeCosts)
    while heap:
        _, edgeLabel = heappop(heap)
        edge = map.edge(edgeLabel)
        lfl = faceLabels[edge.leftFaceLabel()]
        rfl = faceLabels[edge.rightFaceLabel()]
        if lfl != rfl:
            faceLabels.relabel(lfl, rfl)
        else:
            result[edgeLabel] = None

    return result

def regionalMinima(map, mst):
    """Return an array with face labels, where only 'regional minima'
    of the given MST are labeled (i.e. the rest is labeled None).
    See the article 'Fast Implementation of Waterfall Based on Graphs'
    by Marcotegui and Beucher."""

    faceLabels = [None] * map.maxFaceLabel()
    for edge in map.edgeIter():
        edgeCost = mst[edge.label()]
        if edgeCost == None:
            continue
        isMinimum = True
        # FIXME: this finds only local minima - in theory extended
        # local minima should be labeled, too
        for start in edge.dart().alphaOrbit():
            # note that contourDarts is sigmaOrbit in the dual graph:
            for dart in contourDarts(start.leftFace()):
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
    applyFaceClassification(map, faceLabels)
    print "  total waterfall() time: %ss." % (time.clock() - c, )

def dualMap(map, edgeLabels = None, nodePositions = None, midPoints = None,
            onDemandNodeHook = None, initializeMap = True):
    """Compute (a subset of) the dual of a GeoMap.
    `edgeLabels` determines which edges appear in the result.
    If None (default), the complete dual map is returned.
    The labels of the result's nodes will correspond to the source
    map's face labels; edge labels are retained, too.

    `nodePositions` determines the position of each node; be default
    (nodePositions == None), that the nodes are located at the
    centroids of the faces as a simple heuristic.

    `midPoints` determines whether connecting edges should be straight
    connections (i.e. 2 points, midPoints = False) or have an extra
    point in the middle of their dual edges.  Default (midPoints ==
    None): Add midpoints for edges which otherwise do not intersect
    their duals.

    `onDemandNodeHook` can be used to create nodes on demand,
    i.e. when one face should result in several nodes.  (This is
    e.g. used for the infinite face when computing Voronoi maps from
    Delaunay maps.)"""

    result = geomap.GeoMap(map.imageSize())

    if nodePositions is None:
        nodePositions = mapValidFaces(lambda face: geomap.centroid(
            geomap.contourPoly(face.contour())), map)

    nodes = [None] * map.maxFaceLabel()
    for face in map.faceIter(skipInfinite = True):
        l = face.label()
        if nodePositions[l] is not None:
            nodes[l] = result.addNode(nodePositions[l], label = l)

    if edgeLabels is None:
        edgeLabels = map.edgeIter()
    
    for edge in sorted(edgeLabels):
        if not hasattr(edge, "label"):
            edge = map.edge(edge)
        assert not edge.isBridge(), "bridges not supported yet (would need to create self-loops around end node in dual map)"
        sn = nodes[edge.leftFaceLabel()]
        en = nodes[edge.rightFaceLabel()]
        if onDemandNodeHook and (not sn or not en):
            sn, en = onDemandNodeHook(map, edge, result, sn, en)
        if sn and en:
            poly = [sn.position(), en.position()]
            midPoint = midPoints
            if midPoints == None:
                # auto-detect if midPoint is necessary (no intersection):
                clipped = geomap.intersectLine(edge, poly[0], poly[1])
                midPoint = not len(clipped) or (
                    len(clipped) == 1 and len(clipped[0]) == len(edge))
            if midPoint:
                dp = geomap.DartPosition(edge.dart())
                dp.gotoArcLength(edge.length()/2)
                poly.insert(1, dp())
            result.addEdge(sn, en, poly, label = edge.label())
    
    removeIsolatedNodes(result)
    if initializeMap:
        result.initializeMap(False)
    return result

def mst2map(mst, map):
    """Visualize a minimumSpanningTree() result as a GeoMap.
    (See `dualMap()`.)"""
    return dualMap(map, [
        edgeLabel for edgeLabel, cost in enumerate(mst) if cost is not None])

# --------------------------------------------------------------------
#                               tests
# --------------------------------------------------------------------

# if __name__ == '__main__':
    
