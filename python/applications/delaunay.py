# -*- coding: iso-8859-1 -*-

"""(Constrained) Delaunay Triangulation and Chordal Axis Transform

Here's a complete walk-through how to compute the CAT skeleton:

>>> ms = levelcontours.marchingSquares(..., 128)

Note that marchingSquares marks holes/background faces with
OUTER_FACE; depending on your application, you might have to mark
faces whose skeleton you don't need.

Then, you can compute the CDT of all remaining foreground faces:

>>> cdt = delaunay.cdtFromPSLG(ms)

You might want to keep the original CDT map (e.g. needed for
`pruneBySubtendedLength`):

>>> r_cdt = copy.deepcopy(cdt)

The rectified CAT reduces oscillations and leads to much nicer
skeletons:

>>> delaunay.removeWeakChords(r_cdt)
>>> r_cat = delaunay.catMap(r_cdt, rectified = True)

Finally, it makes sense to prune the skeleton, i.e. remove small
barbs:

>>> delaunay.pruneBySubtendedLength(r_cat, cdt, minLength = 10)
"""

__version__ = \
    "$Id$" \
    .split(" ")[2:-2]

import math, sys, hourglass, numpy
#from map import GeoMap, contourPoly
from hourglass import GeoMap, contourPoly
from vigra import Vector2, Vector, dot
from flag_constants import \
     START_NODE_ADDED, END_NODE_ADDED, CONTOUR_SEGMENT, OUTER_FACE, \
     IS_BARB, WEAK_CHORD

try:
    import triangle
except ImportError:
    sys.stderr.write("WARNING: 'triangle' module not found. CDT will be (not very well) faked.\n")
    triangle = None

def _delaunayMapFromData(nodePositions, edgeData, imageSize, sigmaOrbits = None):
    if sigmaOrbits:
        sys.stderr.write(
            "TODO: re-use Delaunay sigma orbits instead of re-sorting!\n")
    edges = [startEnd and
             (startEnd[0], startEnd[1],
              [nodePositions[startEnd[0]], nodePositions[startEnd[1]]])
             for startEnd in edgeData]
    result = GeoMap(nodePositions, edges, imageSize)
    result.initializeMap(initLabelImage = False)
    return result

def _pointInHole(polygon, level = 2):
    sl = hourglass.scanPoly(polygon)
    for y in range(sl.startIndex(), sl.endIndex()):
        inside = 0
        x = 0
        for seg in sl[y]:
            if inside and x < seg.begin:
                return Vector2(x, y)
            x = seg.end
            inside += seg.direction

    result = hourglass.centroid(polygon)
    if polygon.contains(result):
        return result
    
    result = []
    bbox = polygon.boundingBox()
    midPoint = (bbox.begin() + bbox.end())/2
    xRange = numpy.arange(bbox.begin()[0], bbox.end()[0], 1.0/level)
    for y in numpy.arange(bbox.begin()[1], bbox.end()[1], 1.0/level):
        for x in xRange:
            p = Vector2(x, y)
            if polygon.contains(p):
                result.append(((p-midPoint).squaredMagnitude(), p))
    if not result:
        return _pointInHole(polygon, level+1)
    result.sort()
    return result[0][1]

def constrainedDelaunayMap(polygons, imageSize, extraPoints = [],
                           onlyInner = True):
    """constrainedDelaunayMap(polygons, imageSize, extraPoints, onlyInner) -> GeoMap

    Returns a new GeoMap containing a Constrained Delaunay
    Triangulation of all points of the polygons plus the extraPoints
    if given.  The segments of the polygons will be constrained
    segments of the CDT.

    If the optional `onlyInner` parameter is True (default), then all
    edge segments in the 'outside' will be removed.  (This assumes
    that all closed polygons with a negative partialArea() are
    holes.)

    **NOTE**: You probably want to use `cdtFromPSLG` or `faceCDTMap`
    instead, both of which have a much simpler API."""

    assert triangle, """For correct CDT, you need to compile the
    triangle module (vigra/experiments/triangle).  You might want to
    try the home-made fakedConstrainedDelaunayMap instead, which will
    give a correct result if possible, but may just throw an exception
    if it has to give up or even go into an infinite loop."""

    points = list(extraPoints)
    segments = []
    holes = []
    for polygon in polygons:
        l = len(points)
        partPoints = list(polygon)
        partSegments = [(l+i, l+i+1) for i in range(len(partPoints)-1)]
        if polygon[-1] == polygon[0]:
            del partPoints[-1]
            if partPoints[0] in points:
                for i, (s, e) in enumerate(partSegments):
                    partSegments[i] = (s-1, e-1) # partPoints[0] will be deleted
                partSegments[0] = (points.index(partPoints[0]), partSegments[0][1])
                del partPoints[0]
            partSegments[-1] = (partSegments[-1][0], partSegments[0][0])
            if onlyInner and polygon.partialArea() < 0:
                holes.append(_pointInHole(polygon))
        else:
            if partPoints[-1] in points:
                partSegments[-1] = (partSegments[-1][0], points.index(partPoints[-1]))
                del partPoints[-1]
            if partPoints[0] in points:
                for i, (s, e) in enumerate(partSegments):
                    partSegments[i] = (s-1, e-1) # partPoints[0] will be deleted
                partSegments[0] = (points.index(partPoints[0]), partSegments[0][1])
                del partPoints[0]
        points.extend(partPoints)
        segments.extend(partSegments)

    print "- performing Constrained Delaunay Triangulation..."
    print "  (%d points, %s segments, %d holes)" % (
        len(points), len(segments), len(holes))
    nodePositions, edgeData = triangle.constrainedDelaunay(
        points, segments, onlyInner, holes)

    print "- storing result in a GeoMap..."
    result = _delaunayMapFromData(nodePositions, edgeData, imageSize)

    for edge in result.edgeIter():
        if edgeData[edge.label()][2]:
            edge.setFlag(CONTOUR_SEGMENT)

    result.face(0).setFlag(OUTER_FACE)
    for holePoint in holes:
        result.faceAt(holePoint).setFlag(OUTER_FACE)

    return result

def delaunayMap(points, imageSize = (0, 0)):
    """Return a GeoMap containing a Delaunay Triangulation of the given
    points."""
    
    if triangle:
        nodePositions, edges = triangle.delaunay(points)
        sigma = None
    else:
        nodePositions, edges, sigma = hourglass.delaunay(points)
    return _delaunayMapFromData(nodePositions, edges, imageSize,
                                sigma)

def cdtFromPSLG(pslg, onlyInner = False):
    """cdtFromPSLG(pslg, onlyInner) -> GeoMap

    Return a CDT for the given planar straight line graph.  `pslg`
    should be a GeoMap whose node positions are simple input points
    and edges define constraint segments.

    FIXME: ATM, this function may crash if the input contains
    duplicate points (e.g. within the same edge).  This is due to a
    limitation of the triangle module, i.e. J. Shewchuck's code.
    """
#     return constrainedDelaunayMap(
#         list(pslg.edgeIter()), pslg.imageSize(),
#         [node.position() for node in pslg.nodeIter() if node.isIsolated()],
#         onlyInner = False)

    points = []
    segments = []
    holes = []

    nodes = [None] * pslg.maxNodeLabel()
    for node in pslg.nodeIter():
        nodes[node.label()] = len(points)
        points.append(node.position())

    for edge in pslg.edgeIter():
        l = len(points)
        edgeSegments = [(l+i-1, l+i) for i in range(len(edge)-1)]
        edgeSegments[0] = (nodes[edge.startNodeLabel()], edgeSegments[0][1])
        edgeSegments[-1] = (edgeSegments[-1][0], nodes[edge.endNodeLabel()])
        points.extend(edge[1:-1])
#         for i in range(1, len(edge)-1):
#             assert edge[i] not in points, "%s[%d] = ..%s,%s,%s,.. is already present!" % (
#                 edge, i, edge[i-1], edge[i], edge[i+1])
#             points.append(edge[i])
        segments.extend(edgeSegments)

    for face in pslg.faceIter(skipInfinite = True):
        if face.flag(OUTER_FACE):
            # FIXME: re-use available scanlines for face:
            holes.append(_pointInHole(contourPoly(face.contour())))

    print "- performing Constrained Delaunay Triangulation..."
    print "  (%d points, %s segments, %d holes)" % (
        len(points), len(segments), len(holes))
    nodePositions, edgeData = triangle.constrainedDelaunay(
        points, segments, onlyInner, holes)

    print "- storing result in a GeoMap..."
    result = _delaunayMapFromData(nodePositions, edgeData, pslg.imageSize())

    for edge in result.edgeIter():
        if edgeData[edge.label()][2]:
            edge.setFlag(CONTOUR_SEGMENT)

    result.face(0).setFlag(OUTER_FACE)
    for holePoint in holes:
        result.faceAt(holePoint).setFlag(OUTER_FACE)

    return result


def fakedConstrainedDelaunayMap(polygons, imageSize, extraPoints = [],
                                onlyInner = True):

    """See constrainedDelaunayMap, this calculates a DT and throws
    away outer edges retroactively.  This may fail when the DT does
    not contain all constrained segments, which is checked in this
    function and leads to an AssertionError."""

    points = []
    jumpPoints = [0]
    for polygon in polygons:
        points.extend(list(polygon))
        del points[-1]
        jumpPoints.append(len(points))
    
    print "- performing Delaunay Triangulation (%d points)..." % len(points)
    nodePositions, edges, sigma = hourglass.delaunay(points)
    
    print "- storing result in a GeoMap..."
    result = _delaunayMapFromData(nodePositions, edges, imageSize, sigma)

    print "- ex-post marking of contour edges for faked CDT..."
    print "  (keep your fingers crossed that no segment is missing!)"

    edgeSourceDarts = [None] * result.maxEdgeLabel()

    i = 0
    while i < len(jumpPoints) - 1:
        contourStartLabel = jumpPoints[i] + 1 # +1 for conversion of
        contourEndLabel = jumpPoints[i+1] + 1 #    point indices -> node labels
        dart = result.node(contourStartLabel).anchor()
        for nodeLabel in range(contourStartLabel+1, contourEndLabel):
            j = dart.startNode().degree() + 2
            while dart.endNodeLabel() != nodeLabel and j > 0:
                dart.nextSigma()
                j -= 1
            assert j > 0, """Original contour fragment missing in Delauny map!
            (This is a problem of the fakedConstrainedDelaunayMap, try compiling
            the triangle module and using the real constrainedDelaunayMap instead.)"""
            dart.edge().setFlag(CONTOUR_SEGMENT)
            edgeSourceDarts[dart.edgeLabel()] = dart.label()
            dart.nextAlpha()
        j = dart.startNode().degree() + 1
        while dart.endNodeLabel() != contourStartLabel and j > 0:
            dart.nextSigma()
            j -= 1
        assert j > 0, "Original contour fragment missing in Delauny map!"
        dart.edge().setFlag(CONTOUR_SEGMENT)
        edgeSourceDarts[dart.edgeLabel()] = dart.label()
        i += 1

    if onlyInner:
        sys.stdout.write("- reducing delaunay triangulation to inner part...")
        outerFace = [False] * result.maxFaceLabel()
        for edge in result.edgeIter():
            if edge.flag(CONTOUR_SEGMENT):
                dart = result.dart(edgeSourceDarts[edge.label()])
                assert dart.edgeLabel() == edge.label(), str(edge)
                while True:
                    mergeDart = dart.clone().prevSigma()
                    if mergeDart.edge().flag(CONTOUR_SEGMENT):
                        break
                    outerFace[result.mergeFaces(mergeDart).label()] = True
                    #sys.stdout.write(".")
                    #sys.stdout.flush()
        sys.stdout.write("\n")
    
    return result

def faceCDTMap(face, imageSize,
               simplifyEpsilon = None,
               resample = None,
               onlyInner = True):
    """USAGE: dlm = faceCDTMap(face, mapSize)

    `face` should be a GeoMap.Face object, and all its contours will
    be extracted.  `mapSize` is used to initialize the GeoMap with the
    resulting edges.
      
    Optional keyword parameters:

    simplifyEpsilon
      If given, each contour polygon is simplified by calling
      simplifyPolygon with this epsilon as parameter (default None ->
      don't use simplifyPolygon).

    markContour
      If True(default), the point list is expected to be a sorted
      list, and edges between successive entries are marked as contour
      edges (an exception is raised if such a connection is
      missing). A `jumpPoints` list is used to mark multiple
      contours.

    onlyInner
      If True(default), all edges outside (left) of the marked contour
      are removed in a post-processing step.  (This has no effect if
      `markContour` is False.)"""

    polygons = [contourPoly(c) for c in face.contours()]
    if resample:
        polygons = [hourglass.resamplePolygon(p, resample) for p in polygons]
    if simplifyEpsilon != None:
        polygons = [hourglass.simplifyPolygon(p, simplifyEpsilon) for p in polygons]

    if triangle:
        result = constrainedDelaunayMap(polygons, imageSize)
    else:
        result = fakedConstrainedDelaunayMap(polygons, imageSize)

    return result

# --------------------------------------------------------------------

def _oppositeAngle(p1, p2, p3):
    """_oppositeAngle(p1, p2, p3)
    Given the triangle (p1, p2, p3), returns angle at p3 (in radians)"""
    a = p1 - p3
    b = p2 - p3
    c = dot(a, b) / \
        math.sqrt(a.squaredMagnitude()*b.squaredMagnitude())
    if abs(c) < 1:
        return math.acos(c)
    elif abs(c) > 0:
        return 0
    return math.pi

def chordStrength(edge):
    p1 = edge[0]
    p2 = edge[1]
    d = edge.dart()
    p3 = d.nextPhi()[1]
    assert d.nextPhi().endNodeLabel() == edge.startNodeLabel(), \
           "chordStrength: expects Faces to the left and right to be triangles!"
    d = edge.dart().nextAlpha()
    p4 = d.nextPhi()[1]
    assert d.nextPhi().endNodeLabel() == edge.endNodeLabel(), \
           "chordStrength: expects Faces to the left and right to be triangles!"
    return 1 - (_oppositeAngle(p1, p2, p3) + _oppositeAngle(p1, p2, p4)) / math.pi

def calculateChordStrengths(delaunayMap):
    """Calculate the chordStrength of each non-contour edge."""

    result = [None] * delaunayMap.maxEdgeLabel()
    for edge in delaunayMap.edgeIter():
        if edge.flag(CONTOUR_SEGMENT):
            continue
        result[edge.label()] = chordStrength(edge)
    return result

def chordStrengthProfile(delaunayMap, startDart = None, cs = None):
    """Calculate a list of (chordStrength, dartLabel) pairs for each
    dart within the original shape contour for which startDart is an
    anchor."""

    if not startDart:
        startDart = delaunayMap.node(1).anchor()
        while startDart.endNodeLabel() != 2:
            startDart.nextSigma()

    assert startDart.edge().flag(CONTOUR_SEGMENT) and \
           not startDart.leftFace().flag(OUTER_FACE), \
           "expecting startDart to be on the shape's contour," \
           "with the inner triangles on the left"
    
    if not cs:
        cs = calculateChordStrengths(delaunayMap)

    dart = startDart.clone()
    result = []
    while True:
        #print dart
        while not dart.nextSigma().edge().flag(CONTOUR_SEGMENT):
            #print "-", dart
            result.append((cs[dart.edgeLabel()], dart.label()))
        if dart.nextAlpha() == startDart:
            #print "AT BEGIN (%s)" % dart.label()
            break
    return result

def markWeakChords(delaunayMap):
    for edge in delaunayMap.edgeIter():
        if not edge.flag(CONTOUR_SEGMENT):
            edge.setFlag(WEAK_CHORD)
    cs = calculateChordStrengths(delaunayMap)
    for face in delaunayMap.faceIter():
        if not face.flag(OUTER_FACE):
            continue
        for anchor in face.contours():
            csp = chordStrengthProfile(delaunayMap, anchor.nextAlpha(), cs)
            csp.append(csp[0])
            for i in range(len(csp)-1):
                c1 = cmp(csp[i][0], csp[i-1][0])
                c2 = cmp(csp[i][0], csp[i+1][0])
                if c1 + c2 >= 1:
                    delaunayMap.edge(abs(csp[i][1])).setFlag(WEAK_CHORD, False)

def removeWeakChords(delaunayMap):
    result = 0
    markWeakChords(delaunayMap)
    for edge in delaunayMap.edgeIter():
        if edge.flag(WEAK_CHORD):
            delaunayMap.mergeFaces(edge.dart())
            result += 1
    return result

import qt
def chordColors(delaunayMap):
    result = [qt.Qt.red] * delaunayMap.maxEdgeLabel()
    for edge in delaunayMap.edgeIter():
        if edge.flag(CONTOUR_SEGMENT):
            result[edge.label()] = qt.Qt.cyan
        elif edge.flag(WEAK_CHORD):
            result[edge.label()] = qt.Qt.blue
    return result

# --------------------------------------------------------------------

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def baryCenter(*points):
    return sum(points) / len(points)

def oldJunctionNodePosition(p1, p2, p3,
                            joinMiddleThreshold = 1.61):
    "Old, deprecated junction node position defition (ignore this)"
    a = (p2-p1).magnitude()
    b = (p3-p2).magnitude()
    c = (p1-p3).magnitude()
    sv = Vector(a, b, c)
    sv /= sv.magnitude()
    if dot(sv, Vector(1, 1, 1)) < joinMiddleThreshold:
        s = [a,b,c]
        s.sort()
        if s[2] == a:
            nodePos = (p2+p1)/2
        elif s[2] == b:
            nodePos = (p3+p2)/2
        else:
            nodePos = (p1+p3)/2

def _middleOfLongestSide(p1, p2, p3):
    s = [((p2-p1).squaredMagnitude(), (p2+p1)),
         ((p3-p2).squaredMagnitude(), (p3+p2)),
         ((p1-p3).squaredMagnitude(), (p1+p3))]
    s.sort()
    return s[2][1] / 2

def junctionNodePosition(chords):
    """Defition of node position for junction triangles from original
    CAT article (works only for triangles, uses circumcenter if
    possible).  Use this if possible."""
    p1 = chords[0][0]
    p2 = chords[1][0]
    p3 = chords[2][0]
    circumCenter, circumRadius = circumCircle(p1, p2, p3)
    if chords[0].leftFace().contains(circumCenter):
        return circumCenter
    return _middleOfLongestSide(p1, p2, p3)

def rectifiedJunctionNodePosition(chords):
    """Defition of node position for junction triangles from rectified
    CAT article.  Not as good as the original definition for
    triangles, but works for other polygons, too."""
    totalWeights = 0.0
    result = Vector2(0, 0)
    for chord in chords:
        weight = chord.edge().length()
        result += middlePoint(chord) * weight
        totalWeights += weight
    return result / totalWeights

def catMap(delaunayMap,
           rectified = True,
           includeTerminalPositions = False):
    """catMap(delaunayMap,
           rectified = True,
           includeTerminalPositions = False)

    Extract a CAT (chordal axis transform) from a GeoMap object
    containing a Delaunay Triangulation.  Implements the rectified
    version (from L.Prasad's DGCI'05 article), which also works for
    GeoMaps with non-triangular Faces.

    Expects outer faces of the delaunayMap to be marked with the
    OUTER_FACE flag (in order to skip these faces and work only on the
    inner parts of the shapes), and edges which belong to the
    constrained segments of the CDT to be marked with the
    CONTOUR_SEGMENT flag (i.e. an edge is a chord iff
    'not edge.flag(CONTOUR_SEGMENT)')

    Usually, you will do sth. along:

    >>> del1 = delaunay.faceCDTMap(map1.face(173), map1.imageSize())
    >>> cat1 = delaunay.catMap(del1, rectified = False)

    The optional parameter `rectified` has been set to False here to
    use the original definition of junction node positions (using the
    circumcenter if possible) which works only for triangulations (and
    is therefore disabled by default).

    For triangulations, it is also possible to include the opposite
    vertex of the terminal triangles into the skeleton by setting
    `includeTerminalPositions` to True.

    If you want to use the rectified CAT (with weak chords being
    suppressed), you would use removeWeakChords before calling catMap:
    
    >>> del2 = copy.copy(del1) # if you want to keep the original
    >>> removeWeakChords(del2)
    >>> rcat2 = delaunay.catMap(del2)

    The resulting map will have additional attributes that relate the
    skeleton with the original boundary:

    subtendedLengths:
      This is a list mapping edge labels of the skeleton to the
      accumulated lengths of all contour segments within the contours
      of the sleeve- and terminal-faces that comprise the sleeve
      represented by that skeleton edge.

    shapeWidths:
      This is a list mapping edge labels of the skeleton to lists of
      shape widths at each support point of the polygon (may be None
      at the ends).  In theory, this allows to reconstruct the
      original boundary and thus makes the CAT fully reversible (see
      Prasad's papers).

    nodeChordLabels:
      This is a list mapping node labels of the skeleton to lists of
      (sleeveLabel, chordLabel) pairs, where sleeveLabel is the label
      of a skeleton edge and chordLabel is the dart label of the first
      chord from the `delaunayMap` within that sleeve.  (This is
      needed for later corrections of the junction node positions
      after pruning.)"""

    if rectified:
        junctionPos = rectifiedJunctionNodePosition
        assert not includeTerminalPositions, \
               "includeTerminalPositions is not supported for the rectified CAT!"
    else:
        junctionPos = junctionNodePosition

    result = GeoMap(delaunayMap.imageSize())
    result.subtendedLengths = [0] # unused Edge 0
    result.shapeWidths = [None]
    result.nodeChordLabels = []

    faceChords = [None] * delaunayMap.maxFaceLabel()
    boundaryLength = [None] * delaunayMap.maxFaceLabel()
    faceType = [None] * delaunayMap.maxFaceLabel()
    nodeLabel = [None] * delaunayMap.maxFaceLabel()

    for face in delaunayMap.faceIter(skipInfinite = True):
        if face.flag(OUTER_FACE):
            continue

        assert face.holeCount() == 0

        chords = []
        bl = 0.0
        for dart in face.contour().phiOrbit():
            if not dart.edge().flag(CONTOUR_SEGMENT):
                chords.append(dart)
            else:
                bl += dart.edge().length()
        
        boundaryLength[face.label()] = bl
        faceChords[face.label()] = chords
        
        # classify into terminal, sleeve, and junction triangles:
        if len(chords) < 2:
            faceType[face.label()] = "T"
            if chords:
                nodePos = middlePoint(chords[0]) # (may be changed later)
        elif len(chords) == 2:
            faceType[face.label()] = "S"
        else:
            faceType[face.label()] = "J"
            nodePos = junctionPos(chords)

        # add nodes for non-sleeve triangles:
        if faceType[face.label()] != "S" and chords:
            nodeLabel[face.label()] = result.addNode(nodePos).label()
            result.nodeChordLabels.append([])

    for face in delaunayMap.faceIter(skipInfinite = True):
        if face.flag(OUTER_FACE):
            continue
        
        if faceType[face.label()] != "S":
            #print faceType[face.label()] + "-triangle", face.label()
            for dart in faceChords[face.label()]:
                #print "following limb starting with", dart

                startNode = result.node(nodeLabel[face.label()])
                startDart = dart.clone()
                edgePoints = []
                subtendedBoundaryLength = 0.0
                shapeWidth = []

                if faceType[face.label()] == "T":
                    subtendedBoundaryLength += boundaryLength[face.label()]
                    if includeTerminalPositions:
                        # include opposite position
                        startNode.setPosition((dart.clone().nextPhi())[-1])

                while True:
                    edgePoints.append(middlePoint(dart))
                    shapeWidth.append(dart.edge().length())
                    dart.nextAlpha()
                    nextFace = dart.leftFace()
                    if faceType[nextFace.label()] == "S":
                        subtendedBoundaryLength += boundaryLength[nextFace.label()]
                        # continue with opposite dart:
                        if dart == faceChords[nextFace.label()][0]:
                            dart = faceChords[nextFace.label()][1]
                        else:
                            dart = faceChords[nextFace.label()][0]
                    else:
                        faceChords[nextFace.label()].remove(dart)
                        break

                endNode = result.node(nodeLabel[nextFace.label()])

                if faceType[nextFace.label()] == "T":
                    subtendedBoundaryLength += boundaryLength[nextFace.label()]
                    if includeTerminalPositions:
                        endNode.setPosition((dart.clone().nextPhi())[-1])

                flags = 0
                if not edgePoints or edgePoints[0] != startNode.position():
                    edgePoints.insert(0, startNode.position())
                    shapeWidth.insert(0, None)
                    flags |= START_NODE_ADDED
                if edgePoints[-1] != endNode.position():
                    edgePoints.append(endNode.position())
                    shapeWidth.append(None)
                    flags |= END_NODE_ADDED

                if len(edgePoints) < 2:
                    # don't add edges with only 1 point (two adjacent
                    # T-triangles) - may fail for J-faces?!
                    # but J-faces should not have nodes on their borders
                    assert endNode.position() == startNode.position() and \
                           endNode.isIsolated() and startNode.isIsolated()
                    nodeLabel[nextFace.label()] = nodeLabel[face.label()]
                    result.removeIsolatedNode(endNode)
                    continue
                
                sleeve = result.addEdge(
                    startNode, endNode, edgePoints)
                sleeve.setFlag(flags)

                result.subtendedLengths.append(subtendedBoundaryLength)
                result.shapeWidths.append(shapeWidth)
                if faceType[face.label()] == "J":
                    result.nodeChordLabels[startNode.label()].append(
                        (sleeve.label(), startDart.label()))
                if faceType[nextFace.label()] == "J":
                    result.nodeChordLabels[endNode.label()].append(
                        (sleeve.label(), dart.label()))

    result.initializeMap(initLabelImage = False)

    return result

# --------------------------------------------------------------------

def _pruneBarbsInternal(skel):
    count = 0
    for edge in skel.edgeIter():
        if edge.flag(IS_BARB):
            print "pruning", edge
            count += 1
            skel.removeEdge(edge.dart())
    return count

def _leaveCircle(points, dir, center, radius):
    r2 = radius * radius
    i = dir
    try:
        while (points[i] - center).squaredMagnitude() < r2:
            i += dir
    except IndexError:
        i -= dir
    return i

def pruneBarbsByDist(skel, maxDist):
    for edge in skel.edgeIter():
        edge.setFlag(IS_BARB, False)

    maxCutLength = maxDist * math.pi # ;-)
    barbNodeLabels = [None] * skel.maxEdgeLabel()

    for node in skel.nodeIter():
        if not node.hasDegree(1):
            continue

        p = node.position()

        dart = node.anchor()
        while True:
            if (dart[0] - p).magnitude() < maxDist:
                dart.edge().setFlag(IS_BARB, True)
                if (dart[-1] - p).magnitude() > maxCutLength:
                    dart.edge().setFlag(IS_BARB, False)
                    barbNodeLabels[dart.edgeLabel()] = (dart.startNodeLabel(), p)
            elif (dart[-1] - p).magnitude() < maxDist:
                dart.edge().setFlag(IS_BARB, True)
                if (dart[0] - p).magnitude() > maxCutLength:
                    dart.edge().setFlag(IS_BARB, False)
                    barbNodeLabels[dart.edgeLabel()] = (dart.endNodeLabel(), p)
            if dart.nextPhi() == node.anchor():
                break

    result = _pruneBarbsInternal(skel)
    #removeCruft(skel, 3)

    shortenLength = maxCutLength * 2 / 3

    for edge in skel.edgeIter():
        if barbNodeLabels[edge.label()]:
            print "shortening", edge
            barbNodeLabel, endPos = barbNodeLabels[edge.label()]
            if barbNodeLabel == edge.startNodeLabel():
                #i, p = arcLength2Pos(shortenLength, edge)
                i = _leaveCircle(edge, 1, endPos, shortenLength)
                if i < len(edge)-1:
                    skel.splitEdge(edge, i).setFlag(IS_BARB, False)
                    edge.setFlag(IS_BARB, True)
            else:
                #i, p = arcLength2Pos(edge.length()-shortenLength, dart)
                i = _leaveCircle(edge, -1, endPos, shortenLength)
                if i < 0:
                    i += len(edge)
                if i < len(edge)-1:
                    newEdge = skel.splitEdge(edge, i).setFlag(IS_BARB, True)
                    edge.setFlag(IS_BARB, False)

    result += _pruneBarbsInternal(skel)

    return result

def pruneBarbs(skel):
    for edge in skel.edgeIter():
        edge.setFlag(IS_BARB, edge.startNode().hasDegree(1) or \
                      edge.endNode().hasDegree(1))
    return _pruneBarbsInternal(skel)

def pruneByMorphologicalSignificance(skel, ratio = 0.1):
    for edge in skel.edgeIter():
        edge.setFlag(IS_BARB, False)

        if edge.startNode().hasDegree(1):
            edge.setFlag(IS_BARB, True)

            endSide = edge.endSide[0]-edge.endSide[1]
            minDist = ratio * endSide.magnitude()

            endNormal = Vector2(endSide[1], -endSide[0])
            endNormal /= endNormal.magnitude()

            for p in edge:
                if abs(dot(p - edge[-1], endNormal)) > minDist:
                    edge.setFlag(IS_BARB, False)
                    break

            if edge.flag(IS_BARB):
                continue # no need to test other side

        if edge.endNode().hasDegree(1):
            edge.setFlag(IS_BARB, True)

            startSide = edge.startSide[0]-edge.startSide[1]
            minDist = ratio * max(1.0, startSide.magnitude())

            startNormal = Vector2(startSide[1], -startSide[0])
            startNormal /= startNormal.magnitude()

            for p in edge:
                if abs(dot(p - edge[0], startNormal)) > minDist:
                    edge.setFlag(IS_BARB, False)
                    break

    return _pruneBarbsInternal(skel)

def pruneBySubtendedLength(skelMap, delaunayMap = None,
                           minLength = None,
                           ratio = None):
    """pruneBySubtendedLength(skelMap, ...)

    Prunes CAT skeletons by thresholding the boundary length subtended
    by each terminal sleeve.  You can pass either of the following two
    arguments to determine the threshold:

    minLength:
      minimal boundary length in pixels

    ratio:
      fraction of total boundary length (default: 0.01 == 1%)

    Also, you can pass the delaunayMap which the skelMap has been
    extracted from as second argument; this allows this function to
    re-calculate the junction node positions from the remaining
    chords.  (Very useful if and only if using the rectified CAT.)
    Even if no delaunayMap is given, junction positions of junctions
    whose degree went down to two are removed (junction -> sleeve face
    conversion).

    Examples::

      pruneBySubtendedLength(catMap, delMap, minLength = 3)
      pruneBySubtendedLength(catMap, delMap, ratio = 0.02)
      pruneBySubtendedLength(catMap) # use default: 1%"""
    
    changed = True
    result = 0

    if minLength == None:
        totalBL = sum(skelMap.subtendedLengths)
        minLength = (ratio or 0.01) * totalBL

    from heapq import heappush, heappop

    # update nodeChordLabels and subtendedLengths via hooks:
    mergedValues = [None, None, None]
    # (merged edge label, merged subtendedLength, changed endNode)
    def preMergeEdges(dart):
        merged = dart.clone().nextSigma()
        mergedValues[0] = merged.edgeLabel()
        mergedValues[1] = skelMap.subtendedLengths[merged.edgeLabel()]
        mergedValues[2] = merged.endNodeLabel()
        return True
    def postMergeEdges(survivor):
        merged = mergedValues[0]
        chordLabel = skelMap.nodeChordLabels[mergedValues[2]]
        for i, (sleeve, chord) in enumerate(chordLabel):
            if sleeve == merged:
                chordLabel[i] = survivor.label(), chord
        skelMap.subtendedLengths[survivor.label()] += mergedValues[0]
    hooks = skelMap.addMergeEdgesCallbacks(preMergeEdges, postMergeEdges)

    barbs = []
    for edge in skelMap.edgeIter():
        subtendedBoundaryLength = skelMap.subtendedLengths[edge.label()]
        if subtendedBoundaryLength >= minLength:
            continue
        
        dart = edge.dart()
        if dart.startNode().hasMinDegree(2):
            dart.nextAlpha()
        if dart.startNode().hasMinDegree(2):
            continue # no barb (yet?)

        heappush(barbs, (subtendedBoundaryLength, dart))

    while barbs:
        subtendedBoundaryLength, dart = heappop(barbs)
        if not dart.edge() or \
           skelMap.subtendedLengths[dart.edgeLabel()] > subtendedBoundaryLength:
            continue

        # (dart now belongs to a barb, startNode().degree() == 1)
        
        neighbor = dart.clone().nextPhi()
        skelMap.subtendedLengths[neighbor.edgeLabel()] \
            += subtendedBoundaryLength

        # correct junction node position when cutting off sleeve:
        if (delaunayMap and dart.endNode().hasMinDegree(2)) or \
               (not delaunayMap and dart.endNode().hasDegree(3)):

            # update skelMap.nodeChordLabels:
            chordLabels = skelMap.nodeChordLabels[dart.endNodeLabel()]
            assert len(chordLabels) == dart.endNode().degree()
            for i, (sleeve, _) in enumerate(chordLabels):
                if sleeve == dart.edgeLabel():
                    del chordLabels[i]
                    break
            assert sleeve == dart.edgeLabel()
            
            if delaunayMap and dart.endNode().hasMinDegree(3):
                dart.endNode().setPosition(
                    rectifiedJunctionNodePosition(
                    [delaunayMap.dart(dl) for _, dl in chordLabels]))
            else: # junction face -> sleeve face
                remaining = dart.clone().nextPhi()
                re = remaining.edge()
                if remaining.label() < 0:
                    if re.flag(END_NODE_ADDED):
                        re.endNode().setPosition(re[-2])
                        del re[-1] # delete duplicate point
                else:
                    if re.flag(START_NODE_ADDED):
                        re.startNode().setPosition(re[1])
                        del re[1] # delete duplicate point

        skelMap.removeEdge(dart)
        result += 1

        if not neighbor.edge() or neighbor.startNode().hasMinDegree(3):
            continue

        if neighbor.startNode().hasDegree(2):
            survSBL = skelMap.subtendedLengths[neighbor.edgeLabel()] + \
                      skelMap.subtendedLengths[neighbor.nextSigma().edgeLabel()]
            survivor = skelMap.mergeEdges(neighbor)
            skelMap.subtendedLengths[survivor.label()] = survSBL

            if survSBL < minLength:
                dart = survivor.dart()
                if dart.startNode().hasMinDegree(2):
                    dart.nextAlpha()
                if dart.startNode().hasMinDegree(2):
                    continue # no barb (yet?)
                heappush(barbs, (survSBL, dart))
        else:
            if skelMap.subtendedLengths[neighbor.edgeLabel()] < minLength:
                heappush(barbs, (survSBL, neighbor))

    return result

# --------------------------------------------------------------------

def circumCircle(p1, p2, p3):
    """circumCircle(p1, p2, p3) -> (Vector2, float)

    Returns the center and radius of the circumcircle of the given
    three points."""
    
    p1sm = p1.squaredMagnitude()
    x1 = p1[0]
    y1 = p1[1]
    p2sm = p2.squaredMagnitude()
    x2 = p2[0]
    y2 = p2[1]
    p3sm = p3.squaredMagnitude()
    x3 = p3[0]
    y3 = p3[1]
    a = numpy.linalg.det(
        numpy.array([[x1, y1, 1],
                     [x2, y2, 1],
                     [x3, y3, 1]]))
    d = numpy.linalg.det(
        -numpy.array([[p1sm, y1, 1],
                      [p2sm, y2, 1],
                      [p3sm, y3, 1]]))
    e = numpy.linalg.det(
        numpy.array([[p1sm, x1, 1],
                     [p2sm, x2, 1],
                     [p3sm, x3, 1]]))
    f = numpy.linalg.det(
        -numpy.array([[p1sm, x1, y1],
                      [p2sm, x2, y2],
                      [p3sm, x3, y3]]))
    circumCenter = Vector2(-d/(2*a), -e/(2*a))
    try: # FIXME: HACK (DAGM deadline approaching)!!!
        circumRadius = math.sqrt((math.sq(d)+math.sq(e))/(4*math.sq(a))-f/a)
    except ValueError:
        sys.stderr.write("WARNING: Could not calculate triangle circumcircle!\n")
        lengths = [dart.edge().length()
                   for dart in contours[0].phiOrbit()] # FIXME: which `contours`?
        lengths.sort()
        circumRadius = (lengths[1] + lengths[2]) / 4.0
        sys.stderr.write("  Side lengths are %s -> improvised radius = %s\n"
                         % (lengths, circumRadius))
    return circumCenter, circumRadius

def calculateTriangleCircumcircles(delaunayMap):
    """calculateTriangleCircumcircles(delaunayMap) -> list of (Vector2, float)

    Returns a list of circumcircles for each face in the delaunayMap
    (with the face labels as indices).  Entries for faces marked with
    the OUTER_FACE flag are set to None."""
    
    result = [None] * delaunayMap.maxFaceLabel()
    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if triangle.flag(OUTER_FACE):
            continue

        assert triangle.holeCount() == 0, "delaunay triangles should not have holes (face %d)" % triangle.label()
        points = [dart[0] for dart in triangle.contour().phiOrbit()]
        assert len(points) == 3, "triangles should have three points"
        result[triangle.label()] = circumCircle(*points)
    return result

