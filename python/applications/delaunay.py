# -*- coding: iso-8859-1 -*-
_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import math, sys, vigra, hourglass, numpy
from hourglass import Polygon, simplifyPolygon, resamplePolygon
#from map import GeoMap, contourPoly
from hourglass import GeoMap, contourPoly
from vigra import Vector2, Vector, dot

try:
    import triangle
except ImportError:
    triangle = None

CONTOUR_SEGMENT = 256
OUTER_FACE = 1

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
    """constrainedDelaunayMap(polygons, imageSize, extraPoints = [],
                           onlyInner = True)

    Returns a new GeoMap containing a Constrained Delaunay
    Triangulation of all points of the polygons plus the extraPoints
    if given.  The segments of the polygons will be constrained
    segments of the CDT.

    If the optional onlyInner parameter is True (default), then all
    edge segments in the 'outside' will be removed.  (This assumes
    that all closed polygons with a negative partialArea() are
    holes.)"""

    assert triangle, """For correct CDT, you need to compile the
    triangle module (vigra/experiments/triangle).  You might want to
    try the home-made fakedConstrainedDelaunayMap instead, which will
    also give a correct result if possible, and just throws an
    exception if it has to give up."""
    
    points = []
    segments = []
    holes = []
    for polygon in polygons:
        l = len(points)
        partPoints = list(polygon)
        if polygon[-1] == polygon[0]:
            del partPoints[-1]
            segments.append((l+len(partPoints)-1, l))
            if onlyInner and polygon.partialArea() < 0:
                holes.append(_pointInHole(polygon))
        points.extend(partPoints)
        segments.extend([(l+i, l+i+1) for i in range(len(partPoints)-1)])

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

def delaunayMap(points, imageSize):
    """delaunayMap(points, imageSize)

    Returns a GeoMap containing a Delaunay Triangulation of the given
    points."""
    
    if triangle:
        nodePositions, edges = triangle.delaunay(points)
        sigma = None
    else:
        nodePositions, edges, sigma = hourglass.delaunay(points)
    result = _delaunayMapFromData(nodePositions, edges, imageSize,
                                  sigma)

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
        polygons = [resamplePolygon(p, resample) for p in polygons]
    if simplifyEpsilon != None:
        polygons = [simplifyPolygon(p, simplifyEpsilon) for p in polygons]

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
    return math.acos(
        dot(a, b) /
        math.sqrt(a.squaredMagnitude()*b.squaredMagnitude()))

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
    result = [None] * delaunayMap.maxEdgeLabel()
    for edge in delaunayMap.edgeIter():
        if edge.flag(CONTOUR_SEGMENT):
            continue
        result[edge.label()] = chordStrength(edge)
    return result

def chordStrengthProfile(delaunayMap, startIndex = 1, cs = None):
    if not cs:
        cs = calculateChordStrengths(delaunayMap)

    startDart = delaunayMap.node(startIndex).anchor()
    while startDart.endNodeLabel() != 2:
        startDart.nextSigma()
    assert startDart.edge().flag(CONTOUR_SEGMENT), \
           "expecting contour segment between consecutive nodes"

    dart = startDart.clone()
    result = []
    while True:
        #print dart
        while not dart.nextSigma().edge().flag(CONTOUR_SEGMENT):
            #print "-", dart
            result.append(cs[dart.edgeLabel()])
        if dart.nextAlpha() == startDart:
            #print "AT BEGIN (%s)" % dart.label()
            break
    return result

# --------------------------------------------------------------------

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def baryCenter(*points):
    return sum(points) / len(points)

def oldJunctionNodePosition(p1, p2, p3):
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

def middleOfLongestSide(p1, p2, p3):
    s = [((p2-p1).squaredMagnitude(), (p2+p1)),
         ((p3-p2).squaredMagnitude(), (p3+p2)),
         ((p1-p3).squaredMagnitude(), (p1+p3))]
    s.sort()
    return s[2][1] / 2

def junctionNodePosition(triangle):
    c = triangle.contour()
    p1 = c[0]
    p2 = c[1]
    p3 = c.nextPhi()[1]
    circumCenter, circumRadius = circumCircle(p1, p2, p3)
    if triangle.contains(circumCenter):
        return circumCenter
    return middleOfLongestSide(p1, p2, p3)

def catMap(delaunayMap,
           includeTerminalPositions = False,
           joinMiddleThreshold = 1.61):
    """Extract a CAT (chordal axis transform) from a GeoMap object
    containing a Delaunay Triangulation.
    Assumes that all edges have only two points and that all finite
    regions are triangles (such a GeoMap is returned by
    delaunayMap())."""

    result = GeoMap(delaunayMap.imageSize())

    triangleType = [None] * delaunayMap.maxFaceLabel()
    innerDarts = [None] * delaunayMap.maxFaceLabel()
    nodeLabel = [None] * delaunayMap.maxFaceLabel()

    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if triangle.flag(OUTER_FACE):
            continue

        assert triangle.holeCount() == 0
        assert len(list(triangle.contour().phiOrbit())) == 3

        dart1 = triangle.contour()
        dart2 = dart1.clone().nextPhi()
        dart3 = dart2.clone().nextPhi()
        edge1 = dart1.edge()
        edge2 = dart2.edge()
        edge3 = dart3.edge()

        # count number of edges in common with polygon contour:
        contourCount = 0
        innerDarts[triangle.label()] = []
        if edge1.flag(CONTOUR_SEGMENT):
            contourCount += 1
        else:
            innerDarts[triangle.label()].append(dart1)
        if edge2.flag(CONTOUR_SEGMENT):
            contourCount += 1
        else:
            innerDarts[triangle.label()].append(dart2)
        if edge3.flag(CONTOUR_SEGMENT):
            contourCount += 1
        else:
            innerDarts[triangle.label()].append(dart3)

        # classify into sleeve, joint, and terminal triangles:
        if contourCount == 1:
            triangleType[triangle.label()] = "S"
        elif contourCount == 0:
            triangleType[triangle.label()] = "J"
        elif contourCount == 2:
            triangleType[triangle.label()] = "T"
        else:
            triangleType[triangle.label()] = "T?"

        # add nodes for non-sleeve triangles:
        if triangleType[triangle.label()] != "S":
            if triangleType[triangle.label()] == "J":
                nodePos = junctionNodePosition(triangle)
            else:
                # we don't know better yet for terminals:
                nodePos = (dart1[0]+dart2[0]+dart3[0])/3
            nodeLabel[triangle.label()] = result.addNode(nodePos).label()

    edgeTriples = [None]

    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if triangle.flag(OUTER_FACE):
            continue
        
        if triangleType[triangle.label()] != "S":
            #print triangleType[triangle.label()] + "-triangle", triangle.label()
            for dart in innerDarts[triangle.label()]:
                #print "following limb starting with", dart

                startNode = result.node(nodeLabel[triangle.label()])

                edgePoints = []
                if triangleType[triangle.label()] == "T":
                    if not includeTerminalPositions:
                        # correct node position onto this triangle side:
                        startNode.setPosition(middlePoint(dart))
                    else:
                        # include opposite position
                        startNode.setPosition((dart.clone().nextPhi())[-1])
                        edgePoints.append(startNode.position())
                else:
                    # junction position != side -> include in edge geometry
                    edgePoints.append(startNode.position())

                while True:
                    edgePoints.append(middlePoint(dart))
                    dart.nextAlpha()
                    nextTriangle = dart.leftFace()
                    if triangleType[nextTriangle.label()] == "S":
                        if dart == innerDarts[nextTriangle.label()][0]:
                            dart = innerDarts[nextTriangle.label()][1]
                        else:
                            dart = innerDarts[nextTriangle.label()][0]
                    else:
                        innerDarts[nextTriangle.label()].remove(dart)
                        break

                endNode = result.node(nodeLabel[nextTriangle.label()])

                if triangleType[nextTriangle.label()] == "T":
                    if not includeTerminalPositions:
                        endNode.setPosition(middlePoint(dart))
                    else:
                        endNode.setPosition((dart.clone().nextPhi())[-1])
                        edgePoints.append(endNode.position())
                else:
                    edgePoints.append(endNode.position())

                result.addEdge(startNode.label(), endNode.label(), edgePoints)

    result.initializeMap(initLabelImage = False)

    return result

# --------------------------------------------------------------------

def _pruneBarbsInternal(skel):
    count = 0
    for edge in skel.edgeIter():
        if edge.isBarb:
            print "pruning", edge
            count += 1
            dart = edge.dart()
            dart.map().removeBridge(dart)
        else:
            del edge.isBarb
    return count

def pruneBarbsByLength(skel, minLength):
    for edge in skel.edgeIter():
        edge.isBarb = (edge.startNode().degree() == 1 or \
                       edge.endNode().degree() == 1) and \
                      edge.length() < minLength
    return _pruneBarbsInternal(skel)

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
        edge.isBarb = False

    maxCutLength = maxDist * math.pi # ;-)

    for node in skel.nodeIter():
        if node.degree() != 1:
            continue

        p = node.position()

        dart = node.anchor()
        while True:
            if (dart[0] - p).magnitude() < maxDist:
                dart.edge().isBarb = True
                if (dart[-1] - p).magnitude() > maxCutLength:
                    dart.edge().isBarb = False
                    dart.edge().barbNodeLabel = (dart.startNodeLabel(), p)
            elif (dart[-1] - p).magnitude() < maxDist:
                dart.edge().isBarb = True
                if (dart[0] - p).magnitude() > maxCutLength:
                    dart.edge().isBarb = False
                    dart.edge().barbNodeLabel = (dart.endNodeLabel(), p)
            if dart.nextPhi() == node.anchor():
                break

    result = _pruneBarbsInternal(skel)
    #removeCruft(skel, 3)

    shortenLength = maxCutLength * 2 / 3

    for edge in skel.edgeIter():
        edge.isBarb = False
        if hasattr(edge, "barbNodeLabel"):
            print "shortening", edge
            barbNodeLabel, endPos = edge.barbNodeLabel
            if barbNodeLabel == edge.startNodeLabel():
                #i, p = arcLength2Pos(shortenLength, edge)
                i = _leaveCircle(edge, 1, endPos, shortenLength)
                if i < len(edge)-1:
                    splitEdge(edge, i).isBarb = False
                    edge.isBarb = True
            else:
                #i, p = arcLength2Pos(edge.length()-shortenLength, dart)
                i = _leaveCircle(edge, -1, endPos, shortenLength)
                if i < len(edge)-1:
                    newEdge = splitEdge(edge, i).isBarb = True

    result += _pruneBarbsInternal(skel)

    return result

def pruneBarbs(skel):
    for edge in skel.edgeIter():
        edge.isBarb = edge.startNode().degree() == 1 or \
                      edge.endNode().degree() == 1
    return _pruneBarbsInternal(skel)

def pruneByMorphologicalSignificance(skel, ratio = 0.1):
    for edge in skel.edgeIter():
        edge.isBarb = False

        if edge.startNode().degree() == 1:
            edge.isBarb = True

            endSide = edge.endSide[0]-edge.endSide[1]
            minDist = ratio * endSide.magnitude()

            endNormal = Vector2(endSide[1], -endSide[0])
            endNormal /= endNormal.magnitude()

            for p in edge:
                if abs(dot(p - edge[-1], endNormal)) > minDist:
                    edge.isBarb = False
                    break

            if edge.isBarb:
                continue # no need to test other side

        if edge.endNode().degree() == 1:
            edge.isBarb = True

            startSide = edge.startSide[0]-edge.startSide[1]
            minDist = ratio * max(1.0, startSide.magnitude())

            startNormal = Vector2(startSide[1], -startSide[0])
            startNormal /= startNormal.magnitude()

            for p in edge:
                if abs(dot(p - edge[0], startNormal)) > minDist:
                    edge.isBarb = False
                    break

    return _pruneBarbsInternal(skel)

# --------------------------------------------------------------------

def circumCircle(p1, p2, p3):
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
                   for dart in contours[0].phiOrbit()]
        lengths.sort()
        circumRadius = (lengths[1] + lengths[2]) / 4.0
        sys.stderr.write("  Side lengths are %s -> improvised radius = %s\n"
                         % (lengths, triangle.radius))
    return circumCenter, circumRadius

def calculateTriangleCircumcircles(delaunayMap):
    result = [None] * delaunayMap.maxFaceLabel()
    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if triangle.flag(OUTER_FACE):
            continue

        assert triangle.holeCount() == 0, "delaunay triangles should not have holes (face %d)" % triangle.label()
        points = [dart[0] for dart in triangle.contour().phiOrbit()]
        assert len(points) == 3, "triangles should have three points"
        result[triangle.label()] = circumCircle(*points)
    return result

