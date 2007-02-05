# -*- coding: iso-8859-1 -*-
_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import math, sys, vigra, hourglass
from hourglass import Polygon, simplifyPolygon
#from map import GeoMap, contourPoly
from hourglass import GeoMap, contourPoly
from vigra import Vector2, Vector, dot

try:
    import triangle
except ImportError:
    triangle = None

# --------------------------------------------------------------------
#                                 TODO:
# --------------------------------------------------------------------
# * how to pass "outer faces" to the triangle module?
# * how to mark "outer faces" of the CDT maps?
#   (needed by the CAT)
# --------------------------------------------------------------------

CONTOUR_PROTECTION = 2

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

def constrainedDelaunayMap(points, jumpPoints, imageSize,
                           contourProtection = CONTOUR_PROTECTION,
                           onlyInner = True):
    """constrainedDelaunayMap(points, jumpPoints, imageSize,
                           contourProtection = CONTOUR_PROTECTION,
                           onlyInner = True)

    Returns a new GeoMap containing a Constrained Delaunay
    Triangulation of the given points, where jumpPoints determines
    which successive points are not to be connected by constraining
    segments. FIXME: Fix this API. ;-)"""

    assert triangle, """For correct CDT, you need to compile the
    triangle module (vigra/experiments/triangle).  You might want to
    try the home-made fakedConstrainedDelaunayMap instead, which will
    also give a correct result if possible, and just throws an
    exception if it has to give up."""
    
    print "- performing Constrained Delaunay Triangulation (%d points)..." % len(points)
    segments = [(i-1, i) for i in range(len(points)+1)]
    for i, jp in enumerate(jumpPoints[:-1]):
        segments[jp] = (jumpPoints[i+1]-1, jumpPoints[i])
    del segments[-1]

    nodePositions, edgeData = triangle.constrainedDelaunay(
        points, segments, onlyInner)

    print "- storing result in a GeoMap..."
    result = _delaunayMapFromData(nodePositions, edgeData, imageSize)

    for edge in result.edgeIter():
        if edgeData[edge.label()][2]:
            edge.protect(contourProtection)

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

def fakedConstrainedDelaunayMap(points, jumpPoints, imageSize,
                                contourProtection = CONTOUR_PROTECTION,
                                onlyInner = True):

    """See constrainedDelaunayMap, this calculates a DT and throws
    away outer edges retroactively.  This may fail when the DT does
    not contain all constrained segments, which is checked in this
    function and leads to an AssertionError."""
    
    print "- performing Delaunay Triangulation (%d points)..." % len(points)
    nodePositions, edges, sigma = hourglass.delaunay(points)
    
    print "- storing result in a GeoMap..."
    result = _delaunayMapFromData(nodePositions, edges, imageSize, sigma)

    if not contourProtection:
        return result

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
            dart.edge().protect(contourProtection)
            edgeSourceDarts[dart.edgeLabel()] = dart.label()
            dart.nextAlpha()
        j = dart.startNode().degree() + 1
        while dart.endNodeLabel() != contourStartLabel and j > 0:
            dart.nextSigma()
            j -= 1
        assert j > 0, "Original contour fragment missing in Delauny map!"
        dart.edge().protect(contourProtection)
        edgeSourceDarts[dart.edgeLabel()] = dart.label()
        i += 1

    if onlyInner:
        sys.stdout.write("- reducing delaunay triangulation to inner part...")
        outerFace = [False] * result.maxFaceLabel()
        for edge in result.edgeIter():
            if edge.protection() & contourProtection:
                dart = result.dart(edgeSourceDarts[edge.label()])
                assert dart.edgeLabel() == edge.label(), str(edge)
                while True:
                    mergeDart = dart.clone().prevSigma()
                    if mergeDart.edge().protection() & contourProtection:
                        break
                    outerFace[result.mergeFaces(mergeDart).label()] = True
                    sys.stdout.write(".")
                    sys.stdout.flush()
        sys.stdout.write("\n")
    
    return result

def faceCDTMap(face, imageSize, simplifyEpsilon = None,
                contourProtection = CONTOUR_PROTECTION,
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

    if type(face) == tuple:
        points, jumpPoints = face
    elif type(face) in (list, Polygon):
        points = Polygon(face)
        if points[-1] == points[0]:
            del points[-1]
        jumpPoints = [0, len(points)]
    else:
        jumpPoints = [0]

        points = contourPoly(face.contour())
        if simplifyEpsilon != None:
            points = simplifyPolygon(points, simplifyEpsilon, simplifyEpsilon)
        del points[-1]
        jumpPoints.append(len(points))

        for anchor in face.holeContours():
            innerContour = contourPoly(anchor)
            if simplifyEpsilon != None:
                innerContour = simplifyPolygon(
                    innerContour, simplifyEpsilon, simplifyEpsilon)
            points.extend(innerContour)
            del points[-1]
            jumpPoints.append(len(points))

    if triangle:
        return constrainedDelaunayMap(
            points, jumpPoints, imageSize, contourProtection, onlyInner)
    else:
        return fakedConstrainedDelaunayMap(
            points, jumpPoints, imageSize, contourProtection, onlyInner)

# --------------------------------------------------------------------

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def catMap(delaunayMap,
           includeTerminalPositions = False,
           joinMiddleThreshold = 1.61):
    """Extract a CAT (chordal axis transform) from a GeoMap object
    containing a Delaunay Triangulation.
    Assumes that all edges have only two points and that all finite
    regions are triangles (such a GeoMap is returned by
    delaunayMap())."""

    nodePositions = [None]
    sigmaOrbits = [None]
    originalDarts = [None]

    triangleType = [None] * delaunayMap.maxFaceLabel()
    innerDarts = [None] * delaunayMap.maxFaceLabel()
    nodeLabel = [None] * delaunayMap.maxFaceLabel()

    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if hasattr(triangle, "isOutside"):
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
        if edge1.protection() & CONTOUR_PROTECTION:
            contourCount += 1
        else:
            innerDarts[triangle.label()].append(dart1)
        if edge2.protection() & CONTOUR_PROTECTION:
            contourCount += 1
        else:
            innerDarts[triangle.label()].append(dart2)
        if edge3.protection() & CONTOUR_PROTECTION:
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

        if triangleType[triangle.label()] != "S":
            nodeLabel[triangle.label()] = len(nodePositions)
            nodePos = (dart1[0]+dart2[0]+dart3[0])/3
            if triangleType[triangle.label()] == "J":
                a = (dart2[0]-dart1[0]).magnitude()
                b = (dart3[0]-dart2[0]).magnitude()
                c = (dart1[0]-dart3[0]).magnitude()
                sv = Vector(a, b, c)
                sv /= sv.magnitude()
                if dot(sv, Vector(1, 1, 1)) < joinMiddleThreshold:
                    s = [a,b,c]
                    s.sort()
                    if s[2] == a:
                        nodePos = (dart2[0]+dart1[0])/2
                    elif s[2] == b:
                        nodePos = (dart3[0]+dart2[0])/2
                    else:
                        nodePos = (dart1[0]+dart3[0])/2
            nodePositions.append(nodePos)

            originalDarts.append([d.label() for d in innerDarts[triangle.label()]])
            sigmaOrbits.append([None]*len(innerDarts[triangle.label()]))

    edgeTriples = [None]

    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if hasattr(triangle, "isOutside"):
            continue
        
        if triangleType[triangle.label()] != "S":
            #print triangleType[triangle.label()] + "-triangle", triangle.label(), originalDarts[nodeLabel[triangle.label()]]
            for dart in innerDarts[triangle.label()]:
                #print "following limb starting with", dart

                # remember triangle side we started from for later
                # pruning by "morphological ratio":
                startSide = dart[0], dart[1]

                edgeLabel = len(edgeTriples)
                sigmaOrbits[nodeLabel[triangle.label()]][
                    originalDarts[nodeLabel[triangle.label()]].index(dart.label())] \
                    = edgeLabel

                edgePoints = []
                if triangleType[triangle.label()] == "T":
                    if not includeTerminalPositions:
                        # correct node position onto this triangle side:
                        nodePositions[nodeLabel[triangle.label()]] = middlePoint(dart)
                    else:
                        # include opposite position
                        nodePositions[nodeLabel[triangle.label()]] = (
                            dart.clone().nextPhi())[-1]
                        edgePoints.append(nodePositions[nodeLabel[triangle.label()]])
                else:
                    # node position != side -> include in edge geometry
                    edgePoints.append(nodePositions[nodeLabel[triangle.label()]])

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

                if triangleType[nextTriangle.label()] == "T":
                    if not includeTerminalPositions:
                        # correct node position onto this triangle side:
                        nodePositions[nodeLabel[nextTriangle.label()]] = middlePoint(dart)
                    else:
                        # include opposite position
                        nodePositions[nodeLabel[nextTriangle.label()]] = (
                            dart.clone().nextPhi())[-1]
                        edgePoints.append(nodePositions[nodeLabel[nextTriangle.label()]])
                else:
                    edgePoints.append(nodePositions[nodeLabel[nextTriangle.label()]])

                sigmaOrbits[nodeLabel[nextTriangle.label()]][
                    originalDarts[nodeLabel[nextTriangle.label()]].index(dart.label())] \
                    = -edgeLabel

                # see above, for later pruning:
                endSide = dart[0], dart[1]

                edgeTriples.append((
                    nodeLabel[triangle.label()], nodeLabel[nextTriangle.label()], edgePoints,
                    startSide, endSide))

    result = GeoMap(nodePositions, edgeTriples, delaunayMap.imageSize())
    result.initializeMap(initLabelImage = False)
    # FIXME: re-use above sigmaOrbits
    for edge in result.edgeIter():
        edge.startSide = edgeTriples[edge.label()][3]
        edge.endSide = edgeTriples[edge.label()][4]
    return result

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

def calculateTriangleCircumcircles(delaunayMap):
    import Numeric, LinearAlgebra

    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if hasattr(triangle, "isOutside"):
            continue

        contours = triangle.contours()
        assert len(contours) == 1, "delaunay triangles should not have holes (face %d)" % triangle.label()
        points = [dart.startNode().position()
                  for dart in contours[0].phiOrbit()]
        assert len(points) == 3, "triangles should have three points"
        # calculate radius r (formulas from MathWorld):
        p1sm = points[0].squaredMagnitude()
        x1 = points[0][0]
        y1 = points[0][1]
        p2sm = points[1].squaredMagnitude()
        x2 = points[1][0]
        y2 = points[1][1]
        p3sm = points[2].squaredMagnitude()
        x3 = points[2][0]
        y3 = points[2][1]
        a = LinearAlgebra.determinant(
            Numeric.array([[x1, y1, 1],
                           [x2, y2, 1],
                           [x3, y3, 1]]))
        d = LinearAlgebra.determinant(
            -Numeric.array([[p1sm, y1, 1],
                            [p2sm, y2, 1],
                            [p3sm, y3, 1]]))
        e = LinearAlgebra.determinant(
            Numeric.array([[p1sm, x1, 1],
                           [p2sm, x2, 1],
                           [p3sm, x3, 1]]))
        f = LinearAlgebra.determinant(
            -Numeric.array([[p1sm, x1, y1],
                            [p2sm, x2, y2],
                            [p3sm, x3, y3]]))
        triangle.circleCenter = Vector2(-d/(2*a), -e/(2*a))
        try: # FIXME: HACK (DAGM deadline approaching)!!!
            triangle.radius = math.sqrt((math.sq(d)+math.sq(e))/(4*math.sq(a))-f/a)
        except ValueError:
            sys.stderr.write("WARNING: Could not calculate triangle circumcircle!\n")
            lengths = [dart.edge().length()
                       for dart in contours[0].phiOrbit()]
            lengths.sort()
            triangle.radius = (lengths[1] + lengths[2]) / 4.0
            sys.stderr.write("  Side lengths are %s -> improvised radius = %s\n"
                             % (lengths, triangle.radius))
