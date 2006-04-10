import math, sys
from hourglass import Polygon, delaunay, simplifyPolygon
from map import Map, contourPoly, mergeFaces
from vigra import Vector2

def isContourEdge(edge, maxNodeLabel = None):
    "helper function for delaunayMap()"

    dist = abs(edge.startNodeLabel() - edge.endNodeLabel())
    if dist == 1:
        return True
    if maxNodeLabel == None:
        map = edge._map()
        maxNodeLabel = len(map.nodes)-1
        while not map.nodes[maxNodeLabel]:
            maxNodeLabel -= 1
    return (dist == maxNodeLabel - 1)

def delaunayMap(face, size, simplifyEpsilon = None,
                markContour = True, performCleaning = True):
    """USAGE: dlm = delaunayMap(face, mapSize)

    `face` should be a Map.Face object, and all its contours will be
    extracted.  (If a list of points or a Polygon is passed as `face`,
    it is assumed to contain exactly one contour.)  `mapSize` is used
    to initialize the Map with the resulting edges.
      
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

    performCleaning
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

        c = face.contours()
        points = contourPoly(c[0])
        if simplifyEpsilon != None:
            points = simplifyPolygon(points, simplifyEpsilon, simplifyEpsilon)
        del points[-1]
        jumpPoints.append(len(points))

        for anchor in c[1:]:
            innerContour = contourPoly(anchor)
            if simplifyEpsilon != None:
                innerContour = simplifyPolygon(
                    innerContour, simplifyEpsilon, simplifyEpsilon)
            points.extend(innerContour)
            del points[-1]
            jumpPoints.append(len(points))

    print "- performing Delaunay Triangulation (%d points)..." % len(points)
    nodePositions, edges, nodes = delaunay(points)
    #nodePositions = [None] + [node for node in nodePositions if node]

    print "- storing result in a Map..."
    edges = [startEnd and
             (startEnd[0], startEnd[1],
              [nodePositions[startEnd[0]], nodePositions[startEnd[1]]])
             for startEnd in edges]

    result = Map(nodePositions, edges, size,
                 performBorderClosing= False,
                 sigmaOrbits = nodes,
                 skipLabelImage = True)

    if not markContour:
        return result

    for edge in result.edgeIter():
        edge.isContourEdge = False

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
            assert j > 0, "Original contour fragment missing in Delauny map!"
            dart.edge().isContourEdge = dart.label()
            dart.nextAlpha()
        j = dart.startNode().degree() + 1
        while dart.endNodeLabel() != contourStartLabel and j > 0:
            dart.nextSigma()
            j -= 1
        assert j > 0, "Original contour fragment missing in Delauny map!"
        dart.edge().isContourEdge = dart.label()
        i += 1

    if performCleaning:
        cleanOuter(result)
    return result

def cleanOuter(map):
    sys.stdout.write("- reducing delaunay triangulation to inner part")
    for edge in map.edgeIter():
        if edge.isContourEdge:
            dart = map.dart(edge.isContourEdge)
            assert dart.edgeLabel() == edge.label(), str(edge)
            while True:
                mergeDart = dart.clone().prevSigma()
                if mergeDart.edge().isContourEdge:
                    break
                mergeFaces(mergeDart).isOutside = True
                sys.stdout.write(".")
                sys.stdout.flush()
    print

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def catMap(delaunayMap,
           includeTerminalPositions = False,
           joinMiddleThreshold = 1.61):
    """Extract a CAT (chordal axis transform) from a Map object
    containing a Delaunay Triangulation.
    Assumes that all edges have only two points and that all finite
    regions are triangles (such a Map is returned by
    delaunayMap())."""

    nodePositions = [None]
    sigmaOrbits = [None]
    originalDarts = [None]

    for triangle in delaunayMap.faceIter():
        if hasattr(triangle, "isOutside"):
            continue

        c = triangle.contours()
        assert len(c) == 1
        assert len(list(c[0].phiOrbit())) == 3

        dart1 = triangle.contours()[0].clone()
        dart2 = dart1.clone().nextPhi()
        dart3 = dart2.clone().nextPhi()
        edge1 = dart1.edge()
        edge2 = dart2.edge()
        edge3 = dart3.edge()

        # count number of edges in common with polygon contour:
        contourCount = 0
        triangle.innerDarts = []
        if edge1.isContourEdge:
            contourCount += 1
        else:
            triangle.innerDarts.append(dart1)
        if edge2.isContourEdge:
            contourCount += 1
        else:
            triangle.innerDarts.append(dart2)
        if edge3.isContourEdge:
            contourCount += 1
        else:
            triangle.innerDarts.append(dart3)

        # classify into sleeve, joint, and terminal triangles:
        if contourCount == 1:
            triangle.type = "S"
        elif contourCount == 0:
            triangle.type = "J"
        elif contourCount == 2:
            triangle.type = "T"
        else:
            triangle.type = "T?"

        if triangle.type != "S":
            triangle.nodeLabel = len(nodePositions)
            nodePos = (dart1[0]+dart2[0]+dart3[0])/3
            if triangle.type == "J":
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

            originalDarts.append([d.label() for d in triangle.innerDarts])
            sigmaOrbits.append([None]*len(triangle.innerDarts))

    edgeTriples = [None]

    for triangle in delaunayMap.faceIter():
        if hasattr(triangle, "isOutside"):
            continue
        if triangle.type != "S":
            #print triangle.type + "-triangle", triangle.label(), originalDarts[triangle.nodeLabel]
            for dart in triangle.innerDarts:
                #print "following limb starting with", dart

                # remember triangle side we started from for later
                # pruning by "morphological ratio":
                startSide = dart[0], dart[1]

                edgeLabel = len(edgeTriples)
                sigmaOrbits[triangle.nodeLabel][
                    originalDarts[triangle.nodeLabel].index(dart.label())] \
                    = edgeLabel

                edgePoints = []
                if triangle.type == "T":
                    if not includeTerminalPositions:
                        # correct node position onto this triangle side:
                        nodePositions[triangle.nodeLabel] = middlePoint(dart)
                    else:
                        # include opposite position
                        nodePositions[triangle.nodeLabel] = (
                            dart.clone().nextPhi())[-1]
                        edgePoints.append(nodePositions[triangle.nodeLabel])
                else:
                    # node position != side -> include in edge geometry
                    edgePoints.append(nodePositions[triangle.nodeLabel])

                while True:
                    edgePoints.append(middlePoint(dart))
                    dart.nextAlpha()
                    nextTriangle = dart.leftFace()
                    if nextTriangle.type == "S":
                        if dart == nextTriangle.innerDarts[0]:
                            dart = nextTriangle.innerDarts[1]
                        else:
                            dart = nextTriangle.innerDarts[0]
                        del nextTriangle.innerDarts
                    else:
                        nextTriangle.innerDarts.remove(dart)
                        break

                if nextTriangle.type == "T":
                    if not includeTerminalPositions:
                        # correct node position onto this triangle side:
                        nodePositions[nextTriangle.nodeLabel] = middlePoint(dart)
                    else:
                        # include opposite position
                        nodePositions[nextTriangle.nodeLabel] = (
                            dart.clone().nextPhi())[-1]
                        edgePoints.append(nodePositions[nextTriangle.nodeLabel])
                else:
                    edgePoints.append(nodePositions[nextTriangle.nodeLabel])

                sigmaOrbits[nextTriangle.nodeLabel][
                    originalDarts[nextTriangle.nodeLabel].index(dart.label())] \
                    = -edgeLabel

                # see above, for later pruning:
                endSide = dart[0], dart[1]

                edgeTriples.append((
                    triangle.nodeLabel, nextTriangle.nodeLabel, edgePoints,
                    startSide, endSide))

            del triangle.innerDarts

    result = Map(nodePositions, edgeTriples, delaunayMap.imageSize(),
                 sigmaOrbits = sigmaOrbits)
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
            removeBridge(edge.dart())
        else:
            del edge.isBarb
    return count

def pruneBarbsByLength(skel, minLength):
    for edge in skel.edgeIter():
        edge.isBarb = (edge.startNode().degree() == 1 or \
                       edge.endNode().degree() == 1) and \
                      edge.length() < minLength
    return _pruneBarbsInternal(skel)

def leaveCircle(points, dir, center, radius):
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
                i = leaveCircle(edge, 1, endPos, shortenLength)
                if i < len(edge)-1:
                    splitEdge(edge, i).isBarb = False
                    edge.isBarb = True
            else:
                #i, p = arcLength2Pos(edge.length()-shortenLength, dart)
                i = leaveCircle(edge, -1, endPos, shortenLength)
                if i < len(edge)-1:
                    newEdge = splitEdge(edge, i).isBarb = True

    result += _pruneBarbsInternal(skel)

    #return result

#             pruneDart = dart.clone()
#             dart.nextPhi()
#             removeBridge(pruneDart)

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

    it = delaunayMap.faceIter(); it.next() # skip infinite face
    for triangle in it:
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
        triangle.radius = math.sqrt((math.sq(d)+math.sq(e))/(4*math.sq(a))-f/a)
