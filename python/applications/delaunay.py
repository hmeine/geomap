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

def delaunayMap(face, size, performCleaning = True, simplifyEpsilon = None):
    "USAGE: dlm = delaunayMap(polygon, mapSize)"

    c = face.contours()
    jumpPoints = [0]
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
    sys.stdout.flush()
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

            qt.qApp.processEvents()
    print

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def catMap(delaunayMap):
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
                if dot(sv, Vector(1, 1, 1)) < 1.6:
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

                edgeLabel = len(edgeTriples)
                sigmaOrbits[triangle.nodeLabel][
                    originalDarts[triangle.nodeLabel].index(dart.label())] \
                    = edgeLabel

                edgePoints = []
                if triangle.type == "T":
                    nodePositions[triangle.nodeLabel] = middlePoint(dart)
                else:
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
                    nodePositions[nextTriangle.nodeLabel] = middlePoint(dart)
                else:
                    edgePoints.append(nodePositions[nextTriangle.nodeLabel])

                sigmaOrbits[nextTriangle.nodeLabel][
                    originalDarts[nextTriangle.nodeLabel].index(dart.label())] \
                    = -edgeLabel

                edgeTriples.append((
                    triangle.nodeLabel, nextTriangle.nodeLabel, edgePoints))

            del triangle.innerDarts

    return Map(nodePositions, edgeTriples, delaunayMap.imageSize(),
               sigmaOrbits = sigmaOrbits)

def pruneSkeleton(skelMap, minLength):
    count = 0
    for edge in skelMap.edgeIter():
        if (edge.startNode().degree() == 1 or edge.endNode().degree() == 1) and \
               edge.length() < minLength:
            print "pruning", edge
            count += 1
            removeBridge(edge.dart())
    return count
