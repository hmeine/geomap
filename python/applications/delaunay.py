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

def delaunayMap(edge, size):
    "USAGE: dlm = delaunayMap(polygon, mapSize)"
    
    print "- performing Delaunay Triangulation..."
    nodePositions, edges, nodes = delaunay(edge)

    print "- storing result in a Map..."
    edges = [startEnd and
             (startEnd[0], startEnd[1],
              [nodePositions[startEnd[0]], nodePositions[startEnd[1]]])
             for startEnd in edges]
    
    result = Map(nodePositions, edges, size,
                 performBorderClosing= False,
                 sigmaOrbits = nodes,
                 skipLabelImage = True)

    maxNodeLabel = len(edge)-1
    while not result.nodes[maxNodeLabel]:
        # if `edge' is a closed polygon, adding the last point will
        # not have lead to a new node:
        maxNodeLabel -= 1

    # "peel off" outer delaunay edges:
    print "- reducing delaunay triangulation to inner part",
    isClean = False
    while not isClean:
        isClean = True
        # start with outer edge
        dart = result.face(0).contours()[0].clone()
        while True:
            if not isContourEdge(dart):
                while not isContourEdge(dart, maxNodeLabel):
                    sys.stdout.write("."); sys.stdout.flush()
                    startNode = dart.startNode()
                    mergeFaces(dart)
                    dart = startNode.anchor()
                    while dart.leftFaceLabel() != 0:
                        dart.nextSigma()
                isClean = False

            # walk around outer contour:
            if dart.nextPhi() == result.face(0).contours()[0]:
                break

    return result

def middlePoint(twoPointEdge):
    return (twoPointEdge[0] + twoPointEdge[1])/2

def catMap(delaunayMap):
    """Extract a CAT (chordal axis transform) from a Map object
    containing a Delaunay Triangulation.
    Assumes that all edges have only two points and that all finite
    regions are triangles (such a Map is returned by
    delaunayMap())."""

    maxNodeLabel = len(delaunayMap.nodes)-1
    while not delaunayMap.nodes[maxNodeLabel]:
        maxNodeLabel -= 1

    for edge in delaunayMap.edgeIter():
        edge.isContourEdge = isContourEdge(edge, maxNodeLabel)

    nodePositions = [None]
    sigmaOrbits = [None]
    originalDarts = [None]

    faceIt = delaunayMap.faceIter()
    faceIt.next() # skip infinite face
    for triangle in faceIt:
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
            nodePositions.append((dart1[0]+dart2[0]+dart3[0])/3)
            originalDarts.append([d.label() for d in triangle.innerDarts])
            sigmaOrbits.append([None]*len(triangle.innerDarts))

    edgeTriples = [None]

    faceIt = delaunayMap.faceIter()
    faceIt.next() # skip infinite face
    for triangle in faceIt:
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
