def extractMapPoints(map):
    allPoints = [node.position() for node in map.nodeIter()]
    for edge in map.edgeIter():
        if edge.leftFaceLabel() and edge.rightFaceLabel():
            allPoints.extend(list(edge)[1:-1])
    return allPoints

def maxSegmentLength(map):
    result = 0.0
    for edge in map.edgeIter():
        if edge.leftFaceLabel() and edge.rightFaceLabel():
            result = max(result, max(
                [(edge[i+1]-edge[i]).magnitude() for i in range(len(edge)-1)]))
    return result

if not globals().has_key("dm"):
    dm = delaunayMap(allPoints, imageSize, markContour = False)
    print "  %d Delaunay edges, %d triangles" % (dm.edgeCount, dm.faceCount)

    print "- reconstructing triangle circumcircles..."
    calculateTriangleCircumcircles(dm)

# --------------------------------------------------------------------
#                         alpha shape extraction
# --------------------------------------------------------------------

def markAlphaShapes(delaunayMap, alpha, beta):
    if not hasattr(delaunayMap.faceIter(True).next(), "radius"):
        print "- reconstructing triangle circumcircles..."
        calculateTriangleCircumcircles(dm)

    print "- marking triangles with radii < alpha(%s)..." % (alpha, )
    it = delaunayMap.faceIter()
    it.next().mark = False # infinite face
    for triangle in it:
        triangle.mark = triangle.radius < alpha

    print "- marking edges with empty circle radii < alpha(%s)..." % (alpha, )
    for edge in delaunayMap.edgeIter():
        edge.mark = edge.leftFace().mark or edge.rightFace().mark
        if edge.mark:
            continue
        p1 = edge.startNode().position()
        p2 = edge.endNode().position()
        midPoint = (p1 + p2)/2
        radius = edge.length()/2
        edge.mark = radius < alpha
        if edge.mark:
            empty = True
            radius2 = math.sq(radius)
            it = edge.dart().sigmaOrbit(); it.next()
            for dart in it:
                if (dart.endNode().position()-midPoint).squaredMagnitude() < radius2:
                    empty = False
                    break
            it = edge.dart().nextAlpha().sigmaOrbit(); it.next()
            for dart in it:
                if (dart.endNode().position()-midPoint).squaredMagnitude() < radius2:
                    empty = False
                    break
            v0 = p2 - p1
            if edge.leftFaceLabel():
                v1 = edge.leftFace().circleCenter - p1
            else:
                v1 = Vector(v0[1], -v0[0])
            if edge.rightFaceLabel():
                v2 = edge.rightFace().circleCenter - p1
            else:
                v2 = Vector(-v0[1], v0[0])
            empty2 = (v2[1]*v0[0]-v2[0]*v0[1] < 0) != (v1[1]*v0[0]-v1[0]*v0[1] < 0)
            assert empty == empty2, "%s is %s/%s" % (edge, empty, empty2)
            if not empty:
    #             print "  edge %d's circumcircle contains a point, unmarking.." % (
    #                 edge.label(), )
                edge.mark = False
                continue

    # print "unmarking edges whose circles contain points..."
    # c = time.clock()
    # allPoints.sort(lambda p1, p2: cmp(p1[0], p2[0]))
    # radius2 = math.sq(radius)
    # for edge in delaunayMap.edgeIter():
    #     if not edge.mark:
    #         continue
    #     if time.clock() > c + 0.2:
    #         sys.stdout.write("\r  edge %d (length %.2f)... [%.1f%%]" % (edge.label(), edge.length(), 100.*edge.label()/delaunayMap.maxEdgeLabel())); sys.stdout.flush()
    #         c = time.clock()
    #     i = 0
    #     minX = midPoint[0]-radius
    #     while allPoints[i][0] <= minX:
    #         i += 1
    #     maxX = midPoint[0]+radius
    #     while allPoints[i][0] < maxX:
    #         diff = allPoints[i] - midPoint
    #         if abs(diff[1]) >= radius or diff.squaredMagnitude() >= radius2 or \
    #                allPoints[i] in [p1, p2]:
    #             i += 1
    #             continue
    #         print "  edge %d's circumcircle contains a point, unmarking.." % (
    #             edge.label(), )
    #         edge.mark = False
    #         break

    print "  %d/%d edges and %d/%d faces marked." % (
        sum([edge.mark and 1 or 0 for edge in delaunayMap.edgeIter()]), delaunayMap.edgeCount,
        sum([face.mark and 1 or 0 for face in delaunayMap.faceIter()]), delaunayMap.faceCount)

    print "- finding connected components of unlabelled cells..."
    for edge in delaunayMap.edgeIter():
        if not edge.mark:
            edge.componentLabel = None

    for face in delaunayMap.faceIter():
        if not face.mark:
            face.componentLabel = None

    maxSize = 0
    componentCount = 0
    for edge in delaunayMap.edgeIter():
        if edge.mark or edge.componentLabel:
            continue
        componentCount += 1
        boundary = [edge]
        size = 0
        while boundary:
            cell = boundary.pop()
            if hasattr(cell, "leftFace"):
                edge = cell
                if edge.mark or edge.componentLabel:
                    continue
                edge.componentLabel = componentCount
                size += 1
                boundary.append(edge.leftFace())
                boundary.append(edge.rightFace())
            else:
                face = cell
                if face.mark or face.componentLabel:
                    continue
                face.componentLabel = componentCount
                size += 1
                for dart in face.contours()[0].phiOrbit():
                    boundary.append(dart.edge())
        maxSize = max(maxSize, size)

    for face in delaunayMap.faceIter():
        if face.mark or face.componentLabel:
            continue
        

    print "  %s unlabelled components found." % (componentCount, )

    print "- looking for unmarked triangles with radii >= beta (%s)..." % (beta, )

    badComponent = [True] * (componentCount+1)
    for face in delaunayMap.faceIter(skipInfinite = True):
        if face.mark:
            continue
        if face.radius >= beta:
            badComponent[face.componentLabel] = False

    for label in range(1, componentCount+1):
        if badComponent[label]:
            print "  marking connected component %d." % (
                label, )
    for edge in delaunayMap.edgeIter():
        if not edge.mark:
            edge.mark = badComponent[edge.componentLabel]
    for face in delaunayMap.faceIter(skipInfinite = True):
        if not face.mark:
            face.mark = badComponent[face.componentLabel]

# --------------------------------------------------------------------
#                               fig output
# --------------------------------------------------------------------

def ccContour(startFace, condition, close = True):
    startFace.output = True
    contour = list(startFace.contours()[0].phiOrbit())
    i = 0
    while i < len(contour):
        neighbor = contour[i].rightFace()
        if condition(neighbor) and not neighbor.output:
            oldPoint = contour[i][0]
            _ = contour[i].nextAlpha().nextPhi()
            assert oldPoint == contour[i][0]
            contour.insert(i+1, contour[i].clone().nextPhi())
            neighbor.output = True
        else:
            i += 1
    result = Polygon([dart[0] for dart in contour])
    result.append(result[0]) # close poly (for filling)
    return result

def outputMarkedShapes(delaunayMap, fe,
                       regionDepth = 50, edgeDepth = 49):
    # output all cells only once; add flag first
    for edge in delaunayMap.edgeIter():
        edge.output = False

    for face in delaunayMap.faceIter():
        face.output = False

    print "- exporting marked regions as filled polygons..."
    for triangle in delaunayMap.faceIter(skipInfinite = True):
        if not triangle.mark or triangle.output:
            continue
        triangle.output = True
        contour = list(triangle.contours()[0].phiOrbit())
        i = 0
        while i < len(contour):
            contour[i].edge().output = True
            neighbor = contour[i].rightFace()
            if neighbor.mark and not neighbor.output:
                _ = contour[i].nextAlpha().nextPhi()
                contour.insert(i+1, contour[i].clone().nextPhi())
                neighbor.output = True
            else:
                i += 1
        contour = Polygon([dart[0] for dart in contour])
        contour.append(contour[0]) # close poly (for filling)
        i = 2
        while i < len(contour):
            if contour[i] == contour[i-2]:
                del contour[i-2]
                del contour[i-2]
                if i > 2:
                    i -= 1
            else:
                i += 1
        #print "%d points (area %s)" % (len(contour), contour.partialArea())
        pp = fe.addClippedPoly(contour,
                               fillStyle = fig.fillStyleSolid, depth = regionDepth)
        assert len(pp) == 1

    if edgeDepth != None:
        print "- exporting remaining marked edges..."
        for edge in delaunayMap.edgeIter():
            if not edge.mark or edge.output:
                continue
            edge.output = True
            dart = edge.dart()
            poly = Polygon(list(dart))
            drawing = True
            while drawing:
                drawing = False
                _ = dart.nextAlpha()
                next = dart.clone()
                while next.nextSigma() != dart:
                    edge = next.edge()
                    if not edge.mark or edge.output:
                        continue
                    edge.output = True
                    drawing = True
                    poly.append(next[1])
                    dart = next
                    break
            edge.output = True
            fe.addEdge(poly, depth = edgeDepth)

    for edge in dm.edgeIter():
        del edge.output
    for face in dm.faceIter():
        del face.output
