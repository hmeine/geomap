_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import fig, delaunay, sys, math
import map as spmap
from hourglass import Polygon
from math import *

__all__ = ["delaunayMap", "extractMapPoints", "midCrackPoints", "samplingPoints",
           "maxSegmentLength",
           
           "markAlphaShapes", "removeUnmarkedEdges", "alphaBetaMap",
           "findCandidatesForPointCorrection",

           "outputMarkedShapes", "view",

           "findChangeByBisection", "findMinAlpha", "findMaxBeta",

           "alphaShapeThinning1", "alphaShapeThinning"]

# --------------------------------------------------------------------

def delaunayMap(points, imageSize):
    result = delaunay.delaunayMap(points, imageSize, markContour = False)
    print "  %d Delaunay edges, %d triangles" % (result.edgeCount, result.faceCount)

    print "- reconstructing triangle circumcircles..."
    delaunay.calculateTriangleCircumcircles(result)

    return result

def extractMapPoints(map, includeNodes = True):
    if not includeNodes:
        result = []
    else:
        result = [node.position() for node in map.nodeIter()]
    for edge in map.edgeIter():
        if not edge.protection() & spmap.BORDER_PROTECTION:
            result.extend(list(edge)[1:-1])
    return result

def midCrackPoints(img):
    spmap = pixelmap.crackEdgeMap(img)
    return extractMapPoints(spmap, False)

def samplingPoints(img, threshold = 128):
    return [Vector2(pos[0], pos[1]) for pos in img.size()
            if img[pos] < threshold]

def maxSegmentLength(map):
    result = 0.0
    for edge in map.edgeIter():
        if not edge.protection() & spmap.BORDER_PROTECTION:
            result = max(result, max(
                [(edge[i+1]-edge[i]).magnitude() for i in range(len(edge)-1)]))
    return result

# --------------------------------------------------------------------
#                         alpha shape extraction
# --------------------------------------------------------------------

def markAlphaShapes(delaunayMap, alpha, beta = 0.0):
    if not hasattr(delaunayMap.faceIter(True).next(), "radius"):
        print "- reconstructing triangle circumcircles..."
        calculateTriangleCircumcircles(dm)

    # store parameters for convenience:
    delaunayMap.alpha = alpha
    delaunayMap.beta = beta

    print "- marking triangles with radii < alpha(%s)..." % (alpha, )
    it = delaunayMap.faceIter()
    it.next().mark = False # infinite face
    for triangle in it:
        triangle.mark = triangle.radius < alpha

    print "- marking edges with empty circle radii < alpha(%s)..." % (alpha, )
    for edge in delaunayMap.edgeIter():
        assert len(edge) == 2, "markAlphaShapes() expects a delaunay map!"
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
            
#             v0 = p2 - p1
#             if edge.leftFaceLabel():
#                 v1 = edge.leftFace().circleCenter - p1
#             else:
#                 v1 = Vector(v0[1], -v0[0])
#             if edge.rightFaceLabel():
#                 v2 = edge.rightFace().circleCenter - p1
#             else:
#                 v2 = Vector(-v0[1], v0[0])
#             empty2 = (v2[1]*v0[0]-v2[0]*v0[1] < 0) != (v1[1]*v0[0]-v1[0]*v0[1] < 0)

            empty2 = (edge.dart().nextSigma().endNode().position()-midPoint).squaredMagnitude() >= radius2 and (edge.dart().nextAlpha().nextSigma().endNode().position()-midPoint).squaredMagnitude() >= radius2
            
            if empty != empty2:
                sys.stderr.write("WARNING: %s is %s/%s!\n" % (edge, empty, empty2))
            if not empty2:
    #             print "  edge %d's circumcircle contains a point, unmarking.." % (
    #                 edge.label(), )
                edge.mark = False
                continue

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
                for dart in face.contour().phiOrbit():
                    boundary.append(dart.edge())

    for face in delaunayMap.faceIter():
        if face.mark or face.componentLabel:
            continue
        componentCount += 1
        face.componentLabel = componentCount

    print "  %s unlabelled components found." % (componentCount, )

    if not beta:
        return componentCount

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
            componentCount -= 1
    for edge in delaunayMap.edgeIter():
        if not edge.mark:
            edge.mark = badComponent[edge.componentLabel]
    for face in delaunayMap.faceIter(skipInfinite = True):
        if not face.mark:
            face.mark = badComponent[face.componentLabel]

    print "  %s unlabelled components left." % (componentCount, )

    return componentCount

def removeUnmarkedEdges(map, removeInterior = False):
    for edge in map.edgeIter():
        if edge.protection():
            continue
        if not edge.mark:
            removeEdge(edge.dart())
        elif removeInterior and (edge.leftFace().mark and edge.rightFace().mark):
            removeEdge(edge.dart())

def alphaBetaMap(points, imageSize, alpha, beta, removeInteriorEdges = False):
    dm = delaunayMap(points, imageSize)
    markAlphaShapes(dm, alpha, beta)
    removeUnmarkedEdges(dm, removeInteriorEdges)
    return dm

def findCandidatesForPointCorrection(abm):
    mayMove, dontMove = [], []
    for n in abm.nodeIter():
        if n.degree() != 2 or n.isAtBorder():
            dontMove.append(n.position())
        else:
            p = n.position()
            d = n.anchor()
            p0 = d.endNode().position()
            d.nextSigma()
            p1 = d.endNode().position()
            dx, dy = p1 - p0
            orientation = atan2(dy, dx)
            mayMove.append(Edgel(p[0], p[1], 1, orientation))
    return mayMove, dontMove

# --------------------------------------------------------------------
#                               fig output
# --------------------------------------------------------------------

def outputMarkedShapes(delaunayMap, fe, skipInnerEdges = True,
                       regionDepth = 50, edgeDepth = 49,
                       capStyle = fig.capStyleRound, **kwargs):
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
        contour = list(triangle.contour().phiOrbit())
        i = 0
        while i < len(contour):
            contour[i].edge().output = skipInnerEdges
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
        pp = fe.addClippedPoly(contour, depth = regionDepth,
                               fillStyle = fig.fillStyleSolid, capStyle = capStyle,
                               **kwargs)
        if len(pp) > 1:
            sys.stderr.write(
                "############## POLYGON CLIPPING FAILED! ##############\n")

    if edgeDepth != None:
        print "- exporting remaining marked edges..."
        for edge in delaunayMap.edgeIter():
            if not edge.mark or edge.output:
                continue

            dart = edge.dart()
            poly = Polygon(list(dart))
            edge.output = True

            drawing = True
            while drawing:
                drawing = False
                dart.nextAlpha()
                for next in dart.sigmaOrbit():
                    outputEdge = next.edge()
                    if not outputEdge.mark or outputEdge.output:
                        continue

                    drawing = True
                    assert poly[-1] == next[0]
                    if len(outputEdge) == 2:
                        poly.append(next[1])
                    else:
                        poly.extend(Polygon(list(next)[1:]))
                    outputEdge.output = True

                    dart = next
                    break

            # continue in the other direction:
            poly.reverse()
            dart = edge.dart().nextAlpha()
            
            drawing = True
            while drawing:
                drawing = False
                dart.nextAlpha()
                next = dart.clone()
                while next.nextSigma() != dart:
                    outputEdge = next.edge()
                    if not outputEdge.mark or outputEdge.output:
                        continue

                    drawing = True
                    assert poly[-1] == next[0]
                    poly.append(next[1])
                    outputEdge.output = True

                    dart = next
                    break

            fe.addClippedPoly(
                poly, depth = edgeDepth, capStyle = capStyle, **kwargs)

    for edge in delaunayMap.edgeIter():
        del edge.output
    for face in delaunayMap.faceIter():
        del face.output

# --------------------------------------------------------------------

def findChangeByBisection(func, goodParam, badParam, desired = None):
    if not desired:
        desired = func(goodParam)
    param = (goodParam + badParam)/2
    if abs(goodParam - badParam) < 1e-4:
        return goodParam
    current = func(param)
    print "findChangeByBisection: param = %s -> %d." % (
        param, current)
    if current == desired:
        return findChangeByBisection(func, param, badParam, desired)
    else:
        return findChangeByBisection(func, goodParam, param, desired)

def findMinAlpha(dm, goodAlpha, badAlpha, beta = 0.0):
    def countComponents(alpha, dm = dm, beta = beta):
        return markAlphaShapes(dm, alpha, beta)

    return findChangeByBisection(countComponents, goodAlpha, alpha)

def findMaxBeta(dm, alpha, badBeta):
    def countComponents(beta, dm = dm, alpha = alpha):
        return markAlphaShapes(dm, alpha, beta)

    return findChangeByBisection(countComponents, 0.0, badBeta)

# --------------------------------------------------------------------

def alphaShapeThinning1(dm):
    """Old thinning procedure, looking for particular configurations only."""
    
    changed = 0

    for edge in dm.edgeIter():
        if edge.mark == True:
            # (at least one adjacent triangle is marked)
            dart = edge.dart()
            # ensure that we have an unmarked face on the left:
            if dart.leftFace().mark:
                dart.nextAlpha()
            if dart.leftFace().mark or not dart.rightFace().mark:
                continue # no thinning here

            d1 = dart.clone()
            while not d1.nextSigma().edge().mark:
                pass
            if d1.leftFace().mark:
                continue # no thinnable config
            
            d2 = dart.clone().nextAlpha()
            while not d2.prevSigma().edge().mark:
                pass
            if d2.rightFace().mark:
                continue # no thinnable config
            
            edge.mark = False
            assert dart.rightFace().mark
            dart.rightFace().mark = False
            changed += 1

    if changed:
        changed += alphaShapeThinning1(dm)
    
    return changed

from heapq import * # requires Python 2.3+

def alphaShapeThinning(dm):
    """Region-growing based thinning."""

    def isSimple(edge):
        """returns True iff the edge is in the contour of a thick
        alpha shape region"""
        return edge.leftFace().mark != edge.rightFace().mark

    changed = 0
    border = []

    for edge in dm.edgeIter():
        if isSimple(edge):
            heappush(border, (-edge.length(), edge))

    while border:
        _, edge = heappop(border)
        if not isSimple(edge):
            continue

        dart = edge.dart()
        if not dart.leftFace().mark:
            dart.nextAlpha()

        dart.leftFace().mark = False
        edge.mark = False
        changed += 1

        dart.nextPhi()
        if isSimple(dart.edge()):
            heappush(border, (-dart.edge().length(), dart.edge()))

        dart.nextPhi()
        if isSimple(dart.edge()):
            heappush(border, (-dart.edge().length(), dart.edge()))
    
    return changed

# --------------------------------------------------------------------

import os

def view(epsFilename):
    if not os.path.exists(epsFilename) and os.path.exists(epsFilename+".eps"):
        epsFilename = epsFilename+".eps"
    os.system("gv '%s' &" % (epsFilename, ))
