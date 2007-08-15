import sys, copy
import vigra, hourglass, maputils
from vigra import Vector2, Point2D

__all__ = ["levelSetMap", "marchingSquares"]

# TODO:
# 1) *border handling* in followContour() (predictorStep, correctorStep)!
# 3) make followContour() check for both possible intersections
#    and remove the duplicate filtering (cf. levelcontours_both_intersections.diff)
# 4) don't store zc positions in a GeoMap, but use a PositionedMap and
#    cancel out points instead of starting new edges everytime
#    (-> consecutive node and edge labels)
# 5) try larger stepsizes and a good criterion for decreasing h
#    (e.g. the size of the corrector step?!)

def findZeroCrossingsOnGrid(siv, level, minDist = 0.1):
    result = []
    existing = hourglass.PositionedMap()
    minDist2 = minDist*minDist

    def addIntersection(p):
        if not existing(p, minDist2):
            result.append(p)
            existing.insert(p, True)
    
    for y in range(siv.height()-1):
        for x in range(siv.width()-1):
            coeff = siv.coefficients(x, y)

            xPoly = [coeff[k,0] for k in range(coeff.width())]
            xPoly[0] -= level
            try:
                for k in vigra.polynomialRealRoots(xPoly):
                    if k < 0.0 or k >= 1.0:
                        continue
                    addIntersection(Vector2(x+k, y))
            except Exception, e:
                sys.stderr.write("WARNING: no convergence in polynomialRealRoots(%s):\n  %s\n" % (xPoly, e))

            yPoly = [coeff[0,k] for k in range(coeff.height())]
            yPoly[0] -= level
            try:
                for k in vigra.polynomialRealRoots(yPoly):
                    if k < 0.0 or k >= 1.0:
                        continue
                    addIntersection(Vector2(x, y+k))
            except Exception, e:
                sys.stderr.write("WARNING: no convergence in polynomialRealRoots(%s):\n  %s\n" % (yPoly, e))
    return result

# --------------------------------------------------------------------

def tangentDir(siv, pos):
    result = Vector2(-siv.dy(pos[0], pos[1]), siv.dx(pos[0], pos[1]))
    return result / result.magnitude()

def predictorStep(siv, pos, h):
    """predictorStep(siv, pos, h) -> Vector2

    Step distance h from pos in direction of tangent (perpendicular
    to gradient).  Returns None if that point is outside the
    SplineImageView."""
    
    result = pos + h*tangentDir(siv, pos)
    if not siv.isInside(result[0], result[1]):
        return None
    return result

def correctorStep(siv, level, pos, epsilon = 1e-8):
    """Perform corrector step, i.e. perform 1D iterative Newton method in
    direction of gradient in order to return to zero level (with
    accuracy given by epsilon)."""
    
    x, y = pos
    n = Vector2(siv.dx(x, y), siv.dy(x, y))
    n /= n.magnitude()

    for k in range(100):
        value = siv(x, y) - level
        if abs(value) < epsilon:
            break
        
        g = vigra.dot(Vector2(siv.dx(x, y), siv.dy(x, y)), n)
        if not g:
            sys.stderr.write("WARNING: correctorStep: zero gradient!\n")
            break # FIXME: return None instead?
        correction = -value * n / g

        # prevent too large steps (i.e. if norm(g) is small):
        if correction.squaredMagnitude() > 0.25:
            correction /= 20*correction.magnitude()

        x += correction[0]
        y += correction[1]
        
        if not siv.isInside(x, y):
            return # out of range

    return Vector2(x, y)

def predictorCorrectorStep(siv, level, pos, h, epsilon):
    while abs(h) > 1e-6:
        p1 = predictorStep(siv, pos, h)
        p2 = correctorStep(siv, level, p1, epsilon)
        if not p2 or (p2 - p1).squaredMagnitude() > h:
            h /= 2.0
            continue
        h *= 2
        return p2, h
    sys.stderr.write(
        "WARNING: predictorCorrectorStep: not converged at %s!\n" % p2)
    return p2, h

# --------------------------------------------------------------------

def followContour(siv, level, geomap, nodeLabel, h):
    #global pos, ip, poly, startNode, diff, npos, nip, intersection
    startNode = geomap.node(nodeLabel)
    pos = startNode.position()
    ix = int(pos[0])
    iy = int(pos[1])
    poly = [pos]
    while True:
        npos, nh = predictorCorrectorStep(siv, level, pos, h, 1e-6)
        h = max(min(h, nh), 1e-5)
        nix = int(npos[0])
        niy = int(npos[1])
        if nix != ix or niy != iy:
            # determine grid intersection:
            diff = npos - pos
            if nix != ix:
                intersectionX = round(npos[0])
                intersectionY = pos[1]+(intersectionX-pos[0])*diff[1]/diff[0]
            else:
                intersectionY = round(npos[1])
                intersectionX = pos[0]+(intersectionY-pos[1])*diff[0]/diff[1]
            intersection = Vector2(intersectionX, intersectionY)

            # connect to crossed Node:
            endNode = geomap.nearestNode(intersection, 0.01)
            if endNode and endNode.label() == startNode.label():
                if len(poly) > 2:
                    sys.stderr.write("WARNING: ignoring self-crossing from node %d  (poly len: %d)\n" % (startNode.label(), len(poly)))
            elif endNode:
                poly.append(endNode.position())
                geomap.addEdge(startNode, endNode, poly)
                if not endNode.degree() % 2:
                    return
                poly = [endNode.position()]
                startNode = endNode
            else:
                sys.stderr.write("WARNING: level contour crossing grid at %s without intersection Node!\n" % repr(intersection))

            ix = nix
            iy = niy
        
        poly.append(npos)
        pos = npos

def levelSetMap(image, level, sigma = None):
    siv = hasattr(image, "siv") and image.siv or vigra.SplineImageView3(image)
    
    zc = findZeroCrossingsOnGrid(siv, level)
    result = hourglass.GeoMap(zc, [], image.size())

    for node in result.nodeIter():
        if node.isIsolated():
            followContour(siv, level, result, node.label(), 0.1)

    maputils.mergeDegree2Nodes(result)
    #maputils.connectBorderNodes(result, 0.01)

    result.initializeMap()
    return result
    
# --------------------------------------------------------------------

def marchingSquares(image, level, variant = True):
    """Return a new GeoMap with sub-pixel level contours extracted by
    the marching squares method.  (Pixels with values < level are
    separated from pixels >= level.)

    If the image does not have an attribute 'siv', standard linear
    interpolation is used.  If image.siv exists, it should be a
    SplineImageView that is used for the Newton-Raphson method to
    perform another subsequent sub-pixel correction.

    The optional parameter `variant` determines the handling of the
    ambiguous diagonal configuration:

    `variant` = True (default)
      always let the two sampling points above `level` be connected

    `variant` = False
      always let the two opposite sampling points < `level` be connected

    `variant` = SplineImageView(...)
      for each ambiguous configuration, check the midpoint of the
      square; then handle as if variant = (midpoint >= level)"""
    
    connections1 = ((1, 0), (0, 2), (1, 2), (3, 1), (3, 0), (0, 2), (3, 1), (3, 2), (2, 3), (1, 0), (2, 3), (0, 3), (1, 3), (2, 1), (2, 0), (0, 1))
    connections2 = ((1, 0), (0, 2), (1, 2), (3, 1), (3, 0), (0, 1), (3, 2), (3, 2), (2, 3), (1, 3), (2, 0), (0, 3), (1, 3), (2, 1), (2, 0), (0, 1))
    configurations = (0, 0, 1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 16)

    result = hourglass.GeoMap(image.size())
    
    def addNodeX(x, y, ofs):
        return result.addNode((x + ofs, y))

    def addNodeY(x, y, ofs):
        return result.addNode((x, y + ofs))

    def addNodeNewtonRefinementX(x, y, ofs):
        for i in range(100):
            o = -(image.siv(x+ofs, y)-level) / image.siv.dx(x+ofs, y)
            if abs(o) > 0.5:
                o = vigra.sign(o)*0.05
            ofs += o
            if ofs <= 0 or ofs >= 1:
                ofs -= o
                break
            if abs(o) < 1e-4:
                break
        return result.addNode((x+ofs, y))

    def addNodeNewtonRefinementY(x, y, ofs):
        for i in range(100):
            o = -(image.siv(x, y+ofs)-level) / image.siv.dy(x, y+ofs)
            if abs(o) > 0.5:
                o = vigra.sign(o)*0.05
            ofs += o
            if ofs <= 0 or ofs >= 1:
                ofs -= o
                break
            if abs(o) < 1e-4:
                break
        return result.addNode((x, y+ofs))

    if hasattr(image, "siv"):
        addNodeX = addNodeNewtonRefinementX
        addNodeY = addNodeNewtonRefinementY

    hNodes = vigra.GrayImage(image.size())
    for y in range(image.height()):
        for x in range(image.width()-1):
            v1 = image[x,   y]
            v2 = image[x+1, y]
            if (v1 < level) != (v2 < level):
                ofs = (level - v1)/(v2 - v1)
                hNodes[x, y] = addNodeX(x, y, ofs).label()

    vNodes = vigra.GrayImage(image.size())
    for y in range(image.height()-1):
        for x in range(image.width()):
            v1 = image[x, y]
            v2 = image[x, y+1]
            if (v1 < level) != (v2 < level):
                ofs = (level - v1)/(v2 - v1)
                vNodes[x, y] = addNodeY(x, y, ofs).label()

    nodes = (hNodes, vNodes, vNodes, hNodes)
    offsets = (Point2D(0, 0), Point2D(0, 0), Point2D(1, 0), Point2D(0, 1))

    defaultConnections = connections1
    if variant == False:
        defaultConnections = connections2
    if isinstance(variant, bool):
        variant = None

    for y in range(image.height()-1):
        for x in range(image.width()-1):
            config = int(image[x,   y  ] < level)   + \
                     int(image[x+1, y  ] < level)*2 + \
                     int(image[x,   y+1] < level)*4 + \
                     int(image[x+1, y+1] < level)*8
            connections = defaultConnections
            if variant is not None and config in (6, 9):
                if variant(x + 0.5, y + 0.5) < level:
                    connections = connections2
            for s, e in connections[
                configurations[config]:configurations[config+1]]:
                s = result.node(int(nodes[s][offsets[s] + (x, y)]))
                e = result.node(int(nodes[e][offsets[e] + (x, y)]))
                result.addEdge(s, e, [s.position(), e.position()])

    maputils.mergeDegree2Nodes(result) # node suppression
    result = maputils.copyMapContents(result)[0] # compress edge labels
    maputils.connectBorderNodes(result, 0.5)
    result.initializeMap()
    return result

# --------------------------------------------------------------------

# from vigra import addPathFromHere
# addPathFromHere("../evaluation/")
# import edgedetectors

# def levelSetMap(image, level, sigma = None):
#     ed = edgedetectors.EdgeDetector(
#         bi = "Thresholding", s1 = sigma, nonmax = "zerosSubPixel",
#         threshold = level)
#     result, _, _ = ed.computeMap(image)
#     maputils.mergeDegree2Nodes(result)
#     return result
