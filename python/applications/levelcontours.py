import sys, copy
from vigra import Vector2, polynomialRealRoots

def findZeroCrossingsOnGrid(siv):
    result = []
    # FIXME: why is the order of those loops important?!?
    for y in range(siv.height()-1):
        for x in range(siv.width()-1):
            coeff = siv.coefficients(x, y)

            xPoly = [coeff[k,0] for k in range(4)]
            try:
                for k in polynomialRealRoots(xPoly):
                    if k < 0.0 or k >= 1.0:
                        continue
                    result.append(Vector2(x+k, y))
            except Exception, e:
                sys.stderr.write("WARNING: no convergence in polynomialRealRoots(%s):\n  %s\n" % (xPoly, e))

            yPoly = [coeff[0,k] for k in range(4)]
            try:
                for k in polynomialRealRoots(yPoly):
                    if k < 0.0 or k >= 1.0:
                        continue
                    result.append(Vector2(x, y+k))
            except Exception, e:
                sys.stderr.write("WARNING: no convergence in polynomialRealRoots(%s):\n  %s\n" % (yPoly, e))
    
    return result

# --------------------------------------------------------------------

def gradient(siv, pos):
    return Vector2(siv.dx(pos[0], pos[1]), siv.dy(pos[0], pos[1]))

def gradientDir(siv, pos):
    result = gradient(siv, pos)
    return result / result.magnitude()

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

def correctorStep(siv, pos, epsilon = 1e-8, n0 = 3):
    result = copy.copy(pos) # copy

    for k in range(100):
        value = siv(result[0], result[1])
        if abs(value) < epsilon:
            break
        g = gradient(siv, result)
        correction = value / g.squaredMagnitude()
        result -= g * correction
        if not siv.isInside(result[0], result[1]):
            return None, None # out of range
        if (result-pos).squaredMagnitude() > 1.0: # FIXME
            return None, None

    # FIXME: find out why this does not always converge to epsilons
    # around 1e-8, but don't throw away good points (value < 1e-5)!
#     if abs(value) > max(epsilon, 1e-5):
#         return None, None

    contraction = float(k) / n0
    return result, contraction

def predictorCorrectorStep(siv, pos, h, epsilon):
    while abs(h) > 1e-6: # FIXME
        p1 = predictorStep(siv, pos, h)
        if not p1:
            h /= 2.0
            continue
        p2, contraction = correctorStep(siv, p1, epsilon)
        if not p2:
            h /= 2.0
            continue
        h /= min(2.0, max(0.5, pow(contraction, 0.5)))
        if abs(h) > 0.2:
            h = h < 0 and -0.2 or 0.2
        break
    else:
        return None, None
    return p2, h

# --------------------------------------------------------------------

def followContour(siv, geomap, nodeLabel, h):
    pos = geomap.node(nodeLabel).position()
    x = round(pos[0])
    y = round(pos[1])
    poly = [pos]
    while True:
        npos, _ = predictorCorrectorStep(siv, pos, h, 1e-6)
        if not npos:
            return poly
        nx = int(npos[0])
        ny = int(npos[1])
        if nx != x or ny != y:
            # determine grid intersection
            diff = npos - pos
            if nx != x:
                intersectionX = round(npos[0])
                intersectionY = pos[1]+(intersectionX-pos[0])*diff[1]/diff[0]
            else:
                intersectionY = round(npos[1])
                intersectionX = pos[0]+(intersectionY-pos[1])*diff[0]/diff[1]
            intersection = Vector2(intersectionX, intersectionY)

            # connect to crossed Node
            node = geomap.nearestNode(intersection, 0.01)
            if node and node.label() == nodeLabel: # and len(poly) < 2:
                print "coming from node %d to %d, ignoring crossing, poly len: %d" \
                      % (nodeLabel, node.label(), len(poly))                    
                pass
            elif node:
                poly.append(node.position())
                print "added", geomap.addEdge(nodeLabel, node.label(), poly)
                if not node.degree() % 2:
                    return
                poly = [node.position()]
                nodeLabel = node.label()
            else:
                sys.stderr.write("WARNING: level contour crossing grid at %s without intersection Node!\n" % repr(intersection))
            x = nx
            y = ny
        poly.append(npos)
        pos = npos

def haralickConstraint(z, i, x, y, t):
    if i.g2(x, y) < t:
        return False
    gx, gy = i.dx(x, y), i.dy(x, y)
    return (gx**3 * i.dx3(x, y) + \
            3*gx**2 * gy * i.dxxy(x, y) + \
            3*gx * gy**2 * i.dxyy(x, y) + \
            gy**3 * i.dy3(x, y) < 0)

def laplaceConstraint(z, i, x, y, t):
    if i.g2(x, y) < t:
        return False
    gx, gy = -z.dx(x, y), -z.dy(x, y)
    return (gx**3 * i.dx3(x, y) + \
            3*gx**2 * gy * i.dxxy(x, y) + \
            3*gx * gy**2 * i.dxyy(x, y) + \
            gy**3 * i.dy3(x, y) < 0)

def splineConstraint(z, i, x, y, t):
    if i(x,y) < t:
        return False
#      gx, gy, gxx, gxy, gyy = i.dx(x,y), i.dy(x,y), i.dxx(x,y), i.dxy(x,y), i.dyy(x,y)
#      l = gy**2 * gxx - 2*gx*gy*gxy + gx**2 * gyy
#      r = gx**2 * gxx + 2*gx*gy*gxy + gy**2 * gyy
#     return (l > r and l > 0.0)
#      return (l <= 0.0)
    gx, gy, gxx, gxy, gyy = z.dx(x,y), z.dy(x,y), i.dxx(x,y), i.dxy(x,y), i.dyy(x,y)
    return gx**2 * gxx + 2*gx*gy*gxy + gy**2 * gyy <= 0.0

class ZeroEdges:
    """zero crossings of an image, its oriented second derivative, the
    Laplacian, or the height ridge by splines"""
    
    def __init__(self, image, method = "haralick"):
        """method should be one of:
        'direct', 'haralick', 'laplace', or 'splineridge'"""
        
        s = SplineImageView5(image)
        self.i = s
        z = GrayImage(image.size())
        if method is "direct":
            z = image
        elif method is "haralick":
            for x,y in image.size():
                z[x,y] =  s.dx(x,y)**2 * s.dxx(x,y) + \
                          2*s.dx(x,y)*s.dy(x,y)*s.dxy(x,y) + \
                          s.dy(x,y)**2 * s.dyy(x,y)
        elif method is "laplace":
            for x,y in image.size():
                z[x,y] =  s.dxx(x,y) + s.dyy(x,y)
        else:
            for x,y in image.size():
                gx, gy, gxx, gxy, gyy = s.dx(x,y), s.dy(x,y), s.dxx(x,y), s.dxy(x,y), s.dyy(x,y)
                z[x,y] =  gx*gy*(gyy - gxx) + gxy*(gx**2 - gy**2)
        self.z = SplineImageView3(z)
        self.regions = transformImage(z, '\l x: x > 0? 1: x<0? -1: 0')
        self.m = method

    def _addFacetIntersection(self, facetX, facetY, point):
        facetIndex = facetX + 10000*facetY
        if self.facets.has_key(facetIndex):
            self.facets[facetIndex].append(point)
        else:
            self.facets[facetIndex] = [point]

    def points(self, useConstraint = True, t = 1e-7):
        result = []
        self.facets = {}

        if not useConstraint:
            def checkConstraint(z, i, xx, yy, t):
                return True
        elif self.m == "haralick":
            checkConstraint = haralickConstraint
        elif self.m in "laplace":
            checkConstraint = laplaceConstraint
        else:
            checkConstraint = splineConstraint

        for x, y in self.i.size() - Size2D(1, 1):
            c = self.z.coefficients(x, y)
            xPoly = [c[k,0] for k in range(4)]
            yPoly = [c[0,k] for k in range(4)]

            for k in polynomialRealRoots(xPoly):
                if k < 0.0 or k >= 1.0:
                    continue

                p = (x + k, y)
                if not checkConstraint(self.z, self.i, p[0], p[1], t):
                    continue
                result.append(p)

                self._addFacetIntersection(x, y,   p)
                self._addFacetIntersection(x, y-1, p)

            for k in polynomialRealRoots(yPoly):
                if k < 0.0 or k >= 1.0:
                    continue

                p = (x, y + k)
                if not checkConstraint(self.z, self.i, p[0], p[1], t):
                    continue
                result.append(p)

                self._addFacetIntersection(x,   y, p)
                self._addFacetIntersection(x-1, y, p)

        return result

    def edges(self, useConstraint = True, t = 1e-7):
        try:
            f = self.facets
        except:
            self.points(useConstraint, t)
            f = self.facets

        e = []
        for k in f.values():
            if len(k) == 1:
                continue
            elif len(k) == 2:
                e.append(k)
            else:
                xx = reduce(lambda x,y: x+y[0], k, 0.0) / len(k)
                yy = reduce(lambda x,y: x+y[1], k, 0.0) / len(k)
                for kk in k:
                    e.append([(xx,yy),kk])
        return e

# --------------------------------------------------------------------

from operator import setitem
def distinct(l):
    d = {}
    map(setitem, (d,)*len(l), l, [])
    return d.keys()

def levelEdgesMap(edges, imageSize):
    assert len(edges[0]) >= 2 and len(edges[0][0]) == 2, \
           "mapFromZeroEdges() expects output of ZeroEdges.edges() as parameter!"

    nodePositions = [ep[0] for ep in edges]
    nodePositions.extend([ep[1] for ep in edges])

    nodePositions = distinct(nodePositions)
    nodePositions.insert(0, None)

    edgeTriples = [(nodePositions.index(ep[0]), nodePositions.index(ep[1]),
                    [Vector2(*ep[0]), Vector2(*ep[1])])
                   for ep in edges]
    edgeTriples.insert(0, None)

    nodePositions = [np and Vector2(*np) for np in nodePositions]

    result = Map(nodePositions, edgeTriples, imageSize,
                 performEdgeSplits = False)
    removeCruft(result, 2)
    return result

### USAGE EXAMPLE: ###
# ze = ZeroEdges(transformImage(phi, '\l x: x + %s' % level), "direct")
# ee = ze.edges(False)
# levelMap = levelEdgesMap(ee, phi.size())

# --------------------------------------------------------------------

from vigra import addPathFromHere
addPathFromHere("../evaluation/")
import edgedetectors
from maputils import removeCruft

def levelSetMap(image, threshold, sigma = None):
    ed = edgedetectors.EdgeDetector(
        bi = "Thresholding", s1 = sigma, nonmax = "zerosSubPixel",
        threshold = threshold)
    result, _, _ = ed.computeMap(image)
    removeCruft(result, 2)
    return result

__all__ = ["levelSetMap"]
