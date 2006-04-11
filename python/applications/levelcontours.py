def tangentDir(s, x0, y0):
    dx, dy = -s.dy(x0, y0), s.dx(x0, y0)
    norm = sqrt(sq(dx) + sq(dy))
    dx /= norm
    dy /= norm
    return dx, dy

def predictorStep1(s, x0, y0, h):
    dx, dy = -s.dy(x0, y0), s.dx(x0, y0)
    norm = sqrt(sq(dx) + sq(dy))
    dx /= norm
    dy /= norm

    x1 = x0 + h*dx
    y1 = y0 + h*dy
    if not isInside(s, x1, y1):
        return None
    return [x1, y1]

def correctorStep2(s, x0, y0, epsilon = 1e-4, n0 = 3):
    x, y = x0, y0
    for k in range(100):
        value = s(x, y)
        if abs(value) < epsilon:
            contraction = float(k) / n0
            return [x, y], contraction
        dx, dy = s.dx(x, y), s.dy(x, y)
        norm2 = sq(dx) + sq(dy)
        correction = value / norm2
        x -= dx * correction
        y -= dy * correction
        if not isInside(s, x, y):
            return None, None
        diff = sqrt(sq(x - x0) + sq(y - y0))
        if diff > 1.0: # FIXME
            return None, None
    return None, None

def predictorCorrectorStep(s, x, y, h, epsilon):
    while abs(h) > 1e-6: # FIXME
        p1 = predictorStep1(s, x, y, h)
        if not p1:
            h /= 2.0
            continue
        p2, contraction = correctorStep2(s, p1[0], p1[1], epsilon)
        if not p2:
            h /= 2.0
            continue
        # use pow(..., 0.5) in connection with predictorStep1(),
        # use pow(..., 0.33) in connection with predictorStep2(),
        h /= min(2.0, max(0.5, pow(contraction, 0.5)))
        if abs(h) > 0.2:
            h = 0.2
        break
    else:
        return None, None
    return p2, h

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
