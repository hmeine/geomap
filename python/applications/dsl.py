import copy
from vigra import Point2D, Vector2, Rational, sq

def freeman(crackEdge):
    """freeman(crackEdge) -> list

    Returns freeman codes for each polyline segment in the given
    crackEdge polygon.  Uses mathematical orientation, i.e. the code
    meaning is:
    
    0: diff = ( 1,  0)
    1: diff = ( 0,  1)
    2: diff = (-1,  0)
    3: diff = ( 0, -1)

    Effect unspecified if crackEdge is not a crack edge polygon."""
    
    result = [None] * (len(crackEdge)-1)
    it = iter(crackEdge)
    prev = it.next()
    i = 0
    for point in it:
        diff = point - prev
        if diff[0]:
            if diff[0] > 0:
                result[i] = 0
            else:
                result[i] = 2
        elif diff[1] < 0:
            result[i] = 3
        else:
            result[i] = 1
        i += 1
        prev = point
    return result

freeman2Diff8ConnFirstQuadrant    = [Point2D( 1, 0), Point2D( 1,  1)]
freeman2Diff8ConnThirdQuadrant = [Point2D(-1, 0), Point2D(-1, -1)]
freeman2Diff4Conn    = [Point2D( 1, 0), Point2D( 0,  1),
                        Point2D(-1, 0), Point2D( 0, -1)]

def forwardIter(list, index, loop):
    if loop:
        index = index % len(list)
    i = index
    while i < len(list):
        yield list[i]
        i += 1
    if loop:
        i = 0
        while i < index:
            yield list[i]
            i += 1

def backwardIter(list, index, loop):
    if loop:
        index = index % len(list)
    if index < 0:
        index += len(list)
    i = index
    while i >= 0:
        yield list[i]
        i -= 1
    if loop:
        i = len(list) - 1
        while i >= 0:
            yield list[i]
            i -= 1

def searchForwardQuadrant(freemanCodes, pointIndex, closed = True):
    try:
        it = forwardIter(freemanCodes, pointIndex, closed)
        c1 = it.next()
        while True:
            c2 = it.next()
            if c2 != c1:
                return (c1, c2)
    except StopIteration:
        return c1

def searchBackwardQuadrant(freemanCodes, pointIndex, closed = True):
    try:
        it = backwardIter(freemanCodes, pointIndex-1, closed)
        c1 = it.next()
        while True:
            c2 = it.next()
            if c2 != c1:
                return (c1, c2)
    except StopIteration:
        return c1

def searchTangentQuadrant(freemanCodes, pointIndex, closed = True):
    try:
        it1 = forwardIter(freemanCodes, pointIndex, closed)
        it2 = backwardIter(freemanCodes, pointIndex, closed)
        c1 = it2.next()
        while True:
            c2 = it1.next()
            c3 = it2.next()
            if c2 != c1:
                if c3 != c2 and c3 != c1:
                    return c1
                return (c1, c2)
            if c3 != c1:
                return (c1, c3)
    except StopIteration:
        return c1

def straightPoints(poly):
    """Return list of points where forward- and backwards-quadrants
    are different."""
    closed = poly[0] == poly[-1]
    freemanCodes = freeman(poly)
    result = []
    for i in range(len(poly)-1):
        fc = searchForwardQuadrant(freemanCodes, i, closed)
        bc = searchBackwardQuadrant(freemanCodes, i, closed)
        if fc != bc:
            result.append(poly[i])
    return result

class PyDigitalStraightLine(object):
    def __init__(self, a, b, pos):
        self.a = a
        self.b = b
        self.pos = pos
        self.is8Connected = True
    
    def slope(self):
        return Rational(self.a, self.b)
    
    def axisIntercept(self, leaningType = 0):
        """dsl.axisIntercept(leaningType = 0)
        
        leaningType means:
        0: center line
        1: lower leaning line
        2: upper leaning line"""
        
        pos = Rational(self.pos, 1)
        if leaningType == 0:
            pos += Rational(self.width()-1, 2)
        elif leaningType == 1:
            pos += self.width()-1
        return -pos / self.b
    
    def __call__(self, x, y):
        return self.a*x - self.b*y
    
    def __getinitargs__(self):
        return (self.a, self.b, self.pos)
    
    def contains(self, x, y):
        v = self(x, y) - self.pos
        return 0 <= v < self.width()
    
    def width(self):
        if self.is8Connected:
            return max(abs(self.a), abs(self.b))
        else:
            return abs(self.a) + abs(self.b)
    
    def __repr__(self):
        return "PyDigitalStraightLine(%d, %d, %d)" % (self.a, self.b, self.pos)
    
    def addPoint(self, x, y):
        """works only for 8-connected lines in 1st octant"""
        assert self.is8Connected and (0 <= self.a <= self.b)

        v = self(x, y) - self.pos
        width = self.width()
        if 0 <= v < width:
            #print "point already inside:", self, point
            return True # point is within DSL
        
        above = True
        if v == -1:
            #print "point above"
            pass
        elif v == width:
            #print "point below"
            above = False
        else:
            #print point, "cannot be added to", self, self.contains(point)
            return False

        increase = above
        pos = self.pos
        ax = x
        ay = y
        if x < 0:
            ax = -ax
            ay = -ay
            increase = not increase
            pos = 1-self.width()-self.pos

        k = 0
        if increase:
            #print "slope increases:", self, point
            while k < width:
                if (self.a*k-pos) % width == 0:
                    break
                k += 1
            assert k < width
            l = (self.a*k-pos) / width
        else:
            #print "slope decreases:", self, point
            while k < width:
                if (self.a*k-pos-width+1) % width == 0:
                    break
                k += 1
            assert k < width
            l = (self.a*k-pos-width+1) / width

        self.a = ay - l
        self.b = ax - k

        if above:
            # ensure new point is on lower leaning line:
            self.pos = self(x, y)
        else:
            # ensure new point is on upper leaning line:
            self.pos = self(x, y) - self.width() + 1
        
        if not self.contains(x, y):
            #print self, v, point
            assert self.contains(x, y), \
                   "post-condition: addPoint should lead to contains"

        #print "->", self
        return True

    def mirrorX(self):
        self.a = -self.a
    
    def mirrorY(self):
        self.mirrorX()
        self.mirrorXY()
    
    def mirrorXY(self):
        self.pos = 1-self.width()-self.pos
    
    def convertToFourConnected(self):
        assert self.is8Connected, "DSL should not be converted twice"
        self.b = self.b - self.a
        self.is8Connected = False
        return self

import hourglass
from hourglass import LeaningType

def DigitalStraightLine_plotEquation(self, leaningType = LeaningType.CenterLine):
    return ("%s*x + (%s)" % (self.slope(), self.axisIntercept(leaningType))).replace("/", "./")

def DigitalStraightLine_plotItems(self):
    return [Gnuplot.Func(self.plotEquation(), title = "center line", with_ = "l 1"),
            Gnuplot.Func(self.plotEquation(LeaningType.LowerLeaningLine), title = "lower leaning line", with_ = "l 4"),
            Gnuplot.Func(self.plotEquation(LeaningType.UpperLeaningLine), title = "upper leaning line", with_ = "l 5")]

def DigitalStraightLine8__repr__(self):
    return "DigitalStraightLine8(%d, %d, %d)" % (self.a, self.b, self.pos)

def DigitalStraightLine4__repr__(self):
    return "DigitalStraightLine4(%d, %d, %d)" % (self.a, self.b, self.pos)

PyDigitalStraightLine.plotEquation = DigitalStraightLine_plotEquation
PyDigitalStraightLine.plotItems = DigitalStraightLine_plotItems
hourglass.DigitalStraightLine8.plotEquation = DigitalStraightLine_plotEquation
hourglass.DigitalStraightLine8.plotItems = DigitalStraightLine_plotItems
hourglass.DigitalStraightLine8.__repr__ = DigitalStraightLine8__repr__
hourglass.DigitalStraightLine4.plotEquation = DigitalStraightLine_plotEquation
hourglass.DigitalStraightLine4.plotItems = DigitalStraightLine_plotItems
hourglass.DigitalStraightLine4.__repr__ = DigitalStraightLine4__repr__

DigitalStraightLine = hourglass.DigitalStraightLine8
#DigitalStraightLine = PyDigitalStraightLine

def originatingPolyIter(freemanIter, allowed, freeman2Diff):
    cur = Point2D(0, 0)
    count = len(freeman2Diff)
    for fc in freemanIter:
        if not fc in allowed:
            return
        cur += freeman2Diff[fc % count]
        yield copy.copy(cur)

def forwardDSL(freemanCodes, pointIndex, closed, allowed = None):
    assert closed or pointIndex < len(freemanCodes)
    if allowed == None:
        allowed = searchForwardQuadrant(freemanCodes, pointIndex, closed)
    fmi = forwardIter(freemanCodes, pointIndex, closed)
    dsl = DigitalStraightLine(freemanCodes[pointIndex] % 2 and 1 or 0, 1, 0)
    for point in originatingPolyIter(fmi, allowed, freeman2Diff8ConnFirstQuadrant):
        if not dsl.addPoint(*point):
            break
    return dsl, fmi.gi_frame.f_locals["i"]-pointIndex # FIXME: hack & wrong if closed

def backwardDSL(freemanCodes, pointIndex, closed, allowed = None):
    assert closed or pointIndex > 0
    if allowed == None:
        allowed = searchBackwardQuadrant(freemanCodes, pointIndex, closed)
    fmi = backwardIter(freemanCodes, pointIndex-1, closed)
    dsl = DigitalStraightLine(freemanCodes[pointIndex-1] % 2 and 1 or 0, 1, 0)
    for point in originatingPolyIter(fmi, allowed, freeman2Diff8ConnThirdQuadrant):
        if not dsl.addPoint(*point):
            break
    return dsl, -(fmi.gi_frame.f_locals["i"]-pointIndex) # FIXME: hack & wrong if closed

def tangentDSL(freemanCodes, pointIndex, closed, allowed = None):
    assert closed or (pointIndex > 0 and pointIndex < len(freemanCodes))
    if allowed == None:
        allowed = searchTangentQuadrant(freemanCodes, pointIndex, closed)
        if type(allowed) != tuple:
            return DigitalStraightLine(0, 1, 0), None
        print allowed
    if closed:
        pointIndex = pointIndex % len(freemanCodes)
    dsl = DigitalStraightLine(freemanCodes[pointIndex] % 2 and 1 or 0, 1, 0)
    result = dsl
    ffmi = forwardIter(freemanCodes, pointIndex, closed)
    bfmi = backwardIter(freemanCodes, pointIndex-1, closed)
    origin = Point2D(0, 0)
    bopi = originatingPolyIter(bfmi, allowed, freeman2Diff8ConnThirdQuadrant)
    for point1 in originatingPolyIter(ffmi, allowed, freeman2Diff8ConnFirstQuadrant):
        try:
            point2 = bopi.next()
        except StopIteration:
            break
        result = copy.copy(dsl)

        dsl.pos -= dsl(origin[0], origin[1])
        #print "1: adding %s (was %s) to %s" % (point1 - origin, point1, dsl)
        if not dsl.addPoint(point1[0] - origin[0], point1[1] - origin[1]):
            break
        dsl.pos += dsl(origin[0], origin[1])
        origin = point1

        dsl.pos -= dsl(origin[0], origin[1])
        #print "2: adding %s (was %s) to %s" % (point2 - origin, point2, dsl)
        if not dsl.addPoint(point2[0] - origin[0], point2[1] - origin[1]):
            break
        dsl.pos += dsl(origin[0], origin[1])
        origin = point2

        #print "added %s and %s: %s" % (point1, point2, dsl)
    return result, ffmi.gi_frame.f_locals["i"]-pointIndex # FIXME: hack & wrong if closed

def crackEdgeTangent(freemanCodes, pointIndex, closed):
    dsl, ofs = hourglass.tangentDSL(freemanCodes, pointIndex, closed)
    if not ofs:
        return dsl, ofs
    
    dsl = dsl.convertToFourConnected()

    fc1 = freemanCodes[pointIndex - ofs]
    count = len(freemanCodes)
    for crackIndex in range(pointIndex - ofs + 1, pointIndex + ofs):
        fc2 = freemanCodes[crackIndex % count]
        if fc2 != fc1:
            break
    
    q = quadrant(fc1, fc2)
    if q == 0:
        pass
    elif q == 1:
        dsl.mirrorX()
    elif q == 2:
        dsl.mirrorXY()
    else:
        dsl.mirrorY()
    return dsl, ofs

def quadrant(fc1, fc2):
    """quadrant(fc1, fc2) -> 0..3

    Given the two different freeman codes fc1 and fc2, returns the
    number of the respective quadrant."""
    
    if fc2 < fc1:
        fc1, fc2 = fc2, fc1
    if fc1 > 0:
        return fc1
    if fc2 > 1:
        return 3
    return 0

def offset(freemanCodes, pointIndex, closed = True):
    if closed:
        pointIndex = pointIndex % len(freemanCodes)

    dsl, ofs = hourglass.tangentDSL(freemanCodes, pointIndex, closed)
    if not ofs:
        return Vector2(0, 0)

    fc1 = freemanCodes[pointIndex - ofs]
    count = len(freemanCodes)
    for crackIndex in range(pointIndex - ofs + 1, pointIndex + ofs):
        fc2 = freemanCodes[crackIndex % count]
        if fc2 != fc1:
            break

    q = quadrant(fc1, fc2)
    alpha = (2.*dsl.pos+dsl.b-1)/(2*dsl.b)
    if q == 0:
        return Vector2( alpha, -alpha)
    elif q == 1:
        return Vector2(-alpha, -alpha)
    elif q == 2:
        return Vector2(-alpha,  alpha)
    else:
        return Vector2( alpha,  alpha)

def offset2(freemanCodes, pointIndex, closed = True):
    """My own derivation of the necessary offset(), perpendicular to
    the tangent - gives different, but nearly indistinguishable
    results.  The RMSE of the points' radii in the circle example from
    __main__ is 0.4 percent lower, so this is the default / used in
    euclideanPath. ;-)"""
    
    if closed:
        pointIndex = pointIndex % len(freemanCodes)

    dsl, ofs = crackEdgeTangent(freemanCodes, pointIndex, closed)

    cp = Rational(2*dsl.pos+dsl.width()-1, 2*(sq(dsl.a) + sq(dsl.b)))
    return Vector2(float(dsl.a*cp), float(-dsl.b*cp))

from hourglass import Polygon

def euclideanPath(crackPoly, closed = None):
    """Return tangent-driven Euclidean Path for the given `crackPoly`.
    If the polygon is not `closed` (default/None: auto-detect), the
    first and last points will not be changed (no tangent known)."""

    fc = freeman(crackPoly)
    if closed == None:
        closed = crackPoly[-1] == crackPoly[0]
    result = Polygon(crackPoly)
    if closed:
        eo = offset2(fc, 0, closed)
        result[0] += eo
        result[-1] += eo
    for i in range(1, len(result)-1):
        result[i] += offset2(fc, i, closed)
    return result

def crackEdge2EuclideanPath(crackEdge):
    if crackEdge.flag(1): # flag_constants.BORDER_PROTECTION
        return crackEdge
    closed = crackEdge.isLoop() and crackEdge.startNode().hasDegree(2)
    return euclideanPath(crackEdge, closed)

def crackEdges2EuclideanPaths(crackEdgeMap):
    import maputils
    return maputils.copyMap(
        crackEdgeMap, edgeTransform = crackEdge2EuclideanPath)

class DSLExperiment(object):
    def __init__(self, reverse = False):
        g = Gnuplot.Gnuplot()
        g("set xtics 1; set ytics 1; set grid xtics ytics")
        g("set size ratio -1")
        self.g = g
        self.dsl = None
        self.pos = Point2D(0, 0)
        self.points = [self.pos]
        self.code1 = None
        self.code2 = None
        self.freeman2Diff = reverse and freeman2Diff8ConnThirdQuadrant or freeman2Diff8ConnFirstQuadrant
    
    def __call__(self, code, plot = True):
        if self.dsl == None:
            self.dsl = DigitalStraightLine(code % 2 and 1 or 0, 1, 0)
            self.code1 = code
        if code != self.code1:
            if self.code2 == None:
                self.code2 = code
            assert code == self.code2, "more than two different freeman codes passed to DigitalStraightLine!"
        newPos = self.pos + self.freeman2Diff[code % 2]
        if self.dsl.addPoint(*newPos):
            self.pos = newPos
            self.points.append(newPos)
            return self.plot(plot)
        return False
    
    def plot(self, plot = True):
        result = True
        for point in self.points:
            if not self.dsl.contains(point[0], point[1]):
                print "%s lost (no longer in DSL!)" % point
                result = False
        if plot:
            self.g.plot(Gnuplot.Data(self.points,
                                     title = "discrete line", with_ = "lp 3"),
                        *self.dsl.plotItems())
        return result

if __name__ == "__main__":
    dsl = DigitalStraightLine(0, 1, 0)
    dsl.addPoint( 1,  0)
    dsl.addPoint(-1,  0)
    dsl.addPoint(-2, -1)

    from vigra import GrayImage, norm, meshIter
    
    gimg = GrayImage(50, 50)
    for p in meshIter(gimg.size()):
        gimg[p] = norm(p - Point2D(25, 25)) < 20 and 255 or 0

    import pixelmap
    cem = pixelmap.crackEdgeMap(gimg, False)
    crackPoly = cem.edge(2)
    fc = freeman(crackPoly)

    import Gnuplot
    def gpLine(points, with_ = "lines", **kwargs):
        return Gnuplot.Data(points, with_ = with_, **kwargs)

    g = Gnuplot.Gnuplot()
    g("set size ratio -1")

    ep = [p + offset(fc, i) for i, p in enumerate(list(crackPoly)[:-1])]
    ep.append(ep[0])
    ep2 = [p + offset2(fc, i) for i, p in enumerate(list(crackPoly)[:-1])]
    ep2.append(ep2[0])
    g.plot(gpLine(crackPoly), gpLine(ep), gpLine(ep2))

    g2 = Gnuplot.Gnuplot()
    g2.plot([norm(p-Vector2(25,25))-20 for p in ep],
            [norm(p-Vector2(25,25))-20 for p in ep2])

#     tangents = [tangentDSL(fc, i, True)[0] for i in range(len(fc))]
#     g.plot(gpLine(tangentList(ep, 1), "lp"),
#            gpLine([atan2(t.b, t.a)+math.pi/2 for t in tangents], "lp"))

#     e = DSLExperiment(True)
#     for code in [1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1]:
#         if not e(code):
#             break
#         print "code %s -> %s" % (code, e.dsl)

#     index = 43
#     print "debugging %d:" % index
#     dsl, ofs = tangentDSL(fc, index, True)
#     dsl.convertToFourConnected()
#     g.set_range("xrange", (-5, 5))
#     g.set_range("yrange", (-10, 10))
#     g.plot(gpLine(crackPoly + (-crackPoly[index])),
#            Polygon(list(crackPoly)[index-ofs:index+ofs+1]) + (-crackPoly[index]),
#            *(dsl.convertToFourConnected()).plotItems())

#   import mapdisplay
#   d = mapdisplay.MapDisplay(cem, gimg)
