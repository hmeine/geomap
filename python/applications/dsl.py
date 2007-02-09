import copy
from vigra import Point2D, Vector2, Rational

def freeman(crackEdge):
    """Returns freeman codes for each polyline segment in the given
    crackEdge polygon.  Uses mathematical orientation, i.e. the code
    meaning is:
    
    0: diff = ( 1,  0)
    1: diff = ( 0,  1)
    2: diff = (-1,  0)
    3: diff = ( 0, -1)
    """
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

freeman2Diff8Conn    = [Point2D( 1, 0), Point2D( 1,  1)]
freeman2Diff8ConnNeg = [Point2D(-1, 0), Point2D(-1, -1)]
freeman2Diff4Conn    = [Point2D( 1, 0), Point2D( 0,  1)]

def quadrant(c1, c2):
    if c2 < c1:
        c1, c2 = c2, c1
    if c1 > 0:
        return c1
    if c2 > 1:
        return 3
    return 0

def forwardIter(list, index, loop):
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
    i = index
    while i >= 0:
        yield list[i]
        i -= 1
    if loop:
        i = len(list) - 1
        while i >= 0:
            yield list[i]
            i -= 1

def searchForwardQuadrant(freemanCodes, index, closed = True):
    try:
        it = forwardIter(freemanCodes, index, closed)
        c1 = it.next()
        while True:
            c2 = it.next()
            if c2 != c1:
                return (c1, c2)
    except StopIteration:
        return c1

def searchBackwardQuadrant(freemanCodes, index, closed = True):
    try:
        it = backwardIter(freemanCodes, index, closed)
        c1 = it.next() # skip first
        c1 = it.next()
        while True:
            c2 = it.next()
            if c2 != c1:
                return (c1, c2)
    except StopIteration:
        return c1

def searchTangentQuadrant(freemanCodes, index, closed = True):
    try:
        it1 = forwardIter(freemanCodes, index, closed)
        it2 = backwardIter(freemanCodes, index, closed)
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
    
    def plotEquation(self, leaningType = 0):
        return ("%s*x + (%s)" % (self.slope(), self.axisIntercept(leaningType))).replace("/", "./")
    
    def plotItems(self):
        return [Gnuplot.Func(self.plotEquation(), title = "center line"),
                Gnuplot.Func(self.plotEquation(1), title = "lower leaning line"),
                Gnuplot.Func(self.plotEquation(2), title = "upper leaning line")]
    
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
        return "DigitalStraightLine(%d, %d, %d)" % (self.a, self.b, self.pos)
    
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
    return [Gnuplot.Func(self.plotEquation(), title = "center line"),
            Gnuplot.Func(self.plotEquation(LeaningType.LowerLeaningLine), title = "lower leaning line"),
            Gnuplot.Func(self.plotEquation(LeaningType.UpperLeaningLine), title = "upper leaning line")]

def DigitalStraightLine8__repr__(self):
    return "DigitalStraightLine8(%d, %d, %d)" % (self.a, self.b, self.pos)

def DigitalStraightLine4__repr__(self):
    return "DigitalStraightLine4(%d, %d, %d)" % (self.a, self.b, self.pos)

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
    for fc in freemanIter:
        if not fc in allowed:
            return
        cur += freeman2Diff[fc % 2]
        yield copy.copy(cur)

def forwardDSL(freemanCodes, index, closed, allowed = None):
    if allowed == None:
        allowed = searchForwardQuadrant(freemanCodes, index, closed)
    fmi = forwardIter(freemanCodes, index, closed)
    dsl = DigitalStraightLine(freemanCodes[index] % 2 and 1 or 0, 1, 0)
    for point in originatingPolyIter(fmi, allowed, freeman2Diff8Conn):
        if not dsl.addPoint(*point):
            break
    return dsl, fmi.gi_frame.f_locals["i"]-index

def backwardDSL(freemanCodes, index, closed, allowed = None):
    if allowed == None:
        allowed = searchBackwardQuadrant(freemanCodes, index, closed)
    fmi = backwardIter(freemanCodes, index, closed)
    fmi.next()
    dsl = DigitalStraightLine(freemanCodes[index] % 2 and 1 or 0, 1, 0)
    for point in originatingPolyIter(fmi, allowed, freeman2Diff8ConnNeg):
        if not dsl.addPoint(*point):
            break
    return dsl, -(fmi.gi_frame.f_locals["i"]-index)

def tangentDSL(freemanCodes, index, closed, allowed = None):
    if allowed == None:
        allowed = searchTangentQuadrant(freemanCodes, index, closed)
        if type(allowed) != tuple:
            return DigitalStraightLine(0, 1, 0), None
    dsl = DigitalStraightLine(freemanCodes[index] % 2 and 1 or 0, 1, 0)
    result = dsl
    ffmi = forwardIter(freemanCodes, index, closed)
    bfmi = backwardIter(freemanCodes, index, closed)
    bfmi.next()
    origin = Point2D(0, 0)
    bopi = originatingPolyIter(bfmi, allowed, freeman2Diff8ConnNeg)
    for point1 in originatingPolyIter(ffmi, allowed, freeman2Diff8Conn):
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
    return result, ffmi.gi_frame.f_locals["i"]-index

def offset(freemanCodes, index, closed = True):
    dsl, ofs = hourglass.tangentDSL(freemanCodes, index, closed)

    fc1 = freemanCodes[index - ofs]
    for i in range(index - ofs + 1, index + ofs):
        fc2 = freemanCodes[i]
        if fc2 != fc1:
            break

    alpha = (2.*dsl.pos+dsl.b-1)/(2*dsl.b)
    q = quadrant(fc1, fc2)
    if q == 0:
        return Vector2( alpha, -alpha)
    elif q == 1:
        return Vector2(-alpha, -alpha)
    elif q == 2:
        return Vector2(-alpha,  alpha)
    else:
        return Vector2( alpha,  alpha)

def offset2(freemanCodes, index, closed = True):
    dsl, ofs = hourglass.tangentDSL(freemanCodes, index, closed)
    dsl = dsl.convertToFourConnected()

    fc1 = freemanCodes[index - ofs]
    for i in range(index - ofs + 1, index + ofs):
        fc2 = freemanCodes[i]
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

    cp = Rational(2*dsl.pos+dsl.width()-1, 2*(math.sq(dsl.a) + math.sq(dsl.b)))
    return Vector2(float(dsl.a*cp), float(-dsl.b*cp))

import Gnuplot

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
        self.freeman2Diff = reverse and freeman2Diff8ConnNeg or freeman2Diff8Conn
    
    def __call__(self, code):
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
            return self.plot()
        return False
    
    def plot(self):
        result = True
        for point in self.points:
            if not self.dsl.contains(point[0], point[1]):
                print "%s lost (no longer in DSL!)" % point
                result = False
        self.g.plot(self.points, *self.dsl.plotItems())
        return result

if __name__ == "__main__":
    dsl = DigitalStraightLine(0, 1, 0)
    dsl.addPoint( 1,  0)
    dsl.addPoint(-1,  0)
    dsl.addPoint(-2, -1)

    gimg = GrayImage(50, 50)
    for p in gimg.size():
        gimg[p] = norm(p - Point2D(25, 25)) < 20 and 255 or 0

    import pixelmap
    cem = pixelmap.crackEdgeMap(gimg, False)
    crackPoly = cem.edge(2)
    fc = freeman(crackPoly)

    import Gnuplot
    def gpLine(points, with = "lines", **kwargs):
        return Gnuplot.Data(points, with = with, **kwargs)

    g = Gnuplot.Gnuplot()
    g("set size ratio -1")

    ep = [p + offset(fc, i) for i, p in enumerate(list(crackPoly)[:-1])]
    ep.append(ep[0])
    ep2 = [p + offset2(fc, i) for i, p in enumerate(list(crackPoly)[:-1])]
    ep2.append(ep2[0])
    g.plot(gpLine(crackPoly), gpLine(ep), gpLine(ep2))

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
