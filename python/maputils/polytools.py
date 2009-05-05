import vigra, hourglass, math
from vigra import Vector2, dot
from hourglass import Polygon

# FIXME: computing bsd42044 1.5 SPWS (spline order 2)...
# ...
#   File ".../maputils.py", line 473, in addFlowLinesToMap
#     for cp in polytools.clipPoly(points, clipBox):
#   File ".../polytools.py", line 209, in clipPoly
#     assert outside != prevOutside; prevOutside = outside
# AssertionError

def maxDistIter(polygon, maxDist):
    maxDist2 = vigra.sq(maxDist)
    it = iter(polygon)
    prev = it.next()
    yield prev
    for p in it:
        dist2 = (p-prev).squaredMagnitude()
        if dist2 > maxDist2:
            dist = math.sqrt(dist2)
            segmentCount = int(math.ceil(dist / maxDist))
            segment = p - prev
            for i in range(1, segmentCount):
                yield p + segment*i/segmentCount
        yield p
        
# --------------------------------------------------------------------

LEFT = 1
RIGHT = 2
TOP = 4
BOTTOM = 8

def clipPoly(polygon, clipRect, closeAtBorder = None):
    """clipPoly(polygon, clipRect)

    Clips away those parts of polygon which are not in clipRect.
    Returns a list of polygons (since the polygon may leave clipRect,
    enter again, leave, ...).  Polygon segments crossing clipRect's
    borders are cut, such that the resulting polyons get new endpoints
    exactly on the border."""
    
    result = []

#     print "clipPoly(%s..%s)" % (clipRect.begin(), clipRect.end())
#     print list(polygon)

    if closeAtBorder is None:
        closeAtBorder = (polygon[0] == polygon[-1])

    x1, y1 = clipRect.begin()
    x2, y2 = clipRect.end()

    part = None
    startBorder = None
    parts = []

    relPos = None
    for i, p in enumerate(polygon):
        prevRP = relPos
        relPos = 0
        if p[0] < x1:
            relPos |= LEFT
        elif p[0] > x2:
            relPos |= RIGHT
        if p[1] < y1:
            relPos |= TOP
        elif p[1] > y2:
            relPos |= BOTTOM

        if relPos: # outside
            if not i: # incomplete first segment
                continue

            if prevRP & relPos:
                # complete segment outside
                continue

            # calculate leaving intersection
            diff = polygon[i-1] - p
            l = -1.0
            if relPos & LEFT:
                l = max(l, (x1 - p[0]) / diff[0])
                endBorder = LEFT
            if relPos & RIGHT:
                l = max(l, (x2 - p[0]) / diff[0])
                endBorder = RIGHT
            if relPos & TOP:
                nl = (y1 - p[1]) / diff[1]
                if nl > l:
                    l = nl
                    endBorder = TOP
            if relPos & BOTTOM:
                nl = (y2 - p[1]) / diff[1]
                if nl > l:
                    l = nl
                    endBorder = BOTTOM
            ip = p + l * diff

            if prevRP:
                # segment may cross cliprect, calc. start intersection
                pl = 2.0
                if prevRP & LEFT:
                    pl = min(pl, (x1 - p[0]) / diff[0])
                    startBorder = LEFT
                if prevRP & RIGHT:
                    pl = min(pl, (x2 - p[0]) / diff[0])
                    startBorder = RIGHT
                if prevRP & TOP:
                    npl = (y1 - p[1]) / diff[1]
                    if npl < pl:
                        pl = npl
                        startBorder = TOP
                if prevRP & BOTTOM:
                    npl = (y2 - p[1]) / diff[1]
                    if npl < pl:
                        pl = npl
                        startBorder = BOTTOM

                if pl <= l:
                    # we never crossed the clipRect
                    continue

                pip = p + pl * diff
                part = Polygon([pip, ip])
            else:
                part.append(ip)

            if part.length():
                parts.append((startBorder, part, endBorder))
            part = None
            continue

        if not part:
            part = Polygon()
            if i:
                # calculate entering intersection:
                diff = polygon[i-1] - p
                l = 2.0
                if prevRP & LEFT:
                    l = min(l, (x1 - p[0]) / diff[0])
                    startBorder = LEFT
                if prevRP & RIGHT:
                    l = min(l, (x2 - p[0]) / diff[0])
                    startBorder = RIGHT
                if prevRP & TOP:
                    nl = (y1 - p[1]) / diff[1]
                    if nl < l:
                        l = nl
                        startBorder = TOP
                if prevRP & BOTTOM:
                    nl = (y2 - p[1]) / diff[1]
                    if nl < l:
                        l = nl
                        startBorder = BOTTOM
                ip = p + l * diff
                part.append(ip)

        part.append(p)

    if part and part.length():
        parts.append((startBorder, part, None))

    if not parts:
        return []

    if not polygon.closed():
        return [p[1] for p in parts]

    # if polygon[0] (== polygon[-1]) is inside clipRect, we may
    # need to join the first and last part here:
    if parts[0][1][0] == parts[-1][1][-1]:
        assert parts[0][0] is None and parts[-1][-1] is None
        # polygon is entirely within clipRect:
        if len(parts) == 1:
            return [parts[0][1]]
        parts[-1][1].extend(parts[0][1])
        parts[0] = (parts[-1][0], parts[-1][1], parts[0][2])
        del parts[-1]

    if not closeAtBorder:
        return [p[1] for p in parts]

    # compose counterclockwise list of intersection points at clip border:
    sides = (
        ([(-p[1][-1][0], p[1], True ) for p in parts if p[2] == TOP] +
         [(-p[1][ 0][0], p[1], False) for p in parts if p[0] == TOP]),
        ([( p[1][-1][1], p[1], True ) for p in parts if p[2] == LEFT] +
         [( p[1][ 0][1], p[1], False) for p in parts if p[0] == LEFT]),
        ([( p[1][-1][0], p[1], True ) for p in parts if p[2] == BOTTOM] +
         [( p[1][ 0][0], p[1], False) for p in parts if p[0] == BOTTOM]),
        ([(-p[1][-1][1], p[1], True ) for p in parts if p[2] == RIGHT] +
         [(-p[1][ 0][1], p[1], False) for p in parts if p[0] == RIGHT]))

    # counterclockwise list of corner positions:
    corners = (clipRect.begin(),
               clipRect.begin()+(0, clipRect.size()[1]),
               clipRect.end(),
               clipRect.begin()+(clipRect.size()[0], 0))

    isCCW = polygon.partialArea() > 0

    # bookkeeping about mergings (always use the most current polygon)
    merged = {}
    def mergeRoot(poly):
        while True:
            result = merged.get(poly, poly)
            if result is poly:
                break
            poly = result
        return result
    
    lastPoly = None
    prevPoly = None
    prevOutside = None

    for side, end in zip(sides, corners):
        for _, poly, outside in sorted(side):
#            assert outside != prevOutside; prevOutside = outside
            if outside == isCCW:
                prevPoly = poly
            else:
                if prevPoly == None:
                    lastPoly = poly
                    continue
                prevPoly = mergeRoot(prevPoly)
                if prevPoly == poly:
                    poly.append(poly[0])
                    result.append(poly)
                else:
                    prevPoly.extend(poly)
                    merged[poly] = prevPoly
                prevPoly = None

        if prevPoly:
            mergeRoot(prevPoly).append(end)

    if lastPoly:
        lastPoly.append(lastPoly[0])
        if lastPoly.length():
            result.append(lastPoly)

    return result

# --------------------------------------------------------------------

class Line(object):
    def __init__(self, norm, dist):
        self.norm = norm
        self.dist = dist

    def isParallel(self, other):
        return abs(dot(self.norm, other.norm)) == 1.0
    
    def intersect(self, other):
        assert not self.isParallel(other)
        if abs(self.norm[0]) > abs(other.norm[0]):
            a, b = self, other
        else:
            a, b = other, self
        top    = (a.norm/a.norm[0],        a.dist/a.norm[0])
        bottom = (b.norm-top[0]*b.norm[0], b.dist-top[1]*b.norm[0])
        y = bottom[1]/bottom[0][1]
        x = top[1]-y*top[0][1]
        return Vector2(x, y)

    def dir(self):
        """Return direction vector."""
        return Vector2(self.norm[1], -self.norm[0])

    def orthogonalDistance(self, point):
        """Return distance between given point and this line"""
        return dot(point, self.norm) - self.dist

    def point(self, l = 0):
        """Return point on line.  For l == 0 (default), this will be
        the closest point to the origin.  l moves on the line."""
        return self.dist * self.norm + self.dir() * l

class LineSegment(object):
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
    
    def dir(self):
        result = self.p2 - self.p1
        result /= result.magnitude()
        return result
    
    def norm(self):
        d = self.dir()
        return Vector2(-d[1], d[0])
    
    def dist(self):
        """distance to origin"""
        return dot(self.p1, self.norm())
    
    def line(self):
        return Line(self.norm(), self.dist())

def polyLineSegment(poly, index):
    return LineSegment(poly[index % len(poly)], poly[(index+1) % len(poly)])

def polyLineSegments(poly):
    for i in range(len(poly)-1):
        yield polyLineSegment(poly, i)

def shrinkPoly(poly, offset):
    assert poly[0] == poly[-1], "polygon should be closed"
    lines = [seg.line() for seg in polyLineSegments(poly)]
    for line in lines:
        line.dist -= offset
    i = 1
    while i < len(lines):
        if lines[i-1].isParallel(lines[i]):
            del lines[i]
        else:
            i += 1
    if lines[-1].isParallel(lines[0]):
        del lines[-1]
    result = Polygon([lines[i].intersect(lines[(i+1)%len(lines)])
                      for i in range(len(lines))])
    result.append(result[0])
    return result

def rotatePoly(poly, angle):
    """Rotate polygon by the given angle around the origin."""
    
    unitX = Vector2(math.cos(angle), -math.sin(angle))
    unitY = Vector2(math.sin(angle),  math.cos(angle))
    result = Polygon()
    for point in poly:
        result.append(Vector2(dot(point, unitX), dot(point, unitY)))
    return result

def subsetDigitization(poly, shift = None, size = None):
    """Sample poly with a regular grid at integer coordinates starting
    from (0,0) to the given size (which should be a Size2D object)."""

    if size == None:
        size = poly.boundingBox().size()
        size = (int(math.ceil(size[0]))+2,
                int(math.ceil(size[1]))+2)
        if not shift:
            shift = Vector2(0, 0)
        shift = shift + Vector2(1, 1) - poly.boundingBox().begin()
        poly = Polygon(poly + shift)

    result = vigra.GrayImage(size)
    for p in vigra.meshIter(size):
        result[p] = poly.contains(Vector2(p[0], p[1])) and 1 or 0
    return result

# --------------------------------------------------------------------

def smallestBoundingBox(ch):
    """Determine rotated bbox from convex hull"""
    # FIXME: use rotating calipers for O(N) instead of O(N^2)!
    assert ch.closed()
    bboxes = []
    for seg in polyLineSegments(ch):
        line = seg.line()
        norm = line.norm
        dir = line.dir()
        dists = []
        positions = []
        for p in ch:
            dists.append(dot(norm, p))
            positions.append(dot(dir, p))
        l1 = min(positions)
        l2 = max(positions)
        l3 = min(dists)
        l4 = max(dists)
        area = (l2 - l1) * (l4 - l3)
        bboxes.append((area, line, l1, l2, l3, l4))
    bboxes.sort()
    _, line, l1, l2, l3, l4 = bboxes[0]
    p1 = l1 * line.dir() + l3 * line.norm
    p2 = l1 * line.dir() + l4 * line.norm
    p3 = l2 * line.dir() + l4 * line.norm
    p4 = l2 * line.dir() + l3 * line.norm
    return Polygon([p1, p2, p3, p4, p1])

# --------------------------------------------------------------------

if __name__ == "__main__":
    import fig, hourglass
    f = fig.File("cliptest.fig")
    cr = hourglass.BoundingBox((0, 0), (4500, 4500))
    f.layer(1).remove()
    for o in f.findObjects(type = fig.PolylineBase, depth = 42):
        p = Polygon(o.points)
        if o.closed():
            p.append(p[0])
        pp = clipPoly(p, cr)
        for p in pp:
            no = fig.Polygon(p, p[0] == p[-1])
            no.depth = 1
            no.lineWidth = 3
            if no.closed():
                no.fillStyle = fig.fillStyleSolid
                no.fillColor = f.getColor(0.5)
            else:
                no.forwardArrow = fig.Arrow()
            f.append(no)
    f.save(fig2dev = "eps")

