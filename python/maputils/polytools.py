from vigra import Vector2, dot
from hourglass import Polygon

LEFT = 1
RIGHT = 2
TOP = 4
BOTTOM = 8

def _intersectLine(inside, outside, clipRect):
    if outside[1] > clipRect.end()[1]:
        border = BOTTOM
        outside = inside + (outside-inside) * \
                 (clipRect.end()[1]-inside[1])/(outside[1]-inside[1])
    elif outside[1] < clipRect.begin()[1]:
        border = TOP
        outside = inside + (outside-inside) * \
                 (clipRect.begin()[1]-inside[1])/(outside[1]-inside[1])
    if outside[0] > clipRect.end()[0]:
        border = RIGHT
        outside = inside + (outside-inside) * \
               (clipRect.end()[0]-inside[0])/(outside[0]-inside[0])
    elif outside[0] < clipRect.begin()[0]:
        border = LEFT
        outside = inside + (outside-inside) * \
               (clipRect.begin()[0]-inside[0])/(outside[0]-inside[0])
    return border, outside

def clipPoly(polygon, clipRect):
    """clipPoly(polygon, clipRect)

    Clips away those parts of polygon which are not in clipRect.
    Returns a list of polygons (since the polygon may leave clipRect,
    enter again, leave, ...).  Polygon segments crossing clipRect's
    borders are cut, such that the resulting polyons get new endpoints
    exactly on the border."""
    
    result = []

#     print "clipPoly(%s..%s)" % (clipRect.begin(), clipRect.end())
#     print list(polygon)

    def closed(polygon):
        return polygon[0] == polygon[-1]

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
            if not part:
                continue
            # close current part:
            endBorder, ip = _intersectLine(polygon[i-1], p, clipRect)
            part.append(ip)
            parts.append((startBorder, part, endBorder))
            part = None
            continue

        if not part:
            part = Polygon()
            if i:
                startBorder, ip = _intersectLine(p, polygon[i-1], clipRect)
                part.append(ip)

        part.append(p)

    if part:
        parts.append((startBorder, part, None))

    if not closed(polygon):
        result = [p[1] for p in parts]
    else:
        if len(parts) > 1 and (parts[0][1][0] == parts[-1][1][-1]):
            assert parts[0][0] is None and parts[-1][-1] is None
            parts[-1][1].extend(parts[0][1])
            parts[0] = (parts[-1][0], parts[-1][1], parts[0][2])
            del parts[-1]

        isCCW = polygon.partialArea() > 0
        merged = {}
        def mergeRoot(poly):
            while True:
                result = merged.get(poly, poly)
                if result is poly:
                    break
                poly = result
            return result
#             while poly in merged:
#                 poly = merged[poly]
#             return poly
        
        lastPoly = None
        prevPoly = None
        prevOutside = None

        sides = (
            ([(-p[1][-1][0], p[1], True ) for p in parts if p[2] == TOP] +
             [(-p[1][ 0][0], p[1], False) for p in parts if p[0] == TOP]),
            ([( p[1][-1][1], p[1], True ) for p in parts if p[2] == LEFT] +
             [( p[1][ 0][1], p[1], False) for p in parts if p[0] == LEFT]),
            ([( p[1][-1][0], p[1], True ) for p in parts if p[2] == BOTTOM] +
             [( p[1][ 0][0], p[1], False) for p in parts if p[0] == BOTTOM]),
            ([(-p[1][-1][1], p[1], True ) for p in parts if p[2] == RIGHT] +
             [(-p[1][ 0][1], p[1], False) for p in parts if p[0] == RIGHT]))
        corners = (clipRect.begin(),
                   clipRect.begin()+(0, clipRect.size()[1]),
                   clipRect.end(),
                   clipRect.begin()+(clipRect.size()[0], 0))

        for side, end in zip(sides, corners):
            for _, poly, outside in sorted(side):
                assert outside != prevOutside; prevOutside = outside
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

        if prevPoly:
            assert lastPoly
            #prevPoly.extend(lastPoly)
            prevPoly.append(prevPoly[0])
            result.append(prevPoly)
        else:
            assert not lastPoly

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
        return dot(self.p1, self.norm())
    
    def line(self):
        return Line(self.norm(), self.dist())

def polyLineSegment(poly, index):
    return LineSegment(poly[index % len(poly)], poly[(index+1) % len(poly)])

def shrinkPoly(poly, offset):
    assert poly[0] == poly[-1], "polygon should be closed"
    sides = len(poly)-1
    lines = [polyLineSegment(poly, i).line() for i in range(sides)]
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

