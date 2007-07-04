from vigra import Vector2, dot
from hourglass import Polygon

def _intersectLine(inside, outside, clipRect):
    if outside[1] > clipRect.end()[1]:
        outside = inside + (outside-inside) * \
                 (clipRect.end()[1]-inside[1])/(outside[1]-inside[1])
    elif outside[1] < clipRect.begin()[1]:
        outside = inside + (outside-inside) * \
                 (clipRect.begin()[1]-inside[1])/(outside[1]-inside[1])
    if outside[0] > clipRect.end()[0]:
        return inside + (outside-inside) * \
               (clipRect.end()[0]-inside[0])/(outside[0]-inside[0])
    elif outside[0] < clipRect.begin()[0]:
        return inside + (outside-inside) * \
               (clipRect.begin()[0]-inside[0])/(outside[0]-inside[0])
    return outside

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

    i = 0
    while i < len(polygon):
        while i < len(polygon) and not clipRect.contains(polygon[i]):
            i += 1
        if i >= len(polygon):
            break

        p = Polygon()

        if i > 0 and i < len(polygon):
            p.append(_intersectLine(polygon[i], polygon[i-1], clipRect))

        while i < len(polygon) and clipRect.contains(polygon[i]):
            p.append(polygon[i])
            i += 1

        if i < len(polygon):
            p.append(_intersectLine(polygon[i-1], polygon[i], clipRect))

        result.append(p)

    if closed(polygon):
        if len(result) > 1 and (result[0][0] == result[-1][-1]):
            result[-1].extend(result[0])
            del result[0]

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
