def intersectLine(inside, outside, clipRect):
    if outside[1] > clipRect.end()[1]:
        return inside + (outside-inside) * \
               (clipRect.end()[1]-inside[1])/(outside[1]-inside[1])
    if outside[1] < clipRect.begin()[1]:
        return inside + (outside-inside) * \
               (clipRect.begin()[1]-inside[1])/(outside[1]-inside[1])
    if outside[0] > clipRect.end()[0]:
        return inside + (outside-inside) * \
               (clipRect.end()[0]-inside[0])/(outside[0]-inside[0])
    if outside[0] < clipRect.begin()[0]:
        return inside + (outside-inside) * \
               (clipRect.begin()[0]-inside[0])/(outside[0]-inside[0])

def clipPoly(polygon, clipRect):
    result = []

#     print "clipPoly(%s..%s)" % (clipRect.begin(), clipRect.end())
#     print list(polygon)

    i = 0
    while i < len(polygon):
        while i < len(polygon) and not clipRect.contains(polygon[i]):
            i += 1
        if i >= len(polygon):
            break

        p = Polygon()

        if i > 0 and i < len(polygon):
            p.append(intersectLine(polygon[i], polygon[i-1], clipRect))

        while i < len(polygon) and clipRect.contains(polygon[i]):
            p.append(polygon[i])
            i += 1

        if i < len(polygon):
            p.append(intersectLine(polygon[i-1], polygon[i], clipRect))

        result.append(p)

    return result

_qtColor2figColorMapping = {
    Qt.black : fig.colorBlack,
    Qt.blue : fig.colorBlue,
    Qt.green : fig.colorGreen,
    Qt.cyan : fig.colorCyan,
    Qt.red : fig.colorRed,
    Qt.magenta : fig.colorMagenta,
    Qt.yellow : fig.colorYellow,
    Qt.white : fig.colorWhite,
    }

def qtColor2figColor(color):
    return _qtColor2figColorMapping[color]

class FigExporter:
    def __init__(self, scale = 30, roi = None):
        self.f = fig.File()
        self.scale = scale
        self.roi = roi
        self.offset = Vector2(0.5, 0.5)

    def addROIRect(self, roi = None, depth = 1, **attr):
        assert roi or self.roi, "addROIRect(): no ROI given!?"
        if roi == None:
            roi = BoundingBox(self.roi)
        if self.roi:
            roi.moveBy(-self.roi.begin())
        result = fig.PolyBox(roi.begin()[0] * self.scale,
                             roi.begin()[1] * self.scale,
                             roi.end()[0] * self.scale,
                             roi.end()[1] * self.scale)
        result.depth = depth
        self.f.append(result)
        for a in attr:
            setattr(result, a, attr[a])
        return result

    def addBackgroundWithFrame(self, bgImageFilename, **params):
        bgRect = self.addROIRect(**params)

        bgImage = fig.PictureBBox(0, 0, 1, 1, bgImageFilename)
        bgImage.points = list(bgRect.points)
        bgImage.depth = 100
        self.f.append(bgImage)

    def addPointCircles(self, points, radius, returnIndices = False, **attr):
        radius *= self.scale
        if not attr.has_key("fillStyle"):
            attr["fillStyle"] = fig.fillStyleSolid
        result = []
        o = self.offset
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!
        for i, point in enumerate(points):
            if not self.roi.contains(point):
                continue
            p = intPos((point + o) * self.scale)
            dc = fig.Circle(p, radius)
            for a in attr:
                setattr(dc, a, attr[a])
            if returnIndices:
                result.append((i, dc))
            else:
                result.append(dc)
            self.f.append(dc)
        return result

    def addPointOverlay(self, pointOverlay, **attr):
        points = pointOverlay.originalPointlist
        radius = pointOverlay.origRadius
        return self.addPointCircles(points, radius, **attr)

    def addMapNodes(self, map, radius, returnNodes = False, **attr):
        points = [node.position() for node in map.nodeIter()]
        result = self.addPointCircles(
            points, radius, returnIndices = returnNodes, **attr)
        if returnNodes:
            nodes = list(map.nodeIter())
            result = [(nodes[i], circle) for i, circle in result]
        return result

    def addMapEdges(self, map, **attr):
        result = []
        for edge in map.edgeIter():
            if hasattr(edge, "color"):
                parts = self.addClippedPoly(
                    edge, penColor = qtColor2figColor(edge.color), **attr)
            else:
                parts = self.addClippedPoly(edge, **attr)
            result.extend(parts)
        return result

    def addClippedPoly(self, polygon, **attr):
        if not self.roi:
            return [self.addEdge(polygon, **attr)]
        if not self.roi.intersects(polygon.boundingBox()):
            return []
        clipRect = BoundingBox(self.roi)
        clipRect.moveBy(-self.offset)
        if clipRect.contains(polygon.boundingBox()):
            return [self.addEdge(polygon, **attr)]
        return [self.addEdge(polygon, **attr)
                for polygon in clipPoly(polygon, clipRect)]

    def addEdge(self, points, **attr):
        o = self.offset
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!
        pp = simplifyPolygon(
            Polygon([(point + o) * self.scale for point in points]), 0.5)
        fp = fig.Polygon([intPos(v) for v in pp],
                         closed = pp[0] == pp[-1])
        for a in attr:
            setattr(fp, a, attr[a])
        self.f.append(fp)
        return fp
