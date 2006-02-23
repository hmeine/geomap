def _intersectLine(inside, outside, clipRect):
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
    """clipPoly(polygon, clipRect)

    Clips away those parts of polygon which are not in clipRect.
    Returns a list of polygons (since the polygon may leave clipRect,
    enter again, leave, ...).  Polygon segments crossing clipRect's
    borders are cut, such that the resulting polyons get new endpoints
    exactly on the border."""
    
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
            p.append(_intersectLine(polygon[i], polygon[i-1], clipRect))

        while i < len(polygon) and clipRect.contains(polygon[i]):
            p.append(polygon[i])
            i += 1

        if i < len(polygon):
            p.append(_intersectLine(polygon[i-1], polygon[i], clipRect))

        result.append(p)

    return result

_qtColor2figColorMapping = {
    qt.Qt.black : fig.colorBlack,
    qt.Qt.blue : fig.colorBlue,
    qt.Qt.green : fig.colorGreen,
    qt.Qt.cyan : fig.colorCyan,
    qt.Qt.red : fig.colorRed,
    qt.Qt.magenta : fig.colorMagenta,
    qt.Qt.yellow : fig.colorYellow,
    qt.Qt.white : fig.colorWhite,
    }

def qtColor2figColor(color):
    return _qtColor2figColorMapping[color]

class FigExporter:
    """FigExporter objects represent an image range at a given scale,
    and allow adding objects: You can add polygons or points, which
    will be clipped, scaled, positioned, and converted to appropriate
    fig objects and added to an internal fig.File (which is accessible
    via figExporter.f).

    All addFoo() variants allow to set properties of the resulting fig
    objects via keyword arguments, e.g. addROIRect(myroi, depth = 40,
    penColor = fig.colorYellow), cf. documentation of fig.Object ."""
    
    def __init__(self, scale = 30, roi = None):
        """fe = FigExporter(90.0, BoundingBox(10, 10, 50, 50))

        Initializes a FigExporter with the given scale and roi.

        A scale of 450 makes one pixel (source unit) 450 fig units
        large, which equals 1cm.  Default is scale = 30.

        A roi is given as BoundingBox object (default None == no
        clipping) and represents a range in the original image (before
        scaling)."""
        
        self.f = fig.File()
        self.scale = scale
        self.roi = roi
        self.offset = Vector2(0.5, 0.5)

    def addROIRect(self, roi = None, **attr):
        """fe.addROIRect(roi, depth = 85, ...)

        Adds a rectangle around the given roi.
        If roi == None (default), the roi of the FigExporter itself is used.
        The fig.PolyBox object is returned."""
        
        assert roi or self.roi, "addROIRect(): no ROI given!?"
        if roi == None:
            roi = BoundingBox(self.roi)
        if self.roi:
            roi.moveBy(-self.roi.begin())
        result = fig.PolyBox(roi.begin()[0] * self.scale,
                             roi.begin()[1] * self.scale,
                             roi.end()[0] * self.scale,
                             roi.end()[1] * self.scale)
        self.f.append(result)
        for a in attr:
            setattr(result, a, attr[a])
        return result

    def addBackgroundWithFrame(self, bgImageFilename, **params):
        """fe.addROIRect(bgImageFilename, depth = 85, ...)

        Adds a picture object to the fig.File, framed by an additional
        rectangle.  See addROIRect().  If no roi is given (via a
        keyword parameter), the image file is opened (using readImage)
        and its size is used to initialize a BoundingBox positioned at
        the origin."""
        
        if not params.has_key("roi"):
            size = readImage(bgImageFilename).size()
            params["roi"] = BoundingBox(
                Vector2(0, 0), Vector2(size[0], size[1]))

        bgRect = self.addROIRect(**params)

        bgImage = fig.PictureBBox(0, 0, 1, 1, bgImageFilename)
        bgImage.points = list(bgRect.points)
        bgImage.depth = bgRect.depth - 1
        self.f.append(bgImage)

    def addEdge(self, points, simplifyEpsilon = 0.5, **attr):
        """fe.addEdge(points, simplifyEpsilon, ...)

        Adds and returns exactly one fig.Polygon object representing
        the given points.  You will probably want to use
        addClippedPoly() instead.

        If simplifyEpsilon (default: 0.5) is not None, simplifyPolygon
        is called on the *scaled* polygon (i.e. the default is to
        simplify the polygon to 0.5 fig units, which are integer
        anyways)."""
        
        o = self.offset
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!
        pp = Polygon([(point + o) * self.scale for point in points])
        if simplifyEpsilon:
            pp = simplifyPolygon(pp, simplifyEpsilon)
        fp = fig.Polygon([intPos(v) for v in pp],
                         closed = pp[0] == pp[-1])
        for a in attr:
            setattr(fp, a, attr[a])
        self.f.append(fp)
        return fp

    def addClippedPoly(self, polygon, **attr):
        """fe.addClippedPoly(polygon, ...)

        Adds and returns exactly fig.Polygon objects for each part of
        the given polygon which is in the clipping range.  Again, note
        the possibility of setting properties (depth, penColor,
        lineStyle, lineWidth, ...) on all resulting objects via
        keyword arguments (cf. documentation of the FigExporter
        class).

        If simplifyEpsilon (default: 0.5) is not None, simplifyPolygon
        is called on the *scaled* polygon (i.e. the default is to
        simplify the polygon to 0.5 fig units, which are integer
        anyways)."""
        
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

    def addPointCircles(self, points, radius, returnIndices = False, **attr):
        """fe.addPointCircles(points, radius, returnIndices = False,...)

        Marks each point in points with a circle of the given radius
        if it is within the clipping rect.  The list of fig.Circle
        objects is returned.  Again, note the possibility of setting
        properties (depth, penColor, lineStyle, lineWidth, ...) on all
        resulting objects via keyword arguments (cf. documentation of
        the FigExporter class).

        If returnIndices is set to True (default: False), a list of
        (i, c) pairs is returned instead, where c is the fig.Circle
        object, and i is the index of the corresponding position in
        points."""
        
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
        """See addPointCircles(), this function simply takes the
        points and radius from a PointOverlay object for your
        convenience."""
        
        points = pointOverlay.originalPointlist
        radius = pointOverlay.origRadius
        return self.addPointCircles(points, radius, **attr)

    def addMapNodes(self, map, radius, returnNodes = False, **attr):
        """fe.addMapNodes(map, radius, ...)

        See addPointCircles(), this function simply takes the
        positions of all nodes in the given map and marks them with a
        circle of the given radius.

        If the optional parameter returnNodes is set to True, a list
        of (node, circleObject) pairs is returned, similar to the
        returnIndices parameter of addPointCircles()."""
        
        points = [node.position() for node in map.nodeIter()]
        result = self.addPointCircles(
            points, radius, returnIndices = returnNodes, **attr)
        if returnNodes:
            nodes = list(map.nodeIter())
            result = [(nodes[i], circle) for i, circle in result]
        return result

    def addMapEdges(self, map, **attr):
        """fe.addMapEdges(map, ...)

        Adds and returns fig.Polygon for all map edges (or -parts, see
        addClippedPoly)."""
        
        result = []
        for edge in map.edgeIter():
            if hasattr(edge, "color"):
                penColor = edge.color
                if type(penColor) == qt.QColor:
                    penColor = qtColor2figColor(penColor)
                thisattr = dict(attr)
                thisattr["penColor"] = penColor
                parts = self.addClippedPoly(edge, **thisattr)
            else:
                parts = self.addClippedPoly(edge, **attr)
            result.extend(parts)
        return result

    def save(self, filename):
        """Save the resulting XFig file to 'filename' (cf. fig.File.save)."""

        self.f.save(filename)
