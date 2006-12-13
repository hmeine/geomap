import qt, fig

from vigra import Vector2, readImage, Rect2D
from hourglass import BoundingBox, Polygon, simplifyPolygon, intPos
from dartpath import Path

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

import qt, fig
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

def qtColor2figColor(color, figFile):
    try:
        return _qtColor2figColorMapping[color]
    except KeyError:
        return figFile.getColor((color.red(), color.green(), color.blue()))

class FigExporter:
    """FigExporter objects represent an image range at a given scale,
    and allow adding objects: You can add polygons or points, which
    will be clipped, scaled, positioned, and converted to appropriate
    fig objects and added to an internal fig.File (which is accessible
    via figExporter.f).

    All addFoo() variants allow to set properties of the resulting fig
    objects via keyword arguments, e.g. addROIRect(myroi, depth = 40,
    penColor = fig.colorYellow), cf. documentation of fig.Object .

    For a discussion of the scaling/clipping, see documentation of
    __init__()."""
    
    def __init__(self, scale = 30, roi = None, offset = Vector2(0.5, 0.5)):
        """fe = FigExporter(90.0, BoundingBox(Vector2(10, 10), Vector2(50, 50)))

        Initializes a FigExporter with the given scale and roi.

        A scale of 450 makes one pixel (source unit) 450 fig units
        large, which equals 1cm.  Default is scale = 30.

        A roi is given as BoundingBox object (default None == no
        clipping) and represents a range in the original image (before
        scaling).

        An optional parameter offset (default: Vector2(0.5, 0.5)) is
        used to move all edges / points within the roi; this is
        intended to make e.g. Vector2(14.1, 10.0) be a little to the
        right of the *center* of the pixel (14, 10).

        The transformation parameters are stored as properties scale,
        roi, and offset, and these settings can be queried or even
        changed at any time (between adding objects)."""
        
        self.f = fig.File()
        self.scale = scale
        self.roi = roi
        if type(roi) == Rect2D:
            self.roi = BoundingBox(Vector2(*roi.upperLeft()),
                                   Vector2(*roi.lowerRight()))
        self.offset = offset

    def addROIRect(self, roi = None, **attr):
        """fe.addROIRect(roi, depth = 85, ...)

        Adds a rectangle around the given roi (ignoring fe.offset).
        If roi == None (default), the roi of the FigExporter itself is used.
        The fig.PolyBox object is returned."""
        
        assert roi or self.roi, "addROIRect(): no ROI given!?"
        if roi == None:
            roi = self.roi
        if self.roi:
            roi = BoundingBox(roi) # don't modify in-place
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
        """fe.addBackgroundWithFrame(bgImageFilename, depth = 85, ...)

        Adds a picture object to the fig.File, framed by an additional
        rectangle.  See addROIRect().  If no roi is given (via a
        keyword parameter), the image file is opened (using readImage)
        and its size is used to initialize a BoundingBox positioned at
        the origin.

        Returns the pair (bgImage, bgRect) of both added fig objects."""
        
        if not params.has_key("roi"):
            size = readImage(bgImageFilename).size()
            params["roi"] = BoundingBox(
                Vector2(0, 0), Vector2(size[0], size[1]))

        if not params.has_key("depth"):
            params["depth"] = 1
        
        bgRect = self.addROIRect(**params)

        bgImage = fig.PictureBBox(0, 0, 1, 1, bgImageFilename)
        bgImage.points = list(bgRect.points)
        bgImage.depth = 999
        self.f.append(bgImage)

        return bgImage, bgRect

    def addEdge(self, points, simplifyEpsilon = 0.5, **attr):
        """fe.addEdge(points, simplifyEpsilon, ...)

        Adds and returns exactly one fig.Polygon object representing
        the given points.  You will probably want to use
        addClippedPoly() instead.

        If simplifyEpsilon (default: 0.5) is not None, simplifyPolygon
        is called on the *scaled* polygon (i.e. the default is to
        simplify the polygon to 0.5 fig units, which are integer
        anyways)."""
        
        o = self.offset + attr.get('offset', Vector2(0,0))
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
        o = self.offset + attr.get('offset', Vector2(0,0))
        clipRect.moveBy(-o)
        if clipRect.contains(polygon.boundingBox()):
            return [self.addEdge(polygon, **attr)]
        return [self.addEdge(polygon, **attr)
                for polygon in clipPoly(polygon, clipRect)]

    def addPointCircles(self, points, radius, returnIndices = False, **attr):
        """fe.addPointCircles(points, radius, returnIndices = False,...)

        Marks each point in points with a circle of the given radius
        (in pixels) if it is within the clipping rect.  The list of
        fig.Circle objects is returned.  Again, note the possibility
        of setting properties (depth, penColor, lineStyle, lineWidth,
        ...) on all resulting objects via keyword arguments
        (cf. documentation of the FigExporter class).
        
        By default, circles will be filled, but have lineWidth=0. 
        To draw a transparent circle, call:
        
        fi.addPointCircles([Vector2(5.2,5.3)], 2, penColor=fig.colorCyan,fillStyle=fig.fillStyleNone,lineWidth=1)

        If returnIndices is set to True (default: False), a list of
        (i, c) pairs is returned instead, where c is the fig.Circle
        object, and i is the index of the corresponding position in
        points."""
        
        radius *= self.scale
        attr["fillStyle"] = attr.get("fillStyle", fig.fillStyleSolid)
        attr["lineWidth"] = attr.get("lineWidth", 0)
        result = []
        o = self.offset + attr.get('offset', Vector2(0,0))
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!
        for i, point in enumerate(points):
            if self.roi and not self.roi.contains(point+o):
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
        radius = float(pointOverlay.origRadius)
        if not pointOverlay.relativeRadius:
            radius /= pointOverlay.zoom
        if not attr.has_key("fillColor"):
            fillColor = pointOverlay.color
            if type(fillColor) == qt.QColor:
                fillColor = qtColor2figColor(fillColor, self.f)
            attr["fillColor"] = fillColor
        attr["lineWidth"] = attr.get("lineWidth", 0)
        
        return self.addPointCircles(points, radius, **attr)

    def addEdgeOverlay(self, edgeOverlay, **attr):
        """Adds and returns fig.Polygon for all edges (or -parts, see
        addClippedPoly) of the given overlay, using the overlays'
        color."""

        edges = edgeOverlay.originalEdges
        if not attr.has_key("penColor"):
            penColor = edgeOverlay.color
            if type(penColor) == qt.QColor:
                penColor = qtColor2figColor(penColor, self.f)
            attr["penColor"] = penColor
            
        result = []
        for edge in edges:
            if type(edge) == list:
                edge = Polygon(edge)
            parts = self.addClippedPoly(edge, **attr)
            result.extend(parts)
        return result

    def addMapNodes(self, map, radius, returnNodes = False, **attr):
        """fe.addMapNodes(map, radius, ...)

        See addPointCircles(), this function simply takes the
        positions of all nodes in the given map and marks them with a
        circle of the given radius (in pixels).

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

    def addMapEdges(self, map, returnEdges = False, **attr):
        """fe.addMapEdges(map, ...)

        Adds and returns fig.Polygons for all map edges (or -parts,
        see addClippedPoly).  If no penColor is given, only edges with
        a valid 'color' attribute are exported (can be either a fig or
        a Qt color).
        
        For example, to draw only a subregion, and shift the upper left of 
        the region to the origin, call
        
        fi.addMapEdges(map, penColor=fig.colorGreen, \
             offset=Vector2(-13,-13), roi=BoundingBox(Vector2(13,13), Vector2(24,24)))
        
        """

        result = []
        for edge in map.edgeIter():
            if attr.has_key("penColor"):
                parts = self.addClippedPoly(edge, **attr)
            elif hasattr(edge, "color") and edge.color:
                penColor = edge.color
                if type(penColor) == qt.QColor:
                    penColor = qtColor2figColor(penColor, self.f)
                thisattr = dict(attr)
                thisattr["penColor"] = penColor
                parts = self.addClippedPoly(edge, **thisattr)
            else:
                continue # skip invisible edge
            if returnEdges:
                result.extend([(edge, part) for part in parts])
            else:
                result.extend(parts)
        return result

    # FIXME: unfinished!:
    def addMapFaces(self, map, **attr):
        """fe.addMapEdges(map, ...)

        Adds and returns fig.Polygons for all map faces (or -parts,
        see addClippedPoly).  Clipping is experimental, since clipping
        arbitrary closed polygons is much harder than just open ones."""

        assert attr.has_key("fillColor"), \
               "TODO: addMapFaces() does not yet automatically find face colors!"
        
        facePolys = []
        for face in map.faceIter(skipInfinite = True):
            facePolys.append((
                Polygon(list(Path(face.contours()[0].phiOrbit()).points())), face))

        def AreaCompare(c1, c2):
            return -cmp(c1[0].partialArea(), c2[0].partialArea())
        facePolys.sort(AreaCompare)

        result = []
        for poly, face in facePolys:
            if hasattr(face, "color"):
                #fillColor = ...
                thisattr = dict(attr)
                #thisattr["fillColor"] = fillColor
                parts = self.addClippedPoly(poly, **thisattr)
            else:
                parts = self.addClippedPoly(poly, **attr)
            result.extend(parts)
        return result

    def save(self, filename):
        """Save the resulting XFig file to 'filename' (cf. fig.File.save)."""

        return self.f.save(filename)

    def saveEPS(self, basename):
        """Save the resulting XFig file to [basename].{fig,eps} (cf. fig.File.save)."""

        return self.f.saveEPS(basename)

# --------------------------------------------------------------------
#                               USAGE
# --------------------------------------------------------------------

if False:
    # default scale, no ROI:
    fe = figexport.FigExporter()
    # use given scale and ROI:
    fe = figexport.FigExporter(scale = 15, roi = BoundingBox(Rect2D(dm.imageSize())))
    # add background image:
    fe.addBackgroundWithFrame("background.png")
    fe.addMapEdges(someMap) # give optional properties like lineWidth = 4
    someRadius = 0.1 # pixels
    fe.addMapNodes(someMap, someRadius)
    fe.saveEPS("example_basename")
