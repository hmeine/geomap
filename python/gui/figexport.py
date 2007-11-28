import os, sys, qt, fig, math

from vigra import Vector2, readImage, Rect2D, Point2D, Size2D, meshIter
import vigrapyqt
from hourglass import BoundingBox, Polygon, simplifyPolygon, intPos, contourPoly
from dartpath import Path
from polytools import clipPoly
import flag_constants

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

def figRect(rect2d):
    return fig.Rect(rect2d.left(), rect2d.top(),
                    rect2d.right(), rect2d.bottom())

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
        elif type(roi) == Size2D:
            self.roi = BoundingBox(Vector2(0, 0), Vector2(*roi))
        elif type(roi) == tuple:
            self.roi = BoundingBox(*roi)
        self.offset = Vector2(*offset)

    def position2Fig(self, pos):
        """fe.position2Fig(Vector2 pos) -> Vector2

        Maps a pixel position to fig coordinates in the same way in
        which e.g. polygon points are mapped by addEdge()."""
        result = (self.offset + pos)
        if self.roi:
            result -= self.roi.begin()
        result *= self.scale
        return result

    def figPt(self, width):
        """Convert a distance into XFig points (1/80 inches), which
        can be used e.g. for line widths."""
        return width*self.scale*80/self.f.ppi

    def addROIRect(self, roi = None, container = True, **attr):
        """fe.addROIRect(roi, depth = 85, ...)

        Adds a rectangle around the given roi (ignoring fe.offset).
        If roi == None (default), the roi of the FigExporter itself is used.
        The fig.PolyBox object is returned."""
        
        assert roi or self.roi, "addROIRect(): no ROI given!?"
        if roi == None:
            roi = BoundingBox(self.roi)
            roi.moveBy((-0.5, -0.5))
        if container == True:
            container = self.f

        # convert to BoundingBox if necessary and don't modify in-place
        if isinstance(roi, Rect2D):
            roi = BoundingBox(roi)
        else:
            roi = BoundingBox(roi)
            roi.moveBy(self.offset)

        if self.roi and not self.roi.contains(roi):
            sys.stderr.write("WARNING: addROIRect: ROI out of bounds!\n")
            roi.moveBy(-self.offset)
            poly = Polygon([roi.begin(), roi.begin() + (0, roi.size()[1]),
                            roi.end(), roi.begin() + (roi.size()[0], 0),
                            roi.begin()])
            result = self.addClippedPoly(poly,
                                         container = container,
                                         **attr)
            if result:
                assert len(result) == 1
                return result[0]
            return

        if self.roi:
            roi.moveBy(-self.roi.begin())

        result = fig.PolyBox(roi.begin()[0] * self.scale,
                             roi.begin()[1] * self.scale,
                             roi.end()[0] * self.scale,
                             roi.end()[1] * self.scale)
        container.append(result)
        for a in attr:
            setattr(result, a, attr[a])
        return result

    def addBackgroundWithFrame(self, bgImageFilename, container = True, **params):
        """fe.addBackgroundWithFrame(bgImageFilename, depth = 85, ...)

        Adds a picture object to the fig.File, framed by an additional
        rectangle.  See addROIRect() and addImage().  If no roi is
        given (via a keyword parameter), the image file is opened
        (using readImage) and its size is used to initialize a
        BoundingBox positioned at the origin.

        Returns the pair (bgImage, bgRect) of both added fig objects."""
        
        if not params.has_key("roi") and not self.roi:
            size = readImage(bgImageFilename).size()
            params["roi"] = Rect2D(size)
        if container == True:
            container = self.f

        if not params.has_key("depth"):
            params["depth"] = 1
        
        bgRect = self.addROIRect(**params)

        bgImage = fig.PictureBBox(0, 0, 1, 1, bgImageFilename)
        bgImage.points = list(bgRect.points)
        bgImage.depth = 999
        container.append(bgImage)

        return bgImage, bgRect

    def addImage(self, filename, container = True, **params):
        """fe.addImage(filename, depth = 85, ...)

        Add a picture object to the fig.File.  See addROIRect() and
        addBackgroundWithFrame().  If no roi is given (via a keyword
        parameter), the image file is opened (using readImage) and its
        size is used to initialize a BoundingBox positioned at the
        origin.  Returns the fig.PictureBBox object."""
        
        bgImage, bgRect = self.addBackgroundWithFrame(
            filename, container, **params)
        if container == True:
            container = self.f
        container.remove(bgRect)

        return bgImage

    def addPixelRaster(self, rect = None, labels = None, fill = None,
                       container = True, labelType = int, **attr):
        """The default label size looks good with 1cm^2 pixels
        (scale = 450)."""
        
        assert rect or labels or fill
        if rect is None:
            rect = Rect2D((labels or fill).size())
        elif not hasattr(rect, "upperLeft"):
            rect = Rect2D(Size2D(rect[0], rect[1]))
        
        if container == True:
            container = self.f
        result = fig.Compound(container)
        
        for x, y in meshIter(rect):
            boxX1 = (x - rect.left()) * self.scale
            boxY1 = (y - rect.top()) * self.scale
            boxX2 = boxX1 + self.scale
            boxY2 = boxY1 + self.scale

            pixelRect = fig.PolyBox(boxX1, boxY1, boxX2, boxY2)
            for a in attr:
                setattr(pixelRect, a, attr[a])

            if labels:
                textX = (boxX1 + boxX2) / 2
                textY = (boxY1 + boxY2) / 2

                label = fig.Text(Vector2(textX, textY),
                                 str(labelType(labels[x, y])),
                                 alignment = fig.alignCentered)
                label.pos += (0, label.height / 2)

                label.depth = pixelRect.depth - 10
                label.font = fig.fontHelvetica
                result.append(label)

            if fill:
                pixelRect.fillStyle = fig.fillStyleSolid
                pixelRect.fillColor = self.f.getColor(int(fill[x, y]))
                    
            result.append(pixelRect)
        return result

    def addText(self, position, text, container = True, **attr):
        if container == True:
            container = self.f

        obj = fig.Text(self.position2Fig(position), text)
        for a in attr:
            setattr(obj, a, attr[a])
        container.append(obj)
        return obj

    def addEdge(self, points, simplifyEpsilon = 0.5, container = True,
                figClass = fig.Polygon, **attr):
        """fe.addEdge(points, simplifyEpsilon, ...)

        Adds and returns exactly one fig.Polygon object representing
        the given points.  You will probably want to use
        addClippedPoly() instead.

        If simplifyEpsilon (default: 0.5) is not None, simplifyPolygon
        is called on the *scaled* polygon (i.e. the default is to
        simplify the polygon to 0.5 fig units, which are integer
        anyways)."""
        
        if container == True:
            container = self.f
        o = self.offset + attr.get('offset', Vector2(0,0))
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!
        pp = Polygon([(o + point) * self.scale for point in points])
        if simplifyEpsilon:
            pp = simplifyPolygon(pp, simplifyEpsilon)
        fp = figClass([intPos(v) for v in pp], closed = pp[0] == pp[-1])
        for a in attr:
            setattr(fp, a, attr[a])
        container.append(fp)
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

        if "fillColor" in attr and not "fillStyle" in attr:
            attr["fillStyle"] = fig.fillStyleSolid

        # no ROI to clip to?
        if not self.roi:
            return [self.addEdge(polygon, **attr)]

        if type(polygon) != Polygon:
            if not isinstance(polygon, list):
                polygon = Polygon(list(polygon))
            else:
                polygon = Polygon(polygon)

        clipRect = BoundingBox(self.roi)
        o = self.offset + attr.get('offset', Vector2(0,0))
        clipRect.moveBy(-o)

        # handle all-or-none cases:
        if not clipRect.intersects(polygon.boundingBox()):
            return []
        if clipRect.contains(polygon.boundingBox()):
            return [self.addEdge(polygon, **attr)]

        # general case: perform clipping, add parts one-by-one:
        result = [] # fig.Compound(container) - I dont't dare grouping here..
        closeAtBorder = (
            attr.get("fillStyle", fig.fillStyleNone) != fig.fillStyleNone)
        for part in clipPoly(polygon, clipRect, closeAtBorder):
            if part.length(): # don't add zero-length polygons
                result.append(self.addEdge(part, **attr))
        return result

    def addPointCircles(self, points, radius, returnIndices = False, container = True, **attr):
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
        
        if container == True:
            container = self.f

        attr = dict(attr)
        if "fillStyle" not in attr:
            attr["fillStyle"] = fig.fillStyleSolid
        if "lineWidth" not in attr and attr["fillStyle"] != fig.fillStyleNone:
            attr["lineWidth"] = 0

        compound = fig.Compound(container)
        if returnIndices:
            result = []
        else:
            result = compound

        o = self.offset + attr.get('offset', Vector2(0,0))
        o2 = self.roi and o - self.roi.begin() or o
        for i, point in enumerate(points):
            if self.roi:
                radiusOff = Vector2(radius, radius)
                circleBounds = BoundingBox(
                    o + point - radiusOff, o + point + radiusOff)
                if not self.roi.contains(circleBounds):
                    continue # FIXME: add arcs if intersects
            p = intPos((Vector2(point[0], point[1]) + o2) * self.scale)
            dc = fig.Circle(p, radius*self.scale)
            for a in attr:
                setattr(dc, a, attr[a])
            if returnIndices:
                result.append((i, dc))
            compound.append(dc)

        return result

    def addEdgels(self, edgels, length = 0.65, container = True, **attr):
        if container == True:
            container = self.f

        result = fig.Compound(container)

        for edgel in edgels:
            pos = Vector2(edgel.x, edgel.y)
            c =  math.cos(edgel.orientation)*length/2
            s = -math.sin(edgel.orientation)*length/2
            self.addClippedPoly([pos+(c,s), pos-(c,s)],
                                container = result, **attr)

        return result

    def _setOverlayColor(self, overlay, colorAttr, attr):
        """Set color (and possibly line width) from overlay attributes."""
        if colorAttr not in attr:
            color = overlay.color
            if type(color) == qt.QColor:
                color = qtColor2figColor(color, self.f)
            attr[colorAttr] = color
            #print "fetched %s %s from %s" % (colorAttr, color, overlay)

        if colorAttr == "penColor" and \
               hasattr(overlay, "fillColor") and "fillColor" not in attr:
            color = overlay.fillColor
            if type(color) == qt.QColor:
                color = qtColor2figColor(color, self.f)
            if color is not None:
                attr["fillColor"] = color

        if "lineWidth" not in attr:
            if colorAttr == "penColor" and overlay.color is None:
                attr["lineWidth"] = 0
            elif hasattr(overlay, "width"):
                attr["lineWidth"] = overlay.width + 1

    def addPointOverlay(self, pointOverlay, container = True, **attr):
        """See addPointCircles(), this function simply takes the
        points and radius from a PointOverlay object for your
        convenience."""
        
        points = pointOverlay.originalPoints
        radius = float(pointOverlay.origRadius)
        if not pointOverlay.relativeRadius:
            radius /= pointOverlay.zoom

        attr = dict(attr)
        self._setOverlayColor(pointOverlay, "fillColor", attr)
        attr["lineWidth"] = attr.get("lineWidth", 0)
        
        return self.addPointCircles(points, radius, container = container, **attr)

    def addEdgeOverlay(self, edgeOverlay, container = True, **attr):
        """Adds and returns fig.Polygon for all edges (or -parts, see
        addClippedPoly) of the given overlay, using the overlays'
        color."""

        if container == True:
            container = self.f

        edges = edgeOverlay.originalEdges
        attr = dict(attr)
        self._setOverlayColor(edgeOverlay, "penColor", attr)

        result = fig.Compound(container)
        for edge in edges:
            parts = self.addClippedPoly(edge, container = result, **attr)
        return result

    def addCircleOverlay(self, circleOverlay, container = True, **attr):
        """Adds and returns fig.Circle for all circles of the given
        overlay, using the overlays' color and width."""

        if container == True:
            container = self.f

        circles = circleOverlay.originalCircles
        attr = dict(attr)
        self._setOverlayColor(circleOverlay, "penColor", attr)
            
        o = self.offset + attr.get('offset', Vector2(0,0))
        if self.roi:
            o = o - self.roi.begin() # don't modify in-place!

        result = fig.Compound(container)
        for center, radius in circles:
            if self.roi and not self.roi.contains(center+o):
                continue
            p = intPos((Vector2(center[0], center[1]) + o) * self.scale)
            dc = fig.Circle(p, radius * self.scale)
            for a in attr:
                setattr(dc, a, attr[a])
            result.append(dc)

        return result

    def addMapNodes(self, map, radius, returnNodes = False, container = True, **attr):
        """fe.addMapNodes(map, radius, ...)

        See addPointCircles(), this function simply takes the
        positions of all nodes in the given map and marks them with a
        circle of the given radius (in pixels).

        If the optional parameter returnNodes is set to True, a list
        of (node, circleObject) pairs is returned, similar to the
        returnIndices parameter of addPointCircles()."""
        
        points = [node.position() for node in map.nodeIter()]
        result = self.addPointCircles(
            points, radius, returnIndices = returnNodes, container = container, **attr)
        if returnNodes:
            nodes = list(map.nodeIter())
            result = [(nodes[i], circle) for i, circle in result]
        return result

    def addMapEdges(self, map, skipBorder = False,
                    returnEdges = False, container = True, **attr):
        """Adds and returns fig.Polygons for all map edges (or -parts,
        see addClippedPoly).  If no penColor is given, only edges with
        a valid 'color' attribute are exported (can be either a fig or
        a Qt color).
        
        For example, to draw only a subregion, and shift the upper left of
        the region to the origin, call
        
        fi.addMapEdges(map, penColor=fig.colorGreen, \
             offset=Vector2(-13,-13), roi=BoundingBox(Vector2(13,13), Vector2(24,24)))
        
        """

        if container == True:
            container = self.f

        compound = fig.Compound(container)
        if returnEdges:
            result = []
        else:
            result = compound

        for edge in map.edgeIter():
            if skipBorder and edge.flag(flag_constants.BORDER_PROTECTION):
                continue
            parts = self.addClippedPoly(edge, container = compound, **attr)
            if returnEdges:
                result.extend([(edge, part) for part in parts])

        return result

    def addMapFaces(self, geomap, faceMeans = None, similarity = None,
                    returnFaces = False, container = True, **attr):
        """fe.addMapFaces(geomap, faceMeans, ...)

        Adds and returns fig.Polygons for all map faces (or -parts,
        see addClippedPoly).  Clipping is experimental, since clipping
        arbitrary closed polygons is much harder than just open ones."""

        import maputils

        def getGray(face):
            faceColor = faceMeans[face.label()]
            return self.f.gray(int(faceColor))

        def getRGB(face):
            faceColor = faceMeans[face.label()]
            return self.f.getColor(map(int, tuple(faceColor)), similarity)

        if not faceMeans:
            if not hasattr(geomap, "faceMeans"):
                raise ValueError("addMapFaces: need faceMeans for proper coloring")
            faceMeans = geomap.faceMeans

        if hasattr(faceMeans, "bands"):
            getFaceColor = getGray
            if faceMeans.bands() == 3:
                getFaceColor = getRGB
        else:
            getFaceColor = lambda face: faceMeans[face.label()]

        if container == True:
            container = self.f

        attr = dict(attr)
        attr["lineWidth"] = attr.get("lineWidth", 0)
        attr["fillStyle"] = attr.get("fillStyle", fig.fillStyleSolid)

        compound = fig.Compound(container)
        if returnFaces:
            result = []
        else:
            result = compound

        todo = [geomap.face(0)]
        # +1 because the first iteration will not add any objects:
        currentDepth = attr.get("depth", 100) + 1
        while todo:
            thisLayer = todo
            todo = []
            for face in thisLayer:
                if face.area() > 0:
                    color = getFaceColor(face)
                    if color is not None:
                        thisattr = dict(attr)
                        thisattr["fillColor"] = getFaceColor(face)
                        thisattr["depth"] = currentDepth
                        parts = self.addClippedPoly(
                            contourPoly(face.contour()),
                            container = compound, **thisattr)

                        if returnFaces:
                            result.extend([(face, part) for part in parts])

                for anchor in face.holeContours():
                    todo.extend(maputils.holeComponent(anchor))

            currentDepth -= 1

        return result

    def save(self, filename, fig2dev = None):
        """Save the resulting XFig file to 'filename' (cf. fig.File.save)."""

        return self.f.save(filename, fig2dev)

    def saveEPS(self, basename):
        """Save the resulting XFig file to [basename].{fig,eps} (cf. fig.File.save)."""

        return self.save(basename, fig2dev = "eps")

    def savePDF(self, basename):
        """Save the resulting XFig file to [basename].{fig,pdf} (cf. fig.File.save)."""

        return self.save(basename, fig2dev = "pdf")

# --------------------------------------------------------------------

def addStandardOverlay(fe, overlay, **attr):
    # FIXME: str(type(overlay)).contains(...) instead?
    if isinstance(overlay, vigrapyqt.PointOverlay):
        return fe.addPointOverlay(overlay, **attr)
    elif isinstance(overlay, vigrapyqt.EdgeOverlay):
        return fe.addEdgeOverlay(overlay, **attr)
    elif isinstance(overlay, vigrapyqt.CircleOverlay):
        return fe.addCircleOverlay(overlay, **attr)

def _exportOverlays(fe, overlays, overlayHandler, startDepth = 100):
    depth = startDepth
    for overlay in overlays:
        if hasattr(overlay, 'visible') and not overlay.visible:
            continue
        if overlayHandler(fe, overlay, depth = depth) is None:
            if isinstance(overlay, vigrapyqt.OverlayGroup):
                depth = _exportOverlays(
                    fe, overlay.overlays, overlayHandler, startDepth = depth)
            else:
                sys.stderr.write(
                    "exportImageWindow: overlay type %s not handled!\n" % (
                    type(overlay)))
        depth -= 1
    return depth

def visibleROI(imageWindow):
    viewer = hasattr(imageWindow, "viewer") \
             and imageWindow.viewer or imageWindow
    return Rect2D(
        intPos(viewer.toImageCoordinates(0, 0)),
        intPos(viewer.toImageCoordinates(viewer.width(), viewer.height())))

def exportImageWindow(
    w, basepath, roi = None, scale = None,
    bgFilename = None,
    overlayHandler = addStandardOverlay):
    """Exports the currently displayed contents of the ImageWindow `w`
    as XFig file <basepath>.fig.

    `roi` may be a specific ROI to save, or `True` to save the visible
    image region.  The default (None) saves the whole image.

    If `bgFilename` is not given, the displayed background image is
    saved as <basepath>_bg.png and used as background.  The background
    can be disabled entirely by setting `bgFilename` to `False`.

    Overlays are saved using the given overlayHandler.  The default
    handler (addStandardOverlay) can handle the standard vigrapyqt
    point, edge, and circle overlays."""
    
    figFilename = basepath + ".fig"
    pngFilename = basepath + "_bg.png"

    # determine ROI to be saved
    if roi == None:
        roi = Rect2D(w.image.size())
    elif roi == True:
        roi = visibleROI(w)
    elif type(roi) == tuple:
        roi = Rect2D(*roi)
    elif type(roi) == str:
        roi = Rect2D(*fig.parseGeometry(roi))
    
    if bgFilename == None:
        # create .png background
        image, normalize = w.getDisplay()
        image.subImage(roi).write(pngFilename, normalize)
        _, bgFilename = os.path.split(pngFilename)

    # create .fig file
    if scale == None:
        scale = 20*450 / roi.width() # default: 20cm width
        print "auto-adjusted scale to %s." % (scale, )
    fe = FigExporter(scale, roi)
    if bgFilename != False:
        fe.addBackgroundWithFrame(bgFilename, depth = 100, lineWidth = 0)
    else:
        fe.addROIRect(depth = 100, lineWidth = 0)

    _exportOverlays(fe, w.viewer.overlays, overlayHandler)
    fe.save(figFilename)

    return fe

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
