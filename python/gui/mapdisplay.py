_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import fig, figexport, maputils, flag_constants, qt, sys, os, time, tools
from vigra import BYTE, NBYTE, Point2D, Rect2D, Vector2, GrayImage
from vigrapyqt import ImageWindow, EdgeOverlay, PointOverlay, CircleOverlay
from hourglass import simplifyPolygon, intPos, BoundingBox, contourPoly
from maputils import removeCruft, holeComponent
from weakref import ref

# ui-generated base classes:
from displaysettings import DisplaySettings
from dartnavigator import DartNavigatorBase

def findZoomFactor(srcSize, destSize):
    if destSize > srcSize:
        return 1.0/findZoomFactor(destSize, srcSize)
    result = 1.0
    potentialDiff = 0
    while destSize < srcSize:
        srcSize = srcSize * 0.5 # don't modify in-place with *=!
        result *= 0.5
        potentialDiff = potentialDiff * 2 + 1
    while destSize - potentialDiff > srcSize:
        srcSize = srcSize * 2
        result *= 2
    return result

class MapEdges(object):
    def __init__(self, map, color, width = 0):
        self.setMap(map)
        self.color = color
        self.colors = None
        self.width = width
        self._zoom = None
        self._attachedHooks = None

    def setMap(self, map):
        self._map = ref(map)
        self._zoomedEdges = [None] * map.maxEdgeLabel()

    def attachHooks(self):
        map = self._map()
        self._attachedHooks = (
            map.addRemoveBridgeCallbacks(self.preRemoveEdgeHook, None),
            map.addMergeFacesCallbacks(self.preRemoveEdgeHook, None),
            map.addMergeEdgesCallbacks(self.preMergeEdgesHook, self.postMergeEdgesHook))

    def detachHooks(self):
        """Detaches / removes callbacks from the map's hooks.
        Returns True if successful, False if already detached."""
        map = self._map()
        if not map or not self._attachedHooks:
            return False

        for cb in self._attachedHooks:
            cb.disconnect()
        self._attachedHooks = None
        return True

    def changeColor(self, edge, newColor):
        if isinstance(edge, int):
            edge = self._map().edge(edge)
        if not self.colors:
            self.colors = [self.color] * self._map().maxEdgeLabel()
        self.colors[edge.label()] = newColor
        self._updateEdgeROI(edge)

    def _updateEdgeROI(self, edge):
        """FIXME: only used by changeColor so far, similar to setEdgePoints"""
        result = self._getZoomedEdge(edge).boundingRect()
        result.moveBy(self.viewer.x, self.viewer.y)
        self.viewer.update(result)

    def _calculatePoints(self):
        if not self._map():
            return
        c = time.clock()
        if len(self._zoomedEdges) != self._map().maxEdgeLabel():
            self._zoomedEdges = [None] * self._map().maxEdgeLabel()
        if hasattr(self._map(), "edges"):
            # "secret", internal API available? (not for C++ variant ATM)
            edges = self._map().edges
            for label, edge in enumerate(edges):
                if edge:
                    self._calculateEdgePoints(edge.label(), edge)
                else:
                    self._zoomedEdges[label] = None
        else:
            expected = 0
            for edge in self._map().edgeIter():
                label = edge.label()
                while expected < label:
                    self._zoomedEdges[expected] = None
                    expected += 1
                self._calculateEdgePoints(label, edge)
                expected = label + 1
#         sys.stdout.write("MapEdges._calculatePoints(zoom = %s) took %ss.\n" % (
#             self._zoom, time.clock() - c, ))

    def _calculateEdgePoints(self, index, origEdgePoints):
        offset = Vector2(self._zoom / 2.0 - 0.5, self._zoom / 2.0 - 0.5)
        origEdgePoints = (
            simplifyPolygon(origEdgePoints * self._zoom, 0.5)
            + offset).roundToInteger()

        qpa = self._zoomedEdges[index]
        if qpa == None:
            qpa = qt.QPointArray(len(origEdgePoints))
        elif qpa.size() < len(origEdgePoints):
            qpa.resize(len(origEdgePoints))

        points = iter(origEdgePoints)
        oldPos = points.next()
        qpa.setPoint(0, oldPos[0], oldPos[1])

        qpaIndex = 1
        for pos in points:
            if pos != oldPos:
                qpa.setPoint(qpaIndex, pos[0], pos[1])
                qpaIndex += 1
                oldPos = pos

        # extra parameter in c++ would be QGArray.SpeedOptim (1), but
        # is not available here:
        if qpaIndex != qpa.size():
            qpa.resize(qpaIndex)

        self._zoomedEdges[index] = qpa
        return qpa

    def setEdgePoints(self, index, edgePoints):
        if not self._zoomedEdges or index >= len(self._zoomedEdges):
            # FIXME: this happens if new edges were created (e.g. with
            # splitEdge) and we do not know them yet
            return

        if not self.visible:
            # mark zoomed edge as dirty:
            self._zoomedEdges[index] = None
            return

        if self._zoomedEdges[index]:
            updateROI = self._zoomedEdges[index].boundingRect()
        else:
            updateROI = qt.QRect()

        if edgePoints:
            qpa = self._calculateEdgePoints(index, edgePoints)
            updateROI |= qpa.boundingRect()
        else:
            self._zoomedEdges[index] = None

        updateROI.moveBy(self.viewer.x, self.viewer.y)
        self.viewer.update(updateROI)

    def preRemoveEdgeHook(self, dart):
        self.setEdgePoints(dart.edgeLabel(), None)
        return True

    def preMergeEdgesHook(self, dart):
        self.setEdgePoints(
            dart.clone().nextSigma().edgeLabel(), None)
        return True

    def postMergeEdgesHook(self, edge):
        self.setEdgePoints(edge.label(), edge)

    def setZoom(self, zoom):
        if self._map():
            zoom *= findZoomFactor(
                self._map().imageSize()[0],
                self.viewer.pixmap.size().width())
        if self._zoom != zoom:
            self._zoom = zoom
            self._dirty = True
            self._zoomedEdges = [None] * len(self._zoomedEdges)

    def _getZoomedEdge(self, edge):
        """Return the QPointArray for the given GeoMap.Edge object.
        The QPointArrays are cached and will be created on demand if
        you call this function."""
        index = edge.label()
        result = self._zoomedEdges[index]
        if result == None:
            result = self._calculateEdgePoints(index, edge)
        return result
    
    def draw(self, p):
        if not self._map():
            return
        r = p.clipRegion().boundingRect()
        r.moveBy(-self.viewer.x, -self.viewer.y)
        bbox = BoundingBox(Vector2(r.left() / self._zoom - 0.5,
                                   r.top() / self._zoom - 0.5),
                           Vector2(r.right() / self._zoom + 0.5,
                                   r.bottom() / self._zoom + 0.5))
        map = self._map()
        if self.colors:
            try:
                for e in map.edgeIter():
                    if bbox.intersects(e.boundingBox()):
                        edgeColor = self.colors[e.label()]
                        if edgeColor:
                            p.setPen(qt.QPen(edgeColor, self.width))
                            p.drawPolyline(self._getZoomedEdge(e))
            except IndexError, e:
                print e #"IndexError: %d > %d (maxEdgeLabel: %d)!" % (
                    #i, len(self.colors), map.maxEdgeLabel())
        else:
            p.setPen(qt.QPen(self.color, self.width))
            for e in map.edgeIter():
                if bbox.intersects(e.boundingBox()):
                    p.drawPolyline(self._getZoomedEdge(e))

class MapNodes(object):
    def __init__(self, map, color,
                 radius = 1, relativeRadius = False):
        self.setMap(map)
        self.color = color
        self.setRadius(radius, relativeRadius)
        self._zoom = None
        self._attachedHook = None

    def setMap(self, map):
        self._map = ref(map)
        self._qpointlist = None

    def attachHooks(self):
        map = self._map()
        self._attachedHook = map.addRemoveNodeCallback(self.removeNode)

    def detachHooks(self):
        """Detaches / removes callbacks from the map's hooks.
        Returns True if successful, False if already detached."""
        map = self._map()
        if not map or not self._attachedHook:
            return False

        self._attachedHook.disconnect()
        self._attachedHook = None
        return True

    def _calculatePoints(self):
        if not self._map():
            self._qpointlist = []
            return
        c = time.clock()
        if self.relativeRadius:
            self.radius = int(self._zoom * self.origRadius + 0.5)
        d0 = Vector2(0.5 * (self._zoom-1) - self.radius,
                     0.5 * (self._zoom-1) - self.radius)
        w = 2 * self.radius + 1
        self.s = qt.QSize(w, w)
        self._qpointlist = [None] * self._map().maxNodeLabel()
        for node in self._map().nodeIter():
            ip = intPos(node.position() * self._zoom + d0)
            self._qpointlist[node.label()] = qt.QPoint(ip[0], ip[1])
#         sys.stdout.write("MapNodes._calculatePoints(zoom = %s) took %ss.\n" % (
#             self._zoom, time.clock() - c))

    def removeNode(self, node):
        if self._qpointlist:
            nodeLabel = node.label()
            if not self._qpointlist[nodeLabel]:
                sys.stderr.write("WARNING: MapNodes.removeNode(): Node already None!\n")
                return
            if self.visible:
                ur = qt.QRect(self._qpointlist[nodeLabel], self.s)
                ur.moveBy(self.viewer.x, self.viewer.y)
                self.viewer.update(ur)
            self._qpointlist[nodeLabel] = None
        return True

    def setZoom(self, zoom):
        if self._map():
            zoom *= findZoomFactor(
                self._map().imageSize()[0],
                self.viewer.pixmap.size().width())
        if self._zoom != zoom:
            self._zoom = zoom
            self._qpointlist = None

    def setRadius(self, radius, relativeRadius = False):
        """nodeOverlay.setRadius(radius, relativeRadius = False)

        radius is given in
        * display pixels if relativeRadius == False
        * image pixels   if relativeRadius == True
        """
        
        self.origRadius = radius
        self.radius = radius
        self.relativeRadius = relativeRadius
        self._qpointlist = None

    def draw(self, p):
        if not self._qpointlist:
            self._calculatePoints()
        p.setPen(self.color)
        if self.radius == 0:
            for point in self._qpointlist:
                if point: # TODO: boundingRect
                    #p.drawPoint(point)
                    p.drawLine(point, point)
        else:
            p.setBrush(qt.QBrush(self.color))
            for point in self._qpointlist:
                if point: # TODO: boundingRect
                    p.drawEllipse(qt.QRect(point, self.s))

# --------------------------------------------------------------------
#                          DartHighlighter
# --------------------------------------------------------------------

class DartHighlighter(object):
    """The DartHighlighter class is attached to a viewer and a Map and
    is able to highlight any set of darts at a time."""
    
    def __init__(self, map, viewer):
        self._map = ref(map)
        self._viewer = viewer
        self.eo = None
        self.no = None

    def setMap(self, map):
        self._map = ref(map)

    def highlight(self, darts, color = qt.Qt.yellow):
        """highlight(darts)
        Highlight the given darts (can be any iterable returning labels
        or Dart objects)."""

        if darts:
            # hasattr(darts, "__iter__"): would be better, but is True for Edges
            if not isinstance(darts, (list, tuple)):
                darts = [darts]
            dartObjects = []
            for dart in darts:
                if type(dart) == int:
                    dart = self._map().dart(dart)
                elif hasattr(dart, "anchor"): # Nodes
                    dart = dart.anchor()
                elif hasattr(dart, "dart"): # Edges
                    dart = dart.dart()
                dartObjects.append(dart)
            darts = dartObjects

        if self.eo != None:
            self._viewer.removeOverlay(self)
            self.eo = None
            self.no = None

        if darts == None or not len(darts):
            return

        self.eo = EdgeOverlay([dart.edge() for dart in darts], color)
        self.eo.width = 2
        self.no = PointOverlay(
            [dart.startNode().position() for dart in darts], color, 3)
        self.color = color # used in the viewer's RMB menu
        self._viewer.addOverlay(self)

    def setZoom(self, zoom):
        self.eo.viewer = self.viewer
        self.eo.setZoom(zoom)
        self.no.viewer = self.viewer
        self.no.setZoom(zoom)

    def draw(self, p):
        self.eo.draw(p)
        self.no.draw(p)

# --------------------------------------------------------------------
#                             MapDisplay
# --------------------------------------------------------------------

class MapDisplay(DisplaySettings):
    def __init__(self, map, preparedImage = None, immediateShow = True,
                 faceMeans = None):
        DisplaySettings.__init__(self)
        # for backward compatibility:
        if hasattr(preparedImage, "imageSize") and hasattr(map, "width"):
            map, preparedImage = preparedImage, map
        elif preparedImage == None:
            preparedImage = map.labelImage()
            if not preparedImage:
                preparedImage = GrayImage(map.imageSize())
        
        self.preparedImage = preparedImage
        self.map = map
        self.setFaceMeans(faceMeans)
        self._togglingGUI = False
        self._attachedHooks = None

        if not hasattr(preparedImage, "orig"):
            self.backgroundGroup.setEnabled(False)
            self.image = preparedImage
        else:
            self.image = preparedImage.view
            self.displayColoredAction.setEnabled(
                self.preparedImage.orig.bands() == 3)

        self.imageWindow = ImageWindow(self.image, BYTE, self)
        self.imageWindow.label.hide()
        self.setCentralWidget(self.imageWindow)
        self.viewer = self.imageWindow.viewer
        self.viewer.autoZoom(4.0)

        self.tool = None

        self.normalizeStates = [False, False, True, True, False]
        self._backgroundMode = 0
        self.displayOriginalAction.setOn(True)

        self.connect(self.backgroundGroup, qt.SIGNAL("selected(QAction*)"),
                     self.setBackgroundMode)
        self.connect(self.normalizeAction, qt.SIGNAL("toggled(bool)"),
                     self.toggleNormalize)
        self.connect(self.mapCleanupAction, qt.SIGNAL("activated()"),
                     self.cleanupMap)
        self.connect(self.paintbrushAction, qt.SIGNAL("toggled(bool)"),
                     self.activatePaintbrush)
        self.connect(self.navigateAction, qt.SIGNAL("toggled(bool)"),
                     self.activateNavigator)
        self.connect(self.imageWindow, qt.PYSIGNAL("captionChanged"),
                     self.statusMessage)
                     #self.statusBar(), qt.SLOT("message(const QString&)"))

        self.setCaption("Map Display (map.py Version %s)" % (
            sys.modules["map"]._cvsVersion[0][1], ))

        self.edgeOverlay = MapEdges(map, qt.Qt.red)
        self.nodeOverlay = MapNodes(map, qt.Qt.blue, 1)
        self.viewer.addOverlay(self.edgeOverlay)
        self.viewer.addOverlay(self.nodeOverlay)
        self._dh = DartHighlighter(map, self.viewer)

        if immediateShow:
            self.show()

    def __del__(self):
        # delete tool (which may reference the viewer & map)
        self.setTool(None)

    def setMap(self, map):
        attached = self.detachHooks()

        self.setTool(None)
        self.map = map
        self.edgeOverlay.setMap(map)
        self.nodeOverlay.setMap(map)
        if self._backgroundMode == 3:
            self.setImage(self.map.labelImage(), self.normalizeStates[self._backgroundMode] and NBYTE or BYTE)
        self._dh.setMap(map)

        if attached:
            self.attachHooks()

        self.viewer.update()

    def setFaceMeans(self, faceMeans):
        self.faceMeans = faceMeans
        if faceMeans:
            tools.activePathMeasure = faceMeans.faceMeanDiff
        self.displayMeansAction.setEnabled(bool(faceMeans))

    def _adjustSize(self):
        pass # don't change window size out of a sudden

    def statusMessage(self, msg):
        # workaround - don't know why connecting the PYSIGNAL directly
        # does not work (maybe qt.PYSIGNAL("captionChanged(QString)")
        # would work?)
        self.statusBar().message(msg)

    def setBackgroundMode(self, mode):
        if type(mode) != int:
            mode = [self.displayOriginalAction,
                    self.displayColoredAction,
                    self.displayBIAction,
                    self.displayLabelsAction,
                    self.displayMeansAction].index(mode)

        displayImage = None
        if mode == 0:
            displayImage = self.preparedImage.view
        elif mode == 1:
            displayImage = self.preparedImage.orig
        elif mode == 2:
            displayImage = self.preparedImage.bi.gm
        elif mode == 3:
            displayImage = self.map.labelImage()
        elif mode == 4:
            displayImage = self.faceMeans.regionImage(self.map.labelImage())
        else:
            sys.stderr.write("Unknown background mode %d!\n" % mode)
            return

        self.image = displayImage
        self._backgroundMode = mode
        normalize = self.normalizeStates[mode]
        if self.normalizeAction.isOn() != normalize:
            self.normalizeAction.setOn(normalize)
        else:
            self.setImage(self.image, normalize and NBYTE or BYTE)

    def toggleNormalize(self, normalize):
        self.normalizeStates[self._backgroundMode] = normalize
        self.setImage(self.image, normalize and NBYTE or BYTE)

    def attachHooks(self):
        self.nodeOverlay.attachHooks()
        self.edgeOverlay.attachHooks()
        self._attachedHooks = (
            self.map.addMergeFacesCallbacks(None, self._postMergeFacesHook),
            self.map.addRemoveBridgeCallbacks(self._preRemoveBridgeHook,
                                              self._postRemoveBridgeHook))

    def detachHooks(self):
        """Detaches / removes callbacks from the map's hooks.
        Returns True if successful, False if already detached."""
        results = (self.nodeOverlay.detachHooks(),
                   self.edgeOverlay.detachHooks())

        if self._attachedHooks:
            for h in self._attachedHooks:
                h.disconnect()
            self._attachedHooks = None
            
            if results == (True, True):
                return True
        else:
            if results == (False, False):
                return False

    def _redisplayROIImage(self, roi):
        roiImage = self.map.labelImage().subImage(roi)
        if self._backgroundMode > 3:
            roiImage = self.faceMeans.regionImage(roiImage)
        self.viewer.replaceROI(roiImage.toPNM(BYTE),
                               qt.QPoint(*roi.upperLeft()))

    def _postMergeFacesHook(self, survivor):
        if self._backgroundMode < 3:
            return
        self._redisplayROIImage(intPos(survivor.boundingBox()))

    def _preRemoveBridgeHook(self, dart):
        if self._backgroundMode >= 3:
            self._bridgeROI = intPos(dart.edge().boundingBox())
        return True

    def _postRemoveBridgeHook(self, dart):
        if self._backgroundMode < 3:
            return
        self._redisplayROIImage(self._bridgeROI)

    def showEvent(self, e):
        self.attachHooks()
        DisplaySettings.showEvent(self, e)

    def hideEvent(self, e):
        self.detachHooks()
        DisplaySettings.hideEvent(self, e)

    def _mapToOverlays(self):
        self.edgeOverlay._calculatePoints()
        self.nodeOverlay._calculatePoints()
        self.viewer.update()

    def showMarkedEdges(self, colorMarked = qt.Qt.green, colorUnmarked = None,
                        markFlags = flag_constants.ALL_PROTECTION):
        edgeColors = [None] * self.map.maxEdgeLabel()
        for edge in self.map.edgeIter():
            if edge.flag(markFlags):
                edgeColors[edge.label()] = colorMarked
            else:
                edgeColors[edge.label()] = colorUnmarked
        self.edgeOverlay.colors = edgeColors
        self.viewer.update()

    def showAllEdges(self):
        self.edgeOverlay.colors = None
        self.viewer.update()

    def setTool(self, tool):
        """MapDisplay.setTool(tool)

        Deactivates old tool, activates new one if tool != None.
        `tool` can be either a tool object or
        1 for the MapSearcher tool
        2 for the ActivatePaintbrush
        3 for the IntelligentScissors"""

        if self._togglingGUI:
            return # no recursion please
        if self.tool:
            self.tool.disconnectViewer()
            if self.tool.parent() == self:
                self.tool.deleteLater()
            self.tool = None

        if tool == 1:
            self.tool = tools.MapSearcher(self.map, self)
        elif tool == 2:
            self.tool = tools.ActivePaintbrush(self.map, self)
        elif tool == 3:
            if not self.edgeOverlay.colors:
                self.showMarkedEdges(colorUnmarked = qt.Qt.black)
            self.tool = tools.IntelligentScissors(
                self.map, self.edgeOverlay.colors, self)
            tools.activePathMeasure = \
                tools.SimplePathCostMeasure(self.faceMeans.faceMeanDiff)
            self.nodeOverlay.visible = False
        elif hasattr(tool, "disconnectViewer"):
            self.tool = tool
        elif tool != None:
            print "setTool: invalid argument, tool deactivated now."
            print "  give 1 for MapSearcher, 2 for ActivePaintbrush, 3 for IntelligentScissors"
            tool = 0
        self._togglingGUI = True
        self.paintbrushAction.setOn(tool == 2)
        self.navigateAction.setOn(tool == 1)
        self._togglingGUI = False

    def activatePaintbrush(self, onoff = True):
        self.setTool(onoff and 2 or None)

    def activateNavigator(self, onoff = True):
        self.setTool(onoff and 1 or None)

    def cleanupMap(self):
        removeCruft(self.map, 7)

    def navigate(self, dart, center = True):
        if type(dart) == int:
            dart = self.map.dart(dart)
        self.dn = DartNavigator(dart, self)
        self.dn.show()
        if center:
            self.viewer.center = dart[0]
            self.viewer.optimizeUpperLeft()

    def setImage(self, image, pixelFormat = NBYTE):
        if hasattr(image, "orig"):
            self.backgroundGroup.setEnabled(True)
            self.preparedImage = image
            self.displayColoredAction.setEnabled(
                self.preparedImage.orig.bands() == 3)
            image = image.view
        self.image = image
        self.imageWindow.image = self.image
        self.viewer.replaceImage(self.image.toPNM(pixelFormat))

    def highlight(self, darts):
        """highlight(darts)
        Highlight the given darts (can be any iterable returning labels
        or Dart objects)."""
        self._dh.highlight(darts)

    def plotROI(self, roi,
                lineType = 1, pointType = 0,
                g = None):
        import Gnuplot
        if g == None:
            g = Gnuplot.Gnuplot()

        pi = []
        for edge in self.map.edgeIter():
            # FIXME: This was written with the old Rect2D boundingBox() in mind:
            if edge.boundingBox().intersects(roi):
                pi.append(Gnuplot.Data(
                    edge, with="lines %s" % lineType, using="(0.5+$1):(0.5+$2)",
                    title = "Edge %d" % (edge.label())))

        nodes = []
        for node in self.map.nodeIter():
            if roi.contains(intPos(node.position())):
                nodes.append(node.position())
        if not len(nodes):
            sys.stderr.write("WARNING: No nodes found in ROI %s!\n" % roi)
        else:
            pi.append(Gnuplot.Data(
                nodes, with="points pointtype %s" % pointType, using="(0.5+$1):(0.5+$2)"))

        g.set_range("xrange", (roi.left(), roi.right()))
        g.set_range("yrange", (roi.bottom(), roi.top()))
        g.plot(*pi)
        return g, pi

    def savePNG(self, filename, roi):
        image, normalize = self.imageWindow.getDisplay()
        image.subImage(roi).write(filename, normalize)

    def saveFig(self, basepath, geometry = None, scale = None,
                bgFilename = None, faceMeans = False):
        """display.saveFig(basepath,
                           geometry=None, scale=None, bgFilename = None,
                           faceMeans = False)

        Saves an XFig file as <basepath>.fig (and the pixel background
        as <basepath>_bg.png if bgFilename is not given) and returns
        the FigExporter object.

        If geometry is not None, it determines the ROI to be saved.
        geometry can be a BoundingBox or Rect2D object, a string like
        "10,10-40,30" (see fig.parseGeometry) or a tuple usable as
        Rect2D constructor arguments.

        scale is the size of a pixel in fig units (450 = 1cm) and
        defaults to a size resulting in approx. 20cm width of the
        output file.

        If bgFilename is given, it is supposed to be a background
        filename for the picture box in the XFig file, or False if no
        background is desired."""
        
        figFilename = basepath + ".fig"
        pngFilename = basepath + "_bg.png"

        # determine ROI to be saved
        if geometry == None:
            geometry = Rect2D(self.image.size())
        elif type(geometry) == tuple:
            geometry = Rect2D(*geometry)
        elif type(geometry) == str:
            geometry = Rect2D(*fig.parseGeometry(geometry))

        if faceMeans in (None, True):
            faceMeans = self.faceMeans

        if bgFilename == None and not faceMeans:
            # create .png background
            self.savePNG(pngFilename, geometry)
            _, bgFilename = os.path.split(pngFilename)

        # create .fig file
        if scale == None:
            scale = 20*450 / geometry.width() # default: 20cm width
            print "auto-adjusted scale to %s." % (scale, )
        roi = BoundingBox(geometry)
        fe = figexport.FigExporter(scale, roi)
        qtColor2figColor = figexport.qtColor2figColor
        if bgFilename and not faceMeans:
            fe.addBackgroundWithFrame(bgFilename, depth = 100, roi = roi)
        else:
            fe.addROIRect(depth = 100, roi = roi, lineWidth = 0)

        if faceMeans:
            fe.addMapFaces(self.map, faceMeans)
        
        depth = 50
        for overlay in self.viewer.overlays:
            if not overlay.visible:
                continue
            # FIXME: str(type(overlay)).contains(...) instead?
            if type(overlay) in (MapNodes, MapEdges):
                extraZoom = float(overlay._zoom) / self.viewer.scale
                oldScale, oldOffset, oldROI = fe.scale, fe.offset, fe.roi
                fe.scale *= extraZoom
                fe.roi = BoundingBox(fe.roi.begin() / extraZoom,
                                     fe.roi.end() / extraZoom)
                #fe.offset = fe.offset + Vector2(1, 1) * (extraZoom / 2.0 - 0.5)
                if type(overlay) == MapNodes:
                    radius = overlay.origRadius
                    if not overlay.relativeRadius:
                        radius /= float(overlay._zoom)
                    color = qtColor2figColor(overlay.color, fe.f)
                    fe.addMapNodes(self.map, radius,
                                   fillColor = color, lineWidth = 0, depth = depth)
                else:
                    attr = {"depth" : depth}
                    if not overlay.colors:
                        attr["penColor"] = qtColor2figColor(overlay.color, fe.f)
                    if overlay.width:
                        attr["lineWidth"] = overlay.width
                    if overlay.colors:
                        for edge in overlay._map().edgeIter():
                            edgeColor = overlay.colors[edge.label()]
                            if edgeColor:
                                parts = fe.addClippedPoly(edge,
                                    penColor = qtColor2figColor(edgeColor, fe.f),
                                    **attr)
                    else:
                        fe.addMapEdges(overlay._map(), **attr)
                fe.scale, fe.offset, fe.roi = oldScale, oldOffset, oldROI
            elif type(overlay) == PointOverlay:
                fe.addPointOverlay(overlay, depth = depth)
            elif type(overlay) == EdgeOverlay:
                fe.addEdgeOverlay(overlay, depth = depth)
            elif type(overlay) == CircleOverlay:
                fe.addCircleOverlay(overlay, depth = depth)
            elif type(overlay) == ROISelector:
                color = qtColor2figColor(overlay.color, fe.f)
                fe.addROIRect(overlay.roi, penColor = color, depth = depth)
            else:
                sys.stderr.write(
                    "MapDisplay.saveFig: overlay type %s not handled!\n" % (
                    type(overlay)))
            depth -= 1
        fe.save(figFilename)

        return fe

    def saveEPS(self, basepath, *args, **kwargs):
        """display.saveEPS(basepath, geometry=None, scale=None)

        Saves an XFig file as <basepath>.fig (see saveFig()
        documentation for details) and calls fig2dev to create an
        additional <basepath>.eps."""

        fe = self.saveFig(basepath, *args, **kwargs)
        fe.f.fig2dev(lang = "eps")

    def savePDF(self, basepath, *args, **kwargs):
        """display.savePDF(basepath, geometry=None, scale=None)

        Saves an XFig file as <basepath>.fig (see saveFig()
        documentation for details) and calls fig2dev to create an
        additional <basepath>.pdf."""

        fe = self.saveFig(basepath, *args, **kwargs)
        fe.f.fig2dev(lang = "pdf")

# --------------------------------------------------------------------
#                         dart navigation dialog
# --------------------------------------------------------------------

class DartNavigator(DartNavigatorBase):
    def __init__(self, dart, parent, name = None):
        DartNavigatorBase.__init__(self, parent, name)
        self.dart = dart
        self.connect(self.nextPhiButton, qt.SIGNAL("clicked()"),
                     self.nextPhi)
        self.connect(self.prevPhiButton, qt.SIGNAL("clicked()"),
                     self.prevPhi)
        self.connect(self.nextAlphaButton, qt.SIGNAL("clicked()"),
                     self.nextAlpha)
        self.connect(self.nextSigmaButton, qt.SIGNAL("clicked()"),
                     self.nextSigma)
        self.connect(self.prevSigmaButton, qt.SIGNAL("clicked()"),
                     self.prevSigma)

        self.connect(self.continuousCheckBox, qt.SIGNAL("toggled(bool)"),
                     self.toggleContinuous)

        self.timer = qt.QTimer(self)
        self.connect(self.timer, qt.SIGNAL("timeout()"),
                     self.highlightNext)

        self.dh = DartHighlighter(self.parent().map, self.parent().viewer)
        self.updateLabel()

    def closeEvent(self, e):
        self.dh.highlight(None)
        DartNavigatorBase.closeEvent(self, e)
        self.deleteLater() # like qt.Qt.WDestructiveClose ;-)

    def highlightNext(self):
        self.activePerm()
        self.updateLabel()

    def setDart(self, dart):
        self.dart = dart
        self.updateLabel()

    def toggleContinuous(self, onoff):
        if onoff:
            self.timer.start(1200)
        else:
            self.timer.stop()
        self.nextPhiButton.setToggleButton(onoff)
        self.prevPhiButton.setToggleButton(onoff)
        self.nextAlphaButton.setToggleButton(onoff)
        self.nextSigmaButton.setToggleButton(onoff)
        self.prevSigmaButton.setToggleButton(onoff)

    def moveDart(self, perm):
        perm()
        self.updateLabel()
        self.activePerm = perm

    def nextPhi(self):
        self.moveDart(self.dart.nextPhi)

    def prevPhi(self):
        self.moveDart(self.dart.prevPhi)

    def nextAlpha(self):
        self.moveDart(self.dart.nextAlpha)

    def nextSigma(self):
        self.moveDart(self.dart.nextSigma)

    def prevSigma(self):
        self.moveDart(self.dart.prevSigma)

    def updateLabel(self):
        self.dh.highlight(self.dart.label())
        self.dartLabel.setText(
            "%s%s\nStart node: %s\nEnd node: %s\nFaces: %d (left), %d (right)"
            % (self.dart.edge().isLoop() and "Loop" or (
                   self.dart.edge().isBridge() and "Bridge" or "Edge"),
               str(self.dart.edge())[12:-1],
               str(self.dart.startNode())[8:-1], str(self.dart.endNode())[8:-1],
               self.dart.leftFaceLabel(), self.dart.rightFaceLabel()))
        self.setCaption("DartNavigator(%d)" % (self.dart.label(), ))
        self.emit(qt.PYSIGNAL('updateDart'),(self.dart,))

# --------------------------------------------------------------------

import copy

class ROISelector(qt.QObject):
    def __init__(self, parent = None, name = None, imageSize = None,
                 roi = None, viewer = None, color = qt.Qt.yellow, width = 0):
        qt.QObject.__init__(self, parent, name)
        self._painting = False
        self._alwaysVisible = False
        self.roi = roi

        self.color = color
        self.width = width

        if viewer:
            self._viewer = viewer
        else:
            self._viewer = parent.viewer
            if imageSize == None and hasattr(parent, 'image'):
                imageSize = parent.image.size()

        self._validRect = imageSize and Rect2D(imageSize)

        self.connect(self._viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(self._viewer, qt.PYSIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(self._viewer, qt.PYSIGNAL("mouseReleased"),
                     self.mouseReleased)

        self.setVisible(True) # roi != None)

    def setVisible(self, onoff):
        if self._alwaysVisible != onoff:
            self._alwaysVisible = onoff
            if onoff:
                self._viewer.addOverlay(self)
            else:
                self._viewer.removeOverlay(self)

    def setROI(self, roi):
        if roi != self.roi:
            updateRect = self.windowRect()
            self.roi = roi
            if not self.visible:
                return
            updateRect |= self.windowRect()
            self._viewer.update(updateRect)
            self.emit(qt.PYSIGNAL("roiChanged"), (roi, ))

    def _startPainting(self):
        self._painting = True
        self._oldROI = copy.copy(self.roi)
        if not self._alwaysVisible:
            self._viewer.addOverlay(self)

    def _stopPainting(self):
        self._painting = False
        if not self._alwaysVisible:
            self._viewer.removeOverlay(self)

    def mousePressed(self, x, y, button):
        if self._painting and button == qt.Qt.RightButton:
            self._stopPainting()
            self.setROI(self._oldROI)
            return
        if button != qt.Qt.LeftButton:
            return
        if self.roi:
            mousePos = self._viewer.toWindowCoordinates(x, y)
            wr = self.windowRect()
            if (mousePos - wr.topLeft()).manhattanLength() < 9:
                self.startPos = self.roi.lowerRight() - (1,1)
            elif (mousePos - wr.bottomRight()).manhattanLength() < 9:
                self.startPos = self.roi.upperLeft()
            else:
                self.startPos = intPos((x, y))
        else:
            self.startPos = intPos((x, y))
        self.mouseMoved(x, y)
        self._startPainting()

    def mouseMoved(self, x, y):
        if not self._painting: return
        # TODO: update overlay
        x1, y1 = self.startPos
        x, y = intPos((x, y))
        self.setROI(
            Rect2D(min(x1, x), min(y1, y), max(x1, x)+1, max(y1, y)+1))

    def windowRect(self):
        if not self.roi:
            return qt.QRect()
        return qt.QRect(
            self._viewer.toWindowCoordinates(*self.roi.upperLeft()),
            self._viewer.toWindowCoordinates(*self.roi.lowerRight()))

    def mouseReleased(self, x, y, button):
        if self._painting and button == qt.Qt.LeftButton:
            self._stopPainting()
            self.setROI(self.roi & self._validRect)
            self.emit(qt.PYSIGNAL("roiSelected"), (self.roi, ))

    def disconnectViewer(self):
        self.disconnect(self._viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(self._viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(self._viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)
        if self._alwaysVisible:
            self._viewer.removeOverlay(self)

    def setZoom(self, zoom):
        self.zoom = zoom

    def draw(self, p):
        if not self.roi:
            return
        p.setPen(qt.QPen(self.color, self.width))
        p.setBrush(qt.Qt.NoBrush)
        drawRect = self.windowRect()
        # painter is already set up with a shift:
        drawRect.moveBy(-self._viewer.x, -self._viewer.y)
        p.drawRect(drawRect)

# def queryROI(imageWindow):
#     pd = qt.QProgressDialog(imageWindow)
#     pd.setLabelText("Please mark a ROI with drag & drop")
#     rs = ROISelector(imageWindow)
#     qt.QObject.connect(rs, qt.PYSIGNAL("roiSelected"),
#                        pd.accept)
#     if pd.exec_loop() == qt.QDialog.Accepted:
#         return rs.roi
#     return
