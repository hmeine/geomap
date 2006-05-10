_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import figexport, qt, sys, os, time, tools
from vigra import BYTE, NBYTE, Point2D, Rect2D, Vector2
from vigrapyqt import ImageWindow, EdgeOverlay, PointOverlay
from hourglass import simplifyPolygon, intPos
from map import removeCruft, mergeZeroPixelFaces
from weakref import ref

# ui-generated base classes:
from displaysettings import DisplaySettings
from dartnavigator import DartNavigatorBase

def findZoomFactor(srcSize, destSize):
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
        self.width = width
        self.useIndividualColors = False
        self._zoomedEdges = None
        self._zoom = None

    def setMap(self, map):
        self._map = ref(map)
        self._dirty = True

    def _calculatePoints(self):
        if not self._map():
            return
        c = time.clock()
        self._resizeCount = 0
        edges = self._map().edges
        if not self._zoomedEdges or len(self._zoomedEdges) != len(edges):
            self._zoomedEdges = [None] * len(edges)
        for label, edge in enumerate(edges):
            if edge:
                self._calculateEdgePoints(edge._label, edge)
            else:
                self._zoomedEdges[label] = None
        sys.stdout.write("MapEdges._calculatePoints(zoom = %s) took %ss "
                         "(%d undesirable resizes).\n" % (
            self._zoom, time.clock() - c, self._resizeCount))
        self._dirty = False

    def _calculateEdgePoints(self, index, origEdgePoints):
        qpa = self._zoomedEdges[index]

        origEdgePoints = simplifyPolygon(origEdgePoints * self._zoom, 0.5) \
                         .roundToInteger()

        if qpa == None:
            qpa = qt.QPointArray(len(origEdgePoints))
        elif qpa.size() < len(origEdgePoints):
            qpa.resize(len(origEdgePoints))

        d0 = Point2D(int(self._zoom / 2), int(self._zoom / 2))
        points = iter(origEdgePoints)
        oldPos = points.next()
        qpa.setPoint(0, *(oldPos + d0))

        qpaIndex = 1
        for pos in points:
            if pos != oldPos:
                qpa.setPoint(qpaIndex, *(pos + d0))
                qpaIndex += 1
                oldPos = pos

        # extra parameter in c++ would be QGArray.SpeedOptim (1), but
        # is not available here:
        if qpaIndex != qpa.size():
            qpa.resize(qpaIndex)
            self._resizeCount += 1

        self._zoomedEdges[index] = qpa

    def setEdges(self, edges):
        self.originalEdges = edges
        self.setZoom(self._zoom)

    def setEdgePoints(self, index, edgePoints):
        if not self._zoomedEdges:
            return
        if self._zoomedEdges[index]:
            updateROI = self._zoomedEdges[index].boundingRect()
        else:
            updateROI = qt.QRect()
        if edgePoints:
            self._calculateEdgePoints(index, edgePoints)
            updateROI |= self._zoomedEdges[index].boundingRect()
        else:
            self._zoomedEdges[index] = None
        if self.visible:
            updateROI.moveBy(self.viewer.x, self.viewer.y)
            self.viewer.update(updateROI)

    def preRemoveEdgeHook(self, dart):
        self.setEdgePoints(dart.edgeLabel(), None)

    def preMergeEdgesHook(self, dart):
        self.setEdgePoints(
            dart.clone().nextSigma().edgeLabel(), None)

    def postMergeEdgesHook(self, edge):
        self.setEdgePoints(edge._label, edge)

    def setZoom(self, zoom):
        if self._map():
            zoom *= findZoomFactor(
                self._map().imageSize()[0],
                self.viewer.pixmap.size().width())
        if self._zoom != zoom:
            self._zoom = zoom
            self._dirty = True

    def draw(self, p):
        if not self._map():
            return
        if self._dirty:
            self._calculatePoints()
        r = p.clipRegion()
        r.translate(-self.viewer.x, -self.viewer.y)
        edges = self._map().edges
        if self.useIndividualColors:
            for edge, zoomedEdge in map(None, edges, self._zoomedEdges):
                if edge and hasattr(edge, "color") and edge.color:
                    p.setPen(qt.QPen(edge.color, self.width))
                    p.drawPolyline(zoomedEdge)
        else:
            p.setPen(qt.QPen(self.color, self.width))
            for e in self._zoomedEdges:
                if e and (r.contains(e.boundingRect()) or r.isNull()):
                    p.drawPolyline(e)

class MapNodes(object):
    def __init__(self, map, color,
                 radius = 1, relativeRadius = False):
        self.setMap(map)
        self.color = color
        #self.useIndividualColors = False
        self.setRadius(radius, relativeRadius)
        self._zoom = None

    def setMap(self, map):
        self._map = ref(map)
        self._qpointlist = None

    def _calculatePoints(self):
        if not self._map():
            self._qpointlist = []
            return
        c = time.clock()
        if self.relativeRadius:
            self.radius = int(self._zoom * self.origRadius + 0.5)
        d0 = Vector2(0.5 * self._zoom - self.radius, 0.5 * self._zoom - self.radius)
        w = 2 * self.radius + 1
        self.s = qt.QSize(w, w)
        self._qpointlist = [None] * len(self._map().nodes)
        for node in self._map().nodeIter():
            ip = intPos(node.position() * self._zoom + d0)
            self._qpointlist[node._label] = qt.QPoint(ip[0], ip[1])
        sys.stdout.write("MapNodes._calculatePoints(zoom = %s) took %ss.\n" % (
            self._zoom, time.clock() - c))

    def removeNode(self, node):
        if self._qpointlist:
            if not self._qpointlist[node._label]:
                sys.stderr.write("WARNING: MapNodes.removeNode(): Node already None!\n")
                return
            if self.visible:
                ur = qt.QRect(self._qpointlist[node._label], self.s)
                ur.moveBy(self.viewer.x, self.viewer.y)
                self.viewer.update(ur)
            self._qpointlist[node._label] = None

    def setZoom(self, zoom):
        if self._map():
            zoom *= findZoomFactor(
                self._map().imageSize()[0],
                self.viewer.pixmap.size().width())
        if self._zoom != zoom:
            self._zoom = zoom
            self._qpointlist = None

    def setRadius(self, radius, relativeRadius = False):
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

        if type(darts) == int:
            darts = [darts]
        if darts:
            dartObjects = []
            for dart in darts:
                if type(dart) == int:
                    dartObjects.append(self._map().dart(dart))
                else:
                    dartObjects.append(dart)
            darts = dartObjects

        if self.eo != None:
            self._viewer.removeOverlay(self.eo)
            self._viewer.removeOverlay(self.no)
            self.eo = None
            self.no = None

        if darts == None or not len(darts):
            return

        self.eo = EdgeOverlay([dart.edge() for dart in darts], color)
        self.no = PointOverlay(
            [dart.startNode().position() for dart in darts], color, 2)
        self._viewer.addOverlay(self.eo)
        self._viewer.addOverlay(self.no)

# --------------------------------------------------------------------
#                             MapDisplay
# --------------------------------------------------------------------

class MapDisplay(DisplaySettings):
    def __init__(self, map, preparedImage = None, immediateShow = True):
        DisplaySettings.__init__(self)
        # for backward compatibility:
        if hasattr(preparedImage, "imageSize") and hasattr(map, "width"):
            map, preparedImage = preparedImage, map
        elif preparedImage == None:
            preparedImage = map.labelImage
        
        self.preparedImage = preparedImage
        self.map = map
        self._togglingGUI = False

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
            sys.modules["map"]._cvsVersion[0], ))

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
        
        # I wonder why this makes sense; obviously, the callback hooks
        # do not prevent the garbage collector from deleting this
        # object:
        self.detachHooks()

    def setMap(self, map):
        attached = self.detachHooks()
        self.setTool(None)
        self.map = map
        self.edgeOverlay.setMap(map)
        self.nodeOverlay.setMap(map)
        self._dh.setMap(map)
        if attached:
            self.attachHooks()
        self.viewer.update()

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
            displayImage = self.map.labelImage
        elif mode == 4:
            displayImage = faceImage(self.map)
        else:
            sys.stderr.write("Unknown background mode %d!\n" % mode)
            return

        self.image = displayImage
        self._backgroundMode = mode # needed by toggleNormalize()
        normalize = self.normalizeStates[mode]
        if self.normalizeAction.isOn() != normalize:
            self.normalizeAction.setOn(normalize)
        else:
            self.setImage(self.image, normalize and NBYTE or BYTE)

    def toggleNormalize(self, normalize):
        self.normalizeStates[self._backgroundMode] = normalize
        self.setImage(self.image, normalize and NBYTE or BYTE)

    def attachHooks(self):
        map = self.map
        map.removeNodeHooks.append(self.nodeOverlay.removeNode)
        map.preRemoveBridgeHooks.append(self.edgeOverlay.preRemoveEdgeHook)
        map.preMergeFacesHooks.append(self.edgeOverlay.preRemoveEdgeHook)
        map.preMergeEdgesHooks.append(self.edgeOverlay.preMergeEdgesHook)
        map.postMergeEdgesHooks.append(self.edgeOverlay.postMergeEdgesHook)

    def detachHooks(self):
        """Detaches / removes callbacks from the map's hooks.
        Returns True if successful, False if already detached."""
        map = self.map
        #if not map: return # Map already destroyed (when using weakrefs)
        try:
            map.removeNodeHooks.remove(self.nodeOverlay.removeNode)
        except ValueError:
            return False # already detached
        map.preRemoveBridgeHooks.remove(self.edgeOverlay.preRemoveEdgeHook)
        map.preMergeFacesHooks.remove(self.edgeOverlay.preRemoveEdgeHook)
        map.preMergeEdgesHooks.remove(self.edgeOverlay.preMergeEdgesHook)
        map.postMergeEdgesHooks.remove(self.edgeOverlay.postMergeEdgesHook)
        return True

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

    def setTool(self, tool):
        if self._togglingGUI:
            return # no recursion please
        if self.tool:
            self.tool.disconnectViewer()
            self.tool = None

        if tool == 1:
            self.tool = tools.MapSearcher(self.map, self)
        elif tool == 2:
            self.tool = tools.ActivePaintbrush(self.map, self)
        elif tool == 3:
            self.tool = tools.IntelligentScissors(self.map, self)
            for edge in self.map.edgeIter():
                if not hasattr(edge, "color"):
                    edge.color = qt.Qt.black
            self.edgeOverlay.useIndividualColors = True
            self.nodeOverlay.visible = False
        elif type(tool) != int:
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
        removeCruft(self.map, 6)
        mergeZeroPixelFaces(self.map)

    def navigate(self, dart):
        if type(dart) == int:
            dart = self.map.dart(dart)
        self.dn = DartNavigator(dart, self)
        self.dn.show()

    def setImage(self, image, pixelFormat = NBYTE):
        self.image = image
        self.imageWindow.image = image
        self.viewer.replaceImage(image.toPNM(pixelFormat))

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
        self.image.subImage(roi).write(
            filename, self.normalizeStates[self._backgroundMode] and NBYTE or BYTE)

    def saveFig(self, basepath, geometry = None, scale = None):
        figFilename = basepath + ".fig"
        pngFilename = basepath + "_bg.png"
        epsFilename = basepath + ".eps"

        path, figBasename = os.path.split(figFilename)
        path, pngBasename = os.path.split(pngFilename)
        path, epsBasename = os.path.split(epsFilename)

        if geometry == None:
            geometry = Rect2D(self.image.size())
        elif type(geometry) == tuple:
            geometry = Rect2D(*geometry)
        elif type(geometry) == str:
            geometry = Rect2D(*fig.parseGeometry(geometry))

        # create .png background
        self.savePNG(pngFilename, geometry)

        # create .fig file
        if scale == None:
            scale = 20*450 / geometry.width() # default: 20cm width
        roi = BoundingBox(geometry)
        fe = figexport.FigExporter(scale, roi)
        qtColor2figColor = figexport.qtColor2figColor
        fe.addBackgroundWithFrame(pngFilename, depth = 100, roi = roi)
        depth = 50
        for overlay in self.viewer.overlays:
            if not overlay.visible:
                continue
            # FIXME: str(type(overlay)).contains(...) instead?
            if type(overlay) == MapNodes:
                radius = overlay.origRadius
                if not overlay.relativeRadius:
                    radius /= float(overlay._zoom)
                color = qtColor2figColor(overlay.color, fe.f)
                fe.addMapNodes(self.map, radius,
                               fillColor = color, lineWidth = 0, depth = depth)
            elif type(overlay) == MapEdges:
                fe.addMapEdges(self.map, penColor = qtColor2figColor(overlay.color, fe.f),
                               depth = depth)
            elif type(overlay) == PointOverlay:
                fe.addPointOverlay(overlay, depth = depth)
            elif type(overlay) == EdgeOverlay:
                fe.addEdgeOverlay(overlay, depth = depth)
            else:
                sys.stderr.write(
                    "MapDisplay.saveFig: overlay type %s not handled!\n" % (
                    type(overlay)))
            depth -= 1
        fe.save(figFilename)

        # create .eps output
        cin, cout = os.popen4("%sfig2dev -L eps '%s' '%s'" % (
            path and ("cd '%s' && " % path) or "", figBasename, epsBasename))
        cin.close(); print cout.read(),
        del cin, cout

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

        self.dh = DartHighlighter(dart._map(), self.parent().viewer)
        self.updateLabel()

    def closeEvent(self, e):
        self.dh.highlight(None)
        DartNavigatorBase.closeEvent(self, e)

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
               str(self.dart.edge())[9:-1],
               str(self.dart.startNode())[5:-1], str(self.dart.endNode())[5:-1],
               self.dart.leftFaceLabel(), self.dart.rightFaceLabel()))
        self.setCaption("DartNavigator(%d)" % (self.dart.label(), ))
        self.emit(qt.PYSIGNAL('updateDart'),(self.dart,))

# --------------------------------------------------------------------

class ROISelector(qt.QObject):
    def __init__(self, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.painting = False

        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(viewer, qt.PYSIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(viewer, qt.PYSIGNAL("mouseReleased"),
                     self.mouseReleased)

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        self.startPos = Point2D(x, y)
        self.painting = True
        self.mouseMoved(x, y)
        # TODO: add overlay

    def mouseMoved(self, x, y):
        if not self.painting: return
        # TODO: update overlay

    def mouseReleased(self, x, y):
        self.painting = False
        # TODO: remove overlay
        self.roi = Rect2D(self.startPos, Point2D(x, y))

    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)
