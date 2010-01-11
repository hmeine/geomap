import sys, math, time, qt
import fig, figexport
import vigra, vigrapyqt
import maputils, flag_constants, tools, statistics
from vigra import Rect2D, Vector2
from geomap import simplifyPolygon, intPos, BoundingBox, contourPoly
from maputils import removeCruft
from weakref import ref

# ui-generated base classes:
import displaysettings, dartnavigator

assert bool(vigra.NBYTE) and not bool(vigra.BYTE), \
       "this code relies on the fact that `normalize` can be either a bool or a PixelType"

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

class MapFaces(vigrapyqt.Overlay):
    __base = vigrapyqt.Overlay

    # TODO: remove _map (reuse from edgeOverlay)

    def __init__(self, map, edgeOverlay,
                 color = None, width = 0, fillColor = None,
                 flags = None):
        self.__base.__init__(self, color, width, fillColor)
        self.setMap(map)
        self.edgeOverlay = edgeOverlay
        self.flags = flags
        self.colors = None

    def setMap(self, map):
        self._map = ref(map)

    def _faceROI(self, face):
        # FIXME: very inefficient (use face.boundingBox() instead?!)
        return self._getZoomedFace(face)[0].boundingRect()

    def updateFaceROI(self, face):
        """Update the ROI containing the given face(s) in the viewer.
        `face` may be an int, an Face object, or a tuple/list
        of one of the two.  (Passing multiple faces results in only
        one call to viewer.update(ROI).)"""
        
        if isinstance(face, (tuple, list)):
            roi = qt.QRect()
            for face in face:
                roi |= self._faceROI(face)
        else:
            roi = self._faceROI(face)
        roi.moveBy(self.viewer.upperLeft().x(),
                   self.viewer.upperLeft().y())
        self.viewer.update(roi)

    def _getZoomedFace(self, face):
        """Return the QPointArray for the given GeoMap.Face object."""

        zoomedEdge = self.edgeOverlay._getZoomedEdge
        edges = []
        sizes = []
        for anchor in face.contours():
            contourSize = 0
            for dart in anchor.phiOrbit():
                ze = zoomedEdge(dart.edge())
                edges.append((ze, dart.label() < 0))
                contourSize += ze.size()
            sizes.append(contourSize)

        qpa = qt.QPointArray(sum(sizes))
        ti = 0
        for edge, reverse in edges:
            if reverse:
                ti += edge.size()
                for i in range(edge.size()):
                    qpa.setPoint(ti-i-1, edge.at(i))
            else:
                for i in range(edge.size()):
                    qpa.setPoint(ti+i, edge.at(i))
                ti += edge.size()
        return qpa, sizes

    def setZoom(self, zoom):
        pass # zoom is managed via self.edgeOverlay
    
    def draw(self, p):
        if not self._map():
            return

        self._setupPainter(p)

        zoomFactor = self.edgeOverlay._zoom

        #p.drawPolygon(self._qpointarray, True, s, e-s)

        r = p.clipRegion().boundingRect()
        r.moveBy(-self.viewer.upperLeft().x(),
                 -self.viewer.upperLeft().y())
        bbox = BoundingBox(Vector2(r.left() / zoomFactor - 0.5,
                                   r.top() / zoomFactor - 0.5),
                           Vector2(r.right() / zoomFactor + 0.5,
                                   r.bottom() / zoomFactor + 0.5))
        map = self._map()
        if self.colors:
            try:
                for edge in map.edgeIter():
                    edgeColor = self.colors[edge.label()]
                    if edgeColor and bbox.intersects(edge.boundingBox()):
                        p.setPen(qt.QPen(edgeColor, self.width))
                        p.drawPolyline(self._getZoomedEdge(edge))
            except IndexError, e:
                print e #"IndexError: %d > %d (maxEdgeLabel: %d)!" % (
                    #i, len(self.colors), map.maxEdgeLabel())
        elif self.color or self.fillColor:
            for face in map.faceIter():
                if face.flag(self.flags) and \
                   (not face.label() or \
                    bbox.intersects(face.boundingBox())):
                    qpa, sizes = self._getZoomedFace(face)
                    if face.holeCount():
                        index = 0
                        if self.color:
                            for size in sizes:
                                p.drawPolyline(qpa, index, size)
                                index += size
                        if self.fillColor:
                            p.save()
                            p.setPen(qt.Qt.NoPen)
                            p.drawPolygon(qpa)
                            p.restore()
                    else:
                        p.drawPolygon(qpa)

class MapEdges(vigrapyqt.Overlay):
    __slots__ = ("colors", "protectedColor", "protectedWidth",
                 "_map", "_attachedHooks", "_removedEdgeLabel",
                 "_zoom", "_zoomedEdges")

    __base = vigrapyqt.Overlay
    
    def __init__(self, map, color, width = 0,
                 protectedColor = None, protectedWidth = None,
                 name = None):
        self.__base.__init__(self, color, width, name = name)
        self.setMap(map)
        self.colors = None
        self.protectedColor = protectedColor
        self.protectedWidth = protectedWidth
        self._zoom = None
        self._attachedHooks = None

    def setMap(self, map):
        self._map = ref(map)
        self._zoomedEdges = [None] * map.maxEdgeLabel()

    def attachHooks(self):
        map = self._map()
        self._attachedHooks = (
            map.addRemoveBridgeCallbacks(self._preRemoveEdgeHook, self._postRemoveEdgeHook),
            map.addMergeFacesCallbacks(self._preRemoveEdgeHook, self._postRemoveEdgeHook),
            map.addMergeEdgesCallbacks(self._preMergeEdgesHook, self._postMergeEdgesHook))

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
        """Change the color of `edge` to `newColor` and update the
        viewer.  `edge` may be an int, an Edge object, or a tuple/list
        of one of the two.  (Passing multiple edges results in only
        one call to viewer.update(ROI).)"""
        
        if not self.colors:
            self.colors = [self.color] * self._map().maxEdgeLabel()

        if isinstance(edge, (tuple, list)):
            if isinstance(edge[0], int):
                edge = map(self._map().edge, edge)
            updateROI = qt.QRect()
            for edge in edge:
                self.colors[edge.label()] = newColor
                updateROI |= self._edgeROI(edge)
        else:
            if isinstance(edge, int):
                edge = self._map().edge(edge)
            self.colors[edge.label()] = newColor
            updateROI = self._edgeROI(edge)

        self.viewer.update(updateROI)

    def _edgeROI(self, edge):
        return self._getZoomedEdge(edge).boundingRect()

    def updateEdgeROI(self, edge):
        """Update the ROI containing the given edge(s) in the viewer.
        `edge` may be an int, an Edge object, or a tuple/list
        of one of the two.  (Passing multiple edges results in only
        one call to viewer.update(ROI).)"""
        
        if isinstance(edge, (tuple, list)):
            roi = qt.QRect()
            for edge in edge:
                roi |= self._edgeROI(edge)
        else:
            roi = self._edgeROI(edge)
        roi.moveBy(self.viewer.upperLeft().x(),
                   self.viewer.upperLeft().y())
        self.viewer.update(roi)

    def _calculateZoomedEdge(self, edgeLabel, edge):
        offset = Vector2(self._zoom / 2.0 - 0.5, self._zoom / 2.0 - 0.5)
        origEdgePoints = (
            simplifyPolygon(edge * self._zoom, 0.1)
            + offset).roundToInteger()

        qpa = self._zoomedEdges[edgeLabel]
        if qpa == None:
            qpa = qt.QPointArray(len(origEdgePoints))
            self._zoomedEdges[edgeLabel] = qpa
        elif qpa.size() != len(origEdgePoints):
            # extra parameter in c++ would be QGArray.SpeedOptim (1), but
            # is not available here:
            qpa.resize(len(origEdgePoints))

        for i, pos in enumerate(origEdgePoints):
            qpa.setPoint(i, pos[0], pos[1])

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
            # TODO: knowing the bounding box for the update would be
            # enough if the zoomed edge was marked dirty:
            qpa = self._calculateZoomedEdge(index, edgePoints)
            updateROI |= qpa.boundingRect()
        else:
            self._zoomedEdges[index] = None

        updateROI.moveBy(self.viewer.upperLeft().x(),
                         self.viewer.upperLeft().y())
        self.viewer.update(updateROI)

    def _preRemoveEdgeHook(self, dart):
        self._removedEdgeLabel = dart.edgeLabel()
        return True

    def _postRemoveEdgeHook(self, survivor):
        self.setEdgePoints(self._removedEdgeLabel, None)

    def _preMergeEdgesHook(self, dart):
        self._removedEdgeLabel = dart.clone().nextSigma().edgeLabel()
        return True

    def _postMergeEdgesHook(self, edge):
        self.setEdgePoints(self._removedEdgeLabel, None)
        self.setEdgePoints(edge.label(), edge)

    def setZoom(self, zoom):
        if self._map() and self._map().imageSize()[0]:
            zoom *= findZoomFactor(
                self._map().imageSize()[0],
                self.viewer.pixmap.size().width())
        if self._zoom != zoom:
            self._zoom = zoom
            self._zoomedEdges = [None] * len(self._zoomedEdges)

    def _getZoomedEdge(self, edge):
        """Return the QPointArray for the given GeoMap.Edge object.
        The QPointArrays are cached and will be created on demand if
        you call this function."""
        index = edge.label()
        result = self._zoomedEdges[index]
        if result == None:
            result = self._calculateZoomedEdge(index, edge)
        return result
    
    def draw(self, p):
        if not self._map():
            return

        self._setupPainter(p)

        r = p.clipRegion().boundingRect()
        r.moveBy(-self.viewer.upperLeft().x(),
                 -self.viewer.upperLeft().y())
        bbox = BoundingBox(Vector2(r.left() / self._zoom - 0.5,
                                   r.top() / self._zoom - 0.5),
                           Vector2(r.right() / self._zoom + 0.5,
                                   r.bottom() / self._zoom + 0.5))
        map = self._map()
        if self.colors:
            try:
                for edge in map.edgeIter():
                    edgeColor = self.colors[edge.label()]
                    if edgeColor and bbox.intersects(edge.boundingBox()):
                        p.setPen(qt.QPen(edgeColor, self.width))
                        p.drawPolyline(self._getZoomedEdge(edge))
            except IndexError, e:
                print e #"IndexError: %d > %d (maxEdgeLabel: %d)!" % (
                    #i, len(self.colors), map.maxEdgeLabel())
        elif self.color:
            for edge in map.edgeIter():
                if bbox.intersects(edge.boundingBox()):
                    p.drawPolyline(self._getZoomedEdge(edge))

        if self.protectedColor:
            p.setPen(qt.QPen(self.protectedColor,
                             self.protectedWidth or self.width))
            for edge in map.edgeIter():
                if edge.flag(flag_constants.ALL_PROTECTION) and \
                       bbox.intersects(edge.boundingBox()):
                    p.drawPolyline(self._getZoomedEdge(edge))

class MapNodes(vigrapyqt.Overlay):
    __base = vigrapyqt.Overlay
    
    def __init__(self, map, color, radius = 0.2, relativeRadius = True):
        self.__base.__init__(self, color)
        self.setMap(map)
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
        self._elSize = qt.QSize(w, w)
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
                ur = qt.QRect(self._qpointlist[nodeLabel], self._elSize)
                ur.moveBy(self.viewer.upperLeft().x(),
                          self.viewer.upperLeft().y())
                self.viewer.update(ur)
            self._qpointlist[nodeLabel] = None
        return True

    def setZoom(self, zoom):
        if self._map() and self._map().imageSize()[0]:
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

        self._setupPainter(p)

        if self.radius == 0:
            for point in self._qpointlist:
                if point: # TODO: boundingRect
                    #p.drawPoint(point)
                    p.drawLine(point, point)
        else:
            p.setBrush(qt.QBrush(self.color))
            for point in self._qpointlist:
                if point: # TODO: boundingRect
                    p.drawEllipse(qt.QRect(point, self._elSize))

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
                    if not dart.edge():
                        sys.stderr.write("WARNING: Cannot highlight nonexisting Dart %d!\n" % dart.label())
                        continue
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

        self.eo = vigrapyqt.EdgeOverlay([dart.edge() for dart in darts], color)
        self.eo.width = 2
        self.no = vigrapyqt.PointOverlay(
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

def addMapOverlay(fe, overlay, skipBorder = False, **attr):
    qtColor2figColor = figexport.qtColor2figColor

    # FIXME: str(type(overlay)).contains(...) instead?
    if isinstance(overlay, ROISelector):
        color = qtColor2figColor(overlay.color, fe.f)
        return fe.addROIRect(overlay.roi, penColor = color, **attr)
    elif isinstance(overlay, (MapNodes, MapEdges, MapFaces)):
        oldScale, oldOffset, oldROI = fe.scale, fe.offset, fe.roi

        if isinstance(overlay, MapFaces):
            zoom = overlay.edgeOverlay._zoom
        else:
            zoom = overlay._zoom
        extraZoom = float(zoom) / overlay.viewer.zoomFactor()
        fe.scale *= extraZoom
        fe.roi = BoundingBox(fe.roi.begin() / extraZoom,
                             fe.roi.end() / extraZoom)
        map = overlay._map()

        if isinstance(overlay, MapNodes):
            radius = overlay.origRadius
            if not overlay.relativeRadius:
                radius /= float(overlay._zoom)
            color = qtColor2figColor(overlay.color, fe.f)

            result = fe.addMapNodes(map, radius,
                                    fillColor = color, lineWidth = 0, **attr)
        elif isinstance(overlay, MapEdges):
            attr = dict(attr)
            if overlay.width:
                attr["lineWidth"] = overlay.width

            if overlay.colors:
                result = fig.Compound(fe.f)
                for edge in map.edgeIter():
                    edgeColor = overlay.colors[edge.label()]
                    if edgeColor:
                        fe.addClippedPoly(edge,
                            penColor = qtColor2figColor(edgeColor, fe.f),
                            container = result, **attr)
            elif overlay.color:
                result = fe.addMapEdges(
                    map, skipBorder = skipBorder,
                    penColor = qtColor2figColor(overlay.color, fe.f),
                    **attr)
            else:
                result = fig.Compound(fe.f)

            if overlay.protectedColor:
                attr["penColor"] = \
                    qtColor2figColor(overlay.protectedColor, fe.f)
                attr["lineWidth"] = overlay.protectedWidth or overlay.width
                it = skipBorder and maputils.nonBorderEdges(map) \
                     or map.edgeIter()
                for edge in it:
                    if edge.flag(flag_constants.ALL_PROTECTION):
                        fe.addClippedPoly(edge, container = result, **attr)
        else:
            attr = dict(attr)
            if overlay.width:
                attr["lineWidth"] = overlay.width
            attr["penColor"] = qtColor2figColor(overlay.color, fe.f)
            if overlay.fillColor:
                attr["fillColor"] = qtColor2figColor(overlay.fillColor, fe.f)
                attr["fillStyle"] = fig.fillStyleSolid
            result = fig.Compound(fe.f)
            for face in map.faceIter():
                if face.flag(overlay.flags):
                    if face.holeCount:
                        assert not overlay.fillColor, "FIXME: cannot currently export filled polygons with holes"
                    for dart in face.contours():
                        fe.addClippedPoly(contourPoly(dart),
                                          container = result, **attr)

        fe.scale, fe.offset, fe.roi = oldScale, oldOffset, oldROI
        return result
    else:
        return figexport.addStandardOverlay(fe, overlay, **attr)

class MapDisplay(displaysettings.DisplaySettings):
    __base = displaysettings.DisplaySettings

    # actually, this has only documenting effect; since this inherits
    # PyQt widgets, any attribute may be used:
    __slots__ = ("tool", "viewer", "images", "image",
                 "map", "nodeOverlay", "edgeOverlay",
                 "addOverlay", "replaceOverlay", "removeOverlay",
                 "_togglingGUI", "_imageWindow", "_backgroundMode", "_normalizeStates",
                 "_faceMeans", "_attachedHooks",
                 "dn", "_dh")
    
    def __init__(self, map, preparedImage = None, immediateShow = True,
                 faceMeans = None):
        self.__base.__init__(self)
        self.tool = None
        self._togglingGUI = False
        self._faceMeans = None # temp. needed (for _enableImageActions)

        # for backward compatibility:
        if hasattr(preparedImage, "imageSize") and hasattr(map, "width"):
            map, preparedImage = preparedImage, map
        elif preparedImage == None:
            preparedImage = map.labelImage()
            if not preparedImage:
                preparedImage = vigra.GrayImage(map.imageSize())

        self._imageWindow = None # setImage would norm. pass the image on
        if hasattr(preparedImage, "orig"):
            self.images = {
                "original" : preparedImage.view,
                "colored" : preparedImage.orig,
                "bi" : preparedImage.bi.gm,
                }
            preparedImage = preparedImage.view
        else:
            self.images = {}
            # auto-detects role colored/original:
            self.setImage(preparedImage, normalize = False)

        self.image = preparedImage
        self._imageWindow = vigrapyqt.ImageWindow(self.image, vigra.BYTE, self)
        self._imageWindow.label.hide()
        self.setCentralWidget(self._imageWindow)
        self.viewer = self._imageWindow.viewer
        self.viewer.autoZoom()

        # convenience:
        self.addOverlay = self.viewer.addOverlay
        self.replaceOverlay = self.viewer.replaceOverlay
        self.removeOverlay = self.viewer.removeOverlay
        self.overlays = self.viewer.overlays
        self.applyExpression = self._imageWindow.applyExpression

        self.map = map
        if not faceMeans and hasattr(map, "faceMeans"):
            faceMeans = map.faceMeans
        self.setFaceMeans(faceMeans)
        self._attachedHooks = None

        self._normalizeStates = [False, False, True, True, False]
        if self.image.bands() == 3:
            self._backgroundMode = 1
            self.displayColoredAction.setOn(True)
        else:
            self._backgroundMode = 0
            self.displayOriginalAction.setOn(True)
        self._enableImageActions()

        self.connect(self.backgroundGroup, qt.SIGNAL("selected(QAction*)"),
                     self.setBackgroundMode)
        self.connect(self.nodeDisplayAction, qt.SIGNAL("toggled(bool)"),
                     self.toggleNodeDisplay)
        self.connect(self.edgeDisplayAction, qt.SIGNAL("toggled(bool)"),
                     self.toggleEdgeDisplay)
        self.connect(self.normalizeAction, qt.SIGNAL("toggled(bool)"),
                     self.toggleNormalize)
        self.connect(self.mapCleanupAction, qt.SIGNAL("activated()"),
                     self.cleanupMap)
        self.connect(self.paintbrushAction, qt.SIGNAL("toggled(bool)"),
                     self.activatePaintbrush)
        self.connect(self.scissorsAction, qt.SIGNAL("toggled(bool)"),
                     self.activateScissors)
        self.connect(self.navigateAction, qt.SIGNAL("toggled(bool)"),
                     self.activateNavigator)
        self.connect(self._imageWindow, qt.PYSIGNAL("captionChanged"),
                     self.statusMessage)
                     #self.statusBar(), qt.SLOT("message(const QString&)"))

        self.setCaption("Map Display")

        self.edgeOverlay = MapEdges(map, qt.Qt.red,
                                    protectedColor = qt.Qt.green,
                                    protectedWidth = 2)
        self.viewer.addOverlay(self.edgeOverlay)
        self.edgeOverlay.visible = self.edgeDisplayAction.isOn()
        self.nodeOverlay = MapNodes(map, qt.Qt.blue)
        self.viewer.addOverlay(self.nodeOverlay)
        self.nodeOverlay.visible = self.nodeDisplayAction.isOn()
        self._dh = DartHighlighter(map, self.viewer)
        self.dn = None

        if immediateShow:
            self.show()

    def __del__(self):
        # delete tool (which may reference the viewer & map)
        self.setTool(None)

    def _dartNavigatorDestroyed(self, dn):
        """HACK: this is needed because of the deleteLater()..."""
        # FIXME: cannot reliably detect identity of half-destructed
        # object anymore:
        #if dn is self.dn:
        self.dn = None

    def setMap(self, map):
        attached = self.detachHooks()

        self.setTool(None)
        self.map = map
        self.edgeOverlay.setMap(map)
        self.nodeOverlay.setMap(map)
        self._dh.setMap(map)

        updatedDisplayImage = None
        if self._backgroundMode == 3:
            updatedDisplayImage = self.map.labelImage()

        self._faceMeans = None
        if hasattr(map, "faceMeans"):
            self.setFaceMeans(map.faceMeans)
            if self._backgroundMode == 4:
                updatedDisplayImage = \
                    self._faceMeans.regionImage(self.map.labelImage())
        self._enableImageActions()

        if updatedDisplayImage:
            self._setImage(updatedDisplayImage,
                           self._normalizeStates[self._backgroundMode])

        if attached:
            self.attachHooks()

        self.viewer.update()

    def setFaceMeans(self, faceMeans):
        self._faceMeans = faceMeans
        self._enableImageActions()

    def _adjustSize(self):
        pass # don't change window size out of a sudden

    def statusMessage(self, msg):
        # workaround - don't know why connecting the PYSIGNAL directly
        # does not work (maybe qt.PYSIGNAL("captionChanged(QString)")
        # would work?)
        self.statusBar().message(msg)

    def currentRole(self):
        return ("original", "colored", "bi", "labels", "faceMeans")[
            self._backgroundMode]

    def setBackgroundMode(self, mode):
        if type(mode) != int:
            mode = [self.displayOriginalAction,
                    self.displayColoredAction,
                    self.displayBIAction,
                    self.displayLabelsAction,
                    self.displayMeansAction].index(mode)

        displayImage = None
        if mode == 0:
            displayImage = self.images["original"]
        elif mode == 1:
            displayImage = self.images["colored"]
        elif mode == 2:
            displayImage = self.images["bi"]
        elif mode == 3:
            displayImage = self.map.labelImage()
        elif mode == 4:
            displayImage = self._faceMeans.regionImage(self.map.labelImage())
        else:
            sys.stderr.write("Unknown background mode %d!\n" % mode)
            return

        self.image = displayImage
        self._backgroundMode = mode
        normalize = self._normalizeStates[mode]
        if self.normalizeAction.isOn() != normalize:
            self.normalizeAction.setOn(normalize)
        else:
            self._setImage(self.image, normalize)

    def toggleNodeDisplay(self, onoff):
        if self.nodeOverlay.visible != onoff:
            self.nodeOverlay.visible = onoff
            self.viewer.update()

    def toggleEdgeDisplay(self, onoff):
        if self.edgeOverlay.visible != onoff:
            self.edgeOverlay.visible = onoff
            self.viewer.update()

    def toggleNormalize(self, normalize):
        self._normalizeStates[self._backgroundMode] = normalize
        self._setImage(self.image, normalize)

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
        roi &= Rect2D(self.map.imageSize())
        roiImage = self.map.labelImage().subImage(roi)
        if self._backgroundMode > 3:
            roiImage = self._faceMeans.regionImage(roiImage)
        # FIXME: use global normalization here
        self.viewer.replaceROI(roiImage.toPNM(vigra.BYTE),
                               qt.QPoint(*roi.upperLeft()))

    def _postMergeFacesHook(self, survivor):
        if self._backgroundMode < 3:
            return
        if survivor.label():
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
        return self.__base.showEvent(self, e)

    def hideEvent(self, e):
        self.detachHooks()
        return self.__base.hideEvent(self, e)

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

    def setTool(self, tool = None):
        """Deactivates and destroys old tool, activates/sets new one
        if tool != None.  `tool` can be either a tool object or

        1 for the MapSearcher tool
        2 for the ActivatePaintbrush
        3 for the IntelligentScissors
        4 for the SeedSelector"""

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
            if self.edgeOverlay.color == qt.Qt.red:
                self.edgeOverlay.color = qt.Qt.black
            self.nodeOverlay.visible = False
            self.tool = tools.IntelligentScissors(
                self.map, self.edgeOverlay, self)
            if self._faceMeans:
                tools.activeCostMeasure = \
                    statistics.HyperbolicInverse(self._faceMeans.faceMeanDiff)
        elif tool == 4:
            self.tool = tools.SeedSelector(map = self.map, parent = self)
        elif hasattr(tool, "disconnectViewer"):
            self.tool = tool
        elif tool != None:
            print "setTool: invalid argument, tool deactivated now."
            print "  give 1 for MapSearcher, 2 for ActivePaintbrush, 3 for IntelligentScissors, 4 for SeedSelector"
            tool = 0
        self._togglingGUI = True
        self.scissorsAction.setOn(tool == 3)
        self.paintbrushAction.setOn(tool == 2)
        self.navigateAction.setOn(tool == 1)
        self._togglingGUI = False

    def activateScissors(self, onoff = True):
        self.setTool(onoff and 3 or None)

    def activatePaintbrush(self, onoff = True):
        self.setTool(onoff and 2 or None)

    def activateNavigator(self, onoff = True):
        self.setTool(onoff and 1 or None)

    def cleanupMap(self):
        removeCruft(self.map, 7)

    def navigate(self, dart, center = True, costMeasure = None, new = False):
        if type(dart) == int:
            dart = self.map.dart(dart)
        elif hasattr(dart, "anchor"):
            dart = dart.anchor()
        elif hasattr(dart, "contours"):
            dart = list(dart.contours())
        if new or not self.dn:
            self.dn = DartNavigator(dart, costMeasure, self)
            self.connect(self.dn, qt.SIGNAL("destroyed(QObject*)"),
                         self._dartNavigatorDestroyed)
        else:
            if costMeasure:
                self.dn.costMeasure = costMeasure
            self.dn.setDart(dart)
        self.dn.show()
        if center:
            if isinstance(dart, (list, tuple)):
                dart = dart[0]
            self.viewer.center = dart[0]
            self.viewer.optimizeUpperLeft()

    def _enableImageActions(self):
        self.displayBIAction.setEnabled("bi" in self.images)
        self.displayOriginalAction.setEnabled("original" in self.images)
        self.displayColoredAction.setEnabled("colored" in self.images)
        self.displayLabelsAction.setEnabled(True)
        self.displayMeansAction.setEnabled(bool(self._faceMeans))

    def _setImage(self, image, normalize):
        self.image = image
        self._imageWindow.replaceImage(
            self.image, normalize and vigra.NBYTE or vigra.BYTE)

    def setImage(self, image, normalize = None, role = None):
        """Replace displayed background image.  You may pass role as
        one of ('original', 'colored', 'bi') to replace one of the
        predefined image slots (keyboard shortcuts 1-5)."""
        if role == None:
            if image.bands() == 3:
                self.images["original"] = vigra.transformImage(
#                    image, "\l x: RGB2Lab(x)")[0]
                    image, "\l x: norm(x)/%r" % math.sqrt(3))
                role = "colored"
            else:
                role = "original"
        self.images[role] = image
        self._enableImageActions()
        if normalize == None:
            normalize = self._normalizeStates[self._backgroundMode]
        if self._imageWindow and role == self.currentRole():
            self._setImage(image, normalize)

    def highlight(self, darts):
        """highlight(darts)
        Highlight the given darts (can be any iterable returning labels
        or Dart objects)."""
        self._dh.highlight(darts)

    def saveFig(self, basepath, roi = None, scale = None,
                bgFilename = None, faceMeans = False, similarity = None,
                skipBorder = False):
        """Saves an XFig file as <basepath>.fig (and the pixel background
        as <basepath>_bg.png if `bgFilename` is not given) and returns
        the FigExporter object.

        If `roi` is not None, it determines the ROI to be saved.
        `roi` can be a BoundingBox or Rect2D object, a string like
        "10,10-40,30" (see fig.parseGeometry) or a tuple usable as
        Rect2D constructor arguments.  You can also specify `roi` =
        True to export the visible image region.

        scale is the size of a pixel in fig units (450 = 1cm) and
        defaults to a size resulting in approx. 20cm width of the
        output file.

        If bgFilename is given, it is supposed to be a background
        filename for the picture box in the XFig file, or False if no
        background is desired."""

        def addMapOverlayWithoutBorder(*args, **kwargs):
            kwargs["skipBorder"] = True
            return addMapOverlay(*args, **kwargs)

        overlayHandler = \
            skipBorder and addMapOverlayWithoutBorder or addMapOverlay

        if bgFilename == None and faceMeans:
            bgFilename = False

        fe = figexport.exportImageWindow(
            self._imageWindow, basepath, roi = roi, scale = scale,
            bgFilename = bgFilename,
            overlayHandler = overlayHandler)
        
        if faceMeans in (None, True):
            faceMeans = self._faceMeans

        if faceMeans:
            fe.addMapFaces(
                self.map, faceMeans, similarity, depth = 900)
            fe.f.save()
        
        return fe

    def saveEPS(self, basepath, *args, **kwargs):
        """display.saveEPS(basepath, roi = None, scale = None)

        Saves an XFig file as <basepath>.fig (see saveFig()
        documentation for details) and calls fig2dev to create an
        additional <basepath>.eps.

        Returns the final filename (result of calling fig.File.fig2dev)."""

        fe = self.saveFig(basepath, *args, **kwargs)
        return fe.f.fig2dev(lang = "eps")

    def savePDF(self, basepath, *args, **kwargs):
        """display.savePDF(basepath, roi = None, scale = None)

        Saves an XFig file as <basepath>.fig (see saveFig()
        documentation for details) and calls fig2dev to create an
        additional <basepath>.pdf.

        Returns the final filename (result of calling fig.File.fig2dev)."""

        fe = self.saveFig(basepath, *args, **kwargs)
        return fe.f.fig2dev(lang = "pdf")

# --------------------------------------------------------------------
#                         dart navigation dialog
# --------------------------------------------------------------------

class DartNavigator(dartnavigator.DartNavigatorBase):
    __base = dartnavigator.DartNavigatorBase
    
    def __init__(self, dart, costMeasure, parent, name = None):
        self.__base.__init__(self, parent, name)
        self.dart = dart
        self.costMeasure = costMeasure
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
        self.__base.closeEvent(self, e)
        if e.isAccepted():
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
        dartDesc = "Dart %d, length %.1f, partial area %.1f, %d points" % (
            self.dart.label(), self.dart.edge().length(),
            self.dart.partialArea(), len(self.dart))
        if self.costMeasure:
            dartDesc += "\nassociated cost: %s" % self.costMeasure(self.dart)
        self.dartLabel.setText(dartDesc)
        for node, nodeLabel in ((self.dart.startNode(), self.startNodeLabel),
                                (self.dart.endNode(), self.endNodeLabel)):
            nodeLabel.setText(
                "Node %d (deg. %d)\nat %s" % (
                node.label(), node.degree(), node.position()))
        if self.dart.map().mapInitialized():
            leftFace = self.dart.leftFace()
            rightFace = self.dart.rightFace()
            self.faceLabel.setText(
                """Left: %s\nRight: %s""" % (str(leftFace)[8:-1], str(rightFace)[8:-1]))
        self.setCaption("DartNavigator(%d)" % (self.dart.label(), ))
        self.emit(qt.PYSIGNAL('updateDart'),(self.dart,))

# --------------------------------------------------------------------

import copy

class ROISelector(qt.QObject):
    def __init__(self, parent = None, name = None, imageSize = None,
                 roi = None, viewer = None, color = qt.Qt.yellow, width = 0,
                 alwaysVisible = True):
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
        self._viewer.installEventFilter(self)

        self.setVisible(alwaysVisible)

    def eventFilter(self, watched, e):
        if e.type() in (qt.QEvent.KeyPress, qt.QEvent.KeyRelease,
                        qt.QEvent.MouseButtonPress, qt.QEvent.MouseButtonRelease,
                        qt.QEvent.MouseButtonDblClick, qt.QEvent.MouseMove):
            self._keyState = e.stateAfter()
        return False

    def setVisible(self, onoff):
        """Sets flag whether the ROI should be always visible, or only
        during painting."""
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
            self._viewer.toWindowCoordinates(self.roi.left()-0.5, self.roi.top()-0.5),
            self._viewer.toWindowCoordinates(self.roi.right()-0.5, self.roi.bottom()-0.5))

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
        self._viewer.removeEventFilter(self)

    def setZoom(self, zoom):
        self.zoom = zoom

    def draw(self, p):
        if not self.roi:
            return
        p.setPen(qt.QPen(self.color, self.width))
        p.setBrush(qt.Qt.NoBrush)
        drawRect = self.windowRect()
        # painter is already set up with a shift:
        drawRect.moveBy(-self._viewer.upperLeft().x(),
                        -self._viewer.upperLeft().y())
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
