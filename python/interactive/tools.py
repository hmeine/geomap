"""`tools` - module with interactive GeoMap tools:

You will find the following tool classes:
* MapSearcher
* ManualClassifier
* ActivePaintbrush
* IntelligentScissors
"""

_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import sys, qt, math, maputils
from maputils import mergeFacesByLabel, contourDarts
from flag_constants import *
from vigrapyqt import EdgeOverlay, PointOverlay
from vigra import *

__all__ = ["MapSearcher", "ManualClassifier", "ActivePaintbrush", "SeedSelector",
           "IntelligentScissors", "LiveWire"]

# --------------------------------------------------------------------
#                     interaction with image viewer
# --------------------------------------------------------------------

class MapSearcher(qt.QObject):
    def __init__(self, map, display, name = None):
        qt.QObject.__init__(self, display, name)
        self._map = map
        self.display = display
        self.connect(display.viewer, qt.PYSIGNAL("mousePressed"),
                     self.search)

    def search(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        nearestNode = self._map.nearestNode((x, y))
        #sys.stdout.write("Node %d is %.2f from %d/%d\n" % (nearestNode.label(), minDist, x, y))
        if nearestNode.hasMinDegree(1):
            self.display.navigate(nearestNode.anchor(), center = False)

    def disconnectViewer(self):
        self.disconnect(self.display.viewer, qt.PYSIGNAL("mousePressed"),
                        self.search)

class ManualClassifier(qt.QObject):
    classes = [None, False, True]

    def __init__(self, map, foreground, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self._map = map
        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.foreground = foreground
        self.manual = {}

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        face = self._map.faceAt((x, y))
        oldClass = self.foreground[face.label()]
        newClass = self.classes[(self.classes.index(oldClass) + 1) % 3]
        self.foreground[face.label()] = newClass
        self.manual[face.label()] = newClass
        print "manually changed face %d to %s" % (face.label(), newClass)
        self.emit(qt.PYSIGNAL("classChanged"), (face, ))

    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)

class SeedSelector(qt.QObject):
    def __init__(self, map = None, markFlags = SRG_SEED,
                 parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.seeds = []

        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)

        self.overlay = PointOverlay(self.seeds, qt.Qt.cyan, 2)
        viewer.addOverlay(self.overlay)

        self.map = map
        self.markFlags = markFlags

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        self.seeds.append((x, y))

        if self.map and self.markFlags:
            self.map.faceAt((x, y)).setFlag(self.markFlags)
        
        self.overlay.setPoints(self.seeds)
        viewer = self.parent().viewer
        viewer.update()
        
    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        viewer.removeOverlay(self.overlay)

class ActivePaintbrush(qt.QObject):
    def __init__(self, map, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self._map = map
        self._painting = False
        self._currentLabel = None
        self._changed = None
        s = map.imageSize()
        self._mapArea = s.width * s.height

        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(viewer, qt.PYSIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(viewer, qt.PYSIGNAL("mouseReleased"),
                     self.mouseReleased)
        self.connect(viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                     self.mouseDoubleClicked)

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        if not self._map.mapInitialized():
            sys.stderr.write("Paintbrush: Map not initialized. Unable to determine faces.\n")
            return
        self._currentLabel = None
        self._painting = True
        self._path = []
        self._changed = False
        self.mouseMoved(x, y)

    def mouseMoved(self, x, y):
        if not self._painting:
            return

        self._path.append((x, y))

        map = self._map
        otherLabel = map.faceAt((x, y)).label()
        if otherLabel == 0 and map.face(0).area() < -self._mapArea + 1:
            currentLabel = None
            return
        if self._currentLabel == None:
            self._currentLabel = otherLabel
        if self._currentLabel == otherLabel:
            return

        try:
            survivor = mergeFacesByLabel(map, self._currentLabel, otherLabel)
            if survivor:
                self._changed = True
                self._currentLabel = survivor.label()
            else:
                self._currentLabel = otherLabel
        except Exception, e:
            sys.stderr.write("Paintbrush: Merge operation failed. Cancelling paint mode.\n")
            self._painting = False
            raise

    def mouseReleased(self, x, y):
        self._painting = False
        self.emit(qt.PYSIGNAL("paintbrushFinished"), (
            self._map.face(self._currentLabel), ))

    def mouseDoubleClicked(self, x, y):
        face = self._map.faceAt((x, y))
        maputils.protectFace(face, not face.flag(PROTECTED_FACE))
        self.emit(qt.PYSIGNAL("faceProtectionChanged"), (face, ))

    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)
        self.disconnect(viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                        self.mouseDoubleClicked)

# --------------------------------------------------------------------

from heapq import * # requires Python 2.3+

class LiveWire(object):
    """The LiveWire class does not only represent a single live wire
    path, but also performs the complete path search in a dynamic
    programming fashion, i.e. finding the optimal paths to all
    reachable nodes.  You can then immediately switch between desired
    end nodes by calling setEndNodeLabel().

    Edges that are marked with the CURRENT_CONTOUR flag are avoided
    (use this to prevent going the same way back)."""

    def __init__(self, map, costMeasure, startNodeLabel):
        self._map = map
        self._costMeasure = costMeasure
        if hasattr(startNodeLabel, "label"):
            startNodeLabel = startNodeLabel.label()
        self._startNodeLabel = startNodeLabel
        self._endNodeLabel = startNodeLabel

        self._nodePaths = [None] * (self._map.maxNodeLabel() + 1)
        self._nodePaths[self._startNodeLabel] = (0.0, None)

        self._searchBorder = []

        # similar to _expandNode, but the latter starts from the end
        # of a dart, which we do not have yet here:
        for dart in map.node(startNodeLabel).anchor().sigmaOrbit():
            if dart.edge().flag(CURRENT_CONTOUR):
                continue
            heappush(self._searchBorder, (
                costMeasure(dart), dart.label()))

    def startNodeLabel(self):
        """liveWire.startNodeLabel() -> int

        Returns label of the Node at the beginning of the current live
        wire's path."""

        return self._startNodeLabel

    def endNodeLabel(self):
        """liveWire.endNodeLabel() -> int

        Returns label of the Node at the beginning of the current live
        wire's path."""

        return self._endNodeLabel

    def expandBorder(self):
        """liveWire.expandBorder()

        Performs a single step of the dynamic programming for finding
        all optimal paths from the start node.  Returns False iff the
        process finished.  Call this e.g. from an idle loop as long as
        it returns True.

        Internally, picks cheapest path from searchBorder, and if no
        path to its end node is known yet, stores it and calls
        _expandNode() to proceed with its neighbor nodes."""

        if not len(self._searchBorder):
            return False

        path = heappop(self._searchBorder)
        endNodeLabel = self._map.dart(path[1]).endNodeLabel()
        if not self._nodePaths[endNodeLabel]: # or self._nodePaths[endNodeLabel][0] > path[0]
            self._nodePaths[endNodeLabel] = path
            self._expandNode(endNodeLabel)

        return True

    def expandToNode(self, nodeLabel):
        """Expand (via expandBorder()) until a path to the indicated
        node is known."""

        while not self._nodePaths[nodeLabel]:
            if not self.expandBorder():
                return False
        return True

    def expandToCost(self, cost):
        """Expand (via expandBorder()) until all paths < cost are known."""

        while self._searchBorder and self._searchBorder[0][0] < cost:
            self.expandBorder()
        return bool(self._searchBorder)

    def expand(self):
        """Expand completely, until expandBorder() returns False."""

        while self.expandBorder():
            pass

    def _expandNode(self, nodeLabel):
        """Add all neighbors of the given node to the searchBorder."""

        prevPath = self._nodePaths[nodeLabel]
        sigmaOrbit = self._map.dart(-prevPath[1]).sigmaOrbit()
        sigmaOrbit.next() # skip prevPath[1] where we're coming from

        for dart in sigmaOrbit:
            if dart.edge().flag(CURRENT_CONTOUR):
                continue
            heappush(self._searchBorder, (
                prevPath[0] + self._costMeasure(dart), dart.label()))

    def setEndNodeLabel(self, nodeLabel):
        """Try to set the live wire's end node to the given one.
        Returns ``True`` iff successful, i.e. a path to that node is
        already known.  You can then call `pathDarts()` to query the
        darts belonging to that path or `totalCost()` to get the cost
        of that path.."""

        if self._nodePaths[nodeLabel]:
            self._endNodeLabel = nodeLabel
            return True

    def loopPath(self, nodeLabel):
        """Return additional path segment from nodeLabel to
        endNodeLabel.  If the optimal path from startNodeLabel to
        `nodeLabel` (i.e. pathDarts(nodeLabel)) passes endNodeLabel,
        return this (e.g. loop closing) path segment, else return
        None. """
        
        if self._nodePaths[nodeLabel]:
            result = []
            for dart in self.pathDarts(nodeLabel):
                result.append(dart)
                if dart.endNodeLabel() == self._endNodeLabel:
                    return result

    def pathDarts(self, endNodeLabel = None):
        """Generator function returning all darts along the current live
        wire path, ordered and pointing back from the current
        endNodeLabel() to startNodeLabel()."""

        nl = endNodeLabel
        if nl == None:
            nl = self._endNodeLabel
        while nl != self._startNodeLabel:
            dart = self._map.dart(-self._nodePaths[nl][1])
            yield dart
            nl = dart.endNodeLabel()

    def totalCost(self, endNodeLabel = None):
        if endNodeLabel == None:
            endNodeLabel = self._endNodeLabel
        return self._nodePaths[endNodeLabel][0]

class IntelligentScissors(qt.QObject):
    def __init__(self, map, mapEdges, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self._map = map

        self._liveWire = None
        self._loopNodeLabel = None
        self._startNodeLabel = None
        self._loop = None
        self._contour = [] # darts within current (multi-segment) contour
        self._prevContour = None # last finished _contour
        self._seeds = [] # all seeds of all contours (for debugging ATM)
        self._expandTimer = qt.QTimer(self, "expandTimer")
        self._mapEdges = mapEdges
        self.connect(self._expandTimer, qt.SIGNAL("timeout()"),
                     self._expandBorder)
        self.protect = True

        # connect viewer
        self.viewer = parent.viewer
        self.connect(self.viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(self.viewer, qt.PYSIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(self.viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                     self.mouseDoubleClicked)
        self.viewer.installEventFilter(self)
        self.overlayIndex = self.viewer.addOverlay(
            PointOverlay([], qt.Qt.green, 1))

    def disconnectViewer(self):
        self.disconnect(self.viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(self.viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(self.viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                        self.mouseDoubleClicked)
        self.viewer.removeOverlay(self.overlayIndex)
        self.viewer.removeEventFilter(self)

    def eventFilter(self, watched, e):
        if e.type() in (qt.QEvent.KeyPress, qt.QEvent.KeyRelease,
                        qt.QEvent.MouseButtonPress, qt.QEvent.MouseButtonRelease,
                        qt.QEvent.MouseButtonDblClick, qt.QEvent.MouseMove):
            self._keyState = e.stateAfter()
            self.protect = not self._keyState & qt.Qt.ControlButton
        return False

    def startContour(self):
        self._loopNodeLabel = self._startNodeLabel

    def startLiveWire(self):
        """Start a LiveWire at self._startNodeLabel.
        Starts a QTimer for repeated calling of _expandBorder()."""

        self._seeds.append(self._startNodeLabel)

        if not self._liveWire or \
               self._liveWire.startNodeLabel() != self._startNodeLabel:
            self._liveWire = LiveWire(
                self._map, activeCostMeasure, self._startNodeLabel)

        self._expandTimer.start(0)

    def stopLiveWire(self):
        """Stop the LiveWire (i.e., the QTimer)."""

        self._expandTimer.stop()

    def stopCurrentContour(self):
        self.stopLiveWire()

        #updateViewer(self.currentPathBounds)
        self._startNodeLabel = self._liveWire.endNodeLabel()
        self._loopNodeLabel = None # not really needed I think
        self._liveWire = None
        for dart in self._contour:
            dart.edge().setFlag(CURRENT_CONTOUR, False)
        self._prevContour = self._contour
        self._contour = []

    def _expandBorder(self):
        if not self._liveWire.expandBorder():
            self._expandTimer.stop()
            return

    def mousePressed(self, x, y, button):
        """With left mouse button, the live wire is started, with the
        middle mouse button it can be cancelled."""

        if button == qt.Qt.MidButton and self._liveWire:
            return self.stopCurrentContour()

        if button == qt.Qt.LeftButton:
            if not self._liveWire:
                self.startContour()
            else:
                self.stopLiveWire()
                self._protectPath(self._liveWire.pathDarts())
                self._startNodeLabel = self._liveWire.endNodeLabel()
        
        self.viewer.replaceOverlay(
            EdgeOverlay([], qt.Qt.yellow, 2), self.overlayIndex)

        self.startLiveWire()

    def mouseMoved(self, x, y):
        """It the live wire is active, it's end node is set to the
        nearest node and the display (overlay) is updated. Else, the
        nearest node is chosen as start node and highlighted."""

        node = self._map.nearestNode((x, y))
        if not self._liveWire:
            if node.label() != self._startNodeLabel:
                self._startNodeLabel = node.label()
                self.viewer.replaceOverlay(
                    PointOverlay([node.position()], qt.Qt.green, 1),
                    self.overlayIndex)
        else:
            if node.label() != self._liveWire.endNodeLabel():
                if self._liveWire.setEndNodeLabel(node.label()): # TODO: else... (delayed)
                    pathEdges = [dart.edge()
                                 for dart in self._liveWire.pathDarts()]
                    self._loop = self._liveWire.loopPath(self._loopNodeLabel)
                    if self._loop:
                        pathEdges.extend([dart.edge() for dart in self._loop])
                    self.viewer.replaceOverlay(
                        EdgeOverlay(pathEdges, qt.Qt.yellow, 2),
                        self.overlayIndex)

    def _protectPath(self, darts):
        for dart in darts:
            edge = dart.edge()
            edge.setFlag(SCISSOR_PROTECTION | CURRENT_CONTOUR, self.protect)
            self._mapEdges._updateEdgeROI(edge)
            self._contour.append(dart)

    def mouseDoubleClicked(self, x, y, button):
        """With a double left click, the current live wire is fixed
        (by mouseReleased) and becomes inactive."""

        if button == qt.Qt.LeftButton:
            if self._loop:
                self._protectPath(self._loop)
                self._loop = None
                self._seeds.append(self._loopNodeLabel)

            self.stopCurrentContour()

# --------------------------------------------------------------------

# FIXME: This shouldn't be set here:
activeCostMeasure = None
