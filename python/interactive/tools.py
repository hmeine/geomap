"""`tools` - module with interactive GeoMap tools:

You will find three tool classes:
* MapSearcher
* ManualClassifier
* ActivePaintbrush
* IntelligentScissors

tools.activePathMeasure is used to steer the IntelligentScissors tool.
It can be assigned any object which returns a path cost when called
with two arguments:

* The first argument will be a LiveWire object
  representing the beginning of a path, and

* the second argument will be a Dart object (whose startNode is the
  liveWire's end node) which could potentially be added to the path.

``activePathMeasure(liveWire, newDart)`` should return the cost of
combined path (old path plus ``newDart``), and this must be larger
than the old cost ``liveWire.totalCost()``.

Cf. the (documented) class SimplePathCostMeasure which looks at the
path darts independently with a given measure, e.g. use::

  tools.activePathMeasure = SimplePathCostMeasure(faceMeanDiff)

to (locally, dart-wise) depend on the faceMeanDiff measure.

You might also be interested in the ESPathCostMeasure (using Euler
Spirals), or the `LiveWire` helper class (which is the type of the
object passed into the above-mentioned activePathMeasure's __call__
method)."""

_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import sys, qt, math
from maputils import mergeFacesByLabel, contourDarts, protectFace
from flag_constants import *
from vigrapyqt import EdgeOverlay, PointOverlay
from vigra import *

__all__ = ["MapSearcher", "ManualClassifier", "ActivePaintbrush",
           "IntelligentScissors",
           "LiveWire", "SimplePathCostMeasure", "ESPathCostMeasure"]

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
        nearestNode = None
        minDist = None
        for node in self._map.nodeIter():
            dx, dy = abs(node.position()[0]-x), abs(node.position()[1]-y)
            dist = dx*dx + dy*dy
            if minDist == None or minDist > dist:
                minDist = dist
                nearestNode = node
        #sys.stdout.write("Node %d is %.2f from %d/%d\n" % (nearestNode.label(), minDist, x, y))
        if nearestNode.degree() > 0:
            self.display.navigate(nearestNode.anchor())

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
    def __init__(self, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.seeds = []

        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        self.seeds.append((x, y))
        
    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)

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
        self._currentLabel = None
        self._painting = True
        self._changed = False
        self.mouseMoved(x, y)

    def mouseMoved(self, x, y):
        if not self._painting:
            return

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

    def mouseDoubleClicked(self, x, y):
        face = self._map.faceAt((x, y))
        protectFace(face, not face.flag(PROTECTED_FACE))

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
    end nodes by calling setEndNodeLabel()."""
    
    def __init__(self, map, pathCostMeasure, startNodeLabel):
        self._map = map
        self._pathCostMeasure = pathCostMeasure
        self._startNodeLabel = startNodeLabel
        self._endNodeLabel = startNodeLabel

        self._nodePaths = [None] * (self._map.maxNodeLabel() + 1)
        self._nodePaths[self._startNodeLabel] = (0.0, None)

        self._searchBorder = []
        self._expandNode(self._startNodeLabel)

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

    def _expandNode(self, nodeLabel):
        """Add all neighbors of the given node to the searchBorder."""

        prevPath = self._nodePaths[nodeLabel]
        if prevPath[1]:
            comingFrom = -prevPath[1]
            anchor = self._map.dart(comingFrom)
        else:
            comingFrom = None
            anchor = self._map.node(nodeLabel).anchor()

        for dart in anchor.sigmaOrbit():
            if dart.label() == comingFrom or dart.edge().flag(CURRENT_CONTOUR):
                continue
            heappush(self._searchBorder, (
                self._pathCostMeasure(self, dart), dart.label()))

    def setEndNodeLabel(self, nodeLabel):
        """Try to set the live wire's end node to the given one.
        Returns True iff successful, i.e. a path to that node is
        already known.  You can then call pathDarts() to query the
        darts belonging to that path or totalCost() to get the cost
        of that path.."""

        if self._nodePaths[nodeLabel]:
            self._endNodeLabel = nodeLabel
            return True

    def pathDarts(self, endNodeLabel = None):
        """liveWire.pathDarts()

        Generator function returning all darts along the current live
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
    def __init__(self, map, edgeColors, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self._map = map

        self._liveWire = None
        self._startNodeLabel = None
        self._allEdges = []
        self._expandTimer = qt.QTimer(self, "expandTimer")
        self._edgeColors = edgeColors
        self.connect(self._expandTimer, qt.SIGNAL("timeout()"),
                     self._expandBorder)

        # connect viewer
        self.viewer = parent.viewer
        self.connect(self.viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(self.viewer, qt.PYSIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(self.viewer, qt.PYSIGNAL("mouseReleased"),
                     self.mouseReleased)
        self.connect(self.viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                     self.mouseDoubleClicked)
        self.overlayIndex = self.viewer.addOverlay(
            PointOverlay([], qt.Qt.green, 1))

    def disconnectViewer(self):
        self.disconnect(self.viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(self.viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(self.viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)
        self.disconnect(self.viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                        self.mouseDoubleClicked)
        self.viewer.removeOverlay(self.overlayIndex)

    def startSearch(self):
        """Starts a search at self._startNodeLabel.  Initialized the
        search and starts a QTimer repeatedly calling
        _expandBorder()."""

        if not self._liveWire or \
               self._liveWire.startNodeLabel() != self._startNodeLabel:
            self._liveWire = LiveWire(
                self._map, activePathMeasure, self._startNodeLabel)

        self._expandTimer.start(0)

    def stopSearch(self):
        """Stops the current search (i.e., the QTimer)."""

        self._expandTimer.stop()

    def _expandBorder(self):
        if not self._liveWire.expandBorder():
            self.stopSearch()
            return

    def mousePressed(self, x, y, button):
        """With left mouse button, the live wire is started, with the
        middle mouse button it can be cancelled."""

        if button == qt.Qt.MidButton and self._liveWire:
            return self.stopCurrentContour()

        if button == qt.Qt.LeftButton and not self._liveWire:
            self.startSearch()

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
                    self.viewer.replaceOverlay(
                        EdgeOverlay(pathEdges, qt.Qt.yellow, 2),
                        self.overlayIndex)

    def mouseReleased(self, x, y, button):
        """With each left click, fix the current live wire and start a
        new one."""

        if button != qt.Qt.LeftButton or not self._liveWire:
            return

        self.stopSearch()

        pathEdges = [dart.edge()
                     for dart in self._liveWire.pathDarts()]
        for edge in pathEdges:
            edge.setFlag(SCISSOR_PROTECTION | CURRENT_CONTOUR)
            self._edgeColors[edge.label()] = qt.Qt.green
        self._allEdges.extend(pathEdges)

        self.viewer.replaceOverlay(
            EdgeOverlay([], qt.Qt.yellow, 2), self.overlayIndex)

        self._startNodeLabel = self._liveWire.endNodeLabel()

        self.startSearch()

    def mouseDoubleClicked(self, x, y, button):
        """With a double left click, the current live wire is fixed
        (by mouseReleased) and becomes inactive."""

        if button == qt.Qt.LeftButton:
            self.stopCurrentContour()

    def stopCurrentContour(self):
        self.stopSearch()
        #updateViewer(self.currentPathBounds)
        self._startNodeLabel = self._liveWire.endNodeLabel()
        self._liveWire = None
        for edge in self._allEdges:
            edge.setFlag(CURRENT_CONTOUR, False)

# --------------------------------------------------------------------

from hourglass import composeTangentLists
import statistics, saliency

class SimplePathCostMeasure(object):
    """SimplePathCostMeasure: Cost measure for a (e.g. livewire) path;
    this one does not actually evaluate any path continuity
    properties, but adds up the inverse of all darts' removal costs
    given a single dart cost measure.

    E.g. initialize with: SimplePathCostMeasure(minEdgeGradCost)"""

    def __init__(self, singleDartMeasure):
        """initialize with the given edge cost measure"""
        self.measure = singleDartMeasure

    def __call__(self, liveWire, newDart):
        return liveWire.totalCost(newDart.startNodeLabel()) \
               + 1.0 / (1e-4 + self.measure(newDart))

class ESPathCostMeasure(object):
    def __call__(self, liveWire, newDart):
        previousEndNodeLabel = newDart.startNodeLabel()
        
        allTangents = [statistics.dartTangents(newDart.clone().nextAlpha())]
        for dart in liveWire.pathDarts(previousEndNodeLabel):
            allTangents.append(statistics.dartTangents(dart))
        allTangents = composeTangentLists(allTangents)

        # this is really slow (at least with many support points).. :-(
        #allTangents = gaussianConvolveByArcLength(allTangents, 0.25)
        
        error = saliency.fitParabola(allTangents)

        # Don't know what exact formula to use here (tried some
        # without luck); the value should increase with every dart
        # added to the path (the accumulated previous costs are
        # available as liveWire.totalCost()), the livewire takes the
        # path with the *lowest* total cost.

        return max(liveWire.totalCost(previousEndNodeLabel), error)
        #return liveWire.totalCost(previousEndNodeLabel) + error
        # ... + math.exp(-error)

activePathMeasure = SimplePathCostMeasure(statistics.minEdgeGradCost)
