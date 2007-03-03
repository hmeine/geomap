"""`tools` - module with interactive GeoMap tools:

You will find three tool classes:
* MapSearcher
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
from maputils import mergeFacesByLabel
from vigrapyqt import EdgeOverlay, PointOverlay
from vigra import *

__all__ = ["MapSearcher", "ActivePaintbrush", "IntelligentScissors",
           "LiveWire", "SimplePathCostMeasure", "ESPathCostMeasure"]

# --------------------------------------------------------------------
#                     interaction with image viewer
# --------------------------------------------------------------------

class MapSearcher(qt.QObject):
    def __init__(self, map, display, name = None):
        qt.QObject.__init__(self, display, name)
        self.map = map
        self.display = display
        self.connect(display.viewer, qt.PYSIGNAL("mousePressed"),
                     self.search)

    def search(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        nearestNode = None
        minDist = None
        for node in self.map.nodeIter():
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
        self.map = map
        viewer = parent.viewer
        self.connect(viewer, qt.PYSIGNAL("mousePressed"),
                     self.mousePressed)
        self.foreground = foreground
        self.manual = {}

    def mousePressed(self, x, y, button):
        if button != qt.Qt.LeftButton:
            return
        face = self.map.faceAt(Vector2(x, y))
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

class ActivePaintbrush(qt.QObject):
    def __init__(self, map, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.map = map
        self.painting = False
        self.currentLabel = None
        self.changed = None
        s = map.imageSize()
        self.mapArea = s.width * s.height

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
        self.currentLabel = None
        self.painting = True
        self.changed = False
        self.mouseMoved(x, y)

    def mouseMoved(self, x, y):
        if not self.painting: return # comment out to get mouse-over face output

        map = self.map
        otherLabel = map.faceAt(Vector2(x, y)).label()
        if otherLabel == 0 and map.face(0).area() < -self.mapArea + 1:
            currentLabel = None
            return
        if self.currentLabel == None:
            self.currentLabel = otherLabel
        if self.currentLabel == otherLabel:
            return

        if not self.painting:
            print "mouse moved into", map.face(otherLabel)
            self.currentLabel = otherLabel
            return

        try:
            survivor = mergeFacesByLabel(map, self.currentLabel, otherLabel)
            if survivor:
                self.changed = True
                self.currentLabel = survivor.label()
            else:
                self.currentLabel = otherLabel
        except Exception, e:
            sys.stderr.write("Paintbrush: Merge operation failed. Cancelling paint mode.\n")
            self.painting = False
            raise

    def mouseReleased(self, x, y):
        self.painting = False
#         if self.changed:
#             self.parent().nodeOverlay._calculatePoints()
#             self.parent().viewer.update()

    def disconnectViewer(self):
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)

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
            if dart.label() == comingFrom:
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
    def __init__(self, map, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.map = map

        self.liveWire = None
        self.startNodeLabel = None
        self.expandTimer = qt.QTimer(self, "expandTimer")
        self.connect(self.expandTimer, qt.SIGNAL("timeout()"),
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
        viewer = self.parent().viewer
        self.disconnect(viewer, qt.PYSIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(viewer, qt.PYSIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(viewer, qt.PYSIGNAL("mouseReleased"),
                        self.mouseReleased)
        self.disconnect(self.viewer, qt.PYSIGNAL("mouseDoubleClicked"),
                        self.mouseDoubleClicked)
        self.viewer.removeOverlay(self.overlayIndex)

    def startSearch(self):
        """Starts a search at self.startNodeLabel.  Initialized the
        search and starts a QTimer repeatedly calling
        _expandBorder()."""

        if not self.liveWire or \
               self.liveWire.startNodeLabel() != self.startNodeLabel:
            self.liveWire = LiveWire(
                self.map, activePathMeasure, self.startNodeLabel)

        self.expandTimer.start(0)

    def stopSearch(self):
        """Stops the current search (i.e., the QTimer)."""

        self.expandTimer.stop()

    def _expandBorder(self):
        if not self.liveWire.expandBorder():
            self.stopSearch()
            return

    def mousePressed(self, x, y, button):
        """With left mouse button, the live wire is started, with the
        middle mouse button it can be cancelled."""

        if button == qt.Qt.MidButton and self.liveWire:
            self.stopSearch()
            #updateViewer(self.currentPathBounds)
            self.startNodeLabel = self.liveWire.endNodeLabel()
            self.liveWire = None
            return

        if button == qt.Qt.LeftButton and not self.liveWire:
            self.startSearch()

    def mouseMoved(self, x, y):
        """It the live wire is active, it's end node is set to the
        nearest node and the display (overlay) is updated. Else, the
        nearest node is chosen as start node and highlighted."""

        p = Vector2(x, y)
        node = self.map.nearestNode(p)
        if not self.liveWire:
            if node.label() != self.startNodeLabel:
                self.startNodeLabel = node.label()
                self.viewer.replaceOverlay(
                    PointOverlay([node.position()], qt.Qt.green, 1),
                    self.overlayIndex)
        else:
            if node.label() != self.liveWire.endNodeLabel():
                if self.liveWire.setEndNodeLabel(node.label()): # TODO: else... (delayed)
                    pathEdges = [dart.edge() for dart in self.liveWire.pathDarts()]

                    self.viewer.replaceOverlay(
                        EdgeOverlay(pathEdges, qt.Qt.yellow), self.overlayIndex)

    def mouseReleased(self, x, y, button):
        """With each left click, fix the current live wire and start a
        new one."""

        if button != qt.Qt.LeftButton or not self.liveWire:
            return

        self.stopSearch()

        for dart in self.liveWire.pathDarts():
            dart.edge().color = qt.Qt.green
        self.viewer.update()
        self.viewer.replaceOverlay(
            EdgeOverlay([], qt.Qt.yellow), self.overlayIndex)

        self.startNodeLabel = self.liveWire.endNodeLabel()

        if button == qt.Qt.LeftButton:
            self.startSearch()

    def mouseDoubleClicked(self, x, y, button):
        """With a double left click, the current live wire is fixed
        (by mouseReleased) and becomes inactive."""

        if button == qt.Qt.LeftButton:
            self.stopSearch()
            self.startNodeLabel = self.liveWire.endNodeLabel()
            self.liveWire = None

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
