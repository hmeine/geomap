_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

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
        sys.stdout.write("Node %d is %.2f from %d/%d\n" % (nearestNode._label, minDist, x, y))
        if nearestNode.degree() > 0:
            self.display.navigate(nearestNode.anchor())

    def disconnectViewer(self):
        self.disconnect(self.display.viewer, qt.PYSIGNAL("mousePressed"),
                        self.search)

class ActivePaintbrush(qt.QObject):
    def __init__(self, map, parent = None, name = None):
        qt.QObject.__init__(self, parent, name)
        self.map = map
        self.painting = False
        self.currentLabel = None
        self.changed = None
        s = map.labelImage.size()
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
        if x < 0 or x >= map.labelImage.width() or y < 0 or y >= map.labelImage.height():
            return
        otherLabel = int(map.labelImage[(x, y)])
        if otherLabel < 0:
            return
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
                self.currentLabel = survivor._label
            else:
                self.currentLabel = otherLabel
        except CancelOperation:
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
    reachable nodes."""
    
    def __init__(self, map, measure, startNodeLabel):
        self.map = map
        self.measure = measure
        self.startNodeLabel = startNodeLabel
        self.endNodeLabel = startNodeLabel

        self.nodePaths = [None] * (self.map.maxNodeLabel() + 1)
        self.nodePaths[self.startNodeLabel] = (0.0, None)

        self.searchBorder = []
        self.expandNode(self.startNodeLabel)

    def expandBorder(self):
        """Pick cheapest path from searchBorder, and if no path to its
        end node is known yet, store it and call expandNode() to
        proceed with its neighbor nodes."""

        if not len(self.searchBorder):
            return False

        path = heappop(self.searchBorder)
        endNodeLabel = self.map.dart(path[1]).endNodeLabel()
        if not self.nodePaths[endNodeLabel]: # or self.nodePaths[endNodeLabel][0] > path[0]
            self.nodePaths[endNodeLabel] = path
            self.expandNode(endNodeLabel)

        return True

    def expandNode(self, nodeLabel):
        """Add all neighbors of the given node to the searchBorder."""

        prevPath = self.nodePaths[nodeLabel]
        if prevPath[1]:
            comingFrom = -prevPath[1]
            anchor = self.map.dart(comingFrom)
        else:
            comingFrom = None
            anchor = self.map.node(nodeLabel).anchor()

        for dart in anchor.sigmaOrbit():
            if dart.label() == comingFrom:
                continue
            heappush(self.searchBorder, (self.measure(self, dart), dart.label()))

    def setEnd(self, nodeLabel):
        """Try to set the live wire's end node to the given one.
        Returns True iff successful, i.e. a path to that node is
        already known.  You can then call pathDarts() to query the
        darts belonging to that path."""

        if self.nodePaths[nodeLabel]:
            self.endNodeLabel = nodeLabel
            return True

    def pathDarts(self, endNodeLabel = None):
        """Generator function returning all darts along the current
        live wire path (ordered and pointing from end node to start
        node)."""

        nl = endNodeLabel
        if nl == None:
            nl = self.endNodeLabel
        while nl != self.startNodeLabel:
            dart = self.map.dart(-self.nodePaths[nl][1])
            yield dart
            nl = dart.endNodeLabel()

    def totalCost(self, endNodeLabel = None):
        if endNodeLabel == None:
            endNodeLabel = self.endNodeLabel
        return self.nodePaths[endNodeLabel][0]

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
        
        allTangents = [dartTangents(newDart.clone().nextAlpha())]
        for dart in liveWire.pathDarts(previousEndNodeLabel):
            allTangents.append(dartTangents(dart))
        allTangents = composeTangentLists(allTangents)

        # this is really slow (at least with many support points).. :-(
        #allTangents = gaussianConvolveByArcLength(allTangents, 0.25)
        
        error = fitParabola(allTangents)

        # Don't know what exact formula to use here (tried some
        # without luck); the value should increase with every dart
        # added to the path (the accumulated previous costs are
        # available as liveWire.totalCost()), the livewire takes the
        # path with the *lowest* total cost.

        return max(liveWire.totalCost(previousEndNodeLabel), error)
        #return liveWire.totalCost(previousEndNodeLabel) + error
        # ... + math.exp(-error)

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

        if not self.liveWire or self.liveWire.startNodeLabel != self.startNodeLabel:
            self.liveWire = LiveWire(self.map, activePathMeasure, self.startNodeLabel)

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
            self.startNodeLabel = self.liveWire.endNodeLabel
            self.liveWire = None
            return

        if button == qt.Qt.LeftButton and not self.liveWire:
            self.startSearch()

    def mouseMoved(self, x, y):
        """It the live wire is active, it's end node is set to the
        nearest node and the display (overlay) is updated. Else, the
        nearest node is chosen as start node and highlighted."""

        p = Vector2(x, y)
        node = self.map.node(self.map.nodeMap(p))
        if not self.liveWire:
            if node.label() != self.startNodeLabel:
                self.startNodeLabel = node.label()
                self.viewer.replaceOverlay(
                    PointOverlay([node._position], qt.Qt.green, 1),
                    self.overlayIndex)
        else:
            if node.label() != self.liveWire.endNodeLabel:
                if self.liveWire.setEnd(node.label()): # TODO: else... (delayed)
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

        self.startNodeLabel = self.liveWire.endNodeLabel

        if button == qt.Qt.LeftButton:
            self.startSearch()

    def mouseDoubleClicked(self, x, y, button):
        """With a double left click, the current live wire is fixed
        (by mouseReleased) and becomes inactive."""

        if button == qt.Qt.LeftButton:
            self.stopSearch()
            self.startNodeLabel = self.liveWire.endNodeLabel
            self.liveWire = None

