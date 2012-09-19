import sys
from weakref import ref
from PyQt4 import QtCore
import vigra.pyqt

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

    def highlight(self, darts, color = QtCore.Qt.yellow):
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

        self.eo = vigrapyqt4.EdgeOverlay([dart.edge() for dart in darts], color)
        self.eo.width = 2
        self.no = vigrapyqt4.PointOverlay(
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
