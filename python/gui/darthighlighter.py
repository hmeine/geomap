import sys
from weakref import ref
from PyQt4 import QtCore, QtGui
import VigraQt

class DartHighlighter(VigraQt.Overlay):
    """The DartHighlighter class is attached to a viewer and a Map and
    is able to highlight any set of darts at a time."""
    
    def __init__(self, map, viewer):
        VigraQt.Overlay.__init__(self, viewer)
        self._map = ref(map)
        self._viewer = viewer
        self._darts = []
        self._color = QtCore.Qt.yellow
        viewer.addOverlay(self)

    def setMap(self, map):
        self._map = ref(map)

    def highlight(self, darts, color = QtCore.Qt.yellow):
        """highlight(darts)
        Highlight the given darts (can be any iterable returning labels
        or Dart objects)."""

        self._darts = []
        if darts:
            # hasattr(darts, "__iter__"): would be better, but is True for Edges
            if not isinstance(darts, (list, tuple)):
                darts = [darts]
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
                elif hasattr(dart, "contour"): # Faces
                    dart = dart.contour()
                self._darts.append(dart)

        self._viewer.update()

    def draw(self, p, rect):
        if not self._darts:
            return

        p.setPen(self._color)

        p.setBrush(QtGui.QBrush(self._color))
        w = 2 / self._viewer.zoomFactor()
        for dart in self._darts:
            p.drawEllipse(QtCore.QPointF(*dart[0]), w, w)
            poly = QtGui.QPolygonF(len(dart))
            for i, (x, y) in enumerate(dart.edge()):
                poly[i] = QtCore.QPointF(x, y)
            p.drawPolyline(poly)
