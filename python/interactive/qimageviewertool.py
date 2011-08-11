from PyQt4 import QtCore

class QImageViewerTool(QtCore.QObject):
    """Base class for implementing tools that are interested in mouse
    button events of a QImageViewer.  Sub-classes will want to
    re-implement one or more of the following stub methods:

      def mousePressed(self, x, y, button):
          return False

      def mouseDoubleClicked(self, x, y, button):
          return False

      def mouseReleased(self, x, y, button):
          return False

      def mouseMoved(self, x, y):
          pass

    Coordinates will be given in the image coordinate system
    (cf. `connectViewer()` which offers an option for integer
    coordinates).  Button handlers *must* return a boolean in all
    cases; return True iff you handled the event (this will stop event
    propagation) - otherwise, return False, not None."""

    __slots__ = ("_viewer", "_integerCoordinates", "_currentEvent")

    def connectViewer(self, viewer, integerCoordinates = False):
        """Register with the given viewer (installs event filter).  If
        integerCoordinates is set, event coordinates will be given in
        integer pixel coordinates (default: sub-pixel coords)."""
        self._viewer = viewer
        self._integerCoordinates = integerCoordinates
        self._viewer.installEventFilter(self)
        self.connect(self._viewer, QtCore.SIGNAL("mouseOver(int,int)"),
                     self.mouseMoved)

    def disconnectViewer(self):
        self._viewer.removeEventFilter(self)
        self.disconnect(self._viewer, QtCore.SIGNAL("mouseOver(int,int)"),
                        self.mouseMoved)

    def activeModifiers(self):
        """Return QInputEvent.modifiers() of the currently handled event."""
        return self._currentEvent.modifiers()

    def eventFilter(self, watched, event):
        if event.type() in (QtCore.QEvent.MouseButtonPress,
                            QtCore.QEvent.MouseButtonDblClick,
                            QtCore.QEvent.MouseButtonRelease):
            if self._integerCoordinates:
                pixel = self._viewer.imageCoordinate(event.pos())
            else:
                pixel = self._viewer.imageCoordinateF(event.pos())

            self._currentEvent = event
            if event.type() == QtCore.QEvent.MouseButtonPress:
                return self.mousePressed(pixel.x(), pixel.y(), event.button())
            elif event.type() == QtCore.QEvent.MouseButtonDblClick:
                return self.mouseDoubleClicked(pixel.x(), pixel.y(), event.button())
            else:
                return self.mouseReleased(pixel.x(), pixel.y(), event.button())

        return False

    def mousePressed(self, x, y, button):
        return False

    def mouseDoubleClicked(self, x, y, button):
        return False

    def mouseReleased(self, x, y, button):
        return False

    def mouseMoved(self, x, y):
        pass
