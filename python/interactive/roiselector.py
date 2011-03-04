import copy
from PyQt4 import QtCore, QtGui
from vigra import Rect2D
from geomap import intPos

class ROISelector(QtCore.QObject):
    def __init__(self, parent = None, imageSize = None,
                 roi = None, viewer = None, color = QtCore.Qt.yellow, width = 0,
                 alwaysVisible = True):
        QtCore.QObject.__init__(self, parent)
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

        self.connect(self._viewer, QtCore.SIGNAL("mousePressed"),
                     self.mousePressed)
        self.connect(self._viewer, QtCore.SIGNAL("mousePosition"),
                     self.mouseMoved)
        self.connect(self._viewer, QtCore.SIGNAL("mouseReleased"),
                     self.mouseReleased)
        self._viewer.installEventFilter(self)

        self.setVisible(alwaysVisible)

    def eventFilter(self, watched, e):
        if e.type() in (QtCore.QEvent.KeyPress, QtCore.QEvent.KeyRelease,
                        QtCore.QEvent.MouseButtonPress, QtCore.QEvent.MouseButtonRelease,
                        QtCore.QEvent.MouseButtonDblClick, QtCore.QEvent.MouseMove):
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
            self.emit(QtCore.SIGNAL("roiChanged"), (roi, ))

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
        if self._painting and button == QtCore.Qt.RightButton:
            self._stopPainting()
            self.setROI(self._oldROI)
            return
        if button != QtCore.Qt.LeftButton:
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
            return QtCore.QRect()
        return QtCore.QRect(
            self._viewer.toWindowCoordinates(self.roi.left()-0.5, self.roi.top()-0.5),
            self._viewer.toWindowCoordinates(self.roi.right()-0.5, self.roi.bottom()-0.5))

    def mouseReleased(self, x, y, button):
        if self._painting and button == QtCore.Qt.LeftButton:
            self._stopPainting()
            if self.roi is not None:
                self.setROI(self.roi & self._validRect)
                self.emit(QtCore.SIGNAL("roiSelected"), (self.roi, ))

    def disconnectViewer(self):
        self.disconnect(self._viewer, QtCore.SIGNAL("mousePressed"),
                        self.mousePressed)
        self.disconnect(self._viewer, QtCore.SIGNAL("mousePosition"),
                        self.mouseMoved)
        self.disconnect(self._viewer, QtCore.SIGNAL("mouseReleased"),
                        self.mouseReleased)
        if self._alwaysVisible:
            self._viewer.removeOverlay(self)
        self._viewer.removeEventFilter(self)

    def setZoom(self, zoom):
        self.zoom = zoom

    def draw(self, p):
        if not self.roi:
            return
        p.setPen(QtGui.QPen(self.color, self.width))
        p.setBrush(QtCore.Qt.NoBrush)
        drawRect = self.windowRect()
        # painter is already set up with a shift:
        drawRect.translate(-self._viewer.upperLeft().x(),
                           -self._viewer.upperLeft().y())
        p.drawRect(drawRect)

# def queryROI(imageWindow):
#     pd = QtGui.QProgressDialog(imageWindow)
#     pd.setLabelText("Please mark a ROI with drag & drop")
#     rs = ROISelector(imageWindow)
#     QtCore.QObject.connect(rs, QtCore.SIGNAL("roiSelected"),
#                        pd.accept)
#     if pd.exec_loop() == QtGui.QDialog.Accepted:
#         return rs.roi
#     return
