from darthighlighter import DartHighlighter

# ui-generated base classes:
from dartnavigator_ui import Ui_DartNavigator
from PyQt4 import QtCore, QtGui

class DartNavigator(QtGui.QDialog):
    __base = QtGui.QDialog
    
    def __init__(self, dart, costMeasure, parent, name = None):
        self.__base.__init__(self, parent)
        if name:
            self.setObjectName(name)

        self.dart = dart
        self.costMeasure = costMeasure

        self.ui = Ui_DartNavigator()
        self.ui.setupUi(self)
        self.connect(self.ui.nextPhiButton, QtCore.SIGNAL("clicked()"),
                     self.nextPhi)
        self.connect(self.ui.prevPhiButton, QtCore.SIGNAL("clicked()"),
                     self.prevPhi)
        self.connect(self.ui.nextAlphaButton, QtCore.SIGNAL("clicked()"),
                     self.nextAlpha)
        self.connect(self.ui.nextSigmaButton, QtCore.SIGNAL("clicked()"),
                     self.nextSigma)
        self.connect(self.ui.prevSigmaButton, QtCore.SIGNAL("clicked()"),
                     self.prevSigma)

        self.connect(self.ui.continuousCheckBox, QtCore.SIGNAL("toggled(bool)"),
                     self.toggleContinuous)

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self.highlightNext)

        self._darthighlighter = DartHighlighter(parent.map, parent.viewer)
        self.updateLabel()

    def closeEvent(self, e):
        self._darthighlighter.highlight(None)
        self.__base.closeEvent(self, e)
        if e.isAccepted():
            self.deleteLater() # like QtCore.Qt.WDestructiveClose ;-)

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
        self.ui.nextPhiButton.setToggleButton(onoff)
        self.ui.prevPhiButton.setToggleButton(onoff)
        self.ui.nextAlphaButton.setToggleButton(onoff)
        self.ui.nextSigmaButton.setToggleButton(onoff)
        self.ui.prevSigmaButton.setToggleButton(onoff)

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
        self._darthighlighter.highlight(self.dart.label())
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
        self.emit(QtCore.SIGNAL('updateDart'),(self.dart,))
