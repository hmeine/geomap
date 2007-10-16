import copy, qt
import mapdisplay, icons, maputils, statistics, flag_constants

def labelRoot(lut, label):
    result = lut[label]
    if lut[result] != result:
        result = labelRoot(lut, result)
        lut[label] = result
    return result

class PyramidContractionKernel(statistics.DetachableStatistics):
    __slots__ = ("ck", "_mergedFaceLabels")
    
    def __init__(self, map):
        statistics.DetachableStatistics.__init__(self, map)
        self.ck = [None] * map.maxFaceLabel()
        self._attachHooks()

    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeFacesCallbacks(
            self._preMergeFaces, self._postMergeFaces),
            )

    def _preMergeFaces(self, dart):
        self._mergedFaceLabels = (dart.leftFaceLabel(), dart.rightFaceLabel())
        return True

    def _postMergeFaces(self, survivor):
        labeled = (survivor.label(), self._map().faceCount)
        for faceLabel in self._mergedFaceLabels:
            self.ck[faceLabel] = labeled

    def kernelForFaceCount(self, faceCount):
        result = range(len(self.ck))
        for i, lfc in enumerate(self.ck):
            if lfc and lfc[1] >= faceCount:
                result[i] = lfc[0]
        for i, l in enumerate(result):
            result[i] = labelRoot(result, l)
        return result

    def applyToFaceCount(self, map, faceCount):
        ck = self.kernelForFaceCount(faceCount)
        maputils.applyFaceClassification(map, ck)

class PaintbrushStroke(object):
    __slots__ = ("workspace", "faces", "oldValues", "survivorLabel")

    def __init__(self, workspace, survivorLabel):
        self.workspace = workspace
        self.faces = list(workspace.map.faceLabelLUT().merged(survivorLabel))
        self.survivorLabel = survivorLabel
        self.redo()

    def redo(self):
        manualCK = self.workspace._manualCK
        self.oldValues = [manualCK[fl] for fl in self.faces]
        for fl in self.faces:
            manualCK[fl] = self.survivorLabel

    def undo(self):
        manualCK = self.workspace._manualCK
        for fl, ov in zip(self.faces, self.oldValues):
            manualCK[fl] = ov

class Workspace(mapdisplay.MapDisplay):
    """Workspace for region-based segmentation."""
    __base = mapdisplay.MapDisplay

    # actually, this has only documenting effect; since this inherits
    # PyQt widgets, any attribute may be used:
    __slots__ = ("_level0", "_manualCK", "_pyramidCK",
                 "_mapRestartAction", "_levelSlider",
                 "_history")
    
    def __init__(self, level0, originalImage):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0
        self._manualCK = range(level0.maxFaceLabel())
        self._pyramidCK = None
        self._manualCKHooks = ()
        self._history = []

        # needed for backpropagation of protection:
        if not hasattr(self.map, "mergedEdges"):
            self.map.mergedEdges = statistics.MergedEdges(self.map)

        reinitIcon = qt.QPixmap()
        reinitIcon.loadFromData(icons.reinitIconPNGData, "PNG")
        ra = qt.QAction(self, "mapRestartAction")
        ra.setIconSet(qt.QIconSet(reinitIcon))
        ra.setText("Re-start from level 0")
        ra.setMenuText("&Re-start from level 0")
        ra.setToolTip("Re-start with level 0 map")
        ra.addTo(self.Tools)
        self.connect(ra, qt.SIGNAL("activated()"), self.restart)
        self._mapRestartAction = ra

        #sws.imageFrame.addWidget(self.imageFrame)
        self._levelSlider = qt.QSlider(self._imageWindow, "_levelSlider")
        self._levelSlider.setOrientation(qt.QSlider.Horizontal)
        self._levelSlider.setSteps(-1, -10)
        self._levelSlider.setRange(0, level0.faceCount)
        self.connect(self._levelSlider, qt.SIGNAL("valueChanged(int)"),
                     self._levelSliderChanged)
        self._imageWindow._layout.addWidget(self._levelSlider)
        if self.isShown():
            self._levelSlider.show()

    @staticmethod
    def _detachMapStats(map):
        for a in map.__dict__:
            o = getattr(map, a)
            if hasattr(o, "detachHooks"):
                #print "detaching hooks of", o
                o.detachHooks()

    def setTool(self, tool = None):
        self.__base.setTool(self, tool)
        if not self.tool:
            return
        self.connect(self.tool, qt.PYSIGNAL("paintbrushFinished"),
                     self.paintbrushFinished)
        self.connect(self.tool, qt.PYSIGNAL("faceProtectionChanged"),
                     self.faceProtectionChanged)

    def _perform(self, action):
        # FIXME: redo support
        self._history.append(action)

    def undo(self):
        # FIXME: redo support
        self._history[-1].undo()
        del self._history[-1]
        self._pyramidCK = None # FIXME

    def paintbrushFinished(self, survivor):
        self._perform(PaintbrushStroke(self, survivor.label()))
        self._pyramidCK = None

    def faceProtectionChanged(self, face):
        #assert face.map() == self.map
        for dart in maputils.contourDarts(face):
            protected = dart.edge().flag(flag_constants.CONTOUR_PROTECTION)
            for edgeLabel in self.map.mergedEdges[dart.edgeLabel()]:
                self._level0.edge(edgeLabel).setFlag(
                    flag_constants.CONTOUR_PROTECTION, protected)

    def setMap(self, map):
        self._detachMapStats(self.map)
        self.__base.setMap(self, map)
        self._levelSlider.blockSignals(True)
        self._levelSlider.setValue(
            self._levelSlider.maxValue() - map.faceCount)
        self._levelSlider.blockSignals(False)

    def restart(self):
        """Restart with the manually created level."""
        self.setMap(self.manualBaseMap())

    def manualBaseMap(self):
        result = copy.deepcopy(self._level0)
        if not hasattr(result, "mergedEdges"):
            result.mergedEdges = statistics.MergedEdges(result)
        maputils.applyFaceClassification(result, self._manualCK)
        self._levelSlider.blockSignals(True)
        self._levelSlider.setRange(0, result.faceCount)
        self._levelSlider.blockSignals(False)
        return result

    def faceMeans(self, map):
        if not hasattr(map, "faceMeans"):
            img = self.images.get("colored", self.images["original"])
            map.faceMeans = statistics.FaceColorStatistics(map, img)
        return map.faceMeans

    def costMeasure(self, map):
        #return self.faceMeans().faceMeanDiff
        return self.faceMeans(map).faceHomogeneity
        #return self.faceMeans().faceTTest

    def _levelSliderChanged(self, levelIndex):
        self.displayLevel(levelIndex)

    def displayLevel(self, levelIndex):
        if not self._pyramidCK:
            self.automaticRegionMerging()

        faceCount = self._levelSlider.maxValue() - levelIndex

        if faceCount > self.map.faceCount:
            map = self.manualBaseMap()
            self._pyramidCK.applyToFaceCount(map, faceCount)
            self.setMap(map)
        elif faceCount < self.map.faceCount:
            self._pyramidCK.applyToFaceCount(self.map, faceCount)
            
            self._levelSlider.blockSignals(True)
            self._levelSlider.setValue(
                self._levelSlider.maxValue() - self.map.faceCount)
            self._levelSlider.blockSignals(False)

    def automaticRegionMerging(self):
        map = self.manualBaseMap()
        map.pyramidCK = PyramidContractionKernel(map)
        amr = maputils.AutomaticRegionMerger(map, self.costMeasure(map))
        amr.merge()
        self._detachMapStats(map)
        self._pyramidCK = map.pyramidCK
