import copy, qt
import mapdisplay, icons, maputils, statistics

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

class Workspace(mapdisplay.MapDisplay):
    """Workspace for region-based segmentation."""
    __base = mapdisplay.MapDisplay

    # actually, this has only documenting effect; since this inherits
    # PyQt widgets, any attribute may be used:
    __slots__ = ("_level0", "_manualCK", "_pyramidCK",
                 "_mapRestartAction", "_levelSlider")
    
    def __init__(self, level0, originalImage):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0
        self._manualCK = range(level0.maxFaceLabel())
        self._pyramidCK = None
        self._manualCKHooks = ()

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
                print "detaching hooks of", o
                o.detachHooks()

    def setTool(self, tool = None):
        self.__base.setTool(self, tool)
        if not self.tool:
            return
        self.connect(self.tool, qt.PYSIGNAL("paintbrushFinished"),
                     self.paintbrushFinished)
        self.connect(self.tool, qt.PYSIGNAL("protectionChanged"),
                     self.protectionChanged)

    def paintbrushFinished(self, survivor):
        sl = survivor.label()
        for mergedLabel in self.map.faceLabelLUT().merged(sl):
            self._manualCK[mergedLabel] = sl

    def protectionChanged(self, face):
        pass

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
        maputils.applyFaceClassification(result, self._manualCK)
        self._levelSlider.blockSignals(True)
        self._levelSlider.setRange(0, result.faceCount)
        self._levelSlider.blockSignals(False)
        return result

    def faceMeans(self, map):
        if not hasattr(map, "faceMeans"):
            img = self.images.get("colored", self.images["original"])
            map.faceMeans = statistics.FaceColorStatistics(map, img)
            print map.faceMeans.checkConsistency()
        return map.faceMeans

    def costMeasure(self, map):
        #return self.faceMeans().faceMeanDiff
        return self.faceMeans(map).faceHomogeneity
        #return self.faceMeans().faceTTest

    def _levelSliderChanged(self, level):
        self.displayLevel(self._levelSlider.maxValue() - level)

    def displayLevel(self, faceCount):
        if not self._pyramidCK:
            self.automaticRegionMerging()

        ck = self._pyramidCK.kernelForFaceCount(faceCount)
        if faceCount > self.map.faceCount:
            map = self.manualBaseMap()
            maputils.applyFaceClassification(map, ck)
            self.setMap(map)
        else:
            maputils.applyFaceClassification(self.map, ck)
            
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
