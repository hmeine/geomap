import copy, qt
import mapdisplay, icons, maputils, statistics

class Workspace(mapdisplay.MapDisplay):
    """Workspace for region-based segmentation."""
    __base = mapdisplay.MapDisplay

    __slots__ = ("_level0", "_mapRestartAction",
                 "_history", "_levelSlider", "_levelApproachingTimer",
                 "_currentLevelIndex", "_displayLevelIndex")
    
    def __init__(self, level0, originalImage):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0
        self._displayLevelIndex = 0
        self._currentLevelIndex = 0

        reinitIcon = qt.QPixmap()
        reinitIcon.loadFromData(icons.reinitIconPNGData, "PNG")
        ra = qt.QAction(self, "mapRestartAction")
        ra.setIconSet(qt.QIconSet(reinitIcon))
        ra.setText("Re-start from level 0")
        ra.setMenuText("&Re-start from level 0")
        ra.setToolTip("Re-start with level 0 map")
        ra.addTo(self.Tools)
        self.connect(ra, qt.SIGNAL("activated()"), self.restart)
        self.mapRestartAction = ra

        #sws.imageFrame.addWidget(self.imageFrame)
        self._levelSlider = qt.QSlider(self._imageWindow, "_levelSlider")
        self._levelSlider.setOrientation(qt.QSlider.Horizontal)
        self._levelSlider.setEnabled(False)
        self.connect(self._levelSlider, qt.SIGNAL("valueChanged(int)"),
                     self._levelSliderChanged)
        self._imageWindow._layout.addWidget(self._levelSlider)
        if self.isShown():
            self._levelSlider.show()

        self._levelApproachingTimer = qt.QTimer(self)
        self.connect(self._levelApproachingTimer, qt.SIGNAL("timeout()"),
                     self._approachLevel)

    def _detachMapStats(self):
        for a in self.map.__dict__:
            o = getattr(self.map, a)
            if hasattr(o, "detachHooks"):
                print "detaching hooks of", o
                o.detachHooks()

    def restart(self):
        """Restart with level 0/1."""

        self._detachMapStats()
        self.setMap(copy.deepcopy(self._level0))
        self._currentLevelIndex = 0

    def faceMeans(self):
        if not self._faceMeans:
            img = self.images.get("colored", self.images["original"])
            self.setFaceMeans(statistics.FaceColorStatistics(self.map, img))
        return self._faceMeans

    def costMeasure(self):
        #return self.faceMeans().faceMeanDiff
        return self.faceMeans().faceHomogeneity
        #return self.faceMeans().faceTTest

    def _levelSliderChanged(self, level):
        self.displayLevel(level)

    def displayLevel(self, levelIndex):
        if self._displayLevelIndex != levelIndex:
            self.viewer.setUpdatesEnabled(False)
            self.viewer.setEnabled(False) # prevent manual actions
            self._displayLevelIndex = levelIndex
            self._levelApproachingTimer.start(0)

    def _stopApproaching(self):
        self._levelApproachingTimer.stop()
        self.viewer.setUpdatesEnabled(True)
        self.viewer.setEnabled(True)
        self.viewer.update()
        if self._levelSlider.value() != self._currentLevelIndex:
            self._levelSlider.setValue(self._currentLevelIndex)

    def _approachLevel(self):
        if self._currentLevelIndex > self._displayLevelIndex:
            self.restart()
        elif self._currentLevelIndex < self._displayLevelIndex:
            targetIndex = min(self._displayLevelIndex,
                              self._currentLevelIndex + 25)
            try:
                self._history[self._currentLevelIndex:targetIndex].replay(
                    self.map, verbose = False)
                self._currentLevelIndex = targetIndex # FIXME: inc. via callback
            except RuntimeError, e:
                print e
                self._stopApproaching()
        else:
            self._stopApproaching()

    def automaticRegionMerging(self):
        self.viewer.setUpdatesEnabled(False)

        # FIXME: history is not region-based
        self._history = maputils.LiveHistory(self.map)
        amr = maputils.AutomaticRegionMerger(self.map, self.costMeasure())
        amr.merge()
        self._currentLevelIndex = len(self._history)
        self._levelSlider.setRange(0, len(self._history))
        self._levelSlider.setEnabled(True)

        self.viewer.setUpdatesEnabled(True)
        self.viewer.update()
