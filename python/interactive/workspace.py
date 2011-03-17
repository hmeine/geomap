#!/usr/bin/env python
import copy, sys
from PyQt4 import QtCore, QtGui
import vigra, geomap
import mapdisplay, icons, maputils, statistics, flag_constants, tools
import progress

# FIXME: how to adjust the slider range when switching between slider modes?

# FIXME: let recomputeAutomaticLevels() be a no-op if current level is bottom
# (maybe introduce invalidateAuto... and/or change displayLevel(), too?!)

# FIXME: complextest.py works, but slider-to-the-left acts very strangely!

# TODO: add API for adding cost measures (obviously needs factories)

def labelRoot(lut, label):
    result = lut[label]
    if lut[result] != result:
        result = labelRoot(lut, result)
        lut[label] = result
    return result

class PyramidContractionKernel(statistics.DetachableStatistics):
    """Represents a merge tree annotated with the number of remaining
    faces.  Thus, it can be used to compute a contraction kernel that
    encodes all merges up to a specific pyramid level / face count."""
    
    __base = statistics.DetachableStatistics
    __slots__ = ("ck", "_mergedFaceLabels")
    
    def __init__(self, map):
        self.__base.__init__(self, map)
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

# --------------------------------------------------------------------
#                        Undo-able Action Classes
# --------------------------------------------------------------------

class PaintbrushStroke(object):
    __slots__ = ("workspace", "survivorLabel",
                 "faces", "mergeCount", "oldValues")

    def __init__(self, workspace, survivorLabel):
        self.workspace = workspace
        self.survivorLabel = survivorLabel
        self.faces = list(workspace.map.faceLabelLUT().merged(survivorLabel))
        self.mergeCount = len(dict.fromkeys(
            [workspace._manualCK[faceLabel] for faceLabel in self.faces])) - 1
        if self.mergeCount:
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

        # TODO: add number of undone merges (faceCountOffset)
        self.workspace.recomputeAutomaticLevels()

    def __str__(self):
        return "Paintbrush stroke"

class FaceProtection(object):
    __slots__ = ("workspace", "faceLabel", "protected")

    def __init__(self, workspace, faceLabel):
        self.workspace = workspace
        self.faceLabel = faceLabel
        self.protected = workspace.map.face(faceLabel).flag(
            flag_constants.PROTECTED_FACE)
        self.redo()

    def protect(self, protect):
        if protect:
            self.workspace._protectedFaces.append(self.faceLabel)
        else:
            self.workspace._protectedFaces.remove(self.faceLabel)
        map = self.workspace.map
        face = map.face(map.faceLabelLUT()[self.faceLabel])
        maputils.protectFace(face, protect)

    def redo(self):
        self.protect(self.protected)

    def undo(self):
        self.protect(not self.protected)

    def __str__(self):
        if self.protected:
            return "Face protected"
        else:
            return "Face protection removed"

class ScissorsProtection(object):
    __slots__ = ("workspace", "contour")

    def __init__(self, workspace, contour):
        self.workspace = workspace
        self.contour = []
        #self.closed # see TODO below
        
        level0 = workspace._level0
        mergedEdges = workspace.map.mergedEdges
        for dart in contour:
            for edgeLabel in mergedEdges[dart.edgeLabel()]:
                edge = level0.edge(edgeLabel)
                p = edge.flag(flag_constants.SCISSOR_PROTECTION)
                self.contour.append((edgeLabel, p))
                if not p:
                    edge.setFlag(flag_constants.SCISSOR_PROTECTION)

    def redo(self):
        level0 = self.workspace._level0
        for edgeLabel, _ in self.contour:
            edge = level0.edge(edgeLabel)
            edge.setFlag(flag_constants.SCISSOR_PROTECTION)
        # FIXME: if all protected edges still exist, protect them
        # and skip recomputeAutomaticLevels():
        self.workspace.recomputeAutomaticLevels()

    def undo(self):
        level0 = self.workspace._level0
        for edgeLabel, p in self.contour:
            edge = level0.edge(edgeLabel)
            edge.setFlag(flag_constants.SCISSOR_PROTECTION, p)
        # FIXME: remove protection flag from edges, or is
        # recomputeAutomaticLevels always needed?  (I think so.)
        self.workspace.recomputeAutomaticLevels()

    def __str__(self):
        return "Scissors" # TODO: path length, closedness?

# --------------------------------------------------------------------
#                          Main Workspace Class
# --------------------------------------------------------------------

class Workspace(mapdisplay.MapDisplay):
    """Workspace for region-based segmentation.

    The _levelSlider has the levelIndex as value.  This is 0 for the
    pyramid's bottom (where faceCount ==
    self.manualBaseMap().faceCount), and is maximal for the apex.
    (The number of faces in the apex is estimated in manualBaseMap, in
    order to set the slider range appropriately.)
    """
    __base = mapdisplay.MapDisplay

    # actually, this has only documenting effect; since this inherits
    # PyQt widgets, any attribute may be used:
    __slots__ = ("_level0", "_manualCK", "_seeds",
                 "_pyramidCK", "activeCostMeasure", "dynamicCosts",
                 "_mapRestartAction", "_levelSlider",
                 "_manualBaseMapFaceCount", "_estimatedApexFaceCount",
                 "_history", "_activeTool", "colorSpace")
    
    def __init__(self, level0, originalImage, bi = None):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0
        self._manualCK = range(level0.maxFaceLabel())
        self._seeds = None
        self._protectedFaces = []
        self._pyramidCK = None
        self.activeCostMeasure = 1 # faceHomogeneity
        self.dynamicCosts = True
        self._history = []
        self._activeTool = None
        self.colorSpace = ""
        self._estimatedApexFaceCount = 2
        self._sliderMode = 0
        self._sliderCosts = None

        # needed for backpropagation of protection:
        if not hasattr(self.map, "mergedEdges"):
            self.map.mergedEdges = statistics.MergedEdges(self.map)

        reinitIcon = QtGui.QPixmap()
        reinitIcon.loadFromData(icons.reinitIconPNGData, "PNG")
        ra = QtGui.QAction(self)
        ra.setIcon(QtGui.QIcon(reinitIcon))
        ra.setText("&Re-start from level 0")
        ra.setToolTip("Re-start with original level 0 map")
        self.ui.Tools.addAction(ra)
        self.connect(ra, QtCore.SIGNAL("activated()"), self.restart)
        self._mapRestartAction = ra

        self.connect(self.ui.undoAction, QtCore.SIGNAL("activated()"), self.undo)

        # set up HBox with cost measure options:
        automaticOptions = QtGui.QWidget(self._imageWindow)
        l = QtGui.QHBoxLayout(automaticOptions)
        l.setMargin(2)
        l.setSpacing(6)

        cml = QtGui.QLabel("Merge &cost measure:", automaticOptions)
        l.addWidget(cml)
        cmChooser = QtGui.QComboBox(automaticOptions)
        for i, name in enumerate(self.costMeasureNames):
            cmChooser.insertItem(i, name)
        cmChooser.setCurrentIndex(self.activeCostMeasure)
        self.connect(cmChooser, QtCore.SIGNAL("activated(int)"), self.setCostMeasure)
        l.addWidget(cmChooser)
        cml.setBuddy(cmChooser)
        self._cmChooser = cmChooser
        
        csl = QtGui.QLabel("C&olor space:", automaticOptions)
        l.addWidget(csl)
        csChooser = QtGui.QComboBox(automaticOptions)
        for i, name in enumerate(self.colorSpaceNames):
            csChooser.insertItem(i, name)
        self.connect(csChooser, QtCore.SIGNAL("activated(int)"), self.setColorSpace)
        l.addWidget(csChooser)
        csl.setBuddy(csChooser)
        self._csChooser = csChooser

        self.dynamicCheckBox = QtGui.QCheckBox("&Dynamic costs", automaticOptions)
        l.addWidget(self.dynamicCheckBox)
        self.dynamicCheckBox.setChecked(self.dynamicCosts)
        self.connect(self.dynamicCheckBox, QtCore.SIGNAL("toggled(bool)"),
                     self.setDynamicCosts)
        
        l.addItem(QtGui.QSpacerItem(
            10, 1, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        self._imageWindow._layout.insertWidget(0, automaticOptions)

        # set up HBox with level slider:
        sliderBox = QtGui.QWidget(self._imageWindow)
        l = QtGui.QHBoxLayout(sliderBox)
        l.setMargin(2)
        l.setSpacing(6)

        self.sliderModeChooser = QtGui.QComboBox(sliderBox)
        self.sliderModeChooser.insertItem(0, "Steps:")
        self.sliderModeChooser.insertItem(1, "Steps*:")
        self.sliderModeChooser.insertItem(2, "Cost:")
        self.connect(self.sliderModeChooser, QtCore.SIGNAL("activated(int)"),
                     self.setSliderMode)
        l.addWidget(self.sliderModeChooser)

        self._levelSlider = QtGui.QSlider(sliderBox)
        self._levelSlider.setOrientation(QtCore.Qt.Horizontal)
        #self._levelSlider.setSteps(-1, -10)
        self._levelSlider.setRange(
            0, level0.faceCount - self._estimatedApexFaceCount)
        self.connect(self._levelSlider, QtCore.SIGNAL("valueChanged(int)"),
                     self._levelSliderChanged)
        l.addWidget(self._levelSlider)

        self._imageWindow._layout.insertWidget(1, sliderBox)

        # (finished setting up widgets)
        if self.isVisible():
            automaticOptions.show()
            sliderBox.show()

        if bi is not None:
            self.setImage(bi, role = "bi")

    def setTool(self, tool = None):
        self._activeTool = tool
        self.__base.setTool(self, tool)
        if not self.tool:
            return
        if isinstance(self.tool, tools.IntelligentScissors):
            tools.activeCostMeasure = \
                statistics.HyperbolicInverse(self.costMeasure(self.map))
        self.connect(self.tool, QtCore.SIGNAL("paintbrushFinished"),
                     self.paintbrushFinished)
        self.connect(self.tool, QtCore.SIGNAL("faceProtectionChanged"),
                     self.faceProtectionChanged)
        self.connect(self.tool, QtCore.SIGNAL("contourFinished"),
                     self.scissorsFinished)

        if hasattr(self.tool, "setSeeds"): # SeedSelector?
            if self._seeds is not None:
                self.tool.setSeeds(list(self._seeds))
            self.connect(self.tool, QtCore.SIGNAL("seedAdded"),
                         self.seedAdded)
            self.connect(self.tool, QtCore.SIGNAL("seedRemoved"),
                         self.seedRemoved)

    def _perform(self, action):
        # FIXME: redo support
        self._history.append(action)

    def undo(self):
        # FIXME: redo support
        self._history[-1].undo()
        del self._history[-1]

    def paintbrushFinished(self, survivor):
        s = PaintbrushStroke(self, survivor.label())
        if not s.mergeCount:
            return
        self._perform(s)
        self._pyramidCK = None
        # FIXME: in theory, we would need recomputeAutomaticLevels here
        # (the performed operations could differ, if the manually
        # merged faces result in significantly differing costs)

    def faceProtectionChanged(self, face):
        #assert face.map() == self.map
        self._perform(FaceProtection(self, face.label()))

        bbox = face.boundingBox()
        updateRect = QtCore.QRect(
            self.viewer.windowCoordinate(*bbox.begin()),
            self.viewer.windowCoordinate(*bbox.end()))
        lw = self.edgeOverlay.width + 1
        updateRect.adjust(-lw, -lw, lw + 1, lw + 1)
        self.viewer.update(updateRect)

    def scissorsFinished(self, contour):
        self._perform(ScissorsProtection(self, contour))

    def seedAdded(self, pos):
        if self._seeds is None:
            self._seeds = []
        self._seeds.append(pos)
        self.recomputeAutomaticLevels(faceCountOffset = 1)

    def seedRemoved(self, pos):
        self._seeds.remove(pos)
        if not self._seeds:
            self._seeds = None # switch back to ARM
        self.recomputeAutomaticLevels(faceCountOffset = -1)

    def setMap(self, map):
        """Make the given map (which must belong to the current
        pyramid) the currently displayed level."""
        maputils.detachMapStats(self.map)

        activeTool = self._activeTool
        self.__base.setMap(self, map)
        if activeTool:
            self.setTool(activeTool)

        self._updateLevelSlider()

    def restart(self):
        """Restart with the manually created level."""
        self.setMap(self.manualBaseMap())

    def manualBaseMap(self):
        p = progress.StatusMessage("restoring pyramid bottom")
        result = copy.deepcopy(self._level0)
        if not hasattr(result, "mergedEdges"):
            result.mergedEdges = statistics.MergedEdges(result)
        #p = progress.StatusMessage("  applying manual changes")
        maputils.applyFaceClassification(result, self._manualCK)
        self._manualBaseMapFaceCount = result.faceCount

        #p = progress.StatusMessage("  applying face protection + seeds")
        self._estimatedApexFaceCount = 2 # infinite + one remaining finite
        faceLabelLUT = result.faceLabelLUT()
        for faceLabel in self._protectedFaces:
            face = result.face(faceLabelLUT[faceLabel])
            if not face.flag(flag_constants.PROTECTED_FACE):
                maputils.protectFace(face)
                self._estimatedApexFaceCount += 1

        if self._seeds:
            self._estimatedApexFaceCount -= 1 # one seed is for the remaining finite
            for pos in self._seeds:
                face = result.faceAt(pos)
                if face.flag(flag_constants.PROTECTED_FACE):
                    continue
                if not face.flag(flag_constants.SRG_SEED):
                    face.setFlag(flag_constants.SRG_SEED)
                    self._estimatedApexFaceCount += 1

        self._levelSlider.blockSignals(True)
        self._levelSlider.setRange(
            0, result.faceCount - self._estimatedApexFaceCount)
        self._levelSlider.blockSignals(False)
        p.finish()
        return result

    def faceMeans(self, map):
        attr = "faceMeans" + self.colorSpace
        if not hasattr(map, attr):
            img = self.images.get("colored", self.images["original"])
            if self.colorSpace:
                if img.bands() != 3:
                    sys.stderr.write("WARNING: Color space '%s' requested but ignored; image has only %d band!\n" % (img.bands(), ))
                else:
                    rgbNames = self.colorSpaceNames[-3:]
                    if self.colorSpace in rgbNames:
                        img = vigra.GrayImage(img[rgbNames.index(self.colorSpace)])
                    else:
                        img = vigra.transformImage(
                            img, "\l x: RGB2%s(x)" % self.colorSpace)
            stats = statistics.FaceColorStatistics(map, img)
            stats.image = img
            setattr(map, attr, stats)

            # let MapDisplay make use of attr 'faceMeans' if applicable:
            # FIXME - 'map' != self.map here, so the action stays disabled..
            self._enableImageActions()
        return getattr(map, attr)

    def edgeGradients(self, map, quantiles = False):
        if quantiles and hasattr(map, "egs") and not map.egs.supportsQuantile():
            # FIXME: this should only be done if map.egs is from us:
            map.egs.detachHooks()
            del map.egs
        if not hasattr(map, "egs"):
            bi = self.images["bi"]
            if hasattr(bi, "siv"):
                biSiv = bi.siv
            else:
                biSiv = vigra.SplineImageView5(bi)
            Functor = quantiles \
                      and geomap.QuantileStatistics \
                      or geomap.PolylineStatistics
            map.egs = statistics.EdgeGradientStatistics(
                map, biSiv, Functor = Functor)
        return map.egs

    def setDynamicCosts(self, dc):
        if dc != self.dynamicCosts:
            self.dynamicCosts = dc
            self.recomputeAutomaticLevels()
            self.dynamicCheckBox.setChecked(dc)

    def setColorSpace(self, colorSpace):
        if isinstance(colorSpace, int):
            colorSpace = self.colorSpaceNames[colorSpace]
        else:
            assert colorSpace in self.colorSpaceNames, \
                   "invalid colorSpace '%s'" % colorSpace
        if colorSpace == "RGB":
            colorSpace = ""
        if self.colorSpace != colorSpace:
            self.colorSpace = colorSpace
            self.recomputeAutomaticLevels()
            # FIXME: for scissors:?
#             tools.activeCostMeasure = \
#                 statistics.HyperbolicInverse(self.costMeasure(self.map))
            self._csChooser.setCurrentIndex(
                self.colorSpaceNames.index(colorSpace or "RGB"))

    colorSpaceNames = ["RGB", "RGBPrime", "Luv", "Lab", "red", "green", "blue"]

    def setCostMeasure(self, index):
        if isinstance(index, str):
            index = self.costMeasureNames.index(index)
        elif index < 0 or index > len(self.costMeasureNames):
            sys.stderr.write("ERROR: Invalid cost measure %d, give\n" % index)
            for i, name in enumerate(self.costMeasureNames):
                sys.stderr.write("  %d for %s\n" % (i, name))
            return
        if index in (4,5) and not hasattr(self._level0, "wsStats"):
            sys.stderr.write("ERROR: cost measure %d not possible -- no WatershedStatistics available!\n" % index)
            return
        if index == 5 and not hasattr(self._level0, "wsBasinStats"):
            sys.stderr.write("ERROR: cost measure %d not possible -- no WatershedBasinStatistics available!\n" % index)
            return
        if index in range(6, 10) and "bi" not in self.images:
            sys.stderr.write("ERROR: cost measure %d not possible -- no boundary indicator image available!\n" % index)
            return
        #if index == 8 and not foo.supportsQuantile():
        if self.activeCostMeasure != index:
            self.activeCostMeasure = index
            self.recomputeAutomaticLevels()
            tools.activeCostMeasure = \
                statistics.HyperbolicInverse(self.costMeasure(self.map))
            self._cmChooser.setCurrentIndex(index)

    costMeasureNames = ["face mean difference",
                        "face homogeneity",
                        "face t-test",
                        "face brightness",
                        "pass value",
                        "watershed dynamics",
                        "minimum gradient magnitude",
                        "average gradient magnitude",
                        "median gradient magnitude",
                        "maximum gradient magnitude",
                        "isoperimetric quotient of survivor",
                        "contour length"]

    def costMeasure(self, map):
        """Instantiate and return the currently chosen type of cost
        measure for the given map."""

        if self.activeCostMeasure == 0:
            return self.faceMeans(map).faceMeanDiff
        elif self.activeCostMeasure == 1:
            return self.faceMeans(map).faceHomogeneity
        elif self.activeCostMeasure == 2:
            return self.faceMeans(map).faceTTest
        elif self.activeCostMeasure == 3:
            def brightness(dart, fm = self.faceMeans(map)):
                return max(vigra.norm(fm[dart.leftFaceLabel()]),
                           vigra.norm(fm[dart.rightFaceLabel()]))
            return brightness
        elif self.activeCostMeasure == 4:
            return map.wsStats.dartPassValue
        elif self.activeCostMeasure == 5:
            return map.wsStats.dynamics
        elif self.activeCostMeasure == 6:
            return self.edgeGradients(map).dartMin
        elif self.activeCostMeasure == 7:
            return self.edgeGradients(map).average
        elif self.activeCostMeasure == 8:
            return self.edgeGradients(map, True).quantile(0.5)
        elif self.activeCostMeasure == 9:
            return self.edgeGradients(map).max
        elif self.activeCostMeasure == 10:
            return statistics.mergedIsoperimetricQuotient
        elif self.activeCostMeasure == 11:
            return statistics.mergedContourLength
        else:
            raise RuntimeError("Wrong cost measure (%s)" % self.activeCostMeasure)

    def _levelSliderChanged(self, sliderValue):
        if self._sliderMode == 0:
            self.displayLevel(levelIndex = sliderValue)
            return

        p = float(sliderValue) / self._levelSlider.maximum()
        if self._sliderCosts:
            self.displayLevel(cost = self._sliderCosts[-1] * p)
        else:
            self.displayLevel(
                levelIndex = int(p**0.2 * self._levelSlider.maximum()))

    def _updateLevelSlider(self):
        self._levelSlider.blockSignals(True)
        if self._sliderMode == 0:
            # FIXME: why not level0.faceCount - map.faceCount?
            self._levelSlider.setValue(
                self._levelSlider.maximum() -
                (self.map.faceCount - self._estimatedApexFaceCount))
        else:
            pass # FIXME
        self._levelSlider.blockSignals(False)

    def recomputeAutomaticLevels(self, faceCountOffset = 0):
        self.displayLevel(faceCount = self.map.faceCount + faceCountOffset,
                          force = True)

    def displayedLevelIndex(self):
        # FIXME: quick hack, this API is needed, impl unsure:
        return self._levelSlider.value()

    def displayLevel(self, levelIndex = None, faceCount = None,
                     cost = None, force = False):
        """Display the specified level with the desired number of faces.

        If `faceCount` is not given, `levelIndex` must be given and
        faceCount is calculated accordingly (from the slider range).

        If the displayed map has more faces than `faceCount`, and a
        self._pyramidCK is available, it is reduced incrementally (and
        nothing happens if the number of faces is already the desired
        one)."""

        assert (levelIndex is not None) + \
               (faceCount is not None) + \
               (cost is not None) == 1, \
               "displayLevel: give exactly one of levelIndex, faceCount, or cost!"
        if faceCount is None:
            if cost is not None:
                # FIXME: bisect?
                for i, c in enumerate(self._costLog):
                    if c > cost:
                        levelIndex = max(0, i-1)
                        break

            faceCount = (self._levelSlider.maximum() -
                         (levelIndex - self._estimatedApexFaceCount))

        if force:
            self._pyramidCK = None # force recomputation of pyramid CK

#         print "trying to reach faceCount %d, having %d now..." % (
#             faceCount, self.map.faceCount)

        if not self._pyramidCK or faceCount > self.map.faceCount:
            map = self.manualBaseMap()
            if faceCount < map.faceCount:
                self.pyramidCK().applyToFaceCount(map, faceCount)
            self.setMap(map)
        elif faceCount < self.map.faceCount:
            self.pyramidCK().applyToFaceCount(self.map, faceCount)
            self._updateLevelSlider()

    def pyramidCK(self):
        if not self._pyramidCK:
            self.startAutomaticMethod()
        return self._pyramidCK

    def setSliderMode(self, sm):
        self._sliderMode = sm
        if sm != 2:
            self._sliderCosts = None
            self._updateLevelSlider()
            return

        if not self._pyramidCK:
            self.startAutomaticMethod()
        self._sliderCosts = self._costLog
#         self._sliderCosts = []
#         for i, c in enumerate(self._costLog):
#             self._sliderCosts.append(c) # if ...
        self._updateLevelSlider()

    def startAutomaticMethod(self):
        map = self.manualBaseMap()
        map.pyramidCK = PyramidContractionKernel(map)

        if self._seeds is None:
            methodName = "automatic region merging"
            am = maputils.AutomaticRegionMerger(
                map, self.costMeasure(map), updateNeighborHood = self.dynamicCosts)
        else:
            methodName = "seeded region growing"
            am = maputils.SeededRegionGrowing(
                map, self.costMeasure(map), dynamic = self.dynamicCosts)

        am._costLog = []
        self._costLog = am._costLog
        stepsTotal = map.faceCount - self._estimatedApexFaceCount
        p = progress.ProgressHook(
            progress.StatusMessage(methodName)) \
            .rangeTicker(stepsTotal / 50)
        while am.mergeSteps(50):
            p()

        maputils.detachMapStats(map)
        self._pyramidCK = map.pyramidCK

def main(filename, biScale = 1.6, saddleThreshold = 0.2):
    import bi_utils
    img = vigra.readImage(filename)
    gm, grad = bi_utils.gaussianGradient(img, biScale)
    wsm = maputils.subpixelWatershedMap(
        gm, saddleThreshold = saddleThreshold)

    return Workspace(wsm, img, bi = gm)

if __name__ == "__main__":
    import sys
    #filename = "../../../Testimages/lenna_original_color.png"
    filename = "../../../Testimages/blox.png"
    biScale = 1.6
    saddleThreshold = 0.2
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    if len(sys.argv) > 2:
        biScale = float(sys.argv[2])

    app = QtGui.QApplication(sys.argv)
    w = main(filename, biScale)
    app.connect(app, QtCore.SIGNAL("lastWindowClosed()"), app, QtCore.SLOT("quit()"))
    app.exec_()
