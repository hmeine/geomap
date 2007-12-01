import copy, qt
import vigra, mapdisplay, icons, maputils, statistics, flag_constants, tools
import progress

def labelRoot(lut, label):
    result = lut[label]
    if lut[result] != result:
        result = labelRoot(lut, result)
        lut[label] = result
    return result

class PyramidContractionKernel(statistics.DetachableStatistics):
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

class PaintbrushStroke(object):
    __slots__ = ("workspace", "survivorLabel", "faces", "oldValues")

    def __init__(self, workspace, survivorLabel):
        self.workspace = workspace
        self.survivorLabel = survivorLabel
        self.faces = list(workspace.map.faceLabelLUT().merged(survivorLabel))
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

class FaceProtection(object):
    __slots__ = ("workspace", "dartLabel", "protected")

    def __init__(self, workspace, dartLabel):
        self.workspace = workspace
        self.dartLabel = dartLabel
        self.protected = workspace.map.dart(dartLabel).leftFace().flag(
            flag_constants.PROTECTED_FACE)
        self.redo()

    def protect(self, protect):
        if protect:
            self.workspace._protectedFaceAnchors.append(self.dartLabel)
        else:
            self.workspace._protectedFaceAnchors.remove(self.dartLabel)
        maputils.protectFace(
            self.workspace.map.dart(self.dartLabel).leftFace(), protect)

    def redo(self):
        self.protect(self.protected)

    def undo(self):
        self.protect(not self.protected)

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
                 "_mapRestartAction", "_levelSlider", "_estimatedApexFaceCount",
                 "_history", "_activeTool", "colorSpace")
    
    def __init__(self, level0, originalImage, bi = None):
        self.__base.__init__(self, copy.deepcopy(level0), originalImage)
        self._level0 = level0
        self._manualCK = range(level0.maxFaceLabel())
        self._seeds = None
        self._protectedFaceAnchors = []
        self._pyramidCK = None
        self.activeCostMeasure = 1 # faceHomogeneity
        self.dynamicCosts = True
        self._history = []
        self._activeTool = None
        self.colorSpace = ""

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

        # set up HBox with cost measure options:
        automaticOptions = qt.QWidget(self._imageWindow, "automaticOptions")
        l = qt.QHBoxLayout(automaticOptions, 2, 6)

        cml = qt.QLabel("Merge &cost measure:", automaticOptions)
        l.addWidget(cml)
        self.cmChooser = qt.QComboBox(automaticOptions, "cmChooser")
        for name in self.costMeasureNames:
            self.cmChooser.insertItem(name)
        self.cmChooser.setCurrentItem(self.activeCostMeasure)
        self.connect(self.cmChooser, qt.SIGNAL("activated(int)"), self.setCostMeasure)
        l.addWidget(self.cmChooser)
        cml.setBuddy(self.cmChooser)
        
        csl = qt.QLabel("C&olor space:", automaticOptions)
        l.addWidget(csl)
        csChooser = qt.QComboBox(automaticOptions, "csChooser")
        for name in self.colorSpaceNames:
            csChooser.insertItem(name)
        self.connect(csChooser, qt.SIGNAL("activated(int)"), self.setColorSpace)
        l.addWidget(csChooser)
        csl.setBuddy(csChooser)

        self.dynamicCheckBox = qt.QCheckBox("&Dynamic costs", automaticOptions)
        l.addWidget(self.dynamicCheckBox)
        self.dynamicCheckBox.setChecked(self.dynamicCosts)
        self.connect(self.dynamicCheckBox, qt.SIGNAL("toggled(bool)"),
                     self.setDynamicCosts)
        
        l.addItem(qt.QSpacerItem(
            10, 1, qt.QSizePolicy.Expanding, qt.QSizePolicy.Minimum))
        self._imageWindow._layout.insertWidget(0, automaticOptions)

        #sws.imageFrame.addWidget(self.imageFrame)
        self._estimatedApexFaceCount = 2
        self._levelSlider = qt.QSlider(self._imageWindow, "_levelSlider")
        self._levelSlider.setOrientation(qt.QSlider.Horizontal)
        self._levelSlider.setSteps(-1, -10)
        self._levelSlider.setRange(
            0, level0.faceCount - self._estimatedApexFaceCount)
        self.connect(self._levelSlider, qt.SIGNAL("valueChanged(int)"),
                     self._levelSliderChanged)
        self._imageWindow._layout.insertWidget(1, self._levelSlider)
        if self.isShown():
            automaticOptions.show()
            self._levelSlider.show()

        if bi is not None:
            self.setImage(bi, role = "bi")

    @staticmethod
    def _detachMapStats(map):
        for a in map.__dict__:
            o = getattr(map, a)
            if hasattr(o, "detachHooks"):
                #print "detaching hooks of", o
                o.detachHooks()

    def setTool(self, tool = None):
        if tool:
            assert type(tool) == int, "Workspace needs to remember the tool - only ints allowed!"
        self._activeTool = tool
        self.__base.setTool(self, tool)
        if not self.tool:
            return
        if isinstance(self.tool, tools.IntelligentScissors):
            tools.activeCostMeasure = \
                statistics.HyperbolicInverse(self.costMeasure(self.map))
        self.connect(self.tool, qt.PYSIGNAL("paintbrushFinished"),
                     self.paintbrushFinished)
        self.connect(self.tool, qt.PYSIGNAL("faceProtectionChanged"),
                     self.faceProtectionChanged)
        if hasattr(self.tool, "setSeeds"): # SeedSelector?
            if self._seeds is not None:
                self.tool.setSeeds(list(self._seeds))
            self.connect(self.tool, qt.PYSIGNAL("seedAdded"),
                         self.seedAdded)
            self.connect(self.tool, qt.PYSIGNAL("seedRemoved"),
                         self.seedRemoved)

    def _perform(self, action):
        # FIXME: redo support
        self._history.append(action)

    def undo(self):
        # FIXME: redo support
        self._history[-1].undo()
        del self._history[-1]

    def paintbrushFinished(self, survivor):
        self._perform(PaintbrushStroke(self, survivor.label()))
        self._pyramidCK = None
        # FIXME: in theory, we would need recomputeAutomaticLevels here
        # (the performed operations could differ, if the manually
        # merged faces result in significantly differing costs)

    def faceProtectionChanged(self, face):
        #assert face.map() == self.map
        self._perform(FaceProtection(self, face.contour().label()))

        bbox = face.boundingBox()
        updateRect = qt.QRect(
            self.viewer.toWindowCoordinates(*bbox.begin()),
            self.viewer.toWindowCoordinates(*bbox.end()))
        lw = self.edgeOverlay.width + 1
        updateRect.addCoords(-lw, -lw, lw, lw)
        self.viewer.update(updateRect)

#         for dart in maputils.contourDarts(face):
#             protected = dart.edge().flag(flag_constants.CONTOUR_PROTECTION)
#             for edgeLabel in self.map.mergedEdges[dart.edgeLabel()]:
#                 self._level0.edge(edgeLabel).setFlag(
#                     flag_constants.CONTOUR_PROTECTION, protected)

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
        self._detachMapStats(self.map)

        activeTool = self._activeTool
        self.__base.setMap(self, map)
        if activeTool:
            self.setTool(activeTool)
        
        self._levelSlider.blockSignals(True)
        self._levelSlider.setValue(
            self._levelSlider.maxValue() -
            (self.map.faceCount - self._estimatedApexFaceCount))
        self._levelSlider.blockSignals(False)

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

        #p = progress.StatusMessage("  applying face protection + seeds")
        self._estimatedApexFaceCount = 2 # infinite + one remaining finite
        for dartLabel in self._protectedFaceAnchors:
            face = result.dart(dartLabel).leftFace()
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
        return result

    def faceMeans(self, map):
        attr = "faceMeans" + self.colorSpace
        if not hasattr(map, attr):
            img = self.images.get("colored", self.images["original"])
            if self.colorSpace:
                if img.bands() != 3:
                    sys.stderr.write("WARNING: Color space '%s' requested but ignored; image has only %d band!\n" % (img.bands(), ))
                else:
                    img = vigra.transformImage(
                        img, "\l x: RGB2%s(x)" % self.colorSpace)
            stats = statistics.FaceColorStatistics(map, img)
            stats.image = img
            setattr(map, attr, stats)
        return getattr(map, attr)

    def setDynamicCosts(self, dc):
        if dc != self.dynamicCosts:
            self.dynamicCosts = dc
            self.recomputeAutomaticLevels()
            self.dynamicCheckBox.setChecked(dc)

    def setColorSpace(self, colorSpace):
        if isinstance(colorSpace, int):
            colorSpace = self.colorSpaceNames[colorSpace]
        if colorSpace == "RGB":
            colorSpace = ""
        if self.colorSpace != colorSpace:
            self.colorSpace = colorSpace
            self.recomputeAutomaticLevels()
            # FIXME: for scissors:?
#             tools.activeCostMeasure = \
#                 statistics.HyperbolicInverse(self.costMeasure(self.map))

    colorSpaceNames = ["RGB", "RGBPrime", "Luv", "Lab"]

    def setCostMeasure(self, index):
        if self.activeCostMeasure != index:
            self.activeCostMeasure = index
            self.recomputeAutomaticLevels()
            tools.activeCostMeasure = \
                statistics.HyperbolicInverse(self.costMeasure(self.map))
            self.cmChooser.setCurrentItem(index)

    costMeasureNames = ["face mean difference",
                        "face homogeneity",
                        "face t-test",
                        "face brightness",
                        "isoperimetric quotient of result",
                        "isoperimetric quotient after/before",
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
            return statistics.mergedIsoperimetricQuotient
        elif self.activeCostMeasure == 5:
            return statistics.seedIsoperimetricQuotient
        elif self.activeCostMeasure == 6:
            return statistics.mergedContourLength
        else:
            raise RuntimeError("Wrong cost measure (%s)" % self.activeCostMeasure)

    def _levelSliderChanged(self, levelIndex):
        self.displayLevel(levelIndex = levelIndex)

    def recomputeAutomaticLevels(self, faceCountOffset = 0):
        self.displayLevel(faceCount = self.map.faceCount + faceCountOffset,
                          force = True)

    def displayLevel(self, levelIndex = None, faceCount = None,
                     force = False):
        """Display the specified level with the desired number of faces.

        If `faceCount` is not given, `levelIndex` must be given and
        faceCount is calculated accordingly (from the slider range).

        If the displayed map has more faces than `faceCount`, and a
        self._pyramidCK is available, it is reduced incrementally (and
        nothing happens if the number of faces is already the desired
        one)."""
        
        assert (levelIndex is None) != (faceCount is None), \
               "displayLevel: give exactly one of levelIndex or faceCount!"
        if faceCount is None:
            faceCount = (self._levelSlider.maxValue() -
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
            
            self._levelSlider.blockSignals(True)
            self._levelSlider.setValue(
                self._levelSlider.maxValue() -
                (self.map.faceCount - self._estimatedApexFaceCount))
            self._levelSlider.blockSignals(False)

    def pyramidCK(self):
        if not self._pyramidCK:
            self.startAutomaticMethod()
        return self._pyramidCK

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

        self._detachMapStats(map)
        self._pyramidCK = map.pyramidCK
