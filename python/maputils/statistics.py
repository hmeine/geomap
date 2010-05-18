import math, string, copy, weakref
from vigra import *
import geomap
from geomap import \
     FaceGrayStatistics, FaceRGBStatistics, \
     resamplePolygon, tangentList, tangentListGaussianReflective
import sivtools, flag_constants

class DetachableStatistics(object):
    """Base class for all dynamic statistics.

    Manages a sequence of hooks in an attribute _attachedHooks,
    allowing the user to call `detachHooks` for disconnecting from the
    map.

    Subclasses may provide a method _attachHooks() which should use
    the weak reference self._map() for reconnecting the hooks.
    This then allows proper pickle support (cf. `__getstate__` /
    `__setstate__`, which manage the map and the reconnection of the
    hooks)."""
    
    __slots__ = ["_attachedHooks", "_map"]
    
    def __init__(self, map):
        # prevent cycles if this is an attribute of the map:
        self._map = weakref.ref(map) # only needed for pickle support

    def detachHooks(self):
        for cb in self._attachedHooks:
            cb.disconnect()
        self._attachedHooks = ()

    def __getstate__(self):
        return (self._map(), )

    def __setstate__(self, (map, )): # cf. __init__
        self._map = weakref.ref(map)
        self._attachHooks()

    def __repr__(self):
        map = ", map destroyed"

        if not hasattr(self, "_attachedHooks"):
            hooks = "not attached"
            map = ""
        elif not self._attachedHooks:
            hooks = "no hooks"
        else:
            active = self._attachedHooks[0].connected() and "active" or "detached"
            hooks = "%d %s hooks" % (len(self._attachedHooks), active)

        if self._map():
            addr = id(self._map()) # + 0x100000000L
            map = " for GeoMap @0x%8x" % addr

        return "<%s, %s%s>" % (self.__class__.__name__, hooks, map)

def _combinedMeasure(dart, weightedMeasures):
    cost = 0.0
    for weight, measure in weightedMeasures:
        cost += weight * measure(dart)
    return cost

def combinedMeasure(*weightedMeasures):
    return lambda dart: _combinedMeasure(dart, weightedMeasures)

def trainingData(dartPath, faceMeans, edgeGradients):
    Functor = type(faceMeans.functor(dartPath[0].leftFaceLabel()))
    left = Functor()
    right = Functor()
    grad = geomap.FaceGrayStatistics.Functor()
    leftVariance = Functor()
    rightVariance = Functor()
    for dart in dartPath:
        weight = dart.edge().length()
        left(faceMeans.average(dart.leftFaceLabel()), weight)
        right(faceMeans.average(dart.rightFaceLabel()), weight)
        leftVariance(faceMeans.variance(dart.leftFaceLabel(), True), weight)
        rightVariance(faceMeans.variance(dart.rightFaceLabel(), True), weight)
        grad(edgeGradients.dartAverage(dart), weight)

    leftSigma2 = norm(left.variance(True) + leftVariance.average())
    leftAvg = left.average()
    rightSigma2 = norm(right.variance(True) + rightVariance.average())
    rightAvg = right.average()
    gradSigma2 = grad.variance(True)
    gradAvg = grad.average()

    return leftSigma2, leftAvg, rightSigma2, rightAvg, gradSigma2, gradAvg

def trainedMeasure(dartPath, faceMeans, edgeGradients, weights = (1, 1, 1)):
    leftSigma2, leftAvg, rightSigma2, rightAvg, gradSigma2, gradAvg = \
                trainingData(dartPath, faceMeans, edgeGradients)
    w1, w2, w3 = weights
    weightNorm = 1.0/sum(weights)
    return lambda dart: weightNorm * (
        w1*math.exp(-vigra.sq(faceMeans[dart.leftFaceLabel()]-leftAvg)/leftSigma2) + \
        w2*math.exp(-vigra.sq(faceMeans[dart.rightFaceLabel()]-rightAvg)/rightSigma2) + \
        w3*math.exp(-vigra.sq(edgeGradients.dartAverage(dart))/gradSigma2))

# --------------------------------------------------------------------
#              Region-based Statistics & Cost Measures
# --------------------------------------------------------------------

def clampEpsilon(variance, epsilon):
    if type(variance) == float:
        return max(epsilon, variance)
    for i in range(len(variance)):
        variance[i] = max(epsilon, variance[i])
    return variance

def safeSqrt(variance):
    if type(variance) == float:
        return math.sqrt(max(0.0, variance))
    result = [math.sqrt(max(0.0, variance[i])) for i in range(len(variance))]
    return type(variance)(*result)

class FaceMeanFunctor(object):
    __slots__ = ["pixelCount", "sum", "sum2", "defaultValue"]
    
    def __init__(self, defaultValue):
        self.pixelCount = None
        self.defaultValue = defaultValue

    def average(self):
        if self.pixelCount:
            return self.sum / self.pixelCount
        return self.defaultValue

    def variance(self):
        if self.pixelCount:
            return (self.sum2 / self.pixelCount) - \
                   math.sq(self.sum / self.pixelCount)
        return self.defaultValue

    def stdDeviation(self):
        return safeSqrt(self.variance())

    def __call__(self, value):
        if self.pixelCount:
            self.pixelCount += 1
            self.sum += value
            self.sum2 += math.sq(value)
        else:
            self.pixelCount = 1
            self.sum = value
            self.sum2 = math.sq(value)

    def merge(self, otherStats):
        if not otherStats.pixelCount:
            return
        if self.pixelCount:
            self.pixelCount += otherStats.pixelCount
            self.sum += otherStats.sum
            self.sum2 += otherStats.sum2
        else:
            self.pixelCount = otherStats.pixelCount
            self.sum = otherStats.sum
            self.sum2 = otherStats.sum2

class DynamicFaceStatistics(DetachableStatistics):
    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeFacesCallbacks(self.preMergeFaces, self.postMergeFaces),
            self._map().addAssociatePixelsCallback(self.associatePixels))

from numpy import arange

def sampleRegions(map, image, functors):
    """sampleRegions(map, image, functors)

    For each pixel in the map's labelImage() whose facet belongs
    completely to one region, calls the corresponding
    functors[faceLabel] with the corresponding pixel from the given
    image."""

    class LookupFaceFunctors(object):
        def __init__(self, functors):
            self._functors = functors

        def __call__(self, label, value):
            if label >= 0:
                self._functors[int(label)](value)

    inspectImage(map.labelImage(), image,
                 LookupFaceFunctors(functors))

def superSample(face, level = 2):
    bbox = face.boundingBox()
    xRange = arange(bbox.begin()[0], bbox.end()[0], 1.0/level)
    for y in arange(bbox.begin()[1], bbox.end()[1], 1.0/level):
        for x in xRange:
            pos = Vector2(x, y)
            if face.contains(pos):
                yield pos

# new API: does not touch the Face objects themselves
class _FaceColorStatistics(DynamicFaceStatistics):
    def __init__(self, map, originalImage, minSampleCount = 1,
                 defaultValue = None, SIV = SplineImageView5):
        DynamicFaceStatistics.__init__(self, map)
        self.originalImage = originalImage

        if defaultValue == None:
            defaultValue = 0.0
            if originalImage.bands() > 1:
                # initialize zero Vector of appropriate size:
                defaultValue = originalImage[0, 0]
                for i in range(originalImage.bands()):
                    defaultValue[i] = 0.0

        self._functors = [None] * map.maxFaceLabel()
        for face in map.faceIter():
            self._functors[face.label()] = FaceMeanFunctor(defaultValue)
        self._defaultValue = defaultValue
        self._diffNorm = 255.*math.sqrt(self.bands())
        if originalImage.bands() > 1:
            self.faceMeanDiff = self.faceMeanDiffColor
        else:
            self.faceMeanDiff = self.faceMeanDiffGray

        sampleRegions(map, originalImage, self._functors)

        self._SIV = SIV # class
        self._origSIV = None # instance, lazily-initialized using self._SIV
        self._superSampled = [0] * map.maxFaceLabel()
        if minSampleCount:
            for face in map.faceIter(skipInfinite = True):
                level = 1
                while self._functors[face.label()].pixelCount < minSampleCount and level < 32:
                    level *= 2
                    self.superSample(face, level)
        del self._origSIV # only needed for supersampling here in __init__
        del self._SIV
        self._attachHooks()

    def bands(self):
        return self.originalImage.bands()

    def faceMeanFunctor(self, index):
#         if hasattr(index, "label"):
#             index = index.label()
        return self._functors[index]

    def __getitem__(self, index):
#         if hasattr(index, "label"):
#             index = index.label()
        return self._functors[index].average()

    def superSample(self, face, level = 2):
        if not self._origSIV:
            if self.bands() == 3:
                self._origSIV = sivtools.ThreeBandSIVProxy(self.originalImage, self._SIV)
            elif self.bands() == 2:
                self._origSIV = sivtools.GradientSIVProxy(self.originalImage, self._SIV)
            else:
                self._origSIV = self._SIV(self.originalImage)
        for pos in superSample(face, level):
            self._functors[face.label()](self._origSIV[pos])
        self._superSampled[face.label()] = level

    def preMergeFaces(self, dart):
        ssLeft = self._superSampled[dart.leftFaceLabel()]
        ssRight = self._superSampled[dart.rightFaceLabel()]
        self.ssMerged = 0 # hopefully, the merged face has no supersampling
        if ssLeft and (ssLeft == ssRight):
            self.ssMerged = ssLeft
            ssLeft, ssRight = 0, 0 # can be merged
        if ssLeft and ssRight:
            if ssLeft < ssRight:
                self.ssMerged = ssLeft
                ssLeft = 0 # take stats from face with less supersampling
            else:
                self.ssMerged = ssRight
                ssRight = 0
        if not ssLeft:
            self.mergedStats = copy.copy(self._functors[dart.leftFaceLabel()])
        else:
            self.mergedStats = FaceMeanFunctor(self._defaultValue)
        if not ssRight:
            self.mergedStats.merge(self._functors[dart.rightFaceLabel()])
        return True

    def postMergeFaces(self, survivor):
        self._functors[survivor.label()] = self.mergedStats
        self._superSampled[survivor.label()] = self.ssMerged

    def associatePixels(self, face, positions):
        functor = self._functors[face.label()]
        for pos in positions:
            functor(self.originalImage[pos])

    def regionImage(self, labelImage = None):
        class MeanLookupFunctor(object):
            def __init__(self, faceColorStatistics, default):
                self._fis = faceColorStatistics
                self._default = default

            def __call__(self, label):
                if label >= 0:
                    return self._fis[int(label)]
                return self._default

        if not labelImage:
            labelImage = self._map().labelImage()
        zeroPixel = Pixel(*((0.0, )*self.bands()))
        mlf = MeanLookupFunctor(self, zeroPixel)
        return transformImage(labelImage, mlf)

    def faceMeanDiffGray(self, dart):
        f = self._functors
        m1 = f[dart.leftFaceLabel()].average()
        m2 = f[dart.rightFaceLabel()].average()
        return abs(m1 - m2) / self._diffNorm

    def faceMeanDiffColor(self, dart):
        f = self._functors
        m1 = f[dart.leftFaceLabel()].average()
        m2 = f[dart.rightFaceLabel()].average()
        return (m1 - m2).norm() / self._diffNorm

    def facePoissonLikelyhoodRatio(self, dart):
        f1 = self.faceMeanFunctor(dart.leftFaceLabel())
        f2 = self.faceMeanFunctor(dart.rightFaceLabel())
        m1 = f1.average(); clampEpsilon(m1, 1e-5)
        m2 = f2.average(); clampEpsilon(m2, 1e-5)
        n1 = float(f1.pixelCount)
        n2 = float(f2.pixelCount)
        n = n1 + n2
        m = (n1*m1 + n2*m2)/n
        if self.bands() > 1:
            return sum([n1*m1[i]*math.log(m1[i]) + n2*m2[i]*math.log(m2[i])
                        - n*m[i]*math.log(m[i]) for i in range(self.bands())])
        else:
            return n1*m1*math.log(m1) + n2*m2*math.log(m2) - n*m*math.log(m)

    def faceGaussLikelyhoodRatio(self, dart):
        f1 = self.faceMeanFunctor(dart.leftFaceLabel())
        f2 = self.faceMeanFunctor(dart.rightFaceLabel())
        n1 = float(f1.pixelCount)
        n2 = float(f2.pixelCount)
        n = n1 + n2
        v1 = f1.variance(); clampEpsilon(v1, 1e-5)
        v2 = f2.variance(); clampEpsilon(v2, 1e-5)
        v = ((f1.sum2+f2.sum2) / n) - math.sq((f1.sum+f2.sum) / n)
        clampEpsilon(v, 1e-5)
        if self.bands() > 1:
            return sum([(n*math.log(v[i]) - n1*math.log(v1[i]) - n2*math.log(v2[i]))
                        for i in range(self.bands())]) * 0.5
        else:
            return (n*math.log(v) - n1*math.log(v1) - n2*math.log(v2)) * 0.5

    # FIXME: is this needed?
    def faceStdDevDiff(self, dart):
        f1 = self.faceMeanFunctor(dart.leftFaceLabel())
        f2 = self.faceMeanFunctor(dart.rightFaceLabel())
        return abs(f1.stdDeviation() - f2.stdDeviation())

    def faceHomogeneity2(self, d):
        return math.sq(self.faceMeanDiff(d)) * faceAreaHomogeneity(d)

    def faceHomogeneity(self, d):
        return math.sqrt(self.faceHomogeneity2(d))

FaceGrayStatistics.bands = lambda x: 1
FaceRGBStatistics.bands = lambda x: 3

def FaceColorStatistics(map, originalImage, minSampleCount = 1):
    if originalImage.bands() == 1:
        return FaceGrayStatistics(map, originalImage, minSampleCount)
    elif originalImage.bands() == 3:
        return FaceRGBStatistics(map, originalImage, minSampleCount)
    else:
        return _FaceColorStatistics(map, originalImage, minSampleCount)

def faceAreaHomogeneity(dart):
    a1 = dart.leftFace().area()
    a2 = dart.rightFace().area()
    return (a1 * a2) / (a1 + a2)

# --------------------------------------------------------------------

def _ipqLengths(dart):
    ll = dart.leftFaceLabel()
    rl = dart.rightFaceLabel()

    left = 0.0
    common = 0.0
    right = 0.0

    for d in dart.phiOrbit():
        if d.rightFaceLabel() == rl:
            common += d.edge().length()
        else:
            left += d.edge().length()

    for d in dart.clone().nextAlpha().phiOrbit():
        if d.rightFaceLabel() != ll:
            right += d.edge().length()

    return left, common, right

def mergedContourLength(dart):
    left, common, right = _ipqLengths(dart)
    return (left+right)

def mergedIsoperimetricQuotient(dart):
    left, common, right = _ipqLengths(dart)
    return vigra.sq(left+right) \
           / (4*math.pi*(dart.leftFace().area() + dart.rightFace().area()))

def mergedIsoperimetricQuotient2(dart):
    left, common, right = _ipqLengths(dart)
    after = vigra.sq(left+right) \
            / ((dart.leftFace().area() + dart.rightFace().area()))
    beforeLeft = vigra.sq(left+common) / dart.leftFace().area()
    beforeRight = vigra.sq(right+common) / dart.rightFace().area()
    return after / (beforeLeft + beforeRight)

def seedIsoperimetricQuotient(dart):
    left, common, right = _ipqLengths(dart)
    after = vigra.sq(left+right) \
            / ((dart.leftFace().area() + dart.rightFace().area()))
    if dart.leftFace().flag(flag_constants.SRG_SEED):
        beforeLeft = vigra.sq(left+common) / dart.leftFace().area()
        return after / beforeLeft
    beforeRight = vigra.sq(right+common) / dart.rightFace().area()
    return after / beforeRight

# --------------------------------------------------------------------

class HyperbolicInverse(object):
    """Cost measure for a (e.g. livewire) path, which takes the
    inverse of a darts' removal costs (given a dart cost measure) as
    traversal costs.  To be specific,

       c'(d) := (0.0001 + c(d))^-1 * d.edge().length()"""

    __slots__ = ("measure")

    def __init__(self, singleDartMeasure):
        """initialize with the given edge cost measure"""
        self.measure = singleDartMeasure

    def __call__(self, newDart):
        return newDart.edge().length() / \
               (1e-4 + self.measure(newDart))

# --------------------------------------------------------------------

# FIXME: still old API (see FaceColorStatistics)
class FaceColorHistogram(DynamicFaceStatistics):
    def __init__(self, map, image):
        DynamicFaceStatistics.__init__(self, map)
        self.image = image
        for face in map.faceIter():
            face._colorHistogram = MultiHistogram3()
            face._colorHistogram.init(16,[(0,100),(-86.1813,98.2352),(-107.862,94.4758)])
            face._colorHistogram2 = MultiHistogram2()
            face._colorHistogram2.init(16,[(-86.1813,98.2352),(-107.862,94.4758)])

        class HistInitFunctor(object):
            def __init__(self, map):
                self.map = map

            def __call__(self, label, value):
                if label >= 0:
                    self.map.face(int(label))._colorHistogram.addValue(list(value))
                    self.map.face(int(label))._colorHistogram2.addValue(list(value)[1:])

        inspectImage(map.labelImage(), image, HistInitFunctor(map))

        for face in map.faceIter():
            face._colorHistogram.gaussianSmoothing(2.2)
            face._colorHistogram2.gaussianSmoothing(2.2)
        self._attachHooks()

    def preMergeFaces(self, dart):
        self.mergedHist = dart.leftFace()._colorHistogram
        self.mergedHist.add(dart.rightFace()._colorHistogram)
        self.mergedHist2 = dart.leftFace()._colorHistogram2
        self.mergedHist2.add(dart.rightFace()._colorHistogram2)
        return True

    def postMergeFaces(self, survivor):
        survivor._colorHistogram = self.mergedHist
        survivor._colorHistogram2 = self.mergedHist2

    def associatePixels(self, face, positions):
        for pos in positions:
            face._colorHistogram.addValue(list(self.image[pos]))
            face._colorHistogram2.addValue(list(self.image[pos])[1:])

    def faceHistDiff(dart):
        if (dart.leftFace()._colorHistogram.count()==0) or (dart.rightFace()._colorHistogram.count()==0):
            return 0
        return dart.leftFace()._colorHistogram.diff(dart.rightFace()._colorHistogram)

    def faceHist2Diff(dart):
        if (dart.leftFace()._colorHistogram2.count()==0) or (dart.rightFace()._colorHistogram2.count()==0):
            return 0
        return dart.leftFace()._colorHistogram2.diff(dart.rightFace()._colorHistogram2)

# --------------------------------------------------------------------
#               Edge-based Statistics & Cost Measures
# --------------------------------------------------------------------

class DynamicEdgeStatistics(DetachableStatistics):
    __slots__ = []
    
    def _attachHooks(self):
        self._attachedHooks = (self._map().addMergeEdgesCallbacks(
            self.preMergeEdges, self.postMergeEdges), )

def mergedEdgeCostWeightedByLength(edge1, cost1, edge2, cost2):
    """Combiner function for `StaticEdgeCosts`.  Defines the cost for
    a merged edge as the average of the costs from the merged edges,
    weighted by their lengths."""
    
    l1 = edge1.length()
    l2 = edge2.length()
    return (cost1*l1 + cost2*l2)/(l1+l2)

class StaticEdgeCosts(DetachableStatistics):
    """Initially computes and stores costs for all edges, then serves
    as a cost measure that always returns the initial costs.

    In order to define costs for merged edges, one can pass a
    `combiner` to the constructor.  If this is None (default),
    mergeEdges operations are blocked in order to guarantee sensible
    results (i.e. passing the problem to the calling algorithm).
    If a combiner is given, it is called with four parameters::
    
      mergedCost = combiner(edge1, cost1, edge2, cost2)

    and must return a combined cost.  An example of a suitable
    `combiner` is `mergedEdgeCostWeightedByLength`."""

    __slots__ = ["_costs", "_combiner",
                 "_mergedCost"]

    def __init__(self, map, costMeasure = None, costs = None,
                 skipBorder = True, combiner = None):
        DetachableStatistics.__init__(self, map)
        assert costs is not None or costMeasure, "need static costs either as costs-array or indirectly as a costMeasure"
        if costs is None:
            self._costs = [None] * map.maxEdgeLabel()
        else:
            if combiner:
                # let _costs be a copy (we don't want to modify the
                # caller's data structure with the combined values):
                self._costs = list(costs)
            else:
                self._costs = costs
        self._combiner = combiner
        if costMeasure:
            for edge in map.edgeIter():
                if skipBorder and edge.flag(flag_constants.BORDER_PROTECTION):
                    continue
                self._costs[edge.label()] = costMeasure(edge.dart())
        self._attachHooks()

    def _attachHooks(self):
        if self._combiner:
            self._attachedHooks = (self._map().addMergeEdgesCallbacks(
                self.preMergeEdges, self.postMergeEdges), )
        else:
            self._attachedHooks = (self._map().addMergeEdgesCallbacks(
                self.blockMergeEdges, None), )

    def preMergeEdges(self, dart):
        edge1 = dart.edge()
        edge2 = dart.clone().nextSigma().edge()
        try:
            self._mergedCost = self._combiner(
                edge1, self._costs[edge1.label()],
                edge2, self._costs[edge2.label()])
            return True
        except TypeError: # i.e. combiner cannot handle "None" costs
            return False
    
    def postMergeEdges(self, survivor):
        self._costs[survivor.label()] = self._mergedCost

    def blockMergeEdges(self, dart):
        return False
    
    def __call__(self, dart):
        return self._costs[dart.edgeLabel()]

    def __getstate__(self):
        return DynamicEdgeStatistics.__getstate__(self) + (
            self._costs, )

    def __setstate__(self, (map, tree)):
        DynamicEdgeStatistics.__setstate__(self, (map, ))
        self._costs = tree

# formerly, this was called EdgeMergeTree:
class MergedEdges(DynamicEdgeStatistics):
    """This class manages a list of edges (from the moment of
    attaching, usually from level 0) that have been merged into each
    edge."""

    __slots__ = ["_labelLUT",
                 "_merged"]
    
    def __init__(self, map):
        DynamicEdgeStatistics.__init__(self, map)
        self._labelLUT = geomap.LabelLUT(map.maxEdgeLabel())
        self._attachHooks()
    
    def preMergeEdges(self, dart):
        self._merged = dart.clone().nextSigma().edgeLabel()
        return True
    
    def postMergeEdges(self, survivor):
        self._labelLUT.relabel(self._merged, survivor.label())
    
    def __getitem__(self, edge):
        """Returns list of all edges of original map of which the
        given edge of the current map is composed."""
        if hasattr(edge, "label"):
            edge = edge.label()
        return list(self._labelLUT.merged(edge))

    def __getstate__(self):
        return DynamicEdgeStatistics.__getstate__(self) + (
            self._labelLUT, )

    def __setstate__(self, (map, tree)):
        DynamicEdgeStatistics.__setstate__(self, (map, ))
        self._labelLUT = tree

    def level0Labels(self):
        """Return a list of level0 edges (resp. their labels) the
        currently existing edges are composed of.  (That is, sort of a
        list of surviving boundaries; labels of merged edges appear,
        but removed do not.)"""
        result = []
        for edge in self._map().edgeIter():
            result.extend(self[edge])
        return result

class DynamicEdgeIndices(DetachableStatistics):
    __base = DetachableStatistics
    __slots__ = ["_indices",
                 "_mergedIndices", "_newIndices1", "_newIndices2", "_mergeLabels"]
    
    def __init__(self, map):
        DetachableStatistics.__init__(self, map)
        self._indices = [None] * map.maxEdgeLabel()
        for edge in map.edgeIter():
            self._indices[edge.label()] = []
    
    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeEdgesCallbacks(self.preMergeEdges,
                                               self.postMergeEdges),
            self._map().addSplitEdgeCallbacks(self.preSplitEdge,
                                              self.postSplitEdge))

    def __getstate__(self):
        return self.__base.__getstate__(self) + (self._indices, )

    def __setstate__(self, (map, indices)):
        self.__base.__setstate__(self, (map, ))
        self._indices = indices

    def edgeIndices(self, edge):
        if hasattr(edge, "label"):
            edge = edge.label()
        return list(self._indices[edge])

    def dartIndices(self, dart):
        edgeIndices = self._indices[dart.edgeLabel()]
        if dart.label() > 0:
            return list(edgeIndices)
        mi = len(dart)-1
        return [mi-i for i in edgeIndices]

    def preMergeEdges(self, dart):
        # dart belongs to the surviving edge and starts at the merged
        # node:
        
        edge1 = dart.edge()
        edge2 = dart.clone().nextSigma().edge()
        self._mergeLabels = (edge1.label(), edge2.label())

        if dart.label() < 0:
            # dart.edge() will be simply extend()ed
            el1 = len(edge1) - 1
            self._mergedIndices = self.dartIndices(dart.clone().nextAlpha())
            self._mergedIndices.extend([
                el1 + i for i in self.dartIndices(dart.clone().nextSigma())
                if i])
        else:
            # edge2 will be inserted before dart.edge():
            el2 = len(edge2) - 1
            self._mergedIndices = self.dartIndices(
                dart.clone().nextSigma().nextAlpha())
            self._mergedIndices.extend([
                el2 + i for i in self.dartIndices(dart) if i])

        return True

    def postMergeEdges(self, survivor):
        self._indices[survivor.label()] = self._mergedIndices
        if self._mergeLabels[0] == survivor.label():
            self._indices[self._mergeLabels[1]] = None
        else:
            self._indices[self._mergeLabels[0]] = None

    def preSplitEdge(self, edge, segmentIndex, newPoint):
        oldIndices = self._indices[edge.label()]

        self._newIndices1 = [i for i in oldIndices if i <= segmentIndex]

        offset = -segmentIndex
        if newPoint:
            segmentIndex += 1
        self._newIndices2 = [i+offset for i in oldIndices if i >= segmentIndex]

        return True

    def postSplitEdge(self, oldEdge, newEdge):
        self._indices[oldEdge.label()] = self._newIndices1
        assert len(self._indices) == newEdge.label()
        self._indices.append(self._newIndices2)

class WatershedStatistics(DynamicEdgeIndices):
    __base = DynamicEdgeIndices
    __slots__ = ["_passValues", "_indices", "_gmSiv",
                 "_basinDepth",
                 "_mergedPV"]
    
    def __init__(self, map, flowlines, gmSiv):
        self.__base.__init__(self, map)
        self._passValues = [None] * map.maxEdgeLabel()
        for edge in map.edgeIter():
            if edge.flag(flag_constants.BORDER_PROTECTION):
                continue
            
            flowline = flowlines[edge.label()]
            saddleIndex = flowline[3]

            # flowline tracing might have stopped, in which case the
            # polygon can be modified by connecting to the nearest
            # maximum/node, possibly adding a point and thus shifting
            # the indices:
            if len(edge) != len(flowline[2]) and edge[1] == flowline[2][0]:
                assert flowline[0] <= 0
                saddleIndex += 1
            
            self._passValues[edge.label()] = gmSiv[edge[saddleIndex]]
            self._indices[edge.label()].append(saddleIndex)

        self._gmSiv = gmSiv
        self._basinDepth = None
        self._attachHooks()

    def __getstate__(self):
        result = self.__base.__getstate__(self) + (
            self._passValues, self._basinDepth)
        return result

    def __setstate__(self, s):
        map, indices, passValues = s[:3]
        self.__base.__setstate__(self, (map, indices))
        self._passValues = passValues
        if len(s) > 3:
            self._basinDepth = s[3]
        else:
            self._basinDepth = None
        self._gmSiv = None

    def nonWatershedEdgesAdded(self):
        """wsStats.nonWatershedEdgesAdded()

        Call e.g. after connectBorderNodes, when new non-watershed
        edges have been added."""

        mel = self._map().maxEdgeLabel()
        if len(self._passValues) < mel:
            self._passValues.extend([None] * (mel - len(self._passValues)))
            self._indices.extend([[] for i in
                                  range(mel - len(self._passValues))])

    def preMergeEdges(self, dart):
        if not self.__base.preMergeEdges(self, dart):
            return False
        
        edge1 = dart.edge()
        edge2 = dart.clone().nextSigma().edge()

        self._mergedPV = None

        # do not allow merging of watersheds with border:
        if edge1.flag(flag_constants.BORDER_PROTECTION):
            return edge2.flag(flag_constants.BORDER_PROTECTION)
        if edge2.flag(flag_constants.BORDER_PROTECTION):
            return edge1.flag(flag_constants.BORDER_PROTECTION)

        self._mergedPV = min(self._passValues[edge1.label()],
                             self._passValues[edge2.label()])

        return True

    def postMergeEdges(self, survivor):
        self.__base.postMergeEdges(self, survivor)

        self._passValues[survivor.label()] = self._mergedPV

        if self._mergeLabels[0] == survivor.label():
            self._passValues[self._mergeLabels[1]] = None
        else:
            self._passValues[self._mergeLabels[0]] = None

    def _resetPassValue(self, edge):
        # SIV cannot be pickled:
        assert self._gmSiv, "cannot find correct passvalue after edge splitting without SplineImageView!"
        
        ind = self._indices[edge.label()]
        if ind:
            self._passValues[edge.label()] = min([
                self._gmSiv[edge[i]] for i in ind])
        else:
            self._passValues[edge.label()] = min(
                self._gmSiv[edge[0]], self._gmSiv[edge[-1]])

    def postSplitEdge(self, oldEdge, newEdge):
        DynamicEdgeIndices.postSplitEdge(self, oldEdge, newEdge)
        
        self._resetPassValue(oldEdge)
        self._passValues.append(None)
        self._resetPassValue(newEdge)

    def edgeSaddles(self, edge):
        if hasattr(edge, "label"):
            edgeLabel = edge.label()
        else:
            edgeLabel = edge
            edge = self._map().edge(edgeLabel)
        return [edge[i] for i in self._indices[edgeLabel]]

    def dartPassValue(self, dartOrEdge):
        return self._passValues[abs(dartOrEdge.label())]

    def passValue(self, dart):
        return min([self.dartPassValue(d)
                    for d in commonBoundaryDarts(dart)])

    def setBasins(self, basinStatistics):
        self._basinDepth = basinStatistics._basinDepth

    def dynamics(self, edge):
        """At the moment, `edge` may be an Edge, a Dart, or a label,
        but don't rely on this!  Since this is a cost measure, it may
        only work for Dart objects in the future."""
        if hasattr(edge, "label"):
            edgeLabel = abs(edge.label())
        else:
            edgeLabel = abs(edge)
            edge = self._map().edge(edgeLabel)
        return self._passValues[edgeLabel] - min(
            self._basinDepth[edge.leftFaceLabel()],
            self._basinDepth[edge.rightFaceLabel()])

class WatershedBasinStatistics(DetachableStatistics):
    __base = DetachableStatistics
    __slots__ = ["_basinDepth",
                 "_mergedDepth"]
    
    def __init__(self, map, minima, gmSiv):
        self.__base.__init__(self, map)

        self._basinDepth = [None] * map.maxFaceLabel()
        for mpos in minima:
            face = map.faceAt(mpos)
            if face.label() == 0:
                continue
            depth = gmSiv[mpos]
            fDepth = self._basinDepth[face.label()]

            # multiple minima may happen "legally" through saddle
            # filtering (by proximity / saddleThreshold)!
            # (also, for face 0 without border closing!)
#             if fDepth != None:
#                 sys.stderr.write(
#                     "Face %d (area %s) contains more than one minimum!\n" % (
#                     face.label(), face.area()))

            if fDepth == None or fDepth > depth:
                self._basinDepth[face.label()] = depth

        import maputils

        for face in map.faceIter(skipInfinite = True):
            if self._basinDepth[face.label()] == None:
                isAtBorder = False
                for dart in maputils.neighborFaces(face):
                    if dart.rightFaceLabel() == 0:
                        isAtBorder = True
                        break
                if not isAtBorder:
                    sys.stderr.write(
                        "Face %d (area %s, anchor %d) contains no minimum!\n" % (
                        face.label(), face.area(), face.contour().label()))
                level = 2
                while level < 20: # prevent endless loop
                    level += 1
                    samples = [gmSiv[pos] for pos in superSample(face, level)]
                    if len(samples) > 10:
                        break
                if not isAtBorder:
                    sys.stderr.write(
                        "  (got %d samples at supersampling level %d)\n" % (
                        len(samples), level))
                if samples:
                    self._basinDepth[face.label()] = min(samples)
                    assert self._basinDepth[face.label()] != None, str(samples)
                else:
                    sys.stderr.write("  *** BIG FAT WARNING: faking basin depth! ***\n")
                    self._basinDepth[face.label()] = 0.0

        self._attachHooks()

    def superSample(self, face, level = 2):
        if not self._origSIV:
            if self.bands() == 3:
                self._origSIV = sivtools.ThreeBandSIVProxy(self.originalImage, self._SIV)
            else:
                self._origSIV = self._SIV(self.originalImage)
        for pos in superSample(face, level):
            self._functors[face.label()](self._origSIV[pos])
        self._superSampled[face.label()] = level

    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeFacesCallbacks(self.preMergeFaces,
                                               self.postMergeFaces), )

    def preMergeFaces(self, dart):
        self._mergedDepth = min(
            self._basinDepth[dart.leftFaceLabel()],
            self._basinDepth[dart.rightFaceLabel()])
        return True

    def postMergeFaces(self, survivor):
        self._basinDepth[survivor.label()] = self._mergedDepth

    def __getstate__(self):
        return self.__base.__getstate__(self) + (
            self._basinDepth, )

    def __setstate__(self, (map, basinDepth)):
        self.__base.__setstate__(self, (map, ))
        self._basinDepth = basinDepth

def _makeAttrName(someStr):
    attrTrans = string.maketrans(".+-", "__n")
    return someStr.translate(attrTrans)

class BoundaryIndicatorStatistics(DynamicEdgeStatistics):
    __base = DynamicEdgeStatistics
    __slots__ = ["_functors",
                 "_mergedStats"]
    
    def __init__(self, map):
        DynamicEdgeStatistics.__init__(self, map)
        self._functors = [None] * map.maxEdgeLabel()

    def __getstate__(self):
        return self.__base.__getstate__(self) + (self._functors, )

    def __setstate__(self, (map, functors)):
        self.__base.__setstate__(self, (map, ))
        self._functors = functors

    def preMergeEdges(self, dart):
        self._mergedStats = copy.copy(self._functors[dart.edgeLabel()])
        self._mergedStats.merge(
            self._functors[dart.clone().nextSigma().edgeLabel()])
        return True

    def __getitem__(self, edgeLabel):
        """Return the functor for the given edge label."""
        return self._functors[edgeLabel]

    def postMergeEdges(self, survivor):
        self._functors[survivor.label()] = self._mergedStats

def commonBoundaryDarts(dart):
    rightFaceLabel = dart.rightFaceLabel()
    for d in dart.phiOrbit():
        if d.rightFaceLabel() == rightFaceLabel:
            yield d

class EdgeGradientStatistics(BoundaryIndicatorStatistics):
    """Analyzes e.g. the gradient magnitude on the edge.

    To be specific, manages statistics for each edge on the values
    sampled from a SplineImageView (cf. constructor) at the edge's
    support points.  Each edge has an associated object which manages
    the average value (and possibly a list of sorted values for
    quantile queries).  This functor object can be of any
    user-specified type (see constructor)."""

    __base = BoundaryIndicatorStatistics
    __slots__ = ["_Functor"]
    
    def __init__(self, map, gradSiv, resample = 0.1, tangents = None,
                 Functor = geomap.PolylineStatistics):
        """Init statistics for each edge.  Possible values for Functor
        are: geomap.PolylineStatistics (default) and
        geomap.QuantileStatistics (needed if you want to use
        `quantile` or `dartQuantile`).  The latter needs much more
        memory, that's why it's not the default.

        If `tangents` are given, this replaces the former
        EdgeGradDirDotStatistics.  I.e. instead of sampling a simple
        SplineImageView (gradSiv) at the polyline points, the
        `gradSiv` parameter is assumed to give gradient *vector*
        output (cf. GradientSIVProxy) that is sampled at the tangent
        positions and multiplied with tangent unit vectors.

        Note that the values are not weighted with the polyline
        segment lengths if tangents are given, since usually the
        weighting already happened during the tangents computation
        (and the arclengths between tangent samples would be
        inappropriate here)."""
        
        BoundaryIndicatorStatistics.__init__(self, map)
        self._Functor = Functor

        if tangents and resample:
            sys.stderr.write("EdgeGradientStatistics: cannot resample tangents. resample argument (%s) ignored.\n" % (
                resample == 0.1 and "default = 0.1" or resample, ))
        
        for edge in map.edgeIter():
            if not tangents:
                poly = resample and resamplePolygon(edge, resample) or edge
                self._functors[edge.label()] = Functor(poly, gradSiv)
            else:
                stats = Functor()

                dp = geomap.DartPosition(edge.dart())
                for al, theta in tangents[edge.label()]:
                    dp.gotoArcLength(al)

                    gradDir = gradSiv[dp()]
                    #gradDir /= gradDir.magnitude()

                    segment = Vector2(-math.sin(theta), math.cos(theta))

                    stats(dot(gradDir, segment), 1.0)

                self._functors[edge.label()] = stats

        self._attachHooks()

    def __getstate__(self):
        return self.__base.__getstate__(self) + (self._Functor, )

    def __setstate__(self, state):
        self.__base.__setstate__(self, state[:-1])
        self._Functor = state[-1]

    def dartMin(self, dart):
        """Returns the minimum gradient on the edge."""
        return self[dart.edgeLabel()].min()

    def dartMax(self, dart):
        """Returns the maximum gradient on the edge."""
        return self[dart.edgeLabel()].max()

    def dartQuantile(self, q = 0.5):
        """Returns a specific quantile of the sampled gradients."""
        def specificQuantile(dart):
            return self[dart.edgeLabel()].quantile(q)
        return specificQuantile

    def dartAverage(self, dart):
        """Return the average gradient on the edge (properly weighted
        by polyline segment length)."""
        return self[dart.edgeLabel()].average()

    def combinedStatistics(self, dart):
        """Mostly internal; return a functor for all common edges of
        the faces adjacent to `dart`."""
        result = self._Functor()
        for d in commonBoundaryDarts(dart):
            result.merge(self[d.edgeLabel()])
        return result

    def min(self, dart):
        """Returns the minimum gradient on this common contour."""
        return self.combinedStatistics(dart).min()

    def max(self, dart):
        """Returns the maximum gradient on this common contour."""
        return self.combinedStatistics(dart).max()

    def quantile(self, q):
        """Returns a specific quantile of the sampled gradients on the
        common contour."""
        def specificQuantile(dart):
            return self.combinedStatistics(dart).quantile(q)
        return specificQuantile

    def supportsQuantile(self):
        """Returns True iff the internally used statistics functor supports quantiles"""
        return hasattr(self._Functor, "quantile")

    def average(self, dart):
        """Return the average gradient on this common contour
        (properly weighted by polyline segment length)."""
        return self.combinedStatistics(dart).average()

# USAGE:
# >>> boundaryIndicator = SplineImageView5(...)
# >>> egs = EdgeGradientStatistics(map, boundaryIndicator)
# >>> egs[100].average()
# 1.178475851870056
# >>> egs[100].quantile(0.4)
# 1.1785515546798706

# --------------------------------------------------------------------

def scaleArcLengths(arcLengthList, factor):
    return [(al*factor, v) for al, v in arcLengthList]

import bisect

class PercentPointFunction(list):
    """This represents the PPF (inverse CDF) of a piecewise linear
    function, and can be used for computing quantiles.  Furthermore,
    it can compute its inverse (the CDF) and the PDF.  It is intended
    for arcLengthLists of boundary indicator values."""

    __slots__ = ("valueScale")

    def __init__(self, arcLengthList):
        list.__init__(self)
        if not arcLengthList:
            return

        ss = self.segments(self.splitVert(arcLengthList))
        ss.sort()

        al = 0.0
        v1 = ss[0][0]
        iCDF = [(al, v1)]
        cur = Vector2(v1, 0)
        for v1, l, v2 in ss:
            if v1 > cur[0]:
                al += cur[1]
                iCDF.append((al, v1))
                cur = Vector2(v1, l)
            else:
                cur[1] += l
        if cur[1]:
            al += cur[1]
            iCDF.append((al, v2))

        self.valueScale = iCDF[-1][0]
        self.extend(scaleArcLengths(iCDF, 1.0/self.valueScale))

    @staticmethod
    def splitVert(arcLengthList):
        values = [v for al, v in arcLengthList]
        
        result = list(arcLengthList) # copy
        prevV = None
        for v in values:
            if v == prevV:
                continue
            prevV = v
            i = 1
            while i < len(result):
                if result[i-1][1] < v < result[i][1] or result[i-1][1] > v > result[i][1]:
                    al1, v1 = result[i-1]
                    al2, v2 = result[i]
                    splitAl = al1+(al2-al1)*(v-v1)/(v2-v1)
                    result.insert(i, (splitAl, v))
                i += 1
        
        return result

    @staticmethod
    def segments(arcLengthList):
        """Returns a list of (v1, width, v2) for each piecewise linear
        segment in the given arcLengthList (i.e. length(result) ==
        len(arcLengthList)-1)."""
        
        result = []
        for i in range(1, len(arcLengthList)):
            al1, v1 = arcLengthList[i-1]
            al2, v2 = arcLengthList[i]
            if v2 < v1:
                v1, v2 = v2, v1
            result.append((v1, (al2-al1), v2))
        return result

    def unscaled(self):
        return scaleArcLengths(self, self.valueScale)

    def quantile(self, p = 0.5):
        # FIXME: use bisect, too
        for i in range(1, len(self)):
            if self[i][0] >= p:
                al1, v1 = self[i-1]
                al2, v2 = self[i]
                return v1+(v2-v1)*(p-al1)/(al2-al1)

    def cdf(self):
        return [(b,a) for a,b in self]

    def pdf(self):
        return ProbabiliyDensityFunction(self.cdf())

class ProbabiliyDensityFunction(object):
    """This represents the PDF of a piecewise linear CDF, and allows
    the computation of the probability density of a given continuous
    value."""

    def __init__(self, cdf):
        self.pdf = [(cdf[0][0], 0.0)]
        for i in range(len(cdf)-1):
            v1, cp1 = cdf[i]
            v2, cp2 = cdf[i+1]
            p = (cp2-cp1)/(v2-v1)
            self.pdf.append((v2, p))

    def __call__(self, value):
        i = bisect.bisect_left(self.pdf, (value, 1.0))
        if i == len(self.pdf):
            return 0.0
        return self.pdf[i][1]

    def steps(self):
        """Return step-function, e.g. for plotting."""
        result = []
        for i in range(len(self.pdf)-1):
            v1, p1 = self.pdf[i]
            v2, p2 = self.pdf[i+1]
            result.append((v1, p1))
            result.append((v1, p2))
        result.append((v2, p2))
        result.append((v2, 0.0))
        return result

# --------------------------------------------------------------------

def calcGradAngDisp(grad,n):
    ki = GrayImage(2*n+1,2*n+1)
    for x in range(0,ki.width()):
        for y in range(0,ki.height()):
            ki[x,y] = 1.0
    kernel = customKernel(ki,(n,n))
    mi = convolveImage(transformImage(grad, '\l x: sqrt(sq(x[0])+sq(x[1]))'),kernel)
    gi = convolveImage(grad,kernel)
    gad = GrayImage(grad.size())
    for x in range(0,gad.width()):
        for y in range(0,gad.height()):
            if (mi[x,y]>0):
                gad[x,y] = math.hypot(gi[x,y][0],gi[x,y][1])/mi[x,y]

    return gad

class EdgeGradAngDispStatistics(BoundaryIndicatorStatistics):
    def __init__(self, map, bi, n, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, ("gad_%s_" % n) + bi.name)
        gad = calcGradAngDisp(bi.grad, n)
        gad.siv = eval("SplineImageView%d(gad)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, gad.siv))
        self._attachHooks()

class EdgeMinimumDistance(DynamicEdgeStatistics):
    def __init__(self, map, minima):
        self.attrName = "minDist"
        minimaMap = geomap.PositionedMap()
        for min in minima:
            minimaMap.insert(min, min)
        for edge in map.edgeIter():
            mindist2 = 1e16
            for p in edge:
                near = minimaMap(p, mindist2)
                if near:
                    mindist2 = (near-p).squaredMagnitude()
            setattr(edge, self.attrName, math.sqrt(mindist2))
        self._attachHooks()

    def preMergeEdges(self, dart):
        self.mergedStats = min(getattr(dart.edge(), self.attrName),
                               getattr(dart.clone().nextSigma().edge(), self.attrName))
        return True

    def postMergeEdges(self, survivor):
        setattr(survivor, self.attrName, self.mergedStats)

def calcGradScaleSum(image, steps):
    gss = GrayImage(image.size())
    scale = 0.7
    for i in range(steps*2+1):
        ti = vectorToTensor(gaussianGradientAsVector(image.subImage(0),scale))
        if (image.bands()>1):
            for j in range(1,image.bands()):
                ti += vectorToTensor(gaussianGradientAsVector(image.subImage(j),scale))
        gm = transformImage(tensorTrace(ti),'\l x:sqrt(x)')
        mm = MinMax()
        inspectImage(gm,mm)
        gm = linearRangeMapping(gm,oldRange=(0,mm.max()),newRange=(0,1.0))
        gss += gm
        scale *= 1.41421
    return gss

class EdgeGradScaleSum(BoundaryIndicatorStatistics):
    def __init__(self, map, image, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "gradScaleSum")
        gms = calcGradScaleSum(image,4)
        avgFunc=Average()
        inspectImage(gms,avgFunc)
        avg=avgFunc.average()
        gms/=avg[0]
        gms.siv=eval("SplineImageView%d(gms)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, gms.siv))
        self._attachHooks()

class EdgeSpatialStability(BoundaryIndicatorStatistics):
    def __init__(self, map, image, radius, propexp, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "spatialStability_%s_%s" % (radius, propexp))
        ss = geomap.spatialStabilityImage(image, 4, radius, propexp)
        ss.siv = eval("SplineImageView%d(ss)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, ss.siv))
        self._attachHooks()

def calcGradProd(image, sigma1, sigma2):
    ggv1 = gaussianGradientAsVector(image,sigma1)
    ggv2 = gaussianGradientAsVector(image,sigma2)
    gp = transformImage(ggv1,ggv2,"\l g1,g2: norm(Vector(sqrt(max(g1[0]*g2[0],0)),sqrt(max(g1[1]*g2[1],0))))")
    return gp

def calcColGradProd(image, sigma1, sigma2):
    bandTensors1 = [vectorToTensor(gaussianGradientAsVector(image.subImage(i), sigma1)) for i in range(3)]
    bandTensors2 = [vectorToTensor(gaussianGradientAsVector(image.subImage(i), sigma2)) for i in range(3)]
    colorTensor1 = bandTensors1[0] + bandTensors1[1] + bandTensors1[2]
    colorTensor2 = bandTensors2[0] + bandTensors2[1] + bandTensors2[2]
    cgm1 = tensorTrace(colorTensor1)
    cgm2 = tensorTrace(colorTensor2)
    cg1 = transformImage(tensorEigenRepresentation(colorTensor1), cgm1,
        "\l e, mag: Vector(cos(-e[2])*mag, sin(-e[2])*mag)", {})
    cg2 = transformImage(tensorEigenRepresentation(colorTensor2), cgm2,
        "\l e, mag: Vector(cos(-e[2])*mag, sin(-e[2])*mag)", {})
    cgp = transformImage(cg1,cg2,"\l g1,g2: norm(Vector(sqrt(max(g1[0]*g2[0],0)),sqrt(max(g1[1]*g2[1],0))))")
    return cgp

class EdgeGradProd(BoundaryIndicatorStatistics):
    def __init__(self, map, image, sigma1, sigma2, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "gradProd_%s_%s" % (sigma1, sigma2))
        gp = calcGradProd(image, sigma1, sigma2)
        gp.siv = eval("SplineImageView%d(gp)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, gp.siv))
        self._attachHooks()

class EdgeColGradProd(BoundaryIndicatorStatistics):
    def __init__(self, map, image, sigma1, sigma2, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "colGradProd_%s_%s" % (sigma1, sigma2))
        cgp = calcColGradProd(image, sigma1, sigma2)
        cgp.siv = eval("SplineImageView%d(cgp)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, cgp.siv))
        self._attachHooks()

def calcBTPhaseImage(image,scale):
    filterImage=getBTFilterResponses(image,scale)
    dirImage=transformImage(tensorEigenRepresentation(transformImage(filterImage,"\l x: Vector(sq(x[0])+sq(x[1])+sq(x[3]+x[5]),-x[1]*(x[0]+x[2])-(x[3]+x[5])*(x[4]+x[6]),sq(x[1])+sq(x[2])+sq(x[4]+x[6]))")),"\l x: Vector(cos(-x[2]), sin(-x[2]))")
    quadImage=transformImage(filterImage, dirImage,"\l x, y: Vector(sq(y[0])*x[0]+2.0*y[0]*y[1]*x[1]+sq(y[1])*x[2],y[0]*(x[3]+x[5])+y[1]*(x[4]+x[6]))")
#    phaseImage=transformImage(quadImage,"\l x: abs(abs(atan2(x[1],x[0]))-1.5708)")
    phaseImage=transformImage(quadImage,"\l x: atan2(x[1],x[0])")
    return phaseImage

class EdgePhase(DynamicEdgeStatistics):
    def __init__(self, map, image, scaleList, splineOrder):
        self.scaleList=scaleList
        self.attrNames=[]
        for s in scaleList:
          attrName = _makeAttrName("phase_%s" % (s))
          self.attrNames.append(attrName)
        phaseImages=[calcBTPhaseImage(image,scaleList[i]) for i in range(len(scaleList))]
        phaseDiffImages=[transformImage(phaseImages[i],phaseImages[i+1],"\l x1,x2: 3.1416-abs(abs(x1-x2)-3.1416)") for i in range(len(scaleList)-1)]
        phaseDiffSum=GrayImage(image.size())
        for i in range(len(phaseDiffImages)):
            phaseDiffSum+=phaseDiffImages[i]/abs(scaleList[i+1]-scaleList[i])
        phaseDiffSum/=-len(phaseDiffImages)
        for i in range(len(phaseImages)):
            phaseImages[i]=transformImage(phaseImages[i],"\l x: -abs(abs(x)-1.5708)")
        phaseSivs=[eval("SplineImageView%d(phaseImages[i])" % (splineOrder, )) for i in range(len(scaleList))]
        phaseDiffSumSiv=eval("SplineImageView%d(phaseDiffSum)" % (splineOrder, ))
        for edge in map.edgeIter():
            for i in range(len(scaleList)):
                setattr(edge, self.attrNames[i], EdgeStatistics(edge,phaseSivs[i]))
            edge.phaseDiffSum = EdgeStatistics(edge,phaseDiffSumSiv)
        self._attachHooks()

    def preMergeEdges(self, dart):
        self.mergedStats=[]
        for i in range(len(self.scaleList)):
            stats = copy.copy(getattr(dart.edge(), self.attrNames[i]))
            stats.merge(getattr(dart.clone().nextSigma().edge(), self.attrNames[i]))
            self.mergedStats.append(stats)
        self.mergedStatsDiff = copy.copy(dart.edge().phaseDiffSum)
        self.mergedStatsDiff.merge(dart.clone().nextSigma().edge().phaseDiffSum)
        return True

    def postMergeEdges(self, survivor):
        for i in range(len(self.scaleList)):
            setattr(survivor, self.attrNames[i], self.mergedStats[i])
        survivor.phaseDiffSum = self.mergedStatsDiff

def contAngle(d1,d2,length):
    pl1=list(d1)
    pl2=list(d2)
    l=0
    p1=pl1[1]
    for i in range(2,len(pl1)):
        l+=(pl1[i]-pl1[i-1]).norm()
        p1=pl1[i]
        if l>length:
            break
    l=0
    p2=pl2[1]
    for i in range(2,len(pl2)):
        l+=(pl2[i]-pl2[i-1]).norm()
        p2=pl2[i]
        if l>length:
            break
    v1=(p1-pl1[0])
    v1/=v1.norm()
    v2=(p2-pl2[0])
    v2/=v2.norm()
    return abs(v1[0]*v2[0]+v1[1]*v2[1])

class EdgeContAngle(DetachableStatistics):
    def __init__(self, map, length):
        self.length=length
        self.attrName = _makeAttrName("contAngle_%s" % (self.length))
        for edge in map.edgeIter():
            self.calcContAngle(edge)
        
        self._attachedHooks = (map.addMergeFacesCallbacks(
            None, self.postMergeEdges), )

    def calcContAngle(self, edge):
        d1 = edge.dart()
        d2 = edge.dart()
        d2.nextSigma()
        l1=[]
        while not d1.label() == d2.label():
            l1.append(contAngle(d1,d2,self.length))
            d2.nextSigma()
        d1.nextAlpha()
        d2.nextAlpha()
        d2.nextSigma()
        l2=[]
        while not d1.label() == d2.label():
            l2.append(contAngle(d1,d2,self.length))
            d2.nextSigma()
        if len(l1)==0:
            a1=0
        else:
            a1=max(l1)
        if len(l2)==0:
            a2=0
        else:
            a2=max(l2)
        setattr(edge, self.attrName, min(a1,a2))

    def postMergeEdges(self, survivor):
        self.calcContAngle(survivor)

from geomap import fitLine

class EdgeCurvChange(DetachableStatistics):
    def __init__(self, map):
        for edge in map.edgeIter():
            self.calcCurvChange(edge)
        self._attachedHooks = (map.addMergeFacesCallbacks(
            None, self.postMergeEdges), )

    def calcCurvChange(self, edge):
        # use residuum of line fit on tangents as curvature change definition:
        # (FIXME: well-defined residuum? cf. Svens comments... ;*) )
        edge.curvChange = fitLine(edge.tangents)[2]

    def postMergeEdges(self, survivor):
        self.calcCurvChange(survivor)

from geomap import ParabolaFit

class EdgeCurvChangeLin(DetachableStatistics):
    def __init__(self, map):
        for edge in map.edgeIter():
            self.calcCurvChangeLin(edge)
        self._attachedHooks = (map.addMergeFacesCallbacks(
            None, self.postMergeEdges), )

    def calcCurvChangeLin(self, edge):
        # use residuum of parabola fit on tangents:
        # (FIXME: well-defined residuum? cf. Svens comments... ;*) )
        fit = ParabolaFit()
        fit.addTangentList(edge.tangents)
        edge.curvChangeLin = fit.sumOfSquaredErrors() / len(edge.tangents)

    def postMergeEdges(self, survivor):
        self.calcCurvChangeLin(survivor)

class EdgeRegularity(DetachableStatistics):
    def __init__(self, map, seglength):
        self.seglength = seglength
        self.attrName = _makeAttrName("regularity_%s" % (self.seglength))
        for edge in map.edgeIter():
            self.calcRegularity(edge)

        self._attachedHooks = (map.addMergeFacesCallbacks(
            None, self.postMergeEdges), )

    def calcRegularity(self,edge):
        segLengths=[(edge[i]-edge[i-1]).magnitude() for i in range(1,len(edge))]
        regList=[]
        for i in range(len(segLengths)):
            totalLength=0.0
            endPoint=-1
            for j in range(i,len(segLengths)):
                totalLength+=segLengths[j]
                if totalLength>=self.seglength:
                    endPoint=j+1
                    break
            if endPoint<0:
                break;
            regList.append(((edge[i]-edge[endPoint]).magnitude()/totalLength,totalLength))
        stats = EdgeStatistics()
        if len(regList)==0:
            l=sum(segLengths)
            stats((edge[0]-edge[-1]).magnitude()/l,l)
        else:
            for r in regList:
                stats(r[0],r[1])
        setattr(edge, self.attrName, stats)

    def postMergeEdges(self, survivor):
        self.calcRegularity(survivor)

# --------------------------------------------------------------------

def calculateTangentLists(map, dx = 5, skipPoints = 1):
    """Return a list that for each edge contains the result of running
    `geomap.tangentList(edge, ...)` with the given parameters.
    Special care is taken to ensure that the list will never be empty
    (by reducing the `dx` or finally set `skipPoints` to zero to get at
    least one tangent)."""

    result = [None] * map.maxEdgeLabel()
    badCount = 0
    for edge in map.edgeIter():
        size = len(edge)
        if size >= 2*dx + 2*skipPoints + 1:
            tangents = tangentList(edge, dx, skipPoints)
        else:
            #print "too short:", edge
            badCount += 1
            if size < 3:
                edx, edy = edge[1] - edge[0]
                tangents = [(edge.length()/2, math.atan2(edy, edx))]
            elif size < 3 + 2*skipPoints:
                # we cannot afford skipping points
                tangents = tangentList(edge, 1, 0)
            else:
                # we can skip the end points (if desired), but we have
                # to reduce dx:
                maxDx = (size-2*skipPoints-1)/2
                tangents = tangentList(edge, maxDx, skipPoints)
        assert len(tangents), "all edges should have tangents now"
        result[edge.label()] = tangents

    if badCount:
        print "calculateTangentLists: %d/%d edges were done with " \
              "different parameters (too short)" % (badCount, map.edgeCount)
    return result

def calculateTangentListsGaussianReflective(map, sigma, diff=0.0):
    """calculateTangentListsGaussianReflective(map, sigma, diff=0.0)
    Add 'tangents' property to all edges, containing the result of
    running tangentListGaussianReflective on it with the given
    parameters.  Note that all edges which are too small will have an
    empty 'tangents' list."""
    
    result = [None] * map.maxEdgeLabel()
    for edge in map.edgeIter():
        try:
            result[edge.label()] = tangentListGaussianReflective(edge, sigma, diff)
        except RuntimeError:
            sys.stderr.write("calculateTangentListsGaussianReflective: %s too short!\n" % edge)
            edx, edy = edge[-1] - edge[0]
            result[edge.label()] = [(edge.length()/2, math.atan2(edy, edx))]
    return result

class EdgeTangents(DynamicEdgeStatistics):
    """Stores precomputed tangent lists for each edge and merges then
    when merging edges.  (In the past, EdgeTangents would
    automatically call `calculateTangentLists(*args)` - this has been
    changed in order to make it possible to use any tangent
    calculation method and pass the resulting tangents to the
    constructor instead of *args for the fixed method above.)"""

    __slots__ = ["tangents",
                 "_mergedTangents", "_mergeLabels"]
    
    def __init__(self, map, tangents):
        DynamicEdgeStatistics.__init__(self, map)
        self.tangents = tangents
        self._attachHooks()

    def preMergeEdges(self, dart):
        if dart.label() > 0:
            dart1 = dart.clone().prevPhi()
            dart2 = dart
        else:
            dart1 = dart.clone().nextAlpha()
            dart2 = dart.clone().nextSigma()
        self._mergedTangents = geomap.composeTangentLists([
            self.dartTangents(dart1), self.dartTangents(dart2)])
        self._mergeLabels = (dart1.edgeLabel(), dart2.edgeLabel())
        return True

    def postMergeEdges(self, survivor):
        self.tangents[survivor.label()] = self._mergedTangents
        if self._mergeLabels[0] == survivor.label():
            self.tangents[self._mergeLabels[1]] = None
        else:
            self.tangents[self._mergeLabels[0]] = None

    def __getitem__(self, index):
        return self.tangents[index]

    def dartTangents(self, dart):
        """Return (tangents, dirLength) pair, where dirLength is
        negative if dart.label > 0 (the darts are supposed to be
        reversed for composeTangentLists)"""

        e = dart.edge()
        t = self.tangents[e.label()]
        if dart.label() > 0:
            return (t,  e.length())
        else:
            return (t, -e.length())

# --------------------------------------------------------------------

def gcByArcLength(al):
    """Local measure for good continuation: compares secant directions
    around node.  Range is -1..1, where 1 means perfect continuation."""

    def goodContinuation(dart1, dart2):
        assert dart1.startNode() == dart2.startNode()
        # find one-sided tangents:
        dp1 = geomap.DartPosition(dart1)
        dp1.gotoArcLength(al)
        dp2 = geomap.DartPosition(dart2)
        dp2.gotoArcLength(al)

        p1 = dart2[0] # node position
        t1 = p1 - dp1() # tangents
        t2 = dp2() - p1
        # tangent direction agreement:
        return dot(t1, t2) / (t1.magnitude() * t2.magnitude())

    return goodContinuation
