_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import math, string, copy, weakref
from vigra import *
from hourglass import PositionedMap, EdgeStatistics, \
     FaceGrayStatistics, FaceRGBStatistics, \
     spatialStabilityImage, tangentList, resamplePolygon
import sivtools, flag_constants
from map import arcLengthIter

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

def _combinedMeasure(dart, weightedMeasures):
    cost = 0.0
    for weight, measure in weightedMeasures:
        cost += weight * measure(dart)
    return cost

def combinedMeasure(*weightedMeasures):
    return lambda dart: _combinedMeasure(dart, weightedMeasures)

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
    def __init__(self, map, originalImage,
                 defaultValue = None, minSampleCount = 1,
                 SIV = SplineImageView5):
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

    # FIXME: wrap in sqrt?
    def faceHomogenity(self, d):
        return math.sq(self.faceMeanDiff(d)) * faceAreaHomogenity(d)

    # FIXME: double-check, comment on deficiencies
    def faceTTest(self, dart):
        f1 = self.faceMeanFunctor(dart.leftFaceLabel())
        f2 = self.faceMeanFunctor(dart.rightFaceLabel())
        return norm(f1.average() - f2.average()) / \
           max(math.sqrt(f1.variance() / f1.pixelCount + \
                         f2.variance() / f2.pixelCount), 1e-3)

FaceGrayStatistics.bands = lambda x: 1
FaceRGBStatistics.bands = lambda x: 3

def FaceColorStatistics(map, originalImage, minSampleCount = 1):
    if originalImage.bands() == 1:
        return FaceGrayStatistics(map, originalImage, minSampleCount)
    elif originalImage.bands() == 3:
        return FaceRGBStatistics(map, originalImage, minSampleCount)
    else:
        return _FaceColorStatistics(map, originalImage, minSampleCount,
                                    SIV = sivtools.GradientSIVProxy)

def faceAreaHomogenity(dart):
    a1 = dart.leftFace().area()
    a2 = dart.rightFace().area()
    return (a1 * a2) / (a1 + a2)

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

class EdgeMergeTree(DynamicEdgeStatistics):
    """Actually, this is not a tree but it manages a list of edges
    that have been merged into each edge."""

    __slots__ = ["_tree",
                 "_merged"]
    
    def __init__(self, map):
        DynamicEdgeStatistics.__init__(self, map)
        self._tree = range(map.maxEdgeLabel())
        self._attachHooks()
    
    def preMergeEdges(self, dart):
        self._merged = dart.clone().nextSigma().edgeLabel()
        return True
    
    def postMergeEdges(self, survivor):
        label = survivor.label()
        while True: # search list end
            merged = self._tree[label]
            if merged == label:
                break
            label = merged
        self._tree[label] = self._merged # concatenate lists
    
    def __getitem__(self, edge):
        """Returns list of all edges of original map of which the
        given edge of the current map is composed."""
        if hasattr(edge, "label"):
            edge = edge.label()
        result = [edge]
        while True:
            merged = self._tree[edge]
            if merged == edge:
                break
            result.append(merged)
            edge = merged
        return result

    def __getstate__(self):
        return DynamicEdgeStatistics.__getstate__(self) + (
            self._tree, )

    def __setstate__(self, (map, tree)):
        DynamicEdgeStatistics.__setstate__(self, (map, ))
        self._tree = tree

    def level0Labels(self):
        """Return a list of level0 edges (resp. their labels) the
        currently existing edges are composed of."""
        result = []
        for edge in self._map().edgeIter():
            result.extend(self[edge])
        return result

class DynamicEdgeIndices(DetachableStatistics):
    __slots__ = ["_indices",
                 "_mergedIndices", "_newIndices1", "_newIndices2"]
    
    def __init__(self, map):
        DetachableStatistics.__init__(self, map)
        self._indices = [[] for i in range(map.maxEdgeLabel())]
    
    def _attachHooks(self):
        self._attachedHooks = (
            self._map().addMergeEdgesCallbacks(self.preMergeEdges,
                                               self.postMergeEdges),
            self._map().addSplitEdgeCallbacks(self.preSplitEdge,
                                              self.postSplitEdge))

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
    __slots__ = ["_passValues", "_indices", "_gmSiv",
                 "_basinDepth",
                 "_mergedPV"]
    
    def __init__(self, map, flowlines, gmSiv):
        DynamicEdgeIndices.__init__(self, map)
        self._passValues = [None] * map.maxEdgeLabel()
        for edge in map.edgeIter():
            if edge.flag(flag_constants.BORDER_PROTECTION):
                continue
            
            saddleIndex = flowlines[edge.label()][3]

            # flowline tracing might have stopped, in which case the
            # polygon can be modified by connecting to the nearest
            # maximum/node, possibly adding a point and thus shifting
            # the indices:
            flowline = flowlines[edge.label()]
            if len(edge) != len(flowline[2]) and edge[1] == flowline[2][0]:
                assert flowline[0] <= 0
                saddleIndex += 1
            
            self._passValues[edge.label()] = gmSiv[edge[saddleIndex]]
            self._indices[edge.label()].append(saddleIndex)

        self._gmSiv = gmSiv
        self._attachHooks()

    def __getstate__(self):
        return DynamicEdgeIndices.__getstate__(self) + (
            self._passValues, self._indices)

    def __setstate__(self, (map, passValues, indices)):
        DynamicEdgeIndices.__setstate__(self, (map, ))
        self._passValues = passValues
        self._indices = indices
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
        DynamicEdgeIndices.preMergeEdges(self, dart)

        edge1 = dart.edge()
        edge2 = dart.clone().nextSigma().edge()

        self._mergedPV = None

        # do not allow merging of watersheds with border:
        if edge1.flag(flag_constants.BORDER_PROTECTION):
            return edge2.flag(flag_constants.BORDER_PROTECTION)
        if edge2.flag(flag_constants.BORDER_PROTECTION):
            return edge1.flag(flag_constants.BORDER_PROTECTION)

        if not DynamicEdgeIndices.preMergeEdges(self, dart):
            return False
        
        self._mergedPV = min(self._passValues[edge1.label()],
                             self._passValues[edge2.label()])

        return True

    def postMergeEdges(self, survivor):
        DynamicEdgeIndices.postMergeEdges(self, survivor)

        self._passValues[survivor.label()] = self._mergedPV

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

    def dartPassValue(self, dart):
        return self._passValues[dart.edgeLabel()]

    def passValue(self, dart):
        return min([self.dartPassValue(d)
                    for d in commonBoundaryDarts(dart)])

    def setBasins(self, basinStatistics):
        self._basinDepth = basinStatistics._basinDepth

    def dynamic(self, edge):
        if hasattr(edge, "label"):
            edgeLabel = edge.label()
        else:
            edgeLabel = edge
            edge = self._map().edge(edgeLabel)
        return self._passValues[edgeLabel] - min(
            self._basinDepth[edge.leftFaceLabel()],
            self._basinDepth[edge.rightFaceLabel()])

class WatershedBasinStatistics(DetachableStatistics):
    __slots__ = ["_basinDepth",
                 "_mergedDepth"]
    
    def __init__(self, map, minima, gmSiv):
        DetachableStatistics.__init__(self, map)

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

        for face in map.faceIter(skipInfinite = True):
            if self._basinDepth[face.label()] == None:
                sys.stderr.write(
                    "Face %d (area %s, anchor %d) contains no minimum!\n" % (
                    face.label(), face.area(), face.contour().label()))
                level = 2
                while level < 20: # prevent endless loop
                    level += 1
                    samples = [gmSiv[pos] for pos in superSample(face, level)]
                    if len(samples) > 10:
                        break
                print "  (got %d samples at supersampling level %d)" % (
                    len(samples), level)
                self._basinDepth[face.label()] = min(samples)
                assert self._basinDepth[face.label()] != None, str(samples)

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
        return DetachableStatistics.__getstate__(self) + (
            self._basinDepth, )

    def __setstate__(self, (map, basinDepth)):
        DetachableStatistics.__setstate__(self, (map, ))
        self._basinDepth = basinDepth

def _makeAttrName(someStr):
    attrTrans = string.maketrans(".+-", "__n")
    return someStr.translate(attrTrans)

class BoundaryIndicatorStatistics(DynamicEdgeStatistics):
    __slots__ = ["_functors",
                 "_mergedStats"]
    
    def __init__(self, map):
        DynamicEdgeStatistics.__init__(self, map)
        self._functors = [None] * map.maxEdgeLabel()

    def preMergeEdges(self, dart):
        self._mergedStats = copy.copy(self._functors[dart.edgeLabel()])
        self._mergedStats.merge(
            self._functors[dart.clone().nextSigma().edgeLabel()])
        return True

    def __getitem__(self, index):
#         if hasattr(index, "label"):
#             index = index.label()
        return self._functors[index]

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
    support points.  Each edge has an associated `EdgeStatistics`
    object which manages the average value and a list of sorted values
    for quantile queries."""

    __slots__ = []
    
    def __init__(self, map, gmSiv, resample = 0.1):
        BoundaryIndicatorStatistics.__init__(self, map)
        for edge in map.edgeIter():
            poly = resample and resamplePolygon(edge, resample) or edge
            self._functors[edge.label()] = EdgeStatistics(poly, gmSiv)
        self._attachHooks()

    def __getitem__(self, edgeLabel):
        return self._functors[edgeLabel]

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
        return self[dart.edgeLabel()].average()

    def combinedStatistics(self, dart):
        result = EdgeStatistics()
        for d in commonBoundaryDarts(dart):
            result.merge(self[d.edgeLabel()])
        return result

    def min(self, dart):
        return self.combinedStatistics(dart).min()

    def max(self, dart):
        return self.combinedStatistics(dart).max()

    def quantile(self, q):
        def specificQuantile(dart):
            return self.combinedStatistics(dart).quantile(q)
        return specificQuantile

    def average(self, dart):
        return self.combinedStatistics(dart).average()

# USAGE:
# >>> boundaryIndicator = SplineImageView5(...)
# >>> egs = EdgeGradientStatistics(map, boundaryIndicator)
# >>> egs[100].average()
# 1.178475851870056
# >>> egs[100].quantile(0.4)
# 1.1785515546798706

class EdgeGradDirDotStatistics(BoundaryIndicatorStatistics):
    def __init__(self, map, bi):
        BoundaryIndicatorStatistics.__init__(self, map)

        assert hasattr(map.edgeIter().next(), "tangents"), \
               """Edge does not have 'tangents' attribute!
               Use calculateTangentLists(myMap[, ...]) to initialize tangent lists!"""
        
        for edge in map.edgeIter():
            stats = EdgeStatistics()

            ali = arcLengthIter(edge)
            prevPoint = ali.next()
            curPoint = ali.next()
            for al, theta in edge.tangents:
                # seek to the right edge segment:
                while al > curPoint[0]:
                    prevPoint = curPoint
                    curPoint = ali.next()
                    
                # linearly interpolate to the desired arclength (al):
                pos = prevPoint[1] + \
                      (al - prevPoint[0])/(curPoint[0] - prevPoint[0]) * \
                      (curPoint[1] - prevPoint[1])

                gradDir = bi.grad.siv[pos]
                gradDir /= gradDir.magnitude()

                # FIXME: QuantileStatistics expects segment length,
                # we give always 1.0:
                segment = Vector2(math.cos(theta), math.sin(theta))
                stats(1.0 - abs(dot(gradDir, segment)), 1.0)

            self._functors[edge.label()] = stats

        self._attachHooks()

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
        minimaMap = PositionedMap()
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
        ss = spatialStabilityImage(image, 4, radius, propexp)
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

from hourglass import fitLine

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

from hourglass import ParabolaFit

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
    """calculateTangentLists(map, dx = 5, skipPoints = 1)
    Add 'tangents' property to all edges, containing the result of
    running tangentList() on it with the given parameters.  Special
    care is taken to ensure that the list will never be empty (by
    reducing the dx or finally set skipPoints to zero to get at least
    one tangent)."""
    
    assert not skipPoints > 1, "why that??"
    result = 0
    for edge in map.edgeIter():
        size = len(edge)
        if size >= 2*dx + 2*skipPoints + 1:
            edge.tangents = tangentList(edge, dx, skipPoints)
        else:
            #print "too short:", edge
            result += 1
            if size < 3:
                edx, edy = edge[1] - edge[0]
                edge.tangents = [(edge.length()/2, math.atan2(edy, edx))]
            elif size < 3 + 2*skipPoints:
                # we cannot afford skipping points
                edge.tangents = tangentList(edge, 1, 0)
            else:
                # we can skip the end points (if desired), but we have
                # to reduce dx:
                maxDx = (size-2*skipPoints-1)/2
                edge.tangents = tangentList(edge, maxDx, skipPoints)
        assert len(edge.tangents), "all edges should have tangents now"
    print "calculateTangentLists: %d/%d edges were done with " \
          "different parameters (too short)" % (result, map.edgeCount)
    return result

def calculateTangentListsGaussianReflective(map, sigma, diff=0.0):
    """calculateTangentListsGaussianReflective(map, sigma, diff=0.0)
    Add 'tangents' property to all edges, containing the result of
    running tangentListGaussianReflective on it with the given
    parameters.  Note that all edges which are too small will have an
    empty 'tangents' list."""
    
    for edge in map.edgeIter():
        try:
            edge.tangents = tangentListGaussianReflective(edge, sigma, diff)
        except RuntimeError:
            edge.tangents = []
            print "too short:", edge

def dartTangents(dart):
    """dartTangents(dart)

    Returns (tangents, dirLength) pair, where dirLength is
    negative if dart.label > 0 (the darts are supposed to be
    reversed for composeTangentLists)"""

    e = dart.edge()
    if dart.label() > 0:
        return (e.tangents,  e.length())
    else:
        return (e.tangents, -e.length())

class EdgeTangents(DynamicEdgeStatistics):
    """Calls calculateTangentLists with the constructor parameters (if
    given) and merges tangent lists when merging edges."""
    
    def __init__(self, map, *args):
        """Either make sure to call calculateTangentLists() or one of
        its friends on the map by hand, or pass parameters like this:
        et = EdgeTangents(someMap, 4)"""
        if args:
            calculateTangentLists(map, *args)
        self._attachHooks()

    def preMergeEdges(self, dart):
        import dartpath # prevent cyclic import
        path = dartpath.Path(
            [dart.clone().nextAlpha(), dart.clone().nextSigma()])
        if path[0].label() < 0:
            path.reverse()
        self.newTangents = path.tangents()

    def postMergeEdges(self, survivor):
        survivor.tangents = self.newTangents
