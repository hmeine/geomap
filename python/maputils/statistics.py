_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import math, string, copy
from weakref import ref
from vigra import *
from hourglass import PositionedMap, EdgeStatistics, \
     spatialStabilityImage, tangentList
from map import arcLengthIter

# --------------------------------------------------------------------
#              Region-based Statistics & Cost Measures
# --------------------------------------------------------------------

class FaceStatistics(object):
    def __init__(self, defaultValue):
        self.pixelCount = None
        self.defaultValue = defaultValue

    def average(self):
        if self.pixelCount:
            return self.sum / self.pixelCount
        return self.defaultValue

    def variance(self):
        if self.pixelCount:
            return max(0.0,  (self.sum2 / self.pixelCount) -
                       math.sq(self.sum / self.pixelCount))
        return self.defaultValue

    def stdDeviation(self):
        return math.sqrt(self.variance())

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

class FaceIntensityStatistics(object):
    def __init__(self, map, originalImage, defaultValue = None):
        self.originalImage = originalImage

        if defaultValue == None:
            defaultValue = originalImage[0,0]

        for face in map.faceIter():
            face._faceColor = FaceStatistics(defaultValue)

        class MeansInitFunctor(object):
            def __init__(self, map):
                self.map = map

            def __call__(self, label, value):
                if label >= 0:
                    self.map.face(int(label))._faceColor(value)

        inspectImage(map.labelImage, originalImage, MeansInitFunctor(map))

        self._attachedHooks = (
            map.addMergeFacesCallbacks(self.preMergeFaces, self.postMergeFaces) +
            map.addAssociatePixelsCallback(self.associatePixels))

    def detachHooks(self):
        self._attachedHooks.disconnect()

    def preMergeFaces(self, dart):
        self.mergedStats = copy.copy(dart.leftFace()._faceColor)
        self.mergedStats.merge(dart.rightFace()._faceColor)
        return True

    def postMergeFaces(self, survivor):
        survivor._faceColor = self.mergedStats

    def associatePixels(self, face, positions):
        for pos in positions:
            face._faceColor(self.originalImage[pos])

class FaceColorStatistics(object):
    def __init__(self, map, cielabImage):
        self.cielabImage=cielabImage
        for face in map.faceIter():
            face._faceMeanCIE_L = FaceStatistics(0.0) # FIXME: one property
            face._faceMeanCIE_a = FaceStatistics(0.0)
            face._faceMeanCIE_b = FaceStatistics(0.0)

        class MeansInitFunctor(object):
            def __init__(self, map):
                self.map = map

            def __call__(self, label, value):
                if label >= 0:
                    self.map.face(int(label))._faceMeanCIE_L(value[0])
                    self.map.face(int(label))._faceMeanCIE_a(value[1])
                    self.map.face(int(label))._faceMeanCIE_b(value[2])

        inspectImage(map.labelImage, cielabImage, MeansInitFunctor(map))

        self._attachedHooks = (
            map.addMergeFacesCallbacks(self.preMergeFaces, self.postMergeFaces) +
            map.addAssociatePixelsCallback(self.associatePixels))

    def detachHooks(self):
        self._attachedHooks.disconnect()

    def preMergeFaces(self, dart):
        self.mergedStats_L = copy.copy(dart.leftFace()._faceMeanCIE_L)
        self.mergedStats_L.merge(dart.rightFace()._faceMeanCIE_L)
        self.mergedStats_a = copy.copy(dart.leftFace()._faceMeanCIE_a)
        self.mergedStats_a.merge(dart.rightFace()._faceMeanCIE_a)
        self.mergedStats_b = copy.copy(dart.leftFace()._faceMeanCIE_b)
        self.mergedStats_b.merge(dart.rightFace()._faceMeanCIE_b)
        return True

    def postMergeFaces(self, survivor):
        survivor._faceMeanCIE_L = self.mergedStats_L
        survivor._faceMeanCIE_a = self.mergedStats_a
        survivor._faceMeanCIE_b = self.mergedStats_b

    def associatePixels(self, face, positions):
        for pos in positions:
            face._faceMeanCIE_L(self.cielabImage[pos][0])
            face._faceMeanCIE_a(self.cielabImage[pos][1])
            face._faceMeanCIE_b(self.cielabImage[pos][2])

class FaceColorHistogram(object):
    def __init__(self, map, image):
        self.image=image
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

        inspectImage(map.labelImage, image, HistInitFunctor(map))

        for face in map.faceIter():
            face._colorHistogram.gaussianSmoothing(2.2)
            face._colorHistogram2.gaussianSmoothing(2.2)

        self._attachedHooks = (
            map.addMergeFacesCallbacks(self.preMergeFaces, self.postMergeFaces) +
            map.addAssociatePixelsCallback(self.associatePixels))

    def detachHooks(self):
        self._attachedHooks.disconnect()

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

def faceImage(map):
    class MeanLookupFunctor(object):
        def __init__(self, map):
            self.map = map
            fi = map.faceIter()
            fi.next() # skip region 0
            self.hasColor = type(fi.next()._faceColor.average()) != float
            if not self.hasColor:
                self.default = Pixel(0.0)
            else:
                self.default = Pixel(0, 0, 0)

        def __call__(self, label):
            if label >= 0:
                return self.map.face(int(label))._faceColor.average()
            return self.default

    mlf = MeanLookupFunctor(map)
    return transformImage(map.labelImage, mlf)

def faceMeanDiff(dart):
    return abs(dart.leftFace()._faceColor.average() -
               dart.rightFace()._faceColor.average())/255

def faceMeanDiffPoisson(dart):
    m1=dart.leftFace()._faceColor.average()
    m2=dart.rightFace()._faceColor.average()
    n1=float(dart.leftFace()._faceColor.pixelCount)
    n2=float(dart.rightFace()._faceColor.pixelCount)
    n=n1+n2
    m=(n1*m1+n2*m2)/n
    return n1*m1*math.log(m1)+n2*m2*math.log(m2)-n*m*math.log(m)

def vecNorm(vec):
    return math.sqrt(dot(vec,vec))

def faceMeanCDiff(dart):
    if not dart.leftFaceLabel() or not dart.rightFaceLabel():
        return 10000
    return vecNorm(dart.leftFace()._faceColor.average() -
                   dart.rightFace()._faceColor.average())/441.7

def faceMeanDiff_ab(dart):
    return Vector(dart.leftFace()._faceMeanCIE_a.average() -
                  dart.rightFace()._faceMeanCIE_a.average(),
                  dart.leftFace()._faceMeanCIE_b.average() -
                  dart.rightFace()._faceMeanCIE_b.average()).magnitude()

def faceMeanDiff_Lab(dart):
    return Vector(dart.leftFace()._faceMeanCIE_L.average() -
                  dart.rightFace()._faceMeanCIE_L.average(),
                  dart.leftFace()._faceMeanCIE_a.average() -
                  dart.rightFace()._faceMeanCIE_a.average(),
                  dart.leftFace()._faceMeanCIE_b.average() -
                  dart.rightFace()._faceMeanCIE_b.average()).magnitude()

def faceStdDevDiff(dart):
    return abs(dart.leftFace()._faceColor.stdDeviation() -
               dart.rightFace()._faceColor.stdDeviation())

def faceStdDevDiff_ab(dart):
    return Vector(dart.leftFace()._faceMeanCIE_a.stdDeviation() -
                  dart.rightFace()._faceMeanCIE_a.stdDeviation(),
                  dart.leftFace()._faceMeanCIE_b.stdDeviation() -
                  dart.rightFace()._faceMeanCIE_b.stdDeviation()).magnitude()

def faceStdDevDiff_Lab(dart):
    return Vector(dart.leftFace()._faceMeanCIE_L.stdDeviation() -
                  dart.rightFace()._faceMeanCIE_L.stdDeviation(),
                  dart.leftFace()._faceMeanCIE_a.stdDeviation() -
                  dart.rightFace()._faceMeanCIE_a.stdDeviation(),
                  dart.leftFace()._faceMeanCIE_b.stdDeviation() -
                  dart.rightFace()._faceMeanCIE_b.stdDeviation()).magnitude()

def faceAreaHomogenity(d):
    return (d.leftFace().area() * d.rightFace().area()) / \
           (d.leftFace().area() + d.rightFace().area())

def faceHomogenity(d):
    return math.sq(faceMeanDiff(d)) * faceAreaHomogenity(d)

def faceTTest(dart):
    if (dart.leftFace()._faceColor.pixelCount==0) or (dart.rightFace()._faceColor.pixelCount==0):
        return 0
    return abs(dart.leftFace()._faceColor.average() - dart.rightFace()._faceColor.average()) / \
           max(math.sqrt(dart.leftFace()._faceColor.variance()/dart.leftFace()._faceColor.pixelCount + \
           dart.rightFace()._faceColor.variance()/dart.rightFace()._faceColor.pixelCount),0.001)

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

class DynamicEdgeStatistics(object):
    def __init__(self, map):
        self._attachedHooks = map.addMergeEdgesCallbacks(
            self.preMergeEdges, self.postMergeEdges)

    def detachHooks(self):
        self._attachedHooks.disconnect()

class WatershedStatistics(DynamicEdgeStatistics):
    def __init__(self, map, flowlines, gmSiv):
        DynamicEdgeStatistics.__init__(self, map)
        for edge in map.edgeIter():
            edge._passValue = None
            edge._saddles = []

            saddleIndex = None
            if hasattr(edge, "isSplitResultOf"): # split result?
                edgeLabel = abs(edge.isSplitResultOf[0])
                if edgeLabel < len(flowlines):
                    saddleIndex = flowlines[edgeLabel][3] \
                                  - edge.isSplitResultOf[1] # saddle index offset
            elif edge._label < len(flowlines): # could be a border edge, too
                saddleIndex = flowlines[edge._label][3]
            
            if saddleIndex >= 0 and saddleIndex < len(edge):
                edge._passValue = gmSiv[edge[saddleIndex]]
                edge._saddles.append(saddleIndex)
            else:
                edge._passValue = min(gmSiv[edge[0]], gmSiv[edge[-1]])

    def preMergeEdges(self, dart):
        edge1 = dart.edge()
        edge2 = dart.clone().nextSigma().edge()
        self.mergedSaddles = list(edge1._saddles)
        self.mergedSaddles.extend(edge2._saddles)
        self.mergedPV = min(edge1._passValue, edge2._passValue)
        return True

    def postMergeEdges(self, survivor):
        survivor._passValue = self.mergedPV
        survivor._saddles = self.mergedSaddles

def _makeAttrName(someStr):
    attrTrans = string.maketrans(".+-", "__n")
    return someStr.translate(attrTrans)

class BoundaryIndicatorStatistics(DynamicEdgeStatistics):
    def __init__(self, map, attrName):
        DynamicEdgeStatistics.__init__(self, map)
        self.attrName = _makeAttrName(attrName)

    def preMergeEdges(self, dart):
        self.mergedStats = copy.copy(getattr(dart.edge(), self.attrName))
        self.mergedStats.merge(
            getattr(dart.clone().nextSigma().edge(), self.attrName))
        return True

    def postMergeEdges(self, survivor):
        setattr(survivor, self.attrName, self.mergedStats)

class EdgeGradientStatistics(BoundaryIndicatorStatistics):
    def __init__(self, map, bi):
        BoundaryIndicatorStatistics.__init__(self, map, "mag_" + bi.name)
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, bi.gm.siv))

# USAGE:
# >>> boundaryIndicator = ...
# >>> boundaryIndicator.name = 'grad_1.4'
# >>> EdgeGradientStatistics(map, boundaryIndicator)
# >>> map.edge(100).mag_grad_1_4.average()
# 1.178475851870056
# >>> map.edge(100).mag_grad_1_4.quantile(0.4)
# 1.1785515546798706

class EdgeGradDirDotStatistics(BoundaryIndicatorStatistics):
    def __init__(self, map, bi):
        BoundaryIndicatorStatistics.__init__(self, map, "dirdot_" + bi.name)

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

            setattr(edge, self.attrName, stats)

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

class EdgeMinimumDistance(DynamicEdgeStatistics):
    def __init__(self, map, minima):
        DynamicEdgeStatistics.__init__(self, map)
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

class EdgeSpatialStability(BoundaryIndicatorStatistics):
    def __init__(self, map, image, radius, propexp, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "spatialStability_%s_%s" % (radius, propexp))
        ss = spatialStabilityImage(image, 4, radius, propexp)
        ss.siv = eval("SplineImageView%d(ss)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, ss.siv))

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

class EdgeColGradProd(BoundaryIndicatorStatistics):
    def __init__(self, map, image, sigma1, sigma2, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "colGradProd_%s_%s" % (sigma1, sigma2))
        cgp = calcColGradProd(image, sigma1, sigma2)
        cgp.siv = eval("SplineImageView%d(cgp)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, cgp.siv))

def calcBTPhaseImage(image,scale):
    filterImage=getBTFilterResponses(image,scale)
    dirImage=transformImage(tensorEigenRepresentation(transformImage(filterImage,"\l x: Vector(sq(x[0])+sq(x[1])+sq(x[3]+x[5]),-x[1]*(x[0]+x[2])-(x[3]+x[5])*(x[4]+x[6]),sq(x[1])+sq(x[2])+sq(x[4]+x[6]))")),"\l x: Vector(cos(-x[2]), sin(-x[2]))")
    quadImage=transformImage(filterImage, dirImage,"\l x, y: Vector(sq(y[0])*x[0]+2.0*y[0]*y[1]*x[1]+sq(y[1])*x[2],y[0]*(x[3]+x[5])+y[1]*(x[4]+x[6]))")
#    phaseImage=transformImage(quadImage,"\l x: abs(abs(atan2(x[1],x[0]))-1.5708)")
    phaseImage=transformImage(quadImage,"\l x: atan2(x[1],x[0])")
    return phaseImage

class EdgePhase(DynamicEdgeStatistics):
    def __init__(self, map, image, scaleList, splineOrder):
        DynamicEdgeStatistics.__init__(self, map)
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

class EdgeContAngle(object):
    def __init__(self, map, length):
        self.length=length
        self.attrName = _makeAttrName("contAngle_%s" % (self.length))
        for edge in map.edgeIter():
            self.calcContAngle(edge)
        
        self._attachedHooks = map.addMergeFacesCallbacks(
            None, self.postMergeEdges)

    def detachHooks(self):
        self._attachedHooks.disconnect()

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

class EdgeCurvChange(object):
    def __init__(self, map):
        for edge in map.edgeIter():
            self.calcCurvChange(edge)
        self._attachedHooks = map.addMergeFacesCallbacks(
            None, self.postMergeEdges)

    def detachHooks(self):
        self._attachedHooks.disconnect()

    def calcCurvChange(self, edge):
        # use residuum of line fit on tangents as curvature change definition:
        # (FIXME: well-defined residuum? cf. Svens comments... ;*) )
        edge.curvChange = fitLine(edge.tangents)[2]

    def postMergeEdges(self, survivor):
        self.calcCurvChange(survivor)

from hourglass import ParabolaFit

class EdgeCurvChangeLin(object):
    def __init__(self, map):
        for edge in map.edgeIter():
            self.calcCurvChangeLin(edge)
        self._attachedHooks = map.addMergeFacesCallbacks(
            None, self.postMergeEdges)

    def detachHooks(self):
        self._attachedHooks.disconnect()

    def calcCurvChangeLin(self, edge):
        # use residuum of parabola fit on tangents:
        # (FIXME: well-defined residuum? cf. Svens comments... ;*) )
        fit = ParabolaFit()
        fit.addTangentList(edge.tangents)
        edge.curvChangeLin = fit.sumOfSquaredErrors() / len(edge.tangents)

    def postMergeEdges(self, survivor):
        self.calcCurvChangeLin(survivor)

class EdgeRegularity(object):
    def __init__(self, map, seglength):
        self.seglength=seglength
        self.attrName = _makeAttrName("regularity_%s" % (self.seglength))
        for edge in map.edgeIter():
            self.calcRegularity(edge)

        self._attachedHooks = map.addMergeFacesCallbacks(
            None, self.postMergeEdges)

    def detachHooks(self):
        self._attachedHooks.disconnect()

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

def otherCommonDarts(dart):
    it = dart.clone()
    rightFaceLabel = dart.rightFaceLabel()
    while it.nextPhi() != dart:
        if it.rightFaceLabel() == rightFaceLabel:
            yield it

def totalBoundaryStatistics(dart, attrName):
    result = EdgeStatistics()
    result.merge(getattr(dart.edge(), attrName))
    for d in otherCommonDarts(dart):
        result.merge(getattr(d.edge(), attrName))
    return result

def meanEdgeGradCost(dart):
    return totalBoundaryStatistics(dart, "_meanGradient").average()

def medianEdgeGradCost(dart):
    return totalBoundaryStatistics(dart, "_meanGradient").quantile(0.5)

def minEdgeGradCost(dart):
    result = dart.edge()._passValue
    for d in otherCommonDarts(dart):
        result = min(result, d.edge()._passValue)
    return result

def edgeQuantileCost(attrName, quantile):
    def boundaryStatisticsQuantile(dart, attrName = attrName, quantile = quantile):
        return totalBoundaryStatistics(dart, attrName).quantile(quantile)
    return boundaryStatisticsQuantile

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
        DynamicEdgeStatistics.__init__(self, map)
        if args:
            calculateTangentLists(map, *args)

    def preMergeEdges(self, dart):
        import dartpath # prevent cyclic import
        path = dartpath.Path(
            [dart.clone().nextAlpha(), dart.clone().nextSigma()])
        if path[0].label() < 0:
            path.reverse()
        self.newTangents = path.tangents()

    def postMergeEdges(self, survivor):
        survivor.tangents = self.newTangents
