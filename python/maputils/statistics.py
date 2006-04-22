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
            return (self.sum2 / self.pixelCount) - math.sq(self.sum / self.pixelCount)
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

        map.preMergeFacesHooks.append(self.preMergeFaces)
        map.postMergeFacesHooks.append(self.postMergeFaces)
        map.associatedPixelsHooks.append(self.associatePixels)
        self.map = ref(map) # for detachHooks()

    def detachHooks(self):
        self.map().preMergeFacesHooks.remove(self.preMergeFaces)
        self.map().postMergeFacesHooks.remove(self.postMergeFaces)
        self.map().associatedPixelsHooks.remove(self.associatePixels)

    def preMergeFaces(self, dart):
        self.mergedStats = copy.copy(dart.leftFace()._faceColor)
        self.mergedStats.merge(dart.rightFace()._faceColor)

    def postMergeFaces(self, survivor):
        survivor._faceColor = self.mergedStats

    def associatePixels(self, face, positions):
        for pos in positions:
            face._faceColor(self.originalImage[pos])

class FaceColorStatistics(object):
    def __init__(self, map, cielabImage):
        for face in map.faceIter():
            face._faceMeanCIE_a = FaceStatistics(0.0) # FIXME: one property
            face._faceMeanCIE_b = FaceStatistics(0.0)

        class MeansInitFunctor(object):
            def __init__(self, map):
                self.map = map

            def __call__(self, label, value):
                if label >= 0:
                    self.map.face(int(label))._faceMeanCIE_a(value[1])
                    self.map.face(int(label))._faceMeanCIE_b(value[2])

        inspectImage(map.labelImage, cielabImage, MeansInitFunctor(map))

        map.preMergeFacesHooks.append(self.preMergeFaces)
        map.postMergeFacesHooks.append(self.postMergeFaces)
        self.map = ref(map) # for detachHooks()

    def detachHooks(self):
        self.map().preMergeFacesHooks.remove(self.preMergeFaces)
        self.map().postMergeFacesHooks.remove(self.postMergeFaces)

    def preMergeFaces(self, dart):
        self.mergedStats_a = copy.copy(dart.leftFace()._faceMeanCIE_a)
        self.mergedStats_a.merge(dart.rightFace()._faceMeanCIE_a)
        self.mergedStats_b = copy.copy(dart.leftFace()._faceMeanCIE_b)
        self.mergedStats_b.merge(dart.rightFace()._faceMeanCIE_b)

    def postMergeFaces(self, survivor):
        survivor._faceMeanCIE_a = self.mergedStats_a
        survivor._faceMeanCIE_b = self.mergedStats_b

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

def vecNorm(vec):
    return math.sqrt(dot(vec,vec))

def faceMeanCDiff(dart):
    if not dart.leftFaceLabel() or not dart.rightFaceLabel():
        return 10000
    return vecNorm(dart.leftFace()._faceColor.average() -
                   dart.rightFace()._faceColor.average())/441.7

def faceMeanColDiff(dart):
    return math.hypot(dart.leftFace()._faceMeanCIE_a.average() -
                      dart.rightFace()._faceMeanCIE_a.average(),
                      dart.leftFace()._faceMeanCIE_b.average() -
                      dart.rightFace()._faceMeanCIE_b.average())

def faceStdDevDiff(dart):
    return abs(dart.leftFace()._faceColor.stdDeviation() -
               dart.rightFace()._faceColor.stdDeviation())

def faceStdDevColDiff(dart):
    return math.hypot(dart.leftFace()._faceMeanCIE_a.stdDeviation() -
                      dart.rightFace()._faceMeanCIE_a.stdDeviation(),
                      dart.leftFace()._faceMeanCIE_b.stdDeviation() -
                      dart.rightFace()._faceMeanCIE_b.stdDeviation())

def faceAreaHomogenity(d):
    return (d.leftFace().area() * d.rightFace().area()) / \
           (d.leftFace().area() + d.rightFace().area())

def faceHomogenity(d):
    return math.sq(faceMeanDiff(d)) * faceAreaHomogenity(d)

# --------------------------------------------------------------------
#               Edge-based Statistics & Cost Measures
# --------------------------------------------------------------------

class DynamicEdgeStatistics(object):
    def __init__(self, map):
        map.preMergeEdgesHooks.append(self.preMergeEdges)
        map.postMergeEdgesHooks.append(self.postMergeEdges)
        self.map = ref(map) # for detachHooks()

    def detachHooks(self):
        self.map().preMergeEdgesHooks.remove(self.preMergeEdges)
        self.map().postMergeEdgesHooks.remove(self.postMergeEdges)

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

    def postMergeEdges(self, survivor):
        survivor._passValue = self.mergedPV
        survivor._saddles = self.mergedSaddles

class BoundaryIndicatorStatistics(DynamicEdgeStatistics):
    attrTrans = string.maketrans(".+", "__")

    def __init__(self, map, attrName):
        DynamicEdgeStatistics.__init__(self, map)
        self.attrName = attrName.translate(self.attrTrans)

    def preMergeEdges(self, dart):
        self.mergedStats = copy.copy(getattr(dart.edge(), self.attrName))
        self.mergedStats.merge(
            getattr(dart.clone().nextSigma().edge(), self.attrName))

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
        for edge in map.edgeIter():
            stats = EdgeStatistics()
            s_edge = simplifyPolygon(edge, 0.1)

            for i in range(len(s_edge)-1):
                segment = s_edge[i+1] - s_edge[i]
                p = s_edge[i] + 0.5 * segment
                gradDir = bi.grad.siv[p]
                gradDir /= gradDir.magnitude()
                stats(1.0 - abs(dot(gradDir, segment / segment.magnitude())),
                      segment.magnitude())

            setattr(edge, self.attrName, stats)

class EdgeGradAngDispStatistics(BoundaryIndicatorStatistics):
    def __init__(self, map, bi, splineOrder):
        BoundaryIndicatorStatistics.__init__(self, map, "gad_" + bi.name)
        gad = calcGradAngDisp(bi.grad, 3)
        gad.siv = eval("SplineImageView%d(gad)" % (splineOrder, ))
        for edge in map.edgeIter():
            setattr(edge, self.attrName, EdgeStatistics(edge, gad.siv))

class EdgeMinimumDistance(DynamicEdgeStatistics):
    def __init__(self, map, minima):
        DynamicEdgeStatistics.__init__(self, map)
        minimaMap = PositionedMap()
        for p in minima:
            minimaMap.insert(p, p)
        for edge in map.edgeIter():
            mindist = 1e8
            s_edge = simplifyPolygon(edge, 0.1)
            for p in s_edge:
                d = (minimaMap(p)-p).magnitude()
                if (d < mindist):
                    mindist = d
            edge._minDist = mindist

    def preMergeEdges(self, dart):
        self.mergedMinDist = min(dart.edge()._minDist,
                                 dart.clone().nextSigma().edge()._minDist)

    def postMergeEdges(self, survivor):
        survivor._minDist = self.mergedMinDist

def calcGradScaleSum(image, steps):
    gss = GrayImage(image.size())
    scale = 0.7
    for i in range(0,steps):
        img = gaussianSmoothing(image,scale)
        ti = vectorToTensor(gaussianGradientAsVector(img.subImage(0),1.6))
        if (img.bands()>1):
            for j in range(1,image.bands()):
                ti += vectorToTensor(gaussianGradientAsVector(img.subImage(j),1.6))
        gm = transformImage(tensorTrace(ti),'\l x:sqrt(x)')
        mm = MinMax()
        inspectImage(gm,mm)
        gm=linearRangeMapping(gm,oldRange=(0,mm.max()),newRange=(0,1.0))
        gss += gm
        scale *= 2.0
    return gss

class MeanEdgeGradScaleSum(DynamicEdgeStatistics):
    def __init__(self, map, image, splineOrder):
        gms = calcGradScaleSum(image,5)
        DynamicEdgeStatistics.__init__(self, map)
        avgFunc=Average()
        inspectImage(gms,avgFunc)
        avg=avgFunc.average()
        gms/=avg[0]
        gmsSiv=eval("SplineImageView%d(gms)" % (splineOrder, ))
        for edge in map.edgeIter():
            edge._meanGradScaleSum = EdgeStatistics(edge,gmsSiv)

    def preMergeEdges(self, dart):
        self.mergedStats = copy.copy(dart.edge()._meanGradScaleSum)
        self.mergedStats.merge(dart.clone().nextSigma().edge()._meanGradScaleSum)

    def postMergeEdges(self, survivor):
        survivor._meanGradScaleSum = self.mergedStats

def calcSpatialStabilityImage(image,octaves,radius):
    imglist=[]
    scale=0.7
    for i in range(0,octaves*4+1):
        imglist.append(zeroCrossingImage(image,scale,radius))
        scale*=1.1892
    ssimg=GrayImage(image.size())
    for x in range(0,image.width()):
        for y in range(0,image.height()):
            maximum=0
            current=0
            for i in range(0,len(imglist)):
                if (imglist[i][x,y]>0):
                    current+=1
                    maximum=max(maximum,current)
                else:
                    current=0
            ssimg[x,y]=maximum
    return ssimg

class MeanEdgeSpatialStability(DynamicEdgeStatistics):
    def __init__(self, map, image, splineOrder):
        DynamicEdgeStatistics.__init__(self, map)
        ss = spatialStabilityImage(image, 4, 3)
        ssSiv = eval("SplineImageView%d(ss)" % (splineOrder, ))
        for edge in map.edgeIter():
            edge._meanSpatialStability = EdgeStatistics(edge, ssSiv)

    def preMergeEdges(self, dart):
        self.mergedStats = copy.copy(dart.edge()._meanSpatialStability)
        self.mergedStats.merge(dart.clone().nextSigma().edge()._meanSpatialStability)

    def postMergeEdges(self, survivor):
        survivor._meanSpatialStability = self.mergedStats

def calcGradAngDisp(grad,n):
    ki=GrayImage(2*n+1,2*n+1)
    for x in range(0,ki.width()):
        for y in range(0,ki.height()):
            ki[x,y]=1.0
    kernel=customKernel(ki,(n,n))
    mi=convolveImage(transformImage(grad, '\l x: sqrt(sq(x[0])+sq(x[1]))'),kernel)
    gi=convolveImage(grad,kernel)
    gad=GrayImage(grad.size())
    for x in range(0,gad.width()):
        for y in range(0,gad.height()):
            if (mi[x,y]>0):
                gad[x,y]=math.hypot(gi[x,y][0],gi[x,y][1])/mi[x,y]

    return gad

class MeanEdgeCurvChange(object):
    def __init__(self, map):
        for edge in map.edgeIter():
            self.calcCurvChange(edge)

        map.postMergeEdgesHooks.append(self.postMergeEdges)
        self.map = ref(map) # for detachHooks()

    def calcCurvChange(self,edge):
        edge._meanCurvChange = EdgeStatistics()
        s_edge=simplifyPolygon(edge,0.1)
        for i in range(len(s_edge)-3):
            d1=s_edge[i+2]-s_edge[i+1]
            d2=(s_edge[i]-s_edge[i+1]-s_edge[i+2]+s_edge[i+3])/2.0
            d3=s_edge[i+3]-3.0*s_edge[i+2]+3.0*s_edge[i+1]-s_edge[i]
            if (d1.magnitude()>0):
                cc=((d1[0]*d3[1]-d1[1]*d3[0])*(d1[0]*d1[0]+d1[1]*d1[1])-d2[0]*d2[1]*(d1[0]*d1[0]-d1[1]*d1[1])+d1[0]*d1[1]*(d2[0]*d2[0]+d2[1]*d2[1]))/(d1[0]*d1[0]+d1[1]*d1[1])
                edge._meanCurvChange(cc, d1.magnitude())

    def detachHooks(self):
        self.map().postMergeEdgesHooks.remove(self.postMergeEdges)

    def postMergeEdges(self, survivor):
        self.calcCurvChange(survivor)

class MeanEdgeKRegularity(object):
    def __init__(self, map, image):
        for edge in map.edgeIter():
            self.calcKRegularity(edge,8)

        map.postMergeEdgesHooks.append(self.postMergeEdges)
        self.map = ref(map) # for detachHooks()

    def calcKRegularity(self,edge,k):
        edge._meanKRegularity = EdgeStatistics()
        s_edge=simplifyPolygon(edge,0.1)
        if (len(s_edge)<k+1):
            k=len(s_edge)-1

        for i in range(len(s_edge)-k):
            segment = s_edge[i+k] - s_edge[i]
            length=0
            for j in range(0,k):
                length+=(s_edge[i+j]-s_edge[i+j+1]).magnitude()
            edge._meanKRegularity(segment.magnitude()/length, segment.magnitude())

    def detachHooks(self):
        self.map().postMergeEdgesHooks.remove(self.postMergeEdges)

    def postMergeEdges(self, survivor):
        self.calcKRegularity(survivor,8)

def otherCommonDarts(dart):
    it = dart.clone()
    rightFaceLabel = dart.rightFaceLabel()
    while it.nextPhi() != dart:
        if it.rightFaceLabel() == rightFaceLabel:
            yield it

def totalBoundaryStatistics(dart):
    result = EdgeStatistics()
    result.merge(dart.edge()._meanGradient)
    for d in otherCommonDarts(dart):
        result.merge(d.edge()._meanGradient)
    return result

def meanEdgeGradCost(dart):
    return totalBoundaryStatistics(dart).average()

def medianEdgeGradCost(dart):
    return totalBoundaryStatistics(dart).median()

def polyMedianEdgeGradCost(dart):
    return totalBoundaryStatistics(dart).polyMedian()

def minEdgeGradCost(dart):
    result = dart.edge()._passValue
    for d in otherCommonDarts(dart):
        result = min(result, d.edge()._passValue)
    return result
