from vigra import Vector2, Size2D, labelImage4
from hourglass import Polygon
import maputils

from hourglass import GeoMap, contourPoly
#from map import GeoMap
execfile("testSPWS")

# --------------------------------------------------------------------
# 				flowline map creation / face.contains
# --------------------------------------------------------------------

# The following data contains edges that run out of the image range,
# which ATM leads to overlapping edges after border closing. That
# violates some assumptions and would lead to errors if we actually
# worked with that Map.
map = GeoMap(maxima2, [], Size2D(39, 39))
maputils.addFlowLinesToMap(flowlines2, map)
assert map.checkConsistency(), "graph inconsistent"
map.sortEdgesEventually(stepDist = 0.2, minDist = 0.05)
map.splitParallelEdges()
map.initializeMap()
assert map.checkConsistency(), "map inconsistent"
assert maputils.checkLabelConsistency(map), "map.labelImage() inconsistent"

# merge faces so that survivor has a hole:
hole = map.faceAt((16,17))
assert hole.contains((16,17))
assert contourPoly(hole.contour()).contains((16,17))

dart = hole.contour()
while True:
    while dart.startNode().degree() > 2:
        face = maputils.removeEdge(dart.clone().prevSigma())
    if dart.nextPhi() == hole.contour():
        break

assert face.holeCount() > 0 # should have hole
for p in [(3, 24), (10, 15), (13, 21)]: # in region, but not within hole
    assert face.contains(p)

for p in [(16, 17), (13, 17)]: # in hole
    assert not face.contains(p)
    assert hole.contains(p)

# --------------------------------------------------------------------
# 					  copying the map from above
# --------------------------------------------------------------------

# test copying / init from map data with disconnected contours:
import copy
om = copy.copy(map)
assert om.__getstate__() == map.__getstate__()
assert om.checkConsistency(), "map inconsistent"
assert maputils.checkLabelConsistency(om), "om.labelImage() inconsistent"

# --------------------------------------------------------------------
# 				   simple map, WatershedStatistics
# --------------------------------------------------------------------

from vigra import readImage, resizeImageSplineInterpolation, transformImage, gaussianGradientAsVector, gaussianGradientMagnitude, SplineImageView5
from hourglass import SubPixelWatersheds5

img = readImage("test_blox.png")
img = resizeImageSplineInterpolation(img, Size2D(128, 128))
img.gm = gaussianGradientMagnitude(img, 1.0)
# img.grad = gaussianGradientAsVector(img, 1.0)
# img.gm = transformImage(img.grad, "\l x: norm(x)", {})
img.gm.siv = SplineImageView5(img.gm)
img.spws = SubPixelWatersheds5(img.gm)

import maputils
maxima, flowlines = maputils.subpixelWatershedData(
    img.spws, img.gm.siv, 0.7, perpendicularDistEpsilon = None)

map = GeoMap(maxima, [], img.size())
maputils.addFlowLinesToMap(flowlines, map)

def checkSaddles(map, saddles, flowlines = None):
    print "comparing saddles of %d edges.." % (map.edgeCount, ),
    valid = True
    for edge in map.edgeIter():
        for saddle in map.wsStats.edgeSaddles(edge):
            if not saddle in saddles:
                print "illegal saddle %s in %s" % (saddle, edge)
                valid = False
                if flowlines:
                    originalPoly = flowlines[edge.label()][2]
                    originalSaddleIndex = flowlines[edge.label()][3]
                    originalSaddle = flowlines[edge.label()][4]
                    print "  original saddle at index %d/%d: %s = %s" % (
                        originalSaddleIndex, len(originalPoly),
                        originalSaddle, saddles[originalSaddle])
                print "  indices now: %s" % (
                    map.wsStats._indices[edge.label()], )
    if valid:
        print "all OK"
    return valid

def checkSplittingSaddleHandling(edge, offset):
    saddles = map.wsStats.edgeSaddles(edge.label())
    saddleIndex = map.wsStats._indices[edge.label()][0]
    newEdge = map.splitEdge(edge, saddleIndex + offset)
    if offset >= 0:
        assert map.wsStats.edgeSaddles(edge) == saddles
    if offset <= 0:
        assert map.wsStats.edgeSaddles(newEdge) == saddles
    map.mergeEdges(newEdge.dart())
    assert map.wsStats.edgeSaddles(edge.label()) == saddles
    assert map.wsStats._indices[edge.label()][0] == saddleIndex

from flag_constants import BORDER_PROTECTION

def checkPassValues(map, siv):
    for edge in map.edgeIter():
        if edge.flag(BORDER_PROTECTION):
            continue
        pv = map.wsStats.passValue(edge)
        actualPV = min([siv[p] for p in edge])
        assert pv == actualPV, "%s != %s" % (pv, actualPV) # no epsilon needed ;-)

import statistics
map.wsStats = statistics.WatershedStatistics(map, flowlines, img.gm.siv)
assert checkSaddles(map, img.spws.saddles(), flowlines)
checkPassValues(map, img.gm.siv)

innerEdges = [edge for edge in map.edgeIter()
              if not edge.flag(BORDER_PROTECTION)]
longEdges = [edge for edge in innerEdges
             if len(edge) > 50 and not edge.isLoop()]
longLoops = [edge for edge in innerEdges
             if len(edge) > 50 and edge.isLoop()]

print "checking saddle index handling during splits..."
checkSplittingSaddleHandling(longEdges[0], 5)
checkSplittingSaddleHandling(longEdges[1], -5)
checkSplittingSaddleHandling(longEdges[2], 0)
checkSplittingSaddleHandling(longLoops[0], 5)
checkSplittingSaddleHandling(longLoops[1], -5)
checkSplittingSaddleHandling(longLoops[2], 0)
print "  OK"

assert checkSaddles(map, img.spws.saddles(), flowlines)

print "sigma sorting & splitting parallel edges..."
map.sortEdgesEventually(stepDist = 0.2, minDist = 0.05)
map.splitParallelEdges()

assert checkSaddles(map, img.spws.saddles(), flowlines)
checkPassValues(map, img.gm.siv)

# --------------------------------------------------------------------
# 			complete subpixel watersheds map + statistics
# --------------------------------------------------------------------

spmap = maputils.subpixelWatershedMap(
    maxima, flowlines, img.size(),
    wsStatsSpline = img.gm.siv,
    minima = img.spws.minima())

maputils.removeCruft(spmap, 7)
assert checkSaddles(spmap, img.spws.saddles(), flowlines)
checkPassValues(spmap, img.gm.siv)

# --------------------------------------------------------------------

import pixelmap
lab, count = labelImage4(spmap.labelImage())
cm = pixelmap.crackEdgeMap(lab)
assert cm.faceCount == count + 1

# --------------------------------------------------------------------

map = GeoMap(maxima1, [], Size2D(256, 256))
maputils.addFlowLinesToMap(flowlines1, map)
maputils.connectBorderNodes(map, 0.1)
map.sortEdgesEventually(stepDist = 0.2, minDist = 0.05)
map.initializeMap()
assert map.checkConsistency(), "map inconsistent"
assert maputils.checkLabelConsistency(map), "map.labelImage() inconsistent"

# showMapStats(map)
# bg = readImage("../../../Testimages/blox.gif")
# d = MapDisplay(bg, map)

labels = [map.faceAt(p).label() for p in
          [Vector2(91, 86.4), Vector2(91, 85.8), Vector2(91.4, 85.8)]]
assert len(dict.fromkeys(labels).keys()) == 3, \
       "those nearby positions point to three pairwise different faces"

# --------------------------------------------------------------------

def checkPolygons(map):
    clean = True
    for edge in map.edgeIter():
        l, pa = edge.length(), edge.partialArea()
        p = Polygon(edge)
        p.invalidateProperties()
        if abs(l - p.length()) > 1e-6:
            print "edge %d: length wrong (was %s, is %s)" % (
                edge.label(), l, p.length())
            clean = False
        if abs(pa - p.partialArea()) > 1e-6:
            print "edge %d: partial area wrong (was %s, is %s)" % (
                edge.label(), pa, p.partialArea())
            clean = False
    return clean

# --------------------------------------------------------------------

import random, time, sys
if len(sys.argv) > 1:
    seed = long(sys.argv[1])
else:
    seed = long(time.time())

print "using %d as seed." % (seed, )
random.seed(seed)
#print "random state:", random.getstate()

def mergeEdgesCandidates(map):
    result = []
    for node in map.nodeIter():
        if node.degree() == 2 and not node.anchor().edge().isLoop():
            result.append(node.label())
    return result

def removeBridgeCandidates(map):
    result = []
    for edge in map.edgeIter():
        if edge.leftFaceLabel() == edge.rightFaceLabel():
            result.append(edge.label())
    return result

def mergeFacesCandidates(map):
    result = []
    for edge in map.edgeIter():
        if edge.leftFaceLabel() != edge.rightFaceLabel():
            result.append(edge.label())
    return result

for node in map.nodeIter():
    if node.degree() == 0:
        map.removeIsolatedNode(node)

backup = copy.copy(map)

history = maputils.LiveHistory(map)
try:
  changed = True
  while True:
    assert map.checkConsistency()

    if changed:
        possible = range(3)
        mec = rbc = mfc = None

    if not len(possible):
        break

    operation = random.choice(possible)
    changed = False

    if operation == 0:
        if mec == None:
            mec = mergeEdgesCandidates(map)
        if not len(mec):
            possible.remove(0)
        else:
            nodeLabel = random.choice(mec)
            mec.remove(nodeLabel)
            dart = map.node(nodeLabel).anchor()
            print "removing node %d via %s" % (nodeLabel, dart)
#             mergePos = dart[0]
            survivor = map.mergeEdges(dart)
            if survivor:
#                 assert mergePos in [survivor[i] for i in survivor.mergeIndices], \
#                        "mergeIndices do not point to merge position!"
                changed = True

    if operation == 1:
        if rbc == None:
            rbc = removeBridgeCandidates(map)
        if not len(rbc):
            possible.remove(1)
        else:
            dartLabel = random.choice(rbc)
            rbc.remove(dartLabel)
            dart = map.dart(dartLabel)
            print "removing bridge via %s" % (dart, )
            if map.removeBridge(dart):
                changed = True

    if operation == 2:
        if mfc == None:
            mfc = mergeFacesCandidates(map)
        if not len(mfc):
            possible.remove(2)
        else:
            dartLabel = random.choice(mfc)
            mfc.remove(dartLabel)
            dart = map.dart(dartLabel)
            print "removing edge via %s" % (dart, )
            if map.mergeFaces(dart):
                changed = True

except Exception, e:
    print history.commands()
    raise

maputils.checkLabelConsistency(map)

print "finished successfully, no more candidates for Euler ops:"
print map.nodeCount, "nodes:", list(map.nodeIter())
print map.edgeCount, "edges:", list(map.edgeIter())
print map.faceCount, "faces:", list(map.faceIter())

# --------------------------------------------------------------------

map = backup

history = history.commands()
c = time.clock()
exec history

print "replaying history took %ss." % (time.clock() - c, )
