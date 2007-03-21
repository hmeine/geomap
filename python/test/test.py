from vigra import Vector2, Size2D, labelImage4
from hourglass import Polygon
import maputils

from hourglass import GeoMap
#from map import GeoMap
execfile("testSPWS")

# The following data contains edges that run out of the image range,
# which ATM leads to overlapping edges after border closing. That
# violates some assumptions and would lead to errors if we actually
# worked with that Map.
map = GeoMap(maxima2, [], Size2D(39, 39))
maputils.addFlowLinesToMap(flowlines2, map)
map.sortEdgesEventually(stepDist = 0.2, minDist = 0.05)
map.splitParallelEdges()
map.initializeMap()
assert map.checkConsistency(), "map inconsistent"
assert maputils.checkLabelConsistency(map), "map.labelImage() inconsistent"

# merge faces so that survivor has a hole:
hole = map.faceAt((16,17))
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

# test copying / init from map data with disconnected contours:
import copy
om = copy.copy(map)
assert om.__getstate__() == map.__getstate__()
assert om.checkConsistency(), "map inconsistent"
assert maputils.checkLabelConsistency(om), "om.labelImage() inconsistent"

# --------------------------------------------------------------------

def extractAllSaddles(map):
    result = []
    for edge in map.edgeIter():
        result.extend(map.wsStats.edgeSaddles(edge))
    return result

from vigra import readImage, resizeImageSplineInterpolation, transformImage, gaussianGradientAsVector, SplineImageView5
from hourglass import SubPixelWatersheds5

img = resizeImageSplineInterpolation(readImage("test_blox.png"), Size2D(64, 64))
img.grad = gaussianGradientAsVector(img, 1.6)
img.gm = transformImage(img.grad, "\l x: norm(x)", {})
img.gm.siv = SplineImageView5(img.gm)
img.spws = SubPixelWatersheds5(img.gm)

import maputils
maxima, flowlines = maputils.subpixelWatershedData(
    img.spws,
    img.gm.siv, 0.7,
    perpendicularDistEpsilon = None)
spmap = maputils.subpixelWatershedMap(
    maxima, flowlines, img.size(),

    wsStatsSpline = img.gm.siv,
    minima = img.spws.minima())

print "comparing saddles..",
saddles1 = img.spws.saddles()
for saddle in extractAllSaddles(spmap):
    assert saddle in saddles1
print "all OK"

# --------------------------------------------------------------------

import pixelmap
lab, count = labelImage4(map.labelImage())
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

assert map.faceAt(Vector2(91,  86.4)) == map.face(14)
assert map.faceAt(Vector2(91,  85.8)) == map.face(6)
assert map.faceAt(Vector2(91.4,85.8)) == map.face(1)

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

history = ""
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
            history += "map.mergeEdges(map.dart(%d))\n" % (dart.label(), )
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
            history += "map.removeBridge(map.dart(%d))\n" % (dartLabel, )
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
            history += "map.mergeFaces(map.dart(%d))\n" % (dartLabel, )
            if map.mergeFaces(dart):
                changed = True

except Exception, e:
    print history
    raise

maputils.checkLabelConsistency(map)

print "finished successfully, no more candidates for Euler ops:"
print map.nodeCount, "nodes:", list(map.nodeIter())
print map.edgeCount, "edges:", list(map.edgeIter())
print map.faceCount, "faces:", list(map.faceIter())

# --------------------------------------------------------------------

map = backup

c = time.clock()
exec history

print "replaying history took %ss." % (time.clock() - c, )
