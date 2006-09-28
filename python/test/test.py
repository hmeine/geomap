from vigra import Vector2
from hourglass import Polygon

from map import *
execfile("testSPWS")

# The following data contains edges that run out of the image range,
# which ATM leads to overlapping edges after border closing. That
# violates some assumptions and would lead to errors if we actually
# worked with that Map.
map = Map(maxima2, flowlines2, Size2D(39, 39))
assert checkConsistency(map), "map inconsistent"
assert checkLabelConsistency(map), "map.labelImage inconsistent"
assert len(map.history) == 0

# merge faces so that survivor has a hole:
mergeFaces(map.dart(213))
mergeFaces(map.dart(19))
mergeFaces(map.dart(28))
mergeFaces(map.dart(33))
mergeFaces(map.dart(225))
mergeFaces(map.dart(31))
mergeFaces(map.dart(16))
face = removeBridge(map.dart(18))

assert len(face.contours()) > 1 # should have hole
assert face.contains(Vector2(5,12)) # in region, but not within hole
assert not face.contains(Vector2(12,12)) # in hole

# --------------------------------------------------------------------

# test copying / init from map data with disconnected contours:
om = copy.copy(map)
assert checkConsistency(om), "map inconsistent"
assert checkLabelConsistency(om), "om.labelImage inconsistent"

# --------------------------------------------------------------------

import pixelmap
lab, count = labelImage4(map.labelImage)
cm = pixelmap.crackEdgeMap(lab)
assert cm.faceCount == count + 1

# --------------------------------------------------------------------

map = Map(maxima1, flowlines1, Size2D(256, 256))
assert checkConsistency(map), "map inconsistent"
assert checkLabelConsistency(map), "map.labelImage inconsistent"

# showMapStats(map)
# bg = readImage("../../../Testimages/blox.gif")
# d = MapDisplay(bg, map)

assert map.faceAt(Vector2(91,  86.4)) == map.face(13)
assert map.faceAt(Vector2(91,  85.8)) == map.face(5)
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
                edge._label, l, p.length())
            clean = False
        if abs(pa - p.partialArea()) > 1e-6:
            print "edge %d: partial area wrong (was %s, is %s)" % (
                edge._label, pa, p.partialArea())
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
        if node.degree() == 2 \
               and abs(node._darts[0]) != abs(node._darts[1]):
            result.append(node._label)
    return result

def removeBridgeCandidates(map):
    result = []
    for edge in map.edgeIter():
        if edge.leftFaceLabel() == edge.rightFaceLabel():
            result.append(edge._label)
    return result

def mergeFacesCandidates(map):
    result = []
    for edge in map.edgeIter():
        if edge.leftFaceLabel() != edge.rightFaceLabel():
            result.append(edge._label)
    return result

for node in map.nodeIter():
    if node.degree() == 0:
        node.uninitialize()

history = ""
try:
  changed = True
  while True:
    assert checkConsistency(map)

    if changed:
        possible = range(3)
        mec = rbc = mfc = None

    if not len(possible):
        break

    operation = random.choice(possible)
    changed = False
    try:

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
            history += "mergeEdges(map.dart(%d))\n" % (dart.label(), )
            mergePos = dart[0]
            survivor = mergeEdges(dart)
            assert mergePos in [survivor[i] for i in survivor.mergeIndices], \
                   "mergeIndices do not point to merge position!"
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
            history += "removeBridge(map.dart(%d))\n" % (dartLabel, )
            removeBridge(dart)
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
            history += "mergeFaces(map.dart(%d))\n" % (dartLabel, )
            mergeFaces(dart)
            changed = True

    except CancelOperation:
        print "-> cancelled"
        pass # border can't be removed..

except Exception, e:
    print history
    raise

checkLabelConsistency(map)

print "finished successfully, no more candidates for Euler ops:"
print map.nodeCount, "nodes:", list(map.nodeIter())
print map.edgeCount, "edges:", list(map.edgeIter())
print map.faceCount, "faces:", list(map.faceIter())
