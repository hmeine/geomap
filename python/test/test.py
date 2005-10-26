from vigra import Vector2
from hourglass import Polygon, scanPoly

def scanPoly2List(*args):
    return [[(e.begin, e.direction, e.end) for e in list(s)]
            for s in scanPoly(*args)]

execfile("testPolygons")

assert scanPoly2List(miniPoly, 2, 81) == \
       [[(228, 1, 233)], [(225, 1, 229)]]
assert scanPoly2List(smallPoly, 5, 180) == \
       [[], [], [(66, 0, 68)], [], []]
assert scanPoly2List(openPoly, 11, 188) == \
       [[(177, 1, 178)], [(177, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 2, 179)], [(178, 1, 179)]]
assert scanPoly2List(closedPoly, 20, 104) == \
       [[(86, 0, 88)], [(85, 0, 89)], [(83, 2, 86), (88, -2, 89)], [(82, 2, 84), (87, -2, 89)], [(81, 2, 83), (87, -2, 88)], [(80, 2, 82), (87, -2, 88)], [(79, 2, 81), (87, -2, 88)], [(79, 2, 80), (87, -2, 88)], [(78, 2, 80), (87, -2, 88)], [(78, 2, 79), (87, -2, 88)], [(78, 2, 79), (86, -2, 88)], [(77, 2, 79), (86, -2, 87)], [(77, 2, 79), (86, -2, 87)], [(78, 2, 79), (86, -2, 87)], [(78, 2, 79), (86, -2, 88)], [(78, 2, 79), (87, -2, 88)], [(78, 2, 79), (87, -2, 88)], [(78, 2, 79), (81, -2, 88)], [(78, 0, 82)], []]
assert scanPoly2List(Polygon([
    Vector2(4.5, 0.5), Vector2(0.5, 0.5),
    Vector2(0.5, 4.5), Vector2(2.5, 4.5),
    Vector2(2.5, 2.5), Vector2(4.5, 2.5),
    Vector2(4.5, 0.5)]), 7) == \
    [[], [(1, 2, 1), (5, -2, 5)], [(1, 2, 1), (5, -2, 5)], [(1, 2, 1), (3, -2, 3)], [(1, 2, 1), (3, -2, 3)], [], []]

print "scanPoly tested OK."

# --------------------------------------------------------------------

p1 = Polygon([Vector2(1, 0), Vector2(2, 0)])
p2 = Polygon([Vector2(2, 0), Vector2(3, 0)])
p1.extend(p2)
assert len(p1) == 3, "Polygon composition should prevent duplicate points:\n  %s" % (
    list(p1), )

# --------------------------------------------------------------------

execfile("map.py")
execfile("testSPWS")

map = Map(maxima, [fl and fl[0] or None for fl in flowlines], Size2D(256, 256))

class FaceLookup:
    def __init__(self, map):
        self._map = map
        self.errorCount = 0
        self.errorLabels = []

    def __call__(self, faceLabel):
        if faceLabel < 0:
            return
        faceLabel = int(faceLabel)
        try:
            assert self._map.face(faceLabel) != None
        except:
            self.errorCount += 1
            if not faceLabel in self.errorLabels:
                self.errorLabels.append(faceLabel)

def checkLabelConsistency(map):
    #print "checkLabelConsistency skipped (segfaults)."
    #return
    fl = FaceLookup(map)
    inspectImage(map.labelImage, fl)
    if fl.errorCount:
        sys.stderr.write("labelImage contains %d pixels with unknown faces!\n" % (
            fl.errorCount, ))
        sys.stderr.write("  unknown face labels found: %s\n" % (fl.errorLabels, ))
        if fl.errorCount < 40:
            for p in map.labelImage.size():
                if int(map.labelImage[p]) in fl.errorLabels:
                    print "   label %d at %s" % (int(map.labelImage[p]), p)
    assert fl.errorCount == 0

# execfile("maptest.py")
# showMapStats(map)
# bg = readImage("../../../Testimages/blox.gif")
# d = MapDisplay(bg, map)

assert checkConsistency(map)
checkLabelConsistency(map)

# --------------------------------------------------------------------

import random, time
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
possible = range(3)
try:
  while True:
    assert checkConsistency(map)

    if not len(possible):
        break

    operation = random.choice(possible)
    if operation == 0:
        mec = mergeEdgesCandidates(map)
        if not len(mec):
            possible.remove(0)
        else:
            dart = map.dart(map.node(random.choice(mec))._darts[0])
            print "removing node %d via dart %s" % (dart.startNodeLabel(), dart)
            history += "mergeEdges(map.dart(%d))\n" % (dart.label(), )
            mergeEdges(dart)
            possible = range(3)

    if operation == 1:
        rbc = removeBridgeCandidates(map)
        if not len(rbc):
            possible.remove(1)
        else:
            dart = map.dart(random.choice(rbc))
            print "removing bridge via dart %s" % (dart, )
            history += "removeBridge(map.dart(%d))\n" % (dart.label(), )
            removeBridge(dart)
            possible = range(3)

    if operation == 2:
        mfc = mergeFacesCandidates(map)
        if not len(mfc):
            possible.remove(2)
        else:
            dart = map.dart(random.choice(mfc))
            print "removing edge via dart %s" % (dart, )
            history += "mergeFaces(map.dart(%d))\n" % (dart.label(), )
            mergeFaces(dart)
            possible = range(3)

except Exception, e:
    print history
    raise

checkLabelConsistency(map)

print "finished successfully, no more candidates for Euler ops:"
print map.nodeCount, "nodes:", list(map.nodeIter())
print map.edgeCount, "edges:", list(map.edgeIter())
print map.faceCount, "faces:", list(map.faceIter())
