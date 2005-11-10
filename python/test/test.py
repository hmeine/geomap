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

print "basic scanPoly tested OK."

# --------------------------------------------------------------------

from operator import setitem
def distinct(l):
    d = {}
    map(setitem, (d,)*len(l), l, [])
    return d.keys()

def points(scanlines):
    result = []
    for i, scanline in enumerate(scanlines):
        for segment in scanline:
            for x in range(segment.begin, segment.end):
                result.append((x, i+scanlines.startIndex))
    return result

p1 = Polygon([Vector2(95.5, 49.5), Vector2(97, 50), Vector2(98, 51), Vector2(99, 50), Vector2(100, 50), Vector2(101, 50), Vector2(102, 50), Vector2(103, 50), Vector2(104, 49)])
p2 = Polygon([Vector2(101.5, 46.5), Vector2(103, 48), Vector2(104, 49)])

r1 = points(scanPoly(p1, 20, 40))
r1.extend(points(scanPoly(p2, 60)))
r1 = distinct(r1)

p2.reverse()
p1.extend(p2)
r2 = points(scanPoly(p1, 80))
r2 = distinct(r2)

assert r1 == r2, "composition / reversing changes scanlines?!"

p1.reverse()
r2 = points(scanPoly(p1, 80))
r2 = distinct(r2)

assert r1 == r2, "reversing changes scanlines?!"

print "scanPoly composition/reversability tested OK."

# --------------------------------------------------------------------

p1 = Polygon([Vector2(1, 0), Vector2(2, 0)])
p2 = Polygon([Vector2(2, 0), Vector2(3, 0)])
p1.extend(p2)
assert len(p1) == 3, "Polygon composition should prevent duplicate points:\n  %s" % (
    list(p1), )

# --------------------------------------------------------------------

execfile("map.py")
execfile("testSPWS")

for maxima, flowlines, size in [
    (maxima2, flowlines2, Size2D(39, 39)),
    (maxima1, flowlines1, Size2D(256, 256))]:
    map = Map(maxima, flowlines, size)

    assert checkConsistency(map), "map inconsistent"
    assert checkLabelConsistency(map), "map.labelImage inconsistent"

# execfile("maptest.py")
# showMapStats(map)
# bg = readImage("../../../Testimages/blox.gif")
# d = MapDisplay(bg, map)

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
