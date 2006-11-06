import hourglass

from vigra import Vector2, Size2D, addPathFromHere

def showMapStats(map):
    # TODO: also display _edgeSplitGroups statistics
    pointCount = 0
    totalLength = 0.0
    for edge in map.edgeIter():
        pointCount += len(edge)
        totalLength += edge.length()

    if map.edgeCount == 0:
        print "empty map (no edges, %d nodes)!" % (map.nodeCount, )
        return

    print ("%d nodes, %d edges with a total of %d points " +
          "(mean #: %.2f, mean distance: %.2f), " +
           "%d faces") % (
        map.nodeCount, map.edgeCount, pointCount,
        float(pointCount) / map.edgeCount, totalLength / (pointCount - map.edgeCount),
        map.faceCount)
#     if map.deleted:
#         print "%d edges were deleted (e.g. at image border)." % (
#             len(map.deleted), )
#     if hasattr(map, "unsortable") and map.unsortable:
#         print "%d unsortable groups of edges occured." % (
#             len(map.unsortable), )
#     if hasattr(map, "unembeddableContours") and map.unembeddableContours:
#         print "%d outer contours could not be embedded into " \
#               "their surrounding faces!" % (len(map.unembeddableContours), )

# --------------------------------------------------------------------

print "- creating empty GeoMap..."
tm = hourglass.GeoMap([], [], Size2D(256, 256))

points = hourglass.Vector2Array([Vector2(232.20846246994, 81.488755298170375), Vector2(228.16750125077627, 81.481533365106415), Vector2(224.94552025882538, 81.580691309124461)])

n1 = tm.addNode(points[0])
n2 = tm.addNode(points[-1])
print "- endnodes %d and %d created." % (n1.label(), n2.label())

print "- adding edge..."
edge = tm.addEdge(n1.label(), n2.label(), points)

print "- testing DartPointIter..."
dart = tm.dart(-edge.label())
for i, p in enumerate(dart):
    assert p == dart[i]
    assert p == edge[-1-i]

print "- testing ContourPointIter..."
assert len(list(hourglass.ContourPointIter(dart))) == 2*len(list(dart))-2
closed = list(hourglass.ContourPointIter(dart, True))
assert len(closed) == 2*len(list(dart))-1
assert closed[0] == closed[-1]

print "- testing initContours()..."
tm.initContours()
try:
    tm.initContours()
    assert False, "contours must only be initializable once!"
except RuntimeError:
    pass

# --------------------------------------------------------------------

import random, hourglass
from vigra import *
points = [None]
for i in range(100):
    points.append(Vector2(10*random.random(), 10*random.random()))

print "\n- creating Map from %d random points..." % (len(points)-1, )
gm = hourglass.GeoMap(points, [], Size2D(11, 11))

showMapStats(gm)

for i in range(len(points)-1):
    assert gm.node(i+1).position() == points[i+1]

import triangle
points, edgeData = triangle.delaunay(points[1:])

edgeTuples = [startEnd and
              (startEnd[0], startEnd[1],
               #([points[startEnd[0]], points[startEnd[1]]]))
              hourglass.Vector2Array([points[startEnd[0]], points[startEnd[1]]]))
              for startEnd in edgeData]

print "\n- creating Map from delaunay data (%d edges)..." % (len(edgeTuples)-1, )
gm = hourglass.GeoMap(points, edgeTuples, Size2D(11, 11))

showMapStats(gm)

# addPathFromHere("../subpixelWatersheds")
# import mapdisplay
# d = mapdisplay.MapDisplay(gm)

# --------------------------------------------------------------------

print "\n- unpickling large.map"
import hourglass, pickle
amap = pickle.load(file("large.map"))
#data = list(amap.__getstate__()[:3])
#data[1] = [poly and (poly[0], poly[1], hourglass.Vector2Array(poly[2])) for poly in data[1]]
print "- initializing GeoMap with that data..."
cppmap = hourglass.GeoMap(*amap.__getstate__()[:3])

# --------------------------------------------------------------------
#          run the following from within ../subpixelWatersheds
# --------------------------------------------------------------------

from map import Map, addFlowLinesToMap

print "\n- creating Map from maxima2/flowlines2..."
execfile("testSPWS")
pythonMap = Map(maxima2, flowlines2, Size2D(39, 39),
                performEdgeSplits = False, performBorderClosing = False)

if pythonMap.unsortable:
    sys.stderr.write("### UNSORTABLE HANDLING NOT DONE IN C++ YET, SKIPPING...\n")
else:
    cppmap = hourglass.GeoMap([], [], Size2D(39, 39))
    cppnodes = [node and cppmap.addNode(node) for node in maxima2]
    assert cppnodes[-1].label() == len(maxima2)-1, "Oops, label shift :-("
    addFlowLinesToMap(flowlines2, cppmap)
    cppmap.sortEdgesEventually(0.2, 0.05)
    cppmap.initContours()
    cppmap.embedFaces()

# --------------------------------------------------------------------

print
execfile("maptest.py")
e = Experiment(filename, "grad")
gradientScale = 1.4
threshold = 0.05
e.performEdgeSplits = False
e.performBorderClosing = False
e("map")
checkConsistency(e.map)

def recreateWithCPPMap(pythonMap):
    cppmap = hourglass.GeoMap([], [], pythonMap.imageSize())
    cppnodes = [node and cppmap.addNode(node.position()) for node in pythonMap.nodes]
    for edge in pythonMap.edgeIter():
        cppmap.addEdge(cppnodes[edge.startNodeLabel()].label(),
                       cppnodes[edge.endNodeLabel()].label(), edge,
                       edge.label())
    cppmap.sortEdgesEventually(0.2, 0.05)
    cppmap.initContours()
    cppmap.embedFaces()
    return cppmap

cppmap = recreateWithCPPMap(e.map)

# for edge in e.map.edgeIter():
#     o1 = edge.dart().sigmaOrbit()
#     o2 = cppmap.dart(edge.label()).sigmaOrbit()
#     for d1, d2 in map(None, o1, o2):
#         if d1.label() != d2.label():
#             print edge.startNode()

for node in e.map.nodeIter():
    cd = cppmap.dart(node.anchor().label())
    for pd in node.anchor().sigmaOrbit():
        if cd.label() != pd.label():
            print "DIFFERENCE FOUND:", node
            break
        _ = cd.nextSigma()

e.cleanup()
testHist = e.map.history

print "\n- replaying history (%d entries) on C++ map..." % (len(testHist), )
for op, dl in testHist:
    #print op, dl
    eval("hourglass.%s" % op)(cppmap.dart(dl))
    cppmap.checkConsistency()
