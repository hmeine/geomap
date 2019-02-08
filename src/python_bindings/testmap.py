##########################################################################
#
#                Copyright 2007-2019 by Hans Meine
#
#     Permission is hereby granted, free of charge, to any person
#     obtaining a copy of this software and associated documentation
#     files (the "Software"), to deal in the Software without
#     restriction, including without limitation the rights to use,
#     copy, modify, merge, publish, distribute, sublicense, and/or
#     sell copies of the Software, and to permit persons to whom the
#     Software is furnished to do so, subject to the following
#     conditions:
#
#     The above copyright notice and this permission notice shall be
#     included in all copies or substantial portions of the
#     Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#     OTHER DEALINGS IN THE SOFTWARE.
#
##########################################################################

import geomap

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
tm = geomap.GeoMap([], [], Size2D(256, 256))

points = geomap.Vector2Array([Vector2(232.20846246994, 81.488755298170375), Vector2(228.16750125077627, 81.481533365106415), Vector2(224.94552025882538, 81.580691309124461)])

n1 = tm.addNode(points[0])
n2 = tm.addNode(points[-1])
print "- endnodes %d and %d created." % (n1.label(), n2.label())

print "- adding edge..."
edge = tm.addEdge(n1, n2, points)

print "- testing DartPointIter..."
dart = tm.dart(-edge.label())
for i, p in enumerate(dart):
    assert p == dart[i]
    assert p == edge[-1-i]

print "- testing ContourPointIter..."
assert len(list(geomap.ContourPointIter(dart))) == 2*len(list(dart))-2
closed = list(geomap.ContourPointIter(dart, True))
assert len(closed) == 2*len(list(dart))-1
assert closed[0] == closed[-1]

print "- testing initContours()..."
tm.initializeMap()
try:
    tm.initializeMap()
    assert False, "contours must only be initializable once!"
except RuntimeError:
    pass

# --------------------------------------------------------------------

import random, geomap
from vigra import *
points = [None]
for i in range(100):
    points.append(Vector2(10*random.random(), 10*random.random()))

print "\n- creating Map from %d random points..." % (len(points)-1, )
gm = geomap.GeoMap(points, [], Size2D(11, 11))

showMapStats(gm)

# import mapdisplay
# d = mapdisplay.MapDisplay(gm)

for i in range(1, len(points)):
    assert gm.node(i).position() == points[i]

for node in gm.nodeIter():
    if not node.degree():
        gm.removeIsolatedNode(node) # check MapDisplay callbacks

import triangle
points, edgeData = triangle.delaunay(points[1:])

edgeTuples = [startEnd and
              (startEnd[0], startEnd[1],
               #([points[startEnd[0]], points[startEnd[1]]]))
              geomap.Vector2Array([points[startEnd[0]], points[startEnd[1]]]))
              for startEnd in edgeData]

print "\n- creating Map from delaunay data (%d edges)..." % (len(edgeTuples)-1, )
gm = geomap.GeoMap(points, edgeTuples, Size2D(11, 11))

showMapStats(gm)

# addPathFromHere("../subpixelWatersheds")
# import mapdisplay
# d = mapdisplay.MapDisplay(gm)

# --------------------------------------------------------------------

print "\n- unpickling large.map"
import geomap, pickle
amap = pickle.load(file("large.map"))
#data = list(amap.__getstate__()[:3])
#data[1] = [poly and (poly[0], poly[1], geomap.Vector2Array(poly[2])) for poly in data[1]]
print "- initializing GeoMap with that data..."
cppmap = geomap.GeoMap(*amap.__getstate__()[:3])

# --------------------------------------------------------------------
#          run the following from within ../subpixelWatersheds
# --------------------------------------------------------------------

import map as spmap, maputils

print "\n- creating Map from maxima2/flowlines2..."
execfile("testSPWS")
pythonMap = spmap.GeoMap(maxima2, flowlines2, Size2D(39, 39))
pythonMap.sortEdgesEventually(0.2, 0.05)
pythonMap.initializeMap()

if pythonMap.unsortable:
    sys.stderr.write("### UNSORTABLE HANDLING NOT DONE IN C++ YET, SKIPPING...\n")
else:
    cppmap = geomap.GeoMap([], [], Size2D(39, 39))
    cppnodes = [node and cppmap.addNode(node) for node in maxima2]
    assert cppnodes[-1].label() == len(maxima2)-1, "Oops, label shift :-("
    maputils.addFlowLinesToMap(flowlines2, cppmap)
    cppmap.sortEdgesEventually(0.2, 0.05)
    cppmap.initializeMap()

# --------------------------------------------------------------------

print "\n- checking C++ sigma sorting against python one..."
execfile("maptest.py")
e = Experiment(filename, "grad")
gradientScale = 1.4
threshold = 0.05
e.performEdgeSplits = False
e.performBorderClosing = False
e("map")
checkConsistency(e.map)

def recreateWithCPPMap(pythonMap):
    cppmap = geomap.GeoMap([], [], pythonMap.imageSize())
    cppnodes = [node and cppmap.addNode(node.position()) for node in pythonMap.nodes]
    for edge in pythonMap.edgeIter():
        cppmap.addEdge(cppnodes[edge.startNodeLabel()],
                       cppnodes[edge.endNodeLabel()], edge,
                       edge.label())
    cppmap.sortEdgesEventually(0.2, 0.05)
    cppmap.initializeMap()
    return cppmap

print "- creating C++ map from same data...",
c = time.clock(); sys.stdout.flush()
cppmap = recreateWithCPPMap(e.map)
print " (%ss)" % (time.clock() - c, )

# for edge in e.map.edgeIter():
#     o1 = edge.dart().sigmaOrbit()
#     o2 = cppmap.dart(edge.label()).sigmaOrbit()
#     for d1, d2 in map(None, o1, o2):
#         if d1.label() != d2.label():
#             print edge.startNode()

print "- checking sigma orbits of all %d nodes..." % (e.map.nodeCount, )
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
    eval("geomap.%s" % op)(cppmap.dart(dl))
    cppmap.checkConsistency()
