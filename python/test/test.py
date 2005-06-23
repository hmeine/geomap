execfile("map.py")
execfile("testSPWS")

map = Map(maxima, flowlines, Size2D(256, 256))

def checkConsistency(map):
    nodeCount = 0
    edgeCount = 0
    faceCount = 0

    valid = True
    for node in map.nodeIter():
        nodeCount += 1
        assert map.nodes[node._label] == node
        assert node._map == map
        for a in node._darts:
            if not map.edges[abs(a)]:
                sys.stderr.write("%s contains invalid anchor %d!\n" % (node, a))
                valid = False

    for edge in map.edgeIter():
        edgeCount += 1
        assert map.edges[edge._label] == edge
        assert edge._map == map
        try:
            assert edge._label in map.node(edge._startNodeLabel)._darts
            assert -edge._label in map.node(edge._endNodeLabel)._darts
        except:
            sys.stderr.write("%s has invalid ends!\n  start node: %s\n    end node: %s\n" % (edge, map.node(edge._startNodeLabel), map.node(edge._endNodeLabel)))
            valid = False

    for face in map.faceIter():
        faceCount += 1
        assert map.faces[face._label] == face
        assert face._map == map
        for a in face._anchors:
            if a and not map.edges[a.edgeLabel()]:
                sys.stderr.write("%s contains invalid anchor %d!\n" % (face, a.label()))
                valid = False

    if nodeCount != map.nodeCount:
        sys.stderr.write("Map's nodeCount (%d) is wrong: counted %d!\n" %
                         (map.nodeCount, nodeCount))
        valid = False
    if edgeCount != map.edgeCount:
        sys.stderr.write("Map's edgeCount (%d) is wrong: counted %d!\n" %
                         (map.edgeCount, edgeCount))
        valid = False
    if faceCount != map.faceCount:
        sys.stderr.write("Map's faceCount (%d) is wrong: counted %d!\n" %
                         (map.faceCount, faceCount))
        valid = False

    assert valid

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

checkConsistency(map)
checkLabelConsistency(map)

mergeEdges(map.dart(29)) # create loop with two degree-2-nodes
checkConsistency(map)
mergeEdges(map.dart(43)) # lead to infinite loop
checkConsistency(map)

removeBridge(map.dart(32)) # makes edge 22 a self-loop
try:
    mergeEdges(map.dart(22))
except TypeError:
    print "caught expected exception."
    pass

checkConsistency(map)
checkLabelConsistency(map)

# --------------------------------------------------------------------

points = [Vector(232.20846246994, 81.488755298170375),
          Vector(228.16750125077627, 81.481533365106415),
          Vector(224.94552025882538, 81.580691309124461)]
print scanPoly(iter(points), 2, 81)
assert scanPoly(iter(points), 2, 81) == [
    [(227.5675290883523, 0, 232.20846557617188)],
    [(224.94552612304688, 0, 227.5675290883523)]]

# --------------------------------------------------------------------

# labelImage contains 2 pixels with unknown faces!
#   unknown face labels found: [99999.0]
mergeFaces(map.dart(36))
mergeEdges(map.dart(-28))

checkLabelConsistency(map)

#sys.exit(0)

# --------------------------------------------------------------------

import random, time
if len(sys.argv) > 1:
    seed = long(sys.argv[1])
else:
    seed = time.time()

print "using %d as seed." % (seed, )
random.seed(seed)
print "random state:", random.getstate()

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
    checkConsistency(map)

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
