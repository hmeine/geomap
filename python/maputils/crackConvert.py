import sys, time
from hourglass import GeoMap, crackConnectionImage
from vigra import meshIter
from flag_constants import BORDER_PROTECTION
import maputils
import progress

__all__ = ["crackEdgeMap", "crackEdgeGraph"]

# FIXME: 8-connected regions with holes (see testcase)
# (details: normally, every upper left corner is considered to be a
# node, but a holes' upper left corner may be a 4-junction which has
# explicitly been marked as degree 2 position, such that no start/end
# for the hole exists)

# --------------------------------------------------------------------
#                              FRONTEND
# --------------------------------------------------------------------

def crackEdgeMap(labelImage, initLabelImage = True,
                 eightConnectedRegions = True):
    c = time.clock()
    msg = progress.StatusMessage("- following crack edges")
    result = crackEdgeGraph(
        labelImage, eightConnectedRegions = eightConnectedRegions,
        progressHook = progress.ProgressHook(msg))
    msg.finish()

    sys.stdout.write("  removing deg.2 nodes..."); c1 = time.clock()
    maputils.mergeDegree2Nodes(result)
    sys.stdout.write(" (%ss)\n" % (time.clock()-c1, ))

    sys.stdout.write("  sorting edges..."); c1 = time.clock()
    result.sortEdgesDirectly()
    sys.stdout.write(" (%ss)\n" % (time.clock()-c1, ))

    sys.stdout.write("  initializing faces..."); c1 = time.clock()
    result.initializeMap(initLabelImage)
    sys.stdout.write(" (%ss)\n" % (time.clock()-c1, ))
    
    # mark the border edges:
    assert result.face(0).holeCount() == 1, "infinite face should have exactly one contour, not %d!?" % result.face(0).holeCount()
    for dart in result.face(0).holeContours().next().phiOrbit():
        edge = dart.edge()
        if not edge.leftFaceLabel() or not edge.rightFaceLabel():
            edge.setFlag(BORDER_PROTECTION)

    sys.stdout.write("  done. (%ss)\n" % (time.clock()-c, ))

    return result

# --------------------------------------------------------------------
#                          helper functions
# --------------------------------------------------------------------

from vigra import GrayImage, Vector2, Point2D, Size2D, Rect2D
from hourglass import Polygon

CONN_RIGHT = 1
CONN_DOWN = 2
CONN_LEFT = 4
CONN_UP = 8

connections = [CONN_RIGHT, CONN_UP, CONN_LEFT, CONN_DOWN]

degree = [0] * 16
for conn in connections:
    for i in range(16):
        if i & conn:
            degree[i] += 1

# extend degree LUT for special handling of 8-connected regions:
degree = degree + degree
degree[31] = 2
degree = degree + degree
degree[47] = 2

def pyCrackConnectionImage(labelImage):
    result = GrayImage(labelImage.size()+Size2D(1,1))

    for y in range(labelImage.height()-1):
        for x in range(labelImage.width()-1):
            if labelImage[x, y] != labelImage[x+1, y]:
                result[x+1, y  ] += CONN_DOWN
                result[x+1, y+1] += CONN_UP
            if labelImage[x, y] != labelImage[x, y+1]:
                result[x,   y+1] += CONN_RIGHT
                result[x+1, y+1] += CONN_LEFT
        x = labelImage.width()-1
        if labelImage[x, y] != labelImage[x, y+1]:
            result[x,   y+1] += CONN_RIGHT
            result[x+1, y+1] += CONN_LEFT

    lastRow = labelImage.subImage((0, labelImage.height()-1),
                                  Size2D(labelImage.width(), 1))
    for x in range(labelImage.width()-1):
        if lastRow[x, 0] != lastRow[x+1, 0]:
            result[x+1, labelImage.height()-1] += CONN_DOWN
            result[x+1, labelImage.height()  ] += CONN_UP

    # add border:
    for x in range(labelImage.width()):
        result[x  , 0]                   += CONN_RIGHT
        result[x  , labelImage.height()] += CONN_RIGHT
        result[x+1, 0]                   += CONN_LEFT
        result[x+1, labelImage.height()] += CONN_LEFT

    for y in range(labelImage.height()):
        result[0,                  y  ] += CONN_DOWN
        result[labelImage.width(), y  ] += CONN_DOWN
        result[0,                  y+1] += CONN_UP
        result[labelImage.width(), y+1] += CONN_UP

    return result

_dirOffset = [Size2D(1, 0), Size2D(0, -1), Size2D(-1, 0), Size2D(0, 1)]
_dirVector = [Vector2(1, 0), Vector2(0, -1), Vector2(-1, 0), Vector2(0, 1)]
_turnRight = [3, 0, 1, 2]
_turnLeft  = [1, 2, 3, 0]

# for extracting loops, we add nodes at every upper left corner:
def isNode(connValue):
    return degree[connValue] > 2 or connValue == (CONN_RIGHT | CONN_DOWN)

def followEdge(crackConnectionImage, pos, direction):
    """Follow edge starting at `pos` in `direction` until
    crackConnectionImage[pos] fulfills `isNode`."""
    pos = Point2D(pos[0], pos[1])
    vPos = Vector2(pos[0] - 0.5, pos[1] - 0.5)
    result = Polygon([vPos])

    while True:
        vPos += _dirVector[direction]
        result.append(vPos)
        pos += _dirOffset[direction]

        connection = int(crackConnectionImage[pos])
        if connection >= 16:
            if connection & 16 and direction in (1, 3):
                direction = _turnLeft[direction]
                continue
            elif connection & 32 and direction in (0, 2):
                direction = _turnLeft[direction]
                continue
        elif isNode(connection):
            break

        direction = _turnRight[direction]
        while connection & connections[direction] == 0:
            direction = _turnLeft[direction]

    return result, pos, connections[(direction+2)%4]

def crackEdgeGraph(labelImage, eightConnectedRegions = True,
                   progressHook = None):
    result = GeoMap(labelImage.size())

    cc = crackConnectionImage(labelImage)

    if eightConnectedRegions:
        for y in range(1, cc.height()-1):
          for x in range(1, cc.width()-1):
            if cc[x,y] == 15:
                if labelImage[x,y] == labelImage[x-1,y-1]:
                    cc[x,y] += 16
                if labelImage[x-1,y] == labelImage[x,y-1]:
                    cc[x,y] += 32
                if cc[x,y] == 15+16+32: # crossing regions?
                    if labelImage[x,y-1] > labelImage[x-1,y-1]:
                        cc[x,y] -= 16
                    else:
                        cc[x,y] -= 32
    
    nodeImage = GrayImage(cc.size())

    progressHook = progressHook and progressHook.rangeTicker(cc.height())

    for y in range(cc.height()):
      if progressHook:
          progressHook()
      for x in range(cc.width()):
        nodeConn = int(cc[x, y])
        if isNode(nodeConn):
            startNodeInfo = int(nodeImage[x, y])
            if startNodeInfo:
                startNode = result.node(startNodeInfo >> 4)
            else:
                startNode = result.addNode((x - 0.5, y - 0.5))
                nodeImage[x, y] = startNodeInfo = startNode.label() << 4

            for direction, startConn in enumerate(connections):
                if nodeConn & startConn and not startNodeInfo & startConn:
                    edge, endPos, endConn = followEdge(
                        cc, (x, y), direction)
                    endNodeInfo = int(nodeImage[endPos])
                    if not endNodeInfo:
                        endNode = result.addNode((endPos[0] - 0.5, endPos[1] - 0.5))
                        endNodeInfo = endNode.label() << 4
                    else:
                        assert not endNodeInfo & endConn, "double connection?"
                        endNode = result.node(endNodeInfo >> 4)

                    edge = result.addEdge(startNode, endNode, edge)

                    startNodeInfo |= startConn
                    if edge.isLoop():
                        startNodeInfo |= endConn
                        nodeImage[x, y] = startNodeInfo
                    else:
                        nodeImage[x, y] = startNodeInfo
                        nodeImage[endPos] = endNodeInfo | endConn

    return result

def showDegrees(crackConnectionImage):
    degreeImage = GrayImage(crackConnectionImage.size())
    for p in crackConnectionImage.size():
        degreeImage[p] = degree[int(crackConnectionImage[p])]
    return showImage(degreeImage)

if __name__ == "__main__":
    import unittest

    class CrackEdgeMapTest(unittest.TestCase):
        def testSimpleMap(self):
            g = GrayImage(10, 10)
            for y in range(5, 10):
                g[5,y] = 1
            cem = crackEdgeMap(g)
            self.assertEqual(cem.faceCount, 3) # should be one infinite, one background, and one foreground region

        def test8Connected(self):
            g = GrayImage(10, 10)
            for xy in range(3, 7):
                g[xy,xy] = 1
            cem = crackEdgeMap(g)
            self.assertEqual(cem.faceCount, 3) # should be one infinite, one background, and one foreground region

        def test8ConnectedOpposite(self):
            g = GrayImage(10, 10)
            for xy in range(3, 7):
                g[8-xy,xy] = 1
            cem = crackEdgeMap(g)
            self.assertEqual(cem.faceCount, 3) # should be one infinite, one background, and one foreground region

        def test8ConnectedLoop(self):
            g = GrayImage(10, 10)
            g[4,4] = 1
            g[5,5] = 1
            g[4,6] = 1
            g[3,5] = 1
            cem = crackEdgeMap(g)
            self.assertEqual(cem.faceCount, 4) # as above, plus one hole

        def testHeadlineMap(self):
            from vigra import readImage
            labels = readImage("headline.png")[0]
            cem = crackEdgeMap(labels)
            self.assertEqual(cem.nodeCount, 11)
            self.assertEqual(cem.edgeCount, 13)
            self.assertEqual(cem.faceCount, 11)

    sys.exit(unittest.main())
