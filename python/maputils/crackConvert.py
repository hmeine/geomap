import sys, time, copy
from hourglass import GeoMap, crackConnectionImage, crackEdgeGraph
import vigra
from flag_constants import BORDER_PROTECTION
import maputils
import progress

__all__ = ["crackEdgeMap", "crackEdgeGraph"]

# --------------------------------------------------------------------
#                              FRONTEND
# --------------------------------------------------------------------

def crackEdgeMap(labelImage, initLabelImage = True,
                 eightConnectedRegions = True):
    result = crackEdgeGraph(labelImage, eightConnectedRegions = eightConnectedRegions)
    maputils.mergeDegree2Nodes(result)
    result.sortEdgesDirectly()
    result.initializeMap(initLabelImage)
    
    # mark the border edges:
    assert result.face(0).holeCount() == 1, "infinite face should have exactly one contour, not %d!?" % result.face(0).holeCount()
    for dart in result.face(0).holeContours().next().phiOrbit():
        edge = dart.edge()
        if not edge.leftFaceLabel() or not edge.rightFaceLabel():
            edge.setFlag(BORDER_PROTECTION)

    return result

# --------------------------------------------------------------------
#                          helper functions
# --------------------------------------------------------------------

from vigra import GrayImage, Vector2, Point2D, Size2D
from hourglass import Polygon

CONN_RIGHT = 1
CONN_DOWN = 2
CONN_LEFT = 4
CONN_UP = 8
CONN_ALL4 = 15
_debugConn = {
    1 : "right",
    2 : "down",
    4 : "left",
    8 : "up"}

CONN_DIAG_UPLEFT = 16
CONN_DIAG_UPRIGHT = 32
CONN_DIAG = CONN_DIAG_UPLEFT | CONN_DIAG_UPRIGHT

CONN_NODE = 64
CONN_MAYBE_NODE = 128
CONN_ANYNODE = CONN_NODE | CONN_MAYBE_NODE

connections = [CONN_RIGHT, CONN_UP, CONN_LEFT, CONN_DOWN]

degree = [0] * 16
for conn in connections:
    for i in range(16):
        if i & conn:
            degree[i] += 1

# extend degree LUT for special handling of 8-connected regions:
degree = degree + degree
degree[CONN_ALL4 + CONN_DIAG_UPLEFT] = 2
degree = degree + degree
degree[CONN_ALL4 + CONN_DIAG_UPRIGHT] = 2
degree = degree + degree
degree = degree + degree

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

DIR_EAST   = 0
DIR_NORTH  = 1
DIR_WEST   = 2
DIR_SOUTH  = 3
_dirOffset = [Size2D(1, 0), Size2D(0, -1), Size2D(-1, 0), Size2D(0, 1)]
_dirVector = [Vector2(1, 0), Vector2(0, -1), Vector2(-1, 0), Vector2(0, 1)]
_turnRight = [3, 0, 1, 2]
_turnLeft  = [1, 2, 3, 0]
_debugDir  = ("east", "north", "west", "south")

def followEdge(crackConnectionImage, pos, direction):
    """Follow edge starting at `pos` in `direction` until
    crackConnectionImage[pos] has the CONN_NODE bit set, or until we
    arrive at `pos` again (self-loop).  Any CONN_MAYBE_NODE bits along
    the way are cleared."""
    pos = Point2D(pos[0], pos[1])
    vPos = Vector2(pos[0] - 0.5, pos[1] - 0.5)
    result = Polygon([vPos])

    startPos = copy.copy(pos)
    while True:
        vPos += _dirVector[direction]
        result.append(vPos)
        pos += _dirOffset[direction]

        if pos == startPos:
            break

        connection = int(crackConnectionImage[pos])
        if connection & CONN_DIAG:
            if connection & CONN_DIAG_UPLEFT:
                turnLeft = direction in (DIR_NORTH, DIR_SOUTH)
            else:
                turnLeft = direction in (DIR_EAST, DIR_WEST)

            connection &= ~connections[(direction+2)%4]

            if turnLeft:
                direction = _turnLeft[direction]
            else:
                direction = _turnRight[direction]

            connection &= ~connections[direction]

            if not connection & CONN_ALL4:
                connection &= ~CONN_MAYBE_NODE
                
            crackConnectionImage[pos] = connection
            continue
        elif connection & CONN_NODE:
            break

        if connection & CONN_MAYBE_NODE:
            # we simply pass over it, but we do not want to start a
            # new edge here during further down in the process:
            crackConnectionImage[pos] = connection & ~CONN_MAYBE_NODE

        direction = _turnRight[direction]
        while connection & connections[direction] == 0:
            direction = _turnLeft[direction]

    return result, pos, connections[(direction+2)%4]

def pyCrackEdgeGraph(labelImage, eightConnectedRegions = True,
                   progressHook = None):
    result = GeoMap(labelImage.size())

    cc = crackConnectionImage(labelImage)

    if eightConnectedRegions:
        for y in range(1, cc.height()-1):
          for x in range(1, cc.width()-1):
            if cc[x,y] == CONN_ALL4:
                if labelImage[x,y] == labelImage[x-1,y-1]:
                    cc[x,y] += CONN_DIAG_UPLEFT
                if labelImage[x-1,y] == labelImage[x,y-1]:
                    cc[x,y] += CONN_DIAG_UPRIGHT

                # crossing regions?
                if cc[x,y] == CONN_ALL4 + CONN_DIAG_UPLEFT + CONN_DIAG_UPRIGHT:
                    # preserve connectedness of higher label:
                    if labelImage[x,y-1] > labelImage[x-1,y-1]:
                        cc[x,y] -= CONN_DIAG_UPLEFT
                    else:
                        cc[x,y] -= CONN_DIAG_UPRIGHT

    for y in range(cc.height()):
        for x in range(cc.width()):
            conn = int(cc[x,y])
            if degree[conn] > 2:
                cc[x,y] = conn | CONN_NODE
            elif conn & CONN_ALL4 == (CONN_RIGHT | CONN_DOWN):
                cc[x,y] = conn | CONN_MAYBE_NODE
            if conn & CONN_DIAG:
                cc[x,y] = conn | CONN_MAYBE_NODE
    
    nodeImage = GrayImage(cc.size())
    # nodeImage encoding: each pixel's higher 28 bits encode the
    # (label + 1) of a node that has been inserted into the resulting
    # GeoMap at the corresponding position (+1 because zero is a valid
    # node label), while the lower 4 bits encode the four CONN_
    # directions in which a GeoMap edge is already connected to this
    # node

    progressHook = progressHook and progressHook.rangeTicker(cc.height())

    for startAt in (CONN_NODE, CONN_MAYBE_NODE):
     for y in range(cc.height()):
      if progressHook:
          progressHook()
      for x in range(cc.width()):
        nodeConn = int(cc[x, y])
        if nodeConn & startAt:
            startNodeInfo = int(nodeImage[x, y])
            if startNodeInfo:
                startNode = result.node((startNodeInfo >> 4) - 1)
            else:
                startNode = result.addNode((x - 0.5, y - 0.5))
                nodeImage[x, y] = startNodeInfo = \
                                  (startNode.label() + 1) << 4

            for direction, startConn in enumerate(connections):
                if nodeConn & startConn and not startNodeInfo & startConn:
                    edge, endPos, endConn = followEdge(
                        cc, (x, y), direction)
                    endNodeInfo = int(nodeImage[endPos])
                    if not endNodeInfo:
                        endNode = result.addNode((endPos[0] - 0.5, endPos[1] - 0.5))
                        endNodeInfo = (endNode.label() + 1) << 4
                    else:
                        assert not endNodeInfo & endConn, "double connection?"
                        endNode = result.node((endNodeInfo >> 4) - 1)

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
    return vigra.showImage(degreeImage)

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
            labels = vigra.readImage("headline.png")[0]
            cem = crackEdgeMap(labels)
            self.assertEqual(cem.nodeCount, 11)
            self.assertEqual(cem.edgeCount, 13)
            self.assertEqual(cem.faceCount, 11)

        def testComplexImage(self):
            labels = vigra.labelImageWithBackground8(
                vigra.readImage("crackConvert-test1.png")[0], 0)[0]
            cem = crackEdgeMap(labels)
            self.assertEqual(cem.faceCount, 26)

        def testCrackConnectionImage(self):
            labels = vigra.labelImageWithBackground8(
                vigra.readImage("crackConvert-test1.png")[0], 0)[0]
            cc = crackConnectionImage(labels)
            cc_ref = vigra.readImage("crackConvert-test1cc.png")
            self.assertEqual(
                vigra.inspectImage(cc - cc_ref, vigra.MinMax()).max(), 0)

    sys.exit(unittest.main())
