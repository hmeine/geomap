import copy, sys, time
from hourglass import GeoMap
from vigra import meshIter
import maputils

# --------------------------------------------------------------------
#                              FRONTEND
# --------------------------------------------------------------------

def crackEdgeMap(labelImage, initLabelImage = True):
    sys.stdout.write("- looking for local connections..."); c = time.clock()
    cc = crackConnectionImage(labelImage)
    sys.stdout.write("done. (%ss)\n" % (time.clock()-c, ))

    sys.stdout.write("- following crack edges..."); c = time.clock()
    nodes, edges = crackConnections(cc)
    sys.stdout.write("done. (%ss)\n" % (time.clock()-c, ))

    sys.stdout.write("- creating GeoMap...\n"); c = time.clock()
    result = GeoMap(nodes, edges, labelImage.size())

    maputils.mergeDegree2Nodes(result)

    sys.stdout.write("  sorting edges..."); c1 = time.clock()
    result.sortEdgesDirectly()
    sys.stdout.write(" (%ss)\n" % (time.clock()-c1, ))

    sys.stdout.write("  initializing faces..."); c1 = time.clock()
    result.initializeMap(initLabelImage)
    sys.stdout.write(" (%ss)\n" % (time.clock()-c1, ))
    
    sys.stdout.write("  done. (%ss)\n" % (time.clock()-c, ))

    return result

# --------------------------------------------------------------------
# 						   helper functions
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

def crackConnectionImage(labelImage):
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

# flag whether multiple steps into the same direction should be pruned:
simplifyStraight = False

# for extracting loops, we add nodes at every upper left corner:
def isNode(connValue):
    return degree[connValue] > 2 or connValue == (CONN_RIGHT | CONN_DOWN)

def followEdge(crackConnectionImage, pos, direction):
    pos = Point2D(*pos)
# 	dirName = ["right", "up", "left", "down"]
# 	sys.stderr.write("following crack edge from %s in '%s'-direction.." % (
# 		pos, dirName[direction]))
    vPos = Vector2(*pos) - Vector2(0.5, 0.5)
    result = [copy.copy(vPos)]
    prevDirection = None
    while True:
# 		sys.stderr.write(" [%s]" % direction)
        vPos += _dirVector[direction]
        if not simplifyStraight or prevDirection != direction:
            result.append(copy.copy(vPos))
        else:
            result[-1] = copy.copy(vPos)
        #result.append(Vector2(*pos) - Vector2(0.5, 0.5) + 0.5*dirVector[direction])
        pos += _dirOffset[direction]

        connection = int(crackConnectionImage[pos])
        if isNode(connection):
            break

        prevDirection = direction
        direction = _turnRight[direction]
        while connection & connections[direction] == 0:
# 			sys.stderr.write(" (%s)" % direction)
            direction = _turnLeft[direction]
# 	result.append(Vector2(*pos) - Vector2(0.5, 0.5))
# 	sys.stderr.write(" %d steps.\n" % (len(result)-2, ))
    return result, pos, connections[(direction+2)%4]

def crackConnections(crackConnectionImage):
    edges = [None]
    nodeImage = GrayImage(crackConnectionImage.size())
    nodePositions = [None]
    for startPos in meshIter(crackConnectionImage.size()):
        #startPos = Point2D(*startPos)
        nodeConn = int(crackConnectionImage[startPos])
        if isNode(nodeConn):
            startNode = int(nodeImage[startPos])
            if not startNode:
                nodeImage[startPos] = startNode = len(nodePositions) << 4
                nodePositions.append(Vector2(*startPos) - Vector2(0.5, 0.5))
            for direction, startConn in enumerate(connections):
                if nodeConn & startConn and not startNode & startConn:
                    edge, endPos, endConn = followEdge(
                        crackConnectionImage, startPos, direction)
                    endNode = int(nodeImage[endPos])
                    if not endNode:
                        endNode = len(nodePositions) << 4
                        nodePositions.append(
                            Vector2(*endPos) - Vector2(0.5, 0.5))
                    assert not endNode & endConn, "double connection?"
                    edges.append((startNode >> 4, endNode >> 4, Polygon(edge)))
                    startNode |= startConn
                    if endNode >> 4 == startNode >> 4:
                        startNode |= endConn
                        nodeImage[startPos] = startNode
                    else:
                        nodeImage[startPos] = startNode
                        nodeImage[endPos] = endNode | endConn
    return nodePositions, edges

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

        def testHeadlineMap(self):
            from vigra import readImage
            labels = readImage("headline.png")[0]
            cem = crackEdgeMap(labels)
            self.assertEqual(cem.nodeCount, 11)
            self.assertEqual(cem.edgeCount, 13)
            self.assertEqual(cem.faceCount, 11)

    sys.exit(unittest.main())
