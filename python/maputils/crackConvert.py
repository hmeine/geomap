import copy, sys, time
from hourglass import GeoMap
from map import connectBorderNodes

# --------------------------------------------------------------------
#                              FRONTEND
# --------------------------------------------------------------------

def crackEdgeMap(labelImage):
    sys.stdout.write("- looking for local connections..."); c = time.clock()
    cc = crackConnectionImage(labelImage)
    sys.stdout.write("done. (%ss)\n" % (time.clock()-c, ))

    sys.stdout.write("- following crack edges..."); c = time.clock()
    nodes, edges = crackConnections(cc)
    sys.stdout.write("done. (%ss)\n" % (time.clock()-c, ))

    sys.stdout.write("- creating GeoMap...\n")
    result = GeoMap(nodes, edges, labelImage.size())
    connectBorderNodes(result, 0.1, aroundPixels = True)
    result.sortEdgesDirectly()
    result.initializeMap()
    return result

# --------------------------------------------------------------------
# 						   helper functions
# --------------------------------------------------------------------

from vigra import GrayImage, Vector2, Point2D, Size2D
from hourglass import Polygon

CONN_UP = 1
CONN_DOWN = 2
CONN_LEFT = 4
CONN_RIGHT = 8

connections = [CONN_RIGHT, CONN_UP, CONN_LEFT, CONN_DOWN]

degree = [0, 1, 1, 2, # 0000 0001 0010 0011
          1, 2, 2, 3, # 0100 0101 0110 0111
          1, 2, 2, 3, # 1000 1001 1010 1011
          2, 3, 3, 4, # 1100 1101 1110 1111
          ]

def crackConnectionImage(labelImage):
    result = GrayImage(labelImage.size()+Size2D(1,1))

    for y in range(labelImage.height()-1):
        for x in range(labelImage.width()-1):
            if labelImage[(x, y)] != labelImage[(x+1, y)]:
                result[(x+1, y  )] += CONN_DOWN
                result[(x+1, y+1)] += CONN_UP
            if labelImage[(x, y)] != labelImage[(x, y+1)]:
                result[(x,   y+1)] += CONN_RIGHT
                result[(x+1, y+1)] += CONN_LEFT
        x = labelImage.width()-1
        if labelImage[(x, y)] != labelImage[(x, y+1)]:
            result[(x,   y+1)] += CONN_RIGHT
            result[(x+1, y+1)] += CONN_LEFT

    lastRow = labelImage.subImage((0, labelImage.height()-1),
                                  Size2D(labelImage.width(), 1))
    for x in range(labelImage.width()-1):
        if lastRow[x, 0] != lastRow[x+1, 0]:
            result[x+1, labelImage.height()-1] += CONN_DOWN
            result[x+1, labelImage.height()] += CONN_UP

    return result

_dirOffset = [Size2D(1, 0), Size2D(0, -1), Size2D(-1, 0), Size2D(0, 1)]
_dirVector = [Vector2(1, 0), Vector2(0, -1), Vector2(-1, 0), Vector2(0, 1)]
_turnRight = [3, 0, 1, 2]
_turnLeft  = [1, 2, 3, 0]

# flag whether multiple steps into the same direction should be pruned:
simplifyStraight = False

def followEdge(crackConnectionImage, pos, direction):
    pos = Point2D(*pos)
# 	dirName = ["right", "up", "left", "down"]
# 	sys.stderr.write("following crack edge from %s in '%s'-direction.." % (
# 		pos, dirName[direction]))
    vPos = Vector2(*pos) - Vector2(0.5, 0.5)
    result = [copy.copy(vPos)]
    prevDirection = None
    while True:
        vPos += _dirVector[direction]
        if not simplifyStraight or prevDirection != direction:
            result.append(copy.copy(vPos))
        else:
            result[-1] = copy.copy(vPos)
        #result.append(Vector2(*pos) - Vector2(0.5, 0.5) + 0.5*dirVector[direction])
        pos += _dirOffset[direction]
        if not pos in crackConnectionImage.size():
            break

        connection = int(crackConnectionImage[pos])
        if degree[connection] != 2:
            break

        prevDirection = direction
        direction = _turnRight[direction]
        while connection & connections[direction] == 0:
            direction = _turnLeft[direction]
# 	result.append(Vector2(*pos) - Vector2(0.5, 0.5))
# 	sys.stderr.write(" %d steps.\n" % (len(result)-2, ))
    return result, pos, connections[(direction+2)%4]

def crackConnections(crackConnectionImage):
    edges = [None]
    nodeImage = GrayImage(crackConnectionImage.size())
    nodePositions = [None]
    for startPos in crackConnectionImage.size():
        #startPos = Point2D(*startPos)
        nodeConn = int(crackConnectionImage[startPos])
        if degree[nodeConn] > 2:
            startNode = int(nodeImage[startPos])
            if not startNode:
                nodeImage[startPos] = startNode = len(nodePositions) << 4
                nodePositions.append(Vector2(*startPos) - Vector2(0.5, 0.5))
            for direction, connection in enumerate(connections):
                if nodeConn & connection and not startNode & connection:
                    edge, endPos, endConn = followEdge(
                        crackConnectionImage, startPos, direction)
                    endNode = int(nodeImage[endPos])
                    if not endNode:
                        endNode = len(nodePositions) << 4
                        nodePositions.append(Vector2(*endPos) - Vector2(0.5, 0.5))
                    assert not endNode & endConn, "double connection?"
                    edges.append((startNode >> 4, endNode >> 4, Polygon(edge)))
                    startNode |= connection
                    nodeImage[endPos] = endNode | endConn
            nodeImage[startPos] = startNode
    return nodePositions, edges

def showDegrees(crackConnectionImage):
    degreeImage = GrayImage(crackConnectionImage.size())
    for p in crackConnectionImage.size():
        degreeImage[p] = degree[int(crackConnectionImage[p])]
    return showImage(degreeImage)
