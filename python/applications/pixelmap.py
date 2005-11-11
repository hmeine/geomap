execfile("maptest.py")

sys.path.append('../cellimage')
from cellimage import *

def pixelMapData(geomap, offset, scale):
    nodes = [None] * geomap.maxNodeLabel()
    for node in geomap.nodes:
        ul = node.bounds.upperLeft()
        center = Vector2(*ul) + Vector2(node.bounds.width() - 1,
                                        node.bounds.height() - 1) / 2
        nodes[node.label] = (center+offset) * scale
    edges = [None] * geomap.maxEdgeLabel()
    for edge in geomap.edges:
        points = [(Vector2(*p)+offset) * scale for p in iter(edge.start)]
        startNodeLabel = edge.start.startNodeLabel()
        endNodeLabel = edge.end.startNodeLabel()
        points.insert(0, nodes[startNodeLabel])
        points.append(nodes[endNodeLabel])
        edges[edge.label] = (
            startNodeLabel, endNodeLabel, Polygon(points))
    return nodes, edges

def pixelMap2subPixelMap(geomap, scale = 1.0, offset = Vector2(0, 0),
                         labelImageSize = None):
    nodes, edges = pixelMapData(geomap, offset, scale)
    if labelImageSize == None:
        labelImageSize = geomap.cellImage.size() * scale
    return Map(nodes, edges, labelImageSize,
               performBorderClosing = False, performEdgeSplits = False)

e = Experiment(filename, "grad")
e("img")

crackEdges = 8  # can be 0, 4, or 8

print "- watershed segmentation..."
if crackEdges:
    lab, count = eval('watershedUnionFind'+str(crackEdges))(e.img.bi.gm)
    ce = regionImageToCrackEdgeImage(lab, 0)
else:
    lab, count = watershedSegmentation(e.img.bi.gm, KeepContours)

print "- creating pixel-based GeoMap..."
if crackEdges:
    geomap = GeoMap(ce, 0, CellType.Line)
else:
    geomap = GeoMap(lab, 0, CellType.Vertex)

print "- converting pixel-based GeoMap..."
if crackEdges:
    spmap = pixelMap2subPixelMap(
        geomap, 0.5, labelImageSize = (geomap.cellImage.size()-Size2D(3,3))/2)
else:
    spmap = pixelMap2subPixelMap(geomap, 
                labelImageSize = (geomap.cellImage.size()-Size2D(4,4)))
spmap.faceStats = FaceIntensityStatistics(spmap, e.img.view)
mergeZeroPixelFaces(spmap)

print "- creating display..."
d = MapDisplay(e.img, spmap)
d.setCaption(str(e.img.name+' crackEdges '+str(crackEdges)).replace("_", " "))

if spmap.performEdgeSplits:
    print "*** split results: ***"
    for edge in spmap.edgeIter():
        if hasattr(edge, "isSplitResultOf"):
            print edge
