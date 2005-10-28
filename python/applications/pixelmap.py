execfile("maptest.py")

def pixelMapData(geomap, offset, scale):
    nodes = [None] * geomap.maxNodeLabel()
    for node in geomap.nodes:
        ul = node.bounds.upperLeft()
        center = Vector2(*ul) + Vector2(node.bounds.width() - 1,
                                        node.bounds.height() - 1) / 2
        nodes[node.label] = (center+offset) * scale

    edges = [None] * geomap.maxEdgeLabel()
    for edge in geomap.edges:
        points = [Vector2(*p)+offset for p in iter(edge.start)]
        startNodeLabel = edge.start.startNodeLabel()
        endNodeLabel = edge.end.startNodeLabel()
        points.insert(0, nodes[startNodeLabel])
        points.append(nodes[endNodeLabel])
        edges[edge.label] = (
            startNodeLabel, endNodeLabel, Polygon(points) * scale)

    return nodes, edges

def pixelMap2subPixelMap(geomap, scale = 1.0, offset = Vector2(2, 2)):
    nodes, edges = pixelMapData(geomap, offset, scale)
    return Map(nodes, edges, geomap.cellImage.size() * scale)

from cellimage import *
e = Experiment("../../../Testimages/blox.gif", "grad")
e("img")

# lab, count = watershedUnionFind8(e.img.bi.gm)
# ce = regionImageToCrackEdgeImage(lab, 0)
# geomap = GeoMap(ce, 0, CellType.Line)

lab, count = watershedSegmentation(e.img.bi.gm, KeepContours)
geomap = GeoMap(lab, 0, CellType.Vertex)

# face = geomap.faces[14]
# for c in face.contours:
#     print list(iter(c))

spmap = pixelMap2subPixelMap(geomap)
d = MapDisplay(e.img, spmap)
