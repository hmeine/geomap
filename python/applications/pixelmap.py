from vigra import *
addPathFromHere('../cellimage')
import cellimage, hourglass
from flag_constants import BORDER_PROTECTION
#from map import GeoMap
from hourglass import GeoMap

__all__ = ["pixelMap2subPixelMap", "crackEdgeMap",
           "crackEdges2MidCracks", "cannyEdgeMap", "pixelWatershedMap"]

def pixelMapData(geomap, scale = 1.0, offset = Vector2(0, 0),
                 skipEverySecond = False):
    """pixelMapData(geomap, scale = 1.0, offset = Vector2(0, 0),
                 skipEverySecond = False)

    Extracts node positions and edge geometry from a GeoMap object.
    For nodes, this function simply calculates their center of mass.
    All positions are shifted by the optional offset and then scaled
    with the given factor. Set skipEverySecond to True if the geomap
    contains a crack edge map (otherwise, each resulting edge segment
    will have an additional mid crack point)."""

    nodes = [None] * (geomap.maxNodeLabel() + 1)
    for node in geomap.nodes:
        ul = node.bounds.upperLeft()
        center = Vector2(*ul) + Vector2(node.bounds.width() - 1,
                                        node.bounds.height() - 1) / 2
        nodes[node.label] = (center+offset) * scale

    edges = [None] * (geomap.maxEdgeLabel() + 1)
    for edge in geomap.edges:
        points = [(Vector2(*p)+offset) * scale for p in iter(edge.start)]
        startNodeLabel = edge.start.startNodeLabel()
        endNodeLabel = edge.end.startNodeLabel()
        points.insert(0, nodes[startNodeLabel])
        points.append(nodes[endNodeLabel])
        if skipEverySecond:
            points = [points[i] for i in range(0, len(points), 2)]
        edges[edge.label] = (
            startNodeLabel, endNodeLabel, hourglass.Polygon(points))
    return nodes, edges

def pixelMap2subPixelMap(geomap, scale = 1.0, offset = Vector2(0, 0),
                         labelImageSize = None, skipEverySecond = False):
    """pixelMap2subPixelMap(geomap, scale = 1.0, offset = Vector2(0, 0), labelImageSize = None)

    Uses pixelMapData() to extract the pixel-geomap's geometry and
    returns a new subpixel-GeoMap object initialized with it.  The
    labelImageSize defaults to the (scaled) pixel-based geomap's
    cellImage.size().  See also the documentation of pixelMapData()."""
    
    if labelImageSize == None:
        labelImageSize = geomap.cellImage.size() * scale
    result = GeoMap(imageSize = labelImageSize)

    nodes = [None] * (geomap.maxNodeLabel() + 1)
    for node in geomap.nodes:
        ul = node.bounds.upperLeft()
        center = Vector2(*ul) + Vector2(node.bounds.width() - 1,
                                        node.bounds.height() - 1) / 2
        nodes[node.label] = result.addNode((center+offset) * scale)

    edges = [None] * (geomap.maxEdgeLabel() + 1)
    for edge in geomap.edges:
        points = [(Vector2(p[0], p[1])+offset) * scale for p in iter(edge.start)]
        startNodeLabel = edge.start.startNodeLabel()
        endNodeLabel = edge.end.startNodeLabel()
        points.insert(0, nodes[startNodeLabel].position())
        points.append(nodes[endNodeLabel].position())
        if skipEverySecond:
            points = [points[i] for i in range(0, len(points), 2)]
        result.addEdge(
            nodes[startNodeLabel],
            nodes[endNodeLabel],
            points)

    result.initializeMap()

    # the border closing was done in C++, so we have to mark the
    # border edges manually:
    for edge in result.edgeIter():
        if not edge.leftFaceLabel() or not edge.rightFaceLabel():
            edge.setFlag(BORDER_PROTECTION)

    return result

def crackEdges2MidCracks(subpixelMap):
    """crackEdges2MidCracks(subpixelMap)

    Changes all edge geometry in-place, setting one point on the
    middle of each edge segment (and removes each segment's end
    points, except for the edge ends)."""
    
    for edge in subpixelMap.edgeIter():
        p = hourglass.Polygon()
        p.append(edge[0])
        for i in range(0, len(edge)-1):
            p.append((edge[i]+edge[i+1])/2)
        p.append(edge[-1])
        edge.swap(p) # FIXME: use edge.setGeometry as soon as that's finished

def cannyEdgeImageThinning(img):
    lut = [0]*256
    for k in [183, 222, 123, 237, 219, 111, 189, 246, 220, 115, 205,
              55, 103, 157, 118, 217]:
        lut[k] = 1
    res = GrayImage(img.size())
    res[1:-1,1:-1].copyValues(img[1:-1,1:-1])
    for y in range(1, res.height()-1):
        for x in range(1, res.width()-1):
            if res[x,y]:
                res[x,y] = 1
                continue # no edge pixel
            n = [k for k in res.neighborhood8((x,y))]
            conf = 0
            pt, nt = [], []
            for k in range(8):
                conf |= int(n[k] == 0) << k
                if n[k] == 0 and n[k-1] != 0:
                    nt.append(k)
                if n[k] != 0 and n[k-1] == 0:
                    pt.append(k)
                    lab = n[k]
            if len(pt) != 1:
                if lut[conf]:
                    res[x,y] = 1
                continue
            if nt[0] > pt[0]:
                pt[0] += 8
            if pt[0]-nt[0] >= 3:
                res[x,y] = 1
    return res[1:-1,1:-1].clone()

def cannyEdgeMap(image, scale, thresh):
    """cannyEdgeMap(image, scale, thresh)

    Returns a subpixel-GeoMap object containing thinned canny edges
    obtained from cannyEdgeImage(image, scale, thresh).
    (Internally creates a pixel GeoMap first.)"""
    
    edgeImage = cannyEdgeImage(image, scale, thresh)
    edgeImage = cannyEdgeImageThinning(edgeImage)
    geomap = cellimage.GeoMap(edgeImage, 0, cellimage.CellType.Line)
    spmap = pixelMap2subPixelMap(
        geomap, offset = Vector2(1,1), labelImageSize = image.size())
    return spmap

def crackEdgeMap(labelImage, midCracks = False):
    """crackEdgeMap(labelImage, midCracks = False)

    Returns a subpixel-GeoMap containing crack-edge contours extracted
    from the given labelImage.  If the optional parameter 'midCracks'
    is True(default), the resulting edges consist of the connected
    midpoints of the cracks, not of the crack segments themselves."""

    print "- creating pixel-based GeoMap..."
    ce = regionImageToCrackEdgeImage(transformImage(labelImage, "\l x:x+1"), 0)
    geomap = cellimage.GeoMap(ce, 0, cellimage.CellType.Line)

    print "- converting pixel-based GeoMap..."
    result = pixelMap2subPixelMap(
        geomap, 0.5, labelImageSize = (geomap.cellImage.size()-Size2D(3,3))/2,
        skipEverySecond = True) # source image had explicit cracks
    if midCracks:
        print "  (converting cracks to mid-cracks...)"
        crackEdges2MidCracks(result)
    return result

def pixelWatershedMap(biImage, crackEdges = 4, midCracks = False):
    """pixelWatershedMap(biImage, crackEdges = 4, midCracks = False)

    Performs a watershed segmentation on biImage and returns a
    subpixel-GeoMap containing the resulting contours.  The type of
    watershed segmentation depends on the 'crackEdges' parameter:

    0: 8-connected edges on 4-connected background
    4: crack edges between 4-connected watershed regions
    8: crack edges between 8-connected watershed regions
       (8-connected regions will be separated in the result ATM)

    If midCracks is True, the resulting edges consist of the
    connected midpoints of the cracks, not of the crack segments
    themselves (this parameter is ignored if crackEdges == 0)."""
    
    if crackEdges:
        print "- Union-Find watershed segmentation..."
        lab, count = eval('watershedUnionFind' + str(crackEdges))(biImage)
        return crackEdgeMap(lab, midCracks)

    print "- watershed segmentation..."
    lab, count = watershedSegmentation(biImage, KeepContours)

    print "- creating pixel-based GeoMap..."
    geomap = cellimage.GeoMap(lab, 0, cellimage.CellType.Vertex)

    print "- converting pixel-based GeoMap..."
    return pixelMap2subPixelMap(
        geomap, labelImageSize = (geomap.cellImage.size()-Size2D(4,4)))
