import vigra, geomap, crackConvert
import cellimage
from geomap import Vector2, Size2D, Point2D, Rect2D
from flag_constants import BORDER_PROTECTION

__all__ = ["pixelMap2subPixelMap", "crackEdgeMap",
           "crackEdges2MidCracks", "cannyEdgeMap", "pixelWatershedMap"]

def pixelMap2subPixelMap(pixelMap, scale = 1.0, offset = Vector2(0, 0),
                         imageSize = None,
                         skipEverySecond = False, midCracks = False):

    """Extract node positions and edge geometry from a
    cellimage.GeoMap object and returns a new subpixel-GeoMap object
    initialized with it.  For nodes, this function simply uses the
    center of their bounding box.  All positions are shifted by the
    optional offset and then scaled with the given factor.  The
    imageSize defaults to the (scaled) pixel-based pixelMap's
    cellImage.size().

    Set skipEverySecond to True if the pixelMap contains a crack edge
    map (otherwise, each resulting edge segment will have an
    additional mid crack point).  If you set midCracks to True, the
    edge geometry will include the midpoints of each crack instead of
    the endpoints (skipEverySecond is implied by/ignored if midCracks
    == True)."""

    if imageSize == None:
        imageSize = pixelMap.cellImage.size() * scale
    result = geomap.GeoMap(imageSize = imageSize)

    if midCracks:
        skipEverySecond = False

    nodes = [None] * (pixelMap.maxNodeLabel() + 1)
    for node in pixelMap.nodes:
        ul = node.bounds.upperLeft()
        center = Vector2(ul[0], ul[1]) + \
                 Vector2(node.bounds.width() - 1,
                         node.bounds.height() - 1) / 2
        nodes[node.label] = result.addNode((center+offset) * scale)

    # mark as sorted (sigma order will be copied from source map):
    result.sortEdgesDirectly()

    undesirable = []

    edges = [None] * (pixelMap.maxEdgeLabel() + 1)
    for edge in pixelMap.edges:
        it = iter(edge.start)
        startPos = it.nodePosition()
        points = list(it)
        if midCracks:
            points = points[::2]
        endPos = it.nodePosition()

        points.insert(0, startPos)
        points.append(endPos)

        if skipEverySecond:
            points = points[::2]

        points = [(Vector2(p[0], p[1])+offset) * scale for p in points]

        startNeighbor = nodes[edge.start.startNodeLabel()]
        endNeighbor = nodes[edge.end.startNodeLabel()]

        if startNeighbor.position() != points[0]:
            points.insert(0, startNeighbor.position())
        if endNeighbor.position() != points[-1]:
            points.append(endNeighbor.position())

        # re-use sigma order from source map:
        if not startNeighbor.isIsolated():
            neighbor = cellimage.GeoMap.DartTraverser(edge.start)
            while neighbor.nextSigma().edgeLabel() >= edge.label:
                pass
            startNeighbor = edges[neighbor.edgeLabel()].dart()
            if neighbor == neighbor.edge().end:
                startNeighbor.nextAlpha()

        if not endNeighbor.isIsolated():
            neighbor = cellimage.GeoMap.DartTraverser(edge.end)
            while neighbor.nextSigma().edgeLabel() >= edge.label:
                pass
            endNeighbor = edges[neighbor.edgeLabel()].dart()
            if neighbor == neighbor.edge().end:
                endNeighbor.nextAlpha()
                
        newEdge = result.addEdge(startNeighbor, endNeighbor, points)
        edges[edge.label] = newEdge

        if newEdge.isLoop() and newEdge.partialArea() == 0.0:
            undesirable.append(newEdge.dart())

    for dart in undesirable:
        result.removeEdge(dart)

    result.initializeMap()

    # the border closing was done in C++, so we have to mark the
    # border edges manually:
    assert result.face(0).holeCount() == 1
    for dart in result.face(0).holeContours().next().phiOrbit():
        assert not dart.leftFaceLabel()
        dart.edge().setFlag(BORDER_PROTECTION)

    return result

def crackEdges2MidCracks(subpixelMap):
    """crackEdges2MidCracks(subpixelMap)

    Changes all edge geometry in-place, setting one point on the
    middle of each edge segment (and removes each segment's end
    points, except for the edge ends).

    Note that this is done in a brutal way; the resulting map will
    fail checkConsistency()."""
    
    for edge in subpixelMap.edgeIter():
        p = geomap.Polygon()
        p.append(edge[0])
        for i in range(0, len(edge)-1):
            p.append((edge[i]+edge[i+1])/2)
        p.append(edge[-1])
        edge.swap(p) # FIXME: use edge.setGeometry as soon as that's finished

def cannyEdgeMap(image, scale, thresh):
    """cannyEdgeMap(image, scale, thresh)

    Returns a subpixel-GeoMap object containing thinned canny edges
    obtained from cannyEdgeImage(image, scale, thresh).
    (Internally creates a pixel GeoMap first.)"""
    
    edgeImage = vigra.analysis.cannyEdgeImageWithThinning(image, scale, thresh, 1)
    pixelMap = cellimage.GeoMap(edgeImage, 0, cellimage.CellType.Line)
    spmap = pixelMap2subPixelMap(
        pixelMap, offset = Vector2(1,1), imageSize = image.size())
    return spmap

def crackEdgeMap(labelImage, midCracks = False):
    """crackEdgeMap(labelImage, midCracks = False)

    Returns a subpixel-GeoMap containing crack-edge contours extracted
    from the given labelImage.  If the optional parameter 'midCracks'
    is True(default), the resulting edges consist of the connected
    midpoints of the cracks, not of the crack segments themselves."""

    print "- creating pixel-based GeoMap..."
    ce = vigra.analysis.regionImageToCrackEdgeImage(labelImage + 1, 0)
    pixelMap = cellimage.GeoMap(ce, 0, cellimage.CellType.Line)

    print "- converting pixel-based GeoMap..."
    result = pixelMap2subPixelMap(
        pixelMap, 0.5, imageSize = (pixelMap.cellImage.size()-Size2D(3,3))/2,
        skipEverySecond = True, midCracks = midCracks)
    return result

def thresholdMap(scalarImage, threshold, midCracks = False):
    """Shortcut for calling crackEdgeMap with a thresholded image."""
    return crackEdgeMap(scalarImage >= threshold, midCracks)

def pixelWatershedMap(biImage, crackEdges = 4, midCracks = False):
    """pixelWatershedMap(biImage, crackEdges = 4, midCracks = False)

    Performs a watershed segmentation on biImage and returns a
    subpixel-GeoMap containing the resulting contours.  The type of
    watershed segmentation depends on the 'crackEdges' parameter:

    0: 8-connected edges on 4-connected background       (SRG alg.)
    4: crack edges between 4-connected watershed regions (UF alg.)
    8: crack edges between 8-connected watershed regions (UF alg.)

    If midCracks is True, the resulting edges consist of the
    connected midpoints of the cracks, not of the crack segments
    themselves (this parameter is ignored if crackEdges == 0)."""
    
    if crackEdges:
        print "- Union-Find watershed segmentation..."
        lab, count = getattr(vigra, 'watershedUnionFind' + str(crackEdges))(biImage)
        if not midCracks:
            return crackConvert.crackEdgeMap(
                lab, eightConnectedRegions = (crackEdges == 8))
        return crackEdgeMap(lab, midCracks)

    print "- watershed segmentation..."
    lab, count = vigra.analysis.watershedSegmentation(biImage, vigra.KeepContours)

    print "- creating pixel-based GeoMap..."
    pixelMap = cellimage.GeoMap(lab, 0, cellimage.CellType.Vertex)

    print "- converting pixel-based GeoMap..."
    return pixelMap2subPixelMap(
        pixelMap, imageSize = (pixelMap.cellImage.size()-Size2D(4,4)))

def cellImage2display(cellImage, background = None,
                      nodeColor = (0, 0, 255),
                      edgeColor = (255, 0, 0)):
    assert False, "cellImage2display needs to be ported to vigranumpy"
    if hasattr(cellImage, "cellImage"):
        cellImage = cellImage.cellImage
    if hasattr(cellImage, "serialize"):
        cellImage = cellImage.serialize()
    if background is None:
        background = cellImage / 4
    if background.bands() < 3:
        background = vigra.RGBImage(background)
    if background.size() != cellImage.size():
        nbg = vigra.RGBImage(cellImage.size())
        shift = Point2D((cellImage.width() - background.width())/2,
                        (cellImage.height() - background.height())/2)
        nbg.subImage(Rect2D(shift, background.size())).copyValues(background)
        background = nbg
    return transformImage(
        cellImage, background,
        "\l ci, b: ci %% 4 == 1 ? %r : (ci %% 4 == 2 ? %r : b)" % (
        edgeColor, nodeColor))
