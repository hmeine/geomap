execfile("maptest.py")

if not '../cellimage' in sys.path:
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

def crackEdges2MidCracks(spmap, skipEverySecond = False):
    """Changes all edge geometry in-place, setting one point on the
    middle of each segment. Set skipEverySecond to True if each pixel
    crack is represented with two segments."""
    for edge in spmap.edgeIter():
        p = Polygon()
        step = skipEverySecond and 2 or 1
        p.append(edge[0])
        for i in range(0, len(edge)-1, step):
            p.append((edge[i]+edge[i+step])/2)
        p.append(edge[-1])
        edge.swap(p)

def cannyEdgeImageThinning(img):
    lut = [0]*256
    for k in [183, 222, 123, 237, 219, 111, 189, 246, 220, 115, 205,
              55, 103, 157, 118, 217]:
        lut[k] = 1
    res=GrayImage(img.size())
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

def cannyEdgeMap(i, scale, thresh):
    edgeImage = cannyEdgeImage(i, scale, thresh)
    edgeImage = cannyEdgeImageThinning(edgeImage)
    #lab,count = labelImageWithBackground4(edgeImage)
    geomap = GeoMap(edgeImage, 0, CellType.Line)
    spmap = pixelMap2subPixelMap(geomap, offset = Vector2(1,1),
                                 labelImageSize = i.size())
    return spmap

e = Experiment('kreuzung.png', "grad")
e("img")

# crackEdges can be:
# 0: 8-connected edges on 4-connected background
# 4: crack edges between 4-connected watershed regions
# 8: crack edges between 8-connected watershed regions
#    (still leads to disconnected regions in the map ATM)
crackEdges = 8

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
    crackEdges2MidCracks(spmap, True) # comment out to get original crack edges
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
