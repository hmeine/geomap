kochCos = math.cos(math.pi/3)
kochSin = math.sin(math.pi/3)

def kochIteration(poly):
    result = Polygon()
    for i in range(len(poly)-1):
        segment = poly[i+1]-poly[i]
        smaller = segment/3
        left =  Vector2( smaller[0]*kochCos - smaller[1]*kochSin,
                         smaller[0]*kochSin + smaller[1]*kochCos)
        right = Vector2( smaller[0]*kochCos + smaller[1]*kochSin,
                        -smaller[0]*kochSin + smaller[1]*kochCos)
        p1 = poly[i] + smaller
        p2 = p1 + left
        p3 = p2 + right
        result.append(poly[i])
        result.append(p1)
        result.append(p2)
        result.append(p3)
    result.append(result[0])
    return result

def kochCurve(level = 5):
    result = Polygon()
    p0 = Vector2(-0.5, -math.sqrt(1./12))
    result.append(p0)
    p1 = p0 + Vector2(kochCos, kochSin)
    result.append(p1)
    p2 = p1 + Vector2(kochCos, -kochSin)
    result.append(p2)
    result.append(p0)
    for i in range(level):
        result = kochIteration(result)
    return result

def rotatedPoly(poly, angle):
    """rotatedPoly(poly, angle)

    Rotate polygon by the given angle around the origin."""
    
    unitX = Vector2(math.cos(angle), -math.sin(angle))
    unitY = Vector2(math.sin(angle),  math.cos(angle))
    result = Polygon()
    for point in poly:
        result.append(Vector2(dot(point, unitX), dot(point, unitY)))
    return result

def samplePoly(poly, shift = None, size = None):
    """samplePoly(poly, size)
    Sample poly with a regular grid at integer coordinates starting
    from (0,0) to the given size (which has to be a Size2D object)."""

    if size == None:
        size = poly.boundingBox().size()
        size = Size2D(int(math.ceil(size[0]))+2,
                      int(math.ceil(size[1]))+2)
        if not shift:
            shift = Vector2(0, 0)
        shift = shift + Vector2(1, 1) - poly.boundingBox().begin()
        poly = Polygon(poly + shift)

    result = GrayImage(size)
    for p in size:
        result[p] = poly.contains(Vector2(p[0], p[1])) and 1 or 0
    return result

from cellimage import GeoMap, CellType

def polyCrackMap(poly, shift = None, midCracks = True):
    img = samplePoly(poly, shift)
    ce = regionImageToCrackEdgeImage(transformImage(img, "\l x:x+1"), 0)
    geomap = GeoMap(ce, 0, CellType.Line)
    spmap = pixelMap2subPixelMap(
        geomap, 0.5, labelImageSize = (geomap.cellImage.size()-Size2D(3,3))/2)
    if midCracks:
        crackEdges2MidCracks(spmap, True)
    return spmap

# poly = Polygon((rotatedPoly(kochCurve(5), math.pi/4)+Vector2(0.5, 0.5))*100)
