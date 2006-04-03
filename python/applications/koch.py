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
    result.append(Vector2(0,0))
    result.append(Vector2(kochCos, kochSin))
    result.append(Vector2(kochCos, kochSin) + Vector2(kochCos, -kochSin))
    result.append(Vector2(0,0))
    for i in range(level):
        result = kochIteration(result)
    return result

