##########################################################################
#
#                Copyright 2007-2019 by Hans Meine
#
#     Permission is hereby granted, free of charge, to any person
#     obtaining a copy of this software and associated documentation
#     files (the "Software"), to deal in the Software without
#     restriction, including without limitation the rights to use,
#     copy, modify, merge, publish, distribute, sublicense, and/or
#     sell copies of the Software, and to permit persons to whom the
#     Software is furnished to do so, subject to the following
#     conditions:
#
#     The above copyright notice and this permission notice shall be
#     included in all copies or substantial portions of the
#     Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#     FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#     OTHER DEALINGS IN THE SOFTWARE.
#
##########################################################################

import math
from vigra import Vector2
from geomap import Polygon

__all__ = ("kochCurve", )

_kochCos = math.cos(math.pi/3)
_kochSin = math.sin(math.pi/3)

def _kochIteration(poly):
    result = Polygon()
    for i in range(len(poly)-1):
        segment = poly[i+1]-poly[i]
        smaller = segment/3
        left =  Vector2( smaller[0]*_kochCos - smaller[1]*_kochSin,
                         smaller[0]*_kochSin + smaller[1]*_kochCos)
        right = Vector2( smaller[0]*_kochCos + smaller[1]*_kochSin,
                        -smaller[0]*_kochSin + smaller[1]*_kochCos)
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
    p1 = p0 + Vector2(_kochCos,  _kochSin)
    result.append(p1)
    p2 = p1 + Vector2(_kochCos, -_kochSin)
    result.append(p2)
    result.append(p0)
    for i in range(level):
        result = _kochIteration(result)
    return result

# poly = Polygon((rotatePoly(kochCurve(5), math.pi/4)+Vector2(0.5, 0.5))*100)
