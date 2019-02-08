# -*- coding: iso-8859-1 -*-
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

import vigra

def sivByOrder(order):
    return getattr(vigra.sampling, "SplineImageView%d" % order)

class GradientSIVProxy(object):
    def __init__(self, grad, SIV = 5):
        if type(SIV) == int:
            SIV = sivByOrder(SIV)
        self.gxsiv = SIV(grad[0])
        self.gysiv = SIV(grad[1])

    def __getitem__(self, i):
        return geomap.Vector2(self.gxsiv[i], self.gysiv[i])

class ThreeBandSIVProxy(object):
    def __init__(self, image, SIV = 5):
        if type(SIV) == int:
            SIV = sivByOrder(SIV)
        self.siv0 = SIV(image[0])
        self.siv1 = SIV(image[1])
        self.siv2 = SIV(image[2])

    def __getitem__(self, pos):
        return geomap.Vector3(self.siv0[pos], self.siv1[pos], self.siv2[pos])

class HessianSIVProxy(object):
    def __init__(self, siv):
        self.siv = siv

    def __getitem__(self, pos):
        return geomap.Vector3(self.siv.dxx(pos[0], pos[1]),
                            self.siv.dxy(pos[0], pos[1]),
                            self.siv.dyy(pos[0], pos[1]))

TensorSIVProxy = ThreeBandSIVProxy
