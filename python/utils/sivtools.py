# -*- coding: iso-8859-1 -*-
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
