# -*- coding: iso-8859-1 -*-
_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

import vigra

def sivByOrder(order):
    return getattr(vigra, "SplineImageView%d" % order)

class GradientSIVProxy(object):
    def __init__(self, grad, SIV = 5):
        if type(SIV) == int:
            SIV = sivByOrder(SIV)
        self.gxsiv = SIV(grad[0])
        self.gysiv = SIV(grad[1])

    def __getitem__(self, i):
        return vigra.Vector2(self.gxsiv[i], self.gysiv[i])

class ThreeBandSIVProxy(object):
    def __init__(self, image, SIV = 5):
        if type(SIV) == int:
            SIV = sivByOrder(SIV)
        self.siv0 = SIV(image[0])
        self.siv1 = SIV(image[1])
        self.siv2 = SIV(image[2])

    def __getitem__(self, pos):
        return vigra.Vector(self.siv0[pos], self.siv1[pos], self.siv2[pos])

class HessianSIVProxy(object):
    def __init__(self, siv):
        self.siv = siv
    
    def __getitem__(self, pos):
        return vigra.Vector(self.siv.dxx(*pos),
                            self.siv.dxy(*pos),
                            self.siv.dyy(*pos))

TensorSIVProxy = ThreeBandSIVProxy
