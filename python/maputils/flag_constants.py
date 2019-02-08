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

# --------------------------------------------------------------------
#                             edge flags
# --------------------------------------------------------------------

# note that the highest bits are reserved for the C++ code:
#     static const unsigned int IS_BRIDGE = 0x8000000;

# edge protection:
BORDER_PROTECTION  = 1
SCISSOR_PROTECTION = 2
CONTOUR_PROTECTION = 4 # cf. maputils.protectFace()
CUSTOM_PROTECTION  = 8

# sinology:
COLUMN_PROTECTION = 16

import geomap
ALL_PROTECTION = geomap.GeoMap.Edge.ALL_PROTECTION
assert geomap.GeoMap.Edge.BORDER_PROTECTION == BORDER_PROTECTION
assert ALL_PROTECTION & 31 == 31, "must include or'ed values of the ones above"

# delaunay:
CONTOUR_SEGMENT  = 256
WEAK_CHORD       = 512
IS_BARB          = 1024
START_NODE_ADDED = 8192
END_NODE_ADDED   = 16384

# tools/IntelligentScissors:
CURRENT_CONTOUR = 2048

# alphashapes:
ALPHA_MARK = 4096 # used for both edges and faces, see below

EDGE_USER = 0x100000

# --------------------------------------------------------------------
#                             face flags
# --------------------------------------------------------------------

# note that the highest bits are reserved for the C++ code:
#     enum {
#         BOUNDING_BOX_VALID = 0x80000000U,
#         AREA_VALID         = 0x40000000U,
#         INTERNAL_FLAGS     = 0xf0000000U,
#     };

# delaunay:
OUTER_FACE = 1

# tools/ActivePaintbrush:
PROTECTED_FACE = 2 # flag which leads to CONTOUR_PROTECTION

# alphashapes:
#ALPHA_MARK = 4096 # see above
CORNER_FACE = 32

# maputils:
SRG_SEED = 8
SRG_BORDER = 16

# levelcontours:
BACKGROUND_FACE = 1 # cf. OUTER_FACE
FOREGROUND_FACE = 32

FACE_USER = 0x100000

FACE_NONINTERNAL = 0x0fffffff

if __name__ == "__main__":
    for name in dir():
        value = eval(name)
        if isinstance(value, int):
            print "%20s: %8x (%d)" % (name, value, value)
