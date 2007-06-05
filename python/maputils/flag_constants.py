# --------------------------------------------------------------------
#                             edge flags
# --------------------------------------------------------------------

# note that the highest bits are reserved for the C++ code:
#     static const unsigned int IS_BRIDGE = 0x8000000;

# edge protection:
BORDER_PROTECTION  = 1
SCISSOR_PROTECTION = 2
CONTOUR_PROTECTION = 4 # cf. maputils.protectFace()

# sinology:
COLUMN_PROTECTION = 8

ALL_PROTECTION     = 15 # or'ed values of the ones above

# delaunay:
CONTOUR_SEGMENT = 256
WEAK_CHORD      = 512
IS_BARB         = 1024

# tools/IntelligentScissors:
CURRENT_CONTOUR = 2048

# alphashapes:
ALPHA_MARK = 4096

# --------------------------------------------------------------------
#                             face flags
# --------------------------------------------------------------------

# note that the highest bits are reserved for the C++ code:
#     static const unsigned int BOUNDING_BOX_VALID = 0x8000000;
#     static const unsigned int AREA_VALID         = 0x4000000;

# delaunay:
OUTER_FACE = 1

# tools/ActivePaintbrush:
PROTECTED_FACE = 2 # flag which leads to CONTOUR_PROTECTION

# alphashapes:
ALPHA_MARK = 4

# maputils:
SRG_SEED = 8
SRG_BORDER = 16
