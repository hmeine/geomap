# --------------------------------------------------------------------
# 							  edge flags
# --------------------------------------------------------------------

# edge protection:
BORDER_PROTECTION  = 1
SCISSOR_PROTECTION = 2
CONTOUR_PROTECTION = 4
ALL_PROTECTION     = 7 # or'ed values of the ones above

# delaunay:
CONTOUR_SEGMENT = 256
WEAK_CHORD      = 512
IS_BARB         = 1024

# tools/IntelligentScissors:
CURRENT_CONTOUR = 2048

# --------------------------------------------------------------------
# 							  face flags
# --------------------------------------------------------------------

# delaunay:
OUTER_FACE = 1

# tools/ActivePaintbrush:
PROTECTED_FACE = 2 # flag which leads to CONTOUR_PROTECTION
