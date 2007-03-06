_cvsVersion = "$Id$" \
              .split(" ")[2:-2]

from hourglass import composeTangentLists, Polygon
from statistics import dartTangents

class Path(list):
    def reverse(self):
        """Reverses a path represented as a list of Dart objects.
        (list.reverse() is called and all darts are switched with nextAlpha().)"""
        
        prevNode = None
        for dart in self:
            if prevNode:
                assert dart.startNode() == prevNode
            prevNode = dart.endNode()
            dart.nextAlpha()
        list.reverse(self)

    def tangents(self):
        """Returns a composed tangentList for all darts in the given path."""
        result = []
        for dart in self:
            result.append(dartTangents(dart))
        return composeTangentLists(result)

    def points(self):
        """Returns an iterator over all points in this path.
        (Skipping the first points of all darts except the first,
        since they are supposed to be duplicates.)"""
        
        yield self[0][0]
        for dart in self:
            pit = iter(dart); pit.next() # skip first point
            for point in pit:
                yield point

    def polygon(self):
        """Returns a polygon containing all points of this path."""
        return Polygon(list(self.points()))

    def length(self):
        """Returns length of this path (sum over edge.length())."""
        result = 0.0
        for dart in self:
            result += dart.edge().length()
        return result

    def __getslice__(self, *args):
        return self.__class__(list.__getslice__(self, *args))

def contour(anchor):
    """contour(anchor) -> Path

    Returns a Path object representing the contour starting with the
    Dart anchor.  (syntactic sugar for 'Path(anchor.phiOrbit())')"""

    return Path(anchor.phiOrbit())

# --------------------------------------------------------------------
#                          Path comparison
# --------------------------------------------------------------------

def _pathString(path):
    """encodes a dart sequence as unicode string (see _pathDecode)"""
    result = []
    for d in path:
        if type(d) != int:
            d = d.label()
        if d < 0:
            d += 0x100000
        result.append(unichr(d))
    return u"".join(result)

def _pathDecode(pathString):
    """decodes a unicode string to a dart sequence (see _pathString)"""
    result = []
    for d in pathString:
        d = ord(d)
        if d > 0x080000:
            d -= 0x100000
        result.append(d)
    return result

import difflib
def pathCompare(p1, p2):
    """Uses difflib.SequenceMatcher to return the differences between two paths"""
    sm = difflib.SequenceMatcher(None, _pathString(p1), _pathString(p2))
    return sm.get_opcodes()

#     result = []
#     for tag, i1, i2, j1, j2 in sm.get_opcodes():
#         result.append((tag, p1[i1:i2], p2[j1:j2]))
#     return result

# --------------------------------------------------------------------
#                            path enumeration
# --------------------------------------------------------------------

def allContinuations(startStart, length):
    """allContinuations(startStart, length)

    Returns a list of Path objects representing all possible paths of
    the given length (number of darts) starting with the given dart that
    
    * do not contain an edge with BORDER_PROTECTION and
    
    * do not contain a direct pair of opposite darts
      (i.e. loops are allowed, but no"U-turns")."""

    assert length >= 1, "allContinuations: length must be >= 1 (each Path must at least contain the given startStart)"

    result = []
    _fillContinuations(startStart, length, [], result)
    return result

from flag_constants import BORDER_PROTECTION

def _fillContinuations(prefix, dart, length, allPaths):
    currentPath = Path(prefix)
    currentPath.append(dart.clone())

    if length == 1:
        allPaths.append(currentPath)
        return

    comingFrom = dart.clone().nextAlpha()
    it = comingFrom.sigmaOrbit(); it.next()
    for neighbor in it:
        if not neighbor.edge().flag(BORDER_PROTECTION):
            _fillContinuations(currentPath, neighbor, length-1, allPaths)
