from hourglass import composeTangentLists, Polygon

class Path(list):
    """Represents a dart path.  In addition to being a list of darts,
    it offers some methods to derive information on the whole path
    (tangents(), points(), polygon(), length()) or to modify it
    (reverse(), which also reverses each dart)."""

    def reverse(self):
        """path.reverse()

        Reverses a path represented as a list of Dart objects.
        (list.reverse() is called and all darts are switched with
        nextAlpha().)"""
        
        for dart in self:
            dart.nextAlpha()
        super(Path, self).reverse()

    def tangents(self, dartTangents):
        """path.tangents(dartTangents) -> tangent list

        Returns a composed tangentList for all darts in the given
        path.  Calls dartTangents on each dart and uses
        composeTangentLists to return a single tangent list."""

        result = []
        for dart in self:
            result.append(dartTangents(dart))
        return composeTangentLists(result)

    def points(self):
        """path.points() -> generator function

        Returns an iterator over all points in this path.
        (Skipping the first points of all darts except the first,
        since they are supposed to be duplicates.)"""
        
        yield self[0][0]
        for dart in self:
            pit = iter(dart); pit.next() # skip first point
            for point in pit:
                yield point

    def polygon(self):
        """path.polygon() -> Polygon

        Returns a polygon containing all points of this path."""

        return Polygon(list(self.points()))

    def length(self):
        """path.length() -> float

        Returns length of this path (sum over edge.length())."""

        result = 0.0
        for dart in self:
            result += dart.edge().length()
        return result

    def boundingBox(self):
        """path.boundingBox() -> float

        Returns boundingBox of this path (union of all
        dart.edge().boundingBox()es)."""

        it = iter(self)
        result = it.next().boundingBox()
        for dart in it:
            result |= dart.edge().boundingBox()
        return result

    def __getslice__(self, *args):
        """Re-implement slicing to prevent.. eh.. "slicing". ;-)"""
        return self.__class__(super(Path, self).__getslice__(*args))

    def __add__(self, dart_s):
        """Allows (path + dart) as well as (path1 + path2)."""
        result = self.__class__(self)
        if isinstance(dart_s, list):
            result.extend(dart_s)
        else:
            result.append(dart_s)
        return result

def contour(anchor):
    """contour(anchor) -> Path

    Returns a Path object representing the contour starting with the
    Dart anchor.  (syntactic sugar for 'Path(anchor.phiOrbit())')"""

    return Path(anchor.phiOrbit())

# --------------------------------------------------------------------
#                          Path comparison
# --------------------------------------------------------------------

import difflib

def pathCompare(p1, p2):
    """Uses difflib.SequenceMatcher to return the differences between two paths"""
    sm = difflib.SequenceMatcher(None, _pathString(p1), _pathString(p2))
    return sm.get_opcodes()

#   result = []
#   for tag, i1, i2, j1, j2 in sm.get_opcodes():
#       result.append((tag, p1[i1:i2], p2[j1:j2]))
#   return result

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

# --------------------------------------------------------------------
#                            path enumeration
# --------------------------------------------------------------------

def allContinuations(startDart, length, klass = Path):
    """allContinuations(startDart, length, klass = Path)

    Returns a list of Path objects representing all possible paths of
    the given length (number of darts) starting with the given dart that
    
    * do not contain an edge with BORDER_PROTECTION and
    
    * do not contain a direct pair of opposite darts
      (i.e. loops are allowed, but no "U-turns").

    You can change the class of the returned paths (default: Path)
    with the optional 'klass' argument.  (The given type must support
    copy-construction and append().)"""

    assert length >= 1, "allContinuations: length must be >= 1 (each Path must at least contain the given startDart)"

    result = []
    _fillContinuations(klass(), startDart, length, result)
    return result

from flag_constants import BORDER_PROTECTION

def _fillContinuations(prefix, dart, length, allPaths):
    currentPath = prefix + dart.clone()

    if length == 1:
        allPaths.append(currentPath)
        return

    comingFrom = dart.clone().nextAlpha()
    it = comingFrom.sigmaOrbit(); it.next()
    for neighbor in it:
        if not neighbor.edge().flag(BORDER_PROTECTION):
            _fillContinuations(currentPath, neighbor, length-1, allPaths)
