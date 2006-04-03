from hourglass import composeTangentLists
from maptest import dartTangents

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
