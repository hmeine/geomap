#include "mydebug.hxx"
#include "debugimage.hxx"
#include "crop.hxx"

#include "foureightsegmentation.hxx"
#include "cellconfigurations.hxx"

#include <boost/format.hpp>

#include <iostream>
#include <algorithm>
#include <set>

#if !defined(NDEBUG) && !defined(_MSC_VER)
#  warning Consistency checks will be done after every Euler operation!
#endif

namespace vigra {

namespace cellimage {

#ifndef NDEBUG
inline void validateDart(const GeoMap::DartTraverser &dart)
{
    vigra_precondition(dart.neighborCirculator().center()->type() == CellTypeVertex,
                       "dart is not attached to a node");
    vigra_precondition(dart.startNode().initialized(),
                       "dart's startNode is not valid (initialized())");
    if(!dart.isSingular())
        vigra_precondition(dart.edge().initialized(),
                           "dart's edge is not valid (initialized())");
}
#else
// separate def. to prevent compiler warning due to unused param:
inline void validateDart(const GeoMap::DartTraverser &) {}
#endif

struct FindMaxLabel
{
    CellLabel maxLabel_;

    FindMaxLabel() : maxLabel_(0) {}

    void operator()(const CellPixel &p)
    {
        if(p.label() > maxLabel_)
            maxLabel_ = p.label();
    }

    CellLabel operator()() const
    { return maxLabel_; }
};

GeoMap::GeoMap(const CellImage &importImage)
: initialized_(false)
{
    cellImage.resize(importImage.size());
    copyImage(srcImageRange(importImage), destImage(cellImage));

    cells = cellImage.upperLeft() + Diff2D(2, 2);

    FindMaxLabel findMaxNodeLabel;
    inspectImageIf(srcImageRange(cellImage),
                   maskImage(cellImage, CellTypeEquals<CellTypeVertex>()),
                   findMaxNodeLabel);
    CellLabel maxNodeLabel = findMaxNodeLabel();
    std::cerr << "  found maxNodeLabel: " << maxNodeLabel << "\n";

    FindMaxLabel findMaxEdgeLabel;
    inspectImageIf(srcImageRange(cellImage),
                   maskImage(cellImage, CellTypeEquals<CellTypeLine>()),
                   findMaxEdgeLabel);
    CellLabel maxEdgeLabel = findMaxEdgeLabel();
    std::cerr << "  found maxEdgeLabel: " << maxEdgeLabel << "\n";

    FindMaxLabel findMaxFaceLabel;
    inspectImageIf(srcImageRange(cellImage),
                   maskImage(cellImage, CellTypeEquals<CellTypeRegion>()),
                   findMaxFaceLabel);
    CellLabel maxFaceLabel = findMaxFaceLabel();
    std::cerr << "  found maxFaceLabel: " << maxFaceLabel << "\n";

    nodeCount_ = edgeCount_ = faceCount_ = 0;

//    if(maxNodeLabel+maxEdgeLabel+maxFaceLabel>25000)
//        std::cerr << maxNodeLabel/nodeCount_;

    std::cerr << "  initializing nodeList..\n";
    initNodeList(maxNodeLabel);
    std::cerr << "  initializing edgeList..\n";
    initEdgeList(maxEdgeLabel);

    std::cerr << "  creating temporary contourImage..\n";
    BImage contourImage(cellImage.size());
    contourImage = 1;
    initImageIf(destImageRange(contourImage),
                maskImage(cellImage, CellTypeEquals<CellTypeRegion>()),
                0);
    std::cerr << "  initializing faceList..\n";
    initFaceList(contourImage, maxFaceLabel);

    initialized_ = true;
}

unsigned int GeoMap::findComponentAnchor(
    const FaceInfo &face, const DartTraverser & dart)
{
    unsigned int result = 0;
    const ContourComponents &contours = face.contours;

    if(contours.size() == 1)
        return result;

    if(dart.isSingular())
    {
        // look for startNodeLabel
        for(ConstContourComponentsIterator contour= contours.begin();
            contour != contours.end(); ++contour, ++result)
            if(contour->startNodeLabel() == dart.startNodeLabel())
                return result;
    }
    else
    {
        // look for edgeLabel
        for(ConstContourComponentsIterator contour= contours.begin();
            contour != contours.end(); ++contour, ++result)
            if(!contour->isSingular() &&
               (contour->edgeLabel() == dart.edgeLabel()))
                return result;

        // we have to circulate through all contours now.. :-(
        result = 0;
        for(ConstContourComponentsIterator contour= contours.begin();
            contour != contours.end(); ++contour, ++result)
        {
            DartTraverser dart2(*contour);
            while(dart2.nextPhi() != *contour)
            {
                if(dart2.edgeLabel() == dart.edgeLabel())
                    return result;
            }
        }
    }

    std::cerr << "DART NOT FOUND IN CONTOURS OF FACE " << face.label
              << ": " << dart << "\n";
    result = 1;
    for(ConstContourComponentsIterator contour= contours.begin();
        contour != contours.end(); ++contour, ++result)
    {
        std::cerr << "contour " << result << ":\n";
        DartTraverser dart2(*contour);
        do {
            std::cerr << "  " << dart2 << "\n";
        } while(dart2.nextPhi() != *contour);
    }

    vigra_fail("findComponentAnchor: dart not found in any contour");
    return 0;
}

void GeoMap::removeNodeFromContours(ContourComponents &contours,
                                    CellLabel nodeLabel)
{
    // find contour anchor to be changed if it points to this edge
    for(ContourComponentsIterator contour= contours.begin();
        contour != contours.end(); ++contour)
        if(contour->startNodeLabel() == nodeLabel)
            contour->nextPhi();
}

GeoMap::FaceInfo &GeoMap::removeIsolatedNode(const DartTraverser & dart)
{
    //std::cerr << "removeIsolatedNode(" << dart << ")\n";
    validateDart(dart);
    vigra_precondition(dart.isSingular(),
                       "removeIsolatedNode: node is not singular");

    NodeInfo &node= dart.startNode();
    FaceInfo &face= dart.leftFace();

    face.contours.erase(face.contours.begin() + findComponentAnchor(face, dart));

    for(CellScanIterator it= nodeScanIterator(node.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face.label);

    // updating bounds not necessary since all of the node's neighbors
    // are already face pixels, so the bounds should not change
    face.size += node.size;

    node.uninitialize();
    --nodeCount_;
#ifndef NDEBUG
    try
    {
        checkConsistency();
    }
    catch(vigra::StdException &)
    {
        std::cerr << "OPERATION: removeIsolatedNode(" << dart << ")\n";
        throw;
    }
#endif
    return face;
}

GeoMap::FaceInfo &GeoMap::mergeFaces(const DartTraverser & dart)
{
    //std::cerr << "mergeFaces(" << dart << ")\n";
    validateDart(dart);
    // merge smaller face into larger one:
    DartTraverser removedDart = dart;
    if(dart.leftFace().bounds.area() < dart.rightFace().bounds.area())
        removedDart.nextAlpha();

    EdgeInfo &mergedEdge= removedDart.edge();
    FaceInfo &survivor= removedDart.leftFace();
    FaceInfo &mergedFace= removedDart.rightFace();
    NodeInfo &node1= mergedEdge.start.startNode();
    NodeInfo &node2= mergedEdge.end.startNode();

    vigra_precondition(survivor.label != mergedFace.label,
        "GeoMap::mergeFaces(): dart is singular or edge is a bridge");

    // find indices of contour components to be merged
    const unsigned int contour1(findComponentAnchor(survivor, removedDart));
    const unsigned int contour2(findComponentAnchor(mergedFace, removedDart));

    // re-use an old anchor for the merged contour
    if(survivor.contours[contour1].edgeLabel() == mergedEdge.label)
    {
        survivor.contours[contour1].nextPhi();
        if(survivor.contours[contour1].edgeLabel() == mergedEdge.label)
        {
            survivor.contours[contour1] = mergedFace.contours[contour2];
            if(survivor.contours[contour1].edgeLabel() == mergedEdge.label)
                survivor.contours[contour1].nextPhi();
        }
    }
    vigra_postcondition(
        (survivor.contours[contour1].edgeLabel() != mergedEdge.label)
        || (node1.label == node2.label),
        "did not find an anchor for merged contour!\n");

    // update contours
    for(unsigned int i = 0; i < mergedFace.contours.size(); ++i)
        if(i != contour2)
            survivor.contours.push_back(mergedFace.contours[i]);

    // relabel cells in cellImage:
    for(CellScanIterator it= edgeScanIterator(mergedEdge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, survivor.label);
    for(CellScanIterator it= faceScanIterator(mergedFace.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, survivor.label);

    // turn node anchors if they pointed to the removed edge and
    // there's another left
    if(--node1.degree && node1.anchor.isSingular())
        node1.anchor.carefulNextSigma();
    if(--node2.degree && node2.anchor.isSingular())
        node2.anchor.carefulNextSigma();

    // update bounds and sizes:
    survivor.bounds |= mergedEdge.bounds;
    survivor.size += mergedEdge.size;
    survivor.bounds |= mergedFace.bounds;
    survivor.size += mergedFace.size;

    mergedEdge.uninitialize();
    --edgeCount_;
    mergedFace.uninitialize();
    --faceCount_;
    // FIXME: also update maxFaceLabel / maxEdgeLabel
#ifndef NDEBUG
    try
    {
        checkConsistency();
    }
    catch(vigra::StdException &)
    {
        std::cerr << "OPERATION: mergeFaces(" << dart << ")\n";
        throw;
    }
#endif
    return survivor;
}

GeoMap::FaceInfo &GeoMap::removeBridge(const DartTraverser & dart)
{
    //std::cerr << "removeBridge(" << dart << ")\n";
    validateDart(dart);
    vigra_precondition(!dart.isSingular(), "removeBridge: dart does not point to any edge");
    EdgeInfo &edge= dart.edge();
    FaceInfo &face= dart.leftFace();
    NodeInfo &node1= edge.start.startNode();
    NodeInfo &node2= edge.end.startNode();

    vigra_precondition(face.label == dart.rightFaceLabel(),
                       "GeoMap::removeBridge(): edge is not a bridge");

    // prepare new anchors
    DartTraverser newAnchor1(edge.start);
    newAnchor1.prevSigma();
    DartTraverser newAnchor2(edge.end);
    newAnchor2.prevSigma();

    const unsigned int contourIndex(findComponentAnchor(face, dart));

#ifndef NDEBUG
    vigra_invariant(
        (newAnchor1 != newAnchor2) &&
        (newAnchor1.startNodeLabel() == node1.label) &&
        (newAnchor2.startNodeLabel() == node2.label),
        "GeoMap::removeBridge(): new anchors unusable");

#endif // NDEBUG

    // update anchors
    face.contours[contourIndex] = newAnchor1;
    face.contours.push_back(newAnchor2);

    // relabel cell in cellImage:
    for(CellScanIterator it= edgeScanIterator(edge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face.label);

    // turn node anchors if they pointed to the removed edge and
    // there's another one left
    if(--node1.degree && node1.anchor.isSingular())
        node1.anchor.carefulNextSigma();
    if(--node2.degree && node2.anchor.isSingular())
        node2.anchor.carefulNextSigma();

    // update bounds and size:
    face.bounds |= edge.bounds;
    face.size += edge.size;

    edge.uninitialize();
    --edgeCount_;
#ifndef NDEBUG
    try
    {
        checkConsistency();
    }
    catch(vigra::StdException &)
    {
        std::cerr << "OPERATION: removeBridge(" << dart << ")\n";
        throw;
    }
#endif
    return face;
}

GeoMap::EdgeInfo &GeoMap::mergeEdges(const DartTraverser & dart)
{
    //std::cerr << "mergeEdges(" << dart << ")\n";
    // merge smaller edge (mergedEdge) into larger one (survivor):
    DartTraverser dart1(dart);
    dart1.nextSigma();
    bool firstIsSmaller =
        dart.edge().bounds.area() < dart1.edge().bounds.area();

    if(firstIsSmaller)
        dart1.nextSigma(); // dart1 _temporarily_ points to smaller one

    DartTraverser dart2(dart1);
    EdgeInfo &mergedEdge = dart2.edge();
    NodeInfo &node = dart1.startNode();
    dart1.nextSigma();
    EdgeInfo &survivor = dart1.edge();

    vigra_precondition((firstIsSmaller ? dart2 : dart1) == dart,
        "GeoMap::mergeEdges(): node has degree > 2");
    vigra_precondition(survivor.label != mergedEdge.label,
        "GeoMap::mergeEdges(): node has degree one or is loop");

    // update contours of neighbored faces if necessary
    removeNodeFromContours(dart1.leftFace().contours, dart1.startNodeLabel());
    if(dart1.leftFaceLabel() != dart1.rightFaceLabel())
        removeNodeFromContours(dart1.rightFace().contours, dart1.startNodeLabel());

    // update start, end
    dart1.nextAlpha();
    dart2.nextAlpha();
    survivor.start = dart1;
    survivor.end = dart2;

    // relabel cells in cellImage:
    for(CellScanIterator it= nodeScanIterator(node.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeLine, survivor.label);
    for(CellScanIterator it= edgeScanIterator(mergedEdge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeLine, survivor.label);

    // update bounds:
    survivor.bounds |= node.bounds;
    survivor.size += node.size;
    survivor.bounds |= mergedEdge.bounds;
    survivor.size += mergedEdge.size;

    node.uninitialize();
    --nodeCount_;
    mergedEdge.uninitialize();
    --edgeCount_;
#ifndef NDEBUG
    try
    {
        checkConsistency();
    }
    catch(vigra::StdException &)
    {
        std::cerr << "OPERATION: mergeEdges(" << dart << ")\n";
        throw;
    }
#endif
    return survivor;
}

/********************************************************************/

GeoMap &GeoMap::deepCopy(const GeoMap &other)
{
    cellImage = other.cellImage;
    cells = cellImage.upperLeft() + Diff2D(2,2);

    nodeList_ = other.nodeList_;
    nodeCount_ = other.nodeCount_;
    edgeList_ = other.edgeList_;
    edgeCount_ = other.edgeCount_;
    faceList_ = other.faceList_;
    faceCount_ = other.faceCount_;

    for(NodeIterator it= nodesBegin(); it.inRange(); ++it)
    {
        it->anchor.reparent(this);
    }

    for(EdgeIterator it= edgesBegin(); it.inRange(); ++it)
    {
        it->start.reparent(this);
        it->end.reparent(this);
    }

    for(FaceIterator it= facesBegin(); it.inRange(); ++it)
    {
        for(ContourComponentsIterator contour= it->contours.begin();
            contour != it->contours.end(); ++contour)
        {
            contour->reparent(this);
        }
    }

    initialized_ = other.initialized_; // should be always true

    return *this;
}

/********************************************************************/

void GeoMap::checkConsistency()
{
    bool consistent = true;
    // std::cerr << "checking nodes...\n";
    for(NodeIterator it= nodesBegin(); it != nodesEnd(); ++it)
    {
        // std::cerr << "  node " << it->label << "\r";
        unsigned int realDegree = 0;
        if(it->anchor.neighborCirculator().center()->type() != CellTypeVertex)
        {
            consistent = false;
            std::cerr << "node " << it->label << "'s anchor is broken: "
					  << it->anchor << "\n";
        }
        if(!it->anchor.isSingular())
        {
            DartTraverser dart(it->anchor);
            do { ++realDegree; }
            while(dart.nextSigma() != it->anchor);
        }
        if(realDegree != it->degree)
        {
            consistent = false;
            std::cerr << "node " << it->label << " has degree stored as "
                      << it->degree << " but seems to be of degree "
                      << realDegree << "\n  anchor:"
					  << it->anchor << "\n";
            std::cerr << "  sigma-successors:\n";
            for(DartTraverser dart(it->anchor); dart.nextSigma() != it->anchor; )
            {
                std::cerr << "   " << dart << "\n";
            }
        }

        FindBoundingRectangle actualBounds;
        Point2D coord(0, 0);
        inspectCell(nodeScanIterator(it->label, coord, false), actualBounds);
        Rect2D ab(actualBounds.upperLeft, actualBounds.lowerRight);
        if(ab != it->bounds)
        {
            consistent = false;
            std::cerr << "node " << it->label << "'s bounds are stored as "
                      << it->bounds << " but seem to be " << ab << "\n";
        }
    }

    // std::cerr << "checking edges...\n";
    for(EdgeIterator it= edgesBegin(); it != edgesEnd(); ++it)
    {
        // std::cerr << "  edge " << it->label << "\r";
        FindBoundingRectangle actualBounds;
        Point2D coord(0, 0);
        inspectCell(edgeScanIterator(it->label, coord, false), actualBounds);
        Rect2D ab(actualBounds.upperLeft, actualBounds.lowerRight);
        if(ab != it->bounds)
        {
            consistent = false;
            std::cerr << "edge " << it->label << "'s bounds are stored as "
                      << it->bounds << " but seem to be " << ab << "\n";
        }
        if(it->start.neighborCirculator().center()->type() != CellTypeVertex)
        {
            consistent = false;
            std::cerr << "edge " << it->label << "'s start is broken: "
					  << it->start << "\n";
        }
        if(it->start.edgeLabel() != it->label)
        {
            consistent = false;
            std::cerr << "edge " << it->label << "'s start points to edge "
                      << it->start.edgeLabel() << "\n";
        }
        if(it->end.neighborCirculator().center()->type() != CellTypeVertex)
        {
            consistent = false;
            std::cerr << "edge " << it->label << "'s end is broken: "
					  << it->end << "\n";
        }
        if(it->end.edgeLabel() != it->label)
        {
            consistent = false;
            std::cerr << "edge " << it->label << "'s end points to edge "
                      << it->end.edgeLabel() << "\n";
        }
    }

    // std::cerr << "checking faces...\n";
    for(FaceIterator it= facesBegin(); it != facesEnd(); ++it)
    {
        // std::cerr << "  face " << it->label << "\r";
        FindBoundingRectangle actualBounds;
        Point2D coord(0, 0);
        inspectCell(faceScanIterator(it->label, coord, false), actualBounds);
        Rect2D ab(actualBounds.upperLeft, actualBounds.lowerRight);
        if(ab != it->bounds)
        {
            consistent = false;
            std::cerr << "face " << it->label << "'s bounds are stored as "
                      << it->bounds << " but seem to be " << ab << "\n";
        }

        for(ContourComponentsIterator cc(it->contours.begin());
            cc != it->contours.end(); ++cc)
        {
            if((cc->neighborCirculator().center()->type() != CellTypeVertex)
               || (cc->isSingular() && (cc->startNode().degree > 0)))
            {
                consistent = false;
                std::cerr << "face " << it->label << "'s anchor at position "
                          << (cc - it->contours.begin()) << " is broken: "
						  << *cc << "\n";
            }
            DartTraverser dart(*cc);
            dart.nextSigma();
            dart.prevSigma();
            if(dart != *cc)
            {
                consistent = false;
                std::cerr << "face " << it->label << "'s contours are broken:\n"
                          << "  anchor:  " << *cc << "\n"
						  << "  becomes: " << dart << "\n";
            }
            int maxSteps= 10000;
            do
            {
                if(dart.leftFaceLabel() != it->label)
                {
                    std::cerr << "dart.leftFace() is no longer our face!\n";
                    maxSteps = 0;
                    break;
                }
            }
            while((dart.nextPhi() != *cc) && --maxSteps);
            if(!maxSteps)
            {
                consistent = false;
                std::cerr << "face " << it->label << "'s contours are broken:\n"
                          << "  anchor: " << *cc << "\n"
						  << "  did not return after " << maxSteps << " steps: "
						  << dart << "\n";
            }

//  for(unsigned int i = 1; i < face.contours.size(); ++i)
//      for(unsigned int j = 0; j < i; ++j)
//          if(face.contours[i] == face.contours[j])
//          {
//                 consistent = false;
//                 std::cerr << "face " << face.label << "'s contours contains one dart twice:\n  ";
//              debugDart(face.contours[i]);
//          }
        }
    }

    vigra_postcondition(consistent, "consistency check failed");
}


/********************************************************************/

struct ChooseCellConfiguration
{
    CellType * cornerTypeLut_;
    static CellType preferVertex[6], preferEdge[6];

    ChooseCellConfiguration(CellType cornerType)
    : cornerTypeLut_(cornerType == CellTypeVertex ? preferVertex : preferEdge)
    {}

    CellType operator()(CellType c) const
    {
        return cornerTypeLut_[c];
    }
};

CellType ChooseCellConfiguration::preferVertex[6] = {
    CellTypeRegion, CellTypeLine, CellTypeVertex, CellTypeError, CellTypeVertex, CellTypeError };
CellType ChooseCellConfiguration::preferEdge[6] = {
    CellTypeRegion, CellTypeLine, CellTypeVertex, CellTypeError, CellTypeLine, CellTypeLine };

void GeoMap::initCellImage(BImage & contourImage, CellType cornerType)
{
    CellType cellConf[256];
    std::transform(cellConfigurations, cellConfigurations+256, cellConf,
                   ChooseCellConfiguration(cornerType));

    CellPixel regionPixel(CellTypeRegion, 0);

    BImage::traverser rawLine = contourImage.upperLeft();
    CellImage::traverser cellLine = cellImage.upperLeft();
    for(int y=-2; y < cellImage.height()-2; ++y, ++rawLine.y, ++cellLine.y)
    {
        BImage::traverser raw = rawLine;
        CellImage::traverser cell = cellLine;
        for(int x=-2; x < cellImage.width()-2; ++x, ++raw.x, ++cell.x)
        {
            if(*raw == 0)
            {
                *cell = regionPixel;
            }
            else
            {
                vigra::NeighborhoodCirculator<BImage::traverser, EightNeighborCode>
                    neighbors(raw, EightNeighborCode::SouthEast);
                vigra::NeighborhoodCirculator<BImage::traverser, EightNeighborCode>
                    end = neighbors;
                int conf = 0;
                do
                {
                    conf = (conf << 1) | *neighbors;
                }
                while(--neighbors != end);

                if(cellConf[conf] == CellTypeError)
                {
                    debugImage(crop(srcImageRange(contourImage),
                                    Rect2D(x, y, x+5, y+5)),
                               std::cerr);
                    vigra_precondition(0, (boost::format(
                        "GeoMap::init(): Configuration at (%1%, %2%) must be thinned further "
                        "(found configuration %3%)") % x % y % conf).str());
                }

                cell->setType(cellConf[conf], 0);
            }
        }
    }

    // FIXME: this should be done for all contours, not just the boundary!
    cellImage(1, 1).setType(CellTypeVertex, 0); // FIXME: for faceList_[0].anchor
#ifndef NDEBUG
    std::cerr << "*** RUN-TIME WARNING: slowness results from missing NDEBUG flag during compilation of foureightsegmentation.cxx! ***\n";
#endif
}

/********************************************************************/

CellLabel GeoMap::label0Cells()
{
    BImage nodeImage(cellImage.size());
    BImage::traverser nodes = nodeImage.upperLeft() + Diff2D(2,2);

    for(int y=-2; y < cellImage.height()-2; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-2, y);

        for(int x=-2; x < cellImage.width()-2; ++x, ++cell.x)
        {
            if(cell->type() != CellTypeVertex)
            {
                nodes(x,y) = 0;
            }
            else
            {
                nodes(x,y) = 1;
#if 0
                // test for forbidden configuration
                CellImageEightCirculator n(cell);
                CellImageEightCirculator nend = n;

                do
                {
                    vigra_precondition((n->type() != CellTypeLine || n[1].type() || CellTypeLine),
                                       (boost::format("label0Cells(): Node at (%1%, %2%) has two "
                                                      "incident edgels from the same edge (direction: %3%)")
                                        % x % y % (n - nend)).str());
                }
                while(++n != nend);
#endif
            }
        }
    }

    return labelImageWithBackground(
        srcImageRange(nodeImage),
        destImage(cellImage, LabelWriter<CellTypeVertex>()), true, 0);
}

/********************************************************************/

CellLabel GeoMap::label1Cells(CellLabel maxNodeLabel)
{
    std::vector<bool> nodeProcessed(maxNodeLabel + 1, false);

    CellLabel maxEdgeLabel = 0;

    for(int y=-1; y < cellImage.height()-3; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, y);

        for(int x=-1; x < cellImage.width()-3; ++x, ++cell.x)
        {
            if(cell->type() != CellTypeVertex)
                continue;
            if(nodeProcessed[cell->label()])
                continue;

            nodeProcessed[cell->label()] = true;

            DartTraverser rayAtStart(
                this, CellImageEightCirculator(cell, EightNeighborCode::West));
            DartTraverser rayEnd = rayAtStart;

            do
            {
                if(rayAtStart.edgeLabel() != 0)
                    continue;

                labelEdge(rayAtStart.neighborCirculator(), ++maxEdgeLabel);
            }
            while(rayAtStart.nextSigma() != rayEnd);
        }
    }

    return maxEdgeLabel;
}

/********************************************************************/

CellLabel GeoMap::label2Cells(BImage & contourImage)
{
    // labelImageWithBackground() starts with label 1, so don't
    // include outer border (infinite regions shall have label 0)
    return labelImageWithBackground(
        srcIterRange(contourImage.upperLeft() + Diff2D(1,1),
                     contourImage.lowerRight() - Diff2D(1,1)),
        destIter(cellImage.upperLeft() + Diff2D(1,1),
                 LabelWriter<CellTypeRegion>()),
        false, // eight_neighbors = false -> 4-neighborhood
        1);    // background: contour value
}

/********************************************************************/

void GeoMap::labelCircles(CellLabel & maxNodeLabel, CellLabel & maxEdgeLabel)
{
    for(int y=-1; y < cellImage.height()-3; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, y);

        for(int x=-1; x < cellImage.width()-3; ++x, ++cell.x)
        {
            if(cell->label() != 0)
                continue;

            // found a circle (not labeled by previous steps)

            // mark first point as node
            (*cell) = CellPixel(CellTypeVertex, ++maxNodeLabel);

            CellImageEightCirculator rayAtStart(cell);
            CellImageEightCirculator rayEnd = rayAtStart;

            do
            {
                if(rayAtStart->type() != CellTypeLine)
                    continue;
                if(rayAtStart->label() != 0)
                    continue;

                labelEdge(rayAtStart, ++maxEdgeLabel);
            }
            while(++rayAtStart != rayEnd);
        }
    }
}

/********************************************************************/

void GeoMap::labelEdge(CellImageEightCirculator rayAtStart,
                       CellLabel newLabel)
{
    EdgelIterator edge(rayAtStart);

    // follow the edge and relabel it
    for(; !edge.atEnd(); ++edge)
    {
        edge->setLabel(newLabel, CellTypeLine);
    }
}

/********************************************************************/

struct HoleRemover
{
    const CellPixel &p_;
    mutable unsigned int size_;

    HoleRemover(const CellPixel &p) : p_(p), size_(0) {}

    const CellPixel &operator()(const CellPixel &) const
    {
        ++size_;
        return p_;
    }

    unsigned int size() const
    {
        return size_;
    }
};

void GeoMap::initNodeList(CellLabel maxNodeLabel)
{
    nodeList_.resize(maxNodeLabel + 1);
    std::vector<unsigned int> crackCirculatedAreas(maxNodeLabel + 1, 0);

    for(Point2D pos= Point2D(-1, -1); pos.y < cellImage.height()-3; ++pos.y)
    {
        CellImage::traverser cell = cells + Diff2D(-1, pos.y);

        for(pos.x=-1; pos.x < cellImage.width()-3; ++pos.x, ++cell.x)
        {
            if(cell->type() != CellTypeVertex)
                continue;

            CellLabel nodeLabel = cell->label();
            vigra_precondition(nodeLabel < nodeList_.size(),
                               "nodeList_ must be large enough");

            if(!nodeList_[nodeLabel].initialized())
            {
                nodeList_[nodeLabel].label = nodeLabel;
                ++nodeCount_;

                nodeList_[nodeLabel].degree = 0;
                nodeList_[nodeLabel].anchor = DartTraverser(
                    this, CellImageEightCirculator(cell,
                                                   EightNeighborCode::West));
                if(!nodeList_[nodeLabel].anchor.isSingular())
                {
                    DartTraverser dart(nodeList_[nodeLabel].anchor);
                    do
                    {
                        ++nodeList_[nodeLabel].degree;
                    }
                    while(dart.nextSigma() != nodeList_[nodeLabel].anchor);
                }

                // calculate area from following the outer contour of the node
                CrackContourCirculator<CellImage::traverser> crack(cell);
                CrackContourCirculator<CellImage::traverser> crackend(crack);
                do
                {
                    crackCirculatedAreas[nodeLabel] += crack.diff().x * crack.pos().y -
                                                   crack.diff().y * crack.pos().x;
                }
                while(++crack != crackend);

                crackCirculatedAreas[nodeLabel] /= 2;
            }

            nodeList_[nodeLabel].bounds |= pos;
            ++nodeList_[nodeLabel].size;
        }
    }

    for(NodeIterator node= nodesBegin(); node.inRange(); ++node)
    {
        // methods to calculate the area must yield identical values
        if(crackCirculatedAreas[node->label] != node->size)
        {
//             std::cerr << "GeoMap::initNodeList(): "
//                       << "Node " << node->label << " has a "
//                       << (crackCirculatedAreas[node->label] - node->size)
//                       << "-pixel hole, stuffing..\n";

            CellImage::traverser bul(cells + node->bounds.upperLeft());
            CellImage::traverser blr(cells + node->bounds.lowerRight());

            CellImage::traverser anchor(bul);
            CellPixel nodePixel(CellTypeVertex, node->label);
            while(*anchor != nodePixel)
                ++anchor.x;

            std::set<CellLabel> outerRegions;

            CrackContourCirculator<CellImage::traverser> crack(anchor);
            CrackContourCirculator<CellImage::traverser> crackend(crack);
            do
            {
                if(crack.outerPixel()->type() == CellTypeRegion)
                    outerRegions.insert(crack.outerPixel()->label());
            }
            while(++crack != crackend);

//             std::cerr << "found " << outerRegions.size() << " outer regions:\n";
//             for(std::set<CellLabel>::iterator it= outerRegions.begin();
//                 it != outerRegions.end(); ++it)
//             {
//                 std::cerr << *it << ", ";
//             }
//             std::cerr << "\n";
//             debugImage(crop(srcIterRange(cells, cells), node->bounds),
//                        std::cerr, 4);

            for(bul = cells + node->bounds.upperLeft(); bul.y < blr.y; ++bul.y)
            {
                for(anchor= bul; anchor.x < blr.x-1; ++anchor.x)
                {
                    if(*anchor == nodePixel)
                    {
                        ++anchor.x;
                        if((anchor->type() == CellTypeRegion) &&
                           !outerRegions.count(anchor->label()))
                        {
//                             std::cerr << "found hole label: " << *anchor;
                            HoleRemover holeRemover(nodePixel);
                            transformImageIf(crop(srcIterRange(cells, cells), node->bounds),
                                             crop(maskIter(cells, CellMask(*anchor)), node->bounds),
                                             crop(destIter(cells), node->bounds),
                                             holeRemover);

                            node->size += holeRemover.size();
//                             std::cerr << " (" << holeRemover.size() << " pixels)\n";

                            if(crackCirculatedAreas[node->label] == node->size)
                                goto HolesStuffed; // double break
                        }
                        else
                            --anchor.x;
                    }
                }
            }

        HolesStuffed:
            ;
        }
    }
}

/********************************************************************/

void GeoMap::initEdgeList(CellLabel maxEdgeLabel)
{
    edgeList_.resize(maxEdgeLabel + 1);

    for(NodeIterator n = nodesBegin(); n.inRange(); ++n)
    {
        DartTraverser dart = n->anchor;
        if(dart.isSingular())
            continue;
        DartTraverser dartEnd = dart;

        int debugCount = 0;
        do
        {
            CellLabel label = dart.edgeLabel();
            vigra_precondition(label < edgeList_.size(),
                               "edgeList_ must be large enough");
            if(!edgeList_[label].initialized())
            {
                edgeList_[label].label = label;
                ++edgeCount_;
                edgeList_[label].start = dart;
                edgeList_[label].end = dart;
                edgeList_[label].end.nextAlpha();
                // correct size and bounds will be collected by initFaceList()
                edgeList_[label].size = 0;
            }

            if(++debugCount == 40)
            {
                std::cerr << "node " << n->label << " has degree > 40??\n"
                          << "  bounds: " << n->bounds << "\n";
                Rect2D dr(n->bounds);
                dr |= Point2D(dartEnd.neighborCirculator().base() - cells);
                dr |= Point2D(dart.neighborCirculator().base() - cells);
                std::cerr << "debugging rect " << dr << " including "
                          << Point2D(dartEnd.neighborCirculator().base() - cells) << "\n";
                debugImage(crop(srcIterRange(cells, cells), dr),
                           std::cerr, 4);
                std::cerr << "dartEnd: " << dartEnd << "\n";
                std::cerr << "dart:    " << dart << "\n";
            }
        }
        while(dart.nextSigma() != dartEnd);
    }
}

/********************************************************************/

void GeoMap::initFaceList(BImage & contourImage, CellLabel maxFaceLabel)
{
    faceList_.resize(maxFaceLabel + 1);

    IImage contourLabelImage(cellImage.size());
    contourLabelImage = 0;
    int contourComponentsCount =
        labelImageWithBackground(srcImageRange(contourImage),
                                 destImage(contourLabelImage), true, 0);
    IImage::traverser contourLabel =
        contourLabelImage.upperLeft() + Diff2D(2, 2);

    std::vector<bool> contourProcessed(contourComponentsCount + 1, false);

    // process outer face
    faceList_[0].label= 0;
    ++faceCount_;
    faceList_[0].bounds |= Point2D(-2, -2);
    faceList_[0].bounds |= Point2D(-2, -2) + cellImage.size() - Diff2D(1, 1);
    faceList_[0].size = NumericTraits<unsigned int>::max();
    DartTraverser anchor(this, CellImageEightCirculator(cells + Diff2D(-1, -1),
                                                        EightNeighborCode::West));
    faceList_[0].contours.push_back(anchor.prevSigma());
    contourProcessed[contourLabel(-1, -1)] = true;

    for(Point2D pos= Point2D(-1, -1); pos.y < cellImage.height()-3; ++pos.y)
    {
        pos.x= -1;
        CellImage::traverser cell = cells + pos;
        CellImage::traverser leftNeighbor = cell - Diff2D(1, 0);

        for(; pos.x < cellImage.width()-3; ++pos.x, ++cell.x, ++leftNeighbor.x)
        {
            if(cell->type() != CellTypeRegion)
            {
                if(cell->type() == CellTypeLine)
                {
                    edgeList_[cell->label()].bounds |= pos;
                    edgeList_[cell->label()].size += 1;
                }
                continue;
            }

            CellLabel faceLabel = cell->label();
            vigra_precondition(faceLabel < faceList_.size(),
                               "faceList_ must be large enough");

            if(!faceList_[faceLabel].initialized())
            {
                //std::cerr << "found face " << faceLabel << " at " << pos.x << "," << pos.y << "\n";
                faceList_[faceLabel].label = faceLabel;
                ++faceCount_;

                // determine anchor (either via node or via edge):
                if(leftNeighbor->type() == CellTypeVertex)
                {
                    DartTraverser anchor(this,
                                         CellImageEightCirculator(leftNeighbor));
                    anchor.prevSigma();

                    vigra_invariant(anchor.leftFaceLabel() == faceLabel,
                                    "GeoMap::initFaceList()");

                    faceList_[faceLabel].contours.push_back(anchor);
                }
                else
                {
                    vigra_precondition(leftNeighbor->type() == CellTypeLine,
                                       "leftNeighbor expected to be an edge");

                    CellLabel edgeIndex = leftNeighbor->label();

                    vigra_precondition(edgeList_[edgeIndex].initialized(),
                                       (boost::format("EdgeInfo %1% expected to be initialized") % edgeIndex).str());

                    DartTraverser anchor = edgeList_[edgeIndex].start;
                    if(anchor.leftFaceLabel() != faceLabel)
                        anchor.nextAlpha();

                    vigra_invariant(anchor.leftFaceLabel() == faceLabel,
                                    "GeoMap::initFaceList()");

                    faceList_[faceLabel].contours.push_back(anchor);
                }
            }
            else // faceLabel seen before
            {
                // look for inner contours
                CellImageEightCirculator neighbor(cell);
                CellImageEightCirculator nend = neighbor;

                do
                {
                    int contourIndex = contourLabel[neighbor.base() - cells];
                    if(contourIndex == 0 || contourProcessed[contourIndex])
                        continue;

                    // found an inner contour
                    contourProcessed[contourIndex] = true;

                    // find incident node
                    if(neighbor->type() == CellTypeVertex)
                    {
                        // this is the node
                        CellImageEightCirculator n = neighbor;
                        n.swapCenterNeighbor();

                        DartTraverser anchor(this, n);
                        anchor.prevSigma();

                        vigra_invariant(anchor.leftFaceLabel() == faceLabel,
                                        "GeoMap::initFaceList()");

                        faceList_[faceLabel].contours.push_back(anchor);
                    }
                    else
                    {
                        vigra_precondition(neighbor->type() == CellTypeLine,
                                           "neighbor expected to be an edge");

                        CellLabel edgeIndex = neighbor->label();

                        vigra_precondition(edgeList_[edgeIndex].initialized(),
                                           (boost::format("EdgeInfo %1% should be initialized") % edgeIndex).str());

                        DartTraverser anchor = edgeList_[edgeIndex].start;
                        if(anchor.leftFaceLabel() != faceLabel)
                            anchor.nextAlpha();

                        vigra_invariant(anchor.leftFaceLabel() == faceLabel,
                                        "GeoMap::initFaceList()");

                        faceList_[faceLabel].contours.push_back(anchor);
                    }
                }
                while(++neighbor != nend);
            }

            faceList_[faceLabel].bounds |= pos;
            ++faceList_[faceLabel].size;
        }
    }
}

std::ostream &
operator<<(std::ostream & out,
		   const vigra::cellimage::GeoMap::DartTraverser & d)
{
    static const char * const dirStr[] = {
        " E",
        "NE",
        " N",
        "NW",
        " W",
        "SW",
        " S",
        "SE"
    };

    vigra::Point2D pos(d.neighborCirculator().center() -
                       d.segmentation()->cells);

    out << "DartTraverser["
        << dirStr[(int)d.neighborCirculator().direction()]
        << " from " << pos.x << ", " << pos.y;

    bool singularDart = false;
    if(d.neighborCirculator().center()->type() == CellTypeVertex)
    {
        CellLabel nodeLabel(d.neighborCirculator().center()->label());
        out << " (node " << nodeLabel;
        if(nodeLabel > d.segmentation()->maxNodeLabel())
            out << ", LABEL INVALID!";
        else
        {
            const GeoMap::NodeInfo &node(d.segmentation()->node(nodeLabel));
            if(!node.initialized())
                out << ", UNINITIALIZED";
            out << ", " << node.size << " pixels, deg. " << node.degree;
            singularDart = (node.degree == 0);
        }
        out << ")";
    }
    else
    {
        if(d.neighborCirculator().center()->type() == CellTypeLine)
            out << " (*EDGE* ";
        else
            out << " (*FACE* ";
        out << d.neighborCirculator().center()->label() << ")";
    }

    out << " to "; /*****************************************************/

    if(d.neighborCirculator().base()->type() == CellTypeLine)
    {
        CellLabel edgeLabel(d.neighborCirculator().base()->label());
        out << "edge " << edgeLabel;
        if(edgeLabel > d.segmentation()->maxEdgeLabel())
            out << ", LABEL INVALID!";
        else
        {
            const GeoMap::EdgeInfo &edge(d.segmentation()->edge(edgeLabel));
            if(!edge.initialized())
                out << ", UNINITIALIZED";
            out << ", " << edge.size << " pixels";
        }
    }
    else
    {
        if(d.neighborCirculator().base()->type() == CellTypeVertex)
            out << "*NODE* ";
        else
            out << (singularDart ? "face " : "*FACE* ");
        out << d.neighborCirculator().base()->label();
    }

    out << "]";
    return out;
}

} // namespace cellimage

} // namespace vigra
