#include "foureightsegmentation.hxx"
#include "cellconfigurations.hxx"

#include <iostream>
#include "mydebug.hxx"
#include "debugimage.hxx"

namespace vigra {

namespace cellimage {

unsigned char FourEightSegmentation::findContourComponent(
    const ContourComponents &contours,
    const DartTraverser & dart)
{
    unsigned char result = 0;

    if(contours.size() == 1)
        return result;

    if(dart.isSingular())
    {
        for(ConstContourComponentsIterator contour= contours.begin();
            contour != contours.end(); ++contour, ++result)
            if(contour->startNodeLabel() == dart.startNodeLabel())
                return result;
    }

    for(ConstContourComponentsIterator contour= contours.begin();
        contour != contours.end(); ++contour, ++result)
        if(contour->edgeLabel() == dart.edgeLabel())
            return result;

    // argl, we have to circulate through all contours now..
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

    vigra_fail("findContourComponent: dart not found in any contour!");
    return 0;
}

void FourEightSegmentation::removeNodeFromContours(ContourComponents &contours,
                                                   CellLabel nodeLabel)
{
    // find contour anchor to be changed if it points to this edge
    for(ContourComponentsIterator contour= contours.begin();
        contour != contours.end(); ++contour)
        if(contour->startNodeLabel() == nodeLabel)
            contour->nextPhi();
}

FourEightSegmentation::FaceInfo &FourEightSegmentation::removeIsolatedNode(
    const DartTraverser & dart)
{
    vigra_precondition(dart.isSingular(),
                       "removeIsolatedNode: node is not singular");

    NodeInfo &node= dart.startNode();
    FaceInfo &face= dart.leftFace();

    face.contours.erase(face.contours.begin() +
                        findContourComponent(face.contours, dart));

    for(CellScanIterator it= nodeScanIterator(node.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face.label);

    // update bounds:
    face.bounds |= node.bounds;
    face.size += node.size;

    node.uninitialize();
    --nodeCount_;

    return face;
}

FourEightSegmentation::FaceInfo &FourEightSegmentation::mergeFaces(
    const DartTraverser & dart)
{
    // merge smaller face into larger one:
    DartTraverser removedDart = dart;
    if(dart.leftFace().size < dart.rightFace().size)
        removedDart.nextAlpha();

    EdgeInfo &edge= removedDart.edge();
    FaceInfo &face1= removedDart.leftFace();
    FaceInfo &face2= removedDart.rightFace();
    NodeInfo &node1= edge.start.startNode();
    NodeInfo &node2= edge.end.startNode();

    vigra_precondition(face1.label != face2.label,
                       "FourEightSegmentation::mergeFaces(): edge is a bridge");

    bool removedEdgeIsLoop = node1.label == node2.label;

    // find indices of contour components to be merged
    unsigned char contour1 = findContourComponent(face1.contours, removedDart);
    unsigned char contour2 = findContourComponent(face2.contours, removedDart);

    // re-use an old anchor for the merged contour
    DartTraverser newAnchor(
        face1.contours[contour1].startNodeLabel() <
        face2.contours[contour2].startNodeLabel() ?
        face1.contours[contour1] : face2.contours[contour2]);
    if(!removedEdgeIsLoop && (newAnchor.edgeLabel() == edge.label))
        newAnchor.prevSigma();

    // relabel cells in cellImage:
    for(CellScanIterator it= edgeScanIterator(edge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face1.label);
    for(CellScanIterator it= faceScanIterator(face2.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face1.label);

    // update bounds:
    face1.bounds |= edge.bounds;
    face1.size += edge.size;
    face1.bounds |= face2.bounds;
    face1.size += face2.size;

    // update contours
    for(unsigned int i = 0; i < face2.contours.size(); i++)
        if(i != contour2)
            face1.contours.push_back(face2.contours[i]);
    if(removedEdgeIsLoop)
        newAnchor.recheckSingularity();
    face1.contours[contour1] = newAnchor;

    node1.anchor.recheckSingularity();
    node2.anchor.recheckSingularity();

    edge.uninitialize();
    --edgeCount_;
    face2.uninitialize();
    --faceCount_;

    return face1;
}

FourEightSegmentation::FaceInfo &FourEightSegmentation::removeBridge(
    const DartTraverser & dart)
{
    EdgeInfo &edge= dart.edge();
    FaceInfo &face= dart.leftFace();
    NodeInfo &node1= edge.start.startNode();
    NodeInfo &node2= edge.end.startNode();

    vigra_precondition(face.label == dart.rightFaceLabel(),
                       "FourEightSegmentation::removeBridge(): edge is not a bridge");

    DartTraverser newAnchor1(dart);
    newAnchor1.prevSigma();
    DartTraverser newAnchor2(dart);
    newAnchor2.nextAlpha();
    newAnchor2.prevSigma();

    // find index of contour component to be changed
    unsigned char contourIndex =
        findContourComponent(face.contours, dart);

    // relabel cell in cellImage:
    for(CellScanIterator it= edgeScanIterator(edge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face.label);

    // update bounds:
    face.bounds |= edge.bounds;

    face.contours[contourIndex] = newAnchor1;
    face.contours.push_back(newAnchor2);

    node1.anchor.recheckSingularity();
    node2.anchor.recheckSingularity();

    // FIXME: if(!recheckSingularity()) check if node was anchor of
    // contour, then update entries in all neighbored faces if
    // DartTraverser pointed to the removed edge..

    edge.uninitialize();
    --edgeCount_;

    return face;
}

FourEightSegmentation::EdgeInfo &FourEightSegmentation::mergeEdges(
    const DartTraverser & dart)
{
    // merge smaller edge (edge2) into larger one (edge1):
    DartTraverser dart1(dart);
    dart1.nextSigma();
    bool firstIsSmaller =
        dart.edge().bounds.area() < dart1.edge().bounds.area();

    if(firstIsSmaller)
        dart1.nextSigma(); // dart1 _temporarily_ points to smaller one

    DartTraverser dart2(dart1);
    EdgeInfo &edge2 = dart2.edge();
    NodeInfo &node = dart1.startNode();
    dart1.nextSigma();
    EdgeInfo &edge1 = dart1.edge();

    vigra_precondition((firstIsSmaller ? dart2 : dart1) == dart,
        "FourEightSegmentation::mergeEdges(): node has degree > 2!");

    vigra_precondition(edge1.label != edge2.label,
        "FourEightSegmentation::mergeEdges(): node has degree one or is loop!");

    // update contours of neighbored faces if necessary
    removeNodeFromContours(dart1.leftFace().contours, dart1.startNodeLabel());
    if(dart1.leftFaceLabel() != dart1.rightFaceLabel())
        removeNodeFromContours(dart1.rightFace().contours, dart1.startNodeLabel());

    dart1.nextAlpha();
    dart2.nextAlpha();

    // relabel cells in cellImage:
    for(CellScanIterator it= nodeScanIterator(node.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeLine, edge1.label);
    for(CellScanIterator it= edgeScanIterator(edge2.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeLine, edge1.label);

    // update bounds:
    edge1.bounds |= node.bounds;
    edge1.size += node.size;
    edge1.bounds |= edge2.bounds;
    edge1.size += edge2.size;

    // update start, end
    if(dart1.startNodeLabel() < dart2.startNodeLabel())
    {
        edge1.start = dart1;
        edge1.end = dart2;
    }
    else
    {
        edge1.start = dart2;
        edge1.end = dart1;
    }

    node.uninitialize();
    --nodeCount_;
    edge2.uninitialize();
    --edgeCount_;

    return edge1;
}

void FourEightSegmentation::initCellImage(BImage & contourImage)
{
    BImage::traverser rawLine = contourImage.upperLeft() + Diff2D(1,1);
    CellImage::traverser cellLine = cellImage.upperLeft() + Diff2D(1,1);
    for(int y=-1; y < cellImage.height()-3; ++y, ++rawLine.y, ++cellLine.y)
    {
        BImage::traverser raw = rawLine;
        CellImage::traverser cell = cellLine;
        for(int x=-1; x < cellImage.width()-3; ++x, ++raw.x, ++cell.x)
        {
            if(*raw == 0)
            {
                cell->setType(CellTypeRegion);
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

                if(cellConfigurations[conf] == CellTypeError)
                {
                    debugImage(crop(srcImageRange(contourImage),
                                    Rect2D(x, y, x+5, y+5)),
                               std::cerr);
                    char message[200];
                    sprintf(message, "FourEightSegmentation::init(): "
                            "Configuration at (%d, %d) must be thinned further (found configuration %d)",
                            x, y, conf);

                    vigra_precondition(0, message);
                }

                cell->setType(cellConfigurations[conf]);
            }
        }
    }
}

// -------------------------------------------------------------------

CellLabel FourEightSegmentation::label0Cells()
{
    BImage nodeImage(cellImage.size());
    BImage::traverser nodes = nodeImage.upperLeft() + Diff2D(2,2);

    for(int y=-2; y < cellImage.height()-2; ++y)
    {
        CellImage::traverser cell = cells + Diff2D(-2, y);

        for(int x=-2; x < cellImage.width()-2; ++x, ++cell.x)
        {
            if(cell->type() == CellTypeVertex)
            {
                nodes(x,y) = 1;

                // test for forbidden configuration
                CellImageEightCirculator n(cell);
                CellImageEightCirculator nend = n;

                do
                {
                    if(n->type() == CellTypeLine && n[1].type() == CellTypeLine)
                    {
                        char msg[200];
                        sprintf(msg, "initFourEightSegmentation(): "
                                "Node at (%d, %d) has two incident edgels from the same edge (direction: %d)",
                                x, y, n - nend);
                        vigra_precondition(0, msg);
                    }
                }
                while(++n != nend);
            }
            else
            {
                nodes(x,y) = 0;
            }
        }
    }

    return labelImageWithBackground(
        srcImageRange(nodeImage),
        destImage(cellImage, LabelWriter<CellTypeVertex>()), true, 0);
}

// -------------------------------------------------------------------

CellLabel FourEightSegmentation::label1Cells(CellLabel maxNodeLabel)
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

// -------------------------------------------------------------------

CellLabel FourEightSegmentation::label2Cells(BImage & contourImage)
{
    // labelImageWithBackground() starts with label 1, so don't
    // include outer border (infinite regions shall have label 0)
    return labelImageWithBackground(
        srcIterRange(contourImage.upperLeft() + Diff2D(1,1),
                     contourImage.lowerRight() - Diff2D(1,1),
                     contourImage.accessor()),
        destIter(cellImage.upperLeft() + Diff2D(1,1),
                 LabelWriter<CellTypeRegion>()),
        false, 1);
}

// -------------------------------------------------------------------

void FourEightSegmentation::labelCircles(
    CellLabel & maxNodeLabel, CellLabel & maxEdgeLabel)
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

// -------------------------------------------------------------------

void FourEightSegmentation::labelEdge(CellImageEightCirculator rayAtStart,
                                      CellLabel newLabel)
{
    EdgelIterator edge(rayAtStart);

    // follow the edge and relabel it
    for(; !edge.atEnd(); ++edge)
    {
        edge->setLabel(newLabel, CellTypeLine);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initNodeList(CellLabel maxNodeLabel)
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

            CellLabel index = cell->label();
            vigra_precondition(index < nodeList_.size(),
                               "nodeList_ must be large enough!");

            if(!nodeList_[index].initialized())
            {
                nodeList_[index].label = index;
                ++nodeCount_;

                nodeList_[index].centerX = pos.x;
                nodeList_[index].centerY = pos.y;
                nodeList_[index].size = 1;
                nodeList_[index].anchor = DartTraverser(
                    this, CellImageEightCirculator(cell,
                                                   EightNeighborCode::West));

                // calculate area from following the outer contour of the node
                CrackContourCirculator<CellImage::traverser> crack(cell);
                CrackContourCirculator<CellImage::traverser> crackend(crack);
                do
                {
                    crackCirculatedAreas[index] += crack.diff().x * crack.pos().y -
                                                   crack.diff().y * crack.pos().x;
                }
                while(++crack != crackend);

                crackCirculatedAreas[index] /= 2;
            }
            else
            {
                nodeList_[index].centerX += pos.x;
                nodeList_[index].centerY += pos.y;

                // calculate area from counting the pixels of the node
                ++nodeList_[index].size;
            }
            nodeList_[index].bounds |= pos;
        }
    }

    for(NodeIterator node= nodesBegin(); node.inRange(); ++node)
	{
		node->centerX /= node->size;
        node->centerY /= node->size;

        // methods to calculate the area must yield identical values
        if(crackCirculatedAreas[node->label] != node->size)
        {
            std::cerr << "FourEightSegmentation::initNodeList(): "
                      << "Node " << node->label << " at ("
                      << node->centerX << ", "
                      << node->centerY << ") has a hole:\n";
            Point2D center((int)(node->centerX+.5), (int)(node->centerY+.5));
            debugImage(crop(srcImageRange(cellImage),
                            Rect2D(center.x, center.y, center.x+5, center.y+5)),
                       std::cerr, 4);
        }
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initEdgeList(CellLabel maxEdgeLabel)
{
    edgeList_.resize(maxEdgeLabel + 1);

    NodeIterator n(nodeList_.begin(), nodeList_.end());

    for(; n.inRange(); ++n)
    {
        DartTraverser dart = n->anchor;
        if(dart.isSingular())
            continue;
        DartTraverser dartEnd = dart;

        do
        {
            CellLabel label = dart.edgeLabel();
            vigra_precondition(label < edgeList_.size(),
                               "edgeList_ must be large enough!");
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
        }
        while(dart.nextSigma() != dartEnd);
    }
}

// -------------------------------------------------------------------

void FourEightSegmentation::initFaceList(
    BImage & contourImage, CellLabel maxFaceLabel)
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
    faceList_[0].bounds |= Point2D(-2, -2) + cellImage.size();
    faceList_[0].size = NumericTraits<unsigned int>::max();
    //faceList_[0].anchor = Diff2D(-2, -2); FIXME: are the contour-anchors enough?
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

            CellLabel index = cell->label();
            vigra_precondition(index < faceList_.size(),
                               "faceList_ must be large enough!");

            if(!faceList_[index].initialized())
            {
                //std::cerr << "found face " << index << " at " << pos.x << "," << pos.y << "\n";
                faceList_[index].label = index;
                ++faceCount_;
                //faceList_[index].anchor = pos;
                faceList_[index].bounds |= pos;
                faceList_[index].size = 1;

                // find incident node
                if(leftNeighbor->type() == CellTypeVertex)
                {
                    vigra_precondition(leftNeighbor->type() == CellTypeVertex,
                                       "leftNeighbor expected to be a vertex");

                    DartTraverser anchor(this,
                                         CellImageEightCirculator(leftNeighbor));
                    anchor.prevSigma();

                    vigra_invariant(anchor.leftFaceLabel() == index,
                                    "FourEightSegmentation::initFaceList()");

                    faceList_[index].contours.push_back(anchor);
                }
                else
                {
                    vigra_precondition(leftNeighbor->type() == CellTypeLine,
                                       "leftNeighbor expected to be an edge");

                    CellLabel edgeIndex = leftNeighbor->label();

                    vigra_precondition(edgeList_[edgeIndex].initialized(),
                                       "EdgeInfo expected to be initialized");

                    DartTraverser anchor = edgeList_[edgeIndex].start;
                    if(anchor.leftFaceLabel() != index)
                        anchor.nextAlpha();

                    vigra_invariant(anchor.leftFaceLabel() == index,
                                    "FourEightSegmentation::initFaceList()");

                    faceList_[index].contours.push_back(anchor);
                }
            }
            else
            {
                faceList_[index].bounds |= pos;
                ++faceList_[index].size;

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

                        vigra_invariant(anchor.leftFaceLabel() == index,
                                        "FourEightSegmentation::initFaceList()");

                        faceList_[index].contours.push_back(anchor);
                    }
                    else
                    {
                        vigra_precondition(neighbor->type() == CellTypeLine,
                                           "neighbor expected to be an edge");

                        CellLabel edgeIndex = neighbor->label();

                        vigra_precondition(edgeList_[edgeIndex].initialized(),
                                           "EdgeInfo should be initialized");

                        DartTraverser anchor = edgeList_[edgeIndex].start;
                        if(anchor.leftFaceLabel() != index)
                            anchor.nextAlpha();

                        vigra_invariant(anchor.leftFaceLabel() == index,
                                        "FourEightSegmentation::initFaceList()");

                        faceList_[index].contours.push_back(anchor);
                    }
                }
                while(++neighbor != nend);
            }
        }
    }
}

} // namespace cellimage

} // namespace vigra
