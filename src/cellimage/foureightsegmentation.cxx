#include "foureightsegmentation.hxx"
#include "cellconfigurations.hxx"

#include <iostream>
#include <set>
#include "mydebug.hxx"
#include "debugimage.hxx"

namespace vigra {

namespace cellimage {

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

FourEightSegmentation::FourEightSegmentation(const CellImage &importImage)
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

    if(maxNodeLabel+maxEdgeLabel+maxFaceLabel>25000)
        std::cerr << maxNodeLabel/nodeCount_;

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

unsigned int FourEightSegmentation::findContourComponent(
    const ContourComponents &contours,
    const DartTraverser & dart)
{
    unsigned int result = 0;

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

    //bool removedEdgeIsLoop = node1.label == node2.label;

    // find indices of contour components to be merged
    unsigned int contour1 = findContourComponent(face1.contours, removedDart);
    unsigned int contour2 = findContourComponent(face2.contours, removedDart);

    // re-use an old anchor for the merged contour
    DartTraverser newAnchor(face1.contours[contour1]);
    if(newAnchor.edgeLabel() == edge.label)
    {
        newAnchor.nextPhi();
        if(newAnchor.edgeLabel() == edge.label)
        {
            newAnchor = face2.contours[contour2];
            if(newAnchor.edgeLabel() == edge.label)
            {
                newAnchor.nextPhi();
                //if(newAnchor.edgeLabel() == edge.label)
                //    newAnchor."makeSingular"();
            }
        }
    }

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

    face1.contours[contour1] = newAnchor;

    // turn anchors if they pointed to the removed edge and there's another left
    if(--node1.degree && node1.anchor.isSingular())
        node1.anchor.carefulNextSigma();
    if(--node2.degree && node2.anchor.isSingular())
        node2.anchor.carefulNextSigma();

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
    unsigned int contourIndex =
        findContourComponent(face.contours, dart);

    // relabel cell in cellImage:
    for(CellScanIterator it= edgeScanIterator(edge.label, cells, false);
        it.inRange(); ++it)
        *it= CellPixel(CellTypeRegion, face.label);

    // update bounds:
    face.bounds |= edge.bounds;

    face.contours[contourIndex] = newAnchor1;
    face.contours.push_back(newAnchor2);

    // turn anchors if they pointed to the removed edge and there's another left
    if(--node1.degree && node1.anchor.isSingular())
        node1.anchor.carefulNextSigma();
    if(--node2.degree && node2.anchor.isSingular())
        node2.anchor.carefulNextSigma();

    // FIXME: if(!isSingular()) check if node was anchor of
    // contour, then update entries in all neighbored faces if
    // DartTraverser pointed to the removed edge..

    edge.uninitialize();
    --edgeCount_;

    return face;
}

void debugDart(const FourEightSegmentation::DartTraverser &dart)
{
	const char *dirStr[] = {
		"East",
		"NorthEast",
		"North",
		"NorthWest",
		"West",
		"SouthWest",
		"South",
		"SouthEast"
	};

	vigra::Point2D pos(dart.neighborCirculator().center() -
					   dart.segmentation()->cells);

	std::cerr << "DartTraverser pointing "
			  << dirStr[(int)dart.neighborCirculator().direction()]
			  << " from " << pos.x << ", " << pos.y;
	if(dart.neighborCirculator().center()->type() == CellTypeVertex)
	{
		CellLabel nodeLabel(dart.neighborCirculator().center()->label());
		std::cerr << " (node " << nodeLabel;
		if(nodeLabel > dart.segmentation()->maxNodeLabel())
			std::cerr << ", LABEL INVALID!";
		else
		{
			const FourEightSegmentation::NodeInfo &node =
				dart.segmentation()->node(nodeLabel);
			if(!node.initialized())
				std::cerr << ", UNINITIALIZED";
			std::cerr << ", " << node.size << " pixels, degree " << node.degree;
		}
		std::cerr << ")";
	}
	else
	{
		if(dart.neighborCirculator().center()->type() == CellTypeLine)
			std::cerr << " (CellTypeLine, label ";
		else
			std::cerr << " (CellTypeRegion, label ";
		std::cerr << dart.neighborCirculator().center()->label() << ")";
	}
	std::cerr << " to ";
	if(dart.neighborCirculator().base()->type() == CellTypeLine)
	{
		CellLabel edgeLabel(dart.neighborCirculator().base()->label());
		std::cerr << "edge " << edgeLabel;
		if(edgeLabel > dart.segmentation()->maxEdgeLabel())
			std::cerr << ", LABEL INVALID!";
		else
		{
			const FourEightSegmentation::EdgeInfo &edge =
				dart.segmentation()->edge(edgeLabel);
			if(!edge.initialized())
				std::cerr << ", UNINITIALIZED";
			std::cerr << ", " << edge.size << " pixels";
		}
		std::cerr << "\n";
	}
	else
	{
		if(dart.neighborCirculator().base()->type() == CellTypeVertex)
			std::cerr << "CellTypeVertex pixel, label ";
		else
			std::cerr << "CellTypeRegion pixel, label ";
		std::cerr << dart.neighborCirculator().base()->label() << "\n";
	}
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

/********************************************************************/

FourEightSegmentation &FourEightSegmentation::deepCopy(
    const FourEightSegmentation &other)
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

    return *this;
}

/********************************************************************/

void FourEightSegmentation::initCellImage(BImage & contourImage)
{
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

                if(cellConfigurations[conf] == CellTypeError)
                {
                    debugImage(crop(srcImageRange(contourImage),
                                    Rect2D(x, y, x+5, y+5)),
                               std::cerr);
                    char message[200];
                    snprintf(message, 200, "FourEightSegmentation::init(): "
                            "Configuration at (%d, %d) must be thinned further (found configuration %d)",
                            x, y, conf);

                    vigra_precondition(0, message);
                }

                cell->setType(cellConfigurations[conf], 0);
            }
        }
    }
}

/********************************************************************/

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
                        sprintf(msg, "label0Cells(): "
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

/********************************************************************/

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

    std::cerr << "found maxEdgeLabel " << maxEdgeLabel << "\n";
    return maxEdgeLabel;
}

/********************************************************************/

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

/********************************************************************/

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

/********************************************************************/

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

                nodeList_[index].degree = 0;
                nodeList_[index].size = 1;
                nodeList_[index].anchor = DartTraverser(
                    this, CellImageEightCirculator(cell,
                                                   EightNeighborCode::West));
                if(!nodeList_[index].anchor.isSingular())
                {
                    DartTraverser dart(nodeList_[index].anchor);
                    do
                    {
                        ++nodeList_[index].degree;
                    }
                    while(dart.nextSigma() != nodeList_[index].anchor);
                }

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
                // calculate area from counting the pixels of the node
                ++nodeList_[index].size;
            }
            nodeList_[index].bounds |= pos;
        }
    }

    for(NodeIterator node= nodesBegin(); node.inRange(); ++node)
	{
        // methods to calculate the area must yield identical values
        if(crackCirculatedAreas[node->label] != node->size)
        {
//             std::cerr << "FourEightSegmentation::initNodeList(): "
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

void FourEightSegmentation::initEdgeList(CellLabel maxEdgeLabel)
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
			vigra_precondition(index < edgeList_.size(),
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
                std::cerr << "dartEnd: "; debugDart(dartEnd);
                std::cerr << "dart:    "; debugDart(dart);
            }
        }
        while(dart.nextSigma() != dartEnd);
    }
}

/********************************************************************/

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
