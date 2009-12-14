#include "cppmap.hxx"
#include <vigra/tinyvector.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <algorithm>
#include <cmath>

#ifdef _MSC_VER
inline int isnan(double t) { return _isnan(t); }
#endif

template<class Container>
void removeAll(Container &container,
               const typename Container::value_type &element)
{
    typename Container::iterator it;
    while((it = std::find(container.begin(), container.end(), element))
          != container.end())
        container.erase(it);
}

template<class Container>
void removeOne(Container &container,
               const typename Container::value_type &element)
{
    container.erase(std::find(container.begin(), container.end(), element));
}

/********************************************************************/

double contourArea(const GeoMap::Dart &dart)
{
    double result = 0.0;
    GeoMap::Dart d(dart);
    do
    {
        if(!d.guaranteedEdge()->isBridge())
            result += d.partialArea();
    }
    while(d.nextPhi() != dart);
    return result;
}

double contourLength(const GeoMap::Dart &dart)
{
    double result = 0.0;
    GeoMap::Dart d(dart);
    do
    {
        result += d.edge()->length();
    }
    while(d.nextPhi() != dart);
    return result;
}

double isoperimetricQuotient(const GeoMap::Dart &dart)
{
    double area = 0.0, length = 0.0;

    GeoMap::Dart d(dart);
    do
    {
        if(!d.edge()->isBridge())
            area += d.partialArea();
        length += d.edge()->length();
    }
    while(d.nextPhi() != dart);

    return vigra::squaredNorm(length)/(4*M_PI*area);
}

Polygon contourPoly(const GeoMap::Dart &dart)
{
    Polygon result;
    GeoMap::Dart d(dart);
    do
    {
        if(d.label() < 0)
        {
            Polygon rev(*d.edge());
            rev.reverse();
            result.extend(rev);
        }
        else
            result.extend(*d.edge());
    }
    while(d.nextPhi() != dart);
    return result;
}

/********************************************************************/

template<class InputIterator>
inline InputIterator skipFirst(InputIterator it)
{
    return ++it;
}

template<class InputIterator>
inline InputIterator skipLast(InputIterator it)
{
    return --it;
}

/**
 * atEnd: append 'other' edge? (else prepend)
 * reverse: does the 'other' edge have the opposite direction of this?
 */
void GeoMap::Edge::concatenate(Edge &other, bool atEnd, bool reverse)
{
    if(reverse)
    {
        if(atEnd)
        {
            vigra_assert(this->points_.back() == other.points_.back(),
                         "appending non-matching reversed edge");
            this->points_.insert(
                this->points_.end(),
                skipFirst(other.points_.rbegin()), other.points_.rend());
        }
        else
        {
            vigra_assert(this->points_.front() == other.points_.front(),
                         "appending non-matching reversed edge");
            this->points_.insert(
                this->points_.begin(),
                other.points_.rbegin(), skipLast(other.points_.rend()));
        }
    }
    else
    {
        if(atEnd)
        {
            vigra_assert(this->points_.back() == other.points_.front(),
                         "prepending non-matching reversed edge");
            this->points_.insert(
                this->points_.end(),
                skipFirst(other.points_.begin()), other.points_.end());
        }
        else
        {
            vigra_assert(this->points_.front() == other.points_.back(),
                         "prepending non-matching reversed edge");
            this->points_.insert(
                this->points_.begin(),
                other.points_.begin(), skipLast(other.points_.end()));
        }
    }

    if(lengthValid_)
        this->length_ += other.length();
    if(partialAreaValid_)
    {
        if(reverse)
            this->partialArea_ -= other.partialArea();
        else
            this->partialArea_ += other.partialArea();
    }

    if(boundingBoxValid_)
        this->boundingBox_ |= other.boundingBox();

    if(scanLines_.get())
    {
        if(reverse) // FIXME: call other.scanLines() (might be NULL)
            other.scanLines_->reverse();
        (*scanLines_) += other.scanLines();
    }
}

/********************************************************************/

DartPointIter::DartPointIter(GeoMap::Dart const &dart)
: edge_(dart.guaranteedEdge())
{
    if(dart.label() > 0)
    {
        index_ = 0;
        inc_ = 1;
        end_ = dart.size();
    }
    else
    {
        index_ = dart.size() - 1;
        inc_ = -1;
        end_ = -1;
    }
}

/********************************************************************/

GeoMap::GeoMap(vigra::Size2D imageSize)
: sigmaMappingArray_(101, 0),
  sigmaInverseMappingArray_(101, 0),
  sigmaMapping_(sigmaMappingArray_.begin() + 50),
  sigmaInverseMapping_(sigmaInverseMappingArray_.begin() + 50),
  nodeCount_(0),
  edgeCount_(0),
  faceCount_(0),
  imageSize_(imageSize),
  labelImage_(NULL),
  edgesSorted_(false)
{
    edges_.push_back(NULL_PTR(Edge));
}

GeoMap::GeoMap(const GeoMap &other)
: sigmaMappingArray_(other.sigmaMappingArray_.size(), 0),
  sigmaInverseMappingArray_(other.sigmaInverseMappingArray_.size(), 0),
  sigmaMapping_(sigmaMappingArray_.begin() + sigmaMappingArray_.size()/2),
  sigmaInverseMapping_(sigmaInverseMappingArray_.begin()
                       + sigmaInverseMappingArray_.size()/2),
  nodeCount_(other.nodeCount_),
  edgeCount_(other.edgeCount_),
  faceCount_(other.faceCount_),
  imageSize_(other.imageSize()),
  labelImage_(NULL),
  edgesSorted_(false)
{
    nodes_.resize(other.nodes_.size(), NULL_PTR(GeoMap::Node));
    for(ConstNodeIterator it = other.nodesBegin(); it.inRange(); ++it)
    {
        nodes_[(*it)->label()] =
            GeoMap::NodePtr(new GeoMap::Node(this, **it));
    }

    edges_.resize(other.edges_.size(), NULL_PTR(GeoMap::Edge));
    for(ConstEdgeIterator it = other.edgesBegin(); it.inRange(); ++it)
    {
        edges_[(*it)->label()] =
            GeoMap::EdgePtr(new GeoMap::Edge(this, **it));
    }

    // slightly more efficient than calling setSigmaMapping():
    std::copy(other.sigmaMappingArray_.begin(),
              other.sigmaMappingArray_.end(),
              sigmaMappingArray_.begin());
    std::copy(other.sigmaInverseMappingArray_.begin(),
              other.sigmaInverseMappingArray_.end(),
              sigmaInverseMappingArray_.begin());
    edgesSorted_ = other.edgesSorted_;

    if(other.faces_.size())
    {
        faces_.resize(other.faces_.size(), NULL_PTR(GeoMap::Face));
        for(ConstFaceIterator it = other.facesBegin(); it.inRange(); ++it)
        {
            faces_[(*it)->label()] =
                GeoMap::FacePtr(new GeoMap::Face(this, **it));
        }

        if(other.hasLabelImage())
        {
            labelImage_ = new LabelImage(*other.labelImage_);
            faceLabelLUT_ = other.faceLabelLUT_;
        }
    }
}

GeoMap::~GeoMap()
{
    // make sure the cells don't access this map anymore!
    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
        (*it)->uninitialize();
}

double angleTheta(double dy, double dx); // implemented in polygon.cxx

GeoMap::FacePtr GeoMap::faceAt(const Vector2 &position)
{
    vigra_precondition(mapInitialized(),
        "faceAt() called on graph (mapInitialized() == false)!");

    if(labelImage_)
    {
        GeoMap::LabelImage::difference_type p(detail::intVPos(position));
        if(labelImage_->isInside(p))
        {
            int faceLabel = (*labelImage_)[p];
            if(faceLabel > 0)
                return face(faceLabelLUT_[faceLabel]);
        }
    }

    for(FaceIterator it = finiteFacesBegin(); it.inRange(); ++it)
        if((*it)->contains(position))
            return *it;

    return face(0);
}

GeoMap::ConstFacePtr GeoMap::faceAt(const Vector2 &position) const
{
    vigra_precondition(mapInitialized(),
        "faceAt() called on graph (mapInitialized() == false)!");

    if(labelImage_)
    {
        GeoMap::LabelImage::difference_type p(detail::intVPos(position));
        if(labelImage_->isInside(p))
        {
            int faceLabel = (*labelImage_)[p];
            if(faceLabel > 0)
                return face(faceLabelLUT_[faceLabel]);
        }
    }

    for(ConstFaceIterator it = finiteFacesBegin(); it.inRange(); ++it)
        if((*it)->contains(position))
            return *it;

    return face(0);
}

GeoMap::NodePtr GeoMap::addNode(
    const Vector2 &position)
{
    GeoMap::Node *result = new GeoMap::Node(this, position);
    return node(result->label());
}

GeoMap::NodePtr GeoMap::addNode(
    const Vector2 &position, CellLabel label)
{
    if(label > nodes_.size())
        nodes_.resize(label, NULL_PTR(GeoMap::Node));
    GeoMap::Node *result = new GeoMap::Node(this, position);
    return node(result->label());
}

GeoMap::EdgePtr GeoMap::addEdge(
    const GeoMap::SigmaAnchor &startNeighbor,
    const GeoMap::SigmaAnchor &endNeighbor,
    const Vector2Array &points, CellLabel label)
{
    vigra_precondition(!mapInitialized(),
        "addEdge() called after initializeMap() (updating faces not implemented yet)");

    if(label > edges_.size())
        edges_.resize(label, NULL_PTR(GeoMap::Edge));
    GeoMap::Edge *result = new GeoMap::Edge(
        this, startNeighbor.nodeLabel(),  endNeighbor.nodeLabel(), points);

    if(startNeighbor.isSingular())
    {
        insertSigmaPredecessor(result->startNode()->anchor_, (int)result->label());
        result->startNode()->anchor_ = (int)result->label();
    }
    else
        insertSigmaPredecessor(startNeighbor.dartLabel(), (int)result->label());

    if((endNeighbor == startNeighbor) && result->partialArea() < 0)
    {
        insertSigmaPredecessor((int)result->label(), -(int)result->label());
    }
    else if(endNeighbor.isSingular())
    {
        insertSigmaPredecessor(result->endNode()->anchor_, -(int)result->label());
        result->endNode()->anchor_ = -(int)result->label();
    }
    else
        insertSigmaPredecessor(endNeighbor.dartLabel(), -(int)result->label());

    return edge(result->label());
}

GeoMap::FacePtr GeoMap::removeEdge(GeoMap::Dart &dart)
{
    if(mapInitialized())
    {
        if(dart.edge()->isBridge())
            return removeBridge(dart);
        else
            return mergeFaces(dart);
    }
    else
    {
        // this is just a graph -> detach from nodes & uninitialize()

        detachDart( dart.label());
        detachDart(-dart.label());

        dart.edge()->uninitialize();
    }
    return NULL_PTR(GeoMap::Face);
}

void GeoMap::sortEdgesDirectly()
{
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    typedef std::pair<double, int> DartAngle;

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if((*it)->isIsolated())
            continue;

        std::vector<DartAngle> dartAngles;

        GeoMap::Dart anchor((*it)->anchor()), d(anchor);
        do
        {
            vigra_precondition(
                d.size() >= 2, "cannot measure angle of darts with < 2 points!");
            dartAngles.push_back(
                DartAngle(angleTheta(-d[1][1] + d[0][1],
                                      d[1][0] - d[0][0]),
                          d.label()));
        }
        while(d.nextSigma() != anchor);

        std::sort(dartAngles.begin(), dartAngles.end());

        int predecessor = dartAngles[0].second;
        for(unsigned int i = 1; i < dartAngles.size(); ++i)
        {
            if(dartAngles[i-1].first == dartAngles[i].first)
            {
                std::stringstream s;
                s << "sortEdgesDirectly: edges leave node " << (*it)->label() << " at identical angles!";
                vigra_precondition(false, s.str());
            }
            sigmaMapping_[predecessor] = dartAngles[i].second;
            sigmaInverseMapping_[dartAngles[i].second] = predecessor;
            predecessor = dartAngles[i].second;
        }
        sigmaMapping_[predecessor] = dartAngles[0].second;
        sigmaInverseMapping_[dartAngles[0].second] = predecessor;
    }

    edgesSorted_ = true;
}

typedef std::vector<detail::DartPositionAngle> DartPositionAngles;
typedef DartPositionAngles::iterator DPAI;

template<class Iterator>
void rotateArray(Iterator begin, Iterator newBegin, Iterator end)
{
    typedef std::vector<typename Iterator::value_type> TempArray;
    TempArray temp(begin, end);
    typename Iterator::difference_type pos(newBegin - begin);
    std::copy(temp.begin() + pos, temp.end(), begin);
    std::copy(temp.begin(), temp.begin() + pos, begin + (end - newBegin));
}

inline double normAngle(double diff)
{
    if(diff < -M_PI)
        diff += 2*M_PI;
    if(diff >= M_PI)
        diff -= 2*M_PI;
    return diff;
}

void sortEdgesInternal(const Vector2 &currentPos,
                       double referenceAngle,
                       DPAI dpBegin, DPAI dpEnd,
                       double stepDist2, double minAngle,
                       GeoMap::UnsortableGroups &unsortable,
                       detail::PlannedSplits *splitInfo,
                       bool parallel)
{
    if(dpEnd - dpBegin < 2)
        return;

    bool unsortableState = true;
    for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
    {
        if(!dpi->dp.atEnd())
        {
            unsortableState = false;
            dpi->dp.intersectCircle(currentPos, stepDist2);
            //dpi->dp.leaveCircle(currentPos, stepDist2);
        }

        dpi->absAngle =
            std::atan2(-dpi->dp()[1] + currentPos[1],
                        dpi->dp()[0] - currentPos[0]);

        dpi->angle = normAngle(dpi->absAngle - referenceAngle);
    }

    if(unsortableState)
    {
        GeoMap::UnsortableGroups::value_type unsortableGroup;
        for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
            unsortableGroup.push_back(dpi->dp.dartLabel());
        unsortable.push_back(unsortableGroup);
        return;
    }

    std::sort(dpBegin, dpEnd);

    // handle cyclicity of array first (by rotation if necessary):
    DPAI firstGroupStart = dpEnd;
    if((--firstGroupStart)->angle + minAngle - 2*M_PI > dpBegin->angle)
    {
        // first and last dart are less than minAngle apart,..
        DPAI prev(firstGroupStart);
        while(--firstGroupStart > dpBegin)
        {
            // ..determine last decision point..
            if(prev->angle - firstGroupStart->angle >= minAngle)
            {
                // ..and rotate array to have whole dart group together
                // (start with "prev" == last correct firstGroupStart as newBegin)
                rotateArray(dpBegin, prev, dpEnd);
                break;
            }
            prev = firstGroupStart;
        }
    }

    bool storedSplitPos = false;
    detail::PlannedSplits::size_type storedSplitsOffset = 0;

    // look for groups of parallel edges
    DPAI groupStart = dpBegin,
          groupLast = groupStart, // for convenience; this is always groupEnd - 1
           groupEnd = groupLast + 1;
    for(; true; ++groupLast, ++groupEnd)
    {
        // group ending?
        if((groupEnd == dpEnd) || // last group
           (groupEnd->angle >= groupLast->angle + minAngle)) // decision here
        {
            if(splitInfo && groupEnd != dpEnd && parallel && !storedSplitPos)
            {
                // build group of edges to be split:
                // (since we have been given group of parallel darts and could
                // decide about the order at the current position, but not in
                // the last common position)
                unsigned int splitGroup = splitInfo->splitGroupCount++;
                storedSplitsOffset = splitInfo->size();
                for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi)
                    splitInfo->push_back(dpi->splitPos(splitGroup));
                storedSplitPos = true;
            }

            // recursion needed if > one dart in group:
            if(groupLast != groupStart)
            {
                // determine mean position of dart positions in subgroup:
                Vector2 meanPos(0, 0);
                for(DPAI dpi = groupStart; dpi != groupEnd; ++dpi)
                    meanPos += dpi->commonPos.set(dpi->dp);
                meanPos /= (groupEnd - groupStart);

                // sort parallel subgroup recursively:
                sortEdgesInternal(meanPos, normAngle(
                                      groupStart->absAngle +
                                      normAngle(groupLast->absAngle -
                                                groupStart->absAngle) / 2),
                                  groupStart, groupEnd,
                                  stepDist2, minAngle,
                                  unsortable, splitInfo,
                                  true);
            }

            if(groupEnd == dpEnd)
                break; // loop end

            groupStart = groupEnd;
        }
    }

    if(storedSplitPos)
    {
        // Later, in order to get the merging of the split nodes
        // right, we need to know the sigma order we just found out
        // here.  We store it already here, because it becomes much more
        // complicated after the splitting.

        detail::PlannedSplits::iterator
            storedSplitsBegin(splitInfo->begin() + storedSplitsOffset),
            storedSplitsEnd(storedSplitsBegin + (dpEnd - dpBegin));

        // Unfortunately, this is O(n^2) - I thought long about this,
        // but I cannot see any shortcut.  Maybe I am too tired, but
        // right now I have to get it working.
        int sigmaPos = 0;
        for(DPAI dpi = dpBegin; dpi != dpEnd; ++dpi, ++sigmaPos)
        {
            int dartLabel = dpi->dp.dartLabel();
            for(detail::PlannedSplits::iterator splIt = // wordplay ;-)
                    storedSplitsBegin; splIt != storedSplitsEnd; ++splIt)
            {
                if(splIt->dartLabel == dartLabel)
                {
                    splIt->sigmaPos = sigmaPos;
                    break;
                }
            }
        }
    }
}

void GeoMap::sortEdgesEventually(double stepDist, double minDist,
                                 UnsortableGroups &unsortable,
                                 bool splitEdges)
{
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    double minAngle = std::atan2(minDist, stepDist),
          stepDist2 = vigra::sq(stepDist);

    if(splitEdges)
        splitInfo_ = std::auto_ptr<detail::PlannedSplits>(
            new detail::PlannedSplits());

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        if((*it)->isIsolated())
            continue;

        DartPositionAngles dartPositions;

        GeoMap::Dart anchor((*it)->anchor()), d(anchor);
        do
        {
            dartPositions.push_back(detail::DartPositionAngle(d));
        }
        while(d.nextSigma() != anchor);

        sortEdgesInternal((*it)->position(), 0.0,
                          dartPositions.begin(), dartPositions.end(),
                          stepDist2, minAngle,
                          unsortable, splitInfo_.get(), false);

        int predecessor = dartPositions[0].dp.dartLabel();
        for(unsigned int i = 1; i < dartPositions.size(); ++i)
        {
            int successor = dartPositions[i].dp.dartLabel();
            sigmaMapping_[predecessor] = successor;
            sigmaInverseMapping_[successor] = predecessor;
            predecessor = successor;
        }
        sigmaMapping_[predecessor] = dartPositions[0].dp.dartLabel();
        sigmaInverseMapping_[dartPositions[0].dp.dartLabel()] = predecessor;
    }

    edgesSorted_ = true;
}

/*
 * Represents information needed during merging of parallel darts.
 * The dart is represented by a dartLabel and a "turnLater" flag,
 * since darts are constantly removed during the whole process.  The
 * "turnLater" ensures that we get a valid dart.
 */
struct MergeDart
{
    int dartLabel, sigmaPos;
    bool turnLater;

    MergeDart(int dl, bool tl, int sp)
    : dartLabel(dl), sigmaPos(sp), turnLater(tl)
    {}

        // default constructor - mark as uninitialized through zero dartLabel
    MergeDart()
    : dartLabel(0)
    {}

    bool operator<(MergeDart const &other) const
    {
        return !dartLabel || (other.dartLabel &&
                              sigmaPos < other.sigmaPos);
    }
};

void GeoMap::splitParallelEdges()
{
    vigra_precondition(splitInfo_.get(), "splitParallelEdges(): no planned splits (set splitEdges parameter of sortEdgesEventually?)");

    // find split group positions within flat array [O(N_splits)]
    // (this is used for filling the mergeDarts during splitting)
    std::vector<detail::PlannedSplits::difference_type> groupPositions;
    for(detail::PlannedSplits::iterator it = splitInfo_->begin();
        it != splitInfo_->end(); ++it)
    {
        if(it->splitGroup == groupPositions.size())
            groupPositions.push_back(it - splitInfo_->begin());
    }

    std::sort(splitInfo_->begin(), splitInfo_->end());

    bool hasPreferences = edgePreferences_.get();
    if(hasPreferences)
    {
        vigra_invariant(edgePreferences_->size() == edges_.size(),
                        "edge preferences given, but not exactly one per edge");
        edgePreferences_->resize(edges_.size() + splitInfo_->size());
    }
    else
    {
        edgePreferences_ = std::auto_ptr<EdgePreferences>(
            new EdgePreferences(edges_.size() + splitInfo_->size()));
    }
    EdgePreferences *edgePreferences = edgePreferences_.get();

    // define radius of circle for curvature detection (used if !hasPreferences):
    const double checkSurvivorDist  = 0.5;
    const double checkSurvivorDist2 = checkSurvivorDist*checkSurvivorDist;

    typedef std::vector<MergeDart> MergeDarts;
    MergeDarts mergeDarts(splitInfo_->size());
    for(detail::PlannedSplits::iterator it = splitInfo_->begin();
        it != splitInfo_->end(); ++it)
    {
        CellLabel newEdgeLabel =
            splitEdge(*edge(abs(it->dartLabel)),
                      it->segmentIndex, it->position)->label();

        if(hasPreferences)
        {
            (*edgePreferences)[newEdgeLabel] =
                (*edgePreferences)[abs(it->dartLabel)];
        }
        else
        {
            GeoMap::Dart d(dart(newEdgeLabel));
            Vector2 nodePos(d.startNode()->position());

            // intersect checkSurvivorDist-circle with dart
            DartPosition dp1(d);
            dp1.leaveCircle(nodePos, checkSurvivorDist2);
            d.nextSigma();
            DartPosition dp2(d);
            dp2.leaveCircle(nodePos, checkSurvivorDist2);

            // determine vectors between split node pos. & intersections..
            Vector2
                v1(dp1() - nodePos),
                v2(nodePos - dp2());

            // ..and choose dart with smallest enclosed angle:
            (*edgePreferences)[newEdgeLabel] =
                dot(v1, v2)/(v1.magnitude()*v2.magnitude());
        }

        detail::PlannedSplits::difference_type &pos(
            groupPositions[it->splitGroup]);

        mergeDarts[pos] =
            MergeDart((int)newEdgeLabel, it->dartLabel > 0, it->sigmaPos);

        ++pos;
    }

    vigra_invariant(edgePreferences->size() == edges_.size(),
                    "edge preferences should be exactly one per edge");

    splitInfo_.reset(); // splitting finished, free memory

    // for each split group, merge the resulting nodes:
    MergeDarts::iterator mergeDartsGroupEnd = mergeDarts.begin();
    for(unsigned int i = 0; i < groupPositions.size(); ++i)
    {
        MergeDarts::iterator mergeDartsGroupBegin = mergeDartsGroupEnd;
        mergeDartsGroupEnd = mergeDarts.begin() + (int)groupPositions[i];

        std::sort(mergeDartsGroupBegin, mergeDartsGroupEnd);

        // skip deleted edges (e.g. unsortable)
        while(mergeDartsGroupBegin->dartLabel == 0)
            ++mergeDartsGroupBegin;

        // now search for the best continuation to choose the survivor
        double bestContinuationValue = 0.0;
        int bestContinuationIndex = 0;

        for(MergeDarts::iterator it = mergeDartsGroupBegin;
            it != mergeDartsGroupEnd; ++it)
        {
            double pref = (*edgePreferences)[abs(it->dartLabel)];
            if(pref > bestContinuationValue)
            {
                bestContinuationValue = pref;
                bestContinuationIndex = it - mergeDartsGroupBegin;
            }
        }

        Dart survivor(
            dart(mergeDartsGroupBegin[bestContinuationIndex].dartLabel));
        if(mergeDartsGroupBegin[bestContinuationIndex].turnLater)
            survivor.nextSigma();

        GeoMap::Node &survivingNode(*survivor.startNode());
        int sigmaCenterLabel(survivor.clone().nextSigma().label()),
            sigmaNeighborLabel(sigmaCenterLabel);

        if(bestContinuationIndex)
        {
            // loop from survivor-1 to "begin" of sigma orbit
            MergeDarts::iterator it = mergeDartsGroupBegin + bestContinuationIndex;
            do
            {
                --it;

                // determine mergeDart and relocateDart
                GeoMap::Dart mergeDart(dart(it->dartLabel));
                if(it->turnLater)
                    mergeDart.nextSigma();

                vigra_invariant(mergeDart.startNode()->hasDegree(2),
                                "merge nodes are expected to have degree 2");

                GeoMap::Dart relocateDart(mergeDart);
                relocateDart.nextSigma();

                // checkConsistency(); // no modification should've happened so far

                // re-attach relocateDart to surviving node
                GeoMap::EdgePtr relocateEdge(relocateDart.edge());
                if(relocateDart.label() < 0)
                {
                    relocateEdge->endNodeLabel_ = survivingNode.label();
                    (*relocateEdge)[relocateEdge->size() - 1] =
                        survivingNode.position();
                }
                else
                {
                    relocateEdge->startNodeLabel_ = survivingNode.label();
                    (*relocateEdge)[0] = survivingNode.position();
                }

                detachDart(relocateDart.label());
                insertSigmaPredecessor(sigmaNeighborLabel, relocateDart.label());
                sigmaNeighborLabel = relocateDart.label();

                // propagate flags set on any of the merged edges to the
                // surviving, merged edge
                survivor.edge()->setFlag(mergeDart.edge()->flags());

                // remove mergeDart and its startNode:
                GeoMap::Node &mergedNode(*mergeDart.startNode());
                removeEdge(mergeDart);
                removeIsolatedNode(mergedNode);
            }
            while(it != mergeDartsGroupBegin);
        }

        sigmaNeighborLabel = sigmaMapping_[sigmaCenterLabel];

        for(MergeDarts::iterator it =
                mergeDartsGroupBegin + bestContinuationIndex + 1;
            it != mergeDartsGroupEnd; ++it)
        {
            // determine mergeDart and relocateDart
            GeoMap::Dart mergeDart(dart(it->dartLabel));
            if(it->turnLater)
                mergeDart.nextSigma();

            vigra_invariant(mergeDart.startNode()->hasDegree(2),
                            "merge nodes are expected to have degree 2");

            GeoMap::Dart relocateDart(mergeDart);
            relocateDart.nextSigma();

            // re-attach relocateDart to surviving node
            GeoMap::EdgePtr relocateEdge(relocateDart.edge());
            if(relocateDart.label() < 0)
            {
                relocateEdge->endNodeLabel_ = survivingNode.label();
                (*relocateEdge)[relocateEdge->size() - 1] =
                    survivingNode.position();
            }
            else
            {
                relocateEdge->startNodeLabel_ = survivingNode.label();
                (*relocateEdge)[0] = survivingNode.position();
            }
            detachDart(relocateDart.label());
            insertSigmaPredecessor(sigmaNeighborLabel, relocateDart.label());

            // propagate flags set on any of the merged edges to the
            // surviving, merged edge
            survivor.edge()->setFlag(mergeDart.edge()->flags());

            // remove mergeDart and its startNode:
            GeoMap::Node &mergedNode(*mergeDart.startNode());
            removeEdge(mergeDart);
            removeIsolatedNode(mergedNode);
        }
    }

    edgePreferences_.reset();
}

void GeoMap::setSigmaMapping(SigmaMapping const &sigmaMapping, bool sorted)
{
    vigra_precondition(!mapInitialized(),
        "sigma orbits cannot be changed after initializeMap()");

    vigra_precondition(sigmaMapping.size() >= (2*edges_.size() - 1),
                       "setSigmaMapping: sigmaMapping too small!");
    SigmaMapping::const_iterator sigma(
        sigmaMapping.begin() + (sigmaMapping.size() / 2));

    if(sigmaMappingArray_.size() < 2*edges_.size() - 1)
        resizeSigmaMapping(2*edges_.size() - 1);

    std::copy(sigma - (edges_.size() - 1), sigma + edges_.size(),
              sigmaMapping_ - (edges_.size() - 1));
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        int label = (int)(*it)->label();
        sigmaInverseMapping_[sigmaMapping_[ label]] =  label;
        sigmaInverseMapping_[sigmaMapping_[-label]] = -label;
    }

    edgesSorted_ = sorted;
}

void GeoMap::initializeMap(bool initLabelImage)
{
    vigra_precondition(!mapInitialized(),
                       "initializeMap() called more than once");
    if(!edgesSorted())
        sortEdgesDirectly();

    initContours();
    //std::cerr << faceCount_ << " contours found, embedding...\n";
    embedFaces(initLabelImage);
}

typedef vigra::MultiArray<2, int> LabelImage;

void markEdgeInLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage);

void GeoMap::setHasLabelImage(bool onoff)
{
    if(onoff == hasLabelImage())
        return;

    if(onoff)
    {
        vigra_precondition(imageSize_.area() > 0,
                           "initLabelImage: imageSize must be non-zero!");
        labelImage_ = new LabelImage(
            LabelImage::size_type(imageSize().width(), imageSize().height()), 0);
        faceLabelLUT_.initIdentity(faces_.size());

        for(FaceIterator it = finiteFacesBegin(); it.inRange(); ++it)
        {
            std::auto_ptr<vigra::Scanlines> scanlines =
                (*it)->scanLines();
            fillScannedPoly(*scanlines, (int)(*it)->label(),
                            destMultiArrayRange(*labelImage_));
            (*it)->pixelArea_ = 0;
        }

        for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
            markEdgeInLabelImage((*it)->scanLines(),
                                 *labelImage_);

        // determine pixelArea_:
        for(GeoMap::LabelImage::traverser lrow = labelImage_->traverser_begin();
            lrow != labelImage_->traverser_end(); ++lrow)
        {
            for(GeoMap::LabelImage::traverser::next_type lit = lrow.begin();
                lit != lrow.end(); ++lit)
            {
                int label = *lit;
                if(label >= 0)
                    ++face(label)->pixelArea_;
            }
        }
    }
    else
    {
        delete labelImage_;
        labelImage_ = NULL;
    }
}

void GeoMap::initContours()
{
    new Face(this, Dart(this, 0)); // create infinite face, dart will be ignored

    // fill list of faces with contours, i.e. no face will have a
    // contour after this, but all faces will carry the properties
    // (area, bbox, ...) of the phi orbit:
    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        if((*it)->leftFaceLabel_ == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart( (int)(*it)->label()));
        if((*it)->rightFaceLabel_ == UNINITIALIZED_CELL_LABEL)
            new Face(this, dart(-(int)(*it)->label()));
    }
}

struct AbsAreaCompare
{
        // FIXME: actually, const pointers would suffice:
    bool operator()(GeoMap::FacePtr f1, GeoMap::FacePtr f2) const
    {
        double a1 = f1->area(), a2 = f2->area();
        double absdiff = std::fabs(a1) - std::fabs(a2);
        if(std::fabs(absdiff) < 1e-2 && ((a1 < 0) != (a2 < 0)))
            return (a1 < 0); // for faces with equal area, prefer the exterior one
        return absdiff > 0; // else, prefer face with larger absolute area
    }
};

void GeoMap::embedFaces(bool initLabelImage)
{
    // the result of this function is to transform the result of
    // initContours (i.e. preliminary faces, which are just phi
    // orbits) into the final faces; many "faces" will disappear (if
    // they represent hole contours) and be embedded as holes into
    // their surrounding faces.
    //
    // second, the label image will be set up

    vigra_precondition(!labelImage_,
        "embedFaces() called with already-initialized labelImage");

    if(initLabelImage)
    {
        vigra_precondition(imageSize_.area() > 0,
                           "initLabelImage: imageSize must be non-zero!");
        labelImage_ = new LabelImage(
            LabelImage::size_type(imageSize().width(), imageSize().height()), 0);
        faceLabelLUT_.initIdentity(faces_.size());
    }

    // copy and remove all preliminary contours except the infinite one:
    GeoMap::Faces contours(faces_.begin() + 1, faces_.end());
    std::sort(contours.begin(), contours.end(), AbsAreaCompare());
    std::fill(faces_.begin() + 1, faces_.end(), NULL_PTR(Face));

    for(unsigned int i = 0; i < contours.size(); ++i)
    {
        GeoMap::Face &contour(*contours[i]); // FIXME: const

        GeoMap::Dart anchor(contour.contour());

        bool isExterior = contour.area() <= 0;

        if(!isExterior)
        {
            faces_[contour.label()] = contours[i];

            if(initLabelImage)
            {
                std::auto_ptr<vigra::Scanlines> scanlines =
                    contour.scanLines();
                contour.pixelArea_ =
                    fillScannedPoly(*scanlines, (int)contour.label(),
                                    destMultiArrayRange(*labelImage_));
                // no need for rawAddEdgeToLabelImage here, since we
                // work with darts anyways, and there's no easy way to
                // ensure that the negative counts will not be wrong
                // (esp. also for interior bridges etc.)
                drawScannedPoly(*scanlines, -1,
                                destMultiArrayRange(*labelImage_));
            }
        }
        else
        {
            // contour is a hole, determine parent face
            GeoMap::FacePtr parent = NULL_PTR(GeoMap::Face);

            if(initLabelImage)
            {
                ContourPointIter cpi(anchor);
                while(cpi.inRange())
                {
                    GeoMap::LabelImage::difference_type p(detail::intVPos(*cpi++));
                    if(labelImage_->isInside(p))
                    {
                        int parentLabel = (*labelImage_)[p];
                        if(parentLabel >= 0)
                        {
                            parent = face(parentLabel);
                            break;
                        }
                    }
                }
            }

            if(!parent)
            {
                ContourPointIter cpi(anchor);
                while(cpi.inRange())
                {
                    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
                    {
                        if((*it)->contains(*cpi++))
                        {
                            parent = *it;
                            goto parent_found; // double break
                        }
                    }
                }
            }

            if(!parent)
            {
                parent = face(0);
//                 vigra_postcondition(
//                     parent->contains(anchor[0]),
//                     "contour could not be embedded (parent not found)");
            }

        parent_found:
//             std::cerr << "  embedding contour " << contour.label()
//                       << " in face " << parent->label() << "\n";
            parent->embedContour(anchor);
            contour.uninitialize();
        }
    }

    if(initLabelImage)
    {
        // in the case of holes, the face's pixelArea()s will be wrong:
        for(FaceIterator it = facesBegin(); it.inRange(); ++it)
            (*it)->pixelArea_ = 0;

        // interior bridges may not have been set to negative labels yet:
        for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
            if((*it)->isBridge())
                drawScannedPoly((*it)->scanLines(), -1,
                                destMultiArrayRange(*labelImage_));

        // remove temporary edge markings and fix pixelAreas:
        for(GeoMap::LabelImage::traverser lrow = labelImage_->traverser_begin();
            lrow != labelImage_->traverser_end(); ++lrow)
        {
            for(GeoMap::LabelImage::traverser::next_type lit = lrow.begin();
                lit != lrow.end(); ++lit)
            {
                int label = *lit;
                if(label < 0)
                    *lit = 0;
                else
                    ++face(label)->pixelArea_;
            }
        }

        // redo all edge markings correctly:
        for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
            markEdgeInLabelImage((*it)->scanLines(),
                                 *labelImage_);
    }
}

class LookupNewLabel
{
    const std::vector<CellLabel> &newLabels_;

  public:
    LookupNewLabel(const std::vector<CellLabel> &newLabels)
    : newLabels_(newLabels)
    {}

    int operator()(int label) const
    {
        if(label >= 0)
            return (int)newLabels_[label];
        return label;
    }
};

void GeoMap::changeFaceLabels(
    const std::vector<CellLabel> &newFaceLabels,
    CellLabel maxFaceLabel)
{
    vigra_precondition(mapInitialized(),
        "changeFaceLabels() called on graph (no faces available)");
    vigra_precondition(newFaceLabels.size() == faces_.size(),
        "changeFaceLabels(): 1-to-1 mapping expected (wrong newFaceLabels size)");

    GeoMap::Faces newFaces(maxFaceLabel);
    for(CellLabel l = 0; l < newFaceLabels.size(); ++l)
    {
        if(!faces_[l])
            continue;
        newFaces[newFaceLabels[l]] = faces_[l];
        faces_[l]->label_ = newFaceLabels[l];
    }
    std::swap(faces_, newFaces);

    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        (*it)->leftFaceLabel_ = newFaceLabels[(*it)->leftFaceLabel_];
        (*it)->rightFaceLabel_ = newFaceLabels[(*it)->rightFaceLabel_];
    }

    if(hasLabelImage())
    {
        transformMultiArray(srcMultiArrayRange(*labelImage_, labelAccessor()),
                            destMultiArray(*labelImage_),
                            LookupNewLabel(newFaceLabels));
        faceLabelLUT_.initIdentity(faces_.size());
    }
}

void GeoMap::resizeSigmaMapping(SigmaMapping::size_type newSize)
{
    SigmaMapping
        newSigma(newSize, 0),
        newSigmaInverse(newSize, 0);
    std::copy(sigmaMappingArray_.begin(), sigmaMappingArray_.end(),
              newSigma.begin() + sigmaMappingArray_.size()/2);
    std::copy(sigmaInverseMappingArray_.begin(),
              sigmaInverseMappingArray_.end(),
              newSigmaInverse.begin() + sigmaMappingArray_.size()/2);

    std::swap(sigmaMappingArray_, newSigma);
    sigmaMapping_ = sigmaMappingArray_.begin() + newSize/2;

    std::swap(sigmaInverseMappingArray_, newSigmaInverse);
    sigmaInverseMapping_ = sigmaInverseMappingArray_.begin() + newSize/2;
}

void GeoMap::insertSigmaPredecessor(int successor, int newPredecessor)
{
    while(sigmaMappingArray_.size() < 2*(unsigned)abs(newPredecessor)+1)
        resizeSigmaMapping(2*sigmaMappingArray_.size()-1);

    if(!successor)
    {
        sigmaMapping_[newPredecessor] = newPredecessor;
        sigmaInverseMapping_[newPredecessor] = newPredecessor;
        return;
    }

    int oldPredecessor = sigmaInverseMapping_[successor];
    sigmaMapping_[oldPredecessor] = newPredecessor;
    sigmaMapping_[newPredecessor] = successor;
    sigmaInverseMapping_[successor] = newPredecessor;
    sigmaInverseMapping_[newPredecessor] = oldPredecessor;
}

void GeoMap::detachDart(int dartLabel)
{
    GeoMap::Node &node(*dart(dartLabel).startNode());

    int successor = sigmaMapping_[dartLabel];
    int predecessor = sigmaInverseMapping_[dartLabel];

    sigmaMapping_[predecessor] = successor;
    sigmaInverseMapping_[successor] = predecessor;

    if(node.anchor_ == dartLabel)
    {
        if(successor != dartLabel)
            node.anchor_ = successor;
        else
            node.anchor_ = 0;
    }
}

GeoMap::NodePtr GeoMap::nearestNode(
    const Vector2 &position,
    double maxSquaredDist)
{
    NodeMap::iterator n(
        nodeMap_.nearest(PositionedNodeLabel(position, 0), maxSquaredDist));
    if(n != nodeMap_.end())
        return node(n->second.payload);
    return NULL_PTR(GeoMap::Node);
}

bool GeoMap::checkConsistency()
{
    //std::cerr << "GeoMap[" << this << "].checkConsistency()\n";
    bool result = true;

    if((sigmaMapping_ !=
        sigmaMappingArray_.begin() + (sigmaMappingArray_.size()-1)/2) ||
       (sigmaInverseMapping_ !=
        sigmaInverseMappingArray_.begin() + (sigmaMappingArray_.size()-1)/2) ||
       !(sigmaMappingArray_.size() & 1) ||
       (sigmaMappingArray_.size() != sigmaInverseMappingArray_.size()))
    {
        std::cerr << "  Sigma mapping arrays not correctly setup!\n";
        std::cerr << "      sigma: size " << sigmaMappingArray_.size()
                  << ", center @" << (sigmaMapping_ - sigmaMappingArray_.begin())
                  << "\n";
        std::cerr << "    inverse: size " << sigmaInverseMappingArray_.size()
                  << ", center @" << (
                      sigmaInverseMapping_ - sigmaInverseMappingArray_.begin())
                  << "\n";
        result = false;
    }
    else
    {
        if(sigmaMappingArray_.size() < maxEdgeLabel()*2-1)
        {
            std::cerr << "  Sigma mapping arrays not large enough!\n";
            result = false;
        }

        int dist = sigmaMappingArray_.size()/2;
        for(int label = -dist; label <= dist; ++label)
        {
#if 0
            // this test could probably be gotten rid of:
            if(abs(sigmaMapping_[label]) >= dist ||
               abs(sigmaInverseMapping_[label]) >= dist)
            {
                std::cerr << "  Sigma mapping arrays contain junk (sigma["
                          << label << "] = " << sigmaMapping_[label]
                          << ", sigma^-1[" << label << "] = "
                          << sigmaInverseMapping_[label] << ")!\n";
                result = false;
                continue;
            }
#endif
            if(sigmaMapping_[label] &&
               (CellLabel)abs(label) < maxEdgeLabel() &&
               edges_[abs(label)].get() &&
               sigmaInverseMapping_[sigmaMapping_[label]] != label)
            {
                std::cerr << "  Sigma inverse is not correct ("
                          << label << " -> " << sigmaMapping_[label]
                          << " -> " << sigmaInverseMapping_[sigmaMapping_[label]]
                          << ")!\n";
                result = false;
            }
        }
        if(!result)
            return result; // function may not terminate in this state
    }

    unsigned int actualNodeCount = 0, actualEdgeCount = 0, actualFaceCount = 0,
        connectedComponentsCount = 0, isolatedNodesCount = 0, contourCount = 0;

    for(NodeIterator it = nodesBegin(); it.inRange(); ++it)
    {
        ++actualNodeCount;

        if((*it)->map() != this)
        {
            std::cerr << "  Node " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }

        if((*it)->isIsolated())
        {
            ++isolatedNodesCount;
            continue;
        }

        Dart anchor((*it)->anchor()), dart(anchor);
        do
        {
            if(!dart.edge().get())
            {
                std::cerr << "  Node " << (*it)->label()
                          << " has broken sigma orbit: dart "
                          << dart.label() << " does not exist!\n";

                std::cerr << "    orbit so far: "
                          << anchor.label() << " (anchor)";
                Dart d(anchor);
                while(d != dart)
                {
                    d.nextSigma();
                    std::cerr << " -> " << d.label();
                }
                std::cerr << "\n";

                result = false;
                break;
            }

            if(dart.startNodeLabel() != (*it)->label())
            {
                std::cerr << "  Node " << (*it)->label()
                          << " has broken sigma orbit:\n"
                          << "  contains Dart " << dart.label()
                          << " from node " << dart.startNodeLabel()
                          << " -> " << dart.endNodeLabel() << "\n";
                result = false;
                break;
            }
        }
        while(dart.nextSigma() != anchor);
    }
    if(actualNodeCount != nodeCount())
    {
        std::cerr << "  Node count wrong (" << nodeCount()
                  << ", should be " << actualNodeCount << ")!\n";
            result = false;
    }

    for(EdgeIterator it = edgesBegin(); it.inRange(); ++it)
    {
        ++actualEdgeCount;

        if((*it)->map() != this)
        {
            std::cerr << "  Edge " << (*it)->label() << " has wrong map()!\n";
            result = false;
            break;
        }

        if(mapInitialized())
        {
            if((*it)->isBridge() && (*it)->isLoop())
            {
                std::cerr << "  Edge " << (*it)->label()
                          << " is both loop and bridge!?\n";
                result = false;
            }

            if(!(*it)->leftFace() || !(*it)->rightFace())
            {
                std::cerr << "  Edge " << (*it)->label() << " has invalid faces:\n";
                if(!(*it)->leftFace())
                    std::cerr << "     left face (" << (*it)->leftFaceLabel()
                              << ") does not exist!\n";
                if(!(*it)->rightFace())
                    std::cerr << "    right face (" << (*it)->rightFaceLabel()
                              << ") does not exist!\n";
                result = false;
                continue; // no need to check geometry
            }
        }

        if((*it)->size() < 2)
        {
            std::cerr << "  Edge " << (*it)->label() << " is too short ("
                      << (*it)->size() << " points)!\n";
            result = false;
        }
        else
        {
            if((**it)[0] != (*it)->startNode()->position() ||
                (**it)[(*it)->size()-1] != (*it)->endNode()->position())
            {
                std::cerr << "  Edge " << (*it)->label()
                          << " has non-matching end positions:\n"

                          << "    start node " << (*it)->startNodeLabel() << " at "
                          << (*it)->startNode()->position() << " is "
                          << ((**it)[0] - (*it)->startNode()->position()).magnitude()
                          << " pixels from " << (**it)[0] << "\n"

                          << "    end node " << (*it)->endNodeLabel() << " at "
                          << (*it)->endNode()->position() << " is "
                          << ((**it)[(*it)->size()-1]
                              - (*it)->endNode()->position()).magnitude()
                          << " pixels from " << (**it)[(*it)->size()-1] << "\n";
                result = false;
            }
        }

        Polygon poly((*it)->begin(), (*it)->end());
        if(fabs((*it)->partialArea() - poly.partialArea()) > 1e-4)
        {
            std::cerr << "  Edge " << (*it)->label()
                      << " has partialArea " << (*it)->partialArea()
                      << " instead of " << poly.partialArea() << "\n";
            result = false;
        }
        if(fabs((*it)->length() - poly.length()) > 1e-4)
        {
            std::cerr << "  Edge " << (*it)->label()
                      << " has length " << (*it)->length()
                      << " instead of " << poly.length() << "\n";
            result = false;
        }
    }
    if(actualEdgeCount != edgeCount())
    {
        std::cerr << "  Edge count wrong (" << edgeCount()
                  << ", should be " << actualEdgeCount << ")!\n";
            result = false;
    }

    if(!mapInitialized())
        return result; // cannot test Faces yet (do not exist)

    for(FaceIterator it = facesBegin(); it.inRange(); ++it)
    {
        ++actualFaceCount;
        connectedComponentsCount += (*it)->holeCount();

        // simulate (*it)->contourCount()
        contourCount += (*it)->holeCount();
        if((*it)->label())
            ++contourCount;

        GeoMap::Face &face(**it);
        if(face.map() != this)
        {
            std::cerr << "  Face " << face.label() << " has wrong map()!\n";
            result = false;
            break;
        }

        typedef std::map<int, unsigned int> SeenAnchors;
        SeenAnchors seenAnchors;

        bool outerContour = face.label() > 0;
        unsigned int contourIndex = 0;
        for(GeoMap::Face::ContourIterator ci = face.contoursBegin();
            ci != face.contoursEnd();
            ++ci, ++contourIndex, outerContour = false)
        {
            double area = 0.0;
            int canonicalAnchor(ci->label());

            if(!ci->edge())
            {
                std::cerr << "  Face " << face.label() << " contains an invalid contour anchor:\n    dart " << ci->label()
                          << " belongs to removed edge!\n";
                result = false;
                continue;
            }

            Dart anchor(*ci), dart(anchor);
            do
            {
                if(dart.leftFaceLabel() != face.label())
                {
                    std::cerr << "  Dart " << dart.label()
                              << " has leftFaceLabel() " << dart.leftFaceLabel()
                              << " != " << face.label() << "!\n";
                    result = false;
                    break;
                }

                if(!dart.edge()->isBridge())
                    area += dart.partialArea();

                if(ci->label() < canonicalAnchor)
                    canonicalAnchor = ci->label();
            }
            while(dart.nextPhi() != anchor);

            if(dart != anchor)
                continue; // loop canceled, area not valid

            if((area <= 0) == outerContour)
            {
                std::cerr << "  Face " << face.label() << " contains an "
                          << (outerContour
                              ? "outer anchor with negative area"
                              : "inner anchor with positive area")
                          << ":\n  contour from dart " << ci->label()
                          << " has area " << area << "!\n";
                result = false;
            }

            std::pair<SeenAnchors::iterator, bool>
                seenAnchor(seenAnchors.insert(
                               std::make_pair(canonicalAnchor,
                                              contourIndex)));
            if(!seenAnchor.second)
            {
                std::cerr << "  Face " << face.label()
                          << " contains duplicate anchors at contour indices "
                          << contourIndex << " and "
                          << seenAnchor.first->second << "!\n";
                result = false;
            }
        }
    }
    if(actualFaceCount != faceCount())
    {
        std::cerr << "  Face count wrong (" << faceCount()
                  << ", should be " << actualFaceCount << ")!\n";
            result = false;
    }

    // make the following checks and output more intuitive:
    actualNodeCount -= isolatedNodesCount;

    if(actualNodeCount - actualEdgeCount + actualFaceCount
       - connectedComponentsCount != 1)
    {
        std::cerr << "  Euler-Poincare invariant violated! (N - E + F - C = "
                  << actualNodeCount << " - " << actualEdgeCount << " + "
                  << actualFaceCount << " - " << connectedComponentsCount << " = "
                  << (int)(actualNodeCount - actualEdgeCount + actualFaceCount
                           - connectedComponentsCount)
                  << " != 1, map cannot be planar)\n";
        result = false;
    }

    if(actualNodeCount - actualEdgeCount + contourCount
       - 2*connectedComponentsCount != 0)
    {
        std::cerr << "  Euler-Poincare invariant violated! (N - E + B - 2*C = "
                  << actualNodeCount << " - " << actualEdgeCount << " + "
                  << contourCount << " - 2*" << connectedComponentsCount << " = "
                  << (int)(actualNodeCount - actualEdgeCount + contourCount
                           - 2*connectedComponentsCount)
                  << " != 0, map cannot be planar)\n";
        result = false;
    }

    return result;
}

/********************************************************************/
/*                                                                  */
/*                  Euler operations (and helpers)                  */
/*                                                                  */
/********************************************************************/

bool GeoMap::removeIsolatedNode(GeoMap::Node &node)
{
    if(!removeNodeHook(node))
        return false;

    node.uninitialize();
    return true;
}

void rawAddEdgeToLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage, int diff)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min((int)labelImage.size(1), scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size(0))
                end = labelImage.size(0);

            for(int x = begin; x < end; ++x)
                labelImage[LabelImage::difference_type(x, y)] += diff;
        }
    }
}

void markEdgeInLabelImage(
    const vigra::Scanlines &scanlines, LabelImage &labelImage)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min((int)labelImage.size(1), scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size(0))
                end = labelImage.size(0);

            for(int x = begin; x < end; ++x)
            {
                LabelImage::difference_type pos(x, y);
                int label = labelImage[pos];
                labelImage[pos] = (label >= 0 ? -1 : label-1);
            }
        }
    }
}

void removeEdgeFromLabelImage(
    const vigra::Scanlines &scanlines,
    LabelImage &labelImage,
    LabelImage::value_type substituteLabel,
    PixelList &outputPixels)
{
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min((int)labelImage.size(1), scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const vigra::Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.size(0))
                end = labelImage.size(0);

            for(int x = begin; x < end; ++x)
            {
                LabelImage::reference old(labelImage(x, y));
                if(old != -1)
                {
                    old += 1;
                }
                else
                {
                    old = substituteLabel;
                    outputPixels.push_back(vigra::Point2D(x, y));
                }
            }
        }
    }
}

GeoMap::EdgePtr GeoMap::mergeEdges(const GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "mergeEdges called on removed dart!");
    Dart d1(dart);
    d1.nextSigma();
    vigra_precondition(d1.edgeLabel() != dart.edgeLabel(),
                       "mergeEdges called on self-loop or startNode's degree = 1!");

    Dart d2(d1);
    d1.nextSigma();
    vigra_precondition(d1 == dart,
                       "mergeEdges cannot remove node with degree > 2!");

    vigra_assert((d1.leftFaceLabel() == d2.rightFaceLabel()) &&
                 (d2.leftFaceLabel() == d1.rightFaceLabel()),
                 "mergeEdges: broken map (left/rightFaceLabel)");

    if(d1.label() > 0 && d2.label() < 0)
        std::swap(d1, d2); // minimize number of reverse()s necessary

    GeoMap::Edge &survivor(*d1.edge());
    GeoMap::Node &mergedNode(*d2.startNode());
    GeoMap::Edge &mergedEdge(*d2.edge());

    if(mapInitialized())
    {
        GeoMap::Face *faces[2];
        faces[0] = &(*d2.leftFace());
        faces[1] = &(*d2.rightFace());
        for(GeoMap::Face **faceIt = faces; faceIt != faces+2; ++faceIt)
        {
            GeoMap::Face::Contours::iterator cEnd = (*faceIt)->anchors_.end();
            for(GeoMap::Face::Contours::iterator it = (*faceIt)->anchors_.begin();
                it != cEnd; ++it)
            {
                if(it->edgeLabel() == d2.edgeLabel())
                {
                    it->nextPhi();
                    break;
                }
            }
        }
    }

    if(!preMergeEdgesHook(d1))
        return NULL_PTR(GeoMap::Edge);
    if(!removeNodeHook(mergedNode))
        return NULL_PTR(GeoMap::Edge);

    if(labelImage_)
    {
        rawAddEdgeToLabelImage(mergedEdge.scanLines(), *labelImage_, 1);
        rawAddEdgeToLabelImage(survivor.scanLines(), *labelImage_, 1);
    }

    if(survivor.startNodeLabel() != mergedNode.label())
    {
        survivor.concatenate(
            mergedEdge, true, // atEnd = true, append mergedEdge (good)
            mergedEdge.startNodeLabel() != mergedNode.label());

        survivor.endNodeLabel_ = d2.endNodeLabel();
    }
    else
    {
        survivor.concatenate(
            mergedEdge, false, // atEnd = false, prepend mergedEdge
            mergedEdge.startNodeLabel() == mergedNode.label());

        survivor.startNodeLabel_ = d2.endNodeLabel();
    }

    if(labelImage_)
    {
        rawAddEdgeToLabelImage(survivor.scanLines(), *labelImage_, -1);
    }

    // replace -d2 with d1 within orbits / anchor:
    int successor = sigmaMapping_[-d2.label()];
    if(successor == -d2.label())
    {
        // no other dart at that end node
        sigmaMapping_[d1.label()] = d1.label();
        sigmaInverseMapping_[d1.label()] = d1.label();
    }
    else
    {
        int predecessor = sigmaInverseMapping_[-d2.label()];
        sigmaMapping_[predecessor] = d1.label();
        sigmaMapping_[d1.label()] = successor;
        sigmaInverseMapping_[successor] = d1.label();
        sigmaInverseMapping_[d1.label()] = predecessor;
    }
    d2.endNode()->anchor_ = d1.label();

    mergedNode.uninitialize();
    mergedEdge.uninitialize();

    postMergeEdgesHook(survivor);

    return this->edge(survivor.label());
}

GeoMap::EdgePtr GeoMap::splitEdge(
    GeoMap::Edge &edge, unsigned int segmentIndex)
{
    return splitEdge(edge, segmentIndex, edge[0], false);
}

GeoMap::EdgePtr GeoMap::splitEdge(
    GeoMap::Edge &edge, unsigned int segmentIndex,
    const Vector2 &newPoint, bool insertPoint)
{
    vigra_precondition(segmentIndex < edge.size() - 1,
                       "splitEdge: invalid segmentIndex");

    preSplitEdgeHook(edge, segmentIndex, newPoint, insertPoint);

    GeoMap::Node
        &newNode(*addNode(insertPoint ? newPoint : edge[segmentIndex])),
        &changedNode(*edge.endNode());

    if(labelImage_)
        rawAddEdgeToLabelImage(edge.scanLines(), *labelImage_, 1);

    if(insertPoint)
    {
        ++segmentIndex;
        edge.insert(edge.begin() + segmentIndex, newPoint);
    }

    GeoMap::Edge *result = new GeoMap::Edge(
        this, newNode.label(), changedNode.label(), edge.split(segmentIndex));
    result->leftFaceLabel_ = edge.leftFaceLabel_;
    result->rightFaceLabel_ = edge.rightFaceLabel_;
    result->flags_ = edge.flags_;

    if(sigmaMappingArray_.size() < 2*result->label()+1)
        resizeSigmaMapping(2*sigmaMappingArray_.size()-1);

    int successor = sigmaMapping_[-(int)edge.label()];
    int predecessor = sigmaInverseMapping_[-(int)edge.label()];

    // setup newNode (fixed degree 2)
    sigmaMapping_[-(int)edge.label()] = (int)result->label();
    sigmaMapping_[(int)result->label()] = -(int)edge.label();
    sigmaInverseMapping_[(int)result->label()] = -(int)edge.label();
    sigmaInverseMapping_[-(int)edge.label()] = (int)result->label();
    newNode.anchor_ = (int)result->label();

    // edge now ends in newNode
    edge.endNodeLabel_ = newNode.label();

    if(successor == -(int)edge.label())
    {
        // changedNode has degree 1, replace -edgeLabel with -resultLabel
        sigmaMapping_[-(int)result->label()] = -(int)result->label();
        sigmaInverseMapping_[-(int)result->label()] = -(int)result->label();
        changedNode.anchor_ = -(int)result->label();
    }
    else
    {
        sigmaMapping_[predecessor] = -(int)result->label();
        sigmaMapping_[-(int)result->label()] = successor;
        sigmaInverseMapping_[successor] = -(int)result->label();
        sigmaInverseMapping_[-(int)result->label()] = predecessor;
        if(changedNode.anchor_ == -(int)edge.label())
            changedNode.anchor_ = -(int)result->label();
    }

    if(labelImage_)
    {
        edge.scanLines_.reset();
        rawAddEdgeToLabelImage(edge.scanLines(), *labelImage_, -1);
        rawAddEdgeToLabelImage(result->scanLines(), *labelImage_, -1);
    }

    postSplitEdgeHook(edge, *result);

    return this->edge(result->label());
}

void GeoMap::associatePixels(GeoMap::Face &face, const PixelList &pixels)
{
    face.pixelArea_ += pixels.size();
//     for(unsigned int i = 0; i < pixels.size(); ++i)
//         pixelBounds_ |= pixels[i];
    associatePixelsHook(face, pixels);
}

GeoMap::FacePtr GeoMap::removeBridge(const GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "removeBridge called on removed dart!");

    GeoMap::Edge &edge(*dart.edge());
    GeoMap::Face &face(*dart.leftFace());
    vigra_precondition(face.label() == dart.rightFace()->label(),
                       "removeBridge needs a bridge dart!");
    GeoMap::Node &node1(*dart.startNode());
    GeoMap::Node &node2(*dart.endNode());
    vigra_precondition(node1.label() != node2.label(),
                       "Inconsistent map: bridge to be removed is also a self-loop!?");

    // COMPLEXITY: depends on callbacks (preRemoveBridgeHook)
    if(!preRemoveBridgeHook(dart))
        return NULL_PTR(GeoMap::Face);

    Dart newAnchor1(dart), newAnchor2(dart);
    newAnchor1.prevSigma();
    newAnchor2.nextAlpha().prevSigma();
    // COMPLEXITY: depends on number of contours in face
    Face::AnchorIterator contourPos = face.findComponentAnchor(dart);

    // remove both darts from both sigma orbits:
    detachDart( dart.label());
    detachDart(-dart.label());

    // COMPLEXITY: depends on face anchors.erase/push_back, may be O(contour count)
    if(newAnchor1.edgeLabel() == dart.edgeLabel())
    {
        removeIsolatedNode(*newAnchor1.startNode());
        if(newAnchor2.edgeLabel() == dart.edgeLabel())
        {
            removeIsolatedNode(*newAnchor2.startNode());
            face.anchors_.erase(contourPos);
        }
        else
            *contourPos = newAnchor2;
    }
    else if(newAnchor2.edgeLabel() == dart.edgeLabel())
    {
        removeIsolatedNode(*newAnchor2.startNode());
        *contourPos = newAnchor1;
    }
    else
    {
        if(contourPos == face.contoursBegin() and face.label())
        {
            // determine outer anchor, swap if necessary:
            // COMPLEXITY: depends on number of darts in contours
            if(contourArea(newAnchor1) < contourArea(newAnchor2))
                std::swap(newAnchor1, newAnchor2);
        }
        *contourPos = newAnchor1;
        face.anchors_.push_back(newAnchor2);
    }

    // COMPLEXITY: depends on number of pixel facets crossed by the bridge
    PixelList associatedPixels;
    if(labelImage_)
        removeEdgeFromLabelImage(
            edge.scanLines(), *labelImage_, face.label(), associatedPixels);

    edge.uninitialize();

    // COMPLEXITY: depends on callbacks (postRemoveBridgeHook)
    postRemoveBridgeHook(face);

    // COMPLEXITY: depends on callbacks (associatePixelsHook)
    if(associatedPixels.size())
        associatePixels(face, associatedPixels);

    return this->face(face.label());
}

GeoMap::FacePtr GeoMap::mergeFaces(const GeoMap::Dart &dart)
{
    vigra_precondition(dart.edge(),
                       "mergeFaces called on removed dart!");

    GeoMap::Dart removedDart(dart);

    if(dart.leftFace()->area() < dart.rightFace()->area())
        removedDart.nextAlpha();
    if(!removedDart.rightFaceLabel()) // face 0 shall stay face 0
        removedDart.nextAlpha();

    GeoMap::Edge &mergedEdge(*removedDart.edge());
    GeoMap::Face &survivor(*removedDart.leftFace());
    GeoMap::Face &mergedFace(*removedDart.rightFace());
    GeoMap::Node &node1(*removedDart.startNode());

    vigra_precondition(survivor.label() != mergedFace.label(),
                       "mergeFaces(): dart belongs to a bridge!");

    // COMPLEXITY: depends on callbacks (preMergeFacesHook)
    if(!preMergeFacesHook(removedDart))
        return NULL_PTR(GeoMap::Face);

    // remember bounding box of merged face for later updating
    GeoMap::Face::BoundingBox mergedBBox;
    if(survivor.flag(GeoMap::Face::BOUNDING_BOX_VALID))
        mergedBBox = mergedFace.boundingBox();

    // relabel contour's leftFaceLabel
    // COMPLEXITY: depends on number of darts in mergedFace's contours
    for(Face::ContourIterator it = mergedFace.contoursBegin();
        it != mergedFace.contoursEnd(); ++it)
    {
        GeoMap::Dart d(*it);
        while(d.nextPhi().leftFaceLabel() != survivor.label())
            d.internalLeftFaceLabel() = survivor.label();
    }

    // COMPLEXITY: depends on number of contours in F1 + F2
    Face::AnchorIterator contour1 = survivor.findComponentAnchor(removedDart);
    Face::AnchorIterator contour2 = mergedFace.findComponentAnchor(
        GeoMap::Dart(removedDart).nextAlpha());

    // re-use an old anchor for the merged contour
    if(contour1->edgeLabel() == mergedEdge.label())
    {
        contour1->nextPhi();
        if(contour1->edgeLabel() == mergedEdge.label())
        {
            *contour1 = *contour2;
            if(contour1->edgeLabel() == mergedEdge.label())
            {
                contour1->nextPhi();
                if(contour1->edgeLabel() == mergedEdge.label())
                {
                    vigra_assert(node1.label() == removedDart.endNodeLabel(),
                                 "special-case: merging a self-loop");
                    // results in an isolated node:
                    survivor.anchors_.erase(contour1);
                }
            }
        }
    }

    // append all remaining anchors to survivor's list:
    if(mergedFace.anchors_.size() > 1)
    {
#ifndef ANCHOR_LISTS // toggle according to anchors_ being std::vector/std::list
        survivor.anchors_.insert(survivor.anchors_.end(),
                                 mergedFace.anchors_.begin(), contour2);
        survivor.anchors_.insert(survivor.anchors_.end(),
                                 ++contour2, mergedFace.anchors_.end());
#else
        mergedFace.anchors_.erase(contour2);
        survivor.anchors_.splice(survivor.anchors_.end(), mergedFace.anchors_);
#endif
    }

    // relabel region in image
    PixelList associatedPixels;
    if(labelImage_)
    {
//         relabelImage(map.labelImage.subImage(mergedFace.pixelBounds_),
//                      mergedFace.label(), survivor.label())
        // COMPLEXITY: depends on number of faces prev. merged into mergedFace
        faceLabelLUT_.relabel(mergedFace.label(), survivor.label());

        // COMPLEXITY: depends on number of pixel facets crossed by mergedEdge
        removeEdgeFromLabelImage(
            mergedEdge.scanLines(),
            *labelImage_, survivor.label(), associatedPixels);

//         survivor.pixelBounds_ |= mergedFace.pixelBounds_;
    }

    // remove both darts from both sigma orbits:
    detachDart( removedDart.label());
    detachDart(-removedDart.label());

    // remove singular nodes
    if(node1.isIsolated())
    {
        vigra_assert(removedDart.endNodeLabel() == node1.label(),
                     "mergeFaces can only create isolated nodes from self-loops");
        removeIsolatedNode(node1);
    }

    if(survivor.flag(GeoMap::Face::AREA_VALID))
        survivor.area_ += mergedFace.area();
    survivor.pixelArea_ += mergedFace.pixelArea_;

    if(survivor.flag(GeoMap::Face::BOUNDING_BOX_VALID))
        survivor.boundingBox_ |= mergedBBox;

    mergedEdge.uninitialize();
    mergedFace.uninitialize();

    // COMPLEXITY: depends on callbacks (postMergeFacesHook)
    postMergeFacesHook(survivor);

    // COMPLEXITY: depends on callbacks (associatePixelsHook)
    if(associatedPixels.size())
        associatePixels(survivor, associatedPixels);

    return this->face(survivor.label());
}

/********************************************************************/

void GeoMap::Node::setPosition(const Vector2 &p)
{
    vigra_precondition(initialized(), "setPosition() of uninitialized node!");
    map_->nodeMap_.erase(
        map_->nodeMap_.nearest(PositionedNodeLabel(position_, label_),
                               vigra::NumericTraits<double>::epsilon()));
    position_ = p;

    if(!isIsolated())
    {
        GeoMap::Dart d(map_, anchor_);
        do
        {
            if(d.label() > 0)
            {
                (*map_->edge(d.label()))[ 0] = p;
            }
            else
            {
                GeoMap::Edge &edge(*map_->edge(-d.label()));
                edge[edge.size()-1] = p;
            }
        }
        while(d.nextSigma().label() != anchor_);
    }

    map_->nodeMap_.insert(PositionedNodeLabel(p, label_));
}

GeoMap::Face::AnchorIterator
GeoMap::Face::findComponentAnchor(const GeoMap::Dart &dart)
{
    for(AnchorIterator it = anchors_.begin(); it != anchors_.end(); ++it)
        if(*it == dart)
            return it;

    for(AnchorIterator it = anchors_.begin(); it != anchors_.end(); ++it)
    {
        GeoMap::Dart d(*it);
        while(d.nextPhi() != *it)
            if(d == dart)
                return it;
    }

    vigra_fail("findComponentAnchor failed: dart not found in face contours!");
    return anchors_.begin(); // never reached
}

void GeoMap::Face::embedContour(const Dart &anchor)
{
    anchors_.push_back(anchor);

    Dart dart(anchor);
    for(; dart.leftFaceLabel() != label_; dart.nextPhi())
        dart.internalLeftFaceLabel() = label_;

    // don't try to do this on the fly (special bridge handling needed):
    if(flag(AREA_VALID))
        area_ += contourArea(dart);

    vigra_postcondition(dart == anchor,
                        "contour labeled partially?!");
}
