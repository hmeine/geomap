#ifndef CELLSTATISTICS_HXX
#define CELLSTATISTICS_HXX

#include "foureightsegmentation.hxx"
#include "statisticfunctor.hxx"
#include "mydebug.hxx"

#include <vigra/stdimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/tinyvector.hxx>

#include <algorithm>
#include <iostream>
#include <list>
#include <vector>

typedef vigra::FImage OriginalImage;
typedef vigra::FImage GradientImage;

typedef vigra::TinyVector<float, 2> Float2D;
typedef vigra::BasicImage<Float2D> EdgeDirImage;
typedef vigra::BasicImage<bool> EdgelImage;

struct SegmentationData
{
    OriginalImage preparedOriginal_;
    OriginalImage smoothedOriginal_;
    GradientImage gradientMagnitude_;
    EdgeDirImage  edgeDirection_;
    EdgeDirImage  edgeDirGradient_;

    bool       doEdgeRethinning;
    bool       doRemoveDegree2Nodes;
    bool       doMergeCompleteRegions;
    EdgelImage edgelImage;
    bool       doAutoProtection;
    double     autoProtectionThreshold;

    SegmentationData()
    {
        doEdgeRethinning = true;
        doRemoveDegree2Nodes = true;
        doMergeCompleteRegions = true;
        doAutoProtection = false;
        autoProtectionThreshold = 0.3;
    }
};

namespace std { void swap(SegmentationData &a, SegmentationData &b); }

/********************************************************************/

enum { EPBorder = 1, EPCanny = 2, EPManual = 4 };

class EdgeProtection
{
  public:
    typedef unsigned short FlagType;

  protected:
    std::vector<FlagType> edgeProtection_;

  public:
    EdgeProtection(unsigned int size)
    : edgeProtection_(size, 0)
    {}

    bool protectEdge(vigra::cellimage::CellLabel edgeLabel,
                     FlagType flag, bool protect = true)
    {
        if((bool)(edgeProtection_[edgeLabel] & flag) != protect)
        {
            edgeProtection_[edgeLabel] ^= flag;
            return true;
        }
        return false;
    }

    bool edgeProtected(vigra::cellimage::CellLabel edgeLabel) const
    {
        return (bool)edgeProtection_[edgeLabel];
    }

    bool sameProtection(vigra::cellimage::CellLabel edgeLabel1,
                        vigra::cellimage::CellLabel edgeLabel2) const
    {
        return edgeProtection_[edgeLabel1] == edgeProtection_[edgeLabel2];
    }

    bool edgeProtected(vigra::cellimage::CellLabel edgeLabel,
                       FlagType flag) const
    {
        return (bool)(edgeProtection_[edgeLabel] & flag);
    }
};

/********************************************************************/

struct CountEdgels
{
    unsigned int count;

    CountEdgels(): count(0) {}
    void operator()(bool e) { if(e) ++count; }
    unsigned int operator()() const { return count; }
};

/********************************************************************/

struct CellStatistics
{
    typedef vigra::cellimage::GeoMap
        Segmentation;
    typedef StatisticFunctor<OriginalImage::PixelType>
        //typedef vigra::FindAverage<OriginalImage::PixelType>
        FaceStatistics;
    typedef vigra::FindAverage<GradientImage::PixelType>
        EdgeStatistics;

        // members storing statistics
    std::vector<FaceStatistics>
        faceStatistics_;
    std::vector<EdgeStatistics>
        edgeStatistics_;
    std::vector<vigra::cellimage::CellLabel>
        mergedEdges_;
    std::vector<unsigned int>
        edgels_; // number of Canny edgels on an edges' pixels
    std::vector<vigra::TinyVector<float, 2> >
        nodeCenters_;

        // members storing source data
    SegmentationData *segmentationData;
    EdgeProtection   *edgeProtection;

	vigra::Rect2D     segDataBounds;

    static std::vector<Float2D> configurationDirections_;

    CellStatistics(const Segmentation &initialSegmentation,
                   SegmentationData *segmentationData,
                   EdgeProtection *edgeProtection);

  protected:
    mutable vigra::Rect2D lastChanges_;

        // temporary data for storage between preXXX- and postXXX-calls
    FaceStatistics              tempFaceStatistics_;
    EdgeStatistics              tempEdgeStatistics_;
    vigra::cellimage::CellLabel node1Label_, node2Label_;
    vigra::cellimage::CellLabel edge1Label_, edge2Label_;
        // FIXME: vector would be beneficial, since the ScanIterator
        // already delivers sorted values! just push_back
    std::list<vigra::Point2D> nodePoints_;

    inline void preRemoveEdge(const Segmentation::DartTraverser &dart);
    inline void postRemoveEdge(Segmentation::FaceInfo &face);

  public:
    inline void preRemoveIsolatedNode(const Segmentation::DartTraverser &);
    inline void postRemoveIsolatedNode(Segmentation::FaceInfo &face);

    inline void preRemoveBridge(const Segmentation::DartTraverser &dart);
    inline void postRemoveBridge(Segmentation::FaceInfo &face);

    inline void preMergeFaces(const Segmentation::DartTraverser &dart);
    inline void postMergeFaces(Segmentation::FaceInfo &face);

    inline void preMergeEdges(const Segmentation::DartTraverser &dart);
    inline void postMergeEdges(Segmentation::EdgeInfo &edge);

    bool protectEdge(vigra::cellimage::CellLabel edgeLabel,
                     EdgeProtection::FlagType flag, bool protect = true) const
    {
        bool result =
            edgeProtection->protectEdge(edgeLabel, flag, protect);
//         if(result && protect)
//         {
//             edgeProtection->lastValidLevel =
//                 std::min(edgeProtection->lastValidLevel, edgeRemoveLevels_[edgeLabel]);
//         }
        if(mergedEdges_[edgeLabel] != edgeLabel)
            result =
                protectEdge(mergedEdges_[edgeLabel], flag, protect) || result;
        return result;
    }

    bool edgeProtected(vigra::cellimage::CellLabel edgeLabel) const
    {
        if(edgeProtection->edgeProtected(edgeLabel))
            return true;
        else if(mergedEdges_[edgeLabel] == edgeLabel)
            return false;
        else
            return edgeProtected(mergedEdges_[edgeLabel]);
    }

    const vigra::Rect2D initAutoProtection(const Segmentation &seg) const
    {
        vigra::Rect2D result;

        for(Segmentation::ConstEdgeIterator it = seg.edgesBegin();
            it.inRange(); ++it)
        {
            CountEdgels edgelCount;
            inspectCell(seg.edgeScanIterator(
                            it->label, segmentationData->edgelImage.upperLeft()),
                        edgelCount);
            if(protectEdge(it->label, EPCanny,
                           (edgelCount() > it->size *
                            segmentationData->autoProtectionThreshold)))
                result |= it->bounds;
        }

        return result;
    }

    const vigra::Rect2D & lastChanges() const
    {
        return lastChanges_;
    }

    vigra::Rect2D clearLastChanges() const
    {
        vigra::Rect2D result(lastChanges_);
        lastChanges_.setLowerRight(lastChanges_.upperLeft());
        return result;
    }

    CellStatistics &operator=(const CellStatistics &other);
};

/********************************************************************/

void nodeRethinning(vigra::cellimage::GeoMap &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel nodeLabel);

void edgeRethinning(vigra::cellimage::GeoMap &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel edgeLabel,
					const vigra::Rect2D &rethinRange);

/********************************************************************/

inline void CellStatistics::preRemoveIsolatedNode(
    const Segmentation::DartTraverser &)
{}

inline void CellStatistics::postRemoveIsolatedNode(
    Segmentation::FaceInfo &face)
{
    lastChanges_ |= face.bounds;
}

inline void CellStatistics::preRemoveEdge(const Segmentation::DartTraverser &dart)
{
    tempFaceStatistics_ = faceStatistics_[dart.leftFaceLabel()];
    inspectCell(dart.segmentation()->edgeScanIterator
                (dart.edgeLabel(), segmentationData->preparedOriginal_.upperLeft()),
                tempFaceStatistics_);

    node1Label_ = dart.edge().start.startNodeLabel();
    node2Label_ = dart.edge().end.startNodeLabel();
}

inline void CellStatistics::postRemoveEdge(Segmentation::FaceInfo &face)
{
    faceStatistics_[face.label]= tempFaceStatistics_;
    lastChanges_ |= face.bounds;

    /*nodeRethinning(*dart.segmentation(),
      segmentationData->gradientMagnitude_,  node1Label_);
      if(node1Label_ != node2Label_)
      nodeRethinning(*dart.segmentation(),
      segmentationData->gradientMagnitude_, node2Label_);*/
}

inline void CellStatistics::preRemoveBridge(
    const Segmentation::DartTraverser &dart)
{
    preRemoveEdge(dart);
}

inline void CellStatistics::postRemoveBridge(
    Segmentation::FaceInfo &face)
{
    postRemoveEdge(face);
}

inline void CellStatistics::preMergeFaces(const Segmentation::DartTraverser &dart)
{
    // initialize tempFaceStatistics_ with leftFace and edge:
    preRemoveEdge(dart); // also sets node1Label_ and node2Label_

    tempFaceStatistics_(faceStatistics_[dart.rightFaceLabel()]);
}

inline void CellStatistics::postMergeFaces(Segmentation::FaceInfo &face)
{
    postRemoveEdge(face);
}

inline void CellStatistics::preMergeEdges(const Segmentation::DartTraverser &dart)
{
    Segmentation::DartTraverser d(dart);
    d.nextSigma();
    edge1Label_ = dart.edgeLabel();
    edge2Label_ = d.edgeLabel();

    if(segmentationData->doEdgeRethinning)
    {
        Segmentation::ScanIterator<vigra::Point2D>::type
            nodePosScanner(dart.segmentation()->nodeScanIterator
                            (d.startNodeLabel(), vigra::Point2D()));
        nodePoints_.clear();
        for(; nodePosScanner.inRange(); ++nodePosScanner)
            nodePoints_.push_back(*nodePosScanner);
    }

    if(!dart.leftFaceLabel() || !dart.rightFaceLabel())
    {
        tempEdgeStatistics_.reset();
        tempEdgeStatistics_(
            vigra::NumericTraits<GradientImage::PixelType>::max());
        return;
    }
    tempEdgeStatistics_ = edgeStatistics_[edge1Label_];
    tempEdgeStatistics_(edgeStatistics_[edge2Label_]);
    inspectCell(dart.segmentation()->nodeScanIterator
                (d.startNodeLabel(), segmentationData->gradientMagnitude_.upperLeft()),
                tempEdgeStatistics_);
}

// basically, relabels like "*pixel = newCell;"
// additionally, updates bounding boxes of old and new pixel owners
inline void
reassociate(vigra::cellimage::GeoMap &geoMap,
            const vigra::cellimage::CellImage::traverser pixel,
            const vigra::cellimage::CellPixel &newCell)
{
    using namespace vigra::cellimage;

    CellPixel oldCell(*pixel);

    vigra::Point2D pos(pixel - geoMap.cells);

    // get CellInfo for current owner of the pixel
    GeoMap::CellInfo &oldCI(
        oldCell.type() == CellTypeLine ?
        (GeoMap::CellInfo &)geoMap.edge(oldCell.label()) : (
            oldCell.type() == CellTypeVertex ?
            (GeoMap::CellInfo &)geoMap.node(oldCell.label()) :
            (GeoMap::CellInfo &)geoMap.face(oldCell.label())));

//     std::cerr << "reassociate(" << pos << ", " << oldCell << "->"
//               << newCell << "), old bbox: " << oldCI.bounds << "\n";

    // relabel pixel
    *pixel = newCell;

    // now, check if the pixel was on one of the four sides of the
    // bounding box, and if so, scan that row/column and possibly
    // shrink the bounding box:

    if(pos.x == oldCI.bounds.left())
    {
        // check first column of the bounding box..
        if(std::find((geoMap.cells + oldCI.bounds.upperLeft()).columnIterator(),
                     (geoMap.cells + vigra::Diff2D(pos.x, oldCI.bounds.bottom())).columnIterator(),
                     oldCell)
           == (geoMap.cells + vigra::Diff2D(pos.x, oldCI.bounds.bottom())).columnIterator())
        {
            // removed the only pixel from that column -> change bbox
            oldCI.bounds.setUpperLeft(vigra::Point2D(oldCI.bounds.left() + 1,
                                                     oldCI.bounds.top()));
        }
    }

    if(pos.x == oldCI.bounds.right() - 1)
    {
        // check last column of the bounding box..
        if(std::find((geoMap.cells + vigra::Diff2D(pos.x, oldCI.bounds.top())).columnIterator(),
                     (geoMap.cells + vigra::Diff2D(pos.x, oldCI.bounds.bottom())).columnIterator(),
                     oldCell)
           == (geoMap.cells + vigra::Diff2D(pos.x, oldCI.bounds.bottom())).columnIterator())
        {
            // removed the only pixel from that column -> change bbox
            oldCI.bounds.setLowerRight(vigra::Point2D(oldCI.bounds.right() - 1,
                                                      oldCI.bounds.bottom()));
        }
    }

    if(pos.y == oldCI.bounds.top())
    {
        // check first row of the bounding box..
        if(std::find((geoMap.cells + oldCI.bounds.upperLeft()).rowIterator(),
                     (geoMap.cells + vigra::Diff2D(oldCI.bounds.right(), pos.y)).rowIterator(),
                     oldCell)
           == (geoMap.cells + vigra::Diff2D(oldCI.bounds.right(), pos.y)).rowIterator())
        {
            // removed the only pixel from that row -> change bbox
            oldCI.bounds.setUpperLeft(vigra::Point2D(oldCI.bounds.left(),
                                                     oldCI.bounds.top() + 1));
        }
    }

    if(pos.y == oldCI.bounds.bottom() - 1)
    {
        // check last row of the bounding box..
        if(std::find((geoMap.cells + vigra::Diff2D(oldCI.bounds.left(), pos.y)).rowIterator(),
                     (geoMap.cells + vigra::Diff2D(oldCI.bounds.right(), pos.y)).rowIterator(),
                     oldCell)
           == (geoMap.cells + vigra::Diff2D(oldCI.bounds.right(), pos.y)).rowIterator())
        {
            // removed the only pixel from that row -> change bbox
            oldCI.bounds.setLowerRight(vigra::Point2D(oldCI.bounds.right(),
                                                      oldCI.bounds.bottom() - 1));
        }
    }

    // add pixel pos to bounding box of new cell

    (newCell.type() == CellTypeLine ?
     (GeoMap::CellInfo &)geoMap.edge(newCell.label()) : (
         newCell.type() == CellTypeVertex ?
         (GeoMap::CellInfo &)geoMap.node(newCell.label()) :
         (GeoMap::CellInfo &)geoMap.face(newCell.label())))
                     .bounds |= pos;
}

inline void CellStatistics::postMergeEdges(Segmentation::EdgeInfo &edge)
{
    edgeStatistics_[edge.label] = tempEdgeStatistics_;
    lastChanges_ |= edge.bounds;

    // build tree of mergedEdges_
    vigra::cellimage::CellLabel targetLabel = edge.label;
    while(mergedEdges_[targetLabel] != targetLabel)
        targetLabel = mergedEdges_[targetLabel];
    mergedEdges_[targetLabel] =
        (edge.label == edge1Label_) ? edge2Label_ : edge1Label_;

    // to rethin, we need gradient/edgeness information for all edgels:
    if(segDataBounds.contains(edge.bounds))
    {
        if(segmentationData->doEdgeRethinning)
        {
            vigra::Rect2D rethinRange(edge.bounds);
            rethinRange.addBorder(1);
            rethinRange &= segDataBounds;

            edgeRethinning(*edge.start.segmentation(),
                           segmentationData->gradientMagnitude_, edge.label,
                           rethinRange);
            
            /* This loop checks all pixels that belonged to a node
             * previously and relocates them if their neighbor pixel
             * has a higher gradient.
             *
             * There are two cases, three without reflection:
             *    #o#  #o.  / #.. \   # - edge pixel
             *    .#.  .##  | o#. |   . - region pixel
             *    ...  ...  \ .#. /   o - possible new location
             *                            for center pixel
             *
             * The inner loop circulates over all four diagonal
             * directions and looks for the first edge pixel, then
             * checks for the three configurations by looking for the
             * second edge pixels at circulator offsets 2/3 for the
             * left configurations or -3 for the reflected. If the
             * configuration is found, and the position marked 'o' has
             * a higher gradient, the three new neighbor pixels are
             * checked to make sure that the topology would not be
             * changed. If all conditions hold, the center edge pixel
             * is relocated.
             */
            for(std::list<vigra::Point2D>::const_iterator it =
                    nodePoints_.begin(); it != nodePoints_.end(); ++it)
            {
                vigra::cellimage::CellImageEightCirculator
                    cellImageCirc(edge.start.segmentation()->cells + (*it),
                                  vigra::EightNeighborCode::NorthWest),
                    cellImageEnd(cellImageCirc);

                if(cellImageCirc.center()->type() != vigra::cellimage::CellTypeLine)
                    continue;

                GradientImage::traverser gradient(
                    segmentationData->gradientMagnitude_.upperLeft() + (*it));
                GradientImage::value_type gradientValue(*gradient);

                vigra::NeighborhoodCirculator<
                    GradientImage::traverser, vigra::EightNeighborCode>
                    gradCirc(gradient, vigra::EightNeighborCode::NorthWest);

                do
                {
                    if(*cellImageCirc == *cellImageCirc.center())
                    {
                        // one diagonal edge continuation found, check
                        // the above-mentioned configurations:
                        if(((cellImageCirc[2] == *cellImageCirc.center()) ||
                            (cellImageCirc[3] == *cellImageCirc.center()))
                           && (gradCirc[1] > gradientValue))
                        {
                            // one of the two first configurations found:
                            ++cellImageCirc;
                            
                            // check the three new neighbor pixels
                            vigra::cellimage::CellImageEightCirculator
                                newNeighborCirc(cellImageCirc);
                            newNeighborCirc.moveCenterToNeighbor();
                            if((*newNeighborCirc == *cellImageCirc) &&
                               (newNeighborCirc[1] == *cellImageCirc) &&
                               (newNeighborCirc[-1] == *cellImageCirc))
                            {
                                // "turn around" to get the region pixel
                                // with which to fill in the center:
                                vigra::cellimage::CellPixel
                                    regionPixel(cellImageCirc[-4]);

                                // everything's allright, re-associate
                                // the pixels:
                                reassociate(*edge.start.segmentation(),
                                            cellImageCirc.base(),
                                            *cellImageCirc.center());
                                reassociate(*edge.start.segmentation(),
                                            cellImageCirc.center(),
                                            regionPixel);

                                // mark the two edge continuation
                                // pixels as new candidates for
                                // edgeRethinning:
                                --cellImageCirc;
                                nodePoints_.push_back(*it + cellImageCirc.diff());
                                cellImageCirc += 2;
                                if(*newNeighborCirc == *cellImageCirc)
                                    ++cellImageCirc;
                                nodePoints_.push_back(*it + cellImageCirc.diff());
                                break;
                            } else {
                                // this is our loop variable! (and its inc'ed by 2)
                                --cellImageCirc;
                            }
                        } else if((cellImageCirc[-3] == *cellImageCirc.center())
                                  && (gradCirc[-11] > gradientValue))
                        {
                            // reflected configuration found, further
                            // comments see above (only offsets changed)
                            --cellImageCirc;
                                
                            vigra::cellimage::CellImageEightCirculator
                                newNeighborCirc(cellImageCirc);
                            newNeighborCirc.moveCenterToNeighbor();
                            if((*newNeighborCirc == *cellImageCirc) &&
                               (newNeighborCirc[1] == *cellImageCirc) &&
                               (newNeighborCirc[-1] == *cellImageCirc))
                            {
                                vigra::cellimage::CellPixel
                                    regionPixel(cellImageCirc[-4]);
                                reassociate(*edge.start.segmentation(),
                                            cellImageCirc.base(),
                                            *cellImageCirc.center());
                                reassociate(*edge.start.segmentation(),
                                            cellImageCirc.center(),
                                            regionPixel);

                                ++cellImageCirc;
                                nodePoints_.push_back(*it + cellImageCirc.diff());
                                cellImageCirc -= 3;
                                nodePoints_.push_back(*it + cellImageCirc.diff());
                                break;
                            } else {
                                ++cellImageCirc;
                            }
                        }
                    }

                    gradCirc += 2;
                    cellImageCirc += 2;
                }
                while(cellImageCirc != cellImageEnd);
            }
        }

        if(segmentationData->doAutoProtection)
        {
            CountEdgels edgelCount;
            inspectCell(edge.start.segmentation()->edgeScanIterator
                        (edge.label, segmentationData->edgelImage.upperLeft()),
                        edgelCount);
            edgels_[edge.label] = edgelCount();
            protectEdge(edge.label, EPCanny,
                        (edgelCount() > edge.size *
                         segmentationData->autoProtectionThreshold));
        }
    }
}

#endif // CELLSTATISTICS_HXX
