#ifndef CELLSTATISTICS_HXX
#define CELLSTATISTICS_HXX

#include "foureightsegmentation.hxx"
#include "statisticfunctor.hxx"

#include <vigra/stdimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/tinyvector.hxx>

#include <vector>
#include <iostream>

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
        // the following members would fit better in a separate
        // statistics object which is unique-per-pyramid instead of
        // per-level
    bool       doEdgeRethinning;
    bool       doRemoveDegree2Nodes;
    EdgelImage edgelImage;
    bool       doAutoProtection;
    double     autoProtectionThreshold;
};

namespace std { void swap(SegmentationData &a, SegmentationData &b); }

struct CountEdgels
{
    unsigned int count;

    CountEdgels(): count(0) {}
    void operator()(bool e) { if(e) ++count; }
    unsigned int operator()() const { return count; }
};

/********************************************************************/

enum { EPBorder = 1, EPCanny = 2, EPManual = 4 };

class EdgeProtection
{
  public:
    typedef unsigned short FlagType;

  protected:
    std::vector<FlagType> edgeProtection_;
    unsigned int          timestamp_;

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
            ++timestamp_;
            return true;
        }
        return false;
    }

    unsigned int timestamp() const
    {
        return timestamp_;
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

struct CellStatistics
{
    typedef vigra::cellimage::GeoMap
        Segmentation;
//typedef StatisticFunctor<OriginalImage::PixelType>
    typedef vigra::FindAverage<OriginalImage::PixelType>
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

    static std::vector<Float2D> configurationDirections_;

    CellStatistics(const Segmentation &initialSegmentation,
                   SegmentationData *segmentationData,
                   EdgeProtection *edgeProtection);

        // members storing source data
    SegmentationData *segmentationData;
    EdgeProtection   *edgeProtection;

	vigra::Rect2D     segDataBounds;

  protected:
    mutable vigra::Rect2D lastChanges_;

        // temporary data for storage between preXXX- and postXXX-calls
    FaceStatistics tempFaceStatistics_;
    EdgeStatistics tempEdgeStatistics_;
    vigra::cellimage::CellLabel node1Label_, node2Label_;
    vigra::cellimage::CellLabel edge1Label_, edge2Label_;

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
            edgeProtection->protectEdge(edgeLabel, flag, protect);;
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
{
}

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
