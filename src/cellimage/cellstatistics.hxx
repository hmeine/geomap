#ifndef CELLSTATISTICS_HXX
#define CELLSTATISTICS_HXX

#include <foureightsegmentation.hxx>
#include <statisticfunctor.hxx>

#include <vigra/stdimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/tinyvector.hxx>

#include <vector>
#include <iostream>

typedef vigra::FImage OriginalImage;
typedef vigra::FImage GradientImage;

typedef vigra::TinyVector<float, 2> Float2D;
typedef vigra::BasicImage<Float2D> EdgeDirImage;

struct SegmentationData
{
    OriginalImage preparedOriginal_;
    GradientImage gradientMagnitude_;
    EdgeDirImage  edgeDirection_;
    EdgeDirImage  edgeDirGradient_;
};

namespace std { void swap(SegmentationData &a, SegmentationData &b); }

void nodeRethinning(vigra::cellimage::FourEightSegmentation &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel nodeLabel);

void edgeRethinning(vigra::cellimage::FourEightSegmentation &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel edgeLabel);

extern bool doEdgeRethinning;

struct CellStatistics
{
    typedef vigra::cellimage::FourEightSegmentation
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
    std::vector<vigra::TinyVector<float, 2> >
        nodeCenters_;

        // constant information, could even be static?
    std::vector<Float2D>
        configurationDirections_;

    CellStatistics(const Segmentation &initialSegmentation,
                   SegmentationData  *segmentationData);

        // members storing source data
    SegmentationData *segmentationData_;

  protected:
    mutable vigra::Rect2D lastChanges_;

        // temporary data for storage between preXXX- and postXXX-calls
    FaceStatistics tempFaceStatistics_;
    EdgeStatistics tempEdgeStatistics_;
    vigra::cellimage::CellLabel node1Label_, node2Label_;
    vigra::cellimage::CellLabel edge1Label_, edge2Label_;

    void preRemoveEdge(const Segmentation::DartTraverser &dart)
    {
        tempFaceStatistics_ = faceStatistics_[dart.leftFaceLabel()];
        inspectCell(dart.segmentation()->edgeScanIterator
                    (dart.edgeLabel(), segmentationData_->preparedOriginal_.upperLeft()),
                    tempFaceStatistics_);

        node1Label_ = dart.edge().start.startNodeLabel();
        node2Label_ = dart.edge().end.startNodeLabel();
    }
    void postRemoveEdge(Segmentation::FaceInfo &face)
    {
        faceStatistics_[face.label]= tempFaceStatistics_;
        lastChanges_ |= face.bounds;

        /*nodeRethinning(*dart.segmentation(),
                       segmentationData_->gradientMagnitude_,  node1Label_);
        if(node1Label_ != node2Label_)
            nodeRethinning(*dart.segmentation(),
            segmentationData_->gradientMagnitude_, node2Label_);*/
    }

  public:
    void preRemoveIsolatedNode(const Segmentation::DartTraverser &)
    {}
    void postRemoveIsolatedNode(Segmentation::FaceInfo &face)
    {
        lastChanges_ |= face.bounds;
    }

    void preRemoveBridge(const Segmentation::DartTraverser &dart)
    {
        preRemoveEdge(dart);
    }
    void postRemoveBridge(Segmentation::FaceInfo &face)
    {
        postRemoveEdge(face);
    }

    void preMergeFaces(const Segmentation::DartTraverser &dart)
    {
        // initialize tempFaceStatistics_ with leftFace and edge:
        preRemoveBridge(dart); // also sets node1Label_ and node2Label_

        tempFaceStatistics_(faceStatistics_[dart.rightFaceLabel()]);
    }
    void postMergeFaces(Segmentation::FaceInfo &face)
    {
        postRemoveBridge(face);
    }

    void preMergeEdges(const Segmentation::DartTraverser &dart)
    {
        Segmentation::DartTraverser d(dart);
        d.nextSigma();
        edge1Label_ = dart.edgeLabel();
        edge2Label_ = d.edgeLabel();

        if(!dart.leftFaceLabel() || !dart.rightFaceLabel())
        {
            tempEdgeStatistics_.clear();
            tempEdgeStatistics_(
                vigra::NumericTraits<GradientImage::PixelType>::max());
            return;
        }
        tempEdgeStatistics_ = edgeStatistics_[edge1Label_];
        tempEdgeStatistics_(edgeStatistics_[edge2Label_]);
        inspectCell(dart.segmentation()->nodeScanIterator
                    (d.startNodeLabel(), segmentationData_->gradientMagnitude_.upperLeft()),
                    tempEdgeStatistics_);
    }
    void postMergeEdges(Segmentation::EdgeInfo &edge)
    {
        edgeStatistics_[edge.label] = tempEdgeStatistics_;
        lastChanges_ |= edge.bounds;

        if(doEdgeRethinning)
            edgeRethinning(*edge.start.segmentation(),
                           segmentationData_->gradientMagnitude_, edge.label);
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
};

#endif // CELLSTATISTICS_HXX
