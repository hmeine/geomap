#ifndef CELLSTATISTICS_HXX
#define CELLSTATISTICS_HXX

#include <foureightsegmentation.hxx>
#include <statisticfunctor.hxx>

#include <vigra/stdimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vigra/rect2d.hxx>
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

void nodeRethinning(vigra::cellimage::FourEightSegmentation &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel nodeLabel);

void edgeRethinning(vigra::cellimage::FourEightSegmentation &seg,
                    const GradientImage &gradientMagnitude,
                    vigra::cellimage::CellLabel edgeLabel);

struct CellStatistics
{
    typedef vigra::cellimage::FourEightSegmentation
        Segmentation;

        // members storing source data
    SegmentationData *segmentationData_;
    mutable vigra::Rect2D lastChanges_;

        // members storing statistics
    std::vector<OriginalImage::PixelType>
        faceMeans_;
    std::vector<vigra::NumericTraits<OriginalImage::PixelType>::RealPromote>
        faceVariance_; // sigma squared
    std::vector<GradientImage::PixelType>
        meanEdgeGradients_;
    std::vector<vigra::cellimage::CellLabel>
        mergedEdges_;
    std::vector<vigra::TinyVector<float, 2> >
        nodeCenters_;

        // constant information, could even be static?
    std::vector<Float2D>
        configurationDirections_;

    CellStatistics(const Segmentation &initialSegmentation,
                   SegmentationData  *segmentationData);

        // temporary data for storage between preXXX- and postXXX-calls
    StatisticFunctor<GradientImage::PixelType> tempStatistics_;
    vigra::FindAverage<OriginalImage::PixelType> tempAverage_;
    vigra::cellimage::CellLabel node1Label_, node2Label_;
    vigra::cellimage::CellLabel edge1Label_, edge2Label_;

    void preRemoveIsolatedNode(const Segmentation::DartTraverser &)
    {}
    void postRemoveIsolatedNode(const Segmentation::DartTraverser &,
                                Segmentation::FaceInfo &face)
    {
        lastChanges_ |= face.bounds;
    }

    void preRemoveBridge(const Segmentation::DartTraverser &dart)
    {
        // /!\ common code (called from preMergeFaces, too):
        vigra::cellimage::CellLabel faceLabel = dart.leftFaceLabel();
        tempStatistics_.count= dart.segmentation()->face(faceLabel).size;
        tempStatistics_.sum= tempStatistics_.count * faceMeans_[faceLabel];
        tempStatistics_.squareSum=
            tempStatistics_.count * (faceVariance_[faceLabel] +
                                    faceMeans_[faceLabel] * faceMeans_[faceLabel]);

        inspectCell(dart.segmentation()->edgeScanIterator
                    (dart.edgeLabel(), segmentationData_->preparedOriginal_.upperLeft()),
                    tempStatistics_);

        node1Label_ = dart.edge().start.startNodeLabel();
        node2Label_ = dart.edge().end.startNodeLabel();
    }
    void postRemoveBridge(const Segmentation::DartTraverser &,
                          Segmentation::FaceInfo &face)
    {
        faceMeans_[face.label]= tempStatistics_.avg();
        faceVariance_[face.label]= tempStatistics_.var();
        lastChanges_ |= face.bounds;

        /*nodeRethinning(*dart.segmentation(),
                       segmentationData_->gradientMagnitude_, node1Label_);
        nodeRethinning(*dart.segmentation(),
        segmentationData_->gradientMagnitude_, node2Label_);*/
    }

    void preMergeFaces(const Segmentation::DartTraverser &dart)
    {
        // initialize tempStatistics_ with leftFace and edge:
        preRemoveBridge(dart); // also sets node1Label_ and node2Label_

        vigra::cellimage::CellLabel face2Label = dart.rightFaceLabel();
        long count = dart.segmentation()->face(face2Label).size;
        tempStatistics_.count += count;
        tempStatistics_.sum += count * faceMeans_[face2Label];
        tempStatistics_.squareSum +=
            count * (faceVariance_[face2Label] +
                     faceMeans_[face2Label] * faceMeans_[face2Label]);
    }
    void postMergeFaces(const Segmentation::DartTraverser &,
                        Segmentation::FaceInfo &face)
    {
        faceMeans_[face.label]    = tempStatistics_.avg();
        faceVariance_[face.label] = tempStatistics_.var();

        lastChanges_ |= face.bounds;

        /*nodeRethinning(*dart.segmentation(),
                       segmentationData_->gradientMagnitude_,  node1Label_);
        if(node1Label_ != node2Label_)
            nodeRethinning(*dart.segmentation(),
            segmentationData_->gradientMagnitude_, node2Label_);*/
    }

    void preMergeEdges(const Segmentation::DartTraverser &dart)
    {
        Segmentation::DartTraverser d(dart);
        d.nextSigma();
        edge1Label_ = dart.edgeLabel();
        edge2Label_ = d.edgeLabel();

        if(!dart.leftFaceLabel() || !dart.rightFaceLabel())
        {
            tempAverage_.count = 1;
            tempAverage_.sum = vigra::NumericTraits<GradientImage::PixelType>::max();
            return;
        }
        tempAverage_.count = 0;
        tempAverage_.sum = 0.0f;
        inspectCell(dart.segmentation()->edgeScanIterator
                    (edge1Label_, segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage_);
        inspectCell(dart.segmentation()->edgeScanIterator
                    (edge2Label_, segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage_);
        inspectCell(dart.segmentation()->nodeScanIterator
                    (d.startNodeLabel(), segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage_);
    }
    void postMergeEdges(const Segmentation::DartTraverser &dart,
                        Segmentation::EdgeInfo &edge)
    {
        meanEdgeGradients_[edge.label]= tempAverage_();
        lastChanges_ |= edge.bounds;

        edgeRethinning(*dart.segmentation(),
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
