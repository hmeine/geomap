#ifndef CELLSTATISTICS_HXX
#define CELLSTATISTICS_HXX

#include <foureightsegmentation.hxx>
#include <statisticfunctor.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/inspectimage.hxx>
#include <vector>

typedef vigra::FImage OriginalImage;
typedef vigra::FImage GradientImage;

struct SegmentationData
{
    OriginalImage preparedOriginal_;
    GradientImage gradientMagnitude_;
};

struct CellStatistics
{
    typedef vigra::cellimage::FourEightSegmentation Segmentation;

    SegmentationData  *segmentationData_;
    StatisticFunctor<GradientImage::PixelType> tempStatistics;
    vigra::FindAverage<OriginalImage::PixelType> tempAverage;

    std::vector<OriginalImage::PixelType> faceMeans_;
    std::vector<vigra::NumericTraits<OriginalImage::PixelType>::RealPromote>
        faceVariance_; // sigma squared
    std::vector<GradientImage::PixelType> meanEdgeGradients_;

    CellStatistics(const Segmentation &initialSegmentation,
                   SegmentationData  *segmentationData)
    : segmentationData_(segmentationData)
    {
        std::cerr << "initializing face statistics (max face label: "
                  << initialSegmentation.maxFaceLabel() << ")\n";
        faceMeans_.resize(initialSegmentation.maxFaceLabel() + 1);
        faceVariance_.resize(initialSegmentation.maxFaceLabel() + 1);
        faceMeans_[0]= vigra::NumericTraits<OriginalImage::PixelType>::max();
        faceVariance_[0]= 0.0;
        Segmentation::ConstFaceIterator it =
                initialSegmentation.facesBegin();
        for(++it; it.inRange(); ++it)
        {
            std::cerr << "inspecting face " << it->label << "\r";
            StatisticFunctor<OriginalImage::PixelType> faceStatistics;
            inspectCell(initialSegmentation.faceScanIterator
                        (it->label, segmentationData_->preparedOriginal_.upperLeft()),
                        faceStatistics);
            faceMeans_[it->label]= faceStatistics.avg();
            faceVariance_[it->label]= faceStatistics.var();
            if(!faceVariance_[it->label])
            {
                //std::cerr << "setting faceVariance_[" << it->label << "] to epsilon.!!\n";
                faceVariance_[it->label]= vigra::NumericTraits<float>::epsilon();
            }
        }

        std::cerr << "initializing meanEdgeGradients\n";
        meanEdgeGradients_.resize(initialSegmentation.maxEdgeLabel() + 1);
        for(Segmentation::ConstEdgeIterator it =
                initialSegmentation.edgesBegin();
            it.inRange(); ++it)
        {
            const Segmentation::DartTraverser &anchor = it->start;
            if(!anchor.leftFaceLabel() || !anchor.rightFaceLabel())
                meanEdgeGradients_[it->label]=
                    vigra::NumericTraits<GradientImage::PixelType>::max();
            else
            {
                vigra::FindAverage<GradientImage::PixelType> edgeMean;
                inspectCell(initialSegmentation.edgeScanIterator
                            (it->label, segmentationData_->gradientMagnitude_.upperLeft()),
                            edgeMean);
                meanEdgeGradients_[it->label]= edgeMean();
            }
        }
    }

    void preRemoveIsolatedNode(const Segmentation::DartTraverser &)
    {}
    void postRemoveIsolatedNode(const Segmentation::DartTraverser &,
                                Segmentation::FaceInfo &)
    {}

    void preRemoveBridge(const Segmentation::DartTraverser &dart)
    {
        vigra::cellimage::CellLabel faceLabel = dart.leftFaceLabel();
        tempStatistics.count= dart.segmentation()->face(faceLabel).size;
        tempStatistics.sum= tempStatistics.count * faceMeans_[faceLabel];
        tempStatistics.squareSum=
            tempStatistics.count * (faceVariance_[faceLabel] +
                                    faceMeans_[faceLabel] * faceMeans_[faceLabel]);

        inspectCell(dart.segmentation()->edgeScanIterator
                    (dart.edgeLabel(), segmentationData_->preparedOriginal_.upperLeft()),
                    tempStatistics);
    }
    void postRemoveBridge(const Segmentation::DartTraverser &,
                          Segmentation::FaceInfo &face)
    {
        faceMeans_[face.label]= tempStatistics.avg();
        faceVariance_[face.label]= tempStatistics.var();
    }

    void preMergeFaces(const Segmentation::DartTraverser &dart)
    {
        preRemoveBridge(dart);

        vigra::cellimage::CellLabel face2Label = dart.rightFaceLabel();
        long count = dart.segmentation()->face(face2Label).size;
        tempStatistics.count += count;
        tempStatistics.sum += count * faceMeans_[face2Label];
        tempStatistics.squareSum +=
            count * (faceVariance_[face2Label] +
                     faceMeans_[face2Label] * faceMeans_[face2Label]);
    }
    void postMergeFaces(const Segmentation::DartTraverser &,
                        Segmentation::FaceInfo &face)
    {
        faceMeans_[face.label]= tempStatistics.avg();
        faceVariance_[face.label]= tempStatistics.var();
    }

    void preMergeEdges(const Segmentation::DartTraverser &dart)
    {
        if(!dart.leftFaceLabel() || !dart.rightFaceLabel())
        {
            tempAverage.count = 1;
            tempAverage.sum = vigra::NumericTraits<GradientImage::PixelType>::max();
            return;
        }
        Segmentation::DartTraverser d(dart);
        tempAverage.count = 0;
        tempAverage.sum = 0.0f;
        inspectCell(dart.segmentation()->edgeScanIterator
                    (d.edgeLabel(), segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage);
        d.nextSigma();
        inspectCell(dart.segmentation()->edgeScanIterator
                    (d.edgeLabel(), segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage);
        inspectCell(dart.segmentation()->nodeScanIterator
                    (d.startNodeLabel(), segmentationData_->gradientMagnitude_.upperLeft()),
                    tempAverage);
    }
    void postMergeEdges(const Segmentation::DartTraverser &,
                        Segmentation::EdgeInfo &edge)
    {
        meanEdgeGradients_[edge.label]= tempAverage();
    }
};

#endif // CELLSTATISTICS_HXX
