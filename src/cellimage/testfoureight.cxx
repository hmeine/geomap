#include <iostream>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/functorexpression.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/impex.hxx>
#include "seededregiongrowing.hxx"
#include "foureightsegmentation.hxx"
#include "hessematrix.hxx"
#include "gradient.hxx"
#include "findsaddles.hxx"
#include "mydebug.hxx"
#include "debugimage.hxx"
#include "cellimage.hxx"

using namespace vigra;
using namespace vigra::functor;

static const unsigned char imageData[] =
    { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 139,  52,   6,   6,  52, 139, 255, 255, 255, 255, 255,
      255, 255, 255, 221,  52,   0,   0,   0,   0,   0,   0,  52, 221, 255, 255, 255,
      255, 255, 255,  52,   0,   0,   0,   0,   0,   0,   0,   0,  52, 255, 255, 255,
      255, 255, 139,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 139, 255, 255,
      255, 255,  52,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  52, 255, 255,
      255, 255,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6, 255, 255,
      255, 255,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6, 255, 255,
      255, 255,  52,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  52, 255, 255,
      255, 255, 139,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 139, 255, 255,
      255, 255, 255,  52,   0,   0,   0,   0,   0,   0,   0,   0,  52, 255, 255, 255,
      255, 255, 255, 221,  52,   0,   0,   0,   0,   0,   0,  52, 221, 255, 255, 255,
      255, 255, 255, 255, 255, 139,  52,   6,   6,  52, 139, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
      255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

int main(int argc, char ** argv)
{
    try
    {
        vigra::BImage image(16, 16);
        std::copy(imageData, imageData + 256, image.begin());
        int w = image.width();
        int h = image.height();

        FImage grad(w, h);
        gradientMagnitude(srcImageRange(image), destImage(grad));

        IImage seeds(w, h);
        seeds = 0;

        extendedLocalMinima(srcImageRange(grad), destImage(seeds), 1);

        FImage hess(w,h);
        determinantOfHessian(srcImageRange(image), destImage(hess));

        combineTwoImages(srcImageRange(hess), srcImage(seeds), destImage(seeds),
                         ifThenElse(Arg1() < Param(0), Param(0), Arg2()));

        IImage labels(w, h);
        labels = 0;
        int maxFaceLabel =
            labelImageWithBackground(srcImageRange(seeds), destImage(labels),
                                     false, 0);

        // create a statistics functor for region growing
        ArrayOfRegionStatistics<SeedRgDirectValueFunctor<float> >
            gradstat(maxFaceLabel);

        // perform region growing, starting from the minima of the
        // gradient magnitude; as the feature (first input) image
        // contains the gradient magnitude, this calculates the
        // catchment basin of each minimum
        seededRegionGrowing(srcImageRange(grad), srcImage(labels),
                            destImage(labels), gradstat, 100000);

        CellImage::FourEightSegmentation segmentation;

        segmentation.init(srcImageRange(labels));

        if(segmentation.cellImage.width()>50)
        {
            /*exportImage(
                srcImageRange(segmentation.cellImage,
                              CellImage::CellImageTypeAccessor()),
                ImageExportInfo("cellTypeImage.xv"));
                std::cout << "Wrote cellTypeImage.xv" << std::endl;*/

            exportImage(
                srcImageRange(segmentation.cellImage,
                              CellImage::CellImageLabelAccessor()),
                ImageExportInfo("cellLabelImage.xv"));
            std::cout << "Wrote cellLabelImage.xv" << std::endl;
        }
        else
            debugImage(srcImageRange(segmentation.cellImage));

        std::cout << segmentation.nodeCount() << " nodes:" << std::endl;
        for(CellImage::FourEightSegmentation::NodeIterator
                node= segmentation.nodesBegin(); node.inRange(); ++node)
            std::cout << "  " << node->label << ": at "
                      << node->centerX << ","
                      << node->centerY << ", "
                      << node->bounds << std::endl;
        
        std::cout << segmentation.edgeCount() << " edges:" << std::endl;
        for(CellImage::FourEightSegmentation::EdgeIterator
                edge= segmentation.edgesBegin(); edge.inRange(); ++edge)
        {
            vigra::FindAverage<vigra::FImage::PixelType> average;
            if(edge->label>4)
                vigra::CellImage::inspectCell
                    (segmentation.edgeScanIterator(edge->label, image.upperLeft()),
                     average);
            
            std::cout << "  " << edge->label << ": " << edge->bounds
					  << ", average: " << average() << std::endl;
        }
        
        std::cout << segmentation.faceCount() << " faces:" << std::endl;
        for(CellImage::FourEightSegmentation::FaceIterator
                face= segmentation.facesBegin(); face.inRange(); ++face)
        {
            vigra::FindAverage<vigra::FImage::PixelType> average;
            if(face->label>0)
                vigra::CellImage::inspectCell
                    (segmentation.faceScanIterator(face->label, image.upperLeft()),
                     average);
            
            std::cout << "  " << face->label << ": " << face->bounds
                      << ", average: " << average() << std::endl;
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
