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
#include "debugimage.hxx"
#include "cellimage.hxx"
#include "mydebug.hxx"

using namespace vigra;
using namespace vigra::functor;

int main(int argc, char ** argv)
{
    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " infile" << std::endl;
        std::cout << "(supported fomats: " << impexListFormats() << ")" <<
                     std::endl;

        return 1;
    }
    try
    {
        ImageImportInfo info(argv[1]);

        vigra_precondition(info.isGrayscale(), "Unable to process color images");

        int w = info.width();
        int h = info.height();

        // create input image
        FImage in(w, h);
        importImage(info, destImage(in));

        FImage grad(w, h);
        gradientMagnitude(srcImageRange(in), destImage(grad));

        IImage seeds(w, h);
        seeds = 0;

        extendedLocalMinima(srcImageRange(grad), destImage(seeds), 1);

        FImage hess(w,h);
        determinantOfHessian(srcImageRange(in), destImage(hess));

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
        for(int node= 0; node< segmentation.nodeList.size(); ++node)
            if(segmentation.nodeList[node].initialized())
                std::cout << "  " << node << ": at "
                          << segmentation.nodeList[node].centerX << ","
                          << segmentation.nodeList[node].centerY << "  from "
                          << segmentation.nodeList[node].upperLeft << " to "
                          << segmentation.nodeList[node].lowerRight << std::endl;
		
        std::cout << segmentation.edgeCount() << " edges:" << std::endl;
        for(int edge= 0; edge< segmentation.edgeList.size(); ++edge)
            if(segmentation.edgeList[edge].initialized())
			{
				vigra::FindAverage<vigra::FImage::PixelType> average;
				if(edge>4)
					vigra::inspectCell(segmentation.cellScanIterator
									   (segmentation.edgeList[edge], CellTypeLine,
										in.upperLeft()),
									   average);

                std::cout << "  " << edge << ": from "
                          << segmentation.edgeList[edge].upperLeft << " to "
                          << segmentation.edgeList[edge].lowerRight
						  << ", average: " << average() << std::endl;
			}
		
        std::cout << segmentation.faceCount() << " faces:" << std::endl;
        for(int face= 0; face< segmentation.faceList.size(); ++face)
            if(segmentation.faceList[face].initialized())
			{
				vigra::FindAverage<vigra::FImage::PixelType> average;
				if(face>0)
					vigra::inspectCell(segmentation.cellScanIterator
									   (segmentation.faceList[face], CellTypeRegion,
										in.upperLeft()),
									   average);

                std::cout << "  " << face << ": from "
                          << segmentation.faceList[face].upperLeft << " to "
                          << segmentation.faceList[face].lowerRight
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
