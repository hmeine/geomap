/************************************************************************/
/*                                                                      */
/*                Copyright 1998-99 by Ullrich Koethe                   */
/*        Fraunhoferinstitut fuer Graphische Datenverarbeitung,         */
/*                   Institutsteil Rostock, Germany                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is http://www.egd.igd.fhg.de/~ulli/vigra/       */
/*    Please direct questions, bug reports, and contributions to        */
/*                        ulli@egd.igd.fhg.de                           */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/distancetransform.hxx>
#include <vigra/convolution.hxx>
#include <vigra/functorexpression.hxx>
#include <vigra/labelimage.hxx>
#include "seededregiongrowing.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/impex.hxx"
#include "foureightsegmentation.hxx"
#include "hessematrix.hxx"
#include "gradient.hxx"
#include "findsaddles.hxx"
#include "mydebug.hxx"
#include "debugimage.hxx"

using namespace vigra;
using namespace vigra::functor;

void showNode(CellImage::FourEightSegmentation & seg, int label)
{
    CellImage::FourEightSegmentation::NodeAccessor node;
    CellImage::FourEightSegmentation::EdgeAccessor edge;
    CellImage::FourEightSegmentation::NodeAtEndAccessor endnode;
    CellImage::FourEightSegmentation::NodeIterator
        n = seg.findNode(label);

    std::cout << "Label: " << node.label(n) << ", Location: (" <<
		node.x(n) << ", " << node.y(n) << ")" << std::endl;

    CellImage::RayCirculator ray = node.rayCirculator(n);
    CellImage::RayCirculator rend = ray;

    do
    {
        std::cout << "( " << edge.label(ray) << ", " << endnode.label(ray) << ") ";
    }
    while(++ray != rend);

    std::cout << std::endl;
}

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
        int x,y;

        // create input image
        FImage in(w, h);

        importImage(info, destImage(in));

        FImage gradx(w, h), grady(w, h), hess(w,h), grad(w, h);


#if 0
        Kernel1D<float> symgrad;
        symgrad.initSymmetricGradient();

        separableConvolveX(srcImageRange(in), destImage(gradx), kernel1d(symgrad));
        separableConvolveY(srcImageRange(in), destImage(grady), kernel1d(symgrad));
#endif /* #if 0 */



#if 0
        recursiveFirstDerivativeX(srcImageRange(in), destImage(gradx), 0.8);
        recursiveSmoothY(srcImageRange(gradx), destImage(gradx), 0.8);

        recursiveFirstDerivativeY(srcImageRange(in), destImage(grady), 0.8);
        recursiveSmoothX(srcImageRange(grady), destImage(grady), 0.8);

        gradientMagnitude(srcImageRange(gradx), srcImage(grady), destImage(grad));

#endif /* #if 0 */

        gradientMagnitude(srcImageRange(in), destImage(grad));

#if 0
        exportImage(srcImageRange(grad), ImageExportInfo("grad.xv"));
        std::cout << "Wrote grad.xv" << std::endl;
#endif /* #if 0 */


        IImage labels(w, h), seeds(w, h), seeds1(w,h);
        seeds = 0;
        labels = 0;

        extendedLocalMinima(srcImageRange(grad), destImage(seeds), 1);

        determinantOfHessian(srcImageRange(in), destImage(hess));

        combineTwoImages(srcImageRange(hess), srcImage(seeds), destImage(seeds),
                         ifThenElse(Arg1() < Param(0), Param(0), Arg2()));

        float threshold;
        std::cout << "Threshold ? ";
        std::cin >> threshold;

        combineTwoImages(srcImageRange(grad), srcImage(seeds), destImage(seeds1),
                         ifThenElse(Arg1() < Param(threshold), Param(1), Arg2()));

        int raw_region_label =
            labelImageWithBackground(srcImageRange(seeds1), destImage(labels), true, 0);

        ArrayOfRegionStatistics<FindMinMax<float> > minmax(raw_region_label);
        inspectTwoImages(srcImageRange(seeds), srcImage(labels), minmax);

        for(y=0; y<h; ++y)
        {
            for(x=0; x<w; ++x)
            {
                if(labels(x,y) == 0) continue;
                if(minmax[labels(x,y)].max == 0)
                {
                    seeds(x,y) = 0;
                }
                else
                {
                    seeds(x,y) = 1;
                }

                labels(x,y) = 0;
            }
        }

        exportImage(srcImageRange(seeds), ImageExportInfo("seeds.xv"));
        std::cout << "Wrote seeds.xv" << std::endl;

        int maxFaceLabel =
            labelImageWithBackground(srcImageRange(seeds), destImage(labels), false, 0);

        // create a statistics functor for region growing
        ArrayOfRegionStatistics<SeedRgDirectValueFunctor<float> >
			gradstat(maxFaceLabel);

        // perform region growing, starting from the minima of the gradient magnitude;
        // as the feature (first input) image contains the gradient magnitude,
        // this calculates the catchment basin of each minimum
        seededRegionGrowing(srcImageRange(grad), srcImage(labels),
                            destImage(labels), gradstat, 100000);

        // remove regions of size 1
        for(y=0; y<h; ++y)
        {
            for(x=0; x<w; ++x)
            {
                int lab = labels(x,y);
                if(lab == 0) continue;

                int dx[] = {1, 0, -1, 0};
                int dy[] = {0, -1, 0, 1};
                int i;
                for(i=0; i<4; ++i)
                {
                    int xx = x + dx[i];
                    int yy = y + dy[i];

                    if(xx >= 0 && xx < w && yy >= 0 && yy < h)
                    {
                        if(labels(xx,yy) == lab) break;
                    }
                }

                if(i == 4) labels(x,y) = 0;
            }
        }

        // repeat region growing
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

        std::cout <<
			"Nodes: " << segmentation.nodeCount() << std::endl <<
			"Edges: " << segmentation.edgeCount() << std::endl <<
			"Faces: " << segmentation.faceCount() << std::endl;

        // initialize a functor to determine the average gray-value or color
        // for each region (catchment basin) just found
        ArrayOfRegionStatistics<FindAverage<float> >
            averages(maxFaceLabel);

        // calculate the averages
        inspectTwoImages(srcImageRange(in), srcImage(labels), averages);

        FindAverage<float> zero;
        zero(0);
        averages[0] = zero;

        // write the averages into the destination image (the functor 'averages'
        // acts as a look-up table)
        BImage out(w,h);
        transformImage(srcImageRange(labels), destImage(out), averages);


        // mark the watershaed (region boundaries) black

#if 0
        regionImageToEdgeImage(srcImageRange(labels), destImage(in), 0);
#endif /* #if 0 */

//        findSaddles(srcImageRange(in), destImage(out), 255);
        BImage saddle(w,h);
        saddle = 0;
        localMinima(srcImageRange(hess), destImage(saddle), 1);
        combineTwoImages(srcImageRange(saddle), srcImage(hess),
                         destImage(saddle),
						 ifThenElse(Arg2() < Param(-1.0), Arg1(), Param(0)));
        initImageIf(destImageRange(out), maskImage(saddle), 255);

        exportImage(srcImageRange(out), ImageExportInfo("res.xv"));
        std::cout << "Wrote res.xv" << std::endl;
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
