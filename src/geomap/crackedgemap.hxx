#ifndef CRACKEDGEMAP_HXX
#define CRACKEDGEMAP_HXX

#include "cppmap.hxx"
#include "vigra/crackconnections.hxx"
#include <vigra/stdimage.hxx>
#include <vigra/pixelneighborhood.hxx>
#include <memory>

class CrackEdgeMapGenerator
{
  public:
    template <class SrcImageIterator, class SrcAccessor>
    CrackEdgeMapGenerator(
        vigra::triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
        bool eightConnectedRegions = false)
    : crackConnections((src.second - src.first) + vigra::Size2D(1, 1), 0),
      nodeImage(crackConnections.size(), 0),
      result(new GeoMap(vigra::Size2D(src.second - src.first)))
    {
        crackConnectionImage(src, destImage(crackConnections));
        makeCCSymmetric();
        if(eightConnectedRegions)
            markEightConnectedRegions(src.first, src.third);
        markNodes();
        followAllEdgesStartingWith(CONN_NODE);
        followAllEdgesStartingWith(CONN_MAYBE_NODE);
    }

//     template <class SrcImageIterator, class SrcAccessor>
//     CrackEdgeMapGenerator(
//         SrcImageIterator sul,
//         SrcImageIterator slr, SrcAccessor sa)
//     : crackConnections((slr - sul) + vigra::Size2D(1, 1)),
//       nodeImage(crackConnections.size()),
//       result(slr - sul)
//     {
//         crackConnectionImage(sul, slr, sa, crackConnections.upperLeft(), crackConnections.accessor());
//         makeCCSymmetric();
//     }

    enum {
        CONN_RIGHT = 1,
        CONN_DOWN = 2,
        CONN_LEFT = 4,
        CONN_UP = 8,
        CONN_ALL4 = 15,

        CONN_DIAG_UPLEFT = 16,
        CONN_DIAG_UPRIGHT = 32,
        CONN_DIAG = CONN_DIAG_UPLEFT | CONN_DIAG_UPRIGHT,

        CONN_NODE = 64,
        CONN_MAYBE_NODE = 128,
        CONN_ANYNODE = CONN_NODE | CONN_MAYBE_NODE,
    };

    vigra::IImage crackConnections, nodeImage;
    std::auto_ptr<GeoMap> result;

    void makeCCSymmetric();

    template <class SrcImageIterator, class SrcAccessor>
    void markEightConnectedRegions(
        SrcImageIterator sul, SrcAccessor sa);

    void markNodes();

    std::auto_ptr<Vector2Array>
    followEdge(vigra::Point2D &startPos,
               vigra::FourNeighborOffsetCirculator &dir);

    void followAllEdgesStartingWith(int connMask);

    void initializeMap(bool initLabelImage);

    template <class SrcImageIterator, class SrcAccessor>
    void getSourceLabels(vigra::pair<SrcImageIterator, SrcAccessor> src,
                         std::vector<CellLabel> *sourceLabels);
};

template <class SrcImageIterator, class SrcAccessor>
void CrackEdgeMapGenerator::markEightConnectedRegions(
    SrcImageIterator sul, SrcAccessor sa)
{
    vigra::IImage::traverser
        cul = crackConnections.upperLeft(),
        row = crackConnections.upperLeft() + vigra::Diff2D(1, 1),
        end = crackConnections.lowerRight() - vigra::Diff2D(1, 1);

    for(; row.y < end.y; ++row.y)
    {
        vigra::IImage::traverser it = row;
        for(; it.x < end.x; ++it.x)
        {
            int conn(*it);
            if(conn != CONN_ALL4)
                continue;

            if(sa(sul, it - cul) ==
               sa(sul, it - cul - vigra::Diff2D(1, 1)))
                *it |= CONN_DIAG_UPLEFT;
            if(sa(sul, it - cul - vigra::Diff2D(1, 0)) ==
               sa(sul, it - cul - vigra::Diff2D(0, 1)))
                *it |= CONN_DIAG_UPRIGHT;

            // crossing regions?
            if(*it == (CONN_ALL4 | CONN_DIAG_UPLEFT | CONN_DIAG_UPRIGHT))
            {
                // preserve connectedness of higher label:
                if(sa(sul, it - cul - vigra::Diff2D(0, 1)) >
                   sa(sul, it - cul - vigra::Diff2D(1, 1)))
                    *it -= CONN_DIAG_UPLEFT;
                else
                    *it -= CONN_DIAG_UPRIGHT;
            }
        }
    }
}

template <class SrcImageIterator, class SrcAccessor>
void CrackEdgeMapGenerator::getSourceLabels(vigra::pair<SrcImageIterator, SrcAccessor> src,
                                            std::vector<CellLabel> *sourceLabels)
{
    vigra_precondition(result->mapInitialized(),
                       "getSourceLabels: map not yet initialized -> no faces to label");
    vigra_precondition(sourceLabels->size() >= result->maxFaceLabel(),
                       "getSourceLabels: target array not large enough");

    for(GeoMap::FaceIterator it = result->finiteFacesBegin(); it.inRange(); ++it)
    {
        GeoMap::Dart dart((*it)->contour());

        Vector2 p1(dart[0]), p2(dart[1]);

        // candidate pos of integer/pixel position to the left of the
        // current crack (will be adjusted depending on actual crack
        // direction below), cheaply rounded to the lower right:
        vigra::Diff2D leftPos((int)(p1[0] + 1), (int)(p1[1] + 1));

        if(p1[0] == p2[0])
        {
            if(p1[1] < p2[1]) // downward crack
            {
                // above offset already correct
            }
            else              // upward crack, adjust offset
            {
                --leftPos.x;
                --leftPos.y;
            }
        }
        else
        {
            if(p1[0] < p2[0]) // eastward crack, adjust offset
            {
                --leftPos.y;
            }
            else              // westward crack, adjust offset
            {
                --leftPos.x;
            }
        }

        (*sourceLabels)[(*it)->label()] = src.second(src.first, leftPos);
    }
}

void crackEdgesToMidcracks(GeoMap &geomap);

/**
 * Change edge geometry (in-place) from midcracks to marching
 * squares-like level contours.  Use linear interpolation of the `src`
 * image to determine the threshold position.  Leave alone edges with
 * BORDER_PROTECTION set.  Assumes that the `src` image is (at least)
 * the same size as geomap.imageSize().
 */
template<class ITERATOR, class ACCESSOR>
void midcracksToThreshold(GeoMap &geomap,
                          vigra::pair<ITERATOR, ACCESSOR> src,
                          double threshold)
{
    vigra::Size2D size(geomap.imageSize());

    for(GeoMap::EdgeIterator it = geomap.edgesBegin(); it.inRange(); ++it)
    {
        GeoMap::Edge &edge(**it);
        if(edge.flag(GeoMap::Edge::BORDER_PROTECTION))
            continue;

        for(unsigned int i = 0; i < edge.size(); ++i)
        {
            double x = edge[i][0], y = edge[i][1];
            int ix = (int)floor(x), iy = (int)floor(y);
            if(ix < 0 || iy < 0)
                continue;

            // every mid-crack point has exactly one crack coordinate
            // and one pixel grid coordinate; we only want to change
            // the crack coordinate:
            if(x == ix)
            {
                if(ix >= size.x || (iy + 1) >= size.y)
                    continue;

                // move vertically
                typename ITERATOR::value_type
                    v1 = src.second(src.first, vigra::Diff2D(ix, iy)),
                    v2 = src.second(src.first, vigra::Diff2D(ix, iy + 1));
                edge[i][1] = iy + (threshold - v1) / (v2 - v1);
            }
            else
            {
                if((ix + 1) >= size.x || iy >= size.y)
                    continue;

                // move horizontally
                typename ITERATOR::value_type
                    v1 = src.second(src.first, vigra::Diff2D(ix, iy)),
                    v2 = src.second(src.first, vigra::Diff2D(ix + 1, iy));
                edge[i][0] = ix + (threshold - v1) / (v2 - v1);
            }
        }
    }

    if(geomap.mapInitialized())
    {
        for(GeoMap::FaceIterator it = geomap.facesBegin(); it.inRange(); ++it)
        {
            // FIXME: see note about setGeometry above...
            (*it)->setFlag(0xC0000000U, false);
        }
    }
}

#endif // CRACKEDGEMAP_HXX
