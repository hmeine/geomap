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
};

template <class SrcImageIterator, class SrcAccessor>
void CrackEdgeMapGenerator::markEightConnectedRegions(
    SrcImageIterator sul, SrcAccessor sa)
{
    vigra::IImage::traverser
        cul = crackConnections.upperLeft(),
        row = crackConnections.upperLeft() + vigra::Diff2D(1, 1),
        end = crackConnections.lowerRight() - vigra::Diff2D(1, 1);

    // FIXME: use sa
    for(; row.y < end.y; ++row.y)
    {
        vigra::IImage::traverser it = row;
        for(; it.x < end.x; ++it.x)
        {
            int conn(*it);
            if(conn != CONN_ALL4)
                continue;

            if(sul[it - cul] ==
               sul[it - cul - vigra::Diff2D(1, 1)])
                *it |= CONN_DIAG_UPLEFT;
            if(sul[it - cul - vigra::Diff2D(1, 0)] ==
               sul[it - cul - vigra::Diff2D(0, 1)])
                *it |= CONN_DIAG_UPRIGHT;

            // crossing regions?
            if(*it == (CONN_ALL4 | CONN_DIAG_UPLEFT | CONN_DIAG_UPRIGHT))
                // preserve connectedness of higher label:
                if(sul[it - cul - vigra::Diff2D(0, 1)] >
                   sul[it - cul - vigra::Diff2D(1, 1)])
                    *it -= CONN_DIAG_UPLEFT;
                else
                    *it -= CONN_DIAG_UPRIGHT;
        }
    }
}

#endif // CRACKEDGEMAP_HXX
