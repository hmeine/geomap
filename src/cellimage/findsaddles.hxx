#ifndef VIGRA_FINDSADDLES_HXX
#define VIGRA_FINDSADDLES_HXX

#include <vigra/pixelneighborhood.hxx>

namespace vigra {

template <class SIter, class SAcc, class DIter, class DAcc, class Value>
void findSaddles(SIter sul, SIter slr, SAcc src,
                 DIter dul, DAcc dest, Value saddlemarker)
{
    int x, y;
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    typedef typename NumericTraits<typename SAcc::value_type>::Promote TmpType;
    static TmpType zero = NumericTraits<TmpType>::zero();

    ++sul.y;
    ++dul.y;
    for(y=1; y<h-1; ++y, ++sul.y, ++dul.y)
    {
        SIter scur = sul;
        DIter dcur = dul;

        ++scur.x;
        ++dcur.x;
        for(x=1; x<w-1; ++x, ++scur.x, ++dcur.x)
        {
            TmpType grads[8];

            NeighborhoodCirculator<SIter, EightNeighborOffsetCirculator>
                neigh(scur);

            for(int i=0; i<8; ++i, ++neigh)
            {
                grads[i] = *neigh - *scur;
            }

            int up = 0;
            for(int i=0; i<8; ++i)
            {
                if(zero < grads[i]*grads[(i+7)%8]) continue;

                if(zero < grads[i]) ++up;
            }

            if(up > 1) dest.set(saddlemarker, dcur);
        }
    }
}

template <class SIter, class SAcc, class DIter, class DAcc, class Value>
inline
void findSaddles(triple<SIter, SIter, SAcc> src,
                 pair<DIter, DAcc> dest, Value saddlemarker)
{
    findSaddles(src.first, src.second, src.third,
                         dest.first, dest.second, saddlemarker);
}

} // namespace vigra

#endif /* VIGRA_FINDSADDLES_HXX */
