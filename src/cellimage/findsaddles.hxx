/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

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
