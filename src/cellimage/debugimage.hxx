#ifndef VIGRA_DEBUGIMAGE_HXX
#define VIGRA_DEBUGIMAGE_HXX

#include <iostream>

namespace vigra {

template <class SrcIter, class SrcAcc>
void debugImage(SrcIter sul, SrcIter slr, SrcAcc sa,
                std::ostream &outs = std::cerr)
{
    int w= (slr-sul).x, h= (slr-sul).y;
    for(SrcIter row= sul; row.y<slr.y; row.y++)
    {
        for(SrcIter it= row; it.x<slr.x; it.x++)
            outs << sa(it);
        outs << std::endl;
    }
}

template <class SrcIter, class SrcAcc>
void debugImage(triple<SrcIter, SrcIter, SrcAcc> src,
                std::ostream &outs = std::cerr)
{
    debugImage(src.first, src.second, src.third, outs);
}

} // namespace vigra

#endif // VIGRA_DEBUGIMAGE_HXX
