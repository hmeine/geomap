#ifndef VIGRA_DEBUGIMAGE_HXX
#define VIGRA_DEBUGIMAGE_HXX

#include <iostream>

namespace vigra {

template <class SrcIter, class SrcAcc, class ValueType>
void debugImageImpl(SrcIter sul, SrcIter slr, SrcAcc sa,
                    ValueType, std::ostream &outs)
{
    for(SrcIter row= sul; row.y<slr.y; row.y++)
    {
        for(SrcIter it= row; it.x<slr.x; it.x++)
            outs << sa(it);
        outs << std::endl;
    }
}

// partial specialization for chars to output them as int
template <class SrcIter, class SrcAcc>
void debugImageImpl(SrcIter sul, SrcIter slr, SrcAcc sa,
                    unsigned char, std::ostream &outs)
{
    for(SrcIter row= sul; row.y<slr.y; row.y++)
    {
        for(SrcIter it= row; it.x<slr.x; it.x++)
            outs << (int)sa(it);
        outs << std::endl;
    }
}

template <class SrcIter, class SrcAcc>
void debugImage(SrcIter sul, SrcIter slr, SrcAcc sa,
                std::ostream &outs = std::cerr)
{
    if(sul != slr)
        debugImageImpl(sul, slr, sa, *sul, outs);
}

template <class SrcIter, class SrcAcc>
void debugImage(triple<SrcIter, SrcIter, SrcAcc> src,
                std::ostream &outs = std::cerr)
{
    debugImage(src.first, src.second, src.third, outs);
}

} // namespace vigra

#endif // VIGRA_DEBUGIMAGE_HXX
