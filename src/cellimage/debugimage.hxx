#ifndef VIGRA_DEBUGIMAGE_HXX
#define VIGRA_DEBUGIMAGE_HXX

#include <iostream>
#include <iomanip>

namespace vigra {

template <class SrcIter, class SrcAcc, class ValueType>
void debugImageImpl(SrcIter sul, SrcIter slr, SrcAcc sa,
                    ValueType, std::ostream &outs, int width)
{
    for(SrcIter row= sul; row.y<slr.y; row.y++)
    {
        for(SrcIter it= row; it.x<slr.x; it.x++)
            outs << std::setw(width) << sa(it);
        outs << std::endl;
    }
}

// partial specialization for chars to output them as int
template <class SrcIter, class SrcAcc>
void debugImageImpl(SrcIter sul, SrcIter slr, SrcAcc sa,
                    unsigned char, std::ostream &outs, int width)
{
    for(SrcIter row= sul; row.y<slr.y; row.y++)
    {
        for(SrcIter it= row; it.x<slr.x; it.x++)
            outs << std::setw(width) << (int)sa(it);
        outs << std::endl;
    }
}

template <class SrcIter, class SrcAcc>
void debugImage(SrcIter sul, SrcIter slr, SrcAcc sa,
                std::ostream &outs = std::cerr, int width= 1)
{
    if(sul != slr)
        debugImageImpl(sul, slr, sa, *sul, outs, width);
    else
        outs << "(empty image range)\n";
}

template <class SrcIter, class SrcAcc>
void debugImage(triple<SrcIter, SrcIter, SrcAcc> src,
                std::ostream &outs = std::cerr, int width= 1)
{
    debugImage(src.first, src.second, src.third, outs, width);
}

} // namespace vigra

#endif // VIGRA_DEBUGIMAGE_HXX
