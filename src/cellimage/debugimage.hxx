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
