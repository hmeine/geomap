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

#ifndef MYDEBUG_HXX
#define MYDEBUG_HXX

#include <iostream>
#include <iomanip>
#include "cellimage.hxx"
#include <vigra/diff2d.hxx>

namespace vigra {
namespace cellimage {

inline std::ostream &operator <<(std::ostream &s,
                                 const CellPixel &p)
{
    std::streamsize width= s.width();
    switch(p.type())
    {
    case CellTypeRegion:
        s << p.label();
        break;
    case CellTypeLine:
        s << "\033[1;31m" << std::setw(width) << p.label() << "\033[0m";
        break;
    default:
        s << "\033[1;34m" << std::setw(width) << p.label() << "\033[0m";
    }
    return s;
}

inline std::ostream &operator <<(std::ostream &s,
                                 const CellType &p)
{
    switch(p)
    {
    case CellTypeRegion:
        s << " ";
        break;
    case CellTypeLine:
        s << "\033[1;31m*\033[0m";
        break;
    default:
        s << "\033[1;34m*\033[0m";
    }
    return s;
}

}} // namespace vigra::cellimage

#endif // MYDEBUG_HXX
