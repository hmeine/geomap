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
