#ifndef MYDEBUG_HXX
#define MYDEBUG_HXX

#include <iostream>
#include <iomanip>
#include "cellimage.hxx"

std::ostream &operator <<(std::ostream &s, vigra::Diff2D const &d)
{
    s << "(" << d.x << "/" << d.y << ")";
    return s;
}

std::ostream &operator <<(std::ostream &s, vigra::Point2D const &d)
{
    s << "(" << d.x << "," << d.y << ")";
    return s;
}

std::ostream &operator <<(std::ostream &s, vigra::Size2D const &d)
{
    s << "(" << d.x << "x" << d.y << ")";
    return s;
}

std::ostream &operator <<(std::ostream &s, vigra::Rect2D const &r)
{
    s << "[" << r.upperLeft() << " to " << r.lowerRight()
      << " = " << r.size() << "]";
    return s;
}

namespace vigra { namespace CellImage {

std::ostream &operator <<(std::ostream &s,
                          const CellPixel &p)
{
    std::streamsize width= s.width();
    switch(p.type())
    {
    case CellTypeRegion:
        s << "\033[1;34m" << std::setw(width) << p.label() << "\033[0m";
        break;
    case CellTypeLine:
        s << p.label();
        break;
    default:
        s << "\033[1;31m" << std::setw(width) << p.label() << "\033[0m";
    }
    return s;
}

} }

#endif // MYDEBUG_HXX
