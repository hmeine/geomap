#ifndef MYDEBUG_HXX
#define MYDEBUG_HXX

#include <iostream>
#include <iomanip>
#include "cellimage.hxx"
#include <vigra/rect2d.hxx>

inline std::ostream &operator <<(std::ostream &s, vigra::Diff2D const &d)
{
    s << "(" << d.x << "/" << d.y << ")";
    return s;
}

inline std::ostream &operator <<(std::ostream &s, vigra::Point2D const &d)
{
    s << "(" << d.x << "," << d.y << ")";
    return s;
}

inline std::ostream &operator <<(std::ostream &s, vigra::Size2D const &d)
{
    s << "(" << d.x << "x" << d.y << ")";
    return s;
}

inline std::ostream &operator <<(std::ostream &s, vigra::Rect2D const &r)
{
    s << "[" << r.upperLeft() << " to " << r.lowerRight()
      << " = " << r.size() << "]";
    return s;
}

inline std::ostream &operator <<(std::ostream &s,
                                 const vigra::cellimage::CellPixel &p)
{
    std::streamsize width= s.width();
    switch(p.type())
    {
    case vigra::cellimage::CellTypeRegion:
        s << p.label();
        break;
    case vigra::cellimage::CellTypeLine:
        s << "\033[1;31m" << std::setw(width) << p.label() << "\033[0m";
        break;
    default:
        s << "\033[1;34m" << std::setw(width) << p.label() << "\033[0m";
    }
    return s;
}

inline std::ostream &operator <<(std::ostream &s,
                                 const vigra::cellimage::CellType &p)
{
    switch(p)
    {
    case vigra::cellimage::CellTypeRegion:
        s << " ";
        break;
    case vigra::cellimage::CellTypeLine:
        s << "\033[1;31m*\033[0m";
        break;
    default:
        s << "\033[1;34m*\033[0m";
    }
    return s;
}

#endif // MYDEBUG_HXX
