#ifndef MYDEBUG_HXX
#define MYDEBUG_HXX

#include "foureightsegmentation.hxx"

std::ostream &operator <<(std::ostream &s, vigra::Diff2D const &d)
{
    s << "(" << d.x << "," << d.y << ")";
}

std::ostream &operator <<(std::ostream &s,
                          vigra::FourEightSegmentation::CellPixel const &p)
{
    switch(p.type())
    {
    case CellTypeRegion:
        s << "\033[1;34m" << p.label() << "\033[0m";
        break;
    case CellTypeLine:
        s << p.label();
        break;
    default:
        s << "\033[1;31m" << p.label() << "\033[0m";
    }
}

#endif // MYDEBUG_HXX
