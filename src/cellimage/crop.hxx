#ifndef CROP_HXX
#define CROP_HXX

namespace vigra {

    // FIXME: Size2D is a workaround until Iterator::operator- is adopted
template<class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
crop(triple<Iterator, Iterator, Accessor> t, Rect2D r)
{
    //vigra_precondition(Rect2D(Point2D(0, 0), Size2D((t.second - t.first).x, (t.second - t.first).y)).contains(r),
    //                   "crop called with rect which is out of (iterator) bounds!");
    return triple<Iterator, Iterator, Accessor>(t.first + r.upperLeft(),
                                                t.first + r.lowerRight(),
                                                t.third);
}

template<class Iterator, class Accessor>
inline pair<Iterator, Accessor>
crop(pair<Iterator, Accessor> t, Rect2D r)
{
    return pair<Iterator, Accessor>(t.first + r.upperLeft(), t.second);
}

} // namespace vigra

#endif // CROP_HXX
