#ifndef VIGRA_CRACKCONNECTIONS_HXX
#define VIGRA_CRACKCONNECTIONS_HXX

namespace vigra {

// dest. image must be of src size + (1, 1)
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void
crackConnectionImage(SrcImageIterator srcUpperLeft,
                     SrcImageIterator srcLowerRight, SrcAccessor sa,
                     DestImageIterator destUpperLeft, DestAccessor da)
{
    enum { CONN_RIGHT = 1, CONN_DOWN = 2 };

    // first row
    {
        // upper left corner
        SrcImageIterator src(srcUpperLeft);
        DestImageIterator dest(destUpperLeft);
        da.set((int)(CONN_RIGHT | CONN_DOWN), dest);
        
        // middle part of first row
        SrcImageIterator left(src);
        ++src.x, ++dest.x;
        for(; src.x < srcLowerRight.x; ++src.x, ++dest.x, ++left.x)
        {
            int conn = CONN_RIGHT;
            if(sa(left) != sa(src))
                conn |= CONN_DOWN;
            da.set(conn, dest);
        }
        
        // upper right corner
        da.set((int)CONN_DOWN, dest);
    }

    // middle rows
    ++srcUpperLeft.y, ++destUpperLeft.y;
    for(; srcUpperLeft.y < srcLowerRight.y; ++srcUpperLeft.y, ++destUpperLeft.y)
    {
        SrcImageIterator src(srcUpperLeft);
        DestImageIterator dest(destUpperLeft);
        SrcImageIterator up(src); --up.y;

        // leftmost column
        {
            int conn = CONN_DOWN;
            if(sa(up) != sa(src))
                conn |= CONN_RIGHT;
            da.set(conn, dest);
        }

        // middle part
        SrcImageIterator left(src);
        ++src.x, ++dest.x, ++up.x;
        for(; src.x < srcLowerRight.x; ++src.x, ++dest.x, ++left.x, ++up.x)
        {
            int conn = 0;
            if(sa(up) != sa(src))
                conn |= CONN_RIGHT;
            if(sa(left) != sa(src))
                conn |= CONN_DOWN;
            da.set(conn, dest);
        }

        // rightmost column
        da.set((int)CONN_DOWN, dest);
    }

    // last row
    DestImageIterator dest(destUpperLeft);
    SrcImageIterator src(srcUpperLeft);
    for(; src.x < srcLowerRight.x; ++src.x, ++dest.x)
        da.set((int)CONN_RIGHT, dest);
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline
void
crackConnectionImage(
    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
    pair<DestImageIterator, DestAccessor> dest)
{
    crackConnectionImage(src.first, src.second, src.third,
                         dest.first, dest.second);
}

} // namespace vigra

#endif // VIGRA_CRACKCONNECTIONS_HXX
