#ifndef CELLIMAGE_HXX
#define CELLIMAGE_HXX

#include "celltypes.hxx"
#include <vigra/basicimage.hxx>
#include "pixelneighborhood.hxx"

namespace vigra {

namespace CellImage {

// -------------------------------------------------------------------
//                          CellPixel, CellImage
// -------------------------------------------------------------------
struct CellPixel
{
    typedef unsigned int LabelType;

private:
    CellType type_;
    LabelType label_;

public:
    CellPixel() {}
    CellPixel(CellType type, int label = 0)
        : type_(type), label_(label)
    {}

    inline CellType type() const { return type_; }
    inline void setType(CellType type) { type_ = type; }

    inline LabelType label() const { return label_; }
    inline void setLabel(LabelType label) { label_= label; }
    inline void setLabel(LabelType label, CellType) { label_= label; }

    bool operator==(CellPixel const & rhs) const
        { return label_ == rhs.label_ && type_ == rhs.type_; }

    bool operator!=(CellPixel const & rhs) const
        { return label_ != rhs.label_ || type_ != rhs.type_; }
};

typedef BasicImage<CellPixel> CellImage;

typedef vigra::NeighborhoodCirculator<CellImage::Iterator, EightNeighborCode>
    CellImageEightCirculator;

// -------------------------------------------------------------------
//                     CellPixel/CellImage Accessors
// -------------------------------------------------------------------
struct CellImageTypeAccessor
{
    typedef CellType value_type;

    template<class Iterator>
    CellType operator()(const Iterator &it) const
    {
        return it->type();
    }

    template<class Iterator>
    void set(CellType type, const Iterator &it) const
    {
        it->setType(type);
    }
};

struct CellImageLabelAccessor
{
    typedef CellPixel::LabelType value_type;

    template<class Iterator>
    CellPixel::LabelType operator()(const Iterator &it) const
    {
        return it->label();
    }

    template<class Iterator>
    void set(CellPixel::LabelType label, const Iterator &it) const
    {
        it->setLabel(label);
    }
};

template<CellType type>
struct CellImageLabelWriter
{
    typedef CellPixel::LabelType value_type;

    template<class Iterator>
    void set(CellPixel::LabelType label, const Iterator &it) const
    {
        it->setLabel(label, type);
    }
};

// -------------------------------------------------------------------
//                             RelabelFunctor
// -------------------------------------------------------------------
template<class VALUETYPE>
struct RelabelFunctor
{
    typedef VALUETYPE value_type;
    typedef VALUETYPE argument_type;
    typedef VALUETYPE result_type;

    RelabelFunctor(VALUETYPE oldValue, VALUETYPE newValue)
        : oldValue_(oldValue),
          newValue_(newValue)
    {}
    
    VALUETYPE operator()(VALUETYPE value) const
    {
        return (value == oldValue) ? newValue : value;
    }
    
    VALUETYPE oldValue_, newValue_;
};

// -------------------------------------------------------------------
//                              inspectCell
// -------------------------------------------------------------------
// thinking about "typename IteratorTraits<EndIterator>::DefaultAccessor":
// is not needed since we're not implementing srcCellRange here, but
// algorithms.
// srcCellRange can not be implemented that easy, because most VIGRA
// functions expect an ImageIterator, not a std::iterator
template <class EndIterator, class Accessor, class Functor>
void inspectCell(EndIterator endIterator, Accessor a, Functor & f)
{
    for(; endIterator.inRange(); ++endIterator)
        f(a(endIterator));
}

template <class EndIterator, class Functor>
void inspectCell(EndIterator endIterator, Functor & f)
{
    for(; endIterator.inRange(); ++endIterator)
        f(*endIterator);
}

// -------------------------------------------------------------------
//                             transformCell
// -------------------------------------------------------------------
template <class SrcEndIterator, class SrcAccessor,
          class DestEndIterator, class DestAccessor, class Functor>
void transformCell(SrcEndIterator srcEndIterator, SrcAccessor sa,
                   DestEndIterator destEndIterator, DestAccessor da,
                   Functor const & f)
{
    for(; endIterator.inRange(); ++srcEndIterator, ++destEndIterator)
        da.set(f(sa(srcEndIterator)), destEndIterator);
}

template <class SrcEndIterator, class DestEndIterator, class Functor>
void transformCell(SrcEndIterator srcEndIterator,
                   DestEndIterator destEndIterator,
                   Functor const & f)
{
    for(; endIterator.inRange(); ++srcEndIterator, ++destEndIterator)
        *destEndIterator = f(*srcEndIterator);
}

} // namespace CellImage

// -------------------------------------------------------------------
//                                 Rect2D
// -------------------------------------------------------------------
class Rect2D
{
    Diff2D upperLeft_, lowerRight_;

public:
    /// construct a null rectangle
    Rect2D()
    {}

    // construct a rectangle representing the given range
    Rect2D(Diff2D const &upperLeft, Diff2D const &lowerRight)
        : upperLeft_(upperLeft), lowerRight_(lowerRight)
    {}

    // construct a rectangle representing the given range
    Rect2D(int left, int top, int right, int bottom)
        : upperLeft_(left, top), lowerRight_(right, bottom)
    {}

    Diff2D const & upperLeft() const
    {
        return upperLeft_;
    }

    Diff2D const & lowerRight() const
    {
        return lowerRight_;
    }

    int width() const
    {
        return lowerRight_.x - upperLeft_.x;
    }

    int height() const
    {
        return lowerRight_.y - upperLeft_.y;
    }

    bool operator==(Rect2D const &r) const
    {
        return (upperLeft_ == r.upperLeft_) && (lowerRight_ == r.lowerRight_);
    }

    bool operator!=(Rect2D const &r) const
    {
        return (upperLeft_ != r.upperLeft_) || (lowerRight_ != r.lowerRight_);
    }

    /** Return whether this rectangle is considered empty. It is
     * non-empty if both coordinates of the lower right corner are
     * greater than the corresponding coordinate of the upper left
     * corner. Uniting an empty rectangle with something will return
     * the bounding rectangle of the 'something', intersecting with an
     * empty rectangle will yield again an empty rectangle.
     */
    bool isEmpty() const
    {
        return ((lowerRight_.x <= upperLeft_.x) ||
                (lowerRight_.y <= upperLeft_.y));
    }

    /** Return whether this rectangle contains the given point. That
     * is, if the point lies within the valid range of an
     * ImageIterator walking from upperLeft() to lowerRight()
     * (excluding the latter).
     */
    bool contains(Diff2D const &p) const
    {
        return ((p.x >= upperLeft_.x) && (p.y >= upperLeft_.y) &&
                (p.x < lowerRight_.x) && (p.y < lowerRight_.y));
    }

    /** Return whether this rectangle contains the given
     * one. <tt>r1.contains(r2)</tt> returns the same as
     * <tt>r1 == (r1|r2)</tt> (but is of course more
     * efficient). That also means, a rectangle (even an empty one!)
     * contains() any empty rectangle.
     */
    bool contains(Rect2D const &r) const
    {
        return r.isEmpty() ||
            contains(r.upperLeft()) && contains(r.lowerRight());
    }

    /** Return whether this rectangle overlaps with the given
     * one. <tt>r1.intersects(r2)</tt> returns the same as
     * <tt>!(r1&r2).isEmpty()</tt> (but is of course much more
     * efficient).
     */
    bool intersects(Rect2D const &r) const
    {
        return ((r.upperLeft_.x < lowerRight_.x) &&
                (r.lowerRight_.x > upperLeft_.x) &&
                (r.upperLeft_.y < lowerRight_.y) &&
                (r.lowerRight_.y > upperLeft_.y))
            && !r.isEmpty();
    }

    /** Modifies this rectangle by including the given point. The
     * result is the bounding rectangle of the rectangle and the
     * point. If isEmpty returns true, the union will be a rectangle
     * containing only the given point.
     */
    Rect2D &operator|=(Diff2D const &p)
    {
        if(isEmpty())
        {
            upperLeft_ = p;
            lowerRight_ = p + Diff2D(1, 1);
        }
        else
        {
            if(p.x < upperLeft_.x)
                upperLeft_.x = p.x;
            if(p.y < upperLeft_.y)
                upperLeft_.y = p.y;
            if(p.x >= lowerRight_.x)
                lowerRight_.x = p.x + 1;
            if(p.y >= lowerRight_.y)
                lowerRight_.y = p.y + 1;
        }
        return *this;
    }

    /** Modifies this rectangle by uniting it with the given one. The
     * result is the bounding rectangle of both rectangles. If one of
     * the rectangles isEmpty(), the union will be the other one.
     */
    Rect2D &operator|=(Rect2D const &r)
    {
        if(r.isEmpty())
            return *this;
        if(isEmpty())
            return operator=(r);

        if(r.upperLeft_.x < upperLeft_.x)
            upperLeft_.x = r.upperLeft_.x;
        if(r.upperLeft_.y < upperLeft_.y)
            upperLeft_.y = r.upperLeft_.y;
        if(r.lowerRight_.x > lowerRight_.x)
            lowerRight_.x = r.lowerRight_.x;
        if(r.lowerRight_.y > lowerRight_.y)
            lowerRight_.y = r.lowerRight_.y;
        return *this;
    }

    /** Unites this rectangle with the given one. The result is the
     * bounding rectangle of both rectangles.
     */
    Rect2D operator|(Rect2D const &r) const
    {
        Rect2D result(*this);
        result |= r;
        return result;
    }

    /** Modifies this rectangle by intersecting it with the given
     * one. The result is the maximal rectangle contained in both
     * original ones. Intersecting with an empty rectangle will yield
     * again an empty rectangle.
     */
    Rect2D &operator&=(Rect2D const &r)
    {
        if(isEmpty())
            return *this;
        if(r.isEmpty())
            return operator=(r);

        if(r.upperLeft_.x > upperLeft_.x)
            upperLeft_.x = r.upperLeft_.x;
        if(r.upperLeft_.y > upperLeft_.y)
            upperLeft_.y = r.upperLeft_.y;
        if(r.lowerRight_.x < lowerRight_.x)
            lowerRight_.x = r.lowerRight_.x;
        if(r.lowerRight_.y < lowerRight_.y)
            lowerRight_.y = r.lowerRight_.y;
        return *this;
    }

    /** Intersects this rectangle with the given one. The result is
     * the maximal rectangle contained in both original ones.
     * Intersecting with an empty rectangle will yield again an empty
     * rectangle.
     */
    Rect2D operator&(Rect2D const &r) const
    {
        Rect2D result(*this);
        result &= r;
        return result;
    }
};

} // namespace vigra

#endif // CELLIMAGE_HXX
