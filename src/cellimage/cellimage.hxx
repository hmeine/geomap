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

} // namespace vigra

#endif // CELLIMAGE_HXX
