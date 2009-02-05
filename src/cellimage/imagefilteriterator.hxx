#ifndef IMAGEFILTERITERATOR_HXX
#define IMAGEFILTERITERATOR_HXX

#include <iterator>
#include <vigra/iteratortraits.hxx>

namespace vigra {

template<class Predicate, class ImageIterator,
         class Accessor = typename IteratorTraits<ImageIterator>::DefaultAccessor>
class ImageFilterIterator
{
    ImageIterator iter_, lowerRight_;
    Accessor ac_;
    Predicate predicate_;
    int width_;

public:
        /** the iterator's value type
        */
    typedef typename ImageIterator::value_type value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename ImageIterator::reference reference;

        /** the iterator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename ImageIterator::pointer pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    typedef ImageIterator adaptee;

    ImageFilterIterator()
    {}

    ImageFilterIterator(ImageIterator upperLeft, ImageIterator lowerRight,
                        Accessor ac = Accessor(),
                        Predicate predicate = Predicate())
        : iter_(upperLeft), lowerRight_(lowerRight), ac_(ac),
          predicate_(predicate), width_(lowerRight.x - upperLeft.x)
    {
        if((iter_ != lowerRight_) && !predicate_(ac_(iter_)))
            operator++();
    }

    ImageFilterIterator(triple<ImageIterator, ImageIterator, Accessor> src,
                        Predicate predicate = Predicate())
        : iter_(src.first), lowerRight_(src.second), ac_(src.third),
          predicate_(predicate), width_(src.second.x - src.first.x)
    {
        if((iter_ != lowerRight_) && !predicate_(ac_(iter_)))
            operator++();
    }

    ImageFilterIterator & operator++()
    {
        ++iter_.x;
        while((iter_.x != lowerRight_.x) && !predicate_(ac_(iter_)))
            ++iter_.x;

        if(iter_.x == lowerRight_.x)
        {
            iter_.x -= width_ + 1;
            ++iter_.y;

            if(iter_.y != lowerRight_.y)
                operator++();
            else
                iter_ = lowerRight_;
        }
        return *this;
    }

    ImageFilterIterator operator++(int)
    {
        ImageFilterIterator ret(*this);
        operator++();
        return ret;
    }

    /**
     * the opposite of inRange(); true if this iterator is behind the
     * range and should not be dereferenced any more
     */
	bool atEnd() const
	{
		return iter_ == lowerRight_;
	}

    /**
     * the opposite of atEnd(); true if this iterator is dereferencable
     */
	bool inRange() const
	{
		return iter_ != lowerRight_;
	}

	bool operator==(ImageFilterIterator const &other) const
    {
        return iter_ == other.iter_;
    }

    bool operator!=(ImageFilterIterator const &other) const
    {
        return iter_ != other.iter_;
    }

	bool operator==(ImageIterator const &other) const
    {
        return iter_ == other;
    }

    bool operator!=(ImageIterator const &other) const
    {
        return iter_ != other;
    }

    reference operator*() const
    {
        return *iter_;
    }

    pointer operator->() const
    {
        return iter_.operator->();
    }
};

} // namespace vigra

#endif // IMAGEFILTERITERATOR_HXX
