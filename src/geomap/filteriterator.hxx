/************************************************************************/
/*                                                                      */
/*          Copyright 1998-2002 by Hans Meine, Ullrich Koethe           */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_FILTERITERATOR_HXX
#define VIGRA_FILTERITERATOR_HXX

#include <iterator>
#include <vigra/imageiterator.hxx>

namespace vigra {

// -------------------------------------------------------------------
//                             FilterIterator
// -------------------------------------------------------------------
template<class Iterator, class Predicate>
class FilterIterator
{
  protected:
    Iterator iter_, end_;
    Predicate predicate_;

public:
        /** the iterator's value type
        */
    typedef typename std::iterator_traits<Iterator>::value_type value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename std::iterator_traits<Iterator>::reference reference;

        /** the iterator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename std::iterator_traits<Iterator>::pointer pointer;

        /** the iterator tag (forward_iterator_tag)
        */
    typedef std::forward_iterator_tag iterator_category;

    typedef Iterator adaptee;

    FilterIterator(Iterator begin, Iterator end,
                   Predicate predicate = Predicate())
        : iter_(begin), end_(end), predicate_(predicate)
    {
        checkPredicate();
    }

    FilterIterator & operator++()
    {
        ++iter_;
        checkPredicate();
        return *this;
    }

    FilterIterator operator++(int)
    {
        FilterIterator ret(*this);
        operator++();
        return ret;
    }

    /**
     * the opposite of inRange(); true if this iterator is behind the
     * range and should not be dereferenced any more
     */
	bool atEnd() const
	{
		return iter_ == end_;
	}

    /**
     * the opposite of atEnd(); true if this iterator is dereferencable
     */
	bool inRange() const
	{
		return iter_ != end_;
	}

	bool operator==(FilterIterator const &other) const
    {
        return iter_ == other.iter_;
    }

    bool operator!=(FilterIterator const &other) const
    {
        return iter_ != other.iter_;
    }

	bool operator==(Iterator const &other) const
    {
        return iter_ == other;
    }

    bool operator!=(Iterator const &other) const
    {
        return iter_ != other;
    }

    reference operator*() const
    {
        return *iter_;
    }

    pointer operator->() const
    {
        return &(operator*());
    }

    inline void checkPredicate()
    {
        while((iter_ != end_) && !predicate_(*iter_))
            ++iter_;
    }
};

/**
 * More "safe" version of FilterIterator: Checks predicate also before
 * dereferencing or when checking atEnd/inRange.  The downside is that
 * these operations are not const anymore, since they may move the
 * iterator (conceptually, this is wrong so I should probably have
 * made iter_ mutable).
 */
template<class Iterator, class Predicate>
class SafeFilterIterator
: public FilterIterator<Iterator, Predicate>
{
public:
    typedef FilterIterator<Iterator, Predicate> Base;

    SafeFilterIterator(Iterator begin, Iterator end,
                       Predicate predicate = Predicate())
    : Base(begin, end, predicate)
    {}

    typename Base::reference operator*()
    {
        checkPredicate();
        return *iter_;
    }

    typename Base::pointer operator->()
    {
        checkPredicate();
        return &(operator*());
    }

	bool atEnd()
	{
        checkPredicate();
		return Base::atEnd();
	}

	bool inRange()
	{
        checkPredicate();
		return Base::inRange();
	}
};

/** \addtogroup ImageIteratorAdapters
 */
//@{

// -------------------------------------------------------------------
//                          ImageFilterIterator
// -------------------------------------------------------------------

/** \brief Iterator adapter which filters with a given predicate.

    The ImageFilterIterator is a template class which calls a
    predicate functor (first template argument) on each iterator while
    proceeding and only stops if it returns true. It is a forward
    iterator and has to be given the end of the range to stop there
    without calling the predicate on invalid iterator
    positions. Because it knows the end, it fulfills the concept of an
    EndIterator, which means you can call atEnd() or inRange() on it.

    It can be used to scan a labelled region in the image for example:

    \code
    IImage labelImage = ...;

    struct FilterLabel : public std::unary_function<int, bool>
    {
        int label_;

        FilterLabel(int label): label_(label) {}

        bool operator()(int v) const
        {
            return v == label_;
        }
    };

    typedef vigra::ImageFilterIterator<FilterLabel, vigra::IImage::const_traverser>
        LabelFilterIterator

    for(LabelFilterIterator fit(srcImageRange(labelImage)); fit.inRange(); fit++)
    {
        ... // collect some label properties
    }

    \endcode

    <b>\#include</b> "<a href="filteriterator_8hxx-source.html">vigra/filteriterator.hxx</a>"<br>
    Namespace: vigra
*/
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

        /** the underlying adapted ImageIterator
         */
    typedef ImageIterator adaptee;

        /** Default constructor
         */
    ImageFilterIterator()
    {}

        /** Construct an iterator filtering the given image range with
         * the given predicate. Note that the range may be empty, that
         * is, <tt>upperLeft == lowerRight</tt>. Then this iterator
         * can serve as the end iterator for a conventional iterator
         * range end check:
         *
         * \code
         * LabelFilterIterator fit(sul, slr);
         * LabelFilterIterator fitEnd(slr, slr);
         * std::fill(fit, fitEnd, newValue);
         * \endcode
         */
    ImageFilterIterator(ImageIterator upperLeft, ImageIterator lowerRight,
                        Accessor ac = Accessor(),
                        Predicate predicate = Predicate())
        : iter_(upperLeft), lowerRight_(lowerRight), ac_(ac),
          predicate_(predicate), width_(lowerRight.x - upperLeft.x)
    {
        if((iter_ != lowerRight_) && !predicate_(ac_(iter_)))
            operator++();
    }

        /** Construct an iterator filtering the given image range with
         * the given predicate. Variant using argument objects.
         */
    ImageFilterIterator(triple<ImageIterator, ImageIterator, Accessor> src,
                        Predicate predicate = Predicate())
        : iter_(src.first), lowerRight_(src.second), ac_(src.third),
          predicate_(predicate), width_(src.second.x - src.first.x)
    {
        if((iter_ != lowerRight_) && !predicate_(ac_(iter_)))
            operator++();
    }

        /** Move to the next pixel fulfilling the predicate or to the
         * end position given to the constructor (pre-increment).
         */
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

        /** Move to the next pixel fulfilling the predicate or to the
         * end position given to the constructor (post-increment).
         */
    ImageFilterIterator operator++(int)
    {
        ImageFilterIterator ret(*this);
        operator++();
        return ret;
    }

        /**
         * The opposite of inRange(); true if this iterator is behind the
         * range and should not be dereferenced any more.
         */
    bool atEnd() const
    {
        return iter_ == lowerRight_;
    }

        /**
         * The opposite of atEnd(); true if this iterator is dereferencable.
         */
    bool inRange() const
    {
        return iter_ != lowerRight_;
    }

        /** equality
         */
    bool operator==(ImageFilterIterator const &other) const
    {
        return iter_ == other.iter_;
    }

        /** inequality
         */
    bool operator!=(ImageFilterIterator const &other) const
    {
        return iter_ != other.iter_;
    }

        /** equality check with the adaptee type
         */
    bool operator==(ImageIterator const &other) const
    {
        return iter_ == other;
    }

        /** inequality check with the adaptee type
         */
    bool operator!=(ImageIterator const &other) const
    {
        return iter_ != other;
    }

        /** Access pixel at the current coordinate.
         */
    reference operator*() const
    {
        return *iter_;
    }

        /** Access member of pixel at the current coordinate.
         */
    pointer operator->() const
    {
        return iter_.operator->();
    }
};

} // namespace vigra

#endif // VIGRA_FILTERITERATOR_HXX
