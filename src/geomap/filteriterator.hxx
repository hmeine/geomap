/************************************************************************/
/*                                                                      */
/*          Copyright 1998-2002 by Ullrich Koethe, Hans Meine           */
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

namespace vigra {

template<class Iterator, class Predicate>
class FilterIterator
{
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

    FilterIterator()
    {}

    FilterIterator(Iterator begin, Iterator end,
                   Predicate predicate = Predicate())
        : iter_(begin), end_(end), predicate_(predicate)
    {
        if((iter_ != end_) && !predicate_(*iter_))
            operator++();
    }

    FilterIterator & operator++()
    {
        ++iter_;
        while((iter_ != end_) && !predicate_(*iter_))
            ++iter_;
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
        return iter_.operator->();
    }
};

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

#endif // VIGRA_FILTERITERATOR_HXX
