#ifndef FILTERITERATOR_HXX
#define FILTERITERATOR_HXX

#include <iterator>
#include <vigra/iteratortraits.hxx>

namespace vigra {

template<class Iterator, class Predicate>
class FilterIterator
{
    Iterator iter_, end_;
    Predicate predicate_;

public:
        /** the iterator's value type
        */
    typedef typename Iterator::value_type value_type;

        /** the iterator's reference type (return type of <TT>*iter</TT>)
        */
    typedef typename Iterator::reference reference;

        /** the iterator's pointer type (return type of <TT>operator-></TT>)
        */
    typedef typename Iterator::pointer pointer;

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

} // namespace vigra

#endif // FILTERITERATOR_HXX
