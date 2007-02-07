#ifndef STATISTICFUNCTOR_HXX
#define STATISTICFUNCTOR_HXX

#include <vigra/numerictraits.hxx>

/**
 * Collects statistical information by means similar to a simple calculator.
 * That means: Count values, store sum and sum of squares, then calculate
 * other values from that. (Recall: Var(X)=Mean(X^2)-(Mean(X))^2)
 */
template<class ValueType>
class StatisticFunctor
{
public:
	typedef typename vigra::NumericTraits<ValueType>::Promote PromoteType;
	typedef typename vigra::NumericTraits<ValueType>::RealPromote RealPromoteType;

    typedef ValueType argument_type;
    typedef PromoteType result_type; // fake type, for ArrayOfRegionStatistics

	inline PromoteType sqr(const ValueType &x) const
    {
        return x*x;
    }

	StatisticFunctor()
    {
        reset();
    }

	void reset()
	{
        count_     = 0;
        sum_       = vigra::NumericTraits<PromoteType>::zero();
        squareSum_ = vigra::NumericTraits<PromoteType>::zero();
	}

	void operator()(ValueType value)
	{
		++count_;
		sum_       += value;
		squareSum_ += sqr(value);
	}

	void operator()(const StatisticFunctor &other)
	{
		merge(other);
	}

	void merge(const StatisticFunctor &other)
	{
		count_     += other.count_;
		sum_       += other.sum_;
		squareSum_ += other.squareSum_;
	}

	void remove(ValueType value)
	{
		--count_;
		sum_       -= value;
		squareSum_ -= sqr(value);
	}

	RealPromoteType average() const
    {
        return (RealPromoteType)sum_/(double)count_;
    }

	RealPromoteType variance() const
    {
        return (RealPromoteType)squareSum_/(double)count_ - sqr(average());
    }

	RealPromoteType stdDev() const {
        return sqrt(variance());
    }

	long        count_;
	PromoteType sum_;
	PromoteType squareSum_;
};

#endif // STATISTICFUNCTOR_HXX
