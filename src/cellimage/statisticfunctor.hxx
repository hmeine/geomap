#ifndef STATISTICFUNCTOR_HXX
#define STATISTICFUNCTOR_HXX

#include <vigra/numerictraits.hxx>

/**
 * Collects statistical information by means similar to a simple calculator.
 * That means: Count values, store sum and sum of squares, then calculate
 * other values from that. (Recall: Var(X)=Mean(X^2)-(Mean(X))^2)
 */
template<class ValueType>
class StatisticFunctor {
public:
	typedef typename vigra::NumericTraits<ValueType>::Promote PromoteType;
	typedef typename vigra::NumericTraits<ValueType>::RealPromote RealPromoteType;

	long count;
	PromoteType sum;
	RealPromoteType squareSum;

	RealPromoteType sqr(const ValueType &x) const { return x*x; }

	StatisticFunctor() { clear(); }
	void clear()
	{
		ValueType zero= vigra::NumericTraits<ValueType>::zero();
		init(0, zero, zero);
	}
	void init(long newCount, ValueType newAvg, ValueType newStdDev)
	{
		count = newCount;
		sum= count * newAvg;
		squareSum= sqr(newStdDev) + sqr(newAvg);
	}
	void operator()(ValueType value)
	{
		count++;
		sum = sum + value;
		squareSum = squareSum + sqr(value);
	}
	void operator()(const StatisticFunctor &other)
	{
		merge(other);
	}
	void merge(const StatisticFunctor &other)
	{
		count += other.count;
		sum += other.sum;
		squareSum += other.squareSum;
	}
	void remove(ValueType value)
	{
		count--;
		sum = sum - value;
		squareSum = squareSum - sqr(value);
	}
	RealPromoteType avg() const { return (RealPromoteType)(sum/((double)count)); }
	RealPromoteType var() const { return squareSum/((double)count) - sqr(avg()); }
	RealPromoteType stdDev() const { return sqrt(var()); }
};

#endif // STATISTICFUNCTOR_HXX
