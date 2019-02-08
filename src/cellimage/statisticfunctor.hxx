/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

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
