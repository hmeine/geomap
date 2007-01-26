#ifndef VIGRA_DSL_HXX
#define VIGRA_DSL_HXX

#include <vigra/rational.hxx>

namespace vigra {

namespace DSL {

    enum LeaningType {
        CenterLine,
        LowerLeaningLine,
        UpperLeaningLine
    };

};

template<class INTEGER = int, bool EIGHT_CONNECTED = true>
class DigitalStraightLine
{
  public:
    typedef INTEGER Integer;
    typedef Rational<Integer> RationalType;

    bool eightConnected() const
    {
        return EIGHT_CONNECTED;
    }

    DigitalStraightLine(int a, int b, int pos)
    : a_(a), b_(b), pos_(pos)
    {
    }

    Integer a() const { return a_; }
    Integer b() const { return b_; }
    Integer pos() const { return pos_; }

    void setA(Integer a)
    {
        a_ = a;
    }

    void setB(Integer b)
    {
        b_ = b;
    }

    void setPos(Integer pos)
    {
        pos_ = pos;
    }

    bool contains(Integer x, Integer y) const
    {
        Integer v(operator()(x, y) - pos_);
        return (0 <= v) && (v < width());
    }

    Integer width() const
    {
        if(eightConnected())
            return std::max(abs(a_), abs(b_));
        else
            return abs(a_) + abs(b_);
    }

    RationalType slope() const
    {
        return RationalType(a_, b_);
    }

    RationalType axisIntercept(DSL::LeaningType leaningType = DSL::CenterLine) const
    {
        RationalType result(pos_);
        if(leaningType == DSL::CenterLine)
            result += RationalType(width()-1, 2);
        else if(leaningType == DSL::LowerLeaningLine)
            result += width()-1;
        return -result / b_;
    }

    bool addPoint(Integer x, Integer y)
    {
        vigra_precondition(
            eightConnected() && (b_ >= a_) && (a_ >= 0),
            "addPoint() works only for 8-connected lines in the first octant!");

        Integer v(operator()(x, y) - pos_);
        if((0 <= v) && (v < b_))
            return true; // point is already within DSL

        bool above = true;
        if(v != -1)
        {
            if(v == b_)
                above = false;
            else
                // point cannot be added to DSL
                return false;
        }

        bool increaseSlope = above;
        Integer pos = pos_;
        vigra_precondition((x < 0) == (y < 0), "addPoint: invalid point given!");
        if(x < 0)
        {
            increaseSlope = !increaseSlope;
            pos = 1-b_-pos_; // temporarily mirror line at origin
        }

        Integer k = 0, divPos = (increaseSlope ? pos : pos + width() - 1);
        for(; k < b_; ++k)
        {
            if((a_*k-divPos) % b_ == 0)
                break;
        }
        Integer l = (a_*k-divPos) / b_;

        a_ = abs(y) - l;
        b_ = abs(x) - k;

        if(above)
            // ensure new point is on lower leaning line:
            pos_ = operator()(x, y);
        else
            // ensure new point is on upper leaning line:
            pos_ = operator()(x, y) - b_ + 1;

        vigra_postcondition(
            contains(x, y), "addPoint() should lead to contains()!");

        return true;
    }

    DigitalStraightLine<Integer, false>
    convertToFourConnected() const
    {
        if(eightConnected())
            return DigitalStraightLine<Integer, false>(
                a_, b_ - a_, pos_);
        else
            return DigitalStraightLine<Integer, false>(
                a_, b_, pos_);
    }

        // Ulli doesn't like this being public
        // (admittedly, it's not too intuitive):
    Integer operator()(Integer x, Integer y) const
    {
        return a_*x - b_*y;
    }

  protected:
    Integer a_, b_, pos_;
    bool eightConnected_;
};

} // namespace vigra

#endif // VIGRA_DSL_HXX
