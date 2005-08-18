#ifndef VIGRA_POLYGON_HXX
#define VIGRA_POLYGON_HXX

#include <vector>
#include <algorithm>

namespace vigra {

template<class POINT>
class PointArray
{
  protected:
    typedef std::vector<POINT> InternalVector;

  public:
    typedef POINT Point;

    typedef typename InternalVector::value_type             value_type;
    typedef typename InternalVector::reference              reference;
    typedef typename InternalVector::const_reference        const_reference;
    typedef typename InternalVector::pointer                pointer;
    typedef typename InternalVector::const_pointer          const_pointer;
    typedef typename InternalVector::iterator               iterator;
    typedef typename InternalVector::const_iterator         const_iterator;
    typedef typename InternalVector::reverse_iterator       reverse_iterator;
    typedef typename InternalVector::const_reverse_iterator const_reverse_iterator;
    typedef typename InternalVector::size_type              size_type;
    typedef typename InternalVector::difference_type        difference_type;

    PointArray()
    {}

    PointArray(size_type n)
    : points_(n)
    {}

    template <class InputIterator>
    PointArray(InputIterator b, InputIterator e)
    : points_(b, e)
    {}

    const_reference operator[](difference_type index) const
    {
        return points_[index];
    }

    reference operator[](difference_type index)
    {
        return points_[index];
    }

    const_iterator begin() const
    {
        return points_.begin();
    }

    iterator begin()
    {
        return points_.begin();
    }

    const_iterator end() const
    {
        return points_.end();
    }

    iterator end()
    {
        return points_.end();
    }

    const_reverse_iterator rbegin() const
    {
        return points_.rbegin();
    }

    reverse_iterator rbegin()
    {
        return points_.rbegin();
    }

    const_reverse_iterator rend() const
    {
        return points_.rend();
    }

    reverse_iterator rend()
    {
        return points_.rend();
    }

    void push_back(const_reference v)
    {
        points_.push_back(v);
    }

    void erase(iterator pos)
    {
        points_.erase(pos);
    }

    iterator insert(iterator pos, const_reference x)
    {
        points_.insert(pos, x);
    }

    size_type size() const
    {
        return points_.size();
    }

    void swap(PointArray &rhs)
    {
        std::swap(points_, rhs.points_);
    }

    void reverse()
    {
        std::reverse(points_.begin(), points_.end());
    }

    PointArray operator*(double scale) const
    {
        PointArray result(points_.size());
        for(unsigned int i = 0; i < points_.size(); ++i)
            result[i] = points_[i] * scale;
        return result;
    }

    PointArray<Point2D> roundToInteger() const
    {
        PointArray<Point2D> result(points_.size());
        for(unsigned int i = 0; i < points_.size(); ++i)
        {
            result[i].x = (int)(points_[i][0] + 0.5);
            result[i].y = (int)(points_[i][1] + 0.5);
        }
        return result;
    }

  protected:
    InternalVector points_;
};

/********************************************************************/

template<class ITERATOR>
class PointIter
{
  public:
    typedef ITERATOR Iterator;

    PointIter(Iterator begin, Iterator end)
    : begin_(begin),
      end_(end)
    {}

    PointIter __iter__() const
    {
        return *this;
    }

    unsigned int __len__() const
    {
        return end_ - begin_;
    }

    typename Iterator::value_type next()
    {
        if(begin_ == end_)
        {
            PyErr_SetString(PyExc_StopIteration, "");
            boost::python::throw_error_already_set();
        }
        return *(begin_++);
    }

  protected:
    Iterator begin_, end_;
};

/********************************************************************/

template<class POINT>
class Polygon : public PointArray<POINT>
{
  public:
    typedef PointArray<POINT> Base;

    typedef typename Base::Point                  Point;
    typedef typename Base::value_type             value_type;
    typedef typename Base::reference              reference;
    typedef typename Base::const_reference        const_reference;
    typedef typename Base::pointer                pointer;
    typedef typename Base::const_pointer          const_pointer;
    typedef typename Base::iterator               iterator;
    typedef typename Base::const_iterator         const_iterator;
    typedef typename Base::reverse_iterator       reverse_iterator;
    typedef typename Base::const_reverse_iterator const_reverse_iterator;
    typedef typename Base::size_type              size_type;
    typedef typename Base::difference_type        difference_type;

/********************************************************************/
/*                                                                  */
/*             python stuff (should go into a subclass)             */
/*                                                                  */
/********************************************************************/

    Polygon(::boost::python::list l)
    : Base(::boost::python::len(l)),
      lengthValid_(false),
      partialAreaValid_(false)
    {
        for(unsigned int i = 0; i < points_.size(); ++i)
        {
            points_[i] = ::boost::python::extract<POINT>(l[i]);
        }
    }

    void extend(const Polygon &other)
    {
        if(!other.size())
            return;
        if(size())
        {
            if(lengthValid_)
                length_ += (other.points_.front() - points_.back()).magnitude();
            if(partialAreaValid_)
                partialAreaValid_ += (other.points_.front()[0]*points_.back()[1] -
                                      other.points_.front()[1]*points_.back()[0]);
        }
        if(lengthValid_)
            length_ += other.length();
        if(partialAreaValid_)
            partialArea_ += other.partialArea();
        points_.insert(points_.end(), other.begin(), other.end());
    }

/********************************************************************/

    Polygon(Base points)
    : Base(points),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    Polygon()
    : length_(0.0),
      lengthValid_(true),
      partialArea_(0.0),
      partialAreaValid_(true)
    {}

    Polygon(size_type n)
    : Base(n),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

        // Polygon(const Polygon&);

    template <class InputIterator>
    Polygon(InputIterator b, InputIterator e)
    : Base(b, e),
      lengthValid_(false),
      partialAreaValid_(false)
    {}

    void invalidateProperties()
    {
        lengthValid_ = false;
        partialAreaValid_ = false;
    }

    double length() const
    {
        if(!lengthValid_)
        {
            length_ = 0.0;
            for(unsigned int i = 1; i < points_.size(); ++i)
                length_ += (points_[i] - points_[i-1]).magnitude();
            lengthValid_ = true;
        }
        return length_;
    }

    double partialArea() const
    {
        if(!partialAreaValid_)
        {
            partialArea_ = 0.0;
            for(unsigned int i = 1; i < points_.size(); ++i)
                partialArea_ += (points_[i][0]*points_[i-1][1] -
                                 points_[i][1]*points_[i-1][0]);
            partialAreaValid_ = true;
        }
        return partialArea_;
    }

    void push_back(const_reference v)
    {
        if(size())
        {
            if(lengthValid_)
                length_ += (v - points_.back()).magnitude();
            if(partialAreaValid_)
                partialAreaValid_ += (v[0]*points_.back()[1] -
                                      v[1]*points_.back()[0]);
        }
        Base::push_back(v);
    }

    void setPoint(unsigned int pos, const_reference x)
    {
        if(lengthValid_)
        {
            if(pos > 0)
            {
                length_ += (x - points_[pos-1]).magnitude() -
                           (points_[pos] - points_[pos-1]).magnitude();
            }
            if(pos < size() - 1)
            {
                length_ += (x - points_[pos+1]).magnitude() -
                           (points_[pos] - points_[pos+1]).magnitude();
            }
        }
        partialAreaValid_ = false;
        points_[pos] = x;
    }

    void erase(iterator pos)
    {
        invalidateProperties();
        Base::erase(pos);
    }

    iterator insert(iterator pos, const_reference x)
    {
        if(lengthValid_)
        {
            if(pos > begin())
                length_ += (x - pos[-1]).magnitude();
            if(end() - pos >= 1)
            {
                length_ += (x - *pos).magnitude();
                if(pos > begin())
                    length_ -= (*pos - pos[-1]).magnitude();
            }
        }
        partialAreaValid_ = false;
        Base::insert(pos, x);
    }

    Polygon split(unsigned int pos)
    {
        if(pos == 0)
        {
            Polygon result(1);
            result[0] = points_[0];
            swap(result);
            return result;
        }

        Polygon result(begin() + pos, end());
        if(pos > size() / 3)
        {
             // heuristic: when splitting off only a "small part",
             // re-use existing information
            if(lengthValid_)
                length_ -= result.length();
            if(partialAreaValid_)
                partialArea_ -= result.partialArea();
        }
        else
            invalidateProperties();
        points_.erase(begin() + pos + 1, end());
        return result;
    }

    void swap(Polygon &rhs)
    {
        Base::swap(rhs);
        std::swap(length_, rhs.length_);
        std::swap(lengthValid_, rhs.lengthValid_);
        std::swap(partialArea_, rhs.partialArea_);
        std::swap(partialAreaValid_, rhs.partialAreaValid_);
    }

    void reverse()
    {
        Base::reverse();
        if(partialAreaValid_)
            partialArea_ = -partialArea_;
    }

  protected:
    mutable double length_;
    mutable bool lengthValid_;
    mutable double partialArea_;
    mutable bool partialAreaValid_;
};

/********************************************************************/

template<class PointIterator, class TargetArray>
void simplifyPolygonHelper(
    const PointIterator &polyBegin, const PointIterator &polyEnd,
    TargetArray &simple, double epsilon)
{
    if(polyEnd - polyBegin < 3)
        return; // no splitpoint necessary / possible

    PointIterator splitPos(polyEnd), lastPoint(polyEnd);
    --lastPoint;
    double maxDist = epsilon;

    // calculate normal of straight end point connection
    typename TargetArray::value_type
        straight(*lastPoint - *polyBegin),
        normal(straight[1], -straight[0]);

    // search splitpoint
    if(normal.magnitude() > 1e-6)
    {
        normal /= normal.magnitude();

        PointIterator it(polyBegin);
        for(++it; it != lastPoint; ++it)
        {
            double dist = fabs(dot(*it - *polyBegin, normal));
            if(dist > maxDist)
            {
                splitPos = it;
                maxDist = dist;
            }
        }
    }
    else
    {
        // start- and end-points identical?! -> look for most distant point
        PointIterator it(polyBegin);
        for(++it; it != lastPoint; ++it)
        {
            double dist = (*it - *polyBegin).magnitude();
            if(dist > maxDist)
            {
                splitPos = it;
                maxDist = dist;
            }
        }
    }

    if(splitPos != polyEnd)
    {
        simplifyPolygonHelper(polyBegin, splitPos, simple, epsilon);
        simple.push_back(*splitPos);
        simplifyPolygonHelper(splitPos, polyEnd, simple, epsilon);
    }
}

template<class PointArray>
void simplifyPolygon(const PointArray &poly, PointArray &simple, double epsilon)
{
    simple.push_back(poly[0]);
    simplifyPolygonHelper(poly.begin(), poly.end(), simple, epsilon);
    simple.push_back(poly[poly.size()-1]);
}

} // namespace vigra

namespace std
{

template<class T>
void swap(vigra::Polygon<T> &a, vigra::Polygon<T> &b)
{
    a.swap(b);
}

} // namespace std

#endif // VIGRA_POLYGON_HXX
