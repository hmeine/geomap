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

    typedef typename InternalVector::value_type      value_type;
    typedef typename InternalVector::reference       reference;
    typedef typename InternalVector::const_reference const_reference;
    typedef typename InternalVector::pointer         pointer;
    typedef typename InternalVector::const_pointer   const_pointer;
    typedef typename InternalVector::iterator        iterator;
    typedef typename InternalVector::const_iterator  const_iterator;
    typedef typename InternalVector::size_type       size_type;
    typedef typename InternalVector::difference_type difference_type;

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

    void push_back(const_reference v)
    {
        points_.push_back(v);
    }

    void erase(iterator pos)
    {
        points_.erase(pos);
    }

    void erase(iterator from, iterator to)
    {
        points_.erase(from, to);
    }

    iterator insert(iterator pos, const_reference x)
    {
        points_.insert(pos, x);
    }

    template <class InputIterator>
    void insert(iterator pos, InputIterator f, InputIterator l)
    {
        points_.insert(pos, f, l);
    }

    void insert(iterator pos, size_type n, const_reference x)
    {
        points_.insert(pos, n, x);
    }

    size_type size() const { return points_.size(); }

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
        PointArray<Point2D> result;
        for(unsigned int i = 0; i < points_.size(); ++i)
            result[i] = Point2D((int)(points_[i][0] + 0.5),
                                (int)(points_[i][1] + 0.5));
        return result;
    }

  protected:
    InternalVector points_;
};

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

    const typename Iterator::value_type &next()
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

template<class POINT>
class Polygon : public PointArray<POINT>
{
  public:
    typedef PointArray<POINT> Base;

        // TODO: repeat all other typedefs :-(
    typedef typename Base::size_type size_type;

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

    typedef PointIter<typename Base::InternalVector::const_reverse_iterator>
        RevPointIter;
    RevPointIter __reviter__() const
    {
        return RevPointIter(points_.rbegin(), points_.rend());
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

    void reverse()
    {
        Base::reverse();
        if(partialAreaValid_)
            partialArea_ = -partialArea_;
    }

//     void append(typename Base::value_type p)
//     {
//         Base::app
//     }

  protected:
    mutable double length_;
    mutable bool lengthValid_;
    mutable double partialArea_;
    mutable bool partialAreaValid_;
};

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
