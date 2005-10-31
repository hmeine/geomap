#ifndef VIGRA_POLYGON_HXX
#define VIGRA_POLYGON_HXX

#include <vigra/diff2d.hxx>
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
        return points_.insert(pos, x);
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

        Polygon::const_iterator otherBegin(other.begin());
        if(size())
        {
            if(*otherBegin == points_.back())
            {
                // don't copy first pixel
                ++otherBegin;
            }
            else
            {
                if(lengthValid_)
                    length_ +=
                        (other.points_.front() - points_.back()).magnitude();
                if(partialAreaValid_)
                    partialArea_ +=
                        (other.points_.front()[0]*points_.back()[1] -
                         other.points_.front()[1]*points_.back()[0]);
            }
        }
        if(lengthValid_)
            length_ += other.length();
        if(partialAreaValid_)
            partialArea_ += other.partialArea();
        points_.insert(points_.end(), otherBegin, other.end());
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
                partialArea_ += (v[0]*points_.back()[1] -
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
        return Base::insert(pos, x);
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
        simplifyPolygonHelper(polyBegin, splitPos + 1, simple, epsilon);
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

/********************************************************************/

struct ScanlineSegment
{
    int begin, direction, end;
        // enum { Up, Down, Touch } direction;

    ScanlineSegment(int b, int d, int e)
    : begin(b), direction(d), end(e)
    {}

    ScanlineSegment(double b, int d, double e)
    : begin((int)(b+0.5)), direction(d), end((int)(ceil(e + 0.5)))
    {}

    ScanlineSegment()
    {}

        // joins two successive (overlapping) segments,
        // assumes that other.begin >= this->begin
    void joinSuccessor(const ScanlineSegment &other)
    {
        if(other.end > end)
            end = other.end;
        direction += other.direction;
    }
};

struct ScanlineSegmentCompare
{
    bool operator()(const ScanlineSegment &a, const ScanlineSegment &b) const
    {
        return a.begin < b.begin;
    }
};

struct Scanlines
{
    typedef std::vector<ScanlineSegment> Scanline;
    typedef Scanline value_type;

    int startIndex;
    std::vector<Scanline> scanlines;

    Scanlines(int startIndex, unsigned int count)
    : startIndex(startIndex),
      scanlines(count)
    {}

    void append(int line, const ScanlineSegment &seg)
    {
        int index = line - startIndex;
        if(index >= 0 && index < scanlines.size())
            scanlines[line - startIndex].push_back(seg);
    }

    Scanline &operator[](unsigned int index)
    {
        return scanlines[index];
    }

    const Scanline &operator[](unsigned int index) const
    {
        return scanlines[index];
    }

    unsigned int size() const
    {
        return scanlines.size();
    }

    void normalize()
    {
        for(unsigned int i = 0; i < size(); ++i)
        {
            Scanlines::Scanline &scanline((*this)[i]);

            std::sort(scanline.begin(), scanline.end(),
                      ScanlineSegmentCompare());

            for(unsigned int j = 1; j < scanline.size(); )
            {
                if(scanline[j].begin <= scanline[j-1].end)
                {
                    scanline[j-1].joinSuccessor(scanline[j]);
                    scanline.erase(scanline.begin() + j);
                }
                else
                    ++j;
            }
        }
    }
};

template<class Point>
Scanlines *scanPoly(
    const PointArray<Point> &points,
    unsigned int scanLineCount,
    unsigned int startIndex = 0)
{
    Scanlines *result = new Scanlines(startIndex, scanLineCount);

    if(!points.size())
        return result;

    Point prevPoint(points[0]);
    typename Point::value_type
        s(prevPoint[0]),
        e(prevPoint[0]);

    int prevStep = 0; // how did we get into the current segment?
    int prevLine; // the scanline we are currently in
    int prevKLine; // prevLine as Khalimsky coordinate (odd: inter-pixel)

    { // init prevLine and prevKLine from first point
        typename Point::value_type
            py(prevPoint[1] + 0.5);
        prevLine = (int)floor(py);

        prevKLine = 2*prevLine;
        if(py == prevLine)
            --prevKLine;
    }

    for(unsigned int i = 1; i < points.size(); ++i)
    {
        // read current point coordinates
        typename Point::value_type
            x(points[i][0]), y(points[i][1]), ly(y + 0.5);

        // check whether we are in the right scanline
        int line = ((int)floor(ly)), kLine = 2*line;
        if(ly == line)
            --kLine;

        if(kLine != prevKLine)
        {
            int step = (kLine > prevKLine ? 1 : -1);
            for(; prevKLine != kLine; prevKLine += step)
            {
                typename Point::value_type intersectX =
                    prevPoint[0] + (x - prevPoint[0])
                    * (prevLine + 0.5*step - prevPoint[1])
                    / (y - prevPoint[1]);

                if(!(prevKLine & 1))
                {
                    // before leaving the scanline, include the intersection point
                    if(s > intersectX)
                        s = intersectX;
                    if(e < intersectX)
                        e = intersectX;

                    // save scanline
                    result->append(
                        prevKLine/2, ScanlineSegment(s, step + prevStep, e));
                }
                else
                {
                    // initialize next scanline segment
                    s = e = intersectX;

                    prevLine += step; // needed for intersectX calculation
                }
                prevStep = step;
            }
            prevLine = line;
        }

        // we are in the right scanline, include this point
        if(s > x)
            s = x;
        if(e < x)
            e = x;

        // remember point (using (x, y) is possibly faster than points[i])
        prevPoint[0] = x;
        prevPoint[1] = y;
    }

    if(!(prevKLine & 1))
        result->append(prevKLine/2, ScanlineSegment(s, prevStep, e));

    // normalize result (sort, merge overlapping)
    result->normalize();

    return result;
}

template<class DestIterator, class DestAccessor>
unsigned int fillScannedPoly(
    const Scanlines &scanlines,
    typename DestAccessor::value_type value,
    DestIterator dul, DestAccessor a)
{
    bool clean = true;
    unsigned int pixelCount = 0;

    DestIterator row(dul + scanlines.startIndex);
    for(unsigned int i = 0; i < scanlines.size(); ++i, ++row)
    {
        int inside = 0;
        int x = 0;
        typename DestIterator::next_type it(row.begin());
        const Scanlines::Scanline &scanline(scanlines[i]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            if(inside > 0)
            {
                while(x < scanline[j].begin)
                {
                    a.set(value, it);
                    ++pixelCount;
                    ++x, ++it;
                }
            }
            it += scanline[j].end - x;
            x = scanline[j].end;
            inside += scanline[j].direction;
        }
        if(inside)
            clean = false;
    }

    vigra_postcondition(clean, "error in polygon scanlines (not closed?)");
    return pixelCount;
}

template<class DestIterator, class DestAccessor>
unsigned int drawScannedPoly(
    const Scanlines &scanlines,
    typename DestAccessor::value_type value,
    DestIterator dul, DestAccessor a)
{
    unsigned int pixelCount = 0;

    DestIterator row(dul + scanlines.startIndex);
    for(unsigned int i = 0; i < scanlines.size(); ++i, ++row)
    {
        int x = 0;
        typename DestIterator::next_type it(row.begin());
        const Scanlines::Scanline &scanline(scanlines[i]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            it += scanline[j].begin - x;
            x = scanline[j].begin;
            while(x < scanline[j].end)
            {
                a.set(value, it);
                ++pixelCount;
                ++x, ++it;
            }
        }
    }

    return pixelCount;
}

} // namespace vigra

/********************************************************************/

namespace std
{

template<class T>
void swap(vigra::Polygon<T> &a, vigra::Polygon<T> &b)
{
    a.swap(b);
}

} // namespace std

#endif // VIGRA_POLYGON_HXX
