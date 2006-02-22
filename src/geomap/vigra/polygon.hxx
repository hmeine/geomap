#ifndef VIGRA_POLYGON_HXX
#define VIGRA_POLYGON_HXX

#include <vigra/diff2d.hxx>
#include <vector>
#include <algorithm>
#include "box.hxx"

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

/********************************************************************/
/*                                                                  */
/*             python stuff (should go into a subclass)             */
/*                                                                  */
/********************************************************************/

    PointArray(::boost::python::list l)
    : points_(::boost::python::len(l))
    {
        for(unsigned int i = 0; i < points_.size(); ++i)
        {
            points_[i] = ::boost::python::extract<POINT>(l[i]);
        }
    }

/********************************************************************/

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

    PointArray operator+(const Point &offset) const
    {
        PointArray result(points_.size());
        for(unsigned int i = 0; i < points_.size(); ++i)
            result[i] = points_[i] + offset;
        return result;
    }

    PointArray<Point2D> roundToInteger() const
    {
        PointArray<Point2D> result(points_.size());
        for(unsigned int i = 0; i < points_.size(); ++i)
        {
            result[i].x = (int)floor(points_[i][0] + 0.5);
            result[i].y = (int)floor(points_[i][1] + 0.5);
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

        /**
         * Tests whether the given point lies within this polygon.
         * Requires that this polygon is closed.

         * The result of testing points which lie directly on the
         * polylines (or are incident with the support points) is
         * undefined.  (ATM, the implementation uses half-open
         * intervals, so points on the left/top border are included,
         * in contrast to the ones on the right/bottom.)
         */
    bool contains(const_reference point) const
    {
        vigra_precondition(points_[size()-1] == points_[0],
                           "Polygon::contains() requires polygon to be closed!");
        int result = 0;
        bool above = points_[0][1] < point[1];
        for(unsigned int i = 1; i < size(); ++i)
        {
            bool now = points_[i][1] < point[1];
            if(now != above)
            {
                typename Point::value_type intersectX =
                    points_[i-1][0] + 
                    (points_[i][0] - points_[i-1][0]) *
                    (point[1]      - points_[i-1][1]) /
                    (points_[i][1] - points_[i-1][1]);
                if(intersectX < point[0])
                    ++result;

                above = now;
            }
        }
        return result % 2;
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

//     void swap(PointArray &rhs)
//     {
//         Base::swap(rhs);
//         invalidateProperties();
//     }

    void reverse()
    {
        Base::reverse();
        if(partialAreaValid_)
            partialArea_ = -partialArea_;
    }
    
    POINT nearestPoint(const_reference p) const;

  protected:
    mutable double length_;
    mutable bool lengthValid_;
    mutable double partialArea_;
    mutable bool partialAreaValid_;
};

template <class POINT>
POINT Polygon<POINT>::nearestPoint(const_reference p) const
{
    double dist = NumericTraits<double>::max();
    POINT r;
    for(unsigned int k=1; k<this->size(); ++k)
    {
        POINT dp = (*this)[k] - (*this)[k-1];
        POINT dc = p - (*this)[k-1];
        double t = dot(dp, dc);
        if(t != 0.0)
            t /= squaredNorm(dp);
        if(t > 1.0)
        {
            double d = norm((*this)[k]-p);
            if (d < dist)
            {
                dist = d;
                r = (*this)[k];
            }
        }
        else if(t < 0.0)
        {
            double d = norm((*this)[k-1]-p);
            if (d < dist)
            {
                dist = d;
                r = (*this)[k-1];
            }
        }
        else
        {
            POINT pp = (*this)[k-1] + t*dp;
            double d = norm(pp-p);
            if (d < dist)
            {
                dist = d;
                r = pp;
            }
        }
    }
    return r;
}

/********************************************************************/

// assumes that Point is a TinyVector
template<class POINT>
class BBoxPolygon : public Polygon<POINT>
{
  public:
    typedef Polygon<POINT> Base;

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

    typedef Box<typename Point::value_type, 2> BoundingBox;

    BBoxPolygon(typename Base::Base points)
    : Base(points),
      boundingBoxValid_(false)
    {}

    BBoxPolygon(Base poly)
    : Base(poly),
      boundingBoxValid_(false)
    {}

    BBoxPolygon()
    : boundingBoxValid_(false)
    {}

    BBoxPolygon(size_type n)
    : Base(n),
      boundingBoxValid_(false)
    {}

    template <class InputIterator>
    BBoxPolygon(InputIterator b, InputIterator e)
    : Base(b, e),
      boundingBoxValid_(false)
    {}

    bool contains(const_reference point) const
    {
        if(!boundingBox().contains(point))
            return false;
        return Base::contains(point);
    }

    BoundingBox boundingBox() const
    {
        if(!boundingBoxValid_)
        {
            boundingBox_ = BoundingBox();
            for(unsigned int i = 1; i < points_.size(); ++i)
                boundingBox_ |= points_[i];
            boundingBoxValid_ = true;
        }
        return boundingBox_;
    }

    void push_back(const_reference v)
    {
        if(boundingBoxValid_)
            boundingBox_ |= v;
        Base::push_back(v);
    }

    void extend(const BBoxPolygon &other)
    {
        if(boundingBoxValid_)
            boundingBox_ |= other.boundingBox();
        Base::extend(other);
    }

    void setPoint(unsigned int pos, const_reference x)
    {
        if(boundingBoxValid_)
        {
            if((x[0] < points_[pos][0]) &&
               (points_[pos][0] == boundingBox_.end()[0]) ||
               (x[0] > points_[pos][0]) &&
               (points_[pos][0] == boundingBox_.begin()[0]) ||
               (x[1] < points_[pos][1]) &&
               (points_[pos][1] == boundingBox_.end()[1]) ||
               (x[1] > points_[pos][1]) &&
               (points_[pos][1] == boundingBox_.begin()[1]))
                boundingBoxValid_ = false;
        }
        Base::setPoint(pos, x);
    }

    void erase(iterator pos)
    {
        if(boundingBoxValid_ && (
               (*pos)[0] == boundingBox_.begin()[0] ||
               (*pos)[0] == boundingBox_.end()[0] ||
               (*pos)[1] == boundingBox_.begin()[1] ||
               (*pos)[1] == boundingBox_.end()[1]))
            boundingBoxValid_ = false;
        Base::erase(pos);
    }

    iterator insert(iterator pos, const_reference x)
    {
        if(boundingBoxValid_)
            boundingBox_ |= x;
        return Base::insert(pos, x);
    }

    BBoxPolygon split(unsigned int pos)
    {
        BBoxPolygon result;
        Base base(Base::split(pos));
        static_cast<Base &>(result).swap(base);
        result.boundingBoxValid_ = boundingBoxValid_ = false;
        return result;
    }

    void swap(BBoxPolygon &rhs)
    {
        Base::swap(rhs);
        std::swap(boundingBox_, rhs.boundingBox_);
        std::swap(boundingBoxValid_, rhs.boundingBoxValid_);
    }

  protected:
    mutable BoundingBox boundingBox_;
    mutable bool boundingBoxValid_;
};

/********************************************************************/

template<bool useMaxStep, class PointIterator, class TargetArray>
void simplifyPolygonHelper(
    const PointIterator &polyBegin, const PointIterator &polyEnd,
    TargetArray &simple, double epsilon,
    double maxStep = vigra::NumericTraits<double>::max())
{
    if(polyEnd - polyBegin < 3)
        return; // no splitpoint necessary / possible

    PointIterator splitPos(polyEnd), lastPoint(polyEnd);
    --lastPoint;
    double maxDist = useMaxStep ? 0.0 : epsilon;

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

    if(useMaxStep)
    {
        if((maxDist <= epsilon) && (straight.magnitude() <= maxStep))
            return;
    }

    if(splitPos != polyEnd)
    {
        simplifyPolygonHelper<useMaxStep>(
            polyBegin, splitPos + 1, simple, epsilon, maxStep);
        simple.push_back(*splitPos);
        simplifyPolygonHelper<useMaxStep>(
            splitPos, polyEnd, simple, epsilon, maxStep);
    }
}

template<class PointArray>
void simplifyPolygon(
    const PointArray &poly, PointArray &simple, double epsilon)
{
    simple.push_back(poly[0]);
    simplifyPolygonHelper<false>(poly.begin(), poly.end(), simple, epsilon);
    simple.push_back(poly[poly.size()-1]);
}

template<class PointArray>
void simplifyPolygon(
    const PointArray &poly, PointArray &simple, double epsilon, double maxStep)
{
    simple.push_back(poly[0]);
    simplifyPolygonHelper<true>(poly.begin(), poly.end(), simple, epsilon, maxStep);
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
    : begin((int)floor(b+0.5)), direction(d), end((int)(ceil(e + 0.5)))
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

    int startIndex_;
    std::vector<Scanline> scanlines_;

    Scanlines(int startIndex, unsigned int count)
    : startIndex_(startIndex),
      scanlines_(count)
    {}

    void append(int line, const ScanlineSegment &seg)
    {
        int index = line - startIndex_;
        if(index >= 0 && index < (int)scanlines_.size())
            scanlines_[index].push_back(seg);
    }

    int startIndex() const
    {
        return startIndex_;
    }

    int endIndex() const
    {
        return scanlines_.size() - startIndex_;
    }

    Scanline &operator[](unsigned int index)
    {
        return scanlines_[index - startIndex_];
    }

    const Scanline &operator[](unsigned int index) const
    {
        return scanlines_[index - startIndex_];
    }

    unsigned int size() const
    {
        return scanlines_.size();
    }

    void normalize()
    {
        for(unsigned int i = 0; i < size(); ++i)
        {
            Scanlines::Scanline &scanline(scanlines_[i]);

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
    int startIndex = 0)
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
            if(prevKLine & 1)
            {
                // ensure that (prevLine + 0.5*step) will be a correct intersectY:
                prevLine = (prevKLine-step) / 2;
            }

            for(; prevKLine != kLine; prevKLine += step)
            {
                // we are leaving a scanline (actually, possibly only
                // an inter-pixel line if prevKLine is odd)
                typename Point::value_type
                    intersectX = prevPoint[0] + (x - prevPoint[0])
                    * ((prevLine + 0.5*step) - prevPoint[1]) / (y - prevPoint[1]);

                if(!(prevKLine & 1)) // even Khalimsky coordinate?
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
                else // coming from inter-pixel line
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

template<class DestIterator, class SizeType, class DestAccessor>
unsigned int fillScannedPoly(
    const Scanlines &scanlines,
    typename DestAccessor::value_type value,
    DestIterator dul, SizeType ds, DestAccessor a)
{
    bool clean = true;
    unsigned int pixelCount = 0;

    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(ds[1], scanlines.endIndex());

    for(DestIterator row(dul + y); y < endY; ++y, ++row)
    {
        int inside = 0;
        int x = 0;
        typename DestIterator::next_type it(row.begin());
        const Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
            {
                begin = 0;
                if(end < 0)
                    end = 0;
            }
            if(end > ds[0])
            {
                end = ds[0];
                if(begin > ds[0])
                    begin = ds[0];
            }

            if(inside > 0)
            {
                while(x < begin)
                {
                    a.set(value, it);
                    ++pixelCount;
                    ++x, ++it;
                }
            }
            it += end - x;
            x = end;
            inside += scanline[j].direction;
        }
        if(inside)
            clean = false;
    }

    vigra_postcondition(clean, "error in polygon scanlines (not closed?)");
    return pixelCount;
}

template<class DestIterator, class SizeType, class DestAccessor>
unsigned int drawScannedPoly(
    const Scanlines &scanlines,
    typename DestAccessor::value_type value,
    DestIterator dul, SizeType ds, DestAccessor a)
{
    unsigned int pixelCount = 0;

    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(ds[1], scanlines.endIndex());

    for(DestIterator row(dul + y); y < endY; ++y, ++row)
    {
        int x = 0;
        typename DestIterator::next_type it(row.begin());
        const Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > ds[0])
                end = ds[0];

            it += begin - x;
            x = begin;
            while(x < end)
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
