#include "vigra/polygon.hxx"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <vigra/gaussians.hxx>
#include <vigra/pythonimage.hxx>
#include <vigra/linear_algebra.hxx>
#include <vigra/regression.hxx>
#include <cmath>
#include "exporthelpers.hxx"

using namespace vigra;
using namespace boost::python;

class PythonListBackInserter
{
    list & list_;
  public:

    PythonListBackInserter(list & l)
    : list_(l)
    {}

    template <class T>
    void push_back(T const & t)
    {
        list_.append(t);
    }
};

inline Point2D intPos(const Vector2 &p)
{
    return Point2D((int)floor(p[0]+0.5), (int)floor(p[1]+0.5));
}

template<class Box>
Box *createBoxFromRect2D(const Rect2D &r)
{
    return new Box(typename Box::Vector(r.left(), r.top()),
                   typename Box::Vector(r.right(), r.bottom()));
}

template<class Box>
inline Rect2D intPos_Box(const Box &b)
{
    return Rect2D(intPos(b.begin()), intPos(b.end())+Size2D(1,1));
}

template<class Box>
PyObject * Box__repr__(Box const & b)
{
    std::stringstream s;

    s << "<BoundingBox " << b.begin() << ".." << b.end() << ">";

    return PyString_FromString(s.str().c_str());
}

template<class Box>
struct BoxPickleSuite : pickle_suite
{
    static tuple getinitargs(Box const& b)
    {
        return make_tuple(b.begin(), b.end());
    }
};

// index started with zero (commented out) for better python iterations
// (allowing "for .. in .." or "enumerate(..)")
// TODO: a better fix would be an extra __iter__:
const Scanlines::value_type &
Scanlines__getitem__(Scanlines const & s, int i)
{
//     if(i < 0 || i >= s.size())
    if(i < s.startIndex() || i >= s.endIndex())
    {
        PyErr_SetString(PyExc_IndexError,
            "scanline index out of bounds.");
        throw_error_already_set();
    }
//     return s[i+s.startIndex()];
    return s[i];
}

template<class Polygon>
struct PolygonFromPython
{
    typedef Polygon Type;

    PolygonFromPython()
    {
        converter::registry::insert(
            &convertible, &construct, type_id<Type>());
    }

    static void* convertible(PyObject* obj)
    {
        if(!PyList_Check(obj))
            return NULL;
        unsigned int ll(PySequence_Fast_GET_SIZE(obj));
        for(unsigned int i = 0; i < ll; ++i)
        {
            if(!extract<Vector2>(
                PySequence_Fast_GET_ITEM(obj, i)).check())
                return NULL;
        }
        return obj;
    }

    static void construct(PyObject* obj, converter::rvalue_from_python_stage1_data* data)
    {
        Type * const storage = reinterpret_cast<Type *>((
            reinterpret_cast<converter::rvalue_from_python_storage<Type> *>(
                data))->storage.bytes);
        new (storage) Type();
        unsigned int ll(PySequence_Fast_GET_SIZE(obj));
        for(unsigned int i = 0; i < ll; ++i)
        {
            (*storage).push_back(extract<Vector2>(
                PySequence_Fast_GET_ITEM(obj, i))());
        }
        data->convertible = storage;
    }
};

// FIXME: to be usable with make_constructor, this must return a
// pointer type, but this is ugly:
template<class Polygon>
Polygon * createPolygon(boost::python::list l)
{
    return new Polygon(extract<Polygon>(l)());
}

template<class Polygon>
void
Polygon__setitem__(Polygon & p, int i, typename Polygon::value_type v)
{
    checkPythonIndex(i, p.size());
    p.setPoint(i, v);
}

template<class Polygon>
void insert(Polygon & p, int pos, typename Polygon::const_reference x)
{
    // this would be the original python semantics, but I wonder who
    // wants to insert at a non-existent position (e.g. index 42 for a
    // 5-element list)?
//     if(pos >= p.size())
//         p.push_back(x);
    checkPythonIndex(pos, p.size() + 1);
    p.insert(p.begin() + pos, x);
}

template<class Polygon>
void erase(Polygon & p, int pos)
{
    checkPythonIndex(pos, p.size() + 1);
    p.erase(p.begin() + pos);
}

template<class Polygon>
Polygon split(Polygon & p, int pos)
{
    checkPythonIndex(pos, p.size());
    return p.split(pos);
}

template<class Polygon>
std::string Polygon__repr__(Polygon const &polygon)
{
    std::stringstream s;
    s.precision(3);
    s << "<Polygon (" << polygon.size() << " points, length: " << polygon.length()
      << (polygon.closed() ? " px., closed)>" : " px.)>");
    return s.str();
}

template<class Iterator>
void defIter(const char *name)
{
    class_<Iterator>(name, no_init)
        .def("__len__", &Iterator::__len__)
        .def("__iter__", &Iterator::__iter__,
             return_internal_reference<>())
        .def("next", &Iterator::next)
    ;
}

template<class Array>
list arcLengthList(const Array &a)
{
    list result;
    PythonListBackInserter ins(result);
    a.arcLengthList(ins);
    return result;
}

template<class Array>
Array pyConvexHull(const Array &a)
{
    Array result;
    convexHull(a, result);
    return result;
}

template<class Array>
Array simplifyPolygon(const Array &a, double epsilon)
{
    Array result;
    simplifyPolygon(a, result, epsilon);
    return result;
}

template<class Array>
Array simplifyPolygon(const Array &a, double epsilon, double maxStep)
{
    Array result;
    simplifyPolygon(a, result, epsilon, maxStep);
    return result;
}

template<class Array>
Array simplifyPolygonDigitalLine(const Array &a, int connectivity)
{
    Array result;
    simplifyPolygonDigitalLine(a, result, connectivity);
    return result;
}

template<class Array>
Array resamplePolygon(const Array &a, double desiredPointDistance)
{
    Array result;
    resamplePolygon(a, result, desiredPointDistance);
    return result;
}

template<class Array>
Array resamplePolygonLinearInterpolation(const Array &a, double desiredPointDistance)
{
    Array result;
    resamplePolygonLinearInterpolation(a, result, desiredPointDistance);
    return result;
}

template<class Array>
Array resamplePolygonExponentialFilter(const Array &a, double scale, double desiredPointDistance)
{
    Array result;
    resamplePolygonExponentialFilter(a, result, scale, desiredPointDistance);
    return result;
}

template<class Array>
Array resamplePolygonGaussianFilter(const Array &a, double scale, double desiredPointDistance)
{
    Array result;
    resamplePolygonGaussianFilter(a, result, scale, desiredPointDistance);
    return result;
}

template<class Array>
Array polygonSplineControlPoints(const Array &a, int segmentCount)
{
    Array result;
    polygonSplineControlPoints(a, result, segmentCount);
    return result;
}

ScanlinesIter createScanlinesIter(const Scanlines &sl)
{
    return ScanlinesIter(sl);
}

unsigned int pyFillScannedPoly(
    const Scanlines &scanlines,
    PythonImage &targetV,
    GrayValue value)
    //const Pixel &value)
{
    PythonSingleBandImage target(targetV.subImage(0));
    return fillScannedPoly(scanlines, value,
                           target.traverser_begin(),
                           target.size(),
                           StandardValueAccessor<GrayValue>());
}

unsigned int pyDrawScannedPoly(
    const Scanlines &scanlines,
    PythonImage &targetV,
    float value)
{
    PythonSingleBandImage target(targetV.subImage(0));
    return drawScannedPoly(scanlines, value,
                           target.traverser_begin(),
                           target.size(),
                           StandardValueAccessor<GrayValue>());
}

// needed for signature conversion only (bbox poly vs. poly)
template<class Polygon>
Vector2 pyCentroid(const Polygon &p)
{
    return centroid(p);
}

template<class Array>
list curvatureList(const Array &p, int dx = 5, unsigned int skip = 0)
{
    if(p.size() < 2*dx + 2*skip)
    {
        PyErr_SetString(PyExc_ValueError,
            "curvatureList: polygon too small (less than 2*dx + 2*skip points).");
        throw_error_already_set();
    }

    double pos = 0.0;
    for(unsigned int i = 0; i < dx + skip - 1; ++i)
        pos += (p[i+1]-p[i]).magnitude();

    list result;

    for(unsigned int i = skip; i < p.size() - 2*dx - skip; ++i)
    {
        typename Array::value_type
            s1(p[i+  dx] - p[i]),
            s2(p[i+2*dx] - p[i+dx]),
            s3(p[i+2*dx] - p[i]);

        double
            a2 = s1.squaredMagnitude(),
            b2 = s2.squaredMagnitude(),
            c2 = s3.squaredMagnitude(),
            a = sqrt(a2),
            b = sqrt(b2);

        pos += (p[i+dx] - p[i+dx-1]).magnitude();

        static const double eps = 1e-8;
        if(a < eps || b < eps)
            continue;

        double gamma = (a2+b2-c2) / (2*a*b), delta;
        if(gamma > -1.0)
        {
            delta = M_PI - acos(gamma);
            if((s1[0]*s3[1] - s3[0]*s1[1]) < 0)
                result.append(make_tuple(pos, -2*delta / (a+b)));
            else
                result.append(make_tuple(pos,  2*delta / (a+b)));
        }
        else
            result.append(make_tuple(pos, 0.0));
    }

    return result;
}

struct ContinuousDirection
{
    double prevAngle, offset;

    ContinuousDirection(double startAngle = 0.0)
    : prevAngle(startAngle),
      offset(0.0)
    {}

    double operator()(double angle)
    {
        angle += offset;
        double diff = angle - prevAngle;
        if(fabs(diff) > M_PI)
        {
            double delta = -floor((diff + M_PI)/(M_PI*2.0))*M_PI*2.0;
            offset += delta;
            angle  += delta;
        }
        prevAngle = angle;
        return angle;
    }
};

template<class Array>
list tangentList(const Array &p, int dx = 5, unsigned int skip = 0)
{
    if(p.size() < 2*dx + 2*skip)
    {
        PyErr_SetString(PyExc_ValueError,
            "tangentList: polygon too small (less than 2*dx + 2*skip points)");
        throw_error_already_set();
    }

    if(dx < 1)
    {
        PyErr_SetString(PyExc_ValueError,
            "tangentList: dx too small (must be >= 1)");
        throw_error_already_set();
    }

    double pos = 0.0;
    for(unsigned int i = 0; i < dx + skip - 1; ++i)
        pos += (p[i+1]-p[i]).magnitude();

    list result;
    ContinuousDirection makeContinuous;

    for(unsigned int i = skip; i < p.size() - 2*dx - skip; ++i)
    {
        typename Array::value_type
            s(p[i+2*dx] - p[i]);

        pos += (p[i+dx] - p[i+dx-1]).magnitude();

        static const double eps = 1e-14;
        if(s.squaredMagnitude() < eps)
            continue;

        result.append(
            make_tuple(pos, makeContinuous(VIGRA_CSTD::atan2(s[1], s[0]))));
    }

    return result;
}

/**
 * Performs a convolution of the values in arcLengthList with a
 * Gaussian centered at arcLengthPos.  arcLengthList should contain
 * (arcLength, value) pairs, and the arc length is used to sample the
 * Gaussian at positions (arcLengthPos-arcLengthPos). The index i is
 * used as a hint where arcLengthPos is found in the arcLengthList
 * (that is - arcLengthList[i][0] must be < arcLengthPos).
 */
template<class ValueType>
tuple singleGaussianConvolveAtArcLength(
    const list &arcLengthList,
    int i, double arcLengthPos,
    const Gaussian<> &g)
{
    //std::cerr << "singleGaussianConvolveAtArcLength(.., " << i << ", " << arcLengthPos << ");\n";

    ValueType sum(vigra::NumericTraits<ValueType>::zero());

    double lnorm = 0.0;
    for(int j = i; j < len(arcLengthList); ++j)
    {
        tuple posCurv((extract<tuple>(arcLengthList[j])()));
        double dist = extract<double>(posCurv[0])() - arcLengthPos;
        if(dist > g.radius())
            break;
        ValueType curv = extract<ValueType>(posCurv[1])();
        double w(g(dist));
        sum += w*curv;
        lnorm += w;
    }

    double rnorm = 0.0;
    for(int j = i - 1; j >= 0; --j)
    {
        tuple posCurv((extract<tuple>(arcLengthList[j])()));
        double dist = extract<double>(posCurv[0])() - arcLengthPos;
        if(dist < -g.radius())
            break;
        ValueType curv = extract<ValueType>(posCurv[1])();
        double w(g(dist));
        sum += w*curv;
        rnorm += w;
    }

    if(!g.derivativeOrder())
        sum /= (lnorm+rnorm);
    else if(g.derivativeOrder() & 1)
    {
        double disparity = std::fabs(rnorm) - std::fabs(lnorm);
        if(disparity > 0)
            sum += disparity * extract<ValueType>(arcLengthList[0][1])();
        else
            sum += -disparity * extract<ValueType>(arcLengthList[-1][1])();
    }
    else
    {
        sum += -lnorm * extract<ValueType>(arcLengthList[0][1])();
        sum += -rnorm * extract<ValueType>(arcLengthList[-1][1])();
    }

    return make_tuple(arcLengthPos, sum);
}

template<class ValueType>
list gaussianConvolveByArcLengthInternal(
    const list &arcLengthList,
    double sigma, int derivativeOrder = 0)
{
    list result;
    Gaussian<> g(sigma, derivativeOrder);
    for(int i = 0; i < len(arcLengthList); ++i)
    {
        double pos = extract<double>(arcLengthList[i][0])();
        result.append(singleGaussianConvolveAtArcLength<ValueType>(
                          arcLengthList, i, pos, g));
    }
    return result;
}

list gaussianConvolveByArcLength(const list &arcLengthList,
                                 double sigma, int derivativeOrder = 0)
{
    if(!len(arcLengthList))
        return list();
    if(extract<double>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthInternal<double>(
            arcLengthList, sigma, derivativeOrder);
    if(extract<Vector2>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthInternal<Vector2>(
            arcLengthList, sigma, derivativeOrder);
    PyErr_SetString(PyExc_TypeError,
        "gaussianConvolveByArcLength(arcLengthList, ...) can only be applied on lists\n"
        "of pairs containing floats or Vector2s as second values to be processed.");
    throw_error_already_set();
    return list(); // never reached
}

template<class ArcLengthList, class PointList>
typename PointList::value_type
singleGaussianConvolveAtArcLengthReflective(
    const ArcLengthList &arcLengthList,
    const PointList &pointList,
    int i, double arcLengthPos,
    const Gaussian<> &g)
{
    typedef typename PointList::value_type ValueType;

    vigra_precondition(arcLengthList.size() == pointList.size(),
        "singleGaussianConvolveAtArcLengthReflective(): Input lists must have the same size.");

    ValueType sum(vigra::NumericTraits<ValueType>::zero());
    double norm = 0.0;

    int lastIndex = arcLengthList.size() - 1;
    double totalLength = arcLengthList[lastIndex];
    vigra_precondition(totalLength > g.radius(),
        "singleGaussianConvolveAtArcLengthReflective(): Filter longer than polygon.");

    ValueType lastPoint = pointList[lastIndex];
    for(int j = i; true; ++j)
    {
        int k = (j > lastIndex)
                ? 2*lastIndex - j
                : j;
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j > lastIndex)
        {
            pos = 2.0*totalLength - pos;
            point = 2.0*lastPoint - point;
        }
        double diff = pos - arcLengthPos;
        if(diff > g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    ValueType firstPoint = pointList[0];
    for(int j = i - 1; true; --j)
    {
        int k = abs(j);
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j < 0)
        {
            pos = -pos;
            point = 2.0*firstPoint - point;
        }
        double diff = pos - arcLengthPos;
        if(diff < -g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    if(!g.derivativeOrder())
        sum /= norm;
    return sum;
}

template<class ArcLengthList, class PointList>
typename PointList::value_type
singleGaussianConvolveAtArcLengthCyclic(
    const ArcLengthList &arcLengthList,
    const PointList &pointList,
    int i, double arcLengthPos,
    const Gaussian<> &g)
{
    typedef typename PointList::value_type ValueType;

    vigra_precondition(arcLengthList.size() == pointList.size(),
        "singleGaussianConvolveAtArcLengthCyclic(): Input lists must have the same size.");

    ValueType sum(vigra::NumericTraits<ValueType>::zero());
    double norm = 0.0;

    int lastIndex = arcLengthList.size() - 1;
    double totalLength = VIGRA_CSTD::sqrt(squaredNorm(pointList[0] - pointList[lastIndex])) + arcLengthList[lastIndex];
    vigra_precondition(totalLength > g.radius(),
        "singleGaussianConvolveAtArcLengthCyclic(): Filter longer than polygon.");

    for(int j = i; true; ++j)
    {
        int k = (j > lastIndex)
                    ? j - (lastIndex + 1)
                    : j;
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j > lastIndex)
        {
            pos += totalLength;
        }
        double diff = pos - arcLengthPos;
        if(diff > g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    for(int j = i - 1; true; --j)
    {
        int k = (j < 0)
                    ? lastIndex + 1 + j
                    : j;
        double pos = arcLengthList[k];
        ValueType point = pointList[k];
        if(j < 0)
        {
            pos -= totalLength;
        }
        double diff = pos - arcLengthPos;
        if(diff < -g.radius())
            break;
        double w(g(diff));
        sum += w*point;
        norm += w;
    }

    if(!g.derivativeOrder())
        sum /= norm;
    return sum;
}

template<class ValueType>
list gaussianConvolveByArcLengthCyclicInternal(
    const list &arcLengthList,
    double sigma, int derivativeOrder = 0)
{
    // convert the input lists
    unsigned int size = len(arcLengthList);
    ArrayVector<double> arcLengths(size);
    ArrayVector<ValueType> points(size);
    for(unsigned int i = 0; i < size; ++i)
    {
        arcLengths[i] = extract<double>(arcLengthList[i][0])();
        points[i] = extract<ValueType>(arcLengthList[i][1])();
    }

    list result;
    Gaussian<> g(sigma, derivativeOrder);
    for(unsigned int i = 0; i < size; ++i)
    {
        double pos = arcLengths[i];
        result.append(make_tuple(pos, singleGaussianConvolveAtArcLengthCyclic(
                          arcLengths, points, i, pos, g)));
    }
    return result;
}

list gaussianConvolveByArcLengthCyclic(const list &arcLengthList,
                                       double sigma, int derivativeOrder = 0)
{
    if(!len(arcLengthList))
        return list();
    if(extract<double>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthCyclicInternal<double>(
            arcLengthList, sigma, derivativeOrder);
    if(extract<Vector2>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthCyclicInternal<Vector2>(
            arcLengthList, sigma, derivativeOrder);
    PyErr_SetString(PyExc_TypeError,
        "gaussianConvolveByArcLengthCyclic(arcLengthList, ...) can only be applied on lists\n"
        "of pairs containing floats or Vector2s as second values to be processed.");
    throw_error_already_set();
    return list(); // never reached
}

template <class Point>
struct AugmentedPolygonPoint
{
    Point point;
    double arcLength, confidence;

    AugmentedPolygonPoint()
    {}

    AugmentedPolygonPoint(Point const & p, double a, double c)
    : point(p), arcLength(a), confidence(c)
    {}
};

template <class Array1, class Array2, class Array3>
void copyPointsReflectiveBoundaryConditions(Array1 const & points, Array2 const & conf, Array3 & res)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();
    res[size] = APoint(points[0], 0.0, conf[0]);
    double arcLength = 0.0;
    for(unsigned int k=1; k<size; ++k)
    {
        arcLength += norm(points[k] - points[k-1]);
        res[size+k] = APoint(points[k], arcLength, conf[k]);
    }
    Point lastPoint2 = 2.0*points[size-1];
    for(unsigned int k=0; k<size-1; ++k)
    {
        Point reflected = lastPoint2 - points[size-2-k];
        arcLength += norm(res[2*size+k-1].point - reflected);
        res[2*size+k] = APoint(reflected, arcLength, conf[size-2-k]);
    }
    Point firstPoint2 = 2.0*points[0];
    arcLength = 0.0;
    for(unsigned int k=0; k<size-1; ++k)
    {
        Point reflected = firstPoint2 - points[k+1];
        arcLength -= norm(res[size-k].point - reflected);
        res[size-1-k] = APoint(reflected, arcLength, conf[k+1]);
    }
}

template <class Array1, class Array2, class Array3>
void copyPointsCyclicBoundaryConditions(Array1 const & points, Array2 const & conf, Array3 & res)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();
    res[size] = APoint(points[0], 0.0, conf[0]);
    double arcLength = 0.0;
    for(unsigned int k=1; k<size; ++k)
    {
        arcLength += norm(points[k] - points[k-1]);
        res[size+k] = APoint(points[k], arcLength, conf[k]);
    }
    for(unsigned int k=0; k<size-1; ++k)
    {
        arcLength += norm(res[2*size+k-1].point - points[k+1]);
        res[2*size+k] = APoint(points[k+1], arcLength, conf[k+1]);
    }
    arcLength = 0.0;
    for(unsigned int k=0; k<size-1; ++k)
    {
        arcLength -= norm(res[size-k].point - points[size-2-k]);
        res[size-1-k] = APoint(points[size-2-k], arcLength, conf[size-2-k]);
    }
}

template<class Point, class Iterator>
Point estimateFirstDerivative(Iterator points, Iterator end, int index, double dist);

template<class Array1, class Array2>
list tangentListChord(const Array1 &points, const Array2 &confidences, double scale)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        Point r = estimateFirstDerivative<Point>(p, p+size, i, scale);

        result.append(make_tuple(p[i].arcLength, makeContinuous(VIGRA_CSTD::atan2(r[1], r[0]))));
    }

    return result;
}

template<class Array1, class Array2>
list tangentListNormalizedGaussian(const Array1 &points, const Array2 &confidences, double scale)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;
    double s2 = scale*scale;

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        Point  sp = p[i].confidence*p[i].point, spp = Point(0.0, 0.0);
        double sw = p[i].confidence, swp = 0.0;
        for(int k = 1; k < (int)size; ++k)
        {
            double diff = p[i+k].arcLength - p[i].arcLength;
            if(diff > 3.0*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp += diff/s2*w*p[i+k].point;
            sw  += w;
            swp += diff/s2*w;
        }
        for(int k = -1; k > -(int)size; --k)
        {
            double diff = p[i].arcLength - p[i+k].arcLength;
            if(diff > 3.0*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp -= diff/s2*w*p[i+k].point;
            sw  += w;
            swp -= diff/s2*w;
        }

        Point r = (sw*spp - swp*sp) / (sw*sw);

        result.append(make_tuple(p[i].arcLength, makeContinuous(VIGRA_CSTD::atan2(r[1], r[0]))));
    }

    return result;
}

template<class Point, class Iterator>
Point estimateFirstDerivative(Iterator points, Iterator end, int index, double dist)
{
    double a0 = points[index].arcLength, d = 0, o;
    int k, size = end  - points;
    for(k = index+1; k<2*size; ++k)
    {
        d = points[k].arcLength - a0;
        if(d > dist)
            break;
    }
    vigra_invariant(k < 2*size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d - dist) / (points[k].arcLength - points[k-1].arcLength);
    Point p1(o * points[k-1].point + (1.0 - o) * points[k].point);

    for(k = index-1; k>-size; --k)
    {
        d = points[k].arcLength - a0;
        if(d < -dist)
            break;
    }
    vigra_invariant(k>-size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d + dist) / (points[k].arcLength - points[k+1].arcLength);
    Point pm1(o * points[k+1].point + (1.0 - o) * points[k].point);

    return (p1 - pm1) / (2.0 * dist);
}

template<class Point, class Iterator>
Point estimateThirdDerivative(Iterator points, Iterator end, int index, double dist)
{
    double a0 = points[index].arcLength, d = 0, o;
    int k, size = end  - points;
    for(k = index+1; k<2*size; ++k)
    {
        d = points[k].arcLength - a0;
        if(d > dist)
            break;
    }
    vigra_invariant(k < 2*size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d - dist) / (points[k].arcLength - points[k-1].arcLength);
    Point p1(o * points[k-1].point + (1.0 - o) * points[k].point);

    for(; k<2*size; ++k)
    {
        d = points[k].arcLength - a0;
        if(d > 2.0*dist)
            break;
    }
    vigra_invariant(k < 2*size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d - 2.0*dist) / (points[k].arcLength - points[k-1].arcLength);
    Point p2(o * points[k-1].point + (1.0 - o) * points[k].point);

    for(k = index-1; k>-size; --k)
    {
        d = points[k].arcLength - a0;
        if(d < -dist)
            break;
    }
    vigra_invariant(k>-size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d + dist) / (points[k].arcLength - points[k+1].arcLength);
    Point pm1(o * points[k+1].point + (1.0 - o) * points[k].point);

    for(; k>-size; --k)
    {
        d = points[k].arcLength - a0;
        if(d < -dist)
            break;
    }
    vigra_invariant(k>-size,
        "estimateThirdDerivative(): Polygon shorter than requested distance.");
    o = (d + 2.0*dist) / (points[k].arcLength - points[k+1].arcLength);
    Point pm2(o * points[k+1].point + (1.0 - o) * points[k].point);
    return (p2 - 2.0*p1 + 2.0*pm1 - pm2) / (2.0 * dist * dist * dist);
}

template<class Point, class Array>
Point estimateThirdDerivativeQuick(const Array & points, int index, double dist, double averageDist)
{
    int diff = (int)std::floor(dist/averageDist);
    int size = points.size();
    vigra_invariant(index + 2*diff < size && index - 2*diff >= 0,
        "estimateThirdDerivativeQuick(): Polygon shorter than requested distance.");
    return (points[index-2*diff].point - 2.0*points[index-diff].point +
                2.0*points[index+diff].point - points[index+2*diff].point) / (2.0 * dist * dist * dist);
}

template<class Array1, class Array2>
list gaussianOptimalScales(const Array1 &points, const Array2 &confidences,
                double sigmaFilter, double thirdDerivScale, double minScale)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;
    double totalArcLength = p[size-1].arcLength;
    double dist = std::min(thirdDerivScale, totalArcLength / 16.0);

    list result;
    for(unsigned int i = 0; i < size; ++i)
    {
        double q3 = norm(estimateThirdDerivative<Point>(p, p + size, i, dist));
        double s = std::pow(9.0/4.0*sq(sigmaFilter / sq(p[i].confidence*q3)), 1.0/7.0);
        if(s < sq(minScale) + sq(sigmaFilter))
            result.append(minScale);
        else
            result.append(std::sqrt(s - sq(sigmaFilter)));
    }

    return result;
}

template<class Array1, class Array2>
list tangentListNormalizedGaussianOptimal(const Array1 &points, const Array2 &confidences, double maxScale)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;
    double totalArcLength = p[size-1].arcLength;
    double dist = std::min(maxScale, totalArcLength / 16.0);

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        double q3 = norm(estimateThirdDerivative<Point>(p, p + size, i, dist));
        double scale = std::pow(0.75/sqrt(M_PI)/sq(p[i].confidence*q3), 1.0/7.0);
        double s2 = scale*scale;
        Point  sp = p[i].confidence*p[i].point, spp = Point(0.0, 0.0);
        double sw = p[i].confidence, swp = 0.0;
        for(int k = 1; k < (int)size; ++k)
        {
            double diff = p[i+k].arcLength - p[i].arcLength;
            if(diff > 2.5*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp += diff/s2*w*p[i+k].point;
            sw  += w;
            swp += diff/s2*w;
        }
        for(int k = -1; k > -(int)size; --k)
        {
            double diff = p[i].arcLength - p[i+k].arcLength;
            if(diff > 2.5*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp -= diff/s2*w*p[i+k].point;
            sw  += w;
            swp -= diff/s2*w;
        }

        Point r = (sw*spp - swp*sp) / (sw*sw);

        result.append(make_tuple(p[i].arcLength, makeContinuous(VIGRA_CSTD::atan2(r[1], r[0]))));
    }

    return result;
}

template<class Array1, class Array2, class Array3>
list tangentListNormalizedGaussianVariableScale(
        const Array1 &points, const Array2 &confidences, const Array3 & scales)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        double scale = scales[i];
        double s2 = scale*scale;
        Point  sp = p[i].confidence*p[i].point, spp = Point(0.0, 0.0);
        double sw = p[i].confidence, swp = 0.0;
        for(int k = 1; k < (int)size; ++k)
        {
            double diff = p[i+k].arcLength - p[i].arcLength;
            if(diff > 2.5*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp += diff/s2*w*p[i+k].point;
            sw  += w;
            swp += diff/s2*w;
        }
        for(int k = -1; k > -(int)size; --k)
        {
            double diff = p[i].arcLength - p[i+k].arcLength;
            if(diff > 2.5*scale)
                break;
            double w = p[i+k].confidence*VIGRA_CSTD::exp(-diff*diff/2.0/s2);
            sp  += w*p[i+k].point;
            spp -= diff/s2*w*p[i+k].point;
            sw  += w;
            swp -= diff/s2*w;
        }

        Point r = (sw*spp - swp*sp) / (sw*sw);

        result.append(make_tuple(p[i].arcLength, makeContinuous(VIGRA_CSTD::atan2(r[1], r[0]))));
    }

    return result;
}

template<class Array1, class Array2>
list tangentListChordOptimal(const Array1 &points, const Array2 &confidences, double maxScale)
{
    typedef typename Array1::value_type Point;
    typedef AugmentedPolygonPoint<Point> APoint;

    unsigned int size = points.size();

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, confidences, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, confidences, aPoints);

    typename ArrayVector<APoint>::iterator p = aPoints.begin() + size;
    double totalArcLength = p[size-1].arcLength;
    double dist = std::min(maxScale, totalArcLength / 16.0);

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        double q3 = norm(estimateThirdDerivative<Point>(p, p + size, i, dist));
        double scale = std::pow(3.0/p[i].confidence/q3, 1.0/3.0);
        Point  r = estimateFirstDerivative<Point>(p, p + size, i, scale);
        result.append(make_tuple(p[i].arcLength, makeContinuous(VIGRA_CSTD::atan2(r[1], r[0]))));
    }

    return result;
}

list pytangentListChord(list const & pyPoints, list const & pyConf, double scale)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf),
        "tangentListChord(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
    }
    return tangentListChord(points, conf, scale);
}

list pytangentListNormalizedGaussian(list const & pyPoints, list const & pyConf, double scale)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf),
        "tangentListNormalizedGaussian(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
    }
    return tangentListNormalizedGaussian(points, conf, scale);
}

list pygaussianOptimalScales(list const & pyPoints, list const & pyConf,
        double filterScale, double thirdDerivScale, double minScale)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf),
        "gaussianOptimalScales(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
    }
    return gaussianOptimalScales(points, conf, filterScale, thirdDerivScale, minScale);
}

list pytangentListNormalizedGaussianOptimal(list const & pyPoints, list const & pyConf, double maxScale)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf),
        "tangentListNormalizedGaussianOptimal(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
    }
    return tangentListNormalizedGaussianOptimal(points, conf, maxScale);
}

list pytangentListNormalizedGaussianVariableScale(list const & pyPoints, list const & pyConf, list const & pyScales)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf) && size == len(pyScales),
        "tangentListNormalizedGaussianVariableScale(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    ArrayVector<double>  scales(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
        scales[k] = extract<double>(pyScales[k])();
    }
    return tangentListNormalizedGaussianVariableScale(points, conf, scales);
}

list pyThirdDerivativeOfPolygon(list const & pyPoints, double dist)
{
    typedef AugmentedPolygonPoint<Vector2> APoint;
    int size = len(pyPoints);

    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = 1.0;
    }

    ArrayVector<APoint> aPoints(3*size);

    if(points[0] == points[size-1])
        copyPointsCyclicBoundaryConditions(points, conf, aPoints);
    else
        copyPointsReflectiveBoundaryConditions(points, conf, aPoints);

    ArrayVector<APoint>::iterator p = aPoints.begin() + size;

    list result;
    for(int i = 0; i < size; ++i)
    {
        double q3 = norm(estimateThirdDerivative<Vector2>(p, p + size, i, dist));
        result.append(make_tuple(p[i].arcLength, q3));
    }

    return result;
}

list pytangentListChordOptimal(list const & pyPoints, list const & pyConf, double maxScale)
{
    int size = len(pyPoints);
    vigra_precondition(size == len(pyConf),
        "tangentListChordOptimal(): lists must have the same size.");
    ArrayVector<Vector2> points(size);
    ArrayVector<double>  conf(size);
    for(int k=0; k<size; ++k)
    {
        points[k] = extract<Vector2>(pyPoints[k])();
        conf[k] = extract<double>(pyConf[k])();
    }
    return tangentListChordOptimal(points, conf, maxScale);
}

template<class Array>
list tangentListGaussianReflective(const Array &points, double sigma, double diff = 0.0)
{
    typedef typename Array::value_type Point;

    unsigned int size = points.size();

    // compute arc lengths
    ArrayVector<double> arcLengths(size);
    arcLengths[0] = 0.0;
    for(unsigned int i = 1; i < size; ++i)
        arcLengths[i] = norm(points[i] - points[i-1]) + arcLengths[i-1];

    Gaussian<> g(sigma);
    vigra_precondition(arcLengths[size-1] > g.radius(),
        "tangentListGaussianReflective(): Filter longer than polygon.");

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        double pos = arcLengths[i];
        double dist = (diff == 0.0)
                        ? (i == 0)
                            ? (arcLengths[i+1] - pos) / 2.0
                            : (i == size - 1)
                                ? (pos - arcLengths[i-1]) / 2.0
                                : std::min(pos - arcLengths[i-1], arcLengths[i+1] - pos) / 2.0
                        : diff;
        double pos1 = (i == 0)
                          ? pos
                          : pos - dist;
        double pos2 = (i == size - 1)
                          ? pos
                          : pos + dist;
        Point p1(singleGaussianConvolveAtArcLengthReflective(arcLengths, points, i, pos1, g));
        Point p2(singleGaussianConvolveAtArcLengthReflective(arcLengths, points, i, pos2, g));
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];

        result.append(make_tuple(pos, makeContinuous(VIGRA_CSTD::atan2(dy, dx))));
    }

    return result;
}

template<class Array>
list tangentListGaussianCyclic(const Array &points, double sigma, double diff = 0.0)
{
    typedef typename Array::value_type Point;

    unsigned int size = points.size();

    // compute arc lengths
    ArrayVector<double> arcLengths(size);
    arcLengths[0] = 0.0;
    for(unsigned int i = 1; i < size; ++i)
        arcLengths[i] = norm(points[i] - points[i-1]) + arcLengths[i-1];
    double restLength = norm(points[0] - points[size-1]);

    Gaussian<> g(sigma);
    vigra_precondition(arcLengths[size-1] > g.radius(),
        "tangentListGaussianCyclic(): Filter longer than polygon.");

    list result;
    ContinuousDirection makeContinuous;
    for(unsigned int i = 0; i < size; ++i)
    {
        double pos = arcLengths[i];
        double dist = (diff == 0.0)
                        ? (i == 0)
                            ? std::min(restLength, arcLengths[i+1] - pos) / 2.0
                            : (i == size - 1)
                                ? std::min(pos - arcLengths[i-1], restLength) / 2.0
                                : std::min(pos - arcLengths[i-1],
                                           arcLengths[i+1] - pos) / 2.0
                        : diff;
        double pos1 = pos - dist;
        double pos2 = pos + dist;
        Point p1(singleGaussianConvolveAtArcLengthCyclic(arcLengths, points, i, pos1, g));
        Point p2(singleGaussianConvolveAtArcLengthCyclic(arcLengths, points, i, pos2, g));
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];

        result.append(make_tuple(pos, makeContinuous(VIGRA_CSTD::atan2(dy, dx))));
    }

    return result;
}

template<class ValueType>
list equidistantGaussiansInternal(
    const list &arcLengthList, double sigma, double distance)
{
    list result;
    Gaussian<> g(sigma);
    double totalLength = extract<double>(arcLengthList[-1][0]);
    int edgeIndex = 0;
    for(double pos = 0; pos < totalLength; pos += distance)
    {
        while(extract<double>(arcLengthList[edgeIndex][0]) < pos)
            ++edgeIndex;

        result.append(singleGaussianConvolveAtArcLength<ValueType>(
                          arcLengthList, edgeIndex, pos, g));
    }
    return result;
}

list equidistantGaussians(const list &arcLengthList,
                          double sigma, double distance)
{
    if(extract<double>(arcLengthList[0][1]).check())
        return equidistantGaussiansInternal<double>(
            arcLengthList, sigma, distance);
    if(extract<Vector2>(arcLengthList[0][1]).check())
        return equidistantGaussiansInternal<Vector2>(
            arcLengthList, sigma, distance);
    PyErr_SetString(PyExc_TypeError,
        "equidistantGaussians(arcLengthList, ...) can only be applied on lists\n"
        "of pairs containing floats or Vector2s as second values to be processed.");
    throw_error_already_set();
    return list(); // never reached
}

double appendTangentList(
    list &tangents1, double arcLengthOffset,
    const list &tangents2, double offset,
    ContinuousDirection &makeContinuous)
{
    int size = len(tangents2);
    if(offset > 0)
    {
        for(int j = 0; j < size; ++j)
        {
            double arcLength((extract<double>(tangents2[j][0])()));
            double angle((extract<double>(tangents2[j][1])()));
            tangents1.append(make_tuple(arcLength + arcLengthOffset,
                                        makeContinuous(angle)));
        }
        arcLengthOffset += offset;
    }
    else
    {
        arcLengthOffset += -offset;
        for(int j = size - 1; j >= 0; --j)
        {
            double arcLength((extract<double>(tangents2[j][0])()));
            double angle((extract<double>(tangents2[j][1])())+M_PI);
            tangents1.append(make_tuple(arcLengthOffset - arcLength,
                                        makeContinuous(angle)));
        }
    }
    return arcLengthOffset;
}

double appendTangentList(
    list &tangents1, double arcLengthOffset,
    const list &tangents2, double offset)
{
    ContinuousDirection makeContinuous;
    if(len(tangents1))
        makeContinuous.prevAngle = extract<double>(tangents1[-1]);

    return appendTangentList(tangents1, arcLengthOffset,
                             tangents2, offset, makeContinuous);
}

list composeTangentLists(const list &tangentLists)
{
    list result;
    double arcLengthOffset = 0.0;
    ContinuousDirection makeContinuous;
    for(int i = 0; i < len(tangentLists); ++i)
    {
        tuple tllPair((extract<tuple>(tangentLists[i])()));
        list tangentList((extract<list>(tllPair[0])()));
        double length((extract<double>(tllPair[1])()));

        arcLengthOffset = appendTangentList(
            result, arcLengthOffset,
            tangentList, length, makeContinuous);
    }
    return result;
}

struct LineFit : public IncrementalFitLine
{
    void addPoints(const PointArray<Vector2> &a)
    {
        for(unsigned int i = 0; i < a.size(); ++i)
            addPoint(a[i]);
    }

    tuple pyComputeImplictEquation() const
    {
        double a,b,c;
        double result = computeImplictEquation(a, b, c);
        return make_tuple(result, make_tuple(a, b, c));
    }

    tuple pyComputeParametricEquation() const
    {
        Vector2 c, o;
        double result = computeParametricEquation(c, o);
        return make_tuple(result, make_tuple(c, o));
    }
};

struct ParabolaFit
{
    ContinuousDirection makeContinuous;
    double sumAl4, sumAl3, sumAl2, sumAl;
    double sumAl2Theta, sumAlTheta, sumTheta, sumTheta2;
    int count;

    mutable bool fitDirty;
    mutable double p0, p1, p2;

    ParabolaFit()
    {
        sumAl4 = 0.0, sumAl3 = 0.0, sumAl2 = 0.0, sumAl = 0.0;
        sumAl2Theta = 0.0, sumAlTheta = 0.0, sumTheta = 0.0, sumTheta2 = 0.0;
        count = 0;

        fitDirty = true;
    }

    void addTangentList(const list &arcLengthList)
    {
        int size = len(arcLengthList);

        for(int i = 0; i < size; ++i)
        {
            double al = extract<double>(arcLengthList[i][0])();
            double theta = extract<double>(arcLengthList[i][1])();
            theta = makeContinuous(theta);

            sumAl     += al;
            double al2 = al*al;
            sumAl2    += al2;
            double al3 = al2*al;
            sumAl3    += al3;
            double al4 = al3*al;
            sumAl4    += al4;

            sumAl2Theta += theta*al*al;
            sumAlTheta  += theta*al;
            sumTheta    += theta;
            sumTheta2   += theta*theta;
        }

        count += size;
        fitDirty = true;
    }

    double sumOfSquaredErrors0() const
    {
        if(!count)
            return 0.0;

        ensureFit();

        // sum[(theta^2 - (p0*al^2 + p1*al + p2))^2]:
        return sumTheta2 +
            p0*p0*sumAl4 +   p1*p1*sumAl2     + p2*p2*count      +
          2*p0*p1*sumAl3 + 2*p0*p2*sumAl2     - 2*p0*sumAl2Theta +
          2*p1*p2*sumAl  -    2*p1*sumAlTheta - 2*p2*sumTheta;
    }

    double sumOfSquaredErrors(const list &xyList) const
    {
        if(!count)
            return 0.0;

        ensureFit();

        int size = len(xyList);
        double error = 0.0;
        for(int i = 0; i < size; ++i)
        {
            double al    = extract<double>(xyList[i][0])();
            double theta = extract<double>(xyList[i][1])();
            error += squaredNorm((p0*al + p1)*al + p2 - theta);
        }
        return error;
    }

    double meanSquaredError(const list &xyList) const
    {
        return sumOfSquaredErrors(xyList) / count;
    }

        // temporary function, to be removed again:
    double oldCovarianceOfResiduals() const
    {
        if(!count)
            return 0.0;

        vigra::Matrix<double> matrix(3, 3);
        matrix(0, 0) = sumAl4; matrix(0, 1) = sumAl3; matrix(0, 2) = sumAl2;
        matrix(1, 0) = sumAl3; matrix(1, 1) = sumAl2; matrix(1, 2) = sumAl;
        matrix(2, 0) = sumAl2; matrix(2, 1) = sumAl;  matrix(2, 2) = count;

        vigra::Matrix<double> covMatrix(3, 3);
        try
        {
            inverse(matrix, covMatrix);
        }
        catch(vigra::PreconditionViolation)
        {
            return 0.0; // not invertible?!
        }

        double result =
            sumAl4*covMatrix(0,0) +
            sumAl2*covMatrix(1,1) +
            count*covMatrix(2,2);

        return result;
    }

    double totalCovarianceOfResiduals() const
    {
        if(!count)
            return 0.0;

        vigra::Matrix<double> matrix(3, 3);
        matrix(0, 0) = sumAl4; matrix(0, 1) = sumAl3; matrix(0, 2) = sumAl2;
        matrix(1, 0) = sumAl3; matrix(1, 1) = sumAl2; matrix(1, 2) = sumAl;
        matrix(2, 0) = sumAl2; matrix(2, 1) = sumAl;  matrix(2, 2) = count;

        vigra::Matrix<double> covMatrix(3, 3);
        try
        {
            inverse(matrix, covMatrix);
        }
        catch(vigra::PreconditionViolation)
        {
            return 0.0; // not invertible?!
        }

        double result =
            sumAl4*covMatrix(0,0) +
            sumAl3*2*covMatrix(0,1) +
            sumAl2*(2*covMatrix(0,2) + covMatrix(1,1)) +
            sumAl*2*covMatrix(1,2) +
            count*covMatrix(2,2);

        return result;
    }

    double meanCovarianceOfResiduals() const
    {
        if(!count)
            return 0.0;

        return totalCovarianceOfResiduals() / count;
    }

    double testMeanCovarianceOfResiduals(const list &xyList)
    {
        ParabolaFit tmp(*this);
        tmp.addTangentList(xyList);
        return tmp.meanCovarianceOfResiduals();
    }

    tuple parabolaParams() const
    {
        ensureFit();
        return make_tuple(p0, p1, p2);
    }

    double parabola(double x) const
    {
        ensureFit();
        return (p0*x + p1)*x + p2;
    }

  protected:
    void ensureFit() const
    {
        if(fitDirty)
        {
            vigra::Matrix<double> matrix(3, 3);
            matrix(0, 0) = sumAl4; matrix(0, 1) = sumAl3; matrix(0, 2) = sumAl2;
            matrix(1, 0) = sumAl3; matrix(1, 1) = sumAl2; matrix(1, 2) = sumAl;
            matrix(2, 0) = sumAl2; matrix(2, 1) = sumAl;  matrix(2, 2) = count;

            vigra::Matrix<double> ergs(3, 1);
            ergs(0, 0) = sumAl2Theta;
            ergs(1, 0) = sumAlTheta;
            ergs(2, 0) = sumTheta;

            vigra::Matrix<double> result(3, 1);
            linearSolve(matrix, ergs, result);
            p0 = result(0, 0);
            p1 = result(1, 0);
            p2 = result(2, 0);

            fitDirty = false;
        }
    }
};

/********************************************************************/

void markEdgeInLabelImage(
    const Scanlines &scanlines,
    PythonImage &labelVImage)
{
    PythonSingleBandImage labelImage(labelVImage.subImage(0));

    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(labelImage.height(), scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.width())
                end = labelImage.width();

            for(int x = begin; x < end; ++x)
            {
                PythonSingleBandImage::reference old(labelImage(x, y));
                if(old < 0)
                    old -= 1;
                else
                    old = -1;
            }
        }
    }
}

list removeEdgeFromLabelImage(
    const Scanlines &scanlines,
    PythonImage &labelVImage,
    GrayValue substituteLabel)
{
    PythonSingleBandImage labelImage(labelVImage.subImage(0));

    list result;
    // clip to image range vertically:
    int y = std::max(0, scanlines.startIndex()),
     endY = std::min(labelImage.height(), scanlines.endIndex());

    for(; y < endY; ++y)
    {
        const Scanlines::Scanline &scanline(scanlines[y]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            // X range checking
            int begin = scanline[j].begin,
                  end = scanline[j].end;
            if(begin < 0)
                begin = 0;
            if(end > labelImage.width())
                end = labelImage.width();

            for(int x = begin; x < end; ++x)
            {
                PythonSingleBandImage::reference old(labelImage(x, y));
                if(old != -1)
                {
                    old += 1;
                }
                else
                {
                    old = substituteLabel;
                    result.append(Point2D(x, y));
                }
            }
        }
    }
    return result;
}

template<class Array>
struct ArrayPickleSuite : pickle_suite
{
    static tuple getinitargs(Array const& v)
    {
        list l(v);
        return make_tuple(l);
    }
};

template<class Polygon>
list intersectLine(
    const Polygon &polygon, Vector2 const &lineStart, Vector2 const &lineEnd)
{
    list result;

    if(!polygon.size())
        return result;

    Vector2
        lineDir(lineEnd - lineStart),
        lineNormal(lineDir[1], -lineDir[0]);
    lineNormal = lineNormal / lineNormal.magnitude();

    Polygon currentPart;
    bool inside = dot(polygon[0] - lineStart, lineNormal) >= 0;

    double prevDist = 0.0;
    for(unsigned int pointIndex = 0; pointIndex < polygon.size(); ++pointIndex)
    {
        // skip outside points
        double dist = dot(polygon[pointIndex] - lineStart, lineNormal);
        if(dot(polygon[pointIndex] - lineStart, lineNormal) < 0)
        {
            if(inside)
            {
                currentPart.push_back(
                    polygon[pointIndex-1] - (prevDist / (dist - prevDist)) *
                    (polygon[pointIndex] - polygon[pointIndex-1]));
                result.append(currentPart);
                currentPart = Polygon();
                inside = false;
            }
            prevDist = dist;
            continue;
        }

        if(!inside)
        {
            // previous is outside, this is inside:
            currentPart.push_back(
                polygon[pointIndex-1] - (prevDist / (dist - prevDist)) *
                (polygon[pointIndex] - polygon[pointIndex-1]));
            inside = true;
        }

        currentPart.push_back(polygon[pointIndex]);
    }

    if(inside)
        result.append(currentPart);

    return result;
}

double angleTheta(double dy, double dx); // implemented in cppmap.cxx

void defPolygon()
{
    def("angleTheta", &angleTheta,
        "angleTheta(dy, dx)\n"
        "  calculates an efficient substitute value with the same sorting behavior\n"
        "  as the angle returned by math.atan2(dy, dx) (but ranged -2..2).\n"
        "  It does *not* give the angle, and seems to be only around 17% faster.");

    typedef BBoxPolygon<Vector2> PythonPolygon;
    typedef PythonPolygon::Base PythonPolygonBase;

    typedef PythonPolygonBase::Base Vector2Array;
    typedef STLIterWrapper<Vector2Array::const_iterator>
        VectorIter;
    defIter<VectorIter>("Vector2Iter");
    typedef STLIterWrapper<Vector2Array::const_reverse_iterator>
        VectorRevIter;
    defIter<VectorRevIter>("Vector2RevIter");
    class_<Vector2Array>("Vector2Array")
        .def(init<Vector2Array>())
        //.def("__init__", make_constructor(&createPolygon<Vector2Array>))
        .def("reverse", &Vector2Array::reverse)
        .def("interpolate", &Vector2Array::interpolate)
        .def("__len__", &Vector2Array::size)
        .def("__getitem__", &Array__getitem_slice__<Vector2Array>)
        .def("__getitem__", &Array__getitem__<Vector2Array>)
        .def("__setitem__", &Array__setitem__<Vector2Array>)
        .def("__iter__",    &Array__iter__<Vector2Array>,
             with_custodian_and_ward_postcall<0, 1>())
        .def("__reviter__", &Array__reviter__<Vector2Array>,
             with_custodian_and_ward_postcall<0, 1>())
        .def("insert", &insert<Vector2Array>)
        .def(self * double())
        .def(self + Vector2())
        .def("roundToInteger", &Vector2Array::roundToInteger)
    ;
    PolygonFromPython<Vector2Array>();
    register_ptr_to_python< std::auto_ptr<Vector2Array> >();

    typedef PointArray<Point2D> Point2DArray;
    typedef STLIterWrapper<Point2DArray::const_iterator>
        Point2DIter;
    defIter<Point2DIter>("Point2DIter");
    typedef STLIterWrapper<Point2DArray::const_reverse_iterator>
        Point2DRevIter;
    defIter<Point2DRevIter>("Point2DRevIter");
    class_<Point2DArray>("Point2DArray")
        .def(init<Point2DArray>())
        .def("reverse", &Point2DArray::reverse)
        .def("__len__", &Point2DArray::size)
        .def("__getitem__", &Array__getitem_slice__<Point2DArray>)
        .def("__getitem__", &Array__getitem__<Point2DArray>)
        .def("__setitem__", &Array__setitem__<Point2DArray>)
        .def("__iter__",    &Array__iter__<Point2DArray>,
             with_custodian_and_ward_postcall<0, 1>())
        .def("__reviter__", &Array__reviter__<Point2DArray>,
             with_custodian_and_ward_postcall<0, 1>())
        .def("insert", &insert<Point2DArray>)
    ;
    register_ptr_to_python< std::auto_ptr<Point2DArray> >();

    class_<PythonPolygon, bases<Vector2Array> >("Polygon")
        .def(init<Vector2Array>())
        //.def("__init__", make_constructor(&createPolygon<PythonPolygon>))
        .def("__delitem__", &erase<PythonPolygon>)
        .def("insert", &insert<PythonPolygon>)
        .def("append", &PythonPolygon::push_back)
        .def("extend", &PythonPolygon::extend)
        .def("__getitem__", &Array__getitem_slice__<PythonPolygon>)
        .def("__getitem__", &Array__getitem__<PythonPolygon>)
        .def("__setitem__", &Polygon__setitem__<PythonPolygon>)
        .def("split", &split<PythonPolygon>)
        .def("closed", &PythonPolygon::closed,
             "Return True iff poly[-1] == poly[0]")
        .def("length", &PythonPolygon::length)
        .def("partialArea", &PythonPolygon::partialArea)
        .def("boundingBox", &PythonPolygon::boundingBox)
        .def("contains", &PythonPolygon::contains)
        .def("swap", &PythonPolygon::swap)
        .def("reverse", &PythonPolygon::reverse)
        .def("nearestPoint", &PythonPolygon::nearestPoint)
        .def("invalidateProperties", &PythonPolygon::invalidateProperties)
        .def("arcLengthList", &arcLengthList<PythonPolygon>)
        .def("__repr__", &Polygon__repr__<PythonPolygon>)
        .def_pickle(ArrayPickleSuite<PythonPolygon>())
    ;
    PolygonFromPython<PythonPolygon>();
    register_ptr_to_python< std::auto_ptr<PythonPolygon> >();

    typedef PythonPolygon::BoundingBox BoundingBox;
    class_<BoundingBox>("BoundingBox")
        .def(init<BoundingBox>())
        .def(init<const Vector2 &>(args("size")))
        .def(init<const Vector2 &, const Vector2 &>(args("begin", "end")))
        .def("__init__", make_constructor(&createBoxFromRect2D<BoundingBox>))
        .def("begin", (const Vector2 &(BoundingBox::*)() const)
             &BoundingBox::begin, return_internal_reference<>())
        .def("end", (const Vector2 &(BoundingBox::*)() const)
             &BoundingBox::end, return_internal_reference<>())
        .def("setBegin", &BoundingBox::setBegin, arg("begin"))
        .def("setEnd", &BoundingBox::setEnd, arg("end"))
        .def("moveTo", &BoundingBox::moveTo, arg("newBegin"))
        .def("moveBy", &BoundingBox::moveBy, arg("offset"))
        .def("area", &BoundingBox::volume)
        .def("size", &BoundingBox::size)
        .def("setSize", &BoundingBox::setSize)
        .def("addSize", &BoundingBox::addSize)
        .def("addBorder", &BoundingBox::addBorder)
        .def("isEmpty", &BoundingBox::isEmpty)
        .def("contains", (bool (BoundingBox::*)(const Vector2 &) const)
             &BoundingBox::contains)
        .def("contains", (bool (BoundingBox::*)(const BoundingBox &) const)
             &BoundingBox::contains)
        .def("intersects", &BoundingBox::intersects)
        .def("__repr__", &Box__repr__<BoundingBox>) // TODO: __str__?
        .def(self == self)
        .def(self != self)
        .def(self |= Vector2())
        .def(self | Vector2())
        .def(self |= self)
        .def(self | self)
        .def(self &= self)
        .def(self & self)
        .def_pickle(BoxPickleSuite<BoundingBox>())
    ;

    class_<Scanlines>("Scanlines", no_init)
        .def("__len__", &Scanlines::size)
        .def("__getitem__", &Scanlines__getitem__, return_internal_reference<>())
        .def("startIndex", &Scanlines::startIndex)
        .def("endIndex", &Scanlines::endIndex)
        .def("points", &createScanlinesIter,
             with_custodian_and_ward_postcall<0, 1>())
    ;

    RangeIterWrapper<ScanlinesIter>("_ScanlinesIter");

    class_<Scanlines::Scanline>("Scanline", no_init)
        .def("__len__", &Scanlines::Scanline::size)
        .def("__getitem__", &Array__getitem__byref<Scanlines::Scanline>,
             return_internal_reference<>())
    ;

    class_<ScanlineSegment>("ScanlineSegment", no_init)
        .def_readonly("begin", &ScanlineSegment::begin)
        .def_readonly("direction", &ScanlineSegment::direction)
        .def_readonly("end", &ScanlineSegment::end)
        .def(self == self)
    ;

    register_ptr_to_python< std::auto_ptr<Scanlines> >();

    def("scanPoly", (std::auto_ptr<Scanlines>(*)
                     (const PythonPolygon&))&scanPoly,
        args("polygon"));
    def("scanPoly", (std::auto_ptr<Scanlines>(*)
                     (const Vector2Array &, unsigned int, int))&scanPoly,
        (arg("points"), arg("scanLineCount"), arg("startIndex") = 0));
    def("fillScannedPoly", &pyFillScannedPoly);
    def("drawScannedPoly", &pyDrawScannedPoly);
    def("markEdgeInLabelImage", &markEdgeInLabelImage);
    def("removeEdgeFromLabelImage", &removeEdgeFromLabelImage);

    def("convexHull", (PythonPolygon (*)(const PythonPolygon &))&pyConvexHull);

    def("simplifyPolygon",
        (Vector2Array (*)(const Vector2Array &,double))&simplifyPolygon,
        args("points", "perpendicularDistEpsilon"));
    def("simplifyPolygon",
        (PythonPolygon (*)(const PythonPolygon &,double))&simplifyPolygon,
        args("polygon", "perpendicularDistEpsilon"));
    def("simplifyPolygon",
        (Vector2Array (*)(const Vector2Array &,double,double))&simplifyPolygon,
        args("points", "perpendicularDistEpsilon", "maxStep"));
    def("simplifyPolygon",
        (PythonPolygon (*)(const PythonPolygon &,double,double))&simplifyPolygon),
        args("points", "perpendicularDistEpsilon", "maxStep");
    def("simplifyPolygonDigitalLine",
        (Vector2Array (*)(const Vector2Array &,int))&simplifyPolygonDigitalLine),
        args("points", "connectivity");
    def("simplifyPolygonDigitalLine",
        (PythonPolygon (*)(const PythonPolygon &,int))&simplifyPolygonDigitalLine),
        args("points", "connectivity");

    def("resamplePolygon",
        (PythonPolygon (*)(const PythonPolygon &,double))&resamplePolygon,
        args("polygon", "desiredPointDistance"));
    def("resamplePolygonLinearInterpolation",
        (PythonPolygon (*)(const PythonPolygon &,double))&resamplePolygonLinearInterpolation,
        args("polygon", "desiredPointDistance"));
    def("resamplePolygonExponentialFilter",
        (PythonPolygon (*)(const PythonPolygon &,double,double))&resamplePolygonExponentialFilter,
        args("polygon", "scale", "desiredPointDistance"));
    def("resamplePolygonGaussianFilter",
        (PythonPolygon (*)(const PythonPolygon &,double,double))&resamplePolygonGaussianFilter,
        args("polygon", "scale", "desiredPointDistance"));
    def("polygonSplineControlPoints",
        (PythonPolygon (*)(const PythonPolygon &,int))&polygonSplineControlPoints,
        args("polygon", "segmentCount"));

    def("spline3Integral", vigra::detail::spline3Integral<Vector2>);

    def("centroid", &pyCentroid<PythonPolygon>);

    def("curvatureList", &curvatureList<Vector2Array>,
        (arg("pointArray"), arg("dx") = 5, arg("skipPoints") = 1),
        "curvatureList(pointArray, dx = 5, skipPoints = 1)\n"
        "calculates curvatures values for each triangle between point triples\n"
        "with indices (i-dx, i, i+dx), ignoring skipPoints points from both ends.\n"
        "returns a list of (arcLength, curvature) pairs,\n"
        "whose length is len(pointArray) - 2*dx - 2*skipPoints.\n"
        "(It is possible to get a smaller result if some of the triangles have\n"
        "sides smaller than eps=1e-8 and will be skipped.)");

    def("thirdDerivativeOfPolygon", &pyThirdDerivativeOfPolygon,
        (arg("pointArray"), arg("sigma")),
        "thirdDerivativeOfPolygon(pointArray, sigma): "
        "calculates the magnitude of the third derivative at each point using a symmetric difference at distance 'sigma'.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("tangentList", &tangentList<Vector2Array>,
        (arg("pointArray"), arg("dx") = 5, arg("skipPoints") = 1),
        "tangentList(pointArray, dx = 5, skipPoints = 1)\n"
        "calculates tangent angles from each chord between the point with\n"
        "indices (i-dx, i+dx), ignoring skipPoints points from both ends.\n"
        "returns a list of (arcLength, angle) pairs,\n"
        "whose length is len(pointArray) - 2*dx - 2*skipPoints.");
    def("tangentListGaussianReflective", &tangentListGaussianReflective<Vector2Array>,
        (arg("pointArray"), arg("sigma"), arg("diff") = 0.0),
        "tangentListGaussianReflective(pointArray, sigma): calculates tangent angles at each point.\n"
        "It first approximates the curve at +-diff of each point by means of\n"
        "a Gaussian filter with standard deviation sigma and reflective boundary conditions\n"
        "and then approximates the tangent by the chord between these two points.");
    def("tangentListGaussianCyclic", &tangentListGaussianCyclic<Vector2Array>,
        (arg("pointArray"), arg("sigma"), arg("diff") = 0.0),
        "tangentListGaussianCyclic(pointArray, sigma): calculates tangent angles at each point.\n"
        "It first approximates the curve at +-diff of each point by means of\n"
        "a Gaussian filter with standard deviation sigma and cyclic boundary conditions\n"
        "and then approximates the tangent by the chord between these two points.");
    def("tangentListChord", &pytangentListChord,
        (arg("pointArray"), arg("confidenceArray"), arg("sigma")),
        "tangentListChord(pointArray, confidenceArray, sigma): "
        "calculates tangent angles at each point using a symmetric difference at distance 'sigma'.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("tangentListNormalizedGaussian", &pytangentListNormalizedGaussian,
        (arg("pointArray"), arg("confidenceArray"), arg("sigma")),
        "tangentListNormalizedGaussian(pointArray, confidenceArray, sigma): "
        "calculates tangent angles at each point using a normalized Gaussian derivative at scale 'sigma'.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("tangentListChordOptimal", &pytangentListChordOptimal,
        (arg("pointArray"), arg("confidenceArray"), arg("maxScale")),
        "tangentListChordOptimal(pointArray, confidenceArray, maxScale): "
        "calculates tangent angles at each point using a symmetric difference filter.\n"
        "Its optimal scale is determined from the local confidence and third derivative.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("gaussianOptimalScales", &pygaussianOptimalScales,
        (arg("pointArray"), arg("confidenceArray"), arg("filterScale"), arg("thirdDerivScale"), arg("minScale")),
        "gaussianOptimalScales(pointArray, confidenceArray, filterScale, thirdDerivScale, minScale): "
        "calculates the optimal scale for a Gaussian tangent filter at each point\n"
        "using a third derivative at 'thirdDerivScale' and assuming the the polygon comes from.\n"
        "an edge detector at 'filterScale'. The scales returned are boundaed below by 'minScale'.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("tangentListNormalizedGaussianOptimal", &pytangentListNormalizedGaussianOptimal,
        (arg("pointArray"), arg("confidenceArray"), arg("maxScale")),
        "tangentListNormalizedGaussianOptimal(pointArray, confidenceArray, maxScale): "
        "calculates tangent angles at each point using a normalized Gaussian derivative.\n"
        "Its optimal scale is determined from the local confidence and third derivative.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("tangentListNormalizedGaussianVariableScale", &pytangentListNormalizedGaussianVariableScale,
        (arg("pointArray"), arg("confidenceArray"), arg("scales")),
        "tangentListNormalizedGaussianOptimal(pointArray, confidenceArray, scalesArray): "
        "calculates tangent angles at each point using a normalized Gaussian derivative.\n"
        "The local scale is taken from the scalesArray.\n"
        "Cyclic boundary conditions will be used when the first and last point "
        "of 'pointArray' are equal, reflective boundary conditions otherwise.");
    def("gaussianConvolveByArcLength", &gaussianConvolveByArcLength,
        (arg("arcLengthList"), arg("sigma"), arg("derivativeOrder") = 0),
        "gaussianConvolveByArcLength(arcLengthList, sigma, derivativeOrder = 0)\n"
        "performs a convolution with a Gaussian of the specified derivative\n"
        "order (0 by default) for each element in arcLengthList.\n\n"
        "arcLengthList should contain (arcLength, value) tuples, and\n"
        "the Gaussians will be placed and sampled according to the arcLengths.\n"
        "The values to be convolved may be either floating point numbers or\n"
        "Vector2 objects.\n\n"
        "Border treatment: The Gaussians are simply cut at the ends of\n"
        "arcLengthList.  If derivativeOrder == 0, a smoothing is performed\n"
        "and each convolution is normalized.\n"
        "For derivativeOrders > 0, the convolutions cannot be normalized, so\n"
        "the points are assumed to be equidistant.  However, if the derivative\n"
        "order is uneven, a disparity of the sampled, cut Gaussian is detected\n"
        "and the last value is given an appropriate extra weight (\"repeat\"\n"
        "Border treatment).  If the derivative order is even (and > 0), the\n"
        "same is done to normalize each half of the filter to zero.\n");
    def("gaussianConvolveByArcLengthCyclic", &gaussianConvolveByArcLengthCyclic,
        (arg("arcLengthList"), arg("sigma"), arg("derivativeOrder") = 0));
    def("equidistantGaussians", &equidistantGaussians,
        args("edgeArcLengthList", "sigma", "distance"));
    // FIXME: document properly
//     def("appendTangentList", (double(*)(list &, double, const list &, double, ContinuousDirection &))&appendTangentList,
//         args("tangents1", "arcLengthOffset", "tangents2", "offset", "cont"));
    def("appendTangentList", (double(*)(list &, double, const list &, double))&appendTangentList,
        args("tangents1", "arcLengthOffset", "tangents2", "offset"),
        "appendTangentList(tangents1, arcLengthOffset, tangents2, offset)\n\n"

        "Append tangent list tangents2 to existing tangents1.\n"
        "arcLengthOffset is the offset by which the tangents in tangents2 are\n"
        "shifted.  offset is the arc length of the second curve; if it is\n"
        "negative, tangents2 is reversed on-the-fly.  The result of\n"
        "the append operation is that tangents1 has been extended\n"
        "(i.e. modified in-place) and the resulting arcLengthOffset is\n"
        "returned (i.e. arcLengthOffset + abs(offset)).");

    def("composeTangentLists", &composeTangentLists,
        args("tangentLists"),
        "Composes tangent lists of a chain of (directed) edges into a common one.\n"
        "The parameter should be a list of (edge, length) pairs, where lenght is\n"
        "used to increase the arc length offset for the following edges.\n"
        "The length may also be negative, which means that the tangent list is\n"
        "to be reversed (the curve the tangents come from is traversed in the\n"
        "opposite orientation).");

    class_<LineFit>("LineFit")
        .def(init<LineFit>())
        .def("addPoint",
             (void (LineFit::*)(const Vector2 &))&LineFit::addPoint<Vector2>)
        .def("addPoints",
             &LineFit::addPoints)
        .def("computeImplictEquation",
             &LineFit::pyComputeImplictEquation,
             "computeImplictEquation() -> res, (a, b, p)\n\n"
             "Compute parameters such that ::\n\n"
             "  a * x + b * y + p == 0\n\n"
             "holds for all (x,y) on the line.\n"
             "'res' is the residual of the estimate, namely\n"
             "the std.dev. perpendicular to the line.")
        .def("computeParametricEquation",
             &LineFit::pyComputeParametricEquation,
             "computeParametricEquation() -> res, (center, orientation)\n\n"
             "Compute center, orientation such that ::\n\n"
             "  (x,y) = center + t * orientation\n\n"
             "holds for all (x,y) on the line.\n"
             "'res' is the residual of the estimate, namely\n"
             "the std.dev. perpendicular to the line.")
    ;

    class_<ParabolaFit>("ParabolaFit")
        .def(init<ParabolaFit>())
        .def("addTangentList", &ParabolaFit::addTangentList, arg("arcLengthList"),
             "addTangentList(arcLengthList)\n"
             "Adds list of (arcLength/theta) pairs to the fitted data.\n"

             "The theta values (orientation angles) are made a continuation of\n"
             "the last values (by shifting by multiples of PI).")
        .def("sumOfSquaredErrors", &ParabolaFit::sumOfSquaredErrors0)
        .def("sumOfSquaredErrors", &ParabolaFit::sumOfSquaredErrors,
             "FIXME: The angles within this list are NOT shifted/made continuous!")
        .def("meanSquaredError", &ParabolaFit::meanSquaredError,
             "meanSquaredError(xyList)\n"
             "Returns the same as sumOfSquaredErrors() / count.")
        .def("oldCovarianceOfResiduals", &ParabolaFit::oldCovarianceOfResiduals)
        .def("totalCovarianceOfResiduals", &ParabolaFit::totalCovarianceOfResiduals)
        .def("meanCovarianceOfResiduals", &ParabolaFit::meanCovarianceOfResiduals)
        .def("testMeanCovarianceOfResiduals", &ParabolaFit::testMeanCovarianceOfResiduals)
        .def("parabolaParams", &ParabolaFit::parabolaParams,
             "parabolaParams() returns the parameter tuple (a, b, c) of the\n"
             "fitted parabola function a*x^2 + b*x + c.")
        .def("parabola", &ParabolaFit::parabola, args("x"),
             "parabola(x) returns the value of the parabola function at x")
        .def_readonly("count", &ParabolaFit::count,
                      "the number of values included in the fit")
    ;

    def("intPos", &intPos);
    def("intPos", &intPos_Box<BoundingBox>);

    def("intersectLine", &intersectLine<PythonPolygon>);
}
