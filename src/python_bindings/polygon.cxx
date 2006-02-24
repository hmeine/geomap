#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <vigra/gaussians.hxx>
#include <vigra/pythonimage.hxx>
#include <vigra/linear_algebra.hxx>
#include <cmath>
#include "vigra/polygon.hxx"
#include "delaunay.hxx"

using namespace vigra;
using namespace boost::python;

Point2D intPos(const Vector2 &p)
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
PyObject * Box__repr__(Box const & b)
{
    std::stringstream s;

    s << "<BoundingBox " << b.begin() << ".." << b.end() << ">";

    return PyString_FromString(s.str().c_str());
}

double angleTheta(double dy, double dx)
{
    double denom = std::fabs(dx) + std::fabs(dy);
    if(!denom)
        return 0.0;
    double result = dy / denom;
    if(dx < 0)
    {
        result = 2 - result;
        if(dy < 0)
            result = result - 4;
    }
    return result;
}

inline void checkPointIndex(int &i, int size)
{
    if(i < 0)
        i += size;
    if(i < 0 || i >= size)
    {
        PyErr_SetString(PyExc_IndexError,
            "point index out of bounds.");
        throw_error_already_set();
    }
}

template<class Array>
typename Array::value_type
Array__getitem__(Array const & a, int i)
{
    checkPointIndex(i, a.size());
    return a[i];
}

template<class Array>
const typename Array::value_type &
Array__getitem__byref(Array const & a, int i)
{
    checkPointIndex(i, a.size());
    return a[i];
}

const Scanlines::value_type &
Scanlines__getitem__(Scanlines const & s, unsigned int i)
{
    if(i < 0 || i >= s.size())
//     if(i < s.startIndex() || i >= s.endIndex())
    {
        PyErr_SetString(PyExc_IndexError,
            "scanline index out of bounds.");
        throw_error_already_set();
    }
    return s[i+s.startIndex()];
//     return s[i];
}

template<class Array>
void
Array__setitem__(Array & a, int i, typename Array::value_type v)
{
    checkPointIndex(i, a.size());
    a[i] = v;
}

template<class Polygon>
void
Polygon__setitem__(Polygon & p, int i, typename Polygon::value_type v)
{
    checkPointIndex(i, p.size());
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
    checkPointIndex(pos, p.size() + 1);
    p.insert(p.begin() + pos, x);
}

template<class Polygon>
void erase(Polygon & p, int pos)
{
    checkPointIndex(pos, p.size() + 1);
    p.erase(p.begin() + pos);
}

template<class Polygon>
Polygon split(Polygon & p, int pos)
{
    checkPointIndex(pos, p.size());
    return p.split(pos);
}

template<class Iterator>
void defIter(const char *name)
{
    class_<Iterator>(name, no_init)
        .def("__len__", &Iterator::__len__)
        .def("__iter__", &Iterator::__iter__)
        .def("next", &Iterator::next)
    ;
}

template<class Array>
PointIter<typename Array::const_iterator>
__iter__(const Array &a)
{
    return PointIter<typename Array::const_iterator>(
        a.begin(), a.end());
}

template<class Array>
PointIter<typename Array::const_reverse_iterator>
__reviter__(const Array &a)
{
    return PointIter<typename Array::const_reverse_iterator>(
        a.rbegin(), a.rend());
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

template<class Array>
list tangentList(const Array &p, int dx = 5, unsigned int skip = 0)
{
    if(p.size() < 2*dx + 2*skip)
    {
        PyErr_SetString(PyExc_ValueError,
            "tangentList: polygon too small (less than 2*dx + 2*skip points).");
        throw_error_already_set();
    }

    double pos = 0.0, prevSigma = 0.0;
    for(unsigned int i = 0; i < dx + skip - 1; ++i)
        pos += (p[i+1]-p[i]).magnitude();

    list result;

    for(unsigned int i = skip; i < p.size() - 2*dx - skip; ++i)
    {
        typename Array::value_type
            s(p[i+2*dx] - p[i]);

        pos += (p[i+dx] - p[i+dx-1]).magnitude();

        static const double eps = 1e-14;
        if(s.squaredMagnitude() < eps)
            continue;

        double sigma = atan2(s[1], s[0]), diff(sigma - prevSigma);
        if(diff < -M_PI)
            sigma += 2*M_PI;
        else if(diff > M_PI)
            sigma -= 2*M_PI;
        prevSigma = sigma;

        result.append(make_tuple(pos, sigma));
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
    if(extract<double>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthInternal<double>(
            arcLengthList, sigma, derivativeOrder);
    if(extract<Vector2>(arcLengthList[0][1]).check())
        return gaussianConvolveByArcLengthInternal<Vector2>(
            arcLengthList, sigma, derivativeOrder);
    PyErr_SetString(PyExc_TypeError,
        "gaussianConvolveByArcLength: arcLengthList can only be applied on lists\n"
        "of pairs containing floats or Vector2s as second values to be processed.");
    throw_error_already_set();
    return list(); // never reached
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
        "equidistantGaussians: arcLengthList can only be applied on lists\n"
        "of pairs containing floats or Vector2s as second values to be processed.");
    throw_error_already_set();
    return list(); // never reached
}

double fitParabola(const list &xyList)
{
    double sumAl4 = 0.0, sumAl3 = 0.0, sumAl2 = 0.0, sumAl = 0.0;
    double sumAl2Theta = 0.0, sumAlTheta = 0.0, sumTheta = 0.0;
    int count = 0;

    for(int i = 0; i < len(xyList); ++i)
    {
        double al = extract<double>(xyList[i][0])();
        double theta = extract<double>(xyList[i][1])();

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

        ++count;
    }

    vigra::Matrix<double> matrix(3, 3);
    matrix(0, 0) = sumAl4; matrix(0, 1) = sumAl3; matrix(0, 2) = sumAl2;
    matrix(1, 0) = sumAl3; matrix(1, 1) = sumAl2; matrix(1, 2) = sumAl;
    matrix(2, 0) = sumAl2; matrix(2, 1) = sumAl; matrix(2, 2) = count;

    vigra::Matrix<double> ergs(3, 1);
    ergs(0, 0) = sumAl2Theta;
    ergs(1, 0) = sumAlTheta;
    ergs(2, 0) = sumTheta;

    vigra::Matrix<double> result(3, 1);
    linearSolve(matrix, ergs, result);
    double p0 = result(0, 0);
    double p1 = result(1, 0);
    double p2 = result(2, 0);

    double error = 0.0;
    for(int i = 0; i < len(xyList); ++i)
    {
        double al = extract<double>(xyList[i][0])();
        double theta = extract<double>(xyList[i][1])();
        error += squaredNorm(p0*al*al + p1*al + p2 - theta);
    }
    error = sqrt(error / len(xyList));

    return error;
}

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

list sigmaOrbit(const QuadEdge *edge)
{
    const QuadEdge *orig(edge);
    list result;
    do
    {
        if(edge->dest().label())
        {
            int edgeLabel = edge->holderIndex();
            if(edge->isAnchor())
                result.append( edgeLabel);
            else
                result.append(-edgeLabel);
        }
        // this is confusing; I understood nextOrg() was the right
        // one, but apparently that would've been in a clockwise
        // manner (contrary to the documentation)..
        edge = edge->prevOrg();
    }
    while(edge != orig);
    return result;
}

tuple delaunay(const PointArray<Vector2> &points)
{
    // Construct a large surrounding triangle containing all points:
    Vector2 p1(-1e12, -1e8), p2(1e12, -1e8), p3(0.0, 3e16);
    Subdivision mesh(p1, p2, p3);

    list nodePositions, edges, orbits;
    nodePositions.append(object()); // node labels start with 1
    orbits.append(object());

    for(unsigned int i = 0; i < points.size(); ++i)
    {
        if(mesh.insertSite(points[i]) > 0)
            nodePositions.append(points[i]);
//         else
//             nodePositions.append(object());
    }

    for(Subdivision::NodeIterator it = mesh.nodesBegin();
        it != mesh.nodesEnd(); ++it)
    {
        orbits.append(object());
    }

    for(Subdivision::EdgeIterator it = mesh.edgesBegin();
        it != mesh.edgesEnd(); ++it)
    {
        if(!*it)
        {
            edges.append(object());
            continue;
        }

        const QuadEdge *edge((*it)->e);
        const DelaunayNode &o(edge->org()), &d(edge->dest());
        if(!o.label() || !d.label())
        {
            edges.append(object());
            continue;
        }

        edges.append(make_tuple(o.label(), d.label()));

        if(!orbits[o.label()])
            orbits[o.label()] = sigmaOrbit(edge);
        if(!orbits[d.label()])
            orbits[d.label()] = sigmaOrbit(edge->opposite());
    }

    return make_tuple(nodePositions, edges, orbits);
}

void defPolygon()
{
    def("intPos", &intPos);
    def("angleTheta", &angleTheta,
        "angleTheta(dy, dx)\n"
        "  calculates an efficient substitute value with the same sorting behavior\n"
        "  as the angle returned by math.atan2(dy, dx) (but ranged -2..2).\n"
        "  It does *not* give the angle, and seems to be only around 17% faster.");

    typedef BBoxPolygon<Vector2> PythonPolygon;

    typedef PythonPolygon::Base::Base Vector2Array;
    typedef PointIter<Vector2Array::const_iterator>
        VectorIter;
    defIter<VectorIter>("Vector2Iter");
    typedef PointIter<Vector2Array::const_reverse_iterator>
        VectorRevIter;
    defIter<VectorRevIter>("Vector2RevIter");
    class_<Vector2Array>("Vector2Array")
        .def(init<Vector2Array>())
        .def(init<list>())
        .def("reverse", &Vector2Array::reverse)
        .def("__len__", &Vector2Array::size)
        .def("__getitem__", &Array__getitem__<Vector2Array>)
        .def("__setitem__", &Array__setitem__<Vector2Array>)
        .def("__iter__", &__iter__<Vector2Array>)
        .def("__reviter__", &__reviter__<Vector2Array>)
        .def("insert", &insert<Vector2Array>)
        .def(self * double())
        .def(self + Vector2())
        .def("roundToInteger", &Vector2Array::roundToInteger)
    ;

    typedef PointArray<Point2D> Point2DArray;
    typedef PointIter<Point2DArray::const_iterator>
        Point2DIter;
    defIter<Point2DIter>("Point2DIter");
    typedef PointIter<Point2DArray::const_reverse_iterator>
        Point2DRevIter;
    defIter<Point2DRevIter>("Point2DRevIter");
    class_<Point2DArray>("Point2DArray")
        .def(init<Point2DArray>())
        .def("reverse", &Point2DArray::reverse)
        .def("__len__", &Point2DArray::size)
        .def("__getitem__", &Array__getitem__<Point2DArray>)
        .def("__setitem__", &Array__setitem__<Point2DArray>)
        .def("__iter__", &__iter__<Point2DArray>)
        .def("__reviter__", &__reviter__<Point2DArray>)
        .def("insert", &insert<Point2DArray>)
    ;

    class_<PythonPolygon, bases<Vector2Array> >("Polygon")
        .def(init<Vector2Array>())
        .def(init<list>()) // FIXME: use implicitly_convertible if possible
        .def("__delitem__", &erase<PythonPolygon>)
        .def("insert", &insert<PythonPolygon>)
        .def("append", &PythonPolygon::push_back)
        .def("extend", &PythonPolygon::extend)
        .def("__setitem__", &Polygon__setitem__<PythonPolygon>)
        .def("split", &split<PythonPolygon>)
        .def("length", &PythonPolygon::length)
        .def("partialArea", &PythonPolygon::partialArea)
        .def("boundingBox", &PythonPolygon::boundingBox)
        .def("contains", &PythonPolygon::contains)
        .def("swap", &PythonPolygon::swap)
        .def("reverse", &PythonPolygon::reverse)
        .def("nearestPoint", &PythonPolygon::nearestPoint)
        .def("invalidateProperties", &PythonPolygon::invalidateProperties)
        .def_pickle(ArrayPickleSuite<PythonPolygon>())
    ;

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
        .def("moveTo", &BoundingBox::moveTo, arg("newBegin"))
        .def("moveBy", &BoundingBox::moveBy, arg("offset"))
        .def("area", &BoundingBox::volume)
        .def("size", &BoundingBox::size)
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
    ;

    class_<Scanlines>("Scanlines", no_init)
        .def("__len__", &Scanlines::size)
        .def("__getitem__", &Scanlines__getitem__, return_internal_reference<>())
        .def("startIndex", &Scanlines::startIndex)
        .def("endIndex", &Scanlines::endIndex)
    ;

    class_<Scanlines::Scanline>("Scanline", no_init)
        .def("__len__", &Scanlines::Scanline::size)
        .def("__getitem__", &Array__getitem__byref<Scanlines::Scanline>,
             return_internal_reference<>())
    ;

    class_<ScanlineSegment>("ScanlineSegment", no_init)
        .def_readonly("begin", &ScanlineSegment::begin)
        .def_readonly("direction", &ScanlineSegment::direction)
        .def_readonly("end", &ScanlineSegment::end)
    ;

    def("scanPoly", &scanPoly<Vector2>,
        (arg("points"), arg("scanLineCount"), arg("startIndex") = 0),
        return_value_policy<manage_new_object>());
    def("fillScannedPoly", &pyFillScannedPoly);
    def("drawScannedPoly", &pyDrawScannedPoly);
    def("markEdgeInLabelImage", &markEdgeInLabelImage);
    def("removeEdgeFromLabelImage", &removeEdgeFromLabelImage);

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

    def("curvatureList", &curvatureList<Vector2Array>,
        (arg("pointArray"), arg("dx") = 5, arg("skipPoints") = 1),
        "curvatureList(pointArray, dx = 5, skipPoints = 1)\n"
        "calculates curvatures values for each triangle between point triples\n"
        "with indices (i-dx, i, i+dx), ignoring skipPoints points from both ends.\n"
        "returns a list of (arcLength, curvature) pairs,\n"
        "whose length is len(pointArray) - 2*dx - 2*skipPoints.");
    def("tangentList", &tangentList<Vector2Array>,
        (arg("pointArray"), arg("dx") = 5, arg("skipPoints") = 1),
        "tangentList(pointArray, dx = 5, skipPoints = 1)\n"
        "calculates tangent angles from each chord between the point with\n"
        "indices (i-dx, i+dx), ignoring skipPoints points from both ends.\n"
        "returns a list of (arcLength, angle) pairs,\n"
        "whose length is len(pointArray) - 2*dx - 2*skipPoints.");
    def("gaussianConvolveByArcLength", &gaussianConvolveByArcLength,
        (arg("arcLengthList"), arg("sigma"), arg("derivativeOrder") = 0));
    def("equidistantGaussians", &equidistantGaussians,
        args("edgeArcLengthList", "sigma", "distance"));
    // FIXME: document properly

    def("fitParabola", &fitParabola, args("xyList"),
        "fitParabola(xyList)\n"
        "fits a parabola to the (x,y) pairs in xyList and returns the\n"
        "mean squared error.\n");
    def("delaunay", &delaunay);
}
