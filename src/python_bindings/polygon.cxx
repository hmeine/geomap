#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include "vigra/pythonimage.hxx"
#include "vigra/polygon.hxx"

using namespace vigra;
using namespace boost::python;

Point2D intPos(const Vector2 &p)
{
    return Point2D((int)(p[0]+0.5), (int)(p[1]+0.5));
}

inline void checkPointIndex(int &i, unsigned int size)
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

unsigned int pyFillScannedPoly(
    const Scanlines &scanlines,
    PythonImage &targetV,
    GrayValue value)
    //const Pixel &value)
{
    PythonSingleBandImage target(targetV.subImage(0));
    return fillScannedPoly(scanlines, value,
                           target.traverser_begin(),
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
                           StandardValueAccessor<GrayValue>());
}

void markEdgeInLabelImage(
    const Scanlines &scanlines,
    PythonImage &labelVImage)
{
    PythonSingleBandImage labelImage(labelVImage.subImage(0));

    for(unsigned int y = scanlines.startIndex, i = 0;
        i < scanlines.size(); ++y, ++i)
    {
        const Scanlines::Scanline &scanline(scanlines[i]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            for(unsigned int x = scanline[j].begin;
                x < scanline[j].end; ++x)
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
    for(unsigned int y = scanlines.startIndex, i = 0;
        i < scanlines.size(); ++y, ++i)
    {
        const Scanlines::Scanline &scanline(scanlines[i]);
        for(unsigned int j = 0; j < scanline.size(); ++j)
        {
            for(unsigned int x = scanline[j].begin;
                x < scanline[j].end; ++x)
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

void defPolygon()
{
    def("intPos", &intPos);

    typedef Polygon<Vector2> PythonPolygon;

    typedef PythonPolygon::Base VectorArray;
    typedef PointIter<VectorArray::const_iterator>
        VectorIter;
    defIter<VectorIter>("VectorIter");
    typedef PointIter<VectorArray::const_reverse_iterator>
        VectorRevIter;
    defIter<VectorRevIter>("VectorRevIter");
    class_<VectorArray>("VectorArray")
        .def(init<VectorArray>())
        .def("reverse", &VectorArray::reverse)
        .def("__len__", &VectorArray::size)
        .def("__getitem__", &Array__getitem__<VectorArray>)
        .def("__setitem__", &Array__setitem__<VectorArray>)
        .def("__iter__", &__iter__<VectorArray>)
        .def("__reviter__", &__reviter__<VectorArray>)
        .def("insert", &insert<VectorArray>)
        .def(self * double())
        .def("roundToInteger", &VectorArray::roundToInteger)
    ;

    typedef PointArray<Point2D> IPointArray;
    typedef PointIter<IPointArray::const_iterator>
        IPointIter;
    defIter<IPointIter>("IPointIter");
    typedef PointIter<IPointArray::const_reverse_iterator>
        IPointRevIter;
    defIter<IPointRevIter>("IPointRevIter");
    class_<IPointArray>("Point2DArray")
        .def(init<IPointArray>())
        .def("reverse", &IPointArray::reverse)
        .def("__len__", &IPointArray::size)
        .def("__getitem__", &Array__getitem__<IPointArray>)
        .def("__setitem__", &Array__setitem__<IPointArray>)
        .def("__iter__", &__iter__<IPointArray>)
        .def("__reviter__", &__reviter__<IPointArray>)
        .def("insert", &insert<IPointArray>)
    ;

    class_<PythonPolygon, bases<PythonPolygon::Base> >("Polygon")
        .def(init<list>())
        .def(init<VectorArray>())
        .def("__delitem__", &erase<PythonPolygon>)
        .def("insert", &insert<PythonPolygon>)
        .def("append", &PythonPolygon::push_back)
        .def("extend", &PythonPolygon::extend)
        .def("__setitem__", &Polygon__setitem__<PythonPolygon>)
        .def("split", &split<PythonPolygon>)
        .def("length", &PythonPolygon::length)
        .def("partialArea", &PythonPolygon::partialArea)
        .def("reverse", &PythonPolygon::reverse)
        .def("invalidateProperties", &PythonPolygon::invalidateProperties)
        .def_pickle(ArrayPickleSuite<PythonPolygon>())
    ;

    class_<Scanlines>("Scanlines", no_init)
        .def("__len__", &Scanlines::size)
        .def("__getitem__", &Array__getitem__byref<Scanlines>,
             return_internal_reference<>())
        .def_readonly("startIndex", &Scanlines::startIndex)
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
        (VectorArray (*)(const VectorArray &,double))&simplifyPolygon);
    def("simplifyPolygon",
        (PythonPolygon (*)(const PythonPolygon &,double))&simplifyPolygon);
}
