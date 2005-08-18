#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include "vigra/pythonimage.hxx"
#include "vigra/polygon.hxx"

using namespace vigra;
using namespace boost::python;

// class PythonPolygon : public Polygon<Pixel>
// {
//     PythonPolygon(

// };

Point2D intPos(const Pixel &p)
{
    return Point2D((int)(p[0]+0.5f), (int)(p[1]+0.5f));
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

void defPolygon()
{
    def("intPos", &intPos);

    typedef Polygon<Pixel> PythonPolygon;

    typedef PythonPolygon::Base FPointArray;
    typedef PointIter<FPointArray::const_iterator>
        FPointIter;
    defIter<FPointIter>("FPointIter");
    typedef PointIter<FPointArray::const_reverse_iterator>
        FPointRevIter;
    defIter<FPointRevIter>("FPointRevIter");
    class_<FPointArray>("VectorArray")
        .def(init<FPointArray>())
        .def("reverse", &FPointArray::reverse)
        .def("__len__", &FPointArray::size)
        .def("__getitem__", &Array__getitem__<FPointArray>)
        .def("__setitem__", &Array__setitem__<FPointArray>)
        .def("__iter__", &__iter__<FPointArray>)
        .def("__reviter__", &__reviter__<FPointArray>)
        .def("insert", &insert<FPointArray>)
        .def(self * double())
        .def("roundToInteger", &FPointArray::roundToInteger)
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
        .def(init<FPointArray>())
        .def("insert", &insert<PythonPolygon>)
        .def("append", &PythonPolygon::push_back)
        .def("extend", &PythonPolygon::extend)
        .def("__setitem__", &Polygon__setitem__<PythonPolygon>)
        .def("split", &split<PythonPolygon>)
        .def("length", &PythonPolygon::length)
        .def("partialArea", &PythonPolygon::partialArea)
        .def("invalidateProperties", &PythonPolygon::invalidateProperties)
    ;

    def("simplifyPolygon",
        (PythonPolygon (*)(const PythonPolygon &,double))&simplifyPolygon);
    def("simplifyPolygon",
        (FPointArray (*)(const FPointArray &,double))&simplifyPolygon);
}
