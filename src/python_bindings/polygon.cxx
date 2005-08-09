#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "vigra/pythonimage.hxx"
#include "vigra/polygon.hxx"

using namespace vigra;
using namespace boost::python;

// class PythonPolygon : public Polygon<Pixel>
// {
//     PythonPolygon(

// };

typedef Polygon<Pixel> PythonPolygon;

Point2D intPos(const Pixel &p)
{
    return Point2D((int)(p[0]+0.5f), (int)(p[1]+0.5f));
}

template<class Polygon>
typename Polygon::value_type
Polygon__getitem__(Polygon const & p, int i)
{
    if(i < 0)
        i += p.size();
    if(i < 0 || i >= p.size())
    {
        PyErr_SetString(PyExc_IndexError,
            "Polygon.__getitem__(): index out of bounds.");
        throw_error_already_set();
    }
    return p[i];
}

void defPolygon()
{
    def("intPos", &intPos);

    typedef PythonPolygon::Base FPointArray;
    class_<FPointArray>("VectorArray")
        .def(init<FPointArray>())
        .def(vector_indexing_suite<FPointArray>())
        .def("reverse", &FPointArray::reverse)
        .def(self * double())
        .def("roundToInteger", &FPointArray::roundToInteger)
    ;

    typedef PointArray<Point2D> IPointArray;
    class_<IPointArray>("Point2DArray")
        .def(init<IPointArray>())
        .def(vector_indexing_suite<IPointArray>())
    ;

    class_<PythonPolygon, bases<PythonPolygon::Base> >("Polygon")
        .def(init<list>())
        .def(init<FPointArray>())
        .def("__reviter__", &PythonPolygon::__reviter__)
        .def("__getitem__", &Polygon__getitem__<PythonPolygon>)
        .def("length", &PythonPolygon::length)
        .def("partialArea", &PythonPolygon::partialArea)
        .def("invalidateProperties", &PythonPolygon::invalidateProperties)
    ;

    class_<PythonPolygon::RevPointIter>("RevPointIter", no_init)
        .def("__len__", &PythonPolygon::RevPointIter::__len__)
        .def("__iter__", &PythonPolygon::RevPointIter::__iter__)
        .def("next", &PythonPolygon::RevPointIter::next,
             return_internal_reference<>())
    ;
}
