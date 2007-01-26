#include "vigra/dsl.hxx"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

void defDSL()
{
    using namespace boost::python;

    typedef vigra::DigitalStraightLine<int, true> DigitalStraightLine8;

    class_<DigitalStraightLine8>("DigitalStraightLine8", init<int, int, int>())
        .def("eightConnected", &DigitalStraightLine8::eightConnected)
        .def("a", &DigitalStraightLine8::a)
        .def("b", &DigitalStraightLine8::b)
        .def("pos", &DigitalStraightLine8::pos)
        .def("width", &DigitalStraightLine8::width)
        .def("contains", &DigitalStraightLine8::contains)
        .def("slope", &DigitalStraightLine8::slope)
        .def("axisIntercept", &DigitalStraightLine8::axisIntercept)
        .def("addPoint", &DigitalStraightLine8::addPoint)
        .def("convertToFourConnected", &DigitalStraightLine8::convertToFourConnected)
    ;

    typedef vigra::DigitalStraightLine<int, false> DigitalStraightLine4;

    class_<DigitalStraightLine4>("DigitalStraightLine4", init<int, int, int>())
        .def("eightConnected", &DigitalStraightLine4::eightConnected)
        .def("a", &DigitalStraightLine4::a)
        .def("b", &DigitalStraightLine4::b)
        .def("pos", &DigitalStraightLine4::pos)
        .def("width", &DigitalStraightLine4::width)
        .def("contains", &DigitalStraightLine4::contains)
        .def("slope", &DigitalStraightLine4::slope)
        .def("axisIntercept", &DigitalStraightLine4::axisIntercept)
        .def("addPoint", &DigitalStraightLine4::addPoint)
        //.def("convertToFourConnected", &DigitalStraightLine4::convertToFourConnected)
    ;
}
