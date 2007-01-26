#include "vigra/dsl.hxx"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <vector>

using namespace boost::python;
using namespace vigra;

template<class DSL>
struct DSLPickleSuite : pickle_suite
{
    static tuple getinitargs(DSL const& dsl)
    {
        return make_tuple(dsl.a(), dsl.b(), dsl.pos());
    }
};

tuple pyTangentDSL(list pyFreemanCodes, int index, bool closed)
{
    std::vector<unsigned char> freemanCodes(len(pyFreemanCodes));
    for(unsigned int i = 0; i < freemanCodes.size(); ++i)
        freemanCodes[i] = extract<unsigned char>(pyFreemanCodes[i])();
    
    DigitalStraightLine<int, true> result(0, 1, 0);
    int offset = tangentDSL(freemanCodes.begin(), freemanCodes.end(),
                            index, closed, result);
    return make_tuple(result, offset ? object(offset) : object());
}

void defDSL()
{
    enum_<DSL::LeaningType>("LeaningType")
        .value("CenterLine", DSL::CenterLine)
        .value("LowerLeaningLine", DSL::LowerLeaningLine)
        .value("UpperLeaningLine", DSL::UpperLeaningLine)
    ;

    typedef DigitalStraightLine<int, true> DigitalStraightLine8;

    class_<DigitalStraightLine8>("DigitalStraightLine8", init<int, int, int>())
        .def("eightConnected", &DigitalStraightLine8::eightConnected)
        .add_property("a", &DigitalStraightLine8::a, &DigitalStraightLine8::setA)
        .add_property("b", &DigitalStraightLine8::b, &DigitalStraightLine8::setB)
        .add_property("pos", &DigitalStraightLine8::pos, &DigitalStraightLine8::setPos)
        .def("__call__", &DigitalStraightLine8::operator())
        .def("width", &DigitalStraightLine8::width)
        .def("contains", &DigitalStraightLine8::contains)
        .def("slope", &DigitalStraightLine8::slope)
        .def("axisIntercept", &DigitalStraightLine8::axisIntercept,
             arg("leaningType") = DSL::CenterLine)
        .def("addPoint", &DigitalStraightLine8::addPoint)
        .def("convertToFourConnected", &DigitalStraightLine8::convertToFourConnected)
        .def_pickle(DSLPickleSuite<DigitalStraightLine8>())
    ;

    typedef DigitalStraightLine<int, false> DigitalStraightLine4;

    class_<DigitalStraightLine4>("DigitalStraightLine4", init<int, int, int>())
        .def("eightConnected", &DigitalStraightLine4::eightConnected)
        .add_property("a", &DigitalStraightLine4::a, &DigitalStraightLine4::setA)
        .add_property("b", &DigitalStraightLine4::b, &DigitalStraightLine4::setB)
        .add_property("pos", &DigitalStraightLine4::pos, &DigitalStraightLine4::setPos)
        .def("__call__", &DigitalStraightLine4::operator())
        .def("width", &DigitalStraightLine4::width)
        .def("contains", &DigitalStraightLine4::contains)
        .def("slope", &DigitalStraightLine4::slope)
        .def("axisIntercept", &DigitalStraightLine4::axisIntercept,
             arg("leaningType") = DSL::CenterLine)
        .def_pickle(DSLPickleSuite<DigitalStraightLine8>())
    ;

    def("tangentDSL", &pyTangentDSL);
}
