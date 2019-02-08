/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

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
        .def("mirrorX", &DigitalStraightLine8::mirrorX)
        .def("mirrorXY", &DigitalStraightLine8::mirrorXY)
        .def("mirrorY", &DigitalStraightLine8::mirrorY)
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
        .def("mirrorX", &DigitalStraightLine4::mirrorX)
        .def("mirrorXY", &DigitalStraightLine4::mirrorXY)
        .def("mirrorY", &DigitalStraightLine4::mirrorY)
        .def_pickle(DSLPickleSuite<DigitalStraightLine4>())
    ;

    def("tangentDSL", &pyTangentDSL);
}
