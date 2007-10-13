#include "cppmap_utils.hxx"
#include "exporthelpers.hxx"
#include "labellut.hxx"

#include <boost/python.hpp>
namespace bp = boost::python;

struct EdgeProtectionPickleSuite : bp::pickle_suite
{
//     static bp::tuple getinitargs(EdgeProtection &ep)
//     {
//         return bp::make_tuple(
//             bp::handle<>(bp::reference_existing_object::apply<GeoMap *>
//                          ::type()(ep.map())));
//     }

    static bp::tuple getstate(EdgeProtection &ep)
    {
        return bp::make_tuple(ep.map());
//             bp::handle<>(bp::reference_existing_object::apply<GeoMap *>
//                          ::type()));
    }

    static void setstate(EdgeProtection &ep, bp::tuple state)
    {
        ep.attachHooks(bp::extract<boost::shared_ptr<GeoMap> >(state[0])());
    }
};

std::string EdgeProtection__repr__(EdgeProtection const &cb)
{
    std::stringstream s;
    s << "<EdgeProtection, ";
    if(cb.map())
        s << "active>";
    else
        s << "detached>";
    return s.str();
}

/********************************************************************/

#include <vigra/crackconnections.hxx>

vigra::PythonGrayImage
pyCrackConnectionImage(vigra::PythonImage const &labels)
{
    vigra::PythonGrayImage result(labels.size() + vigra::Diff2D(1, 1));
    crackConnectionImage(srcImageRange(labels), destImage(result));

    vigra::PythonGrayImage::traverser
        end = result.lowerRight() - vigra::Diff2D(1, 1),
        row = result.upperLeft();
    for(; row.y < end.y; ++row.y)
    {
        vigra::PythonGrayImage::traverser it = row;
        for(; it.x < end.x; ++it.x)
        {
            if((int)*it & 1)
                it[vigra::Diff2D(1, 0)] = (int)it[vigra::Diff2D(1, 0)] | 4;
            if((int)*it & 2)
                it[vigra::Diff2D(0, 1)] = (int)it[vigra::Diff2D(0, 1)] | 8;
        }
        if((int)*it & 2)
            it[vigra::Diff2D(0, 1)] = (int)it[vigra::Diff2D(0, 1)] | 8;
    }

    for(; row.x < end.x; ++row.x)
        if((int)*row & 1)
            row[vigra::Diff2D(1, 0)] = (int)row[vigra::Diff2D(1, 0)] | 4;

    return result;
}

/********************************************************************/

void defMapUtils()
{
    using namespace boost::python;
    
    {
        scope labelLUT(
            class_<LabelLUT>("LabelLUT", init<unsigned int>())
            .def("initIdentity", &LabelLUT::initIdentity)
            .def("appendOne", &LabelLUT::appendOne)
            .def("__getitem__", &Array__getitem__<LabelLUT>)
            .def("__len__", &LabelLUT::size)
            .def("relabel", &LabelLUT::relabel) // FIXME: check index
            .def("merged", &LabelLUT::mergedBegin)
            );
        
        RangeIterWrapper<LabelLUT::MergedIterator>("_MergedIterator");
    }

    class_<EdgeProtection, boost::noncopyable>(
        "EdgeProtection",
        "Protects GeoMap Edges which have a protection flag set.\n"
        "I.e. all operations that would remove an edge for which\n"
        "edge.flag(ALL_PROTECTION) is !0 will be\n"
        "canceled automatically.",
        init<boost::shared_ptr<GeoMap> >(arg("map")=object()))
        .def("detachHooks", &EdgeProtection::detachHooks)
        .def("map", &EdgeProtection::map)
        .def_pickle(EdgeProtectionPickleSuite())
        .def("__repr__", &EdgeProtection__repr__)
    ;

    register_ptr_to_python<boost::shared_ptr<GeoMap> >();

    def("mergeFacesCompletely", &mergeFacesCompletely,
        (arg("dart"), arg("mergeDegree2Nodes") = true),
        "mergeFacesCompletely(dart, mergeDegree2Nodes = true)\n\n"
        "In contrast to the Euler operation mergeFaces(), this function\n"
        "removes all common edges of the two faces, !only the single\n"
        "edge belonging to dart.\n\n"
        "Furthermore, if the optional parameter mergeDegree2Nodes is\n"
        "true (default), all nodes whose degree is reduced to two will be\n"
        "merged into their surrounding edges.\n\n"
        "Returns the surviving face.");

    def("crackConnectionImage", &pyCrackConnectionImage,
        args("labelImage"),
        "crackConnectionImage(labelImage)\n\n"
        "Tranform a region image into an image with crack connections marked.\n"
        "(Bit 1: connected to the right, bit 2: connected downwards)");
}
