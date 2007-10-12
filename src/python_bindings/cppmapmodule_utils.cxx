#include "cppmap_utils.hxx"

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

void defMapUtils()
{
    using namespace boost::python;

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
}
