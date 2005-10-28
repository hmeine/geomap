#include "cellimage_module.hxx"

using namespace vigra;
using namespace boost::python;
using namespace vigra::cellimage;

void defineNodes()
{
    class_<GeoMap::NodeIterator>("NodeIterator")
        .def("__iter__", (GeoMap::NodeIterator &(*)(GeoMap::NodeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (GeoMap::NodeInfo &(*)(GeoMap::NodeIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the GeoMap object (however, it's
             // lifetime is expected to be long enough):
             return_internal_reference<>());

    class_<NodeListProxy>("NodeList", no_init)
        .def("__len__", &NodeListProxy::__len__)
        .def("__getitem__", &NodeListProxy::__getitem__,
             // this is not really true, see NodeIterator::next
             return_internal_reference<>())
        .def("__iter__", &NodeListProxy::__iter__);
}
