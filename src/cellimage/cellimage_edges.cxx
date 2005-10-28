#include "cellimage_module.hxx"

using namespace vigra;
using namespace boost::python;
using namespace vigra::cellimage;

void defineEdges()
{
    class_<GeoMap::EdgeIterator>("EdgeIterator")
        .def("__iter__", (GeoMap::EdgeIterator &(*)(GeoMap::EdgeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (GeoMap::EdgeInfo &(*)(GeoMap::EdgeIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the GeoMap object (however, it's
             // lifetime is expected to be long enough):
             return_internal_reference<>());

    class_<EdgeListProxy>("EdgeList", no_init)
        .def("__len__", &EdgeListProxy::__len__)
        .def("__getitem__", &EdgeListProxy::__getitem__,
             // this is not really true, see EdgeIterator::next
             return_internal_reference<>())
        .def("__iter__", &EdgeListProxy::__iter__);
}
