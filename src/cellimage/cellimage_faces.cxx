#include "cellimage_module.hxx"

using namespace vigra;
using namespace boost::python;
using namespace vigra::cellimage;

void defineFaces()
{
    class_<GeoMap::FaceIterator>("FaceIterator")
        .def("__iter__", (GeoMap::FaceIterator &(*)(GeoMap::FaceIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (GeoMap::FaceInfo &(*)(GeoMap::FaceIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the GeoMap object (however, it's
             // lifetime is expected to be long enough:
             return_internal_reference<>());

    class_<FaceListProxy>("FaceList", no_init)
        .def("__len__", &FaceListProxy::__len__)
        .def("__getitem__", &FaceListProxy::__getitem__,
             // this is not really true, see FaceIterator::next
             return_internal_reference<>())
        .def("__iter__", &FaceListProxy::__iter__);
}
