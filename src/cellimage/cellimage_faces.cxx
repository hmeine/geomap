#include "cellimage_module.hxx"

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

void defineFaces()
{
    class_<FourEightSegmentation::FaceIterator>("FaceIterator")
        .def("__iter__", (FourEightSegmentation::FaceIterator &(*)(FourEightSegmentation::FaceIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (FourEightSegmentation::FaceInfo &(*)(FourEightSegmentation::FaceIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the FourEightSegmentation object (however, it's
             // lifetime is expected to be long enough:
             return_internal_reference<>());

    class_<FaceListProxy>("FaceList", no_init)
        .def("__len__", &FaceListProxy::__len__)
        .def("__getitem__", &FaceListProxy::__getitem__,
             // this is not really true, see FaceIterator::next
             return_internal_reference<>())
        .def("__iter__", &FaceListProxy::__iter__);
}
