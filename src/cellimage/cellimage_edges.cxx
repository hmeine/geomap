#include "cellimage_module.hxx"

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

void defineEdges()
{
    class_<FourEightSegmentation::EdgeIterator>("EdgeIterator")
        .def("__iter__", (FourEightSegmentation::EdgeIterator &(*)(FourEightSegmentation::EdgeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (FourEightSegmentation::EdgeInfo &(*)(FourEightSegmentation::EdgeIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the FourEightSegmentation object (however, it's
             // lifetime is expected to be long enough):
             return_internal_reference<>());

    class_<EdgeListProxy>("EdgeList", no_init)
        .def("__len__", &EdgeListProxy::__len__)
        .def("__getitem__", &EdgeListProxy::__getitem__,
             // this is not really true, see EdgeIterator::next
             return_internal_reference<>())
        .def("__iter__", &EdgeListProxy::__iter__);
}
