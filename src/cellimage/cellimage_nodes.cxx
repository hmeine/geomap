#include "cellimage_module.hxx"

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

void defineNodes()
{
    class_<FourEightSegmentation::NodeIterator>("NodeIterator")
        .def("__iter__", (FourEightSegmentation::NodeIterator &(*)(FourEightSegmentation::NodeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (FourEightSegmentation::NodeInfo &(*)(FourEightSegmentation::NodeIterator &))&nextIterPos,
             // this is not really true, since the true owner would be
             // the FourEightSegmentation object (however, it's
             // lifetime is expected to be long enough:
             return_internal_reference<>());

    class_<NodeListProxy>("NodeList", no_init)
        .def("__len__", &NodeListProxy::__len__)
        .def("__getitem__", &NodeListProxy::__getitem__,
             // this is not really true, see NodeIterator::next
             return_internal_reference<>())
        .def("__iter__", &NodeListProxy::__iter__);
}
