#include "foureightsegmentation.hxx"
#include <vigra/pythoncore.hxx>

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

void defineDartTraverser()
{
    class_<FourEightSegmentation::DartTraverser>("DartTraverser", no_init)
        .def("nextAlpha", &FourEightSegmentation::DartTraverser::nextAlpha, return_internal_reference<>())
        .def("prevAlpha", &FourEightSegmentation::DartTraverser::prevAlpha, return_internal_reference<>())
        .def("nextPhi", &FourEightSegmentation::DartTraverser::nextPhi, return_internal_reference<>())
        .def("prevPhi", &FourEightSegmentation::DartTraverser::prevPhi, return_internal_reference<>())
        .def("nextSigma", &FourEightSegmentation::DartTraverser::nextSigma, return_internal_reference<>())
        .def("prevSigma", &FourEightSegmentation::DartTraverser::prevSigma, return_internal_reference<>())
        .def("isSingular", &FourEightSegmentation::DartTraverser::isSingular)
		.def("recheckSingularity", &FourEightSegmentation::DartTraverser::recheckSingularity)
        .def("startNodeLabel", &FourEightSegmentation::DartTraverser::startNodeLabel)
        .def("endNodeLabel", &FourEightSegmentation::DartTraverser::endNodeLabel)
        .def("edgeLabel", &FourEightSegmentation::DartTraverser::edgeLabel)
        .def("leftFaceLabel", &FourEightSegmentation::DartTraverser::leftFaceLabel)
        .def("rightFaceLabel", &FourEightSegmentation::DartTraverser::rightFaceLabel)
        .def("startNode", &FourEightSegmentation::DartTraverser::startNode, return_internal_reference<>())
        .def("endNode", &FourEightSegmentation::DartTraverser::endNode, return_internal_reference<>())
        .def("edge", &FourEightSegmentation::DartTraverser::edge, return_internal_reference<>())
        .def("leftFace", &FourEightSegmentation::DartTraverser::leftFace, return_internal_reference<>())
        .def("rightFace", &FourEightSegmentation::DartTraverser::rightFace, return_internal_reference<>())
        .def(self == self)
        .def(self != self);
}

FourEightSegmentation::DartTraverser &contourGetItem(
	std::vector<FourEightSegmentation::DartTraverser> &contours, long index)
{
	if(index >= (long)contours.size())
	{
		PyErr_SetObject(PyExc_IndexError, vigra::ownedPyObject(index));
		python::throw_error_already_set();
	}
	return contours[index];
}

void defineCellInfos()
{
	class_<std::vector<FourEightSegmentation::DartTraverser> >("Contours", no_init)
		.def("__getitem__", &contourGetItem,
             return_internal_reference<>())
		.def("__len__", &std::vector<FourEightSegmentation::DartTraverser>::size);

    class_<FourEightSegmentation::CellInfo>("CellInfo", no_init)
        .def_readonly("label", &FourEightSegmentation::CellInfo::label)
        .def_readwrite("size", &FourEightSegmentation::NodeInfo::size)
        .def("initialized", &FourEightSegmentation::CellInfo::initialized)
        .def("uninitialize", &FourEightSegmentation::CellInfo::uninitialize);

    class_<FourEightSegmentation::NodeInfo,
           bases<FourEightSegmentation::CellInfo> >("NodeInfo", no_init)
        .def_readwrite("anchor", &FourEightSegmentation::NodeInfo::anchor)
        .def_readonly("centerX", &FourEightSegmentation::NodeInfo::centerX)
        .def_readonly("centerY", &FourEightSegmentation::NodeInfo::centerY);

    class_<FourEightSegmentation::EdgeInfo,
           bases<FourEightSegmentation::CellInfo> >("EdgeInfo", no_init)
        .def_readwrite("start", &FourEightSegmentation::EdgeInfo::start)
        .def_readwrite("end", &FourEightSegmentation::EdgeInfo::end);

    class_<FourEightSegmentation::FaceInfo,
		   bases<FourEightSegmentation::CellInfo> >("FaceInfo", no_init)
		.def_readonly("contours", &FourEightSegmentation::FaceInfo::contours);
}
