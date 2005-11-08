#include "foureightsegmentation.hxx"
#include <vigra/pythonutil.hxx>

using namespace vigra;
using namespace boost::python;
using namespace vigra::cellimage;

class DartIterator
{
public:
    DartIterator(const GeoMap::DartTraverser &dart)
    : it_(dart.neighborCirculator()),
      cellsUL_(dart.segmentation()->cells)
    {}

    DartIterator __iter__()
    {
        return *this;
    }

    Diff2D next()
    {
        if(it_.atEnd())
        {
            PyErr_SetString(PyExc_StopIteration, "");
            boost::python::throw_error_already_set();
        }
        Diff2D result(it_.neighborCirculator().base() - cellsUL_);
        ++it_;
        return result;
    }

    EdgelIterator it_;
    CellImage::traverser cellsUL_;
};

DartIterator DartTraverser__iter__(const GeoMap::DartTraverser &dart)
{
    return DartIterator(dart);
}

template <class RET, RET (GeoMap::DartTraverser::* FCT)()>
RET pyTraverserFct(GeoMap::DartTraverser & self)
{
    return (self.*FCT)();
}

template <class RET, RET (GeoMap::DartTraverser::* FCT)() const>
RET pyTraverserConstFct(GeoMap::DartTraverser const & self)
{
    return (self.*FCT)();
}

void defineDartTraverser()
{
    typedef GeoMap::DartTraverser DT;
    
    class_<DT> scope(
        "DartTraverser", init<const DT &>());
#ifndef _MSC_VER
    scope
        .def("nextAlpha", &DT::nextAlpha, return_internal_reference<>());
        .def("prevAlpha", &DT::prevAlpha, return_internal_reference<>())
        .def("nextPhi", &DT::nextPhi, return_internal_reference<>())
        .def("prevPhi", &DT::prevPhi, return_internal_reference<>())
        .def("nextSigma", &DT::nextSigma, return_internal_reference<>())
        .def("prevSigma", &DT::prevSigma, return_internal_reference<>())
        .def("isSingular", &DT::isSingular)
		//.def("carefulNextSigma", &DT::carefulNextSigma)
		//.def("carefulPrevSigma", &DT::carefulPrevSigma)
        .def("startNodeLabel", &DT::startNodeLabel)
        .def("endNodeLabel", &DT::endNodeLabel)
        .def("edgeLabel", &DT::edgeLabel)
        .def("leftFaceLabel", &DT::leftFaceLabel)
        .def("rightFaceLabel", &DT::rightFaceLabel)
        .def("startNode", &DT::startNode, return_internal_reference<>())
        .def("endNode", &DT::endNode, return_internal_reference<>())
        .def("edge", &DT::edge, return_internal_reference<>())
        .def("leftFace", &DT::leftFace, return_internal_reference<>())
        .def("rightFace", &DT::rightFace, return_internal_reference<>())
        .def("__iter__", &DartTraverser__iter__)
        .def(self == self)
        .def(self != self);
#else  // _MSC_VER
    scope
        .def("nextAlpha", &pyTraverserFct<DT&, &DT::nextAlpha>, return_internal_reference<>())
        .def("prevAlpha", &pyTraverserFct<DT&, &DT::prevAlpha>, return_internal_reference<>())
        .def("nextPhi", &pyTraverserFct<DT&, &DT::nextPhi>, return_internal_reference<>())
        .def("prevPhi", &pyTraverserFct<DT&, &DT::prevPhi>, return_internal_reference<>())
        .def("nextSigma", &pyTraverserFct<DT&, &DT::nextSigma>, return_internal_reference<>())
        .def("prevSigma", &pyTraverserFct<DT&, &DT::prevSigma>, return_internal_reference<>())
        .def("isSingular", &pyTraverserConstFct<bool, &DT::isSingular>)
		//.def("carefulNextSigma", &DT::carefulNextSigma)
		//.def("carefulPrevSigma", &DT::carefulPrevSigma)
        .def("startNodeLabel", &pyTraverserConstFct<CellLabel, &DT::startNodeLabel>)
        .def("endNodeLabel", &pyTraverserConstFct<CellLabel, &DT::endNodeLabel>)
        .def("edgeLabel", &pyTraverserConstFct<CellLabel, &DT::edgeLabel>)
        .def("leftFaceLabel", &pyTraverserConstFct<CellLabel, &DT::leftFaceLabel>)
        .def("rightFaceLabel", &pyTraverserConstFct<CellLabel, &DT::rightFaceLabel>)
        .def("startNode", &pyTraverserConstFct<GeoMap::NodeInfo &,&DT::startNode>, return_internal_reference<>())
        .def("endNode", &pyTraverserConstFct<GeoMap::NodeInfo &,&DT::endNode>, return_internal_reference<>())
        .def("edge", &pyTraverserConstFct<GeoMap::EdgeInfo &,&DT::edge>, return_internal_reference<>())
        .def("leftFace", &pyTraverserConstFct<GeoMap::FaceInfo &,&DT::leftFace>, return_internal_reference<>())
        .def("rightFace", &pyTraverserConstFct<GeoMap::FaceInfo &,&DT::rightFace>, return_internal_reference<>())
        .def("__iter__", &DartTraverser__iter__)
        .def(self == self)
        .def(self != self);
#endif

    class_<DartIterator>("DartIterator", init<const DT &>())
        .def("__iter__", &DartIterator::__iter__)
        .def("next", &DartIterator::next);

}

GeoMap::DartTraverser &contourGetItem(
	std::vector<GeoMap::DartTraverser> &contours, long index)
{
	if(index >= (long)contours.size())
	{
		PyErr_SetObject(PyExc_IndexError, incref(object(index).ptr()));
		throw_error_already_set();
	}
	return contours[index];
}

void defineCellInfos()
{
	class_<std::vector<GeoMap::DartTraverser> >("Contours", no_init)
		.def("__getitem__", &contourGetItem,
             return_internal_reference<>())
		.def("__len__", &std::vector<GeoMap::DartTraverser>::size);

    class_<GeoMap::CellInfo>("CellInfo", no_init)
        .def_readonly("label", &GeoMap::CellInfo::label)
        .def_readwrite("bounds", &GeoMap::NodeInfo::bounds)
        .def_readwrite("size", &GeoMap::NodeInfo::size)
        .def("initialized", &GeoMap::CellInfo::initialized)
        .def("uninitialize", &GeoMap::CellInfo::uninitialize);

    class_<GeoMap::NodeInfo,
           bases<GeoMap::CellInfo> >("NodeInfo", no_init)
        .def_readwrite("anchor", &GeoMap::NodeInfo::anchor)
        .def_readwrite("degree", &GeoMap::NodeInfo::degree);

    class_<GeoMap::EdgeInfo,
           bases<GeoMap::CellInfo> >("EdgeInfo", no_init)
        .def_readwrite("start", &GeoMap::EdgeInfo::start)
        .def_readwrite("end", &GeoMap::EdgeInfo::end);

    class_<GeoMap::FaceInfo,
		   bases<GeoMap::CellInfo> >("FaceInfo", no_init)
		.def_readonly("contours", &GeoMap::FaceInfo::contours);
}
