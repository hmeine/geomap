#include "cellimage_module.hxx"
#include "cellstatistics.hxx"
#define protected public
#include "cellpyramid.hxx"

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

typedef vigra::CellPyramid<FourEightSegmentation, CellStatistics>
    FourEightPyramid;

using namespace boost::python;

FourEightPyramid::Operation &historyGetItem(FourEightPyramid::History &history,
											long index)
{
	if((index >= (long)history.size()) || (index < -(long)history.size()))
	{
		PyErr_SetObject(PyExc_IndexError, vigra::ownedPyObject(index));
		python::throw_error_already_set();
	}
	if(index < 0)
		index += history.size();
	return history[index];
}

void definePyramid()
{
	scope fourEightPyramidScope(
	class_<FourEightPyramid>("FourEightPyramid", no_init)
		//init<const FourEightSegmentation&>())
		.def("storeCheckpoint", &FourEightPyramid::storeCheckpoint)
		.def("removeIsolatedNode", &FourEightPyramid::removeIsolatedNode,
			 return_internal_reference<>())
		.def("mergeFaces", &FourEightPyramid::mergeFaces,
			 return_internal_reference<>())
		.def("removeBridge", &FourEightPyramid::removeBridge,
			 return_internal_reference<>())
		.def("mergeEdges", &FourEightPyramid::mergeEdges,
			 return_internal_reference<>())
		.def("removeEdge", &FourEightPyramid::removeEdge,
			 return_internal_reference<>())
		.def("removeEdgeWithEnds", &FourEightPyramid::removeEdgeWithEnds,
			 return_internal_reference<>())
		.def("topLevel", (FourEightPyramid::Level &(FourEightPyramid::*)(void))
			 &FourEightPyramid::topLevel, return_internal_reference<>())
		.def("__get__", (FourEightPyramid::Level *(FourEightPyramid::*)(unsigned int))
			 &FourEightPyramid::getLevel, return_value_policy<manage_new_object>())
		.def("__len__", &FourEightPyramid::levelCount)
		.def("cutAbove", (void(FourEightPyramid::*)(const FourEightPyramid::Level &))
			 &FourEightPyramid::cutAbove)
		.def("cutAbove", (void(FourEightPyramid::*)(unsigned int))
			 &FourEightPyramid::cutAbove)
		.def_readonly("history", &FourEightPyramid::history_));

	class_<FourEightPyramid::Level>("Level", no_init)
		.def("index", &FourEightPyramid::Level::index)
		.def("segmentation", &FourEightPyramid::Level::segmentation,
			 return_internal_reference<>())
		.def("cellStatistics", &FourEightPyramid::Level::cellStatistics,
			 return_internal_reference<>())
		.def("approachLevel", &FourEightPyramid::Level::approachLevel)
		.def("gotoLevel", &FourEightPyramid::Level::gotoLevel)
		.def("cutHead", &FourEightPyramid::Level::cutHead);

	enum_<FourEightPyramid::OperationType>("OperationType")
		.value("RemoveIsolatedNode", FourEightPyramid::RemoveIsolatedNode)
		.value("MergeFaces", FourEightPyramid::MergeFaces)
		.value("RemoveBridge", FourEightPyramid::RemoveBridge)
		.value("MergeEdges", FourEightPyramid::MergeEdges)
		.value("RemoveEdge", FourEightPyramid::RemoveEdge)
		.value("RemoveEdgeWithEnds", FourEightPyramid::RemoveEdgeWithEnds);

	class_<FourEightPyramid::Operation>("Operation", no_init)
		.def_readonly("type", &FourEightPyramid::Operation::type)
		.def_readonly("param", &FourEightPyramid::Operation::param);

	class_<FourEightPyramid::History>("History", no_init)
		.def("__len__", &FourEightPyramid::History::size)
		.def("__getitem__", &historyGetItem,
			 return_internal_reference<>());
}
