#include "cellimage_module.hxx"
#include "cellstatistics.hxx"
#define private public
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
		.def("currentLevel", &FourEightPyramid::currentLevel)
		.def("currentSegmentation", &FourEightPyramid::currentSegmentation,
			 return_internal_reference<>())
		.def("currentCellStatistics", &FourEightPyramid::currentCellStatistics,
			 return_internal_reference<>())
		.def("__getitem__", &FourEightPyramid::segmentation,
			 return_internal_reference<>())
		.def("cutHead", &FourEightPyramid::cutHead)
		.def("__len__", &FourEightPyramid::levelCount)
		.def_readonly("history", &FourEightPyramid::history_));

	enum_<FourEightPyramid::OperationType>("OperationType")
		.value("RemoveIsolatedNode", FourEightPyramid::RemoveIsolatedNode)
		.value("MergeFaces", FourEightPyramid::MergeFaces)
		.value("RemoveBridge", FourEightPyramid::RemoveBridge)
		.value("MergeEdges", FourEightPyramid::MergeEdges)
		.value("RemoveEdge", FourEightPyramid::RemoveEdge)
		.value("RemoveEdgeWithEnds", FourEightPyramid::RemoveEdgeWithEnds);
	
	class_<FourEightPyramid::Operation>("Operation", no_init)
		.def_readonly("type", &FourEightPyramid::Operation::type_)
		.def_readonly("param", &FourEightPyramid::Operation::param_);

	class_<FourEightPyramid::History>("History", no_init)
		.def("__len__", &FourEightPyramid::History::size)
		.def("__getitem__", &historyGetItem,
			 return_internal_reference<>());
}
