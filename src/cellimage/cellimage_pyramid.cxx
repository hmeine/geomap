#include "cellimage_module.hxx"
#include "cellstatistics.hxx"
#include "cellpyramid.hxx"

using namespace vigra;
using namespace python;
using namespace vigra::cellimage;

typedef vigra::CellPyramid<FourEightSegmentation, CellStatistics>
    FourEightPyramid;

using namespace boost::python;

void definePyramid()
{
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
		.def("segmentation", &FourEightPyramid::segmentation,
			 return_internal_reference<>())
		.def("cutHead", &FourEightPyramid::cutHead)
		.def("levelCount", &FourEightPyramid::levelCount);
}
