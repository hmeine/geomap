#include <vigra/vigrapython.hxx>
#include "foureightsegmentation.hxx"

namespace vigra { namespace CellImage {

CellPixel getPixelXY(CellImage const & image, int x, int y)
{
    vigra_precondition(x >= 0 && x < image.width() &&
                       y >= 0 && y < image.height(),
                       "coordinates out of range.");
    return image(x, y);
}

void setPixelXY(CellImage & image, CellPixel const & value, int x, int y)
{
    vigra_precondition(x >= 0 && x < image.width() &&
                       y >= 0 && y < image.height(),
                       "coordinates out of range.");
    image(x, y) = value;
}

CellPixel getPixel(CellImage const & image, Diff2D const & i)
{
    vigra_precondition(i.x >= 0 && i.x < image.width() &&
                       i.y >= 0 && i.y < image.height(),
                       "coordinates out of range.");
    return image[i];
}

void setPixel(CellImage & image, Diff2D const & i, CellPixel const & value)
{
    vigra_precondition(i.x >= 0 && i.x < image.width() &&
                       i.y >= 0 && i.y < image.height(),
                       "coordinates out of range.");
    image[i] = value;
}

struct NodeListProxy
{
    NodeListProxy(FourEightSegmentation *segmentation)
    : segmentation_(segmentation)
    {}
    
    static NodeListProxy create(FourEightSegmentation *segmentation)
    {
        return NodeListProxy(segmentation);
    }

    FourEightSegmentation *segmentation_;

    long __len__() const
    {
        return segmentation_->nodeCount();
    }
    
    FourEightSegmentation::NodeInfo &__getitem__(long index)
    {
        return segmentation_->node(index);
    }

    FourEightSegmentation::NodeIterator __iter__()
    {
        return segmentation_->nodesBegin();
    }
};

void initFourEightSegmentation(FourEightSegmentation &segmentation,
                               SingleBandImage &image)
{
    segmentation.init(srcImageRange(image));
}

} } // namespace vigra::CellImage

template<class T>
T &returnSelf(T &v) 
{
    return v;
}

template<class Iterator>
Iterator &nextIterPos(Iterator &v)
{
    if(!(++v).inRange())
    {
        PyErr_SetString(PyExc_StopIteration,
                        "no more nodes in cellcomplex");
        throw_error_already_set();
    }
    return v;
}

using namespace vigra;
using vigra::CellImage::FourEightSegmentation;
using vigra::CellImage::CellPixel;

BOOST_PYTHON_MODULE_INIT(cellimage)
{
	enum_<CellImage::CellType>("CellType")
		.value("Error", CellImage::CellTypeError)
		.value("Region", CellImage::CellTypeRegion)
        .value("Line", CellImage::CellTypeLine)
        .value("Vertex", CellImage::CellTypeVertex);

    class_<CellPixel>("CellPixel")
        .def(init<CellImage::CellType, CellImage::CellLabel>())
        .add_property("type", &CellPixel::type,
                      &CellPixel::setType)
        .def(self == self);

    class_<CellImage::CellImage>("CellImage")
        .def("width", &CellImage::CellImage::width)
        .def("height", &CellImage::CellImage::height)
        .def("size", &CellImage::CellImage::size)
        .def("__getitem__", &CellImage::getPixel)
        .def("__setitem__", &CellImage::setPixel)
        .def("get", &CellImage::getPixel)
        .def("set", &CellImage::setPixel)
        .def("get", &CellImage::getPixelXY)
        .def("set", &CellImage::setPixelXY);

    //scope fourEightSegmentation =
    class_<FourEightSegmentation>("FourEightSegmentation")
        .def("__init__", &CellImage::initFourEightSegmentation)
        .def("width", &FourEightSegmentation::width)
        .def("height", &FourEightSegmentation::height)
        .def("nodeCount", &FourEightSegmentation::nodeCount)
        .def("maxNodeLabel", &FourEightSegmentation::maxNodeLabel)
        .def("nodes", &CellImage::NodeListProxy::create);

    class_<FourEightSegmentation::NodeIterator>("NodeIterator")
        .def("__iter__", (FourEightSegmentation::NodeIterator &(*)(FourEightSegmentation::NodeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", (FourEightSegmentation::NodeIterator &(*)(FourEightSegmentation::NodeIterator &))&nextIterPos,
             return_internal_reference<>());

    class_<FourEightSegmentation::CellInfo>("CellInfo", no_init)
        .def_readonly("label", &FourEightSegmentation::CellInfo::label)
        .def("initialized", &FourEightSegmentation::CellInfo::initialized)
        .def("uninitialize", &FourEightSegmentation::CellInfo::uninitialize);

    class_<FourEightSegmentation::NodeInfo,
           bases<FourEightSegmentation::CellInfo> >("NodeInfo", no_init)
        .def_readonly("centerX", &FourEightSegmentation::NodeInfo::centerX)
        .def_readonly("centerY", &FourEightSegmentation::NodeInfo::centerY)
        .def_readwrite("size", &FourEightSegmentation::NodeInfo::size);

    class_<CellImage::NodeListProxy>("NodeList", no_init)
        .def("__len__", &CellImage::NodeListProxy::__len__)
        .def("__getitem__", &CellImage::NodeListProxy::__getitem__,
             // this is not really true, since the true owner would be
             // the FourEightSegmentation object (however, it's
             // lifetime is expected to be long enough:
             return_internal_reference<>())
        .def("__iter__", &CellImage::NodeListProxy::__iter__);
}
