#include "cellimage_module.hxx"

using vigra::cellimage::CellImage;
using vigra::cellimage::CellPixel;

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

CellPixel getPixel(CellImage const & image, vigra::Diff2D const & i)
{
    vigra_precondition(i.x >= 0 && i.x < image.width() &&
                       i.y >= 0 && i.y < image.height(),
                       "coordinates out of range.");
    return image[i];
}

void setPixel(CellImage & image, vigra::Diff2D const & i, CellPixel const & value)
{
    vigra_precondition(i.x >= 0 && i.x < image.width() &&
                       i.y >= 0 && i.y < image.height(),
                       "coordinates out of range.");
    image[i] = value;
}

vigra::cellimage::FourEightSegmentation *
createFourEightSegmentation(
    vigra::PythonSingleBandImage &image,
    float boundaryValue,
	vigra::cellimage::CellType cornerType)
{
    return new vigra::cellimage::FourEightSegmentation(
		srcImageRange(image), boundaryValue, cornerType);
}

using namespace python;
using namespace vigra::cellimage;

void definePyramid();

void defineDartTraverser();
void defineCellInfos();
void defineNodes();
void defineEdges();
void defineFaces();

BOOST_PYTHON_MODULE_INIT(cellimage)
{
    enum_<CellType>("CellType")
        .value("Error", CellTypeError)
        .value("Region", CellTypeRegion)
        .value("Line", CellTypeLine)
        .value("Vertex", CellTypeVertex);

    class_<CellPixel>("CellPixel")
        .def(init<CellType, CellLabel>())
        .add_property("type", &CellPixel::type,
                      (void(CellPixel::*)(CellType))&CellPixel::setType)
        .add_property("label", &CellPixel::label,
                      (void(CellPixel::*)(CellLabel))&CellPixel::setLabel)
        .def(self == self);

    class_<CellImage>("CellImage")
        .def("width", &CellImage::width)
        .def("height", &CellImage::height)
        .def("size", &CellImage::size)
        .def("__getitem__", &getPixel)
        .def("__setitem__", &setPixel)
        .def("get", &getPixel)
        .def("set", &setPixel)
        .def("get", &getPixelXY)
        .def("set", &setPixelXY);

    definePyramid();

    def("debugDart", &debugDart);

    def("createFourEightSegmentation", createFourEightSegmentation,
        return_value_policy<manage_new_object>());

    scope fourEightSegmentation(
        class_<FourEightSegmentation>("FourEightSegmentation", no_init)
        .def("maxNodeLabel", &FourEightSegmentation::maxNodeLabel)
        .add_property("nodes", &NodeListProxy::create)
        .def("maxEdgeLabel", &FourEightSegmentation::maxEdgeLabel)
        .add_property("edges", &EdgeListProxy::create)
        .def("maxFaceLabel", &FourEightSegmentation::maxFaceLabel)
        .add_property("faces", &FaceListProxy::create)
        .def("removeIsolatedNode", &FourEightSegmentation::removeIsolatedNode,
             return_internal_reference<>())
        .def("mergeFaces", &FourEightSegmentation::mergeFaces,
             return_internal_reference<>())
        .def("removeBridge", &FourEightSegmentation::removeBridge,
             return_internal_reference<>())
        .def("mergeEdges", &FourEightSegmentation::mergeEdges,
             return_internal_reference<>()));

    defineDartTraverser();
    defineCellInfos();
    defineNodes();
    defineEdges();
    defineFaces();
}
