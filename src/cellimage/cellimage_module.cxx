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

void createFourEightSegmentation(
    vigra::cellimage::FourEightSegmentation *self,
    vigra::SingleBandImage &image)
{
    new (self) vigra::cellimage::FourEightSegmentation(srcImageRange(image));
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
                      &CellPixel::setType)
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

    scope fourEightSegmentation(
		class_<FourEightSegmentation>("FourEightSegmentation", no_init)
        .def("nodeCount", &FourEightSegmentation::nodeCount)
        .def("maxNodeLabel", &FourEightSegmentation::maxNodeLabel)
        .def("nodes", &NodeListProxy::create)
        .def("edgeCount", &FourEightSegmentation::edgeCount)
        .def("maxEdgeLabel", &FourEightSegmentation::maxEdgeLabel)
        .def("edges", &EdgeListProxy::create)
        .def("faceCount", &FourEightSegmentation::faceCount)
        .def("maxFaceLabel", &FourEightSegmentation::maxFaceLabel)
        .def("faces", &FaceListProxy::create)
        .def("removeIsolatedNode", &FourEightSegmentation::removeIsolatedNode,
             return_internal_reference<>())
        .def("mergeFaces", &FourEightSegmentation::mergeFaces,
             return_internal_reference<>())
        .def("removeBridge", &FourEightSegmentation::removeBridge,
             return_internal_reference<>()));

    def("createFourEightSegmentation", createFourEightSegmentation,
        return_value_policy<manage_new_object>());

    defineDartTraverser();
    defineCellInfos();
    defineNodes();
    defineEdges();
    defineFaces();
}
