#include "cellimage_module.hxx"
#include <vigra/pythonimage.hxx>
#include <boost/python/make_constructor.hpp>
#include <sstream>

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

vigra::cellimage::CellImage *
convertCellImage(
    vigra::PythonImage &image)
{
    vigra::cellimage::CellImage *result = new CellImage(image.size());
    vigra::PythonSingleBandImage simage(image.subImage(0));
    transformImage(srcImageRange(simage), destImage(*result),
                   vigra::cellimage::CellPixelSerializer());
    return result;
}   

vigra::cellimage::GeoMap *
createGeoMap(
    vigra::PythonImage &image,
    float boundaryValue,
	vigra::cellimage::CellType cornerType)
{
    vigra::PythonSingleBandImage simage(image.subImage(0));
    return new vigra::cellimage::GeoMap(
		srcImageRange(simage), boundaryValue, cornerType);
}

void validateDart(const vigra::cellimage::GeoMap::DartTraverser &dart)
{
    vigra_precondition(dart.neighborCirculator().center()->type() ==
                       vigra::cellimage::CellTypeVertex,
                       "dart is not attached to a node");
    vigra_precondition(dart.startNode().initialized(),
                       "dart's startNode is not valid (initialized())");
    if(!dart.isSingular())
        vigra_precondition(dart.edge().initialized(),
                           "dart's edge is not valid (initialized())");
}

std::string CellPixel__repr__(const vigra::cellimage::CellPixel p)
{
    std::stringstream s;
    s << "CellPixel(CellType.";
    switch(p.type())
    {
      case vigra::cellimage::CellTypeRegion:
          s << "Region";
          break;
      case vigra::cellimage::CellTypeLine:
          s << "Line";
          break;
      case vigra::cellimage::CellTypeVertex:
          s << "Vertex";
          break;
      case vigra::cellimage::CellTypeVertexOrLine:
          s << "VertexOrLine";
          break;
      case vigra::cellimage::CellTypeErrorOrLine:
          s << "ErrorOrLine";
          break;
      default:
          s << "Error";
    }
    s << ", " << p.label() << ")";
    return s.str();
}

using namespace boost::python;
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

    class_<CellPixel>("CellPixel", no_init)
        .def(init<CellType, CellLabel>())
        .add_property("type", &CellPixel::type,
                      (void(CellPixel::*)(CellType))&CellPixel::setType)
        .add_property("label", &CellPixel::label,
                      (void(CellPixel::*)(CellLabel))&CellPixel::setLabel)
        .def("__repr__", &CellPixel__repr__)
        .def(self == self);

    class_<CellImage>("CellImage", no_init)
        .def("__init__", make_constructor(
                 &convertCellImage, default_call_policies(),
                 args("image")))
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

    def("validateDart", &validateDart);
    def("debugDart", &debugDart);

    scope geoMap(
        class_<GeoMap>("GeoMap", init<const CellImage &>(args("importImage")))
        .def("__init__", make_constructor(
                 &createGeoMap, default_call_policies(),
                 (arg("image"), arg("edgeLabel"),
                  arg("cornerType") = CellTypeVertex)))
        .def("maxNodeLabel", &GeoMap::maxNodeLabel)
        .add_property("nodes", &NodeListProxy::create)
        .def("maxEdgeLabel", &GeoMap::maxEdgeLabel)
        .add_property("edges", &EdgeListProxy::create)
        .def("maxFaceLabel", &GeoMap::maxFaceLabel)
        .add_property("faces", &FaceListProxy::create)
        .def("removeIsolatedNode", &GeoMap::removeIsolatedNode,
             return_internal_reference<>())
        .def("mergeFaces", &GeoMap::mergeFaces,
             return_internal_reference<>())
        .def("removeBridge", &GeoMap::removeBridge,
             return_internal_reference<>())
        .def("mergeEdges", &GeoMap::mergeEdges,
             return_internal_reference<>())
        .def_readonly("cellImage", &GeoMap::cellImage));

    defineDartTraverser();
    defineCellInfos();
    defineNodes();
    defineEdges();
    defineFaces();
}
