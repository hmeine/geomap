#define PY_ARRAY_UNIQUE_SYMBOL geomap_PyArray_API
//#define NO_IMPORT_ARRAY
#include <vigra/numpy_array.hxx>

/********************************************************************/
/*                                                                  */
/*                Python export code & API wrappers                 */
/*                                                                  */
/********************************************************************/

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

void defDiff2D();
void defPolygon();
void defMap();
void defMapStats();
void defMapUtils();
void defDSL();
void defStatistics();
void defVectorConverters();

BOOST_PYTHON_MODULE_INIT(geomap)
{
    import_array();
    defDiff2D();
    defPolygon();
    defMap();
    defMapStats();
    defMapUtils();
    defDSL();
    defStatistics();
    defVectorConverters();
}
