#ifndef VIGRA_POSITIONEDMAP_HXX
#define VIGRA_POSITIONEDMAP_HXX

#include <boost/python.hpp>
#include <vigra/tinyvector.hxx>
#include "map2d.hxx"

class PositionedMap
{
    typedef vigra::TinyVector<double, 2> Vector2;
    typedef vigra::PositionedObject<Vector2, boost::python::object> ElementType;
    typedef vigra::Map2D<ElementType> MapType;

  public:
    void insert(const Vector2 &p, const boost::python::object &o)
    {
        objects_.insert(ElementType(p, o));
    }

    void remove(const Vector2 &p)
    {
        MapType::iterator nearest(
            objects_.nearest(ElementType(p, boost::python::object())));
        if(nearest == objects_.end() || nearest->second.position != p)
        {
            PyErr_SetString(PyExc_KeyError,
                            "PositionedMap.remove(): position not found");
            throw boost::python::error_already_set();
        }
        objects_.erase(nearest);
    }

    boost::python::object __call__(const Vector2 &p, double maxSquaredDist)
    {
        MapType::iterator nearest(
            objects_.nearest(ElementType(p, boost::python::object()),
                             maxSquaredDist));
        if(nearest == objects_.end())
            return boost::python::object();
        return nearest->second.payload;
    }

    MapType::size_type size() const
    {
        return objects_.size();
    }

  protected:
    MapType objects_;
};

#endif // VIGRA_POSITIONEDMAP_HXX
