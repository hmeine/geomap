#ifndef VIGRA_POSITIONEDMAP_HXX
#define VIGRA_POSITIONEDMAP_HXX

#include <boost/python.hpp>
#include <map>

namespace vigra {

template<class Position, class Payload>
class PositionedObject
{
  public:
    typedef typename Position::value_type value_type;

    Position position;
    Payload  payload;

    PositionedObject(const Position &po, const Payload &pl)
    : position(po),
      payload(pl)
    {}

    Position operator-(const PositionedObject &other) const
    {
        return position - other.position;
    }

    value_type operator[](unsigned char index) const
    {
        return position[index];
    }
};

template<class Vector2D>
class Map2D
{
public:
    typedef typename Vector2D::value_type      CoordType;
    typedef std::multimap<CoordType, Vector2D> CoordMap;
    typedef typename CoordMap::iterator        iterator;
    typedef typename CoordMap::const_iterator  const_iterator;
    typedef typename CoordMap::size_type       size_type;

    void insert(const Vector2D &vector2D)
    {
        vectors_.insert(
            typename CoordMap::value_type((vector2D)[0], vector2D));
    }

    void erase(const iterator &it)
    {
        vectors_.erase(it);
    }

        // assumes that Vector2D(x, y) is a valid constructor:
    void insert(CoordType x, CoordType y)
    {
        vectors_.insert(
            typename CoordMap::value_type(x, Vector2D(x, y)));
    }

    template<class Vector2DIterator>
    void fillFrom(const Vector2DIterator &begin, const Vector2DIterator &end)
    {
        vectors_.clear();
        for(Vector2DIterator it = begin; it != end; ++it)
        {
            vectors_.insert(typename CoordMap::value_type((*it)[0], *it));
        }
    }

    const_iterator begin() const
    {
        return vectors_.begin();
    }

    const_iterator end() const
    {
        return vectors_.end();
    }

    iterator begin()
    {
        return vectors_.begin();
    }

    iterator end()
    {
        return vectors_.end();
    }

    typename CoordMap::size_type size() const
    {
        return vectors_.size();
    }

    const_iterator nearest(
        const Vector2D &v, double maxSquaredDist = NumericTraits<double>::max()) const
    {
        return search(vectors_.begin(), vectors_.lower_bound(v[0]), vectors_.end(),
                      v, maxSquaredDist);
    }

    iterator nearest(
        const Vector2D &v, double maxSquaredDist = NumericTraits<double>::max())
    {
        return search(vectors_.begin(), vectors_.lower_bound(v[0]), vectors_.end(),
                      v, maxSquaredDist);
    }

protected:
    template<class ITERATOR>
    static ITERATOR search(
        const ITERATOR &begin, ITERATOR midPos, const ITERATOR &end,
        const Vector2D &v, double maxSquaredDist)
    {
        ITERATOR nearestPos(end);

        for(ITERATOR it = midPos; it != end; ++it)
        {
            if(squaredNorm(it->first - v[0]) > maxSquaredDist)
                break;

            double dist2 = squaredNorm(it->second - v);
            if(dist2 < maxSquaredDist)
            {
                nearestPos = it;
                maxSquaredDist = dist2;
            }
        }

        if(midPos == begin)
            return nearestPos;

        for(ITERATOR it = --midPos; true; --it)
        {
            if(squaredNorm(v[0] - it->first) > maxSquaredDist)
                break;

            double dist2 = squaredNorm(it->second - v);
            if(dist2 < maxSquaredDist)
            {
                nearestPos = it;
                maxSquaredDist = dist2;
            }

            if(it == begin)
                break;
        }

        return nearestPos;
    }

    CoordMap vectors_;
};

} // namespace vigra

/********************************************************************/
/*                                                                  */
/*                       boost::python stuff                        */
/*                                                                  */
/********************************************************************/

#include <vigra/pythonimage.hxx> // for Vector2 typedef

class PositionedMap
{
    typedef vigra::Vector2 Vector2;
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
