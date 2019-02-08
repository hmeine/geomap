/************************************************************************/
/*                                                                      */
/*               Copyright 2007-2019 by Hans Meine                      */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

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
