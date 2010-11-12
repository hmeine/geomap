#pragma warning( disable : 4786 )

#include <boost/python.hpp>

#include <vigra/diff2d.hxx>

#include <string>
#include <sstream>

namespace python = boost::python;

namespace vigra {

#define defineDiff2Str(type) \
std::string diff2Str##type(type const & d) \
{ \
    std::stringstream s; \
    s << #type "(" << d.x << ", " << d.y << ")"; \
    return s.str(); \
}

defineDiff2Str(Diff2D)
defineDiff2Str(Size2D)
defineDiff2Str(Point2D)

std::string rect2D__str__(Rect2D const & r)
{
    std::stringstream s;
    s << "Rect2D<"
      << r.upperLeft() << "-" << r.lowerRight()
      << "=" << r.size() << ">";
    return s.str();
}

std::string rect2D__repr__(Rect2D const & r)
{
    std::stringstream s;
    s << "Rect2D("
      << diff2StrPoint2D(r.upperLeft()) << ", "
      << diff2StrSize2D(r.size()) << ")";
    return s.str();
}

template<class DIFF2D>
struct Diff2DPickleSuite : python::pickle_suite
{
    static python::tuple getinitargs(DIFF2D const& v)
    {
        return python::make_tuple(v.x, v.y);
    }
};

struct Rect2DPickleSuite : python::pickle_suite
{
    static python::tuple getinitargs(Rect2D const& r)
    {
        return python::make_tuple(r.upperLeft(), r.lowerRight());
    }
};

template <class Shape>
struct ShapeFromPython
{
    typedef Shape Type;

    ShapeFromPython()
    {
        python::converter::registry::insert(&convertible, &construct, python::type_id<Type>());
    }

    static void* convertible(PyObject* obj)
    {
        return ((PyTuple_Check(obj) || PyList_Check(obj)) &&
                PySequence_Size(obj) == 2 &&
                python::extract<double>(PySequence_Fast_GET_ITEM(obj, 0)).check() &&
                python::extract<double>(PySequence_Fast_GET_ITEM(obj, 1)).check())
                ? obj
                : 0;
    }

    static void construct(PyObject* obj, python::converter::rvalue_from_python_stage1_data* data)
    {
        double x = python::extract<double>(PySequence_Fast_GET_ITEM(obj, 0))();
        double y = python::extract<double>(PySequence_Fast_GET_ITEM(obj, 1))();
        void* const storage = ((python::converter::rvalue_from_python_storage<Type>*)data)->storage.bytes;
        new (storage) Type((int)(x < 0
                                   ? x - 0.5
                                   : x + 0.5),
                           (int)(y < 0
                                   ? y - 0.5
                                   : y + 0.5));
        data->convertible = storage;
    }
};

int Diff2D__getitem__(Diff2D const & self, int i)
{
    switch(i)
    {
        case 0:
            return self.x;
        case 1:
            return self.y;
        default:
            PyErr_SetString(PyExc_IndexError, "Diff2D.__getitem__(): index out of bounds.");
            python::throw_error_already_set();
            return 0; // unreachable
    }
}

void Diff2D__setitem__(Diff2D & self, int i, int v)
{
    switch(i)
    {
        case 0:
            self.x = v;
            break;
        case 1:
            self.y = v;
            break;
        default:
            PyErr_SetString(PyExc_IndexError, "Diff2D.__setitem__(): index out of bounds.");
            python::throw_error_already_set();
    }
}

int Diff2D__len__(Diff2D const & self)
{
    return 2;
}

class MeshGridIterator
{
    Point2D point_;
    Size2D  stride_;
    Rect2D  rect_;

  public:
    MeshGridIterator(Size2D const & s, Size2D const & stride = Size2D(1,1))
    : point_(-1,0),
      stride_(stride),
      rect_(s)
    {}

    MeshGridIterator(Rect2D const & r, Size2D const & stride = Size2D(1,1))
    : point_(r.left()-1, r.top()),
      stride_(stride),
      rect_(r)
    {}

    MeshGridIterator __iter__() const
    {
        return *this;
    }

    Point2D next()
    {
         point_.x += stride_.x;
         if(point_.x >= rect_.right())
         {
            point_.x = rect_.left();
            point_.y += stride_.y;
         }
         if(point_.y >= rect_.bottom())
         {
            PyErr_SetString(PyExc_StopIteration, "MeshGridIterator.next(): iterator exhausted.");
            python::throw_error_already_set();
         }
         return point_;
    }
};

MeshGridIterator
meshIter_Size2D(Size2D const & self)
{
    return MeshGridIterator(self);
}

MeshGridIterator
meshIter_Rect2D(Rect2D const & self)
{
    return MeshGridIterator(self);
}

Point2D Rect2D__getitem2__(Rect2D const & self, Diff2D const & i)
{
    return self.upperLeft() + i;
}

} // namespace vigra

using namespace vigra;
using namespace python;

void defDiff2D()
{
    class_<MeshGridIterator>("MeshGridIterator", python::no_init)
        .def("__iter__", &MeshGridIterator::__iter__)
        .def("next", &MeshGridIterator::next)
        ;

    class_<Diff2D> diff2d("Diff2D", python::no_init);
    diff2d
        .def_readwrite("x", &Diff2D::x)
        .def_readwrite("y", &Diff2D::y)
        .def("__str__", &diff2StrDiff2D)
        .def("__repr__", &diff2StrDiff2D)
        .def("__getitem__", &Diff2D__getitem__)
        .def("__setitem__", &Diff2D__setitem__)
        .def("__len__", &Diff2D__len__)
        .def("norm", &Diff2D::magnitude)
        .def(self == self)
        .def(self != self)
        .def(self += Diff2D())
        .def(self -= Diff2D())
        .def(self *= double())
        .def(self /= double())
        .def_pickle(Diff2DPickleSuite<Diff2D>())
        ;
    diff2d.attr("magnitude") = diff2d.attr("norm");

    class_<Point2D, bases<Diff2D> >("Point2D")
        .def(init<int, int>())
        .def("__str__", &diff2StrPoint2D)
        .def("__repr__", &diff2StrPoint2D)
        .def(-self)
        .def(self + Diff2D())
        .def(self - Diff2D())
        .def(self - self)
        .def(self * double())
        .def(double() * self)
        .def(self / double())
        .def_pickle(Diff2DPickleSuite<Point2D>())
        ;

    class_<Size2D, bases<Diff2D> >("Size2D")
        .def(init<int, int>())
        .def_readwrite("width", &Diff2D::x)
        .def_readwrite("height", &Diff2D::y)
        .def("area", &Size2D::area)
        .def("__str__", &diff2StrSize2D)
        .def("__repr__", &diff2StrSize2D)
        .def(-self)
        .def(self + Diff2D())
        .def(self + Point2D())
        .def(self - Diff2D())
        .def(self * double())
        .def(double() * self)
        .def(self / double())
        .def_pickle(Diff2DPickleSuite<Size2D>())
        ;

    class_<Rect2D>("Rect2D")
        .def(init<int, int, int, int>())
        .def(init<Point2D const &, Point2D const &>())
        .def(init<Point2D const &, Size2D const &>())
        .def(init<Size2D const &>())
        .def("__str__", &rect2D__str__)
        .def("__repr__", &rect2D__repr__)
        .def("upperLeft", &Rect2D::upperLeft, return_value_policy<copy_const_reference>())
        .def("lowerRight", &Rect2D::lowerRight, return_value_policy<copy_const_reference>())
        .def("setUpperLeft", &Rect2D::setUpperLeft)
        .def("setLowerRight", &Rect2D::setLowerRight)
        .def("moveTo", (void (Rect2D::*)(Point2D const &))&Rect2D::moveTo)
        .def("moveTo", (void (Rect2D::*)(int, int))&Rect2D::moveTo)
        .def("moveBy", (void (Rect2D::*)(Diff2D const &))&Rect2D::moveBy)
        .def("moveBy", (void (Rect2D::*)(int, int))&Rect2D::moveBy)
        .def("left", &Rect2D::left)
        .def("top", &Rect2D::top)
        .def("right", &Rect2D::right)
        .def("bottom", &Rect2D::bottom)
        .def("width", &Rect2D::width)
        .def("height", &Rect2D::height)
        .def("area", &Rect2D::area)
        .def("size", &Rect2D::size)
        .def("setSize", (void (Rect2D::*)(Size2D const &))&Rect2D::setSize)
        .def("setSize", (void (Rect2D::*)(int, int))&Rect2D::setSize)
        .def("addSize", &Rect2D::addSize)
        .def("addBorder", (void (Rect2D::*)(int, int))&Rect2D::addBorder)
        .def("isEmpty", &Rect2D::isEmpty)
        .def("contains", (bool (Rect2D::*)(Point2D const &) const)&Rect2D::contains)
        .def("__contains__", (bool (Rect2D::*)(Point2D const &) const)&Rect2D::contains)
        .def("contains", (bool (Rect2D::*)(Rect2D const &) const)&Rect2D::contains)
        .def("intersects", &Rect2D::intersects)
        .def("__getitem__", &Rect2D__getitem2__)
        .def(self == self)
        .def(self != self)
        .def(self |= Point2D())
        .def(self |= self)
        .def(self &= Point2D())
        .def(self &= self)
        .def(self | Point2D())
        .def(self | self)
        .def(self & Point2D())
        .def(self & self)
        .def(self *= double())
        .def(self *= int())
        .def(self * double())
        .def(self * int())
        .def_pickle(Rect2DPickleSuite())
        ;

    ShapeFromPython<Diff2D>(); // needed for "img[x,y]"
    ShapeFromPython<Size2D>(); // needed for GrayImage((10,10))
    ShapeFromPython<Point2D>(); // for e.g. customKernel(img, (5, 5))
#ifdef VIGRA_PYTHON_HAS_SLICE
    SlicePairFromPython();
#endif

    def("meshIter", &meshIter_Size2D);
    def("meshIter", &meshIter_Rect2D);
}
