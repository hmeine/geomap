#include "python_types.hxx"

#include <vigra/numerictraits.hxx>
#include <vigra/tinyvector.hxx>

#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

#include <sstream>
#include <iomanip>

namespace python = boost::python;

namespace vigra
{

template<class Vector>
PyObject * Vector__str__(Vector const & v)
{
        std::stringstream s;
        if(v.size() == 1)
        {
            s << v[0];
        }
        else if(v.size() > 0)
        {
            s << "(";
            unsigned int i = 0;
            for(; i<v.size()-1; ++i)
                s << v[i] << ", ";
            s << v[i] << ")";
        }
        else
        {
            // should not happen, can only be created from C++ anyways:
            s << "<empty PythonVector>";
        }
        return PyString_FromString(s.str().c_str());
}

template <class Vector>
PyObject * Vector__repr__(Vector const & v)
{
        std::stringstream s;
        s << std::setprecision(16)
          << "Vector" << v.size() << "(";
        unsigned char i = 0;
        for(; i<v.size()-1; ++i)
            s << v[i] << ", ";
        s << v[i] << ")";
        return PyString_FromString(s.str().c_str());
}

template<class Vector>
double Vector__getitem__(Vector const & p, unsigned int i)
{
    if(i >= p.size())
    {
        PyErr_SetString(PyExc_IndexError,
            "PythonVector.__getitem__(): index out of bounds.");
        python::throw_error_already_set();
    }
    return p[i];
}

template<class Vector>
void Vector__setitem__(Vector & p, unsigned int i, double v)
{
    if(i >= p.size())
    {
        PyErr_SetString(PyExc_IndexError,
            "PythonVector.__setitem__(): index out of bounds.");
        python::throw_error_already_set();
    }
    p[i] = v;
}

template<class Vector>
typename vigra::NumericTraits<typename Vector::value_type>::Promote
pyDot(Vector const & l, Vector const & r)
{
    return dot(l, r);
}

template<class Vector>
typename NormTraits<Vector>::NormType
pyNorm(Vector const & v)
{
    return norm(v);
}

using python::self;

template <class Vector>
struct VectorFromPython
{
    typedef Vector Type;

    VectorFromPython()
    {
        python::converter::registry::insert(
            &convertible, &construct, python::type_id<Type>());
    }

    static void* convertible(PyObject* obj)
    {
        PyObject *seq = PySequence_Fast(obj, "vector-from-python conversion");
        bool possible = (seq != NULL);

        if(possible)
        {
            possible = (PySequence_Size(seq) == Type::static_size);
            for(int i = 0; possible && (i < Type::static_size); ++i)
                possible = python::extract<double>(
                    PySequence_Fast_GET_ITEM(seq, i)).check();
            Py_XDECREF(seq);
        }
        else
            PyErr_Clear();

        return possible ? obj : NULL;
    }

    static void construct(PyObject* obj, python::converter::rvalue_from_python_stage1_data* data)
    {
        Type * const storage = reinterpret_cast<Type *>((
            reinterpret_cast<python::converter::rvalue_from_python_storage<Type> *>(
                data))->storage.bytes);

        PyObject *seq = PySequence_Fast(obj, "vector-from-python conversion");
        for(int i = 0; i < Type::static_size; ++i)
        {
            (*storage)[i] = python::extract<double>(
                PySequence_Fast_GET_ITEM(seq, i))();
        }
        Py_XDECREF(seq);

        data->convertible = storage;
    }
};

template <class Vector>
struct VectorPickleSuite : python::pickle_suite
{
    static python::tuple getinitargs(Vector const& v)
    {
        switch(Vector::static_size)
        {
          case 2:
              return python::make_tuple(v[0], v[1]);
          case 3:
              return python::make_tuple(v[0], v[1], v[2]);
          case 4:
              return python::make_tuple(v[0], v[1], v[2], v[3]);
          default:
              ;
        }
        PyErr_SetString(
            PyExc_TypeError,
            "VectorPickleSuite can only pickle Vector2/3/4.");
        python::throw_error_already_set();
        return python::make_tuple(); // never reached
    }
};

template <class Vector>
void defineVectorFunctions(python::class_<Vector> & cl)
{
    typedef typename Vector::value_type Value;
      cl.def(self == self)
        .def(self != self)
        .def(self + self)
        .def(self - self)
        .def(-self)
        .def(self * self)
        .def(self += self)
        .def(self -= self)
        .def(self *= self)
        .def(Value() * self)
        .def(self * Value())
        .def(self / Value())
        .def(self *= Value())
        .def(self /= Value())
        .def("__len__", &Vector::size)
        .def("size", &Vector::size)
        .def("norm", &Vector::magnitude)
        .def("squaredMagnitude", &Vector::squaredMagnitude)
        .def("__getitem__", &Vector__getitem__<Vector>)
        .def("__setitem__", &Vector__setitem__<Vector>)
        .def("__str__", &Vector__str__<Vector>)
        .def("__repr__", &Vector__repr__<Vector>)
        .def_pickle(VectorPickleSuite<Vector>());
    ;
    cl.attr("magnitude") = cl.attr("norm");
    python::def("norm", &pyNorm<Vector>);
    python::def("dot", &pyDot<Vector>);
}

} // namespace vigra

void defVectorConverters()
{
    using namespace vigra;
    using namespace boost::python;

    VectorFromPython<Vector2>();
    VectorFromPython<Vector3>();

    def("norm", &pyNorm<double>);
    scope().attr("magnitude") = scope().attr("norm"); // install alias

    defineVectorFunctions(python::class_<Vector2>("Vector2")
                          .def(python::init<double, double>()));
    defineVectorFunctions(python::class_<Vector3>("Vector3")
                          .def(python::init<double, double, double>()));
}
