#ifndef EXPORTHELPERS_HXX
#define EXPORTHELPERS_HXX

#include <boost/python.hpp>
#include <boost/python/slice.hpp>
#include <memory>
#include "vigra/pythonutil.hxx"

/*
  I find these very convenient, and would like to add them to
  pythonutil.hxx.  However, everything in there currently is in the
  boost::python namespace, which I dislike, since it's not ours.

  Basically, this is a minimal substitute for the
  vector_indexing_suite.
*/

inline void checkPythonIndex(int &i, unsigned int size)
{
    if(i < 0)
        i += size;
    if((unsigned int)i >= size)
    {
        PyErr_SetString(PyExc_IndexError,
                        "index out of bounds.");
        boost::python::throw_error_already_set();
    }
}

template<class Array>
typename Array::value_type
Array__getitem__(Array const & a, int i)
{
    checkPythonIndex(i, a.size());
    return a[i];
}

template<class Array>
std::auto_ptr<Array>
Array__getitem_slice__(Array const & a, boost::python::slice sl)
{
    boost::python::slice::range<typename Array::const_iterator>
        bounds;
    try
    {
#if BOOST_VERSION >= 104800
        bounds = sl.template get_indices<>(a.begin(), a.end());
#else
        bounds = sl.template get_indicies<>(a.begin(), a.end());
#endif
    }
    catch (std::invalid_argument)
    {
        return std::auto_ptr<Array>(new Array());
    }

    if(bounds.step != 1)
    {
        PyErr_SetString(PyExc_IndexError,
                        "No extended slicing supported yet.");
        boost::python::throw_error_already_set();
    }

    return std::auto_ptr<Array>(new Array(bounds.start, bounds.stop+1));
}

template<class Array>
typename Array::value_type &
Array__getitem__byref(Array & a, int i)
{
    checkPythonIndex(i, a.size());
    return a[i];
}

template<class Array>
void
Array__setitem__(Array & a, int i, typename Array::value_type v)
{
    checkPythonIndex(i, a.size());
    a[i] = v;
}

/********************************************************************/

template<class Iterator,
         class CallPolicies = boost::python::default_call_policies>
struct RangeIterWrapper
: boost::python::class_<Iterator>
{
    RangeIterWrapper(const char *name, CallPolicies cp = CallPolicies())
    : boost::python::class_<Iterator>(name, boost::python::no_init)
    {
        def("__iter__", (Iterator &(*)(Iterator &))&returnSelf,
            boost::python::return_internal_reference<>());
        def("next", &nextIterPos, cp);
    }

    static Iterator &returnSelf(Iterator &v)
    {
        return v;
    }

    static typename Iterator::value_type nextIterPos(Iterator &v)
    {
        if(!v.inRange())
        {
            PyErr_SetString(PyExc_StopIteration, "iterator exhausted");
            boost::python::throw_error_already_set();
        }
        return *v++;
    }
};

/********************************************************************/

template<class Copyable>
boost::python::object
generic__copy__(boost::python::object copyable)
{
    namespace bp = boost::python;

    Copyable *newCopyable = new Copyable(bp::extract<const Copyable &>(copyable)());
    bp::object result =
        bp::object(bp::detail::new_reference(bp::managingPyObject(newCopyable)));

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        copyable.attr("__dict__"));

    return result;
}

template<class Copyable>
boost::python::object
generic__deepcopy__(boost::python::object copyable, boost::python::dict memo)
{
    namespace bp = boost::python;

    bp::object copyMod = bp::import("copy");
    bp::object deepcopy = copyMod.attr("deepcopy");

    Copyable *newCopyable = new Copyable(bp::extract<const Copyable &>(copyable)());
    bp::object result =
        bp::object(bp::detail::new_reference(bp::managingPyObject(newCopyable)));

    // (cf. builtin_id() in Python/bltinmodule.c; copyableId must be
    // the same as the value of id(copyable))
    bp::object copyableId(bp::handle<>(PyLong_FromVoidPtr(copyable.ptr())));
    memo[copyableId] = result;

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        deepcopy(bp::extract<bp::dict>(copyable.attr("__dict__"))(), memo));

    return result;
}

#endif // EXPORTHELPERS_HXX
