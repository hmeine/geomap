#ifndef EXPORTHELPERS_HXX
#define EXPORTHELPERS_HXX

#include <boost/python.hpp>

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

#endif // EXPORTHELPERS_HXX
