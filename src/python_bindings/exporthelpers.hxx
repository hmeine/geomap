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


/********************************************************************/

template<class ITERATOR>
class STLIterWrapper
{
  public:
    typedef ITERATOR Iterator;

    STLIterWrapper(Iterator begin, Iterator end)
    : begin_(begin),
      end_(end)
    {}

    STLIterWrapper __iter__() const
    {
        return *this;
    }

    unsigned int __len__() const
    {
        return end_ - begin_;
    }

    typename Iterator::value_type next()
    {
        if(begin_ == end_)
        {
            PyErr_SetString(PyExc_StopIteration, "");
            boost::python::throw_error_already_set();
        }
        return *(begin_++);
    }

  protected:
    Iterator begin_, end_;
};

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
            PyErr_SetString(PyExc_StopIteration, "cells iterator exhausted");
            boost::python::throw_error_already_set();
        }
        return *v++;
    }
};

#endif // EXPORTHELPERS_HXX
