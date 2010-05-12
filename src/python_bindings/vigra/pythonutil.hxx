/************************************************************************/
/*                                                                      */
/*               Copyright 2002-2003 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_PYTHONUTIL_HXX
#define VIGRA_PYTHONUTIL_HXX

#include "boost/version.hpp"
#include "boost/python.hpp"
#include "boost/python/raw_function.hpp"
#include "boost/python/detail/api_placeholder.hpp" // python::len

//#include "boost/mpl/vector.hpp"

#if defined(BOOST_VERSION) && BOOST_VERSION >= 103200
#   define VIGRA_PYTHON_HAS_SLICE
#endif

namespace boost { namespace python {

inline void throw_type_error(char const * message)
{
    PyErr_SetString(PyExc_TypeError, message);
    python::throw_error_already_set();
}

namespace detail {

template <int>
struct ConvertSequenceOfVectors;

template <>
struct ConvertSequenceOfVectors<2>
{
    template <class Sequence>
    static list createImpl(Sequence const & src)
    {
        list dest;
        typedef typename Sequence::const_iterator Iter;
        for(Iter i = src.begin(); i != src.end(); ++i)
        {
            dest.append(make_tuple((*i)[0], (*i)[1]));
        }
        return dest;
    }

    template <class BackInsertable>
    static void extractImpl(list src, BackInsertable & dest)
    {
        typedef typename BackInsertable::value_type Vector;
        typedef typename Vector::value_type Element;

        for(int k = 0; k < len(src); ++k)
        {
            extract<Vector> ev(src[k]);
            extract<tuple> et(src[k]);
            extract<list> el(src[k]);
            if(ev.check())
            {
                dest.push_back(ev());
            }
            else if(et.check())
            {
                tuple t = et();
                if(len(t) != 2)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector2");
                dest.push_back(Vector(extract<Element>(t[0])(), extract<Element>(t[1])()));
            }
            else if(el.check())
            {
                list l = el();
                if(len(l) != 2)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector2");
                dest.push_back(Vector(extract<Element>(l[0])(), extract<Element>(l[1])()));
            }
            else
            {
                throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector2");
            }
        }
    }
};

template <>
struct ConvertSequenceOfVectors<3>
{
    template <class Sequence>
    static list createImpl(Sequence const & src)
    {
        list dest;
        typedef typename Sequence::const_iterator Iter;
        for(Iter i = src.begin(); i != src.end(); ++i)
        {
            dest.append(make_tuple((*i)[0], (*i)[1], (*i)[2]));
        }
        return dest;
    }

    template <class BackInsertable>
    static void extractImpl(list src, BackInsertable & dest)
    {
        typedef typename BackInsertable::value_type Vector;
        typedef typename Vector::value_type Element;

        for(int k = 0; k < len(src); ++k)
        {
            extract<Vector> ev(src[k]);
            if(ev.check())
            {
                dest.push_back(ev());
                continue;
            }

            extract<tuple> et(src[k]);
            if(et.check())
            {
                tuple t = et();
                if(len(t) != 3)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector3");
                dest.push_back(Vector(extract<Element>(t[0])(), extract<Element>(t[1])(), extract<Element>(t[2])()));
                continue;
            }

            extract<list> el(src[k]);
            if(el.check())
            {
                list l = el();
                if(len(l) != 3)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector3");
                dest.push_back(Vector(extract<Element>(l[0])(), extract<Element>(l[1])()), extract<Element>(l[2])());
                continue;
            }

            throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector3");
        }
    }
};

template <>
struct ConvertSequenceOfVectors<4>
{
    template <class Sequence>
    static list createImpl(Sequence const & src)
    {
        list dest;
        typedef typename Sequence::const_iterator Iter;
        for(Iter i = src.begin(); i != src.end(); ++i)
        {
            dest.append(make_tuple((*i)[0], (*i)[1], (*i)[2], (*i)[3]));
        }
        return dest;
    }

    template <class BackInsertable>
    static void extractImpl(list src, BackInsertable & dest)
    {
        typedef typename BackInsertable::value_type Vector;
        typedef typename Vector::value_type Element;

        for(int k = 0; k < len(src); ++k)
        {
            extract<Vector> ev(src[k]);
            extract<tuple> et(src[k]);
            extract<list> el(src[k]);
            if(ev.check())
            {
                dest.push_back(ev());
            }
            else if(et.check())
            {
                tuple t = et();
                if(len(t) != 4)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector4");
                dest.push_back(Vector(extract<Element>(t[0])(), extract<Element>(t[1])(),
                                      extract<Element>(t[2])(), extract<Element>(t[3])()));
            }
            else if(el.check())
            {
                list l = el();
                if(len(l) != 4)
                    throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector4");
                dest.push_back(Vector(extract<Element>(l[0])(), extract<Element>(l[1])(),
                                      extract<Element>(l[2])(), extract<Element>(l[3])()));
            }
            else
            {
                throw_type_error("extractSequenceOfVectors(): source list elements are not compatible to Vector4");
            }
        }
    }
};

} // namespace detail

template <class BackInsertable>
void extractSequenceOfVectors(list src, BackInsertable & dest)
{
    typedef typename BackInsertable::value_type Vector;
    detail::ConvertSequenceOfVectors<Vector::static_size>::extractImpl(src, dest);
}

template <class Sequence>
list createSequenceOfVectors(Sequence const & src)
{
    list dest;
    typedef typename Sequence::const_iterator Iter;
    for(Iter i = src.begin(); i != src.end(); ++i)
    {
        dest.append(*i);
    }
    return dest;
}

template <class Sequence>
list createSequenceOfTuples(Sequence const & src)
{
    typedef typename Sequence::value_type Vector;
    return detail::ConvertSequenceOfVectors<Vector::static_size>::createImpl(src);
}

template<class T>
inline PyObject * managingPyObject(T *p)
{
    return typename manage_new_object::apply<T *>::type()(p);
}

template<class T>
inline PyObject * managingPyObject(std::auto_ptr<T> p)
{
    return managingPyObject(p.release());
}

namespace detail
{

  template <class F>
  struct raw_constructor_dispatcher
  {
      raw_constructor_dispatcher(F f)
     : f(make_constructor(f)) {}

      PyObject* operator()(PyObject* args, PyObject* keywords)
      {
          borrowed_reference_t* ra = borrowed_reference(args);
          object a(ra);
          return incref(
              object(
                  f(
                      object(a[0])
                    , object(a.slice(1, len(a)))
                    , keywords ? dict(borrowed_reference(keywords)) : dict()
                  )
              ).ptr()
          );
      }

   private:
      object f;
  };

} // namespace detail

template <class F>
object raw_constructor(F f, std::size_t min_args = 0)
{
    return detail::make_raw_function(
        objects::py_function(
            detail::raw_constructor_dispatcher<F>(f)
          , mpl::vector2<void, object>()
          , min_args+1
          , (std::numeric_limits<unsigned>::max)()
        )
    );
}

#if BOOST_VERSION < 103400

inline object import(str name)
{
  // should be 'char const *' but older python versions don't use 'const' yet.
  char *n = extract<char *>(name);
  handle<> module(PyImport_ImportModule(n));
  return object(module);
}

#endif

}} // namespace boost::python

#endif /* VIGRA_PYTHONUTIL_HXX */
