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

#ifndef CELLIMAGE_MODULE_HXX
#define CELLIMAGE_MODULE_HXX

#include "foureightsegmentation.hxx"
#include <boost/python.hpp>

template<class T>
T &returnSelf(T &v)
{
    return v;
}

template<class Iterator>
typename Iterator::reference nextIterPos(Iterator &v)
{
    if(!v.inRange())
    {
        PyErr_SetString(PyExc_StopIteration, "cells iterator exhausted");
        boost::python::throw_error_already_set();
    }
    return *v++;
}

struct NodeListProxy
{
    NodeListProxy(vigra::cellimage::GeoMap *segmentation)
    : segmentation_(segmentation)
    {}

    static NodeListProxy create(vigra::cellimage::GeoMap *segmentation)
    {
        return NodeListProxy(segmentation);
    }

    vigra::cellimage::GeoMap *segmentation_;

    long __len__() const
    {
        return segmentation_->nodeCount();
    }

    vigra::cellimage::GeoMap::NodeInfo &__getitem__(long index)
    {
        if(index > (long)segmentation_->maxNodeLabel())
        {
            PyErr_SetObject(
                PyExc_IndexError,
                boost::python::incref(boost::python::object(index).ptr()));
            boost::python::throw_error_already_set();
        }
        return segmentation_->node(index);
    }

    vigra::cellimage::GeoMap::NodeIterator __iter__()
    {
        return segmentation_->nodesBegin();
    }
};

struct EdgeListProxy
{
    EdgeListProxy(vigra::cellimage::GeoMap *segmentation)
    : segmentation_(segmentation)
    {}

    static EdgeListProxy create(vigra::cellimage::GeoMap *segmentation)
    {
        return EdgeListProxy(segmentation);
    }

    vigra::cellimage::GeoMap *segmentation_;

    long __len__() const
    {
        return segmentation_->edgeCount();
    }

    vigra::cellimage::GeoMap::EdgeInfo &__getitem__(long index)
    {
        if(index > (long)segmentation_->maxEdgeLabel())
        {
            PyErr_SetObject(
                PyExc_IndexError,
                boost::python::incref(boost::python::object(index).ptr()));
            boost::python::throw_error_already_set();
        }
        return segmentation_->edge(index);
    }

    vigra::cellimage::GeoMap::EdgeIterator __iter__()
    {
        return segmentation_->edgesBegin();
    }
};

struct FaceListProxy
{
    FaceListProxy(vigra::cellimage::GeoMap *segmentation)
    : segmentation_(segmentation)
    {}

    static FaceListProxy create(vigra::cellimage::GeoMap *segmentation)
    {
        return FaceListProxy(segmentation);
    }

    vigra::cellimage::GeoMap *segmentation_;

    long __len__() const
    {
        return segmentation_->faceCount();
    }

    vigra::cellimage::GeoMap::FaceInfo &__getitem__(long index)
    {
        if(index > (long)segmentation_->maxFaceLabel())
        {
            PyErr_SetObject(
                PyExc_IndexError,
                boost::python::incref(boost::python::object(index).ptr()));
            boost::python::throw_error_already_set();
        }
        return segmentation_->face(index);
    }

    vigra::cellimage::GeoMap::FaceIterator __iter__()
    {
        return segmentation_->facesBegin();
    }
};

#endif // CELLIMAGE_MODULE_HXX
