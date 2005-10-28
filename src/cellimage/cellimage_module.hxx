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
        throw_error_already_set();
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
