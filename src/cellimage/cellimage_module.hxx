#ifndef CELLIMAGE_MODULE_HXX
#define CELLIMAGE_MODULE_HXX

#include <vigra/vigrapython.hxx>
#include "foureightsegmentation.hxx"

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
    NodeListProxy(vigra::cellimage::FourEightSegmentation *segmentation)
    : segmentation_(segmentation)
    {}
    
    static NodeListProxy create(vigra::cellimage::FourEightSegmentation *segmentation)
    {
        return NodeListProxy(segmentation);
    }

    vigra::cellimage::FourEightSegmentation *segmentation_;

    long __len__() const
    {
        return segmentation_->nodeCount();
    }
    
    vigra::cellimage::FourEightSegmentation::NodeInfo &__getitem__(long index)
    {
        return segmentation_->node(index);
    }

    vigra::cellimage::FourEightSegmentation::NodeIterator __iter__()
    {
        return segmentation_->nodesBegin();
    }
};

struct EdgeListProxy
{
    EdgeListProxy(vigra::cellimage::FourEightSegmentation *segmentation)
    : segmentation_(segmentation)
    {}
    
    static EdgeListProxy create(vigra::cellimage::FourEightSegmentation *segmentation)
    {
        return EdgeListProxy(segmentation);
    }

    vigra::cellimage::FourEightSegmentation *segmentation_;

    long __len__() const
    {
        return segmentation_->edgeCount();
    }
    
    vigra::cellimage::FourEightSegmentation::EdgeInfo &__getitem__(long index)
    {
        return segmentation_->edge(index);
    }

    vigra::cellimage::FourEightSegmentation::EdgeIterator __iter__()
    {
        return segmentation_->edgesBegin();
    }
};

struct FaceListProxy
{
    FaceListProxy(vigra::cellimage::FourEightSegmentation *segmentation)
    : segmentation_(segmentation)
    {}
    
    static FaceListProxy create(vigra::cellimage::FourEightSegmentation *segmentation)
    {
        return FaceListProxy(segmentation);
    }

    vigra::cellimage::FourEightSegmentation *segmentation_;

    long __len__() const
    {
        return segmentation_->faceCount();
    }
    
    vigra::cellimage::FourEightSegmentation::FaceInfo &__getitem__(long index)
    {
        return segmentation_->face(index);
    }

    vigra::cellimage::FourEightSegmentation::FaceIterator __iter__()
    {
        return segmentation_->facesBegin();
    }
};

#endif // CELLIMAGE_MODULE_HXX
