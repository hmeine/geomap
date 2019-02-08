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

#include "cellimage_module.hxx"

using namespace vigra;
using namespace boost::python;
using namespace vigra::cellimage;

void defineNodes()
{
    class_<GeoMap::NodeIterator>("NodeIterator", no_init)
        .def("__iter__", (GeoMap::NodeIterator &(*)(GeoMap::NodeIterator &))&returnSelf,
             return_internal_reference<>())
        .def("next", &nextIterPos<GeoMap::NodeIterator>,
             // this is not really true, since the true owner would be
             // the GeoMap object (however, it's
             // lifetime is expected to be long enough):
             return_internal_reference<>());

    class_<NodeListProxy>("NodeList", no_init)
        .def("__len__", &NodeListProxy::__len__)
        .def("__getitem__", &NodeListProxy::__getitem__,
             // this is not really true, see NodeIterator::next
             return_internal_reference<>())
        .def("__iter__", &NodeListProxy::__iter__);
}
