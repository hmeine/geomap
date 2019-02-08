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

#ifndef PYTHON_TYPES_HXX_
#define PYTHON_TYPES_HXX_

#include <vigra/numpy_array.hxx>
#include <vigra/tinyvector.hxx>

namespace vigra
{
  typedef float                 GrayValue;
  typedef TinyVector<double, 2> Vector2;
  typedef TinyVector<double, 3> Vector3;

  typedef NumpyArray<2, Singleband<float>, UnstridedArrayTag> UnstridedNumpyFImage;
  typedef NumpyArray<2, TinyVector<float, 2> >                NumpyVector2FImage;
  typedef NumpyArray<2, Singleband<int> >                     NumpyIImage;
}

#endif /* PYTHON_TYPES_HXX_ */
