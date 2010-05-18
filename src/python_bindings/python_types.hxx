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
