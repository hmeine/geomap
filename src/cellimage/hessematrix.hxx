#ifndef VIGRA_HESSEMATRIX_HXX
#define VIGRA_HESSEMATRIX_HXX

#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/convolution.hxx>
#include <vigra/functorexpression.hxx>

namespace vigra {

template <class SIter, class SAcc, class DIter, class DAcc>
void determinantOfHessian(SIter sul, SIter slr, SAcc src,
                          DIter dul, DAcc dest)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    Kernel1D<double> symgrad;
    symgrad.initSymmetricGradient();

    typedef typename NumericTraits<typename SAcc::value_type>::Promote
        TmpType;

    BasicImage<TmpType> dxx(w,h), dyy(w,h), dxy(w,h), tmp(w,h);

    separableConvolveX(srcIterRange(sul, slr, src), destImage(tmp), kernel1d(symgrad));
    separableConvolveX(srcImageRange(tmp), destImage(dxx), kernel1d(symgrad));

    separableConvolveY(srcImageRange(tmp), destImage(dxy), kernel1d(symgrad));

    separableConvolveY(srcIterRange(sul, slr, src), destImage(tmp), kernel1d(symgrad));
    separableConvolveY(srcImageRange(tmp), destImage(dyy), kernel1d(symgrad));

    combineThreeImages(srcImageRange(dxx), srcImage(dyy), srcImage(dxy),
                       destIter(dul, dest), Arg1()*Arg2() - Arg3()*Arg3());
}

template <class SIter, class SAcc, class DIter, class DAcc>
inline
void determinantOfHessian(triple<SIter, SIter, SAcc> src,
                          pair<DIter, DAcc> dest)
{
    determinantOfHessian(src.first, src.second, src.third,
                         dest.first, dest.second);
}

template <class SIter, class SAcc, class DIter, class DAcc>
void eigenvaluesOfHessian(SIter sul, SIter slr, SAcc src,
                          DIter eul1, DAcc e1, DIter eul2, DAcc e2)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    Kernel1D<double> symgrad;
    symgrad.initSymmetricGradient();

    typedef typename NumericTraits<typename SAcc::value_type>::Promote
        TmpType;

    BasicImage<TmpType> dxx(w,h), dyy(w,h), dxy(w,h), tmp(w,h);

    separableConvolveX(srcIterRange(sul, slr, src), destImage(tmp), kernel1d(symgrad));
    separableConvolveX(srcImageRange(tmp), destImage(dxx), kernel1d(symgrad));

    separableConvolveY(srcImageRange(tmp), destImage(dxy), kernel1d(symgrad));

    separableConvolveY(srcIterRange(sul, slr, src), destImage(tmp), kernel1d(symgrad));
    separableConvolveY(srcImageRange(tmp), destImage(dyy), kernel1d(symgrad));

    int x,y;

    for(y=0; y<h; ++y, ++eul1.y, ++eul2.y)
    {
        DIter ex1 = eul1;
        DIter ex2 = eul2;

        for(x=0; x<w; ++x, ++ex1.x, ++ex2.x)
        {
            TmpType a = (dxx(x,y) + dyy(x,y)) / 2.0;
            TmpType b = (dxx(x,y) - dyy(x,y)) / 2.0;
            TmpType c = sqrt(b*b + dxy(x,y) * dxy(x,y));

            e1.set(a + c, ex1);
            e2.set(a - c, ex2);
        }
    }
}

template <class SIter, class SAcc, class DIter, class DAcc>
inline
void eigenvaluesOfHessian(triple<SIter, SIter, SAcc> src,
                          pair<DIter, DAcc> e1, pair<DIter, DAcc> e2)
{
    eigenvaluesOfHessian(src.first, src.second, src.third,
                         e1.first, e1.second, e2.first, e2.second);
}

} // namespace vigra

#endif /* VIGRA_HESSEMATRIX_HXX */
