#ifndef VIGRA_GRADIENT_HXX
#define VIGRA_GRADIENT_HXX

#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/numerictraits.hxx>
#include <vigra/convolution.hxx>

namespace vigra {

template <class SIter, class SAcc,
          class DIter, class DAcc>
void gradient(SIter sul, SIter slr, SAcc src,
             DIter gul1, DAcc g1, DIter gul2, DAcc g2)
{
    Kernel1D<double> symgrad;
    symgrad.initSymmetricGradient();

    separableConvolveX(srcIterRange(sul, slr, src), destIter(gul1, g1), kernel1d(symgrad));
    separableConvolveY(srcIterRange(sul, slr, src), destIter(gul2, g2), kernel1d(symgrad));
}

template <class SIter, class SAcc,
          class DIter, class DAcc>
inline
void gradient(triple<SIter, SIter, SAcc> src,
              pair<DIter, DAcc> g1, pair<DIter, DAcc> g2)
{
    gradient(src.first, src.second, src.third,
             g1.first, g1.second, g2.first, g2.second);
}

struct GradientMagnitudeFunctor
{
    template <class T>
    T operator()(T const & gx, T const & gy) const
    {
        return sqrt(gx*gx + gy*gy);
    }

    template <class T>
    T operator()(RGBValue<T> const & gx, RGBValue<T> const & gy) const
    {
        return sqrt(gx.red()*gx.red() + gx.green()*gx.green() +
                    gx.blue()*gx.blue() + gy.red()*gy.red() +
                    gy.green()*gy.green() + gy.blue()*gy.blue());
    }
};

template <class SIter, class SAcc,
          class DIter, class DAcc>
inline
void gradientMagnitude(SIter gul1, SIter glr1, SAcc g1, SIter gul2, SAcc g2,
                       DIter dul, DAcc dest)
{
    combineTwoImages(gul1, glr1, g1, gul2, g2, dul, dest, GradientMagnitudeFunctor());
}

template <class SIter, class SAcc,
          class DIter, class DAcc>
inline
void gradientMagnitude(triple<SIter, SIter, SAcc> g1,
                       pair<SIter, SAcc> g2,
                       pair<DIter, DAcc> dest)
{
    combineTwoImages(g1, g2, dest, GradientMagnitudeFunctor());
}

template <class SIter, class SAcc, class DIter, class DAcc>
void gradientMagnitude(SIter sul, SIter slr, SAcc src, DIter dul, DAcc dest)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    typedef typename NumericTraits<typename SAcc::value_type>::Promote
        TmpType;

    BasicImage<TmpType> dx(w,h), dy(w,h);

    gradient(srcIterRange(sul, slr, src), destImage(dx), destImage(dy));

    gradientMagnitude(srcImageRange(dx), srcImage(dy), destIter(dul, dest));
}

template <class SIter, class SAcc, class DIter, class DAcc>
inline
void gradientMagnitude(triple<SIter, SIter, SAcc> src,
                       pair<DIter, DAcc> dest)
{
    gradientMagnitude(src.first, src.second, src.third,
                         dest.first, dest.second);
}

} // namespace vigra

#endif /* VIGRA_GRADIENT_HXX */
