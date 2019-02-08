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

#define PY_ARRAY_UNIQUE_SYMBOL geomap_PyArray_API
#define NO_IMPORT_ARRAY
#include <vigra/numpy_array.hxx>
#include "facestatistics.hxx"
#include "exporthelpers.hxx"
#include <cmath>

namespace bp = boost::python;
using vigra::NumpyFImage;
using vigra::NumpyFRGBImage;
using vigra::Size2D;

template<class TargetType>
struct CastingAccessor
{
    template<class ITERATOR>
    TargetType operator()(ITERATOR const &it) const
    {
        return static_cast<TargetType>(*it);
    }
};

template<class OriginalImage>
class FaceColorStatisticsWrapper
: public bp::class_<FaceColorStatistics<OriginalImage>, boost::noncopyable>
{
  public:
    typedef FaceColorStatistics<OriginalImage> Statistics;
    typedef typename Statistics::Functor StatsFunctor;

    FaceColorStatisticsWrapper(const char *name)
    : bp::class_<Statistics, boost::noncopyable>(name, bp::no_init)
    {
        this->def("__init__", bp::make_constructor(
                      &create,
                      // FaceColorStatistics stores a ref. to originalImage,
                      // actually also to the map, but we need to prevent
                      // cyclic dependencies:
                      bp::default_call_policies(), // FIXME!!
                      // bp::with_custodian_and_ward_postcall<0, 2>(),
                      (bp::arg("map"), bp::arg("originalImage"),
                       bp::arg("minSampleCount") = 1)));

        this->def("__copy__", &generic__copy__<Statistics>);
        this->def("__deepcopy__", &__deepcopy__);

        this->def("pixelCount", &pixelCount);
        this->def("average", &average);
        this->attr("__getitem__") = this->attr("average");
        this->def("variance", &variance);
        this->def("functor", &functor);

        this->def("faceMeanDiff", &Statistics::faceMeanDiff);
        this->def("faceHomogeneity", &Statistics::faceHomogeneity);
        this->def("faceAreaHomogeneity", &Statistics::faceAreaHomogeneity);
#ifdef HAVE_MATH_TOOLKIT
        this->def("faceTTest", &Statistics::faceTTest,
            "Returns the confidence that the distributions of the two\n"
            "adjacent regions do *not* have the same mean, given the\n"
            "samples.  I.e. a value of 0.95 would mean that the confidence\n"
            "of the null hypothesis that the two means are equal is only 5%.");
        this->def("debugTTest", &debugTTest);
#endif

        this->def("regionImage", &regionImage);
        this->def("regionImage", &convertToRegionMeans);

        this->def("attachHooks", &Statistics::attachHooks);
        this->def("detachHooks", &Statistics::detachHooks);
        this->def("map", &Statistics::map);

        this->def("superSampledCount", &Statistics::superSampledCount);
        this->def("minSampleCount", &Statistics::minSampleCount);
        this->def("checkConsistency", &Statistics::checkConsistency);

        bp::scope parent(*this); // Functor shall become a nested class

        bp::class_<StatsFunctor>("Functor")
            .def("pixelCount", &StatsFunctor::count)
            .def("average", &StatsFunctor::average)
            .def("variance", &StatsFunctor::variance)
            .def("__call__",
                 // If you want to see a *really* bad error message,
                 // add a typename before "StatsFunctor const &" in
                 // the next line and compile with GCC 4.1.2... ;-/
                 (void (StatsFunctor::*)(StatsFunctor const &))
                 &StatsFunctor::operator())
            .def("__call__",
                 (void (StatsFunctor::*)(typename StatsFunctor::argument_type const &, double))
                 &StatsFunctor::operator(),
                 bp::args("sample", "weight"))
            .def("__call__",
                 (void (StatsFunctor::*)(typename StatsFunctor::argument_type const &))
                 &StatsFunctor::operator(),
                 bp::args("sample"))
        ;
    }

    // generic__deepcopy__ not applicable - we need to recursively deepcopy the GeoMap
    static boost::python::object
    __deepcopy__(boost::python::object statistics, boost::python::dict memo)
    {
        bp::object copyMod = bp::import("copy");
        bp::object deepcopy = copyMod.attr("deepcopy");

        Statistics *newStatistics = new Statistics(
            bp::extract<const Statistics &>(statistics)());
        bp::object result =
            bp::object(bp::detail::new_reference(bp::managingPyObject(newStatistics)));

        // (cf. builtin_id() in Python/bltinmodule.c; statisticsId must be
        // the same as the value of id(statistics))
        bp::object statisticsId(bp::handle<>(PyLong_FromVoidPtr(statistics.ptr())));
        memo[statisticsId] = result;

        newStatistics->attachHooks(
            bp::extract<boost::shared_ptr<GeoMap> >(
                deepcopy(statistics.attr("map")(), memo)));

        bp::extract<bp::dict>(result.attr("__dict__"))().update(
            deepcopy(bp::extract<bp::dict>(statistics.attr("__dict__"))(), memo));

        return result;
    }

    static inline void
    checkFaceLabel(Statistics const &stats, CellLabel faceLabel)
    {
        if((unsigned int)faceLabel >= stats.size())
        {
            PyErr_SetString(PyExc_IndexError,
                            "face label out of bounds.");
            bp::throw_error_already_set();
        }
        if(!stats[faceLabel])
        {
            PyErr_SetString(PyExc_ValueError,
                            "no information for the given face label.");
            bp::throw_error_already_set();
        }
    }

    static unsigned int
    pixelCount(Statistics const &stats, CellLabel faceLabel)
    {
        checkFaceLabel(stats, faceLabel);
        return stats.pixelCount(faceLabel);
    }

    static typename StatsFunctor::result_type
    average(Statistics const &stats, CellLabel faceLabel)
    {
        checkFaceLabel(stats, faceLabel);
        return stats.average(faceLabel);
    }

    static typename StatsFunctor::result_type
    variance(Statistics const &stats, CellLabel faceLabel, bool unbiased)
    {
        checkFaceLabel(stats, faceLabel);
        return stats.variance(faceLabel, unbiased);
    }

    static bp::object
    functor(Statistics const &stats, CellLabel faceLabel)
    {
        if((unsigned int)faceLabel >= stats.size())
        {
            PyErr_SetString(PyExc_IndexError,
                            "face label out of bounds.");
            bp::throw_error_already_set();
        }
        if(!stats[faceLabel])
            return bp::object();
        return bp::object(*stats[faceLabel]);
    }

#ifdef HAVE_MATH_TOOLKIT
    static bp::dict
    debugTTest(Statistics const &stats, const GeoMap::Dart &dart)
    {
        double
            diff = vigra::norm(average(stats, dart.leftFaceLabel()) -
                               average(stats, dart.rightFaceLabel())),
            sigma1_2 = clampedScalarVariance(
                variance(stats, dart.leftFaceLabel(), true)),
            sigma2_2 = clampedScalarVariance(
                variance(stats, dart.rightFaceLabel(), true)),
            N1 = pixelCount(stats, dart.leftFaceLabel()),
            N2 = pixelCount(stats, dart.rightFaceLabel()),
            dof = N1 + N2 - 2,
            pooledSigma_2 = ((N1-1)*sigma1_2 + (N2-1)*sigma2_2) / dof,
            t = (diff / std::sqrt(pooledSigma_2 * (1.0 / N1 + 1.0 / N2)));

        using namespace boost::math;
        students_t stud(dof);
        bp::dict result;
        result["mu1"] = average(stats, dart.leftFaceLabel());
        result["mu2"] = average(stats, dart.rightFaceLabel());
        result["diff"] = diff;
        result["N1"] = N1;
        result["sigma1_2"] = sigma1_2;
        result["N2"] = N2;
        result["sigma2_2"] = sigma2_2;
        result["pooledSigma_2"] = pooledSigma_2;
        result["dof"] = dof;
        result["t"] = t;
        result["c_tt"] = 2*cdf(stud, t)-1.0;
        result["P(H_0)"] = 2*cdf(complement(stud, t));
        return result;
    }
#endif

    static Statistics *create(
        boost::shared_ptr<GeoMap> map, OriginalImage const &originalImage, int minSampleCount)
    {
        double maxDiffNorm = 255.*std::sqrt((double) OriginalImage::actual_dimension);
        return new Statistics(map, originalImage,
                              maxDiffNorm, minSampleCount);
    }

    static OriginalImage regionImage(const Statistics &stats)
    {
        vigra::TinyVector<long int,2> sizeVector(stats.map()->imageSize());
        OriginalImage result(sizeVector);

        stats.copyRegionImage(destImage(result));

        return result;
    }

    static OriginalImage convertToRegionMeans(
        const Statistics &stats, NumpyFImage labels)
    {

        vigra::TinyVector<long int,2> sizeVector(Size2D(static_cast<int>(labels.shape(0)),
            static_cast<int>(labels.shape(1))));
        OriginalImage result(sizeVector);

        stats.transformRegionImage(
            srcImageRange(labels, CastingAccessor<int>()), destImage(result));

        return result;
    }
};

void defMapStats()
{
    FaceColorStatisticsWrapper<NumpyFImage>("FaceGrayStatistics");
    FaceColorStatisticsWrapper<NumpyFRGBImage>("FaceRGBStatistics");
}
