#include "facestatistics.hxx"
#include "exporthelpers.hxx"
#include <cmath>

namespace bp = boost::python;

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
: bp::class_<FaceColorStatistics<OriginalImage>, boost::noncopyable>
{
  public:
    typedef FaceColorStatistics<OriginalImage> Statistics;
    typedef typename Statistics::Functor StatsFunctor;

    FaceColorStatisticsWrapper(const char *name)
    : bp::class_<Statistics, boost::noncopyable>(name, bp::no_init)
    {
        def("__init__", bp::make_constructor(
                &create,
                // FaceColorStatistics stores a ref. to originalImage,
                // actually also to the map, but we need to prevent
                // cyclic dependencies:
                bp::default_call_policies(), // FIXME!!
                // bp::with_custodian_and_ward_postcall<0, 2>(),
                (bp::arg("map"), bp::arg("originalImage"),
                 bp::arg("minSampleCount") = 1)));

        def("pixelCount", &pixelCount);
        def("average", &average);
        this->attr("__getitem__") = this->attr("average");
        def("variance", &variance);
        def("functor", &functor);

        def("faceMeanDiff", &Statistics::faceMeanDiff);
        def("faceHomogeneity", &Statistics::faceHomogeneity);
        def("faceAreaHomogeneity", &Statistics::faceAreaHomogeneity);
#ifdef HAVE_MATH_TOOLKIT
        def("faceTTest", &Statistics::faceTTest,
            "Returns the confidence that the distributions of the two\n"
            "adjacent regions do *not* have the same mean, given the\n"
            "samples.  I.e. a value of 0.95 would mean that the confidence\n"
            "of the null hypothesis that the two means are equal is only 5%.");
        def("debugTTest", &debugTTest);
#endif

        def("regionImage", &regionImage);
        def("regionImage", &convertToRegionMeans);

        def("detachHooks", &Statistics::detachHooks);

        def("superSampledCount", &Statistics::superSampledCount);
        def("minSampleCount", &Statistics::minSampleCount);
        def("checkConsistency", &Statistics::checkConsistency);

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
                 (void (StatsFunctor::*)(typename StatsFunctor::argument_type const &))
                 &StatsFunctor::operator())
        ;
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
        GeoMap &map, OriginalImage const &originalImage, int minSampleCount)
    {
        double maxDiffNorm = 255.*std::sqrt((double)originalImage.bands());
        return new Statistics(map, originalImage,
                              maxDiffNorm, minSampleCount);
    }

    static OriginalImage regionImage(const Statistics &stats)
    {
        OriginalImage result(stats.map()->imageSize());

        stats.copyRegionImage(destImage(result));

        return result;
    }

    static OriginalImage convertToRegionMeans(
        const Statistics &stats, vigra::PythonSingleBandImage labels)
    {
        OriginalImage result(labels.size());

        stats.transformRegionImage(
            srcImageRange(labels, CastingAccessor<int>()), destImage(result));

        return result;
    }
};

void defMapStats()
{
    FaceColorStatisticsWrapper<vigra::PythonGrayImage>("FaceGrayStatistics");
    FaceColorStatisticsWrapper<vigra::PythonVector3Image>("FaceRGBStatistics");
}
