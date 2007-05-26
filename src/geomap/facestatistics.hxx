#ifndef VIGRA_FACESTATISTICS_HXX
#define VIGRA_FACESTATISTICS_HXX

#include "cppmap.hxx"
#include <vigra/inspectimage.hxx>
#include <vigra/transformimage.hxx>

namespace detail {

template<class FACE_STATS>
class InitFaceFunctors
{
    std::vector<FACE_STATS *> &stats_;

  public:
    InitFaceFunctors(std::vector<FACE_STATS *> &stats)
    : stats_(stats)
    {}

    void operator()(int label, const typename FACE_STATS::argument_type &v) const
    {
        if(label >= 0)
            (*stats_[label])(v);
    }
};

template<class FACE_STATS>
class LookupFaceAverage
{
    const std::vector<FACE_STATS *> &stats_;

  public:
    LookupFaceAverage(const std::vector<FACE_STATS *> &stats)
    : stats_(stats)
    {}

    typename FACE_STATS::result_type
    operator()(int label) const
    {
        if(label >= 0)
            return stats_[label]->average();
        return typename FACE_STATS::result_type();
    }
};

} // namespace detail

/********************************************************************/

template<class OriginalImage>
class FaceColorStatistics : boost::noncopyable
{
  public:
    typedef vigra::FindAverageAndVariance<typename OriginalImage::value_type>
        Functor;

//     template<int SplineOrder>
    FaceColorStatistics(GeoMap &map, const OriginalImage &originalImage,
                        double maxDiffNorm = 1.0, int minSampleCount = 1);
    ~FaceColorStatistics();

    unsigned int
    pixelCount(CellLabel faceLabel) const
    {
        return functors_[faceLabel]->count();
    }

    typename Functor::result_type
    average(CellLabel faceLabel) const
    {
        return functors_[faceLabel]->average();
    }

    typename Functor::result_type
    variance(CellLabel faceLabel, bool unbiased) const
    {
        return functors_[faceLabel]->variance(unbiased);
    }

    double faceMeanDiff(const GeoMap::Dart &dart)
    {
        return vigra::norm(average(dart.leftFaceLabel()) -
                           average(dart.rightFaceLabel()))
            / maxDiffNorm_;
    }

    template<class DEST_ITERATOR, class DEST_ACCESSOR>
    void copyRegionImage(DEST_ITERATOR dul, DEST_ACCESSOR da) const;

    template<class DEST_ITERATOR, class DEST_ACCESSOR>
    void copyRegionImage(std::pair<DEST_ITERATOR, DEST_ACCESSOR> d) const
    {
        copyRegionImage(d.first, d.second);
    }

    const GeoMap *map() const
    {
        return &map_;
    }

  protected:
    GeoMap &map_;
    const OriginalImage &originalImage_;
    std::vector<Functor *> functors_;
    double maxDiffNorm_;
};

template<class OriginalImage>
// template<int SplineOrder>
FaceColorStatistics<OriginalImage>::FaceColorStatistics(
    GeoMap &map, const OriginalImage &originalImage,
    double maxDiffNorm, int minSampleCount)
: map_(map),
  originalImage_(originalImage),
  functors_(map.maxFaceLabel(), NULL),
  maxDiffNorm_(maxDiffNorm)
{
    for(GeoMap::FaceIterator it = map.facesBegin(); it.inRange(); ++it)
        functors_[(*it)->label()] = new Functor();

    detail::InitFaceFunctors<Functor> iff(functors_);
    inspectTwoImages(map.srcLabelRange(),
                     srcImage(originalImage),
                     iff);

    // FIXME: ensure minSampleCount
}

template<class OriginalImage>
FaceColorStatistics<OriginalImage>::~FaceColorStatistics()
{
    for(unsigned int i = 0; i < functors_.size(); ++i)
        delete functors_[i];
}

template<class OriginalImage>
template<class DEST_ITERATOR, class DEST_ACCESSOR>
void FaceColorStatistics<OriginalImage>::copyRegionImage(
    DEST_ITERATOR dul, DEST_ACCESSOR da) const
{
    transformImage(map_.srcLabelRange(), destIter(dul, da),
                   detail::LookupFaceAverage<Functor>(functors_));
}

#endif // VIGRA_FACESTATISTICS_HXX
