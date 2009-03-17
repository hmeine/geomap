#ifndef VIGRA_FACESTATISTICS_HXX
#define VIGRA_FACESTATISTICS_HXX

#include "cppmap.hxx"
#include <vigra/inspectimage.hxx>
#include <vigra/transformimage.hxx>
#include <vigra/splineimageview.hxx>
#include <cmath>

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
class FaceColorStatistics 
: boost::noncopyable
{
  public:
    typedef vigra::FindAverageAndVariance<typename OriginalImage::value_type>
        Functor;

//     template<int SplineOrder>
    FaceColorStatistics(GeoMap &map, const OriginalImage &originalImage,
                        double maxDiffNorm = 1.0, unsigned int minSampleCount = 1);
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

    unsigned int size() const
    {
        return functors_.size();
    }

    const Functor *operator[](unsigned int index) const
    {
        return functors_[index];
    }

    bool preMergeFaces(const GeoMap::Dart &dart)
    {
        if(superSampled_.get())
        {
            unsigned char ssLeft((*superSampled_)[dart.leftFaceLabel()]);
            unsigned char ssRight((*superSampled_)[dart.rightFaceLabel()]);

            if(ssLeft != ssRight)
            {
                if(ssLeft < ssRight)
                {
                    merged_ = *functors_[dart.leftFaceLabel()];
                    mergedSS_ = ssLeft;
                    mergeDecreasesSSCount_ = !ssLeft;
                }
                else
                {
                    merged_ = *functors_[dart.rightFaceLabel()];
                    mergedSS_ = ssRight;
                    mergeDecreasesSSCount_ = !ssRight;
                }
                return true;
            }

            mergedSS_ = ssLeft;
            mergeDecreasesSSCount_ = false;
        }

        merged_ = *functors_[dart.leftFaceLabel()];
        merged_(*functors_[dart.rightFaceLabel()]);
        return true;
    }

    void postMergeFaces(GeoMap::Face &face)
    {
        *functors_[face.label()] = merged_;
        if(superSampled_.get())
        {
            (*superSampled_)[face.label()] = mergedSS_;
            if(mergeDecreasesSSCount_)
                if(!--superSampledCount_)
                    superSampled_.reset();
        }
    }

    void associatePixels(GeoMap::Face &face, const PixelList &pixels)
    {
        Functor &f(*functors_[face.label()]);
        for(PixelList::const_iterator it = pixels.begin();
            it != pixels.end(); ++it)
        {
            f(originalImage_[*it]);
        }
    }

    double faceMeanDiff(const GeoMap::Dart &dart)
    {
        return vigra::norm(average(dart.leftFaceLabel()) -
                           average(dart.rightFaceLabel()))
            / maxDiffNorm_;
    }

    double faceAreaHomogenity(const GeoMap::Dart &dart)
    {
        double a1 = dart.leftFace()->area();
        double a2 = dart.rightFace()->area();
        return (a1 * a2) / (a1 + a2);
    }

    double faceHomogenity2(const GeoMap::Dart &dart)
    {
        return vigra::squaredNorm(faceMeanDiff(dart))
            * faceAreaHomogenity(dart);
    }

    double faceHomogenity(const GeoMap::Dart &dart)
    {
        return std::sqrt(faceHomogenity2(dart));
    }

    template<class DEST_ITERATOR, class DEST_ACCESSOR>
    void copyRegionImage(DEST_ITERATOR dul, DEST_ACCESSOR da) const;

    template<class DEST_ITERATOR, class DEST_ACCESSOR>
    void copyRegionImage(std::pair<DEST_ITERATOR, DEST_ACCESSOR> d) const
    {
        copyRegionImage(d.first, d.second);
    }

    template<class SRC_ITERATOR, class SRC_ACCESSOR,
             class DEST_ITERATOR, class DEST_ACCESSOR>
    void transformRegionImage(
        SRC_ITERATOR sul, SRC_ITERATOR slr, SRC_ACCESSOR sa,
        DEST_ITERATOR dul, DEST_ACCESSOR da) const;

    template<class SRC_ITERATOR, class SRC_ACCESSOR,
             class DEST_ITERATOR, class DEST_ACCESSOR>
    void transformRegionImage(
        vigra::triple<SRC_ITERATOR, SRC_ITERATOR, SRC_ACCESSOR> s,
        std::pair<DEST_ITERATOR, DEST_ACCESSOR> d) const
    {
        transformRegionImage(
            s.first, s.second, s.third,
            d.first, d.second);
    }

    const GeoMap *map() const
    {
        return &map_;
    }

  protected:
    void ensureMinSampleCount(unsigned int minSampleCount);

    GeoMap &map_;
    
    const OriginalImage originalImage_;
    std::vector<Functor *> functors_;

    std::auto_ptr<std::vector<unsigned char> > superSampled_;
    unsigned int superSampledCount_;

    Functor merged_;
    unsigned char mergedSS_;
    bool mergeDecreasesSSCount_;

    double maxDiffNorm_;
};

template<class OriginalImage>
// template<int SplineOrder>
FaceColorStatistics<OriginalImage>::FaceColorStatistics(
    GeoMap &map, const OriginalImage &originalImage,
    double maxDiffNorm, unsigned int minSampleCount)
: map_(map),
  originalImage_(originalImage),
  functors_(map.maxFaceLabel(), NULL),
  superSampledCount_(0),
  maxDiffNorm_(maxDiffNorm)
{
    for(GeoMap::FaceIterator it = map.facesBegin(); it.inRange(); ++it)
        functors_[(*it)->label()] = new Functor();

    detail::InitFaceFunctors<Functor> iff(functors_);
    inspectTwoImages(map.srcLabelRange(),
                     srcImage(originalImage),
                     iff);

    if(minSampleCount)
        ensureMinSampleCount(minSampleCount);
}

template<class OriginalImage>
void FaceColorStatistics<OriginalImage>::ensureMinSampleCount(
    unsigned int minSampleCount)
{
    std::auto_ptr<std::vector<unsigned char> > superSampled;
    typedef typename vigra::NumericTraits<
        typename OriginalImage::value_type>::RealPromote FloatPixel;
    std::auto_ptr<vigra::SplineImageView<5, FloatPixel> > siv;

    GeoMap::FaceIterator it = map_.facesBegin();
    for(++it; it.inRange(); ++it)
    {
        if(functors_[(*it)->label()]->count() < minSampleCount)
        {
            if(!superSampled.get())
            {
                superSampled.reset(
                    new std::vector<unsigned char>(map_.maxFaceLabel(), 0));
                siv.reset(
                    new vigra::SplineImageView<5, FloatPixel>(
                        srcImageRange(originalImage_)));
            }

            do
            {
                unsigned char ss(++(*superSampled)[(*it)->label()]);
                double gridDist = 1.0/(1+ss);
                GeoMap::Face::BoundingBox bbox((*it)->boundingBox());
                for(double y = bbox.begin()[1]; y < bbox.end()[1]; y += gridDist)
                    for(double x = bbox.begin()[0]; x < bbox.end()[0]; x += gridDist)
                        if((*it)->contains(Vector2(x, y)))
                            (*functors_[(*it)->label()])((*siv)(x, y));
            }
            while(functors_[(*it)->label()]->count() < minSampleCount);

            ++superSampledCount_;
        }
    }

    if(superSampledCount_)
        superSampled_ = superSampled;
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
    transformRegionImage(
        map_.srcLabelRange(), destIter(dul, da));
}

template<class OriginalImage>
template<class SRC_ITERATOR, class SRC_ACCESSOR,
         class DEST_ITERATOR, class DEST_ACCESSOR>
void FaceColorStatistics<OriginalImage>::transformRegionImage(
    SRC_ITERATOR sul, SRC_ITERATOR slr, SRC_ACCESSOR sa,
    DEST_ITERATOR dul, DEST_ACCESSOR da) const
{
    transformImage(sul, slr, sa, dul, da,
                   detail::LookupFaceAverage<Functor>(functors_));
}

template<class OriginalImage>
class FaceColorStatisticsCallback
: public FaceColorStatistics<OriginalImage>, public vigra::signal::CallbackBase
{
  public:

    FaceColorStatisticsCallback(GeoMap &map, const OriginalImage &originalImage,
                        double maxDiffNorm = 1.0, unsigned int minSampleCount = 1)
    : FaceColorStatistics<OriginalImage>(map, originalImage, maxDiffNorm, minSampleCount)
    {
        map.preMergeFacesHook.connect(this, &FaceColorStatisticsCallback::preMergeFaces);
        map.postMergeFacesHook.connect(this, &FaceColorStatisticsCallback::postMergeFaces);
        map.associatePixelsHook.connect(this, &FaceColorStatisticsCallback::associatePixels);
    }

    void detachHooks()
    {
        this->disconnectAll();
    }
};

#endif // VIGRA_FACESTATISTICS_HXX
