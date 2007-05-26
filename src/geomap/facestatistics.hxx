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

    bool preMergeFaces(const GeoMap::Dart &dart)
    {
//         ssLeft = self._superSampled[dart.leftFaceLabel()]
//         ssRight = self._superSampled[dart.rightFaceLabel()]
//         self.ssMerged = 0 # hopefully, the merged face has no supersampling
//         if ssLeft and (ssLeft == ssRight):
//             self.ssMerged = ssLeft
//             ssLeft, ssRight = 0, 0 # can be merged
//         if ssLeft and ssRight:
//             if ssLeft < ssRight:
//                 self.ssMerged = ssLeft
//                 ssLeft = 0 # take stats from face with less supersampling
//             else:
//                 self.ssMerged = ssRight
//                 ssRight = 0
        merged_ = *functors_[dart.leftFaceLabel()];
        merged_(*functors_[dart.rightFaceLabel()]);
        return true;
    }

    void postMergeFaces(GeoMap::Face &face)
    {
        *functors_[face.label()] = merged_;
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
    GeoMap &map_;
    std::vector<sigc::connection> connections_;
    const OriginalImage &originalImage_;
    std::vector<Functor *> functors_;
    Functor merged_;
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

    connections_.push_back(
        map.preMergeFacesHook.connect(
            sigc::mem_fun(this, &FaceColorStatistics::preMergeFaces)));
    connections_.push_back(
        map.postMergeFacesHook.connect(
            sigc::mem_fun(this, &FaceColorStatistics::postMergeFaces)));
    connections_.push_back(
        map.associatePixelsHook.connect(
            sigc::mem_fun(this, &FaceColorStatistics::associatePixels)));
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

#endif // VIGRA_FACESTATISTICS_HXX
