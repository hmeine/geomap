#ifndef VIGRA_FACESTATISTICS_HXX
#define VIGRA_FACESTATISTICS_HXX

#include "cppmap.hxx"
#include <vigra/inspectimage.hxx>
#include <vigra/transformimage.hxx>
#include <vigra/splineimageview.hxx>
#include <boost/bind.hpp>
#include <cmath>
#include <iostream>

#ifdef HAVE_MATH_TOOLKIT
#include <boost/math/distributions/students_t.hpp>
#endif

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
        {
            vigra_assert(
                (unsigned int)label < stats_.size() && stats_[(unsigned int)label],
                "invalid functor index during initalization");
            (*stats_[(unsigned int)label])(v);
        }
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
        {
            vigra_assert(
                (unsigned int)label < stats_.size() && stats_[(unsigned int)label],
                "invalid functor index in lookup functor");
            return stats_[label]->average();
        }
        return typename FACE_STATS::result_type();
    }
};

inline vigra::Point2D intPos(const Vector2 &p)
{
    return vigra::Point2D((int)floor(p[0]+0.5), (int)floor(p[1]+0.5));
}

inline vigra::Rect2D intPos(const GeoMap::Face::BoundingBox &b)
{
    return vigra::Rect2D(intPos(b.begin()),
                         intPos(b.end())+vigra::Size2D(1,1));
}

} // namespace detail

/********************************************************************/

template<class Scalar>
Scalar clampedScalarVariance(Scalar s)
{
    return std::max(s, (Scalar)(1/12.));
}

template<class Scalar, int Dim>
Scalar clampedScalarVariance(vigra::TinyVector<Scalar, Dim> s)
{
    Scalar result(0.0);
    for(int i = 0; i < Dim; ++i)
        result += clampedScalarVariance(s[i]);
    return result;
}

template<class OriginalImage>
class FaceColorStatistics
{
  public:
    typedef vigra::FindAverageAndVariance<typename OriginalImage::value_type>
        Functor;

//     template<int SplineOrder>
    FaceColorStatistics(boost::shared_ptr<GeoMap> map, const OriginalImage &originalImage,
                        double maxDiffNorm = 1.0, unsigned int minSampleCount = 1);
    FaceColorStatistics(const FaceColorStatistics &other);
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
        vigra_assert(
            dart.leftFaceLabel() < size() && functors_[dart.leftFaceLabel()] &&
            dart.rightFaceLabel() < size() && functors_[dart.rightFaceLabel()],
            "invalid face labels in preMergeFaces");

        label1_ = dart.leftFaceLabel();
        label2_ = dart.rightFaceLabel();

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
                    if(dart.rightFace()->pixelArea())
                        rescanFace(*dart.rightFace(), merged_);
                }
                else
                {
                    merged_ = *functors_[dart.rightFaceLabel()];
                    mergedSS_ = ssRight;
                    if(dart.leftFace()->pixelArea())
                        rescanFace(*dart.leftFace(), merged_);
                }
                return true;
            }
            else if(ssLeft) // ssLeft == ssRight
            {
                mergedSS_ = ssLeft;

                if(dart.leftFace()->pixelArea() + dart.rightFace()->pixelArea()
                   >= minSampleCount_)
                {
                    merged_.reset();
                    rescanFace(*dart.leftFace(), merged_);
                    rescanFace(*dart.rightFace(), merged_);
                    mergedSS_ = 0;

                    return true;
                }
            }
            else
                mergedSS_ = 0;
        }

        merged_ = *functors_[dart.leftFaceLabel()];
        merged_(*functors_[dart.rightFaceLabel()]);
        return true;
    }

    void postMergeFaces(GeoMap::Face &face)
    {
        vigra_assert(
            face.label() < size() && functors_[face.label()],
            "invalid survivor label in postMergeFaces");

        CellLabel mergedLabel = (
            face.label() == label1_ ? label2_ : label1_);

        *functors_[face.label()] = merged_;
        if(superSampled_.get())
        {
            if((*superSampled_)[face.label()] && !mergedSS_)
            {
                if(!--superSampledCount_)
                    superSampled_.reset();
            }
            (*superSampled_)[face.label()] = mergedSS_;

            if((*superSampled_)[mergedLabel])
            {
                if(!--superSampledCount_)
                    superSampled_.reset();
            }
        }

        delete functors_[mergedLabel];
        functors_[mergedLabel] = NULL;

#ifndef NDEBUG
        if(!superSampled_.get() || !(*superSampled_)[face.label()])
            vigra_postcondition(face.pixelArea() == functors_[face.label()]->count(),
                                "wrong pixelcount after merge in faceMeans");
#endif
    }

    void associatePixels(const GeoMap::Face &face, const PixelList &pixels)
    {
        vigra_assert(
            face.label() < size() && functors_[face.label()],
            "invalid survivor label in associatePixels");

        if(superSampled_.get() && (*superSampled_)[face.label()] &&
           face.pixelArea() >= minSampleCount_)
        {
            (*superSampled_)[face.label()] = 0;
            --superSampledCount_;

            functors_[face.label()]->reset();

            if(face.pixelArea() != pixels.size())
            {
                // we need to re-scan due to previous supersampling:
                rescanFace(face, *functors_[face.label()]);
                return;
            }
        }

        Functor &f(*functors_[face.label()]);
        for(PixelList::const_iterator it = pixels.begin();
            it != pixels.end(); ++it)
        {
#ifndef NDEBUG
            if(originalImage_.isInside(
                   typename OriginalImage::difference_type((*it).x, (*it).y)))
#endif
                f(originalImage_((*it).x,(*it).y));
        }

#ifndef NDEBUG
        if(!superSampled_.get() || !(*superSampled_)[face.label()])
            vigra_postcondition(face.pixelArea() == functors_[face.label()]->count(),
                                "wrong pixelcount after pixel association in faceMeans");
#endif
    }

    double faceMeanDiff(const GeoMap::Dart &dart) const
    {
        return vigra::norm(average(dart.leftFaceLabel()) -
                           average(dart.rightFaceLabel()))
            / maxDiffNorm_;
    }

    double faceAreaHomogeneity(const GeoMap::Dart &dart) const
    {
        double a1 = dart.leftFace()->area();
        double a2 = dart.rightFace()->area();
        return (a1 * a2) / (a1 + a2);
    }

    double faceHomogeneity2(const GeoMap::Dart &dart) const
    {
        return vigra::squaredNorm(faceMeanDiff(dart))
            * faceAreaHomogeneity(dart);
    }

    double faceHomogeneity(const GeoMap::Dart &dart) const
    {
        return std::sqrt(faceHomogeneity2(dart));
    }

#ifdef HAVE_MATH_TOOLKIT
    double faceTTest(const GeoMap::Dart &dart) const
    {
        double
            sigma1_2 = clampedScalarVariance(variance(dart.leftFaceLabel(), true)),
            sigma2_2 = clampedScalarVariance(variance(dart.rightFaceLabel(), true)),
            N1 = pixelCount(dart.leftFaceLabel()),
            N2 = pixelCount(dart.rightFaceLabel()),
#ifdef UNEQUAL_VARIANCES
            // Welch's t-test
            sigmaExpr = (sigma1_2 / N1 + sigma2_2 / N2),
            t = (vigra::norm(average(dart.leftFaceLabel()) -
                             average(dart.rightFaceLabel())) /
                 std::sqrt(sigmaExpr)),
            // this formula (Welch-Satterthwaite approximation) does
            // *not* assume equal variances and results in non-integer
            // DOF:
            dof = vigra::sq(sigmaExpr) / (
                vigra::sq(sigma1_2/N1)/(N1-1) + vigra::sq(sigma2_2/N2)/(N2-1));
#else
            // Student's t-test
            dof = N1 + N2 - 2,
            pooledSigma_2 = ((N1-1)*sigma1_2 + (N2-1)*sigma2_2) / dof,
            t = (vigra::norm(average(dart.leftFaceLabel()) -
                             average(dart.rightFaceLabel())) /
                 std::sqrt(pooledSigma_2 * (1.0 / N1 + 1.0 / N2)));
#endif

        using namespace boost::math;
        students_t stud(dof);
        return 2*cdf(stud, t)-1.0;
    }
#endif

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

    const boost::shared_ptr<GeoMap> map() const
    {
        return map_;
    }

    void attachHooks(boost::shared_ptr<GeoMap> map)
    {
        vigra_precondition(!connections_.size(),
            "FaceColorStatistics: attachHooks() called with existing connections");
        map_ = map;
        connections_.push_back(
            map->preMergeFacesHook.connect(
                boost::bind(boost::mem_fn(&FaceColorStatistics::preMergeFaces), this, _1)));
        connections_.push_back(
            map->postMergeFacesHook.connect(
                boost::bind(boost::mem_fn(&FaceColorStatistics::postMergeFaces), this, _1)));
        connections_.push_back(
            map->associatePixelsHook.connect(
                boost::bind(boost::mem_fn(&FaceColorStatistics::associatePixels), this, _1, _2)));
    }

    void detachHooks()
    {
        for(unsigned int i = 0; i < connections_.size(); ++i)
            if(connections_[i].connected())
                connections_[i].disconnect();
        connections_.clear();
    }

    unsigned int minSampleCount() const
    {
        return minSampleCount_;
    }

    unsigned int superSampledCount() const
    {
        return superSampledCount_;
    }

    bool checkConsistency() const
    {
        bool result = true;

        if(!map_)
        {
            std::cerr << "FaceColorStatistics not attached to any GeoMap!\n";
            result = false;
        }

        if(map_ && superSampled_.get())
        {
            unsigned int realSuperSampledCount = 0;
            for(GeoMap::FaceIterator it = map_->facesBegin(); it.inRange(); ++it)
            {
                if((*superSampled_)[(*it)->label()])
                    ++realSuperSampledCount;
                else if(pixelCount((*it)->label()) != (*it)->pixelArea())
                {
                    std::cerr << "Face " << (*it)->label() << " has a wrong number of samples (" << (*it)->label() << " instead of pixelArea() == " << (*it)->pixelArea() << ")\n";
                    result = false;
                }
            }

            if(realSuperSampledCount != superSampledCount_)
            {
                std::cerr << "superSampledCount() is wrong ("
                          << superSampledCount_ << ", should be "
                          << realSuperSampledCount << ")\n";
                result = false;
            }
        }

        if((bool)superSampled_.get() != (bool)superSampledCount_)
        {
            std::cerr << "superSampledCount() is "
                      << superSampledCount_ << " but ss. array is "
                      << superSampled_.get() << "\n";
            result = false;
        }

        return result;
    }

  protected:
    void ensureMinSampleCount(unsigned int minSampleCount);

    void rescanFace(const GeoMap::Face &face, Functor &f)
    {
        vigra::Rect2D faceROI(detail::intPos(face.boundingBox()));
        faceROI &= vigra::Rect2D(vigra::Size2D(originalImage_.shape(0),originalImage_.shape(1)));

        GeoMap::LabelImageIterator
            labUL = map_->labelsUpperLeft(),
            ul = labUL + faceROI.upperLeft(),
            lr = labUL + faceROI.lowerRight();
        GeoMap::LabelImageAccessor lac(map_->labelAccessor());

        for(; ul.y < lr.y; ++ul.y)
        {
            GeoMap::LabelImageIterator lab(ul);
            for(; lab.x < lr.x; ++lab.x)
                if(lac(lab) == (int)face.label())
                    f(originalImage_((lab-labUL).x, (lab-labUL).y));
        }
    }

    boost::shared_ptr<GeoMap> map_;
    std::vector<boost::signals::connection> connections_;
    const OriginalImage originalImage_;
    std::vector<Functor *> functors_;

    unsigned int minSampleCount_;
    std::auto_ptr<std::vector<unsigned char> > superSampled_;
    unsigned int superSampledCount_;

    Functor merged_;
    unsigned char mergedSS_;
    CellLabel label1_, label2_;

    double maxDiffNorm_;
};

template<class OriginalImage>
// template<int SplineOrder>
FaceColorStatistics<OriginalImage>::FaceColorStatistics(
    boost::shared_ptr<GeoMap> map, const OriginalImage &originalImage,
    double maxDiffNorm, unsigned int minSampleCount)
: map_(map),
  originalImage_(originalImage),
  functors_(map->maxFaceLabel(), NULL),
  minSampleCount_(minSampleCount),
  superSampledCount_(0),
  maxDiffNorm_(maxDiffNorm)
{
    for(GeoMap::FaceIterator it = map->facesBegin(); it.inRange(); ++it)
        functors_[(*it)->label()] = new Functor();

    detail::InitFaceFunctors<Functor> iff(functors_);
    inspectTwoImages(map->srcLabelRange(),
                     srcImage(originalImage),
                     iff);

    if(minSampleCount)
        ensureMinSampleCount(minSampleCount);

    attachHooks(map);
}

template<class OriginalImage>
FaceColorStatistics<OriginalImage>::FaceColorStatistics(
    const FaceColorStatistics &other)
: originalImage_(other.originalImage_),
  functors_(other.functors_),
  minSampleCount_(other.minSampleCount_),
  superSampled_(other.superSampled_.get()
                ? (new std::vector<unsigned char>(*other.superSampled_))
                : NULL),
  superSampledCount_(other.superSampledCount_),
  maxDiffNorm_(other.maxDiffNorm_)
{
}

template<class OriginalImage>
void FaceColorStatistics<OriginalImage>::ensureMinSampleCount(
    unsigned int minSampleCount)
{
    std::auto_ptr<std::vector<unsigned char> > superSampled;
    typedef typename vigra::NumericTraits<
        typename OriginalImage::value_type>::RealPromote FloatPixel;
    std::auto_ptr<vigra::SplineImageView<5, FloatPixel> > siv;

    for(GeoMap::FaceIterator it = map_->finiteFacesBegin();
        it.inRange(); ++it)
    {
        if(functors_[(*it)->label()]->count() < minSampleCount)
        {
            if(!superSampled.get())
            {
                superSampled.reset(
                    new std::vector<unsigned char>(map_->maxFaceLabel(), 0));
                siv.reset(
                    new vigra::SplineImageView<5, FloatPixel>(
                        srcImageRange(originalImage_)));
            }

            do
            {
                unsigned char ss(++(*superSampled)[(*it)->label()]);
                if(ss > 16)
                {
                    std::cerr << "WARNING: FaceColorStatistics giving up sampling face " << (*it)->label() << " (area " << (*it)->area() << ")\n";
                    break;
                }
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
        map_->srcLabelRange(), destIter(dul, da));
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
