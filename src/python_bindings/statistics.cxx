#define PY_ARRAY_UNIQUE_SYMBOL geomap_PyArray_API
#define NO_IMPORT_ARRAY
#include <vigra/numpy_array.hxx>
#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <vigra/splineimageview.hxx>


#include "vigra/polygon.hxx"
#include "python_types.hxx"

#include <vector>
#include <algorithm>

using namespace vigra;
using namespace boost::python;

class PolylineStatistics
{
  public:
    PolylineStatistics()
    : weightedSum_(0),
      length_(0),
      min_(vigra::NumericTraits<double>::max()),
      max_(vigra::NumericTraits<double>::min())
    {}

//     PolylineStatistics(const PointArray<Vector2> &poly, object siv)
//     : weightedSum_(0),
//       length_(0),
//       min_(vigra::NumericTraits<double>::max()),
//       max_(vigra::NumericTraits<double>::min())
//     {
//         for(unsigned int i = 0; i < poly.size() - 1; ++i)
//         {
//             Vector2 segment(poly[i+1]-poly[i]);
//             __call__(extract<double>(siv[poly[i] + 0.5*segment])(),
//                      segment.magnitude());
//         }
//     }

    PolylineStatistics(const PointArray<Vector2> &poly,
                       const SplineImageView<5, GrayValue> &siv)
    : weightedSum_(0),
      length_(0),
      min_(vigra::NumericTraits<double>::max()),
      max_(vigra::NumericTraits<double>::min())
    {
        for(unsigned int i = 0; i < poly.size() - 1; ++i)
        {
            Vector2
                segment(poly[i+1]-poly[i]),
                midPos(poly[i] + 0.5*segment);
            __call__(siv(midPos[0], midPos[1]), segment.magnitude());
        }
    }

    void __call__(double value, double length)
    {
        weightedSum_ += value*length;
        if(min_ > value)
            min_ = value;
        if(max_ < value)
            max_ = value;
        length_ += length;
    }

    double average() const
    {
        if(length_)
            return weightedSum_ / length_;
        return 0;
    }

    double min() const
    {
        return min_;
    }

    double max() const
    {
        return max_;
    }

    void merge(const PolylineStatistics &otherStats)
    {
        weightedSum_ += otherStats.weightedSum_;
        if(min_ > otherStats.min_)
            min_ = otherStats.min_;
        if(max_ < otherStats.max_)
            max_ = otherStats.max_;
        length_ += otherStats.length_;
    }

  protected:
    friend class PolylineStatisticsPickleSuite;

    double weightedSum_, length_, min_, max_;
};

class PolylineStatisticsPickleSuite : public boost::python::pickle_suite
{
  public:
    static tuple getstate(const PolylineStatistics& w)
    {
        return make_tuple(w.weightedSum_, w.length_, w.min_, w.max_);
    }

    static void setstate(PolylineStatistics& w, boost::python::tuple state)
    {
        w.weightedSum_ = extract<double>(state[0])();
        w.length_ = extract<double>(state[1])();
        w.min_ = extract<double>(state[2])();
        w.max_ = extract<double>(state[3])();
    }
};

/********************************************************************/

class QuantileStatistics : public PolylineStatistics
{
    typedef std::vector<std::pair<double, double> > Segments;

  public:
    QuantileStatistics()
    {}

    QuantileStatistics(const PointArray<Vector2> &poly,
                       const SplineImageView<5, GrayValue> &siv)
    {
        for(unsigned int i = 0; i < poly.size() - 1; ++i)
        {
            Vector2
                segment(poly[i+1]-poly[i]),
                midPos(poly[i] + 0.5*segment);
            __call__(siv(midPos[0], midPos[1]), segment.magnitude());
        }
    }

    void __call__(double value, double length)
    {
        PolylineStatistics::__call__(value, length);
        segments_.push_back(std::make_pair(value, length));
        sorted_ = false;
    }

    double quantile(double quantile) const
    {
        vigra_precondition(segments_.size() > 0, "empty polygon?");
        ensureOrdering();
        if(quantile == 1.0) // numeric errors may lead to below vigra_fail
            return segments_[segments_.size()-1].first;
        double partialLength = 0;
        for(unsigned int i = 0; i < segments_.size(); ++i)
        {
            partialLength += segments_[i].second;
            if(partialLength >= length_ * quantile)
                return segments_[i].first;
        }
        vigra_fail("quantile > 1.0 ??");
        return 42; // never reached
    }

    void merge(const QuantileStatistics &otherStats)
    {
        PolylineStatistics::merge(otherStats);
        Segments::size_type oldSize(segments_.size());
        segments_.resize(segments_.size() + otherStats.segments_.size());
        std::copy(otherStats.segments_.begin(), otherStats.segments_.end(),
                  segments_.begin() + oldSize);
        sorted_ = false;
    }

  protected:
    friend class QuantileStatisticsPickleSuite;

    void ensureOrdering() const
    {
        if(!sorted_)
        {
            std::sort(segments_.begin(), segments_.end());
            sorted_ = true;
        }
    }

    mutable Segments segments_;
    mutable bool sorted_;
};

class QuantileStatisticsPickleSuite : public boost::python::pickle_suite
{
  public:
    static tuple getstate(const QuantileStatistics& w)
    {
        tuple result = PolylineStatisticsPickleSuite::getstate(w);

        list segments;
        for(unsigned int i = 0; i < w.segments_.size(); ++i)
            segments.append(make_tuple(w.segments_[i].first,
                                       w.segments_[i].second));
        return extract<tuple>(result + make_tuple(segments))();
    }

    static void setstate(QuantileStatistics& w, boost::python::tuple state)
    {
        PolylineStatisticsPickleSuite::setstate(w, state);
        list segments = extract<list>(state[-1])();
        w.segments_.resize(len(segments));
        for(unsigned int i = 0; i < w.segments_.size(); ++i)
        {
            //tuple valueLength(extract<tuple>(segments[i])());
            w.segments_[i].first = extract<double>(segments[i][0])();
            w.segments_[i].second = extract<double>(segments[i][1])();
        }
    }
};

void defStatistics()
{
    class_<PolylineStatistics>("PolylineStatistics")
        .def(init<const PointArray<Vector2> &,
                  const SplineImageView<5, GrayValue> &>())
        .def("__call__", &PolylineStatistics::__call__)
        .def("average", &PolylineStatistics::average)
        .def("min", &PolylineStatistics::min)
        .def("max", &PolylineStatistics::max)
        .def("merge", &PolylineStatistics::merge)
        .def_pickle(PolylineStatisticsPickleSuite())
    ;

    class_<QuantileStatistics, bases<PolylineStatistics> >("QuantileStatistics")
        .def(init<const PointArray<Vector2> &,
                  const SplineImageView<5, GrayValue> &>())
        .def("__call__", &QuantileStatistics::__call__)
        .def("quantile", &QuantileStatistics::quantile)
        .def("merge", &QuantileStatistics::merge)
        .def_pickle(QuantileStatisticsPickleSuite())
    ;
}
