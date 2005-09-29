#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <vigra/pythonvector.hxx>
#include <vigra/splineimageview.hxx>

#include <vigra/pythonimage.hxx>

#include "vigra/polygon.hxx"

#include <vector>
#include <algorithm>

using namespace vigra;
using namespace boost::python;

class PolylineStatistics
{
  public:
    PolylineStatistics()
    : weightedSum_(0),
      length_(0)
    {}

    PolylineStatistics(const PointArray<Vector2> &poly, object siv)
    : weightedSum_(0),
      length_(0)
    {
        for(unsigned int i = 0; i < poly.size() - 1; ++i)
        {
            Vector2 segment(poly[i+1]-poly[i]);
            __call__(extract<double>(siv[poly[i] + 0.5*segment])(),
                     segment.magnitude());
        }
    }

    PolylineStatistics(const PointArray<Vector2> &poly,
                       const SplineImageView<5, GrayValue> &siv)
    : weightedSum_(0),
      length_(0)
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
        length_ += length;
    }

    double average() const
    {
        if(length_)
            return weightedSum_ / length_;
        return 0;
    }

    void merge(const PolylineStatistics &otherStats)
    {
        weightedSum_ += otherStats.weightedSum_;
        length_ += otherStats.length_;
    }

  protected:
    double weightedSum_, length_;
};

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
        ensureOrdering();
        double partialLength = 0;
        for(unsigned int i = 0; i < segments_.size(); ++i)
        {
            partialLength += segments_[i].second;
            if(partialLength >= length_ * quantile)
                return segments_[i].first;
        }
        vigra_fail("quantile > 1.0 ??");
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

    friend class QuantilePickleSuite;

  protected:
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

struct QuantilePickleSuite : boost::python::pickle_suite
{
    static tuple getstate(const QuantileStatistics& w)
    {
        list segments;
        for(unsigned int i = 0; i < w.segments_.size(); ++i)
            segments.append(make_tuple(w.segments_[i].first,
                                       w.segments_[i].second));
        return make_tuple(w.weightedSum_, w.length_, segments);
    }

    static void setstate(QuantileStatistics& w, boost::python::tuple state)
    {
        w.weightedSum_ = extract<double>(state[0])();
        w.length_ = extract<double>(state[1])();
        w.segments_.resize(len(state[2]));
        for(unsigned int i = 0; i < w.segments_.size(); ++i)
        {
            //tuple valueLength(extract<tuple>(state[2][i])());
            w.segments_[i].first = extract<double>(state[2][i][0])();
            w.segments_[i].second = extract<double>(state[2][i][1])();
        }
    }
};

void defStatistics()
{
    class_<PolylineStatistics>("EdgeAverage")
        .def(init<const PointArray<Vector2> &,
                  const SplineImageView<5, GrayValue> &>())
        .def("__call__", &PolylineStatistics::__call__)
        .def("average", &PolylineStatistics::average)
        .def("merge", &PolylineStatistics::merge)
    ;

    class_<QuantileStatistics, bases<PolylineStatistics> >("EdgeStatistics")
        .def(init<const PointArray<Vector2> &,
                  const SplineImageView<5, GrayValue> &>())
        .def("__call__", &QuantileStatistics::__call__)
        .def("quantile", &QuantileStatistics::quantile)
        .def("merge", &QuantileStatistics::merge)
        .def_pickle(QuantilePickleSuite())
    ;
}
