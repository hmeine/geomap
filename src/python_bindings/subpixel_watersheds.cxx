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

#include "vigra/subpixel_watershed.hxx"

#include "python_types.hxx"
#include <boost/python.hpp>
namespace python = boost::python;

using namespace vigra;

template <class SPWSType>
class SPWSWrapper : public SPWSType
{
  public:
    typedef Vector2 Coordinate;

    SPWSWrapper(NumpyFImage const & img)
    : SPWSType(srcImageRange(img))
    {}

    void
    findCriticalPointsIf(NumpyFImage const & mask)
    {
        this->findCriticalPoints(srcImage(mask));
    }

    python::list
    saddles()
    {
        if(this->saddles_.size() == 0)
            this->findCriticalPoints();
        python::list plist;
        plist.append(python::object());
        for(unsigned int i=1; i<this->saddles_.size(); ++i)
        {
            plist.append(Coordinate(this->saddles_[i][0], this->saddles_[i][1]));
        }
        return plist;
    }

    python::list
    maxima()
    {
        if(this->maxima_.size() == 0)
            this->findCriticalPoints();
        python::list plist;
        plist.append(python::object());
        for(unsigned int i=1; i<this->maxima_.size(); ++i)
        {
            plist.append(Coordinate(this->maxima_[i][0], this->maxima_[i][1]));
        }
        return plist;
    }

    python::list
    minima()
    {
        if(this->minima_.size() == 0)
            this->findCriticalPoints();
        python::list plist;
        plist.append(python::object());
        for(unsigned int i=1; i<this->minima_.size(); ++i)
        {
            plist.append(Coordinate(this->minima_[i][0], this->minima_[i][1]));
        }
        return plist;
    }

    python::object
    edge(int index)
    {
        double x = this->saddles_[index][0];
        double y = this->saddles_[index][1];

        typename SPWSType::PointArray
            forwardCurve, backwardCurve;

        int findex = this->flowLine(x, y, true, 1e-4, forwardCurve);
        int bindex = this->flowLine(x, y, false, 1e-4, backwardCurve);

        python::list ppoints;
        for(int i = forwardCurve.size() - 1; i >= 0; --i)
            ppoints.append(forwardCurve[i]);
        for(int i = 1; i < (int)backwardCurve.size(); ++i)
            ppoints.append(backwardCurve[i]);

        return python::make_tuple(findex, bindex, ppoints, forwardCurve.size() - 1);
    }

    python::list
    edges(double threshold)
    {
        if(this->saddles_.size() == 0)
            this->findCriticalPoints();
        python::list plist;
        for(unsigned int i=1; i<this->saddles_.size(); ++i)
        {
            if(this->image_(this->saddles_[i][0], this->saddles_[i][1]) >= threshold)
                plist.append(edge(i));
            else
                plist.append(python::object());
        }
        return plist;
    }

    python::tuple
    debugCP()
    {
        python::list minima, saddles, maxima;
        //typename SPWSType::PointArray minima, sources;
        //std::vector<unsigned int> counts;

        int w = this->image_.width();
        int h = this->image_.height();

        typedef typename SPWSType::SplineImageView::value_type Value;

        double d = 1.0 / this->cpOversampling_;

        // search for critical points
        int percent = -1, lastPercent = -1;
        for(int y=1; y <= h-2; ++y)
        {
            percent = 100 * y / h;
            if(percent != lastPercent)
            {
                std::cerr << "debugCP(): " << percent << "%\r";
                lastPercent = percent;
            }
            for(int x=1; x <= w-2; ++x)
            {
                for(double dy = 0.0; dy < 1.0; dy += d)
                {
                    for(double dx = 0.0; dx < 1.0; dx += d)
                    {
                        double xx, yy;

                        CriticalPoint type = findCriticalPointNewtonMethod(
                            this->image_, x + dx, y + dy,
                            &xx, &yy, this->stepEpsilon_);
                        if(type == Failed)
                            continue;

                        (type == Saddle ? saddles
                         : (type == Minimum ? minima : maxima))
                               .append(python::make_tuple(
                                           Coordinate(xx, yy),
                                           Coordinate(x + dx, y + dy)));
                    }
                }
            }
        }
        std::cerr << "debugCP(): done.\n";

        return python::make_tuple(minima, saddles, maxima);
    }

    python::tuple
    findCriticalPointsInFacet(int x, int y)
    {
        typename SPWSType::PointArray minima, saddles, maxima;
        findCriticalPointsInFacet(x, y, minima, saddles, maxima);
        python::list pminima, pmaxima, psaddles;
        for(unsigned int i=0; i<minima.size(); ++i)
        {
            pminima.append(minima[i]);
        }
        for(unsigned int i=0; i<saddles.size(); ++i)
        {
            psaddles.append(saddles[i]);
        }
        for(unsigned int i=0; i<maxima.size(); ++i)
        {
            pmaxima.append(maxima[i]);
        }
        return python::make_tuple(pminima, psaddles, pmaxima);
    }
};

template <class SPWSType>
python::class_<SPWSWrapper<SPWSType> > &
defSubPixelWS(char const * name)
{
    typedef SPWSWrapper<SPWSType> SPWS;

    static python::class_<SPWS> theclass(
        name, python::init<NumpyFImage const &>());
    theclass
        .def("width", &SPWS::width)
        .def("height", &SPWS::height)
        .def("saddles", &SPWS::saddles)
        .def("maxima", &SPWS::maxima)
        .def("minima", &SPWS::minima)
        .def("debugCP", &SPWS::debugCP)
        .def("edge", &SPWS::edge)
        .def("edges", &SPWS::edges)
        .def("findCriticalPoints", (void(SPWS::*)())&SPWS::findCriticalPoints)
        .def("findCriticalPoints", &SPWS::findCriticalPointsIf)
        .def_readwrite("initialStep", &SPWS::initialStep_)
        .def_readwrite("minCPDist", &SPWS::minCPDist_)
        .def_readwrite("cpOversampling", &SPWS::cpOversampling_)
        //.def("findCriticalPointsInFacet", &SPWS::findCriticalPointsInFacet)
    ;

    return theclass;
}

void defSPWS()
{
    defSubPixelWS<SubPixelWatersheds<SplineImageView<2, GrayValue> > >(
        "SubPixelWatersheds2");
    defSubPixelWS<SubPixelWatersheds<SplineImageView<3, GrayValue> > >(
        "SubPixelWatersheds3");
    defSubPixelWS<SubPixelWatersheds<SplineImageView<5, GrayValue> > >(
        "SubPixelWatersheds5");
}
