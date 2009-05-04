/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2003 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA_SUBPIXEL_WATERSHED_HXX
#define VIGRA_SUBPIXEL_WATERSHED_HXX

#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include "vigra/mathutil.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/polynomial.hxx"
#include "vigra/eigensystem.hxx"
#include "vigra/splineimageview.hxx"
#include "vigra/pixelneighborhood.hxx"
#include "positionedmap.hxx"

namespace vigra {

enum CriticalPoint { Minimum = -1, Saddle, Maximum, Failed = 999 };

template <class IMAGE, class VECTOR>
void findCriticalPoints48Neighborhood(IMAGE const &image,
                        VECTOR *minima, VECTOR *saddles, VECTOR *maxima,
                        bool eightneighborhood)
{
    typedef typename IMAGE::value_type Value;
    typedef typename IMAGE::const_traverser Traverser;
    typedef typename VECTOR::value_type Coordinate;

    int w = image.width();
    int h = image.height();
    int d = eightneighborhood ? 1 : 2;
    int steps = eightneighborhood ? 8 : 4;

    Traverser iy = image.upperLeft() + Diff2D(1,1);

    for(int y=1; y < h-1; ++y, ++iy.y)
    {
        Traverser ix = iy;
        for(int x=1; x < w-1; ++x, ++ix.x)
        {
            NeighborhoodCirculator<Traverser, EightNeighborCode> c(ix);
            Value v = *ix, z = NumericTraits<Value>::zero(), s;
            int i = 0;
            do
            {
                s = sign(v - *c);
                if(s != z)
                    break;
                c += d;
            }
            while(++i < steps);
            if(i == steps)
                continue; // plateau
            int countSignChanges = 0, countZeros = 0;
            for(i = 0, c += d; i < steps; ++i, c += d)
            {
                Value ss = sign(v - *c);
                if(s * ss < z)
                {
                    ++countSignChanges;
                    s = ss;
                }
                else if(ss == z)
                    ++countZeros;
            }

            if(countSignChanges >= 4)
            {
                saddles->push_back(Coordinate(x, y));
            }
            else if(countSignChanges == 0 && countZeros == 0 && s < z)
            {
                minima->push_back(Coordinate(x, y));
            }
            else if(countSignChanges == 0 && countZeros == 0 && s > z)
            {
                maxima->push_back(Coordinate(x, y));
            }
        }
    }
}

template <class IMAGE, class VECTOR>
void findCriticalPointsLinearInterpolation(IMAGE const & image,
                        VECTOR *minima, VECTOR *saddles, VECTOR *maxima)
{
    typedef typename IMAGE::value_type Value;
    typedef typename IMAGE::const_traverser Traverser;
    typedef typename VECTOR::value_type Coordinate;

    int w = image.width();
    int h = image.height();

    Traverser iy = image.upperLeft() + Diff2D(1,1);

    // The linear interpolation case is tricky because there are many possibiliries to
    // create a saddle in the interpolated function:
    //      * saddle configuration in the 4-neighborhood (this is independent of the 8-neighbor values!)
    //      * all 4-neighbors are equal to the center, but the 8-neighbors are in a saddle configuration
    //      * the 4 corners of a facet are in a saddle (checker-board) configuration.
    // It can be shown that other configurations cannot lead to a saddle.
    // Extrema can only be encountered at sampling points, i.e. at facet break points,
    // because the facet function a*x*y + b*x + c*y + d cannot have extrema.
    for(int y=1; y < h-1; ++y, ++iy.y)
    {
        Traverser ix = iy;
        for(int x=1; x < w-1; ++x, ++ix.x)
        {
            NeighborhoodCirculator<Traverser, EightNeighborCode> c(ix);
            Value v = *ix, z = NumericTraits<Value>::zero(), s;
            int i, countSame = 0;
            for(i = 0; i < 4; ++i, c += 2)
            {
                s = sign(v - *c);
                if(s == z)
                    ++countSame;
            }
            switch(countSame)
            {
              case 0:
              {
                int countSignChanges = 0;
                for(i = 0; i < 4; ++i, c += 2)
                {
                    Value ss = sign(v - *c);
                    if(s * ss < z)
                    {
                        ++countSignChanges;
                        s = ss;
                    }
                }
                if(countSignChanges == 4)
                {
                    saddles->push_back(Coordinate(x, y));
                }
                else if(countSignChanges == 0 && s < z)
                {
                    minima->push_back(Coordinate(x, y));
                }
                else if(countSignChanges == 0 && s > z)
                {
                    maxima->push_back(Coordinate(x, y));
                }
                break;
              }
              case 4:
#if 0 // ignore plateaus for the time being
              {
                // 4-neighbor plateau => check for 8-neighbor saddle
                for(i = 0, c += 1; i < 4; ++i, c += 2)
                {
                    s = sign(v - *c);
                    if(s != z)
                        break;
                }
                if(i == 4)
                    continue; // true plateau
                int countSignChanges = 0;
                for(i = 0, c += 2; i < 4; ++i, c += 2)
                {
                    Value ss = sign(v - *c);
                    if(s * ss < z)
                    {
                        ++countSignChanges;
                        s = ss;
                    }
                }
                if(countSignChanges == 4)
                {
                    saddles->push_back(Coordinate(x, y));
                }
                break;
              }
#endif
              case 1:
              case 2:
              case 3:
                  break;// not implemented yet
            }
            if(x == w - 2 || y == h - 2)
                continue; // don't check the border facets
            // check for saddle (checker board configuration) in the lower right facet
            c.turnTo(EightNeighborhood::South);
            Value v1 = *c; ++c;
            Value v2 = *c; ++c;
            Value v3 = *c;
            if((v > v1 && v > v3 && v2 > v1 && v2 > v3) ||
               (v < v1 && v < v3 && v2 < v1 && v2 < v3))
            {
                Value sv = v - v1 + v2 - v3;
                saddles->push_back(Coordinate(x + (v - v1) / sv, y + (v - v3) / sv));
            }
        }
    }
}

template <class IMAGEVIEW>
CriticalPoint
findCriticalPointNewtonMethod(
    IMAGEVIEW const & image, const double x, const double y,
    double *xx, double *yy, const double stepEpsilon)
{
    typedef typename IMAGEVIEW::value_type Value;
    Value zero = NumericTraits<Value>::zero();

    *xx = x;
    *yy = y;
    double sxx, syy, stepEpsilon2 = stepEpsilon * stepEpsilon;
    for(int i=0; i<100; ++i) // do at most 100 iterations
    {
        Value dx = image.dx(*xx, *yy);
        Value dy = image.dy(*xx, *yy);
        Value dxx = image.dxx(*xx, *yy);
        Value dxy = image.dxy(*xx, *yy);
        Value dyy = image.dyy(*xx, *yy);
        Value det = dxx*dyy - dxy*dxy;

        if(det != zero)
        {
            sxx = (dxy*dy - dyy*dx) / det;
            syy = (dxy*dx - dxx*dy) / det;
            *xx += sxx;
            *yy += syy;
            if(!image.isValid(*xx, *yy))
                return Failed; // coordinates out of range
        }
        else
        {
            sxx = syy = 0.0;
        }

        double dist2 = sxx*sxx + syy*syy;
        if(dist2 < stepEpsilon2) // convergence
        {
            // FIXME: really reuse stepEpsilon for this purpose?
            if(*xx < -stepEpsilon || *xx > (double)(image.width())-1.0+stepEpsilon ||
               *yy < -stepEpsilon || *yy > (double)(image.height())-1.0+stepEpsilon)
            {
                return Failed; // coordinates out of range
            }
            if(det == zero)
            {
                if(dx == zero && dy == zero)
                {
                    return Saddle;
                }
                return Failed; // Hessian singular
            }
            else if(det < zero)
            {
                return Saddle;
            }
            else if(dxx + dyy > zero)
            {
                return Minimum;
            }
            else
            {
                return Maximum;
            }
        }
    }
    return Failed;
}

template <class IMAGEVIEW>
CriticalPoint
findCriticalPointNewtonMethod(
    IMAGEVIEW const & image, const double x, const double y,
    double *xx, double *yy, const double stepEpsilon,
     const double squaredSearchRadius)
{
    typedef typename IMAGEVIEW::value_type Value;
    Value zero = NumericTraits<Value>::zero();

    *xx = x;
    *yy = y;
    double sxx, syy, stepEpsilon2 = stepEpsilon * stepEpsilon;
    for(int i=0; i<100; ++i) // do at most 100 iterations
    {
        Value dx = image.dx(*xx, *yy);
        Value dy = image.dy(*xx, *yy);
        Value dxx = image.dxx(*xx, *yy);
        Value dxy = image.dxy(*xx, *yy);
        Value dyy = image.dyy(*xx, *yy);
        Value det = dxx*dyy - dxy*dxy;

        if(det != zero)
        {
            sxx = (dxy*dy - dyy*dx) / det;
            syy = (dxy*dx - dxx*dy) / det;
            *xx += sxx;
            *yy += syy;
            if(!image.isValid(*xx, *yy))
                return Failed; // coordinates out of range
        }
        else
        {
            sxx = syy = 0.0;
        }

        double dist2 = sxx*sxx + syy*syy;
        if(dist2 < stepEpsilon2) // convergence
        {
            // FIXME: really reuse stepEpsilon for this purpose?
            if(*xx < -stepEpsilon || *xx > (double)(image.width())-1.0+stepEpsilon ||
               *yy < -stepEpsilon || *yy > (double)(image.height())-1.0+stepEpsilon)
            {
                return Failed; // coordinates out of range
            }
            if(det == zero)
            {
                if(dx == zero && dy == zero)
                {
                    return Saddle;
                }
                return Failed; // Hessian singular
            }
            else if(det < zero)
            {
                return Saddle;
            }
            else if(dxx + dyy > zero)
            {
                return Minimum;
            }
            else
            {
                return Maximum;
            }
        }

        if(sq(*xx-x) + sq(*yy-y) > squaredSearchRadius)
            return Failed;
    }
    return Failed;
}

template <class IMAGEVIEW>
bool
findDirectedMaximum(IMAGEVIEW const & image,
                    double x, double y,
                    double nx, double ny,
                    double *xx, double *yy,
                    double stepEpsilon, double maxDist)
{
    typedef typename IMAGEVIEW::value_type Value;
    Value zero = NumericTraits<Value>::zero();

    double norm(sqrt(nx*nx + ny*ny));
    if(norm)
    {
        nx /= norm;
        ny /= norm;
    }

    *xx = x;
    *yy = y;
    double sxx, syy, stepEpsilon2 = stepEpsilon * stepEpsilon, maxDist2 = maxDist * maxDist;
    for(int i=0; i<100; ++i) // do at most 100 iterations
    {
        Value d = image.dx(x, y)*nx + image.dy(x, y)*ny;
        Value dxx = image.dxx(*xx, *yy);
        Value dxy = image.dxy(*xx, *yy);
        Value dyy = image.dyy(*xx, *yy);
        Value d2 = dxx*nx*nx + 2*dxy*ny*nx + dyy*ny*ny;

        if(d2 != zero)
        {
            double l = -d / d2;
            sxx = l * nx;
            syy = l * ny;
        }
        else
        {
            sxx = syy = 0.0;
        }
        *xx += sxx;
        *yy += syy;
        if(*xx < 1 || *xx > (double)(image.width()-2) ||
           *yy < 1 || *yy > (double)(image.height()-2))
        {
            return false; // coordinates out of range
        }

        if(sxx*sxx + syy*syy < stepEpsilon2) // convergence
            return (d2 < zero) && (dxx + dyy < zero); // maximum?

        sxx = *xx - x; syy = *yy - y;
        if(sxx*sxx + syy*syy > maxDist2) // too far from start
            return false;
    }
    return Failed;
}

struct CriticalPointHolder
{
    double x,y;
    CriticalPoint type;

    CriticalPointHolder(double xx, double yy, CriticalPoint t)
        : x(xx), y(yy), type(t)
    {}

    double operator[](int i) const
    {
        return (&x)[i];
    }

    double squaredDistance(const CriticalPointHolder &other) const
    {
        return sq(other.x-x) + sq(other.y-y);
    }
};

struct CriticalPointsCompare
{
    double epsilon;
    CriticalPointHolder self;

    CriticalPointsCompare(double eps, CriticalPointHolder s)
        : epsilon(eps), self(s)
    {}

    bool operator()(CriticalPointHolder const & l) const
    {
        return closeAtTolerance(l.x, self.x, epsilon) && closeAtTolerance(l.y, self.y, epsilon);
    }
};

template <class IMAGEVIEW, class VECTOR>
void findCriticalPointsNewtonMethod(
    IMAGEVIEW const & image,
    VECTOR *minima, VECTOR *saddles, VECTOR *maxima,
    double minCPDist, double stepEpsilon, unsigned int oversampling)
{
    int w = image.width();
    int h = image.height();

    typedef typename IMAGEVIEW::value_type Value;
    typedef typename VECTOR::value_type Coordinate;

    double d = 1.0 / oversampling;
    double squareMinCPDist = sq(minCPDist);

    Map2D<Coordinate> points;

    // search for critical points
    int percent = -1, lastPercent = -1;
    //std::cerr << "\n";
    for(int y=0; y <= h-1; ++y)
    {
        percent = 100 * y / h;
        if(percent != lastPercent)
        {
            std::cerr << "findCriticalPointsNewtonMethod(): " << percent << "%\r";
            lastPercent = percent;
        }
        for(int x=0; x <= w-1; ++x)
        {
            for(double dy = 0.0; dy < 1.0; dy += d)
            {
                for(double dx = 0.0; dx < 1.0; dx += d)
                {
                    double xx, yy;
                    CriticalPoint type = findCriticalPointNewtonMethod(
                        image, x + dx, y + dy, &xx, &yy, stepEpsilon, 4.0);
                    if(type == Failed)
                        continue;

                    Coordinate c(xx, yy);
                    if(points.nearest(c, squareMinCPDist) != points.end())
                        continue;

                    if(type == Saddle)
                        saddles->push_back(c);
                    else if(type == Minimum)
                        minima->push_back(c);
                    else
                        maxima->push_back(c);

                    points.insert(c);
                }
            }
        }
    }
    std::cerr << "findCriticalPointsNewtonMethod(): done.\n";
}

template <class IMAGEVIEW, class VECTOR, class MaskIterator, class MaskAccessor>
void findCriticalPointsNewtonMethodIf(IMAGEVIEW const & image,
    pair<MaskIterator, MaskAccessor> mask,
    VECTOR *minima, VECTOR *saddles, VECTOR *maxima,
    double minCPDist, double stepEpsilon, unsigned int oversampling)
{
    int w = image.width();
    int h = image.height();

    typedef typename IMAGEVIEW::value_type Value;
    typedef typename VECTOR::value_type Coordinate;

    double d = 1.0 / oversampling;
    double squareMinCPDist = sq(minCPDist);

    Map2D<Coordinate> points;

    // search for critical points
    int percent = -1, lastPercent = -1;
    //std::cerr << "\n";
    MaskIterator maskRow(mask.first);
    for(int y=0; y <= h-1; ++y, ++maskRow.y)
    {
        percent = 100 * y / h;
        if(percent != lastPercent)
        {
            std::cerr << "findCriticalPointsNewtonMethodIf(): " << percent << "%\r";
            lastPercent = percent;
        }
        MaskIterator mit(maskRow);
        for(int x=0; x <= w-1; ++x, ++mit.x)
        {
            if(!mask.second(mit))
                continue;
            for(double dy = 0.0; dy < 1.0; dy += d)
            {
                for(double dx = 0.0; dx < 1.0; dx += d)
                {
                    double xx, yy;
                    CriticalPoint type = findCriticalPointNewtonMethod(
                        image, x + dx, y + dy, &xx, &yy, stepEpsilon, 4.0);
                    if(type == Failed)
                        continue;

                    Coordinate c(xx, yy);
                    if(points.nearest(c, squareMinCPDist) != points.end())
                        continue;

                    if(type == Saddle)
                        saddles->push_back(c);
                    else if(type == Minimum)
                        minima->push_back(c);
                    else
                        maxima->push_back(c);

                    points.insert(c);
                }
            }
        }
    }
    std::cerr << "findCriticalPointsNewtonMethodIf(): done.\n";
}

template <class T, class VECTOR>
void
findCriticalPointsInFacet(
    SplineImageView<2, T> const & s, double x0, double y0,
    VECTOR *minima, VECTOR *saddles, VECTOR *maxima)
{
    typedef typename VECTOR::value_type PointType;

    x0 = VIGRA_CSTD::floor(x0 + 0.5);
    y0 = VIGRA_CSTD::floor(y0 + 0.5);

    DImage splineCoeffs(3,3);
    s.coefficientArray(x0, y0, splineCoeffs);
    //double j = splineCoeffs(0,0);
    double g = splineCoeffs(1,0);
    double e = splineCoeffs(2,0);
    double h = splineCoeffs(0,1);
    double d = splineCoeffs(1,1);
    double b = splineCoeffs(2,1);
    double f = splineCoeffs(0,2);
    double c = splineCoeffs(1,2);
    double a = splineCoeffs(2,2);

    double polyCoeffs[6];
    polyCoeffs[0] =  4.0*f*f*g - 2.0*d*f*h + c*h*h;
    polyCoeffs[1] = -2.0*d*d*f + 8.0*e*f*f + 8.0*c*f*g - 4.0*b*f*h + 2.0*a*h*h;
    polyCoeffs[2] = -c*d*d - 6.0*b*d*f + 16.0*c*e*f + 4.0*c*c*g +
                     8.0*a*f*g - 2.0*b*c*h + 2.0*a*d*h;
    polyCoeffs[3] = -4.0*b*c*d + 8.0*c*c*e - 4.0*b*b*f + 16.0*a*e*f + 8.0*a*c*g;
    polyCoeffs[4] = -3.0*b*b*c - 2.0*a*b*d + 16.0*a*c*e + 4.0*a*a*g;
    polyCoeffs[5] = -2.0*a*b*b + 8.0*a*a*e;

    double eps = 1.0e-7;
    StaticPolynomial<5, double> px(polyCoeffs, 5, 1e-14);
    ArrayVector<double> rx;
    if(!polynomialRealRoots(px, rx))
        return;

    double xold = -100.0;
    for(unsigned int i=0; i < rx.size(); ++i)
    {
        double x = rx[i];

        // ensure that x is in the current facet,
        // and that a multiple root is only used once (this may be
        // wrong, as perhaps several y's share the same x)
        if(std::abs(x) <= 0.5 && std::abs(x-xold) >= eps)
        {
            double cy = 2*f + 2*c*x + 2*a*x*x;
            if(cy == 0.0)
                continue;
            double y = -(h + d*x + b*x*x) / cy;
            if(std::abs(y) <= 0.5)
            {
                double xx = x + x0;
                double yy = y + y0;
                double hxx = s.dxx(xx,yy);
                double hxy = s.dxy(xx,yy);
                double hyy = s.dyy(xx,yy);
                double t1 = hxx + hyy;
                double t2 = hxx*hyy - hxy*hxy;
                if(t2 > 0.0)
                {
                    if(t1 > 0.0)
                    {
                        minima->push_back(PointType(xx,yy));
                    }
                    else
                    {
                        maxima->push_back(PointType(xx,yy));
                    }
                }
                else
                {
                    saddles->push_back(PointType(xx,yy));
                }
            }
        }
        xold = x;
    }
}

template <class T, class VECTOR>
void
findCriticalPoints(
    SplineImageView<2, T> const & s,
    VECTOR *minima, VECTOR *saddles, VECTOR *maxima)
{
    for(unsigned int y=1; y<s.height()-1; ++y)
    {
        for(unsigned int x=1; x<s.width()-1; ++x)
        {
            findCriticalPointsInFacet(s, x, y, minima, saddles, maxima);
        }
    }
}

namespace detail {

enum RungeKuttaResult { Success, Outside, StepToLarge };

template <class T>
RungeKuttaResult
rungeKuttaInitial(SplineImageView<2, T> const & s,
                  double x0, double y0, double & h, bool forward,
                  double & x, double & y)
{
    double dxx = s.dxx(x0, y0);
    double dxy = s.dxy(x0, y0);
    double dyy = s.dyy(x0, y0);
    double a = 0.5*VIGRA_CSTD::atan2(-2.0*dxy, dxx-dyy);
    double dx = forward
                 ?  h*VIGRA_CSTD::cos(a)
                 : -h*VIGRA_CSTD::cos(a);
    double dy = forward
                 ? -h*VIGRA_CSTD::sin(a)
                 :  h*VIGRA_CSTD::sin(a);
    x = x0 + dx;
    y = y0 + dy;
    return s.isInside(x, y) ? Success : Outside;
}

template <class T>
RungeKuttaResult
rungeKuttaInitialStepSecondOrder(SplineImageView<2, T> const & s,
                  double x0, double y0, double h,
                  double & x, double & y, double dx, double dy)
{
    double x1, x2, y1, y2;
    x1 = x0 + 0.5*h*dx;
    y1 = y0 + 0.5*h*dy;
    if(!s.isInside(x1, y1))
        return Outside;
    x2 = x0 + h*s.dx(x1, y1);
    y2 = y0 + h*s.dy(x1, y1);
    if(!s.isInside(x2, y2))
        return Outside;
    x = x2;
    y = y2;
    return Success;
}

template <class T>
RungeKuttaResult
rungeKuttaStepSecondOrder(SplineImageView<2, T> const & s,
                  double x0, double y0, double h,
                  double & x, double & y)
{
    double x1, x2, y1, y2;
    x1 = x0 + 0.5*h*s.dx(x0, y0);
    y1 = y0 + 0.5*h*s.dy(x0, y0);
    if(!s.isInside(x1, y1))
        return Outside;
    x2 = x0 + h*s.dx(x1, y1);
    y2 = y0 + h*s.dy(x1, y1);
    if(!s.isInside(x2, y2))
        return Outside;
    x = x2;
    y = y2;
    return Success;
}

template <class T>
RungeKuttaResult
rungeKuttaDoubleStepSecondOrder(SplineImageView<2, T> const & s,
                  double x0, double y0, double & h,
                  double & x, double & y, double epsilon)
{
    double x1, x2, y1, y2;
    if(rungeKuttaStepSecondOrder(s, x0, y0, 2.0 * h, x1, y1) == Outside ||
       rungeKuttaStepSecondOrder(s, x0, y0, h, x2, y2) == Outside ||
       rungeKuttaStepSecondOrder(s, x2, y2, h, x2, y2) == Outside)
    {
        h /= 4.0;
        return Outside;
    }
    double dx = x2 - x1;
    double dy = y2 - y1;
    double d = std::max(std::abs(dx), std::abs(dy));
    double hh = d ? VIGRA_CSTD::pow(epsilon / d, 0.33) * h : 2*h;

    if(hh < h / 2.0)
    {
        h = hh;
        return StepToLarge;
    }
    x1 = x2 + dx / 3.0;
    y1 = y2 + dy / 3.0;
    if(s.isInside(x1, y1))
    {
        x = x1;
        y = y1;
    }
    else
    {
        x = x2;
        y = y2;
    }
    h = hh;
    return Success;
}


template <class IMAGEVIEW>
void rungeKuttaStep1(IMAGEVIEW const & image,
                    double x, double y, double s, double & xx, double & yy)
{
    xx = x;
    yy = y;
    if(!image.isInside(x, y))
        return;
    double x1 = s*image.dx(x,y);
    double y1 = s*image.dy(x,y);
    double nx = x + 0.5*x1;
    double ny = y + 0.5*y1;
    if(!image.isInside(nx, ny))
        return;
    double x2 = s*image.dx(nx, ny);
    double y2 = s*image.dy(nx, ny);
    nx = x + 0.5*x2;
    ny = y + 0.5*y2;
    if(!image.isInside(nx, ny))
        return;
    double x3 = s*image.dx(nx, ny);
    double y3 = s*image.dy(nx, ny);
    nx = x + x3;
    ny = y + y3;
    if(!image.isInside(nx, ny))
        return;
    double x4 = s*image.dx(nx, ny);
    double y4 = s*image.dy(nx, ny);
    xx = x + x1 / 6.0 + x2 / 3.0 + x3 / 3.0 + x4 / 6.0;
    yy = y + y1 / 6.0 + y2 / 3.0 + y3 / 3.0 + y4 / 6.0;
}

template <class IMAGEVIEW>
void rungeKuttaStep1b(IMAGEVIEW const & image,
                      double x, double y, double s, double & xx, double & yy)
{
    xx = x;
    yy = y;
    if(!image.isInside(x, y))
        return;
    double a = VIGRA_CSTD::atan2(-image.dy(x,y), image.dx(x,y));
    double dx = VIGRA_CSTD::cos(a);
    double dy = -VIGRA_CSTD::sin(a);
    double x1 = s*dx;
    double y1 = s*dy;
    double nx = x + 0.5*x1;
    double ny = y + 0.5*y1;
    if(!image.isInside(nx, ny))
        return;
    a = VIGRA_CSTD::atan2(-image.dy(nx,ny), image.dx(nx,ny));
    dx = VIGRA_CSTD::cos(a);
    dy = -VIGRA_CSTD::sin(a);
    double x2 = s*dx;
    double y2 = s*dy;
    nx = x + 0.5*x2;
    ny = y + 0.5*y2;
    if(!image.isInside(nx, ny))
        return;
    a = VIGRA_CSTD::atan2(-image.dy(nx,ny), image.dx(nx,ny));
    dx = VIGRA_CSTD::cos(a);
    dy = -VIGRA_CSTD::sin(a);
    double x3 = s*dx;
    double y3 = s*dy;
    nx = x + x3;
    ny = y + y3;
    if(!image.isInside(nx, ny))
        return;
    a = VIGRA_CSTD::atan2(-image.dy(nx,ny), image.dx(nx,ny));
    dx = VIGRA_CSTD::cos(a);
    dy = -VIGRA_CSTD::sin(a);
    double x4 = s*dx;
    double y4 = s*dy;
    xx = x + x1 / 6.0 + x2 / 3.0 + x3 / 3.0 + x4 / 6.0;
    yy = y + y1 / 6.0 + y2 / 3.0 + y3 / 3.0 + y4 / 6.0;
}

template <class IMAGEVIEW>
void rungeKuttaStep2(IMAGEVIEW const & image,
                    double x, double y, double s, double a, double & xx, double & yy)
{
    xx = x;
    yy = y;
    double dx = VIGRA_CSTD::cos(a);
    double dy = -VIGRA_CSTD::sin(a);
    if(!image.isInside(x, y))
        return;
    double a1 = 0.5*atan2(-2.0*image.dxy(x,y), (double)image.dxx(x,y)-image.dyy(x,y));
    double dx1 = VIGRA_CSTD::cos(a1);
    double dy1 = -VIGRA_CSTD::sin(a1);
    if(dx*dx1 + dy*dy1 < 0.0)
    {
        dx1 = -dx1;
        dy1 = -dy1;
    }
    double x1 = s*dx1;
    double y1 = s*dy1;
    double nx = x + 0.5*x1;
    double ny = y + 0.5*y1;
    if(!image.isInside(nx, ny))
        return;
    a1 = 0.5*atan2(-2.0*image.dxy(nx,ny), (double)image.dxx(nx,ny)-image.dyy(nx,ny));
    dx1 = VIGRA_CSTD::cos(a1);
    dy1 = -VIGRA_CSTD::sin(a1);
    if(dx*dx1 + dy*dy1 < 0.0)
    {
        dx1 = -dx1;
        dy1 = -dy1;
    }
    double x2 = s*dx1;
    double y2 = s*dy1;
    nx = x + 0.5*x2;
    ny = y + 0.5*y2;
    if(!image.isInside(nx, ny))
        return;
    a1 = 0.5*atan2(-2.0*image.dxy(nx,ny), (double)image.dxx(nx,ny)-image.dyy(nx,ny));
    dx1 = VIGRA_CSTD::cos(a1);
    dy1 = -VIGRA_CSTD::sin(a1);
    if(dx*dx1 + dy*dy1 < 0.0)
    {
        dx1 = -dx1;
        dy1 = -dy1;
    }
    double x3 = s*dx1;
    double y3 = s*dy1;
    nx = x + x3;
    ny = y + y3;
    if(!image.isInside(nx, ny))
        return;
    a1 = 0.5*atan2(-2.0*image.dxy(nx,ny), (double)image.dxx(nx,ny)-image.dyy(nx,ny));
    dx1 = VIGRA_CSTD::cos(a1);
    dy1 = -VIGRA_CSTD::sin(a1);
    if(dx*dx1 + dy*dy1 < 0.0)
    {
        dx1 = -dx1;
        dy1 = -dy1;
    }
    double x4 = s*dx1;
    double y4 = s*dy1;
    xx = x + x1 / 6.0 + x2 / 3.0 + x3 / 3.0 + x4 / 6.0;
    yy = y + y1 / 6.0 + y2 / 3.0 + y3 / 3.0 + y4 / 6.0;
}

template <class IMAGEVIEW, class PVECTOR, class CVECTOR>
void findEdgelChain1(IMAGEVIEW const & image, IImage const & maxImage,
                    PVECTOR const & maxima, CVECTOR & chain,
                    double x, double y, double epsilon)
{
    static int idx[] = {0, 1, 0, 1};
    static int idy[] = {0, 0, 1, 1};
    int w = image.width();
    int h = image.height();
    double d, dx, dy, x1, x2, x3, y1, y2, y3, mx, my, xx, yy, s, ss, a;
    double ox = x, oy = y;
    // add initial edgel
    dx = image.dx(x, y);
    dy = image.dy(x, y);
    a = VIGRA_CSTD::atan2(-dy, dx);
    chain.push_back(Edgel(x, y, image(x, y), a));

    // find initial step size
    d = VIGRA_CSTD::sqrt(dx*dx+dy*dy);
    s = 0.25 / d;

    while(1)
    {
        // do Runge/Kutta double step
        rungeKuttaStep1(image, x, y, 2.0 * s, x1, y1);
        rungeKuttaStep1(image, x, y, s, x2, y2);
        rungeKuttaStep1(image, x2, y2, s, x3, y3);
        dx = x3 - x1;
        dy = y3 - y1;
        xx = x3 + dx / 15.0;
        yy = y3 + dy / 15.0;
        d = VIGRA_CSTD::sqrt(dx*dx + dy*dy);
        if(d == 0.0)
            d = epsilon;
        ss = VIGRA_CSTD::pow(epsilon / d, 0.2) * s;
        if(ss >= s / 2.0) // accept step
        {
            // check if we landed outside the image
            if(!image.isInside(xx, yy))
            {
                break; // cannot further follow the edge
            }

            // add an edgel
            dx = image.dx(xx, yy);
            dy = image.dy(xx, yy);
            a = VIGRA_CSTD::atan2(-dy, dx);
            chain.push_back(Edgel(xx, yy, image(xx, yy), a));

            if(chain.size() > 1000)
            {
                std::cerr << "terminated chain: " << ox << ' ' << oy << '\n';
                break;
            }

            // check if we arrived at a maximum
            int ix =(int)xx;
            if(ix == w-1)
                --ix;
            int iy =(int)yy;
            if(iy == h-1)
                --iy;
            int i;
            for(i=0; i<4; ++i)
            {
                int index = maxImage(ix+idx[i], iy+idy[i]);
                if(index >= 0)
                {
                    mx = maxima[index][0];
                    my = maxima[index][1];
                    dx = mx - xx;
                    dy = my - yy;
                    d = VIGRA_CSTD::sqrt(dx*dx + dy*dy);
                    if(d < 0.25)
                        break;
                }
            }
            if(i < 4) // we are within 0.25 pixel of a maximum
            {
                a = VIGRA_CSTD::atan2(-dy, dx);
                chain.push_back(Edgel(mx, my, image(mx, my), a));
                break; // edgel chain is complete
            }
            x = xx;
            y = yy;
        }
        s = ss;
    }
}

template <class IMAGEVIEW, class PVECTOR, class CVECTOR>
void findEdgelChain2(IMAGEVIEW const & image, IImage const & maxImage,
                    PVECTOR const & maxima, CVECTOR & chain,
                    double x, double y, double a, double epsilon)
{
    static int idx[] = {0, 1, 0, 1};
    static int idy[] = {0, 0, 1, 1};
    int w = image.width();
    int h = image.height();
    double d, dx, dy, x1, x2, x3, y1, y2, y3, mx, my, xx, yy, s, ss, a1, d1, d2;
    double ox = x, oy = y;
    // add initial edgel
    chain.push_back(Edgel(x, y, image(x, y), a));

    // find initial step size
    s = 0.25;

    while(1)
    {
        // decide if we have to do a step based on first or second derivatives
        d = image.dxx(x,y) * image.dyy(x,y) - sq(image.dxy(x,y));
        d1 = sq(image.dx(x,y)) + sq(image.dy(x,y));
        d2 = sq(image.dxx(x,y) - image.dyy(x,y)) + sq(2.0*image.dxy(x,y));
        if(d1 < d2 / 100.0) // && d <= 0.0)
        {
            // do 2nd derivative step
            rungeKuttaStep2(image, x, y, 2.0 * s, a, x1, y1);
            rungeKuttaStep2(image, x, y, s, a, x2, y2);
            rungeKuttaStep2(image, x2, y2, s, a, x3, y3);
        }
        else
        {
            // do 1st derivative step
            rungeKuttaStep1b(image, x, y, 2.0 * s, x1, y1);
            rungeKuttaStep1b(image, x, y, s, x2, y2);
            rungeKuttaStep1b(image, x2, y2, s, x3, y3);
        }
        dx = x3 - x1;
        dy = y3 - y1;
        xx = x3 + dx / 15.0;
        yy = y3 + dy / 15.0;
        d = VIGRA_CSTD::sqrt(dx*dx + dy*dy);
        if(d < epsilon / 32.0)
            d = epsilon / 32.0;
        ss = VIGRA_CSTD::pow(epsilon / d, 0.2) * s;
        if(ss >= s / 2.0) // accept step
        {
            // check if we landed outside the image
            if(!image.isInside(xx, yy))
            {
                break; // cannot further follow the edge
            }

            // add an edgel
            if(d1 < d2)
            {
                a1 = 0.5*VIGRA_CSTD::atan2(-2.0*image.dxy(xx,yy), (double)image.dxx(xx,yy)-image.dyy(xx,yy));
                if(VIGRA_CSTD::cos(a)*VIGRA_CSTD::cos(a1)+VIGRA_CSTD::sin(a)*VIGRA_CSTD::sin(a1) < 0.0)
                    a1 += M_PI;
            }
            else
            {
                a1 = VIGRA_CSTD::atan2(-image.dy(xx,yy), image.dxx(xx,yy));
            }
            chain.push_back(Edgel(xx, yy, image(xx, yy), a1));

            // check if we're still moving
            if(sq(x-xx)+sq(y-yy) < epsilon*epsilon || chain.size() > 1000)
            {
                std::cerr << "terminated chain: " << ox << ' ' << oy << '\n';
                break;
            }

            // check if we arrived at a maximum
#if 0
            int ix =(int)xx;
            if(ix == w-1)
                --ix;
            int iy =(int)yy;
            if(iy == h-1)
                --iy;
            int i;
            for(i=0; i<4; ++i)
            {
                int index = maxImage(ix+idx[i], iy+idy[i]);
                if(index >= 0)
                {
                    mx = maxima[index][0];
                    my = maxima[index][1];
                    dx = mx - xx;
                    dy = my - yy;
                    d = VIGRA_CSTD::sqrt(dx*dx + dy*dy);
                    if(d < 0.25)
                        break;
                }
            }
            if(i < 4) // we are within 0.25 pixel of a maximum
            {
                a = VIGRA_CSTD::atan2(-dy, dx);
                chain.push_back(Edgel(mx, my, image(mx, my), a));
                break; // edgel chain is complete
            }
#endif
            x = xx;
            y = yy;
            a = a1;
            CriticalPoint type = findCriticalPoint(image, x, y, xx, yy, epsilon);
            if(type != Failed && VIGRA_CSTD::sqrt(sq(xx-x)+sq(yy-y)) < 0.25)
            {
                 // we are within 0.25 pixel of a maximum
                a = VIGRA_CSTD::atan2(y-yy, xx-x);
                chain.push_back(Edgel(xx, yy, image(xx, yy), a));
                break; // edgel chain is complete
            }
        }
        s = ss;
    }
}

} // namespace detail

template <class T, class VECTOR>
void
flowLine(SplineImageView<2, T> const & s,
         VECTOR & c, bool forward, double epsilon)
{
    typedef typename VECTOR::value_type PointType;
    double x0 = c[0][0];
    double y0 = c[0][1];
    double h = 0.25;
    double x, y;

    detail::RungeKuttaResult r = detail::rungeKuttaInitial(s, x0, y0, h, forward, x, y);
    if(r == detail::Outside)
        return;
    c.push_back(PointType(x, y));
    ArrayVector<PointType> mi, sa, ma;
    findCriticalPointsInFacet(s, x, y, &mi, &sa, &ma);

    for(int k=0; k<1000; ++k)
    {
        double xn, yn;
        detail::RungeKuttaResult r =
                    detail::rungeKuttaDoubleStepSecondOrder(s, x, y, h, xn, yn, epsilon);
        if(r == detail::Success)
        {
            c.push_back(PointType(xn, yn));
            // check if near a maximum
            if(!s.sameFacet(x, y, xn, yn))
            {
                mi.clear();
                sa.clear();
                ma.clear();
                findCriticalPointsInFacet(s, xn, yn, &mi, &sa, &ma);
            }
            unsigned int i = 0;
            for(; i<ma.size(); ++i)
            {
                if(std::abs(xn-ma[i][0]) < 0.2 && std::abs(yn-ma[i][1]) < 0.2)
                {
                    c.push_back(ma[0]);
                    break;
                }
            }
            if(i < ma.size())
                break;
            x = xn;
            y = yn;
        }
        else if(h < 1.0e-6)
        {
            break;
        }
    }
}

template <class IMAGEVIEW, class PVECTOR, class CVECTOR>
void findEdgelChains(IMAGEVIEW const & image,
                     PVECTOR const & maxima, PVECTOR const & saddles,
                     CVECTOR & chains, double epsilon)
{
    typedef typename IMAGEVIEW::value_type Value;
    typedef typename CVECTOR::value_type EdgelChain;

    IImage maxImage(image.size());
    maxImage.init(-1);
    for(unsigned int i=0; i<maxima.size(); ++i)
    {
        int x = (int)(maxima[i][0] + 0.5);
        int y = (int)(maxima[i][1] + 0.5);
        maxImage(x,y) = i;
    }

    for(unsigned int i=0; i<saddles.size(); ++i)
    {
        double x = saddles[i][0];
        double y = saddles[i][1];
        double dxx = image.dxx(x, y);
        double dxy = image.dxy(x, y);
        double dyy = image.dyy(x, y);
        double a = 0.5*VIGRA_CSTD::atan2(-2.0*dxy, dxx-dyy);

        EdgelChain forwardChain;
        detail::findEdgelChain2(image, maxImage, maxima, forwardChain,
                       x+0.25*VIGRA_CSTD::cos(a), y-0.25*VIGRA_CSTD::sin(a), a, epsilon);
        EdgelChain backwardChain;
        detail::findEdgelChain2(image, maxImage, maxima, backwardChain,
                       x-0.25*VIGRA_CSTD::cos(a), y+0.25*VIGRA_CSTD::sin(a), a + M_PI, epsilon);

        // join chains
        chains.push_back(EdgelChain());
        EdgelChain & chain = *(chains.end()-1);
        for(int i = backwardChain.size()-1; i >= 0; --i)
        {
            chain.push_back(backwardChain[i]);
        }
        chain.push_back(Edgel(x, y, image(x, y), a));
        for(unsigned int i = 0; i < forwardChain.size(); ++i)
        {
            chain.push_back(forwardChain[i]);
        }
        std::cerr << '*';
    }
    std::cerr << '\n';
}

template <class Poly>
void polynomialRemainder(Poly const & l, Poly const & r, Poly & res)
{
    typedef typename Poly::value_type Coeff;
    res = l;
    int diffOrder = res.order() - r.order();
    while(diffOrder >= 0)
    {
        Coeff c = res[res.order()] / r[r.order()];
        for(unsigned int k = 0; k < r.order(); ++k)
        {
            res[k+diffOrder] -= c*r[k];
        }
        res.setOrder(res.order()-1);
        res.minimizeOrder(0.0);
        diffOrder = res.order() - r.order();
    }
}

/*******************************************************************/
/*                                                                 */
/*                         SubPixelWatersheds                      */
/*                                                                 */
/*******************************************************************/

#define DEBUG 0

template <class SPLINEIMAGEVIEW>
class SubPixelWatersheds
{
  public:
    typedef SPLINEIMAGEVIEW SplineImageView;

    typedef TinyVector<double, 2> PointType;
    typedef ArrayVector<PointType> PointArray;
    enum RungeKuttaResult { Success, Outside, StepToLarge };

    template <class SrcIterator, class SrcAccessor>
    SubPixelWatersheds(SrcIterator ul, SrcIterator lr, SrcAccessor src)
    : image_(ul, lr, src),
      initialStep_(0.1),
      minCPDist_(1e-3),
      stepEpsilon_(1e-4),
      cpOversampling_(2)
    {}

    template <class SrcIterator, class SrcAccessor>
    SubPixelWatersheds(triple<SrcIterator, SrcIterator, SrcAccessor> src)
    : image_(src),
      initialStep_(0.1),
      minCPDist_(1e-3),
      stepEpsilon_(1e-4),
      cpOversampling_(2)
    {}

    int width() const { return image_.width(); }
    int height() const { return image_.height(); }

    void findCriticalPointsInFacet(double x, double y,
                                   PointArray & mi, PointArray & sa, PointArray & ma);
    template <class MaskIterator, class MaskAccessor>
    void findCriticalPoints(pair<MaskIterator, MaskAccessor> mask);
    void findCriticalPoints();
    void updateMaxImage();
    double nearestMaximum(double x, double y, double dx, double dy, int & resindex) const;
    int flowLine(double x, double y, bool forward, double epsilon, PointArray & curve);
    pair<int, int> findEdge(double x, double y, double epsilon, PointArray & edge);
    RungeKuttaResult rungeKuttaStepSecondOrder(double x0, double y0, double h,
                                               double & x, double & y, double dx, double dy);
    RungeKuttaResult rungeKuttaDoubleStepSecondOrder(double x0, double y0, double & h,
                            double & x, double & y, double epsilon, double dx, double dy);

    SplineImageView image_;
    PointArray minima_, saddles_, maxima_;
    IImage maxImage_;
    double initialStep_, minCPDist_, stepEpsilon_;
    unsigned int cpOversampling_;
};

// static int sturmcount, zeroOrder;

template <class SplineImageView>
void
SubPixelWatersheds<SplineImageView>::findCriticalPointsInFacet(double x0, double y0,
                                   PointArray & mi, PointArray & sa, PointArray & ma)
{
#if 0
    x0 = VIGRA_CSTD::floor(x0 + 0.5);
    y0 = VIGRA_CSTD::floor(y0 + 0.5);

    DImage splineCoeffs(3,3);
    image_.coefficientArray(x0, y0, splineCoeffs);
    double j = splineCoeffs(0,0);
    double g = splineCoeffs(1,0);
    double e = splineCoeffs(2,0);
    double h = splineCoeffs(0,1);
    double d = splineCoeffs(1,1);
    double b = splineCoeffs(2,1);
    double f = splineCoeffs(0,2);
    double c = splineCoeffs(1,2);
    double a = splineCoeffs(2,2);

    double eps = 1.0e-7;
    StaticPolynomial<5, double> polys[6];
    polys[5].setOrder(5);
    polys[5].setEpsilon(eps);
    polys[5][0] =  4.0*f*f*g - 2.0*d*f*h + c*h*h;
    polys[5][1] = -2.0*d*d*f + 8.0*e*f*f + 8.0*c*f*g - 4.0*b*f*h + 2.0*a*h*h;
    polys[5][2] = -c*d*d - 6.0*b*d*f + 16.0*c*e*f + 4.0*c*c*g +
                   8.0*a*f*g - 2.0*b*c*h + 2.0*a*d*h;
    polys[5][3] = -4.0*b*c*d + 8.0*c*c*e - 4.0*b*b*f + 16.0*a*e*f + 8.0*a*c*g;
    polys[5][4] = -3.0*b*b*c - 2.0*a*b*d + 16.0*a*c*e + 4.0*a*a*g;
    polys[5][5] = -2.0*a*b*b + 8.0*a*a*e;

    polys[5].minimizeOrder();
    if(polys[5].order() == 0)
    {
        //zeroOrder++;
        std::cerr << x0 << "/" << y0 << " resulted in zero-order poly!\n";
        return;
    }
    // create the Sturm sequence and count the sign changes
    double left = polys[5](-0.5);
    double right = polys[5](0.5);

    // if the interval bounds are zero, we go directly to root finding
    if(0)//std::abs(left) > eps && std::abs(right) > eps)
    {
        polys[4] = polys[5].getDerivative();
        int m;
        for(m = 3; m >= 0; --m)
        {
            polynomialRemainder(polys[m+2], polys[m+1], polys[m]);
            if(polys[m].order() == 0)
                break;
        }
        int leftCount = 0;
        int rightCount = 0;
        double v;
        for(int k = 4; k >= m; --k)
        {
            v = polys[k](-0.5);
            if(v*left < 0.0)
                leftCount++;
            if(v != 0.0)
                left = v;
            v = polys[k](0.5);
            if(v*right < 0.0)
                rightCount++;
            if(v != 0.0)
                right = v;
        }
        if(leftCount == rightCount)
        {
            //++sturmcount;
            return; // interval [-0.5, 0.5] cannot contain a zero
        }
    }

    ArrayVector<double> rx;
    rx.reserve(5);
//     if((x0 == 19) && (y0 == 12))
//     {
//         std::cerr << "Poly order " << polys[5].order() << " coeffs ";
//         for(int k = 0; k<=polys[5].order(); ++k)
//             std::cerr << polys[5][k] << ' ';
//         std::cerr << '\n';
//         polynomialRealRoots(polys[5], rx, true); // no root polishing necessary ?
//     }
//     else
    polynomialRealRoots(polys[5], rx, true); // no root polishing necessary ?
//     if((x0 == 19) && (y0 == 12))
//     {
//         std::cerr << "Solution size " << rx.size() << " coeffs ";
//         for(int k = 0; k<rx.size(); ++k)
//             std::cerr << rx[k] << ' ';
//         std::cerr << '\n';
//     }

    double xold = -100.0;
    for(unsigned int i=0; i < rx.size(); ++i)
    {
        double x = rx[i];

        // ensure that x is in the current facet,
        // and that a multilple root is only used once (this may be
        // wrong, as perhaps several y' share the same x)
        if(std::abs(x) <= 0.5 && std::abs(x-xold) >= eps)
        {
            double cy = 2*f + 2*c*x + 2*a*x*x;
            if(cy == 0.0)
                continue;
            double y = -(h + d*x + b*x*x) / cy;
            if(std::abs(y) <= 0.5)
            {
                double xx = x + x0;
                double yy = y + y0;
                if(image_.g2(xx,yy) > 1e-4)
                    continue;  // ??? this shouldn't happen
                double hxx = image_.dxx(xx,yy);
                double hxy = image_.dxy(xx,yy);
                double hyy = image_.dyy(xx,yy);
                double t1 = hxx + hyy;
                double t2 = hxx*hyy - hxy*hxy;
                if(t2 > 0.0)
                {
                    if(t1 > 0.0)
                    {
                       mi.push_back(PointType(xx,yy));
                    }
                    else
                    {
                        ma.push_back(PointType(xx,yy));
                    }
                }
                else
                {
                    sa.push_back(PointType(xx,yy));
                }
            }
        }
        xold = x;
    }
#endif
}

template <class SplineImageView>
template<class MaskIterator, class MaskAccessor>
void SubPixelWatersheds<SplineImageView>::findCriticalPoints(
    pair<MaskIterator, MaskAccessor> mask)
{
    minima_.clear();
    saddles_.clear();
    maxima_.clear();
    // use 1-based arrays
    minima_.push_back(PointType());
    saddles_.push_back(PointType());
    maxima_.push_back(PointType());

    findCriticalPointsNewtonMethodIf(
        image_, mask, &minima_, &saddles_, &maxima_, minCPDist_, stepEpsilon_, cpOversampling_);
    updateMaxImage();
}

template <class SplineImageView>
void
SubPixelWatersheds<SplineImageView>::findCriticalPoints()
{
    minima_.clear();
    saddles_.clear();
    maxima_.clear();
    // use 1-based arrays
    minima_.push_back(PointType());
    saddles_.push_back(PointType());
    maxima_.push_back(PointType());

//     sturmcount = 0;
//     zeroOrder = 0;
//     for(unsigned int y=1; y<image_.height()-1; ++y)
//     {
//         std::cerr << ".";
//         for(unsigned int x=1; x<image_.width()-1; ++x)
//         {
//             findCriticalPointsInFacet(x, y, minima_, saddles_, maxima_);
//         }
//     }
    findCriticalPointsNewtonMethod(
        image_, &minima_, &saddles_, &maxima_, minCPDist_, stepEpsilon_, cpOversampling_);
    updateMaxImage();
//     std::cerr << "Sturm fired: " << sturmcount << " times\n";
//     std::cerr << "Zero order fired: " << zeroOrder << " times\n";
}

struct PointSort
{
    template <class Point>
    bool operator()(Point const & l, Point const & r) const
    {
        double yl = VIGRA_CSTD::floor(l[1] + 0.5);
        double yr = VIGRA_CSTD::floor(r[1] + 0.5);
        return (yl < yr) || (yl == yr && VIGRA_CSTD::floor(l[0] + 0.5) < VIGRA_CSTD::floor(r[0] + 0.5));
    }
};

template <class SplineImageView>
void
SubPixelWatersheds<SplineImageView>::updateMaxImage()
{
    maxImage_.resize(image_.size());
    maxImage_.init(0);  // means no maximum
    std::sort(maxima_.begin()+1, maxima_.end(), PointSort()); // necessary so that points in the same facet are consecutive
    for(unsigned int i = 1; i<maxima_.size(); ++i)
    {
        int x = (int)VIGRA_CSTD::floor(maxima_[i][0] + 0.5);
        int y = (int)VIGRA_CSTD::floor(maxima_[i][1] + 0.5);

        if(maxImage_(x,y) == 0)
        {
            maxImage_(x,y) = i;  // store index of the subpixel coordinate object
        }
        else
        {
            // negative index denotes several maxima in the same facet
            maxImage_(x,y) = -std::abs(maxImage_(x,y));
        }
    }
}

template <class SplineImageView>
double
SubPixelWatersheds<SplineImageView>::nearestMaximum(double x, double y, double dx, double dy,
                                      int & resindex) const
{
    PointType diff, p(x,y);
    double dist = NumericTraits<double>::max();

    int xi = (int)VIGRA_CSTD::floor(x + 0.5);
    int yi = (int)VIGRA_CSTD::floor(y + 0.5);

    int dxi = (x - xi >= 0.0)
                 ?  1
                 : -1;
    int dyi = (y - yi >= 0.0)
                 ?  1
                 : -1;

    for(int ys=0; ys<2; ++ys)
    {
        for(int xs=0; xs<2; ++xs)
        {
            int xx = xi + xs*dxi;
            int yy = yi + ys*dyi;
            if(!image_.isInside(xx,yy))
                continue;
            int index = maxImage_(xx,yy);
            if(index > 0)
            {
                diff = maxima_[index] - p;
                double sd = diff.magnitude();
                if(sd < dist && diff[0]*dx + diff[1]*dy >= 0)
                {
                    dist = sd;
                    resindex = index;
                }
            }
            else if(index < 0)
            {
                diff = maxima_[-index] - p;
                double sd = diff.magnitude();
                if(sd < dist && diff[0]*dx + diff[1]*dy >= 0)
                {
                    dist = sd;
                    resindex = -index;
                }
                diff = maxima_[-index+1] - p; // assumes that points from the same facet are consecutive
                sd = diff.magnitude();
                if(sd < dist && diff[0]*dx + diff[1]*dy >= 0)
                {
                    dist = sd;
                    resindex = -index+1;
                }
            }
        }
    }

    return dist;
}

template <class SplineImageView>
typename SubPixelWatersheds<SplineImageView>::RungeKuttaResult
SubPixelWatersheds<SplineImageView>::rungeKuttaStepSecondOrder(
                  double x0, double y0, double h,
                  double & x, double & y, double dx, double dy)
{
    double x1, x2, y1, y2;
    x1 = x0 + 0.5*h*dx;
    y1 = y0 + 0.5*h*dy;
    if(!image_.isInside(x1, y1))
        return Outside;
    dx = image_.dx(x1, y1);
    dy = image_.dy(x1, y1);
    double norm = hypot(dx, dy);
    if(!norm) // critical point?
    {
        x = x1;
        y = y1;
        return Success;
    }
    x2 = x0 + h*dx/norm;
    y2 = y0 + h*dy/norm;
    if(!image_.isInside(x2, y2))
        return Outside;
    x = x2;
    y = y2;
    return Success;
}

template <class SplineImageView>
typename SubPixelWatersheds<SplineImageView>::RungeKuttaResult
SubPixelWatersheds<SplineImageView>::rungeKuttaDoubleStepSecondOrder(
                  double x0, double y0, double & h,
                  double & x, double & y, double epsilon, double dx, double dy)
{
    double x1, x2, y1, y2;
    if(rungeKuttaStepSecondOrder(x0, y0, 2.0 * h, x1, y1, dx, dy) == Outside ||
       rungeKuttaStepSecondOrder(x0, y0, h, x2, y2, dx, dy) == Outside)
    {
        h /= 4.0;
        return Outside;
    }
    dx = image_.dx(x2, y2);
    dy = image_.dy(x2, y2);
    double norm = hypot(dx, dy);
    if(!norm) // critical point
    {
        x = x2;
        y = y2;
        return Success;
    }
    dx /= norm;
    dy /= norm;
    if(rungeKuttaStepSecondOrder(x2, y2, h, x2, y2, dx, dy) == Outside)
    {
        h /= 4.0;
        return Outside;
    }

    // check that we don't jump too far (next step shall not go backwards)
    if((x2-x0)*image_.dx(x2, y2) + (y2-y0)*image_.dy(x2, y2) <= 0.0)
    {
        h /= 2.0;
        return StepToLarge;
    }

    // estimate error and desirable step size (hh)
    dx = x2 - x1;
    dy = y2 - y1;
    double d = std::max(std::abs(dx), std::abs(dy));
    double hh = VIGRA_CSTD::pow(epsilon / d, 0.33) * h;

    if(hh < h / 2.0)
    {
        h = hh;
        return StepToLarge;
    }

    x = x2;
    y = y2;
    h = hh;
    return Success;
}

template <class SplineImageView>
int
SubPixelWatersheds<SplineImageView>::flowLine(double x, double y, bool forward, double epsilon,
                                PointArray & curve)
{
    curve.push_back(PointType(x, y));
    double h = initialStep_;
    double dxx = image_.dxx(x, y);
    double dxy = image_.dxy(x, y);
    double dyy = image_.dyy(x, y);
    double a = 0.5*VIGRA_CSTD::atan2(-2.0*dxy, dxx-dyy);
    double dx = forward
                 ?  h*VIGRA_CSTD::cos(a)
                 : -h*VIGRA_CSTD::cos(a);
    double dy = forward
                 ? -h*VIGRA_CSTD::sin(a)
                 :  h*VIGRA_CSTD::sin(a);
    int failReason = 0; // gave up / too many steps
    if(DEBUG) std::cerr << "x, y, dx, dy " << x << ' ' << y << ' ' << dx << ' ' << dy << '\n';
    int index;
    if(nearestMaximum(x, y, dx, dy, index) < initialStep_)
    {
        curve.push_back(maxima_[index]);
        if(DEBUG) std::cerr << "stop index, x, y " << index << ' ' << curve.back()[0] << ' ' << curve.back()[1]<< '\n';
        return index;
    }

    x = x + dx;
    y = y + dy;

    curve.push_back(PointType(x, y));
    dx = image_.dx(x, y);
    dy = image_.dy(x, y);
    double norm = hypot(dx, dy);
    if(norm) // critical point?
    {
        dx /= norm;
        dy /= norm;
    }
if(DEBUG) std::cerr << "x, y, dx, dy " << x << ' ' << y << ' ' << dx << ' ' << dy << '\n';
    // check if near a maximum in direction dx/dy
    if(nearestMaximum(x, y, dx, dy, index) < initialStep_)
    {
        curve.push_back(maxima_[index]);
if(DEBUG) std::cerr << "stop index, x, y " << index << ' ' << curve.back()[0] << ' ' << curve.back()[1]<< '\n';
        return index;
    }

    for(int k = 0; k < 100000; ++k)
    {
        double xn, yn;
        RungeKuttaResult rungeKuttaResult(
            rungeKuttaDoubleStepSecondOrder(x, y, h, xn, yn, epsilon, dx, dy));
        if(rungeKuttaResult == Success)
        {
            x = xn;
            y = yn;
            curve.push_back(PointType(x, y));
            dx = image_.dx(x, y);
            dy = image_.dy(x, y);
            double norm = hypot(dx, dy);
            if(norm) // critical point?
            {
                dx /= norm;
                dy /= norm;
            }
            if(DEBUG) std::cerr << "x, y, dx, dy " << x << ' ' << y
                                << ' ' << dx << ' ' << dy << '\n';
            // check if near a maximum
            if(nearestMaximum(x, y, dx, dy, index) < initialStep_)
            {
                curve.push_back(maxima_[index]);
                if(DEBUG) std::cerr
                    << "stop index, x, y " << index << ' '
                    << curve.back()[0] << ' ' << curve.back()[1]<< '\n';
                return index;
            }
        }
        else if(rungeKuttaResult == Outside)
        {
            if(DEBUG) std::cerr << "outside image\n";
            failReason = -1; // signal "outside image", curve might still be useful?!
        }
        if(h < 1.0e-6)
        {
            if(DEBUG) std::cerr << "give up\n";
            return failReason; // give up
        }
    }
    std::cerr << "flowLine(): too many steps, h = " << h << "\n";
    return failReason-3;
}

template <class SplineImageView>
pair<int, int>
SubPixelWatersheds<SplineImageView>::findEdge(
    double x, double y, double epsilon, PointArray & edge)
{
    PointArray forwardCurve, backwardCurve;

        if(DEBUG) std::cerr << "starting forward\n";
    int findex = flowLine(x, y, true, epsilon, forwardCurve);
        if(DEBUG) std::cerr << "starting backward\n";
    int bindex = flowLine(x, y, false, epsilon, backwardCurve);

    for(int i = forwardCurve.size() - 1; i >= 0; --i)
        edge.push_back(forwardCurve[i]);
    for(int i = 1; i < (int)backwardCurve.size(); ++i)
        edge.push_back(backwardCurve[i]);

    return pair<int, int>(findex, bindex);
}

} // namespace vigra

#endif // VIGRA_SUBPIXEL_WATERSHED_HXX
