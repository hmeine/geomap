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
#include "vigra/mathutil.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/polynomial.hxx"
#include "vigra/splineimageview.hxx"

namespace vigra {

enum CriticalPoint { Minimum = -1, Saddle, Maximum, Failed = 999 };

template <class IMAGEVIEW>
CriticalPoint
findCriticalPointNewtonMethod(IMAGEVIEW const & image, 
                  double x, double y, double & xx, double & yy,
                  double epsilon)
{
    typedef typename IMAGEVIEW::value_type Value;
    Value zero = NumericTraits<Value>::zero();

    xx = x;
    yy = y;
    double sxx, syy;
    for(int i=0; i<10; ++i) // do at most 10 iterations
    {
        Value dx = image.dx(xx, yy);
        Value dy = image.dy(xx, yy);
        Value dxx = image.dxx(xx, yy);
        Value dxy = image.dxy(xx, yy);
        Value dyy = image.dyy(xx, yy);
        Value d = dxx*dyy - dxy*dxy;
        if (d != zero)
        {
            sxx = (dxy*dy - dyy*dx) / d;
            syy = (dxy*dx - dxx*dy) / d;
        }
        else
        {
            sxx = syy = 0.0;
        }
        xx += sxx;
        yy += syy;
        if(!image.isInside(xx, yy))
        {
            return Failed; // coordinates out of range
        }
        double diff = VIGRA_CSTD::sqrt(sxx*sxx + syy*syy);
        if(diff < epsilon) // convergence
        {
            if(d == zero)
            {
                if(dx == zero && dy == zero)
                {
                    return Saddle;
                }
                return Failed;
                
            }
            else if (d < zero)
            {
                return Saddle;
            }
            else if (dxx + dyy > zero)
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

template <class IMAGEVIEW, class VECTOR>
void findCriticalPointsNewtonMethod(IMAGEVIEW const & image, 
                        VECTOR & minima, VECTOR & maxima, VECTOR & saddles,
                        double epsilon)
{
    int w = image.width();
    int h = image.height();
    
    typedef typename IMAGEVIEW::value_type Value;
    Value zero = NumericTraits<Value>::zero();
    typedef typename VECTOR::value_type Coordinate;
    
    for(int y=1; y<h-1; ++y)
    {
        for(int x=1; x<w-1; ++x)
        {
            double xx, yy;
            CriticalPoint type = findCriticalPointNewtonMethod(image, x, y, xx, yy, epsilon);
            if(type == Failed)
                continue;
            if(VIGRA_CSTD::abs(xx - x) > 0.5 || VIGRA_CSTD::abs(yy - y) > 0.5)
            {
                continue; // out of current pixel region
            }
            if(type == Saddle)
            {
                saddles.push_back(Coordinate(xx, yy));
            }
            else if (type == Minimum)
            {
                minima.push_back(Coordinate(xx, yy));
            }
            else
            {
                maxima.push_back(Coordinate(xx, yy));
            }
        }
    }
}

template <class T, class VECTOR>
void 
findCriticalPointsInFacet(
    SplineImageView<2, T> const & s, double x0, double y0,
    VECTOR & minima, VECTOR & saddles, VECTOR & maxima)
{
    typedef typename VECTOR::value_type PointType;
    
    x0 = VIGRA_CSTD::floor(x0 + 0.5);
    y0 = VIGRA_CSTD::floor(y0 + 0.5);

    DImage splineCoeffs(3,3);
    s.coefficientArray(x0, y0, splineCoeffs);
    double j = splineCoeffs(0,0);
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
    Polynomial<double> px(polyCoeffs, 6 , eps);
    px.minimizeOrder();
    std::vector<double> rx;
    polynomialRealRoots(px, rx);
    
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
                double hxx = s.dxx(xx,yy);
                double hxy = s.dxy(xx,yy);
                double hyy = s.dyy(xx,yy);
                double t1 = hxx + hyy;
                double t2 = hxx*hyy - hxy*hxy;
                if(t2 > 0.0)
                {
                    if (t1 > 0.0)
                    {
                        minima.push_back(PointType(xx,yy));
                    }
                    else
                    {
                        maxima.push_back(PointType(xx,yy));
                    }
                }
                else
                {
                    saddles.push_back(PointType(xx,yy));
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
    VECTOR & minima, VECTOR & saddles, VECTOR & maxima)
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
    double hh = VIGRA_CSTD::pow(epsilon / d, 0.33) * h;
    
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
void findEdgelChain1(IMAGEVIEW const & image, IImage const & maximage,
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
        if (ss >= s / 2.0) // accept step
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
                int index = maximage(ix+idx[i], iy+idy[i]);
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
void findEdgelChain2(IMAGEVIEW const & image, IImage const & maximage,
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
        if (ss >= s / 2.0) // accept step
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
                int index = maximage(ix+idx[i], iy+idy[i]);
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
    findCriticalPointsInFacet(s, x, y, mi, sa, ma);

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
                findCriticalPointsInFacet(s, xn, yn, mi, sa, ma);
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

    IImage maximage(image.size());
    maximage.init(-1);
    for(unsigned int i=0; i<maxima.size(); ++i)
    {
        int x = (int)(maxima[i][0] + 0.5);
        int y = (int)(maxima[i][1] + 0.5);
        maximage(x,y) = i;
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
        detail::findEdgelChain2(image, maximage, maxima, forwardChain, 
                       x+0.25*VIGRA_CSTD::cos(a), y-0.25*VIGRA_CSTD::sin(a), a, epsilon);
        EdgelChain backwardChain;
        detail::findEdgelChain2(image, maximage, maxima, backwardChain, 
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

/*******************************************************************/
/*                                                                 */
/*                         SubPixelWatersheds                      */
/*                                                                 */
/*******************************************************************/

#define DEBUG 0

template <class T>
class SubPixelWatersheds
{
  public:
    typedef TinyVector<double, 2> PointType;
    typedef ArrayVector<PointType> PointArray;
    enum RungeKuttaResult { Success, Outside, StepToLarge };
    
    template <class SrcIterator, class SrcAccessor>
    SubPixelWatersheds(SrcIterator ul, SrcIterator lr, SrcAccessor src)
    : image_(ul, lr, src),
      initialStep_(0.1), storeStep_(0.5)
    {}
    
    template <class SrcIterator, class SrcAccessor>
    SubPixelWatersheds(triple<SrcIterator, SrcIterator, SrcAccessor> src)
    : image_(src),
      initialStep_(0.1), storeStep_(0.5)
    {}
    
    int width() const { return image_.width(); }
    int height() const { return image_.height(); }
    
    void findCriticalPointsInFacet(double x, double y, 
                                   PointArray & mi, PointArray & sa, PointArray & ma);
    void findCriticalPoints();
    void createMaximage();
    double nearestMaximum(double x, double y, double dx, double dy, int & resindex) const;
    int flowLine(double x, double y, bool forward, double epsilon, PointArray & curve);
    pair<int, int> findEdge(double x, double y, double epsilon, PointArray & edge);
    int findEdges(double threshold, double epsilon = 1e-4);
    RungeKuttaResult rungeKuttaStepSecondOrder(double x0, double y0, double h,
                                               double & x, double & y, double dx, double dy);
    RungeKuttaResult rungeKuttaDoubleStepSecondOrder(double x0, double y0, double & h,
                            double & x, double & y, double epsilon, double dx, double dy);

    RungeKuttaResult rungeKuttaStepSecondOrder(double x0, double y0, double h,
                                               double & x, double & y)
    {
        return rungeKuttaStepSecondOrder(x0, y0, h, x, y, 
                                         image_.dx(x0, y0), image_.dy(x0, y0));
    }

    RungeKuttaResult rungeKuttaDoubleStepSecondOrder(double x0, double y0, double & h,
                                                     double & x, double & y, double epsilon)
    {
        return rungeKuttaDoubleStepSecondOrder(x0, y0, h, x, y, epsilon, 
                                            image_.dx(x0,y0), image_.dy(x0,y0));
    }
    
    SplineImageView<2, T> image_;
    PointArray minima_, saddles_, maxima_, edges_;
    ArrayVector<triple<int, int, int> > edgeIndices_;
    IImage maximage_;
    double initialStep_, storeStep_;
};

template <class T>
void 
SubPixelWatersheds<T>::findCriticalPointsInFacet(double x0, double y0, 
                                   PointArray & mi, PointArray & sa, PointArray & ma)
{
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

    double polyCoeffs[6];
    polyCoeffs[0] =  4.0*f*f*g - 2.0*d*f*h + c*h*h;
    polyCoeffs[1] = -2.0*d*d*f + 8.0*e*f*f + 8.0*c*f*g - 4.0*b*f*h + 2.0*a*h*h;
    polyCoeffs[2] = -c*d*d - 6.0*b*d*f + 16.0*c*e*f + 4.0*c*c*g + 
                     8.0*a*f*g - 2.0*b*c*h + 2.0*a*d*h;
    polyCoeffs[3] = -4.0*b*c*d + 8.0*c*c*e - 4.0*b*b*f + 16.0*a*e*f + 8.0*a*c*g;
    polyCoeffs[4] = -3.0*b*b*c - 2.0*a*b*d + 16.0*a*c*e + 4.0*a*a*g;
    polyCoeffs[5] = -2.0*a*b*b + 8.0*a*a*e;
    
    double eps = 1.0e-7;
    Polynomial<double> px(polyCoeffs, 6 , eps);
    ArrayVector<double> rx;
#if 0
    std::cerr << "Poly order " << px.order() << " coeffs ";
    for(int k = 0; k<px.order(); ++k)
        std::cerr << px[k] << ' ';
    std::cerr << '\n';
#endif
    polynomialRealRoots(px, rx);
    
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
                    if (t1 > 0.0)
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
}

template <class T>
void 
SubPixelWatersheds<T>::findCriticalPoints()
{
    minima_.clear();
    saddles_.clear();
    maxima_.clear();
    // use 1-based arrays
    minima_.push_back(PointType());
    saddles_.push_back(PointType());
    maxima_.push_back(PointType());
    
    for(unsigned int y=1; y<image_.height()-1; ++y)
    {
        for(unsigned int x=1; x<image_.width()-1; ++x)
        {
            findCriticalPointsInFacet(x, y, minima_, saddles_, maxima_);
        }
    }
    createMaximage();
}

template <class T>
void 
SubPixelWatersheds<T>::createMaximage()
{
    maximage_.resize(image_.size());
    maximage_.init(0);  // means no maximum
    for(unsigned int i = 1; i<maxima_.size(); ++i)
    {
        int x = VIGRA_CSTD::floor(maxima_[i][0] + 0.5);
        int y = VIGRA_CSTD::floor(maxima_[i][1] + 0.5);
        
        if(maximage_(x,y) == 0)
        {
            maximage_(x,y) = i;  // store index of the subpixel coordinate object
        }
        else
        {
            // negative index denotes several maxima in the same facet
            maximage_(x,y) = -std::abs(maximage_(x,y));
        }
    }
}

template <class T>
double
SubPixelWatersheds<T>::nearestMaximum(double x, double y, double dx, double dy,
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
            int index = maximage_(xx,yy);
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
                diff = maxima_[-index+1] - p;
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

template <class T>
typename SubPixelWatersheds<T>::RungeKuttaResult
SubPixelWatersheds<T>::rungeKuttaStepSecondOrder(
                  double x0, double y0, double h,
                  double & x, double & y, double dx, double dy)
{
    double x1, x2, y1, y2;
    x1 = x0 + 0.5*h*dx;
    y1 = y0 + 0.5*h*dy;
    if(!image_.isInside(x1, y1))
        return Outside;
    x2 = x0 + h*image_.dx(x1, y1);
    y2 = y0 + h*image_.dy(x1, y1);
    if(!image_.isInside(x2, y2))
        return Outside;
    x = x2;
    y = y2;
    return Success;
}

template <class T>
typename SubPixelWatersheds<T>::RungeKuttaResult
SubPixelWatersheds<T>::rungeKuttaDoubleStepSecondOrder(
                  double x0, double y0, double & h,
                  double & x, double & y, double epsilon, double dx, double dy)
{
    double x1, x2, y1, y2;
    if(rungeKuttaStepSecondOrder(x0, y0, 2.0 * h, x1, y1, dx, dy) == Outside ||
       rungeKuttaStepSecondOrder(x0, y0, h, x2, y2, dx, dy) == Outside ||
       rungeKuttaStepSecondOrder(x2, y2, h, x2, y2) == Outside)
    {
        h /= 4.0;
        return Outside;
    }
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

template <class T>
int
SubPixelWatersheds<T>::flowLine(double x, double y, bool forward, double epsilon,
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
if(DEBUG) std::cerr << "x, y, dx, dy " << x << ' ' << y << ' ' << dx << ' ' << dy << '\n';
    // check if near a maximum
    if(nearestMaximum(x, y, dx, dy, index) < initialStep_)
    {
        curve.push_back(maxima_[index]);
if(DEBUG) std::cerr << "stop index, x, y " << index << ' ' << curve.back()[0] << ' ' << curve.back()[1]<< '\n';
        return index;
    }

    for(int k=0; k<1000; ++k)
    {
        double xn, yn;
        if(rungeKuttaDoubleStepSecondOrder(x, y, h, xn, yn, epsilon) == Success)
        {
            x = xn;
            y = yn;
            curve.push_back(PointType(x, y));
            dx = image_.dx(x, y);
            dy = image_.dy(x, y);
    if(DEBUG) std::cerr << "x, y, dx, dy " << x << ' ' << y << ' ' << dx << ' ' << dy << '\n';
            // check if near a maximum
            if(nearestMaximum(x, y, dx, dy, index) < initialStep_)
            {
                curve.push_back(maxima_[index]);
        if(DEBUG) std::cerr << "stop index, x, y " << index << ' ' << curve.back()[0] << ' ' << curve.back()[1]<< '\n';
                return index;
            }
        }
        else if(h < 1.0e-6)
        {
        if(DEBUG) std::cerr << "give up\n";
            return 0; // give up
        }
    }
        if(DEBUG) std::cerr << "too many steps\n";
    return 0;
}

template <class T>
pair<int, int>
SubPixelWatersheds<T>::findEdge(double x, double y, double epsilon, PointArray & edge)
{
    PointArray forwardCurve, backwardCurve;
    
        if(DEBUG) std::cerr << "starting forward\n";
    int findex = flowLine(x, y, true, epsilon, forwardCurve);
        if(DEBUG) std::cerr << "starting backward\n";
    int bindex = flowLine(x, y, false, epsilon, backwardCurve);
    
    int i = forwardCurve.size() - 1;
    if(i > 0)
    {
        edge.push_back(forwardCurve[i--]);
        edge.push_back(forwardCurve[i--]);
        for(; i >= 0 ; --i)
        {
            if((forwardCurve[i] - edge.back()).magnitude() >= storeStep_)
                edge.push_back(forwardCurve[i]);
        }
        if(forwardCurve.front() != edge.back())
                edge.push_back(forwardCurve.front());
    }
    else
    {
        edge.push_back(forwardCurve.front());
    }
    for(i = 1; i < (int)backwardCurve.size()-2; ++i)
    {
        if((backwardCurve[i] - edge.back()).magnitude() >= storeStep_)
            edge.push_back(backwardCurve[i]);
    }
    if(i < (int)backwardCurve.size())
        edge.push_back(backwardCurve[i++]);
    if(i < (int)backwardCurve.size())
        edge.push_back(backwardCurve[i++]);
    return pair<int, int>(findex, bindex);
}

template <class T>
int
SubPixelWatersheds<T>::findEdges(double threshold, double epsilon)
{
    edges_.clear();
    edgeIndices_.clear();
    for(int i = 0; i < saddles_.size(); ++i)
    {
        double x = saddles_[i][0];
        double y = saddles_[i][1];
        if(image_(x, y) < threshold)
            continue;
        edges_.push_back(PointArray());
        pair<int, int> ind = findEdge(x, y, epsilon, edges_.back());
        edgeIndices_.push_back(triple<int, int, int>(ind.first, i, ind.second));
    }
    return edges_size();
}

} // namespace vigra

#endif // VIGRA_SUBPIXEL_WATERSHED_HXX