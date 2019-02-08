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

#ifndef CROP_HXX
#define CROP_HXX

namespace vigra {

    // FIXME: Size2D is a workaround until Iterator::operator- is adopted
template<class Iterator, class Accessor>
inline triple<Iterator, Iterator, Accessor>
crop(triple<Iterator, Iterator, Accessor> t, Rect2D r)
{
    //vigra_precondition(Rect2D(Point2D(0, 0), Size2D((t.second - t.first).x, (t.second - t.first).y)).contains(r),
    //                   "crop called with rect which is out of (iterator) bounds!");
    return triple<Iterator, Iterator, Accessor>(t.first + r.upperLeft(),
                                                t.first + r.lowerRight(),
                                                t.third);
}

template<class Iterator, class Accessor>
inline pair<Iterator, Accessor>
crop(pair<Iterator, Accessor> t, Rect2D r)
{
    return pair<Iterator, Accessor>(t.first + r.upperLeft(), t.second);
}

} // namespace vigra

#endif // CROP_HXX
