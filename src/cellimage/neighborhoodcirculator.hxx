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

#ifndef VIGRA_NEIGHBORHOODCIRCULATOR_HXX
#define VIGRA_NEIGHBORHOODCIRCULATOR_HXX

#include <pixelneighborhood.hxx>

namespace vigra {

template <class ImageIterator>
class Neighborhood8Circulator : private ImageIterator
{
public:
    typedef typename ImageIterator::value_type value_type;

    Neighborhood8Circulator(ImageIterator const & it,
                            EightNeighborCoding::Directions d = EightNeighborCoding::East)
        : ImageIterator(it), neighbor(d)
    {
        operator+=(neighbor.diff());
    }

    Neighborhood8Circulator & operator++()
    {
        operator+=(neighbor.relativeNext());
        neighbor++;
        return *this;
    }

    Neighborhood8Circulator operator++(int)
    {
        Neighborhood8Circulator ret(*this);
        operator++();
        return ret;
    }

    Neighborhood8Circulator & operator--()
    {
        operator+=(neighbor.relativePrev());
        neighbor++;
        return *this;
    }

    Neighborhood8Circulator operator--(int)
    {
        Neighborhood8Circulator ret(*this);
        operator--();
        return ret;
    }

    bool operator==(Neighborhood8Circulator const & rhs) const
    {
        return neighbor == rhs.neighbor;
    }

    bool operator!=(Neighborhood8Circulator const & rhs) const
    {
        return neighbor != rhs.neighbor;
    }

    typename ImageIterator::reference operator*() const
    {
        return ImageIterator::operator*();
    }

private:
    EightNeighborCoding neighbor;
};

template <class ImageIterator>
class Neighborhood4Circulator : private ImageIterator
{
public:
    typedef typename ImageIterator::value_type value_type;

    Neighborhood4Circulator(ImageIterator const & it,
                            EightNeighborCoding::Directions d = EightNeighborCoding::East)
        : ImageIterator(it), neighbor(d)
    {
        operator+=(neighbor.diff());
    }

    Neighborhood4Circulator & operator++()
    {
        operator+=(neighbor.relativeNext());
        neighbor++;
        return *this;
    }

    Neighborhood4Circulator operator++(int)
    {
        Neighborhood4Circulator ret(*this);
        operator++();
        return ret;
    }

    Neighborhood4Circulator & operator--()
    {
        operator+=(neighbor.relativePrev());
        --neighbor;
        return *this;
    }

    Neighborhood4Circulator operator--(int)
    {
        Neighborhood4Circulator ret(*this);
        operator--();
        return ret;
    }

    bool operator==(Neighborhood4Circulator const & rhs) const
    {
        return neighbor == rhs.neighbor;
    }

    bool operator!=(Neighborhood4Circulator const & rhs) const
    {
        return neighbor != rhs.neighbor;
    }

    value_type operator*() const
    {
        return ImageIterator::operator*();
    }

private:
    FourNeighborCoding neighbor;
};

} // namespace vigra

#endif /* VIGRA_NEIGHBORHOODCIRCULATOR_HXX */
