#ifndef VIGRA_NEIGHBORHOODCIRCULATOR_HXX
#define VIGRA_NEIGHBORHOODCIRCULATOR_HXX

#include <pixelneighborhood.hxx>

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
        ++neighbor;
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
        --neighbor;
        operator+=(neighbor.diff());
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
    
    value_type & operator*()
    {
        return ImageIterator::operator*();
    }
    
    value_type operator*() const
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
        ++neighbor;
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
    
    value_type & operator*()
    {
        return ImageIterator::operator*();
    }
    
    value_type operator*() const
    {
        return ImageIterator::operator*();
    }
    
  private:
    
    FourNeighborCoding neighbor;
};



#endif /* VIGRA_NEIGHBORHOODCIRCULATOR_HXX */
