#ifndef LABELLUT_HXX
#define LABELLUT_HXX

#include <vector>

class LabelLUT
{
  public:
    typedef unsigned int           LabelType;
    typedef std::vector<LabelType> LUTType; // internal array type
    typedef LUTType::size_type     size_type;
    typedef LabelType              value_type;

    class MergedIterator
    {
        const LUTType &prevMerged_;
        LabelType currentLabel_;
        bool atEnd_;

      public:
        typedef LabelType value_type;
        typedef LabelType reference;
        typedef std::forward_iterator_tag iterator_category;

        MergedIterator(const LUTType &prevMerged, LabelType start)
        : prevMerged_(prevMerged),
          currentLabel_(start),
          atEnd_(false)
        {
        }

        LabelType operator*() const
        {
            return currentLabel_;
        }

        MergedIterator & operator++()
        {
            LabelType next = prevMerged_[currentLabel_];
            if(next == currentLabel_)
                atEnd_ = true;
            else
                currentLabel_ = next;
            return *this;
        }
        
        MergedIterator operator++(int)
        {
            MergedIterator ret(*this);
            operator++();
            return ret;
        }

        bool atEnd() const
        {
            return atEnd_;
        }

        bool inRange() const
        {
            return !atEnd_;
        }
    };

    LabelLUT()
    {}

    LabelLUT(unsigned int size)
    : labelLUT_(size),
      prevMerged_(size)
    {
        initIdentity(size);
    }

    void initIdentity(unsigned int size)
    {
        labelLUT_.resize(size);
        prevMerged_.resize(size);
        for(unsigned int i = 0; i < size; ++i)
        {
            labelLUT_[i] = i;
            prevMerged_[i] = i;
        }
    }

    void appendOne()
    {
        labelLUT_.push_back(labelLUT_.size());
        prevMerged_.push_back(prevMerged_.size());
    }

    LabelType operator[](size_type index) const
    {
        return labelLUT_[index];
    }

    size_type size() const
    {
        return labelLUT_.size();
    }

    void relabel(LabelType from, LabelType to)
    {
        // relabel elements in "from" list:
        LabelType prev;
        for(LabelType fromIt = from; true; fromIt = prev)
        {
            labelLUT_[fromIt] = to;

            prev = prevMerged_[fromIt];
            if(prev == fromIt)
                break;
        }

        // insert from-list at beginning of to-list:
        if(prevMerged_[to] != to)
            prevMerged_[prev] = prevMerged_[to];
        prevMerged_[to] = from;
    }

    MergedIterator mergedBegin(LabelType start) const
    {
        return MergedIterator(prevMerged_, start);
    }

  protected:
    LUTType labelLUT_;
    LUTType prevMerged_;
};

#endif // LABELLUT_HXX
