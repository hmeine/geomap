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
        LabelType prev, toEnd;

        for(toEnd = to; true; toEnd = prev)
        {
            // find end of "to" list (for later concatenation)
            prev = prevMerged_[toEnd];
            if(prev == toEnd)
                break;
        }

        for(LabelType fromIt = from; true; fromIt = prev)
        {
            // relabel elements in "from" list:
            labelLUT_[fromIt] = to;

            prev = prevMerged_[fromIt];
            if(prev == fromIt)
                break;
        }

        // concatenate lists:
        prevMerged_[toEnd] = from;
    }

  protected:
    LUTType labelLUT_;
    LUTType prevMerged_;
};

#endif // LABELLUT_HXX
