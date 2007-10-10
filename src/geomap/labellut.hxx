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

  protected:
    LUTType labelLUT_;
    LUTType prevMerged_;
};

#endif // LABELLUT_HXX
