#ifndef LABELLUT_HXX
#define LABELLUT_HXX

#include <vector>

class LabelLUT
{
  public:
    typedef unsigned int LabelType;
    typedef std::vector<LabelType> LUTType;
    typedef LUTType::size_type size_type;

    LabelLUT()
    {}

    LabelLUT(unsigned int size)
    : labelLUT_(size)
    {
        initIdentity(size);
    }

    void initIdentity(unsigned int size)
    {
        labelLUT_.resize(size);
        for(unsigned int i = 0; i < size; ++i)
            labelLUT_[i] = i;
    }

    void appendOne()
    {
        labelLUT_.push_back(labelLUT_.size());
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
        for(unsigned int i = 0; i < labelLUT_.size(); ++i)
        {
            if(labelLUT_[i] == from)
                labelLUT_[i] = to;
        }
    }

  protected:
    std::vector<LabelType> labelLUT_;

};

#endif // LABELLUT_HXX
