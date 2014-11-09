#include "filteriterator.hxx"
#include <vigra/stdimagefunctions.hxx>
#include <unittest.hxx>
#include <algorithm>
#include <functional>
#include <iostream>

struct FilterEven : public std::unary_function<int, bool>
{
    bool operator()(int v) const
    {
        return v%2 == 0;
    }
};

struct FilterOdd : public std::unary_function<int, bool>
{
    bool operator()(int v) const
    {
        return v%2 != 0;
    }
};

struct ImageFilterIteratorTest
{
    vigra::IImage img;

    ImageFilterIteratorTest()
        : img(5, 4)
    {
        static const int imgData[] =
            { 0, 1, 2, 3, 0,
              1, 2, 4, 0, 1,
              2, 6, 0, 7, 2,
              3, 5, 1, 2, 3 };

        std::copy(imgData, imgData + 5*4, img.begin());
    }

    void testInRangeFilteringEven()
    {
        int evenData[] = { 0, 2, 0, 2, 4, 0, 2, 6, 0, 2, 2 };
        int *even= evenData;
        for(vigra::ImageFilterIterator<FilterEven, vigra::IImage::const_traverser>
                fit(srcImageRange(img));
            fit.inRange(); fit++, even++)
        {
            should(even < evenData + sizeof(evenData) / sizeof(evenData[0]));
            shouldEqual(*fit, *even);
        }
        should(even == evenData + sizeof(evenData) / sizeof(evenData[0]));
    }

    void testAtEndFilteringOdd()
    {
        int oddData[] = { 1, 3, 1, 1, 7, 3, 5, 1, 3 };
        int *odd= oddData;
        for(vigra::ImageFilterIterator<FilterOdd, vigra::IImage::const_traverser>
                fit(img.upperLeft(), img.lowerRight());
            !fit.atEnd(); fit++, odd++)
        {
            should(odd < oddData + sizeof(oddData) / sizeof(oddData[0]));
            shouldEqual(*fit, *odd);
        }
        should(odd == oddData + sizeof(oddData) / sizeof(oddData[0]));
    }
};

struct ImageFilterIteratorTestSuite
: public vigra::test_suite
{
    ImageFilterIteratorTestSuite()
    : vigra::test_suite("ImageFilterIteratorTestSuite")
    {
        add( testCase( &ImageFilterIteratorTest::testInRangeFilteringEven));
        add( testCase( &ImageFilterIteratorTest::testAtEndFilteringOdd));
    }
};

struct FilterIteratorTest
{
    vigra::IImage img;

    FilterIteratorTest()
        : img(5, 4)
    {
        static const int imgData[] =
            { 0, 1, 2, 3, 0,
              1, 2, 4, 0, 1,
              2, 6, 0, 7, 2,
              3, 5, 1, 2, 3 };

        std::copy(imgData, imgData + 5*4, img.begin());
    }

    void testInRangeFilteringEven()
    {
        int evenData[] = { 0, 2, 0, 2, 4, 0, 2, 6, 0, 2, 2 };
        int *even= evenData;
        for(vigra::FilterIterator<vigra::IImage::ScanOrderIterator, FilterEven>
                fit(img.begin(), img.end());
            fit.inRange(); fit++, even++)
        {
            should(even < evenData + sizeof(evenData) / sizeof(evenData[0]));
            shouldEqual(*fit, *even);
        }
        should(even == evenData + sizeof(evenData) / sizeof(evenData[0]));
    }

    void testAtEndFilteringOdd()
    {
        int oddData[] = { 1, 3, 1, 1, 7, 3, 5, 1, 3 };
        int *odd= oddData;
        for(vigra::FilterIterator<vigra::IImage::ScanOrderIterator, FilterOdd>
                fit(img.begin(), img.end());
            !fit.atEnd(); fit++, odd++)
        {
            should(odd < oddData + sizeof(oddData) / sizeof(oddData[0]));
            shouldEqual(*fit, *odd);
        }
        should(odd == oddData + sizeof(oddData) / sizeof(oddData[0]));
    }
};

struct FilterIteratorTestSuite
: public vigra::test_suite
{
    FilterIteratorTestSuite()
    : vigra::test_suite("FilterIteratorTestSuite")
    {
        add( testCase( &FilterIteratorTest::testInRangeFilteringEven));
        add( testCase( &FilterIteratorTest::testAtEndFilteringOdd));
    }
};

int main(int argc, char **argv)
{
    ImageFilterIteratorTestSuite imageFilterIteratorTestSuite;
    int failed = imageFilterIteratorTestSuite.run();
    std::cout << imageFilterIteratorTestSuite.report() << std::endl;

    FilterIteratorTestSuite FilterIteratorTestSuite;
    failed += FilterIteratorTestSuite.run();
    std::cout << FilterIteratorTestSuite.report() << std::endl;

    return failed;
}