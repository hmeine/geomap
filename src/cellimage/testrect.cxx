#include <vigra/diff2d.hxx>
#include <unittest.hxx>
#include <iostream>

using namespace vigra;

struct Point2DTest
{
    Point2D p11;

    Point2DTest()
        : p11(1, 1)
    {
    }

    void testOperations()
    {
        should(-p11 == Point2D(-1, -1));
    }
};

struct Rect2DTest
{
    Rect2D rect1_1;
    Rect2D emptyRect;
    Rect2D bigRect;

    Rect2DTest()
        : rect1_1(Point2D(1, 1), Point2D(2, 2)),
          bigRect(0, 0, 10, 10)
    {
    }

    void testProperties()
    {
        shouldEqual(rect1_1.width(), 1);
        shouldEqual(rect1_1.height(), 1);
        should(!rect1_1.isEmpty());

        shouldEqual(emptyRect.width(), 0);
        shouldEqual(emptyRect.height(), 0);
        should(emptyRect.isEmpty());

        shouldEqual(bigRect.width(), 10);
        shouldEqual(bigRect.height(), 10);
        should(!bigRect.isEmpty());

        should(rect1_1 != emptyRect);
        should(bigRect != emptyRect);
        should(bigRect != rect1_1);

        bigRect = rect1_1;
        should(bigRect == rect1_1);
    }

    void testContains()
    {
        should(!emptyRect.contains(Point2D(0, 0)));
        should(!emptyRect.contains(Point2D(0, 1)));
        should(!emptyRect.contains(Point2D(0, 2)));
        should(!emptyRect.contains(Point2D(1, 0)));
        should(!emptyRect.contains(Point2D(1, 1)));
        should(!emptyRect.contains(Point2D(1, 2)));
        should(!emptyRect.contains(Point2D(2, 0)));
        should(!emptyRect.contains(Point2D(2, 1)));
        should(!emptyRect.contains(Point2D(2, 2)));

        should( emptyRect.contains(emptyRect));
        should(!emptyRect.contains(rect1_1));
        should(!emptyRect.contains(bigRect));

        should(!rect1_1.contains(Point2D(0, 0)));
        should(!rect1_1.contains(Point2D(0, 1)));
        should(!rect1_1.contains(Point2D(0, 2)));
        should(!rect1_1.contains(Point2D(1, 0)));
        should( rect1_1.contains(Point2D(1, 1)));
        should(!rect1_1.contains(Point2D(1, 2)));
        should(!rect1_1.contains(Point2D(2, 0)));
        should(!rect1_1.contains(Point2D(2, 1)));
        should(!rect1_1.contains(Point2D(2, 2)));

        should( rect1_1.contains(emptyRect));
        should( rect1_1.contains(rect1_1));
        should(!rect1_1.contains(bigRect));

        should(bigRect.contains(Point2D(0, 0)));
        should(bigRect.contains(Point2D(0, 1)));
        should(bigRect.contains(Point2D(0, 2)));
        should(bigRect.contains(Point2D(1, 0)));
        should(bigRect.contains(Point2D(1, 1)));
        should(bigRect.contains(Point2D(1, 2)));
        should(bigRect.contains(Point2D(2, 0)));
        should(bigRect.contains(Point2D(2, 1)));
        should(bigRect.contains(Point2D(2, 2)));

        should( bigRect.contains(emptyRect));
        should( bigRect.contains(rect1_1));
        should( bigRect.contains(bigRect));
    }

    void testIntersection()
    {
        should(!emptyRect.intersects(emptyRect));
        should(!emptyRect.intersects(rect1_1));
        should(!emptyRect.intersects(bigRect));
        should(!rect1_1.intersects(emptyRect));
        should( rect1_1.intersects(rect1_1));
        should( rect1_1.intersects(bigRect));
        should(!bigRect.intersects(emptyRect));
        should( bigRect.intersects(rect1_1));
        should( bigRect.intersects(bigRect));

        should(!bigRect.intersects(Rect2D(Point2D(3, -3), Point2D(3, 3))));
        should( bigRect.intersects(Rect2D(Point2D(3, -3), Point2D(4, 3))));
        should( bigRect.intersects(Rect2D(Point2D(3, -3), Point2D(14, 3))));

        should((rect1_1 & emptyRect).isEmpty());
        should(!(rect1_1 & bigRect).isEmpty());
        should((rect1_1 & bigRect) == rect1_1);
    }

    void testUnion()
    {
        should(!(rect1_1 | emptyRect).isEmpty());
        should((rect1_1 | emptyRect) == rect1_1);
        should((rect1_1 | bigRect) == bigRect);
        rect1_1 |= Point2D(3, 3);
        shouldEqual(rect1_1.upperLeft(), Diff2D(1, 1));
        shouldEqual(rect1_1.lowerRight(), Diff2D(4, 4));
    }

    void testSizes()
    {
        shouldEqual(rect1_1.size(), Size2D(1, 1));
        shouldEqual(bigRect.size(), Size2D(10, 10));
        emptyRect.setSize(10, 10);
        should(bigRect == emptyRect);
        emptyRect.addSize(Size2D(-4, -7));
        shouldEqual(emptyRect.size(), Size2D(6, 3));
        emptyRect.setSize(bigRect.size());
        should(bigRect == emptyRect);
    }
};

struct Point2DTestSuite
: public test_suite
{
    Point2DTestSuite()
    : test_suite("Point2DTestSuite")
    {
        add(testCase(&Point2DTest::testOperations));
    }
};

struct Rect2DTestSuite
: public test_suite
{
    Rect2DTestSuite()
    : test_suite("Rect2DTestSuite")
    {
        add(testCase(&Rect2DTest::testProperties));
        add(testCase(&Rect2DTest::testIntersection));
        add(testCase(&Rect2DTest::testUnion));
        add(testCase(&Rect2DTest::testSizes));
    }
};

int main(int argc, char **argv)
{
    Rect2DTestSuite rect2dTestSuite;
    int failed = rect2dTestSuite.run();
    std::cout << rect2dTestSuite.report() << std::endl;
    Point2DTestSuite point2dTestSuite;
    failed += point2dTestSuite.run();
    std::cout << point2dTestSuite.report() << std::endl;
    return failed;
}
