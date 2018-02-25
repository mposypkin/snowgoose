/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/* 
 * File:   testseries.hpp
 * Author: alusov
 *
 * Created on June 22, 2017, 11:36 AM
 */

#include <iostream>
#include <math.h>
#include "series.hpp"
#include "interval/interval_air.hpp"
#include "common/sgerrcheck.hpp"

#define EPS 0.001

using namespace std;
using namespace snowgoose::derhighorder;
using namespace snowgoose::interval;


void Test1()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = x * x;
    SG_ASSERT(z.der(2)==2);
    SG_ASSERT(z.der(3)==0);
    std::cout << "Test1: " <<  z;
}


void Test2()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = x ^ 2;
    SG_ASSERT(z.der(2)==2);
    SG_ASSERT(z.der(3)==0);
    std::cout << "Test3: " <<  z;
}


void Test3()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = (x^2) + (x^3);
    SG_ASSERT(z.der(2)==8);
    SG_ASSERT(z.der(3)==6);
    std::cout << "Test4: " <<  z;
}

void Test4()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = (x^3) - (x^2);
    SG_ASSERT(z.der(2)==4);
    SG_ASSERT(z.der(3)==6);
    std::cout << "Test4: " <<  z;
}

void Test5()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = 3.0 * (x^2);
    SG_ASSERT(z.value()==3);
    SG_ASSERT(z.der(2)==6);
    SG_ASSERT(z.der(3)==0);
    std::cout << "Test5: " <<  z;
}

void Test6()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = (x^3)/(x^2);
    SG_ASSERT(z.value()==1);
    SG_ASSERT(z.der(1)==1);
    SG_ASSERT(z.der(2)==0);
    SG_ASSERT(z.der(3)==0);
    std::cout << "Test6: " <<  z;
}

void Test7()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = x^0.5;
    SG_ASSERT(z.value()==1);
    SG_ASSERT(z.der(1)==0.5);
    SG_ASSERT(z.der(2)==-0.25);
    SG_ASSERT(z.der(3)==3.0/8.0);
    std::cout << "Test7: " <<  z;
}

void Test8()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = 0.5 + x;
    SG_ASSERT(z.value()==1.5);
    SG_ASSERT(z.der(1)==1.0);
    SG_ASSERT(z.der(2)==0.0);
    SG_ASSERT(z.der(3)==0.0);
    std::cout << "Test8: " <<  z;
}

void Test9()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = 0.5 - x;
    SG_ASSERT(z.value()==-0.5);
    SG_ASSERT(z.der(1)==-1.0);
    SG_ASSERT(z.der(2)==0.0);
    SG_ASSERT(z.der(3)==0.0);
    std::cout << "Test9: " <<  z;
}

void Test10()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = 5.0 * x;
    SG_ASSERT(z.value()== 5.0);
    SG_ASSERT(z.der(1)==5.0);
    SG_ASSERT(z.der(2)==0.0);
    SG_ASSERT(z.der(3)==0.0);
    std::cout << "Test10: " <<  z;
}

void Test11()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = 5.0 / x;
    SG_ASSERT(z.value()== 5.0);
    SG_ASSERT(z.der(1)==-5.0);
    SG_ASSERT(z.der(2)==10.0);
    SG_ASSERT(z.der(3)==-30.0);
    std::cout << "Test11: " <<  z;
}

void Test12()
{
    Series<double> x(5.0, 1.0, 3);
    Series<double> z = sqr(x);
    SG_ASSERT(z.value()== 25.0);
    SG_ASSERT(z.der(1)==10.0);
    SG_ASSERT(z.der(2)==2.0);
    SG_ASSERT(z.der(3)==0.0);
    std::cout << "Test12: " <<  z;
}

void Test13()
{
    Series<double> x(1.0, 1.0, 3);
    Series<double> z = sqrt(x);
    SG_ASSERT(z.value()==1);
    SG_ASSERT(z.der(1)==0.5);
    SG_ASSERT(z.der(2)==-0.25);
    SG_ASSERT(z.der(3)==3.0/8.0);
    std::cout << "Test13: " <<  z;
}

void Test14()
{
    Series<double> x(M_PI/2.0, 1.0, 3);
    Series<double> z = sin(x);
    SG_ASSERT(z.value()==1.0);
    SG_ASSERT_NEAR(z.der(1), 0.0, EPS);
    SG_ASSERT_NEAR(z.der(2), -1.0, EPS);
    SG_ASSERT_NEAR(z.der(3), 0.0, EPS);
    std::cout << "Test14: " <<  z;
}

void Test15()
{
    Series<double> x(M_PI, 1.0, 3);
    Series<double> z = cos(x);
    SG_ASSERT(z.value()==-1.0);
    SG_ASSERT_NEAR(z.der(1), 0.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 1.0, EPS);
    SG_ASSERT_NEAR(z.der(3), 0.0, EPS);
    std::cout << "Test15: " <<  z;
}

void Test16()
{
    Series<double> x(M_PI/4.0, 1.0, 3);
    Series<double> z = tg(x);
    SG_ASSERT_NEAR(z.value(), 1.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 2.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 4.0, EPS);
    SG_ASSERT_NEAR(z.der(3), 16.0, EPS);
    std::cout << "Test16: " <<  z;
}

void Test17()
{
    Series<double> x(M_PI/4.0, 1.0, 3);
    Series<double> z = ctg(x);

    SG_ASSERT_NEAR(z.value(), 1.0, EPS);
    SG_ASSERT_NEAR(z.der(1), -2.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 4.0, EPS);
    SG_ASSERT_NEAR(z.der(3), -16.0, EPS);

    std::cout << "Test17: " <<  z;
}

void Test18()
{
    Series<double> x(0.5, 1.0, 2);
    Series<double> z = asin(x);

    SG_ASSERT_NEAR(z.value(), M_PI/6.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 2.0/std::sqrt(3.0), EPS);
    SG_ASSERT_NEAR(z.der(2), 4.0/(3.0 * std::sqrt(3.0)), EPS);

    std::cout << "Test18: " <<  z;
}

void Test19()
{
    Series<double> x(0.5, 1.0, 2);
    Series<double> z = acos(x);

    SG_ASSERT_NEAR(z.value(), M_PI/3.0, EPS);
    SG_ASSERT_NEAR(z.der(1), -2.0/std::sqrt(3.0), EPS);
    SG_ASSERT_NEAR(z.der(2), -4.0/(3.0 * std::sqrt(3.0)), EPS);

    std::cout << "Test19: " <<  z;
}

void Test20()
{
    Series<double> x(1.0, 1.0, 2);
    Series<double> z = atg(x);

    SG_ASSERT_NEAR(z.value(), M_PI/4.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 0.5, EPS);
    SG_ASSERT_NEAR(z.der(2), -0.5, EPS);

    std::cout << "Test20: " <<  z;
}


void Test21()
{
    Series<double> x( M_PI/2.0, 1.0, 2);
    Series<double> z = atg(sin(x));

    SG_ASSERT_NEAR(z.value(), M_PI/4.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 0.0, EPS);
    SG_ASSERT_NEAR(z.der(2), -0.5, EPS);

    std::cout << "Test21: " <<  z;
}

void Test22()
{
    Series<double> x(0.0, 1.0, 2);
    Series<double> z = actg(x);

    SG_ASSERT_NEAR(z.value(), M_PI/2.0, EPS);
    SG_ASSERT_NEAR(z.der(1), -1.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 0.0, EPS);

    std::cout << "Test22: " <<  z;
}

void Test23()
{
    Series<double> x(10.0, 1.0, 2);
    Series<double> z = ln(x);

    SG_ASSERT_NEAR(z.value(), std::log(10), EPS);
    SG_ASSERT_NEAR(z.der(1), 0.1, EPS);
    SG_ASSERT_NEAR(z.der(2), -0.01, EPS);

    std::cout << "Test23: " <<  z;
}

void Test24()
{
    Series<double> x(10.0, 1.0, 2);
    Series<double> z = log(x, 10);

    SG_ASSERT_NEAR(z.value(), 1.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 1.0/(10*std::log(10)), EPS);
    SG_ASSERT_NEAR(z.der(2), -1.0/(100*std::log(10)), EPS);

    std::cout << "Test24: " <<  z;
}

void Test25()
{
    Series<double> x(1.0, 1.0, 2);
    Series<double> z = exp(x);

    SG_ASSERT_NEAR(z.value(), std::exp(1.0), EPS);
    SG_ASSERT_NEAR(z.der(1), std::exp(1.0), EPS);
    SG_ASSERT_NEAR(z.der(2), std::exp(1.0), EPS);

    std::cout << "Test25: " <<  z;
}

void Test26()
{
    Series<double> x(1.0, 1.0, 2);
    Series<double> z = pow(10.0, x);

    SG_ASSERT_NEAR(z.value(), 10.0, EPS);
    SG_ASSERT_NEAR(z.der(1), 10 * std::log(10), EPS);
    SG_ASSERT_NEAR(z.der(2), 10 * std::log(10) * std::log(10), EPS);

    std::cout << "Test26: " <<  z;
}

void Test27()
{
    Series<double> x(-1.0, 1.0, 2);
    Series<double> z = abs(x);

    SG_ASSERT_NEAR(z.value(), std::abs(-1), EPS);
    SG_ASSERT_NEAR(z.der(1), -1.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 0.0, EPS);

    std::cout << "Test27: " <<  z;
}

void Test28()
{
    Series<double> x(-1.0, 1.0, 2);
    Series<double> z = abs(x^2);

    SG_ASSERT_NEAR(z.value(), 1.0, EPS);
    SG_ASSERT_NEAR(z.der(1), -2.0, EPS);
    SG_ASSERT_NEAR(z.der(2), 2.0, EPS);

    std::cout << "Test28: " <<  z;
}



int main(int argc, char** argv) {

    Test1(); 
    Test2();
    Test3();
    Test4();
    Test5();
    Test6();
    Test7();
    Test8();
    Test9();
    Test10();
    Test11();
    Test12();
    Test13();
    Test14();
    Test15();
    Test16();
    Test17();
    Test18();
    Test19();
    Test20();
    Test21();
    Test22();
    Test23();
    Test24();
    Test25();
    Test26();
    Test27();
    Test28();
    return 0;
}
 




