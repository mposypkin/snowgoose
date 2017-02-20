/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testinterval.cpp
 * Author: Alexander Usov
 *
 * Created on January 20, 2017, 9:13 AM
 */
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cstdlib>
#include "interval_air.hpp"

using namespace std;
using namespace snowgoose::interval;

template <class T> 
class A
{
    private:
        T f;
    public:
  A(T t) : f(t), f2(t){}
  protected:
      void Test2() {}
      T f2;
};

template <class T> 
class B : public A<T>
{
    public:
    B(T t) : A<T>(t) {}
    void test() 
    { 
        this->f2 = 9.0; 
        this->Test2(); 
    }
};

int main(int argc, char** argv) {

    Interval<double> a(2.0, 3.0), b(1, 7);
    auto c = sqr(a);
    std::cout << c;
    auto cc = a ^ 3;
    std::cout << cc;
    auto ccc = sin(a);
    std::cout << ccc;
    auto cccc = cos(a);
    std::cout << cccc;
    auto dd = exp(a);
    std::cout << dd;
    auto ddd = sqrt(a);
    std::cout << ddd;
    auto dddd = a / b;
    std::cout << dddd;
    std::cout << a + b;
    std::cout << a - b;
    std::cout << a * b;

    Interval<double> x(0.9, 1.1), y(-0.1, 0.1);
    Interval<double> z = -20.0 * exp(-0.2 * sqrt(0.5 * (sqr(x) + sqr(y)))) - exp(0.5*(cos(2.0 * M_PI * x) + cos(2.0 * M_PI * y))) + 20 + M_E;
    std::cout << z;
    Interval<double> z2 = sin(1.0/x);
    std::cout << z2;
    
    std::cout << abs(Interval<double>(-3.0, 1.0));
    std::cout << ln(Interval<double>(M_E, M_E * M_E));
    std::cout << log(Interval<double>(100.0, 1000.0), 10.0);

    return 0;
}

