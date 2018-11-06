#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <set>
#include <iomanip>
#include "pwlfunc.hpp"
#include "pwlbound.hpp"

using namespace snowgoose::pwl;

void test_pwlbounds()
{
	std::cout << "*****test_pwlbounds*****" << std::endl;
	PwlFunc<double> lower1({ { 1,::sin(1) },{ 2, ::sin(2) },{ 3,::sin(3) } });
	PwlFunc<double> upper1({ { 1 , ::sin(1) },{ 1.83 , 0.54*1.83 + 0.3 },{ 3,::sin(3) } });
	PwlBound<double> left(lower1, upper1, 5);
	PwlFunc<double> lower2({ { 1,0 },{ 1.83, -1.83*1.83 + 1.83 },{ 3,-6 } });
	PwlFunc<double> upper2({ { 1,0 },{ 2,-1 },{ 3,-6 } });
	PwlBound<double> right(lower2, upper2, 5);
	PwlBound<double> rezult = left * right;
	std::cout << "bounds: " << rezult << std::endl;
}

void test(double c)
{
	double t = floor(c / M_PI);
}

void test_composition() {
	std::cout << "*****sin(PI*cos(x))*****" << std::endl;
	PwlBound<double> x(-2 * M_PI, 0, 4);
	PwlBound<double> y = sin(M_PI*cos(x));
	std::cout << "sin(PI*cos(x)) bounds: " << y << std::endl;
}

void test_division() {
	std::cout << "*****sin(x)/(x^2 + 1)*****" << std::endl;
	PwlBound<double> x(-2 * M_PI+ M_PI/4, 2 * M_PI + M_PI / 4, 8);
	PwlBound<double> y = sin(x) / ((x ^ 2) + 1.0);
	std::cout << std::setprecision(3);
	std::cout << "sin(x)/(x^2 + 1) bounds: " << y << std::endl;
}

void test_division2() {
	std::cout << "*****ln(x)/exp(x)*****" << std::endl;
	PwlBound<double> x(0.1, 2);
	PwlBound<double> y = ln(x) / exp(x);
	std::cout << std::setprecision(3);
	std::cout << "ln(x)/exp(x) bounds: " << y << std::endl;
}

void test_multiplication() {
	std::cout << "*****sin(x)*(-x^2 + x)*****" << std::endl;
	PwlBound<double> x(1, 3, 8);
	PwlBound<double> y = sin(x)*(-1.0*(x ^ 2) + x);
	//PwlBound<double> y = sin(x);
	std::cout << std::setprecision(3);
	std::cout << "sin(x)*(-x^2 + x) bounds: " << y << std::endl;
}

void test_exp() {
	std::cout << "*****-exp(x^3 - x^2)*****" << std::endl;
	PwlBound<double> x(0, 2, 2);
	PwlBound<double> y = -1.0 * exp((x ^ 3) - (x ^ 2));
	std::cout << std::setprecision(3);
	std::cout << "-exp(x^3 - x^2) bounds: " << y << std::endl;
}

void test_sqr()
{
	std::cout << "*****sin(x^3 - sqr(x))*****" << std::endl;
	PwlBound<double> x(-1, 1.5, 4);
	PwlBound<double> y = sin((x^3)-sqr(x));
	std::cout << std::setprecision(3);
	std::cout << "sin(x^3 - sqr(x)) bounds: " << y << std::endl;
}

void test_sqrt() {
	std::cout << "*****sqrt(sin(x))*****" << std::endl;
	PwlBound<double> x(0, 3, 4);
	PwlBound<double> y = sqrt(sin(x));
	std::cout << std::setprecision(3);
	std::cout << "sqrt(sin(x)) bounds: " << y << std::endl;
}

void test_asin() {
	std::cout << "*****asin(x)*****" << std::endl;
	PwlBound<double> x(-1, 1, 4);
	PwlBound<double> y = asin(x);
	std::cout << std::setprecision(3);
	std::cout << "asin(x) bounds: " << y << std::endl;
}


void test_asin2() {
	std::cout << "*****asin(asin(x))*****" << std::endl;
	PwlBound<double> x(-0.8, 0.8, 4);
	PwlBound<double> y = asin(asin(x));
	std::cout << std::setprecision(3);
	std::cout << "asin(asin(x)) bounds: " << y << std::endl;
}

void test_cos() {
	std::cout << "*****cos(asin(x))*****" << std::endl;
	PwlBound<double> x(-1, 1, 16);
	PwlBound<double> y = cos(asin(x));
	std::cout << std::setprecision(3);
	std::cout << "cos(asin(x)) bounds: " << y << std::endl;
}

void test_acos() {
	std::cout << "*****sqrt(acos(x))*****" << std::endl;
	PwlBound<double> x(-1, 1, 4);
	PwlBound<double> y = sqrt(acos(x));
	std::cout << std::setprecision(3);
	std::cout << "sqrt(acos(x)) bounds: " << y << std::endl;
}

void test_tg() {
	std::cout << "*****sqr(x)*tg(x)*****" << std::endl;
	PwlBound<double> x(-1.5, 1.5, 4);
	PwlBound<double> y = sqr(x)*tg(x);
	std::cout << std::setprecision(3);
	std::cout << "sqr(x)*tg(x) bounds: " << y << std::endl;
}

void test_tg2() {
	std::cout << "*****tg(x)*****" << std::endl;
	PwlBound<double> x(-4.6, -1.6, 4);
	PwlBound<double> y = tg(x);
	std::cout << std::setprecision(3);
	std::cout << "tg(x) bounds: " << y << std::endl;
}

void test_ctg() {
	std::cout << "*****sqr(x)*ctg(x)*****" << std::endl;
	PwlBound<double> x(0.1, 3, 4);
	PwlBound<double> y = ctg(x);
	std::cout << std::setprecision(3);
	std::cout << "sqr(x)*ctg(x) bounds: " << y << std::endl;
}

void test_ctg2() {
	std::cout << "*****ctg(x)*****" << std::endl;
	PwlBound<double> x(-3, -0.1, 4);
	PwlBound<double> y = ctg(x);
	std::cout << std::setprecision(3);
	std::cout << "ctg(x) bounds: " << y << std::endl;
}

void test_atg() {
	std::cout << "*****atg(x)*****" << std::endl;
	PwlBound<double> x(-3, 3, 6);
	PwlBound<double> y = atg(x);
	std::cout << std::setprecision(3);
	std::cout << "atg(x) bounds: " << y << std::endl;
}

void test_actg() {
	std::cout << "*****actg(x)*****" << std::endl;
	PwlBound<double> x(-3, 3, 6);
	PwlBound<double> y = actg(x);
	std::cout << std::setprecision(3);
	std::cout << "actg(x) bounds: " << y << std::endl;
}

void test_actg2() {
	std::cout << "*****sqr(actg(x))*****" << std::endl;
	PwlBound<double> x(-3, 3, 6);
	PwlBound<double> y = sqr(actg(x));
	std::cout << std::setprecision(3);
	std::cout << "sqr(actg(x)) bounds: " << y << std::endl;
}

void test_ln() {
	std::cout << "*****ln(x)*sin(x)*****" << std::endl;
	PwlBound<double> x(0.1, 3, 6);
	PwlBound<double> y = ln(x)*sin(x);
	std::cout << std::setprecision(3);
	std::cout << "ln(x)*sin(x) bounds: " << y << std::endl;
}

void test_log() {
	std::cout << "*****log(x,5.0)*****" << std::endl;
	PwlBound<double> x(0.1, 3, 6);
	PwlBound<double> y = log(x, 5.0);
	std::cout << std::setprecision(3);
	std::cout << "log(x,5.0) bounds: " << y << std::endl;
}

void test_minmax() {
	PwlBound<double> x(0.1, 5, 5);
	std::cout << std::setprecision(3);
	PwlBound<double> y = min(-1.0*(x ^ 2) + 5.0, (x^2)-5.0);
	std::cout << "min(-1.0*(x ^ 2) + 5.0, (x^2)-5.0) bounds: " << y << std::endl;

	y = max(-1.0*(x ^ 2) + 5.0, (x ^ 2) - 5.0);
	std::cout << "max(-1.0*(x ^ 2) + 5.0, (x^2)-5.0) bounds: " << y << std::endl;
}


int main()
{
	test_minmax();
	test_exp();
	test_log();
	test_ln();
	test_actg2();
	test_actg();
	test_atg();
	test_ctg2();
	test_tg2();
	test_ctg();
	test_tg();
	test_acos();
	test_cos();
	test_asin2();
	test_asin();
	test_sqrt();
	test_sqr();
	test_multiplication();
	test_division();
	test_division2();
	test_composition();

	PwlFunc<double> upper_bound_x3({ {0,0}, {1.5,3.375}, {2,8} });
	PwlFunc<double> upper_bound_x2({ {0,0}, {1.5,2.25}, {2,4} });
	PwlFunc<double> lower_bound_x3({ { 0,0 }, {4.0/3,0}, {2,8} });
	PwlFunc<double> lower_bound_x2({ { 0,0 }, {1,0}, {2,4} });

	PwlFunc<double> lower_bound = lower_bound_x3 - upper_bound_x2;
	PwlFunc<double> upper_bound = upper_bound_x3 - lower_bound_x2;

	std::cout << "The lower bound x^3-X^2 is \n" << lower_bound << std::endl;
	std::cout << "The upper bound x^3-X^2 is \n" << upper_bound << std::endl;

	PwlFunc<double> lower_bound_2 = 2.0 * lower_bound;
	std::cout << "The lower bound 2*(x^3-X^2) is \n" << lower_bound_2 << std::endl;

	PwlFunc<double> upper_bound_sin({ { 1, 0.84 },{ 1.83, 1.29 },{ 3, 0.14 } });
	PwlFunc<double> lower_bound_x_x2({ { 1, 0.0 },{ 1.83, -1.51 },{ 3, -5.99 } });
	std::vector<Parabola<double>> parabols = upper_bound_sin * lower_bound_x_x2;
	for (auto p : parabols)
		std::cout << p;

	PwlFunc<double> a({ { -2 * M_PI, 0},{ -1.5 * M_PI, 1 },{ -M_PI, 0 } });
	PwlFunc<double> b({ { -2 * M_PI, 4 * M_PI * M_PI +1},{ -M_PI, M_PI * M_PI + 1 } });
	std::vector<FracLinFunc<double>> c = a / b;
	for (auto fl : c)
		std::cout << fl;

	PwlFunc<double> in_f({ { -2 * M_PI, 0 }, { -1.5*M_PI, 1 }, { -0.5*M_PI, -1 }, { 0,0 } });
	PwlFunc<double> out_f({ {-1, 1}, {-0.5, 0}, {0.5, 0}, {1,1} });
	PwlFunc<double> comp = out_f(in_f);
	std::cout << "compound function is " << comp << std::endl;


	PwlFunc<double> pf1({ {0,3},{3,10},{9,3},{20,5} });
	PwlFunc<double> pf2({ { 0,5 },{ 6,9 },{ 20,0 } });
	PwlFunc<double> pf3({ { 0,9 },{ 9,1 },{ 20,9 } });
	PwlFunc<double> pf4({ { 0,9 },{ 7,-1 },{ 20,3 } });

	PwlFunc<double> pf_max = max<double>({ pf1, pf2, pf3, pf4 });
	std::cout << "max: " << pf_max << std::endl;
	PwlFunc<double> pf_min = min<double>({ pf1, pf2, pf3, pf4 });
	std::cout << "min: " << pf_min << std::endl;

	return 0;
}
