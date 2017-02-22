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
#include "expr.hpp"
#include "testfuncs.hpp"
#include "interval/interval_air.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace expression;
using namespace snowgoose::interval;

void calcFunc(const std::string& name, Expr<double> expr, const std::vector<double>& vars)
{
	auto result = expr.calc(vars, FuncAlg<double>());	
	std::cout << name << ": " << result << '\n';
}

void calcInterval(const std::string& name, Expr<Interval<double>> expr, const std::vector<Interval<double>>& vars)
{
	auto result = expr.calc(vars, InterEvalAlg<double>());
	std::cout << name << ": " << result;
}

int main(int argc, char** argv) {
	
	auto expr = Ackley1<double>(2);	
	std::cout << "Ackley1: \n" << expr;
	calcFunc("Ackley1", expr, { 0.0, 0.0});
	calcFunc("Ackley1", expr, { 1.0, 1.0 });
	calcFunc("Ackley1", expr, { 2.0, 2.0 });
	
	auto exprInterval = Ackley1<Interval<double>>(2);
	calcInterval("Ackley1 estimation", exprInterval, { { 0.5, 1.5 }, { -0.5, 0.5 } });
	calcInterval("Ackley1 estimation", exprInterval, { { 0.9, 1.1 }, { -0.1, 0.1 } });

	auto brad = Brad<double>();
	calcFunc("Brad", brad, { 0.0824, 1.133, 2.3437 });

	calcFunc("Factorial 5", Factorial<double>(5), { });
	
	auto alpine = Alpine2<double>(1);
	calcFunc("Alpine", alpine, { 7.917 });

	auto alpineInterval = Alpine2<Interval<double>>(1);
	calcInterval("Alpine estimation", alpineInterval, { { 7.5, 8.5 } });

	auto bigg = BiggsExpr2<double>();
	calcFunc("Biggs", bigg, { 1, 10});

	auto biggInterval = BiggsExpr2<Interval<double>>();
	calcInterval("Biggs estimation", biggInterval, { {0, 20}, { 0, 20 } });

	auto bird = Bird<double>();
	calcFunc("Bird", bird, { 4.70104, 3.15294 });
	calcFunc("Bird second minimum", bird, { -1.58214, -3.13024});

	auto birdInterval = Bird<Interval<double>>();
	calcInterval("Bird estimation", birdInterval, { { -2 * M_PI, 2 * M_PI },{ -2 * M_PI, 2 * M_PI } });

	auto boxBetts = BoxBettsQuadraticSum<double>();
	calcFunc("BoxBetts", boxBetts, { 1.0, 10.0, 1.0 });

	auto boxBettsOnterval = BoxBettsQuadraticSum<Interval<double>>();
	calcInterval("BoxBetts estimation", boxBettsOnterval, { { 0.9, 1.2 }, {9, 11.2}, { 0.9, 1.2 } });
	
	auto brown = Brown<double>(4);
	calcFunc("Brown", brown, { 0, 0, 0, 0 });
	
	auto hartman3 = Hartman3<double>();
	calcFunc("Hartman3", hartman3, { 0.114614, 0.555649, 0.852547 });
	auto hartman6 = Hartman6<double>();
	calcFunc("Hartman6", hartman6, {0.201690, 0.150011, 0.476874, 0.275332, 0.311652, 0.657301});

	auto schewefel = Schewefel<Interval<double>>(2, 2);
	calcInterval("Schewefel", schewefel, { { 2.0, 3.0 },{ 2.0, 3.0 } });

	IL<double> list = { { 1.0, 4.0 },{ 2.0, 3.0 } };
	auto minInterval = min(list);

	auto last = Ackley4<double>(2);
	calcFunc("Last", last, { 1.479252, -0.739807});
	std::cout << last << '\n';


    return 0;
}
