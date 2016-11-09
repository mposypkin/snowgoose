/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testprojecttosimplex.cpp
 * Author: mposypkin
 *
 * Created on November 9, 2016, 4:01 PM
 */

#include <iostream>
#include "projecttosimplex.hpp"
/*
 * 
 */
int main(int argc, char** argv) {

    const double C = 3.;
    const int N = 3;
    snowgoose::ProjectToSimplex<double> prj(C);
    double x[N] = {0, 0, 0};
    prj.project(N, x);
    std::cout << "x = " << snowgoose::VecUtils::vecPrint(N, x) << "\n";
    return 0;
}

