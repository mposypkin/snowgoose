/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

/*
* File:   keys.hpp
* Author: alusov
*
* Created on February 23, 2017, 5:42 PM
*/

#ifndef KEYS_HPP
#define KEYS_HPP

#include <string>

struct Keys
{
	const std::string Ackley1 = "Ackley1";
	const std::string Ackley2 = "Ackley2";
	const std::string Ackley3 = "Ackley3";
	const std::string Ackley4 = "Ackley4";
	const std::string Adjiman = "Adjiman";
	const std::string Alpine1 = "Alpine1";
	const std::string Alpine2 = "Alpine2";
	const std::string Brad = "Brad";
	const std::string BartelsConn = "BartelsConn";
	const std::string Beale = "Beale";
	const std::string BiggsEXP2 = "BiggsEXP2";
	const std::string BiggsEXP3 = "BiggsEXP3";
	const std::string BiggsEXP4 = "BiggsEXP4";
	const std::string BiggsEXP5 = "BiggsEXP5";
	const std::string BiggsEXP6 = "BiggsEXP6";
	const std::string Bird = "Bird";

	const std::string desc = "description";
	const std::string anyDim = "anyDim";
	const std::string dim = "dim";
	const std::string globMinY = "globMinY";
	const std::string bounds = "bounds";
	const std::string lb = "a";
	const std::string rb = "b";
	const std::string globMin = "globMin";
	const std::string x = "x";
};

static const Keys K;

#endif /* KEYS_HPP */
