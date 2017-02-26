/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/

/*
* File:   descfunc.hpp
* Author: alusov
*
* Created on February 23, 2017, 5:59 PM
*/

#ifndef DESCFUNC_HPP
#define DESCFUNC_HPP

#include <string>
#include <vector>
#include <utility> 
#include <map>
#include "keys.hpp"

namespace OPTITEST
{

	class X
	{
	public:
		X(const std::vector<double> &point, bool anyDim) : globPoint(point) {}
		X(){}

		double operator[](int i)
		{
			if (anyDim)
				return globPoint[0];
			else
				return globPoint[i];
		}
	private:
		std::vector<double> globPoint;
		bool anyDim;
	};
	
	class Bound
	{
	public:
		Bound(std::pair<double,double> bound) : bound(bound) {}
		Bound(const std::vector<std::pair<double,double>> &bounds) : bounds(bounds) {}
		Bound(){}

	std::pair<double,double> operator[](int i)
	{
		if(bounds.size()==0)
			return bound;
		else
			return bounds[i];
	}
	private:
		std::vector<std::pair<double,double>> bounds;
		std::pair<double,double> bound;
	};
	
	class DescFunc
	{
	public:
		DescFunc(const std::string &desc, int dim, const std::vector<std::pair<double, double>> &vBound, const std::vector<double> &x, double y ) :
			descrirtion(desc), dimension(dim), bound(vBound), x(x, false), y(y)  {}

		DescFunc(const std::string &desc, const std::pair<double, double> &bound, double x, double y) :
			descrirtion(desc), dimension(0), bound(bound), x({ x }, true), y(y) {}

		DescFunc() {}

		std::string descrirtion;
		int dimension;
		Bound bound;
		X x;
		double y;	
	};

	class DescFuncs
	{
	public:
		DescFuncs() 
		{
			dfs[Key::Ackley1] = DescFunc("Ackley 1", { -35, 35 }, 0.0, 0.0 );
			dfs[Key::Ackley2] = DescFunc("Ackley 2", 2, { { -35, 35 }, { -35, 35 } }, { 0.0, 0.0 }, 0.0 );
		}
		DescFunc operator[](Key key)
		{
			if(dfs.find(key)== dfs.end())
				throw std::invalid_argument("Invalid Key.");
			return dfs[key];
		}
	private:
		std::map<Key, DescFunc> dfs;
	};
}

#endif /* DESCFUNC_HPP */
