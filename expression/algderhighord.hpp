#ifndef ALGDERHIGHORD__HPP
#define ALGDERHIGHORD__HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>

#include "algorithm.hpp"
#include "derhighorder/series.hpp"
//#include "derivatives/intervalder.hpp"

//using namespace snowgoose::interval;
using namespace snowgoose::derhighorder;

namespace snowgoose {
namespace expression {

	template<class T>
	class SeriesAlg : public Algorithm<Series<T>>
	{
    	private:
        T m_v;
        int m_order;     
	public:
        	SeriesAlg(T v, int order) : m_v(v), m_order(order) {}
		Series<T> Plus(const Series<T>& left, const Series<T>& right) const { return left + right; }
		Series<T> Minus(const Series<T>& left, const Series<T>& right) const { return left - right; }
		Series<T> Mul(const Series<T>& left, const Series<T>& right) const { return left * right; }
		Series<T> Div(const Series<T>& left, const Series<T>& right) const { return left / right; }
		Series<T> Sin(const Series<T>& t) const { return sin(t); }
		Series<T> Cos(const Series<T>& t) const { return cos(t); }
		Series<T> Tan(const Series<T>& t) const { return tg(t);  }
		Series<T> Ctg(const Series<T>& t) const { return ctg(t);  }
		Series<T> ArcCos(const Series<T>& t) const { return acos(t); }
		Series<T> ArcSin(const Series<T>& t) const { return asin(t); }
		Series<T> ArcTan(const Series<T>& t) const { return atg(t); }
		Series<T> ArcCtg(const Series<T>& t) const { return actg(t); }
		Series<T> Exp(const Series<T>& t) const { return exp(t); }
		Series<T> Sqrt(const Series<T>& t) const { return sqrt(t); }
		Series<T> Sqr(const Series<T>& t) const { return sqr(t); }
		Series<T> Pow(const Series<T>& base, int exp) const { return base ^ exp; }
		Series<T> Pow(const Series<T>& base, bool isBaseVar, const Series<T>& exp, bool isExpVar) const { 
		    if(isBaseVar && isExpVar) 
		        throw std::invalid_argument("Invalid argument in DerAlg::Pow. The arguments base and exp depend on var.");
		    if(!isBaseVar)
		        return pow(base.value(), exp);
		    else
		        return base ^ exp.value();
        	};
        	Series<T> PowDouble(const Series<T>& base, double exp) const { return base ^ (T)exp; }
		Series<T> Abs(const Series<T>& t) const { return abs(t); }
		Series<T> Ln(const Series<T>& t) const { return ln(t); }
		Series<T> Log(const Series<T>& t, double base) const { return log(t, base); }
		Series<T> Min(const Series<T>& left, const Series<T>& right) const { return min(left, right); }
		Series<T> Max(const Series<T>& left, const Series<T>& right) const { return max(left, right); }
		Series<T> IfTrue(IntervalBool condition, const Series<T>& left, const Series<T>& right) const { return ifThen(condition, left, right);  }
		Series<T> IfTrue(const ptrCNode<Series<T>>& conditionNode, const ptrNode<Series<T>>& left, const ptrNode<Series<T>>& right, MapIterator &map_iterator) const 
		{
			IntervalBool condition = conditionNode->calcCondition(*this, map_iterator);
			return ifThen(condition, left->calc(*this, map_iterator), right->calc(*this, map_iterator));
		}
		IntervalBool Condition(Conditions condition, const Series<T>& left, const Series<T>& right) const { return cond(condition, left, right); }
		Series<T> CreateVar(int index) const 
		{ 
		    	return Series<T>(m_v, 1.0, m_order);
		}      
		Series<T> CreateConst(double cnst) const 
		{ 
			return Series<T>(cnst, 0.0, m_order); 
		}
	};

  
}
}

#endif
