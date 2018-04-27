#ifndef ALGDERHIGHORD__HPP
#define ALGDERHIGHORD__HPP

#include <iostream> 
#include <memory>
#include <cmath>

#include "algorithm.hpp"
#include "derhighorder/series.hpp"
#include "derhighorder/intervalseries.hpp"

using namespace snowgoose::interval;
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
		        throw std::invalid_argument("Invalid argument in SeriesAlg::Pow. The arguments base and exp depend on var.");
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
		Series<T> UnaryMinus(const Series<T>& t) const { return -1.0 * t; }
	};

template<class T>
	class IntervalSeriesAlg : public Algorithm<IntervalSeries<T>>
	{
    private:
        Interval<T> m_v;   
        int m_order;    
	public:
        	IntervalSeriesAlg(const Interval<T> &v, int order) : m_v(v), m_order(order) {}
		IntervalSeries<T> Plus(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return left + right; }
		IntervalSeries<T> Minus(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return left - right; }
		IntervalSeries<T> Mul(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return left * right; }
		IntervalSeries<T> Div(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return left / right; }
		IntervalSeries<T> Sin(const IntervalSeries<T>& t) const { return sin(t); }
		IntervalSeries<T> Cos(const IntervalSeries<T>& t) const { return cos(t); }
		IntervalSeries<T> Tan(const IntervalSeries<T>& t) const { return tg(t);  }
		IntervalSeries<T> Ctg(const IntervalSeries<T>& t) const { return ctg(t);  }
		IntervalSeries<T> ArcCos(const IntervalSeries<T>& t) const { return acos(t); }
		IntervalSeries<T> ArcSin(const IntervalSeries<T>& t) const { return asin(t); }
		IntervalSeries<T> ArcTan(const IntervalSeries<T>& t) const { return atg(t); }
		IntervalSeries<T> ArcCtg(const IntervalSeries<T>& t) const { return actg(t); }
		IntervalSeries<T> Exp(const IntervalSeries<T>& t) const { return exp(t); }
		IntervalSeries<T> Sqrt(const IntervalSeries<T>& t) const { return sqrt(t); }
		IntervalSeries<T> Sqr(const IntervalSeries<T>& t) const { return sqr(t); }
		IntervalSeries<T> Pow(const IntervalSeries<T>& base, int exp) const { return base ^ exp; }
		IntervalSeries<T> Pow(const IntervalSeries<T>& base, bool isBaseVar, const IntervalSeries<T>& exp, bool isExpVar) const { 
		    if(isBaseVar && isExpVar) 
		        throw std::invalid_argument("Invalid argument in IntervalSeriesAlg::Pow. The arguments base and exp depend on var.");
		    if(!isBaseVar)
		        return pow(base.value(), exp);
		    else
		        return base ^ exp.value();
        	}
        	IntervalSeries<T> PowDouble(const IntervalSeries<T>& base, double exp) const { return base ^ Interval<T>(exp); }
		IntervalSeries<T> Abs(const IntervalSeries<T>& t) const { return abs(t); }
		IntervalSeries<T> Ln(const IntervalSeries<T>& t) const { return ln(t); }
		IntervalSeries<T> Log(const IntervalSeries<T>& t, double base) const { return log(t, base); }
		IntervalSeries<T> Min(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return min(left, right); }
		IntervalSeries<T> Max(const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return max(left, right); }
		IntervalSeries<T> IfTrue(IntervalBool condition, const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return ifThen(condition, left, right);  }
		IntervalSeries<T> IfTrue(const ptrCNode<IntervalSeries<T>>& conditionNode, const ptrNode<IntervalSeries<T>>& left, const ptrNode<IntervalSeries<T>>& right, MapIterator &map_iterator) const 
		{
			IntervalBool condition = conditionNode->calcCondition(*this, map_iterator);
			auto var = std::dynamic_pointer_cast<Var<IntervalSeries<T>>>(left);
			auto cnst = std::dynamic_pointer_cast<Const<IntervalSeries<T>>>(right);
			if(condition == IntervalBool::Intermadiate && var && cnst)
			{
				Conditions cond = conditionNode->getCondition();
				Interval<T> v(m_v);
				if(cond == Conditions::More || cond == Conditions::MoreEqual)					
					v = Interval<T>(cnst->getConst(), v.rb());
				else
					v = Interval<T>(v.lb(), cnst->getConst());
				IntervalSeriesAlg<T> alg(v, m_order);
				return ifThen(condition, left->calc(alg, map_iterator), right->calc(alg, map_iterator));
			}
			else
				return ifThen(condition, left->calc(*this, map_iterator), right->calc(*this, map_iterator));
		}
		IntervalBool Condition(Conditions condition, const IntervalSeries<T>& left, const IntervalSeries<T>& right) const { return cond(condition, left, right); }
		IntervalSeries<T> CreateVar(int index) const 
		{ 
		    return IntervalSeries<T>(m_v, 1.0, m_order);
		}      
		IntervalSeries<T> CreateConst(double cnst) const { return IntervalSeries<T>(cnst, 0.0, m_order); }
		vPtrAlg<IntervalSeries<T>> GetNewAlgorithm(Conditions cond, int index, double cnst) const
		{ 
			Interval<T> leftInterval(m_v);
			Interval<T> rightInterval(m_v);
			if(cond == Conditions::More || cond == Conditions::MoreEqual) {	
				leftInterval = Interval<T>(cnst, m_v.rb());				
				rightInterval = Interval<T>(m_v.lb(), cnst);
			}
			else { // Less or LessEqual
				leftInterval = Interval<T>(m_v.lb(), cnst);				
				rightInterval = Interval<T>(cnst, m_v.rb());
			}
			return vPtrAlg<IntervalSeries<T>>({ ptrAlg<IntervalSeries<T>>(new IntervalSeriesAlg<T>(leftInterval, m_order)), ptrAlg<IntervalSeries<T>>(new IntervalSeriesAlg<T>(rightInterval, m_order))});
		}
		IntervalSeries<T> UnaryMinus(const IntervalSeries<T>& t) const { return -1.0*t; }
	};
  
}
}

#endif
