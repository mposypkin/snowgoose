#ifndef ALGPWL__HPP
#define ALGPWL__HPP

#include <iostream> 
#include <memory>
#include <cmath>

#include "algorithm.hpp"
#include "pwl/pwlbound.hpp"
#include "pwl/pwlfunc.hpp"
#include "interval/enums.h"

using namespace snowgoose::pwl;

namespace snowgoose {
namespace expression {

template<class T>
	class PwlBoundAlg : public Algorithm<PwlBound<T>>
	{
    private:
        T m_a; 
		T m_b;   
        int m_steps;    
	public:
        PwlBoundAlg(T a, T b, int steps) : m_a(a), m_b(b), m_steps(steps) {}
		PwlBound<T> Plus(const PwlBound<T>& left, const PwlBound<T>& right) const { return left + right; }
		PwlBound<T> Minus(const PwlBound<T>& left, const PwlBound<T>& right) const { return left - right; }
		PwlBound<T> Mul(const PwlBound<T>& left, const PwlBound<T>& right) const { return left * right; }
		PwlBound<T> Div(const PwlBound<T>& left, const PwlBound<T>& right) const { return left / right; }
		PwlBound<T> Sin(const PwlBound<T>& t) const { return sin(t); }
		PwlBound<T> Cos(const PwlBound<T>& t) const { return cos(t); }
		PwlBound<T> Tan(const PwlBound<T>& t) const { return tg(t);  }
		PwlBound<T> Ctg(const PwlBound<T>& t) const { return ctg(t);  }
		PwlBound<T> ArcCos(const PwlBound<T>& t) const { return acos(t); }
		PwlBound<T> ArcSin(const PwlBound<T>& t) const { return asin(t); }
		PwlBound<T> ArcTan(const PwlBound<T>& t) const { return atg(t); }
		PwlBound<T> ArcCtg(const PwlBound<T>& t) const { return actg(t); }
		PwlBound<T> Exp(const PwlBound<T>& t) const { return exp(t); }
		PwlBound<T> Sqrt(const PwlBound<T>& t) const { return sqrt(t); }
		PwlBound<T> Sqr(const PwlBound<T>& t) const { return sqr(t); }
		PwlBound<T> Pow(const PwlBound<T>& base, int exp) const { return base ^ exp; }
		PwlBound<T> Pow(const PwlBound<T>& base, bool isBaseVar, const PwlBound<T>& exp, bool isExpVar) const { 
	        throw std::invalid_argument("Exception in Pow(const PwlBound<T>& base, bool isBaseVar, const PwlBound<T>& exp, bool isExpVar). The operation is not defined.");
    	}
        PwlBound<T> PowDouble(const PwlBound<T>& base, double exp) const { return base ^ (T)exp; }
		PwlBound<T> Abs(const PwlBound<T>& t) const { return abs(t); }
		PwlBound<T> Ln(const PwlBound<T>& t) const { return ln(t); }
		PwlBound<T> Log(const PwlBound<T>& t, double base) const { return log(t, (T)base); }
		PwlBound<T> Min(const PwlBound<T>& left, const PwlBound<T>& right) const { return min(left, right); }
		PwlBound<T> Max(const PwlBound<T>& left, const PwlBound<T>& right) const { return max(left, right); }
		PwlBound<T> IfTrue(IntervalBool condition, const PwlBound<T>& left, const PwlBound<T>& right) const { return ifThen(condition, left, right);  }
		IntervalBool Condition(Conditions condition, const PwlBound<T>& left, const PwlBound<T>& right) const
		{
			switch (condition)
			{
			case Conditions::More:
				return left > right;
				break;
			case Conditions::Less:
				return left < right;
				break;
			case Conditions::LessEqual:
				return left <= right;
				break;
			case Conditions::MoreEqual:
				return left >= right;
				break;
			}
			throw std::invalid_argument("Exception in Condition(Conditions condition, const PwlBound<T>& left, const PwlBound<T>& right). Invalid condition.");
		}
		PwlBound<T> CreateVar(int index) const 
		{ 
		    return PwlBound<T>(m_a, m_b, m_steps);
		}      
		PwlBound<T> CreateConst(double cnst) const 
		{
			PwlFunc<T> lower({ { m_a, (T)cnst },{ m_b, (T)cnst } });
			PwlFunc<T> upper({ { m_a, (T)cnst },{ m_b, (T)cnst } });
			return PwlBound<T>(lower, upper, m_steps); 
		}
		vPtrAlg<PwlBound<T>> GetNewAlgorithm(Conditions cond, int index, double cnst) const
		{ 
			if(cond == Conditions::More || cond == Conditions::MoreEqual) {	
				return vPtrAlg<PwlBound<T>>({ ptrAlg<PwlBound<T>>(new PwlBoundAlg<T>((T)cnst, m_b, m_steps)), ptrAlg<PwlBound<T>>(new PwlBoundAlg<T>(m_a, (T)cnst, m_steps))});
			}
			else { // Less or LessEqual
				return vPtrAlg<PwlBound<T>>({ ptrAlg<PwlBound<T>>(new PwlBoundAlg<T>(m_a, (T)cnst, m_steps)), ptrAlg<PwlBound<T>>(new PwlBoundAlg<T>((T)cnst, m_b, m_steps))});
			}	
		}
		PwlBound<T> UnaryMinus(const PwlBound<T>& t) const { return -1.0*t; }
	};
  
}
}

#endif
