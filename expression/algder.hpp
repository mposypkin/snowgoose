#ifndef ALGDER__HPP
#define ALGDER__HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>

#include "algorithm.hpp"
#include "derivatives/valder.hpp"
#include "derivatives/intervalder.hpp"



using namespace snowgoose::interval;
using namespace snowgoose::derivative;

namespace snowgoose {
namespace expression {

	template<class T>
	class DerAlg : public Algorithm<T>
	{
    private:
        using T2 = typename T::value_type;//double or Interval<double>
        std::vector<T2> m_v;      
	public:
        DerAlg(const std::vector<T2> &v) : m_v(v) {}
		T Plus(const T& left, const T& right) const { return left + right; }
		T Minus(const T& left, const T& right) const { return left - right; }
		T Mul(const T& left, const T& right) const { return left * right; }
		T Div(const T& left, const T& right) const { return left / right; }
		T Sin(const T& t) const { return sin(t); }
		T Cos(const T& t) const { return cos(t); }
		T Tan(const T& t) const { return tg(t);  }
		T Ctg(const T& t) const { return ctg(t);  }
		T ArcCos(const T& t) const { return acos(t); }
		T ArcSin(const T& t) const { return asin(t); }
		T ArcTan(const T& t) const { return atg(t); }
		T ArcCtg(const T& t) const { return actg(t); }
		T Exp(const T& t) const { return exp(t); }
		T Sqrt(const T& t) const { return sqrt(t); }
		T Sqr(const T& t) const { return sqr(t); }
		T Pow(const T& base, int exp) const { return base ^ exp; }
		T Pow(const T& base, bool isBaseVar, const T& exp, bool isExpVar) const { 
            if(isBaseVar && isExpVar) 
                throw std::invalid_argument("Invalid argument in DerAlg::Pow. The arguments base and exp depend on var.");
            if(!isBaseVar)
                return pow(base.value(), exp);
            else
                return base ^ exp.value();
        };
        T PowDouble(const T& base, double exp) const { return base ^ (typename T::value_type)exp; }
		T Abs(const T& t) const { return abs(t); }
		T Ln(const T& t) const { return ln(t); }
		T Log(const T& t, double base) const { return log(t, base); }
		T Min(const T& left, const T& right) const { return min(left, right); }
		T Max(const T& left, const T& right) const { return max(left, right); }
		T IfTrue(IntervalBool condition, const T& left, const T& right) const { return ifThen(condition, left, right);  }
		IntervalBool Condition(Conditions condition, const T& left, const T& right) const { return cond(condition, left, right); }
        T CreateVar(int index) const 
        { 
            std::vector<T2> grad(m_v.size(), 0.0);
            grad[index] = 1.0;
            return T(m_v[index], Grad<T2>(grad));
        }      
        T CreateConst(double cnst) const { return T(cnst, Grad<T2>(m_v.size(), 0.0)); }
	};
    
    template<class T=double>
    using ValDerAlg = DerAlg<ValDer<T>>;
    
    template<class T=double>
    using IntervalDerAlg = DerAlg<IntervalDer<T>>;
    
}
}

#endif /* ALGDER__HPP */
