#ifndef ALGDER__HPP
#define ALGDER__HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>

#include "interval/interval_air.hpp"
#include "interval/enums.h"
#include "algorithm.hpp"
#include "derivatives/valder.hpp"
#include "derivatives/grad.hpp"


using namespace snowgoose::interval;
using namespace snowgoose::derivative;

namespace snowgoose {
namespace expression {

	template<class T>
	class DerAlg : public Algorithm<ValDer<T>>
	{
    private:
        std::vector<T> m_v;
	public:
        DerAlg(const std::vector<T> &v) : m_v(v) {}
		ValDer<T> Plus(const ValDer<T>& left, const ValDer<T>& right) const { return left + right; }
		ValDer<T> Minus(const ValDer<T>& left, const ValDer<T>& right) const { return left - right; }
		ValDer<T> Mul(const ValDer<T>& left, const ValDer<T>& right) const { return left * right; }
		ValDer<T> Div(const ValDer<T>& left, const ValDer<T>& right) const { return left / right; }
		ValDer<T> Sin(const ValDer<T>& t) const { return sin(t); };
		ValDer<T> Cos(const ValDer<T>& t) const { return cos(t); };
		ValDer<T> Tan(const ValDer<T>& t) const { return tg(t);  };
		ValDer<T> Ctg(const ValDer<T>& t) const { throw "Not implemented yet";  };
		ValDer<T> ArcCos(const ValDer<T>& t) const { return acos(t); };
		ValDer<T> ArcSin(const ValDer<T>& t) const { return asin(t); };
		ValDer<T> ArcTan(const ValDer<T>& t) const { throw "Not implemented yet"; };
		ValDer<T> ArcCtg(const ValDer<T>& t) const { throw "Not implemented yet"; };
		ValDer<T> Exp(const ValDer<T>& t) const { return exp(t); };
		ValDer<T> Sqrt(const ValDer<T>& t) const { return sqrt(t); };
		ValDer<T> Sqr(const ValDer<T>& t) const { return sqr(t); };
		ValDer<T> Pow(const ValDer<T>& base, int exp) const { return base ^ exp; };
		ValDer<T> Pow(const ValDer<T>& base, const ValDer<T>& exp) const { throw "Not implemented yet"; };
		ValDer<T> Abs(const ValDer<T>& t) const { throw "Not implemented yet";  };
		ValDer<T> Ln(const ValDer<T>& t) const { return ln(t); };
		ValDer<T> Log(const ValDer<T>& t, double base) const { throw "Not implemented yet"; };
		ValDer<T> Min(const ValDer<T>& left, const ValDer<T>& right) const { throw "Not implemented yet"; }
		ValDer<T> Max(const ValDer<T>& left, const ValDer<T>& right) const { throw "Not implemented yet"; };
		ValDer<T> IfTrue(IntervalBool condition, const ValDer<T>& left, const ValDer<T>& right) const { throw "Not implemented yet";  }
		IntervalBool Condition(Conditions condition, const ValDer<T>& left, const ValDer<T>& right) const { throw "Not implemented yet"; }
        ValDer<T> CreateVar(int index) const 
        { 
            std::vector<T> grad(m_v.size(), 0.0);
            grad[index] = 1.0;
            return ValDer<T>(m_v[index], Grad<T>(grad)); 
        }      
        ValDer<T> CreateConst(double cnst) const { return ValDer<T>(cnst, Grad<T>(m_v.size(), 0.0)); };
	};
    
    template<class T=double>
	class ValDerAlg : public DerAlg<T>
	{
	public:
        ValDerAlg(const std::vector<T> &v) : DerAlg<T>(v) {}
    };
    
    template<class T=double>
	class IntervalDerAlg : public DerAlg<Interval<T>>
	{
	public:
        IntervalDerAlg(const std::vector<Interval<T>> &v) : DerAlg<Interval<T>>(v) {}
    };
}
}

#endif /* ALGORITHM__HPP */
