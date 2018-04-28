#ifndef ALGSYMDIF_HPP
#define ALGSYMDIF_HPP

#include <iostream> 
#include <memory>
#include <cmath>
#include <vector>
#include "node.hpp"
#include "expr.hpp"
#include "algorithm.hpp"
#include "interval/interval_air.hpp"


using namespace snowgoose::interval;

namespace snowgoose {
namespace expression {


template<class T> 
	class DiffAlg : public virtual Algorithm<T>
	{
	private:
		int m_order;
	public:
		DiffAlg(int order) : m_order(order) {}
                int order() const { return m_order; }
	};

template<class T> 
	class DiffPointAlg : public FuncAlg<T>, public DiffAlg<T>
	{
	public:
		DiffPointAlg(T point, int order) : FuncAlg<T>(std::vector<T>(1, point)), DiffAlg<T>(order) {}
	};

template<class T> class DiffIntervalAlg : public InterEvalAlg<T>, public DiffAlg<Interval<T>>
{
	public:
		DiffIntervalAlg(const Interval<T> &interval, int order) : InterEvalAlg<T>(std::vector<Interval<T>>(1, interval)), DiffAlg<Interval<T>>(order) {}
};

template<class T>
	class SymDifAlg : public HandleNodeAlg<T>
	{    
	public:
        	SymDifAlg(){}
		Expr<T> varNode(int index) const
                {
                	return 1.0;
                }
                Expr<T> constNode(double val) const
                {
			return 0.0;
                }
                Expr<T> plusNode(const Expr<T> &l, const Expr<T> &r) const
                {
			return l.handle(*this) + r.handle(*this);
                }
                Expr<T> minusNode(const Expr<T> &l, const Expr<T> &r) const
                {
			return l.handle(*this) - r.handle(*this);
                }
                Expr<T> mulNode(const Expr<T> &l, const Expr<T> &r) const
                {
			return l.handle(*this) * r.copy() + l.copy() * r.handle(*this);
                }
                Expr<T> divNode(const Expr<T> &l, const Expr<T> &r) const
                {
                	return (l.handle(*this) * r.copy() - l.copy() * r.handle(*this))/sqr(r.copy());
                }
                Expr<T> sinNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return cos(expr.copy());
			else
				return cos(expr.copy()) * dexpr;
                }
                Expr<T> cosNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return -sin(expr.copy());
			else
				return -sin(expr.copy()) * dexpr;			
                }
 		Expr<T> tgNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return 1.0 + sqr(tg(expr.copy()));
			else
				return (1.0 + sqr(tg(expr.copy()))) * dexpr;			
                }
 		Expr<T> ctgNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return -(1.0 + sqr(ctg(expr.copy())));
			else
				return -(1.0 + sqr(ctg(expr.copy()))) * dexpr;			
                }
                Expr<T> asinNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return 1.0/sqrt(1.0-sqr(expr.copy()));
			else
				return dexpr/sqrt(1.0-sqr(expr.copy()));
                }
                Expr<T> acosNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return -1.0/sqrt(1.0-sqr(expr.copy()));
			else
				return -dexpr/sqrt(1.0-sqr(expr.copy()));			
                }
 		Expr<T> atgNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return 1.0/(1.0 + sqr(expr.copy()));
			else
				return dexpr/(1.0 + sqr(expr.copy()));			
                }
 		Expr<T> actgNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return -1.0/(1.0 + sqr(expr.copy()));
			else
				return -dexpr/(1.0 + sqr(expr.copy()));			
                }

		Expr<T> expNode(const Expr<T> &expr) const
		{		                
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return exp(expr.copy());
			else
				return exp(expr.copy()) * dexpr;
		}
		Expr<T> sqrtNode(const Expr<T> &expr) const
		{		                
			return expr.handle(*this)/(2.0 * sqrt(expr.copy()));
		}

                Expr<T> sqrNode(const Expr<T> &expr) const
                {
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return 2.0 * expr.copy();
			else
				return 2.0 * expr.copy() * dexpr;		
                }

		Expr<T> powIntNode(const Expr<T> &expr, int exponent) const
		{
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return (double)exponent * (expr.copy()^(exponent-1));
			else
				return (double)exponent * (expr.copy()^(exponent-1)) * dexpr;
		}

		Expr<T> powNode(const Expr<T> &expr, double exponent) const
		{
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return exponent * (expr.copy()^(exponent-1.0));
			else
				return exponent * (expr.copy()^(exponent-1.0)) * dexpr;
		}

		Expr<T> powExprNode(const Expr<T> &l, const Expr<T> &r) const
		{
			throw std::invalid_argument("SymDifAlg lib. PowExprNode is not supported function.");
		}

		Expr<T> absNode(const Expr<T> &expr) const
		{
			Expr<T> dexpr = expr.handle(*this);
			if(dexpr == 1.0)
				return abs(expr.copy())/expr.copy();
			else
				return dexpr * abs(expr.copy())/expr.copy();
		}

		Expr<T> unaryMinusNode(const Expr<T> &expr) const
		{
			return -expr.handle(*this);
		}

		Expr<T> lnNode(const Expr<T> &expr) const
		{
			return expr.handle(*this)/expr.copy();
		}

		Expr<T> logNode(const Expr<T> &expr, double base) const
		{
			return expr.handle(*this)/(std::log(base) * expr.copy());
		}

		Expr<T> minNode(const Expr<T> &l, const Expr<T> &r) const
		{
			return ifThen(l.copy() < r.copy(), l.handle(*this), r.handle(*this)); 
		}

		Expr<T> maxNode(const Expr<T> &l, const Expr<T> &r) const
		{
			return ifThen(l.copy() > r.copy(), l.handle(*this), r.handle(*this)); 
		}

		Expr<T> ifTrueNode(const Expr<T> &cond, const Expr<T> &l, const Expr<T> &r) const
		{
			return ifThen(cond.copy(), l.handle(*this), r.handle(*this)); 
		}
		Expr<T> iteratorNode() const
		{
			return 0.0;
		}
		Expr<T> indexNode() const
		{
			return 1.0;
		}
		Expr<T> calcIndexNode() const
		{
			return 1.0;
		}
		Expr<T> exprIndexNode() const
		{
			return 1.0;
		}


	};

}
}

#endif /* ALGORITHM__HPP */
