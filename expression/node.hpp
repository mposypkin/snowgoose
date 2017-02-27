#ifndef NODE__HPP
#define NODE__HPP

#include <iostream> 
#include <memory>
#include <string>
#include <vector>
#include "algorithm.hpp"
#include "interval/interval_air.hpp"
#include "interval/enums.h"
#include "expr.hpp"


namespace snowgoose {
namespace expression {

	template<class T> class Node;
	template<class T> class ConditionNode;
	template<typename T> using ptrNode = std::shared_ptr<Node<T>>;
	template<typename T> using ptrCNode = std::shared_ptr<ConditionNode<T>>;
	template<typename T> using vPtrNode = std::vector<ptrNode<T>>;

	template<class T>
	class Node
	{
	public:
		Node(const vPtrNode<T> &childs) : m_childs(childs) {}
		Node() {}
		virtual T calc(const std::vector<T> &v, const Algorithm<T> &) = 0;
		friend std::ostream& operator<<(std::ostream & out, const Node<T>& v) { return v.prn(out);}
		virtual std::ostream& prn(std::ostream & out) const = 0;
	protected:
		vPtrNode<T> m_childs;
	};	

	template<class T>
	class Var : public Node<T>
	{
	private:
		int index;
	public:
		Var(int i) : index(i) {}
		T calc(const std::vector<T> &v, const Algorithm<T> &alg) { return v[index]; };
		std::ostream& prn(std::ostream & out) const { return out << "x[" << index << "]"; };
	};


	template <class T>
	class Const : public Node<T>
	{
	private:
		T m_const;
	public:
		Const(double value) : m_const(value) {}
		T calc(const std::vector<T> &v, const Algorithm<T> &) { return m_const; }
		std::ostream& prn(std::ostream & out) const { return out << m_const; }
	};

	template <class T>
	class Plus : public Node<T>
	{
	public:
		Plus(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Plus(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << " + " << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Minus : public Node<T>
	{
	public:
		Minus(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Minus(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); }
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << " - " << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Mul : public Node<T>
	{
	public:
		Mul(const ptrNode<T> &left,const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Mul(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); }
		std::ostream& prn(std::ostream & out) const { return out << *this->m_childs[0] << " * " << *this->m_childs[1]; }
	};

	template <class T>
	class Div : public Node<T>
	{
	public:
		Div(const ptrNode<T> &left, const ptrNode<T> &right) : Node<T>({ left, right }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Div(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "(" << *this->m_childs[0] << "/" << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Sin : public Node<T>
	{
	public:
		Sin(const ptrNode<T> &node) : Node<T>({node}) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Sin(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "sin(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Cos : public Node<T>
	{
	public:
		Cos(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Cos(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "cos(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Tg : public Node<T>
	{
	public:
		Tg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Tan(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "tg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Ctg : public Node<T>
	{
	public:
		Ctg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Ctg(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "ctg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcCos : public Node<T>
	{
	public:
		ArcCos(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.ArcCos(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "acos(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcSin : public Node<T>
	{
	public:
		ArcSin(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.ArcSin(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "asin(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcTg : public Node<T>
	{
	public:
		ArcTg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.ArcTan(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "atg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class ArcCtg : public Node<T>
	{
	public:
		ArcCtg(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.ArcCtg(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "actg(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Exp : public Node<T>
	{
	public:
		Exp(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Exp(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "exp(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Sqrt : public Node<T>
	{
	public:
		Sqrt(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Sqrt(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "sqrt(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Sqr : public Node<T>
	{
	public:
		Sqr(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Sqr(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "sqr(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class PowInt : public Node<T>
	{
	private:
		int exponent;
	public:
		PowInt(const ptrNode<T> &node, int exp) : Node<T>({ node }), exponent(exp) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { 
			return alg.Pow(this->m_childs[0]->calc(v, alg), exponent); 
		};
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << exponent << ")"; }
	};

	template <class T>
	class Pow : public Node<T>
	{
	private:
		T exponent;
	public:
		Pow(const ptrNode<T> &node, double exp) : Node<T>({ node }), exponent(exp) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) {
			return alg.Pow(this->m_childs[0]->calc(v, alg), exponent);
		};
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << exponent << ")"; }
	};

	template <class T>
	class PowExpr : public Node<T>
	{
	public:
		PowExpr(const ptrNode<T> &base, const ptrNode<T> &exp) : Node<T>({ base, exp }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Pow(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "pow(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Abs : public Node<T>
	{
	public:
		Abs(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Abs(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "abs(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Ln : public Node<T>
	{
	public:
		Ln(const ptrNode<T> &node) : Node<T>({ node }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Ln(this->m_childs[0]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "ln(" << *this->m_childs[0] << ")"; }
	};

	template <class T>
	class Log : public Node<T>
	{
	private:
		double base;
	public:
		Log(const ptrNode<T> &node, double b) : Node<T>({ node }), base(b) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Log(this->m_childs[0]->calc(v, alg), base); };
		std::ostream& prn(std::ostream & out) const { return out << "log(" << *this->m_childs[0] << "," << base << ")"; }
	};

	template <class T>
	class Min : public Node<T>
	{
	public:
		Min(const ptrNode<T> &lv, const ptrNode<T> &rv) : Node<T>({ lv,  rv }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Min(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "min(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class Max : public Node<T>
	{
	public:
		Max(const ptrNode<T> &lv, const ptrNode<T> &rv) : Node<T>({ lv,  rv }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) { return alg.Max(this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg)); };
		std::ostream& prn(std::ostream & out) const { return out << "max(" << *this->m_childs[0] << "," << *this->m_childs[1] << ")"; }
	};

	template <class T>
	class ConditionNode : public Node<T>
	{
	private:
		Conditions condition;
		IntervalBool result;
	public:
		ConditionNode(Conditions cond, const ptrNode<T> &lv, const ptrNode<T> &rv) : condition(cond), Node<T>({ lv,  rv }) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) 
		{ 
			throw "Invalid operation in ConditionNode.";
		}
		IntervalBool calcCondition(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			return alg.Condition(condition, this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg));
		}
		std::ostream& prn(std::ostream & out) const { return out << " " << *this->m_childs[0] << condition << *this->m_childs[1] << " "; }
	};

	template <class T>
	class IfTrue : public Node<T>
	{
	private:
		ptrCNode<T> conditionNode;
	public:
		IfTrue(const ptrCNode<T> &cond, const ptrNode<T> &lv, const ptrNode<T> &rv) : conditionNode(cond), Node<T>({ lv,  rv }) {}

		T calc(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			return alg.IfTrue(conditionNode->calcCondition(v, alg), this->m_childs[0]->calc(v, alg), this->m_childs[1]->calc(v, alg));
		}
		std::ostream& prn(std::ostream & out) const { return out << "ifThen(" << *conditionNode << " , " << *this->m_childs[0] << " , " << *this->m_childs[1] << ")"; }
	};

	template <class T> class Expr;
	template <class T> class IteratorNode;

	/**
	* Iterator is used  for loopSum, loopMul expressions.
	*/
	template <class T> class IteratorNode;
	class Iterator
	{
	public:
		Iterator(int start, int end) : startIterator(new int(start)), endIterator(new int(end)), current(new int(start))
		{
		}
		int Start() const { return *startIterator; }
		int End() const { return *endIterator; }
		int Current() const { return *current; }
		bool CanIterate() const { return *current <= *endIterator; }
		void Next() { if (CanIterate()) (*current)++;}
		void Reset() { *current = *startIterator; }
		template <class T>
		/**
		* Allows to conver iterator to expression
		* @return expression
		*/
		operator Expr<T>() const
		{
			ptrNode<T> pNode(new IteratorNode<T>(*this));
			return Expr<T>(pNode);
		}
	private:
		std::shared_ptr<int> startIterator;
		std::shared_ptr<int> endIterator;
		std::shared_ptr<int> current;
	};

	template <class T>
	class IteratorNode : public Node<T>
	{
		Iterator iterator;
	public:
		IteratorNode(const Iterator &i) : iterator(i) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			return iterator.Current();
		}
		std::ostream& prn(std::ostream & out) const { return out << "i"; }
	};

	template <class T>
	class Index : public Node<T>
	{
		Iterator iterator;
	public:
		Index(const Iterator &i) : iterator(i) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			int index = iterator.Current();
			return v[index];
		}
		std::ostream& prn(std::ostream & out) const { return out << "x[i]"; };
	};

	template <class T>
	class CalcIndex : public Node<T>
	{
		std::function<int()> func;
	public:
		CalcIndex(const std::function<int()> &f) : func(f) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			return v[func()];
		}
		std::ostream& prn(std::ostream & out) const { return out << "x[calc i]"; };
	};

	template <class T>
	class ExprIndex : public Node<T>
	{
		Expr<T> expr;
	public:
		ExprIndex(const Expr<T>& e) : expr(e) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg)
		{
			size_t index = (size_t)((double)expr.calc(v, alg));
			return v[index];
		}
		std::ostream& prn(std::ostream & out) const { return out << "x[" << expr << "]"; };
	};

	template <class T>
	class CycleSum : public Node<T>
	{
	private:
		Iterator iterator;
	public:
		CycleSum(const ptrNode<T> &node, const Iterator &i) : Node<T>({ node }), iterator(i) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) {
			T result = T();
			iterator.Reset();
			while (iterator.CanIterate())
			{
				result = alg.Plus(result, this->m_childs[0]->calc(v, alg));
				iterator.Next();
			}
			return result;
		};
		std::ostream& prn(std::ostream & out) const { return out << "loopSum(" << *this->m_childs[0] << ",i)"; }
	};

	template <class T>
	class CycleMul : public Node<T>
	{
	private:
		Iterator iterator;
	public:
		CycleMul(const ptrNode<T> &node, const Iterator &i) : Node<T>({ node }), iterator(i) {}
		T calc(const std::vector<T> &v, const Algorithm<T> & alg) {
			T result = 1.0;
			iterator.Reset();
			if (!iterator.CanIterate())
				return T();
			while (iterator.CanIterate())
			{
				result = alg.Mul(result, this->m_childs[0]->calc(v, alg));
				iterator.Next();
			}
			return result;
		};
		std::ostream& prn(std::ostream & out) const { return out << "loopMul(" << *this->m_childs[0] << ",i)"; }
	};

}}

#endif /* NODE__HPP */

