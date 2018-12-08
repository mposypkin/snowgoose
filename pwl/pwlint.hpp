#ifndef PWLINT_HPP
#define PWLINT_HPP

#include <iostream>
#include "pwlbound.hpp"
#include "interval/interval_air.hpp"
#include "interval/enums.h"
#include "pwlfunc.hpp"


using snowgoose::interval::Interval;

namespace snowgoose {
	namespace pwl {

		template<class T> class PwlInterval
		{
		public:
			PwlInterval(T a, T b, int steps = 5) : m_bound(a,b,steps), m_interval(a,b) {
			}
			PwlInterval operator+(const PwlInterval& bound) const;
			PwlInterval operator-(const PwlInterval& bound) const; 
			PwlInterval operator*(const PwlInterval& bound) const;
			PwlInterval operator/(const PwlInterval& bound) const;
			PwlInterval operator^(int exp) const;
			PwlInterval operator^(T exp) const;

			IntervalBool operator < (const PwlInterval &y) const { return this->m_bound < y.m_bound; }
			IntervalBool operator <= (const PwlInterval &y) const { return this->m_bound <= y.m_bound; }
			IntervalBool operator > (const PwlInterval &y) const { return this->m_bound > y.m_bound; }
			IntervalBool operator >= (const PwlInterval &y) const { return this->m_bound >= y.m_bound;  }
			bool operator == (const PwlInterval &y) const { return this->m_bound == y.m_bound; }

			std::vector<Point<T>> lower_bound() const { return m_bound.m_lower.Points(); }
			std::vector<Point<T>> upper_bound() const { return m_bound.m_upper.Points(); }
			T get_min() const { return m_bound.m_lower.get_min(); }
			T get_max() const { return m_bound.m_upper.get_max(); }
			T get_a() const { return m_bound.m_lower.get_a(); }
			T get_b() const { return m_bound.m_lower.get_b(); }

			template<class T2> friend PwlInterval<T2> operator*(T2 x, const PwlInterval<T2>& y);
			template<class T2> friend PwlInterval<T2> operator*(const PwlInterval<T2>& x, T2 y) { return y*x; }
			template<class T2> friend PwlInterval<T2> operator/(T2 x, const PwlInterval<T2>& y);
			template<class T2> friend PwlInterval<T2> operator/(const PwlInterval<T2>& x, T2 y);
			template<class T2> friend PwlInterval<T2> operator+(T2 x, const PwlInterval<T2>& y);
			template<class T2> friend PwlInterval<T2> operator+(const PwlInterval<T2>& x, T2 y);
			template<class T2> friend PwlInterval<T2> operator-(T2 x, const PwlInterval<T2>& y);
			template<class T2> friend PwlInterval<T2> operator-(const PwlInterval<T2>& x, T2 y);			

			template<class T2> friend PwlInterval<T2> sqr(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> sqrt(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> sin(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> cos(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> asin(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> acos(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> tg(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> ctg(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> atg(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> actg(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> log(const PwlInterval<T2>& bound, T2 base);
			template<class T2> friend PwlInterval<T2> ln(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> exp(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> abs(const PwlInterval<T2>& bound);
			template<class T2> friend PwlInterval<T2> max(const PwlInterval<T2>& bound1, const PwlInterval<T2>& bound2);
			template<class T2> friend PwlInterval<T2> min(const PwlInterval<T2>& bound1, const PwlInterval<T2>& bound2);
			template<class T2> friend PwlInterval<T2> ifThen(IntervalBool ib, const PwlInterval<T2> &x, const PwlInterval<T2> &y);

			template<class T2> friend std::ostream& operator<<(std::ostream & out, const PwlInterval<T2> x);
			template<class T2> friend PwlInterval<T2> get_PwlInterval(const PwlBound<T2>& bound, const Interval<T2>& interval);
			PwlInterval(const PwlBound<T>& bound, const Interval<T>& interval) : m_bound(bound), m_interval(interval) {}
		private:
			
			PwlBound<T> m_bound;
			Interval<T> m_interval;
		};

		template<class T2> PwlInterval<T2> get_PwlInterval(const PwlBound<T2>& bound, const Interval<T2>& interval) {
			T2 pmax = bound.get_max();
			T2 pmin = bound.get_min();
			T2 imax = interval.rb();
			T2 imin = interval.lb();
			Interval<T2> new_interval(std::max(pmin, imin), std::min(pmax, imax));
			PwlFunc<T2> pwl_lower = mo(imin, pmin) ? max(bound.m_lower, PwlFunc<T2>({ { bound.get_a(),  imin },{ bound.get_b(),  imin } })) : bound.m_lower;
			PwlFunc<T2> pwl_upper = ls(imax, pmax) ? min(bound.m_upper, PwlFunc<T2>({ { bound.get_a(),  imax },{ bound.get_b(),  imax } })) : bound.m_upper;
			PwlBound<T2> newbound(pwl_lower, pwl_upper, bound.m_steps);
			return PwlInterval<T2>(newbound, new_interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator+(const PwlInterval<T>& x) const {
			auto bound = m_bound + x.m_bound;
			auto interval = m_interval + x.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator-(const PwlInterval<T>& x) const {
			auto bound = m_bound - x.m_bound;
			auto interval = m_interval - x.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator*(const PwlInterval& x) const {
			auto bound = m_bound * x.m_bound;
			auto interval = m_interval * x.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator/(const PwlInterval& x) const {
			auto bound = m_bound / x.m_bound;
			auto interval = m_interval / x.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator^(int exp) const {
			auto bound = m_bound ^ exp;
			auto interval = m_interval ^ exp;
			return get_PwlInterval(bound, interval);
		}

		template<class T> PwlInterval<T> PwlInterval<T>::operator^(T exp) const {
			auto bound = m_bound ^ exp;
			auto interval = m_interval ^ exp;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator*(T2 x, const PwlInterval<T2>& y) {
			PwlBound<T2> bound = x * y.m_bound;
			Interval<T2> interval = x * y.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator/(T2 x, const PwlInterval<T2>& y) {
			auto bound = x / y.m_bound;
			auto interval = x / y.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator/(const PwlInterval<T2>& x, T2 y) {
			auto bound = x.m_bound / y;
			auto interval = x.m_interval / y;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator+(T2 x, const PwlInterval<T2>& y) {
			auto bound = x + y.m_bound;
			auto interval = x + y.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator+(const PwlInterval<T2>& x, T2 y) {
			auto bound = x.m_bound + y;
			auto interval = x.m_interval + y;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator-(T2 x, const PwlInterval<T2>& y) {
			auto bound = x - y.m_bound;
			auto interval = x - y.m_interval;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> operator-(const PwlInterval<T2>& x, T2 y) {
			auto bound = x.m_bound - y;
			auto interval = x.m_interval - y;
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> sqr(const PwlInterval<T2>& x) {
			auto bound = sqr(x.m_bound);
			auto interval = sqr(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> sqrt(const PwlInterval<T2>& x)
		{
			auto bound = sqrt(x.m_bound);
			auto interval = sqrt(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> asin(const PwlInterval<T2>& x) {
			auto bound = asin(x.m_bound);
			auto interval = asin(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> acos(const PwlInterval<T2>& x) {
			auto bound = acos(x.m_bound);
			auto interval = acos(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> tg(const PwlInterval<T2>& x) {
			auto bound = tg(x.m_bound);
			auto interval = tg(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> ctg(const PwlInterval<T2>& x) {
			auto bound = ctg(x.m_bound);
			auto interval = ctg(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> atg(const PwlInterval<T2>& x) {
			auto bound = atg(x.m_bound);
			auto interval = atg(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> actg(const PwlInterval<T2>& x) {
			auto bound = actg(x.m_bound);
			auto interval = actg(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> log(const PwlInterval<T2>& x, T2 base) {
			auto bound = log(x.m_bound, base);
			auto interval = log(x.m_interval, base);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> ln(const PwlInterval<T2>& x) {
			auto bound = ln(x.m_bound);
			auto interval = ln(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> exp(const PwlInterval<T2>& x) {
			auto bound = exp(x.m_bound);
			auto interval = exp(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> abs(const PwlInterval<T2>& x) {
			auto bound = abs(x.m_bound);
			auto interval = abs(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> cos(const PwlInterval<T2>& x) {
			auto bound = cos(x.m_bound);
			auto interval = cos(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> sin(const PwlInterval<T2>& x) {
			auto bound = sin(x.m_bound);
			auto interval = sin(x.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> max(const PwlInterval<T2>& bound1, const PwlInterval<T2>& bound2) {
			auto bound = max(bound1.m_bound, bound2.m_bound);
			auto interval = snowgoose::interval::max({ bound1.m_interval, bound2.m_interval });
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> min(const PwlInterval<T2>& bound1, const PwlInterval<T2>& bound2) {
			auto bound = min(bound1.m_bound, bound2.m_bound);
			auto interval = snowgoose::interval::min({ bound1.m_interval, bound2.m_interval });
			return get_PwlInterval(bound, interval);
		}

		template<class T2> PwlInterval<T2> ifThen(IntervalBool ib, const PwlInterval<T2> &x, const PwlInterval<T2> &y) {
			auto bound = ifThen(ib, x.m_bound, y.m_bound);
			auto interval = ifThen(ib, x.m_interval, y.m_interval);
			return get_PwlInterval(bound, interval);
		}

		template<class T2> std::ostream& operator<<(std::ostream & out, const PwlInterval<T2> x) {
			out << x.m_bound << '\n';
			out << x.m_interval << '\n';
			return out;
		}
	}
}

#endif
