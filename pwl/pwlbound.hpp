#ifndef PWLBOUNDS_HPP
#define PWLBOUNDS_HPP

#include <iostream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <cmath>
#include "pwlfunc.hpp"
#include "utils.hpp"

void test_pwlbounds();

namespace snowgoose {
	namespace pwl {

		template<class T> class PwlBound
		{
		public:
			/**
			* Constructor
			* @param a is a left bound of x
			* @param b is a right bound of x
			* @param steps are number of steps along x to create a pwl bound for an elementary function 
			*/
			PwlBound(T a, T b, int steps = 5) : m_lower({ {a,a},{ b,b } }), m_upper({ { a,a }, { b,b } }), m_steps(steps){
			}
			/**
			* Adds two pwl bounds
			* @return a pwl bound
			*/
			PwlBound operator+(const PwlBound& bound) const;
			/**
			* Subtracts two pwl bounds
			* @return a pwl bound
			*/
			PwlBound operator-(const PwlBound& bound) const; 
			/**
			* Multiplies pwl bounds
			* @param bound is a right member
			* @return a pwl bound
			*/
			PwlBound operator*(const PwlBound& bound) const;
			/**
			* Divides pwl bounds
			* @param bound is a right member
			* @return pwl bound
			*/
			PwlBound operator/(const PwlBound& bound) const;
			/**
			* Raising to a power
			* @param exp is an integer 
			* @return a pwl bound
			*/
			PwlBound operator^(int exp) const;
			/**
			* Raising to a power
			* @param exp is a real exponent
			* @return a pwl bound
			*/
			PwlBound operator^(T exp) const;
			/**
			* Is a left upper bound less than a right lower one
			* @param y is a right member
			* @return IntervalBool
			*/
			IntervalBool operator < (const PwlBound &y) const { this->m_upper < y.m_lower; }
			/**
			* Is a left upper bound less or equal than a right lower one
			* @param y is a right member
			* @return IntervalBool
			*/
			IntervalBool operator <= (const PwlBound &y) const { this->m_upper <= y.m_lower; }
			/**
			* Is a left lower bound more than a right upper one
			* @param y is a right member
			* @return IntervalBool
			*/
			IntervalBool operator > (const PwlBound &y) const { this->m_lower > y.m_upper; }
			/**
			* Is a left lower bound more or equal than a right upper one
			* @param y is a right member
			* @return IntervalBool
			*/
			IntervalBool operator >= (const PwlBound &y) const { this->m_lower >= y.m_upper;  }
			/**
			* Is a left bound equals to a right one (lower and upper respectively)
			* @param y is right number
			* @return true or false
			*/
			bool operator == (const PwlBound &y) const { this->m_lower == y.m_lower && this->m_upper == y.m_upper; }
			/**
			* Returns a lower pwl bound
			* @return an ordered by x array of points
			*/
			std::vector<Point<T>> lower_bound() const { return m_lower.Points(); }
			/**
			* Returns an upper pwl bound
			* @return an ordered by x array of points
			*/
			std::vector<Point<T>> upper_bound() const { return m_upper.Points(); }
			/**
			* Returns a min value for a lower pwl bound
			* @return a min value
			*/
			T get_min() const { return m_lower.get_min(); }
			/**
			* Returns a max value for an upper pwl bound
			* @return a max value
			*/
			T get_max() const { return m_upper.get_max(); }
			/**
			* Returns left bound of x
			* @return value
			*/
			T get_a() const { return m_lower.get_a(); }
			/**
			* Returns right bound of x
			* @return value
			*/
			T get_b() const { return m_lower.get_b(); }

			template<class T2> friend PwlBound<T2> operator*(T2 x, const PwlBound<T2>& y);
			template<class T2> friend PwlBound<T2> operator*(const PwlBound<T2>& x, T2 y) { return y*x; }
			template<class T2> friend PwlBound<T2> operator/(T2 x, const PwlBound<T2>& y);
			template<class T2> friend PwlBound<T2> operator/(const PwlBound<T2>& x, T2 y);
			template<class T2> friend PwlBound<T2> operator+(T2 x, const PwlBound<T2>& y);
			template<class T2> friend PwlBound<T2> operator+(const PwlBound<T2>& x, T2 y);
			template<class T2> friend PwlBound<T2> operator-(T2 x, const PwlBound<T2>& y);
			template<class T2> friend PwlBound<T2> operator-(const PwlBound<T2>& x, T2 y);
			/**
			* Returns a pwl bound for the sqr function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> sqr(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the sqrt function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> sqrt(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the sin function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> sin(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the cos function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> cos(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the asin function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> asin(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the acos function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> acos(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the tg function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> tg(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the ctg function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> ctg(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the arctg function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> atg(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the arcctg function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> actg(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the log function
			* @param bound is a pwl bound
			* @param base is a base for log function
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> log(const PwlBound<T2>& bound, T2 base);
			/**
			* Returns a pwl bound for the ln function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> ln(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the exp function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> exp(const PwlBound<T2>& bound);
			/**
			* Returns a pwl bound for the abs function
			* @param bound is a pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> abs(const PwlBound<T2>& bound);
			/**
			* Returns a max upper and a max lower pwl bound from two ones
			* @param bound1 is a first pwl bound
			* @param bound2 is a second pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> max(const PwlBound<T2>& bound1, const PwlBound<T2>& bound2);
			/**
			* Returns a min upper and a min lower pwl bound from two ones
			* @param bound1 is a first pwl bound
			* @param bound2 is a second pwl bound
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> min(const PwlBound<T2>& bound1, const PwlBound<T2>& bound2);
			/**
			* Returns one of two bounds depending on the condition
			* @param x is a left pwl bound
			* @param y is a right pwl bound
			* @param ib is a condition
			* @return a pwl bound
			*/
			template<class T2> friend PwlBound<T2> ifThen(IntervalBool ib, const PwlBound<T2> &x, const PwlBound<T2> &y);
			/**
			* Output a pwl bound
			* @param out is an output stream
			* @param x is a pwl bound to print
			* @return output stream
			*/
			template<class T2> friend std::ostream& operator<<(std::ostream & out, const PwlBound<T2> x);
			friend void ::test_pwlbounds();
			template<class T2> friend class PwlBoundAlg;
		private:
			PwlBound(const PwlFunc<T>& lower, const PwlFunc<T>& upper, int steps = 5) : m_lower(lower), m_upper(upper), m_steps(steps) {
			}
			PwlBound get_pwl_bound(T a, T b) const;
			PwlFunc<T> get_pwl(const std::vector <Parabola<T>>& bounds, T len_step, bool is_lower = true) const;
			PwlFunc<T> get_pwl(const std::vector <FracLinFunc<T>>& bounds, T len_step, bool is_lower = true) const;
			PwlFunc<T> get_pwl(const Parabola<T>& bound, T len_step, bool is_lower = true) const;
			PwlFunc<T> get_pwl(const FracLinFunc<T>& bound, T len_step, bool is_lower = true) const;
			PwlBound get_copmosition(const PwlBound& out_bound, const std::set<Segment<T>>& out_seg) const;
			PwlBound<T> mul(const PwlBound<T>& bound1, const PwlBound<T>& bound2, T len_step) const;
			PwlBound<T> div(const PwlBound<T>& bound1, const PwlBound<T>& bound2, T len_step) const;
			template<class T2> friend void get_intersections(const std::set<Point<T2>>& points, T2 y, std::set<Point<T2>>& x);
			int m_steps;
			PwlFunc<T> m_lower;
			PwlFunc<T> m_upper;
		};

		template<class T> PwlBound<T> PwlBound<T>::operator+(const PwlBound<T>& bound) const {
			return PwlBound<T>(this->m_lower + bound.m_lower, this->m_upper + bound.m_upper, this->m_steps);
		}

		template<class T> PwlBound<T> PwlBound<T>::operator-(const PwlBound<T>& bound) const {
			return PwlBound<T>(this->m_lower - bound.m_upper, this->m_upper - bound.m_lower, this->m_steps);
		}

		template<class T> PwlBound<T> PwlBound<T>::mul(const PwlBound<T>& bound1, const PwlBound<T>& bound2, T len_step) const {
			PwlFunc<T> lower;
			PwlFunc<T> upper;
			if (moeq(bound1.m_lower.get_min(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = get_pwl(bound1.m_lower * bound2.m_lower, len_step);
				upper = get_pwl(bound1.m_upper * bound2.m_upper, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_upper.get_max(), 0.0)) {
				lower = get_pwl(bound1.m_upper * bound2.m_lower, len_step);
				upper = get_pwl(bound1.m_lower * bound2.m_upper, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_lower * bound2.m_lower, len_step), get_pwl(bound1.m_upper * bound2.m_lower, len_step) });
				upper = max<T>({ get_pwl(bound1.m_upper * bound2.m_upper, len_step, false), get_pwl(bound1.m_lower * bound2.m_upper, len_step, false) });
			}
			else if (lseq(bound1.m_upper.get_max(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = get_pwl(bound1.m_lower * bound2.m_upper, len_step);
				upper = get_pwl(bound1.m_upper * bound2.m_lower, len_step, false);
			}
			else if (lseq(bound1.m_upper.get_max(), 0.0) && lseq(bound2.m_upper.get_max(), 0.0)) {
				lower = get_pwl(bound1.m_upper * bound2.m_upper, len_step);
				upper = get_pwl(bound1.m_lower * bound2.m_lower, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_upper * bound2.m_upper, len_step), get_pwl(bound1.m_lower * bound2.m_upper, len_step) });
				upper = max<T>({ get_pwl(bound1.m_lower * bound2.m_lower, len_step, false), get_pwl(bound1.m_upper * bound2.m_lower, len_step, false) });
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_lower * bound2.m_lower, len_step), get_pwl(bound1.m_lower * bound2.m_upper, len_step) });
				upper = max<T>({ get_pwl(bound1.m_upper * bound2.m_upper, len_step, false), get_pwl(bound1.m_upper * bound2.m_lower, len_step, false) });
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && lseq(bound2.m_upper.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_upper * bound2.m_upper, len_step), get_pwl(bound1.m_upper * bound2.m_lower, len_step) });
				upper = max<T>({ get_pwl(bound1.m_lower * bound2.m_lower, len_step, false), get_pwl(bound1.m_lower * bound2.m_upper, len_step, false) });
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && lseq(bound2.m_upper.get_min(), 0.0)
				&& lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_upper * bound2.m_lower, len_step), get_pwl(bound1.m_lower * bound2.m_upper, len_step) });
				upper = max<T>({ get_pwl(bound1.m_lower * bound2.m_lower, len_step, false), get_pwl(bound1.m_upper * bound2.m_upper, len_step, false) });
			}
			else {
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::mul(const PwlBound<T>& bound1, const PwlBound<T>& bound2, T len_step). Invalid bounds.");
			}
			return PwlBound<T>(lower, upper, this->m_steps);
		}

		template<class T> PwlBound<T> PwlBound<T>::operator*(const PwlBound& bound) const {
			std::set<Point<T>> x({ {this->get_a(),0.0}, {this->get_b(), 0.0} });
			T y = 0.0;
			get_intersections(this->m_lower.m_points, y, x);
			get_intersections(this->m_upper.m_points, y, x);
			get_intersections(bound.m_lower.m_points, y, x);
			get_intersections(bound.m_upper.m_points, y, x);

			std::set<Segment<T>> in_seg;
			auto last_x = std::prev(x.cend());
			for (auto it = x.cbegin(); it != last_x; ++it) {
				in_seg.insert({ it->x, std::next(it)->x });
			}
			PwlFunc<T> lower;
			PwlFunc<T> upper;
			T len_step = (this->get_b() - this->get_a()) / this->m_steps;
			for (const Segment<T>& seg : in_seg) {
				PwlBound<T> bound1 = this->get_pwl_bound(seg.start, seg.end);
				PwlBound<T> bound2 = bound.get_pwl_bound(seg.start, seg.end);
				PwlBound<T> bound = mul(bound1, bound2, len_step);
				lower.insert(bound.m_lower.m_points);
				upper.insert(bound.m_upper.m_points);
			}
			return PwlBound<T>(lower, upper, this->m_steps);
		}

		template<class T> PwlBound<T> PwlBound<T>::div(const PwlBound<T>& bound1, const PwlBound<T>& bound2, T len_step) const {
			PwlFunc<T> lower;
			PwlFunc<T> upper;
	
			if (moeq(bound1.m_lower.get_min(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = get_pwl(bound1.m_lower / bound2.m_upper, len_step);
				upper = get_pwl(bound1.m_upper / bound2.m_lower, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_upper.get_max(), 0.0)) {
				lower = get_pwl(bound1.m_upper / bound2.m_upper, len_step);
				upper = get_pwl(bound1.m_lower / bound2.m_lower, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator/(const PwlBound& bound2). The bounds are not defined.");
			}
			else if (lseq(bound1.m_upper.get_max(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = get_pwl(bound1.m_lower / bound2.m_lower, len_step);
				upper = get_pwl(bound1.m_upper / bound2.m_upper, len_step, false);
			}
			else if (lseq(bound1.m_upper.get_max(), 0.0) && lseq(bound2.m_upper.get_max(), 0.0)) {
				lower = get_pwl(bound1.m_upper / bound2.m_lower, len_step);
				upper = get_pwl(bound1.m_lower / bound2.m_upper, len_step, false);
			}
			else if (moeq(bound1.m_lower.get_min(), 0.0) && lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator/(const PwlBound& bound2). The bounds are not defined.");
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && moeq(bound2.m_lower.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_lower / bound2.m_lower, len_step), get_pwl(bound1.m_lower / bound2.m_upper, len_step) });
				upper = max<T>({ get_pwl(bound1.m_upper / bound2.m_upper, len_step, false), get_pwl(bound1.m_upper / bound2.m_lower, len_step, false) });
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && lseq(bound2.m_upper.get_min(), 0.0)) {
				lower = min<T>({ get_pwl(bound1.m_upper / bound2.m_upper, len_step), get_pwl(bound1.m_upper / bound2.m_lower, len_step) });
				upper = max<T>({ get_pwl(bound1.m_lower / bound2.m_lower, len_step, false), get_pwl(bound1.m_lower / bound2.m_upper, len_step, false) });
			}
			else if (lseq(bound1.m_lower.get_max(), 0.0) && moeq(bound1.m_upper.get_min(), 0.0) && lseq(bound2.m_upper.get_min(), 0.0)
				&& lseq(bound2.m_lower.get_max(), 0.0) && moeq(bound2.m_upper.get_min(), 0.0)) {
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator/(const PwlBound& bound2). The bounds are not defined.");
			}
			else {
				throw std::invalid_argument("Exception in PwlBound<T>::operator/(const PwlBound& bound2). Invalid bounds.");
			}
			return PwlBound<T>(lower, upper, this->m_steps);
		}

		template<class T> PwlBound<T> PwlBound<T>::operator/(const PwlBound& bound) const {
			std::set<Point<T>> x({ {this->get_a(), 0.0}, {this->get_b(), 0.0} });
			T y = 0.0;
			get_intersections(this->m_lower.m_points, y, x);
			get_intersections(this->m_upper.m_points, y, x);
			get_intersections(bound.m_lower.m_points, y, x);
			get_intersections(bound.m_upper.m_points, y, x);

			std::set<Segment<T>> in_seg;
			auto last_x = std::prev(x.cend());
			for (auto it = x.cbegin(); it != last_x; ++it) {
				in_seg.insert({ it->x, std::next(it)->x });
			}
			PwlFunc<T> lower;
			PwlFunc<T> upper;
			T len_step = (this->get_b() - this->get_a()) / this->m_steps;
			for (const Segment<T>& seg : in_seg) {
				PwlBound<T> bound1 = this->get_pwl_bound(seg.start, seg.end);
				PwlBound<T> bound2 = bound.get_pwl_bound(seg.start, seg.end);
				PwlBound<T> bound_seg = div(bound1, bound2, len_step);
				lower.insert(bound_seg.m_lower.m_points);
				upper.insert(bound_seg.m_upper.m_points);
			}
			return PwlBound<T>(lower, upper, this->m_steps);
		}

		template<class T2, class T3> void pow_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u, T3 exp) {
			pow_convex(c, d, len_step, u, l, exp);
		}

		template<class T2, class T3> void pow_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u, T3 exp) {
			T2 y1 = pow(c, exp);
			T2 y2 = pow(d, exp);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));

			//lower bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = pow(x, exp);
				u.insert({ { x,y } });
			}
			//upper bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = pow(x1, exp - 1);
				T2 z2 = pow(x2, exp - 1);
				T2 z3 = z1 * x1;
				T2 z4 = z2 * x2;
				T2 x = (1 - exp)*(z4 - z3) / (exp*(z1 - z2)); //between current and next abscissa
				T2 y = z3 + exp*z1*(x - x1);
				l.insert({ { x,y } });
			}
		}

		template<class T> PwlBound<T> PwlBound<T>::operator^(int exp) const {
			if (exp >= -1 && exp <= 1)
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator^(int exp). Invalid argument.");

			T c = this->get_min();
			T d = this->get_max();

			PwlFunc<T> l;
			PwlFunc<T> u;
			T len_step = (d - c) / this->m_steps;
			if (exp >= 2) {
				if (exp % 2) {//odd		
					if (ls(c, 0.0) && mo(d, 0.0)) {
						pow_concave(c, 0.0, len_step, l, u, exp);
						pow_convex(0.0, d, len_step, l, u, exp);
					}
					else if (lseq(d, 0.0)) {
						pow_concave(c, d, len_step, l, u, exp);
					}
					else {
						pow_convex(c, d, len_step, l, u, exp);
					}
					return PwlBound<T>(l(this->m_lower), u(this->m_upper), this->m_steps);
				}
				else {//even
					if (ls(c, 0.0) && mo(d, 0.0)) {
						pow_convex(c, 0.0, len_step, l, u, exp);
						pow_convex(0.0, d, len_step, l, u, exp);
						PwlBound<T> out_bound(l, u);
						std::set<Segment<T>> y_seg({ { c, 0, false },{ 0, d, true } });
						return this->get_copmosition(out_bound, y_seg);
					}
					else {
						pow_convex(c, d, len_step, l, u, exp);
						if (lseq(d, 0.0)) {
							return PwlBound<T>(l(this->m_upper), u(this->m_lower), this->m_steps);//decreasing
						}
						else {
							return PwlBound<T>(l(this->m_lower), u(this->m_upper), this->m_steps);//increasing
						}
					}
				}
			}
			else { // <=-2
				if(ls(c, 0.0) && mo(d, 0.0))
					throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator^(int exp). Invalid domain.");
				if (exp % 2) { //odd
					if (lseq(d, 0.0)) {
						pow_concave(c, d, len_step, l, u, exp);
					}
					else {
						pow_convex(c, d, len_step, l, u, exp);
					}
					return PwlBound<T>(l(this->m_upper), u(this->m_lower), this->m_steps);//decreasing
				}
				else {
					pow_convex(c, d, len_step, l, u, exp);
					if (lseq(d, 0.0)) {
						return PwlBound<T>(l(this->m_lower), u(this->m_upper), this->m_steps);//increasing
					}
					else {
						return PwlBound<T>(l(this->m_upper), u(this->m_lower), this->m_steps);//decreasing
					}
				}
			}
		}

		template<class T> PwlBound<T> PwlBound<T>::operator^(T exp) const {
			T c = this->get_min();
			T d = this->get_max();

			if((ls(c, 0.0) && exp < 0) || (lseq(c, 0.0) && exp > 0))
				throw std::invalid_argument("Exception in PwlBound<T> PwlBound<T>::operator^(T exp). Invalid domain.");

			PwlFunc<T> l;
			PwlFunc<T> u;
			T len_step = (d - c) / this->m_steps;

			if (ls(exp, 0.0)) {
				pow_convex(c, d, len_step, l, u, exp);
				return PwlBound<T>(l(this->m_upper), u(this->m_lower), this->m_steps);//decreasing
			}
			else if (moeq(exp, 0.0) && ls(exp, 1.0)) {
				pow_concave(c, d, len_step, l, u, exp);
				return PwlBound<T>(l(this->m_lower), u(this->m_upper), this->m_steps);//increasing
			}
			else {
				pow_convex(c, d, len_step, l, u, exp);
				return PwlBound<T>(l(this->m_lower), u(this->m_upper), this->m_steps);//increasing
			}
		}

		template<class T> PwlBound<T> PwlBound<T>::get_pwl_bound(T a, T b) const {
			return PwlBound<T>(this->m_lower.get_pwl_func(a, b), this->m_upper.get_pwl_func(a, b), this->m_steps);
		}

		template<class T2> PwlFunc<T2> parabola_concave(const Parabola<T2>& bound, T2 len_step, bool is_lower = true) {
			return parabola_convex(bound, len_step, !is_lower);
		}

		template<class T2> PwlFunc<T2> parabola_convex(const Parabola<T2>& bound, T2 len_step, bool is_lower = true) {
			T2 start = bound.start->x;
			T2 end = bound.end->x;
			T2 a = bound.coeff.a;
			T2 b = bound.coeff.b;
			T2 c = bound.coeff.c;
			T2 y1 = bound.start->y;
			T2 y2 = bound.end->y;
			PwlFunc<T2> pwl({ { start, y1 },{ end, y2 } });

			int steps = static_cast<int>(ceill((end - start) / len_step));

			if (is_lower) {
				for (int i = 0; i < steps; i++) {
					T2 x1 = start + i*len_step;//current
					T2 x2 = moeq(x1 + len_step, end) ? end : x1 + len_step;//next
					T2 x = (x1 + x2) / 2.0;
					T2 y = get_val<T2>(a, b, c, x1) + (2 * a*x1 + b)*(x - x1);
					pwl.insert({ { x,y } });
				}
			}
			else {
				for (int i = 1; i < steps; i++) {
					T2 x = start + i*len_step;
					T2 y = get_val<T2>(a, b, c, x);
					pwl.insert({ { x,y } });
				}
			}
			return pwl;

		}

		template<class T> PwlFunc<T> PwlBound<T>::get_pwl(const std::vector <Parabola<T>>& bounds, T len_step, bool is_lower) const {
			std::set<Point<T>> points;
			for (auto& bound : bounds) {
				PwlFunc<T> f = get_pwl(bound, len_step, is_lower);
				points.insert(f.m_points.begin(), f.m_points.end());
			}
			return PwlFunc<T>(points);
		}

		template<class T> PwlFunc<T> PwlBound<T>::get_pwl(const Parabola<T>& p, T len_step, bool is_lower) const {
			std::set<Point<T>> points;
			T a = p.coeff.a;
			T b = p.coeff.b;
			T c = p.coeff.c;
			Point<T> start = *p.start;
			Point<T> end = *p.end;

			if (eq(a, 0.0)) {
				return PwlFunc<T>({ start, end });
			}
			else {
				T top = -b / (2 * a);
				if (mo(top, start.x) && ls(top, end.x)) {
					T y = get_val<T>(a, b, c, top);
					Parabola<T> p1, p2;
					p1.coeff = { a,b,c };
					p1.start = p.start;
					p1.end = std::shared_ptr<Point<T>>(new Point<T>({ top, y }));
					p2.coeff = { a,b,c };
					p2.start = std::shared_ptr<Point<T>>(new Point<T>({ top, y }));
					p2.end = p.end;
					if (mo(a, 0.0)) {
						PwlFunc<T> pwl = parabola_convex(p1, len_step, is_lower);
						pwl.insert(parabola_convex(p2, len_step, is_lower).m_points);
						return pwl;
					}
					else {
						PwlFunc<T> pwl = parabola_concave(p1, len_step, is_lower);
						pwl.insert(parabola_concave(p2, len_step, is_lower).m_points);
						return pwl;
					}
				}
				else {
					if (mo(a, 0.0)) {
						return parabola_convex(p, len_step, is_lower);
					}else{ 
						return parabola_concave(p, len_step, is_lower);
					}
				}
			}
		}

		template<class T2> PwlFunc<T2> frac_lin_concave(const FracLinFunc<T2>& bound, T2 len_step, bool is_lower = true) {
			return frac_lin_convex(bound, len_step, !is_lower);
		}

		template<class T2> PwlFunc<T2> frac_lin_convex(const FracLinFunc<T2>& bound, T2 len_step, bool is_lower = true) {
			T2 start = bound.start->x;
			T2 end = bound.end->x;
			T2 a = bound.coeff.a;
			T2 b = bound.coeff.b;
			T2 c = bound.coeff.c;
			T2 d = bound.coeff.d;
			T2 y1 = bound.start->y;
			T2 y2 = bound.end->y;
			PwlFunc<T2> pwl({ { start, y1 },{ end, y2 } });

			int steps = static_cast<int>(ceill((end - start) / len_step));

			if (is_lower) {
				for (int i = 0; i < steps; i++) {
					T2 x1 = start + i*len_step;//current
					T2 x2 = moeq(x1 + len_step,  end) ? end : x1 + len_step;//next
					T2 e1 = a*d - c*b;
					T2 e2 = c*x1 + d;
					T2 e22 = e2*e2;
					T2 e3 = c*x2 + d;
					T2 e32 = e3*e3;
					T2 e4 = a*x1 + b;
					T2 e5 = a*x2 + b;
					T2 e6 = e1 / e22;
					T2 e7 = e1 / e32;
					T2 e8 = e5 / e3;
					T2 e9 = e4 / e2;
					T2 e10 = x1*e1 / e22;
					T2 e11 = x2*e1 / e32;
					T2 x = (e8 - e9 + e10 - e11) / (e6 - e7); //between current and next abscissa
					T2 y = e9 + e6*(x - x1);
					pwl.insert({ { x,y } });
				}
			}
			else {
				for (int i = 1; i < steps; i++) {
					T2 x = start + i*len_step;
					T2 y = (a*x + b) / (c*x + d);
					pwl.insert({ { x,y } });
				}
			}
			return pwl;
		}

		template<class T> PwlFunc<T> PwlBound<T>::get_pwl(const std::vector <FracLinFunc<T>>& bounds, T len_step, bool is_lower) const {
			std::set<Point<T>> points;
			for (auto& bound : bounds) {
				PwlFunc<T> f = get_pwl(bound, len_step, is_lower);
				points.insert(f.m_points.begin(), f.m_points.end());
			}
			return PwlFunc<T>(points);
		}

		//https://ege-ok.ru/2013/11/12/asimptotyi-grafika-funktsii
		template<class T> PwlFunc<T> PwlBound<T>::get_pwl(const FracLinFunc<T>& bound, T len_step, bool is_lower) const {
			bool is_c_zero = eq(bound.coeff.c, 0.0);
			bool is_d_zero = eq(bound.coeff.d, 0.0);

			if(is_c_zero && is_d_zero)
				throw std::invalid_argument("Exception in PwlFunc<T> PwlBound<T>::get_pwl(const FracLinFunc<T>& bound, bool is_lower = true). Division by 0.0.");
			if (is_c_zero || eq(bound.coeff.a * bound.coeff.d, bound.coeff.b * bound.coeff.c)) { //straight lines match or denominator is constant
				return PwlFunc<T>({ *bound.start,*bound.end });
			}
			T ver_asym_x = -bound.coeff.d / bound.coeff.c;
			T hor_asym_y = bound.coeff.a / bound.coeff.c;

			if (lseq(bound.start->x, ver_asym_x) && moeq(bound.end->x, ver_asym_x)) {
				throw std::invalid_argument("Exception in PwlFunc<T> PwlBound<T>::get_pwl(const FracLinFunc<T>& bound, bool is_lower = true). Invalid domain");
			}
			int quarter = get_quarter<T>({ bound.start->x, bound.start->y }, { ver_asym_x, hor_asym_y });
			if (quarter == 1 || quarter == 3) {
				if (ls(bound.end->x, ver_asym_x)) {
					return frac_lin_concave(bound, len_step, is_lower);
				}
				else {
					return frac_lin_convex(bound, len_step, is_lower);
				}
			}
			else {
				if (ls(bound.end->x, ver_asym_x)) {
					return frac_lin_convex(bound, len_step, is_lower);
				}
				else {
					return frac_lin_concave(bound, len_step, is_lower);
				}
			}
		}


		template<class T2> PwlBound<T2> operator*(T2 x, const PwlBound<T2>& y) {
			T2 a = y.get_a();
			T2 b = y.get_b();
			PwlFunc<T2> lower({ { a,x },{ b,x } }), upper({ { a,x },{ b,x } });
			PwlBound<T2> z(lower, upper, y.m_steps);
			return z * y;
		}

		template<class T2> PwlBound<T2> operator/(T2 x, const PwlBound<T2>& y) {
			T2 a = y.get_a();
			T2 b = y.get_b();
			PwlFunc<T2> lower({ { a,x },{ b,x } }), upper({ { a,x },{ b,x } });
			PwlBound<T2> z(lower, upper, y.m_steps);
			return z / y;
		}

		template<class T2> PwlBound<T2> operator/(const PwlBound<T2>& x, T2 y) {
			T2 a = x.get_a();
			T2 b = x.get_b();
			PwlFunc<T2> lower({ { a,y },{ b,y } }), upper({ { a,y },{ b,y } });
			PwlBound<T2> z(lower, upper, x.m_steps);
			return x/z;
		}

		template<class T2> PwlBound<T2> operator+(T2 x, const PwlBound<T2>& y) {
			return PwlBound<T2>(x + y.m_lower, x + y.m_upper, y.m_steps);
		}

		template<class T2> PwlBound<T2> operator+(const PwlBound<T2>& x, T2 y) {
			return PwlBound<T2>(x.m_lower + y, x.m_upper + y, x.m_steps);
		}

		template<class T2> PwlBound<T2> operator-(T2 x, const PwlBound<T2>& y) {
			return PwlBound<T2>(x - y.m_upper, x - y.m_lower, y.m_steps);
		}

		template<class T2> PwlBound<T2> operator-(const PwlBound<T2>& x, T2 y) {
			return PwlBound<T2>(x.m_lower - y, x.m_upper - y, x.m_steps);
		}

		template<class T2> PwlBound<T2> sqr(const PwlBound<T2>& bound) {
			return bound ^ 2;
		}

		template<class T2> PwlBound<T2> sqrt(const PwlBound<T2>& bound)
		{
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			if (ls(c, 0.0))
				throw std::invalid_argument("Exception in PwlBound<T2> sqrt(const PwlBound<T2>& bound). The function sqrt is not define for negative numbers");

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 y1 = ::sqrt(c);
			T2 y2 = ::sqrt(d);
	
			l.insert({ {c, y1}, {d, y2} });
			eq(c, 0.0) ? u.insert({ { d, y2 } }) : u.insert({ {c, y1}, {d, y2} });

			T2 len_step = (d - c) / bound.m_steps;
			int steps = bound.m_steps;
			//lower bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::sqrt(x);
				l.insert({ { x,y } });
			}
			//upper bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 x = ::sqrt(x1*x2); //between current and next abscissa
				T2 y = (::sqrt(x1) + ::sqrt(x2)) / 2;
				u.insert({ { x,y } });
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> void asin_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = ::asin(c);
			T2 y2 = ::asin(d);
			T2 x_next = mo(c + len_step, 0.0) ? 0.0 : c + len_step;
			T2 y3 = eq(c, -1.0) ? ::asin(x_next) - len_step / ::sqrt(1 - x_next*x_next) : y1;
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y3 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));

			//lower bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::asin(x);
				l.insert({ { x,y } });
			}
			//upper bound
			int i = eq(c, -1.0) ? 1 : 0;
			for (; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = ::sqrt(1.0 - x1*x1);
				T2 z2 = ::sqrt(1.0 - x2*x2);
				T2 x = (z2*(x1 + z1*::asin(x2)) - z1*(x2 + z2*::asin(x1))) / (z2 - z1); //between current and next abscissa
				T2 y = ::asin(x1) + (x - x1) / z1;
				u.insert({ { x,y } });
			}
		}

		template<class T2> void asin_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = ::asin(c);
			T2 y2 = ::asin(d);
			T2 x_before_last = ls(d - len_step, 0.0) ? 0.0 : d - len_step;
			T2 y3 = eq(d, 1.0) ? ::asin(x_before_last) + len_step / ::sqrt(1 - x_before_last*x_before_last) : y2;
			l.insert({ { c, y1 },{ d, y3 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));
			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::asin(x);
				u.insert({ { x,y } });
			}
			//lower bound
			int last = eq(d, 1.0) ? steps - 1 : steps;
			for (int i = 0; i < last; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = ::sqrt(1 - x1*x1);
				T2 z2 = ::sqrt(1 - x2*x2);
				T2 x = (z2*(x1 + z1*::asin(x2)) - z1*(x2 + z2*::asin(x1))) / (z2 - z1); //between current and next abscissa
				T2 y = ::asin(x1) + (x - x1) / z1;
				l.insert({ { x,y } });
			}
		}

		template<class T2> PwlBound<T2> asin(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			if (ls(c, -1.0) || mo(d, 1.0))
				throw std::invalid_argument("Exception in PwlBound<T2> asin(const PwlBound<T2>& bound).The argument is out of this interval[-1, 1].");

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 len_step = (d - c) / bound.m_steps;
			if ( ls(c, 0.0) && mo(d, 0.0) ) {
				asin_concave(c, 0.0, len_step, l, u);
				asin_convex(0.0, d, len_step, l, u);
			}
			else if (lseq(d, 0.0)) {
				asin_concave(c, d, len_step, l, u);
			}
			else{
				asin_convex(c, d, len_step, l, u);
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> PwlBound<T2> acos(const PwlBound<T2>& bound) {
			return M_PI / 2 - asin(bound);
		}

		template<class T2> void tg_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			tg_convex(c, d, len_step, u, l);
		}

		template<class T2> void tg_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = ::tan(c);
			T2 y2 = ::tan(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));
			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::tan(x);
				u.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = ::cos(x1) * ::cos(x1);
				T2 z2 = ::cos(x2) * ::cos(x2);
				T2 x = (z2*(x1 + z1*::tan(x2)) - z1*(x2 + z2*::tan(x1))) / (z2 - z1); //between current and next abscissa
				T2 y = ::tan(x1) + (x - x1) / z1;
				l.insert({ { x,y } });
			}
		}

		template<class T2> PwlBound<T2> tg(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			bool p2 = false;
			int k = (int)ceil(c * M_2_PI);
			int m = (int)floor(d * M_2_PI);
			for (int i = k; i <= m; i++) {
				if ((i % 2) == 1 || (i % 2) == -1)
					p2 = true;
				if (p2)
					break;
			}
			if (p2)
				throw std::invalid_argument("Invalid argument in PwlBound<T2> tg(const PwlBound<T2>& bound) function. Tg is not defined for PI*N/2.");

			PwlFunc<T2> l;
			PwlFunc<T2> u;
	
			T2 len_step = (d - c) / bound.m_steps;
			T2 t = static_cast<int>(floor(c + M_PI_2)  / M_PI) * M_PI;
	
			if (ls(c, t) && mo(d, t)) {
				tg_concave(c, t, len_step, l, u);
				tg_convex(t, d, len_step, l, u);
			}
			else if (lseq(d, t)) {
				tg_concave(c, d, len_step, l, u);
			}
			else {
				tg_convex(c, d, len_step, l, u);
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> void ctg_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			ctg_convex(c, d, len_step, u, l);
		}

		template<class T2> void ctg_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = 1 / ::tan(c);
			T2 y2 = 1 / ::tan(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));
			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = 1 / ::tan(x);
				u.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = ::sin(x1)*::sin(x1);
				T2 z2 = ::sin(x2)*::sin(x2);
				T2 x = (z1*(x2 + z2 / ::tan(x2)) - z2*(x1 + z1 / ::tan(x1))) / (z1 - z2); //between current and next abscissa
				T2 y = 1 / ::tan(x1) + (x - x1) / -z1;
				l.insert({ { x,y } });
			}
		}

		template<class T2> PwlBound<T2> ctg(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			bool p2 = false;
			int k = (int)ceil(c * M_2_PI);
			int m = (int)floor(d * M_2_PI);
			for (int i = k; i <= m; i++) {
				if ((i % 2) == 0)
					p2 = true;
				if (p2)
					break;
			}
			if (p2)
				throw std::invalid_argument("Invalid argument in PwlBound<T2> ctg(const PwlBound<T2>& bound) function. Ctg is not defined for PI*N.");

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 len_step = (d - c) / bound.m_steps;

			T2 t = static_cast<int>(floor(c / M_PI)) * M_PI + M_PI_2;

			if (ls(c, t) && mo(d, t)) {
				ctg_convex(c, t, len_step, l, u);
				ctg_concave(t, d, len_step, l, u);
			}
			else if (lseq(d, t)) {
				ctg_convex(c, d, len_step, l, u);
			}
			else {
				ctg_concave(c, d, len_step, l, u);
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> void atg_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			atg_convex(c, d, len_step, u, l);
		}

		template<class T2> void atg_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = ::atan(c);
			T2 y2 = ::atan(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));
			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::atan(x);
				u.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 z1 = 1 + x1*x1;
				T2 z2 = 1 + x2*x2;
				T2 x = (z2*(x1 + z1*::atan(x2)) - z1*(x2 + z2*::atan(x1))) / (z2 - z1); //between current and next abscissa
				T2 y = ::atan(x1) + (x - x1) / z1;
				l.insert({ { x,y } });
			}
		}

		template<class T2> PwlBound<T2> atg(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 len_step = (d - c) / bound.m_steps;
			if (ls(c, 0.0) && mo(d, 0.0)) {
				atg_convex(c, 0.0, len_step, l, u);
				atg_concave(0.0, d, len_step, l, u);
			}
			else if (moeq(d, 0.0)) {
				atg_convex(c, d, len_step, l, u);
			}
			else {
				atg_concave(c, d, len_step, l, u);
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> PwlBound<T2> actg(const PwlBound<T2>& bound) {
			return M_PI / 2 - atg(bound);
		}

		template<class T2> PwlBound<T2> log(const PwlBound<T2>& bound, T2 base) {
			return 1.0 / ::log(base) * ln(bound);
		}

		template<class T2> PwlBound<T2> ln(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			if (lseq(c, 0.0))
				throw std::invalid_argument("Exception in PwlBound<T2> ln(const PwlBound<T2>& bound). Invalid domain.");

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 y1 = ::log(c);
			T2 y2 = ::log(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });

			T2 len_step = (d - c) / bound.m_steps;
			int steps = bound.m_steps;

			//lower bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::log(x);
				l.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 ln1 = ::log(x1);
				T2 ln2 = ::log(x2);
				T2 x = x1*x2*(ln2 - ln1)/(x2-x1); //between x1 and x2 abscissa
				T2 y = ln1 + x2*(ln2 - ln1)/(x2-x1) -1;
				u.insert({ { x,y } });
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> PwlBound<T2> exp(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 y1 = ::exp(c);
			T2 y2 = ::exp(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });

			T2 len_step = (d - c) / bound.m_steps;
			int steps = bound.m_steps;

			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::exp(x);
				u.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 exp1 = ::exp(x1);
				T2 exp2 = ::exp(x2);
				T2 x = (exp1*x1 - exp2*x2)/(exp1 - exp2)-1; //between x1 and x2 abscissa
				T2 y = exp1 + exp1*(x - x1);
				l.insert({ { x,y } });
			}
			return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
		}

		template<class T2> void get_intersections(const std::set<Point<T2>>& points, T2 y, std::set<Point<T2>>& x) {
			auto last_point = std::prev(points.cend());
			for (auto it = points.cbegin(); it != last_point; ++it) {
				auto next = std::next(it);
				if (moeq(y, std::min(it->y, next->y)) && lseq(y, std::max(it->y, next->y))) {
					T2 ab = std::abs(it->y - y);
					T2 dc = std::abs(next->y - y);
					T2 bc = next->x - it->x;
					T2 new_x = it->x + ab*bc / (dc + ab);
					x.insert({ { new_x, y } });
				}
			}
		}

		template<class T> PwlBound<T> PwlBound<T>::get_copmosition(const PwlBound<T>& out_bound, const std::set<Segment<T>>& out_seg) const {
			if (this->m_lower == this->m_upper)
				return PwlBound<T>(out_bound.m_lower(this->m_lower), out_bound.m_upper(this->m_upper), this->m_steps);
			std::set<Point<T>> x({ {this->get_a(), 0.0}, { this->get_b(), 0.0} });
			auto last_seg = std::prev(out_seg.cend());
			for (auto it = out_seg.cbegin(); it != last_seg; ++it) {
				T y = it->end;
				get_intersections(this->m_lower.m_points, y, x);
				get_intersections(this->m_upper.m_points, y, x);
			}
			std::set<Segment<T>> in_seg;
			auto last_x = std::prev(x.cend());
			for (auto it = x.cbegin(); it != last_x; ++it) {
				in_seg.insert({ it->x, std::next(it)->x });
			}
			PwlFunc<T> l;
			PwlFunc<T> u;

			for (const Segment<T>& seg : in_seg) {
				PwlBound<T> in = this->get_pwl_bound(seg.start, seg.end);
				T c = in.m_lower.get_min();
				T d = in.m_upper.get_max();
				PwlBound<T> out = out_bound.get_pwl_bound(c, d);
				auto layer1 = out_seg.upper_bound({ c, c });
				auto layer2 = out_seg.lower_bound({ d, d });
				if (layer1 == layer2) {
					std::set<Point<T>> lower, upper;
					if (layer1->isIncrease) {
						lower = out.m_lower(in.m_lower).m_points;
						upper = out.m_upper(in.m_upper).m_points;
					}
					else {
						lower = out.m_lower(in.m_upper).m_points;
						upper = out.m_upper(in.m_lower).m_points;
					}
					l.m_points.insert(lower.begin(), lower.end());
					u.m_points.insert(upper.begin(), upper.end());
				}
				else if (std::next(layer1) == layer2) {//layers are adjacent
					if (layer1->isIncrease && !layer2->isIncrease) {//outer func is concave
						std::set<Point<T>> lower = min<T>(out.m_lower(in.m_lower), out.m_lower(in.m_upper)).m_points;
						l.m_points.insert(lower.begin(), lower.end());
						T max_out = out.m_upper.get_max();
						u.m_points.insert({ {seg.start, max_out },{ seg.end, max_out} });
					}
					else {//outer func is convex
						T min_out = out.m_lower.get_min();
						l.m_points.insert({ { seg.start, min_out },{ seg.end, min_out } });
						std::set<Point<T>> upper = max<T>(out.m_upper(in.m_lower), out.m_upper(in.m_upper)).m_points;
						u.m_points.insert(upper.begin(), upper.end());
					}
				}
				else {//layers aren't adjacent
					T min_out = out.m_lower.get_min();
					T max_out = out.m_upper.get_max();
					l.insert({ { seg.start, min_out },{ seg.end, min_out } });
					u.insert({ { seg.start, max_out },{ seg.end, max_out } });
				}
			}
			return PwlBound<T>(l, u, this->m_steps);
		}

		template<class T2> PwlBound<T2> abs(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();

			PwlFunc<T2> l;
			PwlFunc<T2> u;

			T2 y1 = std::abs(c);
			T2 y2 = std::abs(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });

			if (moeq(c, 0.0)) {
				return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
			}
			else if(lseq(d, 0.0))
			{
				return PwlBound<T2>(l(bound.m_upper), u(bound.m_lower), bound.m_steps);
			}
			else {// c < 0 && d > 0
				l.insert({ { 0, 0 } });
				u.insert({ { 0, 0 } });
				PwlBound<T2> out_bound(l, u, bound.m_steps);
				std::set<Segment<T2>> y_seg({ { c, 0, false }, { 0, d, true} });
				return bound.get_copmosition(out_bound, y_seg);
			}
		}

		template<class T2> void cos_convex(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 y1 = ::cos(c);
			T2 y2 = ::cos(d);
			l.insert({ { c, y1 },{ d, y2 } });
			u.insert({ { c, y1 },{ d, y2 } });
			int steps = static_cast<int>(ceill((d - c) / len_step));
			//upper bound
			for (int i = 1; i < steps; i++) {
				T2 x = c + i*len_step;
				T2 y = ::cos(x);
				u.insert({ { x,y } });
			}
			//lower bound
			for (int i = 0; i < steps; i++) {
				T2 x1 = c + i*len_step;
				T2 x2 = moeq(x1 + len_step, d) ? d : x1 + len_step;
				T2 cosx1 = ::cos(x1);
				T2 cosx2 = ::cos(x2);
				T2 sinx1 = ::sin(x1);
				T2 sinx2 = ::sin(x2);
				T2 x = (cosx2 - cosx1 + sinx2 * x2 - sinx1*x1)/(sinx2-sinx1); //between current and next abscissa
				T2 y = cosx1 - sinx1*(x - x1);
				l.insert({ { x,y } });
			}
		}
		template<class T2> void cos_concave(T2 c, T2 d, T2 len_step, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			cos_convex(c, d, len_step, u, l);
		}

		template<class T2> void cos_build_bound(T2 c, T2 d, T2 len_step, bool is_increase, PwlFunc<T2>& l, PwlFunc<T2>& u) {
			T2 t = static_cast<int>(floor(c / M_PI)) * M_PI + M_PI_2;
			if (ls(c, t) && mo(d,t)) {
				if (is_increase) {
					cos_convex(c, t, len_step, l, u);
					cos_concave(t, d, len_step, l, u);
				}
				else {
					cos_concave(c, t, len_step, l, u);
					cos_convex(t, d, len_step, l, u);
				}			
			}
			else if (lseq(d,t)) {
				if (is_increase) {
					cos_convex(c, d, len_step, l, u);
				}
				else {
					cos_concave(c, d, len_step, l, u);
				}
			}
			else {
				if (is_increase) {
					cos_concave(c, d, len_step, l, u);
				}
				else {
					cos_convex(c, d, len_step, l, u);
				}
			}
		}

		template<class T2> PwlBound<T2> cos(const PwlBound<T2>& bound) {
			T2 c = bound.get_min();
			T2 d = bound.get_max();
			int start = (int)floor(c / M_PI);
			int end = (int)ceil(d / M_PI);

			PwlFunc<T2> l;
			PwlFunc<T2> u;
			T2 len_step = (d - c) / bound.m_steps;
	
			if (end - start < 2) { 
				if (start % 2 == 0) {//decresing	
					cos_build_bound(c, d, len_step, false, l, u);
					return PwlBound<T2>(l(bound.m_upper), u(bound.m_lower), bound.m_steps);
				}
				else {//increasing
					cos_build_bound(c, d, len_step, true, l, u);
					return PwlBound<T2>(l(bound.m_lower), u(bound.m_upper), bound.m_steps);
				}
			}
			else {
				std::set<Segment<T2>> y_seg;
				for (int i = start; i < end; i++) {
					int last = end - 1;
					T2 a = i==start ? c : i*M_PI;
					T2 b = i==last ? d : (i + 1)*M_PI;
					bool inc = i % 2 == 0 ? false : true;
					y_seg.insert({ a, b, inc });
					cos_build_bound(a, b, len_step, inc, l, u);
				}
				PwlBound<T2> out_bound(l, u, bound.m_steps);
				return bound.get_copmosition(out_bound, y_seg);
			}
		}

		template<class T2> PwlBound<T2> sin(const PwlBound<T2>& bound) {
			return cos(bound - M_PI_2);
		}

		template<class T2> PwlBound<T2> max(const PwlBound<T2>& bound1, const PwlBound<T2>& bound2) {
			return PwlBound<T2>(max(bound1.m_lower, bound2.m_lower), max(bound1.m_upper, bound2.m_upper), bound1.m_steps);
		}

		template<class T2> PwlBound<T2> min(const PwlBound<T2>& bound1, const PwlBound<T2>& bound2) {
			return PwlBound<T2>(min(bound1.m_lower, bound2.m_lower), min(bound1.m_upper, bound2.m_upper), bound1.m_steps);
		}

		template<class T2> PwlBound<T2> ifThen(IntervalBool ib, const PwlBound<T2> &x, const PwlBound<T2> &y) {
			if (ib == IntervalBool::True)
				return x;
			else if (ib == IntervalBool::False)
				return y;
			else {
				if (eq(x.get_a(), y.get_a()) && eq(x.get_b(), y.get_b())) {
					return PwlBound<T2>(min(x.m_lower, y.m_lower), max(x.m_upper, y.m_upper), x.m_steps);
				}
				PwlFunc<T2> lower(x.m_lower), upper(x.m_upper);
				lower.insert(y.m_lower.m_points);
				upper.insert(y.m_upper.m_points);
				return PwlBound<T2>(lower, upper, x.m_steps);
			}			 
		}

		template<class T2> std::ostream& operator<<(std::ostream & out, const PwlBound<T2> x) {
			out << "lower bound: " << x.m_lower << std::endl;
			out << "upper bound: " << x.m_upper << std::endl;
			return out;
		}
	}
}

#endif
