#ifndef PWLFUNC_HPP
#define PWLFUNC_HPP

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include "utils.hpp"
#include "interval/enums.h"

namespace snowgoose {
	namespace pwl {

		template<class T> struct Parabola;
		template<class T> struct FracLinFunc;
		template<class T> class PwlBound;

		template<class T> class PwlFunc
		{
		public:
			PwlFunc(const std::set<Point<T>>& points) : m_points(points) {
			}
			PwlFunc() {}
			PwlFunc<T> operator+(const PwlFunc<T>& pwl) const;
			PwlFunc<T> operator-(const PwlFunc<T>& pwl) const;
			std::vector<Parabola<T>> operator*(const PwlFunc<T>& pwl) const;
			std::vector<FracLinFunc<T>> operator/(const PwlFunc<T>& pwl) const;
			PwlFunc<T> operator()(const PwlFunc<T>& pwl) const;

			IntervalBool operator < (const PwlFunc<T> &y) const;
			IntervalBool operator <= (const PwlFunc<T> &y) const;
			IntervalBool operator > (const PwlFunc<T> &y) const;
			IntervalBool operator >= (const PwlFunc<T> &y) const;
			bool operator==(const PwlFunc<T>& pwl) const;

			template<class T2> friend PwlFunc<T2> operator*(T2 x, const PwlFunc<T2>& y);
			template<class T2> friend PwlFunc<T2> operator*(const PwlFunc<T2>& x, T2 y);
			template<class T2> friend PwlFunc<T2> operator/(T2 x, const PwlFunc<T2>& y);
			template<class T2> friend PwlFunc<T2> operator/(const PwlFunc<T2>& x, T2 y);
			template<class T2> friend PwlFunc<T2> operator+(T2 x, const PwlFunc<T2>& y);
			template<class T2> friend PwlFunc<T2> operator+(const PwlFunc<T2>& x, T2 y);
			template<class T2> friend PwlFunc<T2> operator-(T2 x, const PwlFunc<T2>& y);
			template<class T2> friend PwlFunc<T2> operator-(const PwlFunc<T2>& x, T2 y);
			template<class T2> friend PwlFunc<T2> max(const std::vector<PwlFunc<T2>>& pwls);
			template<class T2> friend PwlFunc<T2> min(const std::vector<PwlFunc<T2>>& pwls);
			template<class T2> friend PwlFunc<T2> max(const PwlFunc<T2>& pwls1, const PwlFunc<T2>& pwls2);
			template<class T2> friend PwlFunc<T2> min(const PwlFunc<T2>& pwls1, const PwlFunc<T2>& pwls2);
			template<class T2> friend PwlBound<T2> ifThen(IntervalBool ib, const PwlBound<T2> &x, const PwlBound<T2> &y);
			/*returns the domain of pwl function*/
			T get_a() const { return m_points.begin()->x; }
			T get_b() const { return m_points.rbegin()->x; }
			T get_y(T x) const;
			/*returns min of pwl function*/
			T get_min() const;
			/*returns max of pwl function*/
			T get_max() const;
			std::vector<Point<T>> Points() { return std::vector<Point<T>>(m_points.begin(), m_points.end()); }
			PwlFunc<T> get_pwl_func(T a, T b) const;
			template<class T2> friend std::ostream& operator<<(std::ostream & out, const PwlFunc<T2> x);
			template<class T2> friend class PwlBound;
			void insert(std::initializer_list<Point<T>> points) { m_points.insert(points); }
		private:
			void insert(const std::set<Point<T>>& points) { m_points.insert(points.begin(), points.end()); }
			std::set<Point<T>> m_points;
		};

		template<class T> struct Parabola
		{
			abc<T> coeff;
			std::shared_ptr<Point<T>> start;
			std::shared_ptr<Point<T>> end;
		};

		template<class T> struct FracLinFunc
		{
			abcd<T> coeff;
			std::shared_ptr<Point<T>> start;
			std::shared_ptr<Point<T>> end;
		};

		template<class T> struct Segment
		{
			T start;
			T end;
			bool isIncrease;
			bool operator<(const Segment<T>& s) const {
				return ls(end, s.end);
			}
		};


		template<class T2> class CmpPwlFuncLeft {
		public:
			CmpPwlFuncLeft(T2 x1, T2 x2) : m_x1(x1), m_x2(x2) {}
			bool operator()(const std::pair<int, PwlFunc<T2>>& l, const std::pair<int, PwlFunc<T2>> & r)
			{
				PwlFunc<T2> f1 = l.second;
				PwlFunc<T2> f2 = r.second;
				T2 y1 = f1.get_y(m_x1);
				T2 y3 = f2.get_y(m_x1);
				if (eq(y1, y3)) {
					T2 y2 = f1.get_y(m_x2);
					T2 y4 = f2.get_y(m_x2);
					T2 coef12 = (y2 - y1) / (m_x2 - m_x1);
					T2 coef34 = (y4 - y3) / (m_x2 - m_x1);
					return ls(coef12, coef34); //slope
				}
				else
					return ls(y1, y3);
			}
		private:
			T2 m_x1;
			T2 m_x2;
		};

		template<class T2> class CmpPwlFuncRight {
		public:
			CmpPwlFuncRight(T2 x1, T2 x2) : m_x1(x1), m_x2(x2) {}
			bool operator()(const std::pair<int, PwlFunc<T2>>& l, const std::pair<int, PwlFunc<T2>> & r)
			{
				PwlFunc<T2> f1 = l.second;
				PwlFunc<T2> f2 = r.second;
				T2 y2 = f1.get_y(m_x2);
				T2 y4 = f2.get_y(m_x2);
				if (eq(y2, y4)) {
					T2 y1 = f1.get_y(m_x1);
					T2 y3 = f2.get_y(m_x1);
					T2 coef12 = (y2 - y1) / (m_x2 - m_x1);
					T2 coef34 = (y4 - y3) / (m_x2 - m_x1);
					return mo(coef12, coef34);//slope
				}
				else
					return ls(y2, y4);
			}
		private:
			T2 m_x1;
			T2 m_x2;
		};


		template<class T> PwlFunc<T> PwlFunc<T>::operator+(const PwlFunc<T>& pwl) const {
			std::set<Point<T>> all_x(m_points);
			all_x.insert(pwl.m_points.begin(), pwl.m_points.end());

			std::set<Point<T>> points;
			for (const Point<T>& point : all_x)
				points.insert({ point.x, this->get_y(point.x) + pwl.get_y(point.x) });

			return PwlFunc<T>(points);
		}

		template<class T> PwlFunc<T> PwlFunc<T>::operator-(const PwlFunc<T>& pwl) const {
			std::set<Point<T>> all_x(m_points);
			all_x.insert(pwl.m_points.begin(), pwl.m_points.end());

			std::set<Point<T>> points;
			for (const Point<T>& point : all_x)
				points.insert({ point.x, this->get_y(point.x) - pwl.get_y(point.x) });

			return PwlFunc<T>(points);
		}

		template<class T> std::vector<Parabola<T>> PwlFunc<T>::operator*(const PwlFunc<T>& pwl) const {
			std::set<Point<T>> all_x(m_points);
			all_x.insert(pwl.m_points.begin(), pwl.m_points.end());
			std::vector<Parabola<T>> parabols;

			auto last = std::prev(all_x.cend());
			std::shared_ptr<Point<T>> start;
			for (auto it = all_x.cbegin(); it != last; ++it) {
				auto next = std::next(it);
				Parabola<T> parabola;
				T y1 = this->get_y(it->x);
				T y2 = this->get_y(next->x);
				T y3 = pwl.get_y(it->x);
				T y4 = pwl.get_y(next->x);
				parabola.coeff = get_abc(it->x, y1, next->x, y2, y3, y4);
				parabola.start = start ? start : std::shared_ptr<Point<T>>(new Point<T>({ it->x, y1*y3 }));
				parabola.end = std::shared_ptr<Point<T>>(new Point<T>({ next->x, y2*y4 }));
				start = parabola.end;
				parabols.push_back(parabola);
			}
			return parabols;
		}

		template<class T> std::vector<FracLinFunc<T>> PwlFunc<T>::operator/(const PwlFunc<T>& pwl) const {
			std::set<Point<T>> all_x(m_points);
			all_x.insert(pwl.m_points.begin(), pwl.m_points.end());
			std::vector<FracLinFunc<T>> funcs;

			auto last = std::prev(all_x.cend());
			std::shared_ptr<Point<T>> start;
			for (auto it = all_x.cbegin(); it != last; ++it) {
				auto next = std::next(it);
				FracLinFunc<T> func;
				T y1 = this->get_y(it->x);
				T y2 = this->get_y(next->x);
				T y3 = pwl.get_y(it->x);
				T y4 = pwl.get_y(next->x);
				if (eq(y3, 0.0) || eq(y4, 0.0))
					throw std::invalid_argument("Exception in std::vector<FracLinFunc<T>> PwlFunc<T>::operator/(const PwlFunc<T>& pwl). Division by 0.0");
				func.coeff = get_abcd(it->x, y1, next->x, y2, y3, y4);
				func.start = start ? start : std::shared_ptr<Point<T>>(new Point<T>({ it->x, y1 / y3 }));
				func.end = std::shared_ptr<Point<T>>(new Point<T>({ next->x, y2 / y4 }));
				start = func.end;
				funcs.push_back(func);
			}
			return funcs;
		}

		template<class T> PwlFunc<T> PwlFunc<T>::operator()(const PwlFunc<T>& pwl) const {
			std::set<Point<T>> all_x(pwl.m_points);
			for (const Point<T>& out_point : m_points) {
				auto last = std::prev(pwl.m_points.cend());
				for (auto it = pwl.m_points.cbegin(); it != last; ++it) {
					auto next = std::next(it);
					if (out_point.x > std::min(it->y, next->y) && out_point.x < std::max(it->y, next->y)) { // stay between
						T x = get_intersect(it->x, it->y, next->x, next->y, out_point.x);
						all_x.insert({ x, out_point.x });
					}
				}
			}
			std::set<Point<T>> compound;
			for (const Point<T>& in_point : all_x) {
				compound.insert({ in_point.x, this->get_y(in_point.y) });
			}
			return PwlFunc<T>(compound);
		}

		template<class T> IntervalBool PwlFunc<T>::operator < (const PwlFunc<T> &y) const {
			return ls(this->get_max(), y.get_min()) ? IntervalBool::True : mo(this->get_min(), y.get_max()) ? IntervalBool::False : IntervalBool::Intermadiate;
		}

		template<class T> IntervalBool PwlFunc<T>::operator <= (const PwlFunc<T> &y) const {
			return lseq(this->get_max(), y.get_min()) ? IntervalBool::True : moeq(this->get_min(), y.get_max()) ? IntervalBool::False : IntervalBool::Intermadiate;
		}

		template<class T> IntervalBool PwlFunc<T>::operator > (const PwlFunc<T> &y) const {
			return mo(this->get_min(), y.get_max()) ? IntervalBool::True : ls(this->get_max(), y.get_min()) ? IntervalBool::False : IntervalBool::Intermadiate;
		}

		template<class T> IntervalBool PwlFunc<T>::operator >= (const PwlFunc<T> &y) const {
			return moeq(this->get_min(), y.get_max()) ? IntervalBool::True : lseq(this->get_max(), y.get_min()) ? IntervalBool::False : IntervalBool::Intermadiate;
		}

		template<class T> inline bool PwlFunc<T>::operator==(const PwlFunc<T>& pwl) const{
			return this->m_points == pwl.m_points;
		}


		template<class T> T PwlFunc<T>::get_y(T x) const {
			auto it = m_points.find({ x, 0 });
			if (it != m_points.end()) {
				return it->y;
			}
			else {
				auto p2 = m_points.upper_bound({ x, 0.0 });
				if (p2 == m_points.begin() || p2 == m_points.end())
					throw std::invalid_argument("PwlFunc<T>::get_y(x). Invalid value x.");
				auto p1 = p2;
				--p1;
				return pwl::get_y(p1->x, p1->y, p2->x, p2->y, x);
			}
		}

		template<class T> T PwlFunc<T>::get_min() const {
			return (*std::min_element(m_points.begin(), m_points.end(), [](const Point<T>& p1, const Point<T>& p2) { return ls(p1.y, p2.y); })).y;
		}

		template<class T> T PwlFunc<T>::get_max() const {
			return (*std::max_element(m_points.begin(), m_points.end(), [](const Point<T>& p1, const Point<T>& p2) { return ls(p1.y, p2.y); })).y;
		}

		template<class T2> PwlFunc<T2> operator*(T2 x, const PwlFunc<T2>& f) {
			std::set<Point<T2>> points;
			for_each(f.m_points.begin(), f.m_points.end(), [&points, x](const Point<T2>& point) { points.insert({ point.x, point.y * x }); });
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> operator*(const PwlFunc<T2>& x, T2 y) {
			return y*x;
		}

		template<class T2> PwlFunc<T2> operator/(T2 x, const PwlFunc<T2>& f) {
			std::set<Point<T2>> points;
			for (const Point<T2>& point : f.m_points) {
				if (eq(point.y, 0.0)) {
					throw std::invalid_argument("Exception in PwlFunc<T2> operator/(T2 x, const PwlFunc<T2>& f). Division by 0");
				}
				points.insert({ point.x,  x / point.y });
			}
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> operator/(const PwlFunc<T2>& f, T2 y) {
			if (eq(y, 0.0)) {
				throw std::invalid_argument("Exception in PwlFunc<T2> operator/(const PwlFunc<T2>& x, T2 y). Division by 0");
			}
			std::set<Point<T2>> points;
			for_each(f.m_points.begin(), f.m_points.end(), [&points, y](const Point<T2>& point) { points.insert({ point.x, point.y / y }); });
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> operator+(T2 x, const PwlFunc<T2>& f) {
			std::set<Point<T2>> points;
			for_each(f.m_points.begin(), f.m_points.end(), [&points, x](const Point<T2>& point) { points.insert({ point.x, point.y + x }); });
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> operator+(const PwlFunc<T2>& f, T2 y) {
			return y + f;
		}

		template<class T2> PwlFunc<T2> operator-(T2 x, const PwlFunc<T2>& f) {
			std::set<Point<T2>> points;
			for_each(f.m_points.begin(), f.m_points.end(), [&points, x](const Point<T2>& point) { points.insert({ point.x, x - point.y }); });
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> operator-(const PwlFunc<T2>& f, T2 y) {
			std::set<Point<T2>> points;
			for_each(f.m_points.begin(), f.m_points.end(), [&points, y](const Point<T2>& point) { points.insert({ point.x, point.y - y }); });
			return PwlFunc<T2>(points);
		}

		template<class T2> PwlFunc<T2> max(const std::vector<PwlFunc<T2>>& pwls)
		{
			std::set<Point<T2>> all_x;
			std::unordered_map<int, PwlFunc<T2>> map;
			int i = 0;
			for (auto& f : pwls) {
				all_x.insert(f.m_points.begin(), f.m_points.end());
				map[i++] = f;
			}

			std::set<Point<T2>> rezult;
			auto last = std::prev(all_x.cend());
			for (auto it = all_x.cbegin(); it != last; ++it) {
				auto next = std::next(it);
				std::unordered_map<int, PwlFunc<T2>> temp(map);

				std::pair<int, PwlFunc<T2>> max_curr = *std::max_element(temp.begin(), temp.end(), CmpPwlFuncLeft<T2>(it->x, next->x));
				rezult.insert({ it->x, max_curr.second.get_y(it->x) });
				std::pair<int, PwlFunc<T2>> max_last = *std::max_element(temp.begin(), temp.end(), CmpPwlFuncRight<T2>(it->x, next->x));
				rezult.insert({ next->x , max_last.second.get_y(next->x) });

				Point<T2> most_left = *next;
				while (max_curr.first != max_last.first) {
					auto max_next = max_curr;
					for (auto &item : temp) {
						Point<T2> p;
						if (get_intersect(it->x, max_curr.second.get_y(it->x), next->x, max_curr.second.get_y(next->x), item.second.get_y(it->x), item.second.get_y(next->x), p)) {
							if (ls(p.x, most_left.x)) {
								max_next = item;
								most_left = p;
							}
						}
					}
					rezult.insert(most_left);
					temp.erase(max_next.first);
					max_curr = max_next;
					most_left = *next;
				}
			}
			return rezult;
		}

		template<class T2> PwlFunc<T2> min(const std::vector<PwlFunc<T2>>& pwls)
		{
			std::set<Point<T2>> all_x;
			std::unordered_map<int, PwlFunc<T2>> map;
			int i = 0;
			for (auto& f : pwls) {
				all_x.insert(f.m_points.begin(), f.m_points.end());
				map[i++] = f;
			}

			std::set<Point<T2>> rezult;
			auto last = std::prev(all_x.cend());
			for (auto it = all_x.cbegin(); it != last; ++it) {
				auto next = std::next(it);
				std::unordered_map<int, PwlFunc<T2>> temp(map);

				std::pair<int, PwlFunc<T2>> min_curr = *std::min_element(temp.begin(), temp.end(), CmpPwlFuncLeft<T2>(it->x, next->x));
				rezult.insert({ it->x, min_curr.second.get_y(it->x) });
				std::pair<int, PwlFunc<T2>> min_last = *std::min_element(temp.begin(), temp.end(), CmpPwlFuncRight<T2>(it->x, next->x));
				rezult.insert({ next->x , min_last.second.get_y(next->x) });
				temp.erase(min_curr.first);

				Point<T2> most_left = *next;
				while (min_curr.first != min_last.first) {
					auto min_next = min_curr;
					for (auto & item : temp) {
						Point<T2> p;
						if (get_intersect(it->x, min_curr.second.get_y(it->x), next->x, min_curr.second.get_y(next->x), item.second.get_y(it->x), item.second.get_y(next->x), p)) {
							if (ls(p.x, most_left.x)) {
								min_next = item;
								most_left = p;
							}
						}
					}
					rezult.insert(most_left);
					temp.erase(min_next.first);
					min_curr = min_next;
					most_left = *next;
				}
			}
			return rezult;
		}

		template<class T2> PwlFunc<T2> max(const PwlFunc<T2>& pwls1, const PwlFunc<T2>& pwls2) {
			std::vector<PwlFunc<T2>> v({ pwls1, pwls2 });
			return max(v);
		}
		template<class T2> PwlFunc<T2> min(const PwlFunc<T2>& pwls1, const PwlFunc<T2>& pwls2) {
			std::vector<PwlFunc<T2>> v({ pwls1, pwls2 });
			return min(v);
		}

		template<class T> PwlFunc<T> PwlFunc<T>::get_pwl_func(T a, T b) const {
			if (ls(a, this->get_a()) || mo(b, this->get_b()))
				throw std::invalid_argument("Exception in PwlFunc<T> PwlFunc<T>::get_pwl_func(T a, T b). Arguments are out of range");
			std::set<Point<T>> points;
			points.insert({ {a, this->get_y(a) }, {b, this->get_y(b)} });
			auto start = m_points.upper_bound({ a, 0.0 });
			auto end = m_points.lower_bound({ b, 0.0 });
			points.insert(start, end);
			return PwlFunc<T>(points);
		}

		template<class T2> std::ostream& operator<<(std::ostream & out, const PwlFunc<T2> x) {
			for (const Point<T2>& point : x.m_points) {
				out << "(" << point.x << "; " << point.y << ") ";
			}
			out << std::endl;
			return out;
		}

		template<class T> std::ostream& operator<<(std::ostream & out, const Parabola<T>& p)
		{
			out << "x1 = " << p.start->x << " y1 = " << p.start->y << " x2 = " << p.end->x << " y2 = " << p.end->y << " a = " << p.coeff.a << " b = " << p.coeff.b << " c = " << p.coeff.c << std::endl;
			return out;
		}

		template<class T> std::ostream& operator<<(std::ostream & out, const FracLinFunc<T>& p)
		{
			out << "x1 = " << p.start->x << " y1 = " << p.start->y << " x2 = " << p.end->x << " y2 = " << p.end->y << " a = " << p.coeff.a << " b = " << p.coeff.b << " c = " << p.coeff.c << " d = " << p.coeff.d << std::endl;
			return out;
		}
	}
}

#endif
		



