#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <stdexcept>
#include <cmath>
#define ACCURACY 0.00001

namespace snowgoose {
	namespace pwl {

		template<class T> bool eq(T x, T y) {
			return std::abs(x - y) <= ACCURACY;
		}

		template<class T> bool ls(T x, T y) {
			return y - x > ACCURACY;
		}

		template<class T> bool mo(T x, T y) {
			return x - y > ACCURACY;
		}

		template<class T> bool lseq(T x, T y) {
			return x - y <= ACCURACY;
		}

		template<class T> bool moeq(T x, T y) {
			return y - x <= ACCURACY;
		}

		/*A point with x,y coordinates*/
		template<class T> struct Point
		{
			T x;
			T y;
			bool operator<(const Point<T>& p) const {
				return ls(x, p.x);
			}
			bool operator==(const Point<T>& p) const {
				return eq(x, p.x);
			}
		};

		/*a,b,c coefficient of a parabola*/
		template<class T> struct abc
		{
			T a;
			T b;
			T c;
		};

		/*a,b,c coefficient of a fractional-linear function*/
		template<class T> struct abcd
		{
			T a;
			T b;
			T c;
			T d;
		};

		template <class T> T get_y(T x1, T y1, T x2, T y2, T x)
		{
			return y1 + (y2 - y1)*(x - x1) / (x2 - x1);
		}

		/*y=ax+b; y'=cx+d; yy'=acx^2+(ad+bc)x+bd*/
		template <class T> abc<T> get_abc(T x1, T y1, T x2, T y2, T y3, T y4)
		{
			T a = eq(y2, y1) ? 0.0 : (y2 - y1) / (x2 - x1);
			T b = y1 - a*x1;
			T c = eq(y4, y3) ? 0.0 : (y4 - y3) / (x2 - x1);
			T d = y3 - c*x1;
			return{ a*c, a*d + b*c, b*d };
		}

		template <class T> abcd<T> get_abcd(T x1, T y1, T x2, T y2, T y3, T y4)
		{
			T a = eq(y2, y1) ? 0.0 : (y2 - y1) / (x2 - x1);
			T b = y1 - a*x1;
			T c = eq(y4, y3) ? 0.0 : (y4 - y3) / (x2 - x1);
			T d = y3 - c*x1;

			return{ a, b , c, d };
		}

		template <class T> T get_intersect(T x1, T y1, T x2, T y2, T c)
		{
			T a = (y2 - y1) / (x2 - x1);
			T b = y1 - a*x1;
			return (c - b) / a;
		}

		template <class T> bool get_intersect(T x1, T y1, T x2, T y2, T y3, T y4, Point<T>& p)
		{
			if ((mo(y1, y3) && ls(y2, y4)) || (ls(y1, y3) && mo(y2, y4))) {
				T a = (y2 - y1) / (x2 - x1);
				T b = y1 - a*x1;
				T c = (y4 - y3) / (x2 - x1);
				T d = y3 - c*x1;
				p.x = (d - b) / (a - c);
				p.y = a*p.x + b;
				return true;
			}
			else
				return false;
		}

		template <class T> T get_val(T a, T b, T c, T x) {
			return a*x*x + b*x + c;
		}

		template<class T> int get_quarter(const Point<T> & p, const Point<T>& c) {
			if (mo(p.x, c.x) && mo(p.y, c.y)) {
				return 1;
			}
			else if (mo(p.x, c.x) && ls(p.y, c.y)) {
				return 2;
			}
			else if (ls(p.x, c.x) && ls(p.y, c.y)) {
				return 3;
			}
			else if (ls(p.x, c.x) && mo(p.y, c.y)) {
				return 4;
			}
			else {
				throw std::invalid_argument("Exception in int get_quarter(Point<T> p, Point<T> c). Invalid point.");
			}
		}
	}
}

#endif
