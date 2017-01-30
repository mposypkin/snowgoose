/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   interval_air.hpp
 * Author: alusov
 *
 * Created on January 28, 2017, 12:07 PM
 */

#ifndef INTERVAL_AIR_HPP
#define INTERVAL_AIR_HPP
#include <iostream> 
#include <limits>
#include <math.h>
#include  <exception>
#include <unordered_map>
#include <memory>
#include "common/utilmacro.hpp"
#include "common/sgerrcheck.hpp"
#include "calcvalues.hpp"

namespace snowgoose { namespace interval {

    template<class T> class Interval
    {
        public:
            Interval(const T &min, const T &max)
            {
                this->m_lb = min;
                this->m_rb = max;
            }

            Interval(const T &t)
            {
                this->m_lb = t;
                this->m_rb = t;
            }

            Interval operator+(const Interval &y) const
            {
                Interval z(m_lb + y.m_lb, m_rb + y.m_rb);
                return z;
            }

            Interval operator-(const Interval &y) const
            {
                Interval z(m_lb - y.m_rb, m_rb - y.m_lb);
                return z;
            }

            Interval operator^(int n) const
            {
                T t1 = 1.0, t2 = 1.0, lb, rb;
                for (int i = 1; i <= n; i++) {
                    t1 *= m_lb;
                    t2 *= m_rb;
                }
                if (n % 2) {
                    lb = t1;
                    rb = t2;
                } else {
                    if ((m_lb <= 0.0) && (m_rb >= 0.0))
                        lb = 0.0;
                    else
                        lb = SGMIN(t1, t2);
                    rb = SGMAX(t1, t2);
                }
                return Interval(lb, rb);
            }

            Interval operator*(const Interval &y) const
            {
                T t1, t2, t3, t4, lb, rb;
                t1 = m_lb * y.m_lb;
                lb = t1;
                rb = t1;

                t2 = m_lb * y.m_rb;
                lb = SGMIN(lb, t2);
                rb = SGMAX(rb, t2);

                t3 = m_rb * y.m_lb;
                lb = SGMIN(lb, t3);
                rb = SGMAX(rb, t3);

                t4 = m_rb * y.m_rb;
                lb = SGMIN(lb, t4);
                rb = SGMAX(rb, t4); 

                return Interval(lb, rb);
            }

            Interval operator/(const Interval &y) const
            {
                if(y.m_lb <= 0.0 && y.m_rb >=0.0)
                    throw std::invalid_argument("Invalid operation. Interval is divided by interval that includes 0.0.");

                T t1, t2, t3, t4, lb, rb;
                t1 = m_lb / y.m_lb;
                lb = t1;
                rb = t1;

                t2 = m_lb/y.m_rb;
                lb = SGMIN(lb, t2);
                rb = SGMAX(rb, t2);

                t3 = m_rb/y.m_lb;
                lb = SGMIN(lb, t3);
                rb = SGMAX(rb, t3);

                t4 = m_rb/y.m_rb;
                lb = SGMIN(lb, t4);
                rb = SGMAX(rb, t4); 

                return Interval(lb, rb);
            }

            friend Interval operator*(T t, const Interval &y)
            {
                return Interval(t) * y;
            }
            friend Interval operator+(T t, const Interval &y)
            {
                return Interval(t) + y;
            }
            friend Interval operator-(T t, const Interval &y)
            {
                return Interval(t) - y;
            }
            friend Interval operator/(T t, const Interval &y)
            {
                return Interval(t) / y;
            }
            friend Interval sqr(const Interval<T> &x)
            {
                return x ^ 2;
            }
            friend Interval sin(const Interval<T> &x)
            {
                int k, l;
                bool p2 = false;
                bool p32 = false;
                T lb, rb;

                k = (int)ceil(x.m_lb * M_2_PI);
                l = (int)floor(x.m_rb * M_2_PI);

                for (int i = k; i <= l; i++) {
                    int test = i % 4;
                    if ((i % 4) == 1 || (i % 4) == -3)
                        p2 = true;
                    else if ((i % 4) == 3 || (i % 4) == -1)
                        p32 = true;
                    if (p2 && p32)
                        break;
                }
                if (p2)
                    rb = 1.0;
                else
                    rb = SGMAX(::sin(x.m_lb), ::sin(x.m_rb));
                if (p32)
                    lb = -1.0;
                else
                    lb = SGMIN(::sin(x.m_lb), ::sin(x.m_rb));

                return Interval(lb, rb);
            }
            friend Interval cos(const Interval<T> &x)
            {
                T t1, t2;
                t1 = x.m_lb + M_PI_2;
                t2 = x.m_rb + M_PI_2;
                Interval t(t1, t2);
                return sin(t);
            }
            friend Interval exp(const Interval<T> &x)
            {
                T t1 = ::exp(x.m_lb);
                T t2 = ::exp(x.m_rb);
                T lb = SGMIN(t1, t2);
                T rb = SGMAX(t1, t2);
                return Interval(lb, rb);
            }
            friend Interval sqrt(const Interval<T> &x)
            {
                if (x.m_lb < 0.0 || x.m_rb < 0.0)
                    throw std::invalid_argument("The function sqrt is not define for negative numbers");
                T a = ::sqrt(x.m_lb);
                T b = ::sqrt(x.m_rb);
                T lb = SGMIN(a, b);
                T rb = SGMAX(a, b);
                return Interval(lb, rb);
            }
            friend Interval abs(const Interval<T> &x)
            {
                T t1 = ::abs(x.m_lb);
                T t2 = ::abs(x.m_rb);
                T lb = (SGMIN(x.m_lb, x.m_rb) <= 0.0 && SGMAX(x.m_lb, x.m_rb) >= 0) ? 0.0 : SGMIN(t1, t2);
                T rb = SGMAX(t1, t2);
                return Interval(lb, rb);
            }
            
            friend Interval ln(const Interval<T> &x)
            {
                if (x.m_lb < 0.0 || x.m_rb < 0.0)
                    throw std::invalid_argument("The function ln is not define for negative numbers");
                T a = ::log(x.m_lb);
                T b = ::log(x.m_rb);
                T lb = SGMIN(a, b);
                T rb = SGMAX(a, b);
                return Interval(lb, rb);
            }

            friend Interval log(const Interval<T> &x, T base)
            {
                if (x.m_lb <= 0.0 || x.m_rb <= 0.0)
                    throw std::invalid_argument("The function log is not define for negative numbers and 0.0");
                if(base <= 0.0 || base == 1.0)
                    throw std::invalid_argument("The function log is not define for negative base and if base equals 1.0");

                T a = ::log(x.m_lb)/::log(base);
                T b = ::log(x.m_rb)/::log(base);
                T lb = SGMIN(a, b);
                T rb = SGMAX(a, b);
                return Interval(lb, rb);
            }
            
            friend std::ostream& operator<<(std::ostream & out, const Interval x)
            {
                return std::cout << "lb: " << x.m_lb << " rb: " << x.m_rb << '\n';
            }
            
            T lb()const {return m_lb;}
            T rb()const {return m_rb;}

        private:
            T m_lb;
            T m_rb;
    };            
        
    }
}


#endif /* INTERVAL_AIR_HPP */

