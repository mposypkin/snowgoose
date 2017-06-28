/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   valder.hpp
 * Author: alusov
 *
 * Created on June 22, 2017, 11:36 AM
 */

#ifndef INTERVALDER_HPP
#define INTERVALDER_HPP

#include <iostream> 
#include <cmath>
#include "interval/interval_air.hpp"
#include "grad.hpp"

using namespace snowgoose::interval;

namespace snowgoose {
  namespace derivative {
    template<class T> class IntervalDer
    {
        public:
            IntervalDer(const Interval<T> &val, const Grad<Interval<T>> &der) : m_val(val), m_der(der)
            {
            }
            IntervalDer(T val, const Grad<Interval<T>> &der) : m_val(val), m_der(der)
            {
            }

            IntervalDer operator+(const IntervalDer &y) const;
            IntervalDer operator-(const IntervalDer &y) const;
            IntervalDer operator*(const IntervalDer &y) const;
            IntervalDer operator/(const IntervalDer &y) const;
            
            IntervalDer operator^(int exp) const;
            
            template<class T2> friend IntervalDer<T2> operator+(T2 t, const IntervalDer<T2> &y);
            template<class T2> friend IntervalDer<T2> operator-(T2 t, const IntervalDer<T2> &y);
            template<class T2> friend IntervalDer<T2> operator*(T2 t, const IntervalDer<T2> &y);
            template<class T2> friend IntervalDer<T2> operator/(T2 t, const IntervalDer<T2> &y);
            
            template<class T2> friend IntervalDer<T2> operator+(const IntervalDer<T2> &y, T2 t);
            template<class T2> friend IntervalDer<T2> operator-(const IntervalDer<T2> &y, T2 t);
            template<class T2> friend IntervalDer<T2> operator*(const IntervalDer<T2> &y, T2 t);
            template<class T2> friend IntervalDer<T2> operator/(const IntervalDer<T2> &y, T2 t);
                       
            template<class T2> friend IntervalDer<T2>sqr(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2>sqrt(const IntervalDer<T2> &x);
            
            template<class T2> friend IntervalDer<T2>sin(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2>cos(const IntervalDer<T2>&x);
            
            template<class T2> friend IntervalDer<T2>tg(const IntervalDer<T2> &x);
            
            template<class T2> friend IntervalDer<T2>asin(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2>acos(const IntervalDer<T2> &x);
            
            template<class T2> friend IntervalDer<T2>ln(const IntervalDer<T2> &x);
            template<class T2> friend IntervalDer<T2>exp(const IntervalDer<T2> &x);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const IntervalDer<T2> &x);
            
        private:
            Interval<T> m_val;
            Grad<Interval<T>> m_der;
    };
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator+(const IntervalDer &y) const
    {
        return IntervalDer(m_val + y.m_val, m_der + y.m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator-(const IntervalDer &y) const
    {
        return IntervalDer(m_val - y.m_val, m_der - y.m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator*(const IntervalDer &y) const
    {
        return IntervalDer(m_val * y.m_val, m_val * y.m_der + y.m_val * m_der);
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator/(const IntervalDer &y) const
    {
        return IntervalDer(m_val / y.m_val, (m_der * y.m_val - m_val * y.m_der)/(y.m_val * y.m_val) );
    }
    
    template<class T> IntervalDer<T> IntervalDer<T>::operator^(int exp) const
    {
        return IntervalDer(m_val^ exp, exp * (m_val ^ (exp - 1)) * m_der);//error control is needed
    }
    
    template<class T2> IntervalDer<T2>operator+(T2 t, const IntervalDer<T2>&y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) + y;
    }
    
    template<class T2> IntervalDer<T2>operator-(T2 t, const IntervalDer<T2>&y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) - y;
    }
    
    template<class T2> IntervalDer<T2>operator*(T2 t, const IntervalDer<T2>&y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) * y;
    }
    
    template<class T2> IntervalDer<T2>operator/(T2 t, const IntervalDer<T2>&y)
    {
        return IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0)) / y;
    }
    
    template<class T2> IntervalDer<T2>operator+(const IntervalDer<T2>&y, T2 t)
    {
        return t+y;
    }
    
    template<class T2> IntervalDer<T2>operator-(const IntervalDer<T2>&y, T2 t)
    {
        return y-IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0));
    }
    
    template<class T2> IntervalDer<T2>operator*(const IntervalDer<T2>&y, T2 t)
    {
        return t*y;
    }
    
    template<class T2> IntervalDer<T2>operator/(const IntervalDer<T2>&y, T2 t)
    {
        return y/IntervalDer<T2>(t, Grad<Interval<T2>>(y.m_der.size(), 0.0));
    }
    
    template<class T2> IntervalDer<T2>sqr(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(x.m_val * x.m_val, 2.0 * x.m_val * x.m_der);
    }
    
    template<class T2> IntervalDer<T2>sqrt(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::sqrt(x.m_val), x.m_der/(2.0 * interval::sqrt(x.m_val)));//error control is needed
    }
    
    template<class T2> IntervalDer<T2>sin(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::sin(x.m_val), interval::cos(x.m_val) * x.m_der);
    }
    
    template<class T2> IntervalDer<T2>cos(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::cos(x.m_val), -1.0 * interval::sin(x.m_val) * x.m_der);
    }
      
    template<class T2> IntervalDer<T2>tg(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::tg(x.m_val), x.m_der/(interval::cos(x.m_val) * interval::cos(x.m_val)));//error control is needed
    }
    
    template<class T2> IntervalDer<T2>asin(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::asin(x.m_val), x.m_der/interval::sqrt(1- x.m_val*x.m_val ));//error control is needed
    }
    
    template<class T2> IntervalDer<T2>acos(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::acos(x.m_val), -x.m_der/interval::sqrt(1- x.m_val*x.m_val ));//error control is needed
    }
    
    template<class T2> IntervalDer<T2>ln(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::ln(x.m_val), x.m_der/x.m_val);//error control is needed
    }
    
    template<class T2> IntervalDer<T2>exp(const IntervalDer<T2>&x)
    {
        return IntervalDer<T2>(interval::exp(x.m_val), x.m_der * interval::exp(x.m_val));//optimization
    }
    
    template<class T2> std::ostream& operator<<(std::ostream & out, const IntervalDer<T2> &x)
    {
        return std::cout << "val: " << x.m_val << " der: " << x.m_der << '\n';
    }
   
  }
}

#endif /* VALDER_HPP */

