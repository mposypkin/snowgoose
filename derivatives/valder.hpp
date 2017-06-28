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

#ifndef VALDER_HPP
#define VALDER_HPP

#include <iostream> 
#include <cmath>
#include "interval/interval_air.hpp"
#include "grad.hpp"
#include "helpfunc.hpp"

using namespace snowgoose::interval;

namespace snowgoose {
  namespace derivative {
    template<class T> class ValDer
    {
        public:
            ValDer(T val, const Grad<T> &der) : m_val(val), m_der(der)
            {
            }

            ValDer operator+(const ValDer &y) const;
            ValDer operator-(const ValDer &y) const;
            ValDer operator*(const ValDer &y) const;
            ValDer operator/(const ValDer &y) const;            
            ValDer operator^(int exp) const;
            
            template<class T2, class T3> friend ValDer<T2> operator+(T3 t, const ValDer<T2> &y);
            template<class T2, class T3> friend ValDer<T2> operator-(T3 t, const ValDer<T2> &y);
            template<class T2, class T3> friend ValDer<T2> operator*(T3 t, const ValDer<T2> &y);
            template<class T2, class T3> friend ValDer<T2> operator/(T3 t, const ValDer<T2> &y);           
            template<class T2, class T3> friend ValDer<T2> operator+(const ValDer<T2> &y, T3 t);
            template<class T2, class T3> friend ValDer<T2> operator-(const ValDer<T2> &y, T3 t);
            template<class T2, class T3> friend ValDer<T2> operator*(const ValDer<T2> &y, T3 t);
            template<class T2, class T3> friend ValDer<T2> operator/(const ValDer<T2> &y, T3 t);
                       
            template<class T2> friend ValDer<T2>sqr(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2>sqrt(const ValDer<T2> &x);            
            template<class T2> friend ValDer<T2>sin(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2>cos(const ValDer<T2>&x);            
            template<class T2> friend ValDer<T2>tg(const ValDer<T2> &x);           
            template<class T2> friend ValDer<T2>asin(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2>acos(const ValDer<T2> &x);           
            template<class T2> friend ValDer<T2>ln(const ValDer<T2> &x);
            template<class T2> friend ValDer<T2>exp(const ValDer<T2> &x);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const ValDer<T2> &x);
            
        private:
            T m_val;
            Grad<T> m_der;
    };
    
    template<class T> ValDer<T> ValDer<T>::operator+(const ValDer &y) const
    {
        return ValDer(m_val + y.m_val, m_der + y.m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator-(const ValDer &y) const
    {
        return ValDer(m_val - y.m_val, m_der - y.m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator*(const ValDer &y) const
    {
        return ValDer(m_val * y.m_val, m_val * y.m_der + y.m_val * m_der);
    }
    
    template<class T> ValDer<T> ValDer<T>::operator/(const ValDer &y) const
    {
        return ValDer(m_val / y.m_val, (m_der * y.m_val - m_val * y.m_der)/(y.m_val * y.m_val) );
    }
    
    template<class T> ValDer<T> ValDer<T>::operator^(int exp) const
    {
        return ValDer(pow_(m_val, exp), (double)exp * pow_(m_val, exp - 1.0) * m_der);//error control is needed
    }
    
    template<class T2, class T3> ValDer<T2>operator+(T3 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) + y;
    }
    
    template<class T2, class T3> ValDer<T2>operator-(T3 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) - y;
    }
    
    template<class T2, class T3> ValDer<T2>operator*(T3 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) * y;
    }
    
    template<class T2, class T3> ValDer<T2>operator/(T3 t, const ValDer<T2>&y)
    {
        return ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0)) / y;;
    }
    
    template<class T2, class T3> ValDer<T2>operator+(const ValDer<T2>&y, T3 t)
    {
        return t+y;
    }
    
    template<class T2, class T3> ValDer<T2>operator-(const ValDer<T2>&y, T3 t)
    {
        return y-ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0));
    }
    
    template<class T2, class T3> ValDer<T2>operator*(const ValDer<T2>&y, T3 t)
    {
        return t*y;
    }
    
    template<class T2, class T3> ValDer<T2>operator/(const ValDer<T2>&y, T3 t)
    {
        return y/ValDer<T2>(t, Grad<T2>(y.m_der.size(), 0.0));
    }
    
    template<class T2> ValDer<T2>sqr(const ValDer<T2>&x)
    {
        return ValDer<T2>(x.m_val * x.m_val, 2.0 * x.m_val * x.m_der);
    }
    
    template<class T2> ValDer<T2>sqrt(const ValDer<T2>&x)
    {
        return ValDer<T2>(sqrt_(x.m_val), x.m_der/(2.0 * sqrt_(x.m_val)));//error control is needed
    }
    
    template<class T2> ValDer<T2>sin(const ValDer<T2>&x)
    {
        return ValDer<T2>(sin_(x.m_val), cos_(x.m_val) * x.m_der);
    }
    
    template<class T2> ValDer<T2>cos(const ValDer<T2>&x)
    {
        return ValDer<T2>(cos_(x.m_val), -1.0 * sin_(x.m_val) * x.m_der);
    }
      
    template<class T2> ValDer<T2>tg(const ValDer<T2>&x)
    {
        return ValDer<T2>(tg_(x.m_val), x.m_der/(cos_(x.m_val) * cos_(x.m_val)));//error control is needed
    }
    
    template<class T2> ValDer<T2>asin(const ValDer<T2>&x)
    {
        return ValDer<T2>(asin_(x.m_val), x.m_der/sqrt_(1.0 - x.m_val*x.m_val ));//error control is needed
    }
    
    template<class T2> ValDer<T2>acos(const ValDer<T2>&x)
    {
        return ValDer<T2>(acos_(x.m_val), -1.0 * x.m_der/sqrt_(1.0 - x.m_val*x.m_val ));//error control is needed
    }
    
    template<class T2> ValDer<T2>ln(const ValDer<T2>&x)
    {
        return ValDer<T2>(ln_(x.m_val), x.m_der/x.m_val);//error control is needed
    }
    
    template<class T2> ValDer<T2>exp(const ValDer<T2>&x)
    {
        return ValDer<T2>(exp_(x.m_val), x.m_der * exp_(x.m_val));//optimization
    }
    
    template<class T2> std::ostream& operator<<(std::ostream & out, const ValDer<T2> &x)
    {
        return std::cout << "val: " << x.m_val << " der: " << x.m_der << '\n';
    }
   
  }
}

#endif /* VALDER_HPP */

