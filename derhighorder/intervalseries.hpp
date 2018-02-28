/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   series.hpp
 * Author: alusov
 *
 * Created on June 22, 2017, 11:36 AM
 */

#ifndef INTERVAL_SERIES_HPP
#define INTERVAL_SERIES_HPP

#include <iostream> 
#include <cmath>
#include <vector>
#include "matrixod.hpp"
#include "interval/enums.h"
#include "interval/interval_air.hpp"

using namespace snowgoose::interval;

namespace snowgoose {
  namespace derhighorder {
    template<class T> class IntervalSeries
    {
        public:
            IntervalSeries(const Interval<T> &val, const Interval<T> &der, int order) : m_coef(order + 1, 0.0)
            {
		m_coef[0] = val;
                m_coef[1] = der;
	    }
            IntervalSeries(){}
            using value_type = Interval<T>;
            IntervalSeries operator+(const IntervalSeries &y) const;
            IntervalSeries operator-(const IntervalSeries &y) const;
            IntervalSeries operator*(const IntervalSeries &y) const;
            IntervalSeries operator/(const IntervalSeries &y) const;            
            IntervalSeries operator^(int exp) const;
            IntervalSeries operator^(const Interval<T> &exp) const;
            
            template<class T2, class T3> friend IntervalSeries<T2> operator+(T3 t, const IntervalSeries<T2> &y);
            template<class T2, class T3> friend IntervalSeries<T2> operator-(T3 t, const IntervalSeries<T2> &y);
            template<class T2, class T3> friend IntervalSeries<T2> operator*(T3 t, const IntervalSeries<T2> &y);
            template<class T2, class T3> friend IntervalSeries<T2> operator/(T3 t, const IntervalSeries<T2> &y);           
            template<class T2, class T3> friend IntervalSeries<T2> operator+(const IntervalSeries<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalSeries<T2> operator-(const IntervalSeries<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalSeries<T2> operator*(const IntervalSeries<T2> &y, T3 t);
            template<class T2, class T3> friend IntervalSeries<T2> operator/(const IntervalSeries<T2> &y, T3 t);
                       
            template<class T2> friend IntervalSeries<T2> sqr(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> sqrt(const IntervalSeries<T2> &x);            
            template<class T2> friend IntervalSeries<T2> sin(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> cos(const IntervalSeries<T2>&x);            
            template<class T2> friend IntervalSeries<T2> tg(const IntervalSeries<T2> &x);  
            template<class T2> friend IntervalSeries<T2> ctg(const IntervalSeries<T2> &x); 
            template<class T2> friend IntervalSeries<T2> asin(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> acos(const IntervalSeries<T2> &x); 
            template<class T2> friend IntervalSeries<T2> atg(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> actg(const IntervalSeries<T2> &x); 
            template<class T2> friend IntervalSeries<T2> ln(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> log(const IntervalSeries<T2> &x, double base);
            template<class T2> friend IntervalSeries<T2> exp(const IntervalSeries<T2> &x);
            template<class T2, class T3> friend IntervalSeries<T2> pow(const T3 &base, const IntervalSeries<T2> &exp);
            template<class T2> friend IntervalSeries<T2> abs(const IntervalSeries<T2> &x);
            template<class T2> friend IntervalSeries<T2> min(const IntervalSeries<T2> &x, const IntervalSeries<T2> &y);
            template<class T2> friend IntervalSeries<T2> max(const IntervalSeries<T2> &x, const IntervalSeries<T2> &y);
            template<class T2> friend IntervalSeries<T2> ifThen(IntervalBool condition, const IntervalSeries<T2> &x, const IntervalSeries<T2> &y);
            template<class T2> friend IntervalBool cond(Conditions condition, const IntervalSeries<T2> &x, const IntervalSeries<T2> &y);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const IntervalSeries<T2> &x);
            Interval<T> value() const { return m_coef.item(0); }
            Interval<T> der(int order) const { return m_coef.item(order) * (double)fact(order); }
            explicit operator T(){ return m_coef.item(0); }
        private:
            IntervalSeries(const MatrixOneDim<Interval<T>> &coef) : m_coef(coef)
            {
            }	
            int fact(int n) const { return (n == 1 || n == 0) ? 1 : fact(n - 1) * n; }
            MatrixOneDim<Interval<T>> m_coef;
    };
    
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator+(const IntervalSeries &y) const
    {
        return IntervalSeries(m_coef + y.m_coef);
    } 
    
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator-(const IntervalSeries &y) const
    {
        return IntervalSeries(m_coef - y.m_coef);
    }
    
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator*(const IntervalSeries &y) const
    {
	const MatrixOneDim<Interval<T>>& l = m_coef;
        const MatrixOneDim<Interval<T>>& r = y.m_coef;
	MatrixOneDim<Interval<T>> h(l.size(), 0.0);

        h[0] = l.item(0) * r.item(0);
        for(int i=1; i < l.size(); i++)
	   h[i] = l.subMatrix(0, i) * r.subMatrix(0, i).reverse();
        return IntervalSeries(h);
    }
    
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator/(const IntervalSeries &y) const
    {
	const MatrixOneDim<Interval<T>>& l = m_coef;
        const MatrixOneDim<Interval<T>>& r = y.m_coef;
	MatrixOneDim<Interval<T>> h(l.size(), 0.0);

        if ((r.item(0) == 0.0) == IntervalBool::True)
            throw std::invalid_argument("IntervalSeries library. Invalid operation. Divide by 0.0.");
	
	h[0]=l.item(0)/r.item(0);
	for(int i=1; i < l.size(); i++)
	   h[i]= (1.0/r.item(0)) * (l.item(i) - h.subMatrix(0, i-1) * r.subMatrix(1, i).reverse());

        return IntervalSeries(h);
    }
    
    // h(x) = u(x)^r
    // h'(x) u(x) = r u'(x) h(x)
    // z(x) = h'(x) u(x) = r u'(x) h(x)
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator^(int exp) const
    {
        
        const MatrixOneDim<Interval<T>>& u = m_coef;
	MatrixOneDim<Interval<T>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T>>  z(u.size(), 0.0);

        h[0] = u.item(0) ^ exp;
	z[0] = (double)exp * u.item(1) * h.item(0);

	for(int i=1; i < u.size(); i++)
	{
	   z[i] =  ((double)exp/i) * ( MatrixOneDim<Interval<T>>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;        
	   h[i] =  i==1 ? z.item(i)/u.item(0) : z.item(i)/u.item(0) - 1.0/((double)i*u.item(0)) * ( MatrixOneDim<Interval<T>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse() ) ;
	}	
        return IntervalSeries<T>(h);
    }

    // h(x) = u(x)^r
    // h'(x) u(x) = r u'(x) h(x)
    // z(x) = h'(x) u(x) = r u'(x) h(x)
    template<class T> IntervalSeries<T> IntervalSeries<T>::operator^(const Interval<T> &exp) const
    {
        const MatrixOneDim<Interval<T>>& u = m_coef;
	MatrixOneDim<Interval<T>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T>>  z(u.size(), 0.0);

        h[0] = u.item(0) ^ exp;
	z[0] = exp * u.item(1) * h.item(0);

	for(int i=1; i < u.size(); i++)
	{
	   z[i] =  (exp/i) * ( MatrixOneDim<Interval<T>>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	   h[i] =  i==1 ? z[i]/u.item(0) : z[i]/u.item(0) - 1.0/((double)i*u.item(0)) * ( MatrixOneDim<Interval<T>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse() ) ;
	}	
        return IntervalSeries<T>(h);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator+(T3 t, const IntervalSeries<T2>&y)
    {
        const MatrixOneDim<Interval<T2>>& u = y.m_coef;
	MatrixOneDim<Interval<T2>> h(u.size(), 0.0);
	h[0]=t;
        return IntervalSeries<T2>(h+u);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator-(T3 t, const IntervalSeries<T2>&y)
    {
        const MatrixOneDim<Interval<T2>>& u = y.m_coef;
	MatrixOneDim<Interval<T2>> h(u.size(), 0.0);
	h[0]=t;
        return IntervalSeries<T2>(h-u);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator*(T3 t, const IntervalSeries<T2>&y)
    {
        return IntervalSeries<T2>(t * y.m_coef);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator/(T3 t, const IntervalSeries<T2>&y)
    {
	const MatrixOneDim<Interval<T2>>& u = y.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
        h[0] = t/u.item(0);

        for(int i=1; i < u.size(); i++)
	   h[i] = -(1.0/u.item(0)) * h.subMatrix(0, i-1) * u.subMatrix(1, i).reverse();
        return IntervalSeries<T2>(h);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator+(const IntervalSeries<T2>&y, T3 t)
    {
        const MatrixOneDim<Interval<T2>>& u = y.m_coef;
	MatrixOneDim<Interval<T2>> h(u.size(), 0.0);
	h[0]=t;
        return IntervalSeries<T2>(u+h);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator-(const IntervalSeries<T2>&y, T3 t)
    {
        const MatrixOneDim<Interval<T2>>& u = y.m_coef;
	MatrixOneDim<Interval<T2>> h(u.size(), 0.0);
	h[0]=t;
        return IntervalSeries<T2>(u-h);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator*(const IntervalSeries<T2>&y, T3 t)
    {
        return IntervalSeries<T2>(t * y.m_coef);
    }
    
    template<class T2, class T3> IntervalSeries<T2>operator/(const IntervalSeries<T2>&y, T3 t)
    {
        return IntervalSeries<T2>((1.0/t) * y.m_coef);
    }
    
    template<class T2> IntervalSeries<T2> sqr(const IntervalSeries<T2>&x)
    {
	const MatrixOneDim<Interval<T2>>& l = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(l.size(), 0.0);

        for(int i=0; i < l.size(); i++)
	   h[i] = l.subMatrix(0, i) * l.subMatrix(0, i).reverse();
        return IntervalSeries<T2>(h);
    }
    
    template<class T2> IntervalSeries<T2> sqrt(const IntervalSeries<T2>& x)
    {
	const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	if ((u.item(0) < 0.0) == IntervalBool::True)
            throw std::invalid_argument("The function Series<T>::sqrt is not defined for negative numbers"); 
        return x^0.5;
    }

    // s(x) = sin(u(x))
    // c(x) = cos(u(x))
    // s'(x) = u'(x) c(x)    
    // c'(x) = −u'(x) s(x),
    // s_k = (1/k) [1u_1, 2u_2, . . . , (k − 1)*u_k−1, k*u_k] [c_k−1, c_k−2, . . . , c_1, c_0]
    // c_k = (-1/k) [1u_1, 2u_2, . . . , (k − 1)*u_k−1, k*u_k] [s_k−1, s_k−2, . . . , s_1, s_0]
    template<class T2> IntervalSeries<T2> sin(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  sn(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  cs(u.size(), 0.0);

        sn[0] = interval::sin(u.item(0));
	cs[0] = interval::cos(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   cs[i] =  (-1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * sn.subMatrix(0, i-1).reverse() ) ;
	   sn[i] =  (1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * cs.subMatrix(0, i-1).reverse() ) ;
	}	
        return IntervalSeries<T2>(sn);
    }
    
    template<class T2> IntervalSeries<T2> cos(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  sn(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  cs(u.size(), 0.0);

        sn[0] = interval::sin(u.item(0));
	cs[0] = interval::cos(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   sn[i] =  (1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * cs.subMatrix(0, i-1).reverse() );
	   cs[i] =  (-1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * sn.subMatrix(0, i-1).reverse() );	   
	}	
        return IntervalSeries<T2>(cs);
    }

    // h(x) = tg(u(x))
    // v(x) = 1+ h(x)*h(x)
    // h'(x) = u'(x) v(x)
    // v'(x) = 2h'(x) h(x)  
    template<class T2> IntervalSeries<T2> tg(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

        h[0] = interval::tg(u.item(0));
	v[0] = 1.0 + interval::tg(u.item(0))*interval::tg(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   h[i] =  (1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * v.subMatrix(0, i-1).reverse() ) ;
	   v[i] =  (2.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(h.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	}	
        return IntervalSeries<T2>(h);
    }
    // h(x) = ctg(u(x))
    // v(x) = 1 + h(x)*h(x)
    // h'(x) = -u'(x) v(x)
    // v'(x) = 2h'(x) h(x) 
    template<class T2> IntervalSeries<T2> ctg(const IntervalSeries<T2>&x)
    {
	const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

        h[0] = interval::ctg(u.item(0));
	v[0] = 1.0 + interval::ctg(u.item(0)) * interval::ctg(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   h[i] =  (-1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * v.subMatrix(0, i-1).reverse() ) ;
	   v[i] =  (2.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(h.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	}	
        return IntervalSeries<T2>(h);
    }
    
    // h(x) = arcsin(u(x))
    // v(x) = sqrt(1 − sqr(u(x)))
    // u'(x) = h'(x) v(x)
    // v'(x) = −h'(x) u(x)
    template<class T2> IntervalSeries<T2> asin(const IntervalSeries<T2>&x)
    {       
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
        if ((u.item(0) < -1.0 ) == IntervalBool::True || (u.item(0) > 1.0 ) == IntervalBool::True)
            throw std::invalid_argument("Invalid argument in IntervalSeries::asin. The argument is out of this interval [-1,1]."); 

	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

        h[0] = interval::asin(u.item(0));
	v[0] = interval::sqrt(1.0 - u.item(0)*u.item(0));
        h[1] = u.item(1)/v.item(0);
        v[1] = -h.item(1) * u.item(0);

	for(int i=2; i < u.size(); i++)
	{	   
	   h[i] =  (1.0/v.item(0)) * (u.item(i) - (1.0/i) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;
           v[i] =  (-1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(h.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	}	
        return IntervalSeries<T2>(h);
    }
    // h(x) = arccos(u(x))
    // v(x) = sqrt(1 − sqr(u(x)))
    // u'(x) = -h'(x) v(x)
    // v'(x) = h'(x) u(x)
    template<class T2> IntervalSeries<T2> acos(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
        if ((u.item(0) < -1.0 ) == IntervalBool::True || (u.item(0) > 1.0 ) == IntervalBool::True)
            throw std::invalid_argument("Invalid argument in IntervalSeries::acos. The argument is out of this interval [-1,1]."); 
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

        h[0] = interval::acos(u.item(0));
	v[0] = interval::sqrt(1.0 - u.item(0)*u.item(0));
        h[1] = -u.item(1)/v.item(0);
        v[1] = h.item(1)*u.item(0);

	for(int i=2; i < u.size(); i++)
	{	   
	   h[i] =  (-1.0/v.item(0)) * (u.item(i) + (1.0/i) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;
           v[i] =  (1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1.0, i).mulItems(h.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	}	
        return IntervalSeries<T2>(h);
    }
    
    // h(x) = arctan(u(x)) 
    // v(x) = 1+sqr(u(x)) 
    // u'(x) = h'(x) v(x) 
    // v'(x) = 2u'(x) u(x)
    template<class T2> IntervalSeries<T2> atg(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

	v[0] = 1.0 + u.item(0)*u.item(0); 
        h[0] = interval::atg(u.item(0));
	v[1] = 2.0 * u.item(1)*u.item(0);      
        h[1] = u.item(1)/v.item(0);
        	

	for(int i=2; i < u.size(); i++)
	{
	   v[i] =  2.0 * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	   h[i] =  (1.0/v.item(0)) * (u.item(i) - (1.0/i) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;        	   
	}	
        return IntervalSeries<T2>(h);            
    }
    
    // h(x) = arcctan(u(x)) 
    // v(x) = 1+sqr(u(x)) 
    // u'(x) = -h'(x) v(x) 
    // v'(x) = 2u'(x) u(x)
    template<class T2> IntervalSeries<T2> actg(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);
	MatrixOneDim<Interval<T2>>  v(u.size(), 0.0);

	v[0] = 1.0 + u.item(0)*u.item(0); 
        h[0] = M_PI_2 - interval::atg(u.item(0));
	v[1] = 2.0 * u.item(1)*u.item(0);      
        h[1] = -u.item(1)/v.item(0);
        	

	for(int i=2; i < u.size(); i++)
	{
	   v[i] =  2.0 * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	   h[i] =  (-1.0/v.item(0)) * (u.item(i) + (1.0/i) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;        	   
	}
        return IntervalSeries<T2>(h); 
    }
    // h(x) = ln(u(x))
    // u'(x) = h'(x) u(x). 
    template<class T2> IntervalSeries<T2> ln(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);

        if ((u.item(0) <= 0.0) == IntervalBool::True )
            throw std::invalid_argument("The function IntervalSeries<T>::ln is not define for negative numbers and 0.0");

        h[0] = interval::ln(u.item(0));
	h[1] = u.item(1)/u.item(0);

	for(int i=2; i < u.size(); i++)
	   h[i]= u.item(i)/u.item(0) - 1.0/((double)i*u.item(0)) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse();
	
        return IntervalSeries<T2>(h);
    }
    
    // h(x) = log_b(u(x)) = ln(u(x))/ln(b)
    // u'(x) = ln(b) * h'(x) * u(x). 
    template<class T2> IntervalSeries<T2> log(const IntervalSeries<T2> &x, double base)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);

        if ((u.item(0) <= 0.0) == IntervalBool::True )
            throw std::invalid_argument("The function IntervalSeries<T>::log is not define for negative numbers and 0.0");

        h[0] = interval::log(u.item(0), base);
	h[1] = u.item(1)/(u.item(0) * std::log(base));

	for(int i=2; i < u.size(); i++)
	   h[i]= u.item(i)/(std::log(base)*u.item(0)) - 1.0/((double)i*u.item(0)) * MatrixOneDim<Interval<T2>>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse();
	
        return IntervalSeries<T2>(h);
    }
    
    // h(x)=exp(u(x))
    // h'(x)=u'(x)h(x)
    // h_k = (1/k)[1*u_1,..,k*u_k][h_(k-1),..,h_0]
    template<class T2> IntervalSeries<T2> exp(const IntervalSeries<T2>&x)
    {
        const MatrixOneDim<Interval<T2>>& u = x.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);

        h[0] = interval::exp(u.item(0));

	for(int i=1; i < u.size(); i++)
	   h[i] =  (1.0/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	
        return IntervalSeries<T2>(h);
    }
    
    template<class T2, class T3> IntervalSeries<T2> pow(const T3 &base, const IntervalSeries<T2> &exp)
    {
        const MatrixOneDim<Interval<T2>>& u = exp.m_coef;
	MatrixOneDim<Interval<T2>>  h(u.size(), 0.0);

        h[0] = base ^ u.item(0);

	for(int i=1; i < u.size(); i++)
	   h[i] =  (interval::ln(base)/i) * ( MatrixOneDim<Interval<T2>>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	
        return IntervalSeries<T2>(h);
    }
    
    template<class T2> IntervalSeries<T2> abs(const IntervalSeries<T2>&x)
    {
        if((x.m_coef.item(0)==0.0) == IntervalBool::True )
            throw std::invalid_argument("Invalid argument in IntervalSeries::get_abs. There isn't derivation at zero.");
        
        IntervalBool rez = x.m_coef.item(0) < 0.0;
        if(rez==IntervalBool::True)
            return IntervalSeries<T2>(-1.0 * x.m_coef);
        else if(rez==IntervalBool::False)
            return IntervalSeries<T2>(x.m_coef);
        else //IntervalBool::Intermediate case
            throw std::invalid_argument("Invalid argument in IntervalSeries::get_abs. There isn't derivation at zero.");
    }
    
    template<class T2> IntervalSeries<T2> min(const IntervalSeries<T2> &x, const IntervalSeries<T2> &y)
    {
        auto condition = (x.m_coef.item(0) <= y.m_coef.item(0));
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
        throw std::invalid_argument("IntervalSeries::min operation is not defined.");
    }
    
    template<class T2> IntervalSeries<T2> max(const IntervalSeries<T2> &x, const IntervalSeries<T2> &y)
    {
        auto condition = (x.m_coef.item(0) > y.m_coef.item(0));
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
        throw std::invalid_argument("IntervalSeries::max operation is not defined.");
    }
            
    template<class T2> IntervalSeries<T2> ifThen(IntervalBool condition, const IntervalSeries<T2> &x, const IntervalSeries<T2> &y)
    {
        if(condition==IntervalBool::True) 
            return x;
        else if(condition==IntervalBool::False) 
            return y;
        else
            throw std::invalid_argument("IntervalSeries::ifThen operation is not defined.");
    }
    
    template<class T2> IntervalBool cond(Conditions condition, const IntervalSeries<T2> &x, const IntervalSeries<T2> &y)
    {
        auto left = x.m_coef.item(0);
        auto right = y.m_coef.item(0);
        
        switch (condition)
        {
            case Conditions::More:
                return left > right;
                break;
            case Conditions::Less:
                return left < right;
                break;
            case Conditions::LessEqual:
                return left <= right;
                break;
            case Conditions::MoreEqual:
                return left >= right;
                break;
        }
        throw std::invalid_argument("Invalid condition in IntervalSeries::Condition.");
    }
        
    template<class T2> std::ostream& operator<<(std::ostream & out, const IntervalSeries<T2> &x)
    {
        std::cout << " val: " << x.m_coef.item(0);
        std::cout << " der: ";
        for(int i=1; i < x.m_coef.size(); i++)
	   std::cout << x.der(i) << ' '; 
        std::cout << '\n';
        return out;
    }
   
  }
}

#endif /* INTERVAL_SERIES_HPP */

