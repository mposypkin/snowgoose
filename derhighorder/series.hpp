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

#ifndef SERIES_HPP
#define SERIES_HPP

#include <iostream> 
#include <cmath>
#include <vector>
#include "matrixod.hpp"
#include "interval/enums.h"


namespace snowgoose {
  namespace derhighorder {
    template<class T> class Series
    {
        public:
            Series(T val, T der, int order) : m_coef(order + 1, 0.0)
            {
		m_coef[0] = val;
                m_coef[1] = der;
	    }
            Series(){}
            using value_type = T;
            Series operator+(const Series &y) const;
            Series operator-(const Series &y) const;
            Series operator*(const Series &y) const;
            Series operator/(const Series &y) const;            
            Series operator^(int exp) const;
            Series operator^(T exp) const;
            
            template<class T2> friend Series<T2> operator+(T2 t, const Series<T2> &y);
            template<class T2> friend Series<T2> operator-(T2 t, const Series<T2> &y);
            template<class T2> friend Series<T2> operator*(T2 t, const Series<T2> &y);
            template<class T2> friend Series<T2> operator/(T2 t, const Series<T2> &y);           
            template<class T2> friend Series<T2> operator+(const Series<T2> &y, T2 t);
            template<class T2> friend Series<T2> operator-(const Series<T2> &y, T2 t);
            template<class T2> friend Series<T2> operator*(const Series<T2> &y, T2 t);
            template<class T2> friend Series<T2> operator/(const Series<T2> &y, T2 t);
                       
            template<class T2> friend Series<T2> sqr(const Series<T2> &x);
            template<class T2> friend Series<T2> sqrt(const Series<T2> &x);            
            template<class T2> friend Series<T2> sin(const Series<T2> &x);
            template<class T2> friend Series<T2> cos(const Series<T2>&x);            
            template<class T2> friend Series<T2> tg(const Series<T2> &x);  
            template<class T2> friend Series<T2> ctg(const Series<T2> &x); 
            template<class T2> friend Series<T2> asin(const Series<T2> &x);
            template<class T2> friend Series<T2> acos(const Series<T2> &x); 
            template<class T2> friend Series<T2> atg(const Series<T2> &x);
            template<class T2> friend Series<T2> actg(const Series<T2> &x); 
            template<class T2> friend Series<T2> ln(const Series<T2> &x);
            template<class T2> friend Series<T2> log(const Series<T2> &x, double base);
            template<class T2> friend Series<T2> exp(const Series<T2> &x);
            template<class T2> friend Series<T2> pow(const T2 &base, const Series<T2> &exp);
            template<class T2> friend Series<T2> abs(const Series<T2> &x);
            template<class T2> friend Series<T2> min(const Series<T2> &x, const Series<T2> &y);
            template<class T2> friend Series<T2> max(const Series<T2> &x, const Series<T2> &y);
            template<class T2> friend Series<T2> ifThen(IntervalBool condition, const Series<T2> &x, const Series<T2> &y);
            template<class T2> friend IntervalBool cond(Conditions condition, const Series<T2> &x, const Series<T2> &y);
            
            template<class T2> friend std::ostream& operator<<(std::ostream & out, const Series<T2> &x);
            T value() const { return m_coef.item(0); }
            T der(int order) const { return m_coef.item(order) * fact(order); }
            explicit operator T(){ return m_coef.item(0); }
        private:
            Series(const MatrixOneDim<T> &coef) : m_coef(coef)
            {
            }	
            int fact(int n) const { return (n == 1 || n == 0) ? 1 : fact(n - 1) * n; }
            MatrixOneDim<T> m_coef;
    };
    
    template<class T> Series<T> Series<T>::operator+(const Series &y) const
    {
        return Series(m_coef + y.m_coef);
    } 
    
    template<class T> Series<T> Series<T>::operator-(const Series &y) const
    {
        return Series(m_coef - y.m_coef);
    }
    
    template<class T> Series<T> Series<T>::operator*(const Series &y) const
    {
	const MatrixOneDim<T>& l = m_coef;
        const MatrixOneDim<T>& r = y.m_coef;
	MatrixOneDim<T> h(l.size(), 0.0);

        h[0] = l.item(0) * r.item(0);
        for(int i=1; i < l.size(); i++)
	   h[i] = l.subMatrix(0, i) * r.subMatrix(0, i).reverse();
        return Series(h);
    }
    
    template<class T> Series<T> Series<T>::operator/(const Series &y) const
    {
	const MatrixOneDim<T>& l = m_coef;
        const MatrixOneDim<T>& r = y.m_coef;
	MatrixOneDim<T> h(l.size(), 0.0);

        if (r.item(0) == 0.0)
            throw std::invalid_argument("Invalid operation. Divide by 0.0.");
	
	h[0]=l.item(0)/r.item(0);
	for(int i=1; i < l.size(); i++)
	   h[i]= (1.0/r.item(0)) * (l.item(i) - h.subMatrix(0, i-1) * r.subMatrix(1, i).reverse());

        return Series(h);
    }
    
    template<class T> Series<T> Series<T>::operator^(int exp) const
    {

	const MatrixOneDim<T>& l = m_coef;
	MatrixOneDim<T>  h = l;
	for(int j=1; j < exp; j++)	
	{
	   MatrixOneDim<T> t = h; 
	   for(int i=0; i < l.size(); i++)   	
              h[i] = t.subMatrix(0, i) * l.subMatrix(0, i).reverse();
	}
        return Series<T>(h);
    }

    // h(x) = u(x)^r
    // h'(x) u(x) = r u'(x) h(x)
    // z(x) = h'(x) u(x) = r u'(x) h(x)
    template<class T> Series<T> Series<T>::operator^(T exp) const
    {
        const MatrixOneDim<T>& u = m_coef;
	MatrixOneDim<T>  h(u.size(), 0.0);
	MatrixOneDim<T>  z(u.size(), 0.0);

        if (u.item(0) == 0.0)
            return Series<T>(h);

        h[0] = std::pow(u.item(0), exp);
	z[0] = exp * u.item(1) * h.item(0);

	for(int i=1; i < u.size(); i++)
	{
	   z[i] =  (exp/i) * ( MatrixOneDim<T>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	   h[i] =  i==1 ? z[i]/u.item(0) : z[i]/u.item(0) - 1.0/(i*u.item(0)) * ( MatrixOneDim<T>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse() ) ;
	}	
        return Series<T>(h);
    }
    
    template<class T2> Series<T2>operator+(T2 t, const Series<T2>&y)
    {
        const MatrixOneDim<T2>& u = y.m_coef;
	MatrixOneDim<T2> h(u.size(), 0.0);
	h[0]=t;
        return Series<T2>(h+u);
    }
    
    template<class T2> Series<T2>operator-(T2 t, const Series<T2>&y)
    {
        const MatrixOneDim<T2>& u = y.m_coef;
	MatrixOneDim<T2> h(u.size(), 0.0);
	h[0]=t;
        return Series<T2>(h-u);
    }
    
    template<class T2> Series<T2>operator*(T2 t, const Series<T2>&y)
    {
        return Series<T2>(t * y.m_coef);
    }
    
    template<class T2> Series<T2>operator/(T2 t, const Series<T2>&y)
    {
	const MatrixOneDim<T2>& u = y.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);

        if (u.item(0) == 0.0)
            throw std::invalid_argument("Invalid operation. Divide by 0.0.");

        h[0] = t/u.item(0);

        for(int i=1; i < u.size(); i++)
	   h[i] = -(1.0/u.item(0)) * h.subMatrix(0, i-1) * u.subMatrix(1, i).reverse();
        return Series<T2>(h);
    }
    
    template<class T2> Series<T2>operator+(const Series<T2>&y, T2 t)
    {
        const MatrixOneDim<T2>& u = y.m_coef;
	MatrixOneDim<T2> h(u.size(), 0.0);
	h[0]=t;
        return Series<T2>(u+h);
    }
    
    template<class T2> Series<T2>operator-(const Series<T2>&y, T2 t)
    {
        const MatrixOneDim<T2>& u = y.m_coef;
	MatrixOneDim<T2> h(u.size(), 0.0);
	h[0]=t;
        return Series<T2>(u-h);
    }
    
    template<class T2> Series<T2>operator*(const Series<T2>&y, T2 t)
    {
        return Series<T2>(t * y.m_coef);
    }
    
    template<class T2> Series<T2>operator/(const Series<T2>&y, T2 t)
    {
        if (t == 0.0)
            throw std::invalid_argument("Invalid operation. Divide by 0.0.");
        return Series<T2>((1.0/t) * y.m_coef);
    }
    
    template<class T2> Series<T2> sqr(const Series<T2>&x)
    {
	const MatrixOneDim<T2>& l = x.m_coef;
	MatrixOneDim<T2>  h(l.size(), 0.0);

        for(int i=0; i < l.size(); i++)
	   h[i] = l.subMatrix(0, i) * l.subMatrix(0, i).reverse();
        return Series<T2>(h);
    }
    
    template<class T2> Series<T2> sqrt(const Series<T2>& x)
    { 
        const MatrixOneDim<T2>& u = x.m_coef;

        if (u.item(0) <= 0.0)
            throw std::invalid_argument("Invalid argument in Series<T>::sqrt. There isn't derivation for negative numbers and  zero.");
	
        return x^0.5;
    }

    // s(x) = sin(u(x))
    // c(x) = cos(u(x))
    // s'(x) = u'(x) c(x)    
    // c'(x) = −u'(x) s(x),
    // s_k = (1/k) [1u_1, 2u_2, . . . , (k − 1)*u_k−1, k*u_k] [c_k−1, c_k−2, . . . , c_1, c_0]
    // c_k = (-1/k) [1u_1, 2u_2, . . . , (k − 1)*u_k−1, k*u_k] [s_k−1, s_k−2, . . . , s_1, s_0]
    template<class T2> Series<T2> sin(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  sn(u.size(), 0.0);
	MatrixOneDim<T2>  cs(u.size(), 0.0);

        sn[0] = std::sin(u.item(0));
	cs[0] = std::cos(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   cs[i] =  (-1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * sn.subMatrix(0, i-1).reverse() ) ;
	   sn[i] =  (1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * cs.subMatrix(0, i-1).reverse() ) ;
	}	
        return Series<T2>(sn);
    }
    
    template<class T2> Series<T2> cos(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  sn(u.size(), 0.0);
	MatrixOneDim<T2>  cs(u.size(), 0.0);

        sn[0] = std::sin(u.item(0));
	cs[0] = std::cos(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   sn[i] =  (1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * cs.subMatrix(0, i-1).reverse() );
	   cs[i] =  (-1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * sn.subMatrix(0, i-1).reverse() );	   
	}	
        return Series<T2>(cs);
    }

    // h(x) = tg(u(x))
    // v(x) = 1+ h(x)*h(x)
    // h'(x) = u'(x) v(x)
    // v'(x) = 2h'(x) h(x)  
    template<class T2> Series<T2> tg(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

        h[0] = std::tan(u.item(0));
	v[0] = 1.0 + std::tan(u.item(0))*std::tan(u.item(0));

	for(int i=1; i < u.size(); i++)
	{
	   h[i] =  (1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * v.subMatrix(0, i-1).reverse() ) ;
	   v[i] =  (2.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(h.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	}	
        return Series<T2>(h);
    }
    // h(x) = tg(u(x))
    // v(x) = 1 + h(x)*h(x)
    // h'(x) = -u'(x) v(x)
    // v'(x) = 2h'(x) h(x) 
    template<class T2> Series<T2> ctg(const Series<T2>&x)
    {
	const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

        h[0] = 1.0/std::tan(u.item(0));
	v[0] = 1.0 + 1.0/(std::tan(u.item(0)) * std::tan(u.item(0)));

	for(int i=1; i < u.size(); i++)
	{
	   h[i] =  (-1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * v.subMatrix(0, i-1).reverse() ) ;
	   v[i] =  (2.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(h.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	}	
        return Series<T2>(h);
    }
    
    // h(x) = arcsin(u(x))
    // v(x) = sqrt(1 − sqr(u(x)))
    // u'(x) = h'(x) v(x)
    // v'(x) = −h'(x) u(x)
    template<class T2> Series<T2> asin(const Series<T2>&x)
    {       
        const MatrixOneDim<T2>& u = x.m_coef;
        if(u.item(0) < -1.0 || u.item(0) > 1.0)
            throw std::invalid_argument("Invalid argument in Series::asin. The argument is out of this interval [-1,1]."); 
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

        h[0] = std::asin(u.item(0));
	v[0] = std::sqrt(1.0 - u.item(0)*u.item(0));
        
        if(v.item(0)==0.0)
		throw std::invalid_argument("Invalid argument in Series::asin. Derivation does not exist.");
        h[1] = u.item(1)/v.item(0);
        v[1] = -h.item(1) * u.item(0);

	for(int i=2; i < u.size(); i++)
	{	   
	   h[i] =  (1.0/v.item(0)) * (u.item(i) - (1.0/i) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;
           v[i] =  (-1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(h.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	}	
        return Series<T2>(h);
    }
    // h(x) = arccos(u(x))
    // v(x) = sqrt(1 − sqr(u(x)))
    // u'(x) = -h'(x) v(x)
    // v'(x) = h'(x) u(x)
    template<class T2> Series<T2> acos(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
        if(u.item(0) < -1.0 || u.item(0) > 1.0)
            throw std::invalid_argument("Invalid argument in Series::acos. The argument is out of this interval [-1,1]."); 
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

        h[0] = std::acos(u.item(0));
	v[0] = std::sqrt(1.0 - u.item(0)*u.item(0));
        if(v.item(0)==0.0)
		throw std::invalid_argument("Invalid argument in Series::acos. Derivation does not exist.");
        h[1] = -u.item(1)/v.item(0);
        v[1] = h.item(1)*u.item(0);

	for(int i=2; i < u.size(); i++)
	{	   
	   h[i] =  (-1.0/v.item(0)) * (u.item(i) + (1.0/i) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;
           v[i] =  (1.0/i) * ( MatrixOneDim<T2>::seq(1.0, i).mulItems(h.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	}	
        return Series<T2>(h);
    }
    
    // h(x) = arctan(u(x)) 
    // v(x) = 1+sqr(u(x)) 
    // u'(x) = h'(x) v(x) 
    // v'(x) = 2u'(x) u(x)
    template<class T2> Series<T2> atg(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

	v[0] = 1.0 + u.item(0)*u.item(0); 
        h[0] = std::atan(u.item(0));
	v[1] = 2.0 * u.item(1)*u.item(0);      
        h[1] = u.item(1)/v.item(0);
        	

	for(int i=2; i < u.size(); i++)
	{
	   v[i] =  2.0 * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	   h[i] =  (1.0/v.item(0)) * (u.item(i) - (1.0/i) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;        	   
	}	
        return Series<T2>(h);            
    }
    
    // h(x) = arcctan(u(x)) 
    // v(x) = 1+sqr(u(x)) 
    // u'(x) = -h'(x) v(x) 
    // v'(x) = 2u'(x) u(x)
    template<class T2> Series<T2> actg(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);
	MatrixOneDim<T2>  v(u.size(), 0.0);

	v[0] = 1.0 + u.item(0)*u.item(0); 
        h[0] = M_PI_2 - std::atan(u.item(0));
	v[1] = 2.0 * u.item(1)*u.item(0);      
        h[1] = -u.item(1)/v.item(0);
        	

	for(int i=2; i < u.size(); i++)
	{
	   v[i] =  2.0 * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * u.subMatrix(0, i-1).reverse() ) ;	   
	   h[i] =  (-1.0/v.item(0)) * (u.item(i) + (1.0/i) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * v.subMatrix(1, i-1).reverse() ) ;        	   
	}
        return Series<T2>(h); 
    }
    // h(x) = ln(u(x))
    // u'(x) = h'(x) u(x). 
    template<class T2> Series<T2> ln(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);

        if (u.item(0) <= 0.0)
            throw std::invalid_argument("The function Series<T>::ln is not define for negative numbers and 0.0");

        h[0] = std::log(u.item(0));
	h[1] = u.item(1)/u.item(0);

	for(int i=2; i < u.size(); i++)
	   h[i]= u.item(i)/u.item(0) - 1.0/(i*u.item(0)) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse();
	
        return Series<T2>(h);
    }
    
    // h(x) = log_b(u(x)) = ln(u(x))/ln(b)
    // u'(x) = ln(b) * h'(x) * u(x). 
    template<class T2> Series<T2> log(const Series<T2> &x, double base)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);

        if (u.item(0) <= 0.0)
            throw std::invalid_argument("The function Series<T>::log is not define for negative numbers and 0.0");

        h[0] = std::log(u.item(0))/std::log(base);
	h[1] = u.item(1)/(u.item(0) * std::log(base));

	for(int i=2; i < u.size(); i++)
	   h[i]= u.item(i)/(std::log(base)*u.item(0)) - 1.0/(i*u.item(0)) * MatrixOneDim<T2>::seq(1, i-1).mulItems(h.subMatrix(1, i-1)) * u.subMatrix(1, i-1).reverse();
	
        return Series<T2>(h);
    }
    
    // h(x)=exp(u(x))
    // h'(x)=u'(x)h(x)
    // h_k = (1/k)[1*u_1,..,k*u_k][h_(k-1),..,h_0]
    template<class T2> Series<T2> exp(const Series<T2>&x)
    {
        const MatrixOneDim<T2>& u = x.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);

        h[0] = std::exp(u.item(0));

	for(int i=1; i < u.size(); i++)
	   h[i] =  (1.0/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	
        return Series<T2>(h);
    }
    
    template<class T2> Series<T2> pow(const T2 &base, const Series<T2> &exp)
    {
        const MatrixOneDim<T2>& u = exp.m_coef;
	MatrixOneDim<T2>  h(u.size(), 0.0);

        h[0] = std::pow(base, u.item(0));

	for(int i=1; i < u.size(); i++)
	   h[i] =  (std::log(base)/i) * ( MatrixOneDim<T2>::seq(1, i).mulItems(u.subMatrix(1, i)) * h.subMatrix(0, i-1).reverse() ) ;
	
        return Series<T2>(h);
    }
    
    template<class T2> Series<T2> abs(const Series<T2>&x)
    {
        if(x.m_coef.item(0)==0.0)
            throw std::invalid_argument("Invalid argument in Series::abs. There isn't derivation at zero.");
        else if(x.m_coef.item(0) < 0.0)
            return Series<T2>(-1.0 * x.m_coef);
        else
            return Series<T2>(x.m_coef);
    }
    
    template<class T2> Series<T2> min(const Series<T2> &x, const Series<T2> &y)
    {
        return (x.m_coef.item(0) < y.m_coef.item(0)) ? x : y;
    }
    
    template<class T2> Series<T2> max(const Series<T2> &x, const Series<T2> &y)
    {
        return (x.m_coef.item(0) > y.m_coef.item(0)) ? x : y;
    }
            
    template<class T2> Series<T2> ifThen(IntervalBool condition, const Series<T2> &x, const Series<T2> &y)
    {
        return (condition == IntervalBool::True)? x : y;
    }
    
    template<class T2> IntervalBool cond(Conditions condition, const Series<T2> &x, const Series<T2> &y)
    {
        auto left = x.m_coef.item(0);
        auto right = y.m_coef.item(0);
        switch (condition)
        {
            case Conditions::More :
                return left > right ? IntervalBool::True : IntervalBool::False;
                break;
            case Conditions::Less :
                return left < right ? IntervalBool::True : IntervalBool::False;
                break;
            case Conditions::LessEqual :
                return left <= right ? IntervalBool::True : IntervalBool::False;
                break;
            case Conditions::MoreEqual :
                return left >= right ? IntervalBool::True : IntervalBool::False;
                break;
        }
        throw std::invalid_argument("Invalid condition in Series::Condition.");
    }
        
    template<class T2> std::ostream& operator<<(std::ostream & out, const Series<T2> &x)
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

#endif /* SERIES_HPP */

