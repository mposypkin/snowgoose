/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   interval_.hpp
 * Author: Alexander Usov
 *
 * Created on January 19, 2017, 10:51 AM
 */

#ifndef INTERVAL__HPP
#define INTERVAL__HPP

#include <iostream> 
#include <limits>
#include <math.h>
#include <exception>
#include <unordered_map>
#include <vector>
#include <set>
#include <memory>
#include <functional>
#include "common/utilmacro.hpp"
#include "common/sgerrcheck.hpp"
#include "calcvalues.hpp"

namespace snowgoose { namespace interval {

    template <class T> class SubInterval;
    template <class T> class Value;

    template<class T> class Interval
    {
        public:
            Interval(T lb, T rb)
            {
                SubInterval<T> si(lb, rb);
                m_SubIntervals.insert(si);
            }

            Interval(T t)
            {
                SubInterval<T> si(t, t);
                m_SubIntervals.insert(si);
            }

            Interval operator+(const Interval &y) const
            {
                return DoOperation(y, [](const SubInterval<T>& lsi, const SubInterval<T>& rsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi + rsi); });
            }

            Interval operator-(const Interval &y) const
            {
                return DoOperation(y, [](const SubInterval<T>& lsi, const SubInterval<T>& rsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi - rsi); });
            }

            Interval operator*(const Interval &y) const
            {
                return DoOperation(y, [](const SubInterval<T>& lsi, const SubInterval<T>& rsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi * rsi); });
            }

            Interval operator/(const Interval &y) const
            {
                return DoOperation(y, [](const SubInterval<T>& lsi, const SubInterval<T>& rsi, std::multiset<SubInterval<T>>& ms) { auto v = lsi/rsi;  ms.insert(v.begin(), v.end()); });
            }
            
            Interval operator^(int n) const
            {
                return DoOperation([&](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi ^ n); });
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
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.sin()); });
            }
            friend Interval cos(const Interval<T> &x)
            {
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.cos()); });
            }
            friend Interval exp(const Interval<T> &x)
            {
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.exp()); });
            }
            friend Interval sqrt(const Interval<T> &x)
            {
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.sqrt()); });
            }
            friend Interval abs(const Interval<T> &x)
            {
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.abs()); });
            }
            friend Interval ln(const Interval<T> &x)
            {
                return x.DoOperation([](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.ln()); });
            }
            friend Interval log(const Interval<T> &x, T base)
            {
                return x.DoOperation([&](const SubInterval<T>& lsi, std::multiset<SubInterval<T>>& ms) { ms.insert(lsi.log(base)); });
            }
            friend std::ostream& operator<<(std::ostream &out, const Interval<T> &x)
            {
                out << "Interval begins\n";
                x.DoTraverse([&](const SubInterval<T>& si) {out << si; });
                return out << "Interval ends\n";
            }
            
            T lb()const {return m_SubIntervals.begin()->lb();}
            T rb()const {return m_SubIntervals.begin()->rb();}

        private:

            Interval(const std::multiset<SubInterval<T>>& subIntervals) : m_SubIntervals(subIntervals)
            {
            }

            Interval DoOperation(const Interval& y, std::function<void(const SubInterval<T>&, const SubInterval<T>&, std::multiset<SubInterval<T>>&)> f) const
            {
                std::multiset<SubInterval<T>> pvNew;
                for (auto lsi : m_SubIntervals)
                    for (auto rsi : y.m_SubIntervals)
                        f(lsi, rsi, pvNew);
                Interval z(pvNew);
                z.Union();
                return z;
            }

            Interval DoOperation(std::function<void(const SubInterval<T>&, std::multiset<SubInterval<T>>&)> f) const
            {
                std::multiset<SubInterval<T>> pvNew;
                for (auto lsi : m_SubIntervals)
                    f(lsi, pvNew);
                Interval z(pvNew);
                z.Union();
                return z;
            }

            void DoTraverse(std::function<void(const SubInterval<T>&)> f) const
            {
                for (auto si : m_SubIntervals)
                    f(si);
            }

            void Union()
            {
                std::multiset<SubInterval<T>> msNew;
                if(!m_SubIntervals.empty())
                    msNew.insert(*m_SubIntervals.begin());

                for (auto it = ++m_SubIntervals.begin(); it != m_SubIntervals.end(); it++)
                {
                    auto itNew = --msNew.end();
                    if (itNew->RightBound() >= it->LeftBound()) {
                        SubInterval<T> si = *itNew;
                        si.Union(*it);
                        msNew.erase(itNew);
                        msNew.insert(si);
                    }
                    else
                        msNew.insert(*it);
                }
                if (m_SubIntervals.size() > msNew.size()) {
                    m_SubIntervals.clear();
                    m_SubIntervals.insert(msNew.begin(), msNew.end());
                }					
            }

            std::multiset<SubInterval<T>> m_SubIntervals;
    };    
    

    template <class T> class SubInterval
    {
        protected:
            Value<T> m_lb;
            Value<T> m_rb;
        public:
            SubInterval(const Value<T> &lb, const Value<T> &rb) : m_lb(lb), m_rb(rb) 
            {
                if (lb < 0.0 && rb == 0.0)
                    m_rb.m_isMinusZero = true;
            }
            const Value<T>& LeftBound() const {return m_lb;}
            const Value<T>& RightBound() const {return m_rb;}
            SubInterval operator+(const SubInterval &y) const {
                auto lb = m_lb + y.m_lb;
                auto rb = m_rb + y.m_rb;
                return SubInterval(lb, rb);
            }

            SubInterval operator-(const SubInterval &y) const {
                auto lb = m_lb - y.m_rb;
                auto rb = m_rb - y.m_lb;
                return SubInterval(lb, rb);
            }
            SubInterval operator*(const SubInterval &y) const {
                Value<T> t1 = m_lb * y.m_lb;
                Value<T> lb = t1;
                Value<T> rb = t1;

                Value<T> t2 = m_lb * y.m_rb;
                lb = SGMIN(lb, t2);
                rb = SGMAX(rb, t2);

                Value<T> t3 = m_rb * y.m_lb;
                lb = SGMIN(lb, t3);
                rb = SGMAX(rb, t3);

                Value<T> t4 = m_rb * y.m_rb;
                lb = SGMIN(lb, t4);
                rb = SGMAX(rb, t4); 

                return SubInterval<T>(lb, rb);
            }

            std::vector<SubInterval> operator/(const SubInterval &y) const {
                if ((y.m_lb <= 0.0 && y.m_rb >= 0.0) && (m_lb <= 0.0 && m_rb >= 0.0))
                    throw std::invalid_argument("Invalid operation. Interval that includes 0.0 is divided by interval that includes 0.0.");
                
                std::vector<SubInterval> vsi;
                if(y.m_lb < 0.0 && y.m_rb > 0.0)
                {
                    auto si1 = div(SubInterval(y.m_lb, 0.0));
                    auto si2 = div(SubInterval(0.0, y.m_rb));
                    vsi.push_back(si1);
                    vsi.push_back(si2);
                    return vsi;
                }
                else
                {
                    auto si = div(y);
                    vsi.push_back(si);
                    return vsi;
                }
            }

            SubInterval operator^(int n) const
            {
                Value<T> t1 = 1.0, t2 = 1.0, lb, rb;
                for (int i = 1; i <= n; i++) {
                    t1 = t1 * m_lb;
                    t2 = t2 * m_rb;
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
                return SubInterval(lb, rb);
            }

            SubInterval sin() const
            {
                if (m_lb.GetType() == TypeValue::MinusInf || m_lb.GetType() == TypeValue::PlusInf || m_rb.GetType() == TypeValue::MinusInf || m_rb.GetType() == TypeValue::PlusInf)
                    return SubInterval(-1, 1);

                int k, l;
                bool p2 = false;
                bool p32 = false;
                Value<T> lb, rb;

                k = (int)::ceil(m_lb.GetValue() * M_2_PI);
                l = (int)::floor(m_rb.GetValue() * M_2_PI);

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
                    rb = SGMAX(::sin(m_lb.GetValue()), ::sin(m_rb.GetValue()));
                if (p32)
                    lb = -1.0;
                else
                    lb = SGMIN(::sin(m_lb.GetValue()), ::sin(m_rb.GetValue()));

                return SubInterval(lb, rb);
            }

            SubInterval cos() const
            {
                Value<T> t1 = m_lb.GetValue() + M_PI_2;
                Value<T> t2 = m_rb.GetValue() + M_PI_2;
                SubInterval t(t1, t2);
                return t.sin();
            }

            SubInterval exp() const
            {
                Value<T> t1 = m_lb.exp();
                Value<T> t2 = m_rb.exp();
                Value<T> lb = SGMIN(t1, t2);
                Value<T> rb = SGMAX(t1, t2);
                return SubInterval(lb, rb);
            }

            SubInterval sqrt() const
            {
                Value<T> a = m_lb.sqrt();
                Value<T> b = m_rb.sqrt();
                Value<T> lb = SGMIN(a, b);
                Value<T> rb = SGMAX(a, b);
                return SubInterval(lb, rb);
            }
            
            SubInterval abs() const
            {
                Value<T> t1 = m_lb.abs();
                Value<T> t2 = m_rb.abs();
                Value<T> lb = (SGMIN(m_lb, m_rb) <= 0.0 && SGMAX(m_lb, m_rb) >=0) ? 0.0 : SGMIN(t1, t2);
                Value<T> rb = SGMAX(t1, t2);
                return SubInterval(lb, rb);
            }
            
            SubInterval ln() const
            {
                Value<T> a = m_lb.ln();
                Value<T> b = m_rb.ln();
                Value<T> lb = SGMIN(a, b);
                Value<T> rb = SGMAX(a, b);
                return SubInterval(lb, rb);
            }

            SubInterval log(T base) const
            {
                Value<T> a = m_lb.log(base);
                Value<T> b = m_rb.log(base);
                Value<T> lb = SGMIN(a, b);
                Value<T> rb = SGMAX(a, b);
                return SubInterval(lb, rb);
            }

            bool operator<(const SubInterval &pv)const { return m_lb < pv.m_lb; }

            SubInterval& Union(const SubInterval& pv)
            {
                if (RightBound() >= pv.LeftBound())
                {
                        if (pv.RightBound() >= RightBound())
                                m_rb = pv.RightBound();
                }
                return *this;
            }

            friend std::ostream& operator<<(std::ostream &out, const SubInterval &x)
            {
                return out << "lb: " << x.m_lb << " rb: " << x.m_rb << '\n';
            }
            
            T lb()const {return m_lb.GetValue();}
            T rb()const {return m_rb.GetValue();}

        private:
            SubInterval div(const SubInterval &y) const
            {
                Value<T> t1 = m_lb / y.m_lb;
                Value<T> lb = t1;
                Value<T> rb = t1;

                Value<T> t2 = m_lb / y.m_rb;
                lb = SGMIN(lb, t2);
                rb = SGMAX(rb, t2);

                Value<T> t3 = m_rb / y.m_lb;
                lb = SGMIN(lb, t3);
                rb = SGMAX(rb, t3);

                Value<T> t4 = m_rb / y.m_rb;
                lb = SGMIN(lb, t4);
                rb = SGMAX(rb, t4);
                return SubInterval(lb, rb);
            }
        };              

        template<class T> class Value
        {                
            public:
                Value(T value = 0.0) : m_value(value), m_type(TypeValue::Real), m_isMinusZero(false)
                {                   
                }

                Value(TypeValue typeValue) : m_type(typeValue), m_isMinusZero(false)
                {
                    if(typeValue == TypeValue::PlusInf)
                        m_value = std::numeric_limits<T>::infinity();
                    else if (typeValue == TypeValue::MinusInf)
                        m_value = -std::numeric_limits<T>::infinity();
                    else
                        throw std::invalid_argument("Invalid constructor initialization for class Value<T>.");             
                }

                T GetValue() const { return m_value;}
                TypeValue GetType() const { return m_type;}
                bool IsMinusZero() const { return m_isMinusZero; }
                          
                Value<T> operator+(const Value<T> &v)const { return  calc.CalValues(m_type, v.GetType())->sum(*this, v);}
                Value<T> operator-(const Value<T> &v)const { return  calc.CalValues(m_type, v.GetType())->sub(*this, v);}
                Value<T> operator*(const Value<T> &v)const { return calc.CalValues(m_type, v.GetType())->mul(*this,v); }
                Value<T> operator/(const Value<T> &v)const { return calc.CalValues(m_type, v.GetType())->div(*this, v); }
                Value<T> exp() const { return calc.CalValue(m_type)->exp(*this); }
                Value<T> abs() const { return calc.CalValue(m_type)->abs(*this); }
                Value<T> sqrt() const { return calc.CalValue(m_type)->sqrt(*this); }
                Value<T> ln() const { return calc.CalValue(m_type)->ln(*this); }
                Value<T> log(T base) const { return calc.CalValue(m_type)->log(*this, base); }

                bool operator<(const Value<T> &v)const {return m_value < v.m_value; }
                bool operator>(const Value<T> &v)const {return m_value > v.m_value; }
                bool operator<=(const Value<T> &v)const {return m_value <= v.m_value; }
                bool operator>=(const Value<T> &v)const {return m_value >= v.m_value; }
                bool operator==(const Value<T> &v)const {return m_value == v.m_value; }

                friend std::ostream& operator<<(std::ostream &out, const Value<T> &v){ return v.calc.CalValue(v.GetType())->prn(out, v);}
                friend class SubInterval<T>;
                
            private:
                T m_value;
                TypeValue m_type;
                bool m_isMinusZero;
                static CalcFactory<T> calc;
        };

        template<class T> CalcFactory<T> Value<T>::calc;      
        
    }
}

#endif /* INTERVAL__HPP */

