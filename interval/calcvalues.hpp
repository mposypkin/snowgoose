/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   calcvalues.hpp
 * Author: Alexander Usov
 *
 * Created on January 28, 2017, 11:45 AM
 */

#ifndef CALCVALUES_HPP
#define CALCVALUES_HPP

#include <exception>
#include <unordered_map>
#include <math.h>


namespace snowgoose {
    namespace interval {

        enum TypeValue
        {
            Real = 1,
            PlusInf = 2,
            MinusInf = 4
        };      
           
        template<class T> class Value;
        
        template<class T> class CalcValues
        {
            public:
                virtual Value<T> sum(const Value<T>& lv, const Value<T>& rv) = 0;
                virtual Value<T> sub(const Value<T>& lv, const Value<T>& rv) = 0;
                virtual Value<T> mul(const Value<T>& lv, const Value<T>& rv) = 0;
                virtual Value<T> div(const Value<T>& lv, const Value<T>& rv) = 0;
                virtual ~CalcValues(){};
        };
        
        template<class T> class CalcValue
        {
            public:
                virtual Value<T> exp(const Value<T>& v) = 0;
                virtual Value<T> abs(const Value<T>& v) = 0;
                virtual Value<T> sqrt(const Value<T>& v) = 0;
                virtual Value<T> ln(const Value<T>& v) = 0;
                virtual Value<T> log(const Value<T>& v, T base) = 0;
                virtual std::ostream& prn(std::ostream &out, const Value<T> &v) = 0;
                virtual ~CalcValue(){};
        };

        template<class T> class RealValue : public CalcValue<T>
        {
        public:
            Value<T> exp(const Value<T>& v) { return ::exp(v.GetValue()); }
            Value<T> abs(const Value<T>& v) { return ::abs(v.GetValue()); }
            std::ostream& prn(std::ostream &out, const Value<T> &v) { return out << v.GetValue(); }

            Value<T> sqrt(const Value<T>& v) { 
                if (v.GetValue() < 0)
                    std::invalid_argument("The function sqrt is not defined for negative numbers.");
                return ::sqrt(v.GetValue());
            };
            Value<T> ln(const Value<T>& v){
                if (v.GetValue() < 0)
                    throw std::invalid_argument("The function ln is not define for negative numbers");
                return ::log(v.GetValue());
            }
            Value<T> log(const Value<T>& v, T base) {
                if (v.GetValue() < 0)
                    throw std::invalid_argument("The function log is not define for negative numbers");
                if (base <= 0.0 || base == 1.0)
                    throw std::invalid_argument("The function log is not define for negative base and if base equals 1.0");
                if (v.GetValue() == 0.0)
                {
                    if (base < 1.0)
                        return Value<T>(TypeValue::MinusInf);
                    else
                        return Value<T>(TypeValue::PlusInf);
                }
                else
                    return ::log(v.GetValue()) / ::log(base);
            }

        };

        template<class T> class PlusInfValue : public CalcValue<T>
        {
            virtual Value<T> exp(const Value<T>& v) { return Value<T>(TypeValue::PlusInf); }
            virtual Value<T> abs(const Value<T>& v) { return Value<T>(TypeValue::PlusInf); }
            virtual Value<T> sqrt(const Value<T>& v) { return Value<T>(TypeValue::PlusInf); }
            virtual Value<T> ln(const Value<T>& v) { return Value<T>(TypeValue::PlusInf); }
            virtual Value<T> log(const Value<T>& v, T base) { return Value<T>(TypeValue::PlusInf); }
            std::ostream& prn(std::ostream &out, const Value<T> &v) { return out << "inf+"; }
        };

        template<class T> class MinusInfValue : public CalcValue<T>
        {
            virtual Value<T> exp(const Value<T>& v) { return Value<T>(0.0); }
            virtual Value<T> abs(const Value<T>& v) { return Value<T>(TypeValue::PlusInf); }
            virtual Value<T> sqrt(const Value<T>& v) { throw std::invalid_argument("The function sqrt is not defined for negative numbers."); }
            virtual Value<T> ln(const Value<T>& v) { throw std::invalid_argument("The function ln is not define for negative numbers"); }
            virtual Value<T> log(const Value<T>& v, T base) { throw std::invalid_argument("The function log is not define for negative numbers"); }
            std::ostream& prn(std::ostream &out, const Value<T> &v) { return out << "inf-"; }
        };
        
        template<class T> class LValueRValue : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv,const Value<T>& rv){return Value<T>(lv.GetValue() + rv.GetValue());}
                Value<T> sub(const Value<T>& lv,const Value<T>& rv){return Value<T>(lv.GetValue() - rv.GetValue());}
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){return Value<T>(lv.GetValue() * rv.GetValue());}
                Value<T> div(const Value<T>& lv, const Value<T>& rv) {
                    if (rv.GetValue() == 0.0)
                        return ((lv.GetValue() < 0.0 && rv.IsMinusZero()) || (lv.GetValue() > 0.0 && rv.IsMinusZero() == false)) ? Value<T>(TypeValue::PlusInf) : Value<T>(TypeValue::MinusInf);
                    return Value<T>(lv.GetValue() / rv.GetValue()); 
                }
        };

        template<class T> class LValueRPlusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::PlusInf);}
                Value<T> sub(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::MinusInf);}
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){
                    return lv.GetValue()==0.0 ? Value<T>(0.0) : lv.GetValue() < 0.0 ?  Value<T>(TypeValue::MinusInf) : Value<T>(TypeValue::PlusInf);
                }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { return Value<T>(0.0); }
        };
        
        template<class T> class LValueRMinusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::MinusInf); }
                Value<T> sub(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::PlusInf); }
                Value<T> mul(const Value<T>& lv, const Value<T>& rv) { 
                    return lv.GetValue() == 0.0 ? Value<T>(0.0) : lv.GetValue() < 0 ? Value<T>(TypeValue::PlusInf) : Value<T>(TypeValue::MinusInf); 
                }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { return Value<T>(0.0); }
        };

        template<class T> class LPlusInfRValue : public LValueRPlusInf<T>
        {
            public:
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){
                    return rv.GetValue() == 0.0 ? Value<T>(0.0) : rv.GetValue() < 0 ?  Value<T>(TypeValue::MinusInf) : Value<T>(TypeValue::PlusInf);
                }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { return Value<T>(rv < 0 ? TypeValue::MinusInf : TypeValue::PlusInf); }
        };
        
        template<class T> class LPlusInfRPlusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::PlusInf);}
                Value<T> sub(const Value<T>& lv,const Value<T>& rv){throw std::invalid_argument("Invalid operation. +Inf and +Inf were subtracted.");}
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::PlusInf);}
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. +Inf and +Inf were divided."); }
        };

        template<class T> class LPlusInfRMinusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. +Inf and -Inf were summed."); }
                Value<T> sub(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::PlusInf); }
                Value<T> mul(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::MinusInf); }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. +Inf and -Inf were divided."); }
        };

        template<class T> class LMinusInfRValue : public LValueRMinusInf<T>
        {
            public:
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){
                    return rv.GetValue() == 0.0 ? Value<T>(0.0) : rv.GetValue() < 0 ?  Value<T>(TypeValue::PlusInf) : Value<T>(TypeValue::MinusInf);
                }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { return Value<T>(rv < 0 ? TypeValue::PlusInf : TypeValue::MinusInf); }
        };
        
        template<class T> class LMinusInfRPlusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv,const Value<T>& rv){throw std::invalid_argument("Invalid operation. -Inf and +Inf were summed.");}
                Value<T> sub(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::MinusInf);}
                Value<T> mul(const Value<T>& lv,const Value<T>& rv){return Value<T>(TypeValue::MinusInf);}
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. -Inf and +Inf were divided."); }
        };

        template<class T> class LMinusInfRMinusInf : public CalcValues<T>
        {
            public:
                Value<T> sum(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::MinusInf); }
                Value<T> sub(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. -Inf and -Inf were subtracted."); }
                Value<T> mul(const Value<T>& lv, const Value<T>& rv) { return Value<T>(TypeValue::PlusInf); }
                Value<T> div(const Value<T>& lv, const Value<T>& rv) { throw std::invalid_argument("Invalid operation. -Inf and -Inf were divided."); }
        };

        
        template<class T> class CalcFactory
        {
            private:
                static const int l = 8;
                std::unordered_map<int, CalcValues<T>*> m_map = { 
                    { (TypeValue::Real << l) | TypeValue::Real , new LValueRValue<T>() },
                    { (TypeValue::Real << l) | TypeValue::PlusInf , new LValueRPlusInf<T>() },
                    { (TypeValue::Real << l) | TypeValue::MinusInf, new LValueRMinusInf<T>() },
                    { (TypeValue::PlusInf << l) | TypeValue::Real, new LPlusInfRValue<T>() },
                    { (TypeValue::PlusInf << l) | TypeValue::PlusInf, new LPlusInfRPlusInf<T>() },
                    { (TypeValue::PlusInf << l) | TypeValue::MinusInf, new LPlusInfRMinusInf<T>() },
                    { (TypeValue::MinusInf << l) | TypeValue::Real, new LMinusInfRValue<T>() },
                    { (TypeValue::MinusInf << l) | TypeValue::PlusInf, new LMinusInfRPlusInf<T>() },
                    { (TypeValue::MinusInf << l) | TypeValue::MinusInf, new LMinusInfRMinusInf<T>() },
                };
                std::unordered_map<int, CalcValue<T>*> m_cv = {
                    { TypeValue::Real , new RealValue<T>() },
                    { TypeValue::PlusInf , new PlusInfValue<T>() },
                    { TypeValue::MinusInf , new MinusInfValue<T>() }
                };
            public:
                CalcFactory() {}
                CalcValues<T>* CalValues(TypeValue ltv, TypeValue rtv)
                {
                    return m_map[(ltv << l) | rtv];
                }
                CalcValue<T>* CalValue(TypeValue tv)
                {
                    return m_cv[tv];
                }
                ~CalcFactory()
                {
                    for (auto v : m_map) {
                        delete v.second;
                    }
                    for (auto v : m_cv) {
                        delete v.second;
                    }
                }
        };		
    }
}

#endif /* CALCVALUES_HPP */

