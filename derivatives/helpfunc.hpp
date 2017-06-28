/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   helpfunc.hpp
 * Author: alusov
 *
 * Created on June 27, 2017, 1:22 PM
 */

#ifndef HELPFUNC_HPP
#define HELPFUNC_HPP

#include <cmath>
#include "interval/interval_air.hpp"

using namespace snowgoose::interval;

namespace snowgoose {
  namespace derivative {
      
    template<class T> T cos_(const T &t)
    {
        return std::cos(t);
    }

    template<class T> Interval<T> cos_(const Interval<T> &t)
    {
        return interval::cos(t);
    }
    
    template<class T> T sin_(const T &t)
    {
        return std::sin(t);
    }

    template<class T> Interval<T> sin_(const Interval<T> &t)
    {
        return interval::sin(t);
    }
    
    template<class T> T acos_(const T &t)
    {
        return std::acos(t);
    }

    template<class T> Interval<T> acos_(const Interval<T> &t)
    {
        return interval::acos(t);
    }
    
    template<class T> T asin_(const T &t)
    {
        return std::asin(t);
    }

    template<class T> Interval<T> asin_(const Interval<T> &t)
    {
        return interval::asin(t);
    }
    
    template<class T> T tg_(const T &t)
    {
        return std::tan(t);
    }

    template<class T> Interval<T> tg_(const Interval<T> &t)
    {
        return interval::tg(t);
    }
    
    template<class T> T sqrt_(const T &t)
    {
        return std::sqrt(t);
    }

    template<class T> Interval<T> sqrt_(const Interval<T> &t)
    {
        return interval::sqrt(t);
    }
    
    template<class T> T pow_(const T &t, int exp)
    {
        return std::pow(t, exp);
    }
    
    template<class T> Interval<T> pow_(const Interval<T> &t, int exp)
    {
        return t ^ exp;
    }

    template<class T> Interval<T> sqrt_(const Interval<T> &t, int exp)
    {
        return t ^ exp;
    }
    
    template<class T> T ln_(const T &t)
    {
        return std::log(t);
    }

    template<class T> Interval<T> ln_(const Interval<T> &t)
    {
        return interval::ln(t);
    }
    
    template<class T> T exp_(const T &t)
    {
        return std::exp(t);
    }

    template<class T> Interval<T> exp_(const Interval<T> &t)
    {
        return interval::exp(t);
    }
    
  }
}

#endif /* HELPFUNC_HPP */

