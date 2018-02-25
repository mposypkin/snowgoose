/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   matrix.hpp
 * Author: alusov
 *
 * Created on June 24, 2017, 8:46 PM
 */

#ifndef MATRIX_ONE_DIM_HPP
#define MATRIX_ONE_DIM_HPP

#include <iostream>
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <memory>

namespace snowgoose {
  namespace derhighorder {

        template <class  T> class MatrixOneDim
        {
            public:
                MatrixOneDim(const std::initializer_list<T> &lst);		
                MatrixOneDim(const std::vector<T> &matrix) : m_matrix(matrix) {}
                MatrixOneDim(std::size_t size, const T &t) : m_matrix(size, t) {}
		MatrixOneDim(){}
                T& operator[](std::size_t i);
                T item(std::size_t i) const { return m_matrix[i];}
                MatrixOneDim operator+(const MatrixOneDim &y) const;
                MatrixOneDim operator-(const MatrixOneDim &y) const;
                MatrixOneDim operator*(const T &y) const;
		T operator*(const MatrixOneDim &y) const;
		MatrixOneDim mulItems(const MatrixOneDim &x) const;
                template<class T2, class T3> friend MatrixOneDim<T2> operator*(const T3 &x, const MatrixOneDim<T2> &y);
                MatrixOneDim operator/(const T &y) const;
                template<class T2> friend std::ostream& operator<<(std::ostream & out, const MatrixOneDim<T2> &y);
                std::size_t size() const { return m_matrix.size(); }
		MatrixOneDim reverse() const { std::vector<T> matrix(m_matrix); std::reverse(matrix.begin(), matrix.end()); return MatrixOneDim(matrix); }
		MatrixOneDim subMatrix(std::size_t first, std::size_t last) const { return MatrixOneDim(std::vector<T>( m_matrix.begin() + first, m_matrix.begin() + last + 1) ); }
		static MatrixOneDim<T> seq(std::size_t first, std::size_t last); 
            private:
                std::vector<T> m_matrix;   
        };

        template<class T> MatrixOneDim<T>::MatrixOneDim(const std::initializer_list<T> &lst) : m_matrix(lst)
        {
        }
       
        template<class T> T& MatrixOneDim<T>::operator[](std::size_t i)
        {
            return m_matrix[i];
        }
        
        template<class T> MatrixOneDim<T> MatrixOneDim<T>::operator+(const MatrixOneDim &y) const 
        {
            std::size_t sz = m_matrix.size();
            std::vector<T> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i] = m_matrix[i]+y.m_matrix[i];
            return MatrixOneDim<T>(matrix);
        }

        template<class T> MatrixOneDim<T> MatrixOneDim<T>::operator-(const MatrixOneDim &y) const
        {
            std::size_t sz = m_matrix.size();
            std::vector<T> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i] = m_matrix[i]-y.m_matrix[i];
            return MatrixOneDim<T>(matrix);
        }

        template<class T> MatrixOneDim<T> MatrixOneDim<T>::operator*(const T &y) const
        {
            std::size_t sz = m_matrix.size();
            std::vector<T> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i]=m_matrix[i]*y;
            return matrix;
        }

	template<class T> T MatrixOneDim<T>::operator*(const MatrixOneDim &y) const
	{
            std::size_t sz = m_matrix.size();
	    if(sz != y.m_matrix.size())
		throw std::invalid_argument("Invalid operation. Sizes of matrixs are nor equal.");	    
            T t = 0.0;
            for(int i=0; i < sz; i++)
            { 
                t += m_matrix[i] * y.m_matrix[i];
            }
            return t;		
	}

	template<class T> MatrixOneDim<T> MatrixOneDim<T>::mulItems(const MatrixOneDim &x) const
	{
            std::size_t sz = x.m_matrix.size();
	    if(sz != m_matrix.size())
		throw std::invalid_argument("Invalid operation. Sizes of matrixs are nor equal.");

	    std::vector<T> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i] = m_matrix[i] * x.m_matrix[i];

            return MatrixOneDim<T>(matrix);	    
	}

        template<class T2, class T3> MatrixOneDim<T2> operator*(const T3 &x, const MatrixOneDim<T2> &y)
        {
            std::size_t sz = y.m_matrix.size();
            std::vector<T2> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i]=y.m_matrix[i]*x;
            return matrix;
        }

        template<class T> MatrixOneDim<T> MatrixOneDim<T>::operator/(const T &y) const //error control
        {
            std::size_t sz = m_matrix.size();
            std::vector<T> matrix(sz, 0.0);
            for(int i=0; i < sz; i++)
                matrix[i]=m_matrix[i]/y;
            return matrix;
        }

        template<class T2> std::ostream& operator<<(std::ostream & out, const MatrixOneDim<T2> &y)
        {
            std::size_t sz = y.m_matrix.size();
            for(int i=0; i< sz; i++)
                std::cout << y.m_matrix[i] << ' ';
            return out;
        }

        template<class T> MatrixOneDim<T> MatrixOneDim<T>::seq(std::size_t first, std::size_t last)
        {
            std::vector<T> matrix;
            for(int i = 0; i < last; i++)
		matrix.push_back(first + i);
	    return MatrixOneDim<T>(matrix);
        }

    }
}


#endif 

