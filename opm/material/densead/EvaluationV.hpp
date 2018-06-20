// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Representation of an evaluation of a function and its derivatives w.r.t. a set
 *        of variables in the localized OPM automatic differentiation (AD) framework.
 */
#ifndef OPM_DENSEAD_EVALUATION_V_HPP
#define OPM_DENSEAD_EVALUATION_V_HPP

#include "MathV.hpp"

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/version.hh>

#include <vector>
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <algorithm>

namespace Opm {
namespace DenseAd {

/*!
 * \brief Represents a function evaluation and its derivatives w.r.t. a fixed set of
 *        variables.
 */
template <class ValueT>
class EvaluationV
{
public:
    //! field type
    typedef ValueT ValueType;

    //! number of derivatives
    int size = 0;

protected:
    //! length of internal data vector
    int length_ = 0;

    //! position index for value
    static constexpr int valuepos_ = 0;
    //! start index for derivatives
    static constexpr int dstart_ = 1;
    //! end+1 index for derivatives
    int dend_ = 0;

public:
    //! default constructor
    // TODO: maybe we should make a default usable one?
    EvaluationV() = default;

    // TODO: this is rather dangerous to have
    // put here just help to make ReturnEval_(EvaluationV, double) return EvaluationV.
    // we should not use it in the calculation at all.
    EvaluationV(const ValueType& value)
     : size(0)
     , length_(size + 1)
     , dend_(length_)
     , data_(length_, 0.0)
    {
        throw std::runtime_error(" You should not convert a Scalar to EvaluationV without the number of variables at all! ");
        setValue(value);
    }


    //! when resizing, all the data will be cleaned
    void resize(const int numVars) {
        size = numVars;
        length_ = size + 1;
        dend_ = length_;
        data_.resize(length_, 0.0);
    }

    //! copy other function evaluation
    EvaluationV(const EvaluationV& other) = default;

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    template <class RhsValueType>
    EvaluationV(const RhsValueType& c, const int numVars)
    {
        resize(numVars);
        setValue( c );
        // clearDerivatives();
        Valgrind::CheckDefined( data_ );
    }

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    template <class RhsValueType>
    EvaluationV(const RhsValueType& c, const int varPos, const int numVars)
    {
        resize(numVars);

        // The variable position must be in represented by the given variable descriptor
        assert(0 <= varPos && varPos < size);

        setValue( c );
        // clearDerivatives();

        data_[varPos + dstart_] = 1.0;
        Valgrind::CheckDefined(data_);
    }

    // set all derivatives to zero
    void clearDerivatives()
    {
        for (int i = dstart_; i < dend_; ++i) {
            data_[i] = 0.0;
        }
    }

    // create a function evaluation for a "naked" depending variable (i.e., f(x) = x)
    template <class RhsValueType>
    static EvaluationV createVariable(const RhsValueType& value, const int varPos, const int numVars)
    {
        // copy function value and set all derivatives to 0, except for the variable
        // which is represented by the value (which is set to 1.0)
        return EvaluationV( value, varPos, numVars );
    }

    // "evaluate" a constant function (i.e. a function that does not depend on the set of
    // relevant variables, f(x) = c).
    template <class RhsValueType>
    static EvaluationV createConstant(const RhsValueType& value, const int numVars)
    {
        return EvaluationV( value, numVars );
    }

    // print the value and the derivatives of the function evaluation
    void print(std::ostream& os = std::cout) const
    {
        // print value
        os << "v: " << value() << " / d:";

        // print derivatives
        for (int varIdx = 0; varIdx < size; ++varIdx) {
            os << " " << derivative(varIdx);
        }
    }

    // copy all derivatives from other
    void copyDerivatives(const EvaluationV& other)
    {
        assert(size == other.size);

        for (int i = dstart_; i < dend_; ++i) {
            data_[i] = other.data_[i];
        }
    }


    // add value and derivatives from other to this values and derivatives
    EvaluationV& operator+=(const EvaluationV& other)
    {
        assert(size == other.size);

        for (int i = 0; i < length_; ++i) {
            data_[i] += other.data_[i];
        }

        return *this;
    }

    // add value from other to this values
    template <class RhsValueType>
    EvaluationV& operator+=(const RhsValueType& other)
    {
        // value is added, derivatives stay the same
        data_[valuepos_] += other;

        return *this;
    }

    // subtract other's value and derivatives from this values
    EvaluationV& operator-=(const EvaluationV& other)
    {
        assert(size == other.size);

        for (int i = 0; i < length_; ++i) {
            data_[i] -= other.data_[i];
        }

        return *this;
    }

    // subtract other's value from this values
    template <class RhsValueType>
    EvaluationV& operator-=(const RhsValueType& other)
    {
        // for constants, values are subtracted, derivatives stay the same
        data_[ valuepos_ ] -= other;

        return *this;
    }

    // multiply values and apply chain rule to derivatives: (u*v)' = (v'u + u'v)
    EvaluationV& operator*=(const EvaluationV& other)
    {
        assert(size == other.size);

        // while the values are multiplied, the derivatives follow the product rule,
        // i.e., (u*v)' = (v'u + u'v).
        const ValueType u = this->value();
        const ValueType v = other.value();

        // value
        data_[valuepos_] *= v ;

        //  derivatives
        for (int i = dstart_; i < dend_; ++i) {
            data_[i] = data_[i] * v + other.data_[i] * u;
        }

        return *this;
    }

    // m(c*u)' = c*u'
    template <class RhsValueType>
    EvaluationV& operator*=(const RhsValueType& other)
    {
        for (int i = 0; i < length_; ++i) {
            data_[i] *= other;
        }

        return *this;
    }

    // m(u*v)' = (vu' - uv')/v^2
    EvaluationV& operator/=(const EvaluationV& other)
    {
        assert(size == other.size);
        // values are divided, derivatives follow the rule for division, i.e., (u/v)' = (v'u -
        // u'v)/v^2.
        ValueType& u = data_[ valuepos_ ];
        const ValueType& v = other.value();
        for (int idx = dstart_; idx < dend_; ++idx) {
            const ValueType& uPrime = data_[idx];
            const ValueType& vPrime = other.data_[idx];

            data_[idx] = (v*uPrime - u*vPrime)/(v*v);
        }
        u /= v;

        return *this;
    }

    // divide value and derivatives by value of other
    template <class RhsValueType>
    EvaluationV& operator/=(const RhsValueType& other)
    {
        const ValueType tmp = 1.0/other;

        for (int i = 0; i < length_; ++i) {
            data_[i] *= tmp;
        }

        return *this;
    }

    // add two evaluation objects
    EvaluationV operator+(const EvaluationV& other) const
    {
        EvaluationV result(*this);

        result += other;

        return result;
    }

    // add constant to this object
    template <class RhsValueType>
    EvaluationV operator+(const RhsValueType& other) const
    {
        EvaluationV result(*this);

        result += other;

        return result;
    }

    // subtract two evaluation objects
    EvaluationV operator-(const EvaluationV& other) const
    {
        EvaluationV result(*this);

        result -= other;

        return result;
    }

    // subtract constant from evaluation object
    template <class RhsValueType>
    EvaluationV operator-(const RhsValueType& other) const
    {
        EvaluationV result(*this);

        result -= other;

        return result;
    }

    // negation (unary minus) operator
    EvaluationV operator-() const
    {
        EvaluationV result(*this);

        // set value and derivatives to negative
        for (int i = 0; i < length_; ++i) {
            result.data_[i] = - data_[i];
        }

        return result;
    }

    EvaluationV operator*(const EvaluationV& other) const
    {
        EvaluationV result(*this);

        result *= other;

        return result;
    }

    template <class RhsValueType>
    EvaluationV operator*(const RhsValueType& other) const
    {
        EvaluationV result(*this);

        result *= other;

        return result;
    }

    EvaluationV operator/(const EvaluationV& other) const
    {
        EvaluationV result(*this);

        result /= other;

        return result;
    }

    template <class RhsValueType>
    EvaluationV operator/(const RhsValueType& other) const
    {
        EvaluationV result(*this);

        result /= other;

        return result;
    }

    template <class RhsValueType>
    EvaluationV& operator=(const RhsValueType& other)
    {
        setValue( other );
        clearDerivatives();

        return *this;
    }

    // copy assignment from evaluation
    EvaluationV& operator=(const EvaluationV& other) = default;

    template <class RhsValueType>
    bool operator==(const RhsValueType& other) const
    { return value() == other; }

    bool operator==(const EvaluationV& other) const
    {
        if (size != other.size) {
            return false;
        }

        for (int idx = 0; idx < length_; ++idx) {
            if (data_[idx] != other.data_[idx]) {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const EvaluationV& other) const
    { return !operator==(other); }

    template <class RhsValueType>
    bool operator!=(const RhsValueType& other) const
    { return !operator==(other); }


    template <class RhsValueType>
    bool operator>(RhsValueType other) const
    { return value() > other; }

    bool operator>(const EvaluationV& other) const
    { return value() > other.value(); }

    template <class RhsValueType>
    bool operator<(RhsValueType other) const
    { return value() < other; }

    bool operator<(const EvaluationV& other) const
    { return value() < other.value(); }

    template <class RhsValueType>
    bool operator>=(RhsValueType other) const
    { return value() >= other; }

    bool operator>=(const EvaluationV& other) const
    { return value() >= other.value(); }

    template <class RhsValueType>
    bool operator<=(RhsValueType other) const
    { return value() <= other; }

    bool operator<=(const EvaluationV& other) const
    { return value() <= other.value(); }

    // return value of variable
    const ValueType& value() const
    { return data_[valuepos_]; }

    // set value of variable
    template <class RhsValueType>
    void setValue(const RhsValueType& val)
    { data_[valuepos_] = val; }

    // return varIdx'th derivative
    const ValueType& derivative(const int varIdx) const
    {
        assert(0 <= varIdx && varIdx < size);

        return data_[dstart_ + varIdx];
    }

    // set derivative at position varIdx
    void setDerivative(const int varIdx, const ValueType& derVal)
    {
        assert(0 <= varIdx && varIdx < size);

        data_[dstart_ + varIdx] = derVal;
    }

private:
    std::vector<ValueT> data_;
};

// the generic operators are only required for the unspecialized case
template <class RhsValueType, class ValueType>
bool operator<(const RhsValueType& a, const EvaluationV<ValueType>& b)
{ return b > a; }

template <class RhsValueType, class ValueType>
bool operator>(const RhsValueType& a, const EvaluationV<ValueType>& b)
{ return b < a; }

template <class RhsValueType, class ValueType>
bool operator<=(const RhsValueType& a, const EvaluationV<ValueType>& b)
{ return b >= a; }

template <class RhsValueType, class ValueType>
bool operator>=(const RhsValueType& a, const EvaluationV<ValueType>& b)
{ return b <= a; }

template <class RhsValueType, class ValueType>
bool operator!=(const RhsValueType& a, const EvaluationV<ValueType>& b)
{ return a != b.value(); }

template <class RhsValueType, class ValueType>
EvaluationV<ValueType> operator+(const RhsValueType& a, const EvaluationV<ValueType>& b)
{
    EvaluationV<ValueType> result(b);
    result += a;
    return result;
}

template <class RhsValueType, class ValueType>
EvaluationV<ValueType> operator-(const RhsValueType& a, const EvaluationV<ValueType>& b)
{
    EvaluationV<ValueType> result(a, b.size);
    result -= b;
    return result;
}

template <class RhsValueType, class ValueType>
EvaluationV<ValueType> operator/(const RhsValueType& a, const EvaluationV<ValueType>& b)
{
    EvaluationV<ValueType> tmp(a, b.size);
    tmp /= b;
    return tmp;
}

template <class RhsValueType, class ValueType>
EvaluationV<ValueType> operator*(const RhsValueType& a, const EvaluationV<ValueType>& b)
{
    EvaluationV<ValueType> result(b);
    result *= a;
    return result;
}

template <class ValueType>
std::ostream& operator<<(std::ostream& os, const EvaluationV<ValueType>& eval)
{
    os << eval.value();
    return os;
}
} } // namespace DenseAd, Opm

// In Dune 2.3, the EvaluationV.hpp header must be included before the fmatrix.hh
// header. Dune 2.4+ does not suffer from this because of some c++-foo.
//
// for those who are wondering: in C++ function templates cannot be partially
// specialized, and function argument overloads must be known _before_ they are used. The
// latter is what we do for the 'Dune::fvmeta::absreal()' function.
//
// consider the following test program:
//
// double foo(double i)
// { return i; }
//
// void bar()
// { std::cout << foo(0) << "\n"; }
//
// int foo(int i)
// { return i + 1; }
//
// void foobar()
// { std::cout << foo(0) << "\n"; }
//
// this will print '0' for bar() and '1' for foobar()...
#if !(DUNE_VERSION_NEWER(DUNE_COMMON, 2,4))

namespace Opm {
namespace DenseAd {
template <class ValueType>
EvaluationV<ValueType> abs(const EvaluationV<ValueType>&);
}}

namespace std {
template <class ValueType>
const Opm::DenseAd::EvaluationV<ValueType> abs(const Opm::DenseAd::EvaluationV<ValueType>& x)
{ return Opm::DenseAd::abs(x); }

} // namespace std

#if defined DUNE_DENSEMATRIX_HH
#warning \
 "Due to some C++ peculiarity regarding function overloads, the 'EvaluationV.hpp'" \
 "header file must be included before Dune's 'densematrix.hh' for Dune < 2.4. " \
 "(If Evaluations are to be used in conjunction with a dense matrix.)"
#endif

#endif

// this makes the Dune matrix/vector classes happy...
#include <dune/common/ftraits.hh>

namespace Dune {
template <class ValueType>
struct FieldTraits<Opm::DenseAd::EvaluationV<ValueType> >
{
public:
    typedef Opm::DenseAd::EvaluationV<ValueType> field_type;
    // setting real_type to field_type here potentially leads to slightly worse
    // performance, but at least it makes things compile.
    typedef field_type real_type;
};

} // namespace Dune

#endif // OPM_DENSEAD_EVALUATION_V_HPP
