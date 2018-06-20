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
 * \brief A number of commonly used algebraic functions for the localized OPM automatic
 *        differentiation (AD) framework with variable number of variables.
 *
 * This file provides AD variants of the the most commonly used functions of the <cmath>
 * header file.
 */
#ifndef OPM_LOCAL_AD_MATH_V_HPP
#define OPM_LOCAL_AD_MATH_V_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <cassert>

namespace Opm {
namespace DenseAd {
// forward declaration of the Evaluation template class
template <class ValueT>
class EvaluationV;

// provide some algebraic functions
// TODO: we might need to put the numVars back here?
template <class ValueType>
EvaluationV<ValueType> abs(const EvaluationV<ValueType>& x)
{ return (x > 0.0)?x:-x; }

template <class ValueType>
EvaluationV<ValueType> min(const EvaluationV<ValueType>& x1,
                           const EvaluationV<ValueType>& x2)
{ return (x1 < x2)?x1:x2; }

template <class Arg1ValueType, class ValueType>
EvaluationV<ValueType> min(const Arg1ValueType& x1,
                           const EvaluationV<ValueType>& x2)
{ return (x1 < x2)?x1:x2; }

template <class ValueType, class Arg2ValueType>
EvaluationV<ValueType> min(const EvaluationV<ValueType>& x1,
                           const Arg2ValueType& x2)
{ return min(x2, x1); }

template <class ValueType>
EvaluationV<ValueType> max(const EvaluationV<ValueType>& x1,
                           const EvaluationV<ValueType>& x2)
{ return (x1 > x2)?x1:x2; }

template <class Arg1ValueType, class ValueType>
EvaluationV<ValueType> max(const Arg1ValueType& x1,
                           const EvaluationV<ValueType>& x2)
{ return (x1 > x2)?x1:x2; }

template <class ValueType, class Arg2ValueType>
EvaluationV<ValueType> max(const EvaluationV<ValueType>& x1,
                           const Arg2ValueType& x2)
{ return max(x2, x1); }

template <class ValueType>
EvaluationV<ValueType> tan(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    const ValueType& tmp = ValueTypeToolbox::tan(x.value());
    result.setValue(tmp);

    // derivatives use the chain rule
    const ValueType& df_dx = 1 + tmp*tmp;
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> atan(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::atan(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1/(1 + x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> atan2(const EvaluationV<ValueType>& x,
                             const EvaluationV<ValueType>& y)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    assert(x.size == y.size);

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::atan2(x.value(), y.value()));

    // derivatives use the chain rule
    const ValueType& alpha = 1/(1 + (x.value()*x.value())/(y.value()*y.value()));
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx) {
        result.setDerivative(curVarIdx,
                             alpha/(y.value()*y.value())
                             *(x.derivative(curVarIdx)*y.value() - x.value()*y.derivative(curVarIdx)));
    }

    return result;
}

template <class ValueType>
EvaluationV<ValueType> atan2(const EvaluationV<ValueType>& x,
                             const ValueType& y)
{
    // TODO: not totally sure it is correct
    const EvaluationV<ValueType> yeval = EvaluationV<ValueType>::createConstant(y, x.size);
    return atan2(x, yeval);
}

template <class ValueType>
EvaluationV<ValueType> sin(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::sin(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = ValueTypeToolbox::cos(x.value());
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> asin(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::asin(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1.0/ValueTypeToolbox::sqrt(1 - x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> cos(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::cos(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = -ValueTypeToolbox::sin(x.value());
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> acos(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::acos(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = - 1.0/ValueTypeToolbox::sqrt(1 - x.value()*x.value());
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

template <class ValueType>
EvaluationV<ValueType> sqrt(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    const ValueType& sqrt_x = ValueTypeToolbox::sqrt(x.value());
    result.setValue(sqrt_x);

    // derivatives use the chain rule
    ValueType df_dx = 0.5/sqrt_x;
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx) {
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));
    }

    return result;
}

template <class ValueType>
EvaluationV<ValueType> exp(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    const ValueType& exp_x = ValueTypeToolbox::exp(x.value());
    result.setValue(exp_x);

    // derivatives use the chain rule
    const ValueType& df_dx = exp_x;
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

// exponentiation of arbitrary base with a fixed constant
template <class ValueType, class ExpType>
EvaluationV<ValueType> pow(const EvaluationV<ValueType>& base,
                          const ExpType& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., base.size);

    const ValueType& pow_x = ValueTypeToolbox::pow(base.value(), exp);
    result.setValue(pow_x);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result = 0.0;
    }
    else {
        // derivatives use the chain rule
        const ValueType& df_dx = pow_x/base.value()*exp;
        for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
            result.setDerivative(curVarIdx, df_dx*base.derivative(curVarIdx));
    }

    return result;
}

// exponentiation of constant base with an arbitrary exponent
template <class BaseType, class ValueType>
EvaluationV<ValueType> pow(const BaseType& base,
                          const EvaluationV<ValueType>& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., exp.size);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result.setValue(0.0);
    }
    else {
        const ValueType& lnBase = ValueTypeToolbox::log(base);
        result.setValue(ValueTypeToolbox::exp(lnBase*exp.value()));

        // derivatives use the chain rule
        const ValueType& df_dx = lnBase*result.value();
        for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
            result.setDerivative(curVarIdx, df_dx*exp.derivative(curVarIdx));
    }

    return result;
}

// this is the most expensive power function. Computationally it is pretty expensive, so
// one of the above two variants above should be preferred if possible.
template <class ValueType>
EvaluationV<ValueType> pow(const EvaluationV<ValueType>& base,
                          const EvaluationV<ValueType>& exp)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    assert(base.size == exp.size);

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., base.size);

    if (base == 0.0) {
        // we special case the base 0 case because 0.0 is in the valid range of the
        // base but the generic code leads to NaNs.
        result.setValue(0.0);
    }
    else {
        ValueType valuePow = ValueTypeToolbox::pow(base.value(), exp.value());
        result.setValue(valuePow);

        // use the chain rule for the derivatives. since both, the base and the exponent can
        // potentially depend on the variable set, calculating these is quite elaborate...
        const ValueType& f = base.value();
        const ValueType& g = exp.value();
        const ValueType& logF = ValueTypeToolbox::log(f);
        for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx) {
            const ValueType& fPrime = base.derivative(curVarIdx);
            const ValueType& gPrime = exp.derivative(curVarIdx);
            result.setDerivative(curVarIdx, (g*fPrime/f + logF*gPrime) * valuePow);
        }
    }

    return result;
}

template <class ValueType>
EvaluationV<ValueType> log(const EvaluationV<ValueType>& x)
{
    typedef MathToolbox<ValueType> ValueTypeToolbox;

    EvaluationV<ValueType> result = EvaluationV<ValueType>::createConstant(0., x.size);

    result.setValue(ValueTypeToolbox::log(x.value()));

    // derivatives use the chain rule
    const ValueType& df_dx = 1/x.value();
    for (int curVarIdx = 0; curVarIdx < result.size; ++curVarIdx)
        result.setDerivative(curVarIdx, df_dx*x.derivative(curVarIdx));

    return result;
}

} // namespace DenseAd

// a kind of traits class for the automatic differentiation case. (The toolbox for the
// scalar case is provided by the MathToolbox.hpp header file.)
template <class ValueT>
struct MathToolbox<Opm::DenseAd::EvaluationV<ValueT> >
{
private:
public:
    typedef ValueT ValueType;
    typedef Opm::MathToolbox<ValueType> InnerToolbox;
    typedef typename InnerToolbox::Scalar Scalar;
    typedef Opm::DenseAd::EvaluationV<ValueType> EvaluationV;

    static ValueType value(const EvaluationV& eval)
    { return eval.value(); }

    static decltype(InnerToolbox::scalarValue(0.0)) scalarValue(const EvaluationV& eval)
    { return InnerToolbox::scalarValue(eval.value()); }

    static EvaluationV createConstant(ValueType value, const int numVars)
    { return EvaluationV::createConstant(value, numVars); }

    static EvaluationV createVariable(ValueType value, const int varIdx, const int numVars)
    { return EvaluationV::createVariable(value, varIdx, numVars); }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<EvaluationV, LhsEval>::value,
                                   LhsEval>::type
    decay(const EvaluationV& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_same<EvaluationV, LhsEval>::value,
                                   LhsEval>::type
    decay(const EvaluationV&& eval)
    { return eval; }

    template <class LhsEval>
    static typename std::enable_if<std::is_floating_point<LhsEval>::value,
                                   LhsEval>::type
    decay(const EvaluationV& eval)
    { return eval.value(); }

    // comparison
    static bool isSame(const EvaluationV& a, const EvaluationV& b, Scalar tolerance)
    {
        typedef MathToolbox<ValueType> ValueTypeToolbox;

        if (a.size != b.size)
            return false;

        // make sure that the value of the evaluation is identical
        if (!ValueTypeToolbox::isSame(a.value(), b.value(), tolerance))
            return false;

        // make sure that the derivatives are identical
        for (int curVarIdx = 0; curVarIdx < a.size; ++curVarIdx)
            if (!ValueTypeToolbox::isSame(a.derivative(curVarIdx), b.derivative(curVarIdx), tolerance))
                return false;

        return true;
    }

    // arithmetic functions
    template <class Arg1Eval, class Arg2Eval>
    static EvaluationV max(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return Opm::DenseAd::max(arg1, arg2); }

    /* static EvaluationV max(const EvaluationV& arg1, const ValueType& arg2)
    { return Opm::DenseAd::max(arg1, arg2); } */

    template <class Arg1Eval, class Arg2Eval>
    static EvaluationV min(const Arg1Eval& arg1, const Arg2Eval& arg2)
    { return Opm::DenseAd::min(arg1, arg2); }

    /* static EvaluationV min(const EvaluationV& arg1, const ValueType& arg2)
    { return Opm::DenseAd::min(arg1, arg2); } */

    static EvaluationV abs(const EvaluationV& arg)
    { return Opm::DenseAd::abs(arg); }

    static EvaluationV tan(const EvaluationV& arg)
    { return Opm::DenseAd::tan(arg); }

    static EvaluationV atan(const EvaluationV& arg)
    { return Opm::DenseAd::atan(arg); }

    static EvaluationV atan2(const EvaluationV& arg1, const EvaluationV& arg2)
    { return Opm::DenseAd::atan2(arg1, arg2); }

    static EvaluationV atan2(const EvaluationV& arg1, const ValueType& arg2)
    { return Opm::DenseAd::atan2(arg1, arg2); }

    static EvaluationV sin(const EvaluationV& arg)
    { return Opm::DenseAd::sin(arg); }

    static EvaluationV asin(const EvaluationV& arg)
    { return Opm::DenseAd::asin(arg); }

    static EvaluationV cos(const EvaluationV& arg)
    { return Opm::DenseAd::cos(arg); }

    static EvaluationV acos(const EvaluationV& arg)
    { return Opm::DenseAd::acos(arg); }

    static EvaluationV sqrt(const EvaluationV& arg)
    { return Opm::DenseAd::sqrt(arg); }

    static EvaluationV exp(const EvaluationV& arg)
    { return Opm::DenseAd::exp(arg); }

    static EvaluationV log(const EvaluationV& arg)
    { return Opm::DenseAd::log(arg); }

    template <class RhsValueType>
    static EvaluationV pow(const EvaluationV& arg1, const RhsValueType& arg2)
    { return Opm::DenseAd::pow(arg1, arg2); }

    template <class RhsValueType>
    static EvaluationV pow(const RhsValueType& arg1, const EvaluationV& arg2)
    { return Opm::DenseAd::pow(arg1, arg2); }

    static EvaluationV pow(const EvaluationV& arg1, const EvaluationV& arg2)
    { return Opm::DenseAd::pow(arg1, arg2); }

    static bool isfinite(const EvaluationV& arg)
    {
        if (!InnerToolbox::isfinite(arg.value()))
            return false;

        for (int i = 0; i < arg.size; ++i)
            if (!InnerToolbox::isfinite(arg.derivative(i)))
                return false;

        return true;
    }

    static bool isnan(const EvaluationV& arg)
    {
        if (InnerToolbox::isnan(arg.value()))
            return true;

        for (int i = 0; i < arg.size; ++i)
            if (InnerToolbox::isnan(arg.derivative(i)))
                return true;

        return false;
    }
};

}

#endif
