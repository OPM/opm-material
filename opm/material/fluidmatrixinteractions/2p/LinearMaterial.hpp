// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Copyright (C) 2010 by Philipp Nuske                                     *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::LinearMaterial
 */
#ifndef OPM_LINEAR_MATERIAL_HH
#define OPM_LINEAR_MATERIAL_HH

#include "LinearMaterialParams.hpp"

#include <algorithm>

namespace Opm {
/*!
 * \ingroup fluidmatrixinteractionslaws
 *
 * \brief Linear capillary pressure and relative permeability <->
 *        saturation relations
 *
 * The entry pressure is reached at \f$ \overline S_w = 1\f$, the maximum
 * capillary pressure is observed at \f$ \overline S_w = 0\f$.
 *
 * For general info see EffToAbsLaw
 *
 * \see LinearMaterialParams
 */
template <class ScalarT, class ParamsT = LinearMaterialParams<ScalarT> >
class LinearMaterial
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param Swe Effective saturation of the wetting phase
     *            \f$\overline{S}_w\f$ conversion from absolute
     *            saturation happened in EffToAbsLaw.
     * \param params A object that provides the appropriate coefficients
     *               for the respective law.
     *
     * \return Capillary pressure calculated by linear constitutive
     *         relation.
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        return (1 - Swe)*(params.maxPC() - params.entryPC()) + params.entryPC();
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
         S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC Capillary pressure \f$p_C\f$
     * \param params A container object that is populated with the
     *               appropriate coefficients for the respective law.
     *
     * \return Effective wetting phase saturation calculated as
     *         inverse of the linear constitutive relation.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        return 1 - (pC - params.entryPC())/(params.maxPC() - params.entryPC());
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t. the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     \f]
     *
     * \param Swe Effective saturation of the wetting phase
     *            \f$\overline{S}_w\f$ conversion from absolute
     *            saturation happened in EffToAbsLaw.
     * \param params A object that provides the appropriate coefficients
     *               for the respective law.
     *
     * \return Partial derivative of \f$p_c\f$ w.r.t. effective
     *         saturation according to linear material relation.
     */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        return - (params.maxPC() - params.entryPC());
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     *
     * \param pC Capillary pressure \f$p_C\f$ [Pa]
     * \param params A object that provides the appropriate coefficients
     *               for the respective law.
     *
     * \return Partial derivative of effective saturation
     *           w.r.t. \f$p_c\f$ according to linear relation.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return - 1/(params.maxPC() - params.entryPC());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params A object that provides the appropriate coefficients
     *               for the respective law.
     * \param Swe Effective saturation of the wetting phase
     *            \f$\overline{S}_w\f$. The conversion from absolute
     *            saturation is done by EffToAbsLaw.
     *
     * \return Relative permeability of the wetting phase calculated
     *         as linear relation.
     */
    static Scalar krw(const Params &params, Scalar Swe)
    {
        return std::max(std::min(Swe,1.0),0.0);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     *
     * \param params A object that provides the appropriate coefficients
     *               for the respective law.
     * \param Swe Effective saturation of the wetting phase
     *            \f$\overline{S}_w\f$. The conversion from absolute
     *            saturation is done by EffToAbsLaw.
     *
     * \return Relative permeability of the non-wetting phase
     *         calculated as linear relation.
     */
    static Scalar krn(const Params &params, Scalar Swe)
    {
        Scalar Sne = 1 - Swe;
        return std::max(std::min(Sne,1.0),0.0);
    }
};
}

#endif
