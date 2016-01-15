// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser
  Copyright (C) 2012 by Bernd Flemisch

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
*/
/*!
 * \file
 * \copydoc Opm::PengRobinsonMixture
 */
#ifndef OPM_PENG_ROBINSON_MIXTURE_HPP
#define OPM_PENG_ROBINSON_MIXTURE_HPP

#include "PengRobinson.hpp"

#include <opm/material/Constants.hpp>

#include <iostream>

namespace Opm {
/*!
 * \brief Implements the Peng-Robinson equation of state for a
 *        mixture.
 */
template <class Scalar, class StaticParameters>
class PengRobinsonMixture
{
    enum { numComponents = StaticParameters::numComponents };
    typedef Opm::PengRobinson<Scalar> PengRobinson;

    // this class cannot be instantiated!
    PengRobinsonMixture() {}

    // the ideal gas constant
    static const Scalar R;

    // the u and w parameters as given by the Peng-Robinson EOS
    static const Scalar u;
    static const Scalar w;

public:
    /*!
     * \brief Computes molar volumes where the Peng-Robinson EOS is
     *        true.
     *
     * \return Number of solutions.
     */
    template <class MutableParams, class FluidState>
    static int computeMolarVolumes(Scalar *Vm,
                                   const MutableParams &params,
                                   unsigned phaseIdx,
                                   const FluidState &fs)
    {
        return PengRobinson::computeMolarVolumes(Vm, params, phaseIdx, fs);
    }

    /*!
     * \brief Returns the fugacity coefficient of an individual
     *        component in the phase.
     *
     * The fugacity coefficient \f$\phi_i\f$ of a component \f$i\f$ is
     * defined as
     * \f[
     f_i = \phi_i x_i \;,
     \f]
     * where \f$f_i\f$ is the component's fugacity and \f$x_i\f$ is
     * the component's mole fraction.
     *
     * See:
     *
      * R. Reid, et al.: The Properties of Gases and Liquids,
      * 4th edition, McGraw-Hill, 1987, pp. 42-44, 143-145
      */
    template <class FluidState, class Params>
    static Scalar computeFugacityCoefficient(const FluidState &fs,
                                             const Params &params,
                                             unsigned phaseIdx,
                                             unsigned compIdx)
    {
        // note that we normalize the component mole fractions, so
        // that their sum is 100%. This increases numerical stability
        // considerably if the fluid state is not physical.
        Scalar Vm = params.molarVolume(phaseIdx);

        // Calculate b_i / b
        Scalar bi_b = params.bPure(phaseIdx, compIdx) / params.b(phaseIdx);

        // Calculate the compressibility factor
        Scalar RT = R*fs.temperature(phaseIdx);
        Scalar p = fs.pressure(phaseIdx); // molar volume in [bar]
        Scalar Z = p*Vm/RT; // compressibility factor

        // Calculate A^* and B^* (see: Reid, p. 42)
        Scalar Astar = params.a(phaseIdx)*p/(RT*RT);
        Scalar Bstar = params.b(phaseIdx)*p/(RT);

        // calculate delta_i (see: Reid, p. 145)
        Scalar sumMoleFractions = 0.0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx)
            sumMoleFractions += fs.moleFraction(phaseIdx, compJIdx);
        Scalar deltai = 2*std::sqrt(params.aPure(phaseIdx, compIdx))/params.a(phaseIdx);
        Scalar tmp = 0;
        for (unsigned compJIdx = 0; compJIdx < numComponents; ++compJIdx) {
            tmp +=
                fs.moleFraction(phaseIdx, compJIdx)
                / sumMoleFractions
                * std::sqrt(params.aPure(phaseIdx, compJIdx))
                * (1.0 - StaticParameters::interactionCoefficient(compIdx, compJIdx));
        };
        deltai *= tmp;

        Scalar base =
            (2*Z + Bstar*(u + std::sqrt(u*u - 4*w))) /
            (2*Z + Bstar*(u - std::sqrt(u*u - 4*w)));
        Scalar expo =  Astar/(Bstar*std::sqrt(u*u - 4*w))*(bi_b - deltai);

        Scalar fugCoeff =
            std::exp(bi_b*(Z - 1))/std::max(Scalar(1e-9), Z - Bstar) *
            std::pow(base, expo);

        ////////
        // limit the fugacity coefficient to a reasonable range:
        //
        // on one side, we want the mole fraction to be at
        // least 10^-3 if the fugacity is at the current pressure
        //
        fugCoeff = std::min(Scalar(1e10), fugCoeff);
        //
        // on the other hand, if the mole fraction of the component is 100%, we want the
        // fugacity to be at least 10^-3 Pa
        //
        fugCoeff = std::max(Scalar(1e-10), fugCoeff);
        ///////////

        return fugCoeff;
    }

};

template <class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::R = Opm::Constants<Scalar>::R;
template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::u = 2.0;
template<class Scalar, class StaticParameters>
const Scalar PengRobinsonMixture<Scalar, StaticParameters>::w = -1.0;

} // namespace Opm

#endif
