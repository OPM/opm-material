// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
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
 * \copydoc Opm::TwoPAdapter
 */
#ifndef OPM_MP_2P_ADAPTER_HH
#define OPM_MP_2P_ADAPTER_HH

#include <algorithm>

namespace Opm {
/*!
 * \ingroup fluidmatrixinteractionslaws
 * \brief Provides an adapter class to use two-phase material laws
 *        with the generic M-phase API.
 */
template <int wPhaseIdxT, class TwoPLaw >
class TwoPAdapter
{
public:
    typedef typename TwoPLaw::Params Params;
    enum { numPhases = 2 };
    enum { wPhaseIdx = wPhaseIdxT };
    enum { nPhaseIdx = (wPhaseIdx == 0)?1:0 };

    /*!
     * \brief The capillary pressure-saturation curve.
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &fluidState)
    {
        // wetting phase does not get anything added
        values[wPhaseIdx] = 0;

        // non-wetting phase gets the capillary pressure added
        values[nPhaseIdx] = TwoPLaw::pC(params, fluidState.saturation(wPhaseIdx));
   }

    /*!
     * \brief The inverse of the capillary pressure-saturation curve.
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &fluidState)
    {
        // wetting phase does not get anything added
        values[wPhaseIdx] = TwoPLaw::Sw(params,
                                        fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx));
        values[nPhaseIdx] = 1.0 - values[wPhaseIdx];
    }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &fluidState)
    {
        values[wPhaseIdx] = TwoPLaw::krw(params, fluidState.saturation(wPhaseIdx));
        values[nPhaseIdx] = TwoPLaw::krn(params, fluidState.saturation(wPhaseIdx));
    }
};
}

#endif
