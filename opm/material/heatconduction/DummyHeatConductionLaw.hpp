// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 * \copydoc Opm::DummyHeatConductionLaw
 */
#ifndef OPM_DUMMY_HEATCONDUCTION_LAW_HH
#define OPM_DUMMY_HEATCONDUCTION_LAW_HH

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{
/*!
 * \ingroup material
 *
 * \brief Implements a dummy law for heat conduction to which isothermal models
 *        can fall back to.
 *
 * If any method of this law is called, it throws std::logic_error.
 */
template <class ScalarT>
class DummyHeatConductionLaw
{
public:
    typedef int Params;
    typedef ScalarT Scalar;

    /*!
     * \brief Given a fluid state, return the effective heat conductivity [W/m^2 / (K/m)] of the porous
     *        medium.
     *
     * If this method is called an exception is thrown at run time.
     */
    template <class FluidState>
    static Scalar heatConductivity(const Params &params,
                                   const FluidState &fluidState)
    {
        OPM_THROW(std::logic_error,
                   "No heat conduction law specified!");
    }
};
}

#endif
